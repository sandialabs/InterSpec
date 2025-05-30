/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov.

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "InterSpec_config.h"

#include <memory>
#include <string>
#include <vector>
#include <cassert>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <functional>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <Wt/WText>
#include <Wt/WServer>
#include <Wt/WIOService>
#include <Wt/WApplication>
#include <Wt/WContainerWidget>


#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/BatchGuiInputFile.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"

#include "InterSpec/ReferenceLineInfo.h"

#if ( USE_REL_ACT_TOOL )
#include "InterSpec/RelActCalc.h"
#endif

#include "rapidxml/rapidxml.hpp"

using namespace Wt;
using namespace std;



BatchGuiInputFile::BatchGuiInputFile( const std::string display_name,
                                      const std::string path_to_file,
                                      const bool should_cleanup,
                                      Wt::WContainerWidget *parent )
: Wt::WContainerWidget( parent ),
  m_filename( path_to_file ),
  m_display_name( display_name ),
  m_should_cleanup( should_cleanup ),
  m_txt( nullptr ),
  m_remove_self_request_signal( this )
{
  addStyleClass( "BatchGuiInputFile" );

  InterSpec * const interspec = InterSpec::instance();
  interspec->useMessageResourceBundle( "BatchGuiInputFile" );

  wApp->useStyleSheet( "InterSpec_resources/BatchGuiInputFile.css" );
  
  m_txt = new WText( display_name, this );
  m_txt->addStyleClass( "FilenameText" );
  m_txt->setToolTip( "Full path of disk file: '" + path_to_file + "'" );

  WContainerWidget *close_icon = new WContainerWidget( this );
  close_icon->addStyleClass( "closeicon-wtdefault" );
  close_icon->clicked().connect( boost::bind( &BatchGuiInputFile::requestRemoveSelf, this ) );
  close_icon->clicked().preventPropagation();
}// BatchGuiInputFile constructor

BatchGuiInputFile::~BatchGuiInputFile()
{
  if( m_should_cleanup )
  {
    const bool success = SpecUtils::remove_file( m_filename );
    if( !success )
      cerr << "BatchGuiDialog::BatchGuiInputFile: Warning, could not delete file '" << m_filename << "'" << endl;
  }
}//~BatchGuiInputFile()

void BatchGuiInputFile::requestRemoveSelf()
{
  m_remove_self_request_signal.emit( this );
}// void requestRemoveSelf()

const std::string &BatchGuiInputFile::display_name() const
{
  return m_display_name;
}

const std::string &BatchGuiInputFile::path_to_file() const
{
  return m_filename;
}

Wt::Signal<BatchGuiInputFile *> &BatchGuiInputFile::remove_self_request()
{
  return m_remove_self_request_signal;
}


BatchGuiInputSpectrumFile::BatchGuiInputSpectrumFile( const std::string display_name,
                                                      const std::string path_to_file,
                                                      const bool should_cleanup,
                                                      Wt::WContainerWidget *parent )
: Wt::WContainerWidget( parent ),
  m_filename( path_to_file ),
  m_display_name( display_name ),
  m_should_cleanup( should_cleanup ),
  m_preview_container( nullptr ),
  m_preview_created( false ),
  m_spec_meas( nullptr ),
  m_is_peaks_csv( false ),
  m_preview_created_signal( this ),
  m_remove_self_request_signal( this )
{
  addStyleClass( "BatchGuiInputSpectrumFile" );

  InterSpec * const interspec = InterSpec::instance();
  interspec->useMessageResourceBundle( "BatchGuiInputFile" );
  
  wApp->useStyleSheet( "InterSpec_resources/BatchGuiInputFile.css" );
  
  m_preview_container = new Wt::WContainerWidget( this );
  m_preview_container->addStyleClass( "InputFilePreview" );

  WText *filename_text = new WText( display_name, this );
  filename_text->addStyleClass( "FilenameText" );
  filename_text->setToolTip( "Full path of disk file: '" + path_to_file + "'" );

  WContainerWidget *close_icon = new WContainerWidget( this );
  close_icon->addStyleClass( "closeicon-wtdefault" );
  close_icon->clicked().connect( boost::bind( &BatchGuiInputSpectrumFile::requestRemoveSelf, this ) );
  close_icon->clicked().preventPropagation();

  const string sessionid = wApp->sessionId();
  std::shared_ptr<int> status_ptr = make_shared<int>( 0 );
  std::shared_ptr<SpecMeas> spec_meas = make_shared<SpecMeas>();
  boost::function<void( void )> updateGuiCallback =
    wApp->bind( boost::bind( &BatchGuiInputSpectrumFile::set_spectrum, this, spec_meas, status_ptr ) );

  WServer::instance()->ioService().boost::asio::io_service::post(
    [updateGuiCallback, spec_meas, status_ptr, display_name, path_to_file, sessionid]()
    {
      const bool success = spec_meas->load_file( path_to_file, SpecUtils::ParserType::Auto, display_name );
      if( success )
      {
        *status_ptr = 1;
        spec_meas->set_filename( SpecUtils::filename( display_name ) );
      } else
      {
        // If a peak CSV file, set the status to 2, otherwise 3.
        bool is_csv = false;

        // We will just check the first line of the file to see if it is a peaks CSV file.
        //  The full parsing function, `PeakModel::csv_to_candidate_fit_peaks(...)` requires
        //  a SpecUtils::Measurement, that we dont have here.
        ifstream in_file( path_to_file );
        if( in_file )
        {
          string line;
          while( line.empty() && SpecUtils::safe_get_line( in_file, line, 8 * 1024 ) )
          {
          }

          is_csv = SpecUtils::icontains( line, "Centroid" ) && SpecUtils::icontains( line, "Net_Area" ) &&
                   SpecUtils::icontains( line, "FWHM" );
        }// if( in_file )

        *status_ptr = is_csv ? 2 : 3;
      }// if( parsed as spectrum file ) / else

      WServer::instance()->post( sessionid, updateGuiCallback );
    } );
}// BatchGuiInputSpectrumFile constructor


void BatchGuiInputSpectrumFile::set_spectrum( std::shared_ptr<SpecMeas> spec_meas, std::shared_ptr<int> status_ptr )
{
  if( status_ptr && ( ( *status_ptr ) == 2 ) )
  {
    m_is_peaks_csv = true;

    WText *preview = new WText( m_preview_container );
    preview->setText( WString::tr( "bgw-peaks-csv-txt" ) );
    preview->setStyleClass( "PeaksCsvFile" );
  } else if( !status_ptr || !spec_meas || ( spec_meas->num_measurements() == 0 ) || ( ( *status_ptr ) != 1 ) )
  {
    WText *preview = new WText( m_preview_container );
    preview->setText( WString::tr( "bgw-not-spec-preview" ) );
    preview->setStyleClass( "NotSpectrumFile" );
  } else
  {
    shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = nullptr;
    vector<shared_ptr<const ReferenceLineInfo>> reflines;

    const bool compact = true;
    const int width_px = 120, height_px = 70;
    //const shared_ptr<const ColorTheme> theme = InterSpec::instance()->getColorTheme();

    // Find most likely foreground.  If passthrough, sum everything
    shared_ptr<const SpecUtils::Measurement> preview_meas;

    const std::set<int> &samples = spec_meas->sample_numbers();
    const vector<string> &det_names = spec_meas->detector_names();
    if( spec_meas->passthrough() )
    {
      try
      {
        preview_meas = spec_meas->sum_measurements( samples, det_names, nullptr );
      } catch( std::exception & )
      {
        WText *preview = new WText( m_preview_container );
        preview->setText( WString::tr( "bgw-passthrough-sum-error" ) );
        return;
      }
    } else
    {
      // If not first passthrough, we'll take the first marked Foreground, Unknown,
      //  Background, sample number, in that order
      shared_ptr<const SpecUtils::Measurement> first_back, first_unk, first_fore;

      for( int sample_num : spec_meas->sample_numbers() )
      {
        shared_ptr<const SpecUtils::Measurement> sample_meas;
        try
        {
          if( det_names.size() == 1 )
            sample_meas = spec_meas->measurement( sample_num, det_names.front() );
          else
            sample_meas = spec_meas->sum_measurements( { sample_num }, det_names, nullptr );
        } catch( std::exception &e )
        {}// try / catch

        if( !sample_meas )
          continue;

        switch( sample_meas->source_type() )
        {
        case SpecUtils::SourceType::IntrinsicActivity:
        case SpecUtils::SourceType::Calibration:
          break;
        case SpecUtils::SourceType::Background:
          if( !first_back )
            first_back = sample_meas;
          break;
        case SpecUtils::SourceType::Foreground:
          if( !first_fore )
            first_fore = sample_meas;
          break;
        case SpecUtils::SourceType::Unknown:
          if( !first_unk )
            first_unk = sample_meas;
          break;
        }// switch( sample_meas->source_type() )

        if( first_fore )
          break;
      }// for( int sample_num : spec_meas->sample_numbers() )

      if( first_fore )
        preview_meas = first_fore;
      else if( first_unk )
        preview_meas = first_unk;
      else
        preview_meas = first_back;
    }// if( spec_meas->passthrough() ) / else

    if( preview_meas )
    {
      D3SpectrumDisplayDiv *spec = new D3SpectrumDisplayDiv( m_preview_container );
      spec->clicked().preventPropagation();
      spec->setThumbnailMode();
      spec->setData( preview_meas, false );
      spec->resize( WLength( 100, WLength::Percentage ), WLength( 100, WLength::Percentage ) );

      // We dont currently need to explicitly set the color theme, as all color theme styling for
      //   D3SpectrumDisplayDiv is globablly applied...
      // InterSpec *interspec = InterSpec::instance();
      // spec->applyColorTheme( interspec->getColorTheme() );
      // interspec->colorThemeChanged().connect( boost::bind( &D3SpectrumDisplayDiv::applyColorTheme, spec,
      // boost::placeholders::_1 ) );
    } else
    {
      WText *preview = new WText( m_preview_container );
      preview->setText( WString::tr( "bgw-preview-error" ) );
      preview->setStyleClass( "PreviewError" );
    }

    m_spec_meas = spec_meas;
    m_preview_created = true;
  }// if( !status_ptr || !spec_meas || (spec_meas->num_measurements() == 0) || ((*status_ptr) != 1) )

  m_preview_created_signal.emit( !!m_spec_meas );

  wApp->triggerUpdate();
}// void set_spectrum()

Wt::Signal<bool> &BatchGuiInputSpectrumFile::preview_created_signal()
{
  return m_preview_created_signal;
}

BatchGuiInputSpectrumFile::~BatchGuiInputSpectrumFile()
{
  if( m_should_cleanup )
  {
    const bool success = SpecUtils::remove_file( m_filename );
    if( !success )
      cerr << "BatchGuiDialog::BatchGuiInputSpectrumFile: Warning, could not delete file '" << m_filename << "'"
           << endl;
  }
}//~BatchGuiInputSpectrumFile()

std::shared_ptr<SpecMeas> BatchGuiInputSpectrumFile::spec_meas() const
{
  return m_spec_meas;
}

bool BatchGuiInputSpectrumFile::is_peaks_csv() const
{
  return m_is_peaks_csv;
}

void BatchGuiInputSpectrumFile::requestRemoveSelf()
{
  m_remove_self_request_signal.emit( this );
}// void requestRemoveSelf()

const std::string &BatchGuiInputSpectrumFile::display_name() const
{
  return m_display_name;
}

const std::string &BatchGuiInputSpectrumFile::path_to_file() const
{
  return m_filename;
}

Wt::Signal<BatchGuiInputSpectrumFile *> &BatchGuiInputSpectrumFile::remove_self_request()
{
  return m_remove_self_request_signal;
}
