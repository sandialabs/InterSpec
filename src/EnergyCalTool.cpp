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

#include <map>
#include <string>
#include <vector>
#include <ctime>
#include <memory>
#include <iostream>

#include <Wt/WMenu>
#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WAnchor>
#include <Wt/WResource>
#include <Wt/WCheckBox>
#include <Wt/WFileUpload>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WItemDelegate>
#include <Wt/WStackedWidget>
#include <Wt/WContainerWidget>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/EnergyCalTool.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/EnergyCalGraphical.h"
#include "InterSpec/EnergyCalUndoRedo.h"
#include "InterSpec/EnergyCalDevPairWidget.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/EnergyCalAddActions.h"
#include "InterSpec/IsotopeSelectionAids.h"


using namespace std;
using namespace Wt;

using EnergyCalImp::meas_old_new_cal_t;
using EnergyCalImp::meas_old_new_peaks_t;
using EnergyCalImp::EnergyCalUndoRedoSentry;

namespace
{
  // The EnergyCalUndoRedoSentry (and its do_undo_or_redo machinery) used to live here; it
  //  is now shared with EnergyCalMultiFile, in EnergyCalUndoRedo.{h,cpp}
  // The DevPair / DeviationPairDisplay widgets (and the index_compare_assend sort helper) used
  //  to live here too; they are now shared with EnergyCalMultiFile, in EnergyCalDevPairWidget.{h,cpp}
}//namespace


namespace EnergyCalImp
{

// DevPair and DeviationPairDisplay are defined in EnergyCalDevPairWidget.{h,cpp}

class CALpDownloadResource : public Wt::WResource
{
  Wt::WApplication *m_app;
  InterSpec *m_interspec;
  EnergyCalTool *m_tool;
  
public:
  CALpDownloadResource( EnergyCalTool *tool, InterSpec *viewer, WObject* parent = nullptr )
  : WResource( parent ), m_app( WApplication::instance() ), m_interspec( viewer ), m_tool( tool )
  {
    assert( m_app );
    assert( m_tool );
    assert( m_interspec );
  }
  
  virtual ~CALpDownloadResource()
  {
    beingDeleted();
  }
  
  virtual void handleRequest( const Wt::Http::Request &request, Wt::Http::Response &response )
  {
    assert( m_app );
    assert( m_interspec );
    
    try
    {
      WApplication::UpdateLock lock( m_app );
      
      if( !lock )
        throw std::runtime_error( "Error grabbing application lock to from CALpDownloadResource resource." );
  
      const SpecUtils::SpectrumType type = m_tool->typeOfCurrentlyShowingCoefficients();
      shared_ptr<SpecMeas> meas = m_interspec->measurment( type );
      if( !meas )
        throw std::runtime_error( "Error getting spectrum file currently being shown." );
      
      string filename = meas->filename();
      if( filename.empty() )
        filename = "energy_calibration";
      const string orig_extension = SpecUtils::file_extension(filename);
      if( orig_extension.size() && (orig_extension.size() < filename.size()) )
        filename = filename.substr(0,filename.size() - orig_extension.size());
      filename += ".CALp";
      
      //Remove bad filename characters
      const string notallowed = "\\/:?\"<>|*";
      for( auto it = begin(filename) ; it < end(filename) ; ++it )
      {
        if( notallowed.find(*it) != string::npos )
          *it = ' ';
      }
      
      suggestFileName( filename, WResource::Attachment );
      response.setMimeType( "application/octet-stream" );
      
      // First loop over visible measurements, adding calibrations for each new detector name,
      //  then loop over all measurements to pick up the rest.  This is because there may be multiple
      //  calibrations for a single detector, but we want to prefer the ones in the currently
      //  displayed spectra, and not deal with the complexity of handling multiple calibrations per
      //  detector
      set<string> dets_so_far;
      
      //Note that disp_detectors has both gamma and neutron detectors, but we only care about gamma
      const vector<string> &gamma_detectors = meas->gamma_detector_names();
      const vector<string> &detectors = meas->detector_names();
      const vector<string> disp_detectors = m_interspec->detectorsToDisplay(type);
      const vector<string> &neut_dets = meas->neutron_detector_names();
      
      for( const int sample : m_interspec->displayedSamples(type) )
      {
        for( const string det : disp_detectors )
        {
          if( dets_so_far.count(det) )
            continue;
          
          const shared_ptr<const SpecUtils::Measurement> m = meas->measurement( sample, det );
          shared_ptr<const SpecUtils::EnergyCalibration> cal = m ? m->energy_calibration() : nullptr;
          
          if( !cal || !cal->valid() || (cal->num_channels() < 3) )
          {
            // We'll assume that a detector that only has neutrons, will always only have neutrons...
            //  TODO: this may not actually be that good of an assumptions; re-evaluate later.
            if( std::find(begin(neut_dets), end(neut_dets), det) != end(neut_dets) )
              dets_so_far.insert( det );
            
            continue;
          }//if( energy cal is not valid )
          
          // Dont write the detector name if its unambiguous
          const string detname = (gamma_detectors.size() == 1) ? string() : det;
          
          if( SpecUtils::write_CALp_file(response.out(), cal, detname) )
          {
            dets_so_far.insert( det );
          }else
          {
            log("error") << "Error writing CALp file to WResource output stream.";
          }
        }//
        
        if( disp_detectors.size() == dets_so_far.size() )
          break;
      }//for( const int sample : m_interspec->displayedSamples(type) )
      
      if( disp_detectors.size() != dets_so_far.size() )
      {
        // Note that we are going over all samples and detectors here - being super thorough - we
        //  could probably tighten this up a lot.
        
        for( const int sample : meas->sample_numbers() )
        {
          for( const string &det : detectors )
          {
            if( dets_so_far.count(det) )
              continue;
            
            const shared_ptr<const SpecUtils::Measurement> m = meas->measurement( sample, det );
            shared_ptr<const SpecUtils::EnergyCalibration> cal = m ? m->energy_calibration() : nullptr;
            
            if( !cal || !cal->valid() || (cal->num_channels() < 3) )
            {
              continue;
            }//if( energy cal is not valid )
            
            if( SpecUtils::write_CALp_file(response.out(), cal, det) )
            {
              dets_so_far.insert( det );
            }else
            {
              log("error") << "Error writing CALp file to WResource output stream (2).";
            }
          }//for( const string &det : gamma_dets )
        }//for( const int sample : meas->sample_numbers() )
      }//if( disp_detectors.size() != dets_so_far.size() )
    }catch( std::exception &e )
    {
      log("error") << "Error handling request for CalFileDownloadResource: " << e.what();
      response.out() << "Error creating CALp file: " << e.what()
                     << "\n\nPlease report to InterSpec@sandia.gov.";
      
      //passMessage( "Error getting spectrum file currently being shown", WarningWidget::WarningMsgHigh );
      
      response.setStatus(500);
      assert( 0 );
    }//try / catch
  }//void handleRequest(...)
};//class CalFileDownloadResource


// DevPair and DeviationPairDisplay implementations moved to EnergyCalDevPairWidget.cpp


/** The label (and tooltip) for a coefficient row; for lower-channel-energy calibrations the two
 rows are the {offset, gain} adjustments relative to the sessions original channel energies (see
 #EnergyCal::adjust_lower_channel_energy_cal), so get their own labels/explanations.
 */
WString coef_disp_label_txt( const size_t order, const bool lower_channel )
{
  if( lower_channel )
    return (order == 0) ? WString::tr("ect-offset-label") : WString::tr("ect-lce-gain-label");

  switch( order )
  {
    case 0:  return WString::tr("ect-offset-label");
    case 1:  return WString::tr("ect-linear-label");
    case 2:  return WString::tr("ect-quad-label");
    case 3:  return WString::tr("ect-cubic-label");
  }//switch( order )

  return WString::tr("ect-nth-label").arg( static_cast<int>(order) );
}//coef_disp_label_txt(...)


WString coef_disp_tooltip_txt( const size_t order, const bool lower_channel )
{
  if( lower_channel )
    return (order == 0) ? WString::tr("ect-tt-lce-offset") : WString::tr("ect-tt-lce-gain");
  return WString();
}//coef_disp_tooltip_txt(...)


class CoefDisplay : public WContainerWidget
{
public:
  const size_t m_order;
  WLabel *m_label;
  WCheckBox *m_fit;
  NativeFloatSpinBox *m_value;

  CoefDisplay( const size_t order, WContainerWidget *parent = nullptr )
  : WContainerWidget( parent ),
    m_order( order ),
    m_label( nullptr ),
    m_fit( nullptr ),
    m_value( nullptr )
  {
    addStyleClass( "CoefDisplay" );

    m_label = new WLabel( coef_disp_label_txt(order, false), this );
    m_label->addStyleClass( "CoefLabel" );

    m_value = new NativeFloatSpinBox( this );
    m_value->setSpinnerHidden( true );

    m_fit = new WCheckBox( "Fit", this );
    m_fit->addStyleClass( "CoefFit CbNoLineBreak" );
  }//CoefDisplay

  /** Sets the rows label and tooltip, only updating the DOM if they actually changed (a display
   can flip between polynomial/FRF and lower-channel-energy labeling when the file changes).
   */
  void setLabeling( const bool lower_channel )
  {
    const WString txt = coef_disp_label_txt( m_order, lower_channel );
    if( m_label->text() != txt )
      m_label->setText( txt );

    const WString tooltip = coef_disp_tooltip_txt( m_order, lower_channel );
    if( m_label->toolTip() != tooltip )
      m_label->setToolTip( tooltip );
  }//setLabeling(...)
};//class CoefDisplay

class CalDisplay : public WContainerWidget
{
  static const size_t sm_min_coef_display ;
  
  EnergyCalTool *m_tool;
  const SpecUtils::SpectrumType m_cal_type;
  const std::string m_det_name;
  
  WText *m_type;
  WText *m_convertMsg;
  WContainerWidget *m_coefficients;
  DeviationPairDisplay *m_devPairs;
  
  //I cant decide if I like hiding deviation pair widget when there is none, or not
#define HIDE_EMPTY_DEV_PAIRS 0
#if( HIDE_EMPTY_DEV_PAIRS )
  WPushButton *m_addPairs;
#endif
  
#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
  WPushButton *m_fitCoeffs;
#endif
  
#if( IMP_CALp_BTN_NEAR_COEFS )
  WPushButton *m_downloadCALp;
  WPushButton *m_uploadCALp;
#endif
  
  shared_ptr<const SpecUtils::EnergyCalibration> m_cal;
  
public:
  CalDisplay( EnergyCalTool *tool,
             const SpecUtils::SpectrumType type,
             const std::string &detname,
             const bool isWideLayout,
             WContainerWidget *parent = nullptr )
  : WContainerWidget( parent ),
   m_tool( tool ),
   m_cal_type( type ),
   m_det_name( detname ),
   m_type( nullptr ),
   m_convertMsg( nullptr ),
   m_coefficients( nullptr ),
   m_devPairs( nullptr )
#if( HIDE_EMPTY_DEV_PAIRS )
   , m_addPairs( nullptr )
#endif
#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
  , m_fitCoeffs( nullptr )
#endif
#if( IMP_CALp_BTN_NEAR_COEFS )
  , m_downloadCALp( nullptr )
  , m_uploadCALp( nullptr )
#endif
  {
    addStyleClass( "CalDisplay" );
    
    WGridLayout *layout = new WGridLayout( this );
    layout->setContentsMargins( 0, 0, 0, 0 );
    layout->setVerticalSpacing( 0 );
    layout->setHorizontalSpacing( 0 );
    
    WContainerWidget *coefDiv = new WContainerWidget();
    coefDiv->addStyleClass( "CoefCol" );
    layout->addWidget( coefDiv, 0, 0 );
    
    m_type = new WText( "&nbsp;", coefDiv );
    m_type->setInline( false );
    m_type->addStyleClass( "Wt-itemview Wt-header Wt-label CalType" );
    
    m_coefficients = new WContainerWidget( coefDiv );
    m_coefficients->addStyleClass( "CoefContent" );
    
#if( IMP_COEF_FIT_BTN_NEAR_COEFS || IMP_CALp_BTN_NEAR_COEFS )
    WContainerWidget *btndiv = new WContainerWidget();
    btndiv->addStyleClass( "CalCoefsBtnDiv" );
  
#if( IMP_CALp_BTN_NEAR_COEFS )
    
    WResource *csv = m_model->peakCsvResource();
#if( BUILD_AS_OSX_APP || IOS )
    m_downloadCALp = new WAnchor( WLink(m_tool->calpResources()), btndiv );
    m_downloadCALp->setTarget( AnchorTarget::TargetNewWindow );
    m_downloadCALp->setStyleClass( "LinkBtn DownloadLink" );
#else
    m_downloadCALp = new WPushButton( btndiv );
    m_downloadCALp->setIcon( "InterSpec_resources/images/download_small.svg" );
    m_downloadCALp->setLink( WLink( m_tool->calpResources() ) );
    m_downloadCALp->setLinkTarget( Wt::TargetNewWindow );
    m_downloadCALp->setStyleClass( "LinkBtn DownloadBtn CALp" );
    
#endif //#if( BUILD_AS_OSX_APP || IOS ) / #else

    m_downloadCALp->setText( "CALp" );
    
    m_uploadCALp = new WPushButton( btndiv );
    m_uploadCALp->setIcon( "InterSpec_resources/images/upload_small.svg" );
    m_uploadCALp->setStyleClass( "LinkBtn UploadBtn CALp" );
    m_uploadCALp->clicked().connect( m_tool, &EnergyCalTool::handleRequestToUploadCALp );
#endif //#if( IMP_CALp_BTN_NEAR_COEFS )
    
    WContainerWidget *spacer = new WContainerWidget( btndiv );
    spacer->addStyleClass( "Spacer" );
    
#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
    m_fitCoeffs = new WPushButton( WString::tr("ect-fit-coeff-btn"), btndiv );
    m_fitCoeffs->addStyleClass( "CalCoefFitBtn" );
#endif
    
    layout->addWidget( btndiv, 1, 0 );
    layout->setRowStretch( 0, 1 );
#endif //#if( IMP_COEF_FIT_BTN_NEAR_COEFS || IMP_CALp_BTN_NEAR_COEFS )
    
    m_devPairs = new DeviationPairDisplay();
    
#if( HIDE_EMPTY_DEV_PAIRS )
    //For files with multiple detectors, the "Add dev. pairs" buttons doesnt show up right
    // for the detectors not currently showing - I guess should toggle dev pairs for all detectors.
    m_devPairs->setHidden( true );
    m_addPairs = new WPushButton( "Add dev. pairs" );
    m_addPairs->addStyleClass( "LinkBtn" );
    //m_addPairs->setIcon( "InterSpec_resources/images/plus_min_white.svg" );
    m_addPairs->setHidden( true );
    m_addPairs->clicked().connect( this, &CalDisplay::showDevPairs );
#endif
    
    if( isWideLayout )
    {
#if( IMP_COEF_FIT_BTN_NEAR_COEFS || IMP_CALp_BTN_NEAR_COEFS )
      layout->addWidget( m_devPairs, 0, 1, 2, 1 );
#else
      layout->addWidget( m_devPairs, 0, 1 );
#endif
      
#if( HIDE_EMPTY_DEV_PAIRS )
      layout->addWidget( m_addPairs, 1, 0, AlignmentFlag::AlignRight );
      layout->setRowStretch( 0, 1 );
#endif
    }else
    {
#if( HIDE_EMPTY_DEV_PAIRS )
      layout->addWidget( m_addPairs, layout->rowCount(), 0, AlignmentFlag::AlignCenter );
      layout->addWidget( m_devPairs, layout->rowCount(), 0 );
#else
      layout->addWidget( m_devPairs, layout->rowCount(), 0 );
#endif
      m_devPairs->setHeight( 100 );
    }
    
    m_devPairs->changed().connect( boost::bind( &EnergyCalTool::userChangedDeviationPair, m_tool, this,
                                               boost::placeholders::_1 ) );
  }//CalDisplay( constructor )
  
  SpecUtils::SpectrumType spectrumType() const { return m_cal_type; }
  const std::string &detectorName() const { return m_det_name; }
  
    
  //void setDeviationPairMsg( const std::string &msg )
  //{
  //  m_devPairs->setMsg( msg );
  //}//void setDeviationPairMsg( const std::string &msg )
  
  
  shared_ptr<const SpecUtils::EnergyCalibration> lastSetCalibration()
  {
    return m_cal;
  }
  
  /// @param fitfor The order coefficients that should be set checked to fit for
  void setFitFor( const set<size_t> &fitfor )
  {
    for( auto w : m_coefficients->children() )
    {
      auto ww = dynamic_cast<const CoefDisplay *>( w );
      assert( ww );
      if( ww )
        ww->m_fit->setChecked( fitfor.count(ww->m_order) );
    }//for( auto w : m_coefficients->children() )
  }//void setFitFor( const set<size_t> &fitfor )
  
#if( HIDE_EMPTY_DEV_PAIRS )
  void showDevPairs()
  {
    m_addPairs->setHidden( true );
    m_devPairs->setHidden( false, WAnimation(WAnimation::Fade, WAnimation::Linear, 200) );
  }
#endif
  
  
  /// @returns The order coefficents that are checked to be fit for
  set<size_t> fitForCoefficents() const
  {
    set<size_t> coeffs;
    
    for( auto w : m_coefficients->children() )
    {
      auto ww = dynamic_cast<const CoefDisplay *>( w );
      assert( ww );
      if( !ww )
        continue;
      
      if( ww->m_fit->isChecked() )
        coeffs.insert( ww->m_order );
    }//for( auto w : m_coefficients->children() )
    
    return coeffs;
  }//vector<bool> fitForCoefficents() const
  
  
  vector<float> displayedCoefficents()
  {
    vector<float> coeffs;
  
    const auto existing = m_coefficients->children();
    for( size_t i = 0; i < existing.size(); ++i )
    {
      if( !existing[i] )
        continue;  //shouldnt happen, but JIC
      
      auto ww = dynamic_cast<const CoefDisplay *>( existing[i] );
      assert( ww );
      if( !ww )
        continue;
    
      //NativeFloatSpinBox avoids round-off errors, and we will rely on this here.
      const float dispvalue = ww->m_value->value();
      coeffs.push_back( dispvalue );
    }//for( size_t i = 0; i < existing.size(); ++i )
    
    //remove trailing zeros
    while( !coeffs.empty() && coeffs.back()==0.0f )
      coeffs.resize( coeffs.size() - 1 );
    
    return coeffs;
  }//vector<float> displayedCoefficents() const
  
  
  std::vector< std::pair<float,float> > displayedDeviationPairs() const
  {
    return m_devPairs->deviationPairs();
  }
  
  void setDeviationPairsInvalid()
  {
    m_devPairs->setInvalidValues();
  }
  
  void setDeviationPairsValid()
  {
    m_devPairs->setValidValues();
  }
  
  void updateToGui( const shared_ptr<const SpecUtils::EnergyCalibration> &cal )
  {
#if( HIDE_EMPTY_DEV_PAIRS )
    const bool hadCal = !!m_cal;
#endif
    
    m_cal = cal;
    
#if( IMP_COEF_FIT_BTN_NEAR_COEFS || IMP_CALp_BTN_NEAR_COEFS )
    bool fitCoeffsVisible = false;
    const auto type = m_cal ? m_cal->type() : SpecUtils::EnergyCalType::InvalidEquationType;
    
    switch( type )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::FullRangeFraction:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      case SpecUtils::EnergyCalType::LowerChannelEdge:
        fitCoeffsVisible = true;
        break;

      case SpecUtils::EnergyCalType::InvalidEquationType:
        fitCoeffsVisible = false;
        break;
    }//switch( m_cal->type() )
    
#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
    if( m_fitCoeffs )
      m_fitCoeffs->setHidden( !fitCoeffsVisible );
#endif
    
#if( IMP_CALp_BTN_NEAR_COEFS )
    if( m_downloadCALp )
      m_downloadCALp->setHidden( (type == SpecUtils::EnergyCalType::InvalidEquationType) );
#endif
#endif //IMP_COEF_FIT_BTN_NEAR_COEFS || IMP_CALp_BTN_NEAR_COEFS
    
    
    if( !m_cal )
    {
      m_type->setText( "No Calibration" );
      m_coefficients->clear();
      m_devPairs->setDeviationPairs( {} );
#if( HIDE_EMPTY_DEV_PAIRS )
      m_devPairs->setHidden( true );
      m_addPairs->setHidden( true );
#endif
      return;
    }//if( !m_cal )
    
    WString typetxt;
    switch( m_cal->type() )
    {
      case SpecUtils::EnergyCalType::LowerChannelEdge:    typetxt = WString::tr("ect-cal-type-lower-channel"); break;
      case SpecUtils::EnergyCalType::InvalidEquationType: typetxt = WString::tr("ect-cal-type-not-defined");    break;
      case SpecUtils::EnergyCalType::Polynomial:          typetxt = WString::tr("ect-cal-type-polynomial");     break;
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
                                                          typetxt = WString::tr("ect-cal-type-default-poly");  break;
      case SpecUtils::EnergyCalType::FullRangeFraction:   typetxt = WString::tr("ect-cal-type-frf");            break;
    }//switch( m_cal->type() )
    
    m_type->setText( typetxt );
    
    switch( m_cal->type() )
    {
      case SpecUtils::EnergyCalType::InvalidEquationType:
        m_coefficients->clear();
        m_devPairs->setDeviationPairs( {} );
        if( !m_devPairs->isHidden() )
          m_devPairs->hide();

        if( !m_convertMsg )
        {
          if( auto p = dynamic_cast<WContainerWidget *>( m_type->parent() ) )
          {
            m_convertMsg = new WText( WString::tr("ect-convert-to-poly-msg"), p );
            m_convertMsg->addStyleClass( "ConvertToPolyMsg" );
            m_convertMsg->setInline( false );
          }
        }//if( !m_convertMsg )

#if( HIDE_EMPTY_DEV_PAIRS )
        m_addPairs->setHidden( true );
        m_devPairs->setHidden( true );
#endif
        return;

      case SpecUtils::EnergyCalType::LowerChannelEdge:
        // Lower channel energy calibrations dont have deviation pairs, but their {offset, gain}
        //  adjustments (relative to the sessions original channel energies) are editable
        if( !m_devPairs->isHidden() )
          m_devPairs->hide();
#if( HIDE_EMPTY_DEV_PAIRS )
        m_addPairs->setHidden( true );
#endif
        if( m_convertMsg )
        {
          delete m_convertMsg;
          m_convertMsg = nullptr;
        }
        break;

      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::FullRangeFraction:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
#if( !HIDE_EMPTY_DEV_PAIRS )
        if( m_devPairs->isHidden() )
          m_devPairs->show();
#endif
        if( m_convertMsg )
        {
          delete m_convertMsg;
          m_convertMsg = nullptr;
        }
        break;
    }//switch( m_cal->type() )

    const bool lower_channel = (m_cal->type() == SpecUtils::EnergyCalType::LowerChannelEdge);

    vector<pair<float,float>> devpairs;
    vector<float> coeffs;

    if( lower_channel )
    {
      // The two displayed "coefficients" are the cumulative {offset, gain} adjustments applied
      //  this session, relative to the original channel energies
      const LowerChanCalOriginal info = m_tool->lowerChannelOriginal( m_cal );
      coeffs = { static_cast<float>(info.offset), static_cast<float>(info.gain) };
    }else
    {
      devpairs = m_cal->deviation_pairs();
      coeffs = m_cal->coefficients();
    }
    
#if( HIDE_EMPTY_DEV_PAIRS )
    // Once deviation pairs are showing, we will leave them showing, even if there are no deviation
    //  pairs; this is because if you click to add deviation pairs, dont add any, then adjust the
    //  the gain, the deviation pair display getting hidden, causing a layout update, is really
    //  jarring.
    const auto anim = hadCal ? WAnimation(WAnimation::Fade, WAnimation::Linear, 200) : WAnimation{};
    if( m_devPairs->isHidden() && !devpairs.empty() )
      m_devPairs->setHidden( false, anim );
    m_addPairs->setHidden( !m_devPairs->isHidden(), anim );
#endif
    
    m_devPairs->setDeviationPairs( devpairs );
    //m_devPairs->changed().connect( std::bind( &EnergyCalTool::userChangedDeviationPair, m_tool, this) );
    
    const size_t num_coef_disp = lower_channel ? coeffs.size()
                                               : std::max( coeffs.size(), sm_min_coef_display );
    vector<CoefDisplay *> coef_disps( num_coef_disp, nullptr );

    size_t coefnum = 0;
    const auto existing = m_coefficients->children();
    for( size_t i = 0; i < existing.size(); ++i )
    {
      auto ww = dynamic_cast<CoefDisplay *>( existing[i] );
      assert( ww || !existing[i] );

      if( ww )
      {
        if( coefnum >= num_coef_disp )
        {
          m_coefficients->removeWidget( existing[i] );
          delete existing[i];
        }else
        {
          assert( coefnum < coef_disps.size() );
          coef_disps[coefnum] = ww;
          ww->setLabeling( lower_channel );
          const float value = (coefnum < coeffs.size()) ? coeffs[coefnum] : 0.0f;
          // Only write the value if it actually changed: NativeFloatSpinBox::setValue always
          //  re-formats and re-sets the DOM text, which would lose the cursor position, and
          //  clobber the formatting, of the value the user just typed in.
          if( ww->m_value->value() != value )
            ww->m_value->setValue( value );
        }

        ++coefnum;
      }//if( ww )
    }//for( size_t i = 0; i < existing.size(); ++i )

    for( ; coefnum < num_coef_disp; ++coefnum )
    {
      CoefDisplay *disp = new CoefDisplay( coefnum, m_coefficients );
      coef_disps[coefnum] = disp;
      disp->setLabeling( lower_channel );
      const float value = (coefnum < coeffs.size()) ? coeffs[coefnum] : 0.0f;
      disp->m_value->setValue( value );
      disp->m_fit->setChecked( (coefnum < 2) );
      
      disp->m_fit->changed().connect( m_tool, &EnergyCalTool::updateFitButtonStatus );
      
      /* Note: if the user uses the up.down arrows in a NativeFloatSpinBox to change values, things
               get all messed up (new values get set via c++ messing  up current values, or the
               valueChanged() callback gets called like 10 times per second, causing changes faster
               than everything can keep up, and just generally poor working), so for now I have
               disabled these spinners via #NativeFloatSpinBox::setSpinnerHidden()
       */
      disp->m_value->valueChanged().connect( boost::bind(&EnergyCalTool::userChangedCoefficient, m_tool, coefnum, this) );
    }
    
    
    //Set the step size to move the upper range of energy by about 1 keV per step
    // Set up the little tick/spin/whatever boxes
    /*
     //The NativeFloatSpinBox::setSpinnerHidden() call is currently removing the spin-box up/down
     //  arrow, so we wont set the step-size, as on Firefox if we fit for a value, then it will turn
     //  red if the new value doesnt hit on the step size
    for( size_t i = 0; i < coef_disps.size(); ++i )
    {
      CoefDisplay *disp = coef_disps[i];
      assert( disp );
      
      float stepsize = 1.0f;
      if( m_cal && m_cal->num_channels() > 4 )
      {
        const size_t nchannel = m_cal->num_channels();
        switch( m_cal->type() )
        {
          case SpecUtils::EnergyCalType::Polynomial:
          case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
            stepsize = 1.0f / std::pow(nchannel,i);
            break;
            
          case SpecUtils::EnergyCalType::FullRangeFraction:
            stepsize = 1.0;
            break;
            
          case SpecUtils::EnergyCalType::InvalidEquationType:
          case SpecUtils::EnergyCalType::LowerChannelEdge:
            stepsize = 0.0f;
            break;
        }//switch( m_coeffEquationType )
      }//if( valid calibration )
      
      disp->m_value->setSingleStep( stepsize );
    }//for( int i = 0; i < sm_numCoefs; ++i )
     */
    
  }//updateToGui(...)

#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
  void setFitButtonEnabled( const bool canFitCeofs )
  {
    if( !m_cal )
      return;
      
    if( m_fitCoeffs->isEnabled() != canFitCeofs )
      m_fitCoeffs->setEnabled( canFitCeofs );
  }
  
  Wt::EventSignal<Wt::WMouseEvent> &doFitCoeffs()
  {
    return m_fitCoeffs->clicked();
  }
#endif
  
};//class CalDisplay

const size_t CalDisplay::sm_min_coef_display = 4;


}//namespace


EnergyCalTool::EnergyCalTool( InterSpec *viewer, PeakModel *peakModel, WContainerWidget *parent )
: WContainerWidget( parent ),
  m_interspec( viewer ),
  m_peakModel( peakModel ),
  m_calpResource( new EnergyCalImp::CALpDownloadResource(this, viewer, this) ),
  m_tallLayoutContent( nullptr ),
  m_peakTable( nullptr ),
  m_specTypeMenu( nullptr ),
  m_specTypeMenuStack( nullptr ),
  m_detectorMenu{ nullptr },
  m_calInfoDisplayStack( nullptr ),
  m_noCalTxt( nullptr ),
  m_moreActionsColumn( nullptr ),
  m_applyToColumn( nullptr ),
  m_detColumn( nullptr ),
  m_detColLayout( nullptr ),
  m_calColumn( nullptr ),
  m_peakTableColumn( nullptr ),
  m_layout( nullptr ),
  m_applyToCbs{ nullptr },
  m_moreActions{ nullptr },
#if( !IMP_COEF_FIT_BTN_NEAR_COEFS )
  m_fitCalBtn( nullptr ),
#endif
#if( !IMP_CALp_BTN_NEAR_COEFS )
  m_downloadCALp( nullptr ),
  m_uploadCALp( nullptr ),
#endif
  m_lastGraphicalRecal( 0 ),
  m_lastGraphicalRecalType( EnergyCalGraphicalConfirm::NumRecalTypes ),
  m_lastGraphicalRecalEnergy( -999.0f ),
  m_graphicalRecal( nullptr ),
  m_addActionWindow( nullptr )
{
  wApp->useStyleSheet( "InterSpec_resources/EnergyCalTool.css" );
  
  assert( viewer );
  viewer->useMessageResourceBundle( "EnergyCalTool" );
    
  addStyleClass( "EnergyCalTool" );
  
  initWidgets( EnergyCalTool::LayoutType::Wide );
}


void EnergyCalTool::initWidgets( EnergyCalTool::LayoutType layoutType )
{
  const bool wide = (layoutType == LayoutType::Wide);
  
  if( (wide && !m_tallLayoutContent && m_layout) || (!wide && m_tallLayoutContent) )
    return;
  
  if( m_graphicalRecal )
  {
    AuxWindow::deleteAuxWindow( m_graphicalRecal );
    m_graphicalRecal = nullptr;
  }//if( m_graphicalRecal )
  
  
  
  if( wide )
  {
    removeStyleClass( "TallEnergyCal" );
    if( m_tallLayoutContent )
      delete m_tallLayoutContent;
    m_tallLayoutContent = nullptr;
    
    m_layout = new WGridLayout( this );
  }else
  {
    addStyleClass( "TallEnergyCal" );
    if( m_layout )
      delete m_layout;
    
    m_tallLayoutContent = new WContainerWidget( this );
    m_layout = new WGridLayout( m_tallLayoutContent );
  }//if( wide ) / else
  
  
  // \TODO: null out all the other memener variables
  m_peakTable = nullptr;
  m_specTypeMenu = nullptr;
  m_specTypeMenuStack = nullptr;
  for( auto &m : m_detectorMenu )
    m = nullptr;
  m_calInfoDisplayStack = nullptr;
  m_noCalTxt = nullptr;
  m_moreActionsColumn = nullptr;
  m_applyToColumn = nullptr;
  m_detColumn = nullptr;
  m_detColLayout = nullptr;
  m_calColumn = nullptr;
  m_peakTableColumn = nullptr;
  for( auto &m : m_applyToCbs )
    m = nullptr;
  for( auto &m : m_moreActions )
    m = nullptr;
#if( !IMP_COEF_FIT_BTN_NEAR_COEFS )
  m_fitCalBtn = nullptr;
#endif
#if( !IMP_CALp_BTN_NEAR_COEFS )
  m_downloadCALp = nullptr;
  m_uploadCALp = nullptr;
#endif
  
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", m_interspec );
  
  m_layout->setContentsMargins( 0, 0, 0, 0 );
  m_layout->setVerticalSpacing( 0 );
  m_layout->setHorizontalSpacing( 0 );
  
  m_noCalTxt = new WText( WString::tr("ect-no-spec") );
  m_noCalTxt->addStyleClass( "NoCalContentTxt" );
  if( wide )
    m_layout->addWidget( m_noCalTxt, 0, 0, AlignmentFlag::AlignCenter | AlignmentFlag::AlignMiddle );
  else
    m_layout->addWidget( m_noCalTxt, 0, 0, 1, 2, AlignmentFlag::AlignCenter | AlignmentFlag::AlignMiddle );
  
  //Create the more actions column...
  m_moreActionsColumn = new WContainerWidget();
  m_moreActionsColumn->addStyleClass( "ToolTabTitledColumn MoreActionCol" );
  if( wide )
    m_layout->addWidget( m_moreActionsColumn, 0, 1 );
  else
    m_layout->addWidget( m_moreActionsColumn, 3, 0 );
  
  WGridLayout *collayout = new WGridLayout( m_moreActionsColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );
  
  WText *header = new WText( WString::tr("ect-more-act") );
  header->addStyleClass( "ToolTabColumnTitle" );
  collayout->addWidget( header, 0, 0 );
  
  //We will put the apply-to list inside a div so we can style consistently with other rows
  // (a <ul> element doesnt accept same css as <div>, apparently).
  WContainerWidget *moreActionsDiv = new WContainerWidget();
  moreActionsDiv->addStyleClass( "ToolTabTitledColumnContent MoreActionsMenuContent" );
  collayout->addWidget( moreActionsDiv, 1, 0 );
  collayout->setRowStretch( 1, 1 );
  
  WContainerWidget *moreActionsList = new WContainerWidget( moreActionsDiv );
  moreActionsList->addStyleClass( "MoreActionsMenuList" );
  moreActionsList->setList( true );
  
  
  for( MoreActionsIndex index = static_cast<MoreActionsIndex>(0);
      index < MoreActionsIndex::NumMoreActionsIndex;
      index = MoreActionsIndex(static_cast<int>(index) + 1) )
  {
    const char *label = "", *tooltip = nullptr;
    switch( index )
    {
      case MoreActionsIndex::Linearize:
        label = "ect-linearize";
        tooltip = "ect-tt-linearize";
        break;
        
      case MoreActionsIndex::Truncate:
        label = "ect-truncate";
        tooltip = "ect-tt-truncate";
        break;
        
      case MoreActionsIndex::CombineChannels:
        label = "ect-combine";
        tooltip = "ect-tt-combine";
        break;
        
      case MoreActionsIndex::ConvertToFrf:
        label = "ect-to-frf";
        tooltip = "ect-tt-to-frf";
        break;
        
      case MoreActionsIndex::ConvertToPoly:
        label = "ect-to-poly";
        tooltip = "ect-tt-to-poly";
        break;
        
      case MoreActionsIndex::MultipleFilesCal:
        label = "ect-multi-file";
        tooltip = "ect-tt-multi-file";
        break;
        
      case MoreActionsIndex::NumMoreActionsIndex:
        assert(0);
        break;
    }//switch( index )
    
    WContainerWidget *holder = new WContainerWidget( moreActionsList );
    m_moreActions[static_cast<int>(index)] = new WAnchor( WLink(), WString::tr(label), holder );
    m_moreActions[static_cast<int>(index)]->clicked().connect( boost::bind(&EnergyCalTool::moreActionBtnClicked, this, index) );
    
    assert( tooltip );
    if( tooltip )
      HelpSystem::attachToolTipOn( holder, WString::tr(tooltip), showToolTips );
  }//for( loop over more actions )
  
  WContainerWidget *btndiv = new WContainerWidget();
  btndiv->addStyleClass( "BtmBtnDiv" );
  collayout->addWidget( btndiv, 2, 0 );
  
  auto helpBtn = new WContainerWidget( btndiv );
  helpBtn->addStyleClass( "Wt-icon ContentHelpBtn" );
  helpBtn->clicked().connect( boost::bind( &HelpSystem::createHelpWindow, "energy-calibration" ) );

  
#if( !IMP_CALp_BTN_NEAR_COEFS )
  m_uploadCALp = new WPushButton( btndiv );
  m_uploadCALp->setIcon( "InterSpec_resources/images/upload_small.svg" );
  m_uploadCALp->setStyleClass( "LinkBtn UploadBtn CALp" );
  m_uploadCALp->clicked().connect( this, &EnergyCalTool::handleRequestToUploadCALp );
  
#if( BUILD_AS_OSX_APP || IOS )
  m_downloadCALp = new WAnchor( WLink(m_calpResource), btndiv );
  m_downloadCALp->setTarget( AnchorTarget::TargetNewWindow );
  m_downloadCALp->setStyleClass( "LinkBtn DownloadLink CALp" );
#else
  m_downloadCALp = new WPushButton( btndiv );
  m_downloadCALp->setIcon( "InterSpec_resources/images/download_small.svg" );
  m_downloadCALp->setLink( WLink( m_calpResource ) );
  m_downloadCALp->setLinkTarget( Wt::TargetNewWindow );
  m_downloadCALp->setStyleClass( "LinkBtn DownloadBtn CALp" );
  
#endif
  m_downloadCALp->setText( "CALp" );
  
  m_downloadCALp->clicked().connect( std::bind([this](){
    m_interspec->logMessage( WString::tr("ect-export-CALp-msg"), WarningWidget::WarningMsgInfo );
  }) );
  HelpSystem::attachToolTipOn( m_downloadCALp, WString::tr("ect-tt-CALp"), showToolTips );
  
  m_downloadCALp->setHidden( true );
  m_uploadCALp->setHidden( true );
#endif // !IMP_CALp_BTN_NEAR_COEFS
  
  
#if( !IMP_COEF_FIT_BTN_NEAR_COEFS )
  m_fitCalBtn = new WPushButton( WString::tr("ect-fit-coeff-btn"), btndiv );
  m_fitCalBtn->addStyleClass( "FitCoefBtn" );
  m_fitCalBtn->clicked().connect( this, &EnergyCalTool::fitCoefficients );
  m_fitCalBtn->setDisabled( true );
  HelpSystem::attachToolTipOn( m_fitCalBtn, WString::tr("ect-tt-fit-coeff-btn"), showToolTips );
#endif // !IMP_COEF_FIT_BTN_NEAR_COEFS
  
  // Create the "Apply To" column that determines what to apply changes to
  m_applyToColumn = new WContainerWidget();
  m_applyToColumn->addStyleClass( "ToolTabTitledColumn ApplyToCol" );
  if( wide )
    m_layout->addWidget( m_applyToColumn, 0, 2 );
  else
    m_layout->addWidget( m_applyToColumn, 3, 1 );
  
  collayout = new WGridLayout( m_applyToColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );
  
  header = new WText( WString::tr("ect-apply-changes-to") );
  header->addStyleClass( "ToolTabColumnTitle" );
  collayout->addWidget( header, 0, 0 );
  
  //We will put the apply-to list inside a div so we can style consistently with other rows
  // (a <ul> element doesnt accept same css as <div>, apparently).
  WContainerWidget *applyToDiv = new WContainerWidget();
  applyToDiv->addStyleClass( "ToolTabTitledColumnContent ApplyToMenuContent" );
  collayout->addWidget( applyToDiv, 1, 0 );
  
  WContainerWidget *applyToList = new WContainerWidget( applyToDiv );
  applyToList->addStyleClass( "ApplyToMenuList" );
  applyToList->setList( true );
  
  
  for( ApplyToCbIndex index = static_cast<ApplyToCbIndex>(0);
      index < ApplyToCbIndex::NumApplyToCbIndex;
      index = ApplyToCbIndex(index+1) )
  {
    const char *label = "";
    switch( index )
    {
      case ApplyToCbIndex::ApplyToForeground:         label = "Foreground";       break;
      case ApplyToCbIndex::ApplyToBackground:         label = "Background";       break;
      case ApplyToCbIndex::ApplyToSecondary:          label = "Secondary";        break;
      case ApplyToCbIndex::ApplyToDisplayedDetectors: label = "ect-disp-dets";    break;
      case ApplyToCbIndex::ApplyToAllDetectors:       label = "ect-all-dets";     break;
      case ApplyToCbIndex::ApplyToDisplayedSamples:   label = "ect-disp-samples"; break;
      case ApplyToCbIndex::ApplyToAllSamples:         label = "ect-all-samples";  break;
      case ApplyToCbIndex::NumApplyToCbIndex:
        assert( 0 );
        break;
    }//switch( index )
    
    WContainerWidget *item = new WContainerWidget( applyToList );
    item->addStyleClass( "ApplyToItem" );
    auto cb = new WCheckBox( WString::tr(label), item );
    cb->setWordWrap( false );
    cb->addStyleClass( "ApplyToItem CbNoLineBreak" );
    
    
    switch( index )
    {
      case ApplyToCbIndex::ApplyToForeground:
      case ApplyToCbIndex::ApplyToBackground:
      case ApplyToCbIndex::ApplyToSecondary:
      case ApplyToCbIndex::ApplyToAllDetectors:
      case ApplyToCbIndex::ApplyToAllSamples:
        cb->setChecked( true );
        break;
        
      case ApplyToCbIndex::ApplyToDisplayedDetectors:
      case ApplyToCbIndex::ApplyToDisplayedSamples:
      case ApplyToCbIndex::NumApplyToCbIndex:
        cb->setChecked( false );
        break;
    }//switch( index )
    
    cb->checked().connect( boost::bind( &EnergyCalTool::applyToCbChanged, this, index ) );
    cb->unChecked().connect( boost::bind( &EnergyCalTool::applyToCbChanged, this, index ) );
    
    m_applyToCbs[index] = cb;
  }//for( loop over ApplyToCbIndex )
  
  WAnimation animation(Wt::WAnimation::Fade, Wt::WAnimation::Linear, 200);

  
  // Create the "Coefficients" column that show the polynomial/FRF coefficents.
  m_calColumn = new WContainerWidget();
  m_calColumn->addStyleClass( "ToolTabTitledColumn CoefColumn" );
  if( wide )
    m_layout->addWidget( m_calColumn, 0, 3 );
  else
    m_layout->addWidget( m_calColumn, 1, 0, 1, 2 );
  
  
  collayout = new WGridLayout( m_calColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );
  
  if( wide )
    collayout->setColumnStretch( 1, 1 );
  
  header = new WText( WString::tr("ect-calib-coeffs") );
  header->addStyleClass( "ToolTabColumnTitle" );
  
  collayout->addWidget( header, 0, 0, 1, 2 );
  //collayout->addWidget( m_calInfoDisplayStack, 1, 0 );
  
  
  // Create the "Detector" column that determines which coefficients to show
  m_detColumn = new WContainerWidget();
  m_detColumn->addStyleClass( "DetCol" );
  collayout->addWidget( m_detColumn, 1, 0 );
  
  m_detColLayout = new WGridLayout( m_detColumn );
  m_detColLayout->setContentsMargins( 0, 0, 0, 0 );
  m_detColLayout->setVerticalSpacing( 0 );
  m_detColLayout->setHorizontalSpacing( 0 );
  
  auto detheader = new WText( WString::tr("Detector") );
  detheader->setInline( false );
  detheader->addStyleClass( "DetHdr Wt-itemview Wt-header Wt-label" );
  //detheader->resize( WLength::Auto, WLength(20,WLength::Unit::Pixel) );
  //collayout->addWidget( detheader, 0, 0 );
  m_detColLayout->addWidget( detheader, 0, 0  );
  m_detColLayout->setRowStretch( 2, 1 );
  
  // Create the "Cal Peaks" table
  m_peakTableColumn = new WContainerWidget();
  m_peakTableColumn->addStyleClass( "ToolTabTitledColumn PeakTableCol" );
  if( wide )
  {
    m_layout->addWidget( m_peakTableColumn, 0, 4 );
    m_layout->setColumnStretch( 4, 1 );
  }else
  {
    m_layout->addWidget( m_peakTableColumn, 2, 0, 1, 2 );
  }

  collayout = new WGridLayout( m_peakTableColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  if( wide )
    collayout->setRowStretch( 1, 1 );
  
  header = new WText( WString::tr("ect-cal-peaks") );
  header->addStyleClass( "ToolTabColumnTitle" );
  collayout->addWidget( header, 0, 0 );
  
  m_peakTable = new RowStretchTreeView();
  m_peakTable->addStyleClass( "ToolTabTitledColumnContent PeakTable" );
  collayout->addWidget( m_peakTable, 1, 0 );
  collayout->setRowStretch( 1, 1 );
  
  m_peakTable->setRootIsDecorated( false ); //makes the tree look like a table! :)
  m_peakTable->setModel( m_peakModel );
  const int numModelCol = m_peakModel->columnCount();
  for( int col = 0; col < numModelCol; ++col )
    m_peakTable->setColumnHidden( col, true );
  
  m_peakTable->setSortingEnabled( true );
  m_peakTable->setAlternatingRowColors( true );
  m_peakTable->setSelectable( true );
  m_peakTable->setSelectionMode( SingleSelection );
  m_peakTable->setEditTriggers( WAbstractItemView::SingleClicked
                               | WAbstractItemView::DoubleClicked );
  
  m_peakTable->setColumnHidden( PeakModel::kUseForCalibration, false );
  m_peakTable->setColumnHidden( PeakModel::kMean, false );
  m_peakTable->setColumnHidden( PeakModel::kIsotope, false );
  m_peakTable->setColumnHidden( PeakModel::kPhotoPeakEnergy, false );
  m_peakTable->setColumnHidden( PeakModel::kDifference, false );
  
  
  m_peakTable->setColumnWidth( PeakModel::kUseForCalibration, WLength(3.7, WLength::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kMean, WLength(4.5, WLength::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kIsotope, WLength(4.5, WLength::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kPhotoPeakEnergy, WLength(6.25, WLength::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kDifference, WLength(5, WLength::FontEm) );
  
  
  
  WItemDelegate *dblDelagate = new WItemDelegate( m_peakTable );
  dblDelagate->setTextFormat( "%.2f" );
  m_peakTable->setItemDelegateForColumn( PeakModel::kMean, dblDelagate );
  
  PhotopeakDelegate *nuclideDelegate = new PhotopeakDelegate( PhotopeakDelegate::NuclideDelegate, true, m_peakTable );
  m_peakTable->setItemDelegateForColumn( PeakModel::kIsotope, nuclideDelegate );
  
  PhotopeakDelegate *photopeakDelegate = new PhotopeakDelegate( PhotopeakDelegate::GammaEnergyDelegate, true, m_peakTable );
  m_peakTable->setItemDelegateForColumn( PeakModel::kPhotoPeakEnergy, photopeakDelegate );
  
  m_peakModel->dataChanged().connect( this, &EnergyCalTool::updateFitButtonStatus );
  m_peakModel->rowsRemoved().connect( this, &EnergyCalTool::updateFitButtonStatus );
  m_peakModel->rowsInserted().connect( this, &EnergyCalTool::updateFitButtonStatus );
  m_peakModel->layoutChanged().connect( this, &EnergyCalTool::updateFitButtonStatus );
  
  m_interspec->displayedSpectrumChanged().connect(
              boost::bind( &EnergyCalTool::displayedSpectrumChanged,
                           this, boost::placeholders::_1, boost::placeholders::_2,
                          boost::placeholders::_3, boost::placeholders::_4 ) );
  
  m_renderFlags |= EnergyCalToolRenderFlags::FullGuiUpdate;
  scheduleRender();
}//void initWidgets( EnergyCalTool::LayoutType layout )


void EnergyCalTool::setWideLayout()
{
  initWidgets( EnergyCalTool::LayoutType::Wide );
}//void setWideLayout()


void EnergyCalTool::setTallLayout()
{
  initWidgets( EnergyCalTool::LayoutType::Tall );
}//void setTallLayout()


EnergyCalTool::~EnergyCalTool()
{
}
  

set<string> EnergyCalTool::gammaDetectorsForDisplayedSamples( const SpecUtils::SpectrumType type )
{
  //We want the names of just the detectors that have gamma calibration information, of the
  //  currently displayed samples.
  //  We will assume the first Measurement for a given named detector will have gamma data if
  //  any of the Measurements from that detector will.
  //  \TODO: evaluate if that is true.
  
  auto meas = m_interspec->measurment( type );
  
  if( !meas )
    return {};
  
  const vector<string> &detnames = meas->gamma_detector_names();
  const vector<string> displayedDets = m_interspec->detectorsToDisplay(type);
  const set<int> &samples = m_interspec->displayedSamples(type);
  
  set<string> detectors, nongammadets;
  
  for( const int sample : samples )
  {
    for( const string &name : detnames )
    {
      if( detectors.count(name) || nongammadets.count(name) )
        continue;
      
      auto m = meas->measurement( sample, name );
      if( m && (m->num_gamma_channels() > 4) )
        detectors.insert( name );
      else if( m && !detectors.count(name) )
        nongammadets.insert(name);
    }//for( const string &name : detnames )
    
    if( (detectors.size() + nongammadets.size()) == detnames.size() )
      break;
  }//for( const int sample : samples )
  
  return detectors;
}//set<string> gammaDetectorsForDisplayedSamples( const SpecUtils::SpectrumType type )


#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
vector<EnergyCalImp::CalDisplay *> EnergyCalTool::calDisplays()
{
  vector<EnergyCalImp::CalDisplay *> answer;
  for( int i = 0; i < 3; ++i )
  {
    WMenu *detMenu = m_detectorMenu[i];
    
    if( !detMenu )
      continue;
    
    WStackedWidget *stack = detMenu->contentsStack();
    
    if( !stack )
      continue;
    
    for( WWidget *w : stack->children() )
    {
      EnergyCalImp::CalDisplay *caldisp = dynamic_cast<EnergyCalImp::CalDisplay *>( w );
      assert( caldisp );
      if( caldisp )
        answer.push_back( caldisp );
    }
  }//for( int i = 0; i < 3; ++i )
  
  return answer;
}//vector<EnergyCalImp::CalDisplay *> calDisplays()
#endif





vector<MeasToApplyCoefChangeTo> EnergyCalTool::measurementsToApplyCoeffChangeTo()
{
  std::vector<MeasToApplyCoefChangeTo> answer;
  
  //Lets loop over spectrum types (Foreground, Background, Secondary), and decide if we should apply
  //  changes to that file, and if so, decide which sample numbers/detectors
  const SpecUtils::SpectrumType spectypes[3] = {
    SpecUtils::SpectrumType::Foreground,
    SpecUtils::SpectrumType::Background,
    SpecUtils::SpectrumType::SecondForeground
  };
  
  for( const auto spectype : spectypes )
  {
    auto meas = m_interspec->measurment( spectype );
    if( !meas )
      continue;
    
    switch( spectype )
    {
      case SpecUtils::SpectrumType::Foreground:
        if( !m_applyToCbs[ApplyToCbIndex::ApplyToForeground]->isChecked() )
          continue;
        break;
        
      case SpecUtils::SpectrumType::SecondForeground:
        if( !m_applyToCbs[ApplyToCbIndex::ApplyToSecondary]->isChecked() )
          continue;
        break;
        
      case SpecUtils::SpectrumType::Background:
        if( !m_applyToCbs[ApplyToCbIndex::ApplyToBackground]->isChecked() )
          continue;
        break;
    }//switch( spectype )
    
    //If we're here, we should apply changes to this spectrum type.
    
    //It could be that the background SpecFile is the same as the foreground, so check if we
    //  already have an entry for this file in answer, and if so, use it.  We dont want duplicate
    //  entries for a SpecFile, since consumers assume each measurement is visited only once (they
    //  would otherwise move peaks multiple times, or fail to look up already-updated calibrations).
    //  Note: samples displayed by a spectrum type that is NOT selected get subtracted back out
    //  below, after this loop; the remaining compromise is that *detector* scopes of multiple
    //  selected types showing the same file get unioned (handling that exactly would need
    //  per-(samples,detectors) entries, and consumers hardened against visiting a measurement
    //  twice, which isnt worth the complexity for how rare the situation is).
    MeasToApplyCoefChangeTo *changes = nullptr;
    for( size_t i = 0; !changes && i < answer.size(); ++i )
      changes = (answer[i].meas == meas) ? &(answer[i]) : changes;
    
    if( !changes )
    {
      //We havent seen this SpecFile yet, create an entry in answer for it
      answer.emplace_back();
      changes = &(answer.back()); //C++14 returns a reference to the emplaced object.
      changes->meas = meas;
    }
    
    assert( changes );
    
    //m_applyToCbs[ApplyToCbIndex::ApplyToDisplayedDetectors]->parent()->isHidden();
    const set<string> detectors = gammaDetectorsForDisplayedSamples(spectype);
    const vector<string> displayed_dets = m_interspec->detectorsToDisplay(spectype);
    
    bool displayingAllDets = true;
    for( const auto &det : detectors )
    {
      if( std::find(begin(displayed_dets), end(displayed_dets), det) == end(displayed_dets) )
        displayingAllDets = false;
    }
    
    if( displayingAllDets )
    {
      //We will insert detector names since different SpecType from the same file could have
      //  different detector names available.
      // \TODO: if InterSpec class is upgraded to select detector by SpecType, we will have to
      //        upgrade this part of the code.
      for( const auto &det : detectors )
        changes->detectors.insert( det );
    }else
    {
      const auto applyToAllCb = m_applyToCbs[ApplyToCbIndex::ApplyToAllDetectors];
      const auto applyToDisplayedCb = m_applyToCbs[ApplyToCbIndex::ApplyToDisplayedDetectors];
      
      const bool toAll = (applyToAllCb->parent()
                          && !applyToAllCb->parent()->isHidden()
                          && applyToAllCb->isChecked());
      const bool toDisplayed = (applyToDisplayedCb->parent()
                                && !applyToDisplayedCb->parent()->isHidden()
                                && applyToDisplayedCb->isChecked());

      // The two checkboxes are kept mutually exclusive by #applyToCbChanged, but both can read
      //  un-checked while their row is hidden - in that case (or any inconsistency) we fall back
      //  to applying to all detectors.
      assert( (toAll != toDisplayed)
              || applyToAllCb->parent()->isHidden() || applyToDisplayedCb->parent()->isHidden() );

      if( toAll || (toAll == toDisplayed) )
      {
        for( const auto &det : detectors )
          changes->detectors.insert( det );
      }else
      {
        for( const auto &dispdet : displayed_dets )
        {
          if( detectors.count(dispdet) )
            changes->detectors.insert( dispdet );
        }
      }//if( apply to all detectors ) / else ( only displayed detectors )
    }//if( displayingAllDets ) / else
    
    
    const auto toAllSamplesCb = m_applyToCbs[ApplyToCbIndex::ApplyToAllSamples];
    const bool onlyDispSamples = (toAllSamplesCb->parent()
                                  && !toAllSamplesCb->parent()->isHidden()
                                  && !toAllSamplesCb->isChecked());
    
    if( onlyDispSamples )
    {
      const set<int> &displayed_samples = m_interspec->displayedSamples(spectype);
      for( const int sample : displayed_samples )
        changes->sample_numbers.insert( sample );
    }else
    {
      changes->sample_numbers = meas->sample_numbers();
    }
  }//for( const auto spectype : spectypes )


  // Dont let a change bleed into the display of an un-selected spectrum type, when the same file
  //  is displayed as multiple types.  E.g., with the foreground and background being different
  //  sample numbers of the same file, and only "Background" selected, the samples the foreground
  //  is displaying should be left alone (even when the sample scope is "all samples").
  //  Samples displayed by both a selected and an un-selected type are kept - the change is to the
  //  underlying measurements, so the un-selected display cant help but follow.
  const auto type_is_checked = [this]( const SpecUtils::SpectrumType spectype ) -> bool {
    switch( spectype )
    {
      case SpecUtils::SpectrumType::Foreground:
        return m_applyToCbs[ApplyToCbIndex::ApplyToForeground]->isChecked();
      case SpecUtils::SpectrumType::SecondForeground:
        return m_applyToCbs[ApplyToCbIndex::ApplyToSecondary]->isChecked();
      case SpecUtils::SpectrumType::Background:
        return m_applyToCbs[ApplyToCbIndex::ApplyToBackground]->isChecked();
    }//switch( spectype )

    assert( 0 );
    return false;
  };

  for( MeasToApplyCoefChangeTo &change : answer )
  {
    set<int> checked_disp, unchecked_disp;
    for( const auto spectype : spectypes )
    {
      if( m_interspec->measurment(spectype) != change.meas )
        continue;

      const set<int> &disp = m_interspec->displayedSamples(spectype);
      if( type_is_checked(spectype) )
        checked_disp.insert( begin(disp), end(disp) );
      else
        unchecked_disp.insert( begin(disp), end(disp) );
    }//for( const auto spectype : spectypes )

    for( const int sample : unchecked_disp )
    {
      if( !checked_disp.count(sample) )
        change.sample_numbers.erase( sample );
    }
  }//for( MeasToApplyCoefChangeTo &change : answer )


  return answer;
}//std::vector<MeasToApplyCoefChangeTo> measurementsToApplyCoeffChangeTo()


LowerChanCalOriginal EnergyCalTool::lowerChannelOriginal(
                        const std::shared_ptr<const SpecUtils::EnergyCalibration> &cal ) const
{
  assert( cal && cal->valid() && (cal->type() == SpecUtils::EnergyCalType::LowerChannelEdge) );

  const auto pos = m_lowerChanOrigCals.find( cal );
  if( pos != end(m_lowerChanOrigCals) )
  {
    assert( pos->second.original );
    return pos->second;
  }

  //This calibration hasnt been adjusted this session - it IS the original
  LowerChanCalOriginal answer;
  answer.original = cal;
  return answer;
}//lowerChannelOriginal(...)


void EnergyCalTool::registerLowerChannelAdjustment(
                        const std::shared_ptr<const SpecUtils::EnergyCalibration> &new_cal,
                        const LowerChanCalOriginal &info )
{
  assert( new_cal && info.original );
  assert( new_cal->type() == SpecUtils::EnergyCalType::LowerChannelEdge );
  assert( info.original->type() == SpecUtils::EnergyCalType::LowerChannelEdge );

  if( !new_cal || !info.original || (new_cal == info.original) )
    return;

  m_lowerChanOrigCals[new_cal] = info;
}//registerLowerChannelAdjustment(...)


void EnergyCalTool::applyCALpEnergyCal( std::map<std::string,std::shared_ptr<const SpecUtils::EnergyCalibration>> det_to_cal,
                                         const SpecUtils::SpectrumType specfile,
                                         const bool all_detectors, const bool all_samples )
{
  // TODO: add option to not set deviation pairs
  EnergyCalUndoRedoSentry undo_sentry;
  
  set<string> fore_gamma_dets;
  const shared_ptr<SpecMeas> measurment = m_interspec->measurment( specfile );
  const shared_ptr<const SpecUtils::Measurement> disp_spec = m_interspec->displayedHistogram( specfile );
  const shared_ptr<const SpecUtils::EnergyCalibration> old_disp_cal
                                           = disp_spec ? disp_spec->energy_calibration() : nullptr;
  const set<int> &disp_samples = m_interspec->displayedSamples( specfile );
  const vector<string> disp_detectors = m_interspec->detectorsToDisplay( specfile );
  
  const char * const desc = SpecUtils::descriptionText(specfile);
  
  if( !measurment || !disp_spec || !old_disp_cal || !old_disp_cal->valid() )
    throw runtime_error( WString::tr("ect-CALp-no-meas").arg(desc).toUTF8() );
  
  MeasToApplyCoefChangeTo tochange;
  tochange.meas = measurment;
  if( all_samples )
    tochange.sample_numbers = measurment->sample_numbers();
  else
    tochange.sample_numbers = disp_samples;
  
  const vector<string> &gamma_det_names = measurment->gamma_detector_names();
  if( all_detectors )
  {
    tochange.detectors.insert( begin(gamma_det_names), end(gamma_det_names) );
  }else
  {
    for( const string &det : disp_detectors )
    {
      const auto pos = std::find( begin(gamma_det_names), end(gamma_det_names), det );
      if( pos != end(gamma_det_names) )
        tochange.detectors.insert( det );
    }
  }//if( all_detectors ) / else
  
  
  if( tochange.detectors.empty() )
    throw runtime_error( WString::tr("ect-CALp-not-applic-det").arg(desc).toUTF8() );
  
  if( tochange.sample_numbers.empty() )
    throw runtime_error( WString::tr("ect-CALp-not-applic-sample").arg(desc).toUTF8() );
  
  
  if( det_to_cal.size() == 1 )
  {
    const shared_ptr<const SpecUtils::EnergyCalibration> new_cal = det_to_cal.begin()->second;
    
    if( !new_cal || !new_cal->valid() )
      throw runtime_error( WString::tr("ect-CALp-invalid-input").toUTF8() );
    
    // Its possible things could work out if there is a different number of channels, but this
    //  doesnt make much sense, so we'll require it, at least for now.
    if( new_cal->num_channels() != old_disp_cal->num_channels() )
    {
      throw runtime_error( WString::tr("ect-CALp-num-channel-mismatch")
                          .arg( static_cast<int>(new_cal->num_channels()) )
                          .arg( static_cast<int>(old_disp_cal->num_channels()) )
                          .toUTF8() );
    }
    
    if( tochange.detectors.size() == 1 )
    {
      setEnergyCal( new_cal, tochange, true );
    }else
    {
      // Check to see if all detectors are using the same energy calibration, and if so, set the
      //  calibration, otherwise do the applyCalChange(...)
      set<shared_ptr<const SpecUtils::EnergyCalibration>> old_energy_cals;
      for( const string &det : tochange.detectors )
      {
        for( const int sample : tochange.sample_numbers )
        {
          const auto m = tochange.meas->measurement( sample, det );
          const shared_ptr<const SpecUtils::EnergyCalibration> cal = m ? m->energy_calibration() : nullptr;
          if( cal && cal->valid() && (cal->num_channels() >= 3) )
            old_energy_cals.insert( m->energy_calibration() );
        }//for( loop over samples )
      }//for( loop over detector names )
      
      if( old_energy_cals.size() > 1 )
      {
        // A single-calibration CALp is being applied to data whose detectors currently use
        //  multiple calibrations: we propagate the *coefficient* change to each detector
        //  (relative to the displayed calibration), and deliberately DO NOT apply the CALp
        //  deviation pairs (there is no consistent way to combine them with each detectors own
        //  pairs).  The workaround, if the deviation pairs are wanted, is to display a single
        //  detector at a time and apply the CALp to each.  We let the user know the deviation
        //  pairs got dropped, so the data loss isnt silent.
        if( !new_cal->deviation_pairs().empty() )
          m_interspec->logMessage( WString::tr("ect-calp-dev-pairs-dropped"), 2 );

        applyCalChange( old_disp_cal, new_cal, { tochange }, false );
      }else
      {
        setEnergyCal( new_cal, tochange, true );
      }
    }//if( we are changing a single detector ) / else
  }else
  {
    // We will apply the input calibration on a detector-by-detector basis, requiring the
    //  input to have info for every relevant detector.
    //
    //  We will only adjust peak positions if the detectors energy calibration is the
    //   `suggested_sum_energy_calibration`, and even then, only once.
    //   (note: not trusting disp_spec->energy_calibration() to be same as
    //          `suggested_sum_energy_calibration` - but I think we could).
    bool adjusted_peaks = false;
    shared_ptr<const SpecUtils::EnergyCalibration> display_cal;
    try
    {
      display_cal = tochange.meas->suggested_sum_energy_calibration( disp_samples, disp_detectors );
    }catch( std::exception & )
    {
    }
    
    if( !display_cal )
      display_cal = old_disp_cal;
    
    string missing_dets, extra_dets;
    for( const auto &det : tochange.detectors )
    {
      const auto pos = det_to_cal.find(det);
      
      if( (pos == end(det_to_cal)) || !pos->second || !pos->second->valid() )
        missing_dets += (missing_dets.empty() ? "" : ", ") + det;
    }
    
    for( const auto &det_cal : det_to_cal )
    {
      if( !tochange.detectors.count(det_cal.first) )
        extra_dets += (extra_dets.empty() ? "" : ", ") + det_cal.first;
    }
    
    // If we have extra detectors, no problem - if we dont have calibration for some detectors,
    //  its a bigger problem.
    if( missing_dets.size() )
    {
      WString msg = WString::tr("ect-CALp-missing-det").arg( missing_dets );
      if( extra_dets.size() )
        msg.arg( WString::tr("ect-CALp-extra-det").arg( extra_dets ) );
      else
        msg.arg( "" );
      
      throw runtime_error( msg.toUTF8() );
    }//if( missing_dets.size() )
    
    const set<string> detectors_to_apply = tochange.detectors;
    for( const string &det : detectors_to_apply )
    {
      const auto pos = det_to_cal.find(det);
      
      // TODO: we can probably do a better job of _trying_ something...
      if( pos == end(det_to_cal) )
      {
        if( det.empty() )
          throw runtime_error( WString::tr("ect-CALp-named-vs-unamed").toUTF8() );
        throw runtime_error( WString::tr("ect-CALp-no-cal-for-det").arg(det).toUTF8() );
      }//if( pos == end(det_to_cal) )
        
      const shared_ptr<const SpecUtils::EnergyCalibration> new_cal = pos->second;
      assert( new_cal && new_cal->valid() );
      if( !new_cal || !new_cal->valid() )
        throw runtime_error( "unexpected logic error" );
      
      shared_ptr<const SpecUtils::EnergyCalibration> old_cal;
      for( const int sample : tochange.sample_numbers )
      {
        const auto m = measurment->measurement( sample, det );
        const auto cal = m ? m->energy_calibration() : nullptr;
        if( cal && cal->valid() )
        {
          old_cal = cal;
          break;
        }
      }
      
      if( !old_cal )
        throw runtime_error( WString::tr("ect-CALp-no-prev-for-det").arg(det).toUTF8() );
      
      if( old_cal->num_channels() != new_cal->num_channels() )
        throw runtime_error( WString::tr("ect-CALp-nchannel-mismatch-det")
                               .arg( static_cast<int>(new_cal->num_channels()) )
                               .arg( static_cast<int>(old_cal->num_channels()) )
                               .arg(det).toUTF8() );
      
      MeasToApplyCoefChangeTo det_specific_tochange = tochange;
      
      det_specific_tochange.detectors.clear();
      det_specific_tochange.detectors.insert( det );
      
      // We will adjust the peak energies, only if the calibrations being changed contain the
      //  display energy calibration.
      const bool adjust_peaks = !adjusted_peaks && (display_cal == old_cal);
      adjusted_peaks |= adjust_peaks;
      
      setEnergyCal( new_cal, det_specific_tochange, adjust_peaks );
    }//for( const string &det : detectors_to_apply )
  }//if( det_to_cal.size() == 1 ) / else
  
  
  m_interspec->refreshDisplayedCharts();
  refreshGuiFromFiles();
}//void applyCALpEnergyCal(...)


SpecUtils::SpectrumType EnergyCalTool::typeOfCurrentlyShowingCoefficients() const
{
  // Prefer the currently showing CalDisplay, which knows its own spectrum type; this matters
  //  when every file has a single detector, and all the displays live in the foreground menu
  //  labeled "Foreground"/"Background"/"Secondary" (so the spectrum-type menu index alone would
  //  always say foreground).
  if( m_calInfoDisplayStack )
  {
    auto display = dynamic_cast<const EnergyCalImp::CalDisplay *>( m_calInfoDisplayStack->currentWidget() );
    if( display )
      return display->spectrumType();
  }//if( m_calInfoDisplayStack )

  assert( m_specTypeMenu );

  const SpecUtils::SpectrumType spectypes[3] = {
    SpecUtils::SpectrumType::Foreground,
    SpecUtils::SpectrumType::Background,
    SpecUtils::SpectrumType::SecondForeground
  };

  const int selectedType = m_specTypeMenu->currentIndex();
  if( selectedType < 0 || selectedType > 2 )
  {
    assert( 0 );  //shouldnt ever happen, right?
    throw runtime_error( "EnergyCalTool::typeOfCurrentlyShowingCoefficients(): invalid spec type" );
  }

  return spectypes[selectedType];
}//SpecUtils::SpectrumType typeOfCurrentlyShowingCoefficients() const

/*
 //Commented out because it is unused/untested
std::string EnergyCalTool::detectorNameOfCurrentlyShowingCoefficients() const
{
  assert( m_specTypeMenu );
  
  const SpecUtils::SpectrumType spectypes[3] = {
    SpecUtils::SpectrumType::Foreground,
    SpecUtils::SpectrumType::Background,
    SpecUtils::SpectrumType::SecondForeground
  };
  
  const int selectedType = m_specTypeMenu->currentIndex();
  if( selectedType < 0 || selectedType > 2 )
  {
    assert( 0 );  //shouldnt ever happen, right?
    return;
  }
  
  WMenu *detMenu = m_detectorMenu[selectedType];
  assert( detMenu );
  
  WMenuItem *detitem = detMenu->currentItem();
  if( !detitem || !detMenu->count() )
  {
    assert( 0 );
    return "";
  }
  
  return detitem->text().toUTF8();
}//std::string detectorNameOfCurrentlyShowingCoefficients() const
*/


EnergyCalImp::CALpDownloadResource *EnergyCalTool::calpResources()
{
  return m_calpResource;
}



void EnergyCalTool::handleRequestToUploadCALp()
{
  SimpleDialog *dialog = new SimpleDialog();
  WPushButton *closeButton = dialog->addButton( "Cancel" );
  WGridLayout *stretcher = new WGridLayout();
  stretcher->setContentsMargins( 0, 0, 0, 0 );
  dialog->contents()->setLayout( stretcher );
  dialog->contents()->setOverflow( WContainerWidget::Overflow::OverflowVisible,
                                  Wt::Horizontal | Wt::Vertical );
  WText *title = new WText( WString::tr("ect-import-CALp") );
  title->addStyleClass( "title" );
  stretcher->addWidget( title, 0, 0 );
  
  WText *t = new WText( WString::tr("ect-select-CALp") );
  stretcher->addWidget( t, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
  t->setTextAlignment( Wt::AlignCenter );
  
  
  WFileUpload *upload = new WFileUpload();
  upload->fileTooLarge().connect( std::bind( [=](){
    dialog->contents()->clear();
    dialog->footer()->clear();
    
    WPushButton *closeButton = dialog->addButton( WString::tr("Close") );
    WGridLayout *stretcher = new WGridLayout();
    stretcher->setContentsMargins( 0, 0, 0, 0 );
    dialog->contents()->setLayout( stretcher );
    WText *title = new WText( WString::tr("ect-upload-CALp-to-large") );
    title->addStyleClass( "title" );
    stretcher->addWidget( title, 0, 0 );
  }) );
  
  upload->changed().connect( upload, &WFileUpload::upload );
  upload->uploaded().connect( std::bind( [dialog,upload](){
    InterSpec *interspec = InterSpec::instance();
    SpecMeasManager *measmn = interspec ? interspec->fileManager() : nullptr;
    
    assert( measmn );
    if( !measmn )
      return;
    
    const string calp_spool_path = upload->spoolFileName();
    const size_t calp_file_size = SpecUtils::file_size( calp_spool_path );
    if( calp_file_size > 10 * 1024 * 1024 )
    {
      dialog->contents()->clear();
      dialog->footer()->clear();

      WPushButton *closeButton = dialog->addButton( WString::tr("Close") );
      WGridLayout *stretcher = new WGridLayout();
      stretcher->setContentsMargins( 0, 0, 0, 0 );
      dialog->contents()->setLayout( stretcher );
      WText *title = new WText( WString::tr("ect-upload-CALp-to-large") );
      title->addStyleClass( "title" );
      stretcher->addWidget( title, 0, 0 );

      return;
    }

#ifdef _WIN32
    const std::wstring wcalp_path = SpecUtils::convert_from_utf8_to_utf16( calp_spool_path );
    ifstream input( wcalp_path.c_str(), ios::in | ios::binary );
#else
    ifstream input( calp_spool_path.c_str(), ios::in | ios::binary );
#endif

    if( !measmn->handleCALpFile( input, dialog, true ) )
    {
      dialog->contents()->clear();
      dialog->footer()->clear();
      
      WPushButton *closeButton = dialog->addButton( WString::tr("Close") );
      WGridLayout *stretcher = new WGridLayout();
      stretcher->setContentsMargins( 0, 0, 0, 0 );
      dialog->contents()->setLayout( stretcher );
      WText *title = new WText( WString::tr("ect-invalid-CALp") );
      title->addStyleClass( "title" );
      stretcher->addWidget( title, 0, 0 );
      
      return;
    }//if( was not a valid CALp file )
    
    //wApp->doJavaScript( "$('.Wt-dialogcover').hide();" ); // JIC
    //dialog->done( Wt::WDialog::DialogCode::Accepted );
  } ) );
  
  stretcher->addWidget( upload, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
  
  InterSpec *interspec = InterSpec::instance();
  if( interspec && !interspec->isPhone() )
  {
    t = new WText( WString::tr("ect-CALp-drag-n-drop-note") );
    stretcher->addWidget( t, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
    t->setTextAlignment( Wt::AlignCenter );
  }
  
  /*
   //In case we want to use AuxWindow instead of SimpleDialog
   AuxWindow *window = new AuxWindow( "Import CALp file",
   (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
   | AuxWindowProperties::PhoneNotFullScreen
   | AuxWindowProperties::DisableCollapse
   | AuxWindowProperties::SetCloseable) );
   
   //...
   
   window->rejectWhenEscapePressed();
   window->show();
   window->resizeToFitOnScreen();
   window->centerWindow();
   
   WPushButton *close = window->addCloseButtonToFooter( "Cancel" );
   close->clicked().connect( boost::bind( &AuxWindow::hide, window ) );
   
   window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
   
   // TODO: add link to relevant section of documentation
   //AuxWindow::addHelpInFooter( window->footer(), "energy-cal-CALp" );
   */
}//void handleRequestToUploadCALp();



void EnergyCalTool::applyCalChange( std::shared_ptr<const SpecUtils::EnergyCalibration> disp_prev_cal,
                                    std::shared_ptr<const SpecUtils::EnergyCalibration> new_disp_cal,
                                    const vector<MeasToApplyCoefChangeTo> &changemeas,
                                    const bool isOffsetOnly )
{
  // We get here when:
  //  - we manually change a coefficiecnt on the GUI
  //  - we fit calibration coefficients
  //  - sometimes when we load in a CALp file
  
  using namespace SpecUtils;
  
  assert( disp_prev_cal && disp_prev_cal->valid() );
  assert( new_disp_cal && new_disp_cal->valid() );
  
  EnergyCalUndoRedoSentry undo_sentry;
  
  const auto forgrnd = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  const auto backgrnd = m_interspec->measurment(SpecUtils::SpectrumType::Background);
  const auto secgrnd = m_interspec->measurment(SpecUtils::SpectrumType::SecondForeground);
  
  const set<int> &foresamples = m_interspec->displayedSamples( SpectrumType::Foreground );
  const set<int> &backSamples = m_interspec->displayedSamples( SpecUtils::SpectrumType::Background );
  const set<int> &secoSamples = m_interspec->displayedSamples( SpecUtils::SpectrumType::SecondForeground );
  
  
  // Create a cache of modified calibration both to save time/memory, but also keep it so previous
  //  samples that share a energy calibration will continue to do so (if possible based on what user
  //  wanted calibration applied to).  Also, we wont set any new calibrations until we know all
  //  updated calibrations and peaks are valid
  //Note: we could take this oppritunity to share calibration across SpecFile objects by not just
  //      comparing pointers, but also the actual EnergyCalibration object.  But for now we'll
  //      skip this to avoid trouble, and it isnt clear that it would actually be overall beneficial
  map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> old_to_new_cals;
  
  // We will store updated peaks and not set any of them until we know all the energy calibrations
  //  and peak shifts were successfully done.
  map<shared_ptr<deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_peaks;
  map<shared_ptr<const deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_hint_peaks;
  
  // const vector<MeasToApplyCoefChangeTo> changemeas = measurementsToApplyCoeffChangeTo();
  
  //We will loop over the changes to apply twice.  Once to calculate new calibrations, and make sure
  //  they are valid, then a second time to actually set them.  If a new calibration is invalid,
  //  an exception will be thrown so we will catch that.
  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    assert( change.meas );
    
    //string dbgmsg = "For '" + change.meas->filename() + "' will apply changes to Detectors: {";
    //for( auto iter = begin(change.detectors); iter != end(change.detectors); ++iter )
    //  dbgmsg += (iter==begin(change.detectors) ? "" : ",") + (*iter);
    //dbgmsg += "} and Samples: {";
    //for( auto iter = begin(change.sample_numbers); iter != end(change.sample_numbers); ++iter )
    //  dbgmsg += (iter==begin(change.sample_numbers) ? "" : ",") + std::to_string(*iter);
    //dbgmsg += "}";
    //cout << dbgmsg << endl;
    //wApp->log("app:debug") << dbgmsg;
    
    try
    {
      for( const int sample : change.sample_numbers )
      {
        for( const string &detname : change.detectors )
        {
          auto m = change.meas->measurement( sample, detname );
          if( !m || m->num_gamma_channels() <= 4 )
            continue;
          
          const auto meas_old_cal = m->energy_calibration();
          assert( meas_old_cal );
          
          if( !meas_old_cal || !meas_old_cal->valid() )
            continue;
          
          //If we have already computed the new calibration for a EnergyCalibration object, lets not
          //  re-due it.
          if( old_to_new_cals.count(meas_old_cal) )
            continue;
          
          shared_ptr<const SpecUtils::EnergyCalibration> new_meas_cal;
          if( meas_old_cal == disp_prev_cal )
          {
            new_meas_cal = new_disp_cal;
          }else if( meas_old_cal->type() == SpecUtils::EnergyCalType::LowerChannelEdge )
          {
            // Apply the SAME linear energy remap the user made to the displayed calibration, to
            //  this lower-channel target - expressed relative to the targets own (session-tracked)
            //  original channel energies, so its adjustment stays relative to its original.
            //  The displayed change maps energies as E_new = A*E_prev + B; composing with the
            //  targets current cumulative (off_t, gain_t) [E = off_t + gain_t*E_orig] gives
            //    off_t_new = A*off_t + B,   gain_t_new = A*gain_t.
            bool have_remap = false;
            double remap_a = 1.0, remap_b = 0.0;

            if( disp_prev_cal->type() == SpecUtils::EnergyCalType::LowerChannelEdge )
            {
              const LowerChanCalOriginal previnfo = lowerChannelOriginal( disp_prev_cal );
              const LowerChanCalOriginal newinfo = lowerChannelOriginal( new_disp_cal );

              // Only an offset/gain style change if both displayed cals trace to the same original
              //  (e.g., not a CALp file with un-related exact channel energies)
              if( (newinfo.original == previnfo.original) && (previnfo.gain != 0.0) )
              {
                have_remap = true;
                remap_a = newinfo.gain / previnfo.gain;
                remap_b = newinfo.offset - previnfo.offset*remap_a;
              }
            }else if( isOffsetOnly )
            {
              // Displayed cal is polynomial/FRF; for both, coefficient 0 is the keV offset (a
              //  pure energy shift, so A = 1)
              have_remap = true;
              remap_a = 1.0;
              remap_b = new_disp_cal->coefficients()[0] - disp_prev_cal->coefficients()[0];
            }

            if( have_remap )
            {
              LowerChanCalOriginal info = lowerChannelOriginal( meas_old_cal );
              info.offset = remap_a*info.offset + remap_b;
              info.gain = remap_a*info.gain;

              const auto adjusted = EnergyCal::adjust_lower_channel_energy_cal( info.original,
                                                                        info.offset, info.gain );
              registerLowerChannelAdjustment( adjusted, info );
              new_meas_cal = adjusted;
            }else
            {
              // A general calibration change: the target gets re-derived by channel mapping, and
              //  becomes the new nominal (it is deliberately NOT registered as an offset/gain
              //  adjustment of its original)
              new_meas_cal = EnergyCal::propogate_energy_cal_change( disp_prev_cal, new_disp_cal,
                                                                     meas_old_cal );
            }
          }else if( isOffsetOnly && (disp_prev_cal->type() != SpecUtils::EnergyCalType::LowerChannelEdge) )
          {
            const vector<float> &new_disp_coefs = new_disp_cal->coefficients();
            const vector<float> &prev_disp_coefs = disp_prev_cal->coefficients();
            const vector<pair<float,float>> &dev_pairs = meas_old_cal->deviation_pairs();

            vector<float> new_coefs = meas_old_cal->coefficients();
            new_coefs[0] += (new_disp_coefs[0] - prev_disp_coefs[0]);

            auto cal = make_shared<EnergyCalibration>();
            switch( meas_old_cal->type() )
            {
              case SpecUtils::EnergyCalType::Polynomial:
              case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
                cal->set_polynomial( meas_old_cal->num_channels(), new_coefs, dev_pairs );
                break;

              case SpecUtils::EnergyCalType::FullRangeFraction:
                cal->set_full_range_fraction( meas_old_cal->num_channels(), new_coefs, dev_pairs );
                break;

              case SpecUtils::EnergyCalType::LowerChannelEdge:
              case SpecUtils::EnergyCalType::InvalidEquationType:
                assert( 0 );  //lower channel energy handled in its own branch above
                break;
            }//switch( meas_old_cal->type() )

            new_meas_cal = cal;
          }else
          {
            // Note: an offset-only change of a displayed lower-channel-energy calibration, onto a
            //  polynomial/FRF target, also takes this path (the generalized
            //  propogate_energy_cal_change handles lower-channel displayed calibrations)
            new_meas_cal = EnergyCal::propogate_energy_cal_change( disp_prev_cal, new_disp_cal, meas_old_cal );
          }
          assert( new_meas_cal && new_meas_cal->valid() );
          old_to_new_cals[meas_old_cal] = new_meas_cal;
        }//for( const string &detname : change.detectors )
      }//for( loop over sample numbers )
    }catch( std::exception &e )
    {
      WString msg = WString::tr("ect-change-made-invalid");
      if( (backgrnd && (backgrnd != change.meas)) || (secgrnd && (secgrnd != change.meas)) )
      {
        if( change.meas == forgrnd )
          msg.arg( WString::tr("ect-for-the-fore") );
        else if( change.meas == backgrnd )
          msg.arg( WString::tr("ect-for-the-back") );
        else if( change.meas == secgrnd )
          msg.arg( WString::tr("ect-for-the-sec") );
        else
          msg.arg( "" );
      }else
      {
        msg.arg( "" );
      }//if( it is necessary to say which spectrum had the error )
      
      msg.arg( e.what() );
      
      throw runtime_error( msg.toUTF8() );
    }//try catch
    
    
    // Now go through and translate the peaks, but we wont actually update them to the SpecMeas
    //  until we know we can update all the peaks
    
    const set<set<int>> samples_with_peaks = change.meas->sampleNumsWithPeaks();
    const set<set<int>> samples_with_hint_peaks = change.meas->sampleNumsWithAutomatedSearchPeaks();
    
    // We may not be updating all samples, so we will only update peaks who are owned
    //  by sample numbers that are all in the samples being updated.
    set<set<int>> peaksamples, hintPeakSamples;
    
    for( const set<int> &samples : samples_with_peaks )
    {
      bool all_samples = true;
      
      // Check if the peaks sample numbers are all getting re-calibrated
      for( auto sample_num_iter = begin(samples);
          all_samples && (sample_num_iter != end(samples));
          ++sample_num_iter )
      {
        all_samples = (change.sample_numbers.count(*sample_num_iter) != 0u);
      }
      
      if( all_samples )
        peaksamples.insert( samples );
    }//for( const set<int> &samples : samples_with_peaks )
    
    for( const set<int> &samples : samples_with_hint_peaks )
    {
      bool all_samples = true;

      // Check if the peaks sample numbers are all getting re-calibrated
      for( auto sample_num_iter = begin(samples);
          all_samples && (sample_num_iter != end(samples));
          ++sample_num_iter )
      {
        all_samples = (change.sample_numbers.count(*sample_num_iter) != 0u);
      }
      
      if( all_samples )
        hintPeakSamples.insert( samples );
    }//for( const set<int> &samples : samples_with_peaks )
    
    
    // The peaks position (i.e., mean channel number) is determined by
    //  #SpecFile::suggested_sum_energy_calibration, however we may be applying an energy
    //  calibration change to one or two detectors, that arent that detector, in which case we
    //  dont want to move the peaks.
    //  This is why we dont just use change.detectors.
    vector<string> display_detectors;
    if( change.meas == forgrnd )
    {
      display_detectors = m_interspec->detectorsToDisplay(SpecUtils::SpectrumType::Foreground);
    }else if( change.meas == backgrnd )
    {
      display_detectors = m_interspec->detectorsToDisplay(SpecUtils::SpectrumType::Background);
    }else if( change.meas == secgrnd )
    {
      display_detectors = m_interspec->detectorsToDisplay(SpecUtils::SpectrumType::SecondForeground);
    }
    
    if( display_detectors.empty() )
    {
      // We probably shouldnt ever get here, but if we do, it will probably be fine to just use
      //  change.detectors (
      cerr << __func__ <<  ": Apparently no changing For/Back/Sec spectrum - doesnt seem right!" << endl;
      assert(0);
      display_detectors.insert( end(display_detectors), begin(change.detectors), end(change.detectors) );
    }//if( foreground being change ) / else background / else
    
    
    for( const set<int> &samples : peaksamples )
    {
      // Peak sets are stored per sample-number set, which may only partially overlap the samples
      //  being changed; we shift a peak set exactly when its display calibration is one of the
      //  calibrations being changed.  This is a deliberate compromise - a partially-overlapping
      //  peak set doesnt belong to any single calibration, so there is no exact answer.
      auto oldpeaks = change.meas->peaks(samples);
      auto oldcal = change.meas->suggested_sum_energy_calibration( samples, display_detectors );
      
      if( !oldpeaks || oldpeaks->empty() || !oldcal || !oldcal->valid() )
      {
        if( !oldpeaks || !oldcal || !oldcal->valid() )
        {
          //development sanity check, shouldnt normally get here I dont think
          cerr << __func__ << "Failed to get peaks or oldcal!" << endl; //just for sanity check
          assert( 0 );
        }
        continue;
      }
      
      const auto newcal_pos = old_to_new_cals.find(oldcal);
      if( (newcal_pos == end(old_to_new_cals)) || !newcal_pos->second || !newcal_pos->second->valid() )
      {
        //development sanity check, shouldnt normally happen I dont think
        cerr << __func__ << "Failed to get newcal for peaks shift!" << endl;
        assert( 0 );
        continue;
      }
      
      const shared_ptr<const EnergyCalibration> newcal = newcal_pos->second;
      assert( newcal && newcal->valid() );
      
      if( oldcal == newcal )
      {
        cerr << __func__ <<  ": oldcal == newcal - skipping shifting peak" << endl;
        continue;
      }
      
      try
      {
        auto newpeaks = EnergyCal::translatePeaksForCalibrationChange( *oldpeaks, oldcal, newcal );
        updated_peaks[oldpeaks] = newpeaks;
      }catch( std::exception &e )
      {
        string msg = "There was an issue translating peaks for this energy change;"
        " not applying change.  Error: " + string(e.what());
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, msg.c_str() );
#endif
        
        throw runtime_error( msg );
      }//try / catch
    }//for( const set<int> &samples : peaksampels )
    
    
    // Do similar loop to the above, but for hint peaks
    for( const set<int> &samples : hintPeakSamples )
    {
      auto oldHintPeaks = change.meas->automatedSearchPeaks(samples);
      auto oldcal = change.meas->suggested_sum_energy_calibration( samples, display_detectors );
      
      if( !oldHintPeaks || oldHintPeaks->empty() || !oldcal || !oldcal->valid() )
      {
        if( !oldHintPeaks || !oldcal || !oldcal->valid() )
        {
          assert( 0 );  //shouldnt get here
        }
        continue;
      }
      
      const auto newcal_pos = old_to_new_cals.find(oldcal);
      if( (newcal_pos == end(old_to_new_cals)) || !newcal_pos->second || !newcal_pos->second->valid() )
      {
        assert( 0 ); //shouldnt get here
        continue;
      }
      
      const shared_ptr<const EnergyCalibration> newcal = newcal_pos->second;
      assert( newcal && newcal->valid() );
      if( oldcal == newcal )
        continue;
      
      try
      {
        auto newpeaks = EnergyCal::translatePeaksForCalibrationChange( *oldHintPeaks, oldcal, newcal );
        updated_hint_peaks[oldHintPeaks] = newpeaks;
      }catch( std::exception &e )
      {
        string msg = "There was an issue translating hint peaks for this energy change,"
        " still applying, Error: " + string(e.what());
  #if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, msg.c_str() );
  #endif
      }//try / catch
    }//for( const set<int> &samples : hintPeakSamples )
  }//for( const MeasToApplyCoefChangeTo &change : changemeas )
  
  if( old_to_new_cals.find(disp_prev_cal) == end(old_to_new_cals) )
  {
    //Shouldnt ever happen; check is for development
    string msg = "There was an internal error updating energy calibration - energy cal"
    " associated with GUI wasn't updated - energy calibration state is suspect";
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, msg.c_str() );
#endif
    
    m_interspec->logMessage( msg, 3 );
    
#if( PERFORM_DEVELOPER_CHECKS )
    assert( 0 );
#endif
  }//if( old_to_new_cals.find(disp_prev_cal) == end(old_to_new_cals) )
  
  // Now go through and actually set the energy calibrations; they should all be valid and computed,
  //  as should all the shifted peaks.
  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    assert( change.meas );
    
    meas_old_new_cal_t &meas_old_new_cal = undo_sentry.cal_info( change.meas );
    meas_old_new_peaks_t &meas_old_new_peaks = undo_sentry.peak_info( change.meas );
    meas_old_new_peaks_t &meas_old_new_hint_peaks = undo_sentry.hint_peak_info( change.meas );
    
    
    for( const int sample : change.sample_numbers )
    {
      for( const string &detname : change.detectors )
      {
        auto m = change.meas->measurement( sample, detname );
        if( !m || m->num_gamma_channels() <= 4 )
          continue;
        
        const auto measoldcal = m->energy_calibration();
        assert( measoldcal );
        
        auto iter = old_to_new_cals.find( measoldcal );
        if( iter == end(old_to_new_cals) )
        {
          //Shouldnt ever happen
          string msg = "There was an internal error updating energy calibration - precomputed"
          " calibration couldnt be found - energy calibration will not be fully updated";
#if( PERFORM_DEVELOPER_CHECKS )
          log_developer_error( __func__, msg.c_str() );
#endif
          
          m_interspec->logMessage( msg, 3 );
          assert( 0 );
          continue;
        }//if( we havent already computed a new energy cal )
        
        assert( iter->second );
        assert( iter->second->num_channels() == m->num_gamma_channels() );
        
        meas_old_new_cal.emplace_back( m, measoldcal, iter->second );
        
        change.meas->set_energy_calibration( iter->second, m );
      }//for( loop over detector names )
    }//for( loop over sample numbers )
    
    
    //Now actually set the updated peaks
    const set<set<int>> peaksamples = change.meas->sampleNumsWithPeaks();
    
    for( const set<int> &samples : peaksamples )
    {
      auto oldpeaks = change.meas->peaks(samples);
      if( oldpeaks )
      {
        const auto pos = updated_peaks.find(oldpeaks);
        if( pos == end(updated_peaks) )
        {
          if( oldpeaks && !oldpeaks->empty() )
            cerr << "Couldnt find an expected entry in updated_peaks" << endl;
        }else
        {
          meas_old_new_peaks.emplace_back( samples, *oldpeaks, pos->second );
          
          change.meas->setPeaks( pos->second, samples );
          if( m_peakModel && (change.meas == forgrnd) && (samples == foresamples) )
            m_peakModel->setPeakFromSpecMeas( forgrnd, foresamples, SpecUtils::SpectrumType::Foreground );
          else if( m_peakModel && (change.meas == backgrnd) && (samples == backSamples) )
            m_peakModel->setPeakFromSpecMeas( backgrnd, backSamples, SpecUtils::SpectrumType::Background );
          else if( m_peakModel && (change.meas == secgrnd) && (samples == secoSamples) )
            m_peakModel->setPeakFromSpecMeas(secgrnd, secoSamples, SpecUtils::SpectrumType::SecondForeground );
        }//if( pos == end(updated_peaks) ) / else
      }//if( oldpeaks )
    }//for( const set<int> &samples : peaksampels )
    
    
    // Do similar thing for the hint peaks
    const set<set<int>> hintPeakSamples = change.meas->sampleNumsWithAutomatedSearchPeaks();
    for( const set<int> &samples : hintPeakSamples )
    {
      shared_ptr<const SpecMeas::PeakDeque> hintPeaks = change.meas->automatedSearchPeaks(samples);
      assert( hintPeaks && !hintPeaks->empty() );
      
      if( hintPeaks && !hintPeaks->empty() )
      {
        const auto pos = updated_hint_peaks.find( hintPeaks );
        if( pos == end(updated_hint_peaks) )
        {
          if( hintPeaks && !hintPeaks->empty() )
            cerr << "Couldnt find an expected entry in updated_peaks" << endl;
        }else if( !pos->second.empty() )
        {
          meas_old_new_hint_peaks.emplace_back( samples, *hintPeaks, pos->second );
          
          auto peaks = make_shared<deque<shared_ptr<const PeakDef>>>( pos->second );
          change.meas->setAutomatedSearchPeaks( samples, peaks );
        }//if( pos == end(updated_hint_peaks) ) / else
      }//if( hintPeaks )
    }//for( const set<int> &samples : hintPeakSamples )
  }//for( loop over SpecFiles for change )
  
  m_interspec->refreshDisplayedCharts();
  refreshGuiFromFiles();
}//applyCalChange(...)


void EnergyCalTool::setEnergyCal( shared_ptr<const SpecUtils::EnergyCalibration> new_cal,
                   const MeasToApplyCoefChangeTo &changemeas, const bool adjust_peaks )
{
  // We get here from applying CALp to files, sometimes.
  using namespace SpecUtils;
  
  assert( new_cal && new_cal->valid() );
  
  const auto forgrnd = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  const auto backgrnd = m_interspec->measurment(SpecUtils::SpectrumType::Background);
  const auto secgrnd = m_interspec->measurment(SpecUtils::SpectrumType::SecondForeground);
  
  // Currently its pointless to map old energy calibrations to the new_cal being set, but in the
  //  future we might be a little more lenient about matching number of channels or whatever, so
  //  we'll just toss this mechanism in for now.
  map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> old_to_new_cals;
  
  // We will store updated peaks and not set any of them until we know all the energy calibrations
  //  and peak shifts were successfully done.
  map<shared_ptr<deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_peaks;
  map<shared_ptr<const deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_hint_peaks;
  
  
  // const vector<MeasToApplyCoefChangeTo> changemeas = measurementsToApplyCoeffChangeTo();
  
  //First calculate new calibrations, and make sure they are valid, then actually set them.
  assert( changemeas.meas );
  if( !changemeas.meas )
    throw logic_error( "EnergyCalTool::setEnergyCal: nullptr measurement to change passed in." );
  
  SpecUtils::SpectrumType spectype;
  if( changemeas.meas == forgrnd )
    spectype = SpecUtils::SpectrumType::Foreground;
  else if( changemeas.meas == backgrnd )
    spectype = SpecUtils::SpectrumType::Background;
  else if( changemeas.meas == secgrnd )
    spectype = SpecUtils::SpectrumType::SecondForeground;
  else
  {
    assert( 0 );
    throw runtime_error( "EnergyCalTool::setEnergyCal: measurement to change must be foreground,"
                         " background, or secondary spec displayed" );
  }
  
  
  const set<int> &dispsamples = m_interspec->displayedSamples( spectype );
  const vector<string> dispdets = m_interspec->detectorsToDisplay( spectype );
  const shared_ptr<const SpecUtils::EnergyCalibration> display_cal
                       = changemeas.meas->suggested_sum_energy_calibration( dispsamples, dispdets );
  
  try
  {
    for( const int sample : changemeas.sample_numbers )
    {
      for( const string &detname : changemeas.detectors )
      {
        auto m = changemeas.meas->measurement( sample, detname );
        if( !m || m->num_gamma_channels() <= 4 )
          continue;
          
        const auto meas_old_cal = m->energy_calibration();
        assert( meas_old_cal );
          
        if( !meas_old_cal || !meas_old_cal->valid() )
          continue;
          
        // TODO: maybe be a little more lenient on matching number of channels, by possible making
        //       new energy calibrations when possible.
        
        if( meas_old_cal->num_channels() != new_cal->num_channels() )
        {
          throw runtime_error( WString::tr("ect-set-cal-invalid-nchan")
                              .arg(sample)
                              .arg(detname)
                              .arg( static_cast<int>(meas_old_cal->num_channels()) )
                              .arg( static_cast<int>(new_cal->num_channels()) )
                              .toUTF8() );
        }
        //if( display_cal == meas_old_cal )
        //  adjust_peaks = true;
        
        old_to_new_cals[meas_old_cal] = new_cal;
      }//for( const string &detname : change.detectors )
    }//for( loop over sample numbers )
  }catch( std::exception &e )
  {
    throw runtime_error( WString::tr("ect-error-setting-cal").arg(e.what()).toUTF8() );
  }//try catch
  
  // Now go through and translate the peaks, but we wont actually update them to the SpecMeas
  //  until we know we can update all the peaks
  const set<set<int>> peaksamples = adjust_peaks ? changemeas.meas->sampleNumsWithPeaks()
                                                 : set<set<int>>();
    
  // The peaks position (i.e., mean channel number) is determined by
  //  #SpecFile::suggested_sum_energy_calibration, however we may be applying an energy
  //  calibration change to one or two detectors, that arent that detector, in which case we
  //  dont want to move the peaks.
  for( const set<int> &samples : peaksamples )
  {
    // Peak sets are stored per sample-number set, which may only partially overlap the samples
    //  being changed; we shift a peak set exactly when its display calibration is one of the
    //  calibrations being changed.  This is a deliberate compromise - a partially-overlapping
    //  peak set doesnt belong to any single calibration, so there is no exact answer.
    auto oldpeaks = changemeas.meas->peaks(samples);
    auto oldcal = changemeas.meas->suggested_sum_energy_calibration( samples, dispdets );
      
    if( !oldpeaks || oldpeaks->empty() || !oldcal || !oldcal->valid() )
    {
      if( !oldpeaks || !oldcal || !oldcal->valid() )
      {
        // we can get here if we dont have any display detectors
        cerr << __func__ << "Failed to get peaks or oldcal!" << endl; //just for sanity check
      }
      continue;
    }
    
    const auto newcal_pos = old_to_new_cals.find(oldcal);
    if( (newcal_pos == end(old_to_new_cals)) || !newcal_pos->second || !newcal_pos->second->valid() )
    {
      // This peak-sets display calibration isnt one of the calibrations being changed
      continue;
    }
      
    const shared_ptr<const EnergyCalibration> newcal = newcal_pos->second;
    assert( newcal && newcal->valid() );
      
    if( oldcal == newcal )
    {
      assert( 0 );
      continue;
    }
      
    try
    {
      auto newpeaks = EnergyCal::translatePeaksForCalibrationChange( *oldpeaks, oldcal, newcal );
      updated_peaks[oldpeaks] = newpeaks;
    }catch( std::exception &e )
    {
      string msg = "There was an issue translating peaks for this energy change;"
        " not applying change.  Error: " + string(e.what());
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, msg.c_str() );
#endif
        
      throw runtime_error( msg );
    }//try / catch
  }//for( const set<int> &samples : peaksampels )
  
  
  // Now do same thing as for peaks, but for the hint peaks
  const set<set<int>> hintPeakSamples = adjust_peaks
                                            ? changemeas.meas->sampleNumsWithAutomatedSearchPeaks()
                                            : set<set<int>>();
  for( const set<int> &samples : hintPeakSamples )
  {
    auto oldHintPeaks = changemeas.meas->automatedSearchPeaks( samples );
    auto oldcal = changemeas.meas->suggested_sum_energy_calibration( samples, dispdets );
      
    if( !oldHintPeaks || oldHintPeaks->empty() || !oldcal || !oldcal->valid() )
      continue;
    
    const auto newcal_pos = old_to_new_cals.find(oldcal);
    if( (newcal_pos == end(old_to_new_cals)) || !newcal_pos->second || !newcal_pos->second->valid() )
      continue;
      
    const shared_ptr<const EnergyCalibration> newcal = newcal_pos->second;
    assert( newcal && newcal->valid() );
      
    if( !newcal || !newcal->valid() || (oldcal == newcal) )
    {
      assert( 0 );
      continue;
    }
    
    try
    {
      auto newpeaks = EnergyCal::translatePeaksForCalibrationChange( *oldHintPeaks, oldcal, newcal );
      updated_hint_peaks[oldHintPeaks] = newpeaks;
    }catch( std::exception &e )
    {
#if( PERFORM_DEVELOPER_CHECKS )
      string msg = "There was an issue translating hint peaks for this energy change.  Error: "
                   + string(e.what());
      log_developer_error( __func__, msg.c_str() );
#endif
    }//try / catch
  }//for( const set<int> &samples : hintPeakSamples )
  
    
  // For undo/redo, store changes to energy cal, and peaks.
  EnergyCalUndoRedoSentry undo_sentry;
  
  // Now go through and actually set the energy calibrations; they should all be valid and computed,
  //  as should all the shifted peaks
  for( const int sample : changemeas.sample_numbers )
  {
    meas_old_new_cal_t &meas_old_new_cal = undo_sentry.cal_info(changemeas.meas);
    
    for( const string &detname : changemeas.detectors )
    {
      auto m = changemeas.meas->measurement( sample, detname );
      if( !m || m->num_gamma_channels() <= 4 )
        continue;
        
      const auto measoldcal = m->energy_calibration();
      assert( measoldcal );
        
      auto iter = old_to_new_cals.find( measoldcal );
      if( iter == end(old_to_new_cals) )
      {
        //Shouldnt ever happen
        string msg = "There was an internal error updating energy calibration - precomputed"
        " calibration couldnt be found - energy calibration will not be fully updated";
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, msg.c_str() );
#endif
          
        m_interspec->logMessage( msg, 3 );
        assert( 0 );
        continue;
      }//if( we havent already computed a new energy cal )
        
      assert( iter->second );
      assert( iter->second->num_channels() == m->num_gamma_channels() );
        
      meas_old_new_cal.emplace_back( m, measoldcal, iter->second );
      
      changemeas.meas->set_energy_calibration( iter->second, m );
    }//for( loop over detector names )
  }//for( loop over sample numbers )
    
    
  //Now actually set the updated peaks
  const set<int> &foreSamples = m_interspec->displayedSamples( SpecUtils::SpectrumType::Foreground );
  const set<int> &backSamples = m_interspec->displayedSamples( SpecUtils::SpectrumType::Background );
  const set<int> &secoSamples = m_interspec->displayedSamples( SpecUtils::SpectrumType::SecondForeground );
  
  for( const set<int> &samples : peaksamples )
  {
    auto oldpeaks = changemeas.meas->peaks(samples);
    if( oldpeaks )
    {
      const auto pos = updated_peaks.find(oldpeaks);
      if( pos == end(updated_peaks) )
      {
        if( oldpeaks && !oldpeaks->empty() )
          cerr << "Couldnt find an expected entry in updated_peaks" << endl;
      }else
      {
        meas_old_new_peaks_t &meas_old_new_peaks = undo_sentry.peak_info( changemeas.meas );
        meas_old_new_peaks.emplace_back( samples, *oldpeaks, pos->second );
        
        changemeas.meas->setPeaks( pos->second, samples );
        if( m_peakModel && (changemeas.meas == forgrnd) && (samples == foreSamples) )
          m_peakModel->setPeakFromSpecMeas( forgrnd, foreSamples, SpecUtils::SpectrumType::Foreground );
        else if( m_peakModel && (changemeas.meas == backgrnd) && (samples == backSamples) )
          m_peakModel->setPeakFromSpecMeas( backgrnd, backSamples, SpecUtils::SpectrumType::Background );
        else if( m_peakModel && (changemeas.meas == secgrnd) && (samples == secoSamples) )
          m_peakModel->setPeakFromSpecMeas( secgrnd, secoSamples, SpecUtils::SpectrumType::SecondForeground );
      }//if( pos == end(updated_peaks) ) / else
    }//if( oldpeaks )
  }//for( const set<int> &samples : peaksampels )
  
  // And set the updated hint peaks
  for( const set<int> &samples : hintPeakSamples )
  {
    auto hintPeaks = changemeas.meas->automatedSearchPeaks(samples);
    if( hintPeaks && !hintPeaks->empty() )
    {
      const auto pos = updated_hint_peaks.find(hintPeaks);
      if( pos == end(updated_hint_peaks) )
      {
        if( hintPeaks && !hintPeaks->empty() )
          cerr << "Couldnt find an expected entry in updated_hint_peaks" << endl;
      }else
      {
        meas_old_new_peaks_t &meas_old_new_hint_peaks = undo_sentry.hint_peak_info( changemeas.meas );
        meas_old_new_hint_peaks.emplace_back( samples, *hintPeaks, pos->second );
        auto peaks = make_shared<deque<shared_ptr<const PeakDef>>>( pos->second );
        changemeas.meas->setAutomatedSearchPeaks( samples, peaks );
      }//if( pos == end(updated_hint_peaks) ) / else
    }//if( hintPeaks && !hintPeaks->empty() )
  }//for( const set<int> &samples : peaksampels )
}//void setEnergyCal( new_cal, changemeas )


void EnergyCalTool::addDeviationPair( const std::pair<float,float> &new_pair )
{
  // We get here when we graphically (i.e., cntrl+option+drag) add in a deviation pair
  
  using namespace SpecUtils;
  
  const auto forgrnd = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  const auto backgrnd = m_interspec->measurment(SpecUtils::SpectrumType::Background);
  const auto secgrnd = m_interspec->measurment(SpecUtils::SpectrumType::SecondForeground);
  
  const set<int> &foreSamples = m_interspec->displayedSamples( SpectrumType::Foreground );
  const set<int> &backSamples = m_interspec->displayedSamples( SpectrumType::Background );
  const set<int> &secoSamples = m_interspec->displayedSamples( SpectrumType::SecondForeground );
  
  // We will calculate all new energy calibrations, to make sure we actually can, befor setting them
  map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> old_to_new_cals;
  
  // We will also pre-calculate updated peaks
  map<shared_ptr<deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_peaks;
  map<shared_ptr<const deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_hint_peaks;
  
  const vector<MeasToApplyCoefChangeTo> changemeas = measurementsToApplyCoeffChangeTo();
  
  // For undo/redo, store changes to energy cal, and peaks.
  //  Note: we only record into the sentry in the second (commit) loop, so a throw from this first
  //  loop doesnt insert an undo/redo step for changes that were never applied.
  EnergyCalUndoRedoSentry undo_sentry;

  // Do first loop to calculate new calibrations
  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    assert( change.meas );

    try
    {
      for( const int sample : change.sample_numbers )
      {
        for( const string &detname : change.detectors )
        {
          auto m = change.meas->measurement( sample, detname );
          if( !m || m->num_gamma_channels() <= 4 )
            continue;
          
          const auto old_cal = m->energy_calibration();
          assert( old_cal );
          
          if( !m || (m->num_gamma_channels() <= 4) || !old_cal || !old_cal->valid()
             || (old_cal->type() == EnergyCalType::LowerChannelEdge) )
            continue;
          
          //If we have already computed the new calibration for a EnergyCalibration object, lets not
          //  re-due it.
          if( old_to_new_cals.count(old_cal) )
            continue;
          
          const size_t nchannel = old_cal->num_channels();
          const vector<float> &coefs = old_cal->coefficients();
          vector<pair<float,float>> dev_pairs = old_cal->deviation_pairs();
          if( dev_pairs.empty() )
            dev_pairs.push_back( {0.0f, 0.0f} );
          dev_pairs.push_back( new_pair );
          
          std::sort( begin(dev_pairs), end(dev_pairs) );
          
          shared_ptr<EnergyCalibration> new_cal = make_shared<EnergyCalibration>();
          
          switch( old_cal->type() )
          {
            case SpecUtils::EnergyCalType::Polynomial:
            case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
              new_cal->set_polynomial( nchannel, coefs, dev_pairs );
              break;
              
            case SpecUtils::EnergyCalType::FullRangeFraction:
              new_cal->set_full_range_fraction( nchannel, coefs, dev_pairs );
              break;
              
            case SpecUtils::EnergyCalType::LowerChannelEdge:
            case SpecUtils::EnergyCalType::InvalidEquationType:
              assert( 0 );
              break;
          }//switch( meas_old_cal->type() )
          
          assert( new_cal && new_cal->valid() );
          old_to_new_cals[old_cal] = new_cal;
        }//for( const string &detname : change.detectors )
      }//for( loop over sample numbers )
    }catch( std::exception &e )
    {
      WString msg = WString::tr("ect-add-dev-pair-made-invalid");
      if( (backgrnd && (backgrnd != change.meas)) || (secgrnd && (secgrnd != change.meas)) )
      {
        if( change.meas == forgrnd )
          msg.arg( WString::tr("ect-for-the-fore") );
        else if( change.meas == backgrnd )
          msg.arg( WString::tr("ect-for-the-back") );
        else if( change.meas == secgrnd )
          msg.arg( WString::tr("ect-for-the-sec") );
        else
          msg.arg( "" );
      }else
      {
        msg.arg( "" );
      }//if( it is necessary to say which spectrum had the error ) / else
      
      throw runtime_error( msg.arg(e.what()).toUTF8() );
    }//try catch
    
    
    // Now go through and translate the peaks, but we wont actually update them to the SpecMeas
    //  until we know we can update all the peaks
    const set<set<int>> peaksamples = change.meas->sampleNumsWithPeaks();
    const vector<string> detnamesv( begin(change.detectors), end(change.detectors) );
    
    for( const set<int> &samples : peaksamples )
    {
      // Peak sets are stored per sample-number set, which may only partially overlap the samples
      //  being changed; we shift a peak set exactly when its display calibration is one of the
      //  calibrations being changed.  This is a deliberate compromise - a partially-overlapping
      //  peak set doesnt belong to any single calibration, so there is no exact answer.
      auto oldpeaks = change.meas->peaks(samples);
      auto oldcal = change.meas->suggested_sum_energy_calibration( samples, detnamesv );
      
      if( !oldpeaks || oldpeaks->empty()
         || !oldcal || (oldcal->type() == EnergyCalType::InvalidEquationType) )
      {
        if( !oldpeaks || !oldcal || !oldcal->valid() )
          cerr << "Failed to get peaks or oldcal!" << endl; //just for development
        continue;
      }
      
      const auto newcal_pos = old_to_new_cals.find(oldcal);
      if( (newcal_pos == end(old_to_new_cals)) || !newcal_pos->second || !newcal_pos->second->valid() )
      {
        cerr << "Failed to get newcal for peaks shift!" << endl; //just for development, shouldnt happen I dont think
        continue;
      }
      
      const shared_ptr<const EnergyCalibration> newcal = newcal_pos->second;
      assert( newcal && newcal->valid() );
      
      if( oldcal == newcal )
      {
        cerr << __func__ <<  ": oldcal == newcal - skipping shifting peak" << endl;
        continue;
      }
      
      try
      {
        auto newpeaks = EnergyCal::translatePeaksForCalibrationChange( *oldpeaks, oldcal, newcal );
        updated_peaks[oldpeaks] = newpeaks;
      }catch( std::exception &e )
      {
        string msg = "There was an issue translating peaks for this energy change;"
        " not applying change.  Error: " + string(e.what());
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, msg.c_str() );
#endif
        
        throw runtime_error( msg );
      }//try / catch
    }//for( const set<int> &samples : peaksampels )
    
    
    //Now get new hint peaks
    const set<set<int>> hintPeakSamples = change.meas->sampleNumsWithAutomatedSearchPeaks();
    
    for( const set<int> &samples : hintPeakSamples )
    {
      auto oldHintPeaks = change.meas->automatedSearchPeaks(samples);
      auto oldcal = change.meas->suggested_sum_energy_calibration( samples, detnamesv );
      
      if( !oldHintPeaks || oldHintPeaks->empty()
         || !oldcal || (oldcal->type() == EnergyCalType::InvalidEquationType) )
      {
        if( !oldHintPeaks || !oldcal || !oldcal->valid() )
          cerr << "Failed to get peaks or oldcal!" << endl; //just for development
        continue;
      }
      
      const auto newcal_pos = old_to_new_cals.find(oldcal);
      if( (newcal_pos == end(old_to_new_cals)) || !newcal_pos->second || !newcal_pos->second->valid() )
      {
        cerr << "Failed to get newcal for peaks shift!" << endl; //just for development, shouldnt happen I dont think
        continue;
      }
      
      const shared_ptr<const EnergyCalibration> newcal = newcal_pos->second;
      assert( newcal && newcal->valid() );
      if( !newcal || !newcal->valid() || (oldcal == newcal) )
        continue;
      
      try
      {
        auto newpeaks = EnergyCal::translatePeaksForCalibrationChange( *oldHintPeaks, oldcal, newcal );
        updated_hint_peaks[oldHintPeaks] = newpeaks;
      }catch( std::exception &e )
      {
        string msg = "There was an issue translating hint peaks for this energy change/"
        "  Error: " + string(e.what());
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, msg.c_str() );
#endif
      }//try / catch
    }//for( const set<int> &samples : peaksampels )
  }//for( const MeasToApplyCoefChangeTo &change : changemeas )
  
  
  // Now go through and actually set the energy calibrations; they should all be valid and computed,
  //  as should all the shifted peaks.
  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    assert( change.meas );
    meas_old_new_cal_t &meas_old_new_cal = undo_sentry.cal_info( change.meas );
    meas_old_new_peaks_t &meas_old_new_peaks = undo_sentry.peak_info( change.meas );
    meas_old_new_peaks_t &meas_old_new_hint_peaks = undo_sentry.hint_peak_info( change.meas );

    for( const int sample : change.sample_numbers )
    {
      for( const string &detname : change.detectors )
      {
        auto m = change.meas->measurement( sample, detname );
        if( !m || m->num_gamma_channels() <= 4 )
          continue;
        
        const auto measoldcal = m->energy_calibration();
        assert( measoldcal );
        
        auto iter = old_to_new_cals.find( measoldcal );
        if( iter == end(old_to_new_cals) )
        {
          //Shouldnt ever happen
          string msg = "There was an internal error updating energy calibration - precomputed"
          " calibration couldnt be found - energy calibation will not be fully updated";
#if( PERFORM_DEVELOPER_CHECKS )
          log_developer_error( __func__, msg.c_str() );
#endif
          
          m_interspec->logMessage( msg, 3 );
          assert( 0 );
          continue;
        }//if( we havent already computed a new energy cal )
        
        assert( iter->second );
        assert( iter->second->num_channels() == m->num_gamma_channels() );
        
        change.meas->set_energy_calibration( iter->second, m );
        
        meas_old_new_cal.emplace_back( m, measoldcal, iter->second );
      }//for( loop over detector names )
    }//for( loop over sample numbers )
    
    
    //Now actually set the updated peaks
    const set<set<int>> peaksamples = change.meas->sampleNumsWithPeaks();
    
    for( const set<int> &samples : peaksamples )
    {
      auto oldpeaks = change.meas->peaks(samples);
      if( !oldpeaks || oldpeaks->empty() )
        continue;
      
      const auto pos = updated_peaks.find(oldpeaks);
      if( pos == end(updated_peaks) )
      {
        if( oldpeaks && !oldpeaks->empty() )
          cerr << "Couldnt find an expected entry in updated_peaks" << endl;
        continue;
      }

      meas_old_new_peaks.emplace_back( samples, *oldpeaks, pos->second );

      change.meas->setPeaks( pos->second, samples );
      if( m_peakModel && (change.meas == forgrnd) && (samples == foreSamples) )
        m_peakModel->setPeakFromSpecMeas( forgrnd, foreSamples, SpectrumType::Foreground );
      else if( m_peakModel && (change.meas == backgrnd) && (samples == backSamples) )
        m_peakModel->setPeakFromSpecMeas( backgrnd, backSamples, SpecUtils::SpectrumType::Background );
      else if( m_peakModel && (change.meas == secgrnd) && (samples == secoSamples) )
        m_peakModel->setPeakFromSpecMeas( secgrnd, secoSamples, SpecUtils::SpectrumType::SecondForeground );
    }//for( const set<int> &samples : peaksampels )
    
    // Also grab the updated hint peaks
    const set<set<int>> hintPeakSamples = change.meas->sampleNumsWithAutomatedSearchPeaks();
    for( const set<int> &samples : hintPeakSamples )
    {
      auto oldHintPeaks = change.meas->automatedSearchPeaks( samples );
      if( !oldHintPeaks || oldHintPeaks->empty() )
        continue;
      
      const auto pos = updated_hint_peaks.find(oldHintPeaks);
      if( pos == end(updated_hint_peaks) )
      {
        if( oldHintPeaks && !oldHintPeaks->empty() )
          cerr << "Couldnt find an expected entry in updated_hint_peaks" << endl;
      }else
      {
        meas_old_new_hint_peaks.emplace_back( samples, *oldHintPeaks, pos->second );

        auto peaks = make_shared<deque<shared_ptr<const PeakDef>>>( pos->second );
        change.meas->setAutomatedSearchPeaks( samples, peaks );
      }
    }//for( const set<int> &samples : peaksampels )
  }//for( loop over SpecFiles for change )
  

  m_interspec->refreshDisplayedCharts();
  refreshGuiFromFiles();
}//void addDeviationPair( const std::pair<float,float> &new_pair );



void EnergyCalTool::userChangedCoefficient( const size_t coefnum, EnergyCalImp::CalDisplay *display )
{
  using namespace SpecUtils;
  assert( coefnum < 10 );  //for lower channel energy cals, coefnum is 0 (offset) or 1 (gain)
  
  shared_ptr<const EnergyCalibration> disp_prev_cal = display->lastSetCalibration();
  if( !disp_prev_cal )
  {
    cerr << "unexpected error getting updated energy calibration coefficients" << endl;
    m_interspec->logMessage( WString::tr("ect-unexpected-error-prev-coefs"), 2 );
    doRefreshFromFiles();
    return;
  }//if( !disp_prev_cal )
  
  
  {// Begin check to make sure the changed energy cal is actually checked for it to be applied to...
    bool willBeAppliedToDisplay = false;
    const std::string &detname = display->detectorName();
    const SpecUtils::SpectrumType type = display->spectrumType();
    std::shared_ptr<SpecMeas> cal_disp_meas = m_interspec->measurment(type);
    assert( cal_disp_meas );
    
    const vector<MeasToApplyCoefChangeTo> applyTo = measurementsToApplyCoeffChangeTo();
    
    for( const MeasToApplyCoefChangeTo &delta : applyTo )
    {
      const shared_ptr<SpecMeas> &meas = delta.meas;
      if( meas != cal_disp_meas )
        continue;
      
      const set<string> &detectors = delta.detectors;
      if( !detectors.count(detname) )
        continue;
      
      // Actually, I think if we're here, we're probably good, but we'll check a little deeper to
      //   make sure check in EnergyCalTool::applyCalChange will be satisfied
      const set<int> &samples = delta.sample_numbers;
      for( const int sample : samples )
      {
        auto m = meas->measurement( sample, detname);
        willBeAppliedToDisplay = (m && (m->energy_calibration() == disp_prev_cal));
        if( willBeAppliedToDisplay )
          break;
      }//for( const int sample : samples )
    
      if( willBeAppliedToDisplay )
        break;
    }//for( loop over changes )
    
    if( !willBeAppliedToDisplay )
    {
      m_interspec->logMessage( WString::tr("ect-changed-cal-not-selected"), 2 );
      doRefreshFromFiles();
      return;
    }
  }// End check to make sure the changed energy cal is actually checked for it to be applied to...
  
  
  m_lastGraphicalRecal = 0;
  m_lastGraphicalRecalType = EnergyCalGraphicalConfirm::NumRecalTypes;
  m_lastGraphicalRecalEnergy = -999.0f;
  
  vector<float> dispcoefs = display->displayedCoefficents();
  if( dispcoefs.size() <= coefnum )
    dispcoefs.resize( coefnum+1, 0.0f );

  shared_ptr<const EnergyCalibration> new_disp_cal;
  try
  {
    switch( disp_prev_cal->type() )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      case SpecUtils::EnergyCalType::FullRangeFraction:
      {
        vector<float> prev_disp_coefs = disp_prev_cal->coefficients();
        if( prev_disp_coefs.size() <= coefnum )
          prev_disp_coefs.resize( coefnum+1, 0.0f );

        vector<float> new_disp_coefs = prev_disp_coefs;
        new_disp_coefs[coefnum] = dispcoefs[coefnum];

        const size_t dispnchannel = disp_prev_cal->num_channels();
        const auto &disp_dev_pairs = disp_prev_cal->deviation_pairs();

        auto cal = make_shared<EnergyCalibration>();
        if( disp_prev_cal->type() == SpecUtils::EnergyCalType::FullRangeFraction )
          cal->set_full_range_fraction( dispnchannel, new_disp_coefs, disp_dev_pairs );
        else
          cal->set_polynomial( dispnchannel, new_disp_coefs, disp_dev_pairs );

        new_disp_cal = cal;
        break;
      }//case polynomial or full range fraction

      case SpecUtils::EnergyCalType::LowerChannelEdge:
      {
        // The displayed {offset, gain} are the cumulative adjustments, relative to the sessions
        //  ORIGINAL channel energies - so we always re-adjust from the original.
        if( coefnum >= 2 )
          throw runtime_error( "Unexpected lower channel energy coefficient number" );

        dispcoefs.resize( 2, 0.0f );

        LowerChanCalOriginal info = lowerChannelOriginal( disp_prev_cal );
        info.offset = dispcoefs[0];
        info.gain = dispcoefs[1];

        const auto adjusted = EnergyCal::adjust_lower_channel_energy_cal( info.original,
                                                                      info.offset, info.gain );
        registerLowerChannelAdjustment( adjusted, info );
        new_disp_cal = adjusted;
        break;
      }//case lower channel energy

      case SpecUtils::EnergyCalType::InvalidEquationType:
        throw runtime_error( "Invalid calibration type changed?  Something is way wack." );
        break;
    }//switch( disp_prev_cal->type() )
  }catch( std::exception &e )
  {
    display->updateToGui( disp_prev_cal );
    m_interspec->logMessage( WString::tr("ect-change-made-invalid").arg("").arg(e.what()), 2 );

    return;
  }//try / catch to create new_disp_cal
  
  assert( new_disp_cal && new_disp_cal->valid() );
  
  try
  {
    const vector<MeasToApplyCoefChangeTo> changemeas = measurementsToApplyCoeffChangeTo();
    applyCalChange( disp_prev_cal, new_disp_cal, changemeas, coefnum==0 );
  }catch( std::exception &e )
  {
    display->updateToGui( disp_prev_cal );
    m_interspec->logMessage( WString::tr("ect-change-made-invalid").arg("").arg(e.what()), 2 );
  }//try / catch
}//userChangedCoefficient(...)


void EnergyCalTool::userChangedDeviationPair( EnergyCalImp::CalDisplay *display, const int fieldTypeChanged )
{
  using namespace SpecUtils;

  const auto field = EnergyCalImp::DeviationPairDisplay::UserFieldChanged( fieldTypeChanged );

  switch( field )
  {
    case EnergyCalImp::DeviationPairDisplay::UserFieldChanged::AddedDeviationPair:
      return;  //hasnt been filled out yet; no need to do anything

    case EnergyCalImp::DeviationPairDisplay::UserFieldChanged::RemovedDeviationPair:
    case EnergyCalImp::DeviationPairDisplay::UserFieldChanged::EnergyChanged:
    case EnergyCalImp::DeviationPairDisplay::UserFieldChanged::OffsetChanged:
      break;
  };//enum UserFieldChanged


  m_lastGraphicalRecal = 0;
  m_lastGraphicalRecalType = EnergyCalGraphicalConfirm::NumRecalTypes;
  m_lastGraphicalRecalEnergy = -999.0f;

  const shared_ptr<SpecMeas> foregrnd = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  const shared_ptr<SpecMeas> backgrnd = m_interspec->measurment(SpecUtils::SpectrumType::Background);
  const shared_ptr<SpecMeas> secogrnd = m_interspec->measurment(SpecUtils::SpectrumType::SecondForeground);

  const set<int> &foreSamples = m_interspec->displayedSamples( SpectrumType::Foreground );
  const set<int> &backSamples = m_interspec->displayedSamples( SpectrumType::Background );
  const set<int> &secoSamples = m_interspec->displayedSamples( SpectrumType::SecondForeground );

  assert( foregrnd );
  if( !foregrnd )
    return;

  const shared_ptr<const EnergyCalibration> old_cal = display->lastSetCalibration();
  const vector<pair<float,float>> old_dev_pairs = old_cal
                                        ? old_cal->deviation_pairs() : vector<pair<float,float>>{};
  const vector<pair<float,float>> new_dev_pairs = display->displayedDeviationPairs();

  // After the user clicks to add a new deviation pair, and then fills in one of the fields, if the
  //  other field is not yet filled in, then the GUI wont insert that deviation pair; so for this
  //  case where the deviation pairs havent yet changed, lets skip doing anything yet
  if( new_dev_pairs.size() == old_dev_pairs.size() )
  {
    bool equal = true;
    for( size_t i = 0; equal && (i < new_dev_pairs.size()); ++i )
      equal = ( (fabs(new_dev_pairs[i].first - old_dev_pairs[i].first) < 1.0E-4)
                && (fabs(new_dev_pairs[i].second - old_dev_pairs[i].second) < 1.0E-4) );
    if( equal )
      return;
  }//if( new_dev_pairs.size() == old_dev_pairs.size() )


  const SpectrumType type = display->spectrumType();
  const std::string &detname = display->detectorName();

  const shared_ptr<SpecMeas> dispmeas = m_interspec->measurment( type );
  if( !dispmeas )  //Shouldnt ever happen
  {
    display->updateToGui( old_cal );
    m_interspec->logMessage( "Internal error retrieving correct measurement", 2 );
    return;
  }

  // The change is applied according to the "Apply Changes To" selections, the same as coefficient
  //  edits are.  Each selected measurements deviation pairs are REPLACED with the displayed list,
  //  keeping that measurements own coefficients - a full-list edit has no well-defined
  //  per-measurement delta, so wholesale replacement is the only consistent semantic.
  const vector<MeasToApplyCoefChangeTo> changemeas = measurementsToApplyCoeffChangeTo();

  {// Begin check the display that was edited is actually selected to have changes applied to it
    bool willBeApplied = false;
    for( size_t i = 0; !willBeApplied && (i < changemeas.size()); ++i )
      willBeApplied = ( (changemeas[i].meas == dispmeas) && changemeas[i].detectors.count(detname) );

    if( !willBeApplied )
    {
      display->updateToGui( old_cal );
      m_interspec->logMessage( WString::tr("ect-changed-cal-not-selected"), 2 );
      return;
    }
  }// End check the display that was edited is actually selected to have changes applied to it

  // We will compute all the updated calibrations, and translated peaks, and not set any of them
  //  until we know they can all be successfully computed.
  map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> old_to_new_cal;
  map<shared_ptr<deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_peaks;
  map<shared_ptr<const deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_hint_peaks;

  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    assert( change.meas );

    try
    {
      for( const int sample : change.sample_numbers )
      {
        for( const string &name : change.detectors )
        {
          const shared_ptr<const Measurement> m = change.meas->measurement( sample, name );
          if( !m || (m->num_gamma_channels() <= 4) )
            continue;

          const shared_ptr<const EnergyCalibration> cal = m->energy_calibration();
          if( !cal || !cal->valid() || (cal->type() == EnergyCalType::LowerChannelEdge) )
            continue;

          //If we have already computed the new calibration for a EnergyCalibration object, lets
          //  not re-due it.
          if( old_to_new_cal.count(cal) )
            continue;

          const size_t nchannel = cal->num_channels();
          const vector<float> &coefficients = cal->coefficients();

          auto new_cal = make_shared<EnergyCalibration>();
          switch( cal->type() )
          {
            case EnergyCalType::Polynomial:
            case EnergyCalType::UnspecifiedUsingDefaultPolynomial:
              new_cal->set_polynomial( nchannel, coefficients, new_dev_pairs );
              break;

            case EnergyCalType::FullRangeFraction:
              new_cal->set_full_range_fraction( nchannel, coefficients, new_dev_pairs );
              break;

            case EnergyCalType::LowerChannelEdge:
            case EnergyCalType::InvalidEquationType:
              assert( 0 );
              break;
          }//switch( cal->type() )

          assert( new_cal->valid() );
          old_to_new_cal[cal] = new_cal;
        }//for( const string &name : change.detectors )
      }//for( const int sample : change.sample_numbers )
    }catch( std::exception &e )
    {
      display->updateToGui( old_cal );
      m_interspec->logMessage( WString::tr("ect-dev-pair-change-invalid").arg(e.what()), 2 );
      return;
    }//try / catch

    // Now go through and translate the peaks, but we wont actually update them to the SpecMeas
    //  until we know we can update all the peaks
    const vector<string> detnamesv( begin(change.detectors), end(change.detectors) );
    const set<set<int>> samplesWithPeaks = change.meas->sampleNumsWithPeaks();

    for( const set<int> &samples : samplesWithPeaks )
    {
      try
      {
        auto dispcal = change.meas->suggested_sum_energy_calibration( samples, detnamesv );

        const auto dispcaliter = old_to_new_cal.find(dispcal);
        if( dispcaliter == end(old_to_new_cal) )
          continue;

        auto oldpeaks = change.meas->peaks(samples);
        if( !oldpeaks || oldpeaks->empty() || updated_peaks.count(oldpeaks) )
          continue;

        auto newpeaks = EnergyCal::translatePeaksForCalibrationChange( *oldpeaks, dispcaliter->first,
                                                                       dispcaliter->second );
        updated_peaks[oldpeaks] = newpeaks;
      }catch( std::exception &e )
      {
        display->setDeviationPairsInvalid();

        string msg = "There was an issue translating peaks for this deviation pair change;"
                     " not applying change.  Error: " + string(e.what());
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, msg.c_str() );
#endif
        m_interspec->logMessage( msg, 2 );

        return;
      }//try / catch
    }//for( const set<int> &samples : samplesWithPeaks )


    const set<set<int>> samplesWithHintPeaks = change.meas->sampleNumsWithAutomatedSearchPeaks();
    for( const set<int> &samples : samplesWithHintPeaks )
    {
      try
      {
        auto dispcal = change.meas->suggested_sum_energy_calibration( samples, detnamesv );

        const auto dispcaliter = old_to_new_cal.find(dispcal);
        if( dispcaliter == end(old_to_new_cal) )
          continue;

        auto oldHintPeaks = change.meas->automatedSearchPeaks(samples);
        if( !oldHintPeaks || oldHintPeaks->empty() || updated_hint_peaks.count(oldHintPeaks) )
          continue;

        auto newpeaks = EnergyCal::translatePeaksForCalibrationChange( *oldHintPeaks,
                                                                dispcaliter->first, dispcaliter->second );
        updated_hint_peaks[oldHintPeaks] = newpeaks;
      }catch( std::exception &e )
      {
#if( PERFORM_DEVELOPER_CHECKS )
        string msg = "There was an issue translating hint peaks for this deviation pair change."
        " Error: " + string(e.what());
        log_developer_error( __func__, msg.c_str() );
#endif
      }//try / catch
    }//for( const set<int> &samples : samplesWithHintPeaks )
  }//for( const MeasToApplyCoefChangeTo &change : changemeas )

  display->setDeviationPairsValid();


  size_t num_updated = 0;

  // Track info for undo/redo; changes are recorded per affected SpecMeas (the undo/redo step
  //  itself is associated with the foreground, as all InterSpec undo/redo is).
  //  Note we only record into the sentry in this commit stage, so an early return above doesnt
  //  insert an undo/redo step for changes that were never applied.
  EnergyCalUndoRedoSentry undo_sentry;

  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    assert( change.meas );
    meas_old_new_cal_t &meas_old_new_cal = undo_sentry.cal_info( change.meas );

    for( const int sample : change.sample_numbers )
    {
      for( const string &name : change.detectors )
      {
        const shared_ptr<const Measurement> m = change.meas->measurement( sample, name );
        if( !m || (m->num_gamma_channels() <= 4) )
          continue;

        const shared_ptr<const EnergyCalibration> cal = m->energy_calibration();
        if( !cal || !cal->valid() || (cal->type() == EnergyCalType::LowerChannelEdge) )
          continue;

        const auto calpos = old_to_new_cal.find(cal);
        if( (calpos == end(old_to_new_cal)) || !calpos->second || !calpos->second->valid() )
        {
#if( PERFORM_DEVELOPER_CHECKS )
          log_developer_error( __func__, "Unexpectedly found invalid calibration in old_to_new_cal" );
#endif
          continue;
        }//if( sanity check that new calibration is valid - should always be )

        meas_old_new_cal.emplace_back( m, cal, calpos->second );

        change.meas->set_energy_calibration( calpos->second, m );
        ++num_updated;
      }//for( const string &name : change.detectors )
    }//for( const int sample : change.sample_numbers )
  }//for( const MeasToApplyCoefChangeTo &change : changemeas )

  if( num_updated == 0 )
  {
    display->updateToGui( old_cal );
    m_interspec->logMessage( WString::tr("ect-set-dev-pair-err"), 2 );
    return;
  }

  // Now actually set the new peaks, and hint peaks
  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    meas_old_new_peaks_t &meas_old_new_peaks = undo_sentry.peak_info( change.meas );
    meas_old_new_peaks_t &meas_old_new_hint_peaks = undo_sentry.hint_peak_info( change.meas );

    const set<set<int>> samplesWithPeaks = change.meas->sampleNumsWithPeaks();
    for( const set<int> &samples : samplesWithPeaks )
    {
      auto oldpeaks = change.meas->peaks(samples);
      if( !oldpeaks || oldpeaks->empty() )
        continue;

      const auto peakpos = updated_peaks.find(oldpeaks);
      if( peakpos == end(updated_peaks) )
        continue;

      meas_old_new_peaks.emplace_back( samples, *oldpeaks, peakpos->second );

      change.meas->setPeaks( peakpos->second, samples );
      if( m_peakModel && (change.meas == foregrnd) && (samples == foreSamples) )
        m_peakModel->setPeakFromSpecMeas( foregrnd, foreSamples, SpecUtils::SpectrumType::Foreground );
      else if( m_peakModel && (change.meas == backgrnd) && (samples == backSamples) )
        m_peakModel->setPeakFromSpecMeas( backgrnd, backSamples, SpecUtils::SpectrumType::Background );
      else if( m_peakModel && (change.meas == secogrnd) && (samples == secoSamples) )
        m_peakModel->setPeakFromSpecMeas( secogrnd, secoSamples, SpecUtils::SpectrumType::SecondForeground );
    }//for( const set<int> &samples : samplesWithPeaks )

    const set<set<int>> samplesWithHintPeaks = change.meas->sampleNumsWithAutomatedSearchPeaks();
    for( const set<int> &samples : samplesWithHintPeaks )
    {
      auto oldHintPeaks = change.meas->automatedSearchPeaks(samples);
      if( !oldHintPeaks || oldHintPeaks->empty() )
        continue;

      const auto peakpos = updated_hint_peaks.find(oldHintPeaks);
      if( peakpos == end(updated_hint_peaks) )
        continue;

      meas_old_new_hint_peaks.emplace_back( samples, *oldHintPeaks, peakpos->second );

      auto peaks = make_shared<deque<shared_ptr<const PeakDef>>>( peakpos->second );
      change.meas->setAutomatedSearchPeaks( samples, peaks );
    }//for( const set<int> &samples : samplesWithHintPeaks )
  }//for( const MeasToApplyCoefChangeTo &change : changemeas )

  // Let the user know if the change may have applied beyond just the displayed spectrum
  const size_t ndets = dispmeas->gamma_detector_names().size();
  const size_t nsamples = dispmeas->sample_numbers().size();

  /// \TODO: keep from issuing this message for ever single change!
  if( (ndets > 1) || (nsamples > 1) || (changemeas.size() > 1) )
    m_interspec->logMessage( WString::tr("ect-dev-pairs-applied-scope"), 1 );

  m_interspec->refreshDisplayedCharts();
  refreshGuiFromFiles();
}//void userChangedDeviationPair( CalDisplay *display )


void EnergyCalTool::displayedSpectrumChanged( const SpecUtils::SpectrumType type,
                                              const std::shared_ptr<SpecMeas> &meas,
                                              const std::set<int> &samples,
                                              const std::vector<std::string> &detectors )
{
  static_assert( static_cast<int>(SpecUtils::SpectrumType::Foreground) == 0, "" );
  static_assert( static_cast<int>(SpecUtils::SpectrumType::SecondForeground) == 1, "" );
  static_assert( static_cast<int>(SpecUtils::SpectrumType::Background) == 2, "" );

  const int index = static_cast<int>( type );
  assert( index >= 0 && index < 3 );
  
  if( meas != m_currentSpecMeas[index] )
  {
    //whole new file
    //cout << "EnergyCalTool::displayedSpectrumChanged: new file" << endl;
    
    //We want to cache original energy calibration, if we havent already
  }else if( samples != m_currentSampleNumbers[index] )
  {
    // Just changed what was displayed
    //cout << "EnergyCalTool::displayedSpectrumChanged: changed sample numbers" << endl;
  }else
  {
    //no change...
    //cout << "EnergyCalTool::displayedSpectrumChanged: same file and sample numbers" << endl;
  }

  m_lastGraphicalRecal = 0;
  m_lastGraphicalRecalType = EnergyCalGraphicalConfirm::NumRecalTypes;
  m_lastGraphicalRecalEnergy = -999.0f;
  
  m_currentSpecMeas[index] = meas;
  m_currentSampleNumbers[index] = samples;
  
  refreshGuiFromFiles();
}//void displayedSpectrumChanged(...)


void EnergyCalTool::setShowNoCalInfo( const bool nocal )
{
  m_noCalTxt->setHidden( !nocal );
  m_moreActionsColumn->setHidden( nocal );
  m_applyToColumn->setHidden( nocal );
  m_detColumn->setHidden( nocal );
  m_calColumn->setHidden( nocal );
  m_peakTableColumn->setHidden( nocal );
}//void setShowNoCalInfo( const bool nocal )


void EnergyCalTool::setWasGraphicalRecal( int type, float energy )
{
  time( &m_lastGraphicalRecal );
  m_lastGraphicalRecalType = type;
  m_lastGraphicalRecalEnergy = energy;
}//void setWasGraphicalRecal( int type, double energy )


void EnergyCalTool::specTypeToDisplayForChanged()
{
  const int selectedType = m_specTypeMenu->currentIndex();
  if( selectedType < 0 || selectedType > 2 )
  {
    assert( 0 );  //shouldnt ever happen, right?
    return;
  }
  
  WMenu *detMenu = m_detectorMenu[selectedType];
  WMenuItem *detitem = detMenu->currentItem();
  if( !detitem && detMenu->count() )
    detMenu->select( 0 );
  if( detitem )
    detMenu->select( detitem );
  
  updateFitButtonStatus();
}//void specTypeToDisplayForChanged();


bool EnergyCalTool::canDoEnergyFit()
{
  // Check that "Apply Changes To" for the "Foreground" is checked, otherwise fitting makes no sense
  if( !m_applyToCbs[ApplyToCbIndex::ApplyToForeground]->isChecked() )
    return false;
  
  // Check if there are any peaks currently showing.
  shared_ptr<const deque<PeakModel::PeakShrdPtr>> peaks = m_peakModel->peaks();
  if( !peaks )
    return false;
  
  size_t nPeaksToUse = 0;
  for( const PeakModel::PeakShrdPtr &p : *peaks )
    nPeaksToUse += (p && p->useForEnergyCalibration());
  
  if( nPeaksToUse < 1 )
    return false;
  
  if( !m_calInfoDisplayStack )
    return false;
  
  // We are actually going to fit the coefficients for the currently showing CalDisplay, so only
  //  consult the checkboxes on that one display.
  auto caldisp = dynamic_cast<EnergyCalImp::CalDisplay *>( m_calInfoDisplayStack->currentWidget() );
  if( !caldisp )
    return false;
  
  auto cal = caldisp->lastSetCalibration();
  if( !cal || !cal->valid() )
    return false;

  // The peaks used for the fit are the foregrounds peaks, so the calibration being shown must
  //  belong to the foreground - or at least be the foregrounds displayed calibration (e.g., when
  //  the same file is loaded as multiple spectrum types).
  if( caldisp->spectrumType() != SpecUtils::SpectrumType::Foreground )
  {
    const shared_ptr<const SpecUtils::Measurement> forehist
                           = m_interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
    if( !forehist || (forehist->energy_calibration() != cal) )
      return false;
  }//if( the display being shown isnt for the foreground )

  switch( cal->type() )
  {
    case SpecUtils::EnergyCalType::InvalidEquationType:
      return false;

    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::FullRangeFraction:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
    case SpecUtils::EnergyCalType::LowerChannelEdge:  //fits {offset, gain} relative to original
      break;
  }//switch( cal->type() )

  const set<size_t> ordersToFit = caldisp->fitForCoefficents();
  if( ordersToFit.empty() || ordersToFit.size() > nPeaksToUse )
    return false;

  return true;
}//bool canDoEnergyFit()


void EnergyCalTool::fitCoefficients()
{
  try
  {
    if( !canDoEnergyFit() )
    {
      m_interspec->logMessage( WString::tr("ect-err-not-enough-peaks"), 2 );
      return;
    }//if( double check we can actually do the fit )
    
    
    // Check if there are any peaks currently showing.
    shared_ptr<const deque<PeakModel::PeakShrdPtr>> peaks = m_peakModel->peaks();
    if( !peaks || peaks->empty() )  //shouldnt ever happen.
      throw runtime_error( WString::tr("ect-no-peaks").toUTF8() );
    
    // We are actually going to fit the coefficients for the currently showing CalDisplay, so only
    //  consult the checkboxes on that one display.
    auto caldisp = dynamic_cast<EnergyCalImp::CalDisplay *>( m_calInfoDisplayStack->currentWidget() );
    if( !caldisp )  //shouldnt ever happen, but JIC
      throw runtime_error( "Unexpected error determining current calibration" );
    
    // The peaks we fit to are always the foregrounds peaks, so make sure the display we take the
    //  calibration from (and that gets updated first) is a foreground display of a displayed
    //  detector; #canDoEnergyFit has already checked the currently showing display is consistent
    //  with the foregrounds calibration.
    string detname = caldisp->detectorName();
    const vector<string> displayed = m_interspec->detectorsToDisplay( SpecUtils::SpectrumType::Foreground );

    if( (caldisp->spectrumType() != SpecUtils::SpectrumType::Foreground)
        || (std::find(begin(displayed), end(displayed), detname) == end(displayed)) )
    {
      EnergyCalImp::CalDisplay *foredisp = nullptr;

      for( int menunum = 0; !foredisp && (menunum < 3); ++menunum )
      {
        WMenu * const detMenu = m_detectorMenu[menunum];
        for( int i = 0; !foredisp && detMenu && (i < detMenu->count()); ++i )
        {
          WMenuItem * const item = detMenu->itemAt(i);
          auto display = item ? dynamic_cast<EnergyCalImp::CalDisplay *>( item->contents() ) : nullptr;
          if( display
              && (display->spectrumType() == SpecUtils::SpectrumType::Foreground)
              && (std::find(begin(displayed), end(displayed), display->detectorName()) != end(displayed)) )
          {
            foredisp = display;
            // Select the found display, so the user sees the calibration that got fit
            m_specTypeMenu->select( menunum );
            item->select();
          }//if( we found a displayed foreground detector )
        }//for( loop over the menus items )
      }//for( loop over the detector menus )

      if( !foredisp )
        throw runtime_error( WString::tr("ect-select-cal-of-disp-det").toUTF8() );

      caldisp = foredisp;
      detname = caldisp->detectorName();
    }//if( the currently showing display isnt a displayed foreground detector )


    auto original_cal = caldisp->lastSetCalibration();

    switch( original_cal->type() )
    {
      case SpecUtils::EnergyCalType::InvalidEquationType:
        throw runtime_error( "Unexpected calibration type from display" ); //shouldnt ever happen, but JIC

      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::FullRangeFraction:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      case SpecUtils::EnergyCalType::LowerChannelEdge:
        break;
    }//switch( cal->type() )
    
    if( original_cal->num_channels() < 5 )
      throw runtime_error( WString::tr("ect-not-enough-channel").toUTF8() );
    
    
    const set<size_t> orders_to_fit = caldisp->fitForCoefficents();
    if( orders_to_fit.empty() )  //shouldnt ever happen
      throw runtime_error( WString::tr("ect-no-coeff-selected").toUTF8() );
    
    
    //TODO: meansFitError will currently contain only values of 1.0, eventually
    //      will contain the error of the fit mean for that peak
    vector<EnergyCal::RecalPeakInfo> peakInfos;
    
    for( const auto &peakptr : *peaks )
    {
      if( !peakptr )  //shouldnt be necassary, but JIC
        continue;
      const PeakDef &peak = *peakptr;
      
      if( !peak.useForEnergyCalibration() )
        continue;
      
      const double wantedEnergy = peak.gammaParticleEnergy();
      
      EnergyCal::RecalPeakInfo peakInfo;
      peakInfo.peakMean = peak.mean();
      // Clamp to peak mean uncertainty to be at least 0.25 keV.  This is an arbitrary decision, but
      //  motivated by not wanting a single peak to way, way, dominate the other peaks, when it is
      //  likely non-linearities in the detector may actually dominate the effects
      peakInfo.peakMeanUncert = max( peak.meanUncert(), 0.25 );
      if( IsInf(peakInfo.peakMeanUncert) || IsNan(peakInfo.peakMeanUncert) )
        peakInfo.peakMeanUncert = 0.5;
      
      peakInfo.photopeakEnergy = wantedEnergy;
      peakInfo.peakMeanBinNumber = original_cal->channel_for_energy( peak.mean() );
      
      peakInfos.push_back( peakInfo );
    }//for( int col = 0; col < numModelCol; ++col )
    
    if( orders_to_fit.size() > peakInfos.size() )
      throw runtime_error( WString::tr("ect-err-not-enough-peaks").toUTF8() );
    
    shared_ptr<const SpecUtils::EnergyCalibration> answer;
    double chi2 = -999;

    if( original_cal->type() == SpecUtils::EnergyCalType::LowerChannelEdge )
    {
      // Fit the {offset, gain} adjustment relative to the sessions original channel energies;
      //  this problem is exactly linear in the fit parameters, so no iterative fallback needed.
      const LowerChanCalOriginal previnfo = lowerChannelOriginal( original_cal );

      vector<bool> fitfor( 2, false );
      for( const size_t order : orders_to_fit )
      {
        if( order > 1 )
          throw runtime_error( "Unexpected coefficient order for lower channel energy cal" );
        fitfor[order] = true;
      }

      vector<float> coefs = { static_cast<float>(previnfo.offset),
                              static_cast<float>(previnfo.gain) };
      vector<float> coef_uncerts;
      chi2 = EnergyCal::fit_energy_cal_lower_channel( peakInfos, previnfo.original, fitfor,
                                                      coefs, coef_uncerts );

      LowerChanCalOriginal info = previnfo;
      info.offset = coefs[0];
      info.gain = coefs[1];

      const auto adjusted = EnergyCal::adjust_lower_channel_energy_cal( info.original,
                                                                        info.offset, info.gain );
      registerLowerChannelAdjustment( adjusted, info );
      answer = adjusted;
    }else
    {
      auto fit_cal = make_shared<SpecUtils::EnergyCalibration>();

      const size_t eqn_order = std::max( original_cal->coefficients().size(), (*orders_to_fit.rbegin()) + 1 );
      const size_t nchannel = original_cal->num_channels();
      const auto &devpairs = original_cal->deviation_pairs();

      try
      {
        vector<bool> fitfor( eqn_order, false );

        for( auto order : orders_to_fit )
        {
          assert( order < fitfor.size() );
          fitfor[order] = true;
        }

        vector<float> coefficent_uncerts;
        vector<float> coefficents = original_cal->coefficients();
        if( coefficents.size() < eqn_order )
          coefficents.resize( eqn_order, 0.0f );

        switch ( original_cal->type() )
        {
          case SpecUtils::EnergyCalType::Polynomial:
          case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
            chi2 = EnergyCal::fit_energy_cal_poly( peakInfos, fitfor, nchannel, devpairs,
                                                  coefficents, coefficent_uncerts );
            fit_cal->set_polynomial( nchannel, coefficents, devpairs );
            break;

          case SpecUtils::EnergyCalType::FullRangeFraction:
            chi2 = EnergyCal::fit_energy_cal_frf( peakInfos, fitfor, nchannel, devpairs,
                                                 coefficents, coefficent_uncerts );
            fit_cal->set_full_range_fraction( nchannel, coefficents, devpairs );
            break;

          case SpecUtils::EnergyCalType::LowerChannelEdge:
          case SpecUtils::EnergyCalType::InvalidEquationType:
            assert( 0 );  //lower channel energy handled in its own branch above
            break;
        }//switch ( original_cal->type() )

      }catch( std::exception &e )
      {
#if( PERFORM_DEVELOPER_CHECKS )
        char buffer[512] = { '\0' };
        snprintf( buffer, sizeof(buffer)-1, "fit_energy_cal_poly threw: %s", e.what() );
        log_developer_error( __func__, buffer );
#endif
      }//try / catch fit for coefficents using least linear squares


      if( !fit_cal->valid() )
      {
        // The linear least squares fit failed; fall back to a Ceres based non-linear fit
        EnergyCal::EnergyCalCeresFitSetup ceres_setup;
        ceres_setup.cal_type = original_cal->type();
        ceres_setup.num_channels = nchannel;
        ceres_setup.fitfor = vector<bool>( eqn_order, false );
        for( auto order : orders_to_fit )
          ceres_setup.fitfor[order] = true;
        ceres_setup.starting_coefs = original_cal->coefficients();
        if( ceres_setup.starting_coefs.size() < eqn_order )
          ceres_setup.starting_coefs.resize( eqn_order, 0.0f );
        ceres_setup.dev_pairs = devpairs;

        const EnergyCal::EnergyCalCeresFitResult ceres_result
                                      = EnergyCal::fit_energy_cal_ceres( peakInfos, ceres_setup );
        chi2 = ceres_result.chi2;

        if( !ceres_result.warning_msg.empty() )
          m_interspec->logMessage( WString::fromUTF8(ceres_result.warning_msg), 3 );

        switch ( original_cal->type() )
        {
          case SpecUtils::EnergyCalType::Polynomial:
          case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
            fit_cal->set_polynomial( nchannel, ceres_result.coefs, devpairs );
            break;

          case SpecUtils::EnergyCalType::FullRangeFraction:
            fit_cal->set_full_range_fraction( nchannel, ceres_result.coefs, devpairs );
            break;

          case SpecUtils::EnergyCalType::LowerChannelEdge:
          case SpecUtils::EnergyCalType::InvalidEquationType:
            assert( 0 );
            break;
        }//switch ( original_cal->type() )
      }//if( !fit_coefs )

      answer = fit_cal;
    }//if( lower channel energy cal ) / else

    if( !answer || !answer->valid() )
      throw runtime_error( WString::tr("ect-fail-min").toUTF8() );

    const vector<MeasToApplyCoefChangeTo> changemeas = measurementsToApplyCoeffChangeTo();
    applyCalChange( original_cal, answer, changemeas, false );
    
    
    //To show Chi2 in the message, uncomment out this next section
    /*
    if( peakInfos.size() > orders_to_fit.size() )
    {
      double dof = peakInfos.size() - 1;
      dof -= orders_to_fit.size();
      dof = (dof < 1) ? 1.0 : dof;
      
      char buffer[64];
      snprintf( buffer, sizeof(buffer), " &chi;&sup2;/dof=%.2g", (chi2/dof) );
      msg += buffer;
    }
    */
    
    m_interspec->logMessage( WString::tr("ect-fit-successful"), 1 );
  }catch( std::exception &e )
  {
    WString msg = WString::tr("ect-fail-fit").arg( e.what() );
    cerr << "EnergyCalTool::fitCoefficients():\n\tCaught: " << msg.toUTF8() << endl;
    m_interspec->logMessage( msg, 3 );
  }//try / catch
}//void fitCoefficients()


void EnergyCalTool::updateFitButtonStatus()
{
  const bool canFit = canDoEnergyFit();
  
#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
  for( EnergyCalImp::CalDisplay *disp : calDisplays() )
    disp->setFitButtonEnabled( canFit );
#else
  if( canFit != m_fitCalBtn->isEnabled() )
    m_fitCalBtn->setDisabled( !canFit );
#endif
}//void updateFitButtonStatus()


#if( !IMP_CALp_BTN_NEAR_COEFS )
void EnergyCalTool::updateCALpButtonsStatus()
{
  shared_ptr<const SpecUtils::Measurement> foregrnd
                           = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  
  const bool showDownload = (foregrnd && foregrnd->energy_calibration()
                              && foregrnd->energy_calibration()->valid());
  m_downloadCALp->setHidden( !showDownload );
  
  const bool showUpload = (foregrnd && foregrnd->num_gamma_channels());
  m_uploadCALp->setHidden( !showUpload );
}//void updateCALpButtonsStatus()
#endif


void EnergyCalTool::displayedSpecChangedCallback( const SpecUtils::SpectrumType,
                                                  const std::shared_ptr<SpecMeas>,
                                                  const std::set<int>,
                                                  const std::vector<std::string> )
{
  /// \TODO: set the various m_applyToCbs if it is a new spectrum being shown.
  /// \TODO: if this is the first time seeing a SpecMeas, cache all of its energy calibration
  ///        information
  
  // \TODO: we could maybe save a little time by inspecting what was changed, but the added
  //        complexity probably isnt worth it, so we'll skip this.
  refreshGuiFromFiles();
}//void displayedSpecChangedCallback(...)


void EnergyCalTool::refreshGuiFromFiles()
{
  m_renderFlags |= EnergyCalToolRenderFlags::FullGuiUpdate;
  scheduleRender();
}//void refreshGuiFromFiles()


void EnergyCalTool::handleGraphicalRecalRequest( double xstart, double xfinish )
{
  try
  {
    auto foreground = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    shared_ptr<const SpecUtils::EnergyCalibration> energycal = foreground
                                                               ? foreground->energy_calibration()
                                                               : nullptr;
    if( !energycal || !energycal->valid()
        || (energycal->type() == SpecUtils::EnergyCalType::InvalidEquationType) )
      return;
  
    if( m_graphicalRecal )
    {
      m_graphicalRecal->setEnergies( xstart, xfinish );
    }else
    {
      m_graphicalRecal = new EnergyCalGraphicalConfirm( xstart, xfinish, this,
                                    m_lastGraphicalRecal,
                                    EnergyCalGraphicalConfirm::RecalTypes(m_lastGraphicalRecalType),
                                    m_lastGraphicalRecalEnergy );
    
      m_graphicalRecal->finished().connect( this, &EnergyCalTool::deleteGraphicalRecalConfirmWindow );
    }
  }catch( std::runtime_error & )
  {
    m_interspec->logMessage( "Internal error doing graphical recal; sorry :(", 3 );
  }
}//void handleGraphicalRecalRequest( double xstart, double xfinish )


void EnergyCalTool::deleteGraphicalRecalConfirmWindow()
{
  if( m_graphicalRecal )
  {
    AuxWindow::deleteAuxWindow( m_graphicalRecal );
    m_graphicalRecal = nullptr;
  }//if( m_graphicalRecal )
  
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", m_interspec );
  if( showToolTips )
  {
    m_interspec->logMessage( WString::tr("ect-del-graphical-msg"), 1 );
  }//if( showToolTips )
}//void deleteGraphicalRecalConfirmWindow()


string EnergyCalTool::applyToSummaryTxt() const
{
  string answer;
  
  if( m_applyToColumn->isHidden() )
    return answer;
  
  for( ApplyToCbIndex index = static_cast<ApplyToCbIndex>(0);
      index < ApplyToCbIndex::NumApplyToCbIndex;
      index = ApplyToCbIndex(index+1) )
  {
    Wt::WCheckBox *cb = m_applyToCbs[index];
    const auto cbparent = cb->parent();
    assert( cbparent );
    assert( dynamic_cast<WContainerWidget *>(cbparent) );
    
    if( cbparent->isHidden() || !cb->isChecked() )
      continue;
    
    if( !answer.empty() )
      answer += ", ";
    
    switch( index )
    {
      case ApplyToCbIndex::ApplyToForeground:         answer += WString::tr("foreground").toUTF8(); break;
      case ApplyToCbIndex::ApplyToBackground:         answer += WString::tr("background").toUTF8(); break;
      case ApplyToCbIndex::ApplyToSecondary:          answer += WString::tr("secondary").toUTF8();  break;
      case ApplyToCbIndex::ApplyToDisplayedDetectors: answer += WString::tr("ect-lc-disp-dets").toUTF8(); break;
      case ApplyToCbIndex::ApplyToAllDetectors:       answer += WString::tr("ect-lc-all-dets").toUTF8(); break;
      case ApplyToCbIndex::ApplyToDisplayedSamples:   answer += WString::tr("ect-lc-disp-samples").toUTF8(); break;
      case ApplyToCbIndex::ApplyToAllSamples:         answer += WString::tr("ect-lc-all-samples").toUTF8(); break;
      case ApplyToCbIndex::NumApplyToCbIndex:
        break;
    }//switch( index )
  }//for( loop over ApplyToCbIndex )
  
  return answer;
}//string applyToSummaryTxt() const


namespace
{
  /** The info needed to create, or update, one CalDisplay entry in the detector menus; see
   #EnergyCalTool::doRefreshFromFiles.
   */
  struct CalDisplayEntryInfo
  {
    SpecUtils::SpectrumType type;
    std::string detname;
    Wt::WString label;
    std::shared_ptr<const SpecUtils::EnergyCalibration> cal;
  };//struct CalDisplayEntryInfo
}//namespace


void EnergyCalTool::createCalDisplayWidgets()
{
  assert( !m_specTypeMenu && !m_specTypeMenuStack && !m_calInfoDisplayStack );
  for( int i = 0; i < 3; ++i )
  {
    assert( !m_detectorMenu[i] );
  }

  //Labels for the vertical menu when you have multiple spectra shown, and at least one of them
  //  has more than one detector.
  const char * const spec_type_labels[3] = {"ect-short-fore","ect-short-back","ect-short-secondary"};

  // Note: we deliberately dont set a transition animation on either of the stacks: in Wt 3.7.1 an
  //  animated WStackedWidget::setCurrentIndex() only animates the previous/next widgets, and does
  //  not reconcile the hidden state of its other children, so adding/removing children while
  //  re-using the stack could leave the newly current widget invisible (this was the cause of the
  //  "calibration coefficients wont show up" issue that previously forced re-creating all of
  //  these widgets on nearly every update).
  m_specTypeMenuStack = new WStackedWidget();
  m_specTypeMenuStack->addStyleClass( "CalSpecStack" );

  m_specTypeMenu = new WMenu( m_specTypeMenuStack );
  m_specTypeMenu->addStyleClass( "CalSpecMenu" );
  m_specTypeMenu->itemSelected().connect( this, &EnergyCalTool::specTypeToDisplayForChanged );
  m_detColLayout->addWidget( m_specTypeMenu, 1, 0 );
  m_detColLayout->addWidget( m_specTypeMenuStack, 2, 0 );

  WGridLayout * const callayout = dynamic_cast<WGridLayout *>( m_calColumn->layout() );
  assert( callayout );

  m_calInfoDisplayStack = new WStackedWidget();
  m_calInfoDisplayStack->addStyleClass( "ToolTabTitledColumnContent CalStack" );
  callayout->addWidget( m_calInfoDisplayStack, 1, 1 );

  for( int i = 0; i < 3; ++i )
  {
    WContainerWidget *detMenuDiv = new WContainerWidget();  //this holds the WMenu for this SpecFile
    detMenuDiv->addStyleClass( "DetMenuDiv" );

    WMenuItem *item = m_specTypeMenu->addItem( WString::tr(spec_type_labels[i]), detMenuDiv,
                                               WMenuItem::LoadPolicy::PreLoading );
    //Fix issue, for Wt 3.3.4 at least, if user doesnt click exactly on the <a> element
    item->clicked().connect( boost::bind(&WMenuItem::select, item) );
    item->setHidden( true );  //#doRefreshFromFiles will un-hide, as needed

    m_detectorMenu[i] = new WMenu( m_calInfoDisplayStack, detMenuDiv );
    m_detectorMenu[i]->addStyleClass( "VerticalNavMenu HeavyNavMenu DetCalMenu" );
    m_detectorMenu[i]->itemSelected().connect( this, &EnergyCalTool::updateFitButtonStatus );
  }//for( int i = 0; i < 3; ++i )
}//void createCalDisplayWidgets()


void EnergyCalTool::doRefreshFromFiles()
{
  //Menu entry labels used when each displayed file has no more than a single detector
  const char * const spec_type_labels_vert[3] = {"Foreground", "Background", "Secondary"};

  const SpecUtils::SpectrumType spectypes[3] = {
    SpecUtils::SpectrumType::Foreground,
    SpecUtils::SpectrumType::Background,
    SpecUtils::SpectrumType::SecondForeground
  };

  shared_ptr<const SpecMeas> specfiles[3];
  for( int i = 0; i < 3; ++i )
    specfiles[i] = m_interspec->measurment( spectypes[i] );

  // Prune expired entries from the lower-channel-energy original-calibration tracking (the
  //  calibrations die when their file is unloaded and the undo/redo history lets go of them;
  //  the adjusted calibration then just becomes the new nominal on any future load, which is
  //  fine)
  for( auto iter = begin(m_lowerChanOrigCals); iter != end(m_lowerChanOrigCals); )
  {
    if( iter->first.expired() )
      iter = m_lowerChanOrigCals.erase( iter );
    else
      ++iter;
  }

  const bool isWide = (m_tallLayoutContent ? false : true);

  // Get the gamma detector names, for the displayed sample numbers, of each spectrum type
  set<string> disp_det_names[3];

  //If each spectrum file has no more than a single displayed detector, then instead of having the
  //  vertical menu display detector names, we will have "Foreground", "Background", "Secondary"
  //  entries, all in the foreground detector menu.
  bool specTypeInForgrndMenu = true;

  int nFilesWithCalInfo = 0;

  for( int i = 0; i < 3; ++i )
  {
    disp_det_names[i] = gammaDetectorsForDisplayedSamples( spectypes[i] );
    specTypeInForgrndMenu = (specTypeInForgrndMenu && (disp_det_names[i].size() <= 1));

    if( !disp_det_names[i].empty() )
      nFilesWithCalInfo += 1;
  }//for( int i = 0; i < 3; ++i )

  // The menus and stacks are created just once for each wide/tall layout (#initWidgets nulls the
  //  pointers on a layout change); after that we only diff their contents against the currently
  //  displayed spectra, updating the CalDisplay widgets in-place, and creating/deleting only the
  //  widgets we need to.  This avoids the jitter, and losing the users cursor/focus, that
  //  re-creating everything on each calibration edit used to cause.
  if( !m_specTypeMenu )
    createCalDisplayWidgets();

  assert( m_specTypeMenu && m_specTypeMenuStack && m_calInfoDisplayStack );

  //Dont show spectype menu (the vertical "For.", "Back", "Sec." menu), if we dont need to
  const bool hideSpecType = ( specTypeInForgrndMenu
                             || (nFilesWithCalInfo < 2)
                             || ( (!specfiles[1] || (specfiles[0]==specfiles[1]))
                                 && (!specfiles[2] || (specfiles[0]==specfiles[2]))) );

  // Figure out the CalDisplay entries each detector menu should end up with.  For each detector
  //  we use the energy calibration of the first displayed sample that has gamma data.  If there
  //  is no foreground, all the menus should become empty.
  bool hasFRFCal = false, hasPolyCal = false, hasLowerChanCal = false;

  vector<CalDisplayEntryInfo> wanted_entries[3];  //indexed by detector menu number
  set<string> specdetnames[3]; //Just the names of gamma detectors with at least 5 channels

  for( int i = 0; specfiles[0] && (i < 3); ++i )
  {
    if( !specfiles[i] || disp_det_names[i].empty() )
      continue;

    specdetnames[i] = disp_det_names[i];

    const SpecUtils::SpectrumType type = spectypes[i];
    const set<int> &samples = m_interspec->displayedSamples(type);
    const int menu_index = (specTypeInForgrndMenu ? 0 : i);

    for( const string &detname : disp_det_names[i] )
    {
      for( const int sample : samples )
      {
        const shared_ptr<const SpecUtils::Measurement> m = specfiles[i]->measurement( sample, detname );
        if( !m || (m->num_gamma_channels() <= 4) )
          continue;

        const shared_ptr<const SpecUtils::EnergyCalibration> energycal = m->energy_calibration();

        CalDisplayEntryInfo info;
        info.type = type;
        info.detname = detname;
        info.label = specTypeInForgrndMenu ? WString::tr(spec_type_labels_vert[i])
                                           : WString::fromUTF8(detname);
        info.cal = energycal;
        wanted_entries[menu_index].push_back( std::move(info) );

        if( energycal )
        {
          switch( energycal->type() )
          {
            case SpecUtils::EnergyCalType::Polynomial:
            case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
              hasPolyCal = true;
              break;

            case SpecUtils::EnergyCalType::FullRangeFraction:
              hasFRFCal = true;
              break;

            case SpecUtils::EnergyCalType::LowerChannelEdge:
              hasLowerChanCal = true;
              break;

            case SpecUtils::EnergyCalType::InvalidEquationType:
              break;
          }//switch( energycal->type() )
        }//if( energycal )

        break;
      }//for( const int sample : samples )
    }//for( const string &detname : disp_det_names[i] )
  }//for( int i = 0; i < 3; ++i )

  // Remember which CalDisplay is currently showing, so we can re-show it (or its re-created
  //  equivalent) after the update.
  bool haveCurrent = false;
  SpecUtils::SpectrumType current_type = SpecUtils::SpectrumType::Foreground;
  string current_det;

  {//begin code-block to get the currently showing CalDisplay
    auto current = dynamic_cast<EnergyCalImp::CalDisplay *>( m_calInfoDisplayStack->currentWidget() );
    if( current )
    {
      haveCurrent = true;
      current_type = current->spectrumType();
      current_det = current->detectorName();
    }
  }//end code-block to get the currently showing CalDisplay

  // We want to preserve which "Fit" check boxes are set for any CalDisplay we delete, but will
  //  then re-create in another menu (i.e., when `specTypeInForgrndMenu` flipped because the
  //  displayed files/detectors changed).
  map< pair<SpecUtils::SpectrumType,string>, set<size_t> > set_fit_for_cbs;

  // First pass: remove the CalDisplays that are no longer wanted in their menu.
  //  Note: with the PreLoading policy a menu items contents (the CalDisplay) is owned by the
  //  (shared) contents stack, not the menu item, so we have to delete the contents explicitly,
  //  after removing the item.
  for( int menu_index = 0; menu_index < 3; ++menu_index )
  {
    WMenu * const detMenu = m_detectorMenu[menu_index];
    assert( detMenu );
    if( !detMenu )
      continue;

    for( int j = detMenu->count() - 1; j >= 0; --j )
    {
      WMenuItem * const item = detMenu->itemAt(j);
      auto display = item ? dynamic_cast<EnergyCalImp::CalDisplay *>( item->contents() ) : nullptr;
      assert( display );
      if( !display )
        continue;

      const pair<SpecUtils::SpectrumType,string> key{ display->spectrumType(), display->detectorName() };

      bool iswanted = false;
      for( size_t k = 0; !iswanted && (k < wanted_entries[menu_index].size()); ++k )
      {
        const CalDisplayEntryInfo &info = wanted_entries[menu_index][k];
        iswanted = ((info.type == key.first) && (info.detname == key.second));
      }

      if( iswanted )
        continue;

      for( int other = 0; other < 3; ++other )
      {
        for( const CalDisplayEntryInfo &info : wanted_entries[other] )
        {
          if( (other != menu_index) && (info.type == key.first) && (info.detname == key.second) )
            set_fit_for_cbs[key] = display->fitForCoefficents();
        }
      }//for( loop over the other menus )

      detMenu->removeItem( item );
      delete item;
      delete display;
    }//for( loop backwards over the menus items )
  }//for( int menu_index = 0; menu_index < 3; ++menu_index )

  // Second pass: update the surviving CalDisplays in-place, and create the missing ones.
  for( int menu_index = 0; menu_index < 3; ++menu_index )
  {
    WMenu * const detMenu = m_detectorMenu[menu_index];
    if( !detMenu )
      continue;

    const vector<CalDisplayEntryInfo> &wanted = wanted_entries[menu_index];

    for( size_t j = 0; j < wanted.size(); ++j )
    {
      const CalDisplayEntryInfo &info = wanted[j];

      WMenuItem *item = nullptr;
      for( int k = 0; !item && (k < detMenu->count()); ++k )
      {
        WMenuItem * const kitem = detMenu->itemAt(k);
        auto display = kitem ? dynamic_cast<EnergyCalImp::CalDisplay *>( kitem->contents() ) : nullptr;
        if( display && (display->spectrumType() == info.type)
            && (display->detectorName() == info.detname) )
          item = kitem;
      }//for( look for an existing menu item for this entry )

      if( item )
      {
        auto display = dynamic_cast<EnergyCalImp::CalDisplay *>( item->contents() );
        assert( display );

        if( item->text() != info.label )
          item->setText( info.label );

        if( display )
          display->updateToGui( info.cal );
      }else
      {
        auto display = new EnergyCalImp::CalDisplay( this, info.type, info.detname, isWide );
        item = detMenu->insertItem( static_cast<int>(j), info.label, display,
                                    WMenuItem::LoadPolicy::PreLoading );
        //Fix issue, for Wt 3.3.4 at least, if user doesnt click exactly on the <a> element
        item->clicked().connect( boost::bind(&WMenuItem::select, item) );

#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
        display->doFitCoeffs().connect( this, &EnergyCalTool::fitCoefficients );
#endif

        display->updateToGui( info.cal );

        const auto fitfor_iter = set_fit_for_cbs.find( {info.type, info.detname} );
        if( fitfor_iter != end(set_fit_for_cbs) )
          display->setFitFor( fitfor_iter->second );
      }//if( item ) / else
    }//for( size_t j = 0; j < wanted.size(); ++j )

#if( PERFORM_DEVELOPER_CHECKS )
    // Check the menu ended up with exactly the wanted entries, in the wanted order (the relative
    //  order of surviving entries always matches, since both old and new orderings derive from
    //  the same sorted detector name sets).
    assert( detMenu->count() == static_cast<int>(wanted.size()) );
    for( int j = 0; j < detMenu->count(); ++j )
    {
      auto display = dynamic_cast<EnergyCalImp::CalDisplay *>( detMenu->itemAt(j)->contents() );
      assert( display );
      assert( !display || (display->spectrumType() == wanted[j].type) );
      assert( !display || (display->detectorName() == wanted[j].detname) );
    }
#endif
  }//for( int menu_index = 0; menu_index < 3; ++menu_index )

#if( !IMP_CALp_BTN_NEAR_COEFS )
  updateCALpButtonsStatus();
#endif

  if( !specfiles[0] )
  {
    setShowNoCalInfo( true );
#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
    for( EnergyCalImp::CalDisplay *disp : calDisplays() )
      disp->setFitButtonEnabled( false );
#else
    m_fitCalBtn->disable();
#endif
    return;
  }//if( !specfiles[0] )

  // Update which spectrum type ("For.", "Back", "Sec.") menu items are showing
  for( int i = 0; i < 3; ++i )
  {
    WMenuItem * const specItem = m_specTypeMenu->itemAt(i);
    assert( specItem );
    if( specItem )
      specItem->setHidden( disp_det_names[i].empty() );
  }//for( int i = 0; i < 3; ++i )

  // Re-select the previously showing CalDisplay (or its re-created equivalent), or else fallback
  //  to the first available display.  We re-select even if the same item still looks current,
  //  since this also re-syncs the (shared) contents stack (e.g., adding the first item to a menu
  //  switches the stack to that widget).
  int select_menu_index = -1;
  WMenuItem *select_item = nullptr;

  for( int i = 0; haveCurrent && !select_item && (i < 3); ++i )
  {
    WMenu * const detMenu = m_detectorMenu[i];
    for( int j = 0; detMenu && !select_item && (j < detMenu->count()); ++j )
    {
      WMenuItem * const item = detMenu->itemAt(j);
      auto display = item ? dynamic_cast<EnergyCalImp::CalDisplay *>( item->contents() ) : nullptr;
      if( display && (display->spectrumType() == current_type)
          && (display->detectorName() == current_det) )
      {
        select_menu_index = i;
        select_item = item;
      }
    }//for( loop over the menus items )
  }//for( loop over detector menus, looking for the previously showing display )

  for( int i = 0; !select_item && (i < 3); ++i )
  {
    if( m_detectorMenu[i] && m_detectorMenu[i]->count() )
    {
      select_menu_index = i;
      select_item = m_detectorMenu[i]->itemAt(0);
    }
  }//for( fallback to the first available display )

  if( select_item )
  {
    assert( (select_menu_index >= 0) && (select_menu_index < 3) );
    m_specTypeMenu->select( select_menu_index );
    m_detectorMenu[select_menu_index]->select( m_detectorMenu[select_menu_index]->indexOf(select_item) );

#if( PERFORM_DEVELOPER_CHECKS )
    assert( m_calInfoDisplayStack->currentWidget() == select_item->contents() );
#endif
  }//if( select_item )

  setShowNoCalInfo( !nFilesWithCalInfo );
  m_specTypeMenu->setHidden( hideSpecType );
  
  bool anyApplyToCbShown = false;
  for( ApplyToCbIndex index = static_cast<ApplyToCbIndex>(0);
      index < ApplyToCbIndex::NumApplyToCbIndex;
      index = ApplyToCbIndex(index+1) )
  {
    Wt::WCheckBox *cb = m_applyToCbs[index];
    const auto cbparent = cb->parent();
    assert( cbparent );
    assert( dynamic_cast<WContainerWidget *>(cbparent) );
    
    bool hideRow = false;
    switch( index )
    {
      case ApplyToCbIndex::ApplyToForeground:
        hideRow = (!specfiles[0] || (!specfiles[1] && !specfiles[2]));
        break;
        
      case ApplyToCbIndex::ApplyToBackground:
        // We'll ignore case where foreground and background is same {SpecFile, Samples, Detectors}
        hideRow = !specfiles[1];
        break;
        
      case ApplyToCbIndex::ApplyToSecondary:
        // We'll ignore case where secondary and background is same {SpecFile, Samples, Detectors}
        hideRow = !specfiles[2];
        break;
        
      case ApplyToCbIndex::ApplyToDisplayedDetectors:
      {
        bool displayingAll = true;
        for( int i = 0; displayingAll && (i < 3); ++i )
        {
          const auto meas = specfiles[i];
          if( !meas )
            continue;
          
          const auto type = spectypes[i];
          const vector<string> displayed = m_interspec->detectorsToDisplay( type );
          for( const auto &name : meas->gamma_detector_names() )
          {
            if( std::find( begin(displayed), end(displayed), name ) == end(displayed) )
              displayingAll = false;
          }//for( const auto &name : meas->gamma_detector_names() )
        }//for( loop over the types of spectrum files )
        
        hideRow = displayingAll;
        if( displayingAll )
          cb->setChecked( false );
        
        break;
      }//case ApplyToCbIndex::ApplyToDisplayedDetectors:
        
      case ApplyToCbIndex::ApplyToAllDetectors:
      {
        static_assert( ApplyToCbIndex::ApplyToDisplayedDetectors
                       < ApplyToCbIndex::ApplyToAllDetectors, "" );
        
        const WCheckBox *dispDetsCb = m_applyToCbs[ApplyToCbIndex::ApplyToDisplayedDetectors];
        assert( dispDetsCb->parent() );
        
        const bool dispDetsHid = dispDetsCb->parent()->isHidden();
        
        hideRow = dispDetsHid;
        if( dispDetsHid )
          cb->setChecked( true );
        
        if( dispDetsCb->isChecked() )
          cb->setChecked( false );
        
        break;
      }//case ApplyToAllDetectors:
        
      case ApplyToCbIndex::ApplyToDisplayedSamples:
      {
        // We want to check that for each displayed unique SpecMeas we are displaying all samples.
        /// \TODO: avoid allocating all the set<int>'s below, and then the looping threw every value
        ///        Can probably easily use a hyristic to avoid most of the time, or change how
        ///        things are tracked to avoid probably all the time (at the cost of adding
        ///        complexity to the code)
        map<shared_ptr<const SpecMeas>,set<int>> undisplayed;
        for( int i = 0; i < 3; ++i )
        {
          const auto meas = specfiles[i];
          if( !meas )
            continue;
          
          if( !undisplayed.count(meas) )
            undisplayed[meas] = meas->sample_numbers();
          
          set<int> &undispsamples = undisplayed[meas];
          for( const int sample : m_interspec->displayedSamples(spectypes[i]) )
            undispsamples.erase( sample );
        }//for( loop over the types of spectrum files )
        
        bool displayingAll = true;
        for( const auto &p : undisplayed )
          displayingAll = (displayingAll && p.second.empty());
        
        hideRow = displayingAll;
        if( displayingAll )
          cb->setChecked( false );
        
        break;
      }//case ApplyToDisplayedSamples:
        
      case ApplyToCbIndex::ApplyToAllSamples:
      {
        static_assert( ApplyToCbIndex::ApplyToDisplayedSamples
                         < ApplyToCbIndex::ApplyToAllSamples, "" );
        
        const WCheckBox *dispSamplesCb = m_applyToCbs[ApplyToCbIndex::ApplyToDisplayedSamples];
        assert( dispSamplesCb->parent() );
        const bool dispSamplesHid = dispSamplesCb->parent()->isHidden();
        
        hideRow = dispSamplesHid;
        if( dispSamplesHid )
          cb->setChecked( true );
        
        if( dispSamplesCb->isChecked() )
          cb->setChecked( false );
        
        break;
      }//case ApplyToCbIndex::ApplyToAllSamples:
        
        
      case ApplyToCbIndex::NumApplyToCbIndex:
        assert( 0 );
        break;
    }//switch( index )
    
    cbparent->setHidden( hideRow );
    
    if( !hideRow )
      anyApplyToCbShown = true;
  }//for( loop over ApplyToCbIndex )
  
  m_applyToColumn->setHidden( !anyApplyToCbShown );
  
  // Show the detector/spectrum-type selector whenever there is more than one CalDisplay to pick
  //  from - multiple detectors, multiple files, or the same file displayed as multiple spectrum
  //  types (whose displayed samples, and hence calibrations, can differ)
  const size_t total_cal_displays = wanted_entries[0].size() + wanted_entries[1].size()
                                    + wanted_entries[2].size();
  const bool hideDetCol = (total_cal_displays < 2);

  m_detColumn->setHidden( hideDetCol );
  
  for( MoreActionsIndex index = MoreActionsIndex(0);
      index < MoreActionsIndex::NumMoreActionsIndex;
      index = MoreActionsIndex(static_cast<int>(index) + 1) )
  {
    Wt::WAnchor *anchor = m_moreActions[static_cast<int>(index)];
    assert( anchor );
    auto aparent = anchor->parent();
    assert( dynamic_cast<WContainerWidget *>(aparent) );
    
    switch( index )
    {
      case MoreActionsIndex::Linearize:
      case MoreActionsIndex::Truncate:
      case MoreActionsIndex::CombineChannels:
        aparent->setHidden( !specfiles[0] );
        break;
        
      case MoreActionsIndex::ConvertToFrf:
        aparent->setHidden( !hasPolyCal );
        break;
        
      case MoreActionsIndex::ConvertToPoly:
        aparent->setHidden( !(hasFRFCal || hasLowerChanCal) );
        break;
        
      case MoreActionsIndex::MultipleFilesCal:
      {
        SpecMeasManager *manager = m_interspec->fileManager();
        SpectraFileModel *fmodel = manager ? manager->model() : nullptr;
        const int nfiles = fmodel ? fmodel->rowCount() : 0;
        int nRecordsWithPeaks = 0;
        for( int row = 0; row < nfiles && (nRecordsWithPeaks < 2); ++row )
        {
          shared_ptr<SpectraFileHeader> header = fmodel->fileHeader( row );
          if( !header )
            continue;
          
          int nsamples = header->numSamples();
          shared_ptr<SpecMeas> meas = header->measurementIfInMemory();
          if( meas )
            nRecordsWithPeaks += meas->sampleNumsWithPeaks().size();
          else
            nRecordsWithPeaks += nsamples; //not worth readin file from disk, so we'll be hopeful
        }
        
        aparent->setHidden( nRecordsWithPeaks < 2 );
        break;
      }//case MoreActionsIndex::MultipleFilesCal:
        
      case MoreActionsIndex::NumMoreActionsIndex:
        break;
    }//switch( index )
  }//for( loop over
  
  // Update the "Fit Coeffs" button to be enabled/disabled
  updateFitButtonStatus();
  
}//void doRefreshFromFiles()



void EnergyCalTool::moreActionBtnClicked( const MoreActionsIndex index )
{
  const vector<MeasToApplyCoefChangeTo> measToChange = measurementsToApplyCoeffChangeTo();

  if( m_addActionWindow )
  {
    AuxWindow::deleteAuxWindow( m_addActionWindow );
    m_addActionWindow = nullptr;
  }
  
  m_addActionWindow = new EnergyCalAddActionsWindow( index, measToChange, this );
  m_addActionWindow->finished().connect( this, &EnergyCalTool::cancelMoreActionWindow );
  
  UndoRedoManager *undoManager = m_interspec->undoRedoManager();
  if( !undoManager )
    return;
  
  auto undo = [](){
    InterSpec *viewer = InterSpec::instance();
    EnergyCalTool *tool = viewer ? viewer->energyCalTool() : nullptr;
    if( tool )
      tool->cancelMoreActionWindow();
  };
  
  auto redo = [index](){
    InterSpec *viewer = InterSpec::instance();
    EnergyCalTool *tool = viewer ? viewer->energyCalTool() : nullptr;
    if( tool )
      tool->moreActionBtnClicked(index);
  };
  
  undoManager->addUndoRedoStep( undo, redo, "Show additional energy cal tool.");
}//void moreActionBtnClicked( const MoreActionsIndex index )


void EnergyCalTool::cancelMoreActionWindow()
{
  if( m_addActionWindow )
  {
    AuxWindow::deleteAuxWindow( m_addActionWindow );
    m_addActionWindow = nullptr;
  }
}//void cancelMoreActionWindow()


void EnergyCalTool::render( Wt::WFlags<Wt::RenderFlag> flags)
{
  //flags.testFlag(RenderFlag::RenderFull) will only be true on initial rending of widget, and
  //  after that only the RenderFlag::RenderUpdate flag will be set
  
  if( flags.testFlag(Wt::RenderFlag::RenderFull)
      || m_renderFlags.testFlag(EnergyCalToolRenderFlags::FullGuiUpdate) )
  {
    doRefreshFromFiles();
    m_renderFlags.clear( EnergyCalToolRenderFlags::FullGuiUpdate );
  }
  
  WContainerWidget::render(flags);
}//void render( Wt::WFlags<Wt::RenderFlag> flags)


void EnergyCalTool::applyToCbChanged( const EnergyCalTool::ApplyToCbIndex index )
{
  // We only get here if the user checked/unchecked a checkbox, or in a undo/redo step.
  assert( index >= 0 && index <= EnergyCalTool::NumApplyToCbIndex );
  WCheckBox *cb = m_applyToCbs[index];
  const bool isChecked = cb->isChecked();
  
  // Grab the starting state of all checkboxed
  bool startingState[ApplyToCbIndex::NumApplyToCbIndex];
  for( ApplyToCbIndex i = ApplyToCbIndex(0);
      i < ApplyToCbIndex::NumApplyToCbIndex;
      i = ApplyToCbIndex(i + 1) )
  {
    startingState[i] = m_applyToCbs[i]->isChecked();
  }
  
  // Assume `index` was actually not what it is now.
  startingState[index] = !startingState[index];
  
  
  switch( index )
  {
    case EnergyCalTool::ApplyToForeground:
      updateFitButtonStatus();
      break;
      
    case EnergyCalTool::ApplyToBackground:
    case EnergyCalTool::ApplyToSecondary:
      break;
    
    case EnergyCalTool::ApplyToDisplayedDetectors:
      m_applyToCbs[EnergyCalTool::ApplyToAllDetectors]->setChecked( !isChecked );
      break;
      
    case EnergyCalTool::ApplyToAllDetectors:
      m_applyToCbs[EnergyCalTool::ApplyToDisplayedDetectors]->setChecked( !isChecked );
      break;
      
    case EnergyCalTool::ApplyToDisplayedSamples:
      m_applyToCbs[EnergyCalTool::ApplyToAllSamples]->setChecked( !isChecked );
      break;
      
    case EnergyCalTool::ApplyToAllSamples:
      m_applyToCbs[EnergyCalTool::ApplyToDisplayedSamples]->setChecked( !isChecked );
      break;
      
    case EnergyCalTool::NumApplyToCbIndex:
      break;
  }//switch( index )
  
  UndoRedoManager *undoManager = m_interspec->undoRedoManager();
  if( !undoManager || undoManager->isInUndoOrRedo() )
    return;
  
  
  bool finalState[ApplyToCbIndex::NumApplyToCbIndex];
  for( ApplyToCbIndex i = ApplyToCbIndex(0);
      i < ApplyToCbIndex::NumApplyToCbIndex;
      i = ApplyToCbIndex(i + 1) )
  {
    finalState[i] = m_applyToCbs[i]->isChecked();
  }
  
  
  auto undo = [index,isChecked,startingState](){
    InterSpec *viewer = InterSpec::instance();
    EnergyCalTool *tool = viewer ? viewer->energyCalTool() : nullptr;
    if( !tool )
      return;
    
    for( ApplyToCbIndex i = ApplyToCbIndex(0);
        i < ApplyToCbIndex::NumApplyToCbIndex;
        i = ApplyToCbIndex(i + 1) )
    {
      tool->m_applyToCbs[i]->setChecked( startingState[i] );
    }
    
    tool->applyToCbChanged( index );
  };
  
  auto redo = [index,isChecked,finalState](){
    InterSpec *viewer = InterSpec::instance();
    EnergyCalTool *tool = viewer ? viewer->energyCalTool() : nullptr;
    if( !tool )
      return;
    
    for( ApplyToCbIndex i = ApplyToCbIndex(0);
        i < ApplyToCbIndex::NumApplyToCbIndex;
        i = ApplyToCbIndex(i + 1) )
    {
      tool->m_applyToCbs[i]->setChecked( finalState[i] );
    }
    
    tool->applyToCbChanged( index );
  };
  
  const string label = cb->text().toUTF8();
  undoManager->addUndoRedoStep( undo, redo, "Change " + label );
}//void applyToCbChanged( const ApplyToCbIndex index )
