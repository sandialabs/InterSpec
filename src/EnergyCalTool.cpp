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
#include <array>
#include <string>
#include <vector>
#include <ctime>
#include <memory>
#include <iostream>
#include <algorithm>
#include <functional>

#include <Wt/Utils.h>
#include <Wt/WMenu.h>
#include <Wt/WText.h>
#include <Wt/WLabel.h>
#include <Wt/WAnchor.h>
#include <Wt/WResource.h>
#include <Wt/WCheckBox.h>
#include <Wt/WFileUpload.h>
#include <Wt/WGridLayout.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>
#include <Wt/Http/Request.h>
#include <Wt/Http/Response.h>
#include <Wt/WItemDelegate.h>
#include <Wt/WStackedWidget.h>
#include <Wt/WContainerWidget.h>

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
#include "InterSpec/WidgetUtils.h"
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
  CALpDownloadResource( EnergyCalTool *tool, InterSpec *viewer )
  : WResource(), m_app( WApplication::instance() ), m_interspec( viewer ), m_tool( tool )
  {
    assert( m_app );
    assert( m_tool );
    assert( m_interspec );

    // Have Wt hold the session lock across handleRequest instead of dropping it for us: by
    //  default `WResource::handle` unlocks before invoking us, so the session thread can be
    //  inside `~WResource`'s `beingDeleted()` - waiting on our use-count - while we block
    //  acquiring the lock it holds.  That is a deadlock.  Free: the `WApplication::UpdateLock`
    //  in handleRequest already spanned the whole body, so this only moves it earlier.
    setTakesUpdateLock( true );
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
      // setTakesUpdateLock(true) in the ctor means Wt already holds the session lock here, so this
      //  is a same-thread no-op (WApplication::UpdateLock detects that and returns immediately).
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
      
      suggestFileName( filename, ContentDisposition::Attachment );
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

  CoefDisplay( const size_t order )
  : WContainerWidget(),
    m_order( order ),
    m_label( nullptr ),
    m_fit( nullptr ),
    m_value( nullptr )
  {
    addStyleClass( "CoefDisplay" );

    m_label = this->addNew<WLabel>( coef_disp_label_txt(order, false) );
    m_label->addStyleClass( "CoefLabel" );

    m_value = this->addNew<NativeFloatSpinBox>();
    m_value->setSpinnerHidden( true );

    m_fit = this->addNew<WCheckBox>( "Fit" );
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
  WText *m_fitHint;
#endif

  /** Whether the per-coefficient "Fit" checkboxes are disabled; they are disabled on displays
   whose checkboxes are never consulted (fitting always uses the foreground displays selections;
   the checkbox states are kept in sync across all displays).
   */
  bool m_fitCbsDisabled;
  
#if( IMP_CALp_BTN_NEAR_COEFS )
  WPushButton *m_downloadCALp;
  WPushButton *m_uploadCALp;
#endif
  
  shared_ptr<const SpecUtils::EnergyCalibration> m_cal;
  
public:
  CalDisplay( EnergyCalTool *tool,
             const SpecUtils::SpectrumType type,
             const std::string &detname,
             const bool isWideLayout )
  : WContainerWidget(),
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
  , m_fitCoeffs( nullptr ),
   m_fitHint( nullptr )
#endif
  , m_fitCbsDisabled( false )
#if( IMP_CALp_BTN_NEAR_COEFS )
  , m_downloadCALp( nullptr )
  , m_uploadCALp( nullptr )
#endif
  {
    addStyleClass( "CalDisplay" );
    
    WGridLayout *layout = setLayout( std::make_unique<WGridLayout>() );
    layout->setContentsMargins( 0, 0, 0, 0 );
    layout->setVerticalSpacing( 0 );
    layout->setHorizontalSpacing( 0 );
    
    WContainerWidget *coefDiv = layout->addWidget( std::make_unique<WContainerWidget>(), 0, 0 );
    coefDiv->addStyleClass( "CoefCol" );
    
    m_type = coefDiv->addNew<WText>( "&nbsp;" );
    m_type->setInline( false );
    m_type->addStyleClass( "Wt-itemview Wt-header Wt-label CalType" );
    
    m_coefficients = coefDiv->addNew<WContainerWidget>();
    m_coefficients->addStyleClass( "CoefContent" );
    
#if( IMP_COEF_FIT_BTN_NEAR_COEFS || IMP_CALp_BTN_NEAR_COEFS )
    WContainerWidget *btndiv = layout->addWidget( std::make_unique<WContainerWidget>(), 1, 0 );
    btndiv->addStyleClass( "CalCoefsBtnDiv" );
  
#if( IMP_CALp_BTN_NEAR_COEFS )

#if( BUILD_AS_OSX_APP || IOS )
    {
      WLink calpLink( m_tool->calpResources() );
      calpLink.setTarget( Wt::LinkTarget::NewWindow );
      m_downloadCALp = btndiv->addNew<WAnchor>( calpLink );
    }
    m_downloadCALp->setStyleClass( "LinkBtn DownloadLink" );
#else
    m_downloadCALp = btndiv->addNew<WPushButton>();
    m_downloadCALp->setIcon( "InterSpec_resources/images/download_small.svg" );
    {
      WLink calpLink( m_tool->calpResources() );
      calpLink.setTarget( Wt::LinkTarget::NewWindow );
      m_downloadCALp->setLink( calpLink );
    }
    m_downloadCALp->setStyleClass( "LinkBtn DownloadBtn CALp" );
    
#endif //#if( BUILD_AS_OSX_APP || IOS ) / #else

    m_downloadCALp->setText( "CALp" );
    
    m_uploadCALp = btndiv->addNew<WPushButton>();
    m_uploadCALp->setIcon( "InterSpec_resources/images/upload_small.svg" );
    m_uploadCALp->setStyleClass( "LinkBtn UploadBtn CALp" );
    m_uploadCALp->clicked().connect( m_tool, &EnergyCalTool::handleRequestToUploadCALp );
#endif //#if( IMP_CALp_BTN_NEAR_COEFS )
    
    WContainerWidget *spacer = btndiv->addNew<WContainerWidget>();
    spacer->addStyleClass( "Spacer" );

#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
    // The hint lives inside the flex spacer, absolutely positioned to fill it, so its text wraps
    //  within the leftover space and can never change the column width or move the button.
    m_fitHint = spacer->addNew<WText>();
    m_fitHint->addStyleClass( "FitCoefHint" );
    m_fitHint->hide();

    m_fitCoeffs = btndiv->addNew<WPushButton>( WString::tr("ect-fit-coeff-btn") );
    m_fitCoeffs->addStyleClass( "CalCoefFitBtn" );

    const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips",
                                                                      InterSpec::instance() );
    HelpSystem::attachToolTipOn( m_fitCoeffs, WString::tr("ect-tt-fit-coeff-btn"), showToolTips );
#endif
    
    layout->setRowStretch( 0, 1 );
#endif //#if( IMP_COEF_FIT_BTN_NEAR_COEFS || IMP_CALp_BTN_NEAR_COEFS )
    
    auto devPairsOwned = std::make_unique<DeviationPairDisplay>();
    m_devPairs = devPairsOwned.get();
    
#if( HIDE_EMPTY_DEV_PAIRS )
    //For files with multiple detectors, the "Add dev. pairs" buttons doesnt show up right
    // for the detectors not currently showing - I guess should toggle dev pairs for all detectors.
    m_devPairs->setHidden( true );
    auto addPairsOwned = std::make_unique<WPushButton>("Add dev. pairs");
    m_addPairs = addPairsOwned.get();
    m_addPairs->addStyleClass( "LinkBtn" );
    //m_addPairs->setIcon( "InterSpec_resources/images/plus_min_white.svg" );
    m_addPairs->setHidden( true );
    m_addPairs->clicked().connect( this, &CalDisplay::showDevPairs );
#endif
    
    if( isWideLayout )
    {
#if( IMP_COEF_FIT_BTN_NEAR_COEFS || IMP_CALp_BTN_NEAR_COEFS )
      layout->addWidget( std::move(devPairsOwned), 0, 1, 2, 1 );
#else
      layout->addWidget( std::move(devPairsOwned), 0, 1 );
#endif
      
#if( HIDE_EMPTY_DEV_PAIRS )
      layout->addWidget( std::move(addPairsOwned), 1, 0, AlignmentFlag::Right );
      layout->setRowStretch( 0, 1 );
#endif
    }else
    {
#if( HIDE_EMPTY_DEV_PAIRS )
      layout->addWidget( std::move(addPairsOwned), layout->rowCount(), 0, AlignmentFlag::Center );
      layout->addWidget( std::move(devPairsOwned), layout->rowCount(), 0 );
#else
      layout->addWidget( std::move(devPairsOwned), layout->rowCount(), 0 );
#endif
      m_devPairs->setHeight( 100 );
    }
    
    m_devPairs->changed().connect( this, [this]( const int a1 ){
      m_tool->userChangedDeviationPair( this, a1 );
    } );
  }//CalDisplay( constructor )
  
  SpecUtils::SpectrumType spectrumType() const { return m_cal_type; }
  const std::string &detectorName() const { return m_det_name; }
  
    
  //void setDeviationPairMsg( const std::string &msg )
  //{
  //  m_devPairs->setMsg( msg );
  //}//void setDeviationPairMsg( const std::string &msg )
  
  
  shared_ptr<const SpecUtils::EnergyCalibration> lastSetCalibration() const
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


  /** Enable/disable the per-coefficient "Fit" checkboxes.  They are disabled on displays whose
   checkboxes are never consulted (fitting always uses the foreground displays selections; the
   checkbox states are kept in sync across all displays).
   */
  void setFitCheckboxesDisabled( const bool disabled )
  {
    m_fitCbsDisabled = disabled;

    const WString tip = disabled ? WString::tr("ect-tt-fit-cb-disabled") : WString();
    for( auto w : m_coefficients->children() )
    {
      auto ww = dynamic_cast<CoefDisplay *>( w );
      assert( ww );
      if( ww && (ww->m_fit->isDisabled() != disabled) )
      {
        ww->m_fit->setDisabled( disabled );
        ww->m_fit->setToolTip( tip );
      }
    }//for( auto w : m_coefficients->children() )
  }//void setFitCheckboxesDisabled( const bool disabled )


#if( HIDE_EMPTY_DEV_PAIRS )
  void showDevPairs()
  {
    m_addPairs->setHidden( true );
    m_devPairs->setHidden( false, WAnimation(AnimationEffect::Fade, TimingFunction::Linear, 200) );
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
    if( m_fitHint && !fitCoeffsVisible && !m_fitHint->isHidden() )
      m_fitHint->hide();
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
            m_convertMsg = p->addNew<WText>( WString::tr("ect-convert-to-poly-msg") );
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
          WidgetUtils::removeWidgetNow( m_convertMsg );
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
          WidgetUtils::removeWidgetNow( m_convertMsg );
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
    const auto anim = hadCal ? WAnimation(AnimationEffect::Fade, TimingFunction::Linear, 200) : WAnimation{};
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
          // Wt4: removeWidget() hands back ownership; the unique_ptr destroys it here.
          const std::unique_ptr<WWidget> doomed = m_coefficients->removeWidget( existing[i] );
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
      CoefDisplay *disp = m_coefficients->addNew<CoefDisplay>( coefnum );
      coef_disps[coefnum] = disp;
      disp->setLabeling( lower_channel );
      const float value = (coefnum < coeffs.size()) ? coeffs[coefnum] : 0.0f;
      disp->m_value->setValue( value );
      disp->m_fit->setChecked( (coefnum < 2) );
      if( m_fitCbsDisabled )
      {
        disp->m_fit->disable();
        disp->m_fit->setToolTip( WString::tr("ect-tt-fit-cb-disabled") );
      }

      disp->m_fit->changed().connect( this, [this](){ m_tool->fitCheckboxChanged( this ); } );
      
      /* Note: if the user uses the up.down arrows in a NativeFloatSpinBox to change values, things
               get all messed up (new values get set via c++ messing  up current values, or the
               valueChanged() callback gets called like 10 times per second, causing changes faster
               than everything can keep up, and just generally poor working), so for now I have
               disabled these spinners via #NativeFloatSpinBox::setSpinnerHidden()
       */
      disp->m_value->valueChanged().connect( this, [this,coefnum](){
        m_tool->userChangedCoefficient( coefnum, this );
      } );
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
  /** Sets the "Fit Coeffs" button enabled state, and when disabled, shows a short hint of why
   fitting isnt possible (pass an empty hint to show nothing).
   */
  void setFitStatus( const bool canFitCoefs, const Wt::WString &hint )
  {
    if( !m_cal )
      return;

    if( m_fitCoeffs->isEnabled() != canFitCoefs )
      m_fitCoeffs->setEnabled( canFitCoefs );

    const bool showHint = (!canFitCoefs && !hint.empty() && !m_fitCoeffs->isHidden());
    if( showHint && (m_fitHint->text() != hint) )
      m_fitHint->setText( hint );
    if( m_fitHint->isHidden() == showHint )
      m_fitHint->setHidden( !showHint );
  }//void setFitStatus( const bool canFitCoefs, const Wt::WString &hint )
  
  Wt::EventSignal<Wt::WMouseEvent> &doFitCoeffs()
  {
    return m_fitCoeffs->clicked();
  }
#endif
  
};//class CalDisplay

const size_t CalDisplay::sm_min_coef_display = 4;


}//namespace


EnergyCalTool::EnergyCalTool( InterSpec *viewer, std::shared_ptr<PeakModel> peakModel )
: WContainerWidget(),
  m_interspec( viewer ),
  m_peakModel( std::move(peakModel) ),
  m_calpResource(),
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
  m_haveRelayoutSelection( false ),
  m_relayoutSelectionType( SpecUtils::SpectrumType::Foreground ),
  m_relayoutSelectionDet(),
  m_addActionWindow( nullptr )
{
  wApp->useStyleSheet( "InterSpec_resources/EnergyCalTool.css" );
  
  assert( viewer );
  viewer->useMessageResourceBundle( "EnergyCalTool" );
    
  addStyleClass( "EnergyCalTool" );

  m_calpResource = std::make_shared<EnergyCalImp::CALpDownloadResource>( this, m_interspec );
  
  initWidgets( EnergyCalTool::LayoutType::Wide );
}


void EnergyCalTool::initWidgets( EnergyCalTool::LayoutType layoutType )
{
  const bool wide = (layoutType == LayoutType::Wide);
  
  if( (wide && !m_tallLayoutContent && m_layout) || (!wide && m_tallLayoutContent) )
    return;
  
  if( m_graphicalRecal )
  {
    AuxWindow::deleteAuxWindow( m_graphicalRecal.get() );
    assert( !m_graphicalRecal );   //observing_ptr auto-clears
  }//if( m_graphicalRecal )

  // Keep the users "Apply Changes To" selections, and the selected calibration display (which
  //  "Indep. Dets." acts on), when switching between wide/tall layouts (i.e., when the tool
  //  is moved between the tool tabs and its own window).
  const ApplyToSelection prevApplyTo = currentApplyToSelection();
  {
    const EnergyCalImp::CalDisplay *current = m_calInfoDisplayStack
          ? dynamic_cast<const EnergyCalImp::CalDisplay *>( m_calInfoDisplayStack->currentWidget() )
          : nullptr;
    m_haveRelayoutSelection = !!current;
    if( current )
    {
      m_relayoutSelectionType = current->spectrumType();
      m_relayoutSelectionDet = current->detectorName();
    }
  }

  if( wide )
  {
    removeStyleClass( "TallEnergyCal" );
    if( m_tallLayoutContent )
    {
      // Wt4: removeWidget() hands back ownership; the unique_ptr destroys it here.
      auto p = removeWidget( m_tallLayoutContent );
      m_tallLayoutContent = nullptr;
    }

    m_layout = setLayout( std::make_unique<WGridLayout>() );
  }else
  {
    addStyleClass( "TallEnergyCal" );

    // We must drop any layout we already have on `this` BEFORE adding m_tallLayoutContent:
    //  WContainerWidget renders only its layout's widgets when a layout is set, so a plain child
    //  added while the wide layout is still installed would never get a DOM element (and the stale
    //  wide widgets would keep rendering, wired to the same handlers).
    //  setLayout(nullptr) -> setLogicalLayout() -> clear(), which destroys the old layout and all
    //  the widgets it managed exactly once.
    setLayout( std::unique_ptr<WLayout>() );
    m_layout = nullptr;

    m_tallLayoutContent = addNew<WContainerWidget>();
    m_layout = m_tallLayoutContent->setLayout( std::make_unique<WGridLayout>() );
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
  
  auto noCalTxtOwned = std::make_unique<WText>(WString::tr("ect-no-spec"));
  m_noCalTxt = noCalTxtOwned.get();
  m_noCalTxt->addStyleClass( "NoCalContentTxt" );
  if( wide )
    m_layout->addWidget( std::move(noCalTxtOwned), 0, 0, AlignmentFlag::Center | AlignmentFlag::Middle );
  else
    m_layout->addWidget( std::move(noCalTxtOwned), 0, 0, 1, 2, AlignmentFlag::Center | AlignmentFlag::Middle );
  
  //Create the more actions column...
  auto oreActionsColumnOwned = std::make_unique<WContainerWidget>();
  m_moreActionsColumn = oreActionsColumnOwned.get();
  m_moreActionsColumn->addStyleClass( "ToolTabTitledColumn MoreActionCol" );
  if( wide )
    m_layout->addWidget( std::move(oreActionsColumnOwned), 0, 1 );
  else
    m_layout->addWidget( std::move(oreActionsColumnOwned), 3, 0 );
  
  WGridLayout *collayout = m_moreActionsColumn->setLayout( std::make_unique<WGridLayout>() );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );
  
  WText *header = collayout->addWidget( std::make_unique<WText>(WString::tr("ect-more-act")), 0, 0 );
  header->addStyleClass( "ToolTabColumnTitle" );
  
  //We will put the apply-to list inside a div so we can style consistently with other rows
  // (a <ul> element doesnt accept same css as <div>, apparently).
  WContainerWidget *moreActionsDiv = collayout->addWidget( std::make_unique<WContainerWidget>(), 1, 0 );
  moreActionsDiv->addStyleClass( "ToolTabTitledColumnContent MoreActionsMenuContent" );
  collayout->setRowStretch( 1, 1 );
  
  WContainerWidget *moreActionsList = moreActionsDiv->addNew<WContainerWidget>();
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

      case MoreActionsIndex::GainMatch:
        label = "ect-gain-match";
        tooltip = "ect-tt-gain-match";
        break;

      case MoreActionsIndex::GuessGain:
        label = "ect-guess-gain";
        tooltip = "ect-tt-guess-gain";
        break;

      case MoreActionsIndex::NumMoreActionsIndex:
        assert(0);
        break;
    }//switch( index )
    
    WContainerWidget *holder = moreActionsList->addNew<WContainerWidget>();
    m_moreActions[static_cast<int>(index)] = holder->addNew<WAnchor>( WLink(), WString::tr(label) );
    m_moreActions[static_cast<int>(index)]->clicked().connect( this, [this,index](){
      moreActionBtnClicked( index );
    } );
    
    assert( tooltip );
    if( tooltip )
      HelpSystem::attachToolTipOn( holder, WString::tr(tooltip), showToolTips );
  }//for( loop over more actions )
  
  WContainerWidget *btndiv = collayout->addWidget( std::make_unique<WContainerWidget>(), 2, 0 );
  btndiv->addStyleClass( "BtmBtnDiv" );
  
  auto helpBtn = btndiv->addNew<WContainerWidget>();
  helpBtn->addStyleClass( "Wt-icon ContentHelpBtn" );
  helpBtn->clicked().connect( this, [](){ HelpSystem::createHelpWindow( "energy-calibration" ); } );

  
#if( !IMP_CALp_BTN_NEAR_COEFS )
  m_uploadCALp = btndiv->addNew<WPushButton>();
  m_uploadCALp->setIcon( "InterSpec_resources/images/upload_small.svg" );
  m_uploadCALp->setStyleClass( "LinkBtn UploadBtn CALp" );
  m_uploadCALp->clicked().connect( this, &EnergyCalTool::handleRequestToUploadCALp );
  
#if( BUILD_AS_OSX_APP || IOS )
  {
    WLink calpLink( m_calpResource );
    calpLink.setTarget( Wt::LinkTarget::NewWindow );
    m_downloadCALp = btndiv->addNew<WAnchor>( calpLink );
  }
  m_downloadCALp->setStyleClass( "LinkBtn DownloadLink CALp" );
#else
  m_downloadCALp = btndiv->addNew<WPushButton>();
  m_downloadCALp->setIcon( "InterSpec_resources/images/download_small.svg" );
  {
    WLink calLink( m_calpResource );
    calLink.setTarget( Wt::LinkTarget::NewWindow );
    m_downloadCALp->setLink( calLink );
  }
  m_downloadCALp->setStyleClass( "LinkBtn DownloadBtn CALp" );
  
#endif
  m_downloadCALp->setText( "CALp" );
  
  m_downloadCALp->clicked().connect( this, [this](){
    m_interspec->logMessage( WString::tr("ect-export-CALp-msg"), WarningWidget::WarningMsgInfo );
  } );
  HelpSystem::attachToolTipOn( m_downloadCALp, WString::tr("ect-tt-CALp"), showToolTips );
  
  m_downloadCALp->setHidden( true );
  m_uploadCALp->setHidden( true );
#endif // !IMP_CALp_BTN_NEAR_COEFS
  
  
#if( !IMP_COEF_FIT_BTN_NEAR_COEFS )
  m_fitCalBtn = btndiv->addNew<WPushButton>( WString::tr("ect-fit-coeff-btn") );
  m_fitCalBtn->addStyleClass( "FitCoefBtn" );
  m_fitCalBtn->clicked().connect( this, &EnergyCalTool::fitCoefficients );
  m_fitCalBtn->setDisabled( true );
  HelpSystem::attachToolTipOn( m_fitCalBtn, WString::tr("ect-tt-fit-coeff-btn"), showToolTips );
#endif // !IMP_COEF_FIT_BTN_NEAR_COEFS
  
  // Create the "Apply To" column that determines what to apply changes to
  auto applyToColumnOwned = std::make_unique<WContainerWidget>();
  m_applyToColumn = applyToColumnOwned.get();
  m_applyToColumn->addStyleClass( "ToolTabTitledColumn ApplyToCol" );
  if( wide )
    m_layout->addWidget( std::move(applyToColumnOwned), 0, 2 );
  else
    m_layout->addWidget( std::move(applyToColumnOwned), 3, 1 );
  
  collayout = m_applyToColumn->setLayout( std::make_unique<WGridLayout>() );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );
  
  auto headerOwned = std::make_unique<WText>(WString::tr("ect-apply-changes-to"));
  header = headerOwned.get();
  header->addStyleClass( "ToolTabColumnTitle" );
  collayout->addWidget( std::move(headerOwned), 0, 0 );
  
  //We will put the apply-to list inside a div so we can style consistently with other rows
  // (a <ul> element doesnt accept same css as <div>, apparently).
  WContainerWidget *applyToDiv = collayout->addWidget( std::make_unique<WContainerWidget>(), 1, 0 );
  applyToDiv->addStyleClass( "ToolTabTitledColumnContent ApplyToMenuContent" );
  
  WContainerWidget *applyToList = applyToDiv->addNew<WContainerWidget>();
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
      case ApplyToCbIndex::ApplyToSelectedDetOnly:    label = "ect-selected-det-only"; break;
      case ApplyToCbIndex::ApplyToDisplayedSamples:   label = "ect-disp-samples"; break;
      case ApplyToCbIndex::ApplyToAllSamples:         label = "ect-all-samples";  break;
      case ApplyToCbIndex::NumApplyToCbIndex:
        assert( 0 );
        break;
    }//switch( index )
    
    WContainerWidget *item = applyToList->addNew<WContainerWidget>();
    item->addStyleClass( "ApplyToItem" );
    auto cb = item->addNew<WCheckBox>( WString::tr(label) );
    cb->setWordWrap( false );
    cb->addStyleClass( "ApplyToItem CbNoLineBreak" );
    cb->setChecked( prevApplyTo.checked[index] );  //defaults, or the previous layouts state
    
    if( index == ApplyToCbIndex::ApplyToSelectedDetOnly )
      HelpSystem::attachToolTipOn( item, WString::tr("ect-tt-selected-det-only"), showToolTips );
    
    cb->checked().connect( this, [this,index](){ applyToCbChanged( index ); } );
    cb->unChecked().connect( this, [this,index](){ applyToCbChanged( index ); } );
    
    m_applyToCbs[index] = cb;
  }//for( loop over ApplyToCbIndex )
  
  WAnimation animation(Wt::AnimationEffect::Fade, Wt::TimingFunction::Linear, 200);

  
  // Create the "Coefficients" column that show the polynomial/FRF coefficents.
  auto calColumnOwned = std::make_unique<WContainerWidget>();
  m_calColumn = calColumnOwned.get();
  m_calColumn->addStyleClass( "ToolTabTitledColumn CoefColumn" );
  if( wide )
    m_layout->addWidget( std::move(calColumnOwned), 0, 3 );
  else
    m_layout->addWidget( std::move(calColumnOwned), 1, 0, 1, 2 );
  
  
  collayout = m_calColumn->setLayout( std::make_unique<WGridLayout>() );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );
  
  if( wide )
    collayout->setColumnStretch( 1, 1 );
  
  auto coefHeaderOwned = std::make_unique<WText>( WString::tr("ect-calib-coeffs") );
  header = coefHeaderOwned.get();
  header->addStyleClass( "ToolTabColumnTitle" );

  collayout->addWidget( std::move(coefHeaderOwned), 0, 0, 1, 2 );
  //collayout->addWidget( m_calInfoDisplayStack, 1, 0 );
  
  
  // Create the "Detector" column that determines which coefficients to show
  auto detColumnOwned = std::make_unique<WContainerWidget>();
  m_detColumn = detColumnOwned.get();
  m_detColumn->addStyleClass( "DetCol" );
  collayout->addWidget( std::move(detColumnOwned), 1, 0 );
  
  m_detColLayout = m_detColumn->setLayout( std::make_unique<WGridLayout>() );
  m_detColLayout->setContentsMargins( 0, 0, 0, 0 );
  m_detColLayout->setVerticalSpacing( 0 );
  m_detColLayout->setHorizontalSpacing( 0 );
  
  auto detheaderOwned = std::make_unique<WText>( WString::tr("Detector") );
  WText *detheader = detheaderOwned.get();
  detheader->setInline( false );
  detheader->addStyleClass( "DetHdr Wt-itemview Wt-header Wt-label" );
  //detheader->resize( WLength::Auto, WLength(20,WLength::Unit::Pixel) );
  //collayout->addWidget( std::move(detheaderOwned), 0, 0 );
  m_detColLayout->addWidget( std::move(detheaderOwned), 0, 0  );
  m_detColLayout->setRowStretch( 2, 1 );
  
  // Create the "Cal Peaks" table
  auto peakTableColumnOwned = std::make_unique<WContainerWidget>();
  m_peakTableColumn = peakTableColumnOwned.get();
  m_peakTableColumn->addStyleClass( "ToolTabTitledColumn PeakTableCol" );
  if( wide )
  {
    m_layout->addWidget( std::move(peakTableColumnOwned), 0, 4 );
    m_layout->setColumnStretch( 4, 1 );
  }else
  {
    m_layout->addWidget( std::move(peakTableColumnOwned), 2, 0, 1, 2 );
  }

  collayout = m_peakTableColumn->setLayout( std::make_unique<WGridLayout>() );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  if( wide )
    collayout->setRowStretch( 1, 1 );
  
  auto peakHeaderOwned = std::make_unique<WText>( WString::tr("ect-cal-peaks") );
  header = peakHeaderOwned.get();
  header->addStyleClass( "ToolTabColumnTitle" );
  collayout->addWidget( std::move(peakHeaderOwned), 0, 0 );
  
  auto peakTableOwned = std::make_unique<RowStretchTreeView>();
  m_peakTable = peakTableOwned.get();
  m_peakTable->addStyleClass( "ToolTabTitledColumnContent PeakTable" );
  collayout->addWidget( std::move(peakTableOwned), 1, 0 );
  collayout->setRowStretch( 1, 1 );
  
  m_peakTable->setRootIsDecorated( false ); //makes the tree look like a table! :)
  m_peakTable->setModel( m_peakModel );
  const int numModelCol = m_peakModel->columnCount();
  for( int col = 0; col < numModelCol; ++col )
    m_peakTable->setColumnHidden( col, true );
  
  m_peakTable->setSortingEnabled( true );
  m_peakTable->setAlternatingRowColors( true );
  m_peakTable->setSelectable( true );
  m_peakTable->setSelectionMode( SelectionMode::Single );
  m_peakTable->setEditTriggers( Wt::EditTrigger::SingleClicked
                               | Wt::EditTrigger::DoubleClicked );
  
  m_peakTable->setColumnHidden( PeakModel::kUseForCalibration, false );
  m_peakTable->setColumnHidden( PeakModel::kMean, false );
  m_peakTable->setColumnHidden( PeakModel::kIsotope, false );
  m_peakTable->setColumnHidden( PeakModel::kPhotoPeakEnergy, false );
  m_peakTable->setColumnHidden( PeakModel::kDifference, false );
  
  
  m_peakTable->setColumnWidth( PeakModel::kUseForCalibration, WLength(3.7, WLength::Unit::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kMean, WLength(4.5, WLength::Unit::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kIsotope, WLength(4.5, WLength::Unit::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kPhotoPeakEnergy, WLength(6.25, WLength::Unit::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kDifference, WLength(5, WLength::Unit::FontEm) );
  
  
  
  auto dblDelagate = std::make_shared<WItemDelegate>();
  dblDelagate->setTextFormat( "%.2f" );
  m_peakTable->setItemDelegateForColumn( PeakModel::kMean, dblDelagate );
  
  auto nuclideDelegate = std::make_shared<PhotopeakDelegate>( PhotopeakDelegate::NuclideDelegate, true );
  m_peakTable->setItemDelegateForColumn( PeakModel::kIsotope, nuclideDelegate );
  
  auto photopeakDelegate = std::make_shared<PhotopeakDelegate>( PhotopeakDelegate::GammaEnergyDelegate, true );
  m_peakTable->setItemDelegateForColumn( PeakModel::kPhotoPeakEnergy, photopeakDelegate );
  
  m_peakModel->dataChanged().connect( this, &EnergyCalTool::updateFitButtonStatus );
  m_peakModel->rowsRemoved().connect( this, &EnergyCalTool::updateFitButtonStatus );
  m_peakModel->rowsInserted().connect( this, &EnergyCalTool::updateFitButtonStatus );
  m_peakModel->layoutChanged().connect( this, &EnergyCalTool::updateFitButtonStatus );
  
  m_interspec->displayedSpectrumChanged().connect( this,
    [this]( const SpecUtils::SpectrumType a1, const std::shared_ptr<SpecMeas> a2,
            const std::set<int> &a3, const std::vector<std::string> &a4 ){
      displayedSpectrumChanged( a1, a2, a3, a4 );
    } );
  
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
  // Both of these are AuxWindows created by `AuxWindow::make(...)`, so `wApp` owns them - they are
  //  not our children and will not go away with us.  They hold a raw `EnergyCalTool*` back-pointer
  //  they dereference from their Accept/Apply handlers (EnergyCalAddActions.cpp, EnergyCalGraphical.cpp),
  //  so leaving them up after we die is a use-after-free waiting on a user click.
  cancelMoreActionWindow();

  //  Note: deliberately not `deleteGraphicalRecalConfirmWindow()` - that also reads a user preference
  //  through, and logs a message to, `m_interspec`.  We run when InterSpec destroys its widget
  //  children, i.e. *after* the `~InterSpec` body, so `m_interspec` is already dangling.
  if( m_graphicalRecal )
  {
    AuxWindow::deleteAuxWindow( m_graphicalRecal.get() );
    assert( !m_graphicalRecal );   //observing_ptr auto-clears
  }
}//EnergyCalTool destructor


set<string> EnergyCalTool::gamma_detectors_for_samples( const std::shared_ptr<const SpecMeas> &meas,
                                                         const std::set<int> &samples )
{
  //We want the names of just the detectors that have gamma calibration information, of the
  //  given samples.
  //  We will assume the first Measurement for a given named detector will have gamma data if
  //  any of the Measurements from that detector will.
  //  \TODO: evaluate if that is true.
  if( !meas )
    return {};
  
  const vector<string> &detnames = meas->gamma_detector_names();
  
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
}//set<string> gamma_detectors_for_samples(...)


set<string> EnergyCalTool::gammaDetectorsForDisplayedSamples( const SpecUtils::SpectrumType type )
{
  return gamma_detectors_for_samples( m_interspec->measurment(type), m_interspec->displayedSamples(type) );
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


EnergyCalImp::CalDisplay *EnergyCalTool::foregroundFitDisplay( int *menunum, Wt::WMenuItem **item )
{
  if( menunum )
    *menunum = -1;
  if( item )
    *item = nullptr;

  const vector<string> displayed
                       = m_interspec->detectorsToDisplay( SpecUtils::SpectrumType::Foreground );

  for( int menuIndex = 0; menuIndex < 3; ++menuIndex )
  {
    WMenu * const detMenu = m_detectorMenu[menuIndex];
    for( int i = 0; detMenu && (i < detMenu->count()); ++i )
    {
      WMenuItem * const thisitem = detMenu->itemAt(i);
      auto display = thisitem ? dynamic_cast<EnergyCalImp::CalDisplay *>( thisitem->contents() )
                              : nullptr;
      if( display
          && (display->spectrumType() == SpecUtils::SpectrumType::Foreground)
          && (std::find(begin(displayed), end(displayed), display->detectorName()) != end(displayed)) )
      {
        if( menunum )
          *menunum = menuIndex;
        if( item )
          *item = thisitem;
        return display;
      }//if( we found a displayed foreground detector )
    }//for( loop over the menus items )
  }//for( loop over the detector menus )

  return nullptr;
}//EnergyCalImp::CalDisplay *foregroundFitDisplay( int *menunum, Wt::WMenuItem **item )





std::array<bool,EnergyCalTool::NumApplyToCbIndex> EnergyCalTool::applicable_apply_to_rows(
                                                              const ApplyToDisplayState &state )
{
  std::array<bool,NumApplyToCbIndex> answer;
  answer.fill( false );
  
  if( !state.meas[0] )
    return answer;  //No calibrations are shown (or editable) without a foreground
  
  answer[ApplyToForeground] = (state.meas[1] || state.meas[2]);
  answer[ApplyToBackground] = !!state.meas[1];
  answer[ApplyToSecondary] = !!state.meas[2];
  
  // "Displayed/All Detectors" only matter if some file has a gamma detector that isnt displayed;
  //  "Indep. Dets." only if some file has more than one gamma detector (i.e., when the
  //  "Detector" column lists detector names); and "Displayed/All Samples" only if some file has
  //  samples that none of the spectrum types its displayed as, are displaying.
  bool displayingAllDets = true, multipleDets = false;
  map<shared_ptr<const SpecMeas>,set<int>> undisplayed;
  
  for( int i = 0; i < 3; ++i )
  {
    const shared_ptr<const SpecMeas> meas = state.meas[i];
    if( !meas )
      continue;
    
    const vector<string> &displayed = state.detectors[i];
    for( const string &name : meas->gamma_detector_names() )
    {
      if( std::find( begin(displayed), end(displayed), name ) == end(displayed) )
        displayingAllDets = false;
    }
    
    multipleDets = (multipleDets || (gamma_detectors_for_samples(meas, state.samples[i]).size() > 1));
    
    if( !undisplayed.count(meas) )
      undisplayed[meas] = meas->sample_numbers();
    
    set<int> &undispsamples = undisplayed[meas];
    for( const int sample : state.samples[i] )
      undispsamples.erase( sample );
  }//for( loop over the types of spectrum files )
  
  bool displayingAllSamples = true;
  for( const auto &p : undisplayed )
    displayingAllSamples = (displayingAllSamples && p.second.empty());
  
  answer[ApplyToDisplayedDetectors] = answer[ApplyToAllDetectors] = !displayingAllDets;
  answer[ApplyToSelectedDetOnly] = multipleDets;
  answer[ApplyToDisplayedSamples] = answer[ApplyToAllSamples] = !displayingAllSamples;
  
  return answer;
}//applicable_apply_to_rows(...)


vector<MeasToApplyCoefChangeTo> EnergyCalTool::apply_to_scope( const ApplyToSelection &selection,
                                                               const ApplyToDisplayState &state )
{
  const std::array<bool,NumApplyToCbIndex> applicable = applicable_apply_to_rows( state );
  
  // A checkbox that isnt applicable (and so isnt shown) is interpreted as being a no-op.
  const bool type_selected[3] = {
    (!applicable[ApplyToForeground] || selection.checked[ApplyToForeground]),
    (!applicable[ApplyToBackground] || selection.checked[ApplyToBackground]),
    (!applicable[ApplyToSecondary] || selection.checked[ApplyToSecondary])
  };
  
  // The Displayed/All checkbox pairs are kept mutually exclusive by #applyToCbChanged; if they are
  //  ever inconsistent, the detectors fall back to all detectors, and the samples are decided by
  //  "All Samples" alone (as was previously done).
  const bool only_displayed_dets = applicable[ApplyToDisplayedDetectors]
                                   && selection.checked[ApplyToDisplayedDetectors]
                                   && !selection.checked[ApplyToAllDetectors];
  const bool only_selected_det = applicable[ApplyToSelectedDetOnly]
                                 && selection.checked[ApplyToSelectedDetOnly];
  const bool only_displayed_samples = applicable[ApplyToAllSamples]
                                      && !selection.checked[ApplyToAllSamples];
  
  vector<MeasToApplyCoefChangeTo> answer;
  
  for( int i = 0; i < 3; ++i )
  {
    const shared_ptr<SpecMeas> &meas = state.meas[i];
    if( !meas || !type_selected[i] )
      continue;
    
    //It could be that the background SpecFile is the same as the foreground, so check if we
    //  already have an entry for this file in answer, and if so, use it.  We dont want duplicate
    //  entries for a SpecFile, since consumers assume each measurement is visited only once (they
    //  would otherwise move peaks multiple times, or fail to look up already-updated calibrations).
    //  Note: samples displayed by a spectrum type that is NOT selected get subtracted back out
    //  below, after this loop; the remaining compromise is that *detector* scopes of multiple
    //  selected types showing the same file get unioned (but since InterSpec uses a single
    //  detector selection for all spectrum types of a file, they are the same anyway).
    MeasToApplyCoefChangeTo *changes = nullptr;
    for( size_t j = 0; !changes && (j < answer.size()); ++j )
      changes = (answer[j].meas == meas) ? &(answer[j]) : nullptr;
    
    if( !changes )
    {
      answer.emplace_back();
      changes = &(answer.back());
      changes->meas = meas;
    }
    
    // Note: we use all the gamma detectors of the file (not just the ones with data in the
    //  displayed samples), so "All Samples" really does get every gamma measurement.
    const vector<string> &gamma_dets = meas->gamma_detector_names();
    const auto is_gamma_det = [&gamma_dets]( const string &name ) -> bool {
      return (std::find( begin(gamma_dets), end(gamma_dets), name ) != end(gamma_dets));
    };
    
    if( only_selected_det )
    {
      // Matched by name, so the same detector in the other selected files gets the change too
      if( is_gamma_det(selection.selected_detector) )
        changes->detectors.insert( selection.selected_detector );
    }else if( only_displayed_dets )
    {
      for( const string &det : state.detectors[i] )
      {
        if( is_gamma_det(det) )
          changes->detectors.insert( det );
      }
    }else
    {
      changes->detectors.insert( begin(gamma_dets), end(gamma_dets) );
    }//if( only_selected_det ) / else if( only_displayed_dets ) / else
    
    if( only_displayed_samples )
      changes->sample_numbers.insert( begin(state.samples[i]), end(state.samples[i]) );
    else
      changes->sample_numbers = meas->sample_numbers();
  }//for( loop over spectrum types )
  
  
  // Dont let a change bleed into the display of an un-selected spectrum type, when the same file
  //  is displayed as multiple types.  E.g., with the foreground and background being different
  //  sample numbers of the same file, and only "Background" selected, the samples the foreground
  //  is displaying should be left alone (even when the sample scope is "all samples").
  //  Samples displayed by both a selected and an un-selected type are kept - the change is to the
  //  underlying measurements, so the un-selected display cant help but follow.
  for( MeasToApplyCoefChangeTo &change : answer )
  {
    set<int> checked_disp, unchecked_disp;
    for( int i = 0; i < 3; ++i )
    {
      if( state.meas[i] != change.meas )
        continue;
      
      const set<int> &disp = state.samples[i];
      set<int> &dest = type_selected[i] ? checked_disp : unchecked_disp;
      dest.insert( begin(disp), end(disp) );
    }//for( int i = 0; i < 3; ++i )
    
    for( const int sample : unchecked_disp )
    {
      if( !checked_disp.count(sample) )
        change.sample_numbers.erase( sample );
    }
  }//for( MeasToApplyCoefChangeTo &change : answer )
  
  // Drop files left with nothing to change (e.g., "Indep. Dets.", and that file doesnt have
  //  the selected detector)
  answer.erase( std::remove_if( begin(answer), end(answer), []( const MeasToApplyCoefChangeTo &c ){
    return c.detectors.empty() || c.sample_numbers.empty();
  } ), end(answer) );
  
  return answer;
}//apply_to_scope(...)


vector<MeasToApplyCoefChangeTo> EnergyCalTool::deviation_pair_edit_scope( ApplyToSelection selection,
                                                    const ApplyToDisplayState &state,
                                                    const std::shared_ptr<const SpecMeas> &edited_meas,
                                                    const std::string &edited_det )
{
  // Deviation pairs are usually unique to a detector, so the detector checkboxes dont apply
  selection.checked[ApplyToCbIndex::ApplyToDisplayedDetectors] = false;
  selection.checked[ApplyToCbIndex::ApplyToAllDetectors] = true;
  selection.checked[ApplyToCbIndex::ApplyToSelectedDetOnly] = false;
  
  vector<MeasToApplyCoefChangeTo> answer = apply_to_scope( selection, state );
  
  // When files have only a single detector, the "Detector" column shows the spectrum type (not the
  //  detector names), and the detectors are equivalent, regardless of their names
  const bool edited_single_det = (edited_meas && (edited_meas->gamma_detector_names().size() == 1));
  
  for( MeasToApplyCoefChangeTo &change : answer )
  {
    const vector<string> &gamma_dets = change.meas->gamma_detector_names();
    if( edited_single_det && (gamma_dets.size() == 1) )
      continue;
    
    const bool has_det = change.detectors.count( edited_det );
    change.detectors.clear();
    if( has_det )
      change.detectors.insert( edited_det );
  }//for( MeasToApplyCoefChangeTo &change : answer )
  
  answer.erase( std::remove_if( begin(answer), end(answer), []( const MeasToApplyCoefChangeTo &c ){
    return c.detectors.empty();
  } ), end(answer) );
  
  return answer;
}//deviation_pair_edit_scope(...)


bool EnergyCalTool::scope_includes_displayed_foreground( const vector<MeasToApplyCoefChangeTo> &changemeas,
                                                         const ApplyToDisplayState &state )
{
  const shared_ptr<SpecMeas> &fore = state.meas[0];
  if( !fore )
    return false;
  
  const vector<string> &disp_dets = state.detectors[0];
  
  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    if( change.meas != fore )
      continue;
    
    for( const int sample : state.samples[0] )
    {
      if( !change.sample_numbers.count(sample) )
        continue;
      
      for( const string &det : change.detectors )
      {
        if( std::find( begin(disp_dets), end(disp_dets), det ) == end(disp_dets) )
          continue;
        
        const shared_ptr<const SpecUtils::Measurement> m = fore->measurement( sample, det );
        if( m && (m->num_gamma_channels() > 4) )
          return true;
      }//for( const string &det : change.detectors )
    }//for( const int sample : state.samples[0] )
  }//for( const MeasToApplyCoefChangeTo &change : changemeas )
  
  return false;
}//scope_includes_displayed_foreground(...)


EnergyCalTool::ApplyToDisplayState EnergyCalTool::currentApplyToDisplayState() const
{
  const SpecUtils::SpectrumType spectypes[3] = {
    SpecUtils::SpectrumType::Foreground,
    SpecUtils::SpectrumType::Background,
    SpecUtils::SpectrumType::SecondForeground
  };
  
  ApplyToDisplayState state;
  for( int i = 0; i < 3; ++i )
  {
    state.meas[i] = m_interspec->measurment( spectypes[i] );
    if( state.meas[i] )
    {
      state.samples[i] = m_interspec->displayedSamples( spectypes[i] );
      state.detectors[i] = m_interspec->detectorsToDisplay( spectypes[i] );
    }
  }//for( int i = 0; i < 3; ++i )
  
  return state;
}//ApplyToDisplayState currentApplyToDisplayState() const


EnergyCalTool::ApplyToSelection EnergyCalTool::currentApplyToSelection( const std::string *edited_det ) const
{
  ApplyToSelection selection;
  for( int i = 0; i < NumApplyToCbIndex; ++i )
  {
    if( m_applyToCbs[i] )
      selection.checked[i] = m_applyToCbs[i]->isChecked();
  }
  
  if( edited_det )
  {
    selection.selected_detector = *edited_det;
  }else
  {
    const EnergyCalImp::CalDisplay *display = m_calInfoDisplayStack
          ? dynamic_cast<const EnergyCalImp::CalDisplay *>( m_calInfoDisplayStack->currentWidget() )
          : nullptr;
    
    if( display )
    {
      selection.selected_detector = display->detectorName();
    }else
    {
      // The displays get created at render time; until then, use the detector that
      //  #doRefreshFromFiles will initially select.
      const set<string> dets = gamma_detectors_for_samples(
                           m_interspec->measurment( SpecUtils::SpectrumType::Foreground ),
                           m_interspec->displayedSamples( SpecUtils::SpectrumType::Foreground ) );
      if( !dets.empty() )
        selection.selected_detector = *begin(dets);
    }//if( display ) / else
  }//if( edited_det ) / else
  
  return selection;
}//ApplyToSelection currentApplyToSelection( const std::string *edited_det ) const


vector<MeasToApplyCoefChangeTo> EnergyCalTool::measurementsToApplyCoeffChangeTo( const std::string *edited_det )
{
  return apply_to_scope( currentApplyToSelection(edited_det), currentApplyToDisplayState() );
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


namespace
{
  /** The changes to set `new_cal` as the calibration of every gamma measurement in `change`;
   throws if any of them have a different number of channels than `new_cal`.
   */
  EnergyCalTool::MeasToNewCalMap absolute_cal_changes( const MeasToApplyCoefChangeTo &change,
                                    const shared_ptr<const SpecUtils::EnergyCalibration> &new_cal )
  {
    assert( change.meas && new_cal && new_cal->valid() );
    
    EnergyCalTool::MeasToNewCalMap answer;
    try
    {
      for( const int sample : change.sample_numbers )
      {
        for( const string &detname : change.detectors )
        {
          const shared_ptr<const SpecUtils::Measurement> m = change.meas->measurement( sample, detname );
          const shared_ptr<const SpecUtils::EnergyCalibration> old_cal = m ? m->energy_calibration() : nullptr;
          if( !m || (m->num_gamma_channels() <= 4) || !old_cal || !old_cal->valid() )
            continue;
          
          // TODO: maybe be a little more lenient on matching number of channels, by possible making
          //       new energy calibrations when possible.
          if( old_cal->num_channels() != new_cal->num_channels() )
          {
            throw runtime_error( WString::tr("ect-set-cal-invalid-nchan")
                                .arg(sample)
                                .arg(detname)
                                .arg( static_cast<int>(old_cal->num_channels()) )
                                .arg( static_cast<int>(new_cal->num_channels()) )
                                .toUTF8() );
          }
          
          answer[m] = new_cal;
        }//for( const string &detname : change.detectors )
      }//for( loop over sample numbers )
    }catch( std::exception &e )
    {
      throw runtime_error( WString::tr("ect-error-setting-cal").arg(e.what()).toUTF8() );
    }//try catch
    
    return answer;
  }//absolute_cal_changes(...)
}//namespace


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
      applyNewCalibrations( {tochange}, absolute_cal_changes( tochange, new_cal ) );
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
        applyNewCalibrations( {tochange}, absolute_cal_changes( tochange, new_cal ) );
      }
    }//if( we are changing a single detector ) / else
  }else
  {
    // We will apply the input calibration on a detector-by-detector basis, requiring the
    //  input to have info for every relevant detector.  All the detectors are changed together, so
    //  the peaks get translated (by #applyNewCalibrations) only once, according to the new
    //  calibration of the detector they are displayed with.
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
    
    MeasToNewCalMap new_cals;
    for( const string &det : tochange.detectors )
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
      
      const MeasToNewCalMap det_cals = absolute_cal_changes( det_specific_tochange, new_cal );
      new_cals.insert( begin(det_cals), end(det_cals) );
    }//for( const string &det : tochange.detectors )
    
    applyNewCalibrations( {tochange}, new_cals );
  }//if( det_to_cal.size() == 1 ) / else
  
  // Note: #applyNewCalibrations refreshed the charts and this tool.
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


std::shared_ptr<EnergyCalImp::CALpDownloadResource> EnergyCalTool::calpResources()
{
  return m_calpResource;
}



void EnergyCalTool::handleRequestToUploadCALp()
{
  SimpleDialog *dialog = SimpleDialog::make<SimpleDialog>();
  WPushButton *closeButton = dialog->addButton( "Cancel", WidgetUtils::ButtonRole::Dismiss );
  WGridLayout *stretcher = dialog->contents()->setLayout( std::make_unique<WGridLayout>() );
  stretcher->setContentsMargins( 0, 0, 0, 0 );
  dialog->contents()->setOverflow( Overflow::Visible );
  WText *title = stretcher->addWidget( std::make_unique<WText>(WString::tr("ect-import-CALp")), 0, 0 );
  title->addStyleClass( "title" );
  
  WText *t = stretcher->addWidget( std::make_unique<WText>(WString::tr("ect-select-CALp")), stretcher->rowCount(), 0, AlignmentFlag::Center | AlignmentFlag::Middle );
  t->setTextAlignment( Wt::AlignmentFlag::Center );
  
  
  WFileUpload *upload = stretcher->addWidget( std::make_unique<WFileUpload>(), stretcher->rowCount(), 0, AlignmentFlag::Center | AlignmentFlag::Middle );
  upload->fileTooLarge().connect( upload, [=](){
    dialog->contents()->clear();
    dialog->footer()->clear();
    
    WPushButton *closeButton = dialog->addButton( WString::tr("Close"), WidgetUtils::ButtonRole::Dismiss );
    WGridLayout *stretcher = dialog->contents()->setLayout( std::make_unique<WGridLayout>() );
    stretcher->setContentsMargins( 0, 0, 0, 0 );
    WText *title = stretcher->addWidget( std::make_unique<WText>(WString::tr("ect-upload-CALp-to-large")), 0, 0 );
    title->addStyleClass( "title" );
  } );
  
  upload->changed().connect( upload, &WFileUpload::upload );
  upload->uploaded().connect( upload, [dialog,upload](){
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

      WPushButton *closeButton = dialog->addButton( WString::tr("Close"), WidgetUtils::ButtonRole::Dismiss );
      WGridLayout *stretcher = dialog->contents()->setLayout( std::make_unique<WGridLayout>() );
      stretcher->setContentsMargins( 0, 0, 0, 0 );
      WText *title = stretcher->addWidget( std::make_unique<WText>(WString::tr("ect-upload-CALp-to-large")), 0, 0 );
      title->addStyleClass( "title" );

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
      
      WPushButton *closeButton = dialog->addButton( WString::tr("Close"), WidgetUtils::ButtonRole::Dismiss );
      WGridLayout *stretcher = dialog->contents()->setLayout( std::make_unique<WGridLayout>() );
      stretcher->setContentsMargins( 0, 0, 0, 0 );
      WText *title = stretcher->addWidget( std::make_unique<WText>(WString::tr("ect-invalid-CALp")), 0, 0 );
      title->addStyleClass( "title" );
      
      return;
    }//if( was not a valid CALp file )
    
    //dialog->done( Wt::WDialog::DialogCode::Accepted );
  } );
  
  
  InterSpec *interspec = InterSpec::instance();
  if( interspec && !interspec->isPhone() )
  {
    auto tOwned = std::make_unique<WText>(WString::tr("ect-CALp-drag-n-drop-note"));
    t = tOwned.get();
    stretcher->addWidget( std::move(tOwned), stretcher->rowCount(), 0, AlignmentFlag::Center | AlignmentFlag::Middle );
    t->setTextAlignment( Wt::AlignmentFlag::Center );
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



namespace
{
  /** Merges entries for the same SpecMeas (taking the union of their samples and detectors), since
   the calibration change functions assume each measurement is only visited once.

   Note: two entries for a file with different samples AND different detectors will then also cover
   the cross combinations; no current caller passes such entries (they are either a scope from
   #EnergyCalTool::apply_to_scope, which has one entry per file, or same-detector entries).
   */
  vector<MeasToApplyCoefChangeTo> merge_change_entries( const vector<MeasToApplyCoefChangeTo> &input )
  {
    vector<MeasToApplyCoefChangeTo> answer;
    for( const MeasToApplyCoefChangeTo &change : input )
    {
      if( !change.meas )
        continue;
      
      MeasToApplyCoefChangeTo *existing = nullptr;
      for( size_t i = 0; !existing && (i < answer.size()); ++i )
        existing = (answer[i].meas == change.meas) ? &(answer[i]) : nullptr;
      
      if( existing )
      {
        existing->sample_numbers.insert( begin(change.sample_numbers), end(change.sample_numbers) );
        existing->detectors.insert( begin(change.detectors), end(change.detectors) );
      }else
      {
        answer.push_back( change );
      }
    }//for( const MeasToApplyCoefChangeTo &change : input )
    
    return answer;
  }//merge_change_entries(...)
  
  
  /** Message context of which spectrum had an error, if it is necessary to say. */
  WString spectrum_desc_for_error( InterSpec *viewer, const shared_ptr<SpecMeas> &meas )
  {
    const auto forgrnd = viewer->measurment(SpecUtils::SpectrumType::Foreground);
    const auto backgrnd = viewer->measurment(SpecUtils::SpectrumType::Background);
    const auto secgrnd = viewer->measurment(SpecUtils::SpectrumType::SecondForeground);
    
    if( (backgrnd && (backgrnd != meas)) || (secgrnd && (secgrnd != meas)) )
    {
      if( meas == forgrnd )
        return WString::tr("ect-for-the-fore");
      if( meas == backgrnd )
        return WString::tr("ect-for-the-back");
      if( meas == secgrnd )
        return WString::tr("ect-for-the-sec");
    }//if( it is necessary to say which spectrum had the error )
    
    return WString();
  }//spectrum_desc_for_error(...)
  
  
  /** Returns a new calibration with the same type, channel count, and coefficients as `cal`, but
   with the given deviation pairs.
   */
  shared_ptr<const SpecUtils::EnergyCalibration> with_dev_pairs(
                                      const shared_ptr<const SpecUtils::EnergyCalibration> &cal,
                                      const vector<pair<float,float>> &dev_pairs )
  {
    assert( cal && cal->valid() );
    
    auto new_cal = make_shared<SpecUtils::EnergyCalibration>();
    switch( cal->type() )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        new_cal->set_polynomial( cal->num_channels(), cal->coefficients(), dev_pairs );
        break;
        
      case SpecUtils::EnergyCalType::FullRangeFraction:
        new_cal->set_full_range_fraction( cal->num_channels(), cal->coefficients(), dev_pairs );
        break;
        
      case SpecUtils::EnergyCalType::LowerChannelEdge:
      case SpecUtils::EnergyCalType::InvalidEquationType:
        throw runtime_error( "Deviation pairs can not be set for this energy calibration type." );
    }//switch( cal->type() )
    
    assert( new_cal->valid() );
    return new_cal;
  }//with_dev_pairs(...)
}//namespace


void EnergyCalTool::applyNewCalibrations( const vector<MeasToApplyCoefChangeTo> &input_changemeas,
                                          const MeasToNewCalMap &new_cals )
{
  using namespace SpecUtils;
  
  const vector<MeasToApplyCoefChangeTo> changemeas = merge_change_entries( input_changemeas );
  
  const SpecUtils::SpectrumType spectypes[3] = {
    SpecUtils::SpectrumType::Foreground,
    SpecUtils::SpectrumType::Background,
    SpecUtils::SpectrumType::SecondForeground
  };
  
  shared_ptr<SpecMeas> specfiles[3];
  for( int i = 0; i < 3; ++i )
    specfiles[i] = m_interspec->measurment( spectypes[i] );
  
  // Calls `f( meas, m, new_cal )` for each measurement in `changemeas` that is getting a new cal
  const auto for_each_change = [&changemeas,&new_cals]( const std::function<void(const shared_ptr<SpecMeas> &,
                                                              const shared_ptr<const Measurement> &,
                                                              const shared_ptr<const EnergyCalibration> &)> &f ){
    for( const MeasToApplyCoefChangeTo &change : changemeas )
    {
      for( const int sample : change.sample_numbers )
      {
        for( const string &detname : change.detectors )
        {
          const shared_ptr<const Measurement> m = change.meas->measurement( sample, detname );
          const auto pos = m ? new_cals.find( m ) : end(new_cals);
          if( (pos != end(new_cals)) && pos->second && (pos->second != m->energy_calibration()) )
            f( change.meas, m, pos->second );
        }//for( loop over detector names )
      }//for( loop over sample numbers )
    }//for( const MeasToApplyCoefChangeTo &change : changemeas )
  };//for_each_change
  
  size_t num_to_change = 0;
  for_each_change( [&num_to_change]( const shared_ptr<SpecMeas> &, const shared_ptr<const Measurement> &m,
                                     const shared_ptr<const EnergyCalibration> &cal ){
    assert( cal->valid() && (cal->num_channels() == m->num_gamma_channels()) );
    if( !cal->valid() || (cal->num_channels() != m->num_gamma_channels()) )
      throw runtime_error( "Invalid energy calibration for '" + m->detector_name() + "' sample "
                           + std::to_string(m->sample_number()) );
    ++num_to_change;
  } );
  
  if( !num_to_change )
    throw runtime_error( WString::tr("ect-changed-cal-not-selected").toUTF8() );
  
  
  // Translate the peaks, but we wont actually update them to the SpecMeas until we know we can
  //  update all the peaks.
  typedef deque<shared_ptr<const PeakDef>> PeakDeque;
  map<shared_ptr<PeakDeque>,PeakDeque> updated_peaks;
  map<shared_ptr<const PeakDeque>,PeakDeque> updated_hint_peaks;
  
  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    // The peaks position (i.e., mean channel number) is determined by
    //  #SpecFile::suggested_sum_energy_calibration of the displayed detectors, however we may be
    //  applying an energy calibration change to one or two detectors, that arent that detector, in
    //  which case we dont want to move the peaks.  This is why we dont just use change.detectors.
    vector<string> display_detectors;
    for( int i = 0; display_detectors.empty() && (i < 3); ++i )
    {
      if( specfiles[i] == change.meas )
        display_detectors = m_interspec->detectorsToDisplay( spectypes[i] );
    }
    
    if( display_detectors.empty() )  //File not displayed (shouldnt really happen)
      display_detectors.insert( end(display_detectors), begin(change.detectors), end(change.detectors) );
    
    // Returns the {old, new} calibration to translate the peaks of `samples` with, or nullptrs if
    //  they shouldnt be translated.  The peaks were fit to the sum of the displayed detectors, so
    //  we only translate them if that whole sum moves: all of their sample numbers are being
    //  updated, and every displayed gamma measurement of those samples is getting a new calibration.
    //  When only some of the summed measurements change (e.g., "Indep. Dets.", or a
    //  deviation pair edit of one detector) there is no right answer, so the peaks are left alone.
    //  The translation is from the summed (displayed) calibration, to the new calibration of a
    //  measurement that currently has it.
    typedef pair<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> CalPair;
    const auto peak_cals = [&]( const set<int> &samples ) -> CalPair {
      for( const int sample : samples )
      {
        if( !change.sample_numbers.count(sample) )
          return {};
      }
      
      shared_ptr<const EnergyCalibration> oldcal;
      try
      {
        oldcal = change.meas->suggested_sum_energy_calibration( samples, display_detectors );
      }catch( std::exception & )
      {
      }
      
      if( !oldcal || !oldcal->valid() )
        return {};
      
      CalPair answer;
      for( const int sample : samples )
      {
        for( const string &detname : display_detectors )
        {
          const shared_ptr<const Measurement> m = change.meas->measurement( sample, detname );
          const shared_ptr<const EnergyCalibration> cal = m ? m->energy_calibration() : nullptr;
          if( !m || (m->num_gamma_channels() <= 4) || !cal || !cal->valid() )
            continue;
          
          const MeasToNewCalMap::const_iterator pos = new_cals.find( m );
          if( (pos == end(new_cals)) || !pos->second || (pos->second == cal) )
            return {};  //Part of the displayed sum isnt changing
          
          if( !answer.second && (cal == oldcal) )
            answer = { oldcal, pos->second };
        }//for( const string &detname : display_detectors )
      }//for( const int sample : samples )
      
      return answer;
    };//peak_cals lambda
    
    for( const set<int> &samples : change.meas->sampleNumsWithPeaks() )
    {
      const shared_ptr<PeakDeque> oldpeaks = change.meas->peaks( samples );
      if( !oldpeaks || oldpeaks->empty() )
        continue;
      
      const CalPair cals = peak_cals( samples );
      if( !cals.second )
        continue;
      
      try
      {
        updated_peaks[oldpeaks] = EnergyCal::translatePeaksForCalibrationChange( *oldpeaks, cals.first, cals.second );
      }catch( std::exception &e )
      {
        const string msg = "There was an issue translating peaks for this energy change;"
                           " not applying change.  Error: " + string(e.what());
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, msg.c_str() );
#endif
        throw runtime_error( msg );
      }//try / catch
    }//for( const set<int> &samples : change.meas->sampleNumsWithPeaks() )
    
    for( const set<int> &samples : change.meas->sampleNumsWithAutomatedSearchPeaks() )
    {
      const shared_ptr<const PeakDeque> oldpeaks = change.meas->automatedSearchPeaks( samples );
      if( !oldpeaks || oldpeaks->empty() )
        continue;
      
      const CalPair cals = peak_cals( samples );
      if( !cals.second )
        continue;
      
      try
      {
        updated_hint_peaks[oldpeaks] = EnergyCal::translatePeaksForCalibrationChange( *oldpeaks, cals.first, cals.second );
      }catch( std::exception &e )
      {
        // Hint peaks are not critical, so we will still apply the change
#if( PERFORM_DEVELOPER_CHECKS )
        const string msg = "There was an issue translating hint peaks for this energy change,"
                           " still applying, Error: " + string(e.what());
        log_developer_error( __func__, msg.c_str() );
#endif
      }//try / catch
    }//for( const set<int> &samples : change.meas->sampleNumsWithAutomatedSearchPeaks() )
  }//for( const MeasToApplyCoefChangeTo &change : changemeas )
  
  
  // Now actually set the energy calibrations and peaks, recording them for undo/redo.
  EnergyCalUndoRedoSentry undo_sentry;
  
  for_each_change( [&undo_sentry]( const shared_ptr<SpecMeas> &meas, const shared_ptr<const Measurement> &m,
                                   const shared_ptr<const EnergyCalibration> &cal ){
    undo_sentry.cal_info( meas ).emplace_back( m, m->energy_calibration(), cal );
    meas->set_energy_calibration( cal, m );
  } );
  
  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    for( const set<int> &samples : change.meas->sampleNumsWithPeaks() )
    {
      const shared_ptr<PeakDeque> oldpeaks = change.meas->peaks( samples );
      const auto pos = oldpeaks ? updated_peaks.find( oldpeaks ) : end(updated_peaks);
      if( pos == end(updated_peaks) )
        continue;
      
      undo_sentry.peak_info( change.meas ).emplace_back( samples, *oldpeaks, pos->second );
      change.meas->setPeaks( pos->second, samples );
      
      for( int i = 0; m_peakModel && (i < 3); ++i )
      {
        if( (change.meas == specfiles[i]) && (samples == m_interspec->displayedSamples(spectypes[i])) )
        {
          m_peakModel->setPeakFromSpecMeas( specfiles[i], samples, spectypes[i] );
          break;
        }
      }//for( loop over spectrum types to update the peak model for )
    }//for( const set<int> &samples : change.meas->sampleNumsWithPeaks() )
    
    for( const set<int> &samples : change.meas->sampleNumsWithAutomatedSearchPeaks() )
    {
      const shared_ptr<const PeakDeque> oldpeaks = change.meas->automatedSearchPeaks( samples );
      const auto pos = oldpeaks ? updated_hint_peaks.find( oldpeaks ) : end(updated_hint_peaks);
      if( (pos == end(updated_hint_peaks)) || pos->second.empty() )
        continue;
      
      undo_sentry.hint_peak_info( change.meas ).emplace_back( samples, *oldpeaks, pos->second );
      change.meas->setAutomatedSearchPeaks( samples, make_shared<PeakDeque>( pos->second ) );
    }//for( const set<int> &samples : change.meas->sampleNumsWithAutomatedSearchPeaks() )
  }//for( const MeasToApplyCoefChangeTo &change : changemeas )
  
  // Any calibration change, other than the graphical one itself (which sets this state after
  //  calling into here), invalidates preserving the last graphical calibration.
  m_lastGraphicalRecal = 0;
  m_lastGraphicalRecalType = EnergyCalGraphicalConfirm::NumRecalTypes;
  m_lastGraphicalRecalEnergy = -999.0f;
  
  m_interspec->refreshDisplayedCharts();
  refreshGuiFromFiles();
}//applyNewCalibrations(...)


void EnergyCalTool::applyCalChange( std::shared_ptr<const SpecUtils::EnergyCalibration> disp_prev_cal,
                                    std::shared_ptr<const SpecUtils::EnergyCalibration> new_disp_cal,
                                    const vector<MeasToApplyCoefChangeTo> &input_changemeas,
                                    const bool isOffsetOnly )
{
  // We get here when:
  //  - we manually change a coefficiecnt on the GUI
  //  - we fit calibration coefficients
  //  - we graphically recalibrate
  //  - sometimes when we load in a CALp file
  
  using namespace SpecUtils;
  
  assert( disp_prev_cal && disp_prev_cal->valid() );
  assert( new_disp_cal && new_disp_cal->valid() );
  
  const vector<MeasToApplyCoefChangeTo> changemeas = merge_change_entries( input_changemeas );
  
  // Create a cache of modified calibration both to save time/memory, but also keep it so previous
  //  samples that share a energy calibration will continue to do so (if possible based on what user
  //  wanted calibration applied to).
  //Note: we could take this oppritunity to share calibration across SpecFile objects by not just
  //      comparing pointers, but also the actual EnergyCalibration object.  But for now we'll
  //      skip this to avoid trouble, and it isnt clear that it would actually be overall beneficial
  map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> old_to_new_cals;
  MeasToNewCalMap new_cals;
  
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
          
          const auto meas_old_cal = m->energy_calibration();
          assert( meas_old_cal );
          
          if( !meas_old_cal || !meas_old_cal->valid() )
            continue;
          
          //If we have already computed the new calibration for a EnergyCalibration object, lets not
          //  re-due it.
          const auto prev_pos = old_to_new_cals.find( meas_old_cal );
          if( prev_pos != end(old_to_new_cals) )
          {
            new_cals[m] = prev_pos->second;
            continue;
          }
          
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
            
            // Both cals must be lower-channel-edge to trace back to a common original.  A change
            //  that switches the calibration *type* (e.g. "Guess Gain", which always applies a
            //  polynomial) has no shared original to remap through, and `lowerChannelOriginal()`
            //  asserts on anything else - so fall through to the non-remap path.
            if( (disp_prev_cal->type() == SpecUtils::EnergyCalType::LowerChannelEdge)
                && new_disp_cal
                && (new_disp_cal->type() == SpecUtils::EnergyCalType::LowerChannelEdge) )
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
          new_cals[m] = new_meas_cal;
        }//for( const string &detname : change.detectors )
      }//for( loop over sample numbers )
    }catch( std::exception &e )
    {
      WString msg = WString::tr("ect-change-made-invalid");
      msg.arg( spectrum_desc_for_error( m_interspec, change.meas ) );
      msg.arg( e.what() );
      
      throw runtime_error( msg.toUTF8() );
    }//try catch
  }//for( const MeasToApplyCoefChangeTo &change : changemeas )
  
  applyNewCalibrations( changemeas, new_cals );
}//applyCalChange(...)


void EnergyCalTool::addDeviationPair( const std::pair<float,float> &new_pair )
{
  addDeviationPair( new_pair, measurementsToApplyCoeffChangeTo() );
}//void addDeviationPair( const std::pair<float,float> &new_pair )


void EnergyCalTool::addDeviationPair( const std::pair<float,float> &new_pair,
                                      const vector<MeasToApplyCoefChangeTo> &input_changemeas )
{
  // We get here when we graphically (i.e., cntrl+option+drag) add in a deviation pair; the pair
  //  is added to each measurements own deviation pairs.
  using namespace SpecUtils;
  
  const vector<MeasToApplyCoefChangeTo> changemeas = merge_change_entries( input_changemeas );
  
  map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> old_to_new_cals;
  MeasToNewCalMap new_cals;
  
  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    assert( change.meas );
    
    try
    {
      for( const int sample : change.sample_numbers )
      {
        for( const string &detname : change.detectors )
        {
          const shared_ptr<const Measurement> m = change.meas->measurement( sample, detname );
          const shared_ptr<const EnergyCalibration> old_cal = m ? m->energy_calibration() : nullptr;
          
          if( !m || (m->num_gamma_channels() <= 4) || !old_cal || !old_cal->valid()
             || (old_cal->type() == EnergyCalType::LowerChannelEdge) )
            continue;
          
          shared_ptr<const EnergyCalibration> &new_cal = old_to_new_cals[old_cal];
          if( !new_cal )
          {
            vector<pair<float,float>> dev_pairs = old_cal->deviation_pairs();
            if( dev_pairs.empty() )
              dev_pairs.push_back( {0.0f, 0.0f} );
            dev_pairs.push_back( new_pair );
            std::sort( begin(dev_pairs), end(dev_pairs) );
            
            new_cal = with_dev_pairs( old_cal, dev_pairs );
          }//if( !new_cal )
          
          new_cals[m] = new_cal;
        }//for( const string &detname : change.detectors )
      }//for( loop over sample numbers )
    }catch( std::exception &e )
    {
      WString msg = WString::tr("ect-add-dev-pair-made-invalid");
      msg.arg( spectrum_desc_for_error( m_interspec, change.meas ) );
      throw runtime_error( msg.arg(e.what()).toUTF8() );
    }//try catch
  }//for( const MeasToApplyCoefChangeTo &change : changemeas )
  
  applyNewCalibrations( changemeas, new_cals );
}//void addDeviationPair( const std::pair<float,float> &new_pair, changemeas );


void EnergyCalTool::applyDeviationPairEdit( const vector<pair<float,float>> &new_dev_pairs,
                                            const vector<MeasToApplyCoefChangeTo> &input_changemeas )
{
  using namespace SpecUtils;
  
  const vector<MeasToApplyCoefChangeTo> changemeas = merge_change_entries( input_changemeas );
  
  map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> old_to_new_cals;
  MeasToNewCalMap new_cals;
  
  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    for( const int sample : change.sample_numbers )
    {
      for( const string &name : change.detectors )
      {
        const shared_ptr<const Measurement> m = change.meas->measurement( sample, name );
        const shared_ptr<const EnergyCalibration> cal = m ? m->energy_calibration() : nullptr;
        if( !m || (m->num_gamma_channels() <= 4) || !cal || !cal->valid()
            || (cal->type() == EnergyCalType::LowerChannelEdge) )
          continue;
        
        shared_ptr<const EnergyCalibration> &new_cal = old_to_new_cals[cal];
        if( !new_cal )
          new_cal = with_dev_pairs( cal, new_dev_pairs );
        new_cals[m] = new_cal;
      }//for( const string &name : change.detectors )
    }//for( const int sample : change.sample_numbers )
  }//for( const MeasToApplyCoefChangeTo &change : changemeas )
  
  applyNewCalibrations( changemeas, new_cals );
}//applyDeviationPairEdit(...)


bool EnergyCalTool::applyingToForeground() const
{
  // Note: equivalent to `applicable_apply_to_rows(...)[ApplyToForeground]`, but this gets called
  //  often (e.g., on every peak model change), so we avoid the work of computing every row
  if( !m_interspec->measurment( SpecUtils::SpectrumType::Foreground ) )
    return false;
  
  const bool haveOther = m_interspec->measurment( SpecUtils::SpectrumType::Background )
                         || m_interspec->measurment( SpecUtils::SpectrumType::SecondForeground );
  return !haveOther || m_applyToCbs[ApplyToCbIndex::ApplyToForeground]->isChecked();
}//bool applyingToForeground() const


bool EnergyCalTool::scopeIncludesDisplayedForeground( const vector<MeasToApplyCoefChangeTo> &changemeas ) const
{
  return scope_includes_displayed_foreground( changemeas, currentApplyToDisplayState() );
}//bool scopeIncludesDisplayedForeground(...)


bool EnergyCalTool::scopeIncludesDisplay( const vector<MeasToApplyCoefChangeTo> &changemeas,
                                          const EnergyCalImp::CalDisplay *display ) const
{
  if( !display )
    return false;
  
  const SpecUtils::SpectrumType type = display->spectrumType();
  const string &detname = display->detectorName();
  const shared_ptr<SpecMeas> meas = m_interspec->measurment( type );
  const shared_ptr<const SpecUtils::EnergyCalibration> cal = display->lastSetCalibration();
  if( !meas || !cal )
    return false;
  
  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    if( (change.meas != meas) || !change.detectors.count(detname) )
      continue;
    
    // The display shows the calibration of the first displayed sample with gamma data, for the
    //  detector, so require that one of the displayed samples, with that calibration, is included
    for( const int sample : m_interspec->displayedSamples(type) )
    {
      if( !change.sample_numbers.count(sample) )
        continue;
      
      const shared_ptr<const SpecUtils::Measurement> m = meas->measurement( sample, detname );
      if( m && (m->energy_calibration() == cal) )
        return true;
    }//for( const int sample : m_interspec->displayedSamples(type) )
  }//for( const MeasToApplyCoefChangeTo &change : changemeas )
  
  return false;
}//bool scopeIncludesDisplay(...)


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
  
  // Make sure the changed energy cal is actually selected for the change to be applied to
  const vector<MeasToApplyCoefChangeTo> changemeas
                                    = measurementsToApplyCoeffChangeTo( &display->detectorName() );
  if( !scopeIncludesDisplay( changemeas, display ) )
  {
    m_interspec->logMessage( WString::tr("ect-changed-cal-not-selected"), 2 );
    doRefreshFromFiles();
    return;
  }//if( the edited display wont have the change applied to it )
  
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
  if( !dispmeas || !old_cal )  //Shouldnt ever happen
  {
    display->updateToGui( old_cal );
    m_interspec->logMessage( "Internal error retrieving correct measurement", 2 );
    return;
  }

  // Deviation pairs are usually unique to a detector, so an edit is only applied to the edited
  //  detector (matched by name in the other selected spectra), regardless of the detector choices
  //  under "Apply Changes To"; the spectrum and sample choices still apply.  Each selected
  //  measurements deviation pairs are REPLACED with the displayed list, keeping that measurements
  //  own coefficients - a full-list edit has no well-defined per-measurement delta.
  const vector<MeasToApplyCoefChangeTo> changemeas
                          = deviation_pair_edit_scope( currentApplyToSelection( &detname ),
                                                       currentApplyToDisplayState(), dispmeas, detname );

  if( !scopeIncludesDisplay( changemeas, display ) )
  {
    display->updateToGui( old_cal );
    m_interspec->logMessage( WString::tr("ect-changed-cal-not-selected"), 2 );
    return;
  }

  m_lastGraphicalRecal = 0;
  m_lastGraphicalRecalType = EnergyCalGraphicalConfirm::NumRecalTypes;
  m_lastGraphicalRecalEnergy = -999.0f;

  try
  {
    applyDeviationPairEdit( new_dev_pairs, changemeas );
  }catch( std::exception &e )
  {
    display->updateToGui( old_cal );
    m_interspec->logMessage( WString::tr("ect-dev-pair-change-invalid").arg(e.what()), 2 );
    return;
  }//try / catch

  display->setDeviationPairsValid();

  // Let the user know if the change may have applied beyond just the displayed spectrum
  /// \TODO: keep from issuing this message for ever single change!
  if( (dispmeas->sample_numbers().size() > 1) || (changemeas.size() > 1) )
    m_interspec->logMessage( WString::tr("ect-dev-pairs-applied-scope").arg(detname), 1 );
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
  
  closeStaleMoreActionWindow();
  
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


EnergyCalTool::CanFitCoefStatus EnergyCalTool::canDoEnergyFit()
{
  // Check that "Apply Changes To" for the "Foreground" is checked, otherwise fitting makes no sense
  if( !applyingToForeground() )
    return CanFitCoefStatus::ForegroundNotApplied;

  // Check if there are any peaks currently showing.
  shared_ptr<const deque<PeakModel::PeakShrdPtr>> peaks = m_peakModel->peaks();
  if( !peaks )
    return CanFitCoefStatus::NoPeaksToUse;

  size_t nPeaksToUse = 0;
  for( const PeakModel::PeakShrdPtr &p : *peaks )
    nPeaksToUse += (p && p->useForEnergyCalibration());

  if( nPeaksToUse < 1 )
    return CanFitCoefStatus::NoPeaksToUse;

  if( !m_calInfoDisplayStack )
    return CanFitCoefStatus::InternalError;

  auto caldisp = dynamic_cast<EnergyCalImp::CalDisplay *>( m_calInfoDisplayStack->currentWidget() );
  if( !caldisp )
    return CanFitCoefStatus::InternalError;

  auto cal = caldisp->lastSetCalibration();
  if( !cal || !cal->valid() )
    return CanFitCoefStatus::InternalError;

  // The peaks used for the fit are the foregrounds peaks, so the calibration being shown must
  //  belong to the foreground - or at least be the foregrounds displayed calibration (e.g., when
  //  the same file is loaded as multiple spectrum types).
  if( caldisp->spectrumType() != SpecUtils::SpectrumType::Foreground )
  {
    const shared_ptr<const SpecUtils::Measurement> forehist
                           = m_interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
    if( !forehist || (forehist->energy_calibration() != cal) )
      return CanFitCoefStatus::NotForegroundCal;
  }//if( the display being shown isnt for the foreground )

  switch( cal->type() )
  {
    case SpecUtils::EnergyCalType::InvalidEquationType:
      return CanFitCoefStatus::InvalidCalType;

    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::FullRangeFraction:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
    case SpecUtils::EnergyCalType::LowerChannelEdge:  //fits {offset, gain} relative to original
      break;
  }//switch( cal->type() )

  // The "Fit" checkboxes actually used for the fit are the foreground fit-source displays (see
  //  #fitCoefficients), which is `caldisp` itself, unless we are showing a non-foreground display
  //  that shares the foregrounds calibration.
  const vector<string> foredets
                       = m_interspec->detectorsToDisplay( SpecUtils::SpectrumType::Foreground );
  if( (caldisp->spectrumType() != SpecUtils::SpectrumType::Foreground)
      || (std::find(begin(foredets), end(foredets), caldisp->detectorName()) == end(foredets)) )
  {
    caldisp = foregroundFitDisplay();
    if( !caldisp )
      return CanFitCoefStatus::InternalError;
  }//if( the shown display isnt a displayed foreground detector )

  const set<size_t> ordersToFit = caldisp->fitForCoefficents();
  if( ordersToFit.empty() )
    return CanFitCoefStatus::NoCoefSelected;

  if( ordersToFit.size() > nPeaksToUse )
    return CanFitCoefStatus::MorePeaksNeeded;

  return CanFitCoefStatus::CanFit;
}//CanFitCoefStatus canDoEnergyFit()


void EnergyCalTool::fitCoefficients()
{
  try
  {
    if( canDoEnergyFit() != CanFitCoefStatus::CanFit )
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
      int menunum = -1;
      WMenuItem *item = nullptr;
      EnergyCalImp::CalDisplay * const foredisp = foregroundFitDisplay( &menunum, &item );

      if( !foredisp )
        throw runtime_error( WString::tr("ect-select-cal-of-disp-det").toUTF8() );

      // Select the found display, so the user sees the calibration that got fit
      m_specTypeMenu->select( menunum );
      item->select();

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

    const vector<MeasToApplyCoefChangeTo> changemeas = measurementsToApplyCoeffChangeTo( &detname );
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

    WString msg = WString::tr("ect-fit-successful");
    try
    {
      // Summarize how much the fit improved things: the mean absolute offset of the used peaks
      //  from their true energies, before and after the fit.
      double pre_dev = 0.0, post_dev = 0.0;
      for( const EnergyCal::RecalPeakInfo &info : peakInfos )
      {
        pre_dev += fabs( info.peakMean - info.photopeakEnergy );
        post_dev += fabs( answer->energy_for_channel(info.peakMeanBinNumber) - info.photopeakEnergy );
      }
      pre_dev /= peakInfos.size();
      post_dev /= peakInfos.size();

      char prebuf[32], postbuf[32];
      snprintf( prebuf, sizeof(prebuf), "%.3g", pre_dev );
      snprintf( postbuf, sizeof(postbuf), "%.3g", post_dev );

      msg = WString::tr("ect-fit-successful-detail")
              .arg( static_cast<int>(peakInfos.size()) )
              .arg( prebuf )
              .arg( postbuf );
    }catch( std::exception & )
    {
      //energy_for_channel can throw for out-of-range channels - just show the plain message
    }

    m_interspec->logMessage( msg, 1 );
  }catch( std::exception &e )
  {
    WString msg = WString::tr("ect-fail-fit").arg( e.what() );
    cerr << "EnergyCalTool::fitCoefficients():\n\tCaught: " << msg.toUTF8() << endl;
    m_interspec->logMessage( msg, 3 );
  }//try / catch
}//void fitCoefficients()


void EnergyCalTool::updateFitButtonStatus()
{
  const CanFitCoefStatus status = canDoEnergyFit();
  const bool canFit = (status == CanFitCoefStatus::CanFit);

#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
  WString hint;
  switch( status )
  {
    case CanFitCoefStatus::CanFit:
      break;

    case CanFitCoefStatus::ForegroundNotApplied:
      hint = WString::tr("ect-fit-hint-fore-deselected");
      break;

    case CanFitCoefStatus::NoPeaksToUse:
      hint = WString::tr("ect-fit-hint-no-peaks");
      break;

    case CanFitCoefStatus::NotForegroundCal:
      hint = WString::tr("ect-fit-hint-not-fore-cal");
      break;

    case CanFitCoefStatus::NoCoefSelected:
      hint = WString::tr("ect-fit-hint-no-coefs");
      break;

    case CanFitCoefStatus::MorePeaksNeeded:
      hint = WString::tr("ect-fit-hint-few-peaks");
      break;

    case CanFitCoefStatus::InvalidCalType:  //the convert-to-polynomial message covers this
    case CanFitCoefStatus::InternalError:   //not user-actionable
      break;
  }//switch( status )

  // The "Fit" checkboxes are only consulted on the foreground displayed-detector display (the
  //  "fit source" - see #fitCoefficients), so disable them everywhere else to avoid implying they
  //  have an independent effect (their states are kept in sync across displays).
  const vector<string> foredets
                       = m_interspec->detectorsToDisplay( SpecUtils::SpectrumType::Foreground );

  for( EnergyCalImp::CalDisplay *disp : calDisplays() )
  {
    const bool fitSource = (disp->spectrumType() == SpecUtils::SpectrumType::Foreground)
        && (std::find(begin(foredets), end(foredets), disp->detectorName()) != end(foredets));
    disp->setFitCheckboxesDisabled( !fitSource );
    disp->setFitStatus( canFit, hint );
  }
#else
  if( canFit != m_fitCalBtn->isEnabled() )
    m_fitCalBtn->setDisabled( !canFit );
#endif
}//void updateFitButtonStatus()


void EnergyCalTool::fitCheckboxChanged( EnergyCalImp::CalDisplay *display )
{
#if( IMP_COEF_FIT_BTN_NEAR_COEFS )
  // The selection of which coefficient orders to fit is a property of the fit itself - not of any
  //  one detector - since the fit result gets propagated to all applicable detectors; so keep the
  //  "Fit" checkboxes in sync (by coefficient order) across all displays.
  if( display )
  {
    const set<size_t> fitfor = display->fitForCoefficents();
    for( EnergyCalImp::CalDisplay *disp : calDisplays() )
    {
      if( disp != display )
        disp->setFitFor( fitfor );
    }
  }//if( display )
#endif

  updateFitButtonStatus();
}//void fitCheckboxChanged( EnergyCalImp::CalDisplay *display )


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
      m_graphicalRecal = AuxWindow::make<EnergyCalGraphicalConfirm>( xstart, xfinish, this,
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
    AuxWindow::deleteAuxWindow( m_graphicalRecal.get() );
    assert( !m_graphicalRecal );   //observing_ptr auto-clears
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
  
  const ApplyToSelection selection = currentApplyToSelection();
  const std::array<bool,NumApplyToCbIndex> applicable
                                        = applicable_apply_to_rows( currentApplyToDisplayState() );
  const bool selectedDetOnly = applicable[ApplyToCbIndex::ApplyToSelectedDetOnly]
                               && selection.checked[ApplyToCbIndex::ApplyToSelectedDetOnly];
  
  for( ApplyToCbIndex index = static_cast<ApplyToCbIndex>(0);
      index < ApplyToCbIndex::NumApplyToCbIndex;
      index = ApplyToCbIndex(index+1) )
  {
    if( !applicable[index] || !selection.checked[index] )
      continue;
    
    if( selectedDetOnly && ((index == ApplyToCbIndex::ApplyToDisplayedDetectors)
                            || (index == ApplyToCbIndex::ApplyToAllDetectors)) )
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
      case ApplyToCbIndex::ApplyToSelectedDetOnly:
        answer += WString::tr("ect-lc-selected-det-only")
                    .arg( Wt::Utils::htmlEncode(selection.selected_detector) ).toUTF8();
        break;
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
  auto specTypeMenuStackOwned = std::make_unique<WStackedWidget>();
  m_specTypeMenuStack = specTypeMenuStackOwned.get();
  m_specTypeMenuStack->addStyleClass( "CalSpecStack" );

  auto specTypeMenuOwned = std::make_unique<WMenu>(m_specTypeMenuStack);
  m_specTypeMenu = specTypeMenuOwned.get();
  m_specTypeMenu->addStyleClass( "CalSpecMenu" );
  m_specTypeMenu->itemSelected().connect( this, &EnergyCalTool::specTypeToDisplayForChanged );
  m_detColLayout->addWidget( std::move(specTypeMenuOwned), 1, 0 );
  m_detColLayout->addWidget( std::move(specTypeMenuStackOwned), 2, 0 );

  WGridLayout * const callayout = dynamic_cast<WGridLayout *>( m_calColumn->layout() );
  assert( callayout );

  auto calInfoDisplayStackOwned = std::make_unique<WStackedWidget>();
  m_calInfoDisplayStack = calInfoDisplayStackOwned.get();
  m_calInfoDisplayStack->addStyleClass( "ToolTabTitledColumnContent CalStack" );
  callayout->addWidget( std::move(calInfoDisplayStackOwned), 1, 1 );

  for( int i = 0; i < 3; ++i )
  {
    auto detMenuDivOwned = std::make_unique<WContainerWidget>();  //this holds the WMenu for this SpecFile
    WContainerWidget *detMenuDiv = detMenuDivOwned.get();
    detMenuDiv->addStyleClass( "DetMenuDiv" );

    WMenuItem *item = m_specTypeMenu->addItem( WString::tr(spec_type_labels[i]),
                                               std::move(detMenuDivOwned), ContentLoading::Eager );
    //Fix issue, for Wt 3.3.4 at least, if user doesnt click exactly on the <a> element
    item->clicked().connect( item, [item](){ item->select(); } );
    item->setHidden( true );  //#doRefreshFromFiles will un-hide, as needed

    m_detectorMenu[i] = detMenuDiv->addNew<WMenu>( m_calInfoDisplayStack );
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
    }else if( m_haveRelayoutSelection )
    {
      // The displays were just re-created for a layout change; re-select what was selected before
      haveCurrent = true;
      current_type = m_relayoutSelectionType;
      current_det = m_relayoutSelectionDet;
    }
    
    m_haveRelayoutSelection = false;
  }//end code-block to get the currently showing CalDisplay

  // We want to preserve which "Fit" check boxes are set for any CalDisplay we delete, but will
  //  then re-create in another menu (i.e., when `specTypeInForgrndMenu` flipped because the
  //  displayed files/detectors changed).
  map< pair<SpecUtils::SpectrumType,string>, set<size_t> > set_fit_for_cbs;

  // First pass: remove the CalDisplays that are no longer wanted in their menu.
  //  Note: Wt4's WMenu::removeItem() pulls the (Eager) contents back out of the shared contents
  //  stack and into the item, then hands us the owning unique_ptr for the item - so dropping that
  //  unique_ptr destroys both the WMenuItem and its CalDisplay.  Deleting either explicitly is a
  //  double-free.
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

      const std::unique_ptr<WMenuItem> doomed = detMenu->removeItem( item );
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
        auto displayOwned = std::make_unique<EnergyCalImp::CalDisplay>( this, info.type, info.detname, isWide );
        auto display = displayOwned.get();
        item = detMenu->insertItem( static_cast<int>(j), info.label, std::move(displayOwned),
                                    ContentLoading::Eager );
        //Fix issue, for Wt 3.3.4 at least, if user doesnt click exactly on the <a> element
        item->clicked().connect( item, [item](){ item->select(); } );

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
      disp->setFitStatus( false, WString() );  //silent - the whole cal column is hidden anyway
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
  
  // Show only the "Apply Changes To" rows that are applicable; the hidden rows are also set to
  //  their no-op values, so when a row re-appears it shows what was actually being applied.
  //  (#apply_to_scope interprets a non-applicable row as a no-op regardless of its state.)
  const std::array<bool,NumApplyToCbIndex> applicable
                                        = applicable_apply_to_rows( currentApplyToDisplayState() );
  
  bool anyApplyToCbShown = false;
  for( ApplyToCbIndex index = static_cast<ApplyToCbIndex>(0);
      index < ApplyToCbIndex::NumApplyToCbIndex;
      index = ApplyToCbIndex(index+1) )
  {
    Wt::WCheckBox *cb = m_applyToCbs[index];
    const auto cbparent = cb->parent();
    assert( cbparent );
    assert( dynamic_cast<WContainerWidget *>(cbparent) );
    
    const bool hideRow = !applicable[index];
    
    switch( index )
    {
      case ApplyToCbIndex::ApplyToForeground:
      case ApplyToCbIndex::ApplyToAllDetectors:
      case ApplyToCbIndex::ApplyToAllSamples:
        if( hideRow )
          cb->setChecked( true );
        break;
        
      case ApplyToCbIndex::ApplyToDisplayedDetectors:
      case ApplyToCbIndex::ApplyToSelectedDetOnly:
      case ApplyToCbIndex::ApplyToDisplayedSamples:
        if( hideRow )
          cb->setChecked( false );
        break;
        
      case ApplyToCbIndex::ApplyToBackground:
      case ApplyToCbIndex::ApplyToSecondary:
        // We'll ignore case where background/secondary is same {SpecFile, Samples, Detectors} as
        //  the foreground; we also leave their state alone so a newly loaded background keeps the
        //  users choice.
        break;
        
      case ApplyToCbIndex::NumApplyToCbIndex:
        assert( 0 );
        break;
    }//switch( index )
    
    // Keep the Displayed/All checkbox pairs mutually exclusive
    static_assert( ApplyToCbIndex::ApplyToDisplayedDetectors < ApplyToCbIndex::ApplyToAllDetectors, "" );
    static_assert( ApplyToCbIndex::ApplyToDisplayedSamples < ApplyToCbIndex::ApplyToAllSamples, "" );
    if( (index == ApplyToCbIndex::ApplyToAllDetectors)
        && m_applyToCbs[ApplyToCbIndex::ApplyToDisplayedDetectors]->isChecked() )
      cb->setChecked( false );
    if( (index == ApplyToCbIndex::ApplyToAllSamples)
        && m_applyToCbs[ApplyToCbIndex::ApplyToDisplayedSamples]->isChecked() )
      cb->setChecked( false );
    
    if( cbparent->isHidden() != hideRow )
      cbparent->setHidden( hideRow );
    
    if( !hideRow )
      anyApplyToCbShown = true;
  }//for( loop over ApplyToCbIndex )
  
  // With "Indep. Dets." checked, the Displayed/All Detectors choice doesnt apply
  const bool selectedDetOnly = m_applyToCbs[ApplyToCbIndex::ApplyToSelectedDetOnly]->isChecked();
  for( const ApplyToCbIndex index : {ApplyToDisplayedDetectors, ApplyToAllDetectors} )
  {
    if( m_applyToCbs[index]->isDisabled() != selectedDetOnly )
      m_applyToCbs[index]->setDisabled( selectedDetOnly );
  }
  
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

      case MoreActionsIndex::GainMatch:
      {
        // Show when the foreground has multiple gamma detectors, or when at least two
        //  distinct files are being displayed (i.e., there is something to match).
        bool show = specfiles[0] && (gammaDetectorsForDisplayedSamples(SpecUtils::SpectrumType::Foreground).size() > 1);
        if( !show && specfiles[0] )
          show = (specfiles[1] && (specfiles[1] != specfiles[0]))
                 || (specfiles[2] && (specfiles[2] != specfiles[0]) && (specfiles[2] != specfiles[1]));
        aparent->setHidden( !show );
        break;
      }//case MoreActionsIndex::GainMatch:

      case MoreActionsIndex::GuessGain:
        // Useful whenever a foreground is loaded - not only for cal-less files, since a file
        //  may also carry a *wrong* calibration.
        aparent->setHidden( !specfiles[0] );
        break;

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
    AuxWindow::deleteAuxWindow( m_addActionWindow.get() );
    assert( !m_addActionWindow );   //observing_ptr auto-clears
  }
  
  m_addActionWindow = AuxWindow::make<EnergyCalAddActionsWindow>( index, measToChange, this );
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
    AuxWindow::deleteAuxWindow( m_addActionWindow.get() );
    assert( !m_addActionWindow );   //observing_ptr auto-clears
  }
}//void cancelMoreActionWindow()


void EnergyCalTool::closeStaleMoreActionWindow()
{
  if( m_addActionWindow && m_addActionWindow->usesApplyToSnapshot() )
    cancelMoreActionWindow();
}//void closeStaleMoreActionWindow()


void EnergyCalTool::render( Wt::WFlags<Wt::RenderFlag> flags)
{
  //flags.test(RenderFlag::RenderFull) will only be true on initial rending of widget, and
  //  after that only the RenderFlag::RenderUpdate flag will be set
  
  if( flags.test(Wt::RenderFlag::Full)
      || m_renderFlags.test(EnergyCalToolRenderFlags::FullGuiUpdate) )
  {
    doRefreshFromFiles();
    m_renderFlags.clear( EnergyCalToolRenderFlags::FullGuiUpdate );
  }
  
  WContainerWidget::render(flags);
}//void render( Wt::WFlags<Wt::RenderFlag> flags)


void EnergyCalTool::applyToCbChanged( const EnergyCalTool::ApplyToCbIndex index )
{
  // We only get here if the user checked/unchecked a checkbox, or in a undo/redo step.
  assert( index >= 0 && index < EnergyCalTool::NumApplyToCbIndex );
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
      
    case EnergyCalTool::ApplyToSelectedDetOnly:
      // The Displayed/All Detectors choice doesnt apply while only changing the selected detector
      m_applyToCbs[EnergyCalTool::ApplyToDisplayedDetectors]->setDisabled( isChecked );
      m_applyToCbs[EnergyCalTool::ApplyToAllDetectors]->setDisabled( isChecked );
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
  
  closeStaleMoreActionWindow();
  
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
