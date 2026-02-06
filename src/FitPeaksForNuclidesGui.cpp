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

#include <set>
#include <string>
#include <vector>
#include <memory>
#include <atomic>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WServer>
#include <Wt/WCheckBox>
#include <Wt/WComboBox>
#include <Wt/WGridLayout>
#include <Wt/WIOService>
#include <Wt/WPushButton>
#include <Wt/WApplication>

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/AnalystChecks.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/FitPeaksForNuclides.h"
#include "InterSpec/ReferenceLinePredef.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/PhysicalUnitsLocalized.h"
#include "InterSpec/FitPeaksForNuclidesGui.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"

using namespace std;
using namespace Wt;


namespace FitPeaksForNuclidesGui
{

SimpleDialog *showAdvancedDialog()
{
  InterSpec *viewer = InterSpec::instance();
  if( !viewer )
    return nullptr;

  viewer->useMessageResourceBundle( "FitPeaksForNuclidesGui" );

  FitPeaksAdvancedDialog *dlg = new FitPeaksAdvancedDialog( WString() );
  dlg->rejectWhenEscapePressed();

  return dlg;
}//SimpleDialog *showAdvancedDialog()


void startFitSources( const bool /*from_advanced_dialog*/ )
{
  InterSpec *viewer = InterSpec::instance();
  if( !viewer )
    return;

  ReferencePhotopeakDisplay *ref_disp = viewer->referenceLinesWidget();
  if( !ref_disp )
    return;

  viewer->useMessageResourceBundle( "FitPeaksForNuclidesGui" );

  std::shared_ptr<const SpecUtils::Measurement> foreground
                                 = viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  if( !foreground )
  {
    passMessage( WString::tr("fpn-no-data"), WarningWidget::WarningMsgInfo );
    return;
  }
  
  WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> options;

  // Collect sources (nuclide/element/reaction) from currently showing + persisted.
  std::vector<RelActCalcAuto::SrcVariant> sources;
  std::vector<RelActCalcAuto::NucInputInfo> base_nucs;
  sources.reserve( 8 );
  base_nucs.reserve( 8 );

  auto add_src = [&]( const RelActCalcAuto::SrcVariant &src,
                      const std::string &age_str ) {
    sources.push_back( src );

    RelActCalcAuto::NucInputInfo info;
    info.source = src;
    info.fit_age = false;
    info.age = 0.0;

    const SandiaDecay::Nuclide *nuc = RelActCalcAuto::nuclide( src );
    if( nuc )
    {
      // If age_str is blank, FitPeaksForNuclides will use PeakDef::defaultDecayTime().
      // Otherwise parse user value (including "HL").
      if( age_str.empty() )
      {
        info.age = -1.0;
      }else
      {
        try
        {
          info.age = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( age_str, nuc->halfLife );
        }catch( std::exception & )
        {
          info.age = -1.0;
        }
      }
    }

    base_nucs.push_back( info );
  };

  auto add_from_info = [&]( const ReferenceLineInfo &info ) {
    if( info.m_validity != ReferenceLineInfo::InputValidity::Valid )
      return;

    switch( info.m_source_type )
    {
      case ReferenceLineInfo::SourceType::Nuclide:
        if( info.m_nuclide )
          add_src( info.m_nuclide, info.m_input.m_age );
        break;

      case ReferenceLineInfo::SourceType::FluorescenceXray:
        if( info.m_element )
          add_src( info.m_element, "" );
        break;

      case ReferenceLineInfo::SourceType::Reaction:
        for( const ReactionGamma::Reaction *rctn : info.m_reactions )
        {
          if( rctn )
            add_src( rctn, "" );
        }
        break;
        
      case ReferenceLineInfo::SourceType::NuclideMixture:
        assert( info.m_nuc_mix );
        if( info.m_nuc_mix )
        {
          for( const ReferenceLinePredef::NucMixComp &comp : info.m_nuc_mix->m_components )
          {
            assert( comp.m_nuclide );
            
            sources.push_back( comp.m_nuclide );
            
            RelActCalcAuto::NucInputInfo input_info;
            input_info.source = comp.m_nuclide;
            input_info.fit_age = false;
            input_info.age = std::max( 0.0, info.m_nuc_mix->m_default_age - comp.m_age_offset );
            
            base_nucs.push_back( input_info );
          }//for( const ReferenceLinePredef::NucMixComp &comp : info.m_nuc_mix->m_components )
        }//if( info.m_nuc_mix )
        break;
        
      case ReferenceLineInfo::SourceType::Background:
        options |= FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaks;
        break;

      default:
        break;
    }//switch( m_source_type )
  };

  const ReferenceLineInfo &current_showing = ref_disp->currentlyShowingNuclide();
  const std::vector<ReferenceLineInfo> &persisted = ref_disp->persistedNuclides();

  add_from_info( current_showing );
  for( const ReferenceLineInfo &info : persisted )
    add_from_info( info );

  if( sources.empty() && !options.testFlag(FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaks) )
  {
    assert( 0 );
    passMessage( "Sorry, cant fit for any of the current sources.", WarningWidget::WarningMsgMedium ); //Skip localization since we dont expect to be here
    return;
  }

  const bool isHPGe = PeakFitUtils::is_likely_high_res( viewer );
  std::shared_ptr<const SpecUtils::Measurement> background
                                 = viewer->displayedHistogram( SpecUtils::SpectrumType::Background );

  std::shared_ptr<SpecMeas> fg_meas = viewer->measurment( SpecUtils::SpectrumType::Foreground );
  std::shared_ptr<const DetectorPeakResponse> drf_input = fg_meas ? fg_meas->detector() : nullptr;

  FitPeaksForNuclides::PeakFitForNuclideConfig config;

  // Let user continue using app while fitting.
  SimpleDialog *wait_dlg = new SimpleDialog( WString::tr("fpn-wait-title"),
                                             WString::tr("fpn-wait-content") );
  wait_dlg->rejectWhenEscapePressed();
  wait_dlg->addButton( WString::tr("Close") );

  // Disable button to prevent repeats while running.
  ref_disp->setFitSourcesButtonEnabled( false );

  boost::function<void(void)> close_wait
                = wApp->bind( boost::bind( &SimpleDialog::accept, wait_dlg ) );

  Wt::WServer *server = Wt::WServer::instance();
  if( !server )
  {
    ref_disp->setFitSourcesButtonEnabled( true );
    return;
  }

  // Create copies of measurements on GUI thread for background work.
  std::shared_ptr<SpecUtils::Measurement> fg_copy = std::make_shared<SpecUtils::Measurement>( *foreground );
  std::shared_ptr<SpecUtils::Measurement> bg_copy = background
                                 ? std::make_shared<SpecUtils::Measurement>( *background ) : nullptr;

  std::weak_ptr<const SpecUtils::Measurement> weak_foreground = foreground;
  const std::string sessionid = wApp->sessionId();

  // Store ReferenceLineInfo objects for later display and color mapping
  std::vector<ReferenceLineInfo> ref_lines_for_display;
  if( current_showing.m_validity == ReferenceLineInfo::InputValidity::Valid )
    ref_lines_for_display.push_back( current_showing );
  for( const ReferenceLineInfo &info : persisted )
  {
    if( info.m_validity == ReferenceLineInfo::InputValidity::Valid )
      ref_lines_for_display.push_back( info );
  }

  std::shared_ptr<FitPeaksForNuclides::PeakFitResult> result = std::make_shared<FitPeaksForNuclides::PeakFitResult>();

  // Stage A: GUI thread - detect peaks using AnalystChecks
  WServer::instance()->post( sessionid, [=](){
    InterSpec *viewer_a = InterSpec::instance();
    if( !wApp || !viewer_a )
    {
      close_wait();
      ReferencePhotopeakDisplay *disp = viewer_a ? viewer_a->referenceLinesWidget() : nullptr;
      if( disp )
        disp->setFitSourcesButtonEnabled( true );
      return;
    }

    std::vector<std::shared_ptr<const PeakDef>> auto_search_peaks;
    std::vector<std::shared_ptr<const PeakDef>> user_peaks;
    try
    {
      AnalystChecks::DetectedPeaksOptions peak_opts;
      peak_opts.specType = SpecUtils::SpectrumType::Foreground;
      peak_opts.nonBackgroundPeaksOnly = false;
      const AnalystChecks::DetectedPeakStatus detected = AnalystChecks::detected_peaks( peak_opts, viewer_a );
      auto_search_peaks = detected.peaks;

      // Capture current user peaks (by shared_ptr) for DoNotUseExistingRois / ExistingPeaksAsFreePeak.
      // These are the exact pointers in PeakModel, preserving identity for later removal.
      PeakModel *pm = viewer_a->peakModel();
      if( pm )
      {
        const std::shared_ptr<const std::deque<PeakModel::PeakShrdPtr>> current_peaks = pm->peaks();
        if( current_peaks )
          user_peaks.assign( current_peaks->begin(), current_peaks->end() );
      }
    }catch( std::exception &e )
    {
      // Failed to detect peaks - show error and return
      close_wait();
      ReferencePhotopeakDisplay *disp = viewer_a ? viewer_a->referenceLinesWidget() : nullptr;
      if( disp )
        disp->setFitSourcesButtonEnabled( true );

      SimpleDialog *err = new SimpleDialog( WString::tr("fpn-error-title"),
                                            WString::tr("fpn-error-content").arg( e.what() ) );
      err->addButton( WString::tr("Okay") );
      wApp->triggerUpdate();
      return;
    }

    // Stage B: Background thread - perform peak fitting
    server->ioService().boost::asio::io_service::post( std::bind( [=](){
      // Note: PeakFitResult does not currently initialize `status`, so set it here.
      result->status = RelActCalcAuto::RelActAutoSolution::Status::Success;
      result->error_message.clear();

      try
      {
        *result = FitPeaksForNuclides::fit_peaks_for_nuclides( auto_search_peaks, fg_copy, base_nucs, user_peaks, bg_copy, drf_input, options, config, isHPGe );
      }catch( std::exception &e )
      {
        result->status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
        result->error_message = std::string("Fit failed: ") + e.what();
      }catch( ... )
      {
        result->status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
        result->error_message = "Fit failed with unknown error";
      }

      // Stage C: GUI thread - handle results
      WServer::instance()->post( sessionid, [=](){
        close_wait();

        InterSpec *viewer_c = InterSpec::instance();
        ReferencePhotopeakDisplay *disp = viewer_c ? viewer_c->referenceLinesWidget() : nullptr;
        if( disp )
          disp->setFitSourcesButtonEnabled( true );

        if( !wApp || !viewer_c )
          return;

        // If a new spectrum has been loaded, ignore these results.
        std::shared_ptr<const SpecUtils::Measurement> current_fg = weak_foreground.lock();
        if( !current_fg || viewer_c->displayedHistogram(SpecUtils::SpectrumType::Foreground) != current_fg )
        {
          passMessage( WString::tr("fpn-stale"), WarningWidget::WarningMsgInfo );
          wApp->triggerUpdate();
          return;
        }

        if( result->status != RelActCalcAuto::RelActAutoSolution::Status::Success )
        {
          SimpleDialog *err = new SimpleDialog( WString::tr("fpn-error-title"),
                                                WString::tr("fpn-error-content").arg( result->error_message ) );
          err->addButton( WString::tr("Okay") );
          wApp->triggerUpdate();
          return;
        }

        const bool auto_accept = UserPreferences::preferenceValue<bool>( "AutoAcceptFitSourcesPeaks", viewer_c );

        // Show warnings via passMessage if auto-accepting, otherwise will show in dialog
        if( auto_accept )
        {
          for( const std::string &w : result->warnings )
            passMessage( w, WarningWidget::WarningMsgLow );
        }

        if( result->observable_peaks.empty() )
        {
          Wt::WString msg = WString::tr("fpn-no-peaks");
          if( !result->warnings.empty() && !auto_accept )
          {
            msg += Wt::WString::fromUTF8( "\n\n" );
            msg += WString::tr("fpn-warnings-header");
            for( const std::string &w : result->warnings )
            {
              msg += Wt::WString::fromUTF8( "\n• " );
              msg += Wt::WString::fromUTF8( w );
            }
          }
          SimpleDialog *dlg = new SimpleDialog( WString::tr("fpn-result-title"), msg );
          dlg->addButton( WString::tr("Okay") );
          wApp->triggerUpdate();
          return;
        }

        if( auto_accept )
        {
          PeakModel *peak_model = viewer_c->peakModel();
          if( peak_model )
          {
            UndoRedoManager::PeakModelChange peak_undo_creator;
            if( !result->original_peaks_to_remove.empty() )
              peak_model->removePeaks( result->original_peaks_to_remove );
            peak_model->addPeaks( result->observable_peaks );
          }
          passMessage( WString::tr("fpn-auto-accepted"), WarningWidget::WarningMsgInfo );
          wApp->triggerUpdate();
          return;
        }

        // Show result dialog with warnings
        SimpleDialog *result_dlg = new SimpleDialog( WString::tr("fpn-result-title"), "" );
        result_dlg->addStyleClass( "FitSourcesResultDialog" );

        Wt::WContainerWidget *contents = result_dlg->contents();
        contents->addStyleClass( "FitSourcesResultContents" );

        // Add summary information about the fit
        const size_t num_peaks = result->observable_peaks.size();
        const size_t num_sources = ref_lines_for_display.size();
        Wt::WString summary = WString::tr("fpn-summary").arg( static_cast<int>(num_peaks) ).arg( static_cast<int>(num_sources) );
        Wt::WText *summary_text = new Wt::WText( summary, contents );
        summary_text->addStyleClass( "FitSourcesSummary" );
        summary_text->setInline( false );

        // Add warnings section if there are any
        if( !result->warnings.empty() )
        {
          Wt::WText *warnings_header = new Wt::WText( WString::tr("fpn-warnings-header"), contents );
          warnings_header->addStyleClass( "fpn-warnings-header" );

          for( const std::string &w : result->warnings )
          {
            Wt::WText *warning = new Wt::WText( Wt::WString::fromUTF8( "• " + w ), contents );
            warning->addStyleClass( "fpn-warning" );
          }

          Wt::WText *spacer = new Wt::WText( "", contents );
          spacer->setInline( false );
        }

        // Size chart dynamically based on rendered dimensions
        int chartw = 420, charth = 260;
        if( viewer_c->renderedWidth() > 500 )
          chartw = std::min( ((3*viewer_c->renderedWidth()/4) - 50), 600 );
        if( viewer_c->renderedHeight() > 400 )
          charth = std::min( viewer_c->renderedHeight()/3, (4*chartw)/7 );
        chartw = std::max( chartw, 300 );
        charth = std::max( charth, 200 );

        D3SpectrumDisplayDiv *spectrum = new D3SpectrumDisplayDiv( contents );
        spectrum->clicked().preventPropagation();
        spectrum->setThumbnailMode();
        spectrum->setMinimumSize( 300, 200 );
        spectrum->resize( chartw, charth );
        spectrum->setData( current_fg, false );

        // Set reference photopeak lines for the sources that were fit
        for( const ReferenceLineInfo &ref_info : ref_lines_for_display )
        {
          if( ref_info.m_validity == ReferenceLineInfo::InputValidity::Valid )
          {
            spectrum->setReferncePhotoPeakLines( ref_info );
            spectrum->persistCurrentReferncePhotoPeakLines();
          }
        }

        PeakModel *preview_model = new PeakModel( spectrum );
        preview_model->setNoSpecMeasBacking();
        preview_model->setForeground( current_fg );
        spectrum->setPeakModel( preview_model );

        // Get existing peaks and add new ones, setting colors based on source
        std::vector<PeakDef> preview_peaks;
        PeakModel *peak_model = viewer_c->peakModel();
        if( peak_model )
          preview_peaks = peak_model->peakVec();

        // Add new peaks and assign colors based on their source
        for( PeakDef &peak : result->observable_peaks )
        {
          // Determine source name from peak
          std::string src_name;
          if( peak.parentNuclide() )
            src_name = peak.parentNuclide()->symbol;
          else if( peak.xrayElement() )
            src_name = peak.xrayElement()->symbol;
          else if( peak.reaction() )
            src_name = peak.reaction()->name();

          // Find matching ReferenceLineInfo to get color
          if( !src_name.empty() && disp )
          {
            for( const ReferenceLineInfo &ref_info : ref_lines_for_display )
            {
              bool matches = false;
              Wt::WColor nuc_mix_comp_color;
              if( ref_info.m_nuclide && peak.parentNuclide() == ref_info.m_nuclide )
                matches = true;
              if( !matches && ref_info.m_element && peak.xrayElement() == ref_info.m_element )
                matches = true;
              if( !matches && peak.reaction() && ref_info.m_reactions.count( peak.reaction() ) > 0 )
                matches = true;
              if( !matches
                 && peak.parentNuclide()
                 && (ref_info.m_source_type == ReferenceLineInfo::SourceType::NuclideMixture)
                 && ref_info.m_nuc_mix )
              {
                for( const ReferenceLinePredef::NucMixComp &comp : ref_info.m_nuc_mix->m_components )
                {
                  if( comp.m_nuclide == peak.parentNuclide() )
                  {
                    matches = true;
                    nuc_mix_comp_color = comp.m_color;
                    break;
                  }
                }
              }

              if( matches )
              {
                // Get source name for color lookup
                std::string ref_src;
                if( ref_info.m_nuclide )
                  ref_src = ref_info.m_nuclide->symbol;
                else if( ref_info.m_element )
                  ref_src = ref_info.m_element->symbol;
                else if( !ref_info.m_reactions.empty() )
                {
                  const ReactionGamma::Reaction *rctn = *ref_info.m_reactions.begin();
                  if( rctn )
                    ref_src = rctn->name();
                }

                // Per-component color from NuclideMixture takes priority
                Wt::WColor color = nuc_mix_comp_color.isDefault() ? ref_info.m_input.m_color : nuc_mix_comp_color;
                if( color.isDefault() && !ref_src.empty() )
                  color = disp->suggestColorForSource( ref_src );
                if( color.isDefault() && !ref_src.empty() )
                {
                  color = disp->nextGenericSourceColor();
                  if( !color.isDefault() )
                    disp->updateColorCacheForSource( ref_src, color );
                }
                if( !color.isDefault() )
                  peak.setLineColor( color );
                break;
              }
            }
          }

          preview_peaks.push_back( peak );
        }

        preview_model->addPeaks( preview_peaks );

        // Checkbox goes below the chart (lowest item)
        Wt::WCheckBox *cb = new Wt::WCheckBox( WString::tr("fpn-auto-accept"), contents );
        UserPreferences::associateWidget( "AutoAcceptFitSourcesPeaks", cb, viewer_c );
        cb->setInline( false );

        Wt::WPushButton *accept_btn = result_dlg->addButton( WString::tr("Accept") );
        result_dlg->addButton( WString::tr("Cancel") );
        accept_btn->setFocus();

        accept_btn->clicked().connect( std::bind( [viewer_c, result](){
          if( viewer_c )
          {
            PeakModel *peak_model = viewer_c->peakModel();
            if( peak_model )
            {
              UndoRedoManager::PeakModelChange peak_undo_creator;
              if( !result->original_peaks_to_remove.empty() )
                peak_model->removePeaks( result->original_peaks_to_remove );
              peak_model->addPeaks( result->observable_peaks );
            }
          }
        } ) );
        // Cancel button uses SimpleDialog's default close behavior.

        wApp->triggerUpdate();
      } );
    } ) );
  } );
}//void startFitSources(...)

}//namespace FitPeaksForNuclidesGui


// FitPeaksAdvancedDialog (global class, not in namespace)
FitPeaksAdvancedDialog::FitPeaksAdvancedDialog( const Wt::WString &title )
  : SimpleDialog( title, "" ), m_widget( nullptr ), m_acceptBtn( nullptr )
{
  addStyleClass( "FitPeaksAdvancedDialog" );
  wApp->useStyleSheet( "InterSpec_resources/FitPeaksForNuclidesGui.css" );

  Wt::WGridLayout *layout = new Wt::WGridLayout();
  layout->setVerticalSpacing( 0 );
  layout->setHorizontalSpacing( 0 );
  layout->setContentsMargins( 0, 0, 0, 0 );
  contents()->setLayout( layout );

  m_widget = new FitPeaksAdvancedWidget();
  layout->addWidget( m_widget, 0, 0 );

  addButton( Wt::WString::tr("fpn-cancel") );
  m_acceptBtn = addButton( Wt::WString::tr("fpn-accept") );
  m_acceptBtn->setEnabled( m_widget->canAccept() );

  m_acceptBtn->clicked().connect( this, &FitPeaksAdvancedDialog::onAcceptClicked );
  m_widget->resultUpdated().connect( this, &FitPeaksAdvancedDialog::onResultUpdated );

  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  if( interspec )
  {
    if( interspec->isPhone() )
    {
      addStyleClass( "FitPeaksAdvancedDialog-phone" );
#if ( IOS )
      addStyleClass( "FitPeaksAdvancedDialog-iphone" );
#endif
    }
    const bool portrait = ( ( interspec->renderedWidth() > 100 ) && ( interspec->renderedHeight() > 100 ) &&
                           ( interspec->renderedWidth() < 480 ) );
    if( portrait )
    {
      addStyleClass( "FitPeaksAdvancedDialog-portrait" );
      m_widget->setWidth( 0.95 * interspec->renderedWidth() - 30 );
      m_widget->setHeight( 0.95 * interspec->renderedHeight() - 90 );
      m_widget->setMinimumSize( 0.95 * interspec->renderedWidth() - 30, 0.95 * interspec->renderedHeight() - 90 );
      m_widget->setMaximumSize( 0.95 * interspec->renderedWidth() - 30, 0.95 * interspec->renderedHeight() - 90 );
    }
    else if( ( interspec->renderedWidth() > 100 ) && ( interspec->renderedHeight() > 50 ) )
    {
      const double dialogWidth = std::min( 750.0, 0.95 * interspec->renderedWidth() );
      const double dialogHeight = std::min( 650.0, 0.95 * interspec->renderedHeight() );
      m_widget->setWidth( dialogWidth - 30 );
      m_widget->setHeight( dialogHeight - 90 );
      m_widget->setMinimumSize( dialogWidth - 30, dialogHeight - 90 );
      m_widget->setMaximumSize( dialogWidth - 30, dialogHeight - 90 );
    }
    else
    {
      m_widget->setWidth( 750 - 30 );
      m_widget->setHeight( 500 - 90 );
      m_widget->setMinimumSize( 750 - 30, 500 - 90 );
      m_widget->setMaximumSize( 750 - 30, 500 - 90 );
    }
  }

  rejectWhenEscapePressed();

  m_acceptBtn->setFocus(true);
}

FitPeaksAdvancedDialog::~FitPeaksAdvancedDialog()
{
}

void FitPeaksAdvancedDialog::onAcceptClicked()
{
  if( m_widget )
    m_widget->acceptResult();
  accept();
}

void FitPeaksAdvancedDialog::onResultUpdated()
{
  m_acceptBtn->setEnabled( m_widget ? m_widget->canAccept() : false );
}


// FitPeaksAdvancedWidget
FitPeaksAdvancedWidget::FitPeaksAdvancedWidget( Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_title( nullptr ),
    m_chart( nullptr ),
    m_chart_peak_model( nullptr ),
    m_status( nullptr ),
    m_warnings_div( nullptr ),
    m_options_div( nullptr ),
    m_opt_dont_use_rois( nullptr ),
    m_opt_existing_as_free( nullptr ),
    m_opt_dont_vary_energy_cal( nullptr ),
    m_opt_dont_refine_energy_cal( nullptr ),
    m_opt_fit_bkgnd_peaks( nullptr ),
    m_opt_fit_bkgnd_dont_use( nullptr ),
    m_opt_roi_min_chi2( nullptr ),
    m_opt_roi_min_peak_sig( nullptr ),
    m_opt_obs_initial_sig( nullptr ),
    m_opt_obs_final_sig( nullptr ),
    m_opt_fwhm_form( nullptr ),
    m_opt_rel_eff_type( nullptr ),
    m_opt_rel_eff_order( nullptr ),
    m_rel_eff_order_row( nullptr ),
    m_render_flags( 0 ),
    m_is_calculating( false ),
    m_calc_number( 0 ),
    m_auto_search_peaks( std::make_shared<std::vector<std::shared_ptr<const PeakDef>>>() )
{
  addStyleClass( "FitPeaksAdvancedWidget" );

  InterSpec *viewer = InterSpec::instance();
  if( !viewer )
    return;

  ReferencePhotopeakDisplay *ref_disp = viewer->referenceLinesWidget();
  if( !ref_disp )
    return;

  std::shared_ptr<const SpecUtils::Measurement> foreground
    = viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  if( !foreground )
    return;

  WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> options;

  auto add_src = [&]( const RelActCalcAuto::SrcVariant &src, const std::string &age_str ) {
    m_sources.push_back( src );
    RelActCalcAuto::NucInputInfo info;
    info.source = src;
    info.fit_age = false;
    info.age = 0.0;
    const SandiaDecay::Nuclide *nuc = RelActCalcAuto::nuclide( src );
    if( nuc )
    {
      if( age_str.empty() )
        info.age = -1.0;
      else
      {
        try
        {
          info.age = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( age_str, nuc->halfLife );
        }
        catch( std::exception & )
        {
          info.age = -1.0;
        }
      }
    }
    m_base_nucs.push_back( info );
  };

  auto add_from_info = [&]( const ReferenceLineInfo &info ) {
    if( info.m_validity != ReferenceLineInfo::InputValidity::Valid )
      return;
    switch( info.m_source_type )
    {
      case ReferenceLineInfo::SourceType::Nuclide:
        if( info.m_nuclide )
          add_src( info.m_nuclide, info.m_input.m_age );
        break;
      case ReferenceLineInfo::SourceType::FluorescenceXray:
        if( info.m_element )
          add_src( info.m_element, "" );
        break;
      case ReferenceLineInfo::SourceType::Reaction:
        for( const ReactionGamma::Reaction *rctn : info.m_reactions )
        {
          if( rctn )
            add_src( rctn, "" );
        }
        break;
      case ReferenceLineInfo::SourceType::NuclideMixture:
        if( info.m_nuc_mix )
        {
          for( const ReferenceLinePredef::NucMixComp &comp : info.m_nuc_mix->m_components )
          {
            if( comp.m_nuclide )
            {
              m_sources.push_back( comp.m_nuclide );
              RelActCalcAuto::NucInputInfo input_info;
              input_info.source = comp.m_nuclide;
              input_info.fit_age = false;
              input_info.age = std::max( 0.0, info.m_nuc_mix->m_default_age - comp.m_age_offset );
              m_base_nucs.push_back( input_info );
            }
          }
        }
        break;
      case ReferenceLineInfo::SourceType::Background:
        options |= FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaks;
        break;
      default:
        break;
    }
  };

  const ReferenceLineInfo &current_showing = ref_disp->currentlyShowingNuclide();
  const std::vector<ReferenceLineInfo> &persisted = ref_disp->persistedNuclides();
  add_from_info( current_showing );
  for( const ReferenceLineInfo &info : persisted )
    add_from_info( info );

  if( m_sources.empty() && !options.testFlag( FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaks ) )
    return;

  m_base_options = options;
  m_is_hpge = PeakFitUtils::is_likely_high_res( viewer );
  m_fg_copy = std::make_shared<SpecUtils::Measurement>( *foreground );
  std::shared_ptr<const SpecUtils::Measurement> background
    = viewer->displayedHistogram( SpecUtils::SpectrumType::Background );
  m_bg_copy = background ? std::make_shared<SpecUtils::Measurement>( *background ) : nullptr;
  std::shared_ptr<SpecMeas> fg_meas = viewer->measurment( SpecUtils::SpectrumType::Foreground );
  m_drf = fg_meas ? fg_meas->detector() : nullptr;

  if( current_showing.m_validity == ReferenceLineInfo::InputValidity::Valid )
    m_ref_lines.push_back( current_showing );
  for( const ReferenceLineInfo &info : persisted )
  {
    if( info.m_validity == ReferenceLineInfo::InputValidity::Valid )
      m_ref_lines.push_back( info );
  }

  PeakModel *pm = viewer->peakModel();
  if( pm )
  {
    std::shared_ptr<const std::deque<PeakModel::PeakShrdPtr>> current_peaks = pm->peaks();
    if( current_peaks )
      m_user_peaks.assign( current_peaks->begin(), current_peaks->end() );
  }

  m_title = new WText( WString::fromUTF8( "Fit peaks for: " + sourceListTitle() ), this );
  m_title->setInline( false );

  int chartw = 420;
  int charth = 260;
  if( viewer->renderedWidth() > 500 )
    chartw = std::min( ( ( 3 * viewer->renderedWidth() / 4 ) - 50 ), 600 );
  if( viewer->renderedHeight() > 400 )
    charth = std::min( viewer->renderedHeight() / 3, ( 4 * chartw ) / 7 );
  chartw = std::max( chartw, 300 );
  charth = std::max( charth, 200 );

  m_chart = new D3SpectrumDisplayDiv( this );
  m_chart->clicked().preventPropagation();
  m_chart->setThumbnailMode();
  m_chart->setMinimumSize( 300, 200 );
  m_chart->setData( m_fg_copy, false );
  for( const ReferenceLineInfo &ref_info : m_ref_lines )
  {
    m_chart->setReferncePhotoPeakLines( ref_info );
    m_chart->persistCurrentReferncePhotoPeakLines();
  }
  m_chart_peak_model = new PeakModel( m_chart );
  m_chart_peak_model->setNoSpecMeasBacking();
  m_chart_peak_model->setForeground( m_fg_copy );
  m_chart->setPeakModel( m_chart_peak_model );

  m_status = new WText( WString::tr("fpn-calculating"), this );
  m_status->addStyleClass( "fpn-status" );
  m_status->addStyleClass( "calculating" );

  m_warnings_div = new WContainerWidget( this );
  m_warnings_div->addStyleClass( "fpn-warnings" );
  m_warnings_div->addStyleClass( "fpn-warnings-computing" );

  m_options_div = new WContainerWidget( this );
  m_options_div->addStyleClass( "fpn-options" );
  buildOptionsFromConfig();

  m_render_flags |= UpdateCalculations;
  scheduleRender();
}

FitPeaksAdvancedWidget::~FitPeaksAdvancedWidget()
{
}

void FitPeaksAdvancedWidget::acceptResult()
{
  if( !canAccept() || !m_current_result )
    return;
  InterSpec *viewer = InterSpec::instance();
  if( !viewer )
    return;
  ReferencePhotopeakDisplay *disp = viewer->referenceLinesWidget();
  PeakModel *peak_model = viewer->peakModel();
  if( !peak_model )
    return;
  UndoRedoManager::PeakModelChange peak_undo_creator;
  if( !m_current_result->original_peaks_to_remove.empty() )
    peak_model->removePeaks( m_current_result->original_peaks_to_remove );
  std::vector<PeakDef> peaks_to_add = m_current_result->observable_peaks;

  std::vector<ReferenceLineInfo> ref_lines_to_use( m_ref_lines.begin(), m_ref_lines.end() );
  bool has_background_ref = false;
  for( const ReferenceLineInfo &r : m_ref_lines )
  {
    if( r.m_source_type == ReferenceLineInfo::SourceType::Background )
    {
      has_background_ref = true;
      break;
    }
  }
  if( m_base_options.testFlag( FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaks )
      && !has_background_ref )
  {
    RefLineInput bg_input;
    bg_input.m_input_txt = "background";
    bg_input.m_age = "";
    bg_input.m_showGammas = true;
    bg_input.m_showXrays = true;
    std::shared_ptr<ReferenceLineInfo> bg_ref = ReferenceLineInfo::generateRefLineInfo( bg_input );
    if( bg_ref && bg_ref->m_validity == ReferenceLineInfo::InputValidity::Valid )
      ref_lines_to_use.push_back( *bg_ref );
  }

  static const std::vector<std::string> norm_nuclide_symbols{ "U238", "Ra226", "U235", "Th232", "K40" };

  for( PeakDef &peak : peaks_to_add )
  {
    std::string src_name;
    if( peak.parentNuclide() )
      src_name = peak.parentNuclide()->symbol;
    else if( peak.xrayElement() )
      src_name = peak.xrayElement()->symbol;
    else if( peak.reaction() )
      src_name = peak.reaction()->name();
    if( src_name.empty() || !disp )
      continue;
    for( const ReferenceLineInfo &ref_info : ref_lines_to_use )
    {
      bool matches = false;
      Wt::WColor nuc_mix_comp_color;
      if( ref_info.m_nuclide && peak.parentNuclide() == ref_info.m_nuclide )
        matches = true;
      if( !matches && ref_info.m_element && peak.xrayElement() == ref_info.m_element )
        matches = true;
      if( !matches && peak.reaction() && ref_info.m_reactions.count( peak.reaction() ) > 0 )
        matches = true;
      if( !matches && peak.parentNuclide()
          && ref_info.m_source_type == ReferenceLineInfo::SourceType::NuclideMixture
          && ref_info.m_nuc_mix )
      {
        for( const ReferenceLinePredef::NucMixComp &comp : ref_info.m_nuc_mix->m_components )
        {
          if( comp.m_nuclide == peak.parentNuclide() )
          {
            matches = true;
            nuc_mix_comp_color = comp.m_color;
            break;
          }
        }
      }
      if( !matches && ref_info.m_source_type == ReferenceLineInfo::SourceType::Background
          && peak.parentNuclide() )
      {
        for( const std::string &sym : norm_nuclide_symbols )
        {
          if( peak.parentNuclide()->symbol == sym )
          {
            matches = true;
            break;
          }
        }
      }
      if( matches )
      {
        std::string ref_src;
        if( ref_info.m_nuclide )
          ref_src = ref_info.m_nuclide->symbol;
        else if( ref_info.m_element )
          ref_src = ref_info.m_element->symbol;
        else if( !ref_info.m_reactions.empty() )
        {
          const ReactionGamma::Reaction *rctn = *ref_info.m_reactions.begin();
          if( rctn )
            ref_src = rctn->name();
        }
        else if( ref_info.m_source_type == ReferenceLineInfo::SourceType::Background
                 && peak.parentNuclide() )
          ref_src = peak.parentNuclide()->symbol;
        Wt::WColor color = nuc_mix_comp_color.isDefault() ? ref_info.m_input.m_color : nuc_mix_comp_color;
        if( color.isDefault() && !ref_src.empty() )
          color = disp->suggestColorForSource( ref_src );
        if( color.isDefault() && !ref_src.empty() )
        {
          color = disp->nextGenericSourceColor();
          if( !color.isDefault() )
            disp->updateColorCacheForSource( ref_src, color );
        }
        if( !color.isDefault() )
          peak.setLineColor( color );
        break;
      }
    }
  }
  peak_model->addPeaks( peaks_to_add );
}

bool FitPeaksAdvancedWidget::canAccept() const
{
  return !m_is_calculating && m_current_result
         && m_current_result->status == RelActCalcAuto::RelActAutoSolution::Status::Success
         && !m_current_result->observable_peaks.empty();
}

Wt::Signal<> &FitPeaksAdvancedWidget::resultUpdated()
{
  return m_resultUpdated;
}

void FitPeaksAdvancedWidget::render( WFlags<RenderFlag> flags )
{
  if( m_render_flags.testFlag( UpdateCalculations ) )
  {
    startComputation();
  }
  m_render_flags = 0;
  WContainerWidget::render( flags );
}

void FitPeaksAdvancedWidget::startComputation()
{
  m_status->setText( WString::tr("fpn-calculating") );
  m_status->removeStyleClass( "calculating" );
  m_status->addStyleClass( "calculating" );
  m_warnings_div->removeStyleClass( "fpn-warnings-success" );
  m_warnings_div->removeStyleClass( "fpn-warnings-has-warnings" );
  m_warnings_div->removeStyleClass( "fpn-warnings-error" );
  m_warnings_div->addStyleClass( "fpn-warnings-computing" );
  while( m_warnings_div->count() > 0 )
    m_warnings_div->removeWidget( m_warnings_div->widget( 0 ) );
  if( m_chart_peak_model )
    m_chart_peak_model->removeAllPeaks();

  if( m_cancel_calc )
    m_cancel_calc->store( true );
  m_is_calculating = true;
  m_calc_number += 1;
  std::shared_ptr<std::atomic_bool> cancel_calc = std::make_shared<std::atomic_bool>( false );
  m_cancel_calc = cancel_calc;
  const size_t calc_number = m_calc_number;
  const std::string sessionid = wApp->sessionId();

  syncConfigFromOptions();
  WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> options = currentOptions();
  FitPeaksForNuclides::PeakFitForNuclideConfig config = currentConfig();

  // Peak detection requires InterSpec (GUI thread); run here and capture for worker
  if( m_auto_search_peaks->empty() )
  {
    InterSpec *viewer = InterSpec::instance();
    if( !viewer )
    {
      m_is_calculating = false;
      std::shared_ptr<std::string> err = std::make_shared<std::string>( "InterSpec instance gone" );
      handleCalcError( err, cancel_calc );
      return;
    }
    try
    {
      AnalystChecks::DetectedPeaksOptions peak_opts;
      peak_opts.specType = SpecUtils::SpectrumType::Foreground;
      peak_opts.nonBackgroundPeaksOnly = false;
      AnalystChecks::DetectedPeakStatus detected = AnalystChecks::detected_peaks( peak_opts, viewer );
      m_auto_search_peaks->assign( detected.peaks.begin(), detected.peaks.end() );
    }
    catch( std::exception &e )
    {
      m_is_calculating = false;
      std::shared_ptr<std::string> err = std::make_shared<std::string>( std::string( "Peak detection failed: " ) + e.what() );
      handleCalcError( err, cancel_calc );
      return;
    }
  }

  std::vector<std::shared_ptr<const PeakDef>> peaks_to_use( m_auto_search_peaks->begin(), m_auto_search_peaks->end() );
  std::shared_ptr<FitPeaksForNuclides::PeakFitResult> result = std::make_shared<FitPeaksForNuclides::PeakFitResult>();
  std::shared_ptr<std::string> error_msg = std::make_shared<std::string>();
  auto gui_update = wApp->bind( boost::bind( &FitPeaksAdvancedWidget::updateFromResult, this, result, cancel_calc, calc_number ) );
  auto gui_error = wApp->bind( boost::bind( &FitPeaksAdvancedWidget::handleCalcError, this, error_msg, cancel_calc ) );

  auto worker = [=](){
    try
    {
      *result = FitPeaksForNuclides::fit_peaks_for_nuclides( peaks_to_use, m_fg_copy, m_base_nucs, m_user_peaks,
                                                            m_bg_copy, m_drf, options, config, m_is_hpge );
      WServer::instance()->post( sessionid, [=](){
        WApplication *app = WApplication::instance();
        if( app )
        {
          gui_update();
          app->triggerUpdate();
        }
      } );
    }
    catch( std::exception &e )
    {
      *error_msg = std::string( "Fit failed: " ) + e.what();
      WServer::instance()->post( sessionid, [=](){
        WApplication *app = WApplication::instance();
        if( app )
        {
          gui_error();
          app->triggerUpdate();
        }
      } );
    }
    catch( ... )
    {
      *error_msg = "Fit failed with unknown error";
      WServer::instance()->post( sessionid, [=](){
        WApplication *app = WApplication::instance();
        if( app )
        {
          gui_error();
          app->triggerUpdate();
        }
      } );
    }
  };

  WServer::instance()->ioService().boost::asio::io_service::post( worker );
}

void FitPeaksAdvancedWidget::updateFromResult( std::shared_ptr<FitPeaksForNuclides::PeakFitResult> result,
                                               std::shared_ptr<std::atomic_bool> cancel_flag,
                                               size_t calc_number )
{
  if( cancel_flag != m_cancel_calc )
    return;
  if( cancel_flag && cancel_flag->load() )
    return;
  if( calc_number != m_calc_number )
    return;
  m_is_calculating = false;
  m_current_result = result;
  m_status->removeStyleClass( "calculating" );

  if( result->status != RelActCalcAuto::RelActAutoSolution::Status::Success )
  {
    m_status->setText( WString::tr("fpn-error-content").arg( result->error_message ) );
    m_warnings_div->show();
    m_warnings_div->removeStyleClass( "fpn-warnings-computing" );
    m_warnings_div->removeStyleClass( "fpn-warnings-success" );
    m_warnings_div->removeStyleClass( "fpn-warnings-has-warnings" );
    m_warnings_div->addStyleClass( "fpn-warnings-error" );
    while( m_warnings_div->count() > 0 )
      m_warnings_div->removeWidget( m_warnings_div->widget( 0 ) );
    WText *err_text = new WText( WString::fromUTF8( result->error_message ), m_warnings_div );
    (void)err_text;
  }
  else
  {
    m_warnings_div->removeStyleClass( "fpn-warnings-computing" );
    m_warnings_div->removeStyleClass( "fpn-warnings-error" );
    const RelActCalcAuto::RelActAutoSolution &sol = result->solution;
    WString summary;
    summary = WString::tr("fpn-summary-chi2").arg( WString::fromUTF8( std::to_string( sol.m_chi2 ) + " / " + std::to_string( sol.m_dof ) ) );
    if( sol.m_rel_activities.size() > 1 )
    {
      double max_act = 0;
      for( const auto &row : sol.m_rel_activities )
      {
        for( const auto &ra : row )
        {
          if( ra.rel_activity > max_act )
            max_act = ra.rel_activity;
        }
      }
      if( max_act > 0 )
      {
        std::string rel_str;
        for( const auto &row : sol.m_rel_activities )
        {
          for( const auto &ra : row )
          {
            if( !rel_str.empty() )
              rel_str += ", ";
            rel_str += ra.name() + "=" + std::to_string( ra.rel_activity / max_act );
          }
        }
        summary += WString::fromUTF8( " " );
        summary += WString::tr("fpn-summary-rel-activities").arg( WString::fromUTF8( rel_str ) );
      }
    }
    const size_t num_rois = result->fit_peaks.size();
    summary += WString::fromUTF8( " " );
    summary += WString::tr("fpn-summary-rois").arg( static_cast<int>( num_rois ) );
    m_status->setText( summary );
    if( result->warnings.empty() )
    {
      m_warnings_div->addStyleClass( "fpn-warnings-success" );
      while( m_warnings_div->count() > 0 )
        m_warnings_div->removeWidget( m_warnings_div->widget( 0 ) );
      WText *ok = new WText( WString::tr("fpn-computation-successful"), m_warnings_div );
      (void)ok;
    }
    else
    {
      m_warnings_div->addStyleClass( "fpn-warnings-has-warnings" );
      while( m_warnings_div->count() > 0 )
        m_warnings_div->removeWidget( m_warnings_div->widget( 0 ) );
      for( const std::string &w : result->warnings )
      {
        WText *wt = new WText( WString::fromUTF8( "• " + w ), m_warnings_div );
        wt->setInline( false );
        (void)wt;
      }
    }

    if( !m_chart_peak_model )
    {
      m_chart_peak_model = new PeakModel( m_chart );
      m_chart_peak_model->setNoSpecMeasBacking();
      m_chart_peak_model->setForeground( m_fg_copy );
      m_chart->setPeakModel( m_chart_peak_model );
    }
    std::vector<PeakDef> preview_peaks = result->observable_peaks;
    ReferencePhotopeakDisplay *disp = InterSpec::instance() ? InterSpec::instance()->referenceLinesWidget() : nullptr;
    for( PeakDef &peak : preview_peaks )
    {
      std::string src_name;
      if( peak.parentNuclide() )
        src_name = peak.parentNuclide()->symbol;
      else if( peak.xrayElement() )
        src_name = peak.xrayElement()->symbol;
      else if( peak.reaction() )
        src_name = peak.reaction()->name();
      if( src_name.empty() || !disp )
        continue;
      for( const ReferenceLineInfo &ref_info : m_ref_lines )
      {
        bool matches = false;
        Wt::WColor nuc_mix_comp_color;
        if( ref_info.m_nuclide && peak.parentNuclide() == ref_info.m_nuclide )
          matches = true;
        if( !matches && ref_info.m_element && peak.xrayElement() == ref_info.m_element )
          matches = true;
        if( !matches && peak.reaction() && ref_info.m_reactions.count( peak.reaction() ) > 0 )
          matches = true;
        if( !matches && peak.parentNuclide()
            && ref_info.m_source_type == ReferenceLineInfo::SourceType::NuclideMixture
            && ref_info.m_nuc_mix )
        {
          for( const ReferenceLinePredef::NucMixComp &comp : ref_info.m_nuc_mix->m_components )
          {
            if( comp.m_nuclide == peak.parentNuclide() )
            {
              matches = true;
              nuc_mix_comp_color = comp.m_color;
              break;
            }
          }
        }
        if( matches )
        {
          std::string ref_src;
          if( ref_info.m_nuclide )
            ref_src = ref_info.m_nuclide->symbol;
          else if( ref_info.m_element )
            ref_src = ref_info.m_element->symbol;
          else if( !ref_info.m_reactions.empty() )
          {
            const ReactionGamma::Reaction *rctn = *ref_info.m_reactions.begin();
            if( rctn )
              ref_src = rctn->name();
          }
          Wt::WColor color = nuc_mix_comp_color.isDefault() ? ref_info.m_input.m_color : nuc_mix_comp_color;
          if( color.isDefault() && !ref_src.empty() )
            color = disp->suggestColorForSource( ref_src );
          if( color.isDefault() && !ref_src.empty() )
          {
            color = disp->nextGenericSourceColor();
            if( !color.isDefault() )
              disp->updateColorCacheForSource( ref_src, color );
          }
          if( !color.isDefault() )
            peak.setLineColor( color );
          break;
        }
      }
    }
    m_chart_peak_model->addPeaks( preview_peaks );

    double low_e = 0, high_e = 0;
    bool first = true;
    for( const PeakDef &p : result->observable_peaks )
    {
      const double m = p.mean();
      if( first )
      {
        low_e = high_e = m;
        first = false;
      }
      else
      {
        if( m < low_e ) low_e = m;
        if( m > high_e ) high_e = m;
      }
    }
    if( !first )
    {
      const double margin = 0.1 * ( high_e - low_e );
      m_chart->setXAxisRange( low_e - margin, high_e + margin );
    }
  }
  m_resultUpdated.emit();
}

void FitPeaksAdvancedWidget::handleCalcError( std::shared_ptr<std::string> error_msg,
                                               std::shared_ptr<std::atomic_bool> cancel_flag )
{
  if( cancel_flag != m_cancel_calc )
    return;
  m_is_calculating = false;
  m_current_result = nullptr;
  m_status->removeStyleClass( "calculating" );
  m_status->setText( WString::tr("fpn-error-content").arg( WString::fromUTF8( *error_msg ) ) );
  m_warnings_div->removeStyleClass( "fpn-warnings-computing" );
  m_warnings_div->removeStyleClass( "fpn-warnings-success" );
  m_warnings_div->removeStyleClass( "fpn-warnings-has-warnings" );
  m_warnings_div->addStyleClass( "fpn-warnings-error" );
  while( m_warnings_div->count() > 0 )
    m_warnings_div->removeWidget( m_warnings_div->widget( 0 ) );
  WText *err_text = new WText( WString::fromUTF8( *error_msg ), m_warnings_div );
  (void)err_text;
  m_resultUpdated.emit();
}

std::string FitPeaksAdvancedWidget::sourceListTitle() const
{
  std::string out;
  for( size_t i = 0; i < m_base_nucs.size(); ++i )
  {
    if( i > 0 )
      out += ", ";
    out += RelActCalcAuto::to_name( m_base_nucs[i].source );
  }
  if( m_base_options.testFlag( FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaks ) )
  {
    if( !out.empty() )
      out += ", ";
    out += "background";
  }
  return out.empty() ? "sources" : out;
}

void FitPeaksAdvancedWidget::buildOptionsFromConfig()
{
  FitPeaksForNuclides::PeakFitForNuclideConfig config;

  InterSpec *viewer = InterSpec::instance();
  const bool show_tool_tips = viewer
    ? UserPreferences::preferenceValue<bool>( "ShowTooltips", viewer )
    : false;

  auto add_checkbox_row = [this, show_tool_tips]( const WString &label, WCheckBox *input,
                                                   const WString &tooltip = WString() ) {
    WContainerWidget *row = new WContainerWidget( m_options_div );
    row->addStyleClass( "fpn-option-row" );
    WLabel *lbl = new WLabel( label, row );
    lbl->addStyleClass( "fpn-option-label" );
    lbl->setBuddy( input );
    WContainerWidget *inp_div = new WContainerWidget( row );
    inp_div->addStyleClass( "fpn-option-input" );
    inp_div->addWidget( input );
    if( !tooltip.empty() )
      HelpSystem::attachToolTipOn( {lbl, input}, tooltip, show_tool_tips );
  };

  auto add_form_row = [this, show_tool_tips]( const WString &label, WFormWidget *input,
                                              const WString &tooltip = WString() ) {
    WContainerWidget *row = new WContainerWidget( m_options_div );
    row->addStyleClass( "fpn-option-row" );
    WLabel *lbl = new WLabel( label, row );
    lbl->addStyleClass( "fpn-option-label" );
    lbl->setBuddy( input );
    WContainerWidget *inp_div = new WContainerWidget( row );
    inp_div->addStyleClass( "fpn-option-input" );
    inp_div->addWidget( input );
    if( !tooltip.empty() )
      HelpSystem::attachToolTipOn( {lbl, input}, tooltip, show_tool_tips );
  };

  const bool has_existing_peaks = !m_user_peaks.empty();
  std::set<std::shared_ptr<const PeakContinuum>> distinct_rois;
  for( const std::shared_ptr<const PeakDef> &p : m_user_peaks )
  {
    if( p && p->continuum() )
      distinct_rois.insert( p->continuum() );
  }
  const bool has_existing_rois = !distinct_rois.empty();

  WContainerWidget *checkboxes_container = new WContainerWidget( m_options_div );
  checkboxes_container->addStyleClass( "fpn-options-checkboxes" );

  if( has_existing_rois )
  {
    m_opt_dont_use_rois = new WCheckBox( WString::tr("fpn-opt-dont-use-existing-rois"), checkboxes_container );
    m_opt_dont_use_rois->setChecked( false );
    m_opt_dont_use_rois->changed().connect( this, &FitPeaksAdvancedWidget::onDontUseRoisChanged );
    HelpSystem::attachToolTipOn( m_opt_dont_use_rois, WString::tr("fpn-opt-tt-dont-use-existing-rois"), show_tool_tips );
  }

  if( has_existing_peaks )
  {
    m_opt_existing_as_free = new WCheckBox( WString::tr("fpn-opt-existing-as-free"), checkboxes_container );
    m_opt_existing_as_free->setChecked( false );
    m_opt_existing_as_free->changed().connect( this, &FitPeaksAdvancedWidget::onExistingAsFreeChanged );
    HelpSystem::attachToolTipOn( m_opt_existing_as_free, WString::tr("fpn-opt-tt-existing-as-free"), show_tool_tips );
  }

  m_opt_dont_vary_energy_cal = new WCheckBox( WString::tr("fpn-opt-dont-vary-energy-cal"), checkboxes_container );
  m_opt_dont_vary_energy_cal->setChecked( false );
  m_opt_dont_vary_energy_cal->changed().connect( this, &FitPeaksAdvancedWidget::scheduleOptionsUpdate );
  m_opt_dont_vary_energy_cal->changed().connect( this, &FitPeaksAdvancedWidget::onDontVaryEnergyCalChanged );
  HelpSystem::attachToolTipOn( m_opt_dont_vary_energy_cal, WString::tr("fpn-opt-tt-dont-vary-energy-cal"), show_tool_tips );

  m_opt_dont_refine_energy_cal = new WCheckBox( WString::tr("fpn-opt-dont-refine-energy-cal"), checkboxes_container );
  m_opt_dont_refine_energy_cal->setChecked( false );
  m_opt_dont_refine_energy_cal->changed().connect( this, &FitPeaksAdvancedWidget::scheduleOptionsUpdate );
  HelpSystem::attachToolTipOn( m_opt_dont_refine_energy_cal, WString::tr("fpn-opt-tt-dont-refine-energy-cal"), show_tool_tips );

  onDontVaryEnergyCalChanged();

  m_opt_fit_bkgnd_peaks = new WCheckBox( WString::tr("fpn-opt-fit-bkgnd-peaks"), checkboxes_container );
  m_opt_fit_bkgnd_peaks->setChecked( m_base_options.testFlag( FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaks ) );
  m_opt_fit_bkgnd_peaks->changed().connect( this, &FitPeaksAdvancedWidget::scheduleOptionsUpdate );
  HelpSystem::attachToolTipOn( m_opt_fit_bkgnd_peaks, WString::tr("fpn-opt-tt-fit-bkgnd-peaks"), show_tool_tips );

  // m_opt_fit_bkgnd_dont_use: not exposed to user for now, may add later
  // m_opt_fit_bkgnd_dont_use = new WCheckBox( WString::tr("fpn-opt-fit-bkgnd-dont-use"), checkboxes_container );
  // m_opt_fit_bkgnd_dont_use->setChecked( false );
  // m_opt_fit_bkgnd_dont_use->changed().connect( this, &FitPeaksAdvancedWidget::scheduleOptionsUpdate );

  NativeFloatSpinBox *roi_chi2 = new NativeFloatSpinBox();
  roi_chi2->setWidth( WLength( 4, WLength::Unit::FontEm ) );
  roi_chi2->setSpinnerHidden( true );
  roi_chi2->setValue( static_cast<float>( config.roi_significance_min_chi2_reduction ) );
  roi_chi2->setRange( 0.1f, 1000.f );
  roi_chi2->valueChanged().connect( this, &FitPeaksAdvancedWidget::scheduleOptionsUpdate );
  add_form_row( WString::tr("fpn-opt-roi-min-chi2-red"), roi_chi2, WString::tr("fpn-opt-tt-roi-min-chi2-red") );
  m_opt_roi_min_chi2 = roi_chi2;

  NativeFloatSpinBox *roi_sig = new NativeFloatSpinBox();
  roi_sig->setWidth( WLength( 4, WLength::Unit::FontEm ) );
  roi_sig->setSpinnerHidden( true );
  roi_sig->setValue( static_cast<float>( config.roi_significance_min_peak_sig ) );
  roi_sig->setRange( 0.1f, 20.f );
  roi_sig->valueChanged().connect( this, &FitPeaksAdvancedWidget::scheduleOptionsUpdate );
  add_form_row( WString::tr("fpn-opt-roi-min-peak-sig"), roi_sig, WString::tr("fpn-opt-tt-roi-min-peak-sig") );
  m_opt_roi_min_peak_sig = roi_sig;

  NativeFloatSpinBox *obs_init = new NativeFloatSpinBox();
  obs_init->setWidth( WLength( 4, WLength::Unit::FontEm ) );
  obs_init->setSpinnerHidden( true );
  obs_init->setValue( static_cast<float>( config.observable_peak_initial_significance_threshold ) );
  obs_init->setRange( 0.1f, 20.f );
  obs_init->valueChanged().connect( this, &FitPeaksAdvancedWidget::scheduleOptionsUpdate );
  add_form_row( WString::tr("fpn-opt-obs-initial-sig"), obs_init, WString::tr("fpn-opt-tt-obs-initial-sig") );
  m_opt_obs_initial_sig = obs_init;

  NativeFloatSpinBox *obs_fin = new NativeFloatSpinBox();
  obs_fin->setWidth( WLength( 4, WLength::Unit::FontEm ) );
  obs_fin->setSpinnerHidden( true );
  obs_fin->setValue( static_cast<float>( config.observable_peak_final_significance_threshold ) );
  obs_fin->setRange( 0.1f, 20.f );
  obs_fin->valueChanged().connect( this, &FitPeaksAdvancedWidget::scheduleOptionsUpdate );
  add_form_row( WString::tr("fpn-opt-obs-final-sig"), obs_fin, WString::tr("fpn-opt-tt-obs-final-sig") );
  m_opt_obs_final_sig = obs_fin;

  m_opt_fwhm_form = new WComboBox();
  for( int i = 0; i <= static_cast<int>( RelActCalcAuto::FwhmForm::NotApplicable ); ++i )
  {
    RelActCalcAuto::FwhmForm f = static_cast<RelActCalcAuto::FwhmForm>( i );
    m_opt_fwhm_form->addItem( RelActCalcAuto::to_str( f ) );
    if( f == config.fwhm_form )
      m_opt_fwhm_form->setCurrentIndex( i );
  }
  m_opt_fwhm_form->changed().connect( this, &FitPeaksAdvancedWidget::scheduleOptionsUpdate );
  add_form_row( WString::tr("fpn-opt-fwhm-form"), m_opt_fwhm_form, WString::tr("fpn-opt-tt-fwhm-form") );

  m_opt_rel_eff_type = new WComboBox();
  for( int i = 0; i <= static_cast<int>( RelActCalc::RelEffEqnForm::FramPhysicalModel ); ++i )
  {
    RelActCalc::RelEffEqnForm f = static_cast<RelActCalc::RelEffEqnForm>( i );
    m_opt_rel_eff_type->addItem( RelActCalc::to_str( f ) );
    if( f == config.rel_eff_eqn_type )
      m_opt_rel_eff_type->setCurrentIndex( i );
  }
  m_opt_rel_eff_type->changed().connect( this, &FitPeaksAdvancedWidget::scheduleOptionsUpdate );
  m_opt_rel_eff_type->changed().connect( this, &FitPeaksAdvancedWidget::onRelEffTypeChanged );
  add_form_row( WString::tr("fpn-opt-rel-eff-eqn-type"), m_opt_rel_eff_type, WString::tr("fpn-opt-tt-rel-eff-eqn-type") );

  NativeFloatSpinBox *order_spin = new NativeFloatSpinBox();
  order_spin->setWidth( WLength( 4, WLength::Unit::FontEm ) );
  order_spin->setSpinnerHidden( true );
  order_spin->setValue( static_cast<float>( config.rel_eff_eqn_order ) );
  order_spin->setSingleStep( 1.f );
  order_spin->setRange( 0.f, 10.f );
  order_spin->setFormatString( "%.0f" );
  order_spin->valueChanged().connect( this, &FitPeaksAdvancedWidget::scheduleOptionsUpdate );
  m_rel_eff_order_row = new WContainerWidget( m_options_div );
  m_rel_eff_order_row->addStyleClass( "fpn-option-row" );
  WLabel *order_lbl = new WLabel( WString::tr("fpn-opt-rel-eff-eqn-order"), m_rel_eff_order_row );
  order_lbl->addStyleClass( "fpn-option-label" );
  order_lbl->setBuddy( order_spin );
  WContainerWidget *order_inp_div = new WContainerWidget( m_rel_eff_order_row );
  order_inp_div->addStyleClass( "fpn-option-input" );
  order_inp_div->addWidget( order_spin );
  HelpSystem::attachToolTipOn( {order_lbl, order_spin}, WString::tr("fpn-opt-tt-rel-eff-eqn-order"), show_tool_tips );
  m_opt_rel_eff_order = order_spin;
  m_rel_eff_order_row->setHidden( config.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel );
}

void FitPeaksAdvancedWidget::syncConfigFromOptions()
{
}

FitPeaksForNuclides::PeakFitForNuclideConfig FitPeaksAdvancedWidget::currentConfig() const
{
  FitPeaksForNuclides::PeakFitForNuclideConfig config;
  if( m_opt_roi_min_chi2 )
    config.roi_significance_min_chi2_reduction = m_opt_roi_min_chi2->value();
  if( m_opt_roi_min_peak_sig )
    config.roi_significance_min_peak_sig = m_opt_roi_min_peak_sig->value();
  if( m_opt_obs_initial_sig )
    config.observable_peak_initial_significance_threshold = m_opt_obs_initial_sig->value();
  if( m_opt_obs_final_sig )
    config.observable_peak_final_significance_threshold = m_opt_obs_final_sig->value();
  if( m_opt_fwhm_form )
    config.fwhm_form = static_cast<RelActCalcAuto::FwhmForm>( m_opt_fwhm_form->currentIndex() );
  if( m_opt_rel_eff_type )
    config.rel_eff_eqn_type = static_cast<RelActCalc::RelEffEqnForm>( m_opt_rel_eff_type->currentIndex() );
  if( m_opt_rel_eff_order )
    config.rel_eff_eqn_order = static_cast<size_t>( m_opt_rel_eff_order->value() );
  return config;
}

Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> FitPeaksAdvancedWidget::currentOptions() const
{
  WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> opts = m_base_options;
  if( m_opt_dont_use_rois && m_opt_dont_use_rois->isChecked() )
    opts |= FitPeaksForNuclides::FitSrcPeaksOptions::DoNotUseExistingRois;
  if( m_opt_existing_as_free && m_opt_existing_as_free->isChecked() )
    opts |= FitPeaksForNuclides::FitSrcPeaksOptions::ExistingPeaksAsFreePeak;
  if( m_opt_dont_vary_energy_cal && m_opt_dont_vary_energy_cal->isChecked() )
    opts |= FitPeaksForNuclides::FitSrcPeaksOptions::DoNotVaryEnergyCal;
  if( m_opt_dont_refine_energy_cal && m_opt_dont_refine_energy_cal->isChecked() && !m_opt_dont_refine_energy_cal->isDisabled() )
    opts |= FitPeaksForNuclides::FitSrcPeaksOptions::DoNotRefineEnergyCal;
  if( m_opt_fit_bkgnd_peaks && m_opt_fit_bkgnd_peaks->isChecked() )
    opts |= FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaks;
  // if( m_opt_fit_bkgnd_dont_use && m_opt_fit_bkgnd_dont_use->isChecked() )
  //   opts |= FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaksDontUse;
  return opts;
}

void FitPeaksAdvancedWidget::scheduleOptionsUpdate()
{
  m_render_flags |= UpdateCalculations;
  scheduleRender();
}

void FitPeaksAdvancedWidget::onDontUseRoisChanged()
{
  if( m_opt_dont_use_rois && m_opt_dont_use_rois->isChecked() && m_opt_existing_as_free )
    m_opt_existing_as_free->setChecked( false );
  scheduleOptionsUpdate();
}

void FitPeaksAdvancedWidget::onExistingAsFreeChanged()
{
  if( m_opt_existing_as_free && m_opt_existing_as_free->isChecked() && m_opt_dont_use_rois )
    m_opt_dont_use_rois->setChecked( false );
  scheduleOptionsUpdate();
}

void FitPeaksAdvancedWidget::onDontVaryEnergyCalChanged()
{
  if( m_opt_dont_refine_energy_cal )
    m_opt_dont_refine_energy_cal->setDisabled( m_opt_dont_vary_energy_cal && m_opt_dont_vary_energy_cal->isChecked() );
}

void FitPeaksAdvancedWidget::onRelEffTypeChanged()
{
  const bool is_physical = ( m_opt_rel_eff_type && m_opt_rel_eff_type->currentIndex()
                             == static_cast<int>( RelActCalc::RelEffEqnForm::FramPhysicalModel ) );
  if( m_rel_eff_order_row )
    m_rel_eff_order_row->setHidden( is_physical );
}
