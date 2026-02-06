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

#include <string>
#include <vector>
#include <memory>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <Wt/WText>
#include <Wt/WServer>
#include <Wt/WCheckBox>
#include <Wt/WIOService>
#include <Wt/WPushButton>
#include <Wt/WApplication>

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/AnalystChecks.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/FitPeaksForNuclides.h"
#include "InterSpec/ReferenceLinePredef.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/PhysicalUnitsLocalized.h"
#include "InterSpec/FitPeaksForNuclidesGui.h"
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

  SimpleDialog *dlg = new SimpleDialog( WString::tr("fpn-adv-title"),
                                        WString::tr("fpn-adv-content") );
  dlg->rejectWhenEscapePressed();

  WPushButton *fit_btn = dlg->addButton( WString::tr("fpn-fit-peaks") );
  dlg->addButton( WString::tr("Cancel") );

  // Note: SimpleDialog will close and delete itself when either button is clicked.
  fit_btn->clicked().connect( boost::bind( &startFitSources, true ) );
  // Cancel button uses SimpleDialog's default close behavior.

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
