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

#include <deque>
#include <string>
#include <memory>
#include <functional>

#include <Wt/WText>
#include <Wt/WTable>
#include <Wt/WPanel>
#include <Wt/WServer>
#include <Wt/WPainter>
#include <Wt/WComboBox>
#include <Wt/WCheckBox>
#include <Wt/WSvgImage>
#include <Wt/WAnimation>
#include <Wt/WIOService>
#include <Wt/WTableCell>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/Chart/WDataSeries>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/MakeFwhmForDrf.h"
#include "InterSpec/RefLineDynamic.h"
#include "InterSpec/PeakInfoDisplay.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/SpectrumDataModel.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/PhysicalUnitsLocalized.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"

#if( USE_REL_ACT_TOOL )
#include "InterSpec/RelEffChart.h"
#endif

using namespace Wt;
using namespace std;

namespace
{
  /** The user action that caused a PeakSelectorWindow to be made.
   This influences how things are presented to the user.
   */
  enum class PeakSelectorWindowReason
  {
    /** User hit the "Search For Peaks" button in the Peak Manager */
    PeakSearch,
    
    /** User hit the "Nuc. from Ref." button in the Peak Manager */
    NuclideId,
    
    /** User loaded a new spectrum from same detector as last one, and had the
     "Ask to Propagate Peaks" option selected.
     
     In this case the 'orig_peaks' passed into the window constructor are the
     candidate peaks to try and fit for.
     */
    PeaksFromPreviousSpectrum,
    
    /** User uploaded a peaks CSV file.
     */
    PeaksFromCsvFile,
  };//enum class PeakSelectorWindowReason
  
  
/** Class representing a dialog that allows the user to select which peaks to
 keep, or nuclides to assign to a peak after an automated peak fitting, nuclide
 ID, or peak propagation function.
 
 Currently this window can be triggered from four different user actions,
 leading to a base of four different ways to display things (there are also some
 dynamic display changes based on data and results), which adds some complexity,
 but maybe at the moment this is more manageable than breaking this class up...
 */
class PeakSelectorWindow : public AuxWindow
{
  InterSpec *m_viewer;
  Wt::WTable *m_table;
  int m_previewChartColumn;
  Wt::WCheckBox *m_keepRefLinePeaksOnly;
  Wt::WCheckBox *m_showAllPeaks;
  
#if( USE_REL_ACT_TOOL )
  Wt::WPanel *m_chartPanel;
  RelEffChart *m_rel_eff_chart;
#endif
  
  const PeakSelectorWindowReason m_reason;
  const vector<PeakDef> m_orig_peaks;
  std::shared_ptr<const SpecUtils::Measurement> m_data;
  const vector<PeakDef> m_final_peaks;
  vector< shared_ptr<const ReferenceLineInfo>> m_displayed;
  vector<WCheckBox *> m_keep_peak_cbs;
  vector<WComboBox *> m_nuc_select_combos;
  vector<WCheckBox *> m_dont_change_nuc_cbs;
  vector<bool> m_peak_updated;
  
  /** For each entry the first peak is a previously existing peak (before either
   the nuclide ID or peak search); the second peak is the new peak (either found
   peak, or old peak with new nuclide assigned).  If both the first and second
   peak are valid, then they ciirespond to eachother (e.g., really similar mean).
   The second peak will typically be valid (but you should check), but the first
   peak wont be for newly fit peaks.
   */
  vector<pair<shared_ptr<PeakDef>,shared_ptr<PeakDef>>> m_old_to_new_peaks;
  
  /** If true when doFinish() is called, revert back to how things were before
     searching for peaks or doing ID.
   */
  bool m_cancelOperation;
  
public:
  /** PeakSelectorWindow constructor.
   
   @param viewer The InterSpec instance to create this dialog for.
   @param reason The source of the user action that triggered making this dialog.
   @param orig_peaks For reasons of PeakSearch, NuclideId, and PeaksFromCsvFile,
          these are the peaks that existed before the fit.
          For PeaksFromPreviousSpectrum, these are the peaks that were fed into
          the fitter.
   @param data The gamma spectrum the peaks are for.
   @param final_peaks All peaks; both previously existing and fit for.
   @param displayed The displayed reference lines.  If provided will add user
          option to the GUI to keep only peaks assigned to one of the reference
          lines.
   */
  PeakSelectorWindow( InterSpec *viewer,
                      const PeakSelectorWindowReason reason,
                      const vector<PeakDef> &orig_peaks,
                      std::shared_ptr<const SpecUtils::Measurement> data,
                      const vector<PeakDef> &final_peaks,
                      const vector<ReferenceLineInfo> &displayed )
  : AuxWindow( "Dummy",
               (AuxWindowProperties::IsModal | AuxWindowProperties::TabletNotFullScreen | AuxWindowProperties::DisableCollapse) ),
    m_viewer( viewer ),
    m_table( nullptr ),
    m_previewChartColumn( -1 ),
    m_keepRefLinePeaksOnly( nullptr ),
    m_showAllPeaks( nullptr ),
#if( USE_REL_ACT_TOOL )
    m_chartPanel( nullptr ),
    m_rel_eff_chart( nullptr ),
#endif
    m_reason( reason ),
    m_orig_peaks( orig_peaks ),
    m_data( data ),
    m_final_peaks( final_peaks ),
    m_displayed(),
    m_cancelOperation( false )
  {
    wApp->useStyleSheet( "InterSpec_resources/PeakSelectorWindow.css" );
   
    if( !m_viewer )
      m_viewer = InterSpec::instance();
    assert( m_viewer );
    if( m_viewer )
      m_viewer->useMessageResourceBundle( "PeakSearchGuiUtils" );
    
    switch( m_reason )
    {
      case PeakSelectorWindowReason::NuclideId:
        setWindowTitle( WString::tr("psw-title-confirm-nuc-ass") );
        break;
      case PeakSelectorWindowReason::PeakSearch:
        setWindowTitle( WString::tr("psw-title-check-search") );
        break;
      case PeakSelectorWindowReason::PeaksFromPreviousSpectrum:
        setWindowTitle( WString::tr("psw-title-prev-spec") );
        break;
      case PeakSelectorWindowReason::PeaksFromCsvFile:
        setWindowTitle( WString::tr("psw-title-csv-import") );
        break;
    }//switch( m_reason )
    
    for( const auto &d : displayed )
      m_displayed.push_back( make_shared<ReferenceLineInfo>(d) );
    
    const int appWidth = viewer->renderedWidth();
    const int appHeight = viewer->renderedHeight();
    if( appWidth > 100 && appHeight > 100 )  //should probably always evaluate to true
      setMaximumSize( 0.95*appWidth, 0.75*appHeight );
    
    //We need to map the initial peaks to the final peaks
    vector<bool> used_new( final_peaks.size(), false );
    for( size_t oldindex = 0; oldindex < orig_peaks.size(); ++oldindex )
    {
      const PeakDef &oldPeak = orig_peaks[oldindex];
      //Just match peaks up by mean for now.
      const double oldMean = oldPeak.mean();
      const double oldSigma = oldPeak.gausPeak() ? oldPeak.sigma() : 0.25*oldPeak.roiWidth();
      
      int nearestNew = -1;
      for( size_t newindex = 0; newindex < final_peaks.size(); ++newindex )
      {
        if( used_new[newindex] )
          continue;
        const PeakDef &newPeak = final_peaks[newindex];
        const double newMean = newPeak.mean();
        const double newSigma = newPeak.gausPeak() ? newPeak.sigma() : 0.25*newPeak.roiWidth();
        if( fabs(newMean - oldMean) < 0.25*(newSigma + oldSigma) )
        {
          if( (nearestNew < 0)
             || (fabs(final_peaks[newindex].mean() - oldMean) < fabs(final_peaks[nearestNew].mean() - oldMean)) )
          {
            nearestNew = static_cast<int>(newindex);
          }
        }
      }
      
      std::shared_ptr<PeakDef> newpeak;
      if( nearestNew >= 0 )
      {
        used_new[nearestNew] = true;
        newpeak = make_shared<PeakDef>( final_peaks[nearestNew] );
      }
      
      m_old_to_new_peaks.emplace_back( make_shared<PeakDef>(oldPeak), newpeak );
    }//for( size_t oldindex = 0; oldindex < orig_peaks.size(); ++oldindex )
    
    for( size_t newindex = 0; newindex < final_peaks.size(); ++newindex )
    {
      if( !used_new[newindex] )
        m_old_to_new_peaks.emplace_back( nullptr, make_shared<PeakDef>(final_peaks[newindex]) );
    }//for( size_t newindex = 0; newindex < final_peaks.size(); ++newindex )
    
    std::sort( begin(m_old_to_new_peaks), end(m_old_to_new_peaks),
              []( const pair<shared_ptr<PeakDef>,shared_ptr<PeakDef>> &lhs, const pair<shared_ptr<PeakDef>,shared_ptr<PeakDef>> &rhs ) -> bool {
                const auto &lhsp = lhs.second ? lhs.second : lhs.first;
                const auto &rhsp = rhs.second ? rhs.second : rhs.first;
                if( !lhsp || !rhsp ) return false;
                return lhsp->mean() < rhsp->mean();
              } );
  

    WText *txt = nullptr;
    
    switch( m_reason )
    {
      case PeakSelectorWindowReason::PeakSearch:
        if( m_displayed.size() )
          txt = new WText( WString::tr("psw-search-instruct-adjust"), contents() );
        else
          txt = new WText( WString::tr("psw-search-instruct-unselect"), contents() );
        break;
        
      case PeakSelectorWindowReason::NuclideId:
        txt = new WText( WString::tr("psw-search-instruct-assign"), contents() );
        break;
        
      case PeakSelectorWindowReason::PeaksFromPreviousSpectrum:
      case PeakSelectorWindowReason::PeaksFromCsvFile:
        txt = new WText( WString::tr("psw-search-instruct-keep"), contents() );
        break;
    }//switch( m_reason )
    
   
    txt->setInline( false );
    
    switch( m_reason )
    {
      case PeakSelectorWindowReason::NuclideId:
        m_showAllPeaks = new WCheckBox( WString::tr("psw-show-all-peaks"), contents() );
        break;
        
      case PeakSelectorWindowReason::PeakSearch:
      case PeakSelectorWindowReason::PeaksFromCsvFile:
      {
        if( !m_displayed.empty() )
        {
          bool anyNonAssignedPeaks = false;
          for( size_t i = 0; !anyNonAssignedPeaks && i < m_old_to_new_peaks.size(); ++i )
            if( m_old_to_new_peaks[i].second )
              anyNonAssignedPeaks = !m_old_to_new_peaks[i].second->hasSourceGammaAssigned();
          
          if( anyNonAssignedPeaks )
          {
            string srcs;
            for( size_t i = 0; i < m_displayed.size(); ++i )
            {
              if( i && ((i+1) == m_displayed.size()) )
                srcs += " or ";
              else if( i )
                srcs += ", ";
              srcs += m_displayed[i]->m_input.m_input_txt;
            }
            WString msg = WString::tr("psw-keep-only-assigned-peaks").arg(srcs);
            m_keepRefLinePeaksOnly = new WCheckBox( msg, contents() );
            m_keepRefLinePeaksOnly->setInline( false );
            m_keepRefLinePeaksOnly->checked().connect( this, &PeakSelectorWindow::keepOnlyRefLinesCbChanged );
            m_keepRefLinePeaksOnly->unChecked().connect( this, &PeakSelectorWindow::keepOnlyRefLinesCbChanged );
            m_keepRefLinePeaksOnly->setMargin( 10, Wt::Top );
            m_keepRefLinePeaksOnly->addStyleClass( "KeepRefLinPeaksOnlyCb CbNoLineBreak" );
          }//if( anyNonAssignedPeaks )
        }//if( we searched for peaks, and there were reference lines )
        
        if( orig_peaks.size() )
          m_showAllPeaks = new WCheckBox( WString::tr("psw-show-prev-peaks"), contents() );
        break;
      }//case PeakSelectorWindowReason::PeakSearch:
        
      case PeakSelectorWindowReason::PeaksFromPreviousSpectrum:
        if( final_peaks.size() != orig_peaks.size() )
          m_showAllPeaks = new WCheckBox( WString::tr("psw-show-missing-peaks"), contents() );
        break;
    }//switch ( m_reason )
    
    if( m_showAllPeaks )
    {
      m_showAllPeaks->addStyleClass( "ShowAllPeaksCb" );
      m_showAllPeaks->setInline( false );
      m_showAllPeaks->setChecked( false );
      m_showAllPeaks->checked().connect( this, &PeakSelectorWindow::showAllPeaksCbChanged );
      m_showAllPeaks->unChecked().connect( this, &PeakSelectorWindow::showAllPeaksCbChanged );
      if( !m_keepRefLinePeaksOnly )
        m_showAllPeaks->setMargin( 10, Wt::Top );
    }//if( m_showAllPeaks )
    
    
#if( USE_REL_ACT_TOOL )
    // Potentially create a rel eff chart to help the user decide about interferences and such
    setupRelEffChart();
#endif
    contents()->setOverflow( WContainerWidget::Overflow::OverflowAuto, Orientation::Vertical );
    m_table = new WTable( contents() );
    m_table->addStyleClass( "PeakSelectorTable" );
    m_table->setHeaderCount( 1, Wt::Horizontal );
    
    m_keep_peak_cbs.resize( m_old_to_new_peaks.size(), nullptr );
    m_nuc_select_combos.resize( m_old_to_new_peaks.size(), nullptr );
    m_dont_change_nuc_cbs.resize( m_old_to_new_peaks.size(), nullptr );
    m_peak_updated.resize( m_old_to_new_peaks.size(), false );
    
    auto nucs_changed = [this]( size_t index ){
      assert( index < m_old_to_new_peaks.size() );
      
      switch( m_reason )
      {
        case PeakSelectorWindowReason::NuclideId:
        case PeakSelectorWindowReason::PeakSearch:
        case PeakSelectorWindowReason::PeaksFromCsvFile:
          break;
          
        case PeakSelectorWindowReason::PeaksFromPreviousSpectrum:
          return true;
      }//switch( m_reason )
      
      shared_ptr<PeakDef> &oldpeak = m_old_to_new_peaks[index].first;
      shared_ptr<PeakDef> &newpeak = m_old_to_new_peaks[index].second;
      
      const bool has_changed = ( (!oldpeak || !newpeak)
                                || (oldpeak->parentNuclide() != newpeak->parentNuclide())
                                || (oldpeak->decayParticle() != newpeak->decayParticle())
                                || (oldpeak->xrayElement() != newpeak->xrayElement())
                                || (oldpeak->xrayElement() && (oldpeak->xrayEnergy() != newpeak->xrayEnergy()))
                                || (oldpeak->reaction() != newpeak->reaction())
                                || (oldpeak->reaction() && (oldpeak->reactionEnergy() != newpeak->reactionEnergy())) );
      return has_changed;
    };//nucs_changed();
    
    
    int keepPeakIndex = 0, peakEnergyIndex = 1, origColumnIndex = 2, newColumnIndex = 3, previewIndex = 4;
    //Check if any reference lines are showing, and if not hide columns 2 and 3
    if( displayed.empty() )
    {
      switch( m_reason )
      {
        case PeakSelectorWindowReason::PeakSearch:
          keepPeakIndex   = 0;
          peakEnergyIndex = 1;
          previewIndex    = 2;
          newColumnIndex = origColumnIndex = -1;
          break;
          
        case PeakSelectorWindowReason::NuclideId:
          //Actually we shouldnt be here!  But lets set indexes rather than handling this potential logic error for now.
          peakEnergyIndex = 0;
          origColumnIndex = 1;
          newColumnIndex  = 2;
          previewIndex    = 3;
          keepPeakIndex   = -1;
          break;
          
        case PeakSelectorWindowReason::PeaksFromPreviousSpectrum:
        case PeakSelectorWindowReason::PeaksFromCsvFile:
          keepPeakIndex   = 0;
          peakEnergyIndex = 1;
          origColumnIndex = -1;
          newColumnIndex  = 2;
          previewIndex    = 3;
          break;
      }//switch( m_reason )
    }else
    {
      switch( m_reason )
      {
        case PeakSelectorWindowReason::PeakSearch:
          keepPeakIndex   = 0;
          peakEnergyIndex = 1;
          origColumnIndex = 2;
          newColumnIndex  = 3;
          previewIndex    = 4;
        break;
          
        case PeakSelectorWindowReason::NuclideId:
          keepPeakIndex   = -1;
          peakEnergyIndex = 0;
          origColumnIndex = 1;
          newColumnIndex  = 2;
          previewIndex    = 3;
        break;
          
        case PeakSelectorWindowReason::PeaksFromPreviousSpectrum:
        case PeakSelectorWindowReason::PeaksFromCsvFile:
          //Actually we shouldnt be here!  But lets set indexes rather than handling this potential logic error for now.
          keepPeakIndex   = 0;
          peakEnergyIndex = 1;
          origColumnIndex = -1;
          newColumnIndex  = 2;
          previewIndex    = 3;
          break;
      }//switch( m_reason )
    }//if( displayed.empty() ) / else

    bool anyPrevPeaksHadNucs = false;
    for( size_t i = 0; i < m_old_to_new_peaks.size(); ++i )
    {
      shared_ptr<PeakDef> &oldpeak = m_old_to_new_peaks[i].first;
      if( nucs_changed(i) && oldpeak )
        anyPrevPeaksHadNucs = (anyPrevPeaksHadNucs || (oldpeak->decayParticle() || oldpeak->xrayElement() || oldpeak->reaction()));
    }//for( size_t i = 0; i < m_old_to_new_peaks.size(); ++i )
    
    if( (origColumnIndex >= 0) && !anyPrevPeaksHadNucs )
    {
      origColumnIndex = -1;
      newColumnIndex = (newColumnIndex>=0) ? newColumnIndex-1 : -1;
      previewIndex = (previewIndex>=0) ? previewIndex-1 : -1;
    }
    
    m_previewChartColumn = previewIndex;
    
    WTableCell *cell = nullptr;
    
    if( keepPeakIndex >= 0 )
    {
      cell = m_table->elementAt(0,keepPeakIndex);
      // We will add in some space so the "Keep Peak" text will have enough room to be next to the
      //  actual check box, and they wont be on separate lines.
      txt = new WText( WString::tr("psw-keep-peak-q"), cell );
      txt->setWordWrap( false );
    }
    
    if( peakEnergyIndex >= 0 )
    {
      cell = m_table->elementAt(0,peakEnergyIndex);
      WString labeltxt = WString("{1}, {2}").arg( WString::tr("Energy") ).arg( WString::tr("FWHM") );
      new WText( labeltxt, cell );
    }
    
    if( origColumnIndex >= 0 )
    {
      cell = m_table->elementAt(0,origColumnIndex);
      new WText( WString::tr("psw-orig-nuc"), cell );
    }
    
    if( newColumnIndex >= 0 )
    {
      cell = m_table->elementAt(0,newColumnIndex);
      new WText( WString::tr("psw-assigned-nuc"), cell );
    }
    
    if( previewIndex >= 0 )
    {
      cell = m_table->elementAt(0,previewIndex);
      new WText( WString::tr("psw-peak-preview"), cell );
    }
    
    
    bool some_nuclides_changed = false, all_peaks_modified = true;
    for( size_t i = 0; i < m_old_to_new_peaks.size(); ++i )
    {
      const int table_row = static_cast<int>(i+1);
      
      const bool has_changed = nucs_changed(i);
      m_peak_updated[i] = has_changed;
      
      some_nuclides_changed = (some_nuclides_changed || has_changed);
      
      if( !has_changed )
      {
        all_peaks_modified = false;
        WTableRow *row = m_table->rowAt(table_row);
        row->hide();
      }//if( !has_changed )
      

      if( (keepPeakIndex>=0) )
      {
        WTableCell *cbcell = m_table->elementAt( table_row, keepPeakIndex );
        cbcell->addStyleClass( "KeepPeakCell" );
        
        if( (m_reason == PeakSelectorWindowReason::PeaksFromPreviousSpectrum)
            && !m_old_to_new_peaks[i].second )
        {
          new WText( WString::tr("psw-not-found"), cbcell );
        }else if( has_changed )
        {
          WCheckBox *cb = new WCheckBox( WString::tr("psw-keep-peak"), cbcell );
          cb->setWordWrap(false);
          cb->addStyleClass( "CbNoLineBreak" );
          cb->setChecked(true);
          cb->checked().connect( this, &PeakSelectorWindow::keepPeakChanged );
          cb->unChecked().connect( this, &PeakSelectorWindow::keepPeakChanged );
          m_keep_peak_cbs[i] = cb;
        }
      }//
      
      if( origColumnIndex >= 0 )
      {
        const string origstr = makeSourceDesciption( m_old_to_new_peaks[i].first );
        WTableCell *origNucCell = m_table->elementAt( table_row, origColumnIndex );
        origNucCell->addStyleClass( "OrigNucCell" );
        new WText( origstr, origNucCell );
      }
      
      if( newColumnIndex >= 0 )
      {
        WTableCell *newNucCell = m_table->elementAt( table_row, newColumnIndex );
        newNucCell->addStyleClass( "NewNucCell" );
        
        if( has_changed
            && m_reason!=PeakSelectorWindowReason::PeaksFromPreviousSpectrum
            && m_reason!=PeakSelectorWindowReason::PeaksFromCsvFile )
        {
          m_nuc_select_combos[i] = new WComboBox(newNucCell);
          m_nuc_select_combos[i]->setInline( false );
          m_nuc_select_combos[i]->changed().connect( boost::bind( &PeakSelectorWindow::nucSelectChanged, this, i ) );
          
          // Only show the "Don't change" checkbox if this is a previously existing peak
          if( m_old_to_new_peaks[i].first )
          {
            m_dont_change_nuc_cbs[i] = new WCheckBox( WString::tr("psw-dont-change"), newNucCell );
            m_dont_change_nuc_cbs[i]->addStyleClass( "DontAssignCb" );
            m_dont_change_nuc_cbs[i]->setInline( false );
            m_dont_change_nuc_cbs[i]->checked().connect( boost::bind( &PeakSelectorWindow::dontChangeNucCbChanged, this, i ) );
            m_dont_change_nuc_cbs[i]->unChecked().connect( boost::bind( &PeakSelectorWindow::dontChangeNucCbChanged, this, i ) );
          }//if( m_old_to_new_peaks[i].first )
        }else
        {
          const string origstr = makeSourceDesciption( m_old_to_new_peaks[i].second );
          new WText( origstr, newNucCell );
        }
      }//if( newColumnIndex >= 0 )
      
      auto finalpeak = m_old_to_new_peaks[i].second;
      auto originalpeak = m_old_to_new_peaks[i].first;
      auto displayPeak = finalpeak ? finalpeak : originalpeak;
      
      if( displayPeak && (peakEnergyIndex >= 0) )
      {
        char buffer[512] = { '\0' };
        
        if( displayPeak->gausPeak() )
        {
          snprintf( buffer, sizeof(buffer), "<div style=\"white-space: nowrap;\">"
                   "mean=%.1f keV</div><div style=\"white-space: nowrap;\">FWHM=%.1f keV</div>",
                   displayPeak->mean(), displayPeak->fwhm() );
        }else
        {
          snprintf( buffer, sizeof(buffer), "<div style=\"white-space: nowrap;\">"
                   "mean=%.1f keV</div><div style=\"white-space: nowrap;\">Data Defined</div>",
                   displayPeak->mean() );
        }
        
        WTableCell *peakene = m_table->elementAt( table_row, peakEnergyIndex );
        peakene->addStyleClass( "PeakEnergyCell" );
        new WText( buffer, peakene );
      }
      
      if( displayPeak && (previewIndex>=0) )
      {
        updatePreviewPlot( i );
      }//if( peak )
    }//for( size_t i = 0; i < peaks.size(); ++i )
    
    if( m_showAllPeaks )
      m_showAllPeaks->setHidden( all_peaks_modified );
    
    populateNuclideSelects();
    
    switch( m_reason )
    {
      case PeakSelectorWindowReason::PeakSearch:
      case PeakSelectorWindowReason::PeaksFromPreviousSpectrum:
      case PeakSelectorWindowReason::PeaksFromCsvFile:
        break;
        
      case PeakSelectorWindowReason::NuclideId:
        if( !some_nuclides_changed )
        {
          txt = new WText( WString("<strong>{1}</strong>").arg( WString::tr("psw-no-nuc-changed") ), contents() );
          txt->setPadding( 10, Wt::Top | Wt::Bottom );
          txt->setTextAlignment( Wt::AlignmentFlag::AlignCenter );
          txt->setInline( false );
        }
        break;
    }//switch( m_reason )

    
    keepPeakChanged();
#if( USE_REL_ACT_TOOL )
    refreshRelEffChart();
#endif
    
    WPushButton *acceptButton = addCloseButtonToFooter( WString::tr("Accept"), true );
    acceptButton->clicked().connect( boost::bind( &AuxWindow::hide, this ) );
    
    acceptButton->clicked().connect( std::bind( [viewer, orig_peaks, final_peaks](){
      
      auto undo = [viewer, orig_peaks](){
        PeakModel *pmodel = viewer->peakModel();
        if( pmodel )
          pmodel->setPeaks( orig_peaks );
      };
      
      auto redo = [viewer, final_peaks](){
        PeakModel *pmodel = viewer->peakModel();
        if( pmodel )
          pmodel->setPeaks( final_peaks );
      };
      
      UndoRedoManager *undoManager = viewer->undoRedoManager();
      if( undoManager )
        undoManager->addUndoRedoStep( undo, redo, "Accept automated peak search." );
    } ) );
                                               
    
    
    WPushButton *cancelButton = nullptr;
    if( viewer->isPhone() )
    {
      cancelButton = new WPushButton( WString::tr("Cancel"), footer() );
      cancelButton->setFloatSide( Wt::Right );
    }else
    {
      cancelButton = addCloseButtonToFooter( WString::tr("Cancel"), true );
    }
    
    cancelButton->clicked().connect( boost::bind( &PeakSelectorWindow::cancelOperation, this ) );
    
    finished().connect( this, &PeakSelectorWindow::doFinish );
    
    show();
    
    resizeToFitOnScreen();
    centerWindow();
    centerWindowHeavyHanded();
    
    
    {// Begin handle undo when the dialog is showing, and the user hasnt accepted it
      auto cancel = wApp->bind( boost::bind( &PeakSelectorWindow::cancelOperation, this ) );
      
      auto undo = [cancel](){
        cancel();
      };
      auto redo = [](){
        //We just wont do anything - we could create a PeakSelectorWindow, but we will
        //  have to avoid undo<-->redo cycles, etc, so it makes the most sense to just
        //  skip redo, and use that effort making other actions better.
      };
      
      UndoRedoManager *undoManager = viewer->undoRedoManager();
      if( undoManager )
        undoManager->addUndoRedoStep( undo, redo, "Reject automated peak search." );
    }// End handle undo when the dialog is showing, and the user hasnt accepted it
  }//PeakSelector constructor
  
  //Only consider normal nuclides, and normal gammas (n xrays, or escape peaks)
  static bool eligible_for_rel_eff_chart( const shared_ptr<PeakDef> &p )
  {
    return (p
            && p->parentNuclide()
            && ((p->sourceGammaType() == PeakDef::NormalGamma)
                || (p->sourceGammaType() == PeakDef::XrayGamma)));
  };//eligible_for_rel_eff_chart(...)
  
#if( USE_REL_ACT_TOOL )
  /** The rel eff chart is close, but not fully debugged */
  void setupRelEffChart()
  {
    //First check to see if there is at least two peaks for a nuclide
    //  This next loop only looks at new peaks, or peaks that may have changed nuclide assignment.
    size_t max_peaks = 0;
    set<const SandiaDecay::Nuclide *> new_peak_nucs;
    for( const auto &old_new : m_old_to_new_peaks )
    {
      // The new peak should be eligible to be in the Rel. Eff. chart, and the nuclide, or at least
      //   the specifically attributed gamma/x-ray, should have changed.
      if( eligible_for_rel_eff_chart(old_new.second)
         && (!old_new.first
             || (old_new.first->decayParticle() != old_new.second->decayParticle()) ) )
      {
        new_peak_nucs.insert( old_new.second->parentNuclide() );
      }
    }//for( const auto &old_new : m_old_to_new_peaks )
    
    map<const SandiaDecay::Nuclide *, size_t> num_peaks;
    for( size_t i = 0; i < m_old_to_new_peaks.size(); ++i )
    {
      auto p = m_old_to_new_peaks[i].second;

      //Only consider normal nuclides, and normal gammas (n xrays, or escape peaks)
      if( eligible_for_rel_eff_chart(p) && new_peak_nucs.count(p->parentNuclide()) )
      {
        assert( (num_peaks.find(p->parentNuclide()) != end(num_peaks)) || (num_peaks[p->parentNuclide()] == 0) );
        size_t &npeaks = num_peaks[p->parentNuclide()]; //will be initialized to zero at first call of num_peaks[]
        max_peaks = std::max( max_peaks, ++npeaks );
      }
    }//for( size_t i = 0; i < m_old_to_new_peaks.size(); ++i )
    
    // We will only create a Relative Efficiency chart only if we have exactly one nuclide
    //  in the newly found peaks with two or more peaks (.
    if( (max_peaks < 2) || (num_peaks.size() != 1) )
      return;
    
    m_chartPanel = new WPanel( contents() );
    m_chartPanel->addStyleClass( "PeakSelectRelEffPanel" );
    m_chartPanel->setCollapsible( true );
    m_chartPanel->setCollapsed( true );
    m_chartPanel->setTitle( WString::tr("psw-rel-eff-plot-title") );
    m_chartPanel->setAnimation( WAnimation(WAnimation::SlideInFromTop, WAnimation::EaseOut, 100) );
    
    WContainerWidget *holder = new WContainerWidget();
    m_chartPanel->setCentralWidget( holder );
    WGridLayout *layout = new WGridLayout();
    holder->setLayout( layout );
    layout->setContentsMargins( 0, 0, 0, 0 );
    
    m_rel_eff_chart = new RelEffChart();
    m_rel_eff_chart->setHeight( 250 );
    m_rel_eff_chart->setYAxisTitle( WString::tr("psw-rel-eff-plot-y-axis-title").toUTF8() );
    layout->addWidget( m_rel_eff_chart, 0, 0 );
  }//void setupRelEffChart()
 
  
  void refreshRelEffChart()
  {
    if( !m_chartPanel || !m_rel_eff_chart )
      return;
    
    set<const SandiaDecay::Nuclide *> new_peak_nucs;
    for( const auto &old_new : m_old_to_new_peaks )
    {
      if( eligible_for_rel_eff_chart(old_new.second) 
         && (!old_new.first || (old_new.first->decayParticle() != old_new.second->decayParticle()) ))
      {
        new_peak_nucs.insert( old_new.second->parentNuclide() );
      }
    }//for( const auto &old_new : m_old_to_new_peaks )
    
    vector<RelActCalcManual::GenericPeakInfo> peaks;
    map<string,pair<double,std::string>> relActsColors;  //maps source --> {rel-act,color}
    map<string,vector<shared_ptr<PeakDef>>> label_to_used_new_peaks;
    for( size_t i = 0; i < m_old_to_new_peaks.size(); ++i )
    {
      shared_ptr<PeakDef> p = m_old_to_new_peaks[i].second;
      
      // If we are doing a peak-search, or assigning nuclides from reference lines, the user may
      //  have already fit some peaks for the showing reference nuclide, in which case `p` will
      //  likely be nullptr (because a new peak wouldnt have been fit for or modified), so in this
      //  case, we'll use the old peak, for plotting the Rel. Eff. line, if the user has the
      //  "Show previous peaks too" checkmark checked.
      //  However, if we are doing a search from previous spectrum, then the "old" peak is
      //  from the other spectrum (slightly normalized to this spectrum), so we definitely
      //  dont want to include the possibility of using it.
      if( !p
         && (m_reason != PeakSelectorWindowReason::PeaksFromPreviousSpectrum)
         && m_showAllPeaks && m_showAllPeaks->isChecked() )
      {
        p = m_old_to_new_peaks[i].first;
      }
      
      if( !m_keep_peak_cbs[i] && (!m_showAllPeaks || !m_showAllPeaks->isChecked()) )
        continue;
      
      if( !eligible_for_rel_eff_chart(p) || !new_peak_nucs.count(p->parentNuclide()) )
        continue;
      
      if( m_keep_peak_cbs[i] && !m_keep_peak_cbs[i]->isChecked() )
        continue;
      
      const double peakSigma = p->gausPeak() ? p->sigma() : 0.25*p->roiWidth();
      
      const SandiaDecay::Nuclide * const nuc = p->parentNuclide();
      const double photopeakEnergy = p->decayParticle()->energy;
      
      const string label = p->parentNuclide()->symbol;
      
      label_to_used_new_peaks[label].push_back( p );
      
      double intensity = 0.0;
      bool gotInfoFromRefLines = false;
      
      for( const shared_ptr<const ReferenceLineInfo> &refline : m_displayed )
      {
        if( refline->m_nuclide == p->parentNuclide() )
        {
          if( !refline->m_input.m_color.isDefault()
             && !relActsColors.count(label) )
          {
            relActsColors[label].first = 1.0; //relative-activity, we'll modify later.
            relActsColors[label].second = refline->m_input.m_color.cssText();
          }
          
          
          for( const ReferenceLineInfo::RefLine &line : refline->m_ref_lines )
          {
            // We only want x-rays and gammas to make rel-eff chart from
            switch( line.m_particle_type )
            {
              case ReferenceLineInfo::RefLine::Particle::Alpha:
              case ReferenceLineInfo::RefLine::Particle::Beta:
                continue;
                break;
                
              case ReferenceLineInfo::RefLine::Particle::Gamma:
              case ReferenceLineInfo::RefLine::Particle::Xray:
                break;
            }//switch( line.m_particle_type )
            
            // I guess we *could* be interested in Rel. Eff. line of sum or escape peaks,
            //  but for the moment, we'll just use normal gammas and x-rays
            switch( line.m_source_type )
            {
              case ReferenceLineInfo::RefLine::RefGammaType::Normal:
              case ReferenceLineInfo::RefLine::RefGammaType::Annihilation:
                break;
              
              case ReferenceLineInfo::RefLine::RefGammaType::SingleEscape:
              case ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape:
              case ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak:
              case ReferenceLineInfo::RefLine::RefGammaType::SumGammaPeak:
                continue;
                break;
            }//switch( line.m_source_type )
            
            
            const double refenergy = line.m_energy;
            const double br = line.m_decay_intensity;
            if( fabs(refenergy - photopeakEnergy) < 1.25*peakSigma )
              intensity += br;
          }//for( const ReferenceLineInfo::RefLine &line : refline->m_ref_lines )
          
          gotInfoFromRefLines = true;
          break;
        }//if( refline->m_nuclide == p->parentNuclide() )
        
        if( gotInfoFromRefLines )
          break;
      }//for( const auto &refline : m_displayed )
      
      if( !gotInfoFromRefLines )
      {
        const double dummy_activity = 0.01*SandiaDecay::curie;
        
        const double age = PeakDef::defaultDecayTime( nuc );
        SandiaDecay::NuclideMixture mix;
        mix.addAgedNuclideByActivity( nuc, dummy_activity, age );
      
        for( const auto gamma : mix.photons(0, SandiaDecay::NuclideMixture::OrderByEnergy ) )
        {
          if( fabs(gamma.energy-photopeakEnergy) < 1.25*peakSigma )
            intensity += gamma.numPerSecond;
        }
      
        if( intensity <= 0.0 )
          continue;
         
        intensity /= dummy_activity;
      }//if( !gotInfoFromRefLines )
      
      if( intensity > 0.0 )
      {
        const double peakMean = p->gausPeak() ? p->mean() : 0.5*(p->upperX() + p->lowerX());
        RelActCalcManual::GenericPeakInfo peak;
        peak.m_mean = peak.m_energy = peakMean;
        peak.m_fwhm = 2.35482 * peakSigma;
        peak.m_counts = p->peakArea();
        peak.m_counts_uncert = (p->peakAreaUncert() > 0.0) ? p->peakAreaUncert() : 0.0;
        peak.m_base_rel_eff_uncert = 0.0;
        peak.m_source_gammas.emplace_back( intensity, label );
        peaks.push_back( peak );
      }//if( intensity > 0.0 )
    }//for( loop over peaks )
    
    
    // We want to make the maximum relative efficiency point 1.0, for each source
    size_t max_num_peaks = 0;
    for( auto &src_to_peaks : label_to_used_new_peaks )
    {
      const string &src = src_to_peaks.first;
      
      auto act_color_pos = relActsColors.find(src);
      if( act_color_pos == end(relActsColors) )
        act_color_pos = relActsColors.insert( make_pair(src,pair<double,string>{1.0,""}) ).first;
      assert( act_color_pos != end(relActsColors) );
      
      // If we arent showing reference, we wont have colors, so lets check if any of the peaks
      //  are colored for this source, and if so use the first peaks color.
      if( act_color_pos->second.second.empty() )
      {
        const vector<shared_ptr<PeakDef>> &peak_defs = src_to_peaks.second;
        for( const shared_ptr<PeakDef> &p : peak_defs )
        {
          if( !p->lineColor().isDefault() )
          {
            act_color_pos->second.second = p->lineColor().cssText();
            break;
          }
        }//for( const auto &p : peak_defs )
      }//if( act_color_pos->second.second.empty() )
      
      double max_eff = 0.0;
      size_t num_src_peaks = 0.0;
      for( const RelActCalcManual::GenericPeakInfo &peak : peaks )
      {
        for( const auto &line : peak.m_source_gammas )
        {
          if( line.m_isotope == src )
          {
            num_src_peaks += 1;
            max_eff = std::max( max_eff, peak.m_counts / line.m_yield );
          }
        }
      }//for( const RelActCalcManual::GenericPeakInfo &peak : peaks )
      
      max_num_peaks = std::max( max_num_peaks, num_src_peaks );
      assert( max_eff > 0.0 );
      
      double &rel_act = act_color_pos->second.first;
      rel_act = (max_eff > 0.0) ? max_eff : 1.0;
    }//for( auto &src_to_peaks : label_to_used_new_peaks )
    
    m_chartPanel->setHidden( (max_num_peaks < 2) );
    
    string relEffEqn = "";
    WString title_chi2_info;
    m_rel_eff_chart->setData( peaks, relActsColors, relEffEqn, title_chi2_info, "" );
  }//void refreshRelEffChart()
#endif // USE_REL_ACT_TOOL
  
  
  void showAllPeaksCbChanged()
  {
    if( !m_showAllPeaks )
      return;
    
    const bool showAll = m_showAllPeaks->isChecked();
    if( showAll )
    {
      for( size_t i = 0; i < m_old_to_new_peaks.size(); ++i )
        m_table->rowAt( static_cast<int>(i+1) )->setHidden( false );
    }else
    {
      for( size_t i = 0; i < m_old_to_new_peaks.size(); ++i )
        m_table->rowAt( static_cast<int>(i+1) )->setHidden( !m_peak_updated[i] );
    }
    
#if( USE_REL_ACT_TOOL )
    refreshRelEffChart();
#endif
  }//void showAllPeaksCbChanged()
  
  
  void updatePreviewPlot( size_t i )
  {
    assert( i < m_old_to_new_peaks.size() );
    if( i >= m_old_to_new_peaks.size() || m_previewChartColumn < 0 )
      return;
    
    auto finalpeak = m_old_to_new_peaks[i].second;
    auto originalpeak = m_old_to_new_peaks[i].first;
    auto displayPeak = finalpeak ? finalpeak : originalpeak;
    
    if( !displayPeak )
      return;
    
    std::shared_ptr<const ColorTheme> theme = m_viewer->getColorTheme();
    
    auto peaks_for_plotting = make_shared< std::deque<std::shared_ptr<const PeakDef> > >();
    for( size_t peak_index = 0; peak_index < m_old_to_new_peaks.size(); ++peak_index )
    {
      const auto &p = m_old_to_new_peaks[peak_index];
      
      if( !p.second )
        continue;
      
      auto new_peak = make_shared<PeakDef>(*p.second);
      
      if( peak_index != i )
      {
        WColor color = new_peak->lineColor();
        if( color.isDefault() && theme )
          color = theme->defaultPeakLine;
        
        color.setRgb( color.red(), color.green(), color.blue(), 50 );
        new_peak->setLineColor( color );
      }//if( peak_index != i )
      
      peaks_for_plotting->push_back( new_peak );
    }//for( const auto &p : m_old_to_new_peaks )
    
    const double lowerx = displayPeak->mean() - 0.75*displayPeak->roiWidth();
    const double upperx = displayPeak->mean() + 0.75*displayPeak->roiWidth();
    
    //TODO: Find 'displayPeak' in 'peaks_for_plotting' by energy, and check
    //  if it shares a ROI with other peaks, and if so, change the other
    //  peak colors, and maybe the displayed width.
    
    
    /*
     // The SVG renderer for Wt 3.7.1 appears to glitch-up sometimes, so instead of using SVG,
     //  we'll display a full SpectrumChart.
    std::shared_ptr<WSvgImage> svg = PeakSearchGuiUtils::renderChartToSvg( m_data, peaks_for_plotting, m_displayed, lowerx, upperx, 225, 125, theme, true );
    if( svg )
    {
      stringstream strm;
      svg->write( strm );
      WTableCell *preview = m_table->elementAt( static_cast<int>(i+1), m_previewChartColumn );
      preview->clear();
      if( !preview->hasStyleClass("PeakPreviewCell") )
        preview->addStyleClass( "PeakPreviewCell" );
      new WText( strm.str(), Wt::XHTMLUnsafeText, preview );
    }
     */
    
    
    SpectrumChart *chart = PeakSearchGuiUtils::createFixedSpectrumDisplay( m_data,
                                                    peaks_for_plotting, m_displayed,
                                                    lowerx, upperx,
                                                    225, 125, theme );
    if( chart )
    {
#if( !DYNAMICALLY_ADJUST_LEFT_CHART_PADDING )
      chart->setPlotAreaPadding( 14, Wt::Left );
#endif
      chart->setPlotAreaPadding( 22, Wt::Bottom );
      chart->setPlotAreaPadding( 0, Wt::Right );
      chart->setPlotAreaPadding( 0, Wt::Top );
      
      WFont labelFont( WFont::Default );
      labelFont.setSize(WFont::Size::XXSmall);
      chart->axis(Chart::XAxis).setLabelFont( labelFont );
      
      // The y-axis label font doesnt seem to be respected (and instead a 10pt font is always used)
      //  Rendering the axis is done in SpectrumChart.cpp, and it isnt clear what the problem is,
      //  so we'll just leave off the y-axis labels for now.  Not a great solution, but its time
      //  to move on to more effective ways to spend time.
      chart->axis(Chart::YAxis).setLabelFont( labelFont ); //doesnt seem to be obeyed
      chart->axis(Chart::YAxis).setLabelFormat( " " );
      
      WTableCell *preview = m_table->elementAt( static_cast<int>(i+1), m_previewChartColumn );
      preview->clear();
      if( !preview->hasStyleClass("PeakPreviewCell") )
        preview->addStyleClass( "PeakPreviewCell" );
      preview->addWidget( chart );
    }//if( previewChart )
  }//void updatePreviewPlot( size_t i )
  
  
  void keepPeakChanged()
  {
    bool anyNonAssignedPeaks = false;
    for( size_t i = 0; !anyNonAssignedPeaks && i < m_old_to_new_peaks.size(); ++i )
    {
      if( m_table->rowAt( static_cast<int>(i+1) )->isHidden() )
        continue;
      
      if( m_old_to_new_peaks[i].second && m_keep_peak_cbs[i] && m_keep_peak_cbs[i]->isChecked() )
        anyNonAssignedPeaks = !m_old_to_new_peaks[i].second->hasSourceGammaAssigned();
    }//for( size_t i = 0; !anyNonAssignedPeaks && i < m_old_to_new_peaks.size(); ++i )
    
    for( size_t i = 0; i < m_old_to_new_peaks.size(); ++i )
    {
      if( !m_keep_peak_cbs[i] )
        continue;
      
      const bool enable = m_keep_peak_cbs[i]->isChecked();
        
      if( m_nuc_select_combos[i] )
        m_nuc_select_combos[i]->setDisabled( !enable );
        
      if( m_dont_change_nuc_cbs[i] )
        m_dont_change_nuc_cbs[i]->setDisabled( !enable );
    }//for( size_t i = 0; i < m_old_to_new_peaks.size(); ++i )
    
    if( m_keepRefLinePeaksOnly )
      m_keepRefLinePeaksOnly->setChecked( !anyNonAssignedPeaks );
    
#if( USE_REL_ACT_TOOL )
    refreshRelEffChart();
#endif
  }//void keepPeakChanged()
  
  
  void keepOnlyRefLinesCbChanged()
  {
    if( !m_keepRefLinePeaksOnly )
      return;
    
    const bool keepOnlyAssigned = m_keepRefLinePeaksOnly->isChecked();
    
    for( size_t i = 0; i < m_old_to_new_peaks.size(); ++i )
    {
      if( m_table->rowAt( static_cast<int>(i+1) )->isHidden()
          || !m_keep_peak_cbs[i] )
        continue;
    
      if( !keepOnlyAssigned )
      {
        m_keep_peak_cbs[i]->setChecked(true);
      }else if( m_old_to_new_peaks[i].second )
      {
        m_keep_peak_cbs[i]->setChecked( m_old_to_new_peaks[i].second->hasSourceGammaAssigned() );
      }
    }//for( size_t i = 0; i < m_old_to_new_peaks.size(); ++i )
    
    keepPeakChanged();
#if( USE_REL_ACT_TOOL )
    refreshRelEffChart();
#endif
  }//void keepOnlyRefLinesCbChanged()
  
  
  
  std::string nuclideDesc( const SandiaDecay::Nuclide *nuc,
                           const SandiaDecay::RadParticle *particle,
                           const PeakDef::SourceGammaType type )
  {
    if( !nuc || !particle )
      return "---";
    
    const ReferenceLineInfo *refline = nullptr;
    
    const char *gamtype = "";
    switch( type )
    {
      case PeakDef::NormalGamma:
      case PeakDef::XrayGamma:
        break;
      case PeakDef::AnnihilationGamma: gamtype = " Annih."; break;
      case PeakDef::DoubleEscapeGamma: gamtype = " D.E.";   break;
      case PeakDef::SingleEscapeGamma: gamtype = " S.E.";   break;
    }//switch( displayPeak->sourceGammaType() )
    
    
    for( const auto &d : m_displayed )
    {
      if( nuc == d->m_nuclide )
      {
        refline = d.get();
      }else
      {
        
        for( const ReferenceLineInfo::RefLine &b : d->m_ref_lines )
        {
          if( nuc == b.m_parent_nuclide )
          {
            refline = d.get();
            break;
          }
        }//for( const ReferenceLineInfo::RefLine &b : d->m_ref_lines )
      }//if ( nuc matches ) / else
      
      if( refline )
        break;
    }//for( const auto &d : m_displayed )
    
    char buffer[256] = { '\0' };
    if( refline )
    {
      double intensity = 0.0;
      for( const ReferenceLineInfo::RefLine &line : refline->m_ref_lines )
      {
        if( fabs(line.m_energy - particle->energy) >= 0.01 )
          continue;
        
        // We are only interested i ngammas and x-rays
        switch( line.m_particle_type )
        {
          case ReferenceLineInfo::RefLine::Particle::Alpha:
          case ReferenceLineInfo::RefLine::Particle::Beta:
            continue;
            break;
            
          case ReferenceLineInfo::RefLine::Particle::Gamma:
          case ReferenceLineInfo::RefLine::Particle::Xray:
            break;
        }//switch( line.m_particle_type )
        
        // I guess we *could* be interested in Rel. Eff. line of sum or escape peaks,
        //  but for the moment, we'll just use normal gammas and x-rays
        switch( line.m_source_type )
        {
          case ReferenceLineInfo::RefLine::RefGammaType::Normal:
          case ReferenceLineInfo::RefLine::RefGammaType::Annihilation:
            break;
            
          case ReferenceLineInfo::RefLine::RefGammaType::SingleEscape:
          case ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape:
          case ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak:
          case ReferenceLineInfo::RefLine::RefGammaType::SumGammaPeak:
            continue;
            break;
        }//switch( line.m_source_type )
        
        intensity += line.m_decay_intensity;
      }//for( loop over ref line particles )
      
      snprintf( buffer, sizeof(buffer), "%s %.1f keV%s, I=%.2g%%",
               nuc->symbol.c_str(), particle->energy, gamtype, 100.0*intensity );
    }else
    {
      const double activity = 1.0E-3*SandiaDecay::curie;
      const double age = PeakDef::defaultDecayTime( nuc );
      
      SandiaDecay::NuclideMixture mix;
      mix.addAgedNuclideByActivity( nuc, activity, age );
      const auto gammas = mix.gammas(0.0, SandiaDecay::NuclideMixture::OrderByEnergy, true );
      double intensity = 0.0;
      for( const auto &g : gammas )
      {
        if( fabs(g.energy - particle->energy) < 0.01 )
          intensity += g.numPerSecond;
      }
      intensity /= activity;
      
      snprintf( buffer, sizeof(buffer), "%s %.1f keV%s, I=%.2g%%",
               nuc->symbol.c_str(), particle->energy, gamtype, 100.0*intensity );
    }//if( source is a reference line ) / else
    
    
    return buffer;
  }//std::string nuclideDesc(...)
  
  
  std::string xrayDesc( const SandiaDecay::Element *el, const float energy )
  {
    if( !el )
      return "---";
    
    const SandiaDecay::EnergyIntensityPair *xrayIntensity = PeakDef::findNearestXray(el, energy);
    
    char buffer[256] = { '\0' };
    if( xrayIntensity && fabs(xrayIntensity->energy - energy) < 0.1 )
    {
      snprintf( buffer, sizeof(buffer), "%s %.1f keV xray, I=%.2g%%)",
               el->symbol.c_str(), energy, 100.0*xrayIntensity->intensity );
    }else
    {
      snprintf( buffer, sizeof(buffer), "%s %.1f keV xray", el->symbol.c_str(), energy );
    }
    
    return buffer;
  }//std::string xrayDesc( const SandiaDecay::Element *el )
  
  
  std::string reactionDesc( const ReactionGamma::Reaction *rctn,
                            const float energy,
                            const PeakDef::SourceGammaType type )
  {
    if( !rctn )
      return "---";
    
    const char *gamtype = "";
    switch( type )
    {
      case PeakDef::NormalGamma: case PeakDef::XrayGamma:    break;
      case PeakDef::AnnihilationGamma: gamtype = " Annih."; break;
      case PeakDef::DoubleEscapeGamma: gamtype = " D.E.";   break;
      case PeakDef::SingleEscapeGamma: gamtype = " S.E.";   break;
    }//switch( displayPeak->sourceGammaType() )
    
    
    char buffer[256] = { '\0' };
    
    double intensity = 0.0;
    for( const auto &g : rctn->gammas )
      if( fabs(g.energy - energy) < 0.1 )
        intensity += g.abundance;
    
    snprintf( buffer, sizeof(buffer), "%s %.1f keV%s, I=%.2g%%",
              rctn->name().c_str(), energy, gamtype, 100.0*intensity );
    
    return buffer;
  }//std::string reactionDesc(...)
  
  
  std::string makeSourceDesciption( const shared_ptr<PeakDef> &p )
  {
    
    if( !p )
      return "---";
    
    if( p->parentNuclide() && p->decayParticle() )
    {
      return nuclideDesc( p->parentNuclide(), p->decayParticle(), p->sourceGammaType() );
    }else if( p->xrayElement() )
    {
      return xrayDesc( p->xrayElement(), p->xrayEnergy() );
    }else if( p->reaction() )
    {
      return reactionDesc( p->reaction(), p->reactionEnergy(), p->sourceGammaType() );
    }
    
    return "---";
  }//std::string makeSourceDesciption( const shared_ptr<PeakDef> &p )
  
  
  void nucSelectChanged( const size_t i )
  {
    assert( i < m_dont_change_nuc_cbs.size() );
    assert( m_nuc_select_combos.size() == m_dont_change_nuc_cbs.size() );
    
    if( (i >= m_nuc_select_combos.size()) || (i >= m_dont_change_nuc_cbs.size()) )  //shouldnt ever happen, but JIC
      return;
    
    //shared_ptr<PeakDef> &oldpeak = m_old_to_new_peaks[i].first;
    shared_ptr<PeakDef> &newpeak = m_old_to_new_peaks[i].second;
    assert( newpeak );
    if( !newpeak )  //shouldnt ever happen!
      return;
    
    WComboBox *combo = m_nuc_select_combos[i];
    assert( combo );
    if( !combo ) //shouldnt ever happen!
      return;
    
    WCheckBox *cb = m_dont_change_nuc_cbs[i];
    const bool keepOriginal = (cb && cb->isChecked());
    
    if( keepOriginal )
    {
      
      
      switch( m_reason )
      {
        case PeakSelectorWindowReason::PeakSearch:
        case PeakSelectorWindowReason::PeaksFromPreviousSpectrum:
        case PeakSelectorWindowReason::PeaksFromCsvFile:
          newpeak->clearSources();
          break;
          
        case PeakSelectorWindowReason::NuclideId:
        {
          shared_ptr<PeakDef> &oldpeak = m_old_to_new_peaks[i].first;
          if( oldpeak )
            newpeak->inheritUserSelectedOptions( *oldpeak, false );
          break;
        }
      }//switch( m_reason )
      
      
      updatePreviewPlot(i);
      keepPeakChanged();
#if( USE_REL_ACT_TOOL )
      refreshRelEffChart();
#endif
      return;
    }//if( keepOriginal )
    
    
    const int index = combo->currentIndex();
    
    if( index <= 0 )
    {
      newpeak->clearSources();
      
      if( m_viewer->colorPeaksBasedOnReferenceLines() )
      {
        for( const auto &l : m_displayed )
        {
          if( l->m_input.m_color == newpeak->lineColor() )
            newpeak->setLineColor( WColor() );
        }//
      }//if( m_viewer->colorPeaksBasedOnReferenceLines() )
      
      updatePreviewPlot(i);
      keepPeakChanged();
#if( USE_REL_ACT_TOOL )
      refreshRelEffChart();
#endif
      return;
    }//if( no source )
    
    string newnucstr = combo->currentText().toUTF8();
    const size_t intens_pos = newnucstr.find( ", I=" );
    if( intens_pos != string::npos )
      newnucstr = newnucstr.substr( 0, intens_pos );

    const auto result = PeakModel::setNuclideXrayReaction( *newpeak, newnucstr, 0.0 );
    
    switch( result )
    {
      case PeakModel::SetGammaSource::NoSourceChange:
        //Shouldnt ever happen, but JIC...
        passMessage( "Trouble making change to nuclide - not applying, sorry!<br />"
                    "Please report error, including selected nuclides to interspec@sandia.gov", 2 );
        populateNuclideSelects();
        break;
        
      case PeakModel::SetGammaSource::FailedSourceChange:
        //Shouldnt ever happen, but JIC...
        passMessage( "Failed changing to nuclide - not applying, sorry!<br />"
                     "Please report error, including selected nuclides to interspec@sandia.gov", 2 );
        populateNuclideSelects();
        break;
        
      case PeakModel::SetGammaSource::SourceChange:
      case PeakModel::SetGammaSource::SourceAndUseChanged:
      {
        //Update color of peak
        if( m_viewer->colorPeaksBasedOnReferenceLines() )
        {
          const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
          
          for( const auto &l : m_displayed )
          {
            //A quick labmda to test if a peak nuclide matches a BackgroundLine
            //  nuclide.
            auto isBackgroundNuc = [&l,&newpeak,db]() -> bool {
              if( !newpeak->parentNuclide() || (l->m_source_type != ReferenceLineInfo::SourceType::Background) )
                return false;
              
              for( const ReferenceLineInfo::RefLine &b : l->m_ref_lines )
              {
                if( b.m_parent_nuclide && (b.m_parent_nuclide == newpeak->parentNuclide()) )
                  return true;
              }
              return false;
            };//isBackgroundNuc lambda
            
            
            if( m_viewer->colorPeaksBasedOnReferenceLines()
               && ((l->m_nuclide && (l->m_nuclide==newpeak->parentNuclide()))
               || (l->m_element && (l->m_element==newpeak->xrayElement()))
               || (newpeak->reaction() && l->m_reactions.count(newpeak->reaction()))
               || isBackgroundNuc()) )
            {
              newpeak->setLineColor( l->m_input.m_color );
            }
          }//for( const auto &l : m_displayed )
        }//if( we should set peak color based on reference lines )
        break;
      }//case: source change successful
    }//switch( result )
    
    updatePreviewPlot(i);
    keepPeakChanged();
#if( USE_REL_ACT_TOOL )
    refreshRelEffChart();
#endif
  }//void nucSelectChanged( const size_t i )
  
  
  void dontChangeNucCbChanged( const size_t i )
  {
    assert( i < m_dont_change_nuc_cbs.size() );
    assert( m_nuc_select_combos.size() == m_dont_change_nuc_cbs.size() );
    if( (i >= m_nuc_select_combos.size()) || (i >= m_dont_change_nuc_cbs.size()) )  //shouldnt ever happen, but JIC
      return;
    
    WComboBox *select = m_nuc_select_combos[i];
    WCheckBox *cb = m_dont_change_nuc_cbs[i];
    assert( select && cb );
    
    //We should only get here if select and cb are valid
    if( select && cb )
    {
      select->setEnabled( !cb->isChecked() );
      
      if( cb->isChecked() 
         && (i < m_old_to_new_peaks.size())
         && m_old_to_new_peaks[i].first 
         && !m_old_to_new_peaks[i].first->hasSourceGammaAssigned() )
      {
        select->setCurrentIndex( 0 );
      }//
    }//if( select && cb )
    
    nucSelectChanged( i );
  }//void dontChangeNucCbChanged( const size_t i )
  
  
  void populateNuclideSelects()
  {
    assert( m_old_to_new_peaks.size() == m_nuc_select_combos.size() );
    
    for( size_t i = 0; i < m_old_to_new_peaks.size(); ++i )
    {
      if( !m_nuc_select_combos[i] )
        continue;
      
      m_nuc_select_combos[i]->clear();
      m_nuc_select_combos[i]->addItem( WString::tr("psd-no-source") );
      
      shared_ptr<PeakDef> &p = m_old_to_new_peaks[i].second;
      if( !p )
        continue;
      shared_ptr<PeakDef> &oldpeak = m_old_to_new_peaks[i].first;
      
      const double peakmean = p->mean();
      double energy = peakmean;
      
      if( p->xrayElement() || p->reaction() )
      {
        energy = p->gammaParticleEnergy();
      }else if( p->parentNuclide() )
      {
        if( p->decayParticle() )
          energy = p->decayParticle()->energy;  // Use this so the 511 or 1022 keV wont be subtracted off for S.E. or D.E.
        else
          energy = p->gammaParticleEnergy();    //Annihilation gammas make it here
      }
      
      
      //This string will store the description of what is currently assigned to
      //  the peak so we can make it the current selection
      //  I'm a little un-easy about how this is determined... but seems to work
      string currentstr;
      
      //If this was a peak-search operation leading to this dialog, and the old
      //  peak did not have any source identified, then dont identify a source
      //  for the new peak either.  If the old peak did have a source, then that
      //  would not have been changed by the ID part of the search anyway.
      switch( m_reason )
      {
        case PeakSelectorWindowReason::PeakSearch:
          if( oldpeak && !oldpeak->hasSourceGammaAssigned() )
            p->clearSources();
          break;
          
        case PeakSelectorWindowReason::NuclideId:
          break;
          
        case PeakSelectorWindowReason::PeaksFromPreviousSpectrum:
          break;
          
        case PeakSelectorWindowReason::PeaksFromCsvFile:
          break;
      }//switch( m_reason )

      
      //We'll sort the select by distance away from current nuclide, or peak mean
      //  - will probably miss some weird edgecase, but what ever for the moment.
      map<double,set<string>> distToDescs;
      for( const auto &ref : m_displayed )
      {
        // Make sure we can assign the peak to the source type.
        switch( ref->m_source_type )
        {
          case ReferenceLineInfo::SourceType::Nuclide:
          case ReferenceLineInfo::SourceType::FluorescenceXray:
          case ReferenceLineInfo::SourceType::Reaction:
          case ReferenceLineInfo::SourceType::Background:
          case ReferenceLineInfo::SourceType::NuclideMixture:
          case ReferenceLineInfo::SourceType::OneOffSrcLines:
          case ReferenceLineInfo::SourceType::FissionRefLines:
            break;
          
          case ReferenceLineInfo::SourceType::CustomEnergy:
          case ReferenceLineInfo::SourceType::None:
            continue;
            break;
        }//switch( ref->m_source_type )
        
        for( const ReferenceLineInfo::RefLine &line : ref->m_ref_lines )
        {
          if( (line.m_particle_type != ReferenceLineInfo::RefLine::Particle::Gamma)
              && (line.m_particle_type != ReferenceLineInfo::RefLine::Particle::Xray) )
            continue;
          
          const double refenergy = line.m_energy;
          const double intensity = line.m_decay_intensity; //Before DRF and Shielding
          
          string label = ref->m_input.m_input_txt;
          if( line.m_parent_nuclide )
          {
            label = line.m_parent_nuclide->symbol;
          }else if( line.m_element )
          {
            label = line.m_element->symbol;
          }else if( line.m_reaction )
          {
            //We wont set label to reaction, because the user may have input an element, but this
            //  particular line might be an isotope of the element.
            if( SpecUtils::icontains(label, "background") )
              label = line.m_reaction->name();
          }
          
          char buffer[128];
          snprintf( buffer, sizeof(buffer), "%s %.1f keV, I=%.2g%%",
                    label.c_str(), refenergy, 100.0*intensity );
          
          const double diff_from_mean = fabs( refenergy - peakmean );
          const double diff_from_assigned = fabs( refenergy - energy );
          
          // If we are assigning from a nuclide, we should be pretty close, but if OneOffSrcLines,
          //  we might be off by even a
          const double allowed_diff = (ref->m_source_type == ReferenceLineInfo::SourceType::OneOffSrcLines) ? 0.05 : 0.001;
          
          if( (diff_from_assigned < allowed_diff) //May need to decrease
             && ((p->parentNuclide() && (line.m_parent_nuclide==p->parentNuclide()))
                 || (p->xrayElement() && (line.m_element==p->xrayElement()))
                 || (p->reaction() && (line.m_reaction ==p->reaction()))
             ))
          {
#if( PERFORM_DEVELOPER_CHECKS )
            if( currentstr.size() )  //shouldnt ever happen, right?
              log_developer_error( __func__, ("PeakSelectorWindow::populateNuclideSelects(): Found mutliple matches of select strings for a nuclide: '" + string(buffer) + "'").c_str() );
#endif
            currentstr = buffer;
          }else if( (ref->m_source_type == ReferenceLineInfo::SourceType::Background)
                   && (diff_from_assigned < 0.1) )
          {
            currentstr = buffer;
          }
          
          //Use the minimum of distance between the reference line and either
          //  the peak mean, or currently assigned source, to sort select
          //  contents by.
          const double diff = min( diff_from_mean, diff_from_assigned );
          
          distToDescs[diff].insert( buffer );
        }//for( size_t i = 0; i < ref.energies.size(); ++i )
      }//for( const auto &ref : m_displayed )
      
      int currentitem = -1;
      size_t nitems = 0;
      const double fwhm = p->fwhm();
      for( auto iter = begin(distToDescs); iter != end(distToDescs); ++iter )
      {
        const double energydiff = iter->first;
        //Put at most 30 items into the select, or if we have 15 items, keep
        //  from wandering >5 FWHM from both the current assigned value and peak
        //  mean
        if( (nitems > 30) || (nitems > 15 && energydiff > 5*fwhm) )
          break;
        
        ++nitems;
        for( const auto &s : iter->second )
        {
          m_nuc_select_combos[i]->addItem( s );
          if( s == currentstr )
            currentitem = m_nuc_select_combos[i]->count() - 1;
        }
      }//for( auto iter = begin(distToDescs); iter != end(distToDescs); ++iter )
      
      m_nuc_select_combos[i]->setCurrentIndex( std::max(0, currentitem) );
    }//for( size_t i = 0; i < m_old_to_new_peaks.size(); ++i )
  }//void populateNuclideSelects()
  
  
  void cancelOperation()
  {
    m_cancelOperation = true;
    
    hide();
  }//void cancelOperation()
  
  
  void doFinish()
  {
    assert( m_old_to_new_peaks.size() == m_keep_peak_cbs.size() );
    assert( m_old_to_new_peaks.size() == m_nuc_select_combos.size() );
    
    for( size_t i = 0; i < m_old_to_new_peaks.size(); ++i )
    {
      if( m_dont_change_nuc_cbs[i] && m_dont_change_nuc_cbs[i]->isChecked() )
      {
        if( m_old_to_new_peaks[i].first )
        {
          switch( m_reason )
          {
            case PeakSelectorWindowReason::NuclideId:
            case PeakSelectorWindowReason::PeakSearch:
              *m_old_to_new_peaks[i].second = *m_old_to_new_peaks[i].first;
            break;
            
            case PeakSelectorWindowReason::PeaksFromPreviousSpectrum:
            case PeakSelectorWindowReason::PeaksFromCsvFile:
              m_old_to_new_peaks[i].second->inheritUserSelectedOptions( *m_old_to_new_peaks[i].first, false );
            break;
          }//switch( m_reason )
        }else
        {
          m_old_to_new_peaks[i].second->clearSources();
          m_old_to_new_peaks[i].second->setLineColor( WColor() );
        }
      }else if( m_nuc_select_combos[i] && (m_nuc_select_combos[i]->currentIndex()<=0)
          && m_old_to_new_peaks[i].second && !m_table->rowAt(static_cast<int>(i+1))->isHidden() )
      {
        m_old_to_new_peaks[i].second->clearSources();
        m_old_to_new_peaks[i].second->setLineColor( WColor() );
      }
      
      if( m_keep_peak_cbs[i] && !m_keep_peak_cbs[i]->isChecked() )
      {
        m_old_to_new_peaks[i].second.reset();
        if( m_reason == PeakSelectorWindowReason::PeaksFromCsvFile )
          m_old_to_new_peaks[i].first.reset();
      }
    }//for( size_t i = 0; i < peaks.size(); ++i )
    
    auto peakModel = m_viewer->peakModel();
    if( peakModel )
    {
      vector<PeakDef> final_peaks;
      for( const auto &i : m_old_to_new_peaks )
      {
        switch( m_reason )
        {
          case PeakSelectorWindowReason::NuclideId:
          case PeakSelectorWindowReason::PeakSearch:
            if( i.second && !m_cancelOperation )
              final_peaks.push_back( *i.second );
            else if( i.first )
              final_peaks.push_back( *i.first );
          break;
          
          case PeakSelectorWindowReason::PeaksFromPreviousSpectrum:
            if( i.second && !m_cancelOperation )
              final_peaks.push_back( *i.second );
          break;
            
          case PeakSelectorWindowReason::PeaksFromCsvFile:
            if( i.first )
              final_peaks.push_back( *i.first );
            else if( i.second && !m_cancelOperation )
              final_peaks.push_back( *i.second );
          break;
        }//switch( m_reason )
      }//for( const auto &i : m_old_to_new_peaks )
      
      peakModel->setPeaks( final_peaks );
    }//if( peakModel )
    
    AuxWindow::deleteAuxWindow( this );
  }//void doFinish()
  
};//class PeakSelectorWindow


  
  
  
//setPeaksFromSearch(): sets the peak model peaks to those passed in; should
//  be called from main event loop thread.  'guiupdater' is intended to be a
//  function that will close the waiting dialog, if it still exists.
//  If 'data' does not match the currently displayed histogram (by pointer
//  comparison), then peaks will not be updated.
void set_peaks_from_search( InterSpec *viewer,
                            vector<ReferenceLineInfo> displayed,
                            std::shared_ptr<std::vector<std::shared_ptr<const PeakDef> > > peaks,
                            std::shared_ptr<const SpecUtils::Measurement> originaldata,
                            std::vector<PeakDef> originalPeaks,
                            boost::function<void(void)> guiupdater )
{
  if( !viewer )
    return;
  
  viewer->useMessageResourceBundle( "PeakSearchGuiUtils" );
  
  auto currentdata = viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  
  if( (currentdata != originaldata) || !peaks )
  {
    passMessage( WString::tr("psgu-peaks-not-updated-no-change-msg"), WarningWidget::WarningMsgHigh );
    viewer->automatedPeakSearchCompleted();
    return;
  }
  

  vector<std::shared_ptr<const PeakDef> > filtered_peaks;
  for( size_t i = 0; i < peaks->size(); ++i )
  {
    const std::shared_ptr<const PeakDef> p = (*peaks)[i];
    if( fabs(p->mean()) > 0.1 && fabs(p->amplitude()) > 0.1 )
      filtered_peaks.push_back( p );
  }
  
  PeakModel *peakModel = viewer->peakModel();
  
  peakModel->setPeaks( filtered_peaks );
  
  if( !guiupdater.empty() )
    guiupdater();
  
  viewer->automatedPeakSearchCompleted();
  
  passMessage( WString::tr("psgu-peaks-updated-msg"), WarningWidget::WarningMsgInfo );
  
  ReferencePhotopeakDisplay *refLineDisplay = viewer->referenceLinesWidget();
  
  //InterSpec::addPeak(...) is what normally associations nuclides with peaks
  vector<string> showingRefLines;
  if( refLineDisplay )
  {
    const ReferenceLineInfo &showingLines = refLineDisplay->currentlyShowingNuclide();
    if( !showingLines.m_input.m_input_txt.empty() )
      showingRefLines.push_back( showingLines.m_input.m_input_txt );
    
    for( const ReferenceLineInfo &line : refLineDisplay->persistedNuclides() )
    {
      if( !line.m_input.m_input_txt.empty() )
        showingRefLines.push_back( line.m_input.m_input_txt );
    }
  }//if( m_referenceNuclideLines )
  
  if( showingRefLines.empty() && viewer->isMobile() )
  {
    auto undo = [peakModel, originalPeaks](){
      peakModel->setPeaks(originalPeaks);
    };
    auto redo = [peakModel, filtered_peaks](){
      peakModel->setPeaks(filtered_peaks);
    };
    UndoRedoManager *undoManager = viewer->undoRedoManager();
    
    if( undoManager )
      undoManager->addUndoRedoStep( undo, redo, "Automated peak search." );
    
    // Add undo/redo point here
    wApp->triggerUpdate();
    return;
  }
  
  vector<PeakDef> result_peaks;
  for( const auto &p : filtered_peaks )
    result_peaks.push_back( *p );
  
  new PeakSelectorWindow( viewer, PeakSelectorWindowReason::PeakSearch, originalPeaks,
                          originaldata, result_peaks, displayed );
  
  wApp->triggerUpdate();
}//void set_peaks_from_search( const vector<PeakDef> &peaks )
}//namespace



namespace PeakSearchGuiUtils
{
  /* Setups a static, fixed chart for display.
   
   You must either put the result into the widget hierarchy, or delete it.
   */
  SpectrumChart *createFixedSpectrumDisplay( std::shared_ptr<const SpecUtils::Measurement> inmeas,
                            std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef> > > peaks,
                            const std::vector<std::shared_ptr<const ReferenceLineInfo>> &displayed,
                            double lowx, double upperx,
                            const int width, const int height,
                            std::shared_ptr<const ColorTheme> theme )
  {
    if( !inmeas )
      return nullptr;
    
    auto meas = std::make_shared<SpecUtils::Measurement>( *inmeas );
    std::shared_ptr<SpecMeas> specmeas = std::make_shared<SpecMeas>();
    specmeas->add_measurement( meas, true );
    if( peaks )
      specmeas->setPeaks( *peaks, specmeas->sample_numbers() );
    
    std::unique_ptr<SpectrumChart> chart( new SpectrumChart() );
    chart->setWidth( width );
    chart->setHeight( height );
    PeakModel *peakmodel = new PeakModel( chart.get() );
    SpectrumDataModel *dataModel = new SpectrumDataModel( chart.get() );
    
    chart->setModel( dataModel );
    chart->setPeakModel( peakmodel );
    peakmodel->setForeground( meas );
    peakmodel->setPeakFromSpecMeas( specmeas, specmeas->sample_numbers(), SpecUtils::SpectrumType::Foreground );
    
    dataModel->setDataHistogram( meas );
    
    const vector<Chart::WDataSeries> series = dataModel->suggestDataSeries();
    chart->setSeries( series );
    
    //chart.enableLegend( true );
    
    chart->axis(Chart::YAxis).setScale( Chart::LogScale );
    
    for( size_t i = displayed.size(); i > 1; --i )
    {
      chart->setReferncePhotoPeakLines( *displayed[i-1] );
      chart->persistCurrentReferncePhotoPeakLines();
    }
    if( !displayed.empty() )
      chart->setReferncePhotoPeakLines( *displayed[0] );
    
    
    if( meas && (fabs(upperx-lowx) < 0.000001) )
    {
      const size_t nchannel = meas->num_gamma_channels();
      lowx = meas->gamma_channel_lower( 0 );
      upperx = meas->gamma_channel_upper( nchannel - 1 );
    }//if( lowx == upperx )
    
    chart->setXAxisRange( lowx, upperx );
    
    const size_t displayednbin = meas->find_gamma_channel( upperx )
                              - meas->find_gamma_channel( lowx );
    const int plotAreaWidth = width 
                              - chart->plotAreaPadding(Left)
                              - chart->plotAreaPadding(Right);
    const float bins_per_pixel = float(displayednbin) / float(plotAreaWidth);
    const int factor = max( static_cast<int>(ceil(bins_per_pixel)), 1 );
    
    dataModel->setRebinFactor( factor );
    
    chart->setAutoYAxisRange();
    
    if( theme )
    {
      dataModel->setForegroundSpectrumColor( theme->foregroundLine );
      dataModel->setBackgroundSpectrumColor( theme->backgroundLine );
      dataModel->setSecondarySpectrumColor( theme->secondaryLine );
      chart->setSeries( dataModel->suggestDataSeries() );
      
      chart->setDefaultPeakColor( theme->defaultPeakLine );
      
      chart->setAxisLineColor( theme->spectrumAxisLines );
      chart->setChartMarginColor( theme->spectrumChartMargins );
      chart->setChartBackgroundColor( theme->spectrumChartBackground );
      chart->setTextColor( theme->spectrumChartText );
    }//if( theme )
    
  
    return chart.release();
  }//createFixedSpectrumDisplay(...)
  
  
  std::shared_ptr<WSvgImage> renderChartToSvg( std::shared_ptr<const SpecUtils::Measurement> inmeas,
                                              std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef> > > peaks,
                                              const std::vector<std::shared_ptr<const ReferenceLineInfo>> &displayed,
                                              double lowx, double upperx,
                                              const int width, const int height,
                                              std::shared_ptr<const ColorTheme> theme,
                                              const bool compact )
  {
    SpectrumChart *chart = createFixedSpectrumDisplay( inmeas, peaks, displayed, lowx, upperx,
                              width, height, theme );
    if( !chart )
      return nullptr;
    
    unique_ptr<SpectrumChart> chart_holder(chart);
    
    if( compact )
    {
      chart->axis(Chart::XAxis).setTitle("");
      chart->axis(Chart::YAxis).setTitle("");
      //chart->setLeftYAxisPadding(double width, <#double height#>)
      chart->setPlotAreaPadding( 35, Wt::Left );
      chart->setPlotAreaPadding( 16, Wt::Bottom );
      chart->setPlotAreaPadding( 0, Wt::Right );
      chart->setPlotAreaPadding( 0, Wt::Top );
      
      WFont labelFont( WFont::Default );
      labelFont.setSize(8);
      chart->axis(Chart::XAxis).setLabelFont( labelFont );

      if( width < 125 ) //125 px arbitrarily chosen as widest to not show y-axis numbers
      {
        chart->setPlotAreaPadding( 6, Wt::Left );
        labelFont.setSize(1);
        chart->axis(Chart::YAxis).setLabelFont( labelFont );
        chart->axis(Chart::YAxis).setTextPen( WPen( WColor(0,0,0,0) ) );
      }else
      {
        chart->axis(Chart::YAxis).setLabelFont( labelFont );
      }
    }//if( compact )
    
    shared_ptr<Wt::WSvgImage> img = make_shared<WSvgImage>( width, height );
    
    WPainter p( img.get() );
    chart->paint( p );
    p.end();
    
    return img;
  }//renderChartToSvg(...)
  
  
void fit_peak_from_double_click( InterSpec *interspec, const double x, const double pixPerKeV,
                                shared_ptr<const DetectorPeakResponse> det,
                                std::string ref_line_name,
                                const SpecUtils::SpectrumType spec_type )
{
  PeakModel *pmodel = interspec ? interspec->peakModel() : nullptr;
  assert( pmodel );
  if( !pmodel )
    return;

  shared_ptr<const SpecUtils::Measurement> data = interspec->displayedHistogram(spec_type);
  if( !data )
    return;

  const bool isHPGe = PeakFitUtils::is_likely_high_res( interspec );

  // Get the appropriate SpecMeas and sample numbers for the spectrum type
  shared_ptr<SpecMeas> meas = interspec->measurment( spec_type );
  const std::set<int> sample_nums = interspec->displayedSamples( spec_type );

  vector<shared_ptr<const PeakDef>> origPeaks, nearbyOrigPeaks;

  // Get peaks from the appropriate SpecMeas for the given spectrum type
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> orig_peaks = meas ? meas->peaks(sample_nums) : nullptr;
  if( orig_peaks )
  {
    for( const std::shared_ptr<const PeakDef> &p : *orig_peaks )
    {
      origPeaks.push_back( p );

      if( ((x >= p->lowerX()) && (x <= p->upperX())) || (fabs(x - p->mean()) < 3.0*p->fwhm()) )
        nearbyOrigPeaks.push_back( p );

      //Avoid fitting a peak in the same area a data defined peak is.
      if( !p->gausPeak() && (x >= p->lowerX()) && (x <= p->upperX()) )
        return;
    }
  }//if( orig_peaks )

  shared_ptr<const deque<shared_ptr<const PeakDef>>> auto_peaks = meas ? meas->automatedSearchPeaks(sample_nums) : nullptr;

  // Get an estimated peak FWHM

  pair< PeakShrdVec, PeakShrdVec > foundPeaks;
  foundPeaks = searchForPeakFromUser( x, pixPerKeV, data, origPeaks, det, auto_peaks, isHPGe );
  
  //cerr << "Found " << foundPeaks.first.size() << " peaks to add, and "
  //     << foundPeaks.second.size() << " peaks to remove" << endl;
  
  if( foundPeaks.first.empty()
      || foundPeaks.second.size() >= foundPeaks.first.size() )
  {
    char msg[256];
    snprintf( msg, sizeof(msg), "Couldn't find peak a peak near %.1f keV", x );
    passMessage( msg, 0 );
    return;
  }//if( foundPeaks.first.empty() )
  
  
  //We want to add all of the previously found peaks back into the model, before
  //  adding the new peak.
  vector<shared_ptr<const PeakDef>> peakstoadd( foundPeaks.first.begin(), foundPeaks.first.end() );
  vector<shared_ptr<const PeakDef>> existingpeaks( foundPeaks.second.begin(), foundPeaks.second.end() );
  
  // For ROIs that will now have multiple peaks, we need to match up the previous peak(s) with the new peaks.
  // We know that the newly added peak wont have a source associated with it.
  // We will remove peaks from `peakstoadd`, and put them into `replacement_peaks`.
  vector<shared_ptr<const PeakDef>> replacement_peaks;

  
  //First add all the new peaks that have a nuclide/reaction/xray associated
  //  with them, since we know they are existing peaks
  for( const PeakModel::PeakShrdPtr &p : foundPeaks.first ) //note, not looping over `peakstoadd`, as we are modifying that vector.
  {
    if( p->parentNuclide() || p->reaction() || p->xrayElement() )
    {
      //find nearest previously existing peak, and add new peak, while removing
      //  old one from existingpeaks
      int nearest = -1;
      double smallesdist = DBL_MAX;
      for( size_t i = 0; i < existingpeaks.size(); ++i )
      {
        const PeakModel::PeakShrdPtr &prev = existingpeaks[i];
        if( prev->parentNuclide() != p->parentNuclide() )
          continue;
        if( prev->reaction() != p->reaction() )
          continue;
        if( prev->xrayElement() != p->xrayElement() )
          continue;
        
        const double thisdif = fabs(p->mean() - prev->mean());
        if( thisdif < smallesdist )
        {
          nearest = static_cast<int>( i );
          smallesdist = thisdif;
        }
      }
      
      if( nearest >= 0 )
      {
        replacement_peaks.push_back( p );
        existingpeaks.erase( existingpeaks.begin() + nearest );
        peakstoadd.erase( std::find(peakstoadd.begin(), peakstoadd.end(), p) );
      }//if( nearest >= 0 )
    }//if( p->parentNuclide() || p->reaction() || p->xrayElement() )
  }//for( const PeakModel::PeakShrdPtr &p : peakstoadd )
  
  //Now go through and add the new versions of the previously existing peaks,
  //  using energy to match the previous to current peak.
  for( const PeakModel::PeakShrdPtr &p : existingpeaks )
  {
    size_t nearest = 0;
    double smallesdist = DBL_MAX;
    for( size_t i = 0; i < peakstoadd.size(); ++i )
    {
      const double thisdif = fabs(p->mean() - peakstoadd[i]->mean());
      if( thisdif < smallesdist )
      {
        nearest = i;
        smallesdist = thisdif;
      }
    }
    
    std::shared_ptr<const PeakDef> peakToAdd = peakstoadd[nearest];
    peakstoadd.erase( begin(peakstoadd) + nearest );
    replacement_peaks.push_back( peakToAdd );
  }//for( const PeakModel::PeakShrdPtr &p : existingpeaks )
  
  
  //We'll check that we didnt mess up the associations between the existing peak an
  //  their respective new version, for all peaks with a source associated with them
  for( const PeakModel::PeakShrdPtr &p : peakstoadd )
  {
    assert( !p->parentNuclide() && !p->reaction() && !p->xrayElement() );
  }
  
  // Now check we added exactly one new peak, and it doesnt yet hav a source associated with it
  assert( (peakstoadd.size() == 1) && !peakstoadd.front()->hasSourceGammaAssigned() );
  

  if( spec_type != SpecUtils::SpectrumType::Foreground )
  {
    auto updated_peaks = make_shared<deque<std::shared_ptr<const PeakDef>>>();
    if( orig_peaks )
      updated_peaks->insert( end(*updated_peaks), begin(*orig_peaks), end(*orig_peaks) );

    // Erase any peaks from `existingpeaks` that are in `updated_peaks`
    updated_peaks->erase( std::remove_if( begin(*updated_peaks), end(*updated_peaks),
      [&foundPeaks](const std::shared_ptr<const PeakDef> &peak) {
        return (std::find( begin(foundPeaks.second), end(foundPeaks.second), peak ) != end(foundPeaks.second));
    } ), end(*updated_peaks) );

    // Insert peaks from `peakstoadd` to `updated_peaks`
    updated_peaks->insert( end(*updated_peaks), begin(replacement_peaks), end(replacement_peaks) );
    std::sort( begin(*updated_peaks), end(*updated_peaks), &PeakDef::lessThanByMeanShrdPtr );
    
    interspec->setPeaks( spec_type, updated_peaks );
  }else
  {
    pmodel->removePeaks( {begin(foundPeaks.second), end(foundPeaks.second)} );
    pmodel->addPeaks( replacement_peaks );
  }//if( spec_type != SpecUtils::SpectrumType::Foreground ) / else
  
  for( shared_ptr<const PeakDef> p : peakstoadd )
  {
    assert( !p->parentNuclide() && !p->reaction() && !p->xrayElement() );
     
    interspec->addPeak( *p, true, spec_type, ref_line_name );
  }
}//void fit_peak_from_double_click(...)

  
void automated_search_for_peaks( InterSpec *viewer,
                                const bool keep_old_peaks )
{
  if( !viewer )
    return; //shouldnt happen
  
  viewer->useMessageResourceBundle( "PeakSearchGuiUtils" );
  
  auto peakModel = viewer->peakModel();
  if( !peakModel )
    return; //shouldnt happen
  
  
  auto dataPtr = viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  
  if( !dataPtr )
  {
    passMessage( "There is no data to search for peaks on.", 1 );
    return;
  }//if( !dataPtr )
  
  //Could do a check to see if the user peaks have changed since the hint peaks
  //  were searched for, and if not, just use those...
  
  
  vector<ReferenceLineInfo> displayed;
  auto refLineDisp = viewer->referenceLinesWidget();
  if( refLineDisp )
  {
    const ReferenceLineInfo &currentNuclide = refLineDisp->currentlyShowingNuclide();
    displayed = refLineDisp->persistedNuclides(); //refLineDisp->showingNuclides()
  
    if( currentNuclide.m_validity == ReferenceLineInfo::InputValidity::Valid )
      displayed.insert( displayed.begin(), currentNuclide );
  }
  
  const bool setColor = viewer->colorPeaksBasedOnReferenceLines();
  
  
  //We should indicate to the user that seraching for peaks will take a while,
  //  but also we want to give them the chance to go past this and keep using
  //  the app.
  const WString title = WString::tr("psgu-search-wait-title");
  const WString content = WString::tr("psgu-search-wait-content");
  SimpleDialog *msg = new SimpleDialog( title, content );
  msg->rejectWhenEscapePressed();
  msg->addButton( WString::tr("Close") );
  
  //Make it so users cant keep clicking the search button
  viewer->automatedPeakSearchStarted();
  
  const auto originalPeaks = peakModel->peakVec();
  auto startingPeaks = peakModel->peaks();
  
  if( !keep_old_peaks )
    startingPeaks.reset();
  
  //using WApplication::bind to call msg->hide() will protect against the user
  //  closing the msg window (which will have deleted it)
  boost::function<void(void)> guiupdater
                          = wApp->bind( boost::bind( &SimpleDialog::accept, msg ) );
  
  //The results of the peak search will be placed into the vector pointed to
  // by searchresults, which is why both 'callback' and below and
  // search_for_peaks_worker(...) get a shared pointer to this vector.
  auto searchresults = std::make_shared< vector<std::shared_ptr<const PeakDef> > >();
  
  //Again, use WApplication::bind to create the call that will set the peaks
  //  when the searching is done.
  boost::function<void(void)> callback
                          = wApp->bind( boost::bind( &set_peaks_from_search,
                                        viewer, displayed, searchresults, dataPtr, originalPeaks, guiupdater ) );
  
  auto foreground = viewer->measurment(SpecUtils::SpectrumType::Foreground);
  shared_ptr<const DetectorPeakResponse> drf = foreground ? foreground->detector() : nullptr;

  const bool isHPGe = PeakFitUtils::is_likely_high_res( viewer );
  
  Wt::WServer *server = Wt::WServer::instance();
  if( !server )
    return;
  
  std::weak_ptr<const SpecUtils::Measurement> weakdata = dataPtr;
  const string seshid = wApp->sessionId();
  
  server->ioService().boost::asio::io_service::post( std::bind( [=](){
    search_for_peaks_worker( weakdata, drf, startingPeaks, displayed, setColor,
                            searchresults, callback, seshid, false, isHPGe );
    
  } ) );
}//void automated_search_for_peaks( InterSpec *interspec, const bool keep_old_peaks )


/** Assigns the passed in peak to a reference line if reasonable.  If this
 assignment would cause one of the already existing peaks to change assignment
 (most likely if two peaks of diff amp, but same nuc, are right next to each
 other) then the peak that should be changed, along with the string of the
 assignment it should be changed to.

 TODO: I dont think the return value is what is claimed, or at least not handled correctly; this should all be improved
 */
std::unique_ptr<std::pair<PeakModel::PeakShrdPtr,std::string>>
      assign_nuc_from_ref_lines( PeakDef &peak,
                                 std::shared_ptr<const deque< PeakModel::PeakShrdPtr > > previouspeaks,
                                 const std::shared_ptr<const SpecUtils::Measurement> &data,
                                 const vector<ReferenceLineInfo> &displayed,
                                 const bool colorPeaksBasedOnReferenceLines,
                                 const bool showingEscapePeakFeature )
{
  std::unique_ptr<std::pair<PeakModel::PeakShrdPtr,std::string>> other_change;

  if( !data || !previouspeaks || displayed.empty() )
    return other_change;


  //There is a fairly common situation (especially for HPGe) where there is a
  //  small peak, next to a much larger peak, where if the user first
  //  identifies the large peak, the correct gamma-ray association gets made,
  //  but then whe the second one is identified, it also gets assigned the
  //  same gamma-ray association as the larger peak, which is incoorect.  To
  //  avoid this, we will use previouspeaks and prevpeak to check, and correct
  //  for this condition.

  double prevPeakDist = DBL_MAX, prevIntensity = DBL_MAX;
  double thisIntensity = DBL_MAX;
  PeakDef::SourceGammaType thisGammaType = PeakDef::SourceGammaType::NormalGamma;

  PeakModel::PeakShrdPtr prevpeak;


  //The efficiency of S.E. and D.E. peaks, relative to F.E. peak, for the 20% Generic GADRAS DRF
  //  included in InterSpec, is given pretty well by the following (energy in keV):
  const auto single_escape_sf = []( const double x ) -> double {
    return std::max( 0.0, (1.8768E-11 *x*x*x) - (9.1467E-08 *x*x) + (2.1565E-04 *x) - 0.16367 );
  };

  const auto double_escape_sf = []( const double x ) -> double {
    return std::max( 0.0, (1.8575E-11 *x*x*x) - (9.0329E-08 *x*x) + (2.1302E-04 *x) - 0.16176 );
  };


  try
  {
    // We will also check if single or double escape peak.
    //  We will only check for escape peaks if HPGe, or user is showing escape peak features,
    //    or if the energy is above 4 MeV (arbitrarily chosen).
    //  We will assume a generic 20% HPGe detector to get S.E. and D.E. efficiencies.
    //  And on top of all those assumptions, we will apply an arbitrary amplitude factor of 0.5
    //    to reactions if we arent showing escape peak features, and 0.2 for xrays
    const float pair_prod_thresh = 1255.0f; //The single_escape_sf and double_escape_sf give negative values below 1255.
    const float always_check_escape_thresh = 4000.0f;
    const float escape_suppression_factor = 0.5;
    const float xray_suppression_factor = 0.2;

    // We will consider escape peaks if it is likely a HPGe detector.
    //  We'll judge this off the current peak width
    vector<std::shared_ptr<const PeakDef>> peakv( 1, make_shared<const PeakDef>(peak) );
    const PeakFitUtils::CoarseResolutionType resolutionType
                              = PeakFitUtils::coarse_resolution_from_peaks( peakv );
    bool isHPGe = (resolutionType == PeakFitUtils::CoarseResolutionType::High);
    if( !isHPGe && data && (data->num_gamma_channels() > 1024) && InterSpec::instance() )
      isHPGe = PeakFitUtils::is_likely_high_res( InterSpec::instance() );


    double mindist = 99999999.9;
    double nearestEnergy = -999.9;

    const double minx = peak.lowerX();
    const double maxx = peak.upperX();

    const SandiaDecay::Nuclide *nuclide = NULL;
    const SandiaDecay::Element *element = NULL;
    const ReactionGamma::Reaction *reaction = NULL;

    const double mean = peak.mean();
    const double sigma = peak.gausPeak() ? peak.sigma() : peak.roiWidth();

    Wt::WColor color;

    for( const ReferenceLineInfo &nuc : displayed )
    {
      // Make sure the source is a type we can potentially assign to a peak
      switch( nuc.m_source_type )
      {
        case ReferenceLineInfo::SourceType::Nuclide:
        case ReferenceLineInfo::SourceType::FluorescenceXray:
        case ReferenceLineInfo::SourceType::Reaction:
        case ReferenceLineInfo::SourceType::Background:
        case ReferenceLineInfo::SourceType::NuclideMixture:
        case ReferenceLineInfo::SourceType::OneOffSrcLines:
        case ReferenceLineInfo::SourceType::FissionRefLines:
          break;

        case ReferenceLineInfo::SourceType::CustomEnergy:
        case ReferenceLineInfo::SourceType::None:
          continue;
          break;
      }//switch( nuc.m_source_type )


      for( const ReferenceLineInfo::RefLine &line : nuc.m_ref_lines )
      {
        double energy = line.m_energy;
        PeakDef::SourceGammaType gammaType = PeakDef::SourceGammaType::NormalGamma;

        switch( line.m_source_type )
        {
          case ReferenceLineInfo::RefLine::RefGammaType::SingleEscape:
            energy += 510.998950;
            gammaType = PeakDef::SourceGammaType::SingleEscapeGamma;
            break;

          case ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape:
            energy += 2.0*510.998950;
            gammaType = PeakDef::SourceGammaType::DoubleEscapeGamma;
            break;

          case ReferenceLineInfo::RefLine::RefGammaType::Annihilation:
            gammaType = PeakDef::SourceGammaType::AnnihilationGamma;
            break;

          case ReferenceLineInfo::RefLine::RefGammaType::Normal:
          case ReferenceLineInfo::RefLine::RefGammaType::SumGammaPeak:
          case ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak:
            gammaType = PeakDef::SourceGammaType::NormalGamma;
            break;
        }//switch( line.m_source_type )


        double expectedPhotoPeakEnergy = energy;
        double intensity = line.m_normalized_intensity;

        switch( line.m_particle_type )
        {
          case ReferenceLineInfo::RefLine::Particle::Gamma:
            break;

          case ReferenceLineInfo::RefLine::Particle::Xray:
            intensity *= xray_suppression_factor;
            gammaType = PeakDef::SourceGammaType::XrayGamma;
            break;

          case ReferenceLineInfo::RefLine::Particle::Alpha:
          case ReferenceLineInfo::RefLine::Particle::Beta:
            continue;
            break;
        }//switch( line.m_particle_type )

        if( IsInf(intensity) || IsNan(intensity) || (intensity < numeric_limits<float>::min()) )
          continue;

        const double delta_e = fabs( mean - energy );
        double dist = (0.25*sigma + delta_e) / intensity;

        if( (line.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::Normal)
           && ( ((isHPGe || showingEscapePeakFeature) && (energy > pair_prod_thresh))
               || (energy > always_check_escape_thresh) ) )
        {
          const double suppress_sf = showingEscapePeakFeature ? 1.0 : escape_suppression_factor;

          double se_abundance = intensity * single_escape_sf(energy) * suppress_sf;
          const double se_delta_e = fabs( mean + 510.9989 - energy );
          const double se_dist = (0.25*sigma + se_delta_e) / se_abundance;

          double de_abundance = intensity * double_escape_sf(energy) * suppress_sf;
          const double de_delta_e = fabs( mean + (2*510.9989) - energy );
          const double de_dist = (0.25*sigma + de_delta_e) / de_abundance;

          if( se_dist < dist )
          {
            dist = se_dist;
            intensity = se_abundance;
            expectedPhotoPeakEnergy = energy - 510.9989;
            gammaType = PeakDef::SourceGammaType::SingleEscapeGamma;
          }

          if( de_dist < dist )
          {
            dist = de_dist;
            intensity = de_abundance;
            expectedPhotoPeakEnergy = energy - 2*510.9989;
            gammaType = PeakDef::SourceGammaType::DoubleEscapeGamma;
          }
        }//if( check for S.E. or D.E. )


        if( (dist < mindist)
           && (expectedPhotoPeakEnergy >= minx) && (expectedPhotoPeakEnergy <= maxx) )
        {
          bool currentlyused = false;
          for( const PeakModel::PeakShrdPtr &pp : *previouspeaks )
          {
            const bool sameNuc = (line.m_parent_nuclide
                                  && pp->decayParticle()
                                  && (pp->parentNuclide() == line.m_parent_nuclide)
                                  && (fabs(pp->decayParticle()->energy - energy) < 0.01)
                                  && (pp->sourceGammaType() == gammaType));
            const bool sameXray = ((line.m_particle_type == ReferenceLineInfo::RefLine::Particle::Xray)
                                   && line.m_element
                                   && (pp->xrayElement() == line.m_element)
                                   && (fabs(pp->xrayEnergy() - energy) < 0.01)
                                   && (pp->sourceGammaType() == gammaType));

            if( sameNuc || sameXray )
            {
              currentlyused = true;
              if( dist < prevPeakDist )
              {
                prevpeak = pp;
                prevPeakDist = dist;
                prevIntensity = intensity;
              }
              break;
            }//if( pp->reaction() == rpp.reaction )
          }//for( const PeakModel::PeakShrdPtr &pp : *previouspeaks )

          if( !currentlyused )
          {
            prevpeak.reset();
            mindist = dist;
            color = line.m_color.isDefault() ? nuc.m_input.m_color : line.m_color;
            nuclide = line.m_parent_nuclide;
            element = line.m_element;
            reaction = line.m_reaction;
            nearestEnergy = energy;
            thisIntensity = intensity;
            thisGammaType = gammaType;
          }//if( !currentlyused )
        }//if( we should possible associate this peak with this line )
      }//for( const double energy : nuc.energies )
    }//if( nuc.nuclide )


    if( nuclide || reaction || element )
    {
      assert( mindist >= 0.0 );

      string src;
      char nuclide_label[128];

      const char *prefix = "", *postfix = "";
      switch( thisGammaType )
      {
        case PeakDef::NormalGamma:       break;
        case PeakDef::AnnihilationGamma: break;
        case PeakDef::SingleEscapeGamma: prefix  = "S.E. "; break;
        case PeakDef::DoubleEscapeGamma: prefix  = "D.E. "; break;
        case PeakDef::XrayGamma:         postfix = " xray"; break;
      }//switch( thisGammaType )


      if( !!prevpeak )
      {
#ifdef _MSC_VER
#pragma message("Need to clean up setting peak source, and also verify S.E., and D.E., are actually all correct")
#else
#warning "Need to clean up setting peak source, and also verify S.E., and D.E., are actually all correct"
#endif

        //There is an already existing peak, whos nuclide/xray/reaction gamma was
        //  "closer" (in terms of the distance metric used) than the
        //  nuclide/xray/reaction actually assigned to this new peak.  We need
        //  to check that we shouldnt swap them, based on relative intensities.
        //
        //  TODO: Havent fully verified this logic when a S.E. or D.E. peak is involved
        const bool prevAmpSmaller = (prevpeak->amplitude()<peak.amplitude());
        const bool prevIntenitySmaller = (prevIntensity < thisIntensity);
        if( prevAmpSmaller != prevIntenitySmaller )
        {
          double prevEnergy = 0.0;
          const SandiaDecay::Nuclide *prevnuclide = prevpeak->parentNuclide();
          const SandiaDecay::Element *prevelement = prevpeak->xrayElement();
          const ReactionGamma::Reaction *prevreaction = prevpeak->reaction();

          if( prevnuclide && prevpeak->parentNuclide() && prevpeak->decayParticle() )
          {
            prevEnergy = prevpeak->decayParticle()->energy;
            src = prevpeak->parentNuclide()->symbol;
          }else if( prevelement && prevpeak->xrayElement() )
          {
            prevEnergy = prevpeak->xrayEnergy();
            src = prevpeak->xrayElement()->symbol;
          }else if( prevreaction && prevpeak->reaction() )
          {
            prevEnergy = prevpeak->reactionEnergy();
            src = prevpeak->reaction()->name();
          }else
          {
            throw std::logic_error( "InterSpec::addPeak(): bad logic "
                                   "checking previous peak gamma assignment" );
          }

          snprintf( nuclide_label, sizeof(nuclide_label),
                   "%s%s%s %.6f keV", prefix, postfix, src.c_str(), nearestEnergy );


          other_change.reset( new std::pair<PeakModel::PeakShrdPtr,std::string>(prevpeak,nuclide_label) );

          nuclide = prevnuclide;
          reaction = prevreaction;
          element = prevelement;
          nearestEnergy = prevEnergy;
        }//if( we need to swap things )
      }//if( !!prevpeak )

      if( nuclide )
      {
        assert( thisGammaType == PeakDef::NormalGamma
               || thisGammaType == PeakDef::SingleEscapeGamma
               || thisGammaType == PeakDef::DoubleEscapeGamma
               || thisGammaType == PeakDef::XrayGamma
               || thisGammaType == PeakDef::AnnihilationGamma );

        src = nuclide->symbol;
      }else if( reaction )
      {
        assert( thisGammaType == PeakDef::NormalGamma
               || thisGammaType == PeakDef::SingleEscapeGamma
               || thisGammaType == PeakDef::DoubleEscapeGamma );

        src = reaction->name();
      }else if( element )
      {
        assert( thisGammaType == PeakDef::XrayGamma );
        src = element->symbol;
      }

      snprintf( nuclide_label, sizeof(nuclide_label),
               "%s%s%s %.6f keV", prefix, src.c_str(), postfix, nearestEnergy );

      // TODO: do we really need to use #PeakModel::setNuclideXrayReaction ?  We should be able to just directly set information
      PeakModel::setNuclideXrayReaction( peak, nuclide_label, -1.0 );

      if( colorPeaksBasedOnReferenceLines /* && peak.lineColor().isDefault() */ )
        peak.setLineColor( color );
    }//if( nuclide || reaction || element )
  }catch( std::exception &e )
  {
    passMessage( "Unexpected error searching for isotope for peak: " + string(e.what()),
                WarningWidget::WarningMsgHigh );
  }

  return other_change;
}//void assignNuclideFromReferenceLines( PeakDef &peak )


  
void assign_peak_nuclides_from_reference_lines( InterSpec *viewer,
                                               const bool only_peaks_with_no_src,
                                               const bool only_current_ref_lines )
{
  if( !viewer )
    return;
  
  //Add assigning nuclides from reference lines in set_peaks_from_search(...)
  //  then create the popup window widget to let the user change things - this
  //  function should also use that same popup window!
  
  const bool showingEscapePeakFeature = viewer->showingFeatureMarker(FeatureMarkerType::EscapePeakMarker);
  auto peakModel = viewer->peakModel();
  auto foreground = viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  const ReferencePhotopeakDisplay *refLineDisp = viewer->referenceLinesWidget();
  
  if( !peakModel || !foreground || !refLineDisp )
    return;
  
  const bool assignColor = viewer->colorPeaksBasedOnReferenceLines();
  
  vector<PeakDef> peaks = peakModel->peakVec();
  vector<PeakDef> orig_peaks = peaks;
  vector<PeakDef> result_peaks;
  
  vector<ReferenceLineInfo> displayed;
  
  {
    if( !only_current_ref_lines )
      displayed = refLineDisp->persistedNuclides();
    
    const ReferenceLineInfo &currentNuclide = refLineDisp->currentlyShowingNuclide();
    if( currentNuclide.m_validity == ReferenceLineInfo::InputValidity::Valid )
      displayed.insert( displayed.begin(), currentNuclide );
  }
  
  vector<shared_ptr<const PeakDef>> resultpeaks;
  for( const PeakDef &p : peaks )
    resultpeaks.push_back( make_shared<PeakDef>(p) );
  
  resultpeaks = assign_srcs_from_ref_lines( foreground, resultpeaks, displayed, assignColor,
                             showingEscapePeakFeature, only_peaks_with_no_src );
  
  for( const auto &foundpeak : resultpeaks )
  {
    PeakDef p = *foundpeak;
    
    //Do a sanity check to make sure line is within a few sigma - we want to
    //  be a little more tight than when a user double-clicks on a peak.
    if( p.parentNuclide() || p.reaction() || p.xrayElement() )
    {
      const double gammaEnergy = p.gammaParticleEnergy();
      const double mean = p.mean();
      const double fwhm = (p.gausPeak() ? p.fwhm() : (p.roiWidth()/4));
      if( (fabs(gammaEnergy - mean) / fwhm) > 2.5 ) //2.5 FWHM arbitrary
      {
        p.clearSources();
      }
    }//if( a source was found )
    
    result_peaks.push_back( p );
  }//for( const auto &foundpeak : resultpeaks )
  
  
  //Hack because I'm short on time
  if( result_peaks.size() != orig_peaks.size() )
  {
    passMessage( "Unexpected issue assigning nuclides - accepting all assignments - sorry :(",
                WarningWidget::WarningMsgHigh );
    return;
  }
  
  auto peak_less_than_by_energy = [](const PeakDef &lhs, const PeakDef &rhs) -> bool {
    return lhs.mean() < rhs.mean();
  };
  
  std::sort( begin(orig_peaks), end(orig_peaks), peak_less_than_by_energy );
  std::sort( begin(result_peaks), end(result_peaks), peak_less_than_by_energy );
  
  for( size_t i = 0; i < orig_peaks.size(); ++i )
  {
    if( fabs(orig_peaks[i].mean() - result_peaks[i].mean()) > 1.0 )
    {
      passMessage( "Unexpected issue assigning nuclides, mismatch in energy - accepting all assignments - sorry :(",
                  WarningWidget::WarningMsgHigh );
      return;
    }
  }//for( size_t i = 0; i < peaks.size(); ++i )
  
  peakModel->setPeaks( result_peaks );

  new PeakSelectorWindow( viewer, PeakSelectorWindowReason::NuclideId, orig_peaks,
                          foreground, result_peaks, displayed );
  
  wApp->triggerUpdate();
}//void assign_peak_nuclides_from_reference_lines()
  
  
void assign_nuclide_from_reference_lines( PeakDef &peak,
                                       PeakModel *peakModel,
                                       const std::shared_ptr<const SpecUtils::Measurement> &data,
                                       const ReferencePhotopeakDisplay *refLineDisp,
                                       const bool setColor,
                                       const bool showingEscapePeakFeature )
{
  if( !data || !refLineDisp || !peakModel )
    return;
  
  const ReferenceLineInfo &currentNuclide = refLineDisp->currentlyShowingNuclide();
  vector<ReferenceLineInfo> displayed = refLineDisp->persistedNuclides();
  
  if( currentNuclide.m_validity == ReferenceLineInfo::InputValidity::Valid )
    displayed.insert( displayed.begin(), currentNuclide );
  
  if( displayed.empty() )
    return;
  
  std::shared_ptr<const deque< PeakModel::PeakShrdPtr > > previouspeaks = peakModel->peaks();
  if( !previouspeaks )  //probably never necessary, but JIC
    previouspeaks = std::make_shared<deque< PeakModel::PeakShrdPtr > >();
  
  unique_ptr<pair<PeakModel::PeakShrdPtr,std::string>> addswap
    = assign_nuc_from_ref_lines( peak, previouspeaks, data, displayed, setColor, showingEscapePeakFeature );
  
  if( addswap )
  {
    WModelIndex prevpeakind = peakModel->indexOfPeak( addswap->first );
    
    if( prevpeakind.isValid() )
    {
      prevpeakind = peakModel->index(prevpeakind.row(), PeakModel::kIsotope);
      peakModel->setData( prevpeakind, WString(addswap->second) );
    }else
    {
      const char *msg = "assign_nuclide_from_reference_lines(): bad logic getting"
                        " previous peak index";
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, msg );
#else
      cerr << msg << endl;
#endif
    }//if( prevpeakind.isValid() ) / else
  }//if( addswap )
}//void assignNuclideFromReferenceLines( PeakDef &peak )



  
void search_for_peaks_worker( std::weak_ptr<const SpecUtils::Measurement> weak_data,
                              shared_ptr<const DetectorPeakResponse> drf,
                               std::shared_ptr<const deque< std::shared_ptr<const PeakDef> > > existingPeaks,
                               const vector<ReferenceLineInfo> displayed,
                               const bool setColor,
                               std::shared_ptr<std::vector<std::shared_ptr<const PeakDef> > > resultpeaks,
                               boost::function<void(void)> callback,
                               const std::string sessionID,
                               const bool singleThread,
                               const bool isHPGe )
{
  Wt::WServer *server = Wt::WServer::instance();
  if( !server )  //shouldnt ever happen,
    return;
  
  std::shared_ptr<const SpecUtils::Measurement> data = weak_data.lock();
  
  if( !data || !resultpeaks )
  {
    server->post( sessionID, callback );
    return;
  }

  auto results = make_shared<vector<shared_ptr<const PeakDef>>>();

  try
  {
    *results = ExperimentalAutomatedPeakSearch::search_for_peaks( data, drf, existingPeaks, singleThread, isHPGe );
    *results = assign_srcs_from_ref_lines( data, *results, displayed, setColor, false, true );
  }catch( std::exception &e )
  {
    string msg = "InterSpec::search_for_peaks_worker(): caught exception: '";
    msg += e.what();
    msg += "'";
    
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, msg.c_str() );
#else
    cerr << msg << endl;
#endif
  }//try / catch

  //assert( 0 );
  //throw runtime_error( "Need to implement search for peaks worker." );

  if( results && !results->empty() )
  {
    try
    {
      //DetectorPeakResponse::ResolutionFnctForm fwhm_type = ...;
      //const int sqrtEqnOrder = ...;
      //auto peaks_deque = make_shared<deque<shared_ptr<const PeakDef>>>( begin(*results), end(*results) );

      //vector<float> result, uncerts;
      //const double chi2 = MakeDrfFit::performResolutionFit( peaks_deque, fwhm_type, sqrtEqnOrder, result, uncerts );
    }catch( std::exception &e )
    {

    }
  }

  server->post( sessionID, [callback,results,resultpeaks](){
    if( resultpeaks )
      *resultpeaks = *results;
    callback();
  } );
}//void search_for_peaks_worker(...)

  
std::vector<std::shared_ptr<const PeakDef>> assign_srcs_from_ref_lines( const std::shared_ptr<const SpecUtils::Measurement> &data,
                                     const std::vector<std::shared_ptr<const PeakDef>> &peaks,
                                     const ReferencePhotopeakDisplay * const refwidget,
                                     const bool setColor,
                                     const bool showingEscapePeakFeature,
                                     const bool only_peaks_with_no_src )
{
  if( !data || !refwidget || peaks.empty() )
    return peaks;
  
  const ReferenceLineInfo &currentNuclide = refwidget->currentlyShowingNuclide();
  vector<ReferenceLineInfo> displayed = refwidget->persistedNuclides();
  if( currentNuclide.m_validity == ReferenceLineInfo::InputValidity::Valid )
    displayed.insert( displayed.begin(), currentNuclide );
  
  if( displayed.empty() )
    return peaks;
  
  return assign_srcs_from_ref_lines( data, peaks, displayed, setColor,
                                    showingEscapePeakFeature, only_peaks_with_no_src );
}//assign_srcs_from_ref_lines(....)
  
  
vector<shared_ptr<const PeakDef>> assign_srcs_from_ref_lines( const std::shared_ptr<const SpecUtils::Measurement> &data,
                                  const std::vector<std::shared_ptr<const PeakDef>> &input_peaks,
                                  const vector<ReferenceLineInfo> &displayed,
                                  const bool setColor,
                                  const bool showingEscapePeakFeature,
                                  const bool only_peaks_with_no_src )
{
  if( input_peaks.empty() || !data || displayed.empty() )
    return input_peaks;
  
  auto answerpeaks = make_shared< deque< shared_ptr<const PeakDef> > >();
  vector< shared_ptr<PeakDef> > unassignedpeaks;
  
  for( const auto &p : input_peaks )
  {
    if( only_peaks_with_no_src && p->hasSourceGammaAssigned() )
      answerpeaks->push_back( p );
    else
      unassignedpeaks.push_back( make_shared<PeakDef>( *p ) );
  }
  
  auto peak_less_than_by_energy = [](const shared_ptr<const PeakDef> &lhs, const shared_ptr<const PeakDef> &rhs) -> bool {
    return lhs->mean() < rhs->mean();
  };
  std::sort( begin(*answerpeaks), end(*answerpeaks), peak_less_than_by_energy );
  
  auto peak_less_than_by_amp = [](const shared_ptr<PeakDef> &lhs, const shared_ptr<PeakDef> &rhs) -> bool {
    return lhs->amplitude() > rhs->amplitude();
  };
  std::sort( begin(unassignedpeaks), end(unassignedpeaks), peak_less_than_by_amp );
  
  
  for( shared_ptr<PeakDef> peak : unassignedpeaks )
  {
    unique_ptr<pair<PeakModel::PeakShrdPtr,string>> addswap
           = assign_nuc_from_ref_lines( *peak, answerpeaks, data, displayed, setColor, showingEscapePeakFeature );
    
    if( addswap )  //The assignment caused a better assignment to be made for a previously existing peak.
    {
      auto pos = std::find( begin(*answerpeaks), end(*answerpeaks), addswap->first );
      if( pos != end(*answerpeaks) )
      {
        shared_ptr<const PeakDef> &oldpeak = *pos;
        auto moddedOldPeak = make_shared<PeakDef>( *oldpeak );
        
        PeakModel::SetGammaSource res = PeakModel::setNuclideXrayReaction( *moddedOldPeak, addswap->second, 1.0 );
        oldpeak = moddedOldPeak;
        if( (res == PeakModel::SetGammaSource::NoSourceChange) || (res == PeakModel::SetGammaSource::FailedSourceChange) )
        {
#if( PERFORM_DEVELOPER_CHECKS )
          log_developer_error( __func__, "A suggested change to a peak did not result in a change - prob shouldnt have happened." );
#endif
        }
      }else
      {
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, "A suggested peak to change was not in the input - logic error!" );
#endif
      }//if( pos != end(*answerpeaks) ) / else
    }//if( addswap )
    
    
    auto sorted_pos = std::lower_bound( begin(*answerpeaks), end(*answerpeaks), peak, peak_less_than_by_energy );
    answerpeaks->insert( sorted_pos, peak );

    assert( std::is_sorted( begin(*answerpeaks), end(*answerpeaks), peak_less_than_by_energy ) );
    //answerpeaks->push_back( peak );
    //std::sort( begin(*answerpeaks), end(*answerpeaks), peak_less_than_by_energy );
  }//for( PeakDef peak : unassignedpeaks )
  
  return vector<shared_ptr<const PeakDef>>( begin(*answerpeaks), end(*answerpeaks) );
}//void assign_srcs_from_ref_lines(...)
  
  
void refit_peaks_from_right_click( InterSpec * const interspec, const double rightClickEnergy, const RefitPeakType type )
{
  if( type == RefitPeakType::WithDrfFwhm )
  {
    PeakSearchGuiUtils::refit_peaks_with_drf_fwhm( interspec, rightClickEnergy );
    return;
  }
  
  try
  {
    PeakModel * const model = interspec->peakModel();
    const shared_ptr<const SpecUtils::Measurement> data = interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    const shared_ptr<const SpecMeas> foreground = interspec->measurment( SpecUtils::SpectrumType::Foreground);
    
    
    if( !model || !data || !foreground ) //shouldnt ever happen
    {
      passMessage( "No data loaded to refit", WarningWidget::WarningMsgInfo );
      return;
    }
    
    //shared_ptr<const PeakDef> peak = interspec->nearestPeak( rightClickEnergy );
    shared_ptr<const PeakDef> peak = model->nearestPeak( rightClickEnergy );
    if( !peak )
    {
      passMessage( "There was no peak to refit", WarningWidget::WarningMsgInfo );
      return;
    }
    
    vector<PeakDef> inputPeak, outputPeak;
    const vector<shared_ptr<const PeakDef>> peaksInRoi = model->peaksSharingRoi( peak );
    const vector<shared_ptr<const PeakDef>> peaksNotInRoi = model->peaksNotSharingRoi( peak );
    
    assert( peaksInRoi.size() >= 1 );
    
    for( const auto &m : peaksInRoi )
      inputPeak.push_back( *m );

    std::sort( inputPeak.begin(), inputPeak.end(), &PeakDef::lessThanByMean );

    if( inputPeak.size() > 1 )
    {
      const shared_ptr<const DetectorPeakResponse> &detector = foreground->detector();
      
      vector<shared_ptr<const PeakDef>> result;
      WFlags<PeakFitLM::PeakFitLMOptions> fit_options;
      
      if( type == RefitPeakType::RoiAgressive )
        fit_options |= PeakFitLM::PeakFitLMOptions::AllPeakFwhmIndependent;

      result = refitPeaksThatShareROI( data, detector, peaksInRoi, fit_options );
      

      if( result.size() == inputPeak.size() )
      {
        for( const auto &m : peaksNotInRoi )
          result.push_back( m );

        std::sort( begin(result), end(result), &PeakDef::lessThanByMeanShrdPtr );
        model->setPeaks( result );
        return;
      }else
      {
        cerr << "refit_peaks_from_right_click was not successful" << endl;
      }//if( result.size() == inputPeak.size() ) / else
    }//if( inputPeak.size() > 1 )
    
    
    //  const double lowE = peak->mean() - 0.1;
    //  const double upE = peak->mean() + 0.1;
    const double lowE = inputPeak.front().mean() - 0.1;
    const double upE = inputPeak.back().mean() + 0.1;
    const double ncausalitysigma = 0.0;
    const double stat_threshold  = 0.0;
    const double hypothesis_threshold = 0.0;
    
    const bool isRefit = true;
    const bool isHPGe = PeakFitUtils::is_likely_high_res(interspec); //Only used if peak range is defined, so shouldnt be needed
    
    
    WFlags<PeakFitLM::PeakFitLMOptions> fit_options;
    
    switch( type )
    {
      case RefitPeakType::Standard:
        // No options to add
        break;
        
      case RefitPeakType::Delicate:
        fit_options |= PeakFitLM::PeakFitLMOptions::SmallRefinementOnly;
        break;
        
      case RefitPeakType::WithDrfFwhm:
        assert( 0 );
        break;
        
      case RefitPeakType::RoiAgressive:
        assert( inputPeak.size() > 1 );
        fit_options |= PeakFitLM::PeakFitLMOptions::AllPeakFwhmIndependent;
        break;
    }//switch( type )
    
    if( type == RefitPeakType::RoiAgressive )
      fit_options |= PeakFitLM::PeakFitLMOptions::AllPeakFwhmIndependent;
    
    outputPeak = fitPeaksInRange( lowE, upE, ncausalitysigma, stat_threshold,
                                 hypothesis_threshold, inputPeak, data, fit_options, isHPGe );
    if( outputPeak.size() != inputPeak.size() )
    {
      WStringStream msg;
      msg << "Failed to refit peak (became insignificant), from "
      << int(inputPeak.size()) << " to " << int(outputPeak.size()) << " peaks";
      passMessage( msg.str(), WarningWidget::WarningMsgInfo );
      return;
    }//if( outputPeak.size() != 1 )
    
    if( inputPeak.size() > 1 )
    {
      for( const shared_ptr<const PeakDef> &p : peaksNotInRoi )
        outputPeak.push_back( *p );

      std::sort( begin(outputPeak), end(outputPeak), &PeakDef::lessThanByMean );
      model->setPeaks( outputPeak );
    }else
    {
      assert( !outputPeak.empty() );
      
      model->updatePeak( peak, outputPeak[0] );
    }//if( inputPeak.size() > 1 )
  }catch( std::exception &e )
  {
    passMessage( "Sorry, error encountered refitting ROI.", WarningWidget::WarningMsgInfo );
    cerr << "Error encountered refitting ROI: " << e.what() << endl;
  }
}//void refit_peaks_from_right_click(...)

  
void refit_peaks_with_drf_fwhm( InterSpec * const interspec, const double rightClickEnergy )
{
  try
  {
    interspec->useMessageResourceBundle( "PeakSearchGuiUtils" );
    
    PeakModel * const model = interspec->peakModel();
    const shared_ptr<const SpecUtils::Measurement> data = interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    const shared_ptr<const SpecMeas> foreground = interspec->measurment( SpecUtils::SpectrumType::Foreground);
    
    
    if( !model || !data || !foreground ) //shouldnt ever happen
    {
      passMessage( "No data loaded to refit peaks", WarningWidget::WarningMsgInfo );
      return;
    }
    
    shared_ptr<const DetectorPeakResponse> drf = foreground->detector();
    if( !drf || !drf->hasResolutionInfo() )
    {
      const WString title = WString::tr("psgu-fwhm-info-needed-title");
      const WString content = WString::tr( !!drf ? "psgu-fwhm-info-needed-has-det-eff" : "psgu-fwhm-info-needed-no-det-eff" );
      
      SimpleDialog *msg = new SimpleDialog( title, content );
      // Setting object name of the SimpleDialog causes a javascript error - so we'll set
      //  the SimpleDialog contents name, and then get the dialog from that.
      //  This is all for undo/redo support, so we dont have to store a pointer to the dialog
      //  somewhere.
      msg->contents()->setObjectName( "AskToFitFwhmDialog" );
      msg->rejectWhenEscapePressed();
      WPushButton *yes_btn = msg->addButton( WString::tr("Yes") );
      WPushButton *no_btn = msg->addButton( WString::tr("No") );
      
      no_btn->clicked().connect( std::bind([interspec,rightClickEnergy](){
        auto undo = [interspec,rightClickEnergy](){
          refit_peaks_with_drf_fwhm( interspec, rightClickEnergy );
        };
        auto redo = [](){
          auto w = dynamic_cast<WContainerWidget *>( wApp->findWidget("AskToFitFwhmDialog") );
          WWidget *p = w ? w->parent() : nullptr;
          WWidget *pp = p ? p->parent() : nullptr;
          WWidget *ppp = pp ? pp->parent() : nullptr;
          SimpleDialog *d = dynamic_cast<SimpleDialog *>( ppp );
            
          if( d )
            d->done(Wt::WDialog::DialogCode::Accepted);
          wApp->doJavaScript( "$('.Wt-dialogcover').hide();" );
        };
        
        UndoRedoManager *undoManager = interspec->undoRedoManager();
        if( undoManager && undoManager->canAddUndoRedoNow() )
          undoManager->addUndoRedoStep( undo, redo, "Cancel refit peak with DRF FWHM." );
      } ) );
      
      
      yes_btn->clicked().connect( std::bind( [interspec,rightClickEnergy](){
        MakeFwhmForDrfWindow *window = interspec->fwhmFromForegroundWindow(true);
        if( window )
        {
          window->tool()->updatedDrf().connect(
            boost::bind( &refit_peaks_with_drf_fwhm, interspec, rightClickEnergy
          ) );
        }//if( window )
        
        auto undo = [interspec,rightClickEnergy](){
          interspec->deleteFwhmFromForegroundWindow();
          refit_peaks_with_drf_fwhm( interspec, rightClickEnergy );
        };
        auto redo = [interspec,rightClickEnergy](){
          auto w = dynamic_cast<WContainerWidget *>( wApp->findWidget("AskToFitFwhmDialog") );
          WWidget *p = w ? w->parent() : nullptr;
          WWidget *pp = p ? p->parent() : nullptr;
          WWidget *ppp = pp ? pp->parent() : nullptr;
          SimpleDialog *d = dynamic_cast<SimpleDialog *>( ppp );
          if( d )
            d->done(Wt::WDialog::DialogCode::Accepted);
          wApp->doJavaScript( "$('.Wt-dialogcover').hide();" );
          
          MakeFwhmForDrfWindow *window = interspec->fwhmFromForegroundWindow(true);
          if( window )
          {
            window->tool()->updatedDrf().connect(
              boost::bind( &refit_peaks_with_drf_fwhm, interspec, rightClickEnergy
            ) );
          }//if( window )
        };
        
        UndoRedoManager *undoManager = interspec->undoRedoManager();
        if( undoManager && undoManager->canAddUndoRedoNow() )
          undoManager->addUndoRedoStep( undo, redo, "Start fit FWHM function." );
      } ) );
      
      return;
    }//if( !drf || !drf->hasResolutionInfo() )
    
    shared_ptr<const PeakDef> peak = model->nearestPeak( rightClickEnergy );
    if( !peak )
    {
      passMessage( "There was no peak to refit with fixed FWHM", WarningWidget::WarningMsgInfo );
      return;
    }
    
    UndoRedoManager::PeakModelChange peak_undo_creator;
    
    vector<PeakDef> inputPeak, outputPeak;
    const vector<shared_ptr<const PeakDef>> peaksInRoi = model->peaksSharingRoi( peak );
    const vector<shared_ptr<const PeakDef>> peaksNotInRoi  = model->peaksNotSharingRoi( peak );
    
    assert( peaksInRoi.size() >= 1 );
    
    for( const auto &m : peaksInRoi )
    {
      PeakDef p = *m;
      
      p.setSigma( drf->peakResolutionSigma(p.mean()) );
      p.setSigmaUncert( 0.0 );
      //p.setSigmaUncert( ... DRF FWHM uncert not implemented ... )
      p.setFitFor( PeakDef::CoefficientType::Sigma, false );
      
      inputPeak.push_back( p );
    }

    
    std::sort( inputPeak.begin(), inputPeak.end(), &PeakDef::lessThanByMean );
    
    if( inputPeak.size() > 1 )
    {
      const shared_ptr<const DetectorPeakResponse> &detector = foreground->detector();
      WFlags<PeakFitLM::PeakFitLMOptions> fit_options;
      //fit_options |= PeakFitLM::PeakFitLMOptions::MediumRefinementOnly;
      
      vector<shared_ptr<const PeakDef>> result = refitPeaksThatShareROI( data, detector, peaksInRoi, fit_options );

      if( result.size() == inputPeak.size() )
      {
        for( const shared_ptr<const PeakDef> &p : peaksNotInRoi )
          result.push_back( p );
        std::sort( begin(result), end(result), &PeakDef::lessThanByMeanShrdPtr );
        model->setPeaks( result );
        return;
      }else
      {
        cerr << "refit_peaks_with_drf_fwhm was not successful" << endl;
      }//if( result.size() == inputPeak.size() ) / else
    }//if( inputPeak.size() > 1 )
    
    
    //  const double lowE = peak->mean() - 0.1;
    //  const double upE = peak->mean() + 0.1;
    const double lowE = inputPeak.front().mean() - 0.1;
    const double upE = inputPeak.back().mean() + 0.1;
    const double ncausalitysigma = 0.0;
    const double stat_threshold  = -1000.0;
    const double hypothesis_threshold = -1000.0;
    
    Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options;
    const bool isHPGe = PeakFitUtils::is_likely_high_res(interspec); //Only used if peak range is defined, so shouldnt be needed
    outputPeak = fitPeaksInRange( lowE, upE, ncausalitysigma, stat_threshold,
                                 hypothesis_threshold, inputPeak, data, fit_options, isHPGe );
    if( outputPeak.size() != inputPeak.size() )
    {
      const WString msg = WString::tr("psgu-fwhm-fit-fail")
                            .arg( static_cast<int>(inputPeak.size()) )
                            .arg( static_cast<int>(outputPeak.size()) );
      
      passMessage( msg, WarningWidget::WarningMsgInfo );
      return;
    }//if( outputPeak.size() != 1 )
    
    if( inputPeak.size() > 1 )
    {
      for( const shared_ptr<const PeakDef> &p : peaksNotInRoi )
        outputPeak.push_back( *p );
      std::sort( begin(outputPeak), end(outputPeak), &PeakDef::lessThanByMean );
      model->setPeaks( outputPeak );
    }else
    {
      assert( !outputPeak.empty() );
      
      model->updatePeak( peak, outputPeak[0] );
    }//if( inputPeak.size() > 1 )
  }catch( std::exception &e )
  {
    passMessage( "Sorry, error encountered refitting ROI with fixed FWHM.", WarningWidget::WarningMsgInfo );
    cerr << "Error encountered refitting ROI: " << e.what() << endl;
  }
}//void refit_peaks_with_drf_fwhm( InterSpec * const interspec, const double rightClickEnergy )

  
std::pair<std::unique_ptr<ReferenceLineInfo>,int> reference_line_near_peak( InterSpec * const interspec,
                                          const PeakDef &peak,
                                          const bool only_nuclide,
                                          const std::string &ref_line_hint )
{
  auto refLineTool = interspec->referenceLinesWidget();
  
  const double mean  = peak.gausPeak() ? peak.mean()
                                       : (((peak.mean() > peak.lowerX()) && (peak.mean() < peak.upperX()))
                                           ? peak.mean()
                                           : 0.5*(peak.lowerX() + peak.upperX()));
  const double sigma = peak.gausPeak() ? peak.sigma() :  0.25f*peak.roiWidth();
  const double lx = peak.lowerX(), ux = peak.upperX();
    
  double largest_w = -9999, best_energy = -1.0f;
  
  vector<ReferenceLineInfo> showingNucs;
  
  string parent = ref_line_hint;
  
  const RefLineDynamic * const dynamicRefLineTool = interspec->refLineDynamic();
  if( !parent.empty() && dynamicRefLineTool && dynamicRefLineTool->isActive() )
  {
    const size_t pos = parent.find(';');
    if( pos != string::npos )
      parent = parent.substr(0,pos);
    SpecUtils::trim(parent);
    
    const shared_ptr<vector<pair<double,ReferenceLineInfo>>> dynamic_ref_lines = dynamicRefLineTool->current_lines();
    if( dynamic_ref_lines )
    {
      for( size_t i = 0; showingNucs.empty() && (i < dynamic_ref_lines->size()); ++i )
      {
        const ReferenceLineInfo &ref_line = (*dynamic_ref_lines)[i].second;
        if( SpecUtils::iequals_ascii(ref_line.m_input.m_input_txt, parent) )
          showingNucs.push_back( ref_line );
      }//
    }//if( dynamic_ref_lines )
  }//if( check kinetic ref lines )
  
  if( refLineTool )
  {
    vector<ReferenceLineInfo> user_lines = refLineTool->showingNuclides();
    if( showingNucs.empty() )
      showingNucs = std::move(user_lines);
    else
      showingNucs.insert( end(showingNucs), begin(user_lines), end(user_lines) );
  }//if( refLineTool )
  
  size_t best_info_index = showingNucs.size();
  size_t best_line_index = 0;
  
  for( size_t info_index = 0; info_index < showingNucs.size(); ++info_index )
  {
    const ReferenceLineInfo &info = showingNucs[info_index];
    for( size_t line_index = 0; line_index < info.m_ref_lines.size(); ++line_index )
    {
      const ReferenceLineInfo::RefLine &l = info.m_ref_lines[line_index];
      
      if( only_nuclide
         && (!l.m_parent_nuclide
             || ((l.m_particle_type != ReferenceLineInfo::RefLine::Particle::Gamma)
                 && (l.m_particle_type != ReferenceLineInfo::RefLine::Particle::Xray) )) )
      {
        // Only consider lines who have a nuclide parent, and are a Gamma or Xray
        continue;
      }
      
      const double dist = fabs( l.m_energy - mean );
      const double w = l.m_normalized_intensity / (0.25*sigma + dist);
      
      if( (w > largest_w) && (l.m_energy < ux) && (l.m_energy > lx) && (dist < 5*sigma) )
      {
        largest_w = w;
        best_energy = l.m_energy;
        
        best_info_index = info_index;
        best_line_index = line_index;
      }//if( best candidate so far )
    }//for( const ReferenceLineInfo::RefLine &line : info.m_ref_lines )
    
    // If we have a source reference line info from the UI, and we found a matching candidate from that set
    //   of lines, we will use that.
    if( !parent.empty() && (best_info_index == info_index) && SpecUtils::iequals_ascii(info.m_input.m_input_txt, parent) )
      break;
  }//for( size_t info_index = 0; info_index < showingNucs.size(); ++info_index )
    
  if( (largest_w > 0) && (best_energy > 0) && (best_info_index < showingNucs.size()) )
  {
    const ReferenceLineInfo &best_info = showingNucs[best_info_index];
    assert( best_line_index < best_info.m_ref_lines.size() );
    return pair<unique_ptr<ReferenceLineInfo>,int>( make_unique<ReferenceLineInfo>(best_info),
                                                   static_cast<int>(best_line_index) );
  }
  
  return pair<unique_ptr<ReferenceLineInfo>,int>( nullptr, -1 );
}//reference_line_near_peak(...)
  

float source_or_reference_line_near_peak_energy( InterSpec * const interspec, const PeakDef &peak,
                                                const std::string &ref_line_hint )
{
  if( !interspec )
    return -999.9f;
  
  try
  {
    return peak.gammaParticleEnergy();
  }catch( std::exception & )
  {
  }
  
  const bool only_nucs = false;
  pair<unique_ptr<ReferenceLineInfo>,int> refline
                              = reference_line_near_peak( interspec, peak, only_nucs, ref_line_hint );
  
  if( !refline.first || refline.second < 0 )
    return -999.9f;
  
  const ReferenceLineInfo::RefLine &line = refline.first->m_ref_lines[refline.second];
  
  return static_cast<float>( line.m_energy );
}//float source_or_reference_line_near_peak_energy( InterSpec * const interspec, const PeakDef &peak )

  
  
tuple<const SandiaDecay::Nuclide *, double, float>
    nuclide_reference_line_near( InterSpec *viewer, const float energy, const std::string &hint_parent )
{
  double sigma = estimate_FWHM_of_foreground( energy );
  
  // For HPGe zoomed out, its really hard to be within the 5-peak-sigmas of the reference
  //  line that `reference_line_near_peak(...)` requires, so we'll also take into account
  //  the number of pixels a peak would show up as.
  //  Note: `SpectrumChartD3.prototype.updateMouseCoordText(...)` basically requires within 10px
  const shared_ptr<const SpecUtils::Measurement> hist
                = viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
      
  const bool isHPGe = PeakFitUtils::is_likely_high_res(viewer);
      
  if( isHPGe )
  {
    double xmin, xmax, ymin, ymax;
    viewer->displayedSpectrumRange( xmin, xmax, ymin, ymax );
    
    const int ww = viewer->renderedWidth();
    
    if( (xmin < xmax) && ((xmax - xmin) > 50.0) && (ww > 200) )
    {
      // Will assume chart plot-area is 100px narrower than total window (unchecked).
      const double keV_per_pixel = (xmax - xmin) / (ww - 100);
      // set peak-sigma to be 2 pixels, giving us to require 10 pixels of accuracy from the user
      sigma = std::max( sigma, 2.0*keV_per_pixel );
    }//if( x-range and render size seem reasonable )
  }//if( isHPGe )
      
  PeakDef tmppeak( energy, std::max(sigma,0.1), 100.0 );
  
  const pair<unique_ptr<ReferenceLineInfo>,int> refline
        = reference_line_near_peak( viewer, tmppeak, true, hint_parent );
  
  if( !refline.first || (refline.second < 0) )
    return tuple<const SandiaDecay::Nuclide *, double, float>( nullptr, 0.0, 0.0f );
  
  assert( refline.second < static_cast<int>(refline.first->m_ref_lines.size()) );
  ReferenceLineInfo::RefLine &line = refline.first->m_ref_lines[refline.second];
  
  const SandiaDecay::Nuclide *nuc = line.m_parent_nuclide;
  if( !nuc )
    return tuple<const SandiaDecay::Nuclide *, double, float>( nullptr, 0.0, 0.0f );
  
  
  double age = -1.0;
  const string &age_str = refline.first->m_input.m_age;
  if( !age_str.empty() )
  {
    try
    {
      age = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( age_str, nuc ? nuc->halfLife : -1.0 );
    }catch( std::exception )
    {
      assert( 0 );
    }
  }//if( !age_str.empty() )
  
  return tuple<const SandiaDecay::Nuclide *, double, float>(line.m_parent_nuclide, age, line.m_energy);
}//nuclide_reference_line_near(...)
  
  
float estimate_FWHM_of_foreground( const float energy )
{
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  if( !viewer )
    return 0.0f;
  
  PeakModel *peakmodel = viewer->peakModel();
  assert( peakmodel );
  if( !peakmodel )
    return 0.0f;
  
  const auto peaks = peakmodel->peakVec();
  const auto user_lb = std::lower_bound( begin(peaks), end(peaks), PeakDef(energy,1.0,1.0), &PeakDef::lessThanByMean );
  const auto user_ub = (user_lb == end(peaks)) ? end(peaks) : user_lb + 1;
  if( (user_lb != end(peaks)) && (user_ub != end(peaks)) && user_lb->gausPeak() && user_ub->gausPeak() )
  {
    //Linearly interpolate between peaks ... should probably upgrade to interpolating based on sqrt(energy) between them.
    const double lower_fwhm = user_lb->fwhm();
    const double upper_fwhm = user_ub->fwhm();
    const double lower_energy = user_lb->mean();
    const double upper_energy = user_ub->mean();
    const double fraction = (energy - lower_energy) / (upper_energy - lower_energy);
      
    return static_cast<float>( lower_fwhm + fraction*(upper_fwhm - lower_fwhm) );
  }
    
  std::shared_ptr<SpecMeas> specmeas = viewer->measurment(SpecUtils::SpectrumType::Foreground);
    
  std::shared_ptr<DetectorPeakResponse> drf = specmeas ? specmeas->detector() : nullptr;
  if( drf && drf->hasResolutionInfo() )
    return drf->peakResolutionFWHM(energy);
    
    
  // Check auto-fit peaks
  const set<int> &dispSamples = viewer->displayedSamples(SpecUtils::SpectrumType::Foreground);
  auto hintPeaks = specmeas->automatedSearchPeaks( dispSamples );
  if( hintPeaks && !hintPeaks->empty() )
  {
    SpecMeas::PeakDeque autopeaks( begin(*hintPeaks), end(*hintPeaks) );
    std::sort( begin(autopeaks), end(autopeaks), &PeakDef::lessThanByMeanShrdPtr ); //just to make sure
    
    auto dummy_peak = make_shared<const PeakDef>(energy,1.0,1.0);
    const auto hint_lb = std::lower_bound( begin(autopeaks), end(autopeaks), dummy_peak, &PeakDef::lessThanByMeanShrdPtr );
    const auto hint_ub = hint_lb==end(autopeaks) ? end(autopeaks) : hint_lb + 1;
    if( (hint_lb != end(autopeaks)) && (hint_ub != end(autopeaks)) && (*hint_lb)->gausPeak() && (*hint_ub)->gausPeak() )
    {
      //Linearly interpolate between peaks ... should probably upgrade to interpolating based on sqrt(energy) between them.
      const double lower_fwhm = (*hint_lb)->fwhm();
      const double upper_fwhm = (*hint_ub)->fwhm();
      const double lower_energy = (*hint_lb)->mean();
      const double upper_energy = (*hint_ub)->mean();
      const double fraction = (energy - lower_energy) / (upper_energy - lower_energy);
      
      return static_cast<float>( lower_fwhm + fraction*(upper_fwhm - lower_fwhm) );
    }
  }//if( hintPeaks && !hintPeaks->empty() )
    
    
  if( !peaks.empty() )
  {
    const PeakDef &refpeak = (user_lb!=end(peaks) ? *user_lb : (peaks.front().mean() > energy) ? peaks.front() : peaks.back());
    const double ref_width = refpeak.gausPeak() ? refpeak.fwhm() : 0.25f*refpeak.roiWidth();
    const double ref_energy = refpeak.gausPeak() ? refpeak.mean() : 0.5*(refpeak.upperX() + refpeak.lowerX());
    
    return static_cast<float>( ref_width*sqrt(energy/ref_energy) );
  }//if( !peaks.empty() )
  
  const shared_ptr<const SpecUtils::Measurement> meas
                              = viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  
  const bool isHPGe = PeakFitUtils::is_likely_high_res( viewer );
  
  try
  {
    const string datadir = InterSpec::staticDataDirectory();
    string drf_dir = SpecUtils::append_path(datadir, "GenericGadrasDetectors/HPGe 40%" );
    
    if( !isHPGe )
      drf_dir = SpecUtils::append_path(datadir, "GenericGadrasDetectors/NaI 1x1" );
    
    drf = make_shared<DetectorPeakResponse>();
    drf->fromGadrasDirectory( drf_dir );
    
    return drf->peakResolutionFWHM(energy);
  }catch(...)
  {
  }
  
  if( !isHPGe )
    return 2.35482f*17.5f*sqrt(energy/661.0f);
  return 2.35482f*0.67f*sqrt(energy/661.0f);
}//float estimate_FWHM_of_foreground( const float energy )
  
  
void refit_peak_with_photopeak_mean( InterSpec * const interspec,
                                    const double rightClickEnergy,
                                    const string &ref_line_hint )
{
  try
  {
    interspec->useMessageResourceBundle( "PeakSearchGuiUtils" );
    
    PeakModel * const model = interspec->peakModel();
    const shared_ptr<const SpecUtils::Measurement> data = interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    const shared_ptr<const SpecMeas> foreground = interspec->measurment( SpecUtils::SpectrumType::Foreground);
    
    shared_ptr<const PeakDef> peak = model->nearestPeak( rightClickEnergy );
    if( !peak )
    {
      passMessage( WString::tr("psgu-no-peak-near-right-click-msg"), WarningWidget::WarningMsgInfo );
      return;
    }
    
    const float energy = source_or_reference_line_near_peak_energy( interspec, *peak, ref_line_hint );
    if( energy < 10 )
    {
      passMessage( WString::tr("psgu-no-ref-line-near-peak-msg"), WarningWidget::WarningMsgInfo );
      return;
    }
  
    UndoRedoManager::PeakModelChange peak_undo_creator;
  
    vector<PeakDef> inputPeak;
    const vector<shared_ptr<const PeakDef>> peaksInRoi = model->peaksSharingRoi( peak );
    const vector<shared_ptr<const PeakDef>> peaksNotInRoi = model->peaksNotSharingRoi( peak );
  
    assert( peaksInRoi.size() >= 1 );
  
    for( const auto &m : peaksInRoi )
    {
      PeakDef p = *m;
    
      if( m == peak )
      {
        p.setMean( energy );
        p.setMeanUncert( 0.0 );
        p.setFitFor( PeakDef::CoefficientType::Mean, false );
        
        //Should we also assign the source gamma to be the same as the reference line?
      }//if( m == peak )
    
      inputPeak.push_back( p );
    }//for( const auto &m : peaksInRoi )

  
    std::sort( begin(inputPeak), end(inputPeak), &PeakDef::lessThanByMean );

    if( inputPeak.size() > 1 )
    {
      const shared_ptr<const DetectorPeakResponse> &detector = foreground->detector();
      WFlags<PeakFitLM::PeakFitLMOptions> fit_options;
      fit_options |= PeakFitLM::PeakFitLMOptions::SmallRefinementOnly; //not sure if this matters - or maybe its better to not have this option
      vector<shared_ptr<const PeakDef>> result = refitPeaksThatShareROI( data, detector, peaksInRoi, fit_options );

      if( result.size() == inputPeak.size() )
      {
        for( const shared_ptr<const PeakDef> &p : peaksNotInRoi )
          result.push_back( p );
        std::sort( begin(result), end(result), &PeakDef::lessThanByMeanShrdPtr );
        model->setPeaks( result );
        return;
      }else
      {
        cerr << "refit_peak_with_photopeak_mean was not successful" << endl;
      }//if( result.size() == inputPeak.size() ) / else
    }//if( inputPeak.size() > 1 )
  
  
    //  const double lowE = peak->mean() - 0.1;
    //  const double upE = peak->mean() + 0.1;
    const double lowE = inputPeak.front().mean() - 0.1;
    const double upE = inputPeak.back().mean() + 0.1;
    const double ncausalitysigma = 0.0;
    const double stat_threshold  = 0;
    const double hypothesis_threshold = 0;
  
    Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options;
    const bool isHPGe = PeakFitUtils::is_likely_high_res(interspec); //Only used if peak range is defined, so shouldnt be needed
    vector<PeakDef> outputPeak = fitPeaksInRange( lowE, upE, ncausalitysigma, stat_threshold,
                               hypothesis_threshold, inputPeak, data, fit_options, isHPGe );
    if( outputPeak.size() != inputPeak.size() )
    {
      WString msg = WString::tr("psgu-failed-refit-after-fixing-mean")
        .arg( static_cast<int>(inputPeak.size()) )
        .arg( static_cast<int>(outputPeak.size()) );
      passMessage( msg, WarningWidget::WarningMsgInfo );
      return;
    }//if( outputPeak.size() != 1 )
  
    if( inputPeak.size() > 1 )
    {
      for( const shared_ptr<const PeakDef> &p : peaksNotInRoi )
        outputPeak.push_back( *p );
      std::sort( begin(outputPeak), end(outputPeak), &PeakDef::lessThanByMean );
      model->setPeaks( outputPeak );
    }else
    {
      assert( !outputPeak.empty() );
    
      model->updatePeak( peak, outputPeak[0] );
    }//if( inputPeak.size() > 1 )
  }catch( std::exception &e )
  {
    passMessage( "Sorry, error encountered refitting ROI after fixing peak mean to"
                " reference photopeak.",
                WarningWidget::WarningMsgInfo );
    cerr << "Error encountered refitting ROI: " << e.what() << endl;
  }//try / catch
}//void refit_peak_with_photopeak_mean( InterSpec * const interspec, const double rightClickEnergy )
  
  
void change_continuum_type_from_right_click( InterSpec * const interspec,
                                            const double rightClickEnergy,
                                            const int continuum_type )
{
  try
  {
    assert( interspec );
    if( interspec )
      interspec->useMessageResourceBundle( "PeakSearchGuiUtils" );
    
    const PeakContinuum::OffsetType type = static_cast<PeakContinuum::OffsetType>( continuum_type );
    
    // Make sure a valid continuum type is passed in, although this check should never fail anyway
    bool valid_offset = false;
    switch( type )
    {
      case PeakContinuum::NoOffset:   case PeakContinuum::Constant:
      case PeakContinuum::Linear:     case PeakContinuum::Quadratic:
      case PeakContinuum::Cubic:      case PeakContinuum::FlatStep:
      case PeakContinuum::LinearStep: case PeakContinuum::BiLinearStep:
      case PeakContinuum::External:
        valid_offset = true;
        break;
    };//enum OffsetType
    
    assert( valid_offset );
    if( !valid_offset )
    {
      // Should never happen, but will leave check in for development
      passMessage( "Unexpected error in InterSpec::handleChangeContinuumTypeFromRightClick"
                   " - invalid continuum type - not proceeding", WarningWidget::WarningMsgHigh );
      return;
    }//if( !valid_offset )
    
    PeakModel * const model = interspec->peakModel();
    const shared_ptr<const SpecUtils::Measurement> data = interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    const shared_ptr<const SpecMeas> foreground = interspec->measurment( SpecUtils::SpectrumType::Foreground);
    
    if( !model || !data || !foreground ) //shouldnt ever happen
    {
      passMessage( "No data loaded to refit", WarningWidget::WarningMsgInfo );
      return;
    }
    
    
    std::shared_ptr<const PeakDef> peak = model->nearestPeak( rightClickEnergy );
    if( !peak || (rightClickEnergy < peak->lowerX()) || (rightClickEnergy > peak->upperX()) )
    {
      // Shouldnt happen, but jic
      passMessage( "There was no ROI at " + std::to_string(rightClickEnergy)
                  + " keV to modify continuum of", WarningWidget::WarningMsgInfo );
      return;
    }//if( !peak )
    
    
    if( peak->type() == PeakDef::DefintionType::DataDefined )
    {
      PeakDef newpeak = *peak;
      newpeak.continuum()->setType(type);
      model->updatePeak( peak, newpeak );
      return;
    }//if( peak->type() == PeakDef::DefintionType::DataDefined )
    
    
    const shared_ptr<const PeakContinuum> oldContinuum = peak->continuum();
    assert( oldContinuum );
    if( oldContinuum->type() == type )
    {
      passMessage( "Continuum is already of type " + WString::tr(PeakContinuum::offset_type_label_tr(type)).toUTF8(),
                  WarningWidget::WarningMsgInfo)
      return;
    }
    
    auto newContinuum = make_shared<PeakContinuum>( *oldContinuum );
    
    switch( type )
    {
      case PeakContinuum::NoOffset:   case PeakContinuum::Constant:
      case PeakContinuum::Linear:     case PeakContinuum::Quadratic:
      case PeakContinuum::Cubic:      case PeakContinuum::FlatStep:
      case PeakContinuum::LinearStep: case PeakContinuum::BiLinearStep:
        newContinuum->setType( type );
      break;
      
      case PeakContinuum::External:
      {
        newContinuum->setType( type );
        shared_ptr<SpecUtils::Measurement> continuumHist = estimateContinuum( data );
        newContinuum->setExternalContinuum( continuumHist );
        break;
      }
    }//case( type )
    
    
    vector<shared_ptr<const PeakDef>> oldPeaksInRoi = model->peaksSharingRoi(peak);
    vector<shared_ptr<const PeakDef>> newCandidatePeaks;
    for( const auto &p : oldPeaksInRoi )
    {
      auto newPeak = make_shared<PeakDef>( *p );
      newPeak->setContinuum( newContinuum );
      newCandidatePeaks.push_back( newPeak );
    }//for( const auto &p : oldPeaksInRoi )
    
    if( newCandidatePeaks.empty() ) //shouldnt happen, but jic
    {
      assert( 0 );
      throw std::runtime_error( "Somehow no input/starting peaks???" );
    }
    
    if( newCandidatePeaks.size() > 1 )
    {
      const shared_ptr<const DetectorPeakResponse> &detector = foreground->detector();
      const vector<shared_ptr<const PeakDef>> result
                                = refitPeaksThatShareROI( data, detector, newCandidatePeaks, WFlags<PeakFitLM::PeakFitLMOptions>(0) );
      
      if( result.size() == newCandidatePeaks.size() )
      {
        vector<PeakDef> newPeaks;
        for( const auto &p : result )
          newPeaks.push_back( *p );
        
        model->updatePeaks( oldPeaksInRoi, newPeaks );
      }else
      {
        WString msg = WString::tr("psgu-failed-change-cont-type")
          .arg( WString::tr(PeakContinuum::offset_type_label_tr(type)) )
          .arg( WString::tr((newCandidatePeaks.size() > 1)
                            ? "psgu-failed-change-cont-plural"
                            : "psgu-failed-change-cont-singular") );
        
        passMessage( msg, WarningWidget::WarningMsgInfo);
      }//if( result.size() == inputPeak.size() ) / else
      
      return;
    }else //if( inputPeak.size() > 1 )
    {
      vector<PeakDef> inputPeak, outputPeak;
      
      for( const auto &p : newCandidatePeaks )
        inputPeak.push_back( *p );
      
      const double lowE = peak->mean() - 0.1;
      const double upE = peak->mean() + 0.1;
      const double ncausalitysigma = 0.0;
      const double stat_threshold  = 0.0;
      const double hypothesis_threshold = 0.0;
      
      // Classifying this as a re-fit will keep the means and widths, from wandering too much.
      const bool isRefit = true;
      Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options;
      //fit_options |= PeakFitLM::PeakFitLMOptions::SmallRefinementOnly; Not sure what makes sense here, so we'll just go with default options
      
      const bool isHPGe = PeakFitUtils::is_likely_high_res(interspec); //Only used if peak range is defined, so shouldnt be needed
      outputPeak = fitPeaksInRange( lowE, upE, ncausalitysigma, stat_threshold,
                                   hypothesis_threshold, inputPeak, data, fit_options, isHPGe );

      if( outputPeak.empty() )
      {
        WString msg = WString::tr("psgu-failed-change-cont-insig-msg")
          .arg( WString::tr(PeakContinuum::offset_type_label_tr(type)) );
        
        passMessage( msg, WarningWidget::WarningMsgInfo );
        return;
      }//if( outputPeak.empty() )
      
      assert( outputPeak.size() == 1 );
      
      if( outputPeak.size() == 1 )
      {
        model->updatePeak( peak, outputPeak[0] );
      }else
      {
        assert( 0 );
        throw std::runtime_error( "Some how fitPeaksInRange didnt return exactly one peak" );
      }//if( inputPeak.size() == 1 ) / else
      
    }//if( inputPeak.size() > 1 ) / else
  }catch( std::exception &e )
  {
    cerr << "Unexpected error changing continuum type: " << e.what() << endl;
    passMessage( "Unexpected error changing continuum type." , WarningWidget::WarningMsgHigh );
  }
}//change_continuum_type_from_right_click(...)

  
void change_skew_type_from_right_click( InterSpec * const interspec,
                                         const double rightClickEnergy,
                                         const int skew_type )
{
  try
  {
    assert( interspec );
    if( interspec )
      interspec->useMessageResourceBundle( "PeakSearchGuiUtils" );
    
    const PeakDef::SkewType type = static_cast<PeakDef::SkewType>( skew_type );
    
    // Make sure a valid continuum type is passed in, although this check should never fail anyway
    bool valid_offset = false;
    switch( type )
    {
      case PeakDef::NoSkew:       case PeakDef::Bortel:
      case PeakDef::GaussExp:     case PeakDef::CrystalBall:
      case PeakDef::ExpGaussExp:  case PeakDef::DoubleSidedCrystalBall:
      case PeakDef::VoigtPlusBortel: case PeakDef::GaussPlusBortel:
      case PeakDef::DoubleBortel:
        valid_offset = true;
        break;
      
      case PeakDef::NumSkewType:
        valid_offset = false;
        break;
    };//enum OffsetType
    
    assert( valid_offset );
    if( !valid_offset )
    {
      // Should never happen, but will leave check in for development
      passMessage( "Unexpected error in InterSpec::handleChangeSkewTypeFromRightClick"
                   " - invalid skew type - not proceeding", WarningWidget::WarningMsgHigh );
      return;
    }//if( !valid_offset )
    
    PeakModel * const model = interspec->peakModel();
    const shared_ptr<const SpecUtils::Measurement> data = interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    const shared_ptr<const SpecMeas> foreground = interspec->measurment( SpecUtils::SpectrumType::Foreground);
    
    if( !model || !data || !foreground ) //shouldnt ever happen
    {
      passMessage( "No data loaded to refit", WarningWidget::WarningMsgInfo );
      return;
    }
    
    
    shared_ptr<const PeakDef> peak = model->nearestPeak( rightClickEnergy );
    if( !peak || (rightClickEnergy < peak->lowerX()) || (rightClickEnergy > peak->upperX()) )
    {
      // Shouldnt happen, but jic
      passMessage( "There was no ROI at " + std::to_string(rightClickEnergy)
                  + " keV to modify continuum of", WarningWidget::WarningMsgInfo );
      return;
    }//if( !peak )
    
    
    if( peak->skewType() == type )
    {
      passMessage( "Skew type is already of type " + string(PeakDef::to_string(type)),
                  WarningWidget::WarningMsgInfo)
      return;
    }//if( peak->skewType() == type )
    
    // Make sure all peaks in ROI are Gaussian defined
    vector<shared_ptr<const PeakDef>> peaks_in_roi = model->peaksSharingRoi( peak );
    for( const auto &p : peaks_in_roi )
    {
      if( peak->type() == PeakDef::DefintionType::DataDefined )
      {
        const WString msg = WString::tr("psgu-change-skew-not-gaus")
          .arg( std::to_string(rightClickEnergy) );
        
        passMessage( msg, WarningWidget::WarningMsgInfo );
        return;
      }//if( peak not Gaussian defined )
    }//for( const auto &p : peaks_in_roi )
    
    // `peak` should be in peaks_in_roi
    assert( std::find(begin(peaks_in_roi), end(peaks_in_roi), peak) != end(peaks_in_roi) );
    
    // Check if any peaks already have this skew type defined, and if so, start with those values
    shared_ptr<const PeakDef> near_skew_peak;
    shared_ptr<const deque<shared_ptr<const PeakDef>>> all_peaks = model->peaks();
    if( all_peaks )
    {
      vector<shared_ptr<const PeakDef>> peaks_by_nearest;
      peaks_by_nearest.insert( end(peaks_by_nearest), begin(*all_peaks), end(*all_peaks) );
      std::sort( begin(peaks_by_nearest), end(peaks_by_nearest),
                [peak]( const shared_ptr<const PeakDef> &lhs,
                       const shared_ptr<const PeakDef> &rhs) -> bool {
        return (fabs(lhs->mean() - peak->mean()) < fabs(rhs->mean() - peak->mean()));
      } );
      
      for( const auto &p : peaks_by_nearest )
      {
        if( p->skewType() == type )
        {
          near_skew_peak = p;
          break;
        }
      }//for( const auto &p : peaks_by_nearest )
    }//if( all_peaks )
    
    
    // Grab the suggested skew starting parameters, but use values from `near_skew_peak`, if valid
    const size_t num_skew_pars = PeakDef::num_skew_parameters( type );
    vector<double> skew_pars( num_skew_pars, 0.0 );
    vector<bool> fit_for_skew_pars( num_skew_pars, true );
    for( size_t skew_par_index = 0; skew_par_index < num_skew_pars; ++skew_par_index )
    {
      const auto ct = PeakDef::CoefficientType(PeakDef::CoefficientType::SkewPar0 + skew_par_index);
      double lower, upper, starting, step;
      const bool use = PeakDef::skew_parameter_range( type, ct, lower, upper, starting, step );
      assert( use );
      if( !use )
        throw std::logic_error("inconsistent skew par def");
      
      double val = starting;
      bool fit_for = true;
      if( near_skew_peak )
      {
        const double peak_val = near_skew_peak->coefficient( ct );
        if( (peak_val >= lower) && (peak_val <= upper) )
        {
          val = peak_val;
          fit_for = near_skew_peak->fitFor(ct);
        }
      }//if( near_skew_peak )
      
      skew_pars[skew_par_index] = val;
      fit_for_skew_pars[skew_par_index] = fit_for;
    }//for( loop over skew parameters )
    

    vector<shared_ptr<const PeakDef>> newCandidatePeaks;
    for( const auto &p : peaks_in_roi )
    {
      auto newPeak = make_shared<PeakDef>( *p );
      newPeak->setSkewType( type );
      for( size_t skew_par_index = 0; skew_par_index < num_skew_pars; ++skew_par_index )
      {
        const auto ct = PeakDef::CoefficientType(PeakDef::CoefficientType::SkewPar0 + skew_par_index);
        newPeak->set_coefficient( skew_pars[skew_par_index], ct );
        newPeak->set_uncertainty( 0.0, ct );
        newPeak->setFitFor( ct, fit_for_skew_pars[skew_par_index] );
      }
      
      newCandidatePeaks.push_back( newPeak );
    }//for( const auto &p : oldPeaksInRoi )
    
    if( newCandidatePeaks.empty() ) //shouldnt happen, but jic
    {
      assert( 0 );
      throw std::runtime_error( "Somehow no input/starting peaks???" );
    }
    
    if( newCandidatePeaks.size() > 1 )
    {
      const shared_ptr<const DetectorPeakResponse> &detector = foreground->detector();
      const vector<shared_ptr<const PeakDef>> result
                                = refitPeaksThatShareROI( data, detector, newCandidatePeaks, WFlags<PeakFitLM::PeakFitLMOptions>(0) );
      
      if( result.size() == newCandidatePeaks.size() )
      {
        vector<PeakDef> newPeaks;
        for( const auto &p : result )
          newPeaks.push_back( *p );
        
        model->updatePeaks( peaks_in_roi, newPeaks );
      }else
      {
        const WString msg = WString::tr("psgu-change-skew-make-insig")
          .arg( PeakDef::to_string(type) )
        .arg( WString::tr( (newCandidatePeaks.size() > 1)
                          ? "psgu-failed-change-cont-plural"
                          : "psgu-failed-change-cont-singular" ) );
        
        passMessage( msg, WarningWidget::WarningMsgInfo);
      }//if( result.size() == inputPeak.size() ) / else
      
      return;
    }else //if( inputPeak.size() > 1 )
    {
      vector<PeakDef> inputPeak, fixedPeaks, outputPeak;
      
      for( const auto &p : newCandidatePeaks )
        inputPeak.push_back( *p );
      
      const double lowE = peak->mean() - 0.1;
      const double upE = peak->mean() + 0.1;
      const double ncausalitysigma = 0.0;
      const double stat_threshold  = 0.0;
      const double hypothesis_threshold = 0.0;
      
      // We'll let the means and widths change significantly
      Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options;
      const bool isHPGe = PeakFitUtils::is_likely_high_res(interspec); //Only used if peak range is defined, so shouldnt be needed
      outputPeak = fitPeaksInRange( lowE, upE, ncausalitysigma, stat_threshold,
                                   hypothesis_threshold, inputPeak, data, fit_options, isHPGe );

      if( outputPeak.empty() )
      {
        const WString msg = WString::tr("psgu-failed-change-skew-insig-msg")
          .arg( PeakDef::to_string(type) );
        
        passMessage( msg, WarningWidget::WarningMsgInfo );
        return;
      }//if( outputPeak.empty() )
      
      assert( outputPeak.size() == 1 );
      
      if( outputPeak.size() == 1 )
      {
        model->updatePeak( peak, outputPeak[0] );
      }else
      {
        assert( 0 );
        throw std::runtime_error( "Some how fitPeaksInRange didnt return exactly one peak" );
      }//if( inputPeak.size() == 1 ) / else
      
    }//if( inputPeak.size() > 1 ) / else
  }catch( std::exception &e )
  {
    cerr << "Unexpected error changing skew type: " << e.what() << endl;
    passMessage( "Unexpected error changing skew type." , WarningWidget::WarningMsgHigh );
  }
}//change_skew_type_from_right_click(...)
  

void fit_template_peaks( InterSpec *interspec, std::shared_ptr<const SpecUtils::Measurement> data,
                         std::vector<PeakDef> input_peaks,
                         std::vector<PeakDef> orig_peaks,
                         const PeakTemplateFitSrc fitsrc, const string sessionid )
{
  if( !data )
  {
    cerr << "fit_template_peaks: data is invalid!" << endl; //prob shouldnt ever happen
    return;
  }
  
  unique_copy_continuum( input_peaks );
  
  vector<PeakDef> candidate_peaks, data_def_peaks;
  
  //The input continuum and peak area could be WAY off from current spectrum,
  // so estimate some sane starting values for them
  for( PeakDef peak : input_peaks )
  {
    if( !peak.gausPeak() )
    {
      auto externalcont = peak.continuum()->externalContinuum();
      
      //Update peak area for new data.  I dont think anything else needs updating
      //  (but not 100% sure ATM)
      double peak_area = data->gamma_integral( peak.lowerX(), peak.upperX() );
      if( externalcont )
        peak_area -= externalcont->gamma_integral( peak.lowerX(), peak.upperX() );
      peak.setPeakArea( peak_area );
      
      data_def_peaks.push_back( peak );
      continue;
    }//if( !peak.gausPeak() )
    
    
    //If the template peak is the same as an already existing peak - skip it
    bool already_has = false;
    const double centroid = peak.mean();
    for( const PeakDef &p : orig_peaks )
      already_has = (already_has || (fabs(p.mean() - centroid) < 1.5*p.sigma()));
    if( already_has )
      continue;
    
    const PeakContinuum::OffsetType type = peak.continuum()->type();
    peak.continuum()->calc_linear_continuum_eqn( data, centroid, peak.lowerX(), peak.upperX(), 2, 2 );
    peak.continuum()->setType( type );
    
    const double lowerx = std::max( peak.lowerX(), peak.mean() - 3*peak.sigma() );
    const double upperx = std::min( peak.upperX(), peak.mean() + 3*peak.sigma() );
    const double dataarea = data->gamma_integral( lowerx, upperx );
    const double contarea = peak.continuum()->offset_integral( lowerx, upperx, data );
    const double peakarea = (dataarea > contarea) ? (dataarea - contarea) : 10.0;
    peak.setPeakArea( peakarea );
    
    candidate_peaks.push_back( peak );
  }//for( auto &peak : input_peaks )
  
  bool isHPGe = false;
  
  {// Begin block to see if HPGe
    const vector<PeakDef> &peaks_for_test = orig_peaks.empty() ? input_peaks : orig_peaks;
    vector<shared_ptr<const PeakDef>> peaks;
    for( const PeakDef &p : peaks_for_test )
      peaks.push_back( make_shared<const PeakDef>( p ) );
    
    const PeakFitUtils::CoarseResolutionType type = PeakFitUtils::coarse_resolution_from_peaks( peaks );
    isHPGe = (type == PeakFitUtils::CoarseResolutionType::High);
  }// End block to see if HPGe
  
  
  Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options;
  //fit_options |= PeakFitLM::PeakFitLMOptions::MediumRefinementOnly; Not sure what makes sense here, so we'll just go with default options
  
  const double x0 = data->gamma_energy_min();
  const double x1 = data->gamma_energy_max();
  const double ncausalitysigma = 0.0;
  const double stat_threshold  = 0.0;
  const double hypothesis_threshold = 0.0;
  
  vector<PeakDef> fitpeaks = fitPeaksInRange( x0, x1, ncausalitysigma,
                                             stat_threshold, hypothesis_threshold,
                                             candidate_peaks, data, fit_options, isHPGe );
  
  //Add back in data_def_peaks and orig_peaks.
  fitpeaks.insert( end(fitpeaks), begin(data_def_peaks), end(data_def_peaks) );
  fitpeaks.insert( end(fitpeaks), begin(orig_peaks), end(orig_peaks) );
  
  //Make sure peaks are sorted, for good measure
  std::sort( begin(fitpeaks), end(fitpeaks), &PeakDef::lessThanByMean );
  
  
  //Dont show user the dialog if no peaks were found, and we were looking for peaks from previous
  //  spectrum
  if( fitpeaks.empty() && (fitsrc == PeakTemplateFitSrc::PreviousSpectrum) )
    return;
  
  WServer::instance()->post( sessionid, [=](){
    
    //Check if a new spectrum has been loaded while the peaks were being fit.
    //  If so, dont try to propogate the fit peaks.
    if( !wApp || (interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground) != data) )
      return;
    
    vector<ReferenceLineInfo> reflines;
    
    PeakSelectorWindowReason reason;
    switch( fitsrc )
    {
      case PeakTemplateFitSrc::CsvFile:
        reason = PeakSelectorWindowReason::PeaksFromCsvFile;
        new PeakSelectorWindow( interspec, reason, orig_peaks, data, fitpeaks, reflines );
        break;
        
      case PeakTemplateFitSrc::PreviousSpectrum:
        reason = PeakSelectorWindowReason::PeaksFromPreviousSpectrum;
        new PeakSelectorWindow( interspec, reason, input_peaks, data, fitpeaks, reflines );
        break;
    }//switch( fitsrc )

    wApp->triggerUpdate();
  } );
}//fit_template_peaks( ... )
  

void prepare_and_add_gadras_peaks( std::shared_ptr<const SpecUtils::Measurement> data,
                                  std::vector<PeakDef> gadras_peaks,
                                  std::vector<PeakDef> orig_peaks,
                                  const std::string sessionid )
{
  //cout << "prepare_and_add_gadras_peaks:" << endl;
  //for( const auto &p : gadras_peaks )
  //  cout << "\t" << p.mean() << ": " << p.amplitude() / data->live_time() << " cps" << endl;
  
  if( !data )
  {
    cerr << "prepare_and_add_gadras_peaks: data is invalid!" << endl; //prob shouldnt ever happen
    return;
  }
  
  unique_copy_continuum( gadras_peaks );
  
  vector<PeakDef> candidate_peaks, data_def_peaks;
  
  for( PeakDef peak : gadras_peaks )
  {
    if( !peak.gausPeak() || (peak.continuum()->type() != PeakContinuum::NoOffset) )
    {
      auto externalcont = peak.continuum()->externalContinuum();
      
      //Update peak area for new data.  I dont think anything else needs updating
      //  (but not 100% sure ATM)
      double peak_area = data->gamma_integral( peak.lowerX(), peak.upperX() );
      if( externalcont )
        peak_area -= externalcont->gamma_integral( peak.lowerX(), peak.upperX() );
      peak.setPeakArea( peak_area );
      
      data_def_peaks.push_back( peak );
      continue;
    }//if( !peak.gausPeak() )
    
    
    //If the template peak is the same as an already existing peak - skip it
    bool already_has = false;
    const double centroid = peak.mean();
    for( const PeakDef &p : orig_peaks )
      already_has = (already_has || (fabs(p.mean() - centroid) < 1.5*p.sigma()));
    if( already_has )
      continue;
    
    peak.continuum()->setType( PeakContinuum::Linear );
    peak.continuum()->calc_linear_continuum_eqn( data, centroid, peak.lowerX(), peak.upperX(), 2, 2 );
    
    peak.setFitFor( PeakDef::GaussAmplitude, false );
    peak.setFitFor( PeakDef::Mean, false );
    peak.setFitFor( PeakDef::Sigma, false );
    
    candidate_peaks.push_back( peak );
  }//for( auto &peak : input_peaks )
  
  sort( begin(candidate_peaks), end(candidate_peaks), PeakDef::lessThanByMean );
  
  // Combine peaks into a ROI, if they are close
  size_t npeaks_in_roi = 0;
  for( size_t i = 1; i < candidate_peaks.size(); ++i )
  {
    PeakDef &prev_peak = candidate_peaks[i-1];
    PeakDef &this_peak = candidate_peaks[i];
    
    const double prev_upper = prev_peak.mean() + 3.0*prev_peak.sigma();
    const double this_lower = this_peak.mean() - 3.0*this_peak.sigma();
    
    if( prev_upper < this_lower )
    {
      npeaks_in_roi = 1;
    }else
    {
      npeaks_in_roi += 1;
      
      const double prev_lower = prev_peak.lowerX();
      const double this_upper = this_peak.upperX();
      
      auto cont = prev_peak.continuum();
      assert( cont );
      this_peak.setContinuum( cont );
      cont->setRange( prev_lower, this_upper );
      cont->calc_linear_continuum_eqn( data, 0.5*(prev_lower + this_upper), prev_lower, this_upper, 2, 2 );
      
      if( npeaks_in_roi < 3 )
        cont->setType( PeakContinuum::Linear );
      else if( npeaks_in_roi < 4 )
        cont->setType( PeakContinuum::Quadratic );
      else
        cont->setType( PeakContinuum::Cubic );
    }//if( we should combine these peaks ) / else
  }//for( size_t i = 1; i < candidate_peaks.size(); ++i )
  
  bool isHPGe = false;
  
  {// Begin block to see if HPGe
    const vector<PeakDef> &peaks_for_test = orig_peaks.empty() ? gadras_peaks : orig_peaks;
    vector<shared_ptr<const PeakDef>> peaks;
    for( const PeakDef &p : peaks_for_test )
      peaks.push_back( make_shared<const PeakDef>( p ) );
    
    const PeakFitUtils::CoarseResolutionType type = PeakFitUtils::coarse_resolution_from_peaks( peaks );
    isHPGe = (type == PeakFitUtils::CoarseResolutionType::High);
  }// End block to see if HPGe
  
  const double x0 = data->gamma_energy_min();
  const double x1 = data->gamma_energy_max();
  const double ncausalitysigma = 3.0;
  const double stat_threshold  = 0.0;
  const double hypothesis_threshold = 0.0;
  Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options;
  
  vector<PeakDef> fitpeaks = fitPeaksInRange( x0, x1, ncausalitysigma,
                                             stat_threshold, hypothesis_threshold,
                                             candidate_peaks, data, fit_options, isHPGe );
  
  //Add back in data_def_peaks and orig_peaks.
  fitpeaks.insert( end(fitpeaks), begin(data_def_peaks), end(data_def_peaks) );
  fitpeaks.insert( end(fitpeaks), begin(orig_peaks), end(orig_peaks) );
  
  //Make sure peaks are sorted, for good measure
  std::sort( begin(fitpeaks), end(fitpeaks), &PeakDef::lessThanByMean );
  
  
  for( PeakDef &peak : fitpeaks )
  {
    peak.setFitFor( PeakDef::GaussAmplitude, true );
    peak.setFitFor( PeakDef::Mean, true );
    peak.setFitFor( PeakDef::Sigma, true );
  }//for( PeakDef &peak : fitpeaks )
  
  
  //Dont show user the dialog if no peaks were found, and we were looking for peaks from previous
  //  spectrum
  if( fitpeaks.empty() )
    return;
  
  WServer::instance()->post( sessionid, [=](){
    
    InterSpec *interspec = InterSpec::instance();
    assert( interspec );
    if( !interspec )
      return;
    
    //Check if a new spectrum has been loaded while the peaks were being fit.
    //  If so, dont try to propogate the fit peaks.
    if( !wApp || (interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground) != data) )
      return;
    
    vector<ReferenceLineInfo> reflines;
    
    const PeakSelectorWindowReason reason = PeakSelectorWindowReason::PeaksFromCsvFile; //PeaksFromPreviousSpectrum
    new PeakSelectorWindow( interspec, reason, orig_peaks, data, fitpeaks, reflines );
    
    wApp->triggerUpdate();
  } );
}//void prepare_and_add_gadras_peaks(...)


void add_peak_from_right_click( InterSpec * const interspec,
                               const double rightClickEnergy,
                               const std::string &ref_line_hint )
{
  if( !interspec )
    return;

  std::shared_ptr<const SpecUtils::Measurement> dataH = interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  std::shared_ptr<const PeakDef> peak = interspec->nearestPeak( rightClickEnergy );
  if( !peak
      || rightClickEnergy < peak->lowerX()
      || rightClickEnergy > peak->upperX()
      || !interspec->measurment(SpecUtils::SpectrumType::Foreground)
      || !dataH )
  {
    passMessage( WString::tr("err-no-roi-to-add-peak"), WarningWidget::WarningMsgInfo );
    return;
  }

  //get all the peaks that belong to this ROI
  typedef std::shared_ptr<const PeakContinuum> ContPtr;
  typedef map<ContPtr, std::vector<PeakDef > > ContToPeakMap;

  const auto origContinumm = peak->continuum();
  const std::shared_ptr<const deque<PeakModel::PeakShrdPtr>> allOrigPeaksDeque = interspec->peakModel()->peaks();
  
  // We need to make a copy of all the shared pointers because we modify the deque that
  //  allOrigPeaksDeque points at.
  const deque<PeakModel::PeakShrdPtr> allOrigPeaks( begin(*allOrigPeaksDeque), end(*allOrigPeaksDeque) );
  
  ContToPeakMap contToPeaks;
  for( const PeakModel::PeakShrdPtr &thispeak : allOrigPeaks )
  {
    if( thispeak )
      contToPeaks[thispeak->continuum()].push_back( *thispeak );
  }
  
  std::vector<PeakDef> origRoiPeaks = contToPeaks[peak->continuum()];
  vector<shared_ptr<const PeakDef>> origRoiPeaksPtr;
  
  //Make sure we dont try to mix data defined and gaussian defined peaks
  for( const PeakDef &p : origRoiPeaks )
  {
    if( !p.gausPeak() )
    {
      passMessage( WString::tr("err-no-add-peak-to-data-roi"), WarningWidget::WarningMsgInfo );
      return;
    }

    origRoiPeaksPtr.push_back( make_shared<const PeakDef>( p ) );
  }

  if( origRoiPeaks.empty() )
    throw runtime_error( "Logic error in add_peak_from_right_click()" );
  
  const double x0 = peak->lowerX();
  const double x1 = peak->upperX();
  const size_t lower_channel = dataH->find_gamma_channel( x0 );
  const size_t upper_channel = dataH->find_gamma_channel( x1 );
  const size_t nbin = ((upper_channel > lower_channel) ? (upper_channel-lower_channel) : size_t(0));
  
  const double startingChi2 = chi2_for_region( origRoiPeaksPtr, dataH, static_cast<int>(lower_channel), static_cast<int>(upper_channel) );


  bool inserted = false;
  vector< std::shared_ptr<PeakDef> > answer;
  for( PeakDef p : origRoiPeaks )
  {
    if( fabs(p.mean() - peak->mean()) < 0.01 )
    {
      inserted = true;
      PeakDef newpeak = p;
      newpeak.setMean( rightClickEnergy );
      newpeak.setSigma( newpeak.sigma() / sqrt(2.0) );
      newpeak.setAmplitude( 0.25*p.amplitude() );
      p.setAmplitude( 0.75*p.amplitude() );
      p.setSigma( p.sigma() / sqrt(2.0) );
      
      if( newpeak.mean() < p.mean() )
      {
        answer.push_back( std::make_shared<PeakDef>( newpeak ) );
        answer.push_back( std::make_shared<PeakDef>( p ) );
      }else
      {
        answer.push_back( std::make_shared<PeakDef>(p)  );
        answer.push_back( std::make_shared<PeakDef>(newpeak) );
      }
    }else
    {
      answer.push_back( std::make_shared<PeakDef>( p ) );
    }
  }

  if( !inserted )
    throw runtime_error( "Logic error 2 in add_peak_from_right_click()" );

  const bool isHPGe = PeakFitUtils::is_likely_high_res( interspec );
  
  const MultiPeakInitialGuessMethod methods[] = { FromInputPeaks, UniformInitialGuess, FromDataInitialGuess };
  
  bool found_better_fit = false;
  double fitChi2 = startingChi2;
  const vector<shared_ptr<PeakDef>> orig_answer = answer;
  for( const MultiPeakInitialGuessMethod method : methods )
  { 
    try
    {  
      vector<shared_ptr<PeakDef>> this_answer = orig_answer;
      
      double this_fitChi2;
      findPeaksInUserRange( x0, x1, int(answer.size()), method, dataH,
                           interspec->measurment(SpecUtils::SpectrumType::Foreground)->detector(), isHPGe, this_answer, this_fitChi2 );
      
      if( this_answer.empty() )
      {
        cerr << "Failed to fit " << answer.size() << " peaks in range [" << x0 << ", " << x1 << "] for method "
        << static_cast<int>(method) << " -- will skip." << endl;
        continue; //We failed to fit the peaks
      }
      
      vector<std::shared_ptr<const PeakDef> > this_answer_const;
      for( const std::shared_ptr<PeakDef> &p : this_answer )
        this_answer_const.push_back( p );

      this_fitChi2 = chi2_for_region( this_answer_const, dataH,
                                          static_cast<int>(lower_channel), static_cast<int>(upper_channel) );

      if( !IsNan(this_fitChi2) && !IsInf(this_fitChi2) && (this_fitChi2 <= fitChi2) )
      {
        fitChi2 = this_fitChi2;
        answer = this_answer;
        found_better_fit = true;
      }
    }catch( std::exception & )
    {
      //Do nothing
    }
  }//for( const MultiPeakInitialGuessMethod method : methods )

  const double dof = (nbin + 3*origRoiPeaks.size() + peak->continuum()->type());
  const double chi2Dof = fitChi2 / dof;
  
  cerr << "rightClickEnergy=" << rightClickEnergy << endl;
  cerr << "PreChi2=" << startingChi2 << ", PostChi2=" << fitChi2 << endl;
  cerr << "Got chi2Dof=" << chi2Dof << " when fitting for new peak - did " 
        << (found_better_fit ? "" : "not ") << "improve fit" << endl;
  
  for( size_t i = 0; i < answer.size(); ++i )
  {
    cerr << "Peak " << i << " at " << answer[i]->mean() << " w/ width="
         << answer[i]->sigma() << ", and amp=" << answer[i]->amplitude() << endl;
  }
  
  //Remove old ROI peaks
  if( !found_better_fit )
  {
    passMessage( WString::tr("err-add-peak-no-improve"), WarningWidget::WarningMsgInfo );
    return;
  }
  
  std::map<std::shared_ptr<PeakDef>,PeakModel::PeakShrdPtr> new_to_orig_peaks;
  for( const PeakModel::PeakShrdPtr &thispeak : allOrigPeaks )
  {
    if( !thispeak || (thispeak->continuum() != origContinumm) )
      continue;
    
    interspec->peakModel()->removePeak( thispeak );
    
    //Associate new peak 
    double dist = 99999.9;
    std::shared_ptr<PeakDef> nearest_new;
    for( std::shared_ptr<PeakDef> newpeak : answer )
    {
      const double d = fabs( newpeak->mean() - thispeak->mean() );
      if( d < dist && !new_to_orig_peaks.count(newpeak) )
      {
        dist = d;
        nearest_new = newpeak;
      }
    }
    if( nearest_new )
    {
      nearest_new->inheritUserSelectedOptions( *thispeak, true );
      new_to_orig_peaks[nearest_new] = thispeak;
    }
  }

  for( size_t i = 0; i < answer.size(); ++i )
  {
    const bool isNew = !new_to_orig_peaks.count(answer[i]);
    interspec->addPeak( *answer[i], isNew, SpecUtils::SpectrumType::Foreground, (isNew ? ref_line_hint : "") );
  }
}//void add_peak_from_right_click()

}//namespace PeakSearchGuiUtils
