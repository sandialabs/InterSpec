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
#include <Wt/WStandardItemModel>
#include <Wt/Chart/WCartesianChart>

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpec.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/PeakInfoDisplay.h"
#include "InterSpec/SpectrumDataModel.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"

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
    PeaksFromPreviousSpectum,
    
    /** User updloaded a peaks CSV file.
     */
    PeaksFromCsvFile,
  };//enum class PeakSelectorWindowReason
  
  
/** Class representing a dialog that allows the user to select which peaks to
 keep, or nuclides to assign to a peak after an automated peak fitting, nuclide
 ID, or peak propagation function.
 
 Currently this window can be triggered from four different user actions,
 leading to a base of four different ways to display things (there are also some
 dynamic display changes based on data and results), which adds some complexity,
 but maybe at the moment this is more manageable than breakign this class up...
 */
class PeakSelectorWindow : public AuxWindow
{
  InterSpec *m_viewer;
  Wt::WTable *m_table;
  int m_previewChartColumn;
  Wt::WCheckBox *m_keepRefLinePeaksOnly;
  Wt::WCheckBox *m_showAllPeaks;
  
  Wt::WPanel *m_chartPanel;
  Wt::Chart::WCartesianChart *m_chart;
  Wt::WStandardItemModel* m_relEffModel;
  
  const PeakSelectorWindowReason m_reason;
  const vector<PeakDef> m_orig_peaks;
  std::shared_ptr<const SpecUtils::Measurement> m_data;
  const vector<PeakDef> m_final_peaks;
  vector< shared_ptr<const ReferenceLineInfo>> m_displayed;
  vector<WCheckBox *> m_keep_peak_cbs;
  vector<WComboBox *> m_nuc_select_combos;
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
          For PeaksFromPreviousSpectum, these are the peaks that were fed into
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
               (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal) | AuxWindowProperties::TabletModal | AuxWindowProperties::DisableCollapse) ),
    m_viewer( viewer ),
    m_table( nullptr ),
    m_previewChartColumn( -1 ),
    m_keepRefLinePeaksOnly( nullptr ),
    m_showAllPeaks( nullptr ),
    m_chartPanel( nullptr ),
    m_chart( nullptr ),
    m_relEffModel( nullptr ),
    m_reason( reason ),
    m_orig_peaks( orig_peaks ),
    m_data( data ),
    m_final_peaks( final_peaks ),
    m_displayed(),
    m_cancelOperation( false )
  {
    wApp->useStyleSheet( "InterSpec_resources/PeakSelectorWindow.css" );
   
    switch( m_reason )
    {
      case PeakSelectorWindowReason::NuclideId:
        setWindowTitle( "Confirm Nuclide Assignment" );
        break;
      case PeakSelectorWindowReason::PeakSearch:
        setWindowTitle( "Check Peak Search Results" );
        break;
      case PeakSelectorWindowReason::PeaksFromPreviousSpectum:
        setWindowTitle( "Peaks to Keep from Previous Spectrum" );
        break;
      case PeakSelectorWindowReason::PeaksFromCsvFile:
        setWindowTitle( "Peaks to import from CSV" );
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
        const double newSigma = newPeak.sigma();
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
          txt = new WText( "Adjust nuclide assignments, or which peaks to keep.", contents() );
        else
          txt = new WText( "Unselect any peaks you would not like to keep.", contents() );
        break;
        
      case PeakSelectorWindowReason::NuclideId:
        txt = new WText( "Check and fix nuclide assignments.", contents() );
        break;
        
      case PeakSelectorWindowReason::PeaksFromPreviousSpectum:
      case PeakSelectorWindowReason::PeaksFromCsvFile:
        txt = new WText( "Select peaks peaks you would like to keep.", contents() );
        break;
    }//switch( m_reason )
    
   
    txt->setInline( false );
    
    switch( m_reason )
    {
      case PeakSelectorWindowReason::NuclideId:
        m_showAllPeaks = new WCheckBox( "Show all peaks", contents() );
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
            string msg = "Keep only new peaks assigned to ";
            for( size_t i = 0; i < m_displayed.size(); ++i )
            {
              if( i && (i+1)==m_displayed.size() )
                msg + " or ";
              else if( i )
                msg + ", ";
              msg += m_displayed[i]->labelTxt;
            }
            m_keepRefLinePeaksOnly = new WCheckBox( msg, contents() );
            m_keepRefLinePeaksOnly->setInline( false );
            m_keepRefLinePeaksOnly->changed().connect( this, &PeakSelectorWindow::keepOnlyRefLinesCbChanged );
            m_keepRefLinePeaksOnly->setMargin( 10, Wt::Top );
            m_keepRefLinePeaksOnly->addStyleClass( "KeepRefLinPeaksOnlyCb" );
          }//if( anyNonAssignedPeaks )
        }//if( we searched for peaks, and there were reference lines )
        
        if( orig_peaks.size() )
          m_showAllPeaks = new WCheckBox( "Show previous peaks too", contents() );
        break;
      }//case PeakSelectorWindowReason::PeakSearch:
        
      case PeakSelectorWindowReason::PeaksFromPreviousSpectum:
        if( final_peaks.size() != orig_peaks.size() )
          m_showAllPeaks = new WCheckBox( "Show absent peaks", contents() );
        break;
    }//switch ( m_reason )
    
    if( m_showAllPeaks )
    {
      m_showAllPeaks->addStyleClass( "ShowAllPeaksCb" );
      m_showAllPeaks->setInline( false );
      m_showAllPeaks->setChecked( false );
      m_showAllPeaks->changed().connect( this, &PeakSelectorWindow::showAllPeaksCbChanged );
      if( !m_keepRefLinePeaksOnly )
        m_showAllPeaks->setMargin( 10, Wt::Top );
    }//if( m_showAllPeaks )
    
    
    // The rel eff chart is close, but not fully debugged, so leaving diabled for now
    //setupRelEffChart();
    
    m_table = new WTable( contents() );
    m_table->addStyleClass( "PeakSelectorTable" );
    m_table->setHeaderCount( 1, Wt::Horizontal );
    
    m_keep_peak_cbs.resize( m_old_to_new_peaks.size(), nullptr );
    m_nuc_select_combos.resize( m_old_to_new_peaks.size(), nullptr );
    m_peak_updated.resize( m_old_to_new_peaks.size(), false );
    
    auto nucs_changed = [this]( size_t index ){
      assert( index < m_old_to_new_peaks.size() );
      
      switch( m_reason )
      {
        case PeakSelectorWindowReason::NuclideId:
        case PeakSelectorWindowReason::PeakSearch:
        case PeakSelectorWindowReason::PeaksFromCsvFile:
          break;
          
        case PeakSelectorWindowReason::PeaksFromPreviousSpectum:
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
          
        case PeakSelectorWindowReason::PeaksFromPreviousSpectum:
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
          
        case PeakSelectorWindowReason::PeaksFromPreviousSpectum:
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
      txt = new WText( "Keep Peak?", cell );
      txt->setWordWrap( false );
    }
    
    if( peakEnergyIndex >= 0 )
    {
      cell = m_table->elementAt(0,peakEnergyIndex);
      txt = new WText( "Peak Energy, FWHM", cell );
    }
    
    if( origColumnIndex >= 0 )
    {
      cell = m_table->elementAt(0,origColumnIndex);
      txt = new WText( "Original Nuclide", cell );
    }
    
    if( newColumnIndex >= 0 )
    {
      cell = m_table->elementAt(0,newColumnIndex);
      txt = new WText( "Assigned Nuclide", cell );
    }
    
    if( previewIndex >= 0 )
    {
      cell = m_table->elementAt(0,previewIndex);
      txt = new WText( "Peak Preview", cell );
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
        
        if( m_reason==PeakSelectorWindowReason::PeaksFromPreviousSpectum
            && !m_old_to_new_peaks[i].second )
        {
          new WText( "Not Found", cbcell );
        }else
        {
          WCheckBox *cb = new WCheckBox( "Keep Peak", cbcell );
          cb->setWordWrap(false);
          cb->setChecked(true);
          cb->changed().connect( this, &PeakSelectorWindow::keepPeakChanged );
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
            && m_reason!=PeakSelectorWindowReason::PeaksFromPreviousSpectum
            && m_reason!=PeakSelectorWindowReason::PeaksFromCsvFile )
        {
          m_nuc_select_combos[i] = new WComboBox(newNucCell);
          m_nuc_select_combos[i]->changed().connect( boost::bind( &PeakSelectorWindow::nucSelectChanged, this, i ) );
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
        snprintf( buffer, sizeof(buffer), "<div style=\"white-space: nowrap;\">"
                     "mean=%.1f keV</div><div style=\"white-space: nowrap;\">FWHM=%.1f keV</div>",
                     displayPeak->mean(), displayPeak->fwhm() );
        
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
      case PeakSelectorWindowReason::PeaksFromPreviousSpectum:
      case PeakSelectorWindowReason::PeaksFromCsvFile:
        break;
        
      case PeakSelectorWindowReason::NuclideId:
        if( !some_nuclides_changed )
        {
          txt = new WText( "<strong>No nuclides changed</strong>", contents() );
          txt->setPadding( 10, Wt::Top | Wt::Bottom );
          txt->setTextAlignment( Wt::AlignmentFlag::AlignCenter );
          txt->setInline( false );
        }
        break;
    }//switch( m_reason )

    
    keepPeakChanged();
    refreshRelEffChart();
    
    WPushButton *acceptButton = addCloseButtonToFooter( "Accept", true );
    acceptButton->clicked().connect( boost::bind( &AuxWindow::hide, this ) );
    
    WPushButton *cancelButton = nullptr;
    if( viewer->isPhone() )
    {
      cancelButton = new WPushButton( "Cancel", footer() );
      cancelButton->setFloatSide( Wt::Right );
    }else
    {
      cancelButton = addCloseButtonToFooter( "Cancel", true );
    }
    
    cancelButton->clicked().connect( boost::bind( &PeakSelectorWindow::cancelOperation, this ) );
    
    finished().connect( this, &PeakSelectorWindow::doFinish );
    
    show();
    
    resizeToFitOnScreen();
    centerWindow();
    centerWindowHeavyHanded();
  }//PeakSelector constructor
  
  /** The rel eff chart is close, but not fully debugged */
  void setupRelEffChart()
  {
    m_chartPanel = new WPanel( contents() );
    m_chartPanel->setCollapsible( true );
    m_chartPanel->setCollapsed( true );
    m_chartPanel->setTitle( "Relative Efficiency Plot" );
    m_chartPanel->setAnimation( WAnimation(WAnimation::SlideInFromTop, WAnimation::EaseOut, 100) );
    
    
    m_chart = new Wt::Chart::WCartesianChart();
    m_chart->setBackground(Wt::WColor(220, 220, 220));
    m_chart->setAutoLayoutEnabled();
    m_chart->setType(Wt::Chart::ScatterPlot);
    m_relEffModel = new WStandardItemModel( m_chart );
    m_chart->setModel( m_relEffModel );
    m_chart->setType( Chart::ScatterPlot );
    m_chart->setXSeriesColumn(0);
    
    
    //We should check the color theme for colors
    auto theme = m_viewer->getColorTheme();
    if( theme )
    {
      if( !theme->spectrumChartText.isDefault() )
      {
        WPen txtpen(theme->spectrumChartText);
        m_chart->setTextPen( txtpen );
        m_chart->axis(Chart::XAxis).setTextPen( txtpen );
        m_chart->axis(Chart::YAxis).setTextPen( txtpen );
        m_chart->axis(Chart::Y2Axis).setTextPen( txtpen );
      }
      
      if( theme->spectrumChartBackground.isDefault() )
        m_chart->setBackground( Wt::NoBrush );
      else
        m_chart->setBackground( WBrush(theme->spectrumChartBackground) );
      
      //From what I can tell, we cant change the legend text color easily, so
      //  we'll just cheat and back the legend background different enough from
      //  black so we can always read the text.  Lame, but whatever.
      m_chart->setLegendStyle( m_chart->legendFont(), m_chart->legendBorder(), WBrush(Wt::WColor(220, 220, 220, 120)) );
      
      if( (theme->spectrumChartMargins.isDefault() && !theme->spectrumChartBackground.isDefault()) )
      {
        //theme->spectrumChartBackground
      }else if( !theme->spectrumChartMargins.isDefault() )
      {
        //theme->spectrumChartMargins
      }
      
      if( !theme->spectrumAxisLines.isDefault() )
      {
        WPen defpen = m_chart->axis(Chart::XAxis).pen();
        defpen.setColor( theme->spectrumAxisLines );
        m_chart->axis(Chart::XAxis).setPen( defpen );
        m_chart->axis(Chart::Y1Axis).setPen( defpen );
        m_chart->axis(Chart::Y2Axis).setPen( defpen );
      }
    }//if( theme )
    
    m_chart->setPlotAreaPadding(5, Wt::Top);
    m_chart->setPlotAreaPadding(55, Wt::Bottom);
    m_chart->setPlotAreaPadding(55, Wt::Left);
    m_chart->setPlotAreaPadding(10, Wt::Right);
    
    m_chart->axis(Wt::Chart::XAxis).setTitle("Energy (keV)");
    m_chart->axis(Wt::Chart::Y1Axis).setTitle("Peak Area / BR");
    
#if( WT_VERSION >= 0x3030400 )
    m_chart->axis(Wt::Chart::Y1Axis).setTitleOrientation( Wt::Vertical );
#endif
    
    m_chart->setLegendEnabled( true );
    m_chart->setLegendLocation(Wt::Chart::LegendInside, Wt::Top, Wt::AlignRight);
    
    WContainerWidget *holder = new WContainerWidget();
    m_chartPanel->setCentralWidget( holder );
    WGridLayout *layout = new WGridLayout();
    holder->setLayout( layout );
    layout->addWidget( m_chart, 0, 0 );
    layout->setContentsMargins( 0, 0, 0, 0 );
    
    m_chart->setHeight( 250 );
  }//void setupRelEffChart()
  
  
  void refreshRelEffChart()
  {
    if( !m_chart )
      return;
    
    m_relEffModel->clear();
    
    set<double> energiesset;
    map<string,WColor> src_to_color;
    map<string, vector<pair<double,double>>> src_to_energy_eff;
    for( size_t i = 0; i < m_old_to_new_peaks.size(); ++i )
    {
      auto p = m_old_to_new_peaks[i].second;
      if( !p )
        p = m_old_to_new_peaks[i].first;
      
      //Only consider normal nuclides, and normal gammas (n xrays, or escape peaks)
      if( !p || !p->parentNuclide() || (p->sourceGammaType()!=PeakDef::NormalGamma) )
        continue;
      
      if( m_keep_peak_cbs[i] && !m_keep_peak_cbs[i]->isChecked() )
        continue;
      
      const double peakSigma = p->sigma();
      const SandiaDecay::Nuclide * const nuc = p->parentNuclide();
      const double photopeakEnergy = p->decayParticle()->energy;
      
      const string label = p->parentNuclide()->symbol;
      
      
      double intensity = 0.0;
      bool gotInfoFromRefLines = false;
      
      for( const auto &refline : m_displayed )
      {
        if( refline->nuclide == p->parentNuclide() )
        {
          double intensity = 0.0;
          for( size_t i = 0; i < refline->energies.size(); ++i )
          {
            const string &parttype = refline->particlestrs[i];
            
            if( !SpecUtils::icontains(parttype,"gamma")
               && !SpecUtils::icontains(parttype,"xray") )
              continue;
            
            const double refenergy = refline->energies[i];
            const double br = refline->intensities[i];
            if( fabs(refenergy-photopeakEnergy) < 1.25*peakSigma )
              intensity += br;
          }
          
          gotInfoFromRefLines = true;
          break;
        }
        
        if( gotInfoFromRefLines )
          break;
      }//for( const auto &refline : m_displayed )
      
      if( !gotInfoFromRefLines )
      {
        src_to_color[label] = p->lineColor();
        
        const double age = PeakDef::defaultDecayTime( nuc );
        SandiaDecay::NuclideMixture mix;
        mix.addAgedNuclideByActivity( nuc, 0.01*SandiaDecay::curie, age );
      
        for( const auto gamma : mix.gammas(0, SandiaDecay::NuclideMixture::OrderByEnergy, true ) )
        {
          if( fabs(gamma.energy-photopeakEnergy) < 1.25*peakSigma )
            intensity += gamma.numPerSecond;
        }
      
        if( intensity <= 0.0 )
          continue;
      }//if( !gotInfoFromRefLines )
      
      energiesset.insert( photopeakEnergy );
      
      if( intensity > 0.0 )
        src_to_energy_eff[label].emplace_back( photopeakEnergy, (p->peakArea()/intensity) );
      
    }//for( loop over peaks )
    
    vector<double> energies( begin(energiesset), end(energiesset) );
    
    m_relEffModel->insertColumn( 0 );
    m_relEffModel->insertRows(0, static_cast<int>(energies.size()) );
    
    for( size_t i = 0; i < energies.size(); ++i )
      m_relEffModel->setData( static_cast<int>(i), 0, energies[i] );
    
    m_chart->setXSeriesColumn(0);
    m_relEffModel->setHeaderData(0, WString("Energy (keV)") );
    
    const vector<Wt::Chart::MarkerType> markers{ Chart::SquareMarker,
      Chart::CircleMarker, Chart::CrossMarker, Chart::XCrossMarker,
      Chart::TriangleMarker,
    };
    
    for( auto &lp : src_to_energy_eff )
    {
      if( lp.second.size() < 2 )
        continue;
      
      double maxintensity = 0.0;
      for( const auto &p : lp.second )
        maxintensity = std::max( maxintensity, p.second );
      
      if( maxintensity <= 0.0 )
      {
        lp.second.clear();
        continue;
      }
      
      const int col = m_relEffModel->columnCount();
      m_relEffModel->insertColumn( col );
      m_relEffModel->setHeaderData( col, WString(lp.first) );
      
      for( auto &p : lp.second )
      {
        p.second /= maxintensity;
        auto pos = std::find( begin(energies), end(energies), p.first );
        if( pos != end(energies) )
          m_relEffModel->setData( pos-begin(energies), col, p.second );
      }
      
      Wt::Chart::WDataSeries series( col, Chart::PointSeries, Chart::Y1Axis );
      series.setMarkerBrush( WBrush(src_to_color[lp.first]) );
      series.setMarkerPen( WPen(src_to_color[lp.first]) );
      series.setMarker( markers[col%markers.size()] );
      m_chart->addSeries( series );
    }//for( auto &lp : src_to_energy_eff )
  }//void refreshRelEffChart()
  
  
  
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
  }//void showAllPeaksCbChanged()
  
  
  void updatePreviewPlot( size_t i )
  {
    if( i >= m_old_to_new_peaks.size() || m_previewChartColumn < 0 )
      return;
    
    
    auto finalpeak = m_old_to_new_peaks[i].second;
    auto originalpeak = m_old_to_new_peaks[i].first;
    auto displayPeak = finalpeak ? finalpeak : originalpeak;
    
    if( !displayPeak )
      return;
    
    auto peaks_for_plotting = make_shared< std::deque<std::shared_ptr<const PeakDef> > >();
    for( const auto &p : m_old_to_new_peaks )
      if( p.second )
        peaks_for_plotting->push_back( make_shared<PeakDef>(*p.second) );
    
    std::shared_ptr<const ColorTheme> theme = m_viewer->getColorTheme();
    
    const double lowerx = displayPeak->mean() - 0.75*displayPeak->roiWidth();
    const double upperx = displayPeak->mean() + 0.75*displayPeak->roiWidth();
    
    //TODO: Find 'displayPeak' in 'peaks_for_plotting' by energy, and check
    //  if it shares a ROI with other peaks, and if so, change the other
    //  peak colors, and maybe the displayed width.
    
    
    std::shared_ptr<WSvgImage> svg = PeakSearchGuiUtils::renderChartToSvg( m_data, peaks_for_plotting, m_displayed, lowerx, upperx, 225, 125, theme, true );
    if( svg )
    {
      stringstream strm;
      svg->write( strm );
      WTableCell *preview = m_table->elementAt( i+1, m_previewChartColumn );
      preview->clear();
      if( !preview->hasStyleClass("PeakPreviewCell") )
        preview->addStyleClass( "PeakPreviewCell" );
      new WText( strm.str(), Wt::XHTMLUnsafeText, preview );
    }
  }//void updatePreviewPlot( size_t i )
  
  
  void keepPeakChanged()
  {
    if( !m_keepRefLinePeaksOnly )
      return;
    
    bool anyNonAssignedPeaks = false;
    for( size_t i = 0; !anyNonAssignedPeaks && i < m_old_to_new_peaks.size(); ++i )
    {
      if( m_table->rowAt( static_cast<int>(i+1) )->isHidden() )
        continue;
      
      if( m_old_to_new_peaks[i].second && m_keep_peak_cbs[i] && m_keep_peak_cbs[i]->isChecked() )
        anyNonAssignedPeaks = !m_old_to_new_peaks[i].second->hasSourceGammaAssigned();
    }
    
    m_keepRefLinePeaksOnly->setChecked( !anyNonAssignedPeaks );
    refreshRelEffChart();
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
    
    refreshRelEffChart();
  }//void keepOnlyRefLinesCbChanged()
  
  
  
  std::string nuclideDesc( const SandiaDecay::Nuclide *nuc,
                           const SandiaDecay::RadParticle *particle,
                           const PeakDef::SourceGammaType type )
  {
    if( !nuc || !particle )
      return "---";
    
    const ReferenceLineInfo *refline = nullptr;
    const BackgroundLine *backgroundLine = nullptr;
    
    const char *gamtype = "";
    switch( type )
    {
      case PeakDef::NormalGamma: case PeakDef::XrayGamma:    break;
      case PeakDef::AnnihilationGamma: gamtype = " Annih."; break;
      case PeakDef::DoubleEscapeGamma: gamtype = " D.E.";   break;
      case PeakDef::SingleEscapeGamma: gamtype = " S.E.";   break;
    }//switch( displayPeak->sourceGammaType() )
    
    
    for( const auto &d : m_displayed )
    {
      if( nuc == d->nuclide )
      {
        refline = d.get();
        break;
      }else
      {
        for( const auto &b : d->backgroundLines )
        {
          bool matches = true;
          const BackgroundLineType type = get<3>(*b);
          switch( type )
          {
            case U238Series:    matches = nuc->symbol == "U238"; break;
            case U235Series:    matches = nuc->symbol == "U235"; break;
            case Th232Series:   matches = nuc->symbol == "Th232"; break;
            case Ra226Series:   matches = nuc->symbol == "Ra226"; break;
            case K40Background: matches = nuc->symbol == "K40"; break;
            case BackgroundXRay:
            case BackgroundReaction:
            case OtherBackground:
              break;
          }//switch( type )
          
          if( matches )
          {
            const float energy = get<0>(*b);
            if( fabs(energy - particle->energy) < 0.1 )
            {
              refline = d.get();
              backgroundLine = b;
              break;
            }
          }
        }//for( const auto &b : d.backgroundLines )
      }//if ( nuc matches ) / else
    }//for( const auto &d : m_displayed )
    
    char buffer[256] = { '\0' };
    if( refline )
    {
      double intensity = 0.0;
      for( size_t i = 0; (i < refline->intensities.size()) && (i < refline->energies.size()) && (i < refline->particlestrs.size()); ++i )
      {
        if( fabs(refline->energies[i] - particle->energy) >= 0.01 )
          continue;
        
        const string &parttype = refline->particlestrs[i];
        if( !SpecUtils::icontains(parttype,"gamma")
           && !SpecUtils::icontains(parttype,"xray") )
          continue;
        
        double sf = 1.0;
        auto sf_iter = refline->particle_sf.find(parttype);
        if( sf_iter != std::end(refline->particle_sf) )
        {
          sf = sf_iter->second;
        }else
        {
          const string msg = "PeakSelectorWindow::nuclideDesc(): What - this shouldnt happen! - coudlnt find SF for '" + parttype + "'";
#if( PERFORM_DEVELOPER_CHECKS )
          log_developer_error( __func__, msg.c_str() );
#else
          cerr << msg << endl;
#endif
        }
        
        intensity += refline->intensities[i] * sf;
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
    if( i >= m_nuc_select_combos.size() )  //shouldnt ever happen, but JIC
      return;
    
    //shared_ptr<PeakDef> &oldpeak = m_old_to_new_peaks[i].first;
    shared_ptr<PeakDef> &newpeak = m_old_to_new_peaks[i].second;
    
    if( !newpeak )  //shouldnt ever happen!
      return;
    
    WComboBox *combo = m_nuc_select_combos[i];
    if( !combo ) //shouldnt ever happen!
      return;
    
    const int index = combo->currentIndex();
    
    if( index <= 0 )
    {
      newpeak->clearSources();
      
      if( m_viewer->colorPeaksBasedOnReferenceLines() )
      {
        for( const auto &l : m_displayed )
        {
          if( l->lineColor == newpeak->lineColor() )
            newpeak->setLineColor( WColor() );
        }//
      }//if( m_viewer->colorPeaksBasedOnReferenceLines() )
      
      updatePreviewPlot(i);
      keepPeakChanged();
      refreshRelEffChart();
      return;
    }//if( no source )
    
    const string newnucstr = combo->currentText().toUTF8();
    
    const auto result = PeakModel::setNuclideXrayReaction( *newpeak, newnucstr, 0.0 );
    
    switch( result )
    {
      case PeakModel::SetGammaSource::NoSourceChange:
        //Shouldnt ever happen, but JIC...
        passMessage( "Trouble making change to nuclide - not applying, sorry!<br />"
                     "Please report error, including selected nuclides to interspec@sandia.gov", "", 2 );
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
              if( !newpeak->parentNuclide() || l->backgroundLines.empty() )
                return false;
              for( const auto b : l->backgroundLines )
              {
                string refnuc;
                const BackgroundLineType type = std::get<3>(*b);
                switch( type )
                {
                  case U238Series:         refnuc = "U238";  break;
                  case U235Series:         refnuc = "U235";  break;
                  case Th232Series:        refnuc = "Th232"; break;
                  case Ra226Series:        refnuc = "Ra226"; break;
                  case K40Background:      refnuc = "K40";   break;
                  case BackgroundXRay:
                  case BackgroundReaction:
                  case OtherBackground:    refnuc = std::get<2>(*b); break;
                }//switch( type )
              
                cout << refnuc << endl;
                if( db->nuclide(refnuc) == newpeak->parentNuclide() )
                  return true;
              }
              return false;
            };//isBackgroundNuc lambda
            
            
            if( m_viewer->colorPeaksBasedOnReferenceLines()
               && ((l->nuclide && (l->nuclide==newpeak->parentNuclide()))
               || (l->element && (l->element==newpeak->xrayElement()))
               || (newpeak->reaction() && l->reactionsTxt.size() && SpecUtils::icontains(l->reactionsTxt, newpeak->reaction()->name()))
               || isBackgroundNuc()) )
            {
              newpeak->setLineColor( l->lineColor );
            }
          }//for( const auto &l : m_displayed )
        }//if( we should set peak color based on reference lines )
        break;
      }//case: source change successful
    }//switch( result )
    
    updatePreviewPlot(i);
    keepPeakChanged();
    refreshRelEffChart();
  }//void nucSelectChanged( const size_t i )
  
  
  void populateNuclideSelects()
  {
    assert( m_old_to_new_peaks.size() == m_nuc_select_combos.size() );
    
    for( size_t i = 0; i < m_old_to_new_peaks.size(); ++i )
    {
      if( !m_nuc_select_combos[i] )
        continue;
      
      m_nuc_select_combos[i]->clear();
      m_nuc_select_combos[i]->addItem( "No Source" );
      
      shared_ptr<PeakDef> &p = m_old_to_new_peaks[i].second;
      if( !p )
        continue;
      shared_ptr<PeakDef> &oldpeak = m_old_to_new_peaks[i].first;
      
      const double peakmean = p->mean();
      double energy = peakmean;
      PeakDef::SourceGammaType gammatype = PeakDef::SourceGammaType::NormalGamma;
      
      if( p->xrayElement() || p->reaction() )
      {
        energy = p->gammaParticleEnergy();
      }else if( p->parentNuclide() && p->decayParticle() )
      {
        gammatype = p->sourceGammaType();
        energy = p->decayParticle()->energy;
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
          
        case PeakSelectorWindowReason::PeaksFromPreviousSpectum:
          break;
          
        case PeakSelectorWindowReason::PeaksFromCsvFile:
          break;
      }//switch( m_reason )

      
      //We'll sort the select by distance away from current nuclide, or peak mean
      //  - will probably miss some weird edgecase, but what ever for the moment.
      map<double,set<string>> distToDescs;
      for( const auto &ref : m_displayed )
      {
        for( size_t i = 0; i < ref->energies.size(); ++i )
        {
          const string &parttype = ref->particlestrs[i];
          //SandiaDecay::to_str(SandiaDecay::BetaParticle)
          
          if( !SpecUtils::icontains(parttype,"gamma")
             && !SpecUtils::icontains(parttype,"xray") )
            continue;
          
          double sf = 1.0;
          auto sf_iter = ref->particle_sf.find(parttype);
          if( sf_iter != std::end(ref->particle_sf) )
          {
            sf = sf_iter->second;
          }else
          {
            const string msg = "PeakSelectorWindow::populateNuclideSelects(): What - this shouldnt happen! - coudlnt find SF for '" + parttype + "'";
#if( PERFORM_DEVELOPER_CHECKS )
            log_developer_error( __func__, msg.c_str() );
#else
            cerr << msg << endl;
#endif
          }
          
          string label = ref->labelTxt;
          const double refenergy = ref->energies[i];
          const double intensity = ref->intensities[i] * sf;
          
          bool backgroundLineMatchedRef = false;
          
          if( !ref->backgroundLines.empty() )
          {
            //Find the background line that cooresponds to `refenergy`
            for( const auto b : ref->backgroundLines )
            {
              const float lineenergy = std::get<0>(*b);
              const string &srcstr = std::get<2>(*b);
              const BackgroundLineType type = std::get<3>(*b);
              
              if( fabs(lineenergy - refenergy) < 0.001 )
              {
                switch( type )
                {
                  case BackgroundLineType::U238Series:         label = "U238";  break;
                  case BackgroundLineType::U235Series:         label = "U235";  break;
                  case BackgroundLineType::Th232Series:        label = "Th232"; break;
                  case BackgroundLineType::Ra226Series:        label = "Ra226"; break;
                  case BackgroundLineType::K40Background:      label = "K40";   break;
                  case BackgroundLineType::BackgroundXRay:
                  case BackgroundLineType::BackgroundReaction:
                  case BackgroundLineType::OtherBackground:    label = srcstr;  break;
                }//switch( type )
                
                backgroundLineMatchedRef = true;
                break;
              }//if( we found the background line with the correct energy )
            }//for( const auto b : ref.backgroundLines )
            
          }//if( !ref.backgroundLines.empty() )
          
          
          char buffer[128];
          snprintf( buffer, sizeof(buffer), "%s %.1f keV, I=%.2g%%",
                    label.c_str(), refenergy, 100.0*intensity );
          
          const double diff_from_mean = fabs( refenergy - peakmean );
          const double diff_from_assigned = fabs( refenergy - energy );
          
          if( (diff_from_assigned < 0.001)
             && ((p->parentNuclide() && (ref->nuclide==p->parentNuclide()))
                 || (p->xrayElement() && (ref->element==p->xrayElement()))
                 || (p->reaction() && ref->reactionsTxt.size() && SpecUtils::icontains(ref->reactionsTxt, p->reaction()->name()))
             ))
          {
#if( PERFORM_DEVELOPER_CHECKS )
            if( currentstr.size() )  //shouldnt ever happen, right?
              log_developer_error( __func__, ("PeakSelectorWindow::populateNuclideSelects(): Found mutliple matches of select strings for a nuclide: '" + string(buffer) + "'").c_str() );
#endif
            currentstr = buffer;
          }else if( backgroundLineMatchedRef && (diff_from_assigned < 0.1) )
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
      if( m_nuc_select_combos[i] && (m_nuc_select_combos[i]->currentIndex()<=0)
          && m_old_to_new_peaks[i].second && !m_table->rowAt(static_cast<int>(i+1))->isHidden() )
      {
        m_old_to_new_peaks[i].second->clearSources();
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
          
          case PeakSelectorWindowReason::PeaksFromPreviousSpectum:
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

  
  
/** Assigns the passed in peak to a reference line if reasonable.  If this
 assignment would cause one of the already existing peaks to change assignment
 (most likely if two peaks of diff amp, but same nuc, are right next to each
 other) then the peak that should be changed, along with the string of the
 assignement it should be changed to.
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
    
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  //There is a fairly common situation (especially for HPGe) where there is a
  //  small peak, next to a much larger peak, where if the user first
  //  identifies the large peak, the correct gamma-ray association gets made,
  //  but then whe the second one is identified, it also gets assigned the
  //  same gamma-ray association as the larger peak, which is incoorect.  To
  //  avoid this, we will use previouspeaks and prevpeak to check, and correct
  //  for this condition.
  
  double prevPeakDist = DBL_MAX, prevIntensity = DBL_MAX;
  double thisIntensity = DBL_MAX;
  
  PeakModel::PeakShrdPtr prevpeak;
  
  
  try
  {
    //const bool isHPGe = (data && data->num_gamma_channels() > 2048);
    
    double mindist = 99999999.9;
    double nearestEnergy = -999.9;
    
    double minx(0.0), maxx(0.0);
    findROIEnergyLimits( minx, maxx, peak, data );
    
    const SandiaDecay::Nuclide *nuclide = NULL;
    const SandiaDecay::Element *element = NULL;
    const ReactionGamma::Reaction *reaction = NULL;
    
    const double mean = peak.mean();
    const double sigma = peak.gausPeak() ? peak.sigma() : peak.roiWidth();
    
    Wt::WColor color;
    
    for( const ReferenceLineInfo &nuc : displayed )
    {
      for( const ReactionGamma::ReactionPhotopeak &rpp : nuc.reactionGammas )
      {
        if( rpp.abundance <= 0.0 )
          continue;
        
        const double delta_e = fabs( mean - rpp.energy );
        const double dist = (0.25*sigma + delta_e) / rpp.abundance;
        
        /*
         //ToDo: implement checking if the S.E. or D.E. of this line is a better
         //      fit than this line, and if so, use it, and assign the peak
         //      as a S.E. or D.E.
        bool is_SE = false, is_DE = false;
        if( isHPGe && showingEscapePeakFeature && rpp.energy > 1122.0f )
        {
          //ToDo: Right now we are assuming the escape peaks will have the same
          //      amplitudes as the full energy peaks, which isnt correct.
          
          const double se_delta_e = fabs( mean + 510.9989 - rpp.energy );
          const double se_dist = (0.25*sigma + se_delta_e) / rpp.abundance;
          
          const double de_delta_e = fabs( mean + (2*510.9989) - rpp.energy );
          const double de_dist = (0.25*sigma + de_delta_e) / rpp.abundance;
          
          if( se_dist < dist )
          {
            is_SE = true;
            dist = se_dist;
          }
        }//if( showingEscapePeakFeature )
        */
        
        if( dist < mindist && rpp.energy >= minx && rpp.energy <= maxx )
        {
          //should check to see if this energy is already assigned to another
          //  peak, if it is, ideally we would want to see which one is most
          //  compatible and swap
          
          bool currentlyused = false;
          for( const PeakModel::PeakShrdPtr &pp : *previouspeaks )
          {
            if( pp->reaction() && (pp->reaction() == rpp.reaction)
                && (fabs(pp->reactionEnergy() - rpp.energy) < 1.0) )
            {
              currentlyused = true;
              if( dist < prevPeakDist )
              {
                prevpeak = pp;
                prevPeakDist = dist;
                prevIntensity = rpp.abundance;
              }
              break;
            }//if( pp->reaction() == rpp.reaction )
          }//for( const PeakModel::PeakShrdPtr &pp : *previouspeaks )
          
          if( !currentlyused )
          {
            prevpeak.reset();
            mindist = dist;
            nuclide = NULL;
            element = NULL;
            color = nuc.lineColor;
            reaction = rpp.reaction;
            nearestEnergy = rpp.energy;
            thisIntensity = rpp.abundance;
          }
        }//if( we should possible associate this peak with this line )
      }//for( const ReactionGamma::ReactionPhotopeak &rpp : nuc.reactionGammas )
      
      if( nuc.nuclide || nuc.element )
      {
        for( size_t i = 0; i < nuc.energies.size(); ++i )
        {
          const double energy = nuc.energies[i];
          double intensity = nuc.intensities.at(i);
          
          bool isXray = !nuc.nuclide;
          if( nuc.nuclide && (energy < 120.0*SandiaDecay::keV) )
          {
            const SandiaDecay::Element *el = db->element( nuc.nuclide->atomicNumber );
            for( const SandiaDecay::EnergyIntensityPair &e : el->xrays )
              isXray |= (e.energy==energy);
          }//if( this could possibly be an xray )
          
          if( isXray )
            intensity *= 0.1;
          
          const double delta_e = fabs( mean - energy );
          const double dist = (0.25*sigma + delta_e) / intensity;
          
          
          if( dist < mindist && energy >= minx && energy <= maxx )
          {
            bool currentlyused = false;
            for( const PeakModel::PeakShrdPtr &pp : *previouspeaks )
            {
              const bool sampleNuc = (nuc.nuclide
                                      && pp->decayParticle()
                                      && (pp->parentNuclide() == nuc.nuclide)
                                      && (fabs(pp->decayParticle()->energy - energy) < 0.01));
              const bool sameXray = (isXray
                                     && nuc.element
                                     && (pp->xrayElement() == nuc.element)
                                     &&  (fabs(pp->xrayEnergy() - energy) < 0.01));
              
              if( sampleNuc || sameXray )
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
              color = nuc.lineColor;
              nuclide = nuc.nuclide;
              element = nuc.element;
              reaction = NULL;
              nearestEnergy = energy;
              thisIntensity = intensity;
              
              if( isXray )
              {
                nuclide = NULL;
                if( !element && nuc.nuclide )
                  element = db->element( nuc.nuclide->atomicNumber );
              }//if( isXray )
            }//if( !currentlyused )
            
          }//if( we should possible associate this peak with this line )
        }//for( const double energy : nuc.energies )
      }//if( nuc.nuclide )
      
      for( const BackgroundLine *line : nuc.backgroundLines )
      {
        const float energy = std::get<0>(*line);
        const float intensity = std::get<1>(*line);
        const string &isotope = std::get<2>(*line);
        const BackgroundLineType type = std::get<3>(*line);
        
        if( intensity < FLT_MIN )
          continue;
        
        const double delta_e = fabs( mean - energy );
        const double dist = (0.25*sigma + delta_e) / intensity;
        
        if( dist < mindist && energy >= minx && energy <= maxx )
        {
          const SandiaDecay::Nuclide *thisnuclide = NULL;
          const SandiaDecay::Element *thiselement = NULL;
          const ReactionGamma::Reaction *thisreaction = NULL;
          
          switch( type )
          {
            case U238Series:      thisnuclide = db->nuclide( "U238" );  break;
            case U235Series:      thisnuclide = db->nuclide( "U235" );  break;
            case Th232Series:     thisnuclide = db->nuclide( "Th232" ); break;
            case Ra226Series:     thisnuclide = db->nuclide( "Ra226" ); break;
            case K40Background:   thisnuclide = db->nuclide( "K40" );   break;
            case OtherBackground:
            {
              //ToDo: deal with single and double escapes, both for backgrounds
              //      and normal sources.
              //auto pos = isotope.find( "S.E." );
              //if( pos == string::npos )
              //  pos = isotope.find( "D.E." );
              //if( pos != string::npos )
              //{
              //  string iso = isotope.substr(0,pos);
              //  SpecUtils::trim( iso );
              //  thisnuclide = db->nuclide( iso );
              //}else
              //{
                thisnuclide = db->nuclide( isotope );
              //}
              break;
            }
            case BackgroundXRay:
              thiselement = db->element( isotope.substr(0,isotope.find(' ')) );
              break;
              
            case BackgroundReaction:
              break;
          }//switch( type )
          
          if( !thisnuclide && !thiselement )
            continue;
          
          bool currentlyused = false;
          for( const PeakModel::PeakShrdPtr &pp : *previouspeaks )
          {
            if( ((thisnuclide && (pp->parentNuclide() == thisnuclide))
                 && (fabs(pp->decayParticle()->energy - energy) < 0.01))
               || (thiselement && (pp->xrayElement() == thiselement)
                   && (fabs(energy - pp->xrayEnergy()) < 0.1)) )
            {
              currentlyused = true;
              if( dist < prevPeakDist )
              {
                prevpeak = pp;
                prevPeakDist = dist;
                prevIntensity = intensity;
              }
              break;
            }//if( this nulcide or xray is already taken )
          }//for( const PeakModel::PeakShrdPtr &pp : *previouspeaks )
          
          if( !currentlyused )
          {
            prevpeak.reset();
            color    = nuc.lineColor;
            nuclide  = thisnuclide;
            element  = thiselement;
            reaction = thisreaction;
            mindist  = dist;
            nearestEnergy = energy;
            thisIntensity = intensity;
          }//if( !currentlyused )
        }//if( we should possible associate this peak with this line )
      }//for( const BackgroundLine &line : nuc.backgroundLines )
    }//for( ReferenceLineInfo & )
    
    if( nuclide || reaction || element )
    {
      string src;
      char nuclide_label[128];
      
      if( !!prevpeak )
      {
        //There is an already existing peak, whos nuclide/xray/reaction gamma was
        //  "closer" (in terms of the distance metric used) than the
        //  nuclide/xray/reaction actually assigned to this new peak.  We need
        //  to check that we shouldnt swap them, based on relative intensities.
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
            throw std::logic_error( "InterSpec::addPeak(): bad logic "
                                   "checking previous peak gamma assignment" );
          
          snprintf( nuclide_label, sizeof(nuclide_label),
                   "%s %.6f keV", src.c_str(), nearestEnergy );
          
          
          other_change.reset( new std::pair<PeakModel::PeakShrdPtr,std::string>(prevpeak,nuclide_label) );
          
          nuclide = prevnuclide;
          reaction = prevreaction;
          element = prevelement;
          nearestEnergy = prevEnergy;
        }//if( we need to swap things )
      }//if( !!prevpeak )
      
      if( nuclide )
        src = nuclide->symbol;
      else if( reaction )
        src = reaction->name();
      else if( element )
        src = element->symbol + " xray";
      
      snprintf( nuclide_label, sizeof(nuclide_label),
               "%s %.6f keV", src.c_str(), nearestEnergy );
      
      PeakModel::setNuclideXrayReaction( peak, nuclide_label, -1.0 );
      
      if( colorPeaksBasedOnReferenceLines && peak.lineColor().isDefault() )
        peak.setLineColor( color );
    }//if( nuclide || reaction || element )
  }catch( std::exception &e )
  {
    passMessage( "Unexpected error searching for isotope for peak: "
                + string(e.what()), "", WarningWidget::WarningMsgHigh );
  }
  
  return other_change;
}//void assignNuclideFromReferenceLines( PeakDef &peak )

  
  
  
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
  
  auto currentdata = viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  
  if( (currentdata != originaldata) || !peaks )
  {
    passMessage( "Peaks not updated, spectrum has changed.",
                "", WarningWidget::WarningMsgHigh );
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
  
  passMessage( "Peaks updated.", "", WarningWidget::WarningMsgInfo );
  
  ReferencePhotopeakDisplay *refLineDisplay = viewer->referenceLinesWidget();
  
  //InterSpec::addPeak(...) is what normally associations nuclides with peaks
  vector<string> showingRefLines;
  if( refLineDisplay )
  {
    const ReferenceLineInfo &showingLines = refLineDisplay->currentlyShowingNuclide();
    if( !showingLines.parentLabel().empty() )
      showingRefLines.push_back( showingLines.parentLabel() );
    
    for( const ReferenceLineInfo &line : refLineDisplay->persistedNuclides() )
      if( !line.parentLabel().empty() )
        showingRefLines.push_back( line.parentLabel() );
  }//if( m_referenceNuclideLines )
  
  if( showingRefLines.empty() && viewer->isMobile() )
  {
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
  
  std::shared_ptr<WSvgImage> renderChartToSvg( std::shared_ptr<const SpecUtils::Measurement> inmeas,
                                              std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef> > > peaks,
                                              const std::vector<std::shared_ptr<const ReferenceLineInfo>> &displayed,
                                              double lowx, double upperx,
                                              const int width, const int height,
                                              std::shared_ptr<const ColorTheme> theme,
                                              const bool compact )
  {
    if( !inmeas )
      return nullptr;
    
    auto meas = std::make_shared<SpecUtils::Measurement>( *inmeas );
    
    std::shared_ptr<Wt::WSvgImage> img
    = std::make_shared<WSvgImage>( width, height );
    std::shared_ptr<SpecMeas> specmeas = std::make_shared<SpecMeas>();
    specmeas->add_measurement( meas, true );
    if( peaks )
      specmeas->setPeaks( *peaks, specmeas->sample_numbers() );
    
    PeakModel peakmodel;
    SpectrumDataModel dataModel;
    SpectrumChart chart;
    
    chart.setModel( &dataModel );
    chart.setPeakModel( &peakmodel );
    peakmodel.setDataModel( &dataModel );
    peakmodel.setPeakFromSpecMeas( specmeas, specmeas->sample_numbers() );
    
    const float liveTime = meas->live_time();
    const float realTime = meas->real_time();
    const double neutrons = meas->neutron_counts_sum();
    dataModel.setDataHistogram( meas, liveTime, realTime, neutrons );
    
    const vector<Chart::WDataSeries> series = dataModel.suggestDataSeries();
    chart.setSeries( series );
    
    //chart.enableLegend( true );
    
    chart.axis(Chart::YAxis).setScale( Chart::LogScale );
    
    for( size_t i = displayed.size(); i > 1; --i )
    {
      chart.setReferncePhotoPeakLines( *displayed[i-1] );
      chart.persistCurrentReferncePhotoPeakLines();
    }
    if( !displayed.empty() )
      chart.setReferncePhotoPeakLines( *displayed[0] );
    
    
    if( meas && (fabs(upperx-lowx) < 0.000001) )
    {
      const size_t nchannel = meas->num_gamma_channels();
      lowx = meas->gamma_channel_lower( 0 );
      upperx = meas->gamma_channel_upper( nchannel - 1 );
    }//if( lowx == upperx )
    
    chart.setXAxisRange( lowx, upperx );
    
    if( compact )
    {
      chart.axis(Chart::XAxis).setTitle("");
      chart.axis(Chart::YAxis).setTitle("");
      //chart.setLeftYAxisPadding(double width, <#double height#>)
      chart.setPlotAreaPadding( 35, Wt::Left );
      chart.setPlotAreaPadding( 16, Wt::Bottom );
      chart.setPlotAreaPadding( 0, Wt::Right );
      chart.setPlotAreaPadding( 0, Wt::Top );
      
      WFont labelFont( WFont::Default );
      labelFont.setSize(8);
      chart.axis(Chart::XAxis).setLabelFont( labelFont );
      chart.axis(Chart::YAxis).setLabelFont( labelFont );
    }//if( compact )
    
    const size_t displayednbin = meas->find_gamma_channel( upperx )
    - meas->find_gamma_channel( lowx );
    const int plotAreaWidth = static_cast<int>( img->width().toPixels() )
    - chart.plotAreaPadding(Left)
    - chart.plotAreaPadding(Right);
    const float bins_per_pixel = float(displayednbin) / float(plotAreaWidth);
    const int factor = max( static_cast<int>(ceil(bins_per_pixel)), 1 );
    
    dataModel.setRebinFactor( factor );
    
    chart.setAutoYAxisRange();
    
    if( theme )
    {
      dataModel.setForegroundSpectrumColor( theme->foregroundLine );
      dataModel.setBackgroundSpectrumColor( theme->backgroundLine );
      dataModel.setSecondarySpectrumColor( theme->secondaryLine );
      chart.setSeries( dataModel.suggestDataSeries() );
      
      chart.setDefaultPeakColor( theme->defaultPeakLine );
      
      chart.setAxisLineColor( theme->spectrumAxisLines );
      chart.setChartMarginColor( theme->spectrumChartMargins );
      chart.setChartBackgroundColor( theme->spectrumChartBackground );
      chart.setTextColor( theme->spectrumChartText );
    }//if( theme )
    
    
    WPainter p( img.get() );
    chart.paint( p );
    p.end();
    
    return img;
  }//renderChartToSvg(...)

void automated_search_for_peaks( InterSpec *viewer,
                                const bool keep_old_peaks )
{
  if( !viewer )
    return; //shouldnt happen
  
  auto peakModel = viewer->peakModel();
  if( !peakModel )
    return; //shouldnt happen
  
  
  auto dataPtr = viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );;
  
  if( !dataPtr )
  {
    passMessage( "There is no data to search for peaks on.", "", 1 );
    return;
  }//if( !dataPtr )
  
  //Could do a check to see if the user peaks have changed since the hint peaks
  //  were searched for, and if not, just use those...
  
  
  vector<ReferenceLineInfo> displayed;
  auto refLineDisp = viewer->referenceLinesWidget();
  if( refLineDisp )
  {
    const ReferenceLineInfo &currentNuclide = refLineDisp->currentlyShowingNuclide();
    displayed = refLineDisp->persistedNuclides();
  
    if( currentNuclide.nuclide || currentNuclide.reactionGammas.size()
       || currentNuclide.element || currentNuclide.backgroundLines.size() )
      displayed.insert( displayed.begin(), currentNuclide );
  }
  
  const bool setColor = viewer->colorPeaksBasedOnReferenceLines();
  
  //We should indicate to the user that seraching for peaks will take a while,
  //  but also we want to give them the chance to go past this and keep using
  //  the app.
  AuxWindow *msg = new AuxWindow( "Just a few moments",
                                 (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal) | AuxWindowProperties::PhoneModal) );
  msg->setClosable( true );
  msg->setModal( true );
  msg->show();
  msg->rejectWhenEscapePressed();
  msg->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, msg ) );
  new WText( "<div style=\"white-space: nowrap;\">Searching for peaks - this may take a bit.</div>"
             "<div style=\"white-space: nowrap;\">You can close this dialog and when the search is done,</div>"
             "you will be notified", Wt::XHTMLText, msg->contents() );
  msg->centerWindow();
  msg->disableCollapse();
  
  WPushButton *button = new WPushButton( "Close", msg->footer() );
  button->clicked().connect( msg, &AuxWindow::hide );
  
  //Make it so users cant keep clicking the search button
  viewer->automatedPeakSearchStarted();
  
  const auto originalPeaks = peakModel->peakVec();
  auto startingPeaks = peakModel->peaks();
  
  if( !keep_old_peaks )
    startingPeaks.reset();
  
  //using WApplication::bind to call msg->hide() will protect against the user
  //  closing the msg window (which will have deleted it)
  boost::function<void(void)> guiupdater
                          = wApp->bind( boost::bind( &AuxWindow::hide, msg ) );
  
  //The results of the peak search will be placed into the vector pointed to
  // by searchresults, which is why both 'callback' and below and
  // search_for_peaks_worker(...) get a shared pointer to this vector.
  auto searchresults = std::make_shared< vector<std::shared_ptr<const PeakDef> > >();
  
  //Again, use WApplication::bind to create the call that will set the peaks
  //  when the searching is done.
  boost::function<void(void)> callback
                          = wApp->bind( boost::bind( &set_peaks_from_search,
                                        viewer, displayed, searchresults, dataPtr, originalPeaks, guiupdater ) );
  
  Wt::WServer *server = Wt::WServer::instance();
  if( !server )
    return;
  
  std::weak_ptr<const SpecUtils::Measurement> weakdata = dataPtr;
  const string seshid = wApp->sessionId();
  
  server->ioService().post( std::bind( [=](){
    search_for_peaks_worker( weakdata, startingPeaks, displayed, setColor,
                            searchresults, callback, seshid, false );
    
  } ) );
}//void automated_search_for_peaks( InterSpec *interspec, const bool keep_old_peaks )

  
  
  
void assign_peak_nuclides_from_reference_lines( InterSpec *viewer )
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
    const ReferenceLineInfo &currentNuclide = refLineDisp->currentlyShowingNuclide();
    displayed = refLineDisp->persistedNuclides();
    
    if( currentNuclide.nuclide || currentNuclide.reactionGammas.size()
       || currentNuclide.element || currentNuclide.backgroundLines.size() )
      displayed.insert( displayed.begin(), currentNuclide );
  }
  
  auto resultpeaks = std::make_shared< std::vector<std::shared_ptr<const PeakDef> > >();
  for( const PeakDef &p : peaks )
    resultpeaks->push_back( make_shared<PeakDef>(p) );
  
  assign_srcs_from_ref_lines( foreground, resultpeaks, displayed, assignColor, showingEscapePeakFeature );
  
  for( const auto &foundpeak : *resultpeaks )
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
  }//for( const auto &foundpeak : *resultpeaks )
  
  
  //Hack because I'm short on time
  if( result_peaks.size() != orig_peaks.size() )
  {
    passMessage( "Unexpected issue assigning nuclides - accepting all assignments - sorry :(", "", WarningWidget::WarningMsgHigh );
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
      passMessage( "Unexpected issue assigning nuclides, mismatch in energy - accepting all assignments - sorry :(", "", WarningWidget::WarningMsgHigh );
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
  
  if( currentNuclide.nuclide || currentNuclide.reactionGammas.size()
     || currentNuclide.element || currentNuclide.backgroundLines.size() )
    displayed.insert( displayed.begin(), currentNuclide );
  
  if( displayed.empty() )
    return;
  
  std::shared_ptr<const deque< PeakModel::PeakShrdPtr > > previouspeaks = peakModel->peaks();
  if( !previouspeaks )  //probably never necassary, but JIC
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
                               std::shared_ptr<const deque< std::shared_ptr<const PeakDef> > > existingPeaks,
                               const vector<ReferenceLineInfo> displayed,
                               const bool setColor,
                               std::shared_ptr<std::vector<std::shared_ptr<const PeakDef> > > resultpeaks,
                               boost::function<void(void)> callback,
                               const std::string sessionID,
                               const bool singleThread )
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
  
  try
  {
    *resultpeaks = ExperimentalAutomatedPeakSearch::search_for_peaks( data, existingPeaks, singleThread );
    
    assign_srcs_from_ref_lines( data, resultpeaks, displayed, setColor, false );
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
  
  
  server->post( sessionID, callback );
}//void search_for_peaks_worker(...)

  
/**
*/
void assign_srcs_from_ref_lines( const std::shared_ptr<const SpecUtils::Measurement> &data,
                                  std::shared_ptr<std::vector<std::shared_ptr<const PeakDef> > > resultpeaks,
                                  const vector<ReferenceLineInfo> &displayed,
                                  const bool setColor,
                                  const bool showingEscapePeakFeature )
{
  if( !resultpeaks || !data )
    return;
  
  auto answerpeaks = make_shared< deque< shared_ptr<const PeakDef> > >();
  vector< shared_ptr<PeakDef> > unassignedpeaks;
  
  for( const auto &p : *resultpeaks )
  {
    if( p->parentNuclide() || p->xrayElement() || p->reaction() )
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
    auto addswap = assign_nuc_from_ref_lines( *peak, answerpeaks, data, displayed, setColor, showingEscapePeakFeature );
    
    if( addswap )  //The assignemnt caused a better assignment to be made for a previously existing peak.
    {
      auto pos = std::find( begin(*answerpeaks), end(*answerpeaks), addswap->first );
      if( pos != end(*answerpeaks) )
      {
        shared_ptr<const PeakDef> &oldpeak = *pos;
        auto moddedOldPeak = make_shared<PeakDef>( *oldpeak );
        
        PeakModel::SetGammaSource res = PeakModel::setNuclideXrayReaction( *moddedOldPeak, addswap->second, 1.0 );
        oldpeak = moddedOldPeak;
        if( res == PeakModel::NoSourceChange )
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
    
    
    //ToDo: use std::lower_bound to insert
    answerpeaks->push_back( peak );
    std::sort( begin(*answerpeaks), end(*answerpeaks), peak_less_than_by_energy );
  }//for( PeakDef peak : unassignedpeaks )
  
  resultpeaks->clear();
  for( auto &p : *answerpeaks )
    resultpeaks->push_back( p );
}//void assign_srcs_from_ref_lines(...)
  
  


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
    peak.continuum()->calc_linear_continuum_eqn( data, peak.lowerX(), peak.upperX(), 1 );
    peak.continuum()->setType( type );
    
    const double lowerx = std::max( peak.lowerX(), peak.mean() - 3*peak.sigma() );
    const double upperx = std::min( peak.upperX(), peak.mean() + 3*peak.sigma() );
    const double dataarea = data->gamma_integral( lowerx, upperx );
    const double contarea = peak.continuum()->offset_integral( lowerx, upperx );
    const double peakarea = (dataarea > contarea) ? (dataarea - contarea) : 10.0;
    peak.setPeakArea( peakarea );
    
    candidate_peaks.push_back( peak );
  }//for( auto &peak : input_peaks )
  
  
  const bool isRefit = true;
  const double x0 = data->gamma_energy_min();
  const double x1 = data->gamma_energy_max();
  const double ncausalitysigma = 0.0;
  const double stat_threshold  = 0.0;
  const double hypothesis_threshold = 0.0;
  
  vector<PeakDef> fitpeaks = fitPeaksInRange( x0, x1, ncausalitysigma,
                                             stat_threshold, hypothesis_threshold,
                                             candidate_peaks, data, orig_peaks, isRefit );
  
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
        reason = PeakSelectorWindowReason::PeaksFromPreviousSpectum;
        new PeakSelectorWindow( interspec, reason, input_peaks, data, fitpeaks, reflines );
        break;
    }//switch( fitsrc )

    wApp->triggerUpdate();
  } );
}//fit_template_peaks( ... )
  
}//namespace PeakSearchGuiUtils
