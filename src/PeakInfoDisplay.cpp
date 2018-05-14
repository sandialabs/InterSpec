/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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

// Disable streamsize <=> size_t warnings in boost
#pragma warning(disable:4244)

#include <memory>
#include <vector>
#include <iostream>
#include <stdexcept>

#include <Wt/WText>
#include <Wt/WImage>
#include <Wt/WLabel>
#include <Wt/WAnchor>
#include <Wt/WServer>
#include <Wt/WLineEdit>
#include <Wt/WTabWidget>
#include <Wt/WTreeView>
#include <Wt/WPushButton>
#include <Wt/WBorderLayout>
#include <Wt/WGridLayout>
#include <Wt/WBorderLayout>
#include <Wt/WApplication>
#include <Wt/WEnvironment>
#include <Wt/WItemDelegate>

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/PeakInfoDisplay.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/RowStretchTreeView.h"
#if ( USE_SPECTRUM_CHART_D3 )
#include "InterSpec/D3SpectrumDisplayDiv.h"
#else
#include "InterSpec/SpectrumDisplayDiv.h"
#endif
#include "InterSpec/IsotopeSelectionAids.h"
#include "InterSpec/ShieldingSourceDisplay.h"


using namespace Wt;
using namespace std;



PeakInfoDisplay::PeakInfoDisplay( InterSpec *viewer,
#if ( USE_SPECTRUM_CHART_D3 )
                                  D3SpectrumDisplayDiv *spectrumDisplayDiv,
#else
                                  SpectrumDisplayDiv *spectrumDisplayDiv,
#endif
                                  PeakModel *peakModel,
                                  Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_model( peakModel ),
    m_viewer( viewer ),
    m_spectrumDisplayDiv( spectrumDisplayDiv ),
    m_infoLayout( NULL ),
    m_infoView( NULL ),
    m_deletePeak( NULL ),
    m_searchForPeaks( NULL )
{
  assert( m_spectrumDisplayDiv );
  addStyleClass( "PeakInfoDisplay" );
  init();
//  setLayoutSizeAware(true);
}//PeakInfoDisplay constructor


PeakInfoDisplay::~PeakInfoDisplay()
{
}//~PeakInfoDisplay()


void PeakInfoDisplay::confirmRemoveAllPeaks()
{
  if( m_model->rowCount() <= 0 )
    return;
  
  AuxWindow *window = new AuxWindow( "Confirmation", true );
  window->disableCollapse();
  window->rejectWhenEscapePressed();
  WText * text = new WText("Erase All Peaks?");
  window->stretcher()->addWidget(text,0,0);
  
  WPushButton *button = new WPushButton("Yes");
  button->addStyleClass("BinIcon");
  window->footer()->addWidget(button);
  button->clicked().connect( window, &AuxWindow::hide );
#if ( USE_SPECTRUM_CHART_D3 )
  button->clicked().connect( boost::bind( &D3SpectrumDisplayDiv::removeAllPeaks, m_spectrumDisplayDiv ) );
#else
  button->clicked().connect( boost::bind( &PeakModel::removeAllPeaks, m_model ) );
#endif
  button->setFloatSide(Wt::Right);
  button = window->addCloseButtonToFooter("No");
  button->clicked().connect( window, &AuxWindow::hide );
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  window->show();
  window->setMinimumSize(200, WLength::Auto);
  window->centerWindow();
}//void confirmRemoveAllPeaks()


void PeakInfoDisplay::assignNuclidesFromRefLines()
{
  m_viewer->assignCurrentPeakNuclideFromReferenceLines();
}//void assignNuclidesFromRefLines()


void PeakInfoDisplay::enablePeakSearchButton( bool enable )
{
  m_searchForPeaks->setEnabled( enable );
}//void enablePeakSearchButton( bool enable = true )



void PeakInfoDisplay::createNewPeak()
{
  double width = 5.0;
  double amplitude = m_spectrumDisplayDiv->yAxisMaximum();

  if( m_model->npeaks() )
  {
    const PeakDef &firstPeak = m_model->peak( 0 );
    amplitude = firstPeak.peakArea();
    width = firstPeak.gausPeak() ? firstPeak.sigma() : 0.25*firstPeak.roiWidth();
  }//if( m_model->npeaks() )

  if( amplitude <= 0.0 )
    amplitude = 1000.0;

  double xmin = m_spectrumDisplayDiv->xAxisMinimum();
  double xmax = m_spectrumDisplayDiv->xAxisMaximum();

  std::shared_ptr<const Measurement> data = m_spectrumDisplayDiv->data();
  if( !!data )
  {
    const size_t nbin = data->num_gamma_channels();
    xmin = std::max( xmin, static_cast<double>(data->gamma_channel_lower(0)) );
    xmax = std::min( xmax, static_cast<double>(data->gamma_channel_upper(nbin-1)) );
  }

  const double energy = ceil( xmin + 0.01*(xmax - xmin) );
  m_viewer->addPeak( PeakDef( energy, width, amplitude ), false );

  passMessage( "Newly created peak assigned initial mean of "
              + std::to_string(energy)
              + " keV; please edit peak.",
              "", WarningWidget::WarningMsgInfo );

  //There is a slight issue here that this new peak will not drawn, since the
  //  chart doesnt know it needs to be repainted..., if we want to avoid this
  //  we can do something like:
//  m_spectrumDisplayDiv->refresh();
  
#if ( USE_SPECTRUM_CHART_D3 )
  m_spectrumDisplayDiv->updateData();
#endif
}//void createNewPeak()


void PeakInfoDisplay::deleteSelectedPeak()
{
  //Basically, for lack of better way, we have to figure out which row is
  //  currently being edited, and then remove the peak corresponding to that
  //  row.
  //  Since we will probably have less than ~20 rows, I'll not bother to find
  //  a better method...
  const int nrow = m_model->rowCount();
  const int ncol = m_model->columnCount();

  m_infoView->closeEditors();
    
  for( int row = 0; row < nrow; ++row )
  {
    for( int col = 0; col < ncol; ++col )
    {
      const WModelIndex index = m_model->index( row, col );
      if( m_infoView->isEditing( index ) || m_infoView->isSelected(index))
      {
        m_model->removeRows( row, 1 );
        disablePeakDelete();
#if ( USE_SPECTRUM_CHART_D3 )
        m_spectrumDisplayDiv->updateData();
#endif
        return;
      }//if( we found the field being edditied
    }//for( loop over columns )
  }//for( loop over rows )
}//void deleteSelectedPeak();


void PeakInfoDisplay::enablePeakDelete( WModelIndex index )
{
  if( m_deletePeak->isEnabled() || !index.isValid() )
    return;

  const WFlags<ItemFlag> flags = m_model->flags( index );
  if( (flags & Wt::ItemIsSelectable) == Wt::ItemIsSelectable )
  {
    m_deletePeak->enable();
  }
}//void enablePeakDelete()


void PeakInfoDisplay::disablePeakDelete()
{
    //Check if any row is selected.  If none are selected, then disable this button.
    
    const int nrow = m_model->rowCount();
    const int ncol = m_model->columnCount();
    
    for( int row = 0; row < nrow; ++row )
    {
        for( int col = 0; col < ncol; ++col )
        {
            const WModelIndex index = m_model->index( row, col );
            if(m_infoView->isSelected(index))
            {
                m_deletePeak->enable();
                return;
            }//if( we found a row selected
        }//for( loop over columns )
    }//for( loop over rows )
    
  m_deletePeak->disable(); //default
    
}//void disablePeakDelete()


void PeakInfoDisplay::init()
{
  assert( !m_infoView );
  assert( !m_infoLayout );

  const bool showToolTipInstantly = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_viewer );
  
  m_infoLayout = new WBorderLayout();
  setLayout(m_infoLayout);
  
  m_infoLayout->setContentsMargins( 1, 1, 1, 1 );
  
  m_infoView = new RowStretchTreeView();
  m_infoView->setRootIsDecorated	(	false); //makes the tree look like a table! :)
  
  if( m_model )
    m_infoView->setModel( m_model );
  else
    cerr << SRC_LOCATION
         << "\n\tPeakInfoDisplay::createInfoTab(): Invalid peak model--"
            "GUI will not function correctly." << endl<< endl;

//  m_infoView->setRowHeight(28);
//  m_infoView->setHeaderHeight(28);
  m_infoView->setSortingEnabled( true );
  m_infoView->setAlternatingRowColors( true );
  m_infoView->setSelectable( true );
  m_infoView->setSelectionMode( SingleSelection );
  setLayout( m_infoLayout );
  
  m_infoView->setColumnHidden( PeakModel::kCandidateIsotopes, true );
  m_infoView->setColumnHidden( PeakModel::kUseForCalibration, true );
  m_infoView->setColumnHidden( PeakModel::kUseForShieldingSourceFit, true );

  m_infoView->setEditTriggers( WAbstractItemView::SingleClicked | WAbstractItemView::DoubleClicked );

  WItemDelegate *dblDelagate = new WItemDelegate( m_infoView );
  dblDelagate->setTextFormat( "%.2f" );
  m_infoView->setItemDelegateForColumn( PeakModel::kMean, dblDelagate );
  m_infoView->setItemDelegateForColumn( PeakModel::kFwhm, dblDelagate );
  m_infoView->setItemDelegateForColumn( PeakModel::kLowerX, dblDelagate );
  m_infoView->setItemDelegateForColumn( PeakModel::kUpperX, dblDelagate );

  dblDelagate = new WItemDelegate( m_infoView );
  dblDelagate->setTextFormat( "%.0f" );
  m_infoView->setItemDelegateForColumn( PeakModel::kAmplitude, dblDelagate );
  m_infoView->setItemDelegateForColumn( PeakModel::kContinuumArea, dblDelagate );

  //tweak column widths
  m_infoView->setColumnWidth( PeakModel::kMean,       WLength(8, WLength::FontEx) );
  m_infoView->setColumnWidth( PeakModel::kFwhm,       WLength(9, WLength::FontEx) );
  m_infoView->setColumnWidth( PeakModel::kAmplitude,  WLength(12, WLength::FontEx) );

  m_infoView->setColumnWidth( PeakModel::kIsotope,    WLength(/*8*/9, WLength::FontEx) );
  //Note 20131211, wcjohns: closeOnBlur was previoulsy set to false, however
  //  there was a rare crash that happened with the model row was deleted while
  //  the editor was still open.  Hopefully closeOnBlur==true will fix this
  //  (appears to).
  const bool closeOnBlur = true;
  PhotopeakDelegate *nuclideDelegate = new PhotopeakDelegate( PhotopeakDelegate::NuclideDelegate, closeOnBlur, m_infoView );
  m_infoView->setItemDelegateForColumn( PeakModel::kIsotope, nuclideDelegate );

  m_infoView->setColumnWidth( PeakModel::kPhotoPeakEnergy, WLength(12 /*18*/, WLength::FontEx) );
  PhotopeakDelegate *photopeakDelegate = new PhotopeakDelegate( PhotopeakDelegate::GammaEnergyDelegate, closeOnBlur, m_infoView );
  m_infoView->setItemDelegateForColumn( PeakModel::kPhotoPeakEnergy, photopeakDelegate );

  m_infoView->setColumnWidth( PeakModel::kUserLabel,  WLength(9, WLength::FontEx) );
  
  m_infoView->setColumnHidden( PeakModel::kHasSkew, true );
//  m_infoView->setColumnWidth( PeakModel::kHasSkew,    WLength(9, WLength::FontEx) );
  
  m_infoView->setColumnHidden( PeakModel::kSkewAmount, true );
//  m_infoView->setColumnWidth( PeakModel::kSkewAmount, WLength(/*10*/13, WLength::FontEx) );
  
  m_infoView->setColumnHidden( PeakModel::kType, true );
//  m_infoView->setColumnWidth( PeakModel::kType,       WLength(12, WLength::FontEx) );
  
  m_infoView->setColumnWidth( PeakModel::kLowerX,     WLength(/*11*/16, WLength::FontEx) );
  m_infoView->setColumnWidth( PeakModel::kUpperX,     WLength(/*11*/16, WLength::FontEx) );
  m_infoView->setColumnWidth( PeakModel::kContinuumArea, WLength(/*13*/16, WLength::FontEx) );
  m_infoView->setColumnWidth( PeakModel::kContinuumType,  WLength(/*11*/16, WLength::FontEx) );
//  m_infoView->setColumnWidth( PeakModel::kNumColumns, WLength(6, WLength::FontEx) );
//  m_infoView->setColumnAlignment( PeakModel::kMean, AlignRight );

//  Selection mode seems to have no effect....
//  m_infoView->setSelectionBehavior( Wt::SelectRows );
//  m_infoView->setSelectionMode( Wt::SingleSelection );
  m_infoView->setEditOptions( WAbstractItemView::SingleEditor );

  m_infoView->clicked().connect( boost::bind( &PeakInfoDisplay::enablePeakDelete, this, _1 ) );
  m_infoView->doubleClicked().connect( boost::bind( &PeakInfoDisplay::enablePeakDelete, this, _1 ) );

  m_infoLayout->addWidget( m_infoView, Wt::WBorderLayout::Center);
  

  //Now add buttons to search/clear peaks
  buttonDiv = new WContainerWidget();
  m_infoLayout->addWidget( buttonDiv, Wt::WBorderLayout::South);

  m_searchForPeaks = new WPushButton( "Search for Peaks", buttonDiv );
  m_searchForPeaks->addStyleClass("FindIcon");
  m_searchForPeaks->setMargin(WLength(7),Wt::Left);
  m_searchForPeaks->setMargin(WLength(3),Wt::Bottom);
  HelpSystem::attachToolTipOn( m_searchForPeaks, "Search for peaks using the automated peak finding "
                              "algorithm.", showToolTipInstantly, HelpSystem::Top  );
  m_searchForPeaks->clicked().connect( boost::bind( &InterSpec::searchForPeaks, m_viewer, true ) );

  
  WPushButton *button = new WPushButton( "Clear all Peaks", buttonDiv );
  HelpSystem::attachToolTipOn( button,"Removes <b>all</b> existing peaks.", showToolTipInstantly, HelpSystem::Top  );
  button->addStyleClass("BinIcon");
        button->setMargin(WLength(2),Wt::Left);
  button->clicked().connect( this, &PeakInfoDisplay::confirmRemoveAllPeaks );

  //"Nuc. from Ref."
  button = new WPushButton( "Nuc. from Ref.", buttonDiv );
  button->addStyleClass("WandIcon");
  button->setMargin(WLength(2),Wt::Left|Wt::Right);
  HelpSystem::attachToolTipOn( button,
                              "Assign peak nuclides from reference lines showing. Only applies to "
                              "peaks which do not already have a nuclide associated "
                              "with them." ,
                              showToolTipInstantly , HelpSystem::Top );
  button->clicked().connect( boost::bind( &PeakInfoDisplay::assignNuclidesFromRefLines, this ) );
  
/*
  button = new WPushButton( "ID Nuclides", buttonDiv );
  button->addStyleClass("WandIcon");
  button->setMargin(WLength(2),Wt::Left|Wt::Right);
  HelpSystem::attachToolTipOn( button,
                              "Guess nuclides responsible for peaks. Only applies to "
                              "peaks which do not already have a nuclide associated "
                              "with them.  Works best once all peaks have been fit for." ,
                              showToolTipInstantly , HelpSystem::Top );
  button->clicked().connect( boost::bind( &InterSpec::guessIsotopesForPeaks, m_viewer, (WApplication *)0 ) );
*/
  
  //If you want to post the bellow, so the ID isnt carried out in the main event loop, use the following
//  boost::function<void ()> guessIsotopeWorker = boost::bind( &InterSpec::guessIsotopesForPeaks, m_viewer, wApp );
//  button->clicked().connect( boost::bind( &WServer::post, WServer::instance(),
//                                           wApp->sessionId(), guessIsotopeWorker,
//                                           boost::function<void ()>() ) );
  
  WLabel *label = new WLabel("Peak: ", buttonDiv);
  label->addStyleClass("buttonSeparator");
  label->setMargin(WLength(10),Wt::Left);
  WPushButton *addPeak = new WPushButton("Add", buttonDiv );

  HelpSystem::attachToolTipOn( addPeak,"Add a new peak for manual editing; peak will "
                              "have initial energy near left side of plot. Peak parameters can be editted by double-clicking on the quantity in the <b>Peak Manager</b> tab.", showToolTipInstantly , HelpSystem::Top );
  
  addPeak->clicked().connect( this, &PeakInfoDisplay::createNewPeak );
  addPeak->addStyleClass( "AddIcon" );

  m_deletePeak = new WPushButton( "Delete", buttonDiv );
  m_deletePeak->setMargin(WLength(2),Wt::Left|Wt::Right);
  HelpSystem::attachToolTipOn( m_deletePeak,"Deletes peak currently being edited.", showToolTipInstantly, HelpSystem::Top  );
  m_deletePeak->clicked().connect( this, &PeakInfoDisplay::deleteSelectedPeak );
  m_deletePeak->setHiddenKeepsGeometry( true );
  m_deletePeak->disable();
   m_deletePeak->addStyleClass( "DeleteIcon" );

  //Whenver a delegate gets closed, lets disable the peak delete button
  set<WAbstractItemDelegate *> uniqueDelegates;
  for( int col = 0; col < m_model->columnCount(); ++col )
  {
    WAbstractItemDelegate *delegate = m_infoView->itemDelegateForColumn( col );
    if( delegate && !uniqueDelegates.count(delegate) )
    {
      delegate->closeEditor().connect( this, &PeakInfoDisplay::disablePeakDelete );
      uniqueDelegates.insert( delegate );
    }
  }//for( int col = 0; col < m_model->columnCount(); ++col )

#if( !ANDROID && !IOS )
  WText *txt = new WText( "Right click on peaks for better editor", buttonDiv );
  txt->addStyleClass( "PeakEditHint" );
  txt->hide();
//  mouseWentOver().connect( txt, &WText::show );
//  mouseWentOut().connect( txt, &WText::hide );
//  doJavaScript( "try{$('#" + txt->id() + "').show();}catch(e){}" );
  mouseWentOver().connect( "function(object, event){try{$('#" + txt->id() + "').show();}catch(e){}}" );
  mouseWentOut().connect( "function(object, event){try{$('#" + txt->id() + "').hide();}catch(e){}}" );
  
  WResource *csv = m_model->peakCsvResource();
  WAnchor *csvAnchor = new WAnchor ( WLink(csv), "CSV", buttonDiv );
//  csvAnchor->setPadding( Wt::Right );
  csvAnchor->addStyleClass("DiskIcon");
  csvAnchor->setTarget( Wt::TargetNewWindow );
  csvAnchor->setFloatSide( Wt::Right );

  HelpSystem::attachToolTipOn( csvAnchor,"Export information about the identified peaks to a "
                              "comma seperated format.", showToolTipInstantly );
#endif
  
}//init()


