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
#include <Wt/WGridLayout>
#include <Wt/WBorderLayout>
#include <Wt/WApplication>
#include <Wt/WEnvironment>
#include <Wt/WItemDelegate>

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/ColorSelect.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/PeakInfoDisplay.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#if( USE_SPECTRUM_CHART_D3 )
#include "InterSpec/D3SpectrumDisplayDiv.h"
#else
#include "InterSpec/SpectrumDisplayDiv.h"
#endif
#include "InterSpec/IsotopeSelectionAids.h"
#include "InterSpec/ShieldingSourceDisplay.h"


using namespace Wt;
using namespace std;

namespace
{
#if( ALLOW_PEAK_COLOR_DELEGATE )
  class ColorDelegate : public Wt::WAbstractItemDelegate
  {
  public:
    class EditWidget : public Wt::WContainerWidget
    {
    public:
      EditWidget( const Wt::WModelIndex& index,
                 const Wt::WFlags<Wt::ViewItemRenderFlag> flags,
                 ColorDelegate *parent )
      {
        m_color = nullptr;
        m_parent = parent;
        
        auto model = dynamic_cast<const PeakModel *>( index.model() );
        m_peakModel = const_cast<PeakModel *>(model);
        
        if( model )
        {
          try
          {
            m_origPeak = model->peak( index );
            if( m_origPeak )
              m_origColor = m_origPeak->lineColor();
          }catch(...){}
        }
        updateRender( index, flags );
      }//EditWidget constructor
      
      virtual ~EditWidget(){}
      
      WString cssColor()
      {
        if( m_color && !m_color->valueText().empty() )
          return m_color->valueText();
        else
          return "none";
      }//cssColor()
      
      void setCssCollor( const WString &s )
      {
        cout << "setCssCollor: " << s.toUTF8() << endl;
        const string ss = s.toUTF8();
        if( m_color )
          m_color->setValueText( (ss.empty() || ss=="none") ? "" : s );
        else
        {
          if( ss.empty() || ss=="none" )
            setAttributeValue( "style", "margin: 5px; background-color: null; border: 1px #e1e1e1" );
          else
            setAttributeValue( "style", "margin: 5px; border: 1px #e1e1e1; background-color: " + ss );
        }
      }//void setCssCollor( WString &s )
      
      void closeEditor( bool save )
      {
        auto index = m_peakModel->indexOfPeak(m_origPeak);
        
        WString color = m_origColor;
        if( index.isValid() && save && m_color && m_color->valueText().toUTF8()!="#ebebeb" )
          color = m_color->valueText();
        
        m_peakModel->setData( index, color );
        if( m_color )
          delete m_color;
        m_color = nullptr;
        updateRender(index,0);
      }
      
      void colorSelected( const std::string &color )
      {
        closeEditor( color!=m_origColor.toUTF8() );
      }
      
      void updateRender( const Wt::WModelIndex &index,
                        const Wt::WFlags<Wt::ViewItemRenderFlag> flags )
      {
        cout << "updateRender: flags=" << flags.value() << endl;
        
        auto model = dynamic_cast<const PeakModel *>( index.model() );
        
        if( model && (flags & RenderEditing) ) //RenderSelected, RenderEditing, RenderFocused, RenderInvalid
        {
          if( !m_color )
          {
            m_color = new ColorSelect(ColorSelect::PrefferNative,this);
            setAttributeValue( "style", "margin: 0px; background-color: null;" );
            m_color->setAttributeValue("style", "height: 15px;" );
            m_color->cssColorChanged().connect( boost::bind( &EditWidget::colorSelected, this, _1 ) );
          }
          try
          {
            const auto &p = model->peak(index);
            m_color->setValueText( p->lineColor().empty() ? "#ebebeb" : p->lineColor().c_str() );
          }catch(...)
          {
          }
        }else
        {
          try
          {
            if( !model )
              throw runtime_error("");
            
            const auto &p = model->peak( index );
            if( p->lineColor().empty() )
              setAttributeValue( "style", "margin: 5px; background-color: null; border: 1px #e1e1e1" );
            else
              setAttributeValue( "style", "margin: 5px; border: 1px #e1e1e1; background-color: " + p->lineColor() );
          }catch(...)
          {
            setAttributeValue( "style", "margin: 5px; border: 1px #e1e1e1; background-color: grey" );
          }
        }
      }//void updateRender(...)
      
    protected:
      Wt::WString m_origColor;
      PeakModel::PeakShrdPtr m_origPeak;
      PeakModel *m_peakModel;
      ColorSelect *m_color;
      ColorDelegate *m_parent;
    };//class EditWidget
    
  public:
    ColorDelegate( Wt::WObject *parent = NULL )
      : WAbstractItemDelegate(parent)
    {
      closeEditor().connect( boost::bind( &ColorDelegate::doCloseEditor, this, _1, _2 ) );
    }
    
    virtual ~ColorDelegate(){}
    
    virtual Wt::WWidget *update( Wt::WWidget *widget,
                                const Wt::WModelIndex &index,
                                Wt::WFlags< Wt::ViewItemRenderFlag > flags )
    {
      EditWidget *edit = dynamic_cast<EditWidget *>(widget);
      if( !edit )
        edit = new EditWidget( index, flags, this );
      else
        edit->updateRender( index, flags );
      
      return edit;
    }//WWidget *update(...)
    
    
    void doCloseEditor( WWidget *editor, bool save )
    {
      EditWidget *edit = dynamic_cast<EditWidget *>(editor);
      if( edit )
        edit->closeEditor(save);
    }//void doCloseEditor( WWidget *editor, bool save )
    
  protected:
    
 
    boost::any editState( Wt::WWidget *editor ) const
    {
      EditWidget *w = dynamic_cast<EditWidget *>(editor);
      if( !w )
      {
        cerr << SRC_LOCATION << "\n\tLogic error - fix me!" << endl;
        return boost::any();
      }//if( !w )
      
      WString csscolor = "none";
      if( w )
        csscolor = w->cssColor();
      
      return boost::any( csscolor );
    }//boost::any editState( WWidget *editor ) const

    
    void setEditState( Wt::WWidget *editor, const boost::any &value ) const
    {
      EditWidget *w = dynamic_cast<EditWidget *>(editor);
      if( !w )
      {
        cerr << SRC_LOCATION << "\n\tLogic error - fix me!" << endl;
        return;
      }
      
      try
      {
        w->setCssCollor( boost::any_cast<WString>(value) );
      }catch(...)
      {
        cerr << SRC_LOCATION << "\n\tPossible Logic error - fix me!" << endl;
      }//try / catch
    }//void setEditState( WWidget *editor, const boost::any& value ) const

    
    void setModelData( const boost::any &editState,
                      Wt::WAbstractItemModel *model,
                      const Wt::WModelIndex &index ) const
    {
      if( model )
        model->setData( index, editState, EditRole );
    }//void setModelData(...)
  };//class ColorDelegate
#endif //ALLOW_PEAK_COLOR_DELEGATE
}//namespace


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
  
  AuxWindow *window = new AuxWindow( "Confirmation",
              (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal) | AuxWindowProperties::TabletModal | AuxWindowProperties::DisableCollapse) );
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
  PeakSearchGuiUtils::assign_peak_nuclides_from_reference_lines( m_viewer );
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
  if( !m_model )
    throw runtime_error( "PeakInfoDisplay must be passed a valid PeakModel" );
  
  assert( !m_infoView );
  assert( !m_infoLayout );

  const bool showToolTipInstantly = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_viewer );
  
  m_infoLayout = new WGridLayout();
  setLayout( m_infoLayout );
  
  m_infoLayout->setContentsMargins( 0, 0, 0, 0 );
  m_infoLayout->setVerticalSpacing( 0 );
  m_infoLayout->setHorizontalSpacing( 0 );
  
  m_infoView = new RowStretchTreeView();
  m_infoView->addStyleClass( "PeakInfoDisplayTable" );
  m_infoView->setRootIsDecorated(	false); //makes the tree look like a table! :)
  
  m_infoView->setModel( m_model );

//  m_infoView->setRowHeight(28);
//  m_infoView->setHeaderHeight(28);
  m_infoView->setSortingEnabled( true );
  m_infoView->setAlternatingRowColors( true );
  m_infoView->setSelectable( true );
  m_infoView->setSelectionMode( SingleSelection );
  
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
  
#if( ALLOW_PEAK_COLOR_DELEGATE )
  m_infoView->setColumnWidth( PeakModel::kPeakLineColor, WLength(6,WLength::FontEx) );
  ColorDelegate *colorDelegate = new ColorDelegate( m_infoView );
  m_infoView->setItemDelegateForColumn( PeakModel::kPeakLineColor, colorDelegate );
#else
  m_infoView->setColumnHidden( PeakModel::kPeakLineColor, true );
#endif
  
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

  m_infoLayout->addWidget( m_infoView, 0, 0 );
  m_infoLayout->setRowStretch( 0, 1 );
  

  //Now add buttons to search/clear peaks
  WContainerWidget *buttonDiv = new WContainerWidget();
  buttonDiv->addStyleClass( "PeakInfoDisplayButtonDiv" );
  m_infoLayout->addWidget( buttonDiv, 1, 0 );

  m_searchForPeaks = new WPushButton( "Search for Peaks", buttonDiv );
  m_searchForPeaks->addStyleClass("FindIcon");
  //m_searchForPeaks->setMargin(WLength(7),Wt::Left);
  //m_searchForPeaks->setMargin(WLength(3),Wt::Bottom);
  HelpSystem::attachToolTipOn( m_searchForPeaks, "Search for peaks using the automated peak finding "
                              "algorithm.", showToolTipInstantly, HelpSystem::Top  );
  m_searchForPeaks->clicked().connect( boost::bind( &PeakSearchGuiUtils::automated_search_for_peaks, m_viewer, true ) );

  
  WPushButton *clearPeaksButton = new WPushButton( "Clear all Peaks", buttonDiv );
  HelpSystem::attachToolTipOn( clearPeaksButton, "Removes <b>all</b> existing peaks.", showToolTipInstantly, HelpSystem::Top  );
  clearPeaksButton->addStyleClass("BinIcon");
  //clearPeaksButton->setMargin(WLength(2),Wt::Left);
  clearPeaksButton->clicked().connect( this, &PeakInfoDisplay::confirmRemoveAllPeaks );
  clearPeaksButton->disable();

  //"Nuc. from Ref."
  WPushButton *nucFromRefButton = new WPushButton( "Nuc. from Ref.", buttonDiv );
  nucFromRefButton->setIcon( "InterSpec_resources/images/assign_white.png" );
  
  //button->setMargin(WLength(2),Wt::Left|Wt::Right);
  HelpSystem::attachToolTipOn( nucFromRefButton,
                              "Assign peak nuclides from reference lines showing. Only applies to "
                              "peaks which do not already have a nuclide associated "
                              "with them." ,
                              showToolTipInstantly , HelpSystem::Top );
  nucFromRefButton->clicked().connect( boost::bind( &PeakInfoDisplay::assignNuclidesFromRefLines, this ) );
  nucFromRefButton->disable();
  
  auto enableDisableNucRef = [this,nucFromRefButton,clearPeaksButton](){
    const bool enable = (m_model->rowCount() > 0);
    clearPeaksButton->setEnabled( enable );
    nucFromRefButton->setEnabled( enable );
    //Should check if any reference lines are showing for nucFromRefButton as well...
  };
  
  m_model->dataChanged().connect( std::bind(enableDisableNucRef) );
  m_model->rowsRemoved().connect( std::bind(enableDisableNucRef) );
  m_model->rowsInserted().connect( std::bind(enableDisableNucRef) );
  m_model->layoutChanged().connect( std::bind(enableDisableNucRef) );
  
/*
  WPushButton *button = new WPushButton( "ID Nuclides", buttonDiv );
  button->setIcon( "InterSpec_resources/images/assign_white.png" );
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
  //addPeak->addStyleClass( "PlusIconWhiteBtn" );
  addPeak->setIcon( "InterSpec_resources/images/plus_min_white.png" );
  
  m_deletePeak = new WPushButton( "Delete", buttonDiv );
  //m_deletePeak->setMargin(WLength(2),Wt::Left|Wt::Right);
  HelpSystem::attachToolTipOn( m_deletePeak,"Deletes peak currently being edited.", showToolTipInstantly, HelpSystem::Top  );
  m_deletePeak->clicked().connect( this, &PeakInfoDisplay::deleteSelectedPeak );
  m_deletePeak->setHiddenKeepsGeometry( true );
  m_deletePeak->disable();
  //m_deletePeak->addStyleClass( "MinusIconWhiteBtn" );
  m_deletePeak->setIcon( "InterSpec_resources/images/minus_min_white.png" );

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
  WPushButton *csvButton = new WPushButton( buttonDiv );
  csvButton->setIcon( "InterSpec_resources/images/download_small.png" );
  csvButton->setLink( WLink(csv) );
  csvButton->setLinkTarget( Wt::TargetNewWindow );
  csvButton->setText( "CSV" );
  csvButton->setStyleClass( "CsvLinkBtn" );
  csvButton->disable();

  auto enableDisableCsv = [this,csvButton](){
    csvButton->setEnabled( m_model->rowCount() > 0 );
  };
  
  m_model->dataChanged().connect( std::bind(enableDisableCsv) );
  m_model->rowsRemoved().connect( std::bind(enableDisableCsv) );
  m_model->rowsInserted().connect( std::bind(enableDisableCsv) );
  m_model->layoutChanged().connect( std::bind(enableDisableCsv) );
  
  HelpSystem::attachToolTipOn( csvButton,"Export information about the identified peaks to a "
                              "comma seperated format.", showToolTipInstantly );
#endif
  
}//init()


