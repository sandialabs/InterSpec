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
#include <Wt/WTable>
#include <Wt/WImage>
#include <Wt/WLabel>
#include <Wt/WAnchor>
#include <Wt/WServer>
#include <Wt/WSvgImage>
#include <Wt/WComboBox>
#include <Wt/WLineEdit>
#include <Wt/WPopupMenu>
#include <Wt/WTableCell>
#include <Wt/WTabWidget>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WApplication>
#include <Wt/WEnvironment>
#include <Wt/WItemDelegate>

#include "SpecUtils/Filesystem.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/ColorSelect.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/PeakInfoDisplay.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/AddNewPeakDialog.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeSelectionAids.h"
#include "InterSpec/ShieldingSourceDisplay.h"



using namespace Wt;
using namespace std;


#if( ANDROID )
// Defined in target/android/android.cpp
extern void android_download_workaround( Wt::WResource *resource, std::string description );
#endif


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
            m_color->cssColorChanged().connect( boost::bind( &EditWidget::colorSelected, this,
                                                            boost::placeholders::_1 ) );
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
      closeEditor().connect( boost::bind( &ColorDelegate::doCloseEditor, this,
                                         boost::placeholders::_1, boost::placeholders::_2 ) );
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
        cerr << "editState():\n\tLogic error - fix me!" << endl;
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
        cerr << "setEditState(...)\n\tLogic error - fix me!" << endl;
        return;
      }
      
      try
      {
        w->setCssCollor( boost::any_cast<WString>(value) );
      }catch(...)
      {
        cerr << "setEditState(...)\n\tPossible Logic error - fix me!" << endl;
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

/** Simple specialization of WItemDelegate that closes the editor on blur, if the user
 has pressed any keys that cause the value to change.
 
 TODO: Should probably totally re-implement WItemDelegate for doubles, and to get the exact
       behavior we would like, rather than hacking on top of WItemDelegate, which isnt great
 */
class ItemDelegate : public Wt::WItemDelegate
{
public:
  ItemDelegate(WObject *parent = 0)
  : WItemDelegate( parent )
  {
  }
   
protected:
  
  // Make sure we it is still safe to access lineEdit
  //  (we get some decodeSignal() in the C++ if you press Esc or Enter while editing - even without
  //   this checkWidgetInDom we dont get any crashes, but this is just to make sure)
  bool checkWidgetInDom( const string &lineEditID ) const
  {
    auto app = wApp;
    if( !app || !app->domRoot()->findById(lineEditID) )
    {
      cerr << "ItemDelegate: Somehow lineEdit disappeared from DOM" << endl;
      return false;
    }
    return true;
  }//checkWidgetInDom(... )
  
  
  void closeOnBlur( WWidget *editor, WLineEdit *lineEdit, const string &lineEditID, shared_ptr<bool> changed ) const
  {
    if( !checkWidgetInDom(lineEditID) )
    {
      closeEditor().emit(editor, false);
      return;
    }
    
    bool save = false;
    if( changed && (*changed) && lineEdit && !lineEdit->text().empty() )
      save = true;

    closeEditor().emit(editor, save);
  }//closeOnBlur
  
  void onKeyDown( WLineEdit *e, const string &eID, shared_ptr<bool> changed,
                  shared_ptr<Wt::Signals::connection> sig, WKeyEvent k ) const
  {
    if( !checkWidgetInDom(eID) )
      return;
    
    assert( changed );
    
    switch( k.key() )
    {
      case Wt::Key_Enter: case Wt::Key_Tab: case Wt::Key_Shift: case Wt::Key_Control:
      case Wt::Key_Alt: case Wt::Key_PageUp: case Wt::Key_PageDown: case Wt::Key_End:
      case Wt::Key_Home: case Wt::Key_Left: case Wt::Key_Up: case Wt::Key_Right:
      case Wt::Key_Down: case Wt::Key_Escape: case Wt::Key_Insert:
      case Wt::Key_F1: case Wt::Key_F2: case Wt::Key_F3: case Wt::Key_F4: case Wt::Key_F5:
      case Wt::Key_F6: case Wt::Key_F7: case Wt::Key_F8: case Wt::Key_F9: case Wt::Key_F10:
      case Wt::Key_F11: case Wt::Key_F12:
        break;
        
        // Only set changed if the key actually made a difference.
      default:
        if( changed)
          (*changed) = true;
        if( sig )
          e->keyWentDown().disconnect( *sig );
        break;
    }//switch( k.key() )
  }//onKeyDown
  
  
  virtual WWidget *createEditor( const WModelIndex& index,
                                WFlags<ViewItemRenderFlag> flags) const
  {
    WWidget *w = WItemDelegate::createEditor(index, flags);
    auto div = dynamic_cast<WContainerWidget *>( w );
    assert( div );
    if( !div )
      return w;
    
    for( WWidget *d : div->children() )
    {
      auto e = dynamic_cast<WLineEdit *>( d );
      if( e )
      {
        auto changed = make_shared<bool>( false );
        auto keydownsig = make_shared<Wt::Signals::connection>();
        e->blurred().connect( boost::bind( &ItemDelegate::closeOnBlur, this, w, e, e->id(), changed ) );
        
        *keydownsig = e->keyWentDown().connect( boost::bind( &ItemDelegate::onKeyDown, this,
                                                e, e->id(), changed, keydownsig, boost::placeholders::_1 ) );
      }//if( e )
    }//
    
    return w;
  }//createEditor(...)
};//class ItemDelegate

  
// See also CopyUrlToClipboard in QrCode.cpp, and CopyFluxDataTextToClipboard in FluxTool.cpp
//  Note that modern browsers ONLY allow copying 'text/plain' and 'text/html' mime types
//  to the pasteboard (for security reasons); e.g., using 'text/csv' silently fails.
//  We copy the data to the pasteboard, in both plain text format, and HTML table format
//  (which pastes into Excel correctly).
WT_DECLARE_WT_MEMBER
(CopyPeakCsvDataToClipboard, Wt::JavaScriptFunction, "CopyPeakCsvDataToClipboard",
function( sender, event, dataOwner, dataName )
{
  const csvData = $(dataOwner).data(dataName);
  const htmlData = $(dataOwner).data(dataName + "Html");
  if( !csvData || !htmlData )
  {
    Wt.emit( $('.specviewer').attr('id'), {name:'miscSignal'},
            'peakCsvCopy-error-No CSV data is available - if this is an error, please report to InterSpec@sandia.gov' );
    return;
  }
  
  // Use ClipboardItem if supported
  if( typeof ClipboardItem !== 'undefined' ) {
    const plainBlob = new Blob([csvData], { type: 'text/plain' });
    const htmlBlob = new Blob([htmlData], { type: 'text/html' });
    const clipboardItem = new ClipboardItem({ 'text/plain': plainBlob, 'text/html': htmlBlob });

    navigator.clipboard.write([clipboardItem]).then(function() {
      console.log( 'Copied peak CSV using ClipboardItem method.' );
      Wt.emit( $('.specviewer').attr('id'), {name:'miscSignal'},
              'peakCsvCopy-success-Copied peaks CSV data to the pasteboard.' );
    }).catch(function(err) {
      console.warn( 'Failed to copy peak CSV using ClipboardItem method: ', err );
      
      Wt.emit( $('.specviewer').attr('id'), {name:'miscSignal'},
                  'peakCsvCopy-error-Unable to copy peak CSV info to the pasteboard - sorry.' );
    });
    
    return;
  }//if( typeof ClipboardItem !== 'undefined' )
  
      
  try
  {
    //This bit of code seems to work on Chrome, but not safari
    let didcopy = 0;
    function listener(e) {
      e.clipboardData.setData("text/html", htmlData);
      e.clipboardData.setData("text/plain", csvData);
      didcopy = 1;
      e.preventDefault();
    }
    document.addEventListener("copy", listener);
    document.execCommand("copy");
    document.removeEventListener("copy", listener);
        
    if( didcopy ){
      console.log( 'Copied peak CSV using exec method.' );
      Wt.emit( $('.specviewer').attr('id'), {name:'miscSignal'}, 'peakCsvCopy-success-Copied CSV data to pasteboard.' );
      return;
    }else{
      console.warn( 'Failed copy as csv to e.clipboardData.setData' );
    }
  }catch(error){
    console.warn( 'Failed to copy rich text to copyboard (e.clipboardData.setData)' );
  }
  
  if( window.clipboardData && window.clipboardData.setData ) {
    const didCopyTxt = window.clipboardData.setData("text/plain", csvData);  // IE
    const didCopyHtml = window.clipboardData.setData("text/html", htmlData);  // IE
    if( didCopyTxt && didCopyHtml ){
      console.log( 'Copied peak CSV using window.clipboardData.setData.' );
      Wt.emit( $('.specviewer').attr('id'), {name:'miscSignal'},
                'peakCsvCopy-success-Copied CSV data to pasteboard as csv data.' );
      return;
    }else{
      console.warn( 'Failed to copy peak CSV using window.clipboardData.setData.' );
    }
  }
  
  if( document.queryCommandSupported && document.queryCommandSupported("copy") ) {
    let temparea = document.createElement("textarea");
    temparea.textContent = csvData;
    temparea.style.position = "fixed";
    document.body.appendChild(temparea);
    temparea.select();
    try {
      const copysuccess = document.execCommand("copy");
      if( copysuccess )
      {
        console.log( 'Copied peak CSV using document.execCommand("copy").' );
        Wt.emit( $('.specviewer').attr('id'), {name:'miscSignal'},
                'peakCsvCopy-success-Copied CSV data to pasteboard as text.' );
        return;
      }else
      {
        console.warn( 'Failed to copy peak CSV using document.execCommand("copy").' );
        Wt.emit( $('.specviewer').attr('id'), {name:'miscSignal'},
                'peakCsvCopy-error-Failed to copy CSV data to pasteboard as text.' );
      }
    } catch( ex ) {
      console.warn( 'Caught exception copying peak CSV using document.execCommand("copy").' );
      Wt.emit( $('.specviewer').attr('id'), {name:'miscSignal'},
              'peakCsvCopy-error-Failed to copy peak CSV info to the pasteboard - sorry.' );
    } finally {
      document.body.removeChild( temparea );
    }
  }//if( document.queryCommandSupported && document.queryCommandSupported("copy") ) {
  
  console.log( 'Failed all methods to copy peak CSV.' );
  Wt.emit( $('.specviewer').attr('id'), {name:'miscSignal'},
              'peakCsvCopy-error-Unable to copy peak CSV info to the pasteboard - sorry.' );
}
);
}//namespace


  
  


PeakInfoDisplay::PeakInfoDisplay( InterSpec *viewer,
                                  D3SpectrumDisplayDiv *spectrumDisplayDiv,
                                  PeakModel *peakModel,
                                  Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_model( peakModel ),
    m_viewer( viewer ),
    m_spectrumDisplayDiv( spectrumDisplayDiv ),
    m_infoLayout( NULL ),
    m_infoView( NULL ),
    m_deletePeak( NULL ),
    m_searchForPeaks( NULL ),
    m_clearPeaksButton( nullptr ),
    m_nucFromRefButton( nullptr ),
    m_peakAddRemoveLabel( nullptr )
{
  assert( m_spectrumDisplayDiv );
  addStyleClass( "PeakInfoDisplay" );
  
  auto app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  app->useMessageResourceBundle( "PeakInfoDisplay" );
  
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
  
  /*
  auto make_dialog = [this](){
    InterSpec *viewer = InterSpec::instance();
    PeakInfoDisplay *display = viewer ? viewer->peakInfoDisplay() : nullptr;
    if( !display )
      return;
    
    SimpleDialog *window = new SimpleDialog( "Erase All Peaks?", "" );
    WPushButton *yes_button = window->addButton( "Yes" );
    WPushButton *no_button = window->addButton( "No" );
    
    yes_button->clicked().connect( boost::bind( &PeakInfoDisplay::removeAllPeaks, display ) );
  };
  
  make_dialog();
   */
  
  
  SimpleDialog *window = new SimpleDialog( WString::tr("pid-dialog-peak-erase"), "" );
  WPushButton *yes_button = window->addButton( WString::tr("Yes") );
  window->addButton( WString::tr("No") );
  
  yes_button->clicked().connect( boost::bind( &PeakInfoDisplay::removeAllPeaks, this ) );
  
  
  // We'll put an "undo" to close the dialog we just made - we could do a redo, but then making the
  //  undo becomes a bit harder... good enough, I guess
  auto closerfcn = wApp->bind( boost::bind( &WDialog::done, window, Wt::WDialog::DialogCode::Accepted ) );
  UndoRedoManager *undoManager = m_viewer ? m_viewer->undoRedoManager() : nullptr;
  if( undoManager )
    undoManager->addUndoRedoStep( closerfcn, nullptr, "Click clear all peaks." );
}//void confirmRemoveAllPeaks()


void PeakInfoDisplay::removeAllPeaks()
{
  if( !m_model )
    return;
  
  auto peaks_ptr = m_model->peaks();
  vector<shared_ptr<const PeakDef>> orig_peaks;
  if( peaks_ptr )
    orig_peaks.insert( end(orig_peaks), begin(*peaks_ptr), end(*peaks_ptr) );
  
  m_model->removeAllPeaks();
  
  auto undo = [orig_peaks](){
    InterSpec *viewer = InterSpec::instance();
    PeakModel *pmodel = viewer ? viewer->peakModel() : nullptr;
    if( pmodel )
      pmodel->setPeaks( orig_peaks );
  };
  auto redo = [](){
    InterSpec *viewer = InterSpec::instance();
    PeakModel *pmodel = viewer ? viewer->peakModel() : nullptr;
    if( pmodel )
      pmodel->removeAllPeaks();
  };
  
  UndoRedoManager *undoManager = m_viewer ? m_viewer->undoRedoManager() : nullptr;
  if( undoManager )
    undoManager->addUndoRedoStep( undo, redo, "Clear all peaks." );
}//void PeakInfoDisplay::removeAllPeaks()


void PeakInfoDisplay::assignNuclidesFromRefLines()
{
  PeakSearchGuiUtils::assign_peak_nuclides_from_reference_lines( m_viewer, true, false );
}//void assignNuclidesFromRefLines()


void PeakInfoDisplay::copyCsvPeakDataToClient()
{
  const auto meas = m_viewer->measurment(SpecUtils::SpectrumType::Foreground);
  const auto data = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  
  const string filename = meas ? meas->filename() : string();
  const deque<shared_ptr<const PeakDef>> &peaks = m_model->sortedPeaks();
  
  stringstream full_csv, no_hdr_csv, compact_csv;
  PeakModel::write_peak_csv( full_csv, filename, PeakModel::PeakCsvType::Full, peaks, data );
  PeakModel::write_peak_csv( no_hdr_csv, filename, PeakModel::PeakCsvType::NoHeader, peaks, data );
  PeakModel::write_peak_csv( compact_csv, filename, PeakModel::PeakCsvType::Compact, peaks, data );
  
  stringstream full_html, no_hdr_html, compact_html;
  PeakModel::write_peak_csv( full_html, filename, PeakModel::PeakCsvType::FullHtml, peaks, data );
  PeakModel::write_peak_csv( no_hdr_html, filename, PeakModel::PeakCsvType::NoHeaderHtml, peaks, data );
  PeakModel::write_peak_csv( compact_html, filename, PeakModel::PeakCsvType::CompactHtml, peaks, data );
  
  
  const string full_csv_str = Wt::WWebWidget::jsStringLiteral(full_csv.str(),'\'');
  const string no_hdr_csv_str = Wt::WWebWidget::jsStringLiteral(no_hdr_csv.str(),'\'');
  const string compact_csv_str = Wt::WWebWidget::jsStringLiteral(compact_csv.str(),'\'');
  
  const string full_html_str = Wt::WWebWidget::jsStringLiteral(full_html.str(),'\'');
  const string no_hdr_html_str = Wt::WWebWidget::jsStringLiteral(no_hdr_html.str(),'\'');
  const string compact_html_str = Wt::WWebWidget::jsStringLiteral(compact_html.str(),'\'');
  
  doJavaScript( "$(" + jsRef() + ")"
  ".data('CsvFullData'," + full_csv_str + ")"
  ".data('CsvFullDataHtml'," + full_html_str + ")"
  ".data('CsvNoHeaderData'," + no_hdr_csv_str + ")"
  ".data('CsvNoHeaderDataHtml'," + no_hdr_html_str + ")"
  ".data('CsvCompactData'," + compact_csv_str + ")"
  ".data('CsvCompactDataHtml'," + compact_html_str + ");");
}//void copyCsvPeakDataToClient()


void PeakInfoDisplay::removeCsvPeakDatafromClient()
{
  // We'll add a 1-second delay, to make sure everything is copied.
  doJavaScript( "setTimeout( function(){ $(" + jsRef() + ").data('CsvFullData',null)"
    ".data('CsvFullDataHtml',null)"
    ".data('CsvNoHeaderData',null)"
    ".data('CsvNoHeaderDataHtml',null)"
    ".data('CsvCompactData',null)"
    ".data('CsvCompactDataHtml',null); }, 1000 );" );
}//void removeCsvPeakDatafromClient()


void PeakInfoDisplay::enablePeakSearchButton( bool enable )
{
  m_searchForPeaks->setEnabled( enable );
}//void enablePeakSearchButton( bool enable = true )


void PeakInfoDisplay::handleChartLeftClick( const double energy )
{
  PeakModel::PeakShrdPtr peak = m_model->nearestPeak(energy);
  if( !peak )
    return;
  
  if( (energy < peak->lowerX()) || (energy > peak->upperX()) )
    return;
  
  Wt::WModelIndex index = m_model->indexOfPeak( peak );
  m_infoView->scrollTo( index, WTreeView::ScrollHint::EnsureVisible );
  m_infoView->select( index );
  
  enablePeakDelete( index );
}//void handleChartLeftClick( const double energy )


#if( InterSpec_PHONE_ROTATE_FOR_TABS )
void PeakInfoDisplay::setNarrowPhoneLayout( const bool narrow )
{
  // m_deletePeak //If mobile, will not have text, and minus_min_black.svg icon
  // addPeak //if mobile, will not have text, and plus_min_black.svg icon
  
  m_searchForPeaks->setText( narrow ? WString() : WString::tr("pid-search-peaks-btn") ); //Has magnifier.png as icon
  m_clearPeaksButton->setText( narrow ? WString() : WString::tr("pid-clear-peaks-btn") ); //By default, no icon
  m_clearPeaksButton->setIcon( narrow ? WLink("InterSpec_resources/images/sweep.svg") : WLink() );
  m_nucFromRefButton->setText( narrow ? WString() : WString::tr("pid-nuc-from-ref-btn") ); //assign_white.png as icon
  m_peakAddRemoveLabel->setHidden( narrow );
}//void setNarrowPhoneLayout( const bool narrow )
#endif


void PeakInfoDisplay::createNewPeak()
{
  /* ToDo:
   - PeakSearchGuiUtils::renderChartToSvg() seems to give a chart with a kinda large right margin; should look at.
   - Associating peak with nuclide.  Deal with?
   - Color?
   */
  auto meas = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  
  if( !meas || meas->num_gamma_channels() < 7 )
  {
    passMessage( WString::tr("pid-err-no-foreground"), WarningWidget::WarningMsgHigh );
    return;
  }//if( we dont have a valid foreground )
  
  float xmin = static_cast<float>( m_spectrumDisplayDiv->xAxisMinimum() );
  float xmax = static_cast<float>( m_spectrumDisplayDiv->xAxisMaximum() );
  
  if( meas )
  {
    const size_t nbin = meas->num_gamma_channels();
    xmin = std::max( xmin, meas->gamma_channel_lower(0) );
    xmax = std::min( xmax, meas->gamma_channel_upper(nbin-1) );
  }//if( meas )
  
  const float initialEnergy = 0.5f*(xmin + xmax);
  
  AddNewPeakDialog *window = new AddNewPeakDialog( initialEnergy );
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  
  // We'll make it so "undo" will close the window.
  //  Not implementing "redo", as it will get complicated, and not super useful
  InterSpec *viewer = InterSpec::instance();
  UndoRedoManager *undoManager = viewer ? viewer->undoRedoManager() : nullptr;
  if( undoManager )
  {
    auto closer = wApp->bind( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
    undoManager->addUndoRedoStep( closer, nullptr , "Open add peak dialog." );
  }//if( undoManager )
}//createNewPeak()


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
    m_deletePeak->show();
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
        m_deletePeak->show();
        return;
      }//if( we found a row selected
    }//for( loop over columns )
  }//for( loop over rows )
    
  m_deletePeak->disable(); //default
  if( dynamic_cast<WImage *>(m_deletePeak) )
    m_deletePeak->hide();
}//void disablePeakDelete()


void PeakInfoDisplay::handleSelectionChanged()
{
  m_infoView->closeEditors();
  
  const WModelIndexSet selected = m_infoView->selectedIndexes();
  if( selected.size() != 1 )
  {
    disablePeakDelete();
    return;
  }
  
  const WModelIndex index = *begin(selected);
  
  const auto peak = m_model->peak( index );
  if( peak )
  {
    m_spectrumDisplayDiv->highlightPeakAtEnergy( peak->mean() );
    enablePeakDelete( index );
  }
}//void handleSelectionChanged()


void PeakInfoDisplay::init()
{
  if( !m_model )
    throw runtime_error( "PeakInfoDisplay must be passed a valid PeakModel" );
  
  assert( !m_infoView );
  assert( !m_infoLayout );

  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_viewer );
  
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
  m_infoView->setColumnHidden( PeakModel::kUseForManualRelEff, true );
  m_infoView->setColumnHidden( PeakModel::kUseForShieldingSourceFit, true );

  m_infoView->setEditTriggers( WAbstractItemView::SingleClicked | WAbstractItemView::DoubleClicked );
  blurred().connect( boost::bind(&RowStretchTreeView::closeEditors, m_infoView, true) );
  
  ItemDelegate *dblDelagate = new ItemDelegate( m_infoView );
  dblDelagate->setTextFormat( "%.2f" );
  m_infoView->setItemDelegateForColumn( PeakModel::kMean, dblDelagate );
  m_infoView->setItemDelegateForColumn( PeakModel::kFwhm, dblDelagate );
  m_infoView->setItemDelegateForColumn( PeakModel::kLowerX, dblDelagate );
  m_infoView->setItemDelegateForColumn( PeakModel::kUpperX, dblDelagate );

  dblDelagate = new ItemDelegate( m_infoView );
  dblDelagate->setTextFormat( "%.0f" );
  //m_infoView->setItemDelegateForColumn( PeakModel::kAmplitude, dblDelagate ); // We'll actually return a string for peak amplitude
  m_infoView->setItemDelegateForColumn( PeakModel::kRoiCounts, dblDelagate );

  //tweak column widths
  m_infoView->setColumnWidth( PeakModel::kMean,       WLength(8, WLength::FontEx) );
  m_infoView->setColumnWidth( PeakModel::kFwhm,       WLength(9, WLength::FontEx) );
  m_infoView->setColumnWidth( PeakModel::kAmplitude,  WLength(12, WLength::FontEx) );
  m_infoView->setColumnWidth( PeakModel::kCps,        WLength(13, WLength::FontEx) );
  
  
  m_infoView->setColumnWidth( PeakModel::kIsotope,    WLength(9, WLength::FontEx) );
  m_infoView->setColumnWidth( PeakModel::kDifference,    WLength(8, WLength::FontEx) );
  //Note 20131211, wcjohns: closeOnBlur was previoulsy set to false, however
  //  there was a rare crash that happened with the model row was deleted while
  //  the editor was still open.  Hopefully closeOnBlur==true will fix this
  //  (appears to).
  const bool closeOnBlur = true;
  PhotopeakDelegate *nuclideDelegate = new PhotopeakDelegate( PhotopeakDelegate::NuclideDelegate, closeOnBlur, m_infoView );
  m_infoView->setItemDelegateForColumn( PeakModel::kIsotope, nuclideDelegate );

  m_infoView->setColumnWidth( PeakModel::kPhotoPeakEnergy, WLength(13 /*18*/, WLength::FontEx) );
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
  
  m_infoView->setColumnWidth( PeakModel::kLowerX,     WLength(/*11*/12, WLength::FontEx) );
  m_infoView->setColumnWidth( PeakModel::kUpperX,     WLength(/*11*/12, WLength::FontEx) );
  m_infoView->setColumnWidth( PeakModel::kRoiCounts,  WLength(8, WLength::FontEx) );
  m_infoView->setColumnWidth( PeakModel::kContinuumType,  WLength(/*11*/12, WLength::FontEx) );
//  m_infoView->setColumnWidth( PeakModel::kNumColumns, WLength(6, WLength::FontEx) );
//  m_infoView->setColumnAlignment( PeakModel::kMean, AlignRight );

//  Selection mode seems to have no effect....
//  m_infoView->setSelectionBehavior( Wt::SelectRows );
//  m_infoView->setSelectionMode( Wt::SingleSelection );
  m_infoView->setEditOptions( WAbstractItemView::SingleEditor );

  //m_infoView->clicked().connect( boost::bind( &PeakInfoDisplay::enablePeakDelete, this,
  //                                           boost::placeholders::_1 ) );
  //m_infoView->doubleClicked().connect( boost::bind( &PeakInfoDisplay::enablePeakDelete, this,
  //                                                 boost::placeholders::_1 ) );

  m_infoView->selectionChanged().connect( this, &PeakInfoDisplay::handleSelectionChanged );
  
  m_infoLayout->addWidget( m_infoView, 0, 0 );
  m_infoLayout->setRowStretch( 0, 1 );
  

  //Now add buttons to search/clear peaks
  WContainerWidget *bottomDiv = new WContainerWidget();
  bottomDiv->addStyleClass( "PeakInfoDisplayBottomDiv" );
  m_infoLayout->addWidget( bottomDiv, 1, 0 );

  
  auto helpBtn = new WContainerWidget( bottomDiv );
  helpBtn->addStyleClass( "Wt-icon ContentHelpBtn PeakInfoHlpBtn" );
  helpBtn->clicked().connect( boost::bind( &HelpSystem::createHelpWindow, "peak-manager" ) );
  
  
  WContainerWidget *buttonsDiv = new WContainerWidget( bottomDiv );
  buttonsDiv->addStyleClass( "PeakInfoDisplayButtonsDiv" );
  
  m_searchForPeaks = new WPushButton( WString::tr("pid-search-peaks-btn"), buttonsDiv );
  m_searchForPeaks->setIcon( "InterSpec_resources/images/magnifier.png" );
  
  HelpSystem::attachToolTipOn( m_searchForPeaks, WString::tr("pid-tt-search-peaks-btn"),
                              showToolTips, HelpSystem::ToolTipPosition::Top  );
  m_searchForPeaks->clicked().connect( boost::bind( &PeakSearchGuiUtils::automated_search_for_peaks, m_viewer, true ) );

  
  m_clearPeaksButton = new WPushButton( WString::tr("pid-clear-peaks-btn"), buttonsDiv );
  HelpSystem::attachToolTipOn( m_clearPeaksButton, WString::tr("pid-tt-clear-peaks-btn"),
                              showToolTips, HelpSystem::ToolTipPosition::Top  );
  
  //m_clearPeaksButton->setMargin(WLength(2),Wt::Left);
  m_clearPeaksButton->clicked().connect( this, &PeakInfoDisplay::confirmRemoveAllPeaks );
  m_clearPeaksButton->disable();

  //"Nuc. from Ref."
  m_nucFromRefButton = new WPushButton( WString::tr("pid-nuc-from-ref-btn"), buttonsDiv );
  m_nucFromRefButton->setIcon( "InterSpec_resources/images/assign_white.png" );
  
  //button->setMargin(WLength(2),Wt::Left|Wt::Right);
  HelpSystem::attachToolTipOn( m_nucFromRefButton, WString::tr("pid-tt-nuc-from-ref-btn"),
                              showToolTips , HelpSystem::ToolTipPosition::Top );
  m_nucFromRefButton->clicked().connect( boost::bind( &PeakInfoDisplay::assignNuclidesFromRefLines, this ) );
  m_nucFromRefButton->disable();
  
  auto enableDisableNucRef = [this](){
    const bool enable = (m_model->rowCount() > 0);
    m_clearPeaksButton->setEnabled( enable );
    m_nucFromRefButton->setEnabled( enable );
    //Should check if any reference lines are showing for m_nucFromRefButton as well...
  };
  
  m_model->dataChanged().connect( std::bind(enableDisableNucRef) );
  m_model->rowsRemoved().connect( std::bind(enableDisableNucRef) );
  m_model->rowsInserted().connect( std::bind(enableDisableNucRef) );
  m_model->layoutChanged().connect( std::bind(enableDisableNucRef) );
  
/*
  WPushButton *button = new WPushButton( "ID Nuclides", buttonsDiv );
  button->setIcon( "InterSpec_resources/images/assign_white.png" );
  button->setMargin(WLength(2),Wt::Left|Wt::Right);
  HelpSystem::attachToolTipOn( button,
                              "Guess nuclides responsible for peaks. Only applies to "
                              "peaks which do not already have a nuclide associated "
                              "with them.  Works best once all peaks have been fit for." ,
 showToolTips , HelpSystem::ToolTipPosition::Top );
  button->clicked().connect( boost::bind( &InterSpec::guessIsotopesForPeaks, m_viewer, (WApplication *)0 ) );
*/
  
  //If you want to post the below, so the ID isnt carried out in the main event loop, use the following
//  boost::function<void ()> guessIsotopeWorker = boost::bind( &InterSpec::guessIsotopesForPeaks, m_viewer, wApp );
//  button->clicked().connect( boost::bind( &WServer::post, WServer::instance(),
//                                           wApp->sessionId(), guessIsotopeWorker,
//                                           boost::function<void ()>() ) );
  
  m_peakAddRemoveLabel = new WLabel( WString("{1}: ").arg( WString::tr("Peak") ), buttonsDiv);
  m_peakAddRemoveLabel->addStyleClass("buttonSeparator");
  m_peakAddRemoveLabel->setMargin(WLength(10),Wt::Left);
  
  if( m_viewer->isMobile() )
  {
    WImage *addPeak = new WImage( WLink("InterSpec_resources/images/plus_min_black.svg"), buttonsDiv );
    addPeak->setAttributeValue( "width", "16" );
    addPeak->setAttributeValue( "height", "16" );
    addPeak->addStyleClass( "WhiteIcon" );
    addPeak->setMargin( 2, Wt::Left );
    addPeak->setMargin( 8, Wt::Right );
    addPeak->clicked().connect( this, &PeakInfoDisplay::createNewPeak );
    
    WImage *delPeak = new WImage( WLink("InterSpec_resources/images/minus_min_black.svg"), buttonsDiv );
    delPeak->setAttributeValue( "width", "16" );
    delPeak->setAttributeValue( "height", "16" );
    delPeak->addStyleClass( "WhiteIcon" );
    delPeak->clicked().connect( this, &PeakInfoDisplay::deleteSelectedPeak );
    m_deletePeak = delPeak;
    m_deletePeak->hide();
  }else
  {
    WPushButton *addPeak = new WPushButton( WString::tr("pid-add-peak-btn"), buttonsDiv );
    HelpSystem::attachToolTipOn( addPeak, WString::tr("pid-tt-add-peak-btn"),
                                showToolTips, HelpSystem::ToolTipPosition::Top );
    addPeak->clicked().connect( this, &PeakInfoDisplay::createNewPeak );
    addPeak->setIcon( "InterSpec_resources/images/plus_min_white.svg" );
    
    WPushButton *delPeak = new WPushButton( WString::tr("Delete"), buttonsDiv );
    HelpSystem::attachToolTipOn( delPeak, WString::tr("pid-tt-del-peaks"),
                                showToolTips, HelpSystem::ToolTipPosition::Top  );
    delPeak->setHiddenKeepsGeometry( true );
    delPeak->clicked().connect( this, &PeakInfoDisplay::deleteSelectedPeak );
    delPeak->setIcon( "InterSpec_resources/images/minus_min_white.png" );
    m_deletePeak = delPeak;
  }//if( mobile ) / else
  
  m_deletePeak->disable();
  
  
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
  if( !m_viewer->isPhone() )
  {
    const char *key = m_viewer->isMobile() ? "pid-better-editor-hint-mobile" : "pid-better-editor-hint";
    WText *txt = new WText( WString::tr(key), bottomDiv );
    txt->addStyleClass( "PeakEditHint" );
    const string show_js( "function(object, event){try{$('#" + txt->id() + "').css('visibility','visible');}catch(e){}}" );
    const string hide_js( "function(object, event){try{$('#" + txt->id() + "').css('visibility','hidden');}catch(e){}}" );
    mouseWentOver().connect( show_js );
    mouseWentOut().connect( hide_js );
    txt->doJavaScript( "(" + hide_js + ")();" );
  }
#endif
  
  WContainerWidget *csvDiv = new WContainerWidget( bottomDiv );
  csvDiv->addStyleClass( "PeakInfoDisplayCsvBtns" );
  
  WPushButton *copyButton = new WPushButton( csvDiv );
  copyButton->setIcon( "InterSpec_resources/images/copy_small.svg" );
  copyButton->setStyleClass( "LinkBtn" );
  copyButton->disable();
  HelpSystem::attachToolTipOn( copyButton, WString::tr("pid-tt-csv-copy"), showToolTips );
  WPopupMenu *copyMenu = nullptr;
  if( m_viewer->isMobile() )
    copyMenu = new WPopupMenu();
  else
    copyMenu = new PopupDivMenu( nullptr, PopupDivMenu::MenuType::TransientMenu);
  copyButton->setMenu( copyMenu );
  copyButton->clicked().connect( this, &PeakInfoDisplay::copyCsvPeakDataToClient );
  copyMenu->aboutToHide().connect( this, &PeakInfoDisplay::removeCsvPeakDatafromClient );
  
  
  LOAD_JAVASCRIPT(wApp, "PeakInfoDisplay.cpp", "PeakInfoDisplay", wtjsCopyPeakCsvDataToClipboard);
  WMenuItem *fullCopyMenuItem = copyMenu->addItem( WString::tr("pid-csv-copy-full") );
  WMenuItem *noHeaderCopyMenuItem = copyMenu->addItem( WString::tr("pid-csv-copy-no-hdr") );
  WMenuItem *compactCopyMenuItem = copyMenu->addItem( WString::tr("pid-csv-copy-compact") );
  
  fullCopyMenuItem->clicked().connect( "function(s,e){ "
    "Wt.WT.CopyPeakCsvDataToClipboard(s,e," + jsRef() + ", 'CsvFullData');"
  "}" );
  noHeaderCopyMenuItem->clicked().connect( "function(s,e){ "
    "Wt.WT.CopyPeakCsvDataToClipboard(s,e," + jsRef() + ", 'CsvNoHeaderData');"
  "}" );
  compactCopyMenuItem->clicked().connect( "function(s,e){ "
    "Wt.WT.CopyPeakCsvDataToClipboard(s,e," + jsRef() + ", 'CsvCompactData');"
  "}" );
  
  // TODO: Maybe add tool tips to each menu item?
  //HelpSystem::attachToolTipOn( fullCopyMenuItem, WString::tr("pid-tt-csv-export"), showToolTips );
  //HelpSystem::attachToolTipOn( noHeaderCopyMenuItem, WString::tr("pid-tt-csv-export"), showToolTips );
  //HelpSystem::attachToolTipOn( compactCopyMenuItem, WString::tr("pid-tt-csv-export"), showToolTips );
  
  
  WResource *csv = m_model->peakCsvResource();
#if( BUILD_AS_OSX_APP || IOS )
  WAnchor *csvButton = new WAnchor( WLink(csv), csvDiv );
  csvButton->setTarget( AnchorTarget::TargetNewWindow );
  csvButton->setStyleClass( "LinkBtn DownloadLink" );
#else
  WPushButton *csvButton = new WPushButton( csvDiv );
  csvButton->setIcon( "InterSpec_resources/images/download_small.svg" );
  csvButton->setLink( WLink(csv) );
  csvButton->setLinkTarget( Wt::TargetNewWindow );
  csvButton->setStyleClass( "LinkBtn DownloadBtn" );
  
#if( ANDROID )
  // Using hacked saving to temporary file in Android, instead of via network download of file.
  csvButton->clicked().connect( std::bind([csv](){
    android_download_workaround(csv, "photopeak_ref.csv");
  }) );
#endif //ANDROID
#endif //#if( BUILD_AS_OSX_APP || IOS ) / else
  
  csvButton->setText( WString::tr("CSV") );
  csvButton->disable();

  auto enableDisableCsv = [this,csvButton,copyButton](){
    const bool enable = (m_model->rowCount() > 0);
    csvButton->setEnabled( enable );
    copyButton->setEnabled( enable );
  };
  
  m_model->dataChanged().connect( std::bind(enableDisableCsv) );
  m_model->rowsRemoved().connect( std::bind(enableDisableCsv) );
  m_model->rowsInserted().connect( std::bind(enableDisableCsv) );
  m_model->layoutChanged().connect( std::bind(enableDisableCsv) );
  
  HelpSystem::attachToolTipOn( csvButton, WString::tr("pid-tt-csv-export"), showToolTips );
}//init()


