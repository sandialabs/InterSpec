#define _USE_MATH_DEFINES

#include "InterSpec_config.h"

#include <cmath>
#include <string>
#include <vector>
#include <numeric>
#include <assert.h>
#include <algorithm>

#include <Wt/WText>
#include <Wt/WTable>
#include <Wt/WLabel>
#include <Wt/WCheckBox>
#include <Wt/WResource>
#include <Wt/WLineEdit>
#include <Wt/WModelIndex>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WButtonGroup>
#include <Wt/WRadioButton>
#include <Wt/WApplication>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WRegExpValidator>
#include <Wt/WContainerWidget>
#include <Wt/WAbstractItemModel>
#include <Wt/WAbstractItemDelegate>

#include "InterSpec/FluxTool.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/HelpSystem.h"
#include "SpecUtils/Filesystem.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/RowStretchTreeView.h"

using namespace std;
using namespace Wt;


#if( FLUX_USE_COPY_TO_CLIPBOARD )
WT_DECLARE_WT_MEMBER
(CopyFluxDataTextToClipboard, Wt::JavaScriptFunction, "CopyFluxDataTextToClipboard",
 function( sender, event, id )
{
  var text = $('#'+id).data('TableData');
  if( !text )
    return false;
  
  var didcopy = 0;
  try
  {
    //This bit of code seems to work on Chrome, but not safari
    function listener(e) {
      e.clipboardData.setData("text/html", text);
      e.clipboardData.setData("text/plain", text);
      console.log( 'I think I copied it...' );
      didcopy = 1;
      e.preventDefault();
    }
    document.addEventListener("copy", listener);
    document.execCommand("copy");
    document.removeEventListener("copy", listener);
  
    if( didcopy )
      return didcopy;
  }catch(error){
    console.warn( 'Failed to copy richtext to copyboard' );
  }
  
  console.log( 'Will try to copy HTML to copyboard' );
  
  //We failed to copy richtext, lets just copy the HTML as text.
  //  ToDo: We could probably try the clipboard API to copy formated text
  //        See https://developer.mozilla.org/en-US/docs/Web/API/Clipboard
  
  if( window.clipboardData && window.clipboardData.setData ) {
    return clipboardData.setData("Text", text);  // IE
  }else if( document.queryCommandSupported && document.queryCommandSupported("copy") ) {
    var temparea = document.createElement("textarea");
    temparea.textContent = text;
    temparea.style.position = "fixed";
    document.body.appendChild(temparea);
    temparea.select();
    try {
      var copysuccess = document.execCommand("copy");
      return copysuccess ? 2 : 0;
    } catch( ex ) {
      console.warn("Copy to clipboard failed.", ex);
      return 0;
    } finally {
      document.body.removeChild( temparea );
    }
  } else {
    return 0;
  }
}
);
#endif //FLUX_USE_COPY_TO_CLIPBOARD


namespace FluxToolImp
{
  struct index_compare_sort
  {
    const std::vector<std::array<double,FluxToolWidget::FluxColumns::FluxNumColumns>> &arr;
    const std::vector<std::string> &names;
    int m_sortColumn;
    Wt::SortOrder m_sortOrder;
    
    index_compare_sort(const std::vector<std::array<double,FluxToolWidget::FluxColumns::FluxNumColumns>> &arr,
                       const std::vector<std::string> &nucnames,
                       int col, Wt::SortOrder order)
      : arr(arr),
        names( nucnames ),
        m_sortColumn( col ),
        m_sortOrder( order )
    {
    }
    
    bool operator()(const size_t a, const size_t b) const
    {
      bool less;
      if( m_sortColumn == FluxToolWidget::FluxColumns::FluxNuclideCol )
        less = names[a] < names[b];
      else
        less = arr[a][m_sortColumn] < arr[b][m_sortColumn];
      return (m_sortOrder==Wt::SortOrder::AscendingOrder ? less : !less);
    }
  };//struct index_compare
  
  
  bool showFluxColumn( FluxToolWidget::DisplayInfoLevel disptype, const FluxToolWidget::FluxColumns col )
  {
    switch( disptype )
    {
      case FluxToolWidget::DisplayInfoLevel::Simple:
        switch( col )
        {
          case FluxToolWidget::FluxColumns::FluxEnergyCol:
          case FluxToolWidget::FluxColumns::FluxNuclideCol:
          case FluxToolWidget::FluxColumns::FluxGammasInto4PiCol:
            return true;
            break;
          
          case FluxToolWidget::FluxColumns::FluxPeakCpsCol:
          case FluxToolWidget::FluxColumns::FluxFluxPerCm2PerSCol:
          case FluxToolWidget::FluxColumns::FluxIntrinsicEffCol:
          case FluxToolWidget::FluxColumns::FluxGeometricEffCol:
          case FluxToolWidget::FluxColumns::FluxFluxOnDetCol:
          case FluxToolWidget::FluxColumns::FluxNumColumns:
            return false;
            break;
        }//switch( col )
        
        break;
        
      case FluxToolWidget::DisplayInfoLevel::Normal:
        switch( col )
        {
          case FluxToolWidget::FluxColumns::FluxEnergyCol:
          case FluxToolWidget::FluxColumns::FluxPeakCpsCol:
          case FluxToolWidget::FluxColumns::FluxFluxPerCm2PerSCol:
          case FluxToolWidget::FluxColumns::FluxGammasInto4PiCol:
            return true;
            break;
            
          case FluxToolWidget::FluxColumns::FluxNuclideCol:
          case FluxToolWidget::FluxColumns::FluxIntrinsicEffCol:
          case FluxToolWidget::FluxColumns::FluxGeometricEffCol:
          case FluxToolWidget::FluxColumns::FluxFluxOnDetCol:
          case FluxToolWidget::FluxColumns::FluxNumColumns:
            return false;
            break;
        }//switch( col )
        
        break;
        
      case FluxToolWidget::DisplayInfoLevel::Extended:
        return (col != FluxToolWidget::FluxColumns::FluxNumColumns); //e.g., always true
        break;
    }//switch( disptype )
    
    
    return false;
  }//bool showFluxColumn( const FluxToolWidget::FluxColumns col )
  
  std::string print_value( const double value, const double uncert )
  {
    int nsigfigs = 6;
    char buffer[256] = { '\0' };
    
    const bool haveUncert = ((uncert > std::numeric_limits<double>::epsilon())
                             && (value > 0.0) && !(std::isinf)(value) && !(std::isnan)(value) );
    if( haveUncert )
    {
      //We'll be conservative and add at least one extra sig fig by using the nsigfigs to be
      //  number of digits after the decimal when printed in scientific notation.
      nsigfigs = static_cast<int>( std::round( std::log10( value / uncert ) + 0.5 ) );
      nsigfigs = std::min( std::max( nsigfigs, 3 ), 9 ); //clamp num sig figs between 3 and 9
    }//if( haveUncert )
    
    snprintf( buffer, sizeof(buffer), "%.*G", nsigfigs, value );
    
    return buffer;
  }//print_value(...)


  std::string print_uncert( const double value, const double uncert, const bool uncertAsPercent )
  {
    const bool haveUncert = ((uncert > std::numeric_limits<double>::epsilon())
                             && (value > 0.0) && !(std::isinf)(value) && !(std::isnan)(value) );
    
    if( !haveUncert )
      return "";
    
    int nsigfigs = 6;
    char buffer[256] = { '\0' };
    
    
    //We'll be conservative and add at least one extra sig fig by using the nsigfigs to be
    //  number of digits after the decimal when printed in scientific notation.
    nsigfigs = static_cast<int>( std::round( std::log10( value / uncert ) + 0.5 ) );
    nsigfigs = std::min( std::max( nsigfigs, 3 ), 9 ); //clamp num sig figs between 3 and 9
    
    if( uncertAsPercent )
    {
      const double percentUncert = 100.0 * uncert / value;
      
      if( percentUncert > 10 )
        snprintf( buffer, sizeof(buffer), "%.1f%%", percentUncert );
      else if( percentUncert > 1 )
        snprintf( buffer, sizeof(buffer), "%.2f%%", percentUncert );
      else
        snprintf( buffer, sizeof(buffer), "%.*G%%", nsigfigs, percentUncert );
    }else
    {
      snprintf( buffer, sizeof(buffer), "%.*G", nsigfigs, uncert );
    }
    
  return buffer;
}//print_uncert(...)


  class FluxModel : public  Wt::WAbstractItemModel
  {
  protected:
    FluxToolWidget *m_fluxtool;
    FluxToolWidget::FluxColumns m_sortColumn;
    Wt::SortOrder m_sortOrder;
    
    //Depending on sort column and order, we will map the rows of
    //  m_fluxtool->m_data to the row in this model
    vector<size_t> m_sort_indices;
    
  public:
    FluxModel( FluxToolWidget *parent )
    : WAbstractItemModel( parent ),
      m_fluxtool( parent ),
      m_sortColumn( FluxToolWidget::FluxColumns::FluxEnergyCol ),
      m_sortOrder( Wt::SortOrder::AscendingOrder )
    {
      parent->tableUpdated().connect( this, &FluxModel::handleFluxToolWidgetUpdated );
    }
    
    
    virtual ~FluxModel()
    {
    }
    
    
    virtual Wt::WFlags<Wt::ItemFlag> flags( const Wt::WModelIndex &index ) const
    {
      return Wt::WFlags<Wt::ItemFlag>(Wt::ItemFlag::ItemIsXHTMLText);
    }
    
    
    virtual WFlags<Wt::HeaderFlag> headerFlags( int section, Wt::Orientation orientation = Wt::Horizontal ) const
    {
      return Wt::WFlags<Wt::HeaderFlag>(Wt::HeaderFlag::HeaderIsXHTMLText);
    }
    
    
    void handleFluxToolWidgetUpdated()
    {
      if( !m_sort_indices.empty() )
      {
        beginRemoveRows( WModelIndex(), 0, static_cast<int>(m_sort_indices.size() - 1) );
        m_sort_indices.clear();
        endRemoveRows();
      }
      
      if( !m_fluxtool->m_data.empty() )
      {
        beginInsertRows( WModelIndex(), 0, static_cast<int>(m_fluxtool->m_data.size() - 1) );
        doSortWork();
        endInsertRows();
      }//
    }//void handleFluxToolWidgetUpdated()
    
    
    virtual boost::any data( const Wt::WModelIndex &index,
                             int role = Wt::DisplayRole ) const
    {
      if( index.parent().isValid() )
        return boost::any();
      
      const int row = index.row();
      const int col = index.column();
      
      if( row < 0 || col < 0
         || col >= FluxToolWidget::FluxColumns::FluxNumColumns
         || row >= static_cast<int>(m_sort_indices.size()) )
        return boost::any();
      
      const size_t realind = m_sort_indices[row];
      
      const auto &nucs = m_fluxtool->m_nucNames;
      const auto &data = m_fluxtool->m_data;
      const auto &uncerts = m_fluxtool->m_uncertainties;
      
      assert( nucs.size() == data.size() );
      assert( uncerts.size() == data.size() );
      assert( realind < data.size() );
      
      switch( role )
      {
        case Wt::DisplayRole:
          if( col == FluxToolWidget::FluxColumns::FluxNuclideCol )
            return nucs[realind].empty() ? boost::any() : boost::any( WString::fromUTF8(nucs[realind]) );
          return boost::any( data[realind][col] );
          
        //case Wt::StyleClassRole: break;
        //case Wt::ToolTipRole: break;
        case Wt::UserRole:
          if( col == FluxToolWidget::FluxColumns::FluxNuclideCol )
            return boost::any();
          return uncerts[realind][col] > std::numeric_limits<double>::epsilon() ? boost::any(uncerts[realind][col]) : boost::any();

        default:
          return boost::any();
      }//switch( role )
      
      return boost::any();
    }//data(...)
    
    
    boost::any headerData( int section,
                          Wt::Orientation orientation = Wt::Horizontal,
                          int role = Wt::DisplayRole ) const
    {
      if( section < 0 || section >= FluxToolWidget::FluxColumns::FluxNumColumns
         || orientation != Wt::Horizontal
         || role != Wt::DisplayRole )
        return boost::any();
      return boost::any( m_fluxtool->m_colnames[section] );
    }
    
    virtual int columnCount( const Wt::WModelIndex &parent
                            = Wt::WModelIndex() ) const
    {
      if( parent.isValid() )
        return 0;
      return FluxToolWidget::FluxColumns::FluxNumColumns;
    }
    
    virtual int rowCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const
    {
      if( parent.isValid() )
        return 0;
      return static_cast<int>( m_sort_indices.size() );
    }
    
    virtual Wt::WModelIndex parent( const Wt::WModelIndex &index) const
    {
      return WModelIndex();
    }
    
    virtual Wt::WModelIndex index( int row, int column,
                                  const Wt::WModelIndex &parent = Wt::WModelIndex() ) const
    {
      return WAbstractItemModel::createIndex(row, column, nullptr);
    }
    
    void doSortWork()
    {
      assert( m_fluxtool->m_data.size() == m_fluxtool->m_nucNames.size() );
      
      m_sort_indices.resize( m_fluxtool->m_data.size() );
      
      for( size_t i = 0; i < m_sort_indices.size(); ++i )
        m_sort_indices[i] = i;
      
      std::stable_sort( m_sort_indices.begin(), m_sort_indices.end(),
                       index_compare_sort( m_fluxtool->m_data, m_fluxtool->m_nucNames, m_sortColumn, m_sortOrder) );
    }//void doSortWork()
    
    virtual void sort( int column, Wt::SortOrder order = Wt::AscendingOrder )
    {
      m_sortOrder = order;
      m_sortColumn = FluxToolWidget::FluxColumns( column );
      
      layoutAboutToBeChanged().emit();
      doSortWork();
      layoutChanged().emit();
    }//sort(...)
  };//class FluxModel(...)
  
  
  class FluxRenderDelegate : public WAbstractItemDelegate
  {
    FluxToolWidget *m_fluxTool;
  public:
    FluxRenderDelegate( FluxToolWidget *parent )
    : WAbstractItemDelegate( parent ),
      m_fluxTool( parent )
    {
    }
    
    
    virtual WWidget *update( WWidget *widget, const WModelIndex& index, WFlags<ViewItemRenderFlag> flags )
    {
      auto oldwidget = dynamic_cast<WText *>( widget );
      
      if( oldwidget )
      {
        oldwidget->setText("");
      }else
      {
        if( widget )
          delete widget;
        oldwidget = new WText();
      }
      
      
      auto model = index.model();
      if( !model )
        return oldwidget;
      
      auto data = index.data();
      auto uncertdata = model->data( index.row(), index.column(), Wt::UserRole );
      
      if( data.empty() )
        return oldwidget;
      
      WString valstr;
      if( index.column() == FluxToolWidget::FluxColumns::FluxNuclideCol )
      {
        valstr = Wt::asString(data);
      }else
      {
        const bool hasUncert = !uncertdata.empty();
        
        double val = 0.0, uncert = 0.0;
        try
        {
          val = boost::any_cast<double>(data);
          if( hasUncert )
            uncert = boost::any_cast<double>(uncertdata);
        }catch(...)
        {
          return oldwidget;
        }
      
        
        if( index.column() == FluxToolWidget::FluxColumns::FluxEnergyCol )
        {
          char buffer[64];
          snprintf( buffer, sizeof(buffer), "%.2f", val );
          valstr = WString::fromUTF8(buffer);
        }else
        {
          if( hasUncert )
          {
            const bool uncertAsPercent = (m_fluxTool
                  && (m_fluxTool->displayInfoLevel() == FluxToolWidget::DisplayInfoLevel::Simple));
            
            const string tmpstr = print_value( val, uncert )
                                  + " &plusmn; "
                                  + print_uncert( val, uncert, uncertAsPercent );
            
            valstr = WString::fromUTF8( tmpstr );
          }else
          {
            valstr = WString::fromUTF8( print_value( val, -1.0 ) );
          }
        }//if( index.column() == FluxToolWidget::FluxColumns::FluxEnergyCol )
      }//if( nucname ) / else
      
      
      oldwidget->setText( valstr );
      return oldwidget;
    }//WWidget *update(...)
  };//class FluxRenderDelegate


  
  class FluxCsvResource : public Wt::WResource
  {
  protected:
    FluxToolWidget *m_fluxtool;
    Wt::WApplication *m_app;
    
  public:
    FluxCsvResource( FluxToolWidget *parent )
    : WResource( parent ),
      m_fluxtool( parent ),
      m_app( WApplication::instance() )
    {
      assert( m_app );
    }
    
    virtual ~FluxCsvResource()
    {
      beingDeleted();
    }
    
    static void data_to_strm( FluxToolWidget *tool, std::ostream &strm, const bool html, const FluxToolWidget::DisplayInfoLevel disptype )
    {
      const string eol_char = html ? "\\n" : "\r\n"; //for windows - could potentially customize this for the users operating system
      
      if( html )
        strm << "<table border=\"1\" cellpadding=\"2\" style=\"border-collapse: collapse\">" << eol_char;
      
      for( FluxToolWidget::FluxColumns col = FluxToolWidget::FluxColumns(0);
          col < FluxToolWidget::FluxColumns::FluxNumColumns; col = FluxToolWidget::FluxColumns(col+1) )
      {
        const WString &colname = tool->m_colnamesCsv[col];
        
        if( !FluxToolImp::showFluxColumn(disptype, col) )
          continue;
        
        if( html )
          strm << (col==0 ? "\\t<tr><th>" : "</th><th>") << colname;
        else
          strm << (col==0 ? "" : ",") << colname;
        
        //No uncertainty on energy, effs, or nuc
        switch( col )
        {
          case FluxToolWidget::FluxEnergyCol:
          case FluxToolWidget::FluxIntrinsicEffCol:
          case FluxToolWidget::FluxGeometricEffCol:
          case FluxToolWidget::FluxNuclideCol:
          case FluxToolWidget::FluxNumColumns:
            break;
            
          case FluxToolWidget::FluxPeakCpsCol:
          case FluxToolWidget::FluxFluxOnDetCol:
          case FluxToolWidget::FluxFluxPerCm2PerSCol:
          case FluxToolWidget::FluxGammasInto4PiCol:
            switch( disptype )
            {
              case FluxToolWidget::DisplayInfoLevel::Simple:
                if( html )
                  strm << "</th><th> Uncert (%)";
                else
                  strm << ", Uncert (%)";
                break;
              
              case FluxToolWidget::DisplayInfoLevel::Normal:
              case FluxToolWidget::DisplayInfoLevel::Extended:
                if( html )
                  strm << "</th><th>" << (colname + " Uncertainty");
                else
                  strm << "," << (colname + " Uncertainty");
                break;
            }//switch( disptype )
            
            break;
        }//switch( col )
      }//for( loop over columns )
      
      if( html )
        strm << "</th></tr>";
      strm << eol_char;
      
      for( size_t row = 0; row < tool->m_data.size(); ++row )
      {
        if( html )
          strm << "\\t<tr>";
        
        for( FluxToolWidget::FluxColumns col = FluxToolWidget::FluxColumns(0);
            col < FluxToolWidget::FluxColumns::FluxNumColumns; col = FluxToolWidget::FluxColumns(col+1) )
        {
          if( !FluxToolImp::showFluxColumn(disptype, col) )
            continue;
          
          const double data = tool->m_data[row][col];
          const double uncert = tool->m_uncertainties[row][col];
          
          
          std::string datastr;
          
          switch( col )
          {
            case FluxToolWidget::FluxEnergyCol:
            {
              char buffer[128] = { '\0' };
              snprintf( buffer, sizeof(buffer), "%.2f", data );
              datastr = buffer;
              break;
            }
              
            case FluxToolWidget::FluxNuclideCol:
              datastr = tool->m_nucNames[row];
              if( !html )
                SpecUtils::ireplace_all( datastr, ",", " " );
              break;
              
            case FluxToolWidget::FluxPeakCpsCol:
            case FluxToolWidget::FluxIntrinsicEffCol:
            case FluxToolWidget::FluxGeometricEffCol:
            case FluxToolWidget::FluxFluxOnDetCol:
            case FluxToolWidget::FluxFluxPerCm2PerSCol:
            case FluxToolWidget::FluxGammasInto4PiCol:
            {
              datastr = print_value( data, uncert );
              break;
            }
              
            case FluxToolWidget::FluxNumColumns: //wont get here, but whatever
              break;
          }//switch( col )
          
          
          if( html )
            strm << (col==0 ? "<td>" : "</td><td>") << datastr;
          else
            strm << (col==0 ? "" : ",") << datastr;
          
          switch( col )
          {
            case FluxToolWidget::FluxEnergyCol:
            case FluxToolWidget::FluxIntrinsicEffCol:
            case FluxToolWidget::FluxGeometricEffCol:
            case FluxToolWidget::FluxNuclideCol:
            case FluxToolWidget::FluxNumColumns:
              break;
              
            case FluxToolWidget::FluxPeakCpsCol:
            case FluxToolWidget::FluxFluxOnDetCol:
            case FluxToolWidget::FluxFluxPerCm2PerSCol:
            case FluxToolWidget::FluxGammasInto4PiCol:
            {
              switch( disptype )
              {
                case FluxToolWidget::DisplayInfoLevel::Simple:
                  strm << (html ? "</td><td>" : ",") << print_uncert( data, uncert, true );
                  break;
                  
                case FluxToolWidget::DisplayInfoLevel::Normal:
                case FluxToolWidget::DisplayInfoLevel::Extended:
                  strm << (html ? "</td><td>" : ",") << print_uncert( data, uncert, false );
                  break;
              }//switch( disptype )
              
              break;
            }//case CPS, FluxOnDet, FluxPerArea, GammasInto 4pi
          }//switch( col )
        }//for( loop over columns )
        
        if( html )
          strm << "</td></tr>";
        strm << eol_char;
      }//for( size_t row = 0; row < m_fluxtool->m_data.size(); ++row )
      
      if( html )
        strm << "</table>";
    }//static void data_to_strm( FluxToolWidget *tool, std::ostream &strm, const bool html )
    
    
  private:
    virtual void handleRequest( const Wt::Http::Request &, Wt::Http::Response &response )
    {
      WApplication::UpdateLock lock( m_app );
      
      if( !lock )
      {
        log("error") << "Failed to WApplication::UpdateLock in FluxRenderDelegate.";
        
        response.out() << "Error grabbing application lock to form FluxRenderDelegate resource; please report to InterSpec@sandia.gov.";
        response.setStatus(500);
        assert( 0 );
        
        return;
      }//if( !lock )
      
      const string eol_char = "\r\n"; //for windows - could potentially cosutomize this for the users operating system
      
      if( !m_fluxtool )
        return;
      
      string filename;
      auto meas = m_fluxtool->m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
      if( meas && !meas->filename().empty() )
      {
        filename = SpecUtils::filename( meas->filename() );
        const string extension = SpecUtils::file_extension( filename );
        filename = filename.substr( 0, filename.size() - extension.size() );
        if( filename.size() )
          filename += "_";
      }
      filename += "flux.csv";
      
      suggestFileName( filename, WResource::Attachment );
      
      FluxToolWidget::DisplayInfoLevel disptype = FluxToolWidget::DisplayInfoLevel::Extended;
      if( m_fluxtool )
        disptype = m_fluxtool->m_displayInfoLevel;
      
      data_to_strm( m_fluxtool, response.out(), false, disptype );
    }//handleRequest(...)
    
  };//class FluxCsvResource
}//namespace


FluxToolWindow::FluxToolWindow( InterSpec *viewer )
: AuxWindow( "Flux Tool",
  (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
   | AuxWindowProperties::SetCloseable
   | AuxWindowProperties::DisableCollapse) ),
  m_fluxTool( nullptr )
{
  assert( viewer );
  rejectWhenEscapePressed( true );
  
  m_fluxTool = new FluxToolWidget( viewer );
  //m_fluxTool->setHeight( WLength(100, WLength::Percentage) );
  
  stretcher()->addWidget( m_fluxTool, 0, 0 );
  
  AuxWindow::addHelpInFooter( footer(), "flux-tool" );
  
  
  WContainerWidget *buttonDiv = footer();
  
  WResource *csv = new FluxToolImp::FluxCsvResource( m_fluxTool );
#if( BUILD_AS_OSX_APP )
  WAnchor *csvButton = new WAnchor( WLink(csv), buttonDiv );
  csvButton->setTarget( AnchorTarget::TargetNewWindow );
  csvButton->setStyleClass( "LinkBtn" );
#else
  WPushButton *csvButton = new WPushButton( buttonDiv );
  csvButton->setIcon( "InterSpec_resources/images/download_small.svg" );
  csvButton->setLink( WLink(csv) );
  csvButton->setLinkTarget( Wt::TargetNewWindow );
  csvButton->setStyleClass( "LinkBtn DownloadBtn" );
#endif
  
  csvButton->setText( "CSV" );
  //csvButton->setAttributeValue( "style", "float: none;" );  //Keep the CSV download to the left side of the close button.  .LinkBtn style class has the button float right..
  csvButton->setFloatSide( Wt::Side::Left );
  
  auto enableDisableCsv = [csvButton,this](){
    csvButton->setDisabled( m_fluxTool->m_data.empty() );
  };
  
  m_fluxTool->m_tableUpdated.connect( std::bind(enableDisableCsv) );
  csvButton->disable();
  
  
  WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( this, &AuxWindow::hide );
  
  finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );
  
  show();
  
  const int screenW = viewer->renderedWidth();
  const int screenH = viewer->renderedHeight();
  const int width = ((screenW < 680) ? screenW : 680);
  const int height = ((screenH < 420) ? screenH : 420);
  resizeWindow( width, height );
  
  resizeToFitOnScreen();
  centerWindowHeavyHanded();
}//FluxToolWindow(...) constrctor


FluxToolWindow::~FluxToolWindow()
{
}


FluxToolWidget::FluxToolWidget( InterSpec *viewer, Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_interspec( viewer ),
    m_detector( nullptr ),
    m_msg( nullptr ),
    m_distance( nullptr ),
    m_table( nullptr ),
#if( FLUX_USE_COPY_TO_CLIPBOARD )
    m_copyBtn( nullptr ),
    m_infoCopied( this, "infocopied", true ),
#endif
    m_needsTableRefresh( true ),
    m_displayInfoLevel( DisplayInfoLevel::Normal ),
    m_tableUpdated( this )
{
  init();
}


FluxToolWidget::~FluxToolWidget()
{
}


Wt::Signal<> &FluxToolWidget::tableUpdated()
{
  return m_tableUpdated;
}

FluxToolWidget::DisplayInfoLevel FluxToolWidget::displayInfoLevel() const
{
  return m_displayInfoLevel;
}

void FluxToolWidget::init()
{
  assert( m_interspec );
  assert( !m_detector );
  
  for( FluxColumns col = FluxColumns(0); col < FluxColumns::FluxNumColumns; col = FluxColumns(col + 1) )
  {
    switch( col )
    {
      case FluxEnergyCol:
        m_colnames[col] = WString::fromUTF8("Energy (keV)");
        m_colnamesCsv[col] = WString::fromUTF8("Energy (keV)");
        break;
      case FluxNuclideCol:
        m_colnames[col] = WString::fromUTF8("Nuclide");
        m_colnamesCsv[col] = WString::fromUTF8("Nuclide");
        break;
      case FluxPeakCpsCol:
        m_colnames[col] = WString::fromUTF8("Peak CPS");
        m_colnamesCsv[col] = WString::fromUTF8("Peak CPS");
        break;
      case FluxIntrinsicEffCol:
        m_colnames[col] = WString::fromUTF8("Intr. Eff.");
        m_colnamesCsv[col] = WString::fromUTF8("Intrinsic Efficiency");
        break;
      case FluxGeometricEffCol:
        m_colnames[col] = WString::fromUTF8("Geom. Eff.");
        m_colnamesCsv[col] = WString::fromUTF8("Geometric Efficiency");
        break;
      case FluxFluxOnDetCol:
        m_colnames[col] = WString::fromUTF8("Flux on Det. (&gamma;/s)");
        m_colnamesCsv[col] = WString::fromUTF8("Flux on Detector (gammas/s)");
        break;
      case FluxFluxPerCm2PerSCol:
        m_colnames[col] = WString::fromUTF8("Flux (&gamma;/cm&sup2;/s)");
        m_colnamesCsv[col] = WString::fromUTF8("Flux (gammas/cm2/s)");
        break;
      case FluxGammasInto4PiCol:
        m_colnames[col] = WString::fromUTF8("&gamma;/4&pi;/s");
        m_colnamesCsv[col] = WString::fromUTF8("gammas/4pi/s");
        break;
      case FluxNumColumns:        break;
    }//switch( col )
  }//for( loop over columns )
  
  
  wApp->useStyleSheet( "InterSpec_resources/FluxTool.css" );
  
  const bool showToolTips = m_interspec ? InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_interspec ) : false;
  
  addStyleClass( "FluxToolWidget" );
  
  WGridLayout *layout = new WGridLayout();
  layout->setContentsMargins( 0, 0, 0, 0 );
  setLayout( layout );
  
  WTable *distDetRow = new WTable();
  distDetRow->addStyleClass( "FluxDistMsgDetTbl" );
#if( FLUX_USE_COPY_TO_CLIPBOARD )
  layout->addWidget( distDetRow, 0, 0, 1, 2 );
#else
  layout->addWidget( distDetRow, 0, 0 );
#endif
  
  SpectraFileModel *specFileModel = m_interspec->fileManager()->model();
  m_detector = new DetectorDisplay( m_interspec, specFileModel );
  m_detector->addStyleClass( "FluxDet" );
  m_interspec->detectorChanged().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  m_interspec->detectorModified().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  distDetRow->elementAt(0,1)->addWidget( m_detector );
  
  auto distCell = distDetRow->elementAt(0,0);
  distCell->addStyleClass( "FluxDistCell" );
  WLabel *label = new WLabel( "Distance:", distCell );
  label->addStyleClass( "FluxDistLabel" );
  
  m_distance = new WLineEdit( "100 cm", distCell );
  m_distance->addStyleClass( "FluxDistanceEnter" );
  label->setBuddy(m_distance);
  WRegExpValidator *validator = new WRegExpValidator( PhysicalUnits::sm_distanceUnitOptionalRegex, this );
  validator->setFlags( Wt::MatchCaseInsensitive );
  m_distance->setValidator( validator );
  HelpSystem::attachToolTipOn( m_distance,
                              "Distance from center of source to face of detector. Number must be"
                              " followed by units; valid units are: meters, m, cm, mm, km, feet,"
                              " ft, ', in, inches, or \".  You may also add multiple distances,"
                              " such as '3ft 4in', or '3.6E-2 m 12 cm' which are equivalent to "
                              " 40inches and 15.6cm respectively.", showToolTips );
  m_distance->changed().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  m_distance->enterPressed().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  
  
  PeakModel *peakmodel = m_interspec->peakModel();
  assert( peakmodel );
  peakmodel->dataChanged().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  peakmodel->rowsRemoved().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  peakmodel->rowsInserted().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  peakmodel->layoutChanged().connect( this, &FluxToolWidget::setTableNeedsUpdating );

  auto msgCell = distDetRow->elementAt(1,0);
  msgCell->setColumnSpan( 2 );
  msgCell->addStyleClass( "FluxMsgCell" );
  m_msg = new WText( "", Wt::XHTMLText, msgCell );
  m_msg->addStyleClass( "FluxMsg" );
  
  
  FluxToolImp::FluxModel *fluxmodel = new FluxToolImp::FluxModel( this );
  m_table = new RowStretchTreeView();
  m_table->setModel( fluxmodel );
  m_table->setRootIsDecorated( false );
  m_table->addStyleClass( "FluxTable" );
  m_table->setAlternatingRowColors( true );
  m_table->sortByColumn( FluxColumns::FluxEnergyCol, Wt::AscendingOrder );
  
  //Setting rows selectable doesnt seem to work... not that it would do us
  //  any good anyway since you probably want to copy the info to the clipboard.
  //m_table->setSelectable( true );
  //m_table->setSelectionMode( Wt::SelectionMode::SingleSelection );
  
  FluxToolImp::FluxRenderDelegate *renderdel = new FluxToolImp::FluxRenderDelegate( this );
  m_table->setItemDelegate( renderdel );
  
  for( FluxColumns col = FluxColumns(0); col < FluxColumns::FluxNumColumns; col = FluxColumns(col + 1) )
    m_table->setSortingEnabled( col, (col!=FluxColumns::FluxGeometricEffCol) );
  
#if( FLUX_USE_COPY_TO_CLIPBOARD )
  layout->addWidget( m_table, 1, 0, 1, 2 );
#else
  layout->addWidget( m_table, 1, 0, 1, 1 );
#endif
  layout->setRowStretch( 1, 1 );
  
  WContainerWidget *buttonBox = new WContainerWidget();
  buttonBox->addStyleClass( "FluxInfoAmount" );
  
  WRadioButton *simpleInfo = new WRadioButton( "Simple", buttonBox );
  WRadioButton *standardInfo = new WRadioButton( "Standard", buttonBox );
  WRadioButton *moreInfo = new WRadioButton( "More", buttonBox );
  
  WButtonGroup *btnGrp = new WButtonGroup( buttonBox );
  btnGrp->addButton( simpleInfo, static_cast<int>(DisplayInfoLevel::Simple) );
  btnGrp->addButton( standardInfo, static_cast<int>(DisplayInfoLevel::Normal) );
  btnGrp->addButton( moreInfo, static_cast<int>(DisplayInfoLevel::Extended) );
  btnGrp->setCheckedButton( standardInfo );

  btnGrp->checkedChanged().connect( std::bind( [this,btnGrp](){
    const auto level = static_cast<DisplayInfoLevel>( btnGrp->checkedId() );
    
    assert( (level == DisplayInfoLevel::Simple)
           || (level == DisplayInfoLevel::Normal)
           || (level == DisplayInfoLevel::Extended) );
    
    setDisplayInfoLevel( level );
  }) );
  
  
#if( FLUX_USE_COPY_TO_CLIPBOARD )
  layout->addWidget( buttonBox, 2, 1, AlignRight | AlignMiddle );
#else
  layout->addWidget( buttonBox, 2, 0, AlignRight );
#endif
  
#if( FLUX_USE_COPY_TO_CLIPBOARD )
  LOAD_JAVASCRIPT(wApp, "FluxTool.cpp", "FluxTool", wtjsCopyFluxDataTextToClipboard );
  
  m_copyBtn = new WPushButton( "Copy To Clipboard" );
  m_copyBtn->clicked().connect( "function(s,e){ "
    "var success = Wt.WT.CopyFluxDataTextToClipboard(s,e,'" + m_copyBtn->id() + "'); "
    "Wt.emit( '" + id() + "', {name:'infocopied', eventObject:e}, success );"
  "}" );
  layout->addWidget( m_copyBtn, 2, 0, AlignLeft | AlignMiddle );
  
  m_infoCopied.connect( boost::bind( &FluxToolWidget::tableCopiedToCliboardCallback, this, _1 ) );
#endif
  
  m_displayInfoLevel = DisplayInfoLevel::Simple; //set to simple so #setDisplayInfoLevel(...) will do its work
  setDisplayInfoLevel( DisplayInfoLevel::Normal );
}//void init()


void FluxToolWidget::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  if( m_needsTableRefresh )
    refreshPeakTable();
  
  WContainerWidget::render( flags );
}//void render( Wt::WFlags<Wt::RenderFlag> flags );


void FluxToolWidget::setTableNeedsUpdating()
{
  m_needsTableRefresh = true;
  scheduleRender();
}//void setTableNeedsUpdating()


void FluxToolWidget::refreshPeakTable()
{
  PeakModel *peakmodel = m_interspec->peakModel();
  
  m_nucNames.clear();
  m_data.clear();
  m_uncertainties.clear();
  
  m_msg->setText( "" );
  
  double distance = 1.0*PhysicalUnits::meter;
  try
  {
    distance = PhysicalUnits::stringToDistance( m_distance->text().toUTF8() );
  }catch(...)
  {
    m_msg->setText( "Invalid Distance" );
    m_tableUpdated.emit();
    return;
  }
  
  auto det = m_detector->detector();
  if( !det || !det->isValid() )
  {
    m_msg->setText( "No Detector Response Function Chosen" );
    m_tableUpdated.emit();
    return;
  }
  
  auto spec = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  if( !spec )
  {
    m_msg->setText( "No foreground spectrum loaded" );
    m_tableUpdated.emit();
    return;
  }
  
  const float live_time = spec->gamma_live_time();
  if( live_time <= 0.0f )
  {
    m_msg->setText( "Invalid foregorund livetime" );
    m_tableUpdated.emit();
    return;
  }
  
  const vector<PeakDef> peaks = peakmodel->peakVec();
  
  const size_t npeaks = peaks.size();
  m_nucNames.resize( npeaks );
  m_data.resize( npeaks );
  m_uncertainties.resize( npeaks );

  for( int i = 0; i < static_cast<int>(npeaks); ++i )
  {
    m_data[i].fill( 0.0 );
    m_uncertainties[i].fill( 0.0 );
    
    const PeakDef &peak = peaks[i];
    
    const double energy = peak.mean();
    
    const double amp = peak.peakArea();  //ToDO: make sure this works for non-Gaussian peaks
    const double ampUncert = peak.peakAreaUncert();
    const double cps = amp / live_time;
    const double cpsUncert = ampUncert / live_time;
    const double intrisic = det->intrinsicEfficiency(energy);
    const double geomEff = det->fractionalSolidAngle( det->detectorDiameter(), distance );
    const double totaleff = det->efficiency(energy, distance );
    
    
    if( peak.parentNuclide() )
      m_nucNames[i] = peak.parentNuclide()->symbol;
    else if( peak.xrayElement() )
      m_nucNames[i] = peak.xrayElement()->symbol;
    else if( peak.reaction() )
      m_nucNames[i] = peak.reaction()->name();
    
    m_data[i][FluxColumns::FluxEnergyCol] = energy;
    m_data[i][FluxColumns::FluxPeakCpsCol] = cps;
    m_uncertainties[i][FluxColumns::FluxPeakCpsCol] = cpsUncert;
    
    m_data[i][FluxColumns::FluxGeometricEffCol] = geomEff;
    m_data[i][FluxColumns::FluxIntrinsicEffCol] = intrisic;
    //ToDo: Check if there is an uncertainty on DRF, and if so include that.
    
    if( totaleff <= 0.0 || intrisic <= 0.0 )
    {
      m_data[i][FluxColumns::FluxFluxOnDetCol]      = std::numeric_limits<double>::infinity();
      m_data[i][FluxColumns::FluxFluxPerCm2PerSCol] = std::numeric_limits<double>::infinity();
      m_data[i][FluxColumns::FluxGammasInto4PiCol]  = std::numeric_limits<double>::infinity();
    }else
    {
      const double fluxOnDet = cps / intrisic;
      const double fluxOnDetUncert = cpsUncert / intrisic;
      
      //gammas into 4pi
      const double gammaInto4pi = cps / totaleff;
      const double gammaInto4piUncert = cpsUncert / totaleff;
      
      //Flux in g/cm2/s
      const double distance_cm = distance / PhysicalUnits::cm;
      const double flux = gammaInto4pi / (4*M_PI*distance_cm*distance_cm);
      const double fluxUncert = gammaInto4piUncert / (4*M_PI*distance_cm*distance_cm);
      
      
      m_data[i][FluxColumns::FluxFluxOnDetCol] = fluxOnDet;
      m_uncertainties[i][FluxColumns::FluxFluxOnDetCol] = fluxOnDetUncert;
      
      m_data[i][FluxColumns::FluxFluxPerCm2PerSCol] = flux;
      m_uncertainties[i][FluxColumns::FluxFluxPerCm2PerSCol] = fluxUncert;
      
      m_data[i][FluxColumns::FluxGammasInto4PiCol] = gammaInto4pi;
      m_uncertainties[i][FluxColumns::FluxGammasInto4PiCol] = gammaInto4piUncert;
    }//if( eff > 0 ) / else
  }//for( const PeakDef &peak : peaks )
  
#if( FLUX_USE_COPY_TO_CLIPBOARD )
  stringstream pastebrdtxt;
  FluxToolImp::FluxCsvResource::data_to_strm( this, pastebrdtxt, true, m_displayInfoLevel );
  m_copyBtn->doJavaScript( "$('#" + m_copyBtn->id() + "').data('TableData','" + pastebrdtxt.str() + "');" );
#endif
  
  m_tableUpdated.emit();
}//void refreshPeakTable()




void FluxToolWidget::setDisplayInfoLevel( const DisplayInfoLevel disptype )
{
  if( disptype == m_displayInfoLevel )
    return;
  
  m_displayInfoLevel = disptype;
  
  setTableNeedsUpdating();
  
  for( FluxColumns col = FluxColumns(0); col < FluxColumns::FluxNumColumns; col = FluxColumns(col + 1) )
  {
    const bool show = FluxToolImp::showFluxColumn(m_displayInfoLevel, col);
    m_table->setColumnHidden( col, !show );
    
    WLength length;
    switch( col )
    {
      case FluxEnergyCol:         length = WLength(7.5, WLength::FontEm); break;
      case FluxNuclideCol:        length = WLength(5.0, WLength::FontEm); break;
      case FluxPeakCpsCol:        length = WLength(7.5, WLength::FontEm); break;
      case FluxIntrinsicEffCol:   length = WLength(6.5, WLength::FontEm); break;
      case FluxGeometricEffCol:   length = WLength(6.5, WLength::FontEm); break;
      case FluxFluxOnDetCol:      length = WLength(7.5, WLength::FontEm); break;
      case FluxFluxPerCm2PerSCol: length = WLength(9.0, WLength::FontEm); break;
      case FluxGammasInto4PiCol:  length = WLength(9.0, WLength::FontEm); break;
      case FluxNumColumns:        break;
    }//switch( col )
      
    m_table->setColumnWidth( col, length);
  }//for( loop over columns )
  
  m_table->refreshColWidthLayout();
}//void setDisplayInfoLevel( const bool minonly )


#if( FLUX_USE_COPY_TO_CLIPBOARD )
void FluxToolWidget::tableCopiedToCliboardCallback( const int copied )
{
  switch( copied )
  {
    case 0:
      passMessage( "Failed to copy to clipboard - maybe a permissions issue - sorry.", "", 3 );
      break;
    
    case 1:
      passMessage( "Copied table to clipboard.", "", 0 );
      break;
    
    case 2:
      passMessage( "Copied table to clipboard as HTML text.", "", 0 );
      break;
      
    default:
      passMessage( "Unknown result of copy command - sorry.", "", 3 );
      break;
  }//switch( copied )
}//void tableCopiedToCliboardCallback( const int copied )
#endif

