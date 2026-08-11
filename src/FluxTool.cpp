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
#include <Wt/WApplication>
#include <Wt/WEnvironment>
#include <Wt/WRadioButton>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WItemDelegate>
#include <Wt/WRegExpValidator>
#include <Wt/WContainerWidget>
#include <Wt/WAbstractItemModel>
#include <Wt/WAbstractItemDelegate>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/FluxTool.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/RowStretchTreeView.h"

#if( USE_QR_CODES )
#include <Wt/Utils>

#include "InterSpec/QrCode.h"
#endif

using namespace std;
using namespace Wt;


#if( FLUX_USE_COPY_TO_CLIPBOARD )

// See also CopyUrlToClipboard in QrCode.cpp
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
  //  ToDo: We could probably try the clipboard API to copy formatted text
  //        See https://developer.mozilla.org/en-US/docs/Web/API/Clipboard
  
  if( window.clipboardData && window.clipboardData.setData ) {
    return window.clipboardData.setData("Text", text);  // IE
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
  
  
  bool showFluxColumn( FluxToolWidget::DisplayInfoLevel disptype,
                      const FluxToolWidget::FluxColumns col,
                      const shared_ptr<const DetectorPeakResponse> &det )
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
          
          case FluxToolWidget::FluxColumns::FluxIntrinsicEffCol:
          case FluxToolWidget::FluxColumns::FluxGeometricEffCol:
          case FluxToolWidget::FluxColumns::FluxPeakCpsCol:
          case FluxToolWidget::FluxColumns::FluxFluxPerCm2PerSCol:
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
            
          case FluxToolWidget::FluxColumns::FluxIntrinsicEffCol:
          case FluxToolWidget::FluxColumns::FluxGeometricEffCol:
          case FluxToolWidget::FluxColumns::FluxNuclideCol:
          case FluxToolWidget::FluxColumns::FluxFluxOnDetCol:
          case FluxToolWidget::FluxColumns::FluxNumColumns:
            return false;
            break;
        }//switch( col )
        
        break;
        
      case FluxToolWidget::DisplayInfoLevel::Extended:
        switch( col )
        {
          case FluxToolWidget::FluxColumns::FluxEnergyCol:
          case FluxToolWidget::FluxColumns::FluxPeakCpsCol:
          case FluxToolWidget::FluxColumns::FluxFluxPerCm2PerSCol:
          case FluxToolWidget::FluxColumns::FluxGammasInto4PiCol:
          case FluxToolWidget::FluxColumns::FluxNuclideCol:
          case FluxToolWidget::FluxColumns::FluxFluxOnDetCol:
            return true;
            break;
            
          case FluxToolWidget::FluxColumns::FluxIntrinsicEffCol:
          case FluxToolWidget::FluxColumns::FluxGeometricEffCol:
            return (!det || !det->isFixedGeometry());
            break;
            
          case FluxToolWidget::FluxColumns::FluxNumColumns:
            return false;
            break;
        }//switch( col )
        
        return true;
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


  /** Returns, for each entry of `peaks`, whether any of `energies` names that peak.

   An energy names the peak whose mean is nearest to it, and then only if that mean is
   within 0.5 FWHM.  The tolerance necessarily comes from the candidate peak, since only
   the mean is remembered (both in the app-URL, and in FluxToolWidget::m_deselectedEnergies).

   Marking is idempotent, so two energies landing on one peak, or an energy sitting between
   the two halves of a doublet, give an answer that doesnt depend on iteration order, and
   hence doesnt flap between table refreshes.

   This is the only place peak/energy matching is done; #FluxToolWidget::refreshPeakTable,
   #FluxToolWidget::handleAppUrl, and #FluxToolWidget::encodeStateToUrl all come through here.
   */
  vector<bool> peaks_matching_energies( const vector<PeakDef> &peaks,
                                        const vector<double> &energies )
  {
    vector<bool> matched( peaks.size(), false );

    for( const double energy : energies )
    {
      size_t nearest = peaks.size();
      double nearest_dist = std::numeric_limits<double>::max();

      for( size_t i = 0; i < peaks.size(); ++i )
      {
        const double dist = fabs( peaks[i].mean() - energy );
        if( dist < nearest_dist )
        {
          nearest_dist = dist;
          nearest = i;
        }
      }//for( loop over peaks, finding nearest )

      if( nearest >= peaks.size() )
        continue;

      // PeakDef::fwhm() throws for data-defined peaks, so use the same fallback as the
      //  rest of the code-base (e.g., InterSpec::searchForSinglePeak(...)).
      const PeakDef &peak = peaks[nearest];
      const double fwhm = peak.gausPeak() ? peak.fwhm() : 0.5*peak.roiWidth();

      if( nearest_dist <= 0.5*fwhm )
        matched[nearest] = true;
    }//for( const double energy : energies )

    return matched;
  }//vector<bool> peaks_matching_energies(...)


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
      if( index.column() == FluxToolWidget::sm_selectColumn )
        return Wt::WFlags<Wt::ItemFlag>(Wt::ItemFlag::ItemIsUserCheckable);

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
      
      reset(); //JIC, to make sure the GUI refreshes
    }//void handleFluxToolWidgetUpdated()
    
    
    virtual boost::any data( const Wt::WModelIndex &index,
                             int role = Wt::DisplayRole ) const
    {
      if( index.parent().isValid() )
        return boost::any();
      
      const int row = index.row();
      const int col = index.column();

      if( row < 0 || col < 0
         || col > FluxToolWidget::sm_selectColumn
         || row >= static_cast<int>(m_sort_indices.size()) )
        return boost::any();

      const size_t realind = m_sort_indices[row];

      const auto &nucs = m_fluxtool->m_nucNames;
      const auto &data = m_fluxtool->m_data;
      const auto &uncerts = m_fluxtool->m_uncertainties;

      assert( nucs.size() == data.size() );
      assert( uncerts.size() == data.size() );
      assert( realind < data.size() );

      // The selection column has no value/uncertainty in m_data; it is checkbox-only.
      //  Note that Wt::WItemDelegate decides whether to render a checkbox at all based on
      //  CheckStateRole being non-empty here - the ItemIsUserCheckable flag only controls
      //  if that checkbox is enabled - so the other columns must keep returning an empty
      //  boost::any for CheckStateRole (via the `default:` case, below).
      if( col == FluxToolWidget::sm_selectColumn )
      {
        if( role != Wt::CheckStateRole )
          return boost::any();

        assert( m_fluxtool->m_rowSelected.size() == data.size() );
        if( realind >= m_fluxtool->m_rowSelected.size() )
          return boost::any();

        return boost::any( static_cast<bool>(m_fluxtool->m_rowSelected[realind]) );
      }//if( col == FluxToolWidget::sm_selectColumn )

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


    /** Handles the user clicking a rows checkbox.

     Wt::WItemDelegate::onCheckedChange sends a bool (rather than a Wt::CheckState), since
     we dont set ItemIsTristate.
     */
    virtual bool setData( const Wt::WModelIndex &index, const boost::any &value,
                          int role = Wt::EditRole )
    {
      if( index.parent().isValid()
         || (index.column() != FluxToolWidget::sm_selectColumn)
         || (role != Wt::CheckStateRole)
         || (index.row() < 0)
         || (index.row() >= static_cast<int>(m_sort_indices.size()))
         || (value.type() != typeid(bool)) )
      {
        return false;
      }

      // Sorting means the displayed row is not the m_data row; m_rowSelected parallels
      //  m_data, so map through m_sort_indices, just like data(...) does.
      const size_t realind = m_sort_indices[index.row()];
      assert( realind < m_fluxtool->m_rowSelected.size() );
      if( realind >= m_fluxtool->m_rowSelected.size() )
        return false;

      m_fluxtool->handleRowSelectionChanged( realind, boost::any_cast<bool>(value) );

      dataChanged().emit( index, index );

      return true;
    }//bool setData(...)


    /** Tells the view every rows checkbox may have changed.

     Used when selection changes from something other than a click on that checkbox (an
     undo/redo, the "select rows" checkbox, or an app-URL), so the view redraws without the
     full model reset() that #handleFluxToolWidgetUpdated does.
     */
    void notifySelectionChanged()
    {
      if( m_sort_indices.empty() )
        return;

      const int lastrow = static_cast<int>(m_sort_indices.size()) - 1;
      dataChanged().emit( index(0, FluxToolWidget::sm_selectColumn),
                          index(lastrow, FluxToolWidget::sm_selectColumn) );
    }//void notifySelectionChanged()
    
    
    boost::any headerData( int section,
                          Wt::Orientation orientation = Wt::Horizontal,
                          int role = Wt::DisplayRole ) const
    {
      if( section < 0 || section > FluxToolWidget::sm_selectColumn
         || orientation != Wt::Horizontal
         || role != Wt::DisplayRole )
        return boost::any();

      if( section == FluxToolWidget::sm_selectColumn )
        return boost::any( m_fluxtool->m_selectColName );

      return boost::any( m_fluxtool->m_colnames[section] );
    }

    virtual int columnCount( const Wt::WModelIndex &parent
                            = Wt::WModelIndex() ) const
    {
      if( parent.isValid() )
        return 0;

      // Always report the selection column, even when rows arent being selected - the view
      //  hides it instead (see FluxToolWidget::setDisplayInfoLevel).  Wt::WAbstractItemView
      //  only truncates/appends its column list on a model reset(), and never re-initializes
      //  the entries that survive, so a column count that changed would have its ColumnInfo
      //  destroyed by the reset() in handleFluxToolWidgetUpdated(), then silently recreated
      //  with default width/sorting/delegate the next time rows became selectable.
      return FluxToolWidget::sm_selectColumn + 1;
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
      // index_compare_sort indexes m_data, which only has the FluxColumns; sorting on
      //  anything else (i.e., the selection column) would read off the end of it.
      //  RowStretchTreeView is told not to offer sorting on that column, but that is set
      //  per-column and defaults to on, so dont rely on it alone.
      assert( (column >= 0) && (column < FluxToolWidget::FluxColumns::FluxNumColumns) );
      if( (column < 0) || (column >= FluxToolWidget::FluxColumns::FluxNumColumns) )
        return;

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
      
      const shared_ptr<const DetectorPeakResponse> det = tool->m_detector->detector();
      
      for( FluxToolWidget::FluxColumns col = FluxToolWidget::FluxColumns(0);
          col < FluxToolWidget::FluxColumns::FluxNumColumns; col = FluxToolWidget::FluxColumns(col+1) )
      {
        const WString &colname = tool->m_colnamesCsv[col];
        
        if( !FluxToolImp::showFluxColumn(disptype, col, det) )
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
        // When the user is picking rows, only the checked ones get output.  m_rowSelected
        //  parallels m_data, both in un-sorted order, so it indexes straight off `row`.
        //  The size check is just belt-and-braces (refreshPeakTable() keeps the two the
        //  same size, including on its early-returns, where both end up empty).
        assert( tool->m_rowSelected.size() == tool->m_data.size() );
        if( tool->m_selectRows
           && (row < tool->m_rowSelected.size())
           && !tool->m_rowSelected[row] )
        {
          continue;
        }

        if( html )
          strm << "\\t<tr>";

        for( FluxToolWidget::FluxColumns col = FluxToolWidget::FluxColumns(0);
            col < FluxToolWidget::FluxColumns::FluxNumColumns; col = FluxToolWidget::FluxColumns(col+1) )
        {
          if( !FluxToolImp::showFluxColumn(disptype, col, det) )
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
: AuxWindow( WString::tr("window-title-flux-tool"),
  (AuxWindowProperties::TabletNotFullScreen
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
  
  WPushButton *closeButton = nullptr;
  if( viewer && viewer->isPhone() )
    closeButton = addCloseButtonToFooter();
  
#if( USE_QR_CODES )
  WPushButton *qr_btn = new WPushButton( footer() );
  qr_btn->setText( WString::tr("QR Code") );
  qr_btn->setIcon( "InterSpec_resources/images/qr-code.svg" );
  qr_btn->setStyleClass( "LinkBtn DownloadBtn DialogFooterQrBtn" );
  qr_btn->clicked().preventPropagation();
  qr_btn->clicked().connect( std::bind( [this](){
    try
    {
      const string url = "interspec://flux/?" + Wt::Utils::urlEncode(m_fluxTool->encodeStateToUrl());
      QrCode::displayTxtAsQrCode( url, WString::tr("ftw-qr-tool-state-title"),
                                 WString::tr("ftw-qr-tool-state-text") );
    }catch( std::exception &e )
    {
      passMessage( WString::tr("app-qr-err").arg(e.what()), WarningWidget::WarningMsgHigh );
    }
  }) );
#endif //USE_QR_CODES
  
  WContainerWidget *buttonDiv = footer();
  
  if( !viewer || !viewer->isPhone() )
    closeButton = addCloseButtonToFooter();
  
  assert( closeButton );
  closeButton->clicked().connect( this, &AuxWindow::hide );
  
  WResource *csv = new FluxToolImp::FluxCsvResource( m_fluxTool );
#if( BUILD_AS_OSX_APP || IOS )
  WAnchor *csvButton = new WAnchor( WLink(csv), buttonDiv );
  csvButton->setTarget( AnchorTarget::TargetNewWindow );
  csvButton->setStyleClass( "LinkBtn DownloadLink" );
#else
  WPushButton *csvButton = new WPushButton( buttonDiv );
  csvButton->setIcon( "InterSpec_resources/images/download_small.svg" );
  csvButton->setLink( WLink(csv) );
  csvButton->setLinkTarget( Wt::TargetNewWindow );
  csvButton->setStyleClass( "LinkBtn DownloadBtn" );
  
#endif
  
  csvButton->setText( WString::tr("CSV") );
  //csvButton->setAttributeValue( "style", "float: none;" );  //Keep the CSV download to the left side of the close button.  .LinkBtn style class has the button float right..
  csvButton->setFloatSide( Wt::Side::Left );
  
  auto enableDisableCsv = [csvButton,this](){
    csvButton->setDisabled( m_fluxTool->m_data.empty() );
  };
  
  m_fluxTool->m_tableUpdated.connect( std::bind(enableDisableCsv) );
  csvButton->disable();
  
  show();
  
  const int screenW = viewer->renderedWidth();
  const int screenH = viewer->renderedHeight();
  const int width = ((screenW < 680) ? screenW : 680);
  const int height = ((screenH < 420) ? screenH : 420);
  resizeWindow( width, height );
  
  resizeToFitOnScreen();
  centerWindowHeavyHanded();
}//FluxToolWindow(...) constructor


FluxToolWindow::~FluxToolWindow()
{
}


void FluxToolWindow::handleAppUrl( const std::string &query_str )
{
  m_fluxTool->handleAppUrl( query_str );
}//void handleAppUrl( std::string query_str )


std::string FluxToolWindow::encodeStateToUrl() const
{
  return m_fluxTool->encodeStateToUrl();
}


FluxToolWidget::FluxToolWidget( InterSpec *viewer, Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_interspec( viewer ),
    m_detector( nullptr ),
    m_narrowLayout( false ),
    m_msg( nullptr ),
    m_distance( nullptr ),
    m_prevDistance(),
    m_table( nullptr ),
    m_fluxModel( nullptr ),
    m_selectRowsCb( nullptr ),
    m_selectRows( false ),
#if( FLUX_USE_COPY_TO_CLIPBOARD )
    m_copyBtn( nullptr ),
    m_infoCopied( this, "infocopied", true ),
#endif
    m_needsTableRefresh( true ),
    m_displayInfoLevel( DisplayInfoLevel::Normal ),
    m_displayLevelButtons( nullptr ),
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


void FluxToolWidget::handleAppUrl( std::string query_str )
{
  // Do we want to add an undo/redo step here?
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  
  const map<string,string> parts = AppUtils::query_str_key_values( query_str );
  
  const auto ver_iter = parts.find( "VER" );
  if( ver_iter == end(parts) )
    Wt::log("warn") << "No 'VER' field in Flux Tool URI.";
  else if( ver_iter->second != "1" && !SpecUtils::starts_with(ver_iter->second, "1.") )
    throw runtime_error( "Can not read Flux Tool URI version '" + ver_iter->second + "'" );
  
  const auto dist_iter = parts.find( "DIST" );
  if( dist_iter == end(parts) )
    throw runtime_error( "URL is invalid Flux Tool URI = no 'DIST' component" );
  
  string dist = dist_iter->second;
  try
  {
    PhysicalUnits::stringToDistance( dist );
  }catch( std::exception & )
  {
    Wt::log("warn") << "FluxToolWidget URI had invalid distance: '" << dist << "'";
    throw std::runtime_error( "Invalid Flux Tool distance specified in URI" );
  }
  
  
  DisplayInfoLevel level = DisplayInfoLevel::Normal;
  string disp_str = "NORMAL";
  const auto disp_iter = parts.find( "DISPLAY" );
  if( disp_iter != end(parts) )
    disp_str = disp_iter->second;
  if( SpecUtils::iequals_ascii(disp_str, "SIMPLE") )
    level = DisplayInfoLevel::Simple;
  else if( SpecUtils::iequals_ascii(disp_str, "NORMAL") )
    level = DisplayInfoLevel::Normal;
  else if( SpecUtils::iequals_ascii(disp_str, "EXTENDED") )
    level = DisplayInfoLevel::Extended;
  else
    Wt::log("warn") << "FluxToolWidget URI had invalid display type: '" << disp_str << "'";
  
  // "SELROWS" and "SEL" are both optional (added in URI version 1.1), so URIs saved by
  //  previous versions - and QR codes already out in the wild - keep working.
  bool selectRows = false;
  const auto selrows_iter = parts.find( "SELROWS" );
  if( selrows_iter != end(parts) )
    selectRows = (selrows_iter->second == "1") || SpecUtils::iequals_ascii(selrows_iter->second, "true");

  // The URI names the *selected* peaks, but we track the un-selected ones (so a peak fit
  //  later defaults to being included), so invert against the peaks currently loaded.  An
  //  absent "SEL" means every peak is selected.
  vector<double> deselected;
  const auto sel_iter = parts.find( "SEL" );
  string sel_str = (sel_iter != end(parts)) ? sel_iter->second : string();
  SpecUtils::trim( sel_str );

  // Keyed off the field being *present*, not off it being non-empty: "SEL=" is how a URI
  //  says no rows are selected (which is what encodeStateToUrl() writes when the user has
  //  un-checked every row).  Treating an empty value as "all selected" would invert the
  //  state on the way back in.
  if( sel_iter != end(parts) )
  {
    // SpecUtils::split drops empty tokens, so "SEL=" gives no energies, hence no peaks
    //  matched, hence every peak de-selected.
    vector<string> fields;
    SpecUtils::split( fields, sel_str, "," );

    vector<double> selected;
    for( string field : fields )
    {
      SpecUtils::trim( field );

      double energy = 0.0;
      if( !SpecUtils::parse_double( field.c_str(), field.size(), energy ) )
        throw runtime_error( "Invalid Flux Tool selected-row energy in URI: '" + field + "'" );

      selected.push_back( energy );
    }//for( string field : fields )

    PeakModel *peakmodel = m_interspec->peakModel();
    const vector<PeakDef> peaks = peakmodel ? peakmodel->peakVec() : vector<PeakDef>{};
    const vector<bool> is_selected = FluxToolImp::peaks_matching_energies( peaks, selected );

    for( size_t i = 0; i < peaks.size(); ++i )
    {
      if( !is_selected[i] )
        deselected.push_back( peaks[i].mean() );
    }
  }//if( sel_iter != end(parts) )

  setDisplayInfoLevel( level, true );
  applySelectionState( selectRows, deselected );
  m_prevDistance = WString::fromUTF8(dist);
  m_distance->setText( m_prevDistance );
  setTableNeedsUpdating();
}//void handleAppUrl( std::string query_str )


std::string FluxToolWidget::encodeStateToUrl() const
{
  // "dist=1.2m&display=low&ver=1"
  
  string dist = m_distance->text().toUTF8();
  try
  {
    PhysicalUnits::stringToDistance(dist);
  }catch( std::exception & )
  {
    dist = "";
  }
  
  string answer = "VER=1.1&DIST=" + dist + "&DISPLAY=";

  switch( displayInfoLevel() )
  {
    case DisplayInfoLevel::Simple:   answer += "SIMPLE";   break;
    case DisplayInfoLevel::Normal:   answer += "NORMAL";   break;
    case DisplayInfoLevel::Extended: answer += "EXTENDED"; break;
  }//switch( displayInfoLevel() )

  if( m_selectRows )
    answer += "&SELROWS=1";

  // List the selected peaks, but only when that isnt all of them - it keeps the URI (and
  //  so the QR code) short in the common case.  Note this is emitted even when
  //  !m_selectRows, since the URI is the only thing that carries state across the window
  //  being closed, and the user gets their rows back if they re-tick "select rows".
  //
  //  Deliberately worked out from the peaks, rather than from m_data, so this is right even
  //  when the table is empty (no DRF chosen, invalid distance, ...).
  if( !m_deselectedEnergies.empty() )
  {
    PeakModel *peakmodel = m_interspec ? m_interspec->peakModel() : nullptr;
    const vector<PeakDef> peaks = peakmodel ? peakmodel->peakVec() : vector<PeakDef>{};
    const vector<bool> deselected = FluxToolImp::peaks_matching_energies( peaks, m_deselectedEnergies );

    const bool anyDeselected = std::any_of( begin(deselected), end(deselected),
                                            []( const bool v ){ return v; } );
    if( anyDeselected )
    {
      answer += "&SEL=";
      char buffer[32] = { '\0' };
      bool first = true;
      for( size_t i = 0; i < peaks.size(); ++i )
      {
        if( deselected[i] )
          continue;

        // "%.2f" to match how the table prints energies; it also keeps the '+' that "%g"
        //  would put in an exponent out of the URI (Wt::Utils::urlEncode escapes '+', as
        //  it does the ',' separator).
        snprintf( buffer, sizeof(buffer), "%.2f", peaks[i].mean() );
        answer += (first ? "" : ",") + string(buffer);
        first = false;
      }//for( loop over peaks )
    }//if( anyDeselected )
  }//if( !m_deselectedEnergies.empty() )

  //SpecUtils::erase_any_character( answer, " \t\n\r" );

  return answer;
}//std::string encodeStateToUrl() const


void FluxToolWidget::init()
{
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  
  assert( m_interspec );
  assert( !m_detector );
  
  if( m_interspec )
    m_interspec->useMessageResourceBundle( "FluxTool" );
    
  if( m_interspec && m_interspec->isPhone() )
  {
    int w = m_interspec->renderedWidth();
    int h = m_interspec->renderedHeight();
    if( w < 100 )
    {
      w = wApp->environment().screenWidth();
      h = wApp->environment().screenHeight();
    }
    m_narrowLayout = ((w > 100) && (w < 500));
  }//if( phone )
  
  for( FluxColumns col = FluxColumns(0); col < FluxColumns::FluxNumColumns; col = FluxColumns(col + 1) )
  {
    switch( col )
    {
      case FluxEnergyCol:
        m_colnames[col] = WString::tr( m_narrowLayout ? "Energy" : "Energy (keV)" );
        m_colnamesCsv[col] = WString::fromUTF8("Energy (keV)");
        break;
      case FluxNuclideCol:
        m_colnames[col] = WString::tr( m_narrowLayout ? "Nuc." : "Nuclide");
        m_colnamesCsv[col] = WString::fromUTF8("Nuclide");
        break;
      case FluxPeakCpsCol:
        m_colnames[col] = WString("{1} {2}").arg(WString::tr("Peak")).arg(WString::tr("CPS"));
        m_colnamesCsv[col] = WString::fromUTF8("Peak CPS");
        break;
      case FluxIntrinsicEffCol:
        m_colnames[col] = WString::tr("ftw-hdr-intrinsic-eff");
        m_colnamesCsv[col] = WString::fromUTF8("Intrinsic Efficiency");
        break;
      case FluxGeometricEffCol:
        m_colnames[col] = WString::tr("ftw-hdr-geom-eff");
        m_colnamesCsv[col] = WString::fromUTF8("Geometric Efficiency");
        break;
      case FluxFluxOnDetCol:
        m_colnames[col] = WString::tr("ftw-hdr-flux-on-det");
        m_colnamesCsv[col] = WString::fromUTF8("Flux on Detector (gammas/s)");
        break;
      case FluxFluxPerCm2PerSCol:
        m_colnames[col] = WString::tr("ftw-hdr-flux-cm");
        m_colnamesCsv[col] = WString::fromUTF8("Flux (gammas/cm2/s)");
        break;
      case FluxGammasInto4PiCol:
        m_colnames[col] = WString::tr("ftw-hdr-flux-4pi");
        m_colnamesCsv[col] = WString::fromUTF8("gammas/4pi/s");
        break;
      case FluxNumColumns:        break;
    }//switch( col )
  }//for( loop over columns )
  
  
  wApp->useStyleSheet( "InterSpec_resources/FluxTool.css" );
  
  const bool showToolTips = m_interspec ? UserPreferences::preferenceValue<bool>( "ShowTooltips", m_interspec ) : false;
  
  addStyleClass( "FluxToolWidget" );
  
  WGridLayout *layout = new WGridLayout();
  layout->setContentsMargins( 0, 0, 0, 0 );
  setLayout( layout );
  
  WTable *distDetRow = new WTable();
  distDetRow->addStyleClass( "FluxDistMsgDetTbl" );
  layout->addWidget( distDetRow, 0, 0, 1, 2 );

  
  auto distCell = distDetRow->elementAt( m_narrowLayout ? 1 : 0, 0 );
  distCell->addStyleClass( "FluxDistCell" );
  WLabel *label = new WLabel( WString::tr("distance-label"), distCell );
  label->addStyleClass( "FluxDistLabel" );
  
  m_prevDistance = "100 cm";
  m_distance = new WLineEdit( m_prevDistance, distCell );
  m_distance->addStyleClass( "FluxDistanceEnter" );
  label->setBuddy(m_distance);
  
  m_distance->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_distance->setAttributeValue( "autocorrect", "off" );
  m_distance->setAttributeValue( "spellcheck", "off" );
#endif
  
  WRegExpValidator *validator = new WRegExpValidator( PhysicalUnits::sm_distanceUnitOptionalRegex, this );
  validator->setFlags( Wt::MatchCaseInsensitive );
  m_distance->setValidator( validator );
  HelpSystem::attachToolTipOn( m_distance, WString::tr("ftw-tt-distance"), showToolTips );
  m_distance->changed().connect( this, &FluxToolWidget::distanceUpdated );
  m_distance->enterPressed().connect( this, &FluxToolWidget::distanceUpdated );
  
  SpectraFileModel *specFileModel = m_interspec->fileManager()->model();
  m_detector = new DetectorDisplay( m_interspec, specFileModel );
  m_detector->addStyleClass( "FluxDet" );
  m_interspec->detectorChanged().connect( boost::bind( &FluxToolWidget::handleDrfChange, this, boost::placeholders::_1 ) );
  m_interspec->detectorModified().connect( boost::bind( &FluxToolWidget::handleDrfChange, this, boost::placeholders::_1 ) );
  
  if( m_narrowLayout )
  {
    auto detCell = distDetRow->elementAt(0,0);
    detCell->setColumnSpan( 2 );
    detCell->addWidget( m_detector );
  }else
  {
    distDetRow->elementAt(0,1)->addWidget( m_detector );
  }
  
  
  PeakModel *peakmodel = m_interspec->peakModel();
  assert( peakmodel );
  peakmodel->dataChanged().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  peakmodel->rowsRemoved().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  peakmodel->rowsInserted().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  peakmodel->layoutChanged().connect( this, &FluxToolWidget::setTableNeedsUpdating );

  auto msgCell = distDetRow->elementAt( m_narrowLayout ? 2 : 1,0);
  msgCell->setColumnSpan( 2 );
  msgCell->addStyleClass( "FluxMsgCell" );
  m_msg = new WText( "", Wt::XHTMLText, msgCell );
  m_msg->addStyleClass( "FluxMsg" );
  
  
  m_selectColName = WString::tr("ftw-hdr-use");

  m_fluxModel = new FluxToolImp::FluxModel( this );
  m_table = new RowStretchTreeView();
  m_table->setModel( m_fluxModel );
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

  // The stock delegate knows how to render a checkbox from CheckStateRole; using it just
  //  for the selection column saves teaching FluxRenderDelegate about checkboxes.
  m_table->setItemDelegateForColumn( sm_selectColumn, new WItemDelegate(this) );

  for( FluxColumns col = FluxColumns(0); col < FluxColumns::FluxNumColumns; col = FluxColumns(col + 1) )
    m_table->setSortingEnabled( col, (col!=FluxColumns::FluxGeometricEffCol) );

  // Sorting is per-column, and defaults to enabled, so the selection column needs this
  //  explicitly - the loop above only covers the data columns.  See FluxModel::sort().
  m_table->setSortingEnabled( sm_selectColumn, false );

  layout->addWidget( m_table, 1, 0, 1, 2 );
  layout->setRowStretch( 1, 1 );

  WContainerWidget *buttonBox = new WContainerWidget();
  buttonBox->addStyleClass( "FluxInfoAmount" );
  
  WRadioButton *simpleInfo = new WRadioButton( WString::tr("ftw-simple"), buttonBox );
  WRadioButton *standardInfo = new WRadioButton( WString::tr("ftw-standard"), buttonBox );
  WRadioButton *moreInfo = new WRadioButton( WString::tr("ftw-more"), buttonBox );
  
  m_displayLevelButtons = new WButtonGroup( buttonBox );
  m_displayLevelButtons->addButton( simpleInfo, static_cast<int>(DisplayInfoLevel::Simple) );
  m_displayLevelButtons->addButton( standardInfo, static_cast<int>(DisplayInfoLevel::Normal) );
  m_displayLevelButtons->addButton( moreInfo, static_cast<int>(DisplayInfoLevel::Extended) );
  m_displayLevelButtons->setCheckedButton( m_narrowLayout ? simpleInfo : standardInfo );

  m_displayLevelButtons->checkedChanged().connect( std::bind( [this](){
    const auto level = static_cast<DisplayInfoLevel>( m_displayLevelButtons->checkedId() );
    
    assert( (level == DisplayInfoLevel::Simple)
           || (level == DisplayInfoLevel::Normal)
           || (level == DisplayInfoLevel::Extended) );
    
    setDisplayInfoLevel( level, false );
  }) );
  
  
  layout->addWidget( buttonBox, 2, 1, AlignRight | AlignMiddle );

  WContainerWidget *outputOptsDiv = new WContainerWidget();
  outputOptsDiv->addStyleClass( "FluxOutputOpts" );
  layout->addWidget( outputOptsDiv, 2, 0, AlignLeft | AlignMiddle );

#if( FLUX_USE_COPY_TO_CLIPBOARD )
  LOAD_JAVASCRIPT(wApp, "FluxTool.cpp", "FluxTool", wtjsCopyFluxDataTextToClipboard );

  m_copyBtn = new WPushButton( WString::tr( m_narrowLayout ? "ftw-copy-btn-narrow" : "ftw-copy-btn"),
                               outputOptsDiv );

  // TODO: "upgrade" to using the InterSpecApp 'miscSignal' directly in CopyFluxDataTextToClipboard, and get rid of m_infoCopied signal handler
  //"Wt.emit( $('.specviewer').attr('id'), {name:'miscSignal'}, 'showMsg-info-' );"

  m_copyBtn->clicked().connect( "function(s,e){ "
    "var success = Wt.WT.CopyFluxDataTextToClipboard(s,e,'" + m_copyBtn->id() + "'); "
    "Wt.emit( '" + id() + "', {name:'infocopied', eventObject:e}, success );"
  "}" );

  m_infoCopied.connect( boost::bind( &FluxToolWidget::tableCopiedToCliboardCallback, this,
                                    boost::placeholders::_1 ) );
#endif

  m_selectRowsCb = new WCheckBox( WString::tr( m_narrowLayout ? "ftw-select-rows-narrow"
                                                              : "ftw-select-rows"),
                                  outputOptsDiv );
  m_selectRowsCb->addStyleClass( "FluxSelectRowsCb" );
  HelpSystem::attachToolTipOn( m_selectRowsCb, WString::tr("ftw-tt-select-rows"), showToolTips );
  m_selectRowsCb->checked().connect( boost::bind( &FluxToolWidget::setSelectRows, this, true ) );
  m_selectRowsCb->unChecked().connect( boost::bind( &FluxToolWidget::setSelectRows, this, false ) );

  setDisplayInfoLevel( m_narrowLayout ? DisplayInfoLevel::Simple : DisplayInfoLevel::Normal, true );
}//void init()


void FluxToolWidget::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  if( m_needsTableRefresh )
    refreshPeakTable();
  
  WContainerWidget::render( flags );
}//void render( Wt::WFlags<Wt::RenderFlag> flags );


void FluxToolWidget::distanceUpdated()
{
  const WString dist = m_distance->text();
  
  if( dist != m_prevDistance )
  {
    const WString prev = m_prevDistance;
    
    UndoRedoManager *undoRedo = UndoRedoManager::instance();
    if( undoRedo && undoRedo->canAddUndoRedoNow() )
    {
      auto undo_redo = [dist, prev]( const bool is_undo ){
        InterSpec *viewer = InterSpec::instance();
        FluxToolWindow *fluxwin = viewer ? viewer->createFluxTool() : nullptr;
        FluxToolWidget *tool = fluxwin ? fluxwin->m_fluxTool : nullptr;
        if( tool )
        {
          tool->m_prevDistance = is_undo ? prev : dist;
          tool->m_distance->setText( tool->m_prevDistance );
          tool->m_needsTableRefresh = true;
          tool->scheduleRender();
        }//if( tool )
      };
      
      auto undo = [=](){ undo_redo(true); };
      auto redo = [=](){ undo_redo(false); };
      undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Change flux tool distance." );
    }//if( add undo/redo )
    
    m_prevDistance = dist;
  }//if( dist != m_prevDistance )
  
  m_needsTableRefresh = true;
  scheduleRender();
}

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

  // Note that m_deselectedEnergies is deliberately *not* touched here - it is keyed by
  //  energy exactly so it can outlive the table.  Several of the early-returns below leave
  //  m_data empty, and re-deriving the users selections from an empty table would throw
  //  them away.
  m_rowSelected.clear();

  m_msg->setText( "" );

  // Shared tail for the cases below where there is nothing to compute.  The table is left
  //  empty, so the clipboard text is refreshed to match, rather than leaving the button
  //  holding the last table that could be computed.
  auto showMessageForEmptyTable = [this]( const char * const msg_key ){
    m_msg->setText( WString::tr(msg_key) );
#if( FLUX_USE_COPY_TO_CLIPBOARD )
    updateCopyToClipboardText();
#endif
    m_tableUpdated.emit();
  };

  double distance = 1.0*PhysicalUnits::meter;
  try
  {
    distance = PhysicalUnits::stringToDistance( m_distance->text().toUTF8() );
  }catch(...)
  {
    showMessageForEmptyTable( "ftw-invalid-dist" );
    return;
  }

  auto det = m_detector->detector();
  if( !det || !det->isValid() )
  {
    showMessageForEmptyTable( "ftw-no-drf" );
    return;
  }

  auto spec = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  if( !spec )
  {
    showMessageForEmptyTable( "ftw-no-foreground" );
    return;
  }

  const float live_time = spec->gamma_live_time();
  if( live_time <= 0.0f )
  {
    showMessageForEmptyTable( "ftw-invalid-livetime" );
    return;
  }

  const vector<PeakDef> peaks = peakmodel->peakVec();
  
  const size_t npeaks = peaks.size();
  m_nucNames.resize( npeaks );
  m_data.resize( npeaks );
  m_uncertainties.resize( npeaks );
  
  const bool fixed_geom = det->isFixedGeometry();

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
    const double intrinsic = det->intrinsicEfficiency(energy);
    const double geomEff = fixed_geom ? 1.0 : det->fractionalSolidAngle( det->detectorDiameter(), distance + det->detectorSetback() );
    const double totaleff = fixed_geom ? intrinsic : det->efficiency(energy, distance );
    
    
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
    m_data[i][FluxColumns::FluxIntrinsicEffCol] = intrinsic;
    
    // TODO: Check if there is an uncertainty on DRF, and if so include that.
    
    if( totaleff <= 0.0 || intrinsic <= 0.0 )
    {
      m_data[i][FluxColumns::FluxFluxOnDetCol]      = std::numeric_limits<double>::infinity();
      m_data[i][FluxColumns::FluxFluxPerCm2PerSCol] = std::numeric_limits<double>::infinity();
      m_data[i][FluxColumns::FluxGammasInto4PiCol]  = std::numeric_limits<double>::infinity();
    }else
    {
      const double fluxOnDet = cps / intrinsic;
      const double fluxOnDetUncert = cpsUncert / intrinsic;
      
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

  syncRowSelectionFromEnergies( peaks );

#if( FLUX_USE_COPY_TO_CLIPBOARD )
  updateCopyToClipboardText();
#endif

  m_tableUpdated.emit();
}//void refreshPeakTable()


void FluxToolWidget::syncRowSelectionFromEnergies( const vector<PeakDef> &peaks )
{
  assert( peaks.size() == m_data.size() );

  m_rowSelected.resize( m_data.size(), true );

  if( m_deselectedEnergies.empty() )
  {
    // Nothing remembered, so every row is selected - which is also the state after the
    //  resize() above, but be explicit since m_rowSelected may have been left over.
    std::fill( begin(m_rowSelected), end(m_rowSelected), true );
    return;
  }

  // deselected is sized to `peaks`, m_rowSelected to m_data; callers pass the peaks m_data
  //  was built from, so they agree - but bound the loop by both, since reading past a
  //  vector<bool> wouldnt crash, it would just silently give the wrong rows.
  const vector<bool> deselected = FluxToolImp::peaks_matching_energies( peaks, m_deselectedEnergies );
  assert( deselected.size() == m_rowSelected.size() );

  const size_t nrows = std::min( deselected.size(), m_rowSelected.size() );
  for( size_t i = 0; i < nrows; ++i )
    m_rowSelected[i] = !deselected[i];
}//void syncRowSelectionFromEnergies( const vector<PeakDef> &peaks )


#if( FLUX_USE_COPY_TO_CLIPBOARD )
void FluxToolWidget::updateCopyToClipboardText()
{
  stringstream pastebrdtxt;
  FluxToolImp::FluxCsvResource::data_to_strm( this, pastebrdtxt, true, m_displayInfoLevel );
  m_copyBtn->doJavaScript( "$('#" + m_copyBtn->id() + "').data('TableData','" + pastebrdtxt.str() + "');" );
}//void updateCopyToClipboardText()
#endif


void FluxToolWidget::handleDrfChange( std::shared_ptr<DetectorPeakResponse> drf )
{
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  
  const auto level = static_cast<DisplayInfoLevel>( m_displayLevelButtons->checkedId() );
  
  assert( (level == DisplayInfoLevel::Simple)
         || (level == DisplayInfoLevel::Normal)
         || (level == DisplayInfoLevel::Extended) );
  
  // The m_interspec->detectorChanged() signal will actually call this function
  //  (i.e., FluxToolWidget::handleDrfChange()) before DetectorDisplay::setDetector, so
  //  the following calls would use the wrong DRF, unless we explicitly set it here.
  m_detector->setDetector( drf );
  
  setDisplayInfoLevel( level, true );
  
  setTableNeedsUpdating();
}//void handleDrfChange()


void FluxToolWidget::setDisplayInfoLevel( const DisplayInfoLevel disptype, const bool force )
{
  if( !force && (disptype == m_displayInfoLevel) )
    return;
  
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( (disptype != m_displayInfoLevel) && undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    const DisplayInfoLevel prev = m_displayInfoLevel;
    auto undo_redo = [prev, disptype]( const bool is_undo ){
      InterSpec *viewer = InterSpec::instance();
      FluxToolWindow *fluxwin = viewer ? viewer->createFluxTool() : nullptr;
      FluxToolWidget *tool = fluxwin ? fluxwin->m_fluxTool : nullptr;
      if( tool )
        tool->setDisplayInfoLevel( is_undo ? prev : disptype, false );
    };
    auto undo = [=](){ undo_redo(true); };
    auto redo = [=](){ undo_redo(false); };
    undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Change flux tool display level." );
  }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
  
  m_displayLevelButtons->setSelectedButtonIndex( static_cast<int>(disptype) );
  
  m_displayInfoLevel = disptype;
  
  setTableNeedsUpdating();
  
  const shared_ptr<const DetectorPeakResponse> det = m_detector->detector();
  
  for( FluxColumns col = FluxColumns(0); col < FluxColumns::FluxNumColumns; col = FluxColumns(col + 1) )
  {
    const bool show = FluxToolImp::showFluxColumn(m_displayInfoLevel, col, det);
    m_table->setColumnHidden( col, !show );
    
    WLength length;
    switch( col )
    {
      case FluxEnergyCol:         
        length = WLength( m_narrowLayout ? 4.25 : 7.5, WLength::FontEm);
        break;
        
      case FluxNuclideCol:        
        length = WLength( m_narrowLayout ? 3.25 : 5.0, WLength::FontEm);
        break;
        
      case FluxPeakCpsCol:        
        length = WLength(7.5, WLength::FontEm);
        break;
        
      case FluxIntrinsicEffCol:   
        length = WLength(6.5, WLength::FontEm);
        break;
        
      case FluxGeometricEffCol:   
        length = WLength(6.5, WLength::FontEm);
        break;
        
      case FluxFluxOnDetCol:      
        length = WLength(7.5, WLength::FontEm);
        break;
        
      case FluxFluxPerCm2PerSCol:
        length = WLength(9.0, WLength::FontEm);
        break;
        
      case FluxGammasInto4PiCol:  
        length = WLength( m_narrowLayout ? 7.0 : 9.0, WLength::FontEm);
        break;
        
      case FluxNumColumns:
        break;
    }//switch( col )
      
    m_table->setColumnWidth( col, length);
  }//for( loop over columns )

  m_table->setColumnHidden( sm_selectColumn, !m_selectRows );

  // Must set the width explicitly: WTreeView::setColumnHidden() calls the (virtual)
  //  setColumnWidth() with the columns current width, which for a column we have never
  //  sized is Wt's 150px default - and RowStretchTreeView takes that as the columns share
  //  of the tables width.  The loop above is what keeps the data columns from doing this.
  m_table->setColumnWidth( sm_selectColumn, WLength(2.5, WLength::FontEm) );

  m_table->refreshColWidthLayout();
}//void setDisplayInfoLevel( const bool minonly )


void FluxToolWidget::setSelectRows( const bool selectRows )
{
  if( selectRows == m_selectRows )
    return;

  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    const bool prev = m_selectRows;
    const vector<double> deselected = m_deselectedEnergies;
    auto undo_redo = [prev, selectRows, deselected]( const bool is_undo ){
      InterSpec *viewer = InterSpec::instance();
      FluxToolWindow *fluxwin = viewer ? viewer->createFluxTool() : nullptr;
      FluxToolWidget *tool = fluxwin ? fluxwin->m_fluxTool : nullptr;
      if( tool )
        tool->applySelectionState( is_undo ? prev : selectRows, deselected );
    };
    auto undo = [=](){ undo_redo(true); };
    auto redo = [=](){ undo_redo(false); };
    undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Change flux tool row selection." );
  }//if( undoRedo && undoRedo->canAddUndoRedoNow() )

  // m_deselectedEnergies is intentionally left alone, so un-ticking and re-ticking gets the
  //  user back the rows they had picked.
  applySelectionState( selectRows, m_deselectedEnergies );
}//void setSelectRows( const bool selectRows, const bool force )


void FluxToolWidget::handleRowSelectionChanged( const size_t dataRow, const bool selected )
{
  PeakModel *peakmodel = m_interspec->peakModel();
  const vector<PeakDef> peaks = peakmodel ? peakmodel->peakVec() : vector<PeakDef>{};

  assert( dataRow < m_rowSelected.size() );
  assert( m_rowSelected.size() == m_data.size() );
  if( (dataRow >= m_rowSelected.size())
     || (m_rowSelected.size() != m_data.size())
     || (peaks.size() != m_data.size()) )
  {
    return;
  }

  if( m_rowSelected[dataRow] == selected )
    return;

  vector<double> deselected = m_deselectedEnergies;

  if( selected )
  {
    // Drop the remembered energy that names this row.  It need not equal the peaks current
    //  mean (the peak may have been re-fit since it was remembered), so match it the same
    //  way refreshPeakTable() does.  Entries that name some other peak - or no peak at all
    //  - are kept, so this set only ever changes through an explicit user action.
    vector<double> keep;
    for( const double energy : deselected )
    {
      const vector<bool> named = FluxToolImp::peaks_matching_energies( peaks, vector<double>{energy} );
      if( !named[dataRow] )
        keep.push_back( energy );
    }
    deselected.swap( keep );
  }else
  {
    deselected.push_back( m_data[dataRow][FluxColumns::FluxEnergyCol] );
  }

  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    const bool selectRows = m_selectRows;
    const vector<double> prev = m_deselectedEnergies;
    const vector<double> current = deselected;
    auto undo_redo = [selectRows, prev, current]( const bool is_undo ){
      InterSpec *viewer = InterSpec::instance();
      FluxToolWindow *fluxwin = viewer ? viewer->createFluxTool() : nullptr;
      FluxToolWidget *tool = fluxwin ? fluxwin->m_fluxTool : nullptr;
      if( tool )
        tool->applySelectionState( selectRows, is_undo ? prev : current );
    };
    auto undo = [=](){ undo_redo(true); };
    auto redo = [=](){ undo_redo(false); };
    undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Change flux tool selected row." );
  }//if( undoRedo && undoRedo->canAddUndoRedoNow() )

  applySelectionState( m_selectRows, deselected );
}//void handleRowSelectionChanged( const size_t dataRow, const bool selected )


void FluxToolWidget::applySelectionState( const bool selectRows, vector<double> deselected )
{
  // Taken by value: setSelectRows() passes m_deselectedEnergies itself, and we both assign
  //  to and then read that member below.
  const bool colVisibilityChanged = (selectRows != m_selectRows);

  m_selectRows = selectRows;
  m_deselectedEnergies = std::move( deselected );

  if( m_selectRowsCb->isChecked() != selectRows )
    m_selectRowsCb->setChecked( selectRows );

  // Re-derive which rows are checked, without rebuilding the table - the peaks and their
  //  flux values havent changed, only which of them the user wants in the output.
  //  A matching size means these are the same peaks m_data was built from: every peak
  //  change routes through setTableNeedsUpdating() (see init()), and Wt renders at the end
  //  of the event, so m_data can only be stale here if a refresh is already scheduled.
  PeakModel *peakmodel = m_interspec->peakModel();
  const vector<PeakDef> peaks = peakmodel ? peakmodel->peakVec() : vector<PeakDef>{};
  if( peaks.size() == m_data.size() )
    syncRowSelectionFromEnergies( peaks );
  else
    setTableNeedsUpdating();  //Table is stale; it will pick the selections up when rebuilt.

  if( colVisibilityChanged )
  {
    // Re-runs the column show/hide + width logic, which is where sm_selectColumn gets shown.
    setDisplayInfoLevel( m_displayInfoLevel, true );
  }

  if( m_fluxModel )
    m_fluxModel->notifySelectionChanged();

#if( FLUX_USE_COPY_TO_CLIPBOARD )
  // Deliberately not emitting m_tableUpdated - that resets the model, re-rendering (and
  //  re-sorting) the whole table, which is far too heavy for a checkbox click.
  updateCopyToClipboardText();
#endif
}//void applySelectionState( const bool selectRows, const vector<double> &deselected )


#if( FLUX_USE_COPY_TO_CLIPBOARD )
void FluxToolWidget::tableCopiedToCliboardCallback( const int copied )
{
  switch( copied )
  {
    case 0:
      passMessage( WString::tr("ftw-err-copy-clipboard"), 3 );
      break;
    
    case 1:
      passMessage( WString::tr("ftw-copied-to-clipboard"), 0 );
      break;
    
    case 2:
      passMessage( WString::tr("ftw-copied-as-html"), 0 );
      break;
      
    default:
      passMessage( WString::tr("ftw-copy-unknown-status"), 3 );
      break;
  }//switch( copied )
}//void tableCopiedToCliboardCallback( const int copied )
#endif

