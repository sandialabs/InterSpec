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
#include <memory>
#include <string>
#include <vector>
#include <numeric>
#include <assert.h>
#include <algorithm>

#include <Wt/WText.h>
#include <Wt/WTable.h>
#include <Wt/WLabel.h>
#include <Wt/WCheckBox.h>
#include <Wt/WResource.h>
#include <Wt/WLineEdit.h>
#include <Wt/WModelIndex.h>
#include <Wt/WGridLayout.h>
#include <Wt/WPushButton.h>
#include <Wt/WButtonGroup.h>
#include <Wt/WApplication.h>
#include <Wt/WEnvironment.h>
#include <Wt/WRadioButton.h>
#include <Wt/Http/Request.h>
#include <Wt/Http/Response.h>
#include <Wt/WRegExpValidator.h>
#include <Wt/WContainerWidget.h>
#include <Wt/cpp17/any.hpp>
#include <Wt/WAbstractItemModel.h>
#include <Wt/WItemDelegate.h>
#include <Wt/WAbstractItemDelegate.h>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/FluxTool.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/FluxToolCalc.h"
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
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/DetectorPeakResponse.h"

#if( USE_QR_CODES )
#include <Wt/Utils.h>

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
  var el = document.getElementById(id);
  var text = el && el._isData ? el._isData.TableData : null;
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
      return (m_sortOrder==Wt::SortOrder::Ascending ? less : !less);
    }
  };//struct index_compare
  
  
  /** Whether a column is shown, in the table and the CSV/clipboard output alike.

   @param bkgSub Whether background peaks are being subtracted (see
          FluxToolWidget::backgroundSubtractionActive).
   */
  bool showFluxColumn( FluxToolWidget::DisplayInfoLevel disptype,
                      const FluxToolWidget::FluxColumns col,
                      const shared_ptr<const DetectorPeakResponse> &det,
                      const bool bkgSub )
  {
    // A fixed-geometry efficiency describes the whole source-detector setup: there is no distance,
    //  and so no flux at the detector, or solid angle, to give.
    const bool fixed_geom = (det && det->isFixedGeometry());
    switch( col )
    {
      case FluxToolWidget::FluxColumns::FluxIntrinsicEffCol:
      case FluxToolWidget::FluxColumns::FluxGeometricEffCol:
      case FluxToolWidget::FluxColumns::FluxFluxOnDetCol:
      case FluxToolWidget::FluxColumns::FluxFluxPerCm2PerSCol:
        if( fixed_geom )
          return false;
        break;

      case FluxToolWidget::FluxColumns::FluxBackgroundCpsCol:
        return bkgSub && (disptype != FluxToolWidget::DisplayInfoLevel::Simple);

      default:
        break;
    }//switch( col )

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
          case FluxToolWidget::FluxColumns::FluxBackgroundCpsCol:
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

          case FluxToolWidget::FluxColumns::FluxBackgroundCpsCol:
            return bkgSub;
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
          case FluxToolWidget::FluxColumns::FluxIntrinsicEffCol:
          case FluxToolWidget::FluxColumns::FluxGeometricEffCol:
            return true;  //Fixed-geometry DRFs were taken care of above
            break;

          case FluxToolWidget::FluxColumns::FluxBackgroundCpsCol:
            return bkgSub;
            
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
    
    // Negative values are possible once background is subtracted, so the size is what matters
    const bool haveUncert = ((uncert > std::numeric_limits<double>::epsilon())
                             && (value != 0.0) && !(std::isinf)(value) && !(std::isnan)(value) );
    if( haveUncert )
    {
      //We'll be conservative and add at least one extra sig fig by using the nsigfigs to be
      //  number of digits after the decimal when printed in scientific notation.
      nsigfigs = static_cast<int>( std::round( std::log10( std::fabs(value) / uncert ) + 0.5 ) );
      nsigfigs = std::min( std::max( nsigfigs, 3 ), 9 ); //clamp num sig figs between 3 and 9
    }//if( haveUncert )
    
    snprintf( buffer, sizeof(buffer), "%.*G", nsigfigs, value );
    
    return buffer;
  }//print_value(...)


  std::string print_uncert( const double value, const double uncert, const bool uncertAsPercent )
  {
    const bool haveUncert = ((uncert > std::numeric_limits<double>::epsilon())
                             && (value != 0.0) && !(std::isinf)(value) && !(std::isnan)(value) );
    
    if( !haveUncert )
      return "";
    
    int nsigfigs = 6;
    char buffer[256] = { '\0' };
    
    
    //We'll be conservative and add at least one extra sig fig by using the nsigfigs to be
    //  number of digits after the decimal when printed in scientific notation.
    nsigfigs = static_cast<int>( std::round( std::log10( std::fabs(value) / uncert ) + 0.5 ) );
    nsigfigs = std::min( std::max( nsigfigs, 3 ), 9 ); //clamp num sig figs between 3 and 9
    
    if( uncertAsPercent )
    {
      const double percentUncert = 100.0 * uncert / std::fabs(value);
      
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


  /** An uncertainty split into its statistical and systematic (detector efficiency) parts, like
   "0.401 (stat) +- 0.112 (syst)", as the exports and the table's tooltips give it.  The systematic
   part is only given when there is one (the note under the table says when the detector states
   none); peak (and background) CPS never has one, and intrinsic efficiency has no statistical part.
   Empty when `value` is zero, or not a finite number.
   */
  std::string print_stat_syst( const FluxToolWidget::FluxColumns col, const double value,
                               const double stat, const double syst, const bool asPercent,
                               const std::string &plus_minus )
  {
    if( (value == 0.0) || (std::isinf)(value) || (std::isnan)(value) )
      return "";

    const auto part = [=]( const double uncert ) -> std::string {
      const std::string txt = print_uncert( value, uncert, asPercent );
      return txt.empty() ? std::string(asPercent ? "0%" : "0") : txt;
    };

    std::string answer;
    if( col != FluxToolWidget::FluxColumns::FluxIntrinsicEffCol )
      answer = part( stat ) + " (stat)";
    if( (col != FluxToolWidget::FluxColumns::FluxPeakCpsCol)
       && (col != FluxToolWidget::FluxColumns::FluxBackgroundCpsCol)
       && (syst > std::numeric_limits<double>::epsilon()) )
      answer += (answer.empty() ? "" : (" " + plus_minus + " ")) + part( syst ) + " (syst)";

    return answer;
  }//print_stat_syst(...)


  /** The energy range the DRF states it is valid for, like "50-3000", for messages; empty if it states
   none (see FluxToolCalc::outside_energy_range).
   */
  std::string drf_energy_range_txt( const shared_ptr<const DetectorPeakResponse> &det )
  {
    if( !det || !(det->upperEnergy() > det->lowerEnergy()) )
      return "";
    char buffer[64] = { '\0' };
    snprintf( buffer, sizeof(buffer), "%.6g-%.6g", det->lowerEnergy(), det->upperEnergy() );
    return buffer;
  }//drf_energy_range_txt(...)


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
    : WAbstractItemModel(),
      m_fluxtool( parent ),
      m_sortColumn( FluxToolWidget::FluxColumns::FluxEnergyCol ),
      m_sortOrder( Wt::SortOrder::Ascending )
    {
      parent->tableUpdated().connect( this, &FluxModel::handleFluxToolWidgetUpdated );
    }
    
    
    virtual ~FluxModel()
    {
    }
    
    
    virtual Wt::WFlags<Wt::ItemFlag> flags( const Wt::WModelIndex &index ) const
    {
      if( index.column() == FluxToolWidget::sm_selectColumn )
        return Wt::WFlags<Wt::ItemFlag>(Wt::ItemFlag::UserCheckable);

      return Wt::WFlags<Wt::ItemFlag>(Wt::ItemFlag::XHTMLText);
    }
    
    
    virtual WFlags<Wt::HeaderFlag> headerFlags( int section, Wt::Orientation orientation = Wt::Orientation::Horizontal ) const
    {
      return Wt::WFlags<Wt::HeaderFlag>(Wt::HeaderFlag::XHTMLText);
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
      
      // Deliberately NO reset() here.  The remove/insert-rows signals above already refresh the
      //  view, and under Wt4 `WAbstractItemModel::reset()` -> `WAbstractItemView::modelReset()`
      //  is literally `setModel(model_)`: it re-connects every model signal (and
      //  RowStretchTreeView::setModel's rowsInserted/rowsRemoved callbacks) a second time, and the
      //  header re-render it forces leaves the *previous* header DOM behind on the client - so the
      //  table visibly grows a duplicate set of column headers and rows on the first info-level or
      //  "select rows" change.  (Under Wt3 this was a harmless "just in case".)
    }//void handleFluxToolWidgetUpdated()
    
    
    virtual Wt::cpp17::any data( const Wt::WModelIndex &index,
                             Wt::ItemDataRole role = Wt::ItemDataRole::Display ) const
    {
      if( index.parent().isValid() )
        return Wt::cpp17::any();
      
      const int row = index.row();
      const int col = index.column();

      if( row < 0 || col < 0
         || col > FluxToolWidget::sm_selectColumn
         || row >= static_cast<int>(m_sort_indices.size()) )
        return Wt::cpp17::any();

      const size_t realind = m_sort_indices[row];

      const auto &nucs = m_fluxtool->m_nucNames;
      const auto &data = m_fluxtool->m_data;
      const auto &stat_uncerts = m_fluxtool->m_statUncerts;
      const auto &syst_uncerts = m_fluxtool->m_systUncerts;

      assert( nucs.size() == data.size() );
      assert( stat_uncerts.size() == data.size() );
      assert( syst_uncerts.size() == data.size() );
      assert( realind < data.size() );

      // The selection column has no value/uncertainty in m_data; it is checkbox-only.
      //  Note that Wt::WItemDelegate decides whether to render a checkbox at all based on
      //  CheckStateRole being non-empty here - the ItemIsUserCheckable flag only controls
      //  if that checkbox is enabled - so the other columns must keep returning an empty
      //  Wt::cpp17::any for CheckStateRole (via the `default:` case, below).
      if( col == FluxToolWidget::sm_selectColumn )
      {
        if( role != Wt::ItemDataRole::Checked )
          return Wt::cpp17::any();

        assert( m_fluxtool->m_rowSelected.size() == data.size() );
        if( realind >= m_fluxtool->m_rowSelected.size() )
          return Wt::cpp17::any();

        return Wt::cpp17::any( static_cast<bool>(m_fluxtool->m_rowSelected[realind]) );
      }//if( col == FluxToolWidget::sm_selectColumn )

      // Wt4's ItemDataRole is a class (not an integer enum), so this cannot be a switch.
      if( role == Wt::ItemDataRole::Display )
      {
        if( col == FluxToolWidget::FluxColumns::FluxNuclideCol )
          return nucs[realind].empty() ? Wt::cpp17::any() : Wt::cpp17::any( WString::fromUTF8(nucs[realind]) );
        return Wt::cpp17::any( data[realind][col] );
      }//if( role == Wt::ItemDataRole::Display )

      if( col == FluxToolWidget::FluxColumns::FluxNuclideCol )
        return Wt::cpp17::any();

      // Rows outside the DRFs stated energy range get their energy flagged
      if( col == FluxToolWidget::FluxColumns::FluxEnergyCol )
      {
        const vector<bool> &out_of_range = m_fluxtool->m_rowOutOfRange;
        if( (realind >= out_of_range.size()) || !out_of_range[realind] )
          return Wt::cpp17::any();

        if( role == Wt::ItemDataRole::StyleClass )
          return Wt::cpp17::any( WString::fromUTF8("FluxOutOfRange") );

        if( role == Wt::ItemDataRole::ToolTip )
        {
          const string range = drf_energy_range_txt( m_fluxtool->m_detector->detector() );
          return Wt::cpp17::any( WString::tr("ftw-tt-outside-range").arg( range ) );
        }

        return Wt::cpp17::any();
      }//if( col == FluxToolWidget::FluxColumns::FluxEnergyCol )

      const double stat = stat_uncerts[realind][col], syst = syst_uncerts[realind][col];
      const double total = std::sqrt( stat*stat + syst*syst );
      const bool has_uncert = (total > std::numeric_limits<double>::epsilon());

      // The combined uncertainty, which the cells show
      if( role == Wt::ItemDataRole::User )
        return has_uncert ? Wt::cpp17::any( total ) : Wt::cpp17::any();

      // The statistical / systematic split of it
      if( (role == Wt::ItemDataRole::ToolTip) && has_uncert )
      {
        const bool asPercent = (m_fluxtool->displayInfoLevel() == FluxToolWidget::DisplayInfoLevel::Simple);
        const FluxToolWidget::FluxColumns fcol = static_cast<FluxToolWidget::FluxColumns>( col );
        const double value = data[realind][col];
        const std::string split = print_stat_syst( fcol, value, stat, syst, asPercent, "\xC2\xB1" );
        if( !split.empty() )
          return Wt::cpp17::any( WString::fromUTF8( print_value( value, total ) + " \xC2\xB1 " + split ) );
      }//if( role == Wt::ItemDataRole::ToolTip )

      return Wt::cpp17::any();
    }//data(...)


    /** Handles the user clicking a rows checkbox.

     Wt::WItemDelegate::onCheckedChange sends a bool (rather than a Wt::CheckState), since
     we dont set ItemIsTristate.
     */
    virtual bool setData( const Wt::WModelIndex &index, const Wt::cpp17::any &value,
                          Wt::ItemDataRole role = Wt::ItemDataRole::Edit )
    {
      if( index.parent().isValid()
         || (index.column() != FluxToolWidget::sm_selectColumn)
         || (role != Wt::ItemDataRole::Checked)
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

      m_fluxtool->handleRowSelectionChanged( realind, Wt::cpp17::any_cast<bool>(value) );

      dataChanged().emit( index, index );

      return true;
    }//bool setData(...)


    /** Tells the view every rows checkbox may have changed.

     Used when selection changes from something other than a click on that checkbox (an
     undo/redo, the "select rows" checkbox, or an app-URL), so the view redraws the checkboxes
     without rebuilding every row the way #handleFluxToolWidgetUpdated does.
     */
    void notifySelectionChanged()
    {
      if( m_sort_indices.empty() )
        return;

      const int lastrow = static_cast<int>(m_sort_indices.size()) - 1;
      dataChanged().emit( index(0, FluxToolWidget::sm_selectColumn),
                          index(lastrow, FluxToolWidget::sm_selectColumn) );
    }//void notifySelectionChanged()
    
    
    Wt::cpp17::any headerData( int section,
                          Wt::Orientation orientation = Wt::Orientation::Horizontal,
                          Wt::ItemDataRole role = Wt::ItemDataRole::Display ) const
    {
      if( section < 0 || section > FluxToolWidget::sm_selectColumn
         || orientation != Wt::Orientation::Horizontal
         || ((role != Wt::ItemDataRole::Display) && (role != Wt::ItemDataRole::ToolTip)) )
        return Wt::cpp17::any();

      if( role == Wt::ItemDataRole::ToolTip )
      {
        // What the intrinsic and geometric columns mean - the intrinsic one is at this distance, so
        //  for a detector with a characterized response it is not its far-field curve.
        switch( section )
        {
          case FluxToolWidget::FluxColumns::FluxIntrinsicEffCol:
            return Wt::cpp17::any( WString::tr("ftw-tt-hdr-intrinsic-eff") );
          case FluxToolWidget::FluxColumns::FluxGeometricEffCol:
            return Wt::cpp17::any( WString::tr("ftw-tt-hdr-geom-eff") );
          case FluxToolWidget::FluxColumns::FluxBackgroundCpsCol:
            return Wt::cpp17::any( WString::tr("ftw-tt-hdr-bkg-cps") );
          case FluxToolWidget::FluxColumns::FluxFluxOnDetCol:
            return Wt::cpp17::any( WString::tr("ftw-tt-hdr-flux-on-det") );
          case FluxToolWidget::FluxColumns::FluxFluxPerCm2PerSCol:
            return Wt::cpp17::any( WString::tr("ftw-tt-hdr-flux-cm") );
          case FluxToolWidget::FluxColumns::FluxGammasInto4PiCol:
          {
            const shared_ptr<const DetectorPeakResponse> det = m_fluxtool->m_detector->detector();
            const bool fixed_geom = (det && det->isFixedGeometry());
            return Wt::cpp17::any( WString::tr( fixed_geom ? "ftw-tt-hdr-fixed-geom-rate" : "ftw-tt-hdr-flux-4pi") );
          }
          default:
            return Wt::cpp17::any();
        }//switch( section )
      }//if( role == ToolTip )

      if( section == FluxToolWidget::sm_selectColumn )
        return Wt::cpp17::any( m_fluxtool->m_selectColName );

      return Wt::cpp17::any( m_fluxtool->m_colnames[section] );
    }

    virtual int columnCount( const Wt::WModelIndex &parent
                            = Wt::WModelIndex() ) const
    {
      if( parent.isValid() )
        return 0;

      // Always report the selection column, even when rows arent being selected - the view
      //  hides it instead (see FluxToolWidget::setDisplayInfoLevel).  Wt::WAbstractItemView
      //  only truncates/appends its column list when the model column count changes, and never
      //  re-initializes the entries that survive, so a varying column count would leave the
      //  ColumnInfo silently recreated with a default width/sorting/delegate the next time rows
      //  became selectable.
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
    
    virtual void sort( int column, Wt::SortOrder order = Wt::SortOrder::Ascending )
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
    : WAbstractItemDelegate(),
      m_fluxTool( parent )
    {
    }


    virtual std::unique_ptr<WWidget> update( WWidget *widget, const WModelIndex& index, WFlags<ViewItemRenderFlag> flags )
    {
      WText *oldwidget = dynamic_cast<WText *>( widget );

      std::unique_ptr<WText> newWidget;
      if( oldwidget )
      {
        oldwidget->setText("");
      }else
      {
        newWidget = std::make_unique<WText>();
        oldwidget = newWidget.get();
      }


      const auto model = index.model();
      if( !model )
        return newWidget;

      const auto data = index.data();
      const auto uncertdata = model->data( index.row(), index.column(), Wt::ItemDataRole::User );

      if( !Wt::cpp17::any_has_value(data) )
        return newWidget;

      WString valstr;
      if( index.column() == FluxToolWidget::FluxColumns::FluxNuclideCol )
      {
        valstr = Wt::asString(data);
      }else
      {
        const bool hasUncert = Wt::cpp17::any_has_value(uncertdata);

        double val = 0.0, uncert = 0.0;
        try
        {
          val = Wt::cpp17::any_cast<double>(data);
          if( hasUncert )
            uncert = Wt::cpp17::any_cast<double>(uncertdata);
        }catch(...)
        {
          return newWidget;
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

      const auto tooltip = model->data( index.row(), index.column(), Wt::ItemDataRole::ToolTip );
      oldwidget->setToolTip( Wt::cpp17::any_has_value(tooltip) ? Wt::asString(tooltip) : WString() );

      const auto style = model->data( index.row(), index.column(), Wt::ItemDataRole::StyleClass );
      oldwidget->toggleStyleClass( "FluxOutOfRange", Wt::cpp17::any_has_value(style) );

      return newWidget;
    }//std::unique_ptr<WWidget> update(...)
  };//class FluxRenderDelegate


  
  class FluxCsvResource : public Wt::WResource
  {
  protected:
    FluxToolWidget *m_fluxtool;
    Wt::WApplication *m_app;
    
  public:
    FluxCsvResource( FluxToolWidget *parent )
    : WResource(),
      m_fluxtool( parent ),
      m_app( WApplication::instance() )
    {
      assert( m_app );

      // Have Wt hold the session lock across handleRequest instead of dropping it for us: by
      //  default `WResource::handle` unlocks before invoking us, so the session thread can be
      //  inside `~WResource`'s `beingDeleted()` - waiting on our use-count - while we block
      //  acquiring the lock it holds.  That is a deadlock.  Free: the `WApplication::UpdateLock`
      //  in handleRequest already spanned the whole body, so this only moves it earlier.
      setTakesUpdateLock( true );
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
      const bool bkgSub = tool->backgroundSubtractionActive();

      // Rows that will be output; when the user is picking rows, only the checked ones.
      //  m_rowSelected parallels m_data, both in un-sorted order, so it indexes straight off the row.
      //  The size check is just belt-and-braces (refreshPeakTable() keeps the two the same size,
      //  including on its early-returns, where both end up empty).
      assert( tool->m_rowSelected.size() == tool->m_data.size() );
      const auto output_row = [tool]( const size_t row ) -> bool {
        return !tool->m_selectRows
               || (row >= tool->m_rowSelected.size())
               || tool->m_rowSelected[row];
      };

      // A "Notes" column is added only when some row output needs one.
      bool any_notes = false;
      for( size_t row = 0; row < tool->m_data.size(); ++row )
        any_notes |= (output_row(row) && (row < tool->m_rowOutOfRange.size()) && tool->m_rowOutOfRange[row]);
      const string range_note = any_notes ? ("outside DRF energy range (" + drf_energy_range_txt(det) + " keV)") : string();
      
      for( FluxToolWidget::FluxColumns col = FluxToolWidget::FluxColumns(0);
          col < FluxToolWidget::FluxColumns::FluxNumColumns; col = FluxToolWidget::FluxColumns(col+1) )
      {
        const WString &colname = tool->m_colnamesCsv[col];
        
        if( !FluxToolImp::showFluxColumn(disptype, col, det, bkgSub) )
          continue;
        
        if( html )
          strm << (col==0 ? "\\t<tr><th>" : "</th><th>") << colname;
        else
          strm << (col==0 ? "" : ",") << colname;
        
        //No uncertainty on energy, geometric eff, or nuc; the intrinsic efficiency's is the DRF's
        switch( col )
        {
          case FluxToolWidget::FluxEnergyCol:
          case FluxToolWidget::FluxGeometricEffCol:
          case FluxToolWidget::FluxNuclideCol:
          case FluxToolWidget::FluxNumColumns:
            break;
            
          case FluxToolWidget::FluxIntrinsicEffCol:
          case FluxToolWidget::FluxPeakCpsCol:
          case FluxToolWidget::FluxBackgroundCpsCol:
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

      if( any_notes )
        strm << (html ? "</th><th>" : ",") << "Notes";
      
      if( html )
        strm << "</th></tr>";
      strm << eol_char;
      
      for( size_t row = 0; row < tool->m_data.size(); ++row )
      {
        if( !output_row(row) )
          continue;

        if( html )
          strm << "\\t<tr>";

        for( FluxToolWidget::FluxColumns col = FluxToolWidget::FluxColumns(0);
            col < FluxToolWidget::FluxColumns::FluxNumColumns; col = FluxToolWidget::FluxColumns(col+1) )
        {
          if( !FluxToolImp::showFluxColumn(disptype, col, det, bkgSub) )
            continue;
          
          const double data = tool->m_data[row][col];
          const double stat = tool->m_statUncerts[row][col];
          const double syst = tool->m_systUncerts[row][col];
          const double uncert = std::sqrt( stat*stat + syst*syst );
          
          
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
            case FluxToolWidget::FluxBackgroundCpsCol:
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
            case FluxToolWidget::FluxGeometricEffCol:
            case FluxToolWidget::FluxNuclideCol:
            case FluxToolWidget::FluxNumColumns:
              break;
              
            case FluxToolWidget::FluxIntrinsicEffCol:
            case FluxToolWidget::FluxPeakCpsCol:
            case FluxToolWidget::FluxBackgroundCpsCol:
            case FluxToolWidget::FluxFluxOnDetCol:
            case FluxToolWidget::FluxFluxPerCm2PerSCol:
            case FluxToolWidget::FluxGammasInto4PiCol:
            {
              // e.g. "0.401 (stat) +- 0.112 (syst)" - percents in the Simple view
              const bool asPercent = (disptype == FluxToolWidget::DisplayInfoLevel::Simple);
              strm << (html ? "</td><td>" : ",")
                   << print_stat_syst( col, data, stat, syst, asPercent, (html ? "&plusmn;" : "+-") );
              
              break;
            }//case CPS, FluxOnDet, FluxPerArea, GammasInto 4pi
          }//switch( col )
        }//for( loop over columns )

        if( any_notes )
        {
          const bool out_of_range = (row < tool->m_rowOutOfRange.size()) && tool->m_rowOutOfRange[row];
          strm << (html ? "</td><td>" : ",") << (out_of_range ? range_note : string());
        }
        
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
      
      suggestFileName( filename, ContentDisposition::Attachment );
      
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

  {
    auto fluxToolOwned = std::make_unique<FluxToolWidget>( viewer );
    m_fluxTool = fluxToolOwned.get();
    stretcher()->addWidget( std::move(fluxToolOwned), 0, 0 );
  }
  //m_fluxTool->setHeight( WLength(100, WLength::Unit::Percentage) );

  AuxWindow::addHelpInFooter( footer(), "flux-tool" );

  WPushButton *closeButton = nullptr;
  if( viewer && viewer->isPhone() )
    closeButton = addCloseButtonToFooter();

#if( USE_QR_CODES )
  WPushButton *qr_btn = footer()->addNew<WPushButton>();
  qr_btn->setText( WString::tr("QR Code") );
  qr_btn->setIcon( "InterSpec_resources/images/qr-code.svg" );
  qr_btn->setStyleClass( "LinkBtn DownloadBtn DialogFooterQrBtn" );
  qr_btn->clicked().preventPropagation();
  qr_btn->clicked().connect( this, [this](){
    try
    {
      const string url = "interspec://flux/?" + Wt::Utils::urlEncode(m_fluxTool->encodeStateToUrl());
      QrCode::displayTxtAsQrCode( url, WString::tr("ftw-qr-tool-state-title"),
                                 WString::tr("ftw-qr-tool-state-text") );
    }catch( std::exception &e )
    {
      passMessage( WString::tr("app-qr-err").arg(e.what()), WarningWidget::WarningMsgHigh );
    }
  } );
#endif //USE_QR_CODES

  WContainerWidget *buttonDiv = footer();

  if( !viewer || !viewer->isPhone() )
    closeButton = addCloseButtonToFooter();

  assert( closeButton );
  closeButton->clicked().connect( this, &AuxWindow::hide );

  auto csv = std::make_shared<FluxToolImp::FluxCsvResource>( m_fluxTool );
#if( BUILD_AS_OSX_APP || IOS )
  WLink csvLink( csv );
  csvLink.setTarget( LinkTarget::NewWindow );
  WAnchor *csvButton = buttonDiv->addNew<WAnchor>( csvLink, WString::tr("CSV") );
  csvButton->setStyleClass( "LinkBtn DownloadLink" );
#else
  WPushButton *csvButton = buttonDiv->addNew<WPushButton>();
  csvButton->setIcon( "InterSpec_resources/images/download_small.svg" );
  WLink csvLink( csv );
  csvLink.setTarget( LinkTarget::NewWindow );
  csvButton->setLink( csvLink );
  csvButton->setStyleClass( "LinkBtn DownloadBtn" );
#endif

  csvButton->setText( WString::tr("CSV") );
  //csvButton->setAttributeValue( "style", "float: none;" );  //Keep the CSV download to the left side of the close button.  .LinkBtn style class has the button float right..
  csvButton->setFloatSide( Wt::Side::Left );
  
  auto enableDisableCsv = [csvButton,this](){
    csvButton->setDisabled( m_fluxTool->m_data.empty() );
  };
  
  m_fluxTool->m_tableUpdated.connect( this, enableDisableCsv );
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


FluxToolWidget::FluxToolWidget( InterSpec *viewer )
  : WContainerWidget(),
    m_interspec( viewer ),
    m_detector( nullptr ),
    m_narrowLayout( false ),
    m_msg( nullptr ),
    m_drfUncertNote( nullptr ),
    m_calcNote( nullptr ),
    m_rangeWarning( nullptr ),
    m_airAttenCb( nullptr ),
    m_attenuateForAir( true ),
    m_backgroundPeakSub( nullptr ),
    m_subtractBackground( false ),
    m_skipBackgroundPeakCheck( false ),
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
    m_tableUpdated()
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

  // URIs before version 1.2 predate the air-attenuation option, so had none applied.
  const bool pre_ver_1_2 = (ver_iter == end(parts))
                           || (ver_iter->second == "1")
                           || (ver_iter->second == "1.0")
                           || (ver_iter->second == "1.1");
  
  const auto dist_iter = parts.find( "DIST" );
  if( dist_iter == end(parts) )
    throw runtime_error( "URL is invalid Flux Tool URI = no 'DIST' component" );
  
  // encodeStateToUrl() writes an empty distance when the distance box held an invalid one (which
  //  for a fixed-geometry DRF isnt used, so may be left that way); keep the current distance then,
  //  so the rest of the state - e.g., when undoing closing the tool - still gets restored.
  string dist = dist_iter->second;
  SpecUtils::trim( dist );
  if( dist.empty() )
  {
    dist = m_distance->text().toUTF8();
  }else
  {
    try
    {
      PhysicalUnits::stringToDistance( dist );
    }catch( std::exception & )
    {
      Wt::log("warn") << "FluxToolWidget URI had invalid distance: '" << dist << "'";
      throw std::runtime_error( "Invalid Flux Tool distance specified in URI" );
    }
  }//if( dist.empty() ) / else
  
  
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

    const std::shared_ptr<PeakModel> peakmodel = m_interspec->peakModel();
    const vector<PeakDef> peaks = peakmodel ? peakmodel->peakVec() : vector<PeakDef>{};
    const vector<bool> is_selected = FluxToolImp::peaks_matching_energies( peaks, selected );

    for( size_t i = 0; i < peaks.size(); ++i )
    {
      if( !is_selected[i] )
        deselected.push_back( peaks[i].mean() );
    }
  }//if( sel_iter != end(parts) )

  // "AIR" and "BKGSUB" were added in URI version 1.2; see the default for an absent "AIR" above.
  bool attenuateForAir = !pre_ver_1_2;
  const auto air_iter = parts.find( "AIR" );
  if( air_iter != end(parts) )
    attenuateForAir = (air_iter->second == "1") || SpecUtils::iequals_ascii(air_iter->second, "true");

  bool subtractBackground = false;
  const auto bkgsub_iter = parts.find( "BKGSUB" );
  if( bkgsub_iter != end(parts) )
    subtractBackground = (bkgsub_iter->second == "1") || SpecUtils::iequals_ascii(bkgsub_iter->second, "true");

  applyAttenuateForAir( attenuateForAir );
  applySubtractBackground( subtractBackground );
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
  
  string answer = "VER=1.2&DIST=" + dist + "&DISPLAY=";

  switch( displayInfoLevel() )
  {
    case DisplayInfoLevel::Simple:   answer += "SIMPLE";   break;
    case DisplayInfoLevel::Normal:   answer += "NORMAL";   break;
    case DisplayInfoLevel::Extended: answer += "EXTENDED"; break;
  }//switch( displayInfoLevel() )

  // Only written when not the default (on), to keep the URI short.
  if( !m_attenuateForAir )
    answer += "&AIR=0";

  if( m_subtractBackground )
    answer += "&BKGSUB=1";

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
    const std::shared_ptr<PeakModel> peakmodel = m_interspec ? m_interspec->peakModel() : nullptr;
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
  
  updateColumnLabels();
  
  
  wApp->useStyleSheet( "InterSpec_resources/FluxTool.css" );
  
  const bool showToolTips = m_interspec ? UserPreferences::preferenceValue<bool>( "ShowTooltips", m_interspec ) : false;
  
  addStyleClass( "FluxToolWidget" );
  
  WGridLayout *layout = setLayout( std::make_unique<WGridLayout>() );
  layout->setContentsMargins( 0, 0, 0, 0 );

  WTable *distDetRow = layout->addWidget( std::make_unique<WTable>(), 0, 0, 1, 2 );
  distDetRow->addStyleClass( "FluxDistMsgDetTbl" );

  
  auto distCell = distDetRow->elementAt( m_narrowLayout ? 1 : 0, 0 );
  distCell->addStyleClass( "FluxDistCell" );
  WLabel *label = distCell->addNew<WLabel>( WString::tr("distance-label") );
  label->addStyleClass( "FluxDistLabel" );

  m_prevDistance = "100 cm";
  m_distance = distCell->addNew<WLineEdit>( m_prevDistance );
  m_distance->addStyleClass( "FluxDistanceEnter" );
  label->setBuddy(m_distance);

  m_distance->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_distance->setAttributeValue( "autocorrect", "off" );
  m_distance->setAttributeValue( "spellcheck", "off" );
#endif

  auto validator = std::make_shared<WRegExpValidator>( PhysicalUnits::sm_distanceUnitOptionalRegex );
  validator->setFlags( Wt::RegExpFlag::MatchCaseInsensitive );
  m_distance->setValidator( validator );
  HelpSystem::attachToolTipOn( m_distance, WString::tr("ftw-tt-distance"), showToolTips );
  m_distance->changed().connect( this, &FluxToolWidget::distanceUpdated );
  m_distance->enterPressed().connect( this, &FluxToolWidget::distanceUpdated );

  SpectraFileModel *specFileModel = m_interspec->fileManager()->model();
  if( m_narrowLayout )
  {
    auto detCell = distDetRow->elementAt(0,0);
    detCell->setColumnSpan( 2 );
    m_detector = detCell->addNew<DetectorDisplay>( m_interspec, specFileModel );
  }else
  {
    m_detector = distDetRow->elementAt(0,1)->addNew<DetectorDisplay>( m_interspec, specFileModel );
  }
  m_detector->addStyleClass( "FluxDet" );
  m_interspec->detectorChanged().connect( this, [this]( std::shared_ptr<DetectorPeakResponse> a1 ){ handleDrfChange( a1 ); } );
  m_interspec->detectorModified().connect( this, [this]( std::shared_ptr<DetectorPeakResponse> a1 ){ handleDrfChange( a1 ); } );

  // The calculation options, and the notes on how the numbers were computed, go under the table (and
  //  under the row of output options right below it)
  WContainerWidget *optsAndNotesDiv = nullptr;
  {
    auto optsAndNotesOwned = std::make_unique<WContainerWidget>();
    optsAndNotesDiv = optsAndNotesOwned.get();
    optsAndNotesDiv->addStyleClass( "FluxCalcOptsAndNotes" );
    layout->addWidget( std::move(optsAndNotesOwned), 3, 0, 1, 2 );
  }

  WContainerWidget *calcOptsDiv = optsAndNotesDiv->addNew<WContainerWidget>();
  calcOptsDiv->addStyleClass( "FluxCalcOpts" );

  m_airAttenCb = calcOptsDiv->addNew<WCheckBox>( WString::tr( m_narrowLayout ? "ftw-cb-air-narrow" : "ftw-cb-air" ) );
  m_airAttenCb->addStyleClass( "CbNoLineBreak FluxAirCb" );
  m_airAttenCb->setChecked( m_attenuateForAir );
  HelpSystem::attachToolTipOn( m_airAttenCb, WString::tr("ftw-tt-air"), showToolTips );
  m_airAttenCb->checked().connect( this, [this](){ setAttenuateForAir( true ); } );
  m_airAttenCb->unChecked().connect( this, [this](){ setAttenuateForAir( false ); } );

  m_backgroundPeakSub = calcOptsDiv->addNew<WCheckBox>( WString::tr( m_narrowLayout ? "ftw-cb-bkg-sub-narrow" : "ftw-cb-bkg-sub" ) );
  m_backgroundPeakSub->addStyleClass( "CbNoLineBreak FluxBkgSubCb" );
  m_backgroundPeakSub->setChecked( m_subtractBackground );
  m_backgroundPeakSub->setHidden( true );  //see updateBackgroundSubVisibility()
  HelpSystem::attachToolTipOn( m_backgroundPeakSub, WString::tr("ftw-tt-bkg-sub"), showToolTips );
  m_backgroundPeakSub->checked().connect( this, [this](){ setSubtractBackground( true ); } );
  m_backgroundPeakSub->unChecked().connect( this, [this](){ setSubtractBackground( false ); } );

  // A new foreground means new peaks, so it is worth offering missing background peaks again.
  m_interspec->displayedSpectrumChanged().connect( this, [this]( SpecUtils::SpectrumType type,
                                                        std::shared_ptr<SpecMeas>, std::set<int>,
                                                        std::vector<std::string> ){
    if( type == SpecUtils::SpectrumType::Foreground )
      m_skipBackgroundPeakCheck = false;
    setTableNeedsUpdating();
  } );
  m_interspec->spectrumScaleFactorChanged().connect( this, [this]( SpecUtils::SpectrumType, double ){
    setTableNeedsUpdating();
  } );
  m_interspec->hintPeaksSet().connect( this, [this]( SpecUtils::SpectrumType type ){
    if( type == SpecUtils::SpectrumType::Background )
      setTableNeedsUpdating();  //can change whether background subtraction is offered
  } );


  const std::shared_ptr<PeakModel> peakmodel = m_interspec->peakModel();
  assert( peakmodel );
  peakmodel->dataChanged().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  peakmodel->rowsRemoved().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  peakmodel->rowsInserted().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  peakmodel->layoutChanged().connect( this, &FluxToolWidget::setTableNeedsUpdating );
  peakmodel->backgroundPeaksChanged().connect( this, &FluxToolWidget::setTableNeedsUpdating );

  // Why the table is empty, if it is, stays prominent at the top
  auto msgCell = distDetRow->elementAt( m_narrowLayout ? 2 : 1, 0 );
  msgCell->setColumnSpan( 2 );
  msgCell->addStyleClass( "FluxMsgCell" );
  m_msg = msgCell->addNew<WText>( WString::fromUTF8(""), Wt::TextFormat::XHTML );
  m_msg->addStyleClass( "FluxMsg" );

  WContainerWidget *notesDiv = optsAndNotesDiv->addNew<WContainerWidget>();
  notesDiv->addStyleClass( "FluxNotes" );
  m_rangeWarning = notesDiv->addNew<WText>( WString::fromUTF8("") );
  m_rangeWarning->addStyleClass( "FluxRangeWarn" );
  m_rangeWarning->setInline( false );
  m_drfUncertNote = notesDiv->addNew<WText>( WString::fromUTF8("") );
  m_drfUncertNote->addStyleClass( "FluxDrfUncertNote" );
  m_drfUncertNote->setInline( false );
  m_calcNote = notesDiv->addNew<WText>( WString::fromUTF8(""), Wt::TextFormat::XHTML );
  m_calcNote->addStyleClass( "FluxDrfUncertNote FluxCalcNote" );
  m_calcNote->setInline( false );
  
  
  m_selectColName = WString::tr("ftw-hdr-use");

  // Wt4 owns models and item delegates through shared_ptr; `m_fluxModel` is just an observer
  //  (the table's setModel() keeps the model alive).
  auto fluxModelOwned = std::make_shared<FluxToolImp::FluxModel>( this );
  m_fluxModel = fluxModelOwned.get();

  {
    auto tbl = std::make_unique<RowStretchTreeView>();
    m_table = tbl.get();
#if( FLUX_USE_COPY_TO_CLIPBOARD )
    layout->addWidget( std::move(tbl), 1, 0, 1, 2 );
#else
    layout->addWidget( std::move(tbl), 1, 0, 1, 1 );
#endif
  }

  m_table->setModel( fluxModelOwned );
  m_table->setRootIsDecorated( false );
  m_table->addStyleClass( "FluxTable" );
  m_table->setAlternatingRowColors( true );
  m_table->sortByColumn( FluxColumns::FluxEnergyCol, Wt::SortOrder::Ascending );

  //Setting rows selectable doesnt seem to work... not that it would do us
  //  any good anyway since you probably want to copy the info to the clipboard.
  //m_table->setSelectable( true );
  //m_table->setSelectionMode( Wt::SelectionMode::SingleSelection );

  auto renderdel = std::make_shared<FluxToolImp::FluxRenderDelegate>( this );
  m_table->setItemDelegate( renderdel );

  // The stock delegate knows how to render a checkbox from CheckStateRole; using it just
  //  for the selection column saves teaching FluxRenderDelegate about checkboxes.
  m_table->setItemDelegateForColumn( sm_selectColumn, std::make_shared<WItemDelegate>() );

  for( FluxColumns col = FluxColumns(0); col < FluxColumns::FluxNumColumns; col = FluxColumns(col + 1) )
    m_table->setSortingEnabled( col, (col!=FluxColumns::FluxGeometricEffCol) );

  // Sorting is per-column, and defaults to enabled, so the selection column needs this
  //  explicitly - the loop above only covers the data columns.  See FluxModel::sort().
  m_table->setSortingEnabled( sm_selectColumn, false );

  layout->setRowStretch( 1, 1 );

  // Immediately under the table: row selection on the left, copying in the middle, and how much
  //  detail to show on the right - wrapping onto more lines if there isnt room (e.g., phones).
  WContainerWidget *outputOptsDiv = nullptr;
  {
    auto outputOptsOwned = std::make_unique<WContainerWidget>();
    outputOptsDiv = outputOptsOwned.get();
    outputOptsDiv->addStyleClass( "FluxOutputOpts" );
    layout->addWidget( std::move(outputOptsOwned), 2, 0, 1, 2 );
  }

  m_selectRowsCb = outputOptsDiv->addNew<WCheckBox>( WString::tr( m_narrowLayout ? "ftw-select-rows-narrow"
                                                              : "ftw-select-rows") );
  m_selectRowsCb->addStyleClass( "CbNoLineBreak FluxSelectRowsCb" );
  HelpSystem::attachToolTipOn( m_selectRowsCb, WString::tr("ftw-tt-select-rows"), showToolTips );
  m_selectRowsCb->checked().connect( this, [this](){ setSelectRows( true ); } );
  m_selectRowsCb->unChecked().connect( this, [this](){ setSelectRows( false ); } );

#if( FLUX_USE_COPY_TO_CLIPBOARD )
  LOAD_JAVASCRIPT(wApp, "FluxTool.cpp", "FluxTool", wtjsCopyFluxDataTextToClipboard );

  m_copyBtn = outputOptsDiv->addNew<WPushButton>( WString::tr("ftw-copy-btn-narrow") );
  m_copyBtn->setIcon( "InterSpec_resources/images/copy_small.svg" );
  m_copyBtn->setStyleClass( "LinkBtn FluxCopyBtn" );
  HelpSystem::attachToolTipOn( m_copyBtn, WString::tr("ftw-tt-copy-btn"), showToolTips );

  // TODO: "upgrade" to using the InterSpecApp 'miscSignal' directly in CopyFluxDataTextToClipboard, and get rid of m_infoCopied signal handler
  //"Wt.emit( $('.specviewer').attr('id'), {name:'miscSignal'}, 'showMsg-info-' );"

  m_copyBtn->clicked().connect( "function(s,e){ "
    "var success = Wt.WT.CopyFluxDataTextToClipboard(s,e,'" + m_copyBtn->id() + "'); "
    "Wt.emit( '" + id() + "', {name:'infocopied', eventObject:e}, success );"
  "}" );

  m_infoCopied.connect( this, [this]( const int copied ){
    tableCopiedToCliboardCallback( copied );
  } );
#endif

  WContainerWidget *buttonBox = outputOptsDiv->addNew<WContainerWidget>();
  buttonBox->addStyleClass( "FluxInfoAmount" );

  WRadioButton *simpleInfo = buttonBox->addNew<WRadioButton>( WString::tr("ftw-simple") );
  WRadioButton *standardInfo = buttonBox->addNew<WRadioButton>( WString::tr("ftw-standard") );
  WRadioButton *moreInfo = buttonBox->addNew<WRadioButton>( WString::tr("ftw-more") );

  m_displayLevelButtons = std::make_shared<WButtonGroup>();
  m_displayLevelButtons->addButton( simpleInfo, static_cast<int>(DisplayInfoLevel::Simple) );
  m_displayLevelButtons->addButton( standardInfo, static_cast<int>(DisplayInfoLevel::Normal) );
  m_displayLevelButtons->addButton( moreInfo, static_cast<int>(DisplayInfoLevel::Extended) );
  m_displayLevelButtons->setCheckedButton( m_narrowLayout ? simpleInfo : standardInfo );

  m_displayLevelButtons->checkedChanged().connect( this, [this](){
    const auto level = static_cast<DisplayInfoLevel>( m_displayLevelButtons->checkedId() );

    assert( (level == DisplayInfoLevel::Simple)
           || (level == DisplayInfoLevel::Normal)
           || (level == DisplayInfoLevel::Extended) );

    setDisplayInfoLevel( level, false );
  } );

  updateForDrf();
  setDisplayInfoLevel( m_narrowLayout ? DisplayInfoLevel::Simple : DisplayInfoLevel::Normal, true );
}//void init()


void FluxToolWidget::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  if( m_needsTableRefresh )
  {
    // Can change which columns are shown, which itself asks for a refresh - that we are doing now.
    updateBackgroundSubVisibility();
    m_needsTableRefresh = false;
    refreshPeakTable();
  }
  
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
  const std::shared_ptr<PeakModel> peakmodel = m_interspec->peakModel();
  
  m_nucNames.clear();
  m_data.clear();
  m_statUncerts.clear();
  m_systUncerts.clear();
  m_rowOutOfRange.clear();

  // Note that m_deselectedEnergies is deliberately *not* touched here - it is keyed by
  //  energy exactly so it can outlive the table.  Several of the early-returns below leave
  //  m_data empty, and re-deriving the users selections from an empty table would throw
  //  them away.
  m_rowSelected.clear();

  m_msg->setText( "" );
  m_drfUncertNote->setText( "" );
  m_calcNote->setText( "" );
  m_rangeWarning->setText( "" );

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

  const shared_ptr<const DetectorPeakResponse> det = m_detector->detector();
  if( !det || !det->isValid() )
  {
    showMessageForEmptyTable( "ftw-no-drf" );
    return;
  }

  // A fixed-geometry efficiency already describes the whole source-detector setup, so the distance
  //  isnt used - or even parsed, so an invalid one doesnt blank the table.
  const bool fixed_geom = det->isFixedGeometry();

  double distance = 1.0*PhysicalUnits::meter;
  if( !fixed_geom )
  {
    try
    {
      distance = PhysicalUnits::stringToDistance( m_distance->text().toUTF8() );
    }catch(...)
    {
      distance = -1.0;
    }

    // e.g., "1e305 km" parses, but overflows to infinity
    if( !(distance >= 0.0) || std::isinf(distance) )
    {
      showMessageForEmptyTable( "ftw-invalid-dist" );
      return;
    }
  }//if( !fixed_geom )

  // The displayed histogram (not the spectrum file), so the live time is that of the samples shown.
  const shared_ptr<const SpecUtils::Measurement> spec = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  if( !spec )
  {
    showMessageForEmptyTable( "ftw-no-foreground" );
    return;
  }

  const float live_time = spec->live_time();
  if( !(live_time > 0.0f) || std::isinf(live_time) )
  {
    showMessageForEmptyTable( "ftw-invalid-livetime" );
    return;
  }

  // Background peaks to subtract, put on the foregrounds footing - the same way the
  //  Activity/Shielding fit does it.
  const bool bkg_sub = backgroundSubtractionActive();
  vector<PeakDef> scaled_back_peaks;
  double back_scale = 0.0, back_live_time_ratio = 0.0;
  size_t num_non_gaus_back = 0;
  if( bkg_sub )
  {
    const shared_ptr<const SpecUtils::Measurement> back = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Background);
    const shared_ptr<const deque<shared_ptr<const PeakDef>>> back_peaks
                                      = peakmodel->peaks( SpecUtils::SpectrumType::Background );
    if( back && back_peaks && !back_peaks->empty() )
    {
      vector<PeakDef> peaks;
      for( const shared_ptr<const PeakDef> &p : *back_peaks )
      {
        if( p )
          peaks.push_back( *p );
      }

      const float back_live_time = back->live_time();
      if( (back_live_time > 0.0f) && !std::isinf(back_live_time) )
        back_live_time_ratio = static_cast<double>(live_time) / back_live_time;

      // The displayed scale factor (the live-time ratio, unless the user has changed it), falling back
      //  to the live-time ratio - so a background with no live time still uses the display scale.
      const double scale = GammaInteractionCalc::background_peak_scale(
                              m_interspec->displayScaleFactor(SpecUtils::SpectrumType::Background),
                              live_time, back_live_time );
      if( (scale > 0.0) && !std::isinf(scale) && !std::isnan(scale) )
      {
        back_scale = scale;
        vector<double> skipped_energies;
        scaled_back_peaks = GammaInteractionCalc::scale_background_peaks( peaks, back_scale, &skipped_energies );
        num_non_gaus_back = skipped_energies.size();
      }
    }//if( we have background peaks )
  }//if( bkg_sub )

  FluxToolCalc::Options options;
  options.distance = distance;
  options.correct_for_air = m_attenuateForAir;
  options.scaled_background_peaks = bkg_sub ? &scaled_back_peaks : nullptr;

  const vector<PeakDef> peaks = peakmodel->peakVec();
  
  const size_t npeaks = peaks.size();
  m_nucNames.resize( npeaks );
  m_data.resize( npeaks );
  m_statUncerts.resize( npeaks );
  m_systUncerts.resize( npeaks );
  m_rowOutOfRange.resize( npeaks, false );

  // The range of detector-efficiency uncertainty folded into the rows, for the note below.
  double min_eff_uncert = std::numeric_limits<double>::infinity(), max_eff_uncert = 0.0;
  bool eff_uncert_from_model_only = false;
  size_t num_back_subtracted = 0;
  vector<string> out_of_range_energies;

  for( size_t i = 0; i < npeaks; ++i )
  {
    m_data[i].fill( 0.0 );
    m_statUncerts[i].fill( 0.0 );
    m_systUncerts[i].fill( 0.0 );
    
    const PeakDef &peak = peaks[i];
    
    if( peak.parentNuclide() )
      m_nucNames[i] = peak.parentNuclide()->symbol;
    else if( peak.xrayElement() )
      m_nucNames[i] = peak.xrayElement()->symbol;
    else if( peak.reaction() )
      m_nucNames[i] = peak.reaction()->name();

    //ToDO: make sure peakArea() works for non-Gaussian peaks
    FluxToolCalc::PeakFlux flux;
    try
    {
      flux = FluxToolCalc::peak_flux( peak, *det, *spec, options );
    }catch( std::exception &e )
    {
      // The inputs were checked above, so this shouldnt happen - but this runs from render(), where an
      //  exception would end the session, so show an empty table instead.
      Wt::log("error") << "FluxToolWidget::refreshPeakTable: " << e.what();
      m_nucNames.clear();
      m_data.clear();
      m_statUncerts.clear();
      m_systUncerts.clear();
      m_rowOutOfRange.clear();
      showMessageForEmptyTable( "ftw-err-calc" );
      return;
    }//try / catch

    min_eff_uncert = std::min( min_eff_uncert, flux.eff_uncert_frac );
    max_eff_uncert = std::max( max_eff_uncert, flux.eff_uncert_frac );
    eff_uncert_from_model_only |= !flux.eff_uncert_states_own;
    num_back_subtracted += flux.background_subtracted ? 1 : 0;

    m_rowOutOfRange[i] = flux.outside_drf_energy_range;
    if( flux.outside_drf_energy_range )
    {
      char buffer[32] = { '\0' };
      snprintf( buffer, sizeof(buffer), "%.2f", flux.energy );
      out_of_range_energies.push_back( buffer );
    }
    
    m_data[i][FluxColumns::FluxEnergyCol] = flux.energy;
    m_data[i][FluxColumns::FluxPeakCpsCol] = flux.cps;
    m_statUncerts[i][FluxColumns::FluxPeakCpsCol] = flux.cps_stat;

    m_data[i][FluxColumns::FluxBackgroundCpsCol] = flux.background_cps;
    m_statUncerts[i][FluxColumns::FluxBackgroundCpsCol] = flux.background_cps_stat;
    
    m_data[i][FluxColumns::FluxGeometricEffCol] = flux.geom_eff;
    m_data[i][FluxColumns::FluxIntrinsicEffCol] = flux.intrinsic_eff;
    m_systUncerts[i][FluxColumns::FluxIntrinsicEffCol] = flux.intrinsic_syst;

    m_data[i][FluxColumns::FluxFluxOnDetCol] = flux.flux_on_det;
    m_statUncerts[i][FluxColumns::FluxFluxOnDetCol] = flux.flux_on_det_stat;
    m_systUncerts[i][FluxColumns::FluxFluxOnDetCol] = flux.flux_on_det_syst;

    m_data[i][FluxColumns::FluxFluxPerCm2PerSCol] = flux.flux_per_cm2;
    m_statUncerts[i][FluxColumns::FluxFluxPerCm2PerSCol] = flux.flux_per_cm2_stat;
    m_systUncerts[i][FluxColumns::FluxFluxPerCm2PerSCol] = flux.flux_per_cm2_syst;

    m_data[i][FluxColumns::FluxGammasInto4PiCol] = flux.gammas_4pi;
    m_statUncerts[i][FluxColumns::FluxGammasInto4PiCol] = flux.gammas_4pi_stat;
    m_systUncerts[i][FluxColumns::FluxGammasInto4PiCol] = flux.gammas_4pi_syst;
  }//for( const PeakDef &peak : peaks )

  // Say whether the uncertainties include the detector efficiency's, and what kind.
  if( npeaks > 0 )
  {
    string range = SpecUtils::printCompact( 100.0*min_eff_uncert, 2 );
    if( range != SpecUtils::printCompact( 100.0*max_eff_uncert, 2 ) )
      range += "-" + SpecUtils::printCompact( 100.0*max_eff_uncert, 2 );

    if( max_eff_uncert <= 0.0 )
      m_drfUncertNote->setText( WString::tr("ftw-drf-uncert-none") );
    else if( eff_uncert_from_model_only )
      m_drfUncertNote->setText( WString::tr("ftw-drf-uncert-model-only").arg( range ) );
    else
      m_drfUncertNote->setText( WString::tr("ftw-drf-uncert-stated").arg( range ) );
  }//if( npeaks > 0 )

  // Short notes on how the numbers were arrived at
  vector<string> notes;
  if( fixed_geom )
  {
    const string &postfix = DetectorPeakResponse::det_eff_geom_type_postfix( det->geometryType() );
    notes.push_back( WString::tr("ftw-note-fixed-geom").arg( postfix ).toUTF8() );
  }//if( fixed_geom )

  if( bkg_sub && (npeaks > 0) )
  {
    // The scale is only worth pointing out when it isnt the live-time ratio
    if( (num_back_subtracted > 0) && !(back_live_time_ratio > 0.0) )
      notes.push_back( WString::tr("ftw-note-bkg-scale-no-lt").arg( SpecUtils::printCompact(back_scale, 4) ).toUTF8() );
    else if( (num_back_subtracted > 0)
       && (fabs(back_scale - back_live_time_ratio) > 1.0E-4*back_live_time_ratio) )
      notes.push_back( WString::tr("ftw-note-bkg-scale").arg( SpecUtils::printCompact(back_scale, 4) )
                                    .arg( SpecUtils::printCompact(back_live_time_ratio, 4) ).toUTF8() );

    if( num_back_subtracted == 0 )
      notes.push_back( WString::tr("ftw-note-bkg-none").toUTF8() );

    if( num_non_gaus_back > 0 )
      notes.push_back( WString::tr("ftw-note-bkg-non-gaus").arg( static_cast<int>(num_non_gaus_back) ).toUTF8() );
  }//if( bkg_sub )

  string note_txt;
  for( const string &note : notes )
    note_txt += (note_txt.empty() ? "" : "<br />") + note;
  m_calcNote->setText( WString::fromUTF8( note_txt ) );

  if( !out_of_range_energies.empty() )
  {
    string energies;
    for( size_t i = 0; i < out_of_range_energies.size(); ++i )
      energies += (i ? ", " : "") + out_of_range_energies[i];
    m_rangeWarning->setText( WString::tr("ftw-warn-outside-range").arg( energies )
                                                         .arg( FluxToolImp::drf_energy_range_txt(det) ) );
  }//if( !out_of_range_energies.empty() )

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
  // Must stash the text on the elements own `_isData` expando: that is what
  //  `CopyFluxDataTextToClipboard` (top of this file) reads.  jQuery's `.data()` writes into
  //  jQuery's internal cache instead, which the reader never looks at.
  //  jsStringLiteral escapes the quotes/newlines the CSV payload contains.
  m_copyBtn->doJavaScript( "var el=document.getElementById('" + m_copyBtn->id() + "');"
                           " el._isData = el._isData || {};"
                           " el._isData.TableData = "
                           + Wt::WWebWidget::jsStringLiteral( pastebrdtxt.str(), '\'' ) + ";" );
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

#if( PERFORM_DEVELOPER_CHECKS )
  // For a DRF with a geometry-modeled response, the geometric efficiency uses the responses crystal
  //  (see FluxToolCalc::crystal_diameter) - note when the DRFs own stated values say otherwise.
  if( drf && drf->ceeloResponse() && !drf->isFixedGeometry() )
  {
    const double crystal_diam = FluxToolCalc::crystal_diameter( *drf );
    const double crystal_setback = FluxToolCalc::crystal_setback( *drf );
    if( (fabs( drf->detectorDiameter() - crystal_diam ) > 0.01*crystal_diam)
       || (fabs( drf->detectorSetback() - crystal_setback ) > 0.01*crystal_diam) )
    {
      const string msg = "Flux Tool: DRF '" + drf->name() + "' states a diameter of "
                         + PhysicalUnits::printToBestLengthUnits( drf->detectorDiameter(), 3 )
                         + " and setback of " + PhysicalUnits::printToBestLengthUnits( drf->detectorSetback(), 3 )
                         + ", but its geometry-modeled response has a crystal diameter of "
                         + PhysicalUnits::printToBestLengthUnits( crystal_diam, 3 ) + " at "
                         + PhysicalUnits::printToBestLengthUnits( crystal_setback, 3 )
                         + " behind the face; the response's values are used for the geometric efficiency.";
      log_developer_error( __func__, msg.c_str() );
    }
  }//if( a geometry-modeled response )
#endif

  updateForDrf();
  setDisplayInfoLevel( level, true );
  
  setTableNeedsUpdating();
}//void handleDrfChange()


void FluxToolWidget::updateForDrf()
{
  const shared_ptr<const DetectorPeakResponse> det = m_detector->detector();
  const bool fixed_geom = (det && det->isFixedGeometry());

  // A fixed-geometry efficiency describes the whole source-detector setup, so has no distance, or air
  //  between source and detector, to account for.
  m_distance->setDisabled( fixed_geom );
  m_airAttenCb->setDisabled( fixed_geom );

  updateColumnLabels();
}//void updateForDrf()


void FluxToolWidget::updateColumnLabels()
{
  const shared_ptr<const DetectorPeakResponse> det = m_detector ? m_detector->detector() : nullptr;
  const bool fixed_geom = (det && det->isFixedGeometry());
  const bool bkg_sub = backgroundSubtractionActive();

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
        if( bkg_sub )
        {
          m_colnames[col] = WString::tr("ftw-hdr-net-cps");
          m_colnamesCsv[col] = WString::fromUTF8("Background-subtracted Peak CPS");
        }else
        {
          m_colnames[col] = WString("{1} {2}").arg(WString::tr("Peak")).arg(WString::tr("CPS"));
          m_colnamesCsv[col] = WString::fromUTF8("Peak CPS");
        }
        break;
      case FluxBackgroundCpsCol:
        m_colnames[col] = WString::tr("ftw-hdr-bkg-cps");
        m_colnamesCsv[col] = WString::fromUTF8("Background Peak CPS (scaled to foreground)");
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
        if( fixed_geom )
        {
          // The sources emission rate - per cm2, m2, or gram, for those geometry types
          const string &postfix = DetectorPeakResponse::det_eff_geom_type_postfix( det->geometryType() );
          m_colnames[col] = WString::tr("ftw-hdr-fixed-geom-rate").arg( postfix );
          m_colnamesCsv[col] = WString::fromUTF8( "gammas/s" + postfix );
        }else
        {
          m_colnames[col] = WString::tr("ftw-hdr-flux-4pi");
          m_colnamesCsv[col] = WString::fromUTF8("gammas/4pi/s");
        }

        // The Simple level has no CPS column to say so
        if( bkg_sub )
          m_colnamesCsv[col] = WString::fromUTF8( "Background-subtracted " + m_colnamesCsv[col].toUTF8() );
        break;
      case FluxNumColumns:        break;
    }//switch( col )
  }//for( loop over columns )

  // Not a model reset, which under Wt4 leaves duplicate headers behind - see handleFluxToolWidgetUpdated
  if( m_fluxModel )
    m_fluxModel->headerDataChanged().emit( Wt::Orientation::Horizontal, 0, sm_selectColumn );
}//void updateColumnLabels()


bool FluxToolWidget::backgroundSubtractionActive() const
{
  return m_subtractBackground && m_backgroundPeakSub && !m_backgroundPeakSub->isHidden();
}//bool backgroundSubtractionActive() const


void FluxToolWidget::updateBackgroundSubVisibility()
{
  bool available = false;

  const shared_ptr<const SpecUtils::Measurement> back_hist
                     = m_interspec->displayedHistogram( SpecUtils::SpectrumType::Background );
  const shared_ptr<PeakModel> peakmodel = m_interspec->peakModel();
  if( back_hist && peakmodel )
  {
    const shared_ptr<const deque<shared_ptr<const PeakDef>>> back_peaks
                                    = peakmodel->peaks( SpecUtils::SpectrumType::Background );

    // Only what the automated search already found - never wait on a search here.
    const shared_ptr<const SpecMeas> back_meas = m_interspec->measurment( SpecUtils::SpectrumType::Background );
    const shared_ptr<const deque<shared_ptr<const PeakDef>>> auto_peaks = back_meas
              ? back_meas->automatedSearchPeaks( m_interspec->displayedSamples(SpecUtils::SpectrumType::Background) )
              : nullptr;

    // Whether any background peak is close enough to `fore` to be subtracted from it; the same
    //  test as GammaInteractionCalc::overlapping_background.
    const auto overlaps = []( const PeakDef &fore, const shared_ptr<const deque<shared_ptr<const PeakDef>>> &back ) -> bool {
      if( !back )
        return false;
      const double sigma = fore.gausPeak() ? fore.sigma() : 0.25*fore.roiWidth();
      for( const shared_ptr<const PeakDef> &b : *back )
      {
        // Only Gaussian background peaks are subtracted (see GammaInteractionCalc::scale_background_peaks)
        if( b && b->gausPeak()
           && (fabs(b->mean() - fore.mean()) < (GammaInteractionCalc::sm_background_peak_nsigma_near*sigma)) )
          return true;
      }
      return false;
    };//overlaps lambda

    for( const PeakDef &fore : peakmodel->peakVec() )
    {
      if( overlaps( fore, back_peaks ) || overlaps( fore, auto_peaks ) )
      {
        available = true;
        break;
      }
    }//for( const PeakDef &fore : peakmodel->peakVec() )
  }//if( back_hist && peakmodel )

  const bool was_active = backgroundSubtractionActive();
  m_backgroundPeakSub->setHidden( !available );

  if( was_active != backgroundSubtractionActive() )
  {
    updateColumnLabels();
    setDisplayInfoLevel( m_displayInfoLevel, true );  //shows/hides the background CPS column
  }
}//void updateBackgroundSubVisibility()


void FluxToolWidget::setAttenuateForAir( const bool attenuate )
{
  if( attenuate == m_attenuateForAir )
    return;

  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    auto undo_redo = [attenuate]( const bool is_undo ){
      InterSpec *viewer = InterSpec::instance();
      FluxToolWindow *fluxwin = viewer ? viewer->createFluxTool() : nullptr;
      FluxToolWidget *tool = fluxwin ? fluxwin->m_fluxTool : nullptr;
      if( tool )
        tool->applyAttenuateForAir( is_undo ? !attenuate : attenuate );
    };
    auto undo = [=](){ undo_redo(true); };
    auto redo = [=](){ undo_redo(false); };
    undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Change flux tool air attenuation." );
  }//if( undoRedo && undoRedo->canAddUndoRedoNow() )

  applyAttenuateForAir( attenuate );
}//void setAttenuateForAir( const bool attenuate )


void FluxToolWidget::applyAttenuateForAir( const bool attenuate )
{
  m_attenuateForAir = attenuate;
  if( m_airAttenCb->isChecked() != attenuate )
    m_airAttenCb->setChecked( attenuate );

  setTableNeedsUpdating();
}//void applyAttenuateForAir( const bool attenuate )


void FluxToolWidget::setSubtractBackground( const bool subtract )
{
  if( subtract == m_subtractBackground )
    return;

  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    auto undo_redo = [subtract]( const bool is_undo ){
      InterSpec *viewer = InterSpec::instance();
      FluxToolWindow *fluxwin = viewer ? viewer->createFluxTool() : nullptr;
      FluxToolWidget *tool = fluxwin ? fluxwin->m_fluxTool : nullptr;
      if( tool )
        tool->applySubtractBackground( is_undo ? !subtract : subtract );
    };
    auto undo = [=](){ undo_redo(true); };
    auto redo = [=](){ undo_redo(false); };
    undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Change flux tool background subtraction." );
  }//if( undoRedo && undoRedo->canAddUndoRedoNow() )

  applySubtractBackground( subtract );

  // Offer to fit background peaks the user hasnt, under foreground peaks - as the
  //  Activity/Shielding fit does.  Adding them is its own undo/redo step (InterSpec::setPeaks).
  //  Only when the backgrounds automated search has already finished: the prompt would otherwise
  //  wait on it, freezing the session.
  const shared_ptr<const SpecMeas> back_meas = m_interspec->measurment( SpecUtils::SpectrumType::Background );
  const bool have_auto_peaks = back_meas
      && back_meas->automatedSearchPeaks( m_interspec->displayedSamples(SpecUtils::SpectrumType::Background) );
  if( subtract && !m_skipBackgroundPeakCheck && have_auto_peaks )
  {
    const std::shared_ptr<PeakModel> peakmodel = m_interspec->peakModel();
    const shared_ptr<const deque<PeakModel::PeakShrdPtr>> fore_peaks = peakmodel ? peakmodel->peaks() : nullptr;
    if( fore_peaks && !fore_peaks->empty() )
    {
      const vector<shared_ptr<const PeakDef>> peaks( begin(*fore_peaks), end(*fore_peaks) );

      // Tied to `this`, so not called if this tool is closed before the user answers.
      const auto on_done = [this]( const PeakSearchGuiUtils::MissingBackPeaksOutcome outcome ){
        if( outcome == PeakSearchGuiUtils::MissingBackPeaksOutcome::Declined )
          m_skipBackgroundPeakCheck = true;
        setTableNeedsUpdating();
      };

      PeakSearchGuiUtils::prompt_to_add_missing_background_peaks( m_interspec, peaks, this, on_done );
    }//if( fore_peaks && !fore_peaks->empty() )
  }//if( subtract && !m_skipBackgroundPeakCheck )
}//void setSubtractBackground( const bool subtract )


void FluxToolWidget::applySubtractBackground( const bool subtract )
{
  const bool was_active = backgroundSubtractionActive();

  m_subtractBackground = subtract;
  if( m_backgroundPeakSub->isChecked() != subtract )
    m_backgroundPeakSub->setChecked( subtract );

  if( was_active != backgroundSubtractionActive() )
  {
    updateColumnLabels();
    setDisplayInfoLevel( m_displayInfoLevel, true );  //shows/hides the background CPS column
  }

  setTableNeedsUpdating();
}//void applySubtractBackground( const bool subtract )


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
  const bool bkgSub = backgroundSubtractionActive();
  
  for( FluxColumns col = FluxColumns(0); col < FluxColumns::FluxNumColumns; col = FluxColumns(col + 1) )
  {
    const bool show = FluxToolImp::showFluxColumn(m_displayInfoLevel, col, det, bkgSub);
    m_table->setColumnHidden( col, !show );
    
    WLength length;
    switch( col )
    {
      case FluxEnergyCol:         
        length = WLength( m_narrowLayout ? 4.25 : 7.5, WLength::Unit::FontEm);
        break;
        
      case FluxNuclideCol:        
        length = WLength( m_narrowLayout ? 3.25 : 5.0, WLength::Unit::FontEm);
        break;
        
      case FluxPeakCpsCol:        
        length = WLength(7.5, WLength::Unit::FontEm);
        break;

      case FluxBackgroundCpsCol:
        length = WLength(7.5, WLength::Unit::FontEm);
        break;
        
      case FluxIntrinsicEffCol:   
        length = WLength(6.5, WLength::Unit::FontEm);
        break;
        
      case FluxGeometricEffCol:   
        length = WLength(6.5, WLength::Unit::FontEm);
        break;
        
      case FluxFluxOnDetCol:      
        length = WLength(7.5, WLength::Unit::FontEm);
        break;
        
      case FluxFluxPerCm2PerSCol:
        length = WLength(9.0, WLength::Unit::FontEm);
        break;
        
      case FluxGammasInto4PiCol:  
        length = WLength( m_narrowLayout ? 7.0 : 9.0, WLength::Unit::FontEm);
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
  m_table->setColumnWidth( sm_selectColumn, WLength(2.5, WLength::Unit::FontEm) );

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
  const std::shared_ptr<PeakModel> peakmodel = m_interspec->peakModel();
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
  const std::shared_ptr<PeakModel> peakmodel = m_interspec->peakModel();
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

