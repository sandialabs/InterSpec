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

#include <Wt/WTreeView>
#include <Wt/WApplication>
#include <Wt/WAbstractItemModel>

#include "InterSpec/RowStretchTreeView.h"

using namespace Wt;
using namespace std;

WT_DECLARE_WT_MEMBER(TreeViewCheckWidth, Wt::JavaScriptFunction, "TreeViewCheckWidth",
  function(el){
    var p = $(el), w = -1, sw = -1;
    if( p.data('noadj') )
      return;
    
    w = p.children().children('.cwidth').not('.Wt-header').children().width();
    p.children().children('.cwidth').not('.Wt-header').each( function(i, val){
      var t = $("<div/>");
      $(val).append(t);
      //get width not including scroll bar, in pxs.  We could maybye just get the child divs width, but appending to be safe to force the recompuation of width...
      sw = Math.max( sw, val.offsetWidth - val.clientWidth );  //get scrollbar width
      w = w > 0 ? Math.min(w,t.width()) : t.width();
      t.remove();
    } );
    var oldw = p.data('RSTVW');
    if( (!oldw || oldw !== w) && w !== 0 ){
      //console.log( 'Will emit w=' + w + ', sw=' + sw );
      p.data('RSTVW',w);
      Wt.emit( el.id, {name: 'widthchanged'}, Math.round(w), Math.round(sw) );
    }
  }
);


RowStretchTreeView::RowStretchTreeView( WContainerWidget *parent )
  : WTreeView( parent ),
    m_rowWidthPx( -1 ),
    m_hasScrollBar( true ),
    m_scrollBarWidth( -1 ),
    m_width( -1 ),
    m_height( -1 ),
    m_rowWidthChanged( this, "widthchanged", true )
{
#if( (WT_VERSION < 0x3030100 || WT_VERSION > 0x3030400) && (WT_VERSION != 0x3070100) )
#warning The RowStretchTreeView JavaScript has only been verified for Wt 3.3.1 through 3.3.4
#endif
  
  setLayoutSizeAware( true );
  m_rowWidthChanged.connect( this, &RowStretchTreeView::widthChanged );
  
  //Could handle users effecting column sizes...
  //Signal<int,WLength> & WTreeView::columnResized();
  
  LOAD_JAVASCRIPT(wApp, "src/RowStretchTreeView.cpp", "RowStretchTreeView", wtjsTreeViewCheckWidth);
  

  //Add a little square in the upper right hand corner, that is the same color
  //  as the header, so if scroll bars appear, that little area doesnt look
  //  weird.  Definitely a bit of a hack, but whatever for now.
  WStringStream backjs;
  backjs << "$(" << jsRef() << ")"
            ".prepend(\"<div style='position: absolute; width: 20px;"
            " height: " << headerHeight().toPixels() << "px; top: 0px;"
            " right: 0px; background-color: #EEEEEE;'></div>\");";
  doJavaScript( backjs.str() );
  //doJavaScript( "$(" + jsRef() + ").css('background-color', '#EEEEEE');" );
  //doJavaScript( "$(" + jsRef() + ").find('.Wt-header.headerrh.cwidth').each( function(i,o){ console.log(o); $(o).css('background-color', '#EEEEEE');} );" );
  
  //  this->setHeaderItemDelegate( new HeaderDelegate( this ) );
};//RowStretchTreeView constructor


void RowStretchTreeView::setColumnWidth( int column, const Wt::WLength &width )
{
  WTreeView::setColumnWidth( column, width );
  
  if( width.isAuto() )
    m_nominalWidth.erase( column );
  else
    m_nominalWidth[column] = width.toPixels();
}//void setColumnWidth( int column, const Wt::WLength &width )


void RowStretchTreeView::rowsAddedCallback( Wt::WModelIndex index, int first, int last )
{
  if( m_scrollBarWidth <= 0 || m_hasScrollBar )
    return;

  //4 on next line is empirically found
  const double nrowsfit = (m_height - headerHeight().toPixels() - 4) / rowHeight().toPixels();

  if( model()->rowCount() < nrowsfit )
    return;
  
  //If we made it here, we will add scroll bars
  const int newwidth = m_rowWidthPx - m_scrollBarWidth;
  doJavaScript( "$(" + jsRef() + ").data('noadj',true);" );
  
  widthChanged( newwidth, m_scrollBarWidth );
  
  WStringStream js;
  js << "$(" << jsRef() << ").data('RSTVW'," << newwidth << ");"
     << "setTimeout( function(){$(" << jsRef() << ").data('noadj',false);"
     << "Wt.WT.TreeViewCheckWidth(" << jsRef() << ");},0);";
  doJavaScript( js.str() );
}//void rowsAddedCallback( Wt::WModelIndex index, int first, int last )


void RowStretchTreeView::rowsRemovedCallback( Wt::WModelIndex index, int first, int last )
{
  if( m_scrollBarWidth <= 0.0 || !m_hasScrollBar )
    return;
  
  //4 on next line is empirically found
  const double nrowsfit = (m_height - headerHeight().toPixels() - 4) / rowHeight().toPixels();
  
  if( model()->rowCount() <= nrowsfit )
  {
    const int newwidth = m_rowWidthPx + m_scrollBarWidth;
    doJavaScript( "$(" + jsRef() + ").data('noadj',true);" );
    widthChanged( newwidth, -1 );
    WStringStream js;
    js << "$(" << jsRef() << ").data('RSTVW'," << newwidth << ");"
       << "setTimeout( function(){$(" << jsRef() << ").data('noadj',false);},0);";
    doJavaScript( js.str() );
  }
}//void rowsRemovedCallback( Wt::WModelIndex index, int first, int last )



void RowStretchTreeView::setModel( Wt::WAbstractItemModel *model )
{
  WTreeView::setModel( model );
  
  if( model )
  {
    model->rowsInserted().connect( this, &RowStretchTreeView::rowsAddedCallback );
    model->rowsRemoved().connect(  this, &RowStretchTreeView::rowsRemovedCallback );
  }//if( model )
  
}//void setModel( Wt::WAbstractItemModel *model )


void RowStretchTreeView::refreshColWidthLayout()
{
  //doJavaScript( "Wt.WT.TreeViewCheckWidth(" + jsRef() + ");" );
  widthChanged( m_rowWidthPx, m_scrollBarWidth );
}//virtual void refreshColWidthLayout();


void RowStretchTreeView::render(	Wt::WFlags<Wt::RenderFlag> flags )
{
  WTreeView::render( flags );
  
  if( flags & RenderFull )
  {
    //Redefine the wtResize so this way adjusting column sizes saves one client
    //  to server round trip.
    //  Note that this JS is only checked to be valid Wt 3.3.1 through 3.3.4,
    //  and only very superficially tested for Wt 3.7.1
    //  The WT_RESIZE_JS member seems to be called even if a CSS flex layout
    //  is controlling the size of the item.
    setJavaScriptMember(WT_RESIZE_JS,
                        "function(self,w,h) {"
#if( WT_VERSION < 0x3030400 )
                        "$(self).data('obj').wtResize();"
#else
                        "self.wtObj.wtResize();"
#endif
                        "Wt.WT.TreeViewCheckWidth(" + jsRef() + ");"
                        "}");
  }//if( flags & RenderFull )
  
  //If we added enough entries such that a scroll bar will now apear, we will
  //  check for and fix this.
  doJavaScript( "Wt.WT.TreeViewCheckWidth(" + jsRef() + ");" );
}//void render(	Wt::WFlags<Wt::RenderFlag> flags )


void RowStretchTreeView::widthChanged( int widthpx, int scrollbarWidth )
{
  m_hasScrollBar = (scrollbarWidth > 0);

  if( m_hasScrollBar && m_scrollBarWidth >= 0 && m_scrollBarWidth != scrollbarWidth )
  {
    cerr << "\n\n\n\nGot a new scroll bar width!!!!\n\n\n\n" << endl;
#if(PERFORM_DEVELOPER_CHECKS)
    char errormsg[1024];
    snprintf( errormsg, sizeof(errormsg),
              "Got a new scroll bar width: %i when I had %i",
              int(m_scrollBarWidth), int(scrollbarWidth) );
    log_developer_error( __func__, errormsg );
#endif
  }
  
  if( m_hasScrollBar && m_scrollBarWidth < 0 )
    m_scrollBarWidth = scrollbarWidth;
  
  m_rowWidthPx = widthpx;
  
  double setWidthTotal = 0.0, nonSetWidthsTotal = 0.0;;
  int nVisibleSetWidth = 0;
  const int ncolumn = model()->columnCount();
  
  typedef map<int,double> RowWidthMap;
  
  for( const RowWidthMap::value_type &vt : m_nominalWidth )
  {
    if( !isColumnHidden(vt.first) )
    {
      ++nVisibleSetWidth;
      setWidthTotal += vt.second;
    }
  }

  for( int column = 0; column < ncolumn; ++column )
  {
    if( !isColumnHidden(column) && !m_nominalWidth.count(column) )
    {
      const WLength w = columnWidth(column);
      if( w.isAuto() )
        nonSetWidthsTotal += 157.0;
      else
        nonSetWidthsTotal += w.toPixels() + 7.0;
    }//if( !isColumnHidden(column) && !m_nominalWidth.count(column) )
  }//for( int column = 0; column < ncolumn; ++column )
  
  const double remaniningWidth = m_rowWidthPx - nonSetWidthsTotal;
  const double multiple = (remaniningWidth - 7.0*nVisibleSetWidth) / setWidthTotal;
  
  const string containerjs = "$(" + jsRef() + ").children().children('.cwidth').not('.Wt-header')";
  
  if( multiple <= 1.0 )
  {
    for( const RowWidthMap::value_type &vt : m_nominalWidth )
    {
      if( !isColumnHidden(vt.first) )
        WTreeView::setColumnWidth( vt.first, vt.second );
    }
    
    doJavaScript( containerjs + ".css('overflow-x','auto');" );
    
    return;
  }//if( multiple <= 1.0 )
  
  
  for( const RowWidthMap::value_type &vt : m_nominalWidth )
  {
    const int col = vt.first;
    const double w = multiple*vt.second;
    const double delta = columnWidth(col).toPixels() - w;
    
    //If the column is within 1 px and is to skinny, dont adjust it (possible
    //  savings on rendering?); but if its too wide, bring it down to size.
    //  (this is probably a pointless optimiztion that doesnt help - consider
    //  removing)
    if( !isColumnHidden(col) && (delta < -1.0 || delta > 0.0) )
      WTreeView::setColumnWidth( col, WLength(w,WLength::Pixel) );
  }//for( const RowWidthMap::value_type &vt : m_nominalWidth )
  
  doJavaScript( containerjs + ".css('overflow-x','hidden');" );
}//void RowStretchTreeView::widthChanged( double widthpx )


void RowStretchTreeView::layoutSizeChanged( int width, int height )
{
  m_width = width;
  m_height = height;
  
  if( m_rowWidthPx > 0 || width < 0 )
    return;
  
  //Below here is all a fallback incase the JavaScript TreeViewCheckWidth
  //  function for some reason doesnt work okay (possibly newer/oler versions of
  //  Wt, or use of Ajax, or something happening with WLayouts that I dont know
  //  about)
  double nonSetWidths = 0.0;
  const int ncolumn = model()->columnCount();
  
  for( int column = 0; column < ncolumn; ++column )
  {
    if( !isColumnHidden(column) && !m_nominalWidth.count(column) )
    {
      const WLength w = columnWidth(column);
      if( w.isAuto() )
        nonSetWidths += 157.0;
      else
        nonSetWidths += w.toPixels() + 7.0;
    }//if( !isColumnHidden(column) && !m_nominalWidth.count(column) )
  }//for( int column = 0; column < ncolumn; ++column )
  
  int nvisibleSetWidth = 0;
  double setWidthTotal = 0.0;
  for( std::map<int,double>::const_iterator iter = m_nominalWidth.begin();
      iter != m_nominalWidth.end(); ++iter )
  {
    if( !isColumnHidden(iter->first) )
    {
      ++nvisibleSetWidth;
      setWidthTotal += iter->second;
    }
  }
  
  const double remaniningWidth = width - nonSetWidths;
  
  //Giving 8px of padding instead of the 7px Wt actually does, as well as
  //  additional scrollbar room, all in order to keep there from being
  //  a scroll bar for the 'x' direction unecassarily.
  //Note: this is kinda a hack - we could figure out if there is a scroll bar
  //      and then account for that, but I'm tired tonight.
  
  int scrollBW = (m_scrollBarWidth > 0 ? m_scrollBarWidth : 15);
  if( !m_hasScrollBar && m_scrollBarWidth > 0 )
    scrollBW = 0;
  
  const double multiple = (remaniningWidth - 8.0*nvisibleSetWidth - scrollBW) / setWidthTotal;
  
  if( multiple <= 1.0 )
  {
    for( std::map<int,double>::const_iterator iter = m_nominalWidth.begin();
        iter != m_nominalWidth.end(); ++iter )
    {
      const int col = iter->first;
      const double w = iter->second;
      if( !isColumnHidden(col) && fabs(columnWidth(col).toPixels() - w) > 0.9 )
        WTreeView::setColumnWidth( col, WLength(w,WLength::Pixel));
    }
    return;
  }//if( remaniningWidth <= 0.0 )
  
  
  for( std::map<int,double>::const_iterator iter = m_nominalWidth.begin();
      iter != m_nominalWidth.end(); ++iter )
  {
    const int col = iter->first;
    const double w = multiple*iter->second;
    
    //we could probably do something a little smarter in how we decide if we
    //  need to adjust the width of the column, in terms of optimizing
    //  renderings, bandwidth, and display
    if( !isColumnHidden(col)  && fabs(columnWidth(col).toPixels() - w) > 0.9 )
      WTreeView::setColumnWidth( col, WLength(w,WLength::Pixel) );
  }
}//void layoutSizeChanged (const int width, const int height)


void RowStretchTreeView::selectRange( const Wt::WModelIndex &first, const Wt::WModelIndex &last )
{
  // Override of Wt::WTreeView::selectRange to fix an infinite loop bug in the original Wt 3.7.1
  // implementation. The original algorithm has no termination condition when backtracking past
  // the root of the tree - when the selection range ends at the last visible item, it would
  // wrap around to the first item and loop forever.
  //
  // We use WAbstractItemView::select() instead of internalSelect() since internalSelect is not
  // accessible. We block the selectionChanged() signal during batch selection.

  struct SignalBlocker
  {
    Wt::Signal<> &m_signal;
    const bool m_was_blocked;
    SignalBlocker( Wt::Signal<> &ref )
      : m_signal( ref ), m_was_blocked( ref.isBlocked() ){ m_signal.setBlocked( true ); }
    ~SignalBlocker(){ m_signal.setBlocked( m_was_blocked ); }
  };//struct SignalBlocker

  SignalBlocker blocker( selectionChanged() );

  // Helper to promote a child index to its parent if the parent is collapsed.
  // This ensures we use visible row positions for range endpoints when children
  // of a collapsed node are selected (e.g., user selects samples, then collapses file).
  auto promoteToVisibleAncestor = [this]( const WModelIndex &idx ) -> WModelIndex {
    if( !idx.isValid() )
      return idx;

    WModelIndex current = idx;

    // Walk up the parent chain, checking if any ancestor is collapsed
    while( current.parent().isValid() )
    {
      const WModelIndex parentCol0 = model()->index( current.parent().row(), 0, current.parent().parent() );
      if( !isExpanded( parentCol0 ) )
        current = current.parent(); // Parent is collapsed - promote to parent level
      else
        break; // Parent is expanded - this level is visible
    }//

    return current;
  };//promoteToVisibleAncestor lambda

  // Promote endpoints to visible ancestors if they're children of collapsed nodes
  const WModelIndex visibleFirst = promoteToVisibleAncestor( first );
  const WModelIndex visibleLast = promoteToVisibleAncestor( last );

  // Helper to compare indices in tree traversal order.
  // Returns true if 'a' comes before 'b' in depth-first traversal.
  auto indexComesBefore = []( const WModelIndex &a, const WModelIndex &b ) -> bool {
    if( !a.isValid() || !b.isValid() )
      return a.isValid(); // valid comes before invalid

    // Build ancestor chains (from index up to root)
    std::vector<WModelIndex> chainA, chainB;
    for( WModelIndex idx = a; idx.isValid(); idx = idx.parent() )
      chainA.push_back( idx );
    for( WModelIndex idx = b; idx.isValid(); idx = idx.parent() )
      chainB.push_back( idx );

    // Compare from root down (reverse order since we built bottom-up)
    const size_t minDepth = std::min( chainA.size(), chainB.size() );
    for( size_t i = 0; i < minDepth; ++i )
    {
      const WModelIndex &idxA = chainA[chainA.size() - 1 - i];
      const WModelIndex &idxB = chainB[chainB.size() - 1 - i];

      if( idxA.row() != idxB.row() )
        return idxA.row() < idxB.row();
      // Same row at this level - if same parent, compare columns
      if( idxA.column() != idxB.column() )
        return idxA.column() < idxB.column();
    }

    // One is ancestor of the other - the shorter chain (ancestor) comes first
    return chainA.size() < chainB.size();
  };

  // Ensure we traverse from earlier to later in tree order (using visible endpoints)
  const WModelIndex &startIdx = indexComesBefore( visibleFirst, visibleLast ) ? visibleFirst : visibleLast;
  const WModelIndex &endIdx = indexComesBefore( visibleFirst, visibleLast ) ? visibleLast : visibleFirst;

  // Determine column range - handle case where first and last columns may be in either order
  const int minCol = std::min( visibleFirst.column(), visibleLast.column() );
  const int maxCol = std::max( visibleFirst.column(), visibleLast.column() );

  // Safety limit to prevent infinite loops from model inconsistencies
  const int maxIterations = 10000;
  int iterations = 0;

  WModelIndex index = startIdx;

  // Helper to check if two indices refer to the same tree position
  // (maybe a bit overkill over operator==, )
  auto sameTreePosition = []( const WModelIndex &a, const WModelIndex &b ) -> bool {
    if( a.row() != b.row() )
      return false;

    // Compare full parent chains for arbitrary tree depth
    WModelIndex pa = a.parent();
    WModelIndex pb = b.parent();

    while( pa.isValid() && pb.isValid() )
    {
      if( pa.row() != pb.row() )
        return false;
      pa = pa.parent();
      pb = pb.parent();
    }

    // Both should be invalid (reached root) at the same time for same position
    return !pa.isValid() && !pb.isValid();
  };

  while( iterations++ < maxIterations )
  {
    // Check if we've reached the target row before selecting
    const bool atLastRow = sameTreePosition( index, endIdx );

    // Select all columns in current row
    for( int c = minCol; c <= maxCol; ++c )
    {
      const WModelIndex ic = model()->index( index.row(), c, index.parent() );
      if( !ic.isValid() )
        continue;

      WAbstractItemView::select( ic, Select );

      // Check both exact match and position match (in case indices differ in representation)
      if( ic == endIdx || sameTreePosition( ic, endIdx ) )
        return;
    }//for( int c = minCol; c <= maxCol; ++c )

    // If we were at the last row but didn't find exact match, we're done
    if( atLastRow )
      return;

    // Additional safety: check if we're at or past the target row (ignoring column)
    // This catches cases where position comparison fails due to index representation differences
    {
      const bool indexParentValid = index.parent().isValid();
      const bool endParentValid = endIdx.parent().isValid();

      if( indexParentValid == endParentValid )
      {
        if( !indexParentValid )
        {
          // Both root level - compare rows directly
          if( index.row() >= endIdx.row() )
            return;
        }else
        {
          // Both children - compare parent rows, then child rows
          const int indexParentRow = index.parent().row();
          const int endParentRow = endIdx.parent().row();
          if( indexParentRow > endParentRow )
            return;
          if( (indexParentRow == endParentRow) && (index.row() >= endIdx.row()) )
            return;
        }//if( !indexParentValid ) / else
      }else if( indexParentValid && !endParentValid )
      {
        // index is a child, endIdx is root - index comes after if parent row >= endIdx row
        if( index.parent().row() >= endIdx.row() )
          return;
      }
      // If index is root and endIdx is child, index always comes before (no action needed)
    }

    // Try to descend into expanded children
    const WModelIndex indexc0 = (index.column() == 0)
                                ? index
                                : model()->index( index.row(), 0, index.parent() );

    if( isExpanded( indexc0 ) && model()->hasChildren( indexc0 ) )
    {
      index = model()->index( 0, minCol, indexc0 );
      if( index.isValid() )
        continue;
    }

    // Find next row: try siblings, then backtrack to parent's sibling
    bool foundNext = false;
    WModelIndex current = index;

    while( !foundNext )
    {
      const WModelIndex parent = current.parent();
      const int nextRow = current.row() + 1;

      if( nextRow < model()->rowCount( parent ) )
      {
        index = model()->index( nextRow, minCol, parent ); // Found a next sibling
        foundNext = true;
      }else if( parent.isValid() )
      {
        current = parent; // Backtrack to parent level
      }
      else
      {
        // Reached root with no more items - we're done
        return;
      }
    }//while( !foundNext )
  }//while( iterations++ < maxIterations )

  // If we hit max iterations, log a warning (developer check)
#if( PERFORM_DEVELOPER_CHECKS )
  log_developer_error( __func__, "selectRange hit max iteration limit - possible model issue" );
#endif
}//void selectRange( const Wt::WModelIndex &first, const Wt::WModelIndex &last )
