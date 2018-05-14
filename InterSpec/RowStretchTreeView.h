#ifndef RowStretchTreeView_h
#define RowStretchTreeView_h
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

#include <map>

#include <Wt/WTreeView>
//#include <Wt/WItemDelegate>


// This class stretchs the column widths to be as cumulatively as wide as the
// table is
// Reference: http://redmine.emweb.be/boards/2/topics/2659

class RowStretchTreeView : public Wt::WTreeView
{
  /**
   Note: DO NOT RESIZE INSTANCES OF THIS CLASS.  YOU WILL GET WT RESIZE ERRORS.
   */

public:
  RowStretchTreeView();
  virtual ~RowStretchTreeView(){}
  
  //Setting the width of a column will ensure it is never narrower
  //  than that, but the column will be stretched by an equal factor
  //  as all the other columns with set widths in order to make the
  //  table fill out 100% of its allocated width.
  virtual void setColumnWidth( int column, const Wt::WLength &width );
  
  //Hi-jack the model calls so we can intercept when rows are being added or
  //  removed, in order to deal with if the scroll bars will appear or disapear.
  virtual void setModel( Wt::WAbstractItemModel *model );
  
protected:
  //render is specialized to check to see if the row widths have changed, and
  //  thus column widths should be adjusted; this check happens both when the
  //  outer width of this widget changes, as well as when ever it is rendered.
  virtual void render(	Wt::WFlags<Wt::RenderFlag> flags );
  
  //widthChanged(): called by m_rowWidthChanged whenever the clientside detects
  //  the inner width of the table area has changed.
  void widthChanged( int widthpx, int scrollbarWidth );
  
  //layoutSizeChanged(): Adjusts column widths, only if m_rowWidthPx is not
  //  valid (essentially a falback incase determining the inner width fails
  //  due to table heirarchy not being as expected).  Will leave 15px plus 1px
  //  for each column with a set width, as spacing on the right hand side to
  //  account for the scrollbar.
  virtual void layoutSizeChanged( int width, int height );
  
  //Called when rows are added to the model.
  void rowsAddedCallback( Wt::WModelIndex index, int first, int last );
  
  //Called when rows are removed from the model.
  void rowsRemovedCallback( Wt::WModelIndex index, int first, int last );
  
protected:
  //m_nominalWidth: width of the columns; if they add up to more that rendered
  //  table width, columns will be as sepcified; if hte table is wider than this
  //  then rows will be scaled relative to these widths to fill up the entire
  //  table.
  std::map<int,double> m_nominalWidth;
  
  //m_rowWidthPx: the width of the rows in the table (doesnt include vertical
  //  scroll bar for this area).  Will be -1 until the table has been rendered
  //  and server notified of the width.
  int m_rowWidthPx;
  
  //m_hasScrollBar: set if the javascript TreeViewCheckWidth function
  //  detected scroll bars are present last time the internal width of the table
  //  changed.
  bool m_hasScrollBar;
  
  //m_scrollBarWidth: The detected scrollbar width in px; will be -1 until the
  //  scroll bars are detected.
  int m_scrollBarWidth;
  
  //m_width: width of the total widget, as set in layoutSizeChanged(...)
  int m_width;
  
  //m_height: hight of the total widget, as set in layoutSizeChanged(...)
  int m_height;
  
  //m_rowWidthChanged: signal to notify server-side that the width of the area
  //  holding the table rows has changed (for instance if scroll bars existence
  //  changed, or table resized).
  Wt::JSignal<int,int> m_rowWidthChanged;
};//class RowStretchTreeView

//class HeaderDelegate: public Wt::WItemDelegate
//{
//public:
//  HeaderDelegate( RowStretchTreeView* tv ) : m_TableView( tv ) {}
//  
//  virtual Wt::WWidget *update( Wt::WWidget *widget, const Wt::WModelIndex& index, Wt::WFlags< Wt::ViewItemRenderFlag > flags )
//  {
//    //    Wt::WInteractWidget *iw = dynamic_cast< Wt::WInteractWidget* >( widget );
//    if (!widget)
//    {
//      widget = Wt::WItemDelegate::update( widget, index, flags );
//      Wt::WInteractWidget *iw = dynamic_cast< Wt::WInteractWidget* >( widget );
//      //    iw->mouseWentUp().connect( this, &HeaderDelegate::ShowHeaderPopup ); // <-- access violation exception here
////      widget->setToolTip("test");
//    }
//    return widget;
//  }
//  
//  
//private:
//  RowStretchTreeView* m_TableView;
//}; //class HeaderDelegate

#endif  //RowStretchTreeView_h
