#ifndef GroupBox_h
#define GroupBox_h
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

#include <Wt/WString.h>
#include <Wt/WContainerWidget.h>


/** A drop-in replacement for WGroupBox that renders as a `<fieldset>` with
 *  a `<legend>`, but does NOT create an internal layout.
 *
 *  In Wt 4.12, WGroupBox creates an internal WVBoxLayout which prevents children
 *  added via addNew()/addWidget() from rendering (they go into children_ but
 *  not into the layout). This class avoids that problem by not setting any
 *  internal layout, while still producing the same HTML structure and
 *  participating in the Wt JS layout system for proper size computation.
 *
 *  Note: do NOT call setLayout() on this class. The legend is a regular child
 *  widget, and a layout would prevent it from rendering. If you need a layout
 *  inside a fieldset, use Wt::WGroupBox with setLayout() instead.
 */
class GroupBox : public Wt::WContainerWidget
{
public:
  GroupBox();
  explicit GroupBox( const Wt::WString &title );

  /** Returns the fieldset DOM element type, same as WGroupBox. */
  virtual Wt::DomElementType domElementType() const override;

  void setTitle( const Wt::WString &title );
  const Wt::WString &title() const;

private:
  void init();
};//class GroupBox

#endif //GroupBox_h
