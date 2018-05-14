#ifndef DbStateBrowser_h
#define DbStateBrowser_h
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

#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpecUser.h"

namespace Wt
{
  class WTreeView;
  class WPushButton;
  class WTreeView;
}//namespace Wt

class DbStateBrowser : public AuxWindow
{
protected:
  
  std::shared_ptr<DataBaseUtils::DbSession> m_session;
  InterSpec  *m_viewer;
  Wt::Dbo::QueryModel< Wt::Dbo::ptr<UserState> >    *m_model;
  Wt::WTreeView      *m_table;
  Wt::WPushButton     *m_loadButton;
public:
  DbStateBrowser( InterSpec *viewer, bool testStatesOnly );
  virtual ~DbStateBrowser();
  void selectionChanged();
  void handleDoubleClicked( Wt::WModelIndex index,Wt::WMouseEvent event );
  void loadSelected();
};//class DbStateBrowser

#endif //DbStateBrowser_h
