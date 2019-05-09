#ifndef MakeDrf_h
#define MakeDrf_h
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

#include <Wt/WContainerWidget>

#include "InterSpec/AuxWindow.h"

class MakeDrf;
class InterSpec;
class MaterialDB;
class MakeDrfChart;

namespace Wt
{
  class WText;
  class WLineEdit;
  class WSuggestionPopup;
}

class MakeDrfWindow : public AuxWindow
{
public:
  MakeDrfWindow( InterSpec *viewer, MaterialDB *materialDB, Wt::WSuggestionPopup *materialSuggest );
  
  virtual ~MakeDrfWindow();
  
protected:
  MakeDrf *m_makeDrf;
};//class MakeDrfWindow


class MakeDrf : public Wt::WContainerWidget
{
public:
  MakeDrf( InterSpec *viewer,
           MaterialDB *materialDB,
           Wt::WSuggestionPopup *materialSuggest,
           Wt::WContainerWidget *parent = nullptr );
  
  virtual ~MakeDrf();
  
protected:
  void handleSourcesUpdates();
  
  InterSpec *m_interspec;
  MaterialDB *m_materialDB;
  Wt::WSuggestionPopup *m_materialSuggest;
  
  MakeDrfChart *m_chart;
  
  Wt::WContainerWidget *m_files;
  
  Wt::WLineEdit *m_detDiameter;
  
  Wt::WText *m_errorMsg;
};//class MakeDrf

#endif //MakeDrf_h
