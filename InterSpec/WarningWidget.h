#ifndef WarningWidget_h
#define WarningWidget_h
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

#include <Wt/WText>
#include <Wt/WContainerWidget>
#include <Wt/WGridLayout>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#if ( USE_SPECTRUM_CHART_D3 )
#include "InterSpec/D3SpectrumDisplayDiv.h"
#else
#include "InterSpec/SpectrumDisplayDiv.h"
#endif
#include "InterSpec/RowStretchTreeView.h"

class InterSpec;
#if ( USE_SPECTRUM_CHART_D3 )
class D3SpectrumDisplayDiv;
#else
class SpectrumDisplayDiv;
#endif

namespace Wt
{
  class WApplication;
  class WStandardItemModel;
}//namespace Wt

class WarningWidget : public Wt::WContainerWidget
{
public:
  enum WarningMsgLevel
  {
    WarningMsgInfo = 0, WarningMsgLow = 1, WarningMsgMedium = 2, WarningMsgHigh = 3 /*Not actually used anywhere*/, WarningMsgSave = 4
  };//enum WarningMsgLevel
  
  static const char *tostr( const WarningMsgLevel level );
  static const char *popupToStr( const WarningMsgLevel level );
  
  static const char *description( const WarningMsgLevel level );
  
public:
  WarningWidget(
#if ( USE_SPECTRUM_CHART_D3 )
                 D3SpectrumDisplayDiv *spectrumDisplayDiv,
#else
                 SpectrumDisplayDiv *spectrumDisplayDiv,
#endif
                 InterSpec *hostViewer,
                 Wt::WContainerWidget *parent = 0 );
  virtual ~WarningWidget();

  void createContent();
  void addMessage( const Wt::WString &msg, const Wt::WString &src, int level );

  
  void setPopupActivity( WarningMsgLevel priority, bool allowed );
  void setActivity( WarningMsgLevel priority, bool allowed );

  bool active( WarningMsgLevel level ) const;

  void resultSelectionChanged();
    
  void clearMessages();
  
protected:
#if ( USE_SPECTRUM_CHART_D3 )
  D3SpectrumDisplayDiv *m_spectrumDisplayDiv;
#else
  SpectrumDisplayDiv *m_spectrumDisplayDiv;
#endif
  
  InterSpec *m_hostViewer;

  int m_totalMessages;

  bool m_popupActive[int(WarningMsgSave)+1];
  bool m_active[int(WarningMsgSave)+1];

  Wt::WApplication *m_app;
  Wt::WGridLayout* m_layout;
  Wt::WStandardItemModel* m_messageModel;
    
  RowStretchTreeView *m_tableView;
  Wt::WText * m_description;
    
}; // class WarningWidget

#endif
