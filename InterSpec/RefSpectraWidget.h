#ifndef RefSpectraWidget_h
#define RefSpectraWidget_h

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

#include <memory>
#include <string>
#include <vector>

#include <Wt/WContainerWidget>

#include "InterSpec/SimpleDialog.h"

namespace Wt 
{
  class WComboBox;
  class WTreeView;
  class WCheckBox;
  class WPushButton;
  class WStackedWidget;
}

namespace SpecUtils{ enum class SpectrumType : int; }

class RefSpectraModel;
class RefSpectraWidget;
class D3SpectrumDisplayDiv;

/** An enum to describe the initial behaviour of the RefSpectraWidget.
*/
enum class RefSpectraInitialBehaviour
{
  /** All ref spectra directories will be shown collapsed, for the user to expand. */
  Default,

  /** The tool will attempt to guess a built-in reference spectrum that most closely
      matches the current foreground spectrum, and open with that spectrum selected. 

      Not implemented yet.
  */
  GuessMatchToForeground
};

class RefSpectraDialog : public SimpleDialog
{
public:
  RefSpectraDialog( const Wt::WString &title = "Reference Spectra" );
  virtual ~RefSpectraDialog();

  RefSpectraWidget *widget() { return m_widget; }

  /** Create a new RefSpectraDialog with the given initial behaviour and spectrum type.
   * 
   *  @param initialBehaviour How to configure the initial view of the tool.
   *  @param type The default spectrum type to load selected file as.
   *  @return A new RefSpectraDialog.
   */
  static RefSpectraDialog *createDialog( const RefSpectraInitialBehaviour initialBehaviour, 
                                         const SpecUtils::SpectrumType type );
private:
  RefSpectraWidget *m_widget;
};//class RefSpectraDialog



class RefSpectraWidget : public Wt::WContainerWidget
{
public:

  RefSpectraWidget( Wt::WContainerWidget *parent = nullptr );
  virtual ~RefSpectraWidget();

  void setLoadSpectrumType( SpecUtils::SpectrumType type );
  void selectSimilarToForeground();

protected:
  void handleSelectionChanged();
#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
  void startAddDirectory();
  void removeCurrentlySelectedDirectory();
  void addDirectory( const std::string &path );
#endif // !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT
  void updatePreview();
  void updateTreeView();
  
private:
  RefSpectraModel *m_treeModel;
  Wt::WTreeView *m_treeView;
  Wt::WStackedWidget *m_stack;
  D3SpectrumDisplayDiv *m_spectrum;
  Wt::WContainerWidget *m_dirInfoContainer;
  Wt::WCheckBox *m_showCurrentForeground;
  Wt::WCheckBox *m_refBackground;
#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
  Wt::WContainerWidget *m_addDirButton;
  Wt::WContainerWidget *m_deleteDirButton;
  Wt::WPushButton *m_showInExplorerButton;
#endif // !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT
  Wt::WComboBox *m_showAsComboBox;
  std::vector<std::string> m_baseDirectories;
  
  void setupUI();
  void createTreeModel();
  void initBaseDirs();
  void showInExplorer();
};

#endif // RefSpectraWidget_h 