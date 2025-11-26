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

#include <chrono>
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
class RowStretchTreeView;
class D3SpectrumDisplayDiv;

/** An enum to describe the initial behaviour of the RefSpectraWidget.
*/
enum class RefSpectraInitialBehaviour
{
  /** All ref spectra directories will be shown collapsed, for the user to expand. */
  Default,

  LastUserSelectedSpectra,

  /** The tool will attempt to guess a built-in reference spectrum that most closely
      matches the current foreground spectrum, and open with that spectrum selected. 

      Not implemented yet.
  */
  //GuessMatchToForeground
};

enum class RefSpectraWidgetSelectionType
{
    None, Directory, File
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

protected:
  void handleSelectionChanged( RefSpectraWidgetSelectionType type );

private:
  RefSpectraWidget *m_widget;
  Wt::WPushButton *m_loadBtn;
};//class RefSpectraDialog



class RefSpectraWidget : public Wt::WContainerWidget
{
public:
  RefSpectraWidget( Wt::WContainerWidget *parent = nullptr );
  virtual ~RefSpectraWidget();

  void setLoadSpectrumType( SpecUtils::SpectrumType type );
  void loadSelectedSpectrum();

  /** Select the last selected path from the user preference name "RefSpectraLastSelection". */
  void selectLastSelectedPath();

  /** Select a spectrum similar to the current foreground spectrum. 
   * 
   *  Not implemented yet.
  */
  //void selectSimilarToForeground();

  Wt::Signal<RefSpectraWidgetSelectionType> &fileSelectionChangedSignal();
protected:
  void handleSelectionChanged();
  void handleCollapsed( const Wt::WModelIndex &index );
  void handleExpanded( const Wt::WModelIndex &index );
  void tryExpandNode( const Wt::WModelIndex &index );
  void handleLayoutAboutToBeChanged();
  void handleLayoutChanged();

  /** Store the last selected path into the user preference name "RefSpectraLastSelection". 
   
   This function is called from the destructor of this widget.
  */
  void storeLastSelectedPath() noexcept;

#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
  void startAddDirectory();
  void removeCurrentlySelectedDirectory();
  void addDirectory( const std::string &path );
#endif // !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT
  void updatePreview();
  void updateTreeView();
  void updateCheckboxVisibility();

  //Some commented out code to help prepair the default ref spectra from the source directory.
  //static void dev_code();
private:
  RefSpectraModel *m_treeModel;
  Wt::WTreeView *m_treeView;
  Wt::WStackedWidget *m_stack;
  D3SpectrumDisplayDiv *m_spectrum;
  Wt::WContainerWidget *m_dirInfoContainer;
  Wt::WCheckBox *m_showCurrentForeground;
  Wt::WCheckBox *m_refBackground;
  Wt::WCheckBox *m_loadBackground;
#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
  Wt::WContainerWidget *m_addDirButton;
  Wt::WContainerWidget *m_deleteDirButton;
  Wt::WPushButton *m_showInExplorerButton;
#endif // !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT
  Wt::WComboBox *m_showAsComboBox;
  std::vector<std::string> m_baseDirectories;
  Wt::Signal<RefSpectraWidgetSelectionType> m_fileSelectionChangedSignal;
  
  /** We need to track last collapsed node so we can auto-expand nodes when they get selected,
   * but avoid expanding nodes that *just* got collapsed (colapsing them causes them to be selected).
   */
  Wt::WModelIndex m_lastCollapsedIndex;

  /** If the user uses the little "+" icon to expand a folder, the node will quickly open then collapse;
   not sure if its something I'm doing to cause this (probably), or Wt, but as a workaround, when we
   get the signal that a node is collapsed, if its `m_lastExpandedIndex` and its been less than
   half a second (arbitrary), then we will re-expand it.
   */
  Wt::WModelIndex m_lastExpandedIndex;
  std::chrono::time_point<std::chrono::steady_clock> m_lastExpandTime;

  void setupUI();
  void createTreeModel();
  void initBaseDirs();
#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
  void showInExplorer();
#endif
};

#endif // RefSpectraWidget_h 
