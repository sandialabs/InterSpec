#ifndef DirectorySelector_h
#define DirectorySelector_h

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

#include <string>
#include <functional>

#include <Wt/WContainerWidget>
#include <Wt/WSignal>

namespace Wt
{
  class WLabel;
  class WText;
  class WLineEdit;
  class WPushButton;
}

#if( defined(IOS) || defined(ANDROID) || BUILD_FOR_WEB_DEPLOYMENT )
static_assert( false, "DirectorySelector is not supported on this platform" );
#endif

/** A widget for selecting directory paths that provides both native file dialogs 
 * (on supported platforms) and text input fallback.
 * 
 * This widget encapsulates directory selection functionality and provides:
 * - Native directory picker on Electron and macOS apps
 * - Text input fallback for other platforms
 * - Path validation
 * - Signal emission when path changes
 */
class DirectorySelector : public Wt::WContainerWidget
{
public:
  /** DirectorySelector constructor */
  DirectorySelector( Wt::WContainerWidget *parent = nullptr );
                     
  virtual ~DirectorySelector();

  /** Set the label text
    Default is "Path:"
   */
  void setLabelTxt( const Wt::WString &labelTxt );

  /** Get the currently selected/entered path
   * 
   * @return The current path as a string
   */
  std::string path() const;
  
  /** Set the current path
   * 
   * @param path The path to set
   * @param validate Whether to validate the path exists as a directory
   */
  void setPath( const std::string &path, bool validate = true );
  
  /** Check if the current path is valid (exists and is a directory)
   * 
   * @return true if path is valid, false otherwise
   */
  bool isPathValid() const;
  
  /** Enable or disable the widget
   * 
   * @param enabled Whether the widget should be enabled
   */
  void setEnabled( bool enabled );
  
  /** Get the signal emitted when the path changes
   * 
   * The signal is emitted with the new path string when:
   * - User selects a path via native dialog
   * - User enters/changes text in the input field
   * - Path is set programmatically via setPath()
   * 
   * @return Reference to the path changed signal
   */
  Wt::Signal<std::string> &pathChanged() { return m_pathChanged; }
  
  /** Get the signal emitted when path validation status changes
   * 
   * The signal is emitted with true when path becomes valid,
   * false when it becomes invalid.
   * 
   * @return Reference to the path validity changed signal
   */
  Wt::Signal<bool> &pathValidityChanged() { return m_pathValidityChanged; }

protected:
  /** Initialize the UI components */
  void setupUI();
  
  /** Validate the current path and update UI accordingly */
  void validatePath();
  
  /** Handle native directory selection (Electron/macOS) */
  void handleNativeDirectorySelection();
  
  /** Handle text input changes */
  void handleTextInputChange();

private:
  std::string m_currentPath;
  bool m_isValid;
  
  // UI components
  Wt::WLabel *m_label;
  
  // Native path selection components (when available)
  Wt::WText *m_pathDisplay;
  Wt::WPushButton *m_selectButton;
  
  // Text input fallback components
  Wt::WLineEdit *m_pathInput;
  
  // Signals
  Wt::Signal<std::string> m_pathChanged;
  Wt::Signal<bool> m_pathValidityChanged;
  
  bool m_useNativeDialog;
};

#endif // DirectorySelector_h 