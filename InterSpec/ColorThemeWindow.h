#ifndef ColorThemeWindow_h
#define ColorThemeWindow_h
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

#include <Wt/WContainerWidget>

#include "InterSpec/AuxWindow.h"


//Forward declarations
namespace Wt
{
  class WMenu;
  class WMenuItem;
  class WPushButton;
}//namespace Wt

class InterSpec;
struct ColorTheme;
class ColorThemeWidget;

/** This widget is a window that displays the default and database saved color
 themes in a menu and allows the user to select between them, to allow editing
 and applying for use within the app.
 A #ColorThemeWidget handles displaying and editing a single color theme, but
 this widget, #ColorThemeWindow, handles displaying available color themes
 (either default in the app, or that the user has saved to a database), allowing
 the user to select the different themes, or if they edit a theme save to the
 database.  It also handles the logic for adding/removing themes from database
 as well as applying the actual theme to the app.
 */
class ColorThemeWindow : public AuxWindow
{
public:
  /** ColorThemeWindow constructor
   \param interspec Pointer to the owning InterSpec instance for this widget.
          Must be valid.
   */
  ColorThemeWindow( InterSpec *interspec );
  virtual ~ColorThemeWindow();
  
  /** Called when the user selects a new theme. */
  void themeSelected(Wt::WMenuItem *item);
  
  /** Saves the passed in color theme to the database, returning a pointer
      loaded from the DB, or nullptr if error
   */
  std::unique_ptr<ColorTheme> saveThemeToDb( const ColorTheme *theme );
  
  /** Returns title/name for current theme. */
  std::string currentThemeTitle();
  
  /** Returns JSON for current theme. */
  std::string currentThemeJson();
  
protected:
  /** Called when the user edits the current theme. */
  void themEditedCallback();
  
  /** Retrieves the users previously saved themes from the database. */
  std::vector<std::unique_ptr<ColorTheme>> userDbThemes();
  
  /** Called when the dialog is closed, and checks to see if the user
   should be prompted to save anything.
   */
  void checkForSavesAndCleanUp();
  void saveAndDelete();
  
  void saveCallback();
  void applyCallback();
  
  void cloneThemeCallback();
  void removeThemeCallback();
  void uploadThemeCallback();
  
  /** To work around Wt (at least <=3.3.4) bug where user had to click on the
   text of the menu item; this fixes to anywhere on the menu item.
   */
  void selectItem( Wt::WMenuItem *item );
  
  void showOrHideApplyButton();
  
private:
  InterSpec        *m_interspec;
  Wt::WMenu        *m_menu;
  ColorThemeWidget *m_edit;
  
  Wt::WContainerWidget *m_removeIcn;
  Wt::WContainerWidget *m_addIcn;
  
  Wt::WPushButton  *m_close;
  Wt::WPushButton  *m_save;
  Wt::WPushButton  *m_apply;
};//class ColorThemeWindow


#endif
