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
#include <Wt/WLabel>
#include <Wt/WCheckBox>
#include <Wt/WLineEdit>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WContainerWidget>

#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/FarmOptionsWindow.h"

using namespace Wt;
using namespace std;


FarmOptionsWindow::FarmOptionsWindow( InterSpec *viewer )
  : AuxWindow( WString::tr("fow-window-title"),
               (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
                | AuxWindowProperties::TabletNotFullScreen
                | AuxWindowProperties::DisableCollapse) ),
    m_interspec( viewer ),
    m_enableFarm( nullptr ),
    m_enableGadras( nullptr ),
    m_gadrasExePath( nullptr ),
    m_synthesizeBackground( nullptr ),
    m_enableRelAct( nullptr ),
    m_enableFram( nullptr ),
    m_framExePath( nullptr ),
    m_framOutputPath( nullptr ),
    m_writeFertilized( nullptr )
{
  viewer->useMessageResourceBundle( "FarmOptionsWindow" );
  wApp->useStyleSheet( "InterSpec_resources/DirectorySelector.css" );
  init();
  loadFromPreferences();
}//FarmOptionsWindow constructor


FarmOptionsWindow::~FarmOptionsWindow()
{
}//~FarmOptionsWindow()


void FarmOptionsWindow::init()
{
  WContainerWidget *content = contents();

  WGridLayout *layout = new WGridLayout();
  content->setLayout( layout );

  int row = 0;

  // Master enable checkbox
  m_enableFarm = new WCheckBox( WString::tr("fow-enable-farm") );
  layout->addWidget( m_enableFarm, row, 0, 1, 2 );
  m_enableFarm->changed().connect( this, &FarmOptionsWindow::handleEnableChanged );
  HelpSystem::attachToolTipOn( m_enableFarm, WString::tr("fow-tt-enable-farm"), true );
  row++;

  // Separator
  WText *sep1 = new WText( "<hr />" );
  sep1->setTextFormat( Wt::XHTMLText );
  layout->addWidget( sep1, row, 0, 1, 2 );
  row++;

  // GADRAS section header
  WText *gadrasHeader = new WText( WString::tr("fow-gadras-header") );
  gadrasHeader->addStyleClass( "FarmSectionHeader" );
  layout->addWidget( gadrasHeader, row, 0, 1, 2 );
  row++;

  // GADRAS enable checkbox
  m_enableGadras = new WCheckBox( WString::tr("fow-enable-gadras") );
  layout->addWidget( m_enableGadras, row, 0, 1, 2 );
  m_enableGadras->changed().connect( this, &FarmOptionsWindow::handleEnableChanged );
  HelpSystem::attachToolTipOn( m_enableGadras, WString::tr("fow-tt-enable-gadras"), true );
  row++;

  // GADRAS Full Spectrum exe path
  WLabel *gadrasLabel = new WLabel( WString::tr("fow-gadras-path") );
  layout->addWidget( gadrasLabel, row, 0 );
  m_gadrasExePath = new WLineEdit();
  m_gadrasExePath->addStyleClass( "DirectorySelectorPathTxt" );
  m_gadrasExePath->setPlaceholderText( WString::tr("fow-gadras-path-placeholder") );
  m_gadrasExePath->changed().connect( this, &FarmOptionsWindow::validatePathInputs );
  m_gadrasExePath->keyPressed().connect( this, &FarmOptionsWindow::validatePathInputs );
  layout->addWidget( m_gadrasExePath, row, 1 );
  row++;

  // Synthesize background checkbox
  m_synthesizeBackground = new WCheckBox( WString::tr("fow-synthesize-background") );
  layout->addWidget( m_synthesizeBackground, row, 0, 1, 2 );
  HelpSystem::attachToolTipOn( m_synthesizeBackground, WString::tr("fow-tt-synthesize-background"), true );
  row++;

  // RelActCalcAuto section header
  WText *relactHeader = new WText( WString::tr("fow-relact-header") );
  relactHeader->addStyleClass( "FarmSectionHeader" );
  layout->addWidget( relactHeader, row, 0, 1, 2 );
  row++;

  // RelActCalcAuto enable checkbox
  m_enableRelAct = new WCheckBox( WString::tr("fow-enable-relact") );
  layout->addWidget( m_enableRelAct, row, 0, 1, 2 );
  HelpSystem::attachToolTipOn( m_enableRelAct, WString::tr("fow-tt-enable-relact"), true );
  row++;

  // FRAM section header
  WText *framHeader = new WText( WString::tr("fow-fram-header") );
  framHeader->addStyleClass( "FarmSectionHeader" );
  layout->addWidget( framHeader, row, 0, 1, 2 );
  row++;

  // FRAM enable checkbox
  m_enableFram = new WCheckBox( WString::tr("fow-enable-fram") );
  layout->addWidget( m_enableFram, row, 0, 1, 2 );
  m_enableFram->changed().connect( this, &FarmOptionsWindow::handleEnableChanged );
  HelpSystem::attachToolTipOn( m_enableFram, WString::tr("fow-tt-enable-fram"), true );
  row++;

  // FRAM exe path
  WLabel *framLabel = new WLabel( WString::tr("fow-fram-path") );
  layout->addWidget( framLabel, row, 0 );
  m_framExePath = new WLineEdit();
  m_framExePath->addStyleClass( "DirectorySelectorPathTxt" );
  m_framExePath->setPlaceholderText( WString::tr("fow-fram-path-placeholder") );
  m_framExePath->changed().connect( this, &FarmOptionsWindow::validatePathInputs );
  m_framExePath->keyPressed().connect( this, &FarmOptionsWindow::validatePathInputs );
  layout->addWidget( m_framExePath, row, 1 );
  row++;

  // FRAM output path
  WLabel *framOutLabel = new WLabel( WString::tr("fow-fram-output-path") );
  layout->addWidget( framOutLabel, row, 0 );
  m_framOutputPath = new WLineEdit();
  m_framOutputPath->addStyleClass( "DirectorySelectorPathTxt" );
  m_framOutputPath->setPlaceholderText( WString::tr("fow-fram-output-placeholder") );
  m_framOutputPath->changed().connect( this, &FarmOptionsWindow::validatePathInputs );
  m_framOutputPath->keyPressed().connect( this, &FarmOptionsWindow::validatePathInputs );
  layout->addWidget( m_framOutputPath, row, 1 );
  row++;

  // Output section header
  WText *outputHeader = new WText( WString::tr("fow-output-header") );
  outputHeader->addStyleClass( "FarmSectionHeader" );
  layout->addWidget( outputHeader, row, 0, 1, 2 );
  row++;

  // Write fertilized N42 checkbox
  m_writeFertilized = new WCheckBox( WString::tr("fow-write-fertilized") );
  layout->addWidget( m_writeFertilized, row, 0, 1, 2 );
  HelpSystem::attachToolTipOn( m_writeFertilized, WString::tr("fow-tt-write-fertilized"), true );
  row++;

  // Footer buttons
  WContainerWidget *buttonDiv = footer();

  WPushButton *cancelBtn = new WPushButton( WString::tr("Cancel"), buttonDiv );
  cancelBtn->clicked().connect( this, &FarmOptionsWindow::hide );

  WPushButton *saveBtn = new WPushButton( WString::tr("Save"), buttonDiv );
  saveBtn->clicked().connect( this, &FarmOptionsWindow::closeAndSave );

  // Initialize enable/disable states
  handleEnableChanged();

  // Allow escape to close
  rejectWhenEscapePressed();

  // Set window size
  setWidth( 500 );
  centerWindow();
}//init()




Farm::FarmOptions FarmOptionsWindow::currentOptions() const
{
  Farm::FarmOptions opts;

  opts.enable_farm_analysis = m_enableFarm->isChecked();
  opts.enable_gadras_rid = m_enableGadras->isChecked();
  opts.gadras_exe_path = m_gadrasExePath->text().toUTF8();
  opts.synthesize_background_if_missing = m_synthesizeBackground->isChecked();
  opts.enable_relact_isotopics = m_enableRelAct->isChecked();
  opts.enable_fram_isotopics = m_enableFram->isChecked();
  opts.fram_exe_path = m_framExePath->text().toUTF8();
  opts.fram_output_path = m_framOutputPath->text().toUTF8();
  opts.write_fertilized_n42 = m_writeFertilized->isChecked();

  return opts;
}//currentOptions()


void FarmOptionsWindow::setOptions( const Farm::FarmOptions &opts )
{
  m_enableFarm->setChecked( opts.enable_farm_analysis );
  m_enableGadras->setChecked( opts.enable_gadras_rid );
  m_gadrasExePath->setText( WString::fromUTF8(opts.gadras_exe_path) );
  m_synthesizeBackground->setChecked( opts.synthesize_background_if_missing );
  m_enableRelAct->setChecked( opts.enable_relact_isotopics );
  m_enableFram->setChecked( opts.enable_fram_isotopics );
  m_framExePath->setText( WString::fromUTF8(opts.fram_exe_path) );
  m_framOutputPath->setText( WString::fromUTF8(opts.fram_output_path) );
  m_writeFertilized->setChecked( opts.write_fertilized_n42 );

  handleEnableChanged();
  validatePathInputs();
}//setOptions()


void FarmOptionsWindow::saveToPreferences()
{
  const Farm::FarmOptions opts = currentOptions();
  const std::string json = opts.toJson();

  UserPreferences::setPreferenceValue<std::string>( "FarmOptions", json, m_interspec );

  m_optionsChanged.emit( opts );
}//saveToPreferences()

void FarmOptionsWindow::closeAndSave()
{
  saveToPreferences();
  hide();
}


void FarmOptionsWindow::loadFromPreferences()
{
  try
  {
    const std::string json = UserPreferences::preferenceValue<std::string>( "FarmOptions", m_interspec );
    const Farm::FarmOptions opts = Farm::FarmOptions::fromJson( json );
    setOptions( opts );
  }catch( std::exception &e )
  {
    // Use defaults on error
    cerr << "Failed to load FARM preferences:" << e.what() << endl;
    setOptions( Farm::FarmOptions() );
  }
}//loadFromPreferences()


void FarmOptionsWindow::handleEnableChanged()
{
  const bool enabled = m_enableFarm->isChecked();

  // GADRAS controls
  m_enableGadras->setEnabled( enabled );
  const bool gadrasEnabled = enabled && m_enableGadras->isChecked();
  m_gadrasExePath->setEnabled( gadrasEnabled );
  m_synthesizeBackground->setEnabled( gadrasEnabled );

  // RelActCalcAuto controls
  m_enableRelAct->setEnabled( enabled );

  // FRAM controls
  m_enableFram->setEnabled( enabled );
  const bool framEnabled = enabled && m_enableFram->isChecked();
  m_framExePath->setEnabled( framEnabled );
  m_framOutputPath->setEnabled( framEnabled );

  // Output controls
  m_writeFertilized->setEnabled( enabled );
}//handleEnableChanged()


void FarmOptionsWindow::validatePathInputs()
{
  // Helper: applies valid/invalid styling to a WLineEdit following the same
  // pattern as DirectorySelector::validatePath().  InvalidPath class is only
  // applied when the field is non-empty (empty == no styling, field is optional).
  auto applyPathStyle = []( Wt::WLineEdit *edit, const bool isValid )
  {
    if( isValid )
    {
      if( edit->hasStyleClass( "DirectorySelectorInvalidPath" ) )
        edit->removeStyleClass( "DirectorySelectorInvalidPath" );
      if( !edit->hasStyleClass( "DirectorySelectorValidPath" ) )
        edit->addStyleClass( "DirectorySelectorValidPath" );
    }else
    {
      if( edit->hasStyleClass( "DirectorySelectorValidPath" ) )
        edit->removeStyleClass( "DirectorySelectorValidPath" );
      if( !edit->hasStyleClass( "DirectorySelectorInvalidPath" ) )
        edit->addStyleClass( "DirectorySelectorInvalidPath" );
    }
  };

  // GADRAS Full Spectrum exe: must be an existing file; empty clears styling (field is optional)
  {
    const std::string path = m_gadrasExePath->text().toUTF8();
    if( path.empty() )
    {
      m_gadrasExePath->removeStyleClass( "DirectorySelectorValidPath" );
      m_gadrasExePath->removeStyleClass( "DirectorySelectorInvalidPath" );
    }else
      applyPathStyle( m_gadrasExePath, SpecUtils::is_file( path ) );
  }

  // FRAM exe: same logic as GADRAS Full Spectrum exe
  {
    const std::string path = m_framExePath->text().toUTF8();
    if( path.empty() )
    {
      m_framExePath->removeStyleClass( "DirectorySelectorValidPath" );
      m_framExePath->removeStyleClass( "DirectorySelectorInvalidPath" );
    }else
      applyPathStyle( m_framExePath, SpecUtils::is_file( path ) );
  }

  // FRAM output path: the target file itself may not exist yet, so we validate
  // that its parent directory exists
  {
    const std::string path = m_framOutputPath->text().toUTF8();
    if( path.empty() )
    {
      m_framOutputPath->removeStyleClass( "DirectorySelectorValidPath" );
      m_framOutputPath->removeStyleClass( "DirectorySelectorInvalidPath" );
    }else
      applyPathStyle( m_framOutputPath, SpecUtils::is_directory( SpecUtils::parent_path( path ) ) );
  }
}//validatePathInputs()
