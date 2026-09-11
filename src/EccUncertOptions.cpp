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

#include <cmath>
#include <memory>
#include <vector>
#include <cstdio>

#include <Wt/WText.h>
#include <Wt/WLabel.h>
#include <Wt/WTable.h>
#include <Wt/WComboBox.h>
#include <Wt/WCheckBox.h>
#include <Wt/WLineEdit.h>
#include <Wt/WTableCell.h>
#include <Wt/WApplication.h>
#include <Wt/WDoubleValidator.h>
#include <Wt/WContainerWidget.h>

#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/EccUncertOptions.h"
#include "InterSpec/DetectorEfficiency.h"

using namespace Wt;
using namespace std;


EccUncertOptions::EccUncertOptions( const vector<float> &energies,
                                    const vector<float> &baselineFrac,
                                    const vector<float> &convergenceFrac )
  : WContainerWidget(),
    m_energies( energies ),
    m_baselineFrac( baselineFrac ),
    m_convergenceFrac( convergenceFrac ),
    m_import( nullptr ),
    m_mode( nullptr ),
    m_corrLenLabel( nullptr ),
    m_corrLen( nullptr ),
    m_exampleTitle( nullptr ),
    m_exampleTable( nullptr ),
    m_changed()
{
  addStyleClass( "EccUncertOptions" );

  WApplication * const wapp = WApplication::instance();
  if( wapp )
    wapp->useStyleSheet( "InterSpec_resources/EccUncertOptions.css" );

  InterSpecApp *app = dynamic_cast<InterSpecApp *>( wapp );
  if( app )
    app->useMessageResourceBundle( "EccUncertOptions" );

  m_import = addNew<WCheckBox>( WString::tr("euo-import-cb") );
  m_import->setChecked( true );
  m_import->changed().connect( this, &EccUncertOptions::handleImportToggled );

  WTable *table = addNew<WTable>();
  table->addStyleClass( "EccUncertOptsTable" );

  {
    table->elementAt( 0, 0 )->addNew<WLabel>( WString::tr("euo-mode-label") );
    m_mode = table->elementAt( 0, 1 )->addNew<WComboBox>();
    m_mode->addItem( WString::tr("euo-mode-fully-corr") );  //FullyCorrelated
    m_mode->addItem( WString::tr("euo-mode-gaussian") );    //Gaussian
    m_mode->addItem( WString::tr("euo-mode-uncorr") );      //Uncorrelated
    m_mode->setCurrentIndex( FullyCorrelated );
    m_mode->activated().connect( this, &EccUncertOptions::handleModeChanged );
    HelpSystem::attachToolTipOn( m_mode, WString::tr("euo-tt-mode"), true );
  }

  {
    m_corrLenLabel = table->elementAt( 1, 0 )->addNew<WLabel>( WString::tr("euo-corrlen-label") );
    m_corrLen = table->elementAt( 1, 1 )->addNew<WLineEdit>();
    m_corrLen->setTextSize( 6 );
    char buf[32] = { '\0' };
    snprintf( buf, sizeof(buf), "%.2f", DetectorEfficiencyUncert::sm_defaultEccCorrLength );
    m_corrLen->setText( WString::fromUTF8(buf) );
    std::shared_ptr<WDoubleValidator> validator = std::make_shared<WDoubleValidator>( 0.0, 1.0E6 );
    validator->setMandatory( true );
    m_corrLen->setValidator( validator );
    m_corrLen->changed().connect( this, &EccUncertOptions::handleModeChanged );
    m_corrLenLabel->setBuddy( m_corrLen );
    HelpSystem::attachToolTipOn( m_corrLen, WString::tr("euo-tt-corrlen"), true );
    HelpSystem::attachToolTipOn( m_corrLenLabel, WString::tr("euo-tt-corrlen"), true );
  }

  m_exampleTitle = addNew<WText>( WString::tr("euo-example-title") );
  m_exampleTitle->addStyleClass( "EccUncertExampleTitle" );

  m_exampleTable = addNew<WTable>();
  m_exampleTable->addStyleClass( "EccUncertExampleTable" );
  m_exampleTable->setHeaderCount( 1 );  // Render the first row as <th> cells.

  handleModeChanged();
}//EccUncertOptions constructor


EccUncertOptions::~EccUncertOptions()
{
}


bool EccUncertOptions::importUncertainties() const
{
  return m_import->isChecked();
}//importUncertainties()


double EccUncertOptions::effectiveCorrLength() const
{
  switch( m_mode->currentIndex() )
  {
    case Uncorrelated:
      return -1.0;

    case Gaussian:
    {
      double value = 0.0;
      const string txt = m_corrLen->text().toUTF8();
      if( (sscanf( txt.c_str(), "%lf", &value ) == 1) && (value > 0.0) )
        return value;
      return DetectorEfficiencyUncert::sm_defaultEccCorrLength;
    }//case Gaussian

    case FullyCorrelated:
    default:
      break;
  }//switch( m_mode->currentIndex() )

  return DetectorEfficiencyUncert::sm_fullyCorrelatedLength;
}//effectiveCorrLength()


shared_ptr<DetectorEfficiencyUncert> EccUncertOptions::buildUncert() const
{
  if( !m_import->isChecked() || (m_energies.size() < 2) )
    return nullptr;

  try
  {
    return DetectorEfficiencyUncert::fromCorrelatedPlusDiagonal( m_energies,
                                        m_baselineFrac, m_convergenceFrac,
                                        effectiveCorrLength() );
  }catch( std::exception & )
  {
    return nullptr;
  }
}//buildUncert()


Wt::Signal<> &EccUncertOptions::changed()
{
  return m_changed;
}


void EccUncertOptions::handleImportToggled()
{
  const bool import = m_import->isChecked();
  m_mode->setEnabled( import );
  m_corrLen->setEnabled( import && (m_mode->currentIndex() == Gaussian) );
  m_exampleTitle->setHidden( !import );
  m_exampleTable->setHidden( !import );

  m_changed.emit();
}//handleImportToggled()


void EccUncertOptions::handleModeChanged()
{
  const bool gaussian = (m_mode->currentIndex() == Gaussian);
  m_corrLenLabel->setHidden( !gaussian );
  m_corrLen->setHidden( !gaussian );
  m_corrLen->setEnabled( m_import->isChecked() && gaussian );

  rebuildExampleTable();

  m_changed.emit();
}//handleModeChanged()


void EccUncertOptions::rebuildExampleTable()
{
  m_exampleTable->clear();

  if( m_energies.size() < 2 )
  {
    m_exampleTitle->hide();
    m_exampleTable->hide();
    return;
  }

  m_exampleTitle->setHidden( !m_import->isChecked() );
  m_exampleTable->setHidden( !m_import->isChecked() );

  const double corrLength = effectiveCorrLength();

  // rho(E1,E2) for the current mode; matches DetectorEfficiencyUncert's kernel.
  auto rho = [corrLength]( const double e1, const double e2 ) -> double {
    if( corrLength <= 0.0 )
      return (e1 == e2) ? 1.0 : 0.0;
    const double dlne = std::log(e1) - std::log(e2);
    return std::exp( -0.5 * std::pow( dlne / corrLength, 2.0 ) );
  };

  // A few representative energy pairs (keV) illustrating the correlation.
  const std::vector<std::pair<float,float>> pairs = {
    { 59.0f, 122.0f },
    { 84.0f, 185.0f },
    { 185.0f, 1001.0f },
    { 356.0f, 661.0f },
  };

  m_exampleTable->elementAt( 0, 0 )->addNew<WText>( WString::tr("euo-example-e1") );
  m_exampleTable->elementAt( 0, 1 )->addNew<WText>( WString::tr("euo-example-e2") );
  m_exampleTable->elementAt( 0, 2 )->addNew<WText>( WString::tr("euo-example-rho") );

  int trow = 1;
  for( const std::pair<float,float> &p : pairs )
  {
    // Skip degenerate pairs (e.g., when there are very few nodes).
    if( p.first == p.second )
      continue;

    char e1buf[32] = { '\0' }, e2buf[32] = { '\0' }, rbuf[32] = { '\0' };
    snprintf( e1buf, sizeof(e1buf), "%.4g", p.first );
    snprintf( e2buf, sizeof(e2buf), "%.4g", p.second );
    snprintf( rbuf, sizeof(rbuf), "%.3f", rho( p.first, p.second ) );

    m_exampleTable->elementAt( trow, 0 )->addNew<WText>( WString::fromUTF8(e1buf) );
    m_exampleTable->elementAt( trow, 1 )->addNew<WText>( WString::fromUTF8(e2buf) );
    m_exampleTable->elementAt( trow, 2 )->addNew<WText>( WString::fromUTF8(rbuf) );
    ++trow;
  }//for( const std::pair<float,float> &p : pairs )
}//rebuildExampleTable()
