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
#include <limits>

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WTable>
#include <Wt/WString>
#include <Wt/WCheckBox>
#include <Wt/WComboBox>
#include <Wt/WSvgImage>
#include <Wt/WTableCell>
#include <Wt/WPushButton>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/AddNewPeakDialog.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"


using namespace std;
using namespace Wt;


const float AddNewPeakDialog::m_minfwhm = 0.05f;
const float AddNewPeakDialog::m_maxfwhm = 450.0f;  //reasonable range of peak widths


AddNewPeakDialog::AddNewPeakDialog( const float initialEnergy, const std::string &ref_line_hint )
: AuxWindow( WString::tr("anpd-dialog-add-peak-title"),
            (AuxWindowProperties::IsModal
             | AuxWindowProperties::TabletNotFullScreen
             | AuxWindowProperties::DisableCollapse) ),
m_renderFlags( 0 ),
m_viewer( InterSpec::instance() ),
m_peakModel( m_viewer->peakModel() ),
m_isPhone( m_viewer->isPhone() ),
m_meas( m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground) ),
m_candidatePeak( nullptr ),
m_energySB( nullptr ),
m_fwhmSB( nullptr ),
m_areaSB( nullptr ),
m_roiLowerSB( nullptr ),
m_roiUpperSB( nullptr ),
m_continuumType( nullptr ),
m_fitBtn( nullptr ),
m_fitEnergy( nullptr ),
m_fitFWHM( nullptr ),
m_fitAmplitude( nullptr ),
m_initial_src_hint( ref_line_hint ),
m_chart( nullptr )
{
  m_viewer->useMessageResourceBundle( "AddNewPeakDialog" );
  
  if( !m_peakModel || !m_meas || m_meas->num_gamma_channels() < 7 )
  {
    passMessage( WString::tr("anpd-err-no-foreground"), WarningWidget::WarningMsgHigh );
    return;
  }//if( we dont have a valid foreground )
  
  const size_t nbin = m_meas->num_gamma_channels();
  
  // Now create a dialog where user can specify the mean, sigma, roi range, and roi polynomial order...
  const float initialFWHM = estimateFWHM(initialEnergy);
  const double initialArea = std::max(10.0,0.2*m_meas->gamma_integral(initialEnergy-initialFWHM, initialEnergy+initialFWHM));
  
  m_candidatePeak = make_shared<PeakDef>(initialEnergy,initialFWHM/2.35482,initialArea);
  
  const bool isHPGe = PeakFitUtils::is_likely_high_res( m_viewer );
  
  size_t roi_lower_channel, roi_upper_channel;
  estimatePeakFitRange( *m_candidatePeak, m_meas, isHPGe, roi_lower_channel, roi_upper_channel );
  //I think we are garunteed the bins to be in range, but we'll enforce this JIC
  roi_upper_channel = std::min( roi_upper_channel, m_meas->num_gamma_channels()-1 );
  
  const double initialRoiLower = m_meas->gamma_channel_lower( roi_lower_channel );
  const double initialRoiUpper = m_meas->gamma_channel_upper( roi_upper_channel );
  m_candidatePeak->continuum()->setRange( initialRoiLower, initialRoiUpper );
  
  const float minEnergy = m_meas->gamma_channel_lower(0);
  const float maxEnergy = m_meas->gamma_channel_upper(nbin-1);
  
  
  WTable *table = new WTable( contents() );
  table->addStyleClass( "AddPeakTbl" );
  table->setHeaderCount( 1, Wt::Orientation::Vertical );
  
  if( !wApp->styleSheet().isDefined("AddPeakTblCell") )
    wApp->styleSheet().addRule( ".AddPeakTbl > tbody > tr > td", "padding-left: 10px; vertical-align: middle; padding-top: 3px", "AddPeakTblCell" );
  
  if( m_isPhone && !wApp->styleSheet().isDefined("AddPeakTblMbl") )
    wApp->styleSheet().addRule( ".AddPeakTbl", "width: 100%; margin: 10px;", "AddPeakTblMbl" );
  
  
  WLabel *label = new WLabel( WString::tr("anpd-peak-mean") );
  table->elementAt(0,0)->addWidget( label );
  
  const WLength inputWidth( 75, WLength::Unit::Pixel );
  
  m_energySB = new NativeFloatSpinBox();
  m_energySB->setWidth( inputWidth );
  table->elementAt(0,1)->addWidget( m_energySB );
  m_energySB->setRange( minEnergy, maxEnergy );
  m_energySB->setSpinnerHidden( true );
  m_energySB->setValue( initialEnergy );
  
  
  label = new WLabel( WString::tr("anpd-peak-fwhm") );
  table->elementAt(1,0)->addWidget( label );
  
  m_fwhmSB = new NativeFloatSpinBox();
  m_fwhmSB->setWidth( inputWidth );
  table->elementAt(1,1)->addWidget( m_fwhmSB );
  m_fwhmSB->setRange( m_meas->gamma_channel_lower(0), m_meas->gamma_channel_upper(nbin-1) );
  m_fwhmSB->setRange( m_minfwhm, m_maxfwhm );
  m_fwhmSB->setSpinnerHidden( true );
  m_fwhmSB->setValue( initialFWHM );
  
  
  label = new WLabel( WString::tr("anpd-peak-amp") );
  table->elementAt(2,0)->addWidget( label );
  m_areaSB = new NativeFloatSpinBox();
  m_areaSB->setWidth( inputWidth );
  table->elementAt(2,1)->addWidget( m_areaSB );
  m_areaSB->setRange( 0, m_meas->gamma_channels_sum(0, nbin-1) );
  m_areaSB->setSpinnerHidden( true );
  m_areaSB->setValue( initialArea );
  
  label = new WLabel( WString::tr("anpd-roi-lower") );
  table->elementAt(3,0)->addWidget( label );
  m_roiLowerSB = new NativeFloatSpinBox();
  m_roiLowerSB->setWidth( inputWidth );
  table->elementAt(3,1)->addWidget( m_roiLowerSB );
  m_roiLowerSB->setRange( minEnergy, maxEnergy );
  m_roiLowerSB->setSpinnerHidden( true );
  m_roiLowerSB->setValue( initialRoiLower );
  
  label = new WLabel( WString::tr("anpd-roi-upper") );
  table->elementAt(4,0)->addWidget( label );
  m_roiUpperSB = new NativeFloatSpinBox();
  m_roiUpperSB->setWidth( inputWidth );
  table->elementAt(4,1)->addWidget( m_roiUpperSB );
  m_roiUpperSB->setRange( minEnergy, maxEnergy );
  m_roiUpperSB->setSpinnerHidden( true );
  m_roiUpperSB->setValue( initialRoiUpper );
  
  
  label = new WLabel( WString::tr("Continuum") );
  table->elementAt(5,0)->addWidget( label );
  m_continuumType = new WComboBox();
  table->elementAt(5,1)->addWidget( m_continuumType );
  //layout->addWidget( m_continuumType, 5, 1 );
  m_continuumType->addItem( WString::tr("pct-none") );
  m_continuumType->addItem( WString::tr("pct-constant") );
  m_continuumType->addItem( WString::tr("pct-linear") );
  m_continuumType->addItem( WString::tr("pct-quadratic") );
  m_continuumType->addItem( WString::tr("pct-cubic") );
  m_continuumType->addItem( WString::tr("pct-flat-step") );
  m_continuumType->addItem( WString::tr("pct-linear-step") );
  m_continuumType->addItem( WString::tr("pct-bilinear-step") );
  m_continuumType->addItem( WString::tr("pct-global") );
  m_continuumType->setCurrentIndex( 2 );
  
  auto fitCell = table->elementAt(6,0);
  fitCell->setColumnSpan( 2 );
  fitCell->setVerticalAlignment( Wt::AlignmentFlag::AlignBottom );
  //layout->addWidget( fitDiv, 0, 6, 1, 2 );
  m_fitBtn = new WPushButton( WString::tr("Fit"), fitCell );
  m_fitEnergy = new WCheckBox( WString::tr("Energy"), fitCell );
  m_fitFWHM = new WCheckBox( WString::tr("FWHM"), fitCell );
  m_fitAmplitude = new WCheckBox( WString::tr("Amp."), fitCell );
  m_associateWithNuclide = new WCheckBox( WString::tr("anpd-associate-with").arg("N/A"), fitCell );
  m_fitEnergy->setMargin(3,Wt::Left);
  m_fitFWHM->setMargin(3,Wt::Left);
  m_fitAmplitude->setMargin(3,Wt::Left);
  m_associateWithNuclide->setMargin(3,Wt::Left);
  
  WFont fitCbFont;
  fitCbFont.setWeight( WFont::Weight::NormalWeight );
  fitCbFont.setSize( WFont::Size::Smaller );
  
  m_fitEnergy->decorationStyle().setFont( fitCbFont );
  m_fitFWHM->decorationStyle().setFont( fitCbFont );
  m_fitAmplitude->decorationStyle().setFont( fitCbFont );
  m_associateWithNuclide->decorationStyle().setFont( fitCbFont );
  m_associateWithNuclide->setHiddenKeepsGeometry( true );
  
  m_fitEnergy->setChecked( false );
  m_fitFWHM->setChecked( true );
  m_fitAmplitude->setChecked( true );
  
  m_chart = new WText( "", Wt::XHTMLUnsafeText );
  m_chart->setInline( false );
  
  const bool narrowLayout = (m_isPhone && (m_viewer->renderedWidth() < m_viewer->renderedHeight()));
  
  if( narrowLayout )
  {
    // Narrow vertical phone layout
    WTableCell *cell = table->elementAt(7,0);
    cell->addWidget( m_chart );
    cell->setRowSpan( 7 );
    cell->setColumnSpan( table->columnCount() );
    cell->setVerticalAlignment( Wt::AlignmentFlag::AlignMiddle );
    cell->setContentAlignment( Wt::AlignmentFlag::AlignCenter );
  }else
  {
    // Normal wide Layout
    WTableCell *cell = table->elementAt(0,2);
    cell->addWidget( m_chart );
    cell->setRowSpan( 7 );
    
    if( m_isPhone )
    {
      cell->setVerticalAlignment( Wt::AlignmentFlag::AlignMiddle );
      cell->setContentAlignment( Wt::AlignmentFlag::AlignCenter );
    }
  }//if( vertical phone ) / else
  
  
  WText *msg = new WText( WString::tr("anpd-tt-more-adv") );
  msg->decorationStyle().setFont( fitCbFont );
  if( m_isPhone )
  {
    WTableCell *cell = table->elementAt(table->rowCount(),0);
    cell->addWidget( msg );
    cell->setColumnSpan( table->columnCount() );
  }else
  {
    footer()->addWidget( msg );
    msg->setFloatSide( Wt::Side::Left );
  }
  
  
  m_energySB->valueChanged().connect( this, &AddNewPeakDialog::meanChanged );
  m_fwhmSB->valueChanged().connect( this, &AddNewPeakDialog::fwhmChanged );
  m_areaSB->valueChanged().connect( this, &AddNewPeakDialog::ampChanged );
  m_roiLowerSB->valueChanged().connect( this, &AddNewPeakDialog::roiRangeChanged );
  m_roiUpperSB->valueChanged().connect( this, &AddNewPeakDialog::roiRangeChanged );
  m_continuumType->activated().connect( this, &AddNewPeakDialog::roiTypeChanged );
  m_fitBtn->clicked().connect( this, &AddNewPeakDialog::doFit );
  m_associateWithNuclide->checked().connect( this, &AddNewPeakDialog::associateWithSrcChanged );
  m_associateWithNuclide->unChecked().connect( this, &AddNewPeakDialog::associateWithSrcChanged );
  
  WPushButton *closeButton = addCloseButtonToFooter(WString::tr("Cancel"));
  closeButton->clicked().connect( this, &AuxWindow::hide );
  WPushButton *doAdd = new WPushButton( WString::tr("Add"), footer() );
  
  doAdd->clicked().connect( std::bind( [this](){
    UndoRedoManager::PeakModelChange peak_undo_creator;
    m_viewer->addPeak( *m_candidatePeak, false );
    this->hide();
  }) );
  
  
  roiTypeChanged();
  m_renderFlags |= RenderActions::UpdatePreview;
  scheduleRender();
  
  rejectWhenEscapePressed();
  
  show();
  resizeToFitOnScreen();
  centerWindowHeavyHanded();
}//AddNewPeakDialog constructor


void AddNewPeakDialog::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  if( m_renderFlags.testFlag(RenderActions::UpdatePreview) )
  {
    updateCandidatePeakPreview();
    updateAssociateWithCb();
  }
  
  AuxWindow::render( flags );
  
  m_renderFlags = 0;
}//void render( Wt::WFlags<Wt::RenderFlag> flags )


float AddNewPeakDialog::estimateFWHM( const float energy )
{
  float fwhm = PeakSearchGuiUtils::estimate_FWHM_of_foreground( energy );
  
  return std::min( m_maxfwhm, std::max(m_minfwhm,fwhm) );
}//float estimateFWHM( const float energy )


void AddNewPeakDialog::updateCandidatePeakPreview()
{
  /** Setting the width makes the AuxWindow 517 px wide
   
   */
  const int ww = m_viewer->renderedWidth();
  const int wh = m_viewer->renderedHeight();
  const bool narrowLayout = (m_isPhone && (m_viewer->renderedWidth() < m_viewer->renderedHeight()));
  
  if( (!narrowLayout && (ww < 350)) || (narrowLayout && ww < 310) )
  {
    m_chart->setText( WString::tr("anpd-screen-to-small") );
    return;
  }
  
  auto meas = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  
  if( !meas )
  {
    m_chart->setText( WString::tr("anpd-err-no-foreground-1") );
    return;
  }
  
  auto peaks = std::make_shared< std::deque<std::shared_ptr<const PeakDef> > >();
  peaks->push_back( m_candidatePeak );
  
  const bool compact = false;
  const int width_px = m_isPhone ? (narrowLayout ? (ww - 20) : (ww - 250)) : std::min(550,ww-250);
  const int height_px = narrowLayout ? std::min(350,wh-225) : std::min(350,wh-50);
  const double roiWidth = m_candidatePeak->upperX() - m_candidatePeak->lowerX();
  const double lower_energy = m_candidatePeak->lowerX() - roiWidth;
  const double upper_energy = m_candidatePeak->upperX() + roiWidth;
  std::shared_ptr<const ColorTheme> theme = m_viewer->getColorTheme();
  const std::vector<std::shared_ptr<const ReferenceLineInfo>> reflines;
  
  std::shared_ptr<Wt::WSvgImage> preview
  = PeakSearchGuiUtils::renderChartToSvg( meas, peaks, reflines, lower_energy,
                                         upper_energy, width_px, height_px, theme, compact );
  
  if( preview )
  {
    stringstream strm;
    preview->write( strm );
    m_chart->setText( strm.str() );
  }else
  {
    m_chart->setText( WString::tr("anpd-err-preview") );
  }
}//void updateCandidatePeakPreview()


void AddNewPeakDialog::updateAssociateWithCb()
{
  if( !m_candidatePeak )
  {
    m_associateWithNuclide->hide();
    return;
  }
  
  
  
  const pair<unique_ptr<ReferenceLineInfo>,int> nearby_ref_line
      = PeakSearchGuiUtils::reference_line_near_peak( m_viewer, *m_candidatePeak, false, m_initial_src_hint );
  
  const ReferenceLineInfo::RefLine *line = nullptr;
  if( nearby_ref_line.first && (nearby_ref_line.second < static_cast<int>(nearby_ref_line.first->m_ref_lines.size())) )
    line = &(nearby_ref_line.first->m_ref_lines[nearby_ref_line.second]);
  
  if( line && (line->m_parent_nuclide || line->m_element || line->m_reaction) )
  {
    string src;
    if( line->m_parent_nuclide )
      src = line->m_parent_nuclide->symbol;
    else if( line->m_element )
      src = line->m_element->symbol;
    else if( line->m_reaction )
      src = line->m_reaction->name();
    
    char buffer[32] = { '\0' };
    snprintf( buffer, sizeof(buffer), " %.1f keV", line->m_energy );
    src += buffer + line->particlestr();
    
    m_associateWithNuclide->setText( WString::tr("anpd-associate-with").arg(src) );
  }else
  {
    m_associateWithNuclide->setText( WString::tr("anpd-associate-with").arg("N/A") );
    m_associateWithNuclide->hide();
    m_candidatePeak->clearSources();
  }
}//void AddNewPeakDialog::updateAssociateWithCb()


void AddNewPeakDialog::roiTypeChanged()
{
  auto meas = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  
  if( !meas )
  {
    m_chart->setText( WString::tr("anpd-err-no-foreground-1") );
    return;
  }
  
  const double roi_lower = m_candidatePeak->lowerX();
  const double roi_upper = m_candidatePeak->upperX();
  const double ref_energy = 0.5*(roi_upper + roi_lower);
  
  const WString txt = m_continuumType->currentText();
  PeakContinuum::OffsetType offset = PeakContinuum::OffsetType::Linear;
  if( txt.key() == "pct-none" )
    offset = PeakContinuum::OffsetType::NoOffset;
  else if( txt.key() == "pct-constant" )
    offset = PeakContinuum::OffsetType::Constant;
  else if( txt.key() == "pct-linear" )
    offset = PeakContinuum::OffsetType::Linear;
  else if( txt.key() == "pct-quadratic" )
    offset = PeakContinuum::OffsetType::Quadratic;
  else if( txt.key() == "pct-cubic" )
    offset = PeakContinuum::OffsetType::Cubic;
  else if( txt.key() == "pct-flat-step" )
    offset = PeakContinuum::OffsetType::FlatStep;
  else if( txt.key() == "pct-linear-step" )
    offset = PeakContinuum::OffsetType::LinearStep;
  else if( txt.key() == "pct-bilinear-step" )
    offset = PeakContinuum::OffsetType::BiLinearStep;
  else if( txt.key() == "pct-global" )
    offset = PeakContinuum::OffsetType::External;
  else
  {
    cerr << "Totally whack yo - unrecognized continuum type: {val: '" << txt.toUTF8() << "', key: '" << txt.key() << "'}" << endl;
    m_continuumType->setCurrentIndex( 2 );
  }
  
  
  switch( offset )
  {
    case PeakContinuum::OffsetType::NoOffset: //None
      m_candidatePeak->continuum()->setType( offset );
      break;
      
    case PeakContinuum::OffsetType::Constant:
    case PeakContinuum::OffsetType::Linear:
    case PeakContinuum::OffsetType::Quadratic:
    case PeakContinuum::OffsetType::Cubic:
    case PeakContinuum::OffsetType::FlatStep:
    case PeakContinuum::OffsetType::LinearStep:
    case PeakContinuum::OffsetType::BiLinearStep:
      m_candidatePeak->continuum()->calc_linear_continuum_eqn( meas, ref_energy, roi_lower, roi_upper, 5, 5 );
      m_candidatePeak->continuum()->setType( offset );
      break;
      
    case PeakContinuum::OffsetType::External:
    {
      auto gcontinuum = m_candidatePeak->continuum()->externalContinuum();
      if( !gcontinuum )
        gcontinuum = estimateContinuum( meas );
      m_candidatePeak->continuum()->setType( PeakContinuum::OffsetType::External );
      m_candidatePeak->continuum()->setExternalContinuum( gcontinuum );
      break;
    }
  }//switch( m_continuumType->currentIndex() )
  
  m_renderFlags |= RenderActions::UpdatePreview;
  scheduleRender();
}//void roiTypeChanged()


void AddNewPeakDialog::meanChanged()
{
  const double oldmean = m_candidatePeak->mean();
  
  if( WValidator::State::Valid != m_energySB->validate() )
    m_energySB->setValue( oldmean );
  
  const double newmean = m_energySB->value();
  const double change = newmean - oldmean;
  
  const double newLowerX = m_candidatePeak->lowerX() + change;
  const double newUpperX = m_candidatePeak->upperX() + change;
  m_roiLowerSB->setValue( newLowerX );
  m_roiUpperSB->setValue( newUpperX );
  
  m_candidatePeak->setMean( newmean );
  m_candidatePeak->continuum()->setRange( newLowerX, newUpperX );
  roiTypeChanged();
}//void meanChanged()


void AddNewPeakDialog::fwhmChanged()
{
  const double oldFWHM = m_candidatePeak->fwhm();
  
  if( WValidator::State::Valid != m_fwhmSB->validate() )
    m_candidatePeak->setSigma( oldFWHM/2.35482 );
  
  const double newFWHM = m_fwhmSB->value();
  m_candidatePeak->setSigma( newFWHM/2.35482 );
  
  m_renderFlags |= RenderActions::UpdatePreview;
  scheduleRender();
}//void fwhmChanged()


void AddNewPeakDialog::ampChanged()
{
  const double oldAmp = m_candidatePeak->amplitude();
  
  if( WValidator::State::Valid != m_areaSB->validate() )
    m_areaSB->setValue( oldAmp );
  
  const double newAmp = m_areaSB->value();
  m_candidatePeak->setAmplitude( newAmp );
  
  m_renderFlags |= RenderActions::UpdatePreview;
  scheduleRender();
}//void ampChanged()


void AddNewPeakDialog::roiRangeChanged()
{
  if( WValidator::State::Valid != m_roiLowerSB->validate() )
    m_roiLowerSB->setValue( m_candidatePeak->lowerX() );
  
  if( WValidator::State::Valid != m_roiUpperSB->validate() )
    m_roiUpperSB->setValue( m_candidatePeak->upperX() );
  
  double roiLower = m_roiLowerSB->value();
  double roiUpper = m_roiUpperSB->value();
  if( roiLower >= roiUpper )
  {
    m_roiLowerSB->setValue( roiUpper );
    m_roiUpperSB->setValue( roiLower );
    std::swap( roiLower, roiUpper );
  }
  
  if( m_candidatePeak->mean() < roiLower )
  {
    m_energySB->setValue( roiLower );
    m_candidatePeak->setMean( roiLower );
  }
  
  if( m_candidatePeak->mean() > roiUpper )
  {
    m_energySB->setValue( roiUpper );
    m_candidatePeak->setMean( roiUpper );
  }
  
  m_candidatePeak->continuum()->setRange( roiLower, roiUpper );
  
  roiTypeChanged();
}//void roiRangeChanged()


void AddNewPeakDialog::associateWithSrcChanged()
{
  setRefLineSrcToPeak();
  
  m_renderFlags |= RenderActions::UpdatePreview;
  scheduleRender();
}//void associateWithSrcChanged()


void AddNewPeakDialog::setRefLineSrcToPeak()
{
  if( !m_candidatePeak )
    return;
  
  const bool asscociate = (m_associateWithNuclide->isVisible() && m_associateWithNuclide->isChecked());
  
  const pair<unique_ptr<ReferenceLineInfo>,int> nearby_ref_line
    = PeakSearchGuiUtils::reference_line_near_peak( m_viewer, *m_candidatePeak, false, m_initial_src_hint );
  
  const ReferenceLineInfo::RefLine *line = nullptr;
  PeakDef::SourceGammaType src_type = PeakDef::SourceGammaType::NormalGamma;
  if( nearby_ref_line.first && (nearby_ref_line.second < static_cast<int>(nearby_ref_line.first->m_ref_lines.size())) )
  {
    line = &(nearby_ref_line.first->m_ref_lines[nearby_ref_line.second]);
    
    switch( line->m_source_type )
    {
      case ReferenceLineInfo::RefLine::RefGammaType::Normal:
        src_type = PeakDef::SourceGammaType::NormalGamma;
        break;
        
      case ReferenceLineInfo::RefLine::RefGammaType::Annihilation:
        src_type = PeakDef::SourceGammaType::AnnihilationGamma;
        break;
        
      case ReferenceLineInfo::RefLine::RefGammaType::SingleEscape:
        src_type = PeakDef::SourceGammaType::SingleEscapeGamma;
        break;
        
      case ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape:
        src_type = PeakDef::SourceGammaType::DoubleEscapeGamma;
        break;
        
      case ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak:
      case ReferenceLineInfo::RefLine::RefGammaType::SumGammaPeak:
        line = nullptr; //Setting this realtion in peaks not supported yet...
        break;
    }//switch( line->m_source_type )
  }//if( nearby_ref_line.first && (nearby_ref_line.second < static_cast<int>(nearby_ref_line.first->m_ref_lines.size())) )
  
  
  if( asscociate && (line && (line->m_parent_nuclide || line->m_element || line->m_reaction)) )
  {
    if( line->m_parent_nuclide )
    {
      int particle_index = 0;
      double nearest_gamma = 999999.9;
      for( size_t i = 0; line->m_transition && (i < line->m_transition->products.size()); ++i )
      {
        if( (line->m_transition->products[i].type == SandiaDecay::ProductType::GammaParticle)
           || (line->m_transition->products[i].type == SandiaDecay::ProductType::XrayParticle) )
        {
          double dist = fabs(line->m_energy - line->m_transition->products[i].energy);
          if( dist < nearest_gamma )
          {
            nearest_gamma = dist;
            particle_index = static_cast<int>( i );
          }//if( dist < nearest_gamma )
        }
      }//for( loop over products )
      
      m_candidatePeak->setNuclearTransition( line->m_parent_nuclide, line->m_transition, particle_index, src_type );
    }else if( line->m_element )
    {
      m_candidatePeak->setXray( line->m_element, line->m_energy );
    }else if( line->m_reaction )
    {
      m_candidatePeak->setReaction( line->m_reaction, line->m_energy, src_type );
    }
    
    string src;
    if( line->m_parent_nuclide )
      src = line->m_parent_nuclide->symbol;
    else if( line->m_element )
      src = line->m_element->symbol;
    else if( line->m_reaction )
      src = line->m_reaction->name();
    
    Wt::WColor color;
    ReferencePhotopeakDisplay * const refTool = m_viewer->referenceLinesWidget();
    if( !src.empty() && refTool )
      color = refTool->suggestColorForSource(src);
    if( color.isDefault() )
      color = line->m_color;
    m_candidatePeak->setLineColor( color );
  }else
  {
    m_candidatePeak->clearSources();
    m_candidatePeak->setLineColor( WColor{} );
  }//if( m_associateWithNuclide->isVisible() && m_associateWithNuclide->isChecked() )
}//void setRefLineSrcToPeak()


void AddNewPeakDialog::doFit()
{
  auto meas = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  
  if( !meas )
    return;
  
  m_candidatePeak->setFitFor( PeakDef::CoefficientType::Mean, m_fitEnergy->isChecked() );
  m_candidatePeak->setFitFor( PeakDef::CoefficientType::Sigma, m_fitFWHM->isChecked() );
  m_candidatePeak->setFitFor( PeakDef::CoefficientType::GaussAmplitude, m_fitAmplitude->isChecked() );
  
  const bool isHPGe = PeakFitUtils::is_likely_high_res( m_viewer ); //Or we could use the current peaks FWHM to estimate this...
  
  vector<PeakDef> input_peaks( 1, *m_candidatePeak ), results;
  const double stat_threshold  = 0.0, hypothesis_threshold = 0.0;
  
  Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options;
  fitPeaks( input_peaks, stat_threshold, hypothesis_threshold, meas, results, fit_options, isHPGe );
  
  if( results.empty() )
  {
    passMessage( WString::tr("anpd-err-fit-failed"), WarningWidget::WarningMsgLow );
  }else
  {
    *m_candidatePeak = results[0];
    m_energySB->setValue( m_candidatePeak->mean() );
    m_fwhmSB->setValue( m_candidatePeak->fwhm() );
    m_areaSB->setValue( m_candidatePeak->amplitude() );
    
    setRefLineSrcToPeak();
  }//if( results.empty() ) / else
  
  m_renderFlags |= RenderActions::UpdatePreview;
  scheduleRender();
};//void doFit()
