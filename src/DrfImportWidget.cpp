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
#include <algorithm>

#include <Wt/Utils.h>
#include <Wt/WText.h>
#include <Wt/WLabel.h>
#include <Wt/WServer.h>
#include <Wt/WCheckBox.h>
#include <Wt/WLineEdit.h>
#include <Wt/WComboBox.h>
#include <Wt/WIOService.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>
#include <Wt/WEnvironment.h>
#include <Wt/WRegExpValidator.h>

#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/DrfChart.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/WidgetUtils.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/DrfImportWidget.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/EccUncertOptions.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/FileDragUploadResource.h"

using namespace std;
using namespace Wt;

using DrfImport::FileKind;
using DrfImport::Status;
using DrfImport::Interpretation;


namespace
{
  const char *kind_key( const FileKind kind )
  {
    switch( kind )
    {
      case FileKind::DrfXml:              return "dsi-kind-drf-xml";
      case FileKind::MakeDrfCsv:          return "dsi-kind-makedrf-csv";
      case FileKind::MultiDrfCsv:         return "dsi-kind-multi-csv";
      case FileKind::GammaQuantCsv:       return "dsi-kind-gq-csv";
      case FileKind::IsocsEcc:            return "dsi-kind-ecc";
      case FileKind::Angle:               return "dsi-kind-angle";
      case FileKind::EfficiencyCsv:       return "dsi-kind-eff-csv";
      case FileKind::GadrasEfficiencyCsv: return "dsi-kind-gadras-eff";
      case FileKind::GadrasDetectorDat:   return "dsi-kind-gadras-dat";
      case FileKind::ParGrid:             return "dsi-kind-par";
      case FileKind::ParDetectorTxt:      return "dsi-kind-detector-txt";
    }//switch( kind )

    return "dsi-kind-drf-xml";
  }//kind_key(...)


  const char *interpretation_key( const Interpretation interp )
  {
    switch( interp )
    {
      case Interpretation::AsIs:              return "dsi-interp-as-is";
      case Interpretation::FarFieldIntrinsic: return "ds-intrinsic-eff";
      case Interpretation::FarFieldAbsolute:  return "ds-abs-eff";
      case Interpretation::FixedTotal:        return "ds-fixed-geom-total";
      case Interpretation::FixedPerCm2:       return "ds-fixed-geom-cm2";
      case Interpretation::FixedPerM2:        return "ds-fixed-geom-m2";
      case Interpretation::FixedPerGram:      return "ds-fixed-geom-gram";
      case Interpretation::GenericDetector:   return "dsi-interp-generic";
    }//switch( interp )

    return "dsi-interp-as-is";
  }//interpretation_key(...)


  // The key describing the companion file `kind` wants; nullptr if it has none.
  const char *companion_key( const FileKind kind )
  {
    switch( kind )
    {
      case FileKind::GadrasEfficiencyCsv: return "dsi-companion-dat";
      case FileKind::GadrasDetectorDat:   return "dsi-companion-eff";
      case FileKind::ParGrid:             return "dsi-companion-txt";
      case FileKind::ParDetectorTxt:      return "dsi-companion-par";
      default:                            break;
    }//switch( kind )

    return nullptr;
  }//companion_key(...)


  // A distance from a line edit (PhysicalUnits), or `fallback` if empty or not a distance.
  double edit_distance( const WLineEdit *edit, const double fallback )
  {
    const string txt = SpecUtils::trim_copy( edit->text().toUTF8() );
    if( txt.empty() )
      return fallback;

    try
    {
      return PhysicalUnits::stringToDistance( txt );
    }catch( std::exception & )
    {
    }

    return fallback;
  }//edit_distance(...)


  WLineEdit *add_distance_row( WContainerWidget *parent, const WString &label,
                               shared_ptr<WValidator> validator, WContainerWidget *&row )
  {
    row = parent->addNew<WContainerWidget>();
    row->addStyleClass( "DrfImportRow" );
    WLabel *lbl = row->addNew<WLabel>( label );
    WLineEdit *edit = row->addNew<WLineEdit>();
    lbl->setBuddy( edit );
    edit->setValidator( validator );
    edit->setAttributeValue( "ondragstart", "return false" );
    edit->setAttributeValue( "autocorrect", "off" );
    edit->setAttributeValue( "spellcheck", "off" );
    edit->setPlaceholderText( "cm" );
    edit->setTextSize( 10 );
    return edit;
  }//add_distance_row(...)
}//namespace


DrfImportWidget::DrfImportWidget( const Host host )
  : WContainerWidget(),
    m_host( host ),
    m_mainUpload( nullptr ),
    m_companionUpload( nullptr ),
    m_mainDrop( nullptr ),
    m_fileList( nullptr ),
    m_companionDrop( nullptr ),
    m_companionTxt( nullptr ),
    m_recordDiv( nullptr ),
    m_recordCombo( nullptr ),
    m_interpDiv( nullptr ),
    m_interpCombo( nullptr ),
    m_interpretations{},
    m_diameterDiv( nullptr ),
    m_diameterEdit( nullptr ),
    m_setbackDiv( nullptr ),
    m_setbackEdit( nullptr ),
    m_distanceDiv( nullptr ),
    m_distanceEdit( nullptr ),
    m_eccUncertHolder( nullptr ),
    m_eccUncert( nullptr ),
    m_nameDiv( nullptr ),
    m_nameEdit( nullptr ),
    m_status( nullptr ),
    m_notes( nullptr ),
    m_characterizeBtn( nullptr ),
    m_primary{},
    m_companion{},
    m_source{},
    m_generation( 0 ),
    m_candidate{},
    m_seed{},
    m_changed{},
    m_characterizeRequested{}
{
  InterSpec * const viewer = InterSpec::instance();
  if( viewer )
    viewer->useMessageResourceBundle( "DrfSelect" );

  WApplication * const app = WApplication::instance();
  app->useStyleSheet( "InterSpec_resources/DrfImportWidget.css" );
  app->require( "InterSpec_resources/BatchGuiWidget.js" );

  addStyleClass( "DrfImportWidget" );

  const bool showToolTips = viewer && UserPreferences::preferenceValue<bool>( "ShowTooltips", viewer );
  const bool isMobile = viewer && viewer->isMobile();

  // Drop area to choose the file; also click (or tap) to browse, where more than one file may be
  //  chosen.  The drop dialog already has its file.
  if( m_host == Host::ImportTab )
  {
    m_mainUpload = make_unique<FileDragUploadResource>();
    m_mainUpload->fileDrop().connect( this, [this]( const string &name, const string &spool ){
      handleUpload( name, spool, false );
    } );

    m_mainDrop = addNew<WContainerWidget>();
    m_mainDrop->addStyleClass( "DrfImportDrop DrfImportMainDrop" );
    m_mainDrop->addNew<WText>( WString::tr( isMobile ? "dsi-drop-txt-mobile" : "dsi-drop-txt" ) );
    m_mainDrop->doJavaScript( "BatchInputDropUploadSetup(" + m_mainDrop->jsRef() + ", '"
                              + m_mainUpload->url() + "', true);" );
    HelpSystem::attachToolTipOn( m_mainDrop, WString::tr("dsi-supported-tt"), showToolTips );
  }//if( m_host == Host::ImportTab )

  m_fileList = addNew<WContainerWidget>();
  m_fileList->addStyleClass( "DrfImportFiles" );

  // The other half of a Detector.dat/Efficiency.csv or .par/DETECTOR.txt pair.
  m_companionUpload = make_unique<FileDragUploadResource>();
  m_companionUpload->fileDrop().connect( this, [this]( const string &name, const string &spool ){
    handleUpload( name, spool, true );
  } );

  m_companionDrop = addNew<WContainerWidget>();
  m_companionDrop->addStyleClass( "DrfImportDrop DrfImportCompanion" );
  m_companionTxt = m_companionDrop->addNew<WText>();
  m_companionDrop->doJavaScript( "BatchInputDropUploadSetup(" + m_companionDrop->jsRef() + ", '"
                                 + m_companionUpload->url() + "');" );
  m_companionDrop->hide();

  string drop_ids = "'" + m_companionDrop->id() + "'";
  if( m_mainDrop )
    drop_ids += ",'" + m_mainDrop->id() + "'";
  doJavaScript( "setupOnDragEnterDom([" + drop_ids + "]);" );

  m_recordDiv = addNew<WContainerWidget>();
  m_recordDiv->addStyleClass( "DrfImportRow" );
  WLabel *label = m_recordDiv->addNew<WLabel>( WString::tr("dsi-record-label") );
  m_recordCombo = m_recordDiv->addNew<WComboBox>();
  label->setBuddy( m_recordCombo );
  m_recordCombo->activated().connect( this, &DrfImportWidget::handleRecordChanged );
  m_recordDiv->hide();

  m_interpDiv = addNew<WContainerWidget>();
  m_interpDiv->addStyleClass( "DrfImportRow" );
  label = m_interpDiv->addNew<WLabel>( WString::tr("dsi-interp-label") );
  m_interpCombo = m_interpDiv->addNew<WComboBox>();
  label->setBuddy( m_interpCombo );
  m_interpCombo->activated().connect( this, &DrfImportWidget::updateCandidate );
  HelpSystem::attachToolTipOn( m_interpCombo, WString::tr("ds-tt-manual-drf-type"), showToolTips );
  m_interpDiv->hide();

  auto distValidator = make_shared<WRegExpValidator>( PhysicalUnits::sm_distanceRegex );
  distValidator->setFlags( Wt::RegExpFlag::MatchCaseInsensitive );

  m_diameterEdit = add_distance_row( this, WString::tr("ds-det-diam"), distValidator, m_diameterDiv );
  m_setbackEdit = add_distance_row( this, WString::tr("ds-det-setback"), distValidator, m_setbackDiv );
  m_distanceEdit = add_distance_row( this, WString::tr("ds-dist-label"), distValidator, m_distanceDiv );
  m_diameterDiv->hide();
  m_setbackDiv->hide();
  m_distanceDiv->hide();

  m_eccUncertHolder = addNew<WContainerWidget>();
  m_eccUncertHolder->hide();

  m_nameDiv = addNew<WContainerWidget>();
  m_nameDiv->addStyleClass( "DrfImportRow" );
  label = m_nameDiv->addNew<WLabel>( WString::tr("ds-name-label") );
  m_nameEdit = m_nameDiv->addNew<WLineEdit>();
  label->setBuddy( m_nameEdit );
  m_nameEdit->setAttributeValue( "ondragstart", "return false" );
  m_nameEdit->setAttributeValue( "autocorrect", "off" );
  m_nameEdit->setAttributeValue( "spellcheck", "off" );
  m_nameEdit->setTextSize( 30 );
  m_nameDiv->hide();

  for( WLineEdit *edit : { m_diameterEdit, m_setbackEdit, m_distanceEdit, m_nameEdit } )
  {
    edit->changed().connect( this, &DrfImportWidget::updateCandidate );
    edit->enterPressed().connect( this, &DrfImportWidget::updateCandidate );
    if( m_host == Host::DropDialog )
      edit->textInput().connect( this, &DrfImportWidget::updateCandidate );
  }

  m_characterizeBtn = addNew<WPushButton>( WString::tr("dsi-characterize") );
  m_characterizeBtn->addStyleClass( "DrfImportCharacterize" );
  m_characterizeBtn->clicked().connect( this, [this](){ m_characterizeRequested.emit(); } );
  m_characterizeBtn->hide();

  m_status = addNew<WText>();
  m_status->addStyleClass( "DrfImportStatus" );
  m_status->setInline( false );

  m_notes = addNew<WContainerWidget>();
  m_notes->addStyleClass( "DrfImportNotes" );
}//DrfImportWidget constructor


DrfImportWidget::~DrfImportWidget()
{
  WApplication * const app = WApplication::instance();
  if( !app || !m_companionDrop )
    return;

  string drop_ids = "'" + m_companionDrop->id() + "'";
  if( m_mainDrop )
    drop_ids += ",'" + m_mainDrop->id() + "'";
  app->doJavaScript( "if(window.removeOnDragEnterDom)removeOnDragEnterDom([" + drop_ids + "]);" );
}//~DrfImportWidget()


Wt::Signal<> &DrfImportWidget::changed()
{
  return m_changed;
}


Wt::Signal<> &DrfImportWidget::characterizeRequested()
{
  return m_characterizeRequested;
}


std::shared_ptr<DetectorPeakResponse> DrfImportWidget::candidate() const
{
  return m_candidate;
}


std::shared_ptr<DetectorPeakResponse> DrfImportWidget::characterizationSeed() const
{
  return m_seed;
}


void DrfImportWidget::setStatus( const Wt::WString &txt, const bool isError )
{
  m_status->setText( txt );
  m_status->toggleStyleClass( "DrfImportError", isError );
}//setStatus(...)


void DrfImportWidget::handleUpload( const std::string &displayName,
                                    const std::string &spoolName,
                                    const bool companionArea )
{
  FileDragUploadResource * const resource = companionArea ? m_companionUpload.get()
                                                          : m_mainUpload.get();

  shared_ptr<const DrfImport::ParsedFile> parsed;
  try
  {
    const shared_ptr<const string> data = DrfImport::readFile( spoolName );
    parsed = DrfImport::parseFile( SpecUtils::filename( displayName ), data, true );
  }catch( std::exception &e )
  {
    setStatus( WString::tr("dsi-parse-error")
                 .arg( Wt::Utils::htmlEncode( SpecUtils::filename(displayName) ) )
                 .arg( Wt::Utils::htmlEncode( string( e.what() ) ) ), true );
  }//try / catch

  if( resource )
    resource->clearSpooledFiles();

  if( parsed )
  {
    // The companion area only takes the other half of the current file.
    if( companionArea
        && (!m_primary || !DrfImport::areCompanions( m_primary->kind, parsed->kind )) )
    {
      setStatus( WString::tr("dsi-companion-wrong")
                   .arg( Wt::Utils::htmlEncode( parsed->displayName ) )
                   .arg( WString::tr( kind_key( parsed->kind ) ) ), true );
    }else
    {
      addFile( parsed );
    }
  }//if( parsed )

  WApplication * const app = WApplication::instance();
  if( app )
    app->triggerUpdate();
}//handleUpload(...)


bool DrfImportWidget::addFile( const std::string &displayName, std::shared_ptr<const std::string> data )
{
  try
  {
    addFile( DrfImport::parseFile( SpecUtils::filename( displayName ), data, true ) );
    return true;
  }catch( std::exception &e )
  {
    setStatus( WString::tr("dsi-parse-error")
                 .arg( Wt::Utils::htmlEncode( SpecUtils::filename(displayName) ) )
                 .arg( Wt::Utils::htmlEncode( string( e.what() ) ) ), true );
  }

  return false;
}//bool addFile(...)


void DrfImportWidget::addFile( std::shared_ptr<const DrfImport::ParsedFile> file )
{
  if( !file )
    return;

  // A file is one half of the current pair, replacing whichever half it is (the two files of a
  //  pair arrive one at a time, in either order); anything else starts a new import.
  if( m_primary && DrfImport::areCompanions( m_primary->kind, file->kind ) )
  {
    m_companion = file;
  }else if( m_primary && m_companion && (m_primary->kind == file->kind) )
  {
    m_primary = file;
  }else
  {
    m_primary = file;
    m_companion.reset();
  }

  startSource( 0 );
}//void addFile( std::shared_ptr<const DrfImport::ParsedFile> file )


void DrfImportWidget::updateFileList()
{
  m_fileList->clear();

  for( const shared_ptr<const DrfImport::ParsedFile> &file : { m_primary, m_companion } )
  {
    if( !file )
      continue;

    WText *txt = m_fileList->addNew<WText>( WString::tr("dsi-file-line")
                                             .arg( Wt::Utils::htmlEncode( file->displayName ) )
                                             .arg( WString::tr( kind_key( file->kind ) ) ) );
    txt->setInline( false );
  }//for( both files )
}//void updateFileList()


void DrfImportWidget::startSource( const size_t record )
{
  ++m_generation;

  updateFileList();

  const shared_ptr<const DrfImport::ParsedFile> primary = m_primary, companion = m_companion;

  if( !DrfImport::isSlow( primary.get(), companion.get() ) )
  {
    setSource( make_shared<const DrfImport::Source>(
                                          DrfImport::makeSource( primary, companion, record ) ) );
    return;
  }

  // Building from a .par grid takes seconds, so do it off the session thread; show the pair as
  //  pending in the meantime.
  m_source.reset();
  m_candidate.reset();
  m_seed.reset();
  const std::initializer_list<WWidget *> options{ m_companionDrop, m_recordDiv, m_interpDiv,
    m_diameterDiv, m_setbackDiv, m_distanceDiv, m_eccUncertHolder, m_nameDiv, m_characterizeBtn };
  for( WWidget *w : options )
    w->hide();
  m_notes->clear();
  setStatus( WString::tr("dsi-building"), false );
  m_changed.emit();

  const uint64_t generation = m_generation;
  const WidgetUtils::WidgetHandle self( this );
  const auto result = make_shared<DrfImport::Source>();

  AppUtils::run_on_ioservice(
    [primary, companion, record, result](){
      *result = DrfImport::makeSource( primary, companion, record );
    },
    [self, generation, result](){
      DrfImportWidget * const widget = self.resolve_as<DrfImportWidget>();
      if( !widget || (widget->m_generation != generation) )
        return;  //the user has since dropped something else

      widget->setSource( result );

      WApplication * const app = WApplication::instance();
      if( app )
        app->triggerUpdate();
    } );
}//void startSource()


void DrfImportWidget::handleRecordChanged()
{
  startSource( static_cast<size_t>( std::max( m_recordCombo->currentIndex(), 0 ) ) );
}//void handleRecordChanged()


void DrfImportWidget::setSource( std::shared_ptr<const DrfImport::Source> src )
{
  m_source = src;
  if( !src )
    return;

  // Companion file: what it adds, and whether it is required
  const char * const companion = (m_primary && !m_companion) ? companion_key( m_primary->kind ) : nullptr;
  m_companionDrop->setHidden( !companion );
  if( companion )
    m_companionTxt->setText( WString::tr( companion ) );

  // Several DRFs in one file, or several DETECTOR.txt records none of which name the .par
  m_recordCombo->clear();
  for( const string &name : src->recordNames )
    m_recordCombo->addItem( WString::fromUTF8( name ) );
  m_recordDiv->setHidden( src->recordNames.size() < 2 );
  if( src->recordNames.size() > 1 )
    m_recordCombo->setCurrentIndex( static_cast<int>( src->record ) );

  m_interpretations = src->interpretations;
  m_interpCombo->clear();
  for( const Interpretation interp : m_interpretations )
    m_interpCombo->addItem( WString::tr( interpretation_key( interp ) ) );
  const auto default_pos = std::find( begin(m_interpretations), end(m_interpretations),
                                      src->defaultInterpretation );
  if( default_pos != end(m_interpretations) )
    m_interpCombo->setCurrentIndex( static_cast<int>( default_pos - begin(m_interpretations) ) );
  m_interpDiv->setHidden( m_interpretations.empty() );

  // The file's own diameter/setback, when it has one; otherwise keep what the user entered
  const shared_ptr<const DetectorPeakResponse> &base = src->base;
  if( base && (base->detectorDiameter() > 0.0f) )
    m_diameterEdit->setText( PhysicalUnits::printToBestLengthUnits( base->detectorDiameter() ) );
  if( base && (base->detectorSetback() > 0.0) )
    m_setbackEdit->setText( PhysicalUnits::printToBestLengthUnits( base->detectorSetback() ) );

  m_eccUncertHolder->clear();
  m_eccUncert = nullptr;
  if( src->uncertEnergies.size() >= 2 )
  {
    m_eccUncert = m_eccUncertHolder->addNew<EccUncertOptions>( src->uncertEnergies,
                                                src->baselineFrac, src->convergenceFrac );
    m_eccUncert->changed().connect( this, &DrfImportWidget::updateCandidate );
  }
  m_eccUncertHolder->setHidden( !m_eccUncert );

  // Common file names (Efficiency.csv, Detector.dat) say little, so date-stamp short ones
  string name = src->name;
  if( src->nameIsFileStem && (name.size() < 15) )
  {
    auto now = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
    WApplication * const app = WApplication::instance();
    if( app )
      now += app->environment().timeZoneOffset();
    name += " " + SpecUtils::to_vax_string( now );
  }
  m_nameEdit->setText( WString::fromUTF8( name ) );
  m_nameDiv->setHidden( (src->status != Status::Ready)
                        && (src->status != Status::NeedsCharacterization) );

  m_characterizeBtn->setHidden( src->status != Status::NeedsCharacterization );

  m_notes->clear();
  for( const string &note : src->notes )
  {
    WText *txt = m_notes->addNew<WText>( WString::fromUTF8( note ), Wt::TextFormat::Plain );
    txt->setInline( false );
  }

  // Credits are HTML by design (as on the "Rel. Eff." tab); WText filters XHTML of anything unsafe.
  if( !src->credits.empty() )
  {
    string html;
    for( const string &credit : src->credits )
      html += "<div>" + credit + "</div>";
    WText *txt = m_notes->addNew<WText>( WString::fromUTF8( html ), Wt::TextFormat::XHTML );
    txt->setInline( false );
  }

  updateCandidate();
}//void setSource(...)


DrfImport::Options DrfImportWidget::currentOptions() const
{
  DrfImport::Options options;

  const int index = m_interpCombo->currentIndex();
  if( (index >= 0) && (index < static_cast<int>( m_interpretations.size() )) )
    options.interpretation = m_interpretations[index];

  options.diameter = edit_distance( m_diameterEdit, 0.0 );
  options.setback = edit_distance( m_setbackEdit, -1.0 );
  options.distance = edit_distance( m_distanceEdit, 0.0 );

  if( m_eccUncert )
  {
    options.setUncert = true;
    options.uncert = m_eccUncert->buildUncert();
  }

  options.name = SpecUtils::trim_copy( m_nameEdit->text().toUTF8() );

  return options;
}//DrfImport::Options currentOptions() const


void DrfImportWidget::updateCandidate()
{
  m_candidate.reset();
  m_seed.reset();

  if( !m_source )
  {
    m_changed.emit();
    return;
  }

  const DrfImport::Options options = currentOptions();

  const bool far_field = (!m_interpretations.empty()
                          && ((options.interpretation == Interpretation::FarFieldIntrinsic)
                              || (options.interpretation == Interpretation::FarFieldAbsolute)));
  m_diameterDiv->setHidden( !far_field );
  m_setbackDiv->setHidden( !far_field );
  m_distanceDiv->setHidden( !far_field || (options.interpretation != Interpretation::FarFieldAbsolute) );

  const DrfImport::Result result = DrfImport::build( *m_source, options );

  switch( result.status )
  {
    case Status::Ready:
      m_candidate = result.drf;
      setStatus( WString(), false );
      break;

    case Status::NeedsCompanion:
      setStatus( WString(), false );  //the companion drop area says what is needed
      break;

    case Status::NeedsCharacterization:
      if( m_source->base )
      {
        m_seed = make_shared<DetectorPeakResponse>( *m_source->base );
        if( !options.name.empty() )
          m_seed->setName( options.name );
      }
      setStatus( WString::tr("dsi-geometry-only"), false );
      break;

    case Status::NeedsDiameter:
      setStatus( WString::tr("dsi-need-diameter"), false );
      break;

    case Status::NeedsDistance:
      setStatus( WString::tr("dsi-need-distance"), false );
      break;

    case Status::Error:
      setStatus( WString::tr("dsi-build-error").arg( Wt::Utils::htmlEncode( result.error ) ), true );
      break;
  }//switch( result.status )

  m_changed.emit();
}//void updateCandidate()


DrfImportWidget *DrfImportWidget::setupDropDialog( SimpleDialog *dialog,
                                     std::shared_ptr<const DrfImport::ParsedFile> file )
{
  InterSpec * const viewer = InterSpec::instance();
  assert( dialog && viewer );
  if( !dialog || !viewer )
    return nullptr;

  viewer->useMessageResourceBundle( "DrfSelect" );

  dialog->addStyleClass( "DrfImportDialog" );

  WContainerWidget * const contents = dialog->contents();

  WText *title = contents->addNew<WText>( WString::tr("dsi-dialog-title") );
  title->addStyleClass( "title" );
  title->setInline( false );

  DrfChart *chart = contents->addNew<DrfChart>();
  chart->addStyleClass( "DrfImportChart" );
  chart->setMinimumSize( 300, 175 );
  int chartw = 350, charth = 200;
  if( viewer->renderedWidth() > 500 )
    chartw = std::min( ((3 * viewer->renderedWidth() / 4) - 50), 500 );
  if( viewer->renderedHeight() > 400 )
    charth = std::min( viewer->renderedHeight() / 4, (4 * chartw) / 7 );
  chart->resize( std::max( chartw, 300 ), std::max( charth, 175 ) );

  DrfImportWidget *widget = contents->addNew<DrfImportWidget>( DrfImportWidget::Host::DropDialog );

  const shared_ptr<SpecMeas> fore = viewer->measurment( SpecUtils::SpectrumType::Foreground );
  const shared_ptr<DetectorPeakResponse> prev = fore ? fore->detector() : nullptr;

  WCheckBox *defaultForSerialNumber = nullptr, *defaultForDetectorModel = nullptr;
  if( fore && !fore->instrument_id().empty() )
  {
    defaultForSerialNumber = contents->addNew<WCheckBox>(
                              WString::tr("ds-use-for-serialnum-cb").arg( fore->instrument_id() ) );
    defaultForSerialNumber->addStyleClass( "CbNoLineBreak DrfImportDefaultCb" );
    defaultForSerialNumber->setInline( false );
  }

  if( fore && ((fore->detector_type() != SpecUtils::DetectorType::Unknown)
               || !fore->instrument_model().empty()) )
  {
    const string model = (fore->detector_type() != SpecUtils::DetectorType::Unknown)
                           ? detectorTypeToString( fore->detector_type() )
                           : fore->instrument_model();
    defaultForDetectorModel = contents->addNew<WCheckBox>(
                                            WString::tr("ds-use-for-model-cb").arg( model ) );
    defaultForDetectorModel->addStyleClass( "CbNoLineBreak DrfImportDefaultCb" );
    defaultForDetectorModel->setInline( false );
  }

  dialog->addButton( WString::tr("Cancel"), WidgetUtils::ButtonRole::Dismiss );
  WPushButton * const further = dialog->addButton( WString::tr("dsi-further-options"),
                                                   WidgetUtils::ButtonRole::Neutral );
  WPushButton * const accept = dialog->addButton( WString::tr("dsi-use-drf"),
                                                  WidgetUtils::ButtonRole::Affirm );

  widget->changed().connect( widget, [widget, chart, accept, further](){
    const shared_ptr<DetectorPeakResponse> drf = widget->candidate();
    const shared_ptr<const DetectorPeakResponse> shown = drf ? drf : widget->characterizationSeed();
    chart->updateChart( (shown && shown->isValid()) ? shown : nullptr );
    if( shown && (shown->upperEnergy() > 10000.0) )
      chart->setXAxisRange( std::max( shown->lowerEnergy(), 0.0 ), 4000.0 );
    accept->setEnabled( !!drf );
    further->setEnabled( !!drf );
  } );

  // Geometry only: Monte Carlo characterization happens in the Modify tool, which applies the
  //  result itself, so this dialog is done.
  widget->characterizeRequested().connect( widget, [widget, dialog](){
    InterSpec * const viewer = InterSpec::instance();
    const shared_ptr<DetectorPeakResponse> seed = widget->characterizationSeed();
    if( viewer && seed && viewer->showDrfModifyWindow( seed ) )
      dialog->done( Wt::DialogCode::Accepted );
  } );

  // Shared by "Use DRF" and "Further options..."; the latter also opens the editor on the result.
  auto accept_drf = [widget, fore, prev, defaultForSerialNumber, defaultForDetectorModel]( const bool open_modify ){
    InterSpec * const viewer = InterSpec::instance();
    const shared_ptr<DetectorPeakResponse> drf = widget->candidate();
    if( !viewer || !drf )
      return;

    shared_ptr<DataBaseUtils::DbSession> sql = viewer->sql();
    const Wt::Dbo::ptr<InterSpecUser> &user = viewer->user();
    DrfSelect::updateLastUsedTimeOrAddToDb( drf, user.id(), sql );
    viewer->detectorChanged().emit( drf );

    for( const pair<WCheckBox *,UseDrfPref::UseDrfType> &cb
          : { make_pair( defaultForSerialNumber, UseDrfPref::UseDrfType::UseDetectorSerialNumber ),
              make_pair( defaultForDetectorModel, UseDrfPref::UseDrfType::UseDetectorModelName ) } )
    {
      if( !fore || !cb.first || !cb.first->isChecked() )
        continue;

      const UseDrfPref::UseDrfType preftype = cb.second;
      WServer::instance()->ioService().boost::asio::io_service::post( std::bind( [=](){
        DrfSelect::setUserPrefferedDetector( drf, sql, user, preftype, fore );
      } ) );
    }//for( both "use as default" checkboxes )

    UndoRedoManager * const undoManager = viewer->undoRedoManager();
    if( undoManager && undoManager->canAddUndoRedoNow() )
    {
      auto undo = [prev](){
        InterSpec * const viewer = InterSpec::instance();
        if( viewer )
          viewer->detectorChanged().emit( prev );
      };

      auto redo = [drf](){
        InterSpec * const viewer = InterSpec::instance();
        if( viewer )
          viewer->detectorChanged().emit( drf );
      };

      undoManager->addUndoRedoStep( undo, redo, "Import DRF from file" );
    }//if( undoManager )

    if( open_modify )
      viewer->showDrfModifyWindow( drf );
  };//accept_drf

  accept->clicked().connect( widget, [accept_drf](){ accept_drf( false ); } );
  further->clicked().connect( widget, [accept_drf](){ accept_drf( true ); } );

  accept->disable();
  further->disable();

  widget->addFile( file );

  return widget;
}//DrfImportWidget *setupDropDialog(...)
