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
#include <ctime>
#include <chrono>
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>

#include <Wt/Utils.h>
#include <Wt/WDate.h>
#include <Wt/WLink.h>
#include <Wt/WText.h>
#include <Wt/WAnchor.h>
#include <Wt/WCheckBox.h>
#include <Wt/WResource.h>
#include <Wt/WEnvironment.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>
#include <Wt/Http/Request.h>
#include <Wt/Http/Response.h>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/DrfChart.h"
#include "InterSpec/DrfAppearance.h"
#include "InterSpec/DrfLibraryWidgets.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace Wt;

using DrfLibrary::Row;
using DrfLibrary::Entry;
using DrfLibrary::Group;
using DrfLibrary::Range;
using DrfLibrary::Origin;
using DrfLibrary::RowType;
using DrfLibrary::Summary;


namespace
{
  string html( const string &text )
  {
    return Wt::Utils::htmlEncode( text );
  }

  string html( const WString &text )
  {
    return Wt::Utils::htmlEncode( text.toUTF8() );
  }

  /** A length in mm, to the precision a crystal dimension is quoted with. */
  string format_mm( const double cm )
  {
    const double mm = 10.0 * cm;
    char buffer[32];
    const bool whole = (mm >= 100.0) || (std::fabs( mm - std::round(mm) ) < 0.05);
    snprintf( buffer, sizeof(buffer), whole ? "%.0f" : "%.1f", mm );
    return buffer;
  }

  /** The DRF XML of whatever a DrfPreviewPanel shows, for the "XML" button. */
  class DrfXmlResource : public Wt::WResource
  {
    DrfPreviewPanel * const m_panel;
    Wt::WApplication * const m_app;

  public:
    explicit DrfXmlResource( DrfPreviewPanel *panel )
      : WResource(), m_panel( panel ), m_app( WApplication::instance() )
    {
      // Hold the session lock across handleRequest (see DrfSelect.cpp for why: otherwise the
      //  session thread can wait in ~WResource on us while we wait for the lock it holds).
      setTakesUpdateLock( true );
    }

    virtual ~DrfXmlResource()
    {
      beingDeleted();
    }

    virtual void handleRequest( const Wt::Http::Request &, Wt::Http::Response &response )
    {
      WApplication::UpdateLock lock( m_app );
      if( !lock || !m_panel )
      {
        response.setStatus( 500 );
        return;
      }

      const shared_ptr<const DetectorPeakResponse> det = m_panel->detector();
      if( !det || !det->isValid() )
      {
        response.setStatus( 404 );
        return;
      }

      string name = det->name();
      name += (name.empty() ? "" : ".");
      name += "drf.xml";
      for( char &c : name )
      {
        if( string("\\/:?\"<>|*").find( c ) != string::npos )
          c = '_';
      }

      suggestFileName( name, ContentDisposition::Attachment );
      response.setMimeType( "application/octet-stream" );

      try
      {
        rapidxml::xml_document<char> doc;
        det->toXml( &doc, &doc );
        rapidxml::print( response.out(), doc, 0 );
      }catch( std::exception &e )
      {
        log("error") << "Failed to write DRF XML: " << e.what();
        response.setStatus( 500 );
      }
    }//handleRequest(...)
  };//class DrfXmlResource


  /** One row of the "What it includes" checklist. */
  struct Capability
  {
    WString label;
    bool ok = false;
    WString detail;
  };

  WString efficiency_detail( const Summary &s )
  {
    using DrfSource = DetectorPeakResponse::DrfSource;

    if( !s.hasEfficiency )
      return WString::tr( "dsl-eff-needs-mc" );
    if( s.importedGrid )
      return WString::tr( "dsl-eff-par" );

    switch( s.source )
    {
      case DrfSource::IsocsEcc:
        return WString::tr( "dsl-eff-ecc" );
      case DrfSource::AngleOutx:
        return WString::tr( "dsl-eff-angle" );
      case DrfSource::DefaultGadrasDrf:
      case DrfSource::UserAddedGadrasDrf:
        return WString::tr( "dsl-eff-gadras" );
      case DrfSource::UserSpecifiedFormulaDrf:
        return WString::tr( "dsl-eff-formula" );
      default:
        break;
    }//switch( s.source )

    if( s.hasResponse )
      return WString::tr( "dsl-eff-mc" );
    if( s.measuredPointCount > 0 )
      return WString::tr( "dsl-eff-fit-points" ).arg( static_cast<int>(s.measuredPointCount) );

    switch( s.efficiencyForm )
    {
      case DetectorPeakResponse::kEnergyEfficiencyPairs:
        return WString::tr( "dsl-eff-table" ).arg( static_cast<int>(s.efficiencyParamCount) );
      case DetectorPeakResponse::kFunctialEfficienyForm:
        return WString::tr( "dsl-eff-formula" );
      case DetectorPeakResponse::kExpOfLogPowerSeries:
        return WString::tr( "dsl-eff-coefs" ).arg( static_cast<int>(s.efficiencyParamCount) );
      case DetectorPeakResponse::kNumEfficiencyFnctForms:
        break;
    }

    return WString::tr( "dsl-eff-included" );
  }//efficiency_detail(...)


  vector<Capability> capabilities( const Summary &s )
  {
    vector<Capability> caps( 6 );

    caps[0].label = WString::tr( "dsl-cap-eff" );
    caps[0].ok = s.hasEfficiency;
    caps[0].detail = efficiency_detail( s );

    caps[1].label = WString::tr( "dsl-cap-fwhm" );
    caps[1].ok = s.hasFwhm && (s.fwhmAt662 > 0.0f);
    if( caps[1].ok )
    {
      char buffer[32];
      snprintf( buffer, sizeof(buffer), (s.fwhmAt662 < 10.0f) ? "%.2f" : "%.0f", s.fwhmAt662 );
      caps[1].detail = WString::tr( "dsl-cap-fwhm-yes" ).arg( buffer );
    }else
    {
      caps[1].detail = WString::tr( "dsl-cap-not-included" );
    }

    caps[2].label = WString::tr( "dsl-cap-geom" );
    caps[2].ok = s.hasGeometry;
    caps[2].detail = WString::tr( s.hasGeometry ? "dsl-cap-geom-yes" : "dsl-cap-geom-no" );

    caps[3].label = WString::tr( "dsl-cap-near" );
    caps[3].ok = (s.range == Range::NearFieldOk);
    switch( s.range )
    {
      case Range::NearFieldOk:
        caps[3].detail = WString::tr( "dsl-cap-near-yes" );
        break;
      case Range::FixedGeometry:
        caps[3].detail = WString::tr( "dsl-cap-near-fixed" );
        break;
      case Range::FarFieldOnly:
      case Range::NeedsCharacterization:
        caps[3].detail = WString::tr( s.hasGeometry ? "dsl-cap-near-mc" : "dsl-cap-near-no" );
        break;
    }//switch( s.range )

    caps[4].label = WString::tr( "dsl-cap-unc" );
    caps[4].ok = s.hasUncertainty && !s.uncertaintyIsAssumed;
    if( caps[4].ok )
      caps[4].detail = WString::tr( "dsl-cap-unc-yes" );
    else if( s.hasUncertainty )
      caps[4].detail = WString::tr( "dsl-cap-unc-assumed" );
    else
      caps[4].detail = WString::tr( "dsl-cap-not-included" );

    caps[5].label = WString::tr( "dsl-cap-total" );
    caps[5].ok = s.hasTotalEfficiency;
    caps[5].detail = WString::tr( s.hasTotalEfficiency ? "dsl-cap-total-yes" : "dsl-cap-not-included" );

    return caps;
  }//capabilities(...)


  /** Sets a segmented-control button's pressed look. */
  void set_pressed( WPushButton *btn, const bool pressed )
  {
    btn->toggleStyleClass( "Pressed", pressed );
    btn->setAttributeValue( "aria-pressed", pressed ? "true" : "false" );
  }
}//namespace


namespace DrfLibraryText
{
WString originLabel( const Origin origin )
{
  switch( origin )
  {
    case Origin::BuiltIn:      return WString::tr( "dsl-origin-builtin" );
    case Origin::Gadras:       return WString::tr( "dsl-origin-gadras" );
    case Origin::Isocs:        return WString::tr( "dsl-origin-isocs" );
    case Origin::Angle:        return WString::tr( "dsl-origin-angle" );
    case Origin::RelEffTable:  return WString::tr( "dsl-origin-reltable" );
    case Origin::UserMade:     return WString::tr( "dsl-origin-usermade" );
    case Origin::ThisSpectrum: return WString::tr( "dsl-origin-spectrum" );
  }
  return WString();
}//originLabel(...)


std::string crystalLabel( const Summary &summary )
{
  const DrfLibrary::Crystal &c = summary.crystal;

  string dims;
  if( c.fromGeometry && (c.diameter_cm > 0.0) && (c.length_cm > 0.0) )
  {
    if( c.isBox && (c.height_cm > 0.0) )
      dims = format_mm( c.diameter_cm ) + " × " + format_mm( c.height_cm ) + " × " + format_mm( c.length_cm ) + " mm";
    else
      dims = format_mm( c.diameter_cm ) + " × " + format_mm( c.length_cm ) + " mm";
  }else if( c.faceDiameter_cm > 0.0 )
  {
    char buffer[32];
    snprintf( buffer, sizeof(buffer), "Ø %.1f cm", c.faceDiameter_cm );
    dims = buffer;
  }

  if( c.material.empty() || dims.empty() )
    return c.material + dims;
  return c.material + " · " + dims;
}//crystalLabel(...)


WString rangeLabel( const Range range )
{
  switch( range )
  {
    case Range::FarFieldOnly:          return WString::tr( "dsl-range-far" );
    case Range::NearFieldOk:           return WString::tr( "dsl-range-near" );
    case Range::FixedGeometry:         return WString::tr( "dsl-range-fixed" );
    case Range::NeedsCharacterization: return WString::tr( "dsl-range-needs-mc" );
  }
  return WString();
}//rangeLabel(...)


WString relEffLabel( const Summary &summary )
{
  if( (summary.relEffPercent <= 0.0) || summary.nameStatesRelEff )
    return WString();

  char buffer[32];
  snprintf( buffer, sizeof(buffer), (summary.relEffPercent >= 10.0) ? "%.0f" : "%.1f", summary.relEffPercent );
  return WString::tr( "dsl-rel-eff" ).arg( buffer );
}//relEffLabel(...)


WString metaLine( const Summary &summary, const Origin origin )
{
  vector<WString> parts;
  const string crystal = crystalLabel( summary );
  if( !crystal.empty() )
    parts.push_back( WString::fromUTF8( crystal ) );
  const WString rel_eff = relEffLabel( summary );
  if( !rel_eff.empty() )
    parts.push_back( rel_eff );
  parts.push_back( originLabel( origin ) );
  parts.push_back( rangeLabel( summary.range ) );

  WString answer;
  for( size_t i = 0; i < parts.size(); ++i )
  {
    if( i )
      answer += WString::fromUTF8( " · " );
    answer += parts[i];
  }
  return answer;
}//metaLine(...)


WString lastUsedLabel( const int64_t utc )
{
  if( utc <= 0 )
    return WString::fromUTF8( "—" );

  // Calendar days in the user's time zone
  std::chrono::minutes offset( 0 );
  const WApplication * const app = WApplication::instance();
  if( app )
    offset = app->environment().timeZoneOffset();

  const int64_t now_utc = static_cast<int64_t>( std::time( nullptr ) );
  const int64_t local_then = utc + 60 * static_cast<int64_t>( offset.count() );
  const int64_t local_now = now_utc + 60 * static_cast<int64_t>( offset.count() );
  const int64_t day_then = (local_then >= 0) ? (local_then / 86400) : ((local_then - 86399) / 86400);
  const int64_t day_now = (local_now >= 0) ? (local_now / 86400) : ((local_now - 86399) / 86400);

  if( day_then == day_now )
    return WString::tr( "dsl-used-today" );
  if( day_then == (day_now - 1) )
    return WString::tr( "dsl-used-yesterday" );

  // 1970-01-01 is day zero
  const WDate then = WDate( 1970, 1, 1 ).addDays( static_cast<int>(day_then) );
  const WDate today = WDate( 1970, 1, 1 ).addDays( static_cast<int>(day_now) );
  WString answer = WString::fromUTF8( std::to_string( then.day() ) + " " );
  answer += WDate::shortMonthName( then.month() );
  if( then.year() != today.year() )
    answer += WString::fromUTF8( " " + std::to_string( then.year() ) );
  return answer;
}//lastUsedLabel(...)


std::string iconHtml( const Summary *summary, const double height_px, const double max_width_px )
{
  if( summary && summary->appearance )
  {
    const string own = summary->appearance->icon_html( height_px, max_width_px );
    if( !own.empty() )
      return own;
  }
  return "<span class=\"DrfIcon DrfIconGeneric\"></span>";
}//iconHtml(...)
}//namespace DrfLibraryText


// ---------------------------------------------------------------------------------------------
//  DrfLibraryList
// ---------------------------------------------------------------------------------------------

DrfLibraryList::DrfLibraryList()
  : WContainerWidget(),
    m_entries(),
    m_filter(),
    m_extraText( nullptr ),
    m_expanded(),
    m_pickedId(),
    m_inUseId(),
    m_suggestedNote(),
    m_checking( false ),
    m_narrow( false ),
    m_numDetectorRows( 0 ),
    m_rows( nullptr ),
    m_pickedRow( nullptr )
{
  addStyleClass( "DrfLibList" );

  WContainerWidget *head = addNew<WContainerWidget>();
  head->addStyleClass( "DrfLibListHead DrfLibGrid" );
  head->setAttributeValue( "aria-hidden", "true" );
  for( const char *key : { "dsl-col-detector", "dsl-col-crystal", "dsl-col-source", "dsl-col-used" } )
    head->addNew<WText>( WString::tr(key), TextFormat::Plain );

  m_rows = addNew<WContainerWidget>();
  m_rows->addStyleClass( "DrfLibListRows" );

  rebuild();
}//DrfLibraryList constructor


void DrfLibraryList::setEntries( std::shared_ptr<const std::vector<DrfLibrary::Entry>> entries )
{
  m_entries = entries;
}


void DrfLibraryList::setFilter( const DrfLibrary::Filter &filter )
{
  m_filter = filter;
}


void DrfLibraryList::setPicked( const std::string &entryId )
{
  if( entryId == m_pickedId )
    return;

  const string previous = m_pickedId;
  m_pickedId = entryId;

  // Picking a detector that has customized copies lists them (see DrfLibrary::buildRows), as does
  //  leaving one whose copies were listed only because it was picked.
  const auto now_row = m_rowButtons.find( entryId );
  const bool now_has_copies = m_parentsWithCopies.count( entryId );
  const bool prev_had_copies = m_parentsWithCopies.count( previous );
  if( now_has_copies || prev_had_copies || (!entryId.empty() && (now_row == end(m_rowButtons))) )
  {
    rebuild();
    return;
  }

  const auto prev_row = m_rowButtons.find( previous );
  if( prev_row != end(m_rowButtons) )
  {
    prev_row->second->removeStyleClass( "Picked" );
    prev_row->second->setAttributeValue( "aria-pressed", "false" );
  }

  m_pickedRow = nullptr;
  if( now_row != end(m_rowButtons) )
  {
    m_pickedRow = now_row->second;
    m_pickedRow->addStyleClass( "Picked" );
    m_pickedRow->setAttributeValue( "aria-pressed", "true" );
  }
}//setPicked(...)


void DrfLibraryList::setInUse( const std::string &entryId )
{
  m_inUseId = entryId;
}


void DrfLibraryList::setSuggestedNote( const Wt::WString &note )
{
  m_suggestedNote = note;
}


void DrfLibraryList::setChecking( const bool checking )
{
  m_checking = checking;
}


void DrfLibraryList::setExtraSearchText( DrfLibrary::ExtraSearchText fcn )
{
  m_extraText = fcn;
}


void DrfLibraryList::setNarrow( const bool narrow )
{
  m_narrow = narrow;
  toggleStyleClass( "Narrow", narrow );
}


void DrfLibraryList::scrollToPicked()
{
  // After this response's DOM changes (the rows may have just been re-created), from the app
  //  rather than this widget, whose own JavaScript can run before its new children exist.
  if( m_pickedRow && wApp )
    wApp->doJavaScript( "setTimeout(function(){const e=" + m_pickedRow->jsRef()
                        + ";if(e)e.scrollIntoView({block:'center'});},0);" );
}


size_t DrfLibraryList::numDetectorRows() const
{
  return m_numDetectorRows;
}


Wt::Signal<std::string> &DrfLibraryList::picked()
{
  return m_picked;
}


Wt::Signal<> &DrfLibraryList::importRequested()
{
  return m_importRequested;
}


Wt::WString DrfLibraryList::subText( const DrfLibrary::Entry &entry, const DrfLibrary::Row &row ) const
{
  using namespace DrfLibrary;

  if( entry.pending )
    return WString::tr( "dsl-reading" );
  if( !entry.summary )
    return WString::tr( "dsl-could-not-read" );  //the error itself is the row's tooltip

  // Why a suggestion is suggested, in the group that says so
  if( row.group == Group::Suggested )
  {
    if( entry.suggestReasons & SuggestForSerialNumber )
      return WString::tr( "dsl-why-serial" );
    if( entry.suggestReasons & SuggestForModel )
      return WString::tr( "dsl-why-model" );
    if( entry.suggestReasons & SuggestInSpectrumFile )
      return WString::tr( "dsl-why-file" );
    if( entry.suggestReasons & SuggestDefaultForType )
      return WString::tr( "dsl-why-default" );
  }//if( a suggestion )

  const string &desc = entry.summary->description;
  if( row.nested )
    return desc.empty() ? WString::tr( "dsl-your-copy" ) : WString::tr( "dsl-your-copy-desc" ).arg( desc );

  if( row.copyOfHidden && m_entries )
  {
    const int parent = DrfLibrary::parentIndex( *m_entries, static_cast<size_t>(&entry - &(*m_entries)[0]) );
    const Entry * const p = (parent >= 0) ? &(*m_entries)[parent] : nullptr;
    if( p && p->summary )
      return WString::tr( "dsl-copy-of" ).arg( p->summary->name );
    return WString::tr( "dsl-copy-of-missing" );
  }//if( a copy listed on its own )

  return WString::fromUTF8( desc );
}//subText(...)


void DrfLibraryList::rebuild()
{
  using namespace DrfLibrary;

  m_rows->clear();
  m_pickedRow = nullptr;
  m_rowButtons.clear();
  m_parentsWithCopies.clear();
  m_numDetectorRows = 0;

  const bool have_entries = (m_entries && !m_entries->empty());
  const vector<Row> rows = have_entries
                           ? buildRows( *m_entries, m_filter, m_expanded, m_pickedId, m_inUseId, m_extraText )
                           : vector<Row>();

  if( rows.empty() )
  {
    WContainerWidget *empty = m_rows->addNew<WContainerWidget>();
    empty->addStyleClass( "DrfLibEmpty" );

    if( m_checking && !have_entries )
    {
      WText *txt = empty->addNew<WText>( WString::tr("dsl-checking"), TextFormat::Plain );
      txt->addStyleClass( "DrfLibChecking" );
      return;
    }

    const string query = SpecUtils::trim_copy( m_filter.query );
    WText *title = empty->addNew<WText>( query.empty() ? WString::tr("dsl-empty-source")
                                                       : WString::tr("dsl-empty-query").arg(query),
                                         TextFormat::Plain );
    title->addStyleClass( "DrfLibEmptyTitle" );
    title->setInline( false );
    WText *hint = empty->addNew<WText>( WString::tr("dsl-empty-hint"), TextFormat::Plain );
    hint->setInline( false );
    WPushButton *import_btn = empty->addNew<WPushButton>( WString::tr("dsl-import-a-detector") );
    import_btn->addStyleClass( "DrfBtn" );
    import_btn->clicked().connect( [this](){ m_importRequested.emit(); } );
    return;
  }//if( rows.empty() )

  for( const Row &row : rows )
  {
    switch( row.type )
    {
      case RowType::GroupHeader:
      {
        WContainerWidget *header = m_rows->addNew<WContainerWidget>();
        header->addStyleClass( "DrfLibGroup" );
        const char *key = "dsl-grp-all";
        switch( row.group )
        {
          case Group::Suggested:  key = "dsl-grp-suggested";  break;
          case Group::Recent:     key = "dsl-grp-recent";     break;
          case Group::Everything: key = "dsl-grp-everything"; break;
          case Group::All:        key = "dsl-grp-all";        break;
        }
        WText *title = header->addNew<WText>( WString::tr(key), TextFormat::Plain );
        title->addStyleClass( "DrfLibGroupTitle" );
        if( (row.group == Group::Suggested) && !m_suggestedNote.empty() )
          header->addNew<WText>( m_suggestedNote, TextFormat::Plain )->addStyleClass( "DrfLibGroupNote" );
        break;
      }//case RowType::GroupHeader:

      case RowType::Detector:
      {
        const Entry &entry = (*m_entries)[row.entry];
        const Summary * const summary = entry.summary.get();
        const bool picked = (entry.id == m_pickedId);
        const bool in_use = (!m_inUseId.empty() && (entry.id == m_inUseId));
        const bool suggested_why = (row.group == Group::Suggested) && entry.suggestReasons;

        const string name = summary ? summary->name : SpecUtils::filename( entry.path );
        const string crystal = summary ? DrfLibraryText::crystalLabel( *summary ) : string();
        const WString range = summary ? DrfLibraryText::rangeLabel( summary->range ) : WString();
        const WString badge = (summary && summary->parentHash) ? WString::tr( "dsl-badge-your-copy" )
                              : (!entry.sourceLabel.empty() ? WString::fromUTF8( entry.sourceLabel )
                                                            : DrfLibraryText::originLabel( entry.origin ));

        // On a phone, the columns fold into the second line: crystal · source · last used
        WString sub = subText( entry, row );
        if( m_narrow && !suggested_why && !row.nested && !row.copyOfHidden && summary )
        {
          sub = WString::fromUTF8( crystal );
          const WString parts[] = { badge, (entry.lastUsedUtc > 0) ? DrfLibraryText::lastUsedLabel( entry.lastUsedUtc ) : WString() };
          for( const WString &part : parts )
          {
            if( part.empty() )
              continue;
            if( !sub.empty() )
              sub += WString::fromUTF8( " · " );
            sub += part;
          }
        }//if( m_narrow )

        string markup = "<span class=\"DrfLibCol DrfLibNameCol\"><span class=\"DrfLibNameLine\">";
        if( row.nested )
          markup += "<span class=\"DrfLibNestArrow\" aria-hidden=\"true\">↳</span>";
        markup += "<span class=\"DrfLibName\">" + html( name ) + "</span>";
        if( in_use )
          markup += "<span class=\"DrfLibInUse\">" + html( WString::tr("dsl-in-use-pill") ) + "</span>";
        markup += "</span><span class=\"DrfLibSub" + string(suggested_why ? " DrfLibWhy" : "") + "\">"
                  + html( sub ) + "</span></span>";
        if( m_narrow )
        {
          markup += "<span class=\"DrfLibChevron\" aria-hidden=\"true\">›</span>";
        }else
        {
          markup += "<span class=\"DrfLibCol\"><span class=\"DrfLibCrystal\">" + html( crystal )
                    + "</span><span class=\"DrfLibRange\">" + html( range ) + "</span></span>";
          markup += "<span class=\"DrfLibCol\"><span class=\"DrfLibBadge"
                    + string(entry.origin == Origin::ThisSpectrum ? " Spectrum" : "") + "\">"
                    + html( badge ) + "</span></span>";
          markup += "<span class=\"DrfLibUsed\">" + html( DrfLibraryText::lastUsedLabel( entry.lastUsedUtc ) ) + "</span>";
        }

        WPushButton *btn = m_rows->addNew<WPushButton>( WString::fromUTF8( markup ), TextFormat::XHTML );
        btn->addStyleClass( m_narrow ? "DrfLibRow Narrow" : "DrfLibRow DrfLibGrid" );
        if( row.nested )
          btn->addStyleClass( "Nested" );
        if( picked )
          btn->addStyleClass( "Picked" );
        if( entry.pending )
        {
          btn->addStyleClass( "Pending" );
          btn->disable();  //nothing to preview until it is read
        }
        if( !entry.summary && !entry.pending && !entry.error.empty() )
          btn->setToolTip( WString::fromUTF8( entry.error ) );
        btn->setAttributeValue( "aria-pressed", picked ? "true" : "false" );

        const string id = entry.id;
        btn->clicked().connect( [this, id](){ m_picked.emit( id ); } );

        if( picked )
          m_pickedRow = btn;
        m_rowButtons[id] = btn;
        if( row.nested )
          m_parentsWithCopies.insert( id );  //picking away from a copy can hide its family
        m_numDetectorRows += 1;
        break;
      }//case RowType::Detector:

      case RowType::CopiesToggle:
      {
        const Entry &parent = (*m_entries)[row.entry];
        const int n = static_cast<int>( row.numCopies );
        const char *key = row.expanded ? ((n == 1) ? "dsl-copies-hide-1" : "dsl-copies-hide-n")
                                       : ((n == 1) ? "dsl-copies-show-1" : "dsl-copies-show-n");
        WPushButton *btn = m_rows->addNew<WPushButton>( WString::tr(key).arg( n ) );
        btn->addStyleClass( "DrfLibCopies" );
        btn->toggleStyleClass( "Open", row.expanded );
        btn->setAttributeValue( "aria-expanded", row.expanded ? "true" : "false" );

        const string parent_id = parent.id;
        const bool expanded = row.expanded;
        m_parentsWithCopies.insert( parent_id );
        btn->clicked().connect( [this, parent_id, expanded](){
          m_expanded[parent_id] = !expanded;
          rebuild();
        } );
        break;
      }//case RowType::CopiesToggle:
    }//switch( row.type )
  }//for( const Row &row : rows )
}//void DrfLibraryList::rebuild()


// ---------------------------------------------------------------------------------------------
//  DrfPreviewPanel
// ---------------------------------------------------------------------------------------------

DrfPreviewPanel::DrfPreviewPanel( const Context context )
  : WContainerWidget(),
    m_context( context ),
    m_summary(),
    m_details(),
    m_drf(),
    m_reference(),
    m_chartShowsFwhm( false ),
    m_content( nullptr ),
    m_empty( nullptr ),
    m_emptyTitle( nullptr ),
    m_emptyText( nullptr ),
    m_tag( nullptr ),
    m_icon( nullptr ),
    m_name( nullptr ),
    m_meta( nullptr ),
    m_kidsNote( nullptr ),
    m_warning( nullptr ),
    m_lineage( nullptr ),
    m_parentLink( nullptr ),
    m_lineageEnd( nullptr ),
    m_effBtn( nullptr ),
    m_fwhmBtn( nullptr ),
    m_chartHolder( nullptr ),
    m_chart( nullptr ),
    m_chartMsg( nullptr ),
    m_legend( nullptr ),
    m_legendSolid( nullptr ),
    m_legendDashed( nullptr ),
    m_legendSolidTxt( nullptr ),
    m_legendDashedTxt( nullptr ),
    m_caps( nullptr ),
    m_where( nullptr ),
    m_actions( nullptr ),
    m_customizeBtn( nullptr ),
    m_qrBtn( nullptr ),
    m_deleteBtn( nullptr ),
    m_useNoneBtn( nullptr ),
    m_xmlBtn( nullptr )
{
  addStyleClass( "DrfPreview" );

  // -- Placeholder when nothing is picked
  m_empty = addNew<WContainerWidget>();
  m_empty->addStyleClass( "DrfPvEmpty" );
  WContainerWidget *empty_icon = m_empty->addNew<WContainerWidget>();
  empty_icon->addStyleClass( "DrfPvEmptyIcon" );
  m_emptyTitle = m_empty->addNew<WText>( "", TextFormat::Plain );
  m_emptyTitle->addStyleClass( "DrfPvEmptyTitle" );
  m_emptyTitle->setInline( false );
  m_emptyText = m_empty->addNew<WText>( "", TextFormat::Plain );
  m_emptyText->setInline( false );

  // -- The preview
  m_content = addNew<WContainerWidget>();
  m_content->addStyleClass( "DrfPvContent" );

  WContainerWidget *head = m_content->addNew<WContainerWidget>();
  head->addStyleClass( "DrfPvHead" );
  m_tag = head->addNew<WText>( "", TextFormat::Plain );
  m_tag->addStyleClass( "DrfPvTag" );

  WContainerWidget *title = head->addNew<WContainerWidget>();
  title->addStyleClass( "DrfPvTitle" );
  m_icon = title->addNew<WText>( "", TextFormat::UnsafeXHTML );  //see DrfAppearance::icon_html
  m_icon->addStyleClass( "DrfPvIcon" );
  m_name = title->addNew<WText>( "", TextFormat::Plain );
  m_name->addStyleClass( "DrfPvName" );
  m_name->addStyleClass( "SelectAllTxt" );

  m_meta = head->addNew<WText>( "", TextFormat::Plain );
  m_meta->addStyleClass( "DrfPvMeta" );
  m_meta->setInline( false );

  m_lineage = head->addNew<WContainerWidget>();
  m_lineage->addStyleClass( "DrfPvLineage" );
  m_lineage->addNew<WText>( WString::tr("dsl-pv-your-copy-of"), TextFormat::Plain );
  m_parentLink = m_lineage->addNew<WPushButton>();
  m_parentLink->addStyleClass( "DrfLinkBtn" );
  m_parentLink->clicked().connect( [this](){ m_parentRequested.emit(); } );
  m_lineageEnd = m_lineage->addNew<WText>( "", TextFormat::Plain );

  m_kidsNote = head->addNew<WText>( "", TextFormat::Plain );
  m_kidsNote->addStyleClass( "DrfPvKids" );
  m_kidsNote->setInline( false );
  m_warning = head->addNew<WText>( "", TextFormat::Plain );
  m_warning->addStyleClass( "DrfPvWarning WarnTxt" );
  m_warning->setInline( false );
  m_warning->hide();

  // Efficiency / FWHM segmented control
  WContainerWidget *seg = m_content->addNew<WContainerWidget>();
  seg->addStyleClass( "DrfSeg" );
  seg->setAttributeValue( "role", "group" );
  seg->setAttributeValue( "aria-label", WString::tr("dsl-chart-shows") );
  m_effBtn = seg->addNew<WPushButton>( WString::tr("dsl-chart-eff") );
  m_fwhmBtn = seg->addNew<WPushButton>( WString::tr("dsl-chart-fwhm") );
  m_effBtn->clicked().connect( [this](){ setChartShowsFwhm( false ); } );
  m_fwhmBtn->clicked().connect( [this](){ setChartShowsFwhm( true ); } );

  m_chartHolder = m_content->addNew<WContainerWidget>();
  m_chartHolder->addStyleClass( "DrfPvChart" );
  m_chart = m_chartHolder->addNew<DrfChart>();
  m_chart->setShowFwhm( false );
  m_chartMsg = m_chartHolder->addNew<WText>( "", TextFormat::Plain );
  m_chartMsg->addStyleClass( "DrfPvChartMsg" );
  m_chartMsg->hide();

  m_legend = m_content->addNew<WContainerWidget>();
  m_legend->addStyleClass( "DrfPvLegend" );
  WContainerWidget * const legend = m_legend;
  m_legendSolid = legend->addNew<WContainerWidget>();
  m_legendSolid->addNew<WContainerWidget>()->addStyleClass( "DrfPvSwatch" );
  m_legendSolidTxt = m_legendSolid->addNew<WText>( "", TextFormat::Plain );
  m_legendDashed = legend->addNew<WContainerWidget>();
  m_legendDashed->addNew<WContainerWidget>()->addStyleClass( "DrfPvSwatch Dashed" );
  m_legendDashedTxt = m_legendDashed->addNew<WText>( "", TextFormat::Plain );

  WContainerWidget *caps = m_content->addNew<WContainerWidget>();
  caps->addStyleClass( "DrfPvCaps" );
  WText *caps_title = caps->addNew<WText>( WString::tr("dsl-includes"), TextFormat::Plain );
  caps_title->addStyleClass( "DrfSectionTitle" );
  caps_title->setInline( false );
  m_caps = caps->addNew<WContainerWidget>();

  m_where = m_content->addNew<WText>( "", TextFormat::XHTML );
  m_where->addStyleClass( "DrfPvWhere" );
  m_where->addStyleClass( "SelectableTxt" );
  m_where->setInline( false );

  m_actions = m_content->addNew<WContainerWidget>();
  m_actions->addStyleClass( "DrfPvActions" );

  m_customizeBtn = m_actions->addNew<WPushButton>( WString::tr("dsl-customize") );
  m_customizeBtn->addStyleClass( "DrfBtn DrfBtnCustomize" );
  m_customizeBtn->clicked().connect( [this](){ m_customizeRequested.emit(); } );

  const shared_ptr<DrfXmlResource> xml_resource = make_shared<DrfXmlResource>( this );
#if( BUILD_AS_OSX_APP || IOS )
  WAnchor *xml_anchor = m_actions->addNew<WAnchor>( WLink( xml_resource ), "XML" );
  xml_anchor->setTarget( AnchorTarget::NewWindow );
  xml_anchor->addStyleClass( "DrfBtn DrfBtnXml" );
  m_xmlBtn = xml_anchor;
#else
  WPushButton *xml_btn = m_actions->addNew<WPushButton>( "XML" );
  xml_btn->setLink( WLink( xml_resource ) );
  xml_btn->addStyleClass( "DrfBtn DrfBtnXml" );
  m_xmlBtn = xml_btn;
#endif

  m_qrBtn = m_actions->addNew<WPushButton>( WString::tr("dsl-qr") );
  m_qrBtn->addStyleClass( "DrfBtn DrfBtnQr" );
  m_qrBtn->clicked().connect( [this](){ m_qrRequested.emit(); } );
#if( !USE_QR_CODES )
  m_qrBtn->hide();
#endif

  m_deleteBtn = m_actions->addNew<WPushButton>( WString::tr("dsl-delete") );
  m_deleteBtn->addStyleClass( "DrfBtn DrfBtnDelete" );
  m_deleteBtn->clicked().connect( [this](){ m_deleteRequested.emit(); } );

  m_useNoneBtn = m_actions->addNew<WPushButton>( WString::tr("dsl-use-none") );
  m_useNoneBtn->addStyleClass( "DrfBtn DrfBtnUseNone" );
  m_useNoneBtn->clicked().connect( [this](){ m_useNoneRequested.emit(); } );
  m_useNoneBtn->hide();

  setChartShowsFwhm( false );
  showNothing();
}//DrfPreviewPanel constructor


DrfPreviewPanel::~DrfPreviewPanel()
{
}


void DrfPreviewPanel::setNarrow( const bool narrow )
{
  toggleStyleClass( "Narrow", narrow );
  const WString label = narrow ? WString::tr("dsl-share") : WString::fromUTF8( "XML" );
  if( WPushButton * const btn = dynamic_cast<WPushButton *>( m_xmlBtn ) )
    btn->setText( label );
  else if( WAnchor * const anchor = dynamic_cast<WAnchor *>( m_xmlBtn ) )
    anchor->setText( label );
}//setNarrow(...)


void DrfPreviewPanel::setContext( const Context context )
{
  m_context = context;
  if( m_summary )
    showSummary( m_summary, m_details );
  else
    showNothing();
}//setContext(...)


void DrfPreviewPanel::showNothing()
{
  m_summary.reset();
  m_drf.reset();
  m_details = Details();
  m_chart->updateChart( nullptr );
  m_chart->setReferenceDetector( nullptr );

  const bool import = (m_context == Context::Import);
  m_emptyTitle->setText( WString::tr( import ? "dsl-pv-empty-import-title" : "dsl-pv-empty-title" ) );
  m_emptyText->setText( WString::tr( import ? "dsl-pv-empty-import-text" : "dsl-pv-empty-text" ) );
  m_empty->show();
  m_content->hide();
}//showNothing()


void DrfPreviewPanel::showSummary( std::shared_ptr<const DrfLibrary::Summary> summary,
                                   const Details &details )
{
  if( !summary )
  {
    showNothing();
    return;
  }

  const bool same_detector = (m_summary && (m_summary->hash == summary->hash));
  m_summary = summary;
  m_details = details;
  if( !same_detector )
    m_drf.reset();

  m_empty->hide();
  m_content->show();

  // Tag
  const bool import = (m_context == Context::Import);
  m_tag->setText( WString::tr( details.inUse ? "dsl-pv-tag-in-use"
                               : (import ? "dsl-pv-tag-import" : "dsl-pv-tag-preview") ) );
  m_tag->toggleStyleClass( "InUse", details.inUse );

  // Icon, name, meta
  m_icon->setText( WString::fromUTF8( DrfLibraryText::iconHtml( summary.get(), 16.0, 30.0 ) ) );
  m_name->setText( WString::fromUTF8( summary->name ) );
  m_meta->setText( DrfLibraryText::metaLine( *summary, details.origin ) );
  if( DrfLibraryText::relEffLabel( *summary ).empty() )
    m_meta->setToolTip( WString() );
  else
    m_meta->setToolTip( WString::tr("dsl-rel-eff-tip") );

  // Lineage
  m_lineage->setHidden( details.parentName.empty() );
  m_parentLink->setText( WString::fromUTF8( details.parentName ) );
  m_lineageEnd->setText( summary->description.empty() ? WString(".")
                                                      : WString::fromUTF8( ": " + summary->description ) );

  const int num_copies = static_cast<int>( details.numCopies );
  m_kidsNote->setHidden( import || (num_copies == 0) );
  if( num_copies )
    m_kidsNote->setText( WString::tr( (num_copies == 1) ? "dsl-pv-kids-1" : "dsl-pv-kids-n" ).arg( num_copies ) );

  m_warning->setText( details.warning );
  m_warning->setHidden( details.warning.empty() );

  // What it includes
  m_caps->clear();
  for( const Capability &cap : capabilities( *summary ) )
  {
    WContainerWidget *row = m_caps->addNew<WContainerWidget>();
    row->addStyleClass( cap.ok ? "DrfCap Yes" : "DrfCap No" );
    row->addNew<WContainerWidget>()->addStyleClass( "DrfCapMark" );
    row->addNew<WText>( cap.label, TextFormat::Plain )->addStyleClass( "DrfCapLabel" );
    row->addNew<WText>( cap.detail, TextFormat::Plain )->addStyleClass( "DrfCapDetail" );
  }

  // Where it came from
  if( details.where.empty() )
  {
    m_where->setText( WString() );
    m_where->hide();
  }else if( details.whereIsPath )
  {
    const string label = html( WString::tr( import ? "dsl-pv-from" : "dsl-pv-found-in" ) );
    m_where->setText( WString::fromUTF8( label + " <span class=\"DrfMono\">" + html( details.where ) + "</span>" ) );
    m_where->setToolTip( WString::fromUTF8( details.where ) );
    m_where->show();
  }else
  {
    m_where->setText( WString::fromUTF8( html( details.where ) ) );
    m_where->setToolTip( WString() );
    m_where->show();
  }

  // Actions: only for Library detectors
  m_actions->setHidden( import );
  m_deleteBtn->setHidden( !details.deletable );
  m_useNoneBtn->setHidden( !details.inUse || !hasStyleClass( "Narrow" ) );

  updateChart();
}//showSummary(...)


void DrfPreviewPanel::setDetector( std::shared_ptr<const DetectorPeakResponse> drf )
{
  m_drf = drf;
  updateChart();
}//setDetector(...)


void DrfPreviewPanel::setReference( std::shared_ptr<const DetectorPeakResponse> drf )
{
  m_reference = drf;
  updateChart();
}//setReference(...)


void DrfPreviewPanel::showError( const Wt::WString &message )
{
  m_drf.reset();
  m_chart->updateChart( nullptr );
  m_chartMsg->setText( message );
  m_chartMsg->show();
}//showError(...)


std::shared_ptr<const DetectorPeakResponse> DrfPreviewPanel::detector() const
{
  return m_drf;
}


Wt::Signal<> &DrfPreviewPanel::customizeRequested()
{
  return m_customizeRequested;
}


Wt::Signal<> &DrfPreviewPanel::qrRequested()
{
  return m_qrRequested;
}


Wt::Signal<> &DrfPreviewPanel::deleteRequested()
{
  return m_deleteRequested;
}


Wt::Signal<> &DrfPreviewPanel::parentRequested()
{
  return m_parentRequested;
}


Wt::Signal<> &DrfPreviewPanel::useNoneRequested()
{
  return m_useNoneRequested;
}


void DrfPreviewPanel::setChartShowsFwhm( const bool fwhm )
{
  m_chartShowsFwhm = fwhm;
  set_pressed( m_effBtn, !fwhm );
  set_pressed( m_fwhmBtn, fwhm );
  m_chart->setShowFwhm( fwhm );
  m_chart->setShowEfficiency( !fwhm );
  m_legend->toggleStyleClass( "Fwhm", fwhm );  //the swatches take the FWHM curve's colour
  updateLegendAndMessage();
}//setChartShowsFwhm(...)


void DrfPreviewPanel::updateChart()
{
  const bool have_drf = !!m_drf;
  m_chart->updateChart( m_drf );

  // The detector in use, dashed - unless it is the one previewed.
  const bool show_reference = m_reference && !m_details.inUse;
  m_chart->setReferenceDetector( show_reference ? m_reference : nullptr );

  m_xmlBtn->setDisabled( !have_drf || !m_drf->isValid() );
  m_qrBtn->setDisabled( !have_drf || !m_drf->isValid() );
  m_customizeBtn->setDisabled( !have_drf );

  updateLegendAndMessage();
}//updateChart()


void DrfPreviewPanel::updateLegendAndMessage()
{
  if( !m_summary )
    return;

  const bool show_reference = m_reference && !m_details.inUse;
  m_legendSolidTxt->setText( WString::tr("dsl-legend-solid").arg( m_summary->name ) );
  m_legendDashed->setHidden( !show_reference );
  if( show_reference )
    m_legendDashedTxt->setText( WString::tr("dsl-legend-dashed").arg( m_reference->name() ) );

  WString msg;
  if( !m_chartShowsFwhm && !m_summary->hasEfficiency )
    msg = WString::tr( "dsl-chart-no-eff" );
  else if( m_chartShowsFwhm && !m_summary->hasFwhm )
    msg = WString::tr( "dsl-chart-no-fwhm" );

  m_chartMsg->setText( msg );
  m_chartMsg->setHidden( msg.empty() );
}//updateLegendAndMessage()


// ---------------------------------------------------------------------------------------------
//  DrfLocationsPage
// ---------------------------------------------------------------------------------------------

DrfLocationsPage::DrfLocationsPage()
  : WContainerWidget(),
    m_autoRows( nullptr ),
    m_userRows( nullptr ),
    m_yours( nullptr )
{
  addStyleClass( "DrfLocations" );

  WPushButton *back = addNew<WPushButton>( WString::tr("dsl-loc-back") );
  back->addStyleClass( "DrfLinkBtn DrfLocBack" );
  back->clicked().connect( [this](){ m_backRequested.emit(); } );

  WText *title = addNew<WText>( WString::tr("dsl-loc-title"), TextFormat::Plain );
  title->addStyleClass( "DrfLocTitle" );
  title->setInline( false );
  WText *desc = addNew<WText>( WString::tr("dsl-loc-desc"), TextFormat::Plain );
  desc->addStyleClass( "DrfLocDesc" );
  desc->setInline( false );

  WText *auto_title = addNew<WText>( WString::tr("dsl-loc-auto"), TextFormat::Plain );
  auto_title->addStyleClass( "DrfSectionTitle" );
  auto_title->setInline( false );
  m_autoRows = addNew<WContainerWidget>();
  m_autoRows->addStyleClass( "DrfLocTable" );

  WText *user_title = addNew<WText>( WString::tr("dsl-loc-user"), TextFormat::Plain );
  user_title->addStyleClass( "DrfSectionTitle" );
  user_title->setInline( false );
  m_userRows = addNew<WContainerWidget>();
  m_userRows->addStyleClass( "DrfLocTable" );

  WContainerWidget *add = addNew<WContainerWidget>();
  add->addStyleClass( "DrfLocAdd" );
  WPushButton *add_folder = add->addNew<WPushButton>( WString::tr("dsl-loc-add-folder") );
  add_folder->addStyleClass( "DrfBtn DrfBtnAddFolder" );
  add_folder->clicked().connect( [this](){ m_addFolderRequested.emit(); } );
  WPushButton *add_file = add->addNew<WPushButton>( WString::tr("dsl-loc-add-file") );
  add_file->addStyleClass( "DrfBtn DrfBtnAddFile" );
  add_file->clicked().connect( [this](){ m_addFileRequested.emit(); } );

  m_yours = addNew<WText>( "", TextFormat::Plain );
  m_yours->addStyleClass( "DrfLocYours" );
  m_yours->setInline( false );
}//DrfLocationsPage constructor


void DrfLocationsPage::setLocations( const std::vector<DrfLibrary::Location> &locations,
                                     const std::vector<DrfLibrary::LocationStatus> &statuses,
                                     const size_t numMadeByYou )
{
  using DrfLibrary::LocationRole;

  m_autoRows->clear();
  m_userRows->clear();

  for( const DrfLibrary::Location &loc : locations )
  {
    const DrfLibrary::LocationStatus *status = nullptr;
    for( const DrfLibrary::LocationStatus &s : statuses )
    {
      if( s.id == loc.id )
        status = &s;
    }

    const char *name_key = nullptr, *kind_key = nullptr;
    switch( loc.role )
    {
      case LocationRole::GenericGadras: name_key = "dsl-loc-generic";  kind_key = "dsl-loc-generic-kind";  break;
      case LocationRole::ExtraGadras:   name_key = "dsl-loc-extra";    kind_key = "dsl-loc-generic-kind";  break;
      case LocationRole::BuiltInTable:  name_key = "dsl-loc-table";    kind_key = "dsl-loc-table-kind";    break;
      case LocationRole::KeptImports:   name_key = "dsl-loc-kept";     kind_key = "dsl-loc-kept-kind";     break;
      case LocationRole::UserData:      name_key = "dsl-loc-userdata"; kind_key = "dsl-loc-userdata-kind"; break;
      case LocationRole::GadrasInstall: name_key = "dsl-loc-gadras";   kind_key = "dsl-loc-install-kind";  break;
      case LocationRole::IsocsInstall:  name_key = "dsl-loc-isocs";    kind_key = "dsl-loc-install-kind";  break;
      case LocationRole::AngleInstall:  name_key = "dsl-loc-angle";    kind_key = "dsl-loc-install-kind";  break;
      case LocationRole::Package:       kind_key = "dsl-loc-package-kind"; break;
      case LocationRole::User:
        kind_key = loc.isFolder ? "dsl-loc-user-folder-kind" : "dsl-loc-user-file-kind";
        break;
    }//switch( loc.role )

    const bool is_user = (loc.role == LocationRole::User);
    const bool is_install = (loc.role == LocationRole::GadrasInstall)
                            || (loc.role == LocationRole::IsocsInstall)
                            || (loc.role == LocationRole::AngleInstall);
    const bool found = status && status->found;

    WContainerWidget *row = (is_user ? m_userRows : m_autoRows)->addNew<WContainerWidget>();
    row->addStyleClass( "DrfLocRow" );
    row->addNew<WContainerWidget>()->addStyleClass( loc.isFolder ? "DrfLocIcon Folder" : "DrfLocIcon File" );

    WContainerWidget *what = row->addNew<WContainerWidget>();
    what->addStyleClass( "DrfLocWhat" );
    const WString name = name_key ? WString::tr( name_key )
                                  : WString::fromUTF8( loc.label.empty() ? SpecUtils::filename( loc.path ) : loc.label );
    WString name_line = name;
    name_line += WString::fromUTF8( " · " );
    name_line += WString::tr( kind_key );
    what->addNew<WText>( name_line, TextFormat::Plain )->addStyleClass( "DrfLocName" );
    WText *path = what->addNew<WText>( WString::fromUTF8( loc.path ), TextFormat::Plain );
    path->addStyleClass( "DrfLocPath DrfMono" );
    path->addStyleClass( "SelectAllTxt" );
    path->setToolTip( WString::fromUTF8( loc.path ) );

    // InterSpec's own folder for kept imports is only made when something is kept.
    const bool not_yet = (loc.role == LocationRole::KeptImports) && status && !found;

    WString status_txt;
    if( !status )
      status_txt = WString::tr( "dsl-loc-checking" );
    else if( not_yet )
      status_txt = WString::tr( "dsl-loc-kept-none" );
    else if( !found )
      status_txt = WString::tr( "dsl-loc-missing" );
    else if( status->numDetectors == 0 )
      status_txt = WString::tr( "dsl-loc-empty" );
    else
      status_txt = WString::tr( (status->numDetectors == 1) ? "dsl-loc-num-1" : "dsl-loc-num-n" )
                     .arg( static_cast<int>(status->numDetectors) );
    if( found && status->numUnreadable )
    {
      status_txt += WString::fromUTF8( " · " );
      status_txt += WString::tr( "dsl-loc-unreadable" ).arg( static_cast<int>(status->numUnreadable) );
    }
    WText *status_w = row->addNew<WText>( status_txt, TextFormat::Plain );
    status_w->addStyleClass( "DrfLocStatus" );
    status_w->toggleStyleClass( "Found", found && status->numDetectors );
    status_w->toggleStyleClass( "Missing", status && !found && !not_yet );

    WCheckBox *shown = row->addNew<WCheckBox>( WString::tr("dsl-loc-show") );
    shown->addStyleClass( "DrfLocShown" );
    shown->setChecked( !loc.hidden );
    shown->setDisabled( status && !found );  //nothing to show or hide
    const string id = loc.id;
    shown->changed().connect( [this, id, shown](){ m_shownChanged.emit( id, shown->isChecked() ); } );

    WContainerWidget *actions = row->addNew<WContainerWidget>();
    actions->addStyleClass( "DrfLocActions" );
    if( found && !loc.builtIn )
    {
      WPushButton *btn = actions->addNew<WPushButton>( WString::tr("dsl-loc-rescan") );
      btn->addStyleClass( "DrfBtn" );
      btn->clicked().connect( [this, id](){ m_rescanRequested.emit( id ); } );
    }
    if( is_install )
    {
      WPushButton *btn = actions->addNew<WPushButton>( WString::tr("dsl-loc-locate") );
      btn->addStyleClass( "DrfBtn" );
      btn->clicked().connect( [this, id](){ m_locateRequested.emit( id ); } );
      btn->setHidden( found );
    }
    if( is_user )
    {
      WPushButton *btn = actions->addNew<WPushButton>( WString::tr("dsl-loc-remove") );
      btn->addStyleClass( "DrfBtn" );
      btn->clicked().connect( [this, id](){ m_removeRequested.emit( id ); } );
    }
  }//for( const DrfLibrary::Location &loc : locations )

  if( m_userRows->count() == 0 )
    m_userRows->addNew<WText>( WString::tr("dsl-loc-none-added"), TextFormat::Plain )->addStyleClass( "DrfLocNone" );

  m_yours->setText( WString::tr("dsl-loc-yours").arg( static_cast<int>(numMadeByYou) ) );
}//setLocations(...)


Wt::Signal<> &DrfLocationsPage::backRequested()
{
  return m_backRequested;
}

Wt::Signal<std::string,bool> &DrfLocationsPage::shownChanged()
{
  return m_shownChanged;
}

Wt::Signal<std::string> &DrfLocationsPage::rescanRequested()
{
  return m_rescanRequested;
}

Wt::Signal<std::string> &DrfLocationsPage::locateRequested()
{
  return m_locateRequested;
}

Wt::Signal<std::string> &DrfLocationsPage::removeRequested()
{
  return m_removeRequested;
}

Wt::Signal<> &DrfLocationsPage::addFolderRequested()
{
  return m_addFolderRequested;
}

Wt::Signal<> &DrfLocationsPage::addFileRequested()
{
  return m_addFileRequested;
}


// ---------------------------------------------------------------------------------------------
//  DrfFoundList
// ---------------------------------------------------------------------------------------------

DrfFoundList::DrfFoundList()
  : WContainerWidget(),
    m_title( nullptr ),
    m_subtitle( nullptr ),
    m_summary( nullptr ),
    m_hint( nullptr ),
    m_rows( nullptr ),
    m_keepRow( nullptr ),
    m_keep( nullptr ),
    m_keepDesc( nullptr )
{
  addStyleClass( "DrfFound" );

  WContainerWidget *head = addNew<WContainerWidget>();
  head->addStyleClass( "DrfFoundHead" );
  head->addNew<WContainerWidget>()->addStyleClass( "DrfLocIcon Folder" );
  WContainerWidget *what = head->addNew<WContainerWidget>();
  what->addStyleClass( "DrfFoundWhat" );
  m_title = what->addNew<WText>( "", TextFormat::Plain );
  m_title->addStyleClass( "DrfFoundTitle" );
  m_title->setInline( false );
  m_subtitle = what->addNew<WText>( "", TextFormat::Plain );
  m_subtitle->addStyleClass( "DrfFoundSub DrfMono" );
  m_subtitle->setInline( false );
  m_summary = head->addNew<WText>( "", TextFormat::Plain );
  m_summary->addStyleClass( "DrfFoundSummary" );
  WPushButton *clear = head->addNew<WPushButton>();
  clear->addStyleClass( "DrfNoticeClose DrfFoundClear" );
  clear->setAttributeValue( "aria-label", WString::tr("dsl-found-clear") );
  clear->clicked().connect( [this](){ m_cleared.emit(); } );

  m_hint = addNew<WText>( WString::tr("dsl-found-hint"), TextFormat::Plain );
  m_hint->addStyleClass( "DrfFormHint" );
  m_hint->setInline( false );

  m_rows = addNew<WContainerWidget>();
  m_rows->addStyleClass( "DrfFoundRows" );

  m_keepRow = addNew<WContainerWidget>();
  m_keepRow->addStyleClass( "DrfFoundKeep" );
  m_keep = m_keepRow->addNew<WCheckBox>( WString::tr("dsl-found-keep") );
  m_keep->setChecked( true );
  m_keepDesc = m_keepRow->addNew<WText>( WString::tr("dsl-found-keep-desc"), TextFormat::Plain );
  m_keepDesc->addStyleClass( "DrfFormHint" );
  m_keepDesc->setInline( false );
}//DrfFoundList constructor


void DrfFoundList::setReading( const Wt::WString &message )
{
  m_title->setText( message );
  m_subtitle->setText( WString() );
  m_summary->setText( WString() );
  m_hint->hide();
  m_rows->clear();
  m_items.clear();
  m_keepRow->hide();
}//setReading(...)


void DrfFoundList::setFound( const Wt::WString &title, const Wt::WString &subtitle,
                             const Wt::WString &summary, const std::vector<Item> &items,
                             const bool canKeep )
{
  m_title->setText( title );
  m_subtitle->setText( subtitle );
  m_summary->setText( summary );
  m_hint->setHidden( items.empty() );
  m_keepRow->setHidden( !canKeep || items.empty() );

  m_rows->clear();
  m_items.clear();
  for( size_t i = 0; i < items.size(); ++i )
  {
    const Item &item = items[i];
    string markup = "<span class=\"DrfFoundRadio\" aria-hidden=\"true\"></span>"
                    "<span class=\"DrfFoundText\"><span class=\"DrfFoundName\">" + html( item.name ) + "</span>"
                    "<span class=\"DrfFoundFmt\">" + html( item.format );
    if( !item.crystal.empty() )
      markup += " · " + html( item.crystal );
    markup += "</span></span><span class=\"DrfFoundTag" + string(item.usable ? "" : " Warn") + "\">"
              + html( item.tag ) + "</span>";

    WPushButton *btn = m_rows->addNew<WPushButton>( WString::fromUTF8( markup ), TextFormat::XHTML );
    btn->addStyleClass( "DrfFoundRow" );
    btn->setAttributeValue( "aria-pressed", "false" );
    btn->clicked().connect( [this, i](){ m_picked.emit( i ); } );
    m_items.push_back( btn );
  }//for( loop over items )
}//setFound(...)


void DrfFoundList::setPicked( const size_t index )
{
  for( size_t i = 0; i < m_items.size(); ++i )
  {
    m_items[i]->toggleStyleClass( "Picked", i == index );
    m_items[i]->setAttributeValue( "aria-pressed", (i == index) ? "true" : "false" );
  }
}//setPicked(...)


void DrfFoundList::setKeepDescription( const Wt::WString &description )
{
  m_keepDesc->setText( description );
}


bool DrfFoundList::keep() const
{
  return !m_keepRow->isHidden() && m_keep->isChecked();
}


Wt::Signal<size_t> &DrfFoundList::picked()
{
  return m_picked;
}


Wt::Signal<> &DrfFoundList::cleared()
{
  return m_cleared;
}
