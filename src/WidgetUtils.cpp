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

#include <regex>
#include <cctype>
#include <memory>
#include <string>
#include <vector>
#include <cassert>
#include <algorithm>

#include <Wt/WDate.h>
#include <Wt/WPanel.h>
#include <Wt/WServer.h>
#include <Wt/WWidget.h>
#include <Wt/WDateEdit.h>
#include <Wt/WFormWidget.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>
#include <Wt/WDateValidator.h>
#include <Wt/WContainerWidget.h>

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/WidgetUtils.h"


namespace WidgetUtils
{

namespace
{
/** Destroy `doomed` on the session thread, on the next event-loop iteration.

 The widget must be destroyed *inside* the posted lambda's body, not by letting the lambda's own
 destructor drop the last reference: `WServer::post()` -> `WServer::schedule(0ms)` ->
 `asio::post(strand_, fn)`, and asio destroys the handler on the **io-service thread** after
 invoking it.  Letting `~WWidget` run there would be a disaster - it happens with the session lock
 released and `WApplication::instance()` null, so `WWebWidget::repaint()` (reached from
 `~WContainerWidget::clear()`) dereferences a null application, and `renderOk()` /
 `EventSignalBase::prepareDestruct()` silently skip un-registering the widget from the renderers
 update map and the applications exposed-signal map, leaving dangling entries.

 Running `reset()` in the body destroys the widget on the session thread under the update lock;
 whatever copies of the (now empty) holder asio later destroys are inert.
 */
void destroy_on_session_thread( std::unique_ptr<Wt::WWidget> doomed )
{
  if( !doomed )
    return;

  Wt::WApplication * const app = Wt::WApplication::instance();
  Wt::WServer * const server = Wt::WServer::instance();

  // No session/server to post to - destroy now, on this (session) thread, rather than leak.
  if( !app || !server )
    return;

  // Park ownership on wApp rather than inside the posted task.  If the task is never delivered -
  //  the session dies between the post and its pickup, so WebController::handleApplicationEvent
  //  returns without running it - asio destroys the undelivered task on an io-service thread, and
  //  a widget owned by that task would be destroyed there.  Parked on wApp instead, the worst case
  //  is that it lives until ~WApplication destroys it, on the session thread.
  Wt::WWidget * const raw = doomed.get();
  app->addChild( std::move(doomed) );

  server->post( app->sessionId(), [raw](){
    Wt::WApplication * const a = Wt::WApplication::instance();
    if( a )
    {
      // Dropping the returned unique_ptr destroys `raw` here, on the session thread, under the
      //  update lock.
      const std::unique_ptr<Wt::WObject> owned = a->removeChild( raw );
    }
  } );
}//void destroy_on_session_thread( std::unique_ptr<Wt::WWidget> )
}//namespace


void destroyLater( std::unique_ptr<Wt::WWidget> widget )
{
  destroy_on_session_thread( std::move(widget) );
}//void destroyLater( std::unique_ptr<Wt::WWidget> )


void removeWidgetNow( Wt::WWidget *child )
{
  if( !child )
    return;

  // Dropping the returned unique_ptr destroys `child` exactly once.
  const std::unique_ptr<Wt::WWidget> doomed = child->removeFromParent();

  // A null return means nobody owned `child` through its widget parent, so nothing was destroyed
  //  (see the header for which parent shapes do this).  Callers of this function are meant to be
  //  replacing a `delete child;`, so that is virtually always a bug at the call site.
  assert( doomed );
}//void removeWidgetNow( Wt::WWidget *child )


void removeWidgetLater( Wt::WWidget *child )
{
  if( !child )
    return;

  std::unique_ptr<Wt::WWidget> doomed = child->removeFromParent();
  assert( doomed );   //see #removeWidgetNow
  destroy_on_session_thread( std::move(doomed) );
}//void removeWidgetLater( Wt::WWidget *child )


void removeWidgetLater( Wt::WContainerWidget *parent, Wt::WWidget *child )
{
  if( !parent || !child )
    return;

  std::unique_ptr<Wt::WWidget> doomed = parent->removeWidget( child );
  assert( doomed );   //see #removeWidgetNow
  destroy_on_session_thread( std::move(doomed) );
}//void removeWidgetLater( Wt::WContainerWidget *parent, Wt::WWidget *child )


WidgetHandle::WidgetHandle( Wt::WWidget *widget )
 : m_id( widget ? widget->id() : std::string() )
{
}


Wt::WWidget *WidgetHandle::resolve() const
{
  if( m_id.empty() )
    return nullptr;

  Wt::WApplication * const app = Wt::WApplication::instance();
  if( !app )
    return nullptr;

  // WApplication::findById searches domRoot_ and domRoot2_.  It reaches widgets inside dialogs
  //  (WCompositeWidget::findById forwards to its implementation, and WTemplate::iterateChildren
  //  walks the bound widgets) and inside layouts/WStackedWidgets (WContainerWidget::iterateChildren
  //  walks children_ *and* layout_->iterateWidgets).
  return app->findById( m_id );
}


void trackSessionDialog( Wt::WDialog *dialog )
{
  InterSpec * const viewer = InterSpec::instance();
  if( viewer && dialog )
    viewer->trackToolDialog( dialog );
}//void trackSessionDialog( Wt::WDialog *dialog )


void applyButtonRole( Wt::WPushButton *button, ButtonRole role )
{
  if( !button )
    return;

  switch( role )
  {
    case ButtonRole::Affirm:      button->addStyleClass( "DialogBtnAffirm" );      break;
    case ButtonRole::Dismiss:     button->addStyleClass( "DialogBtnDismiss" );     break;
    case ButtonRole::Neutral:     button->addStyleClass( "DialogBtnNeutral" );     break;
    case ButtonRole::Destructive: button->addStyleClass( "DialogBtnDestructive" ); break;
  }//switch( role )
}//void applyButtonRole( Wt::WPushButton *button, ButtonRole role )


void setInputValueError( Wt::WFormWidget *w, const bool error )
{
  if( !w || (w->hasStyleClass("InputValueError") == error) )
    return;

  w->toggleStyleClass( "InputValueError", error );
  w->setAttributeValue( "aria-invalid", error ? "true" : "false" );
}//void setInputValueError( Wt::WFormWidget *w, const bool error )


void revealAndFocus( Wt::WFormWidget *w )
{
  if( !w )
    return;

  for( Wt::WWidget *p = w->parent(); p; p = p->parent() )
  {
    Wt::WPanel * const panel = dynamic_cast<Wt::WPanel *>( p );
    if( panel && panel->isCollapsible() && panel->isCollapsed() )
      panel->setCollapsed( false );
  }//for( loop over ancestors )

  w->setFocus( true );
  w->doJavaScript( "setTimeout(function(){"
                     "var e=" + w->jsRef() + ";"
                     "if(e&&e.scrollIntoView)e.scrollIntoView({block:'center',inline:'nearest'});"
                   "},0);" );
}//void revealAndFocus( Wt::WFormWidget *w )


namespace
{
  /** The (1-based) month a word names, or 0.  English names and abbreviations of at least three
   letters ("Jan", "janu", "September", "Sept"), plus the session locale's names when there is one.
   */
  int month_from_name( std::string word )
  {
    SpecUtils::to_lower_ascii( word );
    if( word.size() < 3 )
      return 0;

    static const char * const english[12] = { "january", "february", "march", "april", "may",
      "june", "july", "august", "september", "october", "november", "december" };
    for( int i = 0; i < 12; ++i )
    {
      const std::string name = english[i];
      if( (word.size() <= name.size()) && (name.compare( 0, word.size(), word ) == 0) )
        return i + 1;
    }

    // The localized names come from the message bundles, so only when there is a session.
    if( Wt::WApplication::instance() )
    {
      for( int i = 1; i <= 12; ++i )
      {
        for( std::string name : { Wt::WDate::longMonthName(i).toUTF8(), Wt::WDate::shortMonthName(i).toUTF8() } )
        {
          SpecUtils::to_lower_ascii( name );
          if( name == word )
            return i;
        }
      }
    }//if( Wt::WApplication::instance() )

    return 0;
  }//int month_from_name( std::string word )


  /** The order the year ('y'), month ('M') and day ('d') appear in a WDate format, e.g. "yMd" for
   "yyyy-MM-dd", or an empty string if the format does not have each of them once. */
  std::string format_field_order( const std::string &format )
  {
    std::string order;
    for( size_t i = 0; i < format.size(); ++i )
    {
      const char c = format[i];
      if( (c != 'y') && (c != 'M') && (c != 'd') )
        continue;
      if( (i > 0) && (format[i-1] == c) )
        continue;  //rest of the same field
      if( order.find(c) != std::string::npos )
        return "";
      order += c;
    }//for( loop over format characters )

    return (order.size() == 3) ? order : std::string();
  }//std::string format_field_order( const std::string &format )
}//namespace


DateParse parseLenientDate( const std::string &input, const std::string &format, Wt::WDate &result )
{
  std::string text = input;
  SpecUtils::trim( text );
  if( text.empty() )
    return DateParse::Empty;

  // Ignore a trailing time - a pasted timestamp ("2019-01-05T13:24:00", "2019-01-05 13:24")
  {
    const std::regex time_re( R"(^(.+?)(?:[T ]\d{1,2}:\d{2}(?::\d{2}(?:\.\d*)?)?\s*(?:[AaPp][Mm])?\s*(?:Z|[+-]\d{2}:?\d{2})?)$)" );
    std::smatch mtch;
    if( std::regex_match( text, mtch, time_re ) )
    {
      text = mtch[1].str();
      SpecUtils::trim( text );
    }
  }

  auto make_date = [&result]( const int year, const int month, const int day ) -> DateParse {
    // Checked before constructing the WDate, which logs a warning for an out-of-range field
    if( (year < 1000) || (year > 9999) || (month < 1) || (month > 12) || (day < 1) || (day > 31) )
      return DateParse::Invalid;
    const Wt::WDate date( year, month, day );
    if( !date.isValid() )
      return DateParse::Invalid;
    result = date;
    return DateParse::Valid;
  };//make_date

  // 1) Exactly the expected format
  if( !format.empty() )
  {
    const Wt::WDate exact = Wt::WDate::fromString( Wt::WString::fromUTF8(text), Wt::WString::fromUTF8(format) );
    if( exact.isValid() )
    {
      result = exact;
      return DateParse::Valid;
    }
  }//if( !format.empty() )

  std::smatch mtch;

  // 2) "yyyyMMdd"
  if( std::regex_match( text, mtch, std::regex( R"(^(\d{4})(\d{2})(\d{2})$)" ) ) )
    return make_date( std::stoi(mtch[1].str()), std::stoi(mtch[2].str()), std::stoi(mtch[3].str()) );

  // 3) Three groups of digits
  if( std::regex_match( text, mtch, std::regex( R"(^(\d{1,4})\s*[-/. ]\s*(\d{1,4})\s*[-/. ]\s*(\d{1,4})$)" ) ) )
  {
    const std::string a = mtch[1].str(), b = mtch[2].str(), c = mtch[3].str();
    const int va = std::stoi(a), vb = std::stoi(b), vc = std::stoi(c);

    // Year first is unambiguous - nobody writes year-day-month.
    if( a.size() == 4 )
      return ((b.size() <= 2) && (c.size() <= 2)) ? make_date( va, vb, vc ) : DateParse::Invalid;

    if( c.size() != 4 )
      return ((b.size() <= 2) && (a.size() <= 2)) ? DateParse::Ambiguous : DateParse::Invalid;

    // Year last: only the field's own order is assumed; "1/5/2019" in a "yyyy-MM-dd" field could be
    //  either January 5th or May 1st.
    const std::string order = format_field_order( format );
    if( order == "Mdy" )
      return make_date( vc, va, vb );
    if( order == "dMy" )
      return make_date( vc, vb, va );
    return DateParse::Ambiguous;
  }//if( three groups of digits )

  // 4) A month name, a 4-digit year, and a day - in any order
  {
    std::vector<std::string> tokens;
    SpecUtils::split( tokens, text, "-/., \t" );
    tokens.erase( std::remove_if( begin(tokens), end(tokens),
                                  []( const std::string &t ){ return t.empty(); } ), end(tokens) );

    if( tokens.size() == 3 )
    {
      int year = -1, month = -1, day = -1;
      bool ok = true;
      for( std::string t : tokens )
      {
        // "5th", "1st" - drop an English ordinal suffix from a day
        if( (t.size() > 2) && std::isdigit(static_cast<unsigned char>(t[0]))
           && std::regex_match( t, std::regex( R"(^\d{1,2}(st|nd|rd|th)$)", std::regex::icase ) ) )
          t = t.substr( 0, t.size() - 2 );

        const bool all_digits = (t.find_first_not_of("0123456789") == std::string::npos);
        if( all_digits && (t.size() == 4) && (year < 0) )
          year = std::stoi( t );
        else if( all_digits && (t.size() <= 2) && (day < 0) )
          day = std::stoi( t );
        else if( !all_digits && (month < 0) && (month_from_name(t) > 0) )
          month = month_from_name( t );
        else
          ok = false;
      }//for( std::string t : tokens )

      if( ok && (year > 0) && (month > 0) && (day > 0) )
        return make_date( year, month, day );

      // e.g. "5 Jan 19" - a month name, but a 2-digit year
      if( (month > 0) && (year < 0) )
        return DateParse::Ambiguous;
    }//if( tokens.size() == 3 )
  }

  return DateParse::Invalid;
}//DateParse parseLenientDate(...)


void setupDateEdit( Wt::WDateEdit *edit, const Wt::WString &format )
{
  if( !edit )
    return;

  edit->setFormat( format );

  std::string placeholder = format.toUTF8();
  SpecUtils::to_lower_ascii( placeholder );
  edit->setPlaceholderText( Wt::WString::fromUTF8( placeholder ) );

  edit->setBottom( Wt::WDate( 1900, 1, 1 ) );
  edit->setTop( Wt::WDate( 2200, 12, 31 ) );
}//void setupDateEdit( Wt::WDateEdit *edit, const Wt::WString &format )


DateParse normalizeDateEdit( Wt::WDateEdit *edit )
{
  if( !edit )
    return DateParse::Invalid;

  Wt::WDate date;
  const DateParse status = parseLenientDate( edit->text().toUTF8(), edit->format().toUTF8(), date );
  if( (status == DateParse::Valid) && (date.toString( edit->format() ) != edit->text()) )
    edit->setDate( date );  //re-writes the text in the edits format, which re-validates

  return status;
}//DateParse normalizeDateEdit( Wt::WDateEdit *edit )

}//namespace WidgetUtils
