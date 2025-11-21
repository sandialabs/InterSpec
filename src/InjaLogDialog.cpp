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

#include "InterSpec/InjaLogDialog.h"

#include <string>
#include <memory>
#include <sstream>

#include <Wt/WText>
#include <Wt/Utils>
#include <Wt/WComboBox>
#include <Wt/WResource>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WContainerWidget>
#include <Wt/WMemoryResource>

#include <nlohmann/json.hpp>

#include "external_libs/SpecUtils/3rdparty/inja/inja.hpp"

#include "InterSpec/InterSpec.h"
#include "InterSpec/BatchInfoLog.h"


namespace
{
  /** Helper resource class to serve rendered HTML content via a URL. */
  class InjaReportResource : public Wt::WResource
  {
  private:
    const std::string m_html_content;

  public:
    InjaReportResource( std::string &&html_content, Wt::WObject *parent = nullptr )
    : Wt::WResource( parent ), m_html_content( std::move( html_content ) )
    {
      setDispositionType( DispositionType::Inline );
    }

    virtual ~InjaReportResource()
    {
      beingDeleted();
    }

    virtual void handleRequest( const Wt::Http::Request &request, Wt::Http::Response &response ) override
    {
      response.setMimeType( "text/html; charset=utf-8" );
      response.addHeader( "Cache-Control", "no-cache, no-store, must-revalidate" );
      response.addHeader( "Pragma", "no-cache" );
      response.addHeader( "Expires", "0" );

      response.out() << m_html_content;
    }
  };//class InjaReportResource

}//namespace


const char * const InjaLogDialog::sm_txt_log_pre_wrapper = "<html><body><pre>";
const char * const InjaLogDialog::sm_txt_log_post_wrapper = "</pre></body></html>";


class InjaLogResource : public Wt::WResource
{
  // Resource class for safely downloading rendered HTML content
public:
  InjaLogResource( InjaLogDialog *parent )
  : Wt::WResource( parent ),
  m_dialog( parent ),
  m_app( Wt::WApplication::instance() )
  {
    assert( m_app );
  }
  
  virtual ~InjaLogResource()
  {
    beingDeleted();
  }

private:
  InjaLogDialog *m_dialog;
  Wt::WApplication *m_app;

  virtual void handleRequest( const Wt::Http::Request &request,
                              Wt::Http::Response &response )
  {
    Wt::WApplication::UpdateLock lock( m_app );
    
    if( !lock )
    {
      Wt::log("error") << "Failed to WApplication::UpdateLock in InjaLogResource.";
      response.out() << "Error grabbing application lock to form InjaLogResource resource; please report to InterSpec@sandia.gov.";
      response.setStatus( 500 );
      assert( 0 );
      return;
    }//if( !lock )
    
    if( !m_dialog )
    {
      response.out() << "Error: InjaLogDialog is null";
      response.setStatus( 500 );
      return;
    }
    
    // Get the filename and content for download
    const std::string filename = m_dialog->current_suggested_name();
    const std::string &content = m_dialog->current_content();

    suggestFileName( filename, Wt::WResource::Attachment );

    // Determine MIME type from file extension
    std::string mime_type = "text/html; charset=utf-8";
    if( filename.size() >= 4 && filename.compare( filename.size() - 4, 4, ".txt" ) == 0 )
      mime_type = "text/plain; charset=utf-8";

    response.setMimeType( mime_type );
    response.out() << content;
  }//InjaLogResource::handleRequest(...)
};//class InjaLogResource



InjaLogDialog::InjaLogDialog( const Wt::WString &title,
                               const nlohmann::json &data,
                               std::vector<std::tuple<Wt::WString, std::string, LogType, std::function<std::string(inja::Environment&, const nlohmann::json&)>>> template_options )
: SimpleDialog( title ),
  m_download_resource( nullptr ),
  m_data( data ),
  m_template_options( std::move( template_options ) ),
  m_toolbar( nullptr ),
  m_template_selector( nullptr ),
  m_download_button( nullptr ),
  m_iframe_holder( nullptr ),
  m_current_resource( nullptr ),
  m_current_content()
{
  addStyleClass( "InjaLogDialog" );

  Wt::WApplication *app = Wt::WApplication::instance();
  assert( app );
  app->useStyleSheet( "InterSpec_resources/InjaLogDialog.css" );
  
  // We will load the dialog size rule directly here because the InjaLogDialog will load before InjaLogDialog.css
  //  which will lead to the dialog being at most 50% the screen width - at least for the first InjaLogDialog created.
  //  But loading the rule directly here will avoid this issue.
  static const std::string chart_css_rule_size_name = "InjaLogDialog-size";
  Wt::WCssStyleSheet &style = app->styleSheet();
  if( !style.isDefined(chart_css_rule_size_name) )
    style.addRule( ".simple-dialog.InjaLogDialog", "max-width: 95vw; max-height: 95vh; width: 95vw; height: 95vw;", chart_css_rule_size_name );
  
  InterSpec *interspec = InterSpec::instance();
  if( interspec )
    interspec->useMessageResourceBundle( "InjaLogDialog" );

  // Create iframe holder
  m_iframe_holder = new Wt::WText( contents() );
  m_iframe_holder->addStyleClass( "InjaLogDialogIFrameHolder" );
  m_iframe_holder->setWidth( Wt::WLength( 100.0, Wt::WLength::Percentage ) );
  m_iframe_holder->setHeight( Wt::WLength( 100.0, Wt::WLength::Percentage ) );

  
  // Create toolbar for template selector and download button
  m_toolbar = new Wt::WContainerWidget( contents() );
  m_toolbar->addStyleClass( "InjaLogToolbar" );

  // Only show template selector if we have multiple templates
  if( m_template_options.size() > 1 )
  {
    m_template_selector = new Wt::WComboBox( m_toolbar );
    m_template_selector->addStyleClass( "InjaLogTemplateSelector" );
    for( size_t i = 0; i < m_template_options.size(); ++i )
    {
      m_template_selector->addItem( std::get<0>(m_template_options[i]) );
    }

    m_template_selector->setCurrentIndex( 0 );
    m_template_selector->changed().connect( this, &InjaLogDialog::handleTemplateChange );
  }

  // Create download resource
  m_download_resource = new InjaLogResource( this );

  // Add download button
  m_download_button = new Wt::WPushButton( m_toolbar );
  m_download_button->setIcon( "InterSpec_resources/images/download_small.svg" );
  m_download_button->setLink( Wt::WLink( static_cast<Wt::WResource*>(m_download_resource) ) );
  m_download_button->setLinkTarget( Wt::TargetNewWindow );
  m_download_button->setStyleClass( "LinkBtn DownloadBtn InjaLogDownloadBtn" );

  // Set dialog size
  resize( Wt::WLength( 95.0, Wt::WLength::Percentage ), Wt::WLength( 95.0, Wt::WLength::Percentage ) );

  // Add standard close button
  addButton( Wt::WString::tr( "Close" ) );

  // Render initial content
  updateDisplay();
}


InjaLogDialog::~InjaLogDialog()
{
  // m_current_resource will be deleted automatically as it's parented to this dialog
}


const std::string &InjaLogDialog::current_content() const
{
  return m_current_content;
}

void InjaLogDialog::setToolbarVisible( const bool show )
{
  if( m_toolbar )
  {
    if( show )
      m_toolbar->show();
    else
      m_toolbar->hide();
  }
}


const std::string InjaLogDialog::current_suggested_name() const
{
  // Get base filename from JSON data (the foreground spectrum filename)
  std::string base_filename = "report";

  try
  {
    if( m_data.contains("Filename") && m_data["Filename"].is_string() )
    {
      std::string spec_filename = m_data["Filename"].get<std::string>();

      // Remove file extension from spectrum filename
      const size_t dot_pos = spec_filename.find_last_of( '.' );
      if( dot_pos != std::string::npos )
        spec_filename = spec_filename.substr( 0, dot_pos );

      if( !spec_filename.empty() )
        base_filename = spec_filename;
    }
  } catch( ... )
  {
    // If we can't get filename from JSON, use default
  }

  // Get filename suffix from template options
  std::string filename_suffix = "_report.html";  // Default fallback

  if( m_template_selector )
  {
    const int current_index = m_template_selector->currentIndex();
    if( current_index >= 0 && static_cast<size_t>(current_index) < m_template_options.size() )
    {
      filename_suffix = std::get<1>( m_template_options[current_index] );
    }
  }
  else if( m_template_options.size() == 1 )
  {
    // If there's only one template and no selector, get its suffix
    filename_suffix = std::get<1>( m_template_options[0] );
  }

  // Append suffix to base filename
  std::string filename = base_filename + filename_suffix;

  //Remove bad filename characters
  const std::string notallowed = "\\/:?\"<>|*";
  for( auto it = begin(filename) ; it < end(filename) ; ++it )
  {
    if( notallowed.find(*it) != std::string::npos )
      *it = '_';
  }
  
  return filename;
}


void InjaLogDialog::updateDisplay()
{
  if( m_template_options.empty() )
  {
    const std::string error_msg = "<p style='color: red;'>" + Wt::WString::tr( "ild-no-templates" ).toUTF8() + "</p>";
    m_iframe_holder->setText( error_msg );
    m_iframe_holder->setTextFormat( Wt::TextFormat::XHTMLUnsafeText );
    return;
  }

  // Get the current template index
  int template_index = 0;
  if( m_template_selector )
    template_index = m_template_selector->currentIndex();

  if( template_index < 0 || static_cast<size_t>(template_index) >= m_template_options.size() )
    template_index = 0;

  try
  {
    // Create inja environment with default setup
    const std::string tmplt_dir = BatchInfoLog::default_template_dir();
    inja::Environment env = BatchInfoLog::get_default_inja_env( tmplt_dir );

    // Call the template rendering function
    const auto &render_func = std::get<3>( m_template_options[template_index] );
    m_current_content = render_func( env, m_data );

    // Get the log type
    const LogType log_type = std::get<2>( m_template_options[template_index] );

    // Prepare content for display
    std::string display_content;
    if( log_type == LogType::Text )
    {
      // For text content, escape HTML and wrap in <pre> tags
      const std::string escaped_content = Wt::Utils::htmlEncode( m_current_content, Wt::Utils::EncodeNewLines );
      display_content = std::string( sm_txt_log_pre_wrapper ) + escaped_content + sm_txt_log_post_wrapper;
    }else
    {
      // For HTML content, use as-is
      display_content = m_current_content;
    }

    // Clean up old resource if it exists
    if( m_current_resource )
    {
      delete m_current_resource;
      m_current_resource = nullptr;
    }

    // Create new resource and get its URL
    m_current_resource = new InjaReportResource( std::move( display_content ), this );
    const std::string resource_url = m_current_resource->url();

    // Create iframe HTML
    const std::string iframe_html = "<iframe "
                                    "style=\"width: 100%; height: 100%; border: none;\" "
                                    "sandbox=\"allow-scripts\" "
                                    "src=\"" +
                                    resource_url +
                                    "\" "
                                    "onerror=\"console.error('Iframe failed to load')\"> "
                                    "</iframe>";

    m_iframe_holder->setText( iframe_html );
    m_iframe_holder->setTextFormat( Wt::TextFormat::XHTMLUnsafeText );

  } catch( nlohmann::json::parse_error &e )
  {
    const std::string error_msg = "<p>" + Wt::WString::tr( "ild-json-parse-error" ).arg( std::string( e.what() ) ).toUTF8() + "</p>";
    m_iframe_holder->setText( error_msg );
    m_iframe_holder->setTextFormat( Wt::TextFormat::XHTMLUnsafeText );
    m_iframe_holder->removeStyleClass( "InjaLogDialogIFrameHolder" );
    m_iframe_holder->addStyleClass( "InjaLogDialogJsonError" );

  } catch( inja::InjaError &e )
  {
    const std::string error_msg = "<p>" +
                                  Wt::WString::tr( "ild-inja-error" )
                                    .arg( e.message )
                                    .arg( std::to_string( e.location.line ) )
                                    .arg( std::to_string( e.location.column ) )
                                    .toUTF8() +
                                  "</p>";

    m_iframe_holder->setText( error_msg );
    m_iframe_holder->setTextFormat( Wt::TextFormat::XHTMLUnsafeText );
    m_iframe_holder->removeStyleClass( "InjaLogDialogIFrameHolder" );
    m_iframe_holder->addStyleClass( "InjaLogDialogInjaError" );

  } catch( std::exception &e )
  {
    const std::string error_msg = "<p>" + Wt::WString::tr( "ild-misc-error" ).arg( std::string( e.what() ) ).toUTF8() + "</p>";
    m_iframe_holder->setText( error_msg );
    m_iframe_holder->setTextFormat( Wt::TextFormat::XHTMLUnsafeText );
    m_iframe_holder->removeStyleClass( "InjaLogDialogIFrameHolder" );
    m_iframe_holder->addStyleClass( "InjaLogDialogMiscError" );
  }
}


void InjaLogDialog::handleTemplateChange()
{
  updateDisplay();
}
