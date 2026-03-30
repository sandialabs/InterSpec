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

#include <Wt/Utils.h>
#include <Wt/WLink.h>
#include <Wt/WText.h>
#include <Wt/WImage.h>
#include <Wt/WLabel.h>
#include <Wt/WAnchor.h>
#include <Wt/WCheckBox.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>
#include <Wt/WMemoryResource.h>
#include <Wt/WContainerWidget.h>

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/MultimediaDisplay.h"

using namespace std;
using namespace Wt;

#if( ANDROID )
// Defined in target/android/android.cpp
extern void android_download_workaround( Wt::WResource *resource, std::string description );
#endif

namespace
{

//Class to display a DetectorAnalysis; just a first go at it
//  Some sort of model and table, or something might be a better implementation.
class MultimediaDisplay : public WContainerWidget
{
  std::shared_ptr<const SpecMeas> m_meas;
  size_t m_current_index;
  shared_ptr<WMemoryResource> m_resource;
  WImage *m_image;
  WText *m_error;
  WContainerWidget *m_add_info;
  WText *m_remark;
  WText *m_description;
  WText *m_time;
  WContainerWidget *m_nav;
  WPushButton *m_prev;
  WText *m_pos_txt;
  WPushButton *m_next;
#if( BUILD_AS_OSX_APP || IOS )
  WAnchor *m_download;
#else
  WPushButton *m_download;
#endif
  
public:
  MultimediaDisplay()
  : WContainerWidget(),
  m_meas{ nullptr },
  m_current_index( 0 ),
  m_resource( nullptr ),
  m_image( nullptr ),
  m_error( nullptr ),
  m_add_info( nullptr ),
  m_remark( nullptr ),
  m_description( nullptr ),
  m_time( nullptr ),
  m_nav( nullptr ),
  m_prev( nullptr ),
  m_pos_txt( nullptr ),
  m_next( nullptr ),
  m_download( nullptr )
  {
    wApp->useStyleSheet( "InterSpec_resources/MultimediaDisplay.css" );
    
    InterSpec::instance()->useMessageResourceBundle( "MultimediaDisplay" );
    
    addStyleClass( "MultimediaDisplay" );
  
    m_resource = make_shared<WMemoryResource>();

    m_image = addNew<WImage>( WLink(m_resource) );
    m_image->addStyleClass( "SpecImage" );
    m_image->setHidden( true );

    m_error = addNew<WText>();
    m_error->setStyleClass( "NoContentTxt" );

    m_add_info = addNew<WContainerWidget>();
    m_add_info->addStyleClass( "AddInfo" );

    m_remark = m_add_info->addNew<WText>();
    m_remark->addStyleClass( WString::tr("Remark") );

    m_description = m_add_info->addNew<WText>();
    m_description->addStyleClass( WString::tr("Description") );

    m_time = m_add_info->addNew<WText>();
    m_time->addStyleClass( "Time" );

    m_nav = addNew<WContainerWidget>();
    m_nav->addStyleClass( "Nav" );

    m_prev = m_nav->addNew<WPushButton>( WString::tr("smmd-prev-btn") );
    m_prev->addStyleClass( "Prev" );
    m_prev->clicked().connect( this, &MultimediaDisplay::prevIndex );

    m_pos_txt = m_nav->addNew<WText>();
    m_pos_txt->addStyleClass( "NavPos" );

    m_next = m_nav->addNew<WPushButton>( WString::tr("smmd-next-btn") );
    m_next->addStyleClass( "Next" );
    m_next->clicked().connect( this, &MultimediaDisplay::nextIndex );


    WContainerWidget *footer = addNew<WContainerWidget>();
    footer->addStyleClass( "PrefAndDownload" );

    WCheckBox *cb = footer->addNew<WCheckBox>( WString::tr("smmd-auto-show-images-cb") );
    cb->setToolTip( WString::tr("smmd-tt-auto-show-images") );
    cb->addStyleClass( "CbNoLineBreak" );
    m_next->setFocus();
    
    InterSpec *interspec = InterSpec::instance();
    assert( interspec );
    if( interspec )
      UserPreferences::associateWidget( "AutoShowSpecMultimedia", cb, interspec );
    
#if( BUILD_AS_OSX_APP || IOS )
    m_download = footer->addNew<WAnchor>( WLink(m_resource) );
    m_download->setTarget( AnchorTarget::TargetNewWindow );
    m_download->setStyleClass( "LinkBtn DownloadLink" );
    m_download->setTarget( AnchorTarget::TargetNewWindow ); //TargetDownload
    m_download->setText( "Save..." );
#else
    m_download = footer->addNew<WPushButton>();
    m_download->setIcon( "InterSpec_resources/images/download_small.svg" );
    {
      WLink lnk( m_resource );
      lnk.setTarget( LinkTarget::NewWindow );
      m_download->setLink( lnk );
    }
    m_download->setStyleClass( "LinkBtn DownloadBtn" );

#if( ANDROID )
    // Using hacked saving to temporary file in Android, instead of via network download of file.
    m_download->clicked().connect( std::bind([this](){
      android_download_workaround(m_resource.get(), "image_from_spec_file");
    }) );
#endif //ANDROID
#endif
    
    m_download->setToolTip( WString::tr("smmd-tt-export-image-file"), Wt::TextFormat::Plain );
    
    setIndex( 0 );
  }//constructor
  
  void nextIndex()
  {
    if( !m_meas || (m_meas->multimedia_data().size() < 2) )
    {
      setIndex( 0 );
      return;
    }
    
    const size_t index = ((m_current_index + 1) % m_meas->multimedia_data().size());
    setIndex( index );
    m_next->setFocus();
  }//void nextIndex()
  
  
  void prevIndex()
  {
    if( !m_meas || (m_meas->multimedia_data().size() < 2) )
    {
      setIndex( 0 );
      return;
    }
    
    size_t index = m_current_index > 0 ? m_current_index - 1 : m_meas->multimedia_data().size() - 1;
    setIndex( index );
    m_prev->setFocus();
  }//void prevIndex()
  
  
  void setIndex( size_t index )
  {
    vector<shared_ptr<const SpecUtils::MultimediaData>> all_data;
    if( m_meas )
      all_data = m_meas->multimedia_data();
    
    auto showErrorMsg = [this]( const WString &msg ){
      m_image->hide();
      m_resource->setData( {} );
      m_add_info->hide();
      m_nav->hide();
      m_error->show();
      m_error->setText( msg );
      m_download->hide();
      m_nav->setHidden( !m_meas || (m_meas->multimedia_data().size() < 2) );
    };//showErrorMsg lamda
    
    if( all_data.empty() )
    {
      showErrorMsg( WString::tr("smmd-err-no-content") );
      return;
    }//if( all_data.empty() )
  
  
    index = (index % all_data.size());
    m_current_index = index;
    
    shared_ptr<const SpecUtils::MultimediaData> data = all_data[index];
    
    if( !data || (data->data_.size() < 25) )
    {
      if( data && data->file_uri_.empty() )
        showErrorMsg(  WString::tr("smmd-err-file-is-uri").arg(data->file_uri_) );
      else
        showErrorMsg( WString::tr("smmd-err-no-multimedia-src") );
      
      return;
    }//if( !data || (data->data_.size() < 25) )
    
    string data_str( begin(data->data_), end(data->data_) );
      
    switch( data->data_encoding_ )
    {
      case SpecUtils::MultimediaData::EncodingType::BinaryUTF8:
        data_str.clear();
        break;
          
      case SpecUtils::MultimediaData::EncodingType::BinaryHex:
        data_str = Wt::Utils::hexDecode(data_str);
        break;
          
      case SpecUtils::MultimediaData::EncodingType::BinaryBase64:
        data_str = Wt::Utils::base64Decode(data_str);
        break;
    }//switch( data->data_encoding_ )
      
    vector<unsigned char> data_uchar( (const unsigned char *)data_str.c_str(),
                                      (const unsigned char *)(data_str.c_str() + data_str.size()) );
      
    const string mime = Wt::Utils::guessImageMimeTypeData(data_uchar);
      
    if( data_str.empty() )
    {
      showErrorMsg( WString::tr("smmd-err-encoding-not-supported") );
      return;
    }
      
    if( mime.empty() )
    {
      showErrorMsg( WString::tr("smmd-err-not-image") );
      return;
    }
      
    m_resource->setData( data_uchar );
    m_resource->setMimeType( mime );
    
    string filename = "image_" + to_string(m_current_index + 1) + "_of_" + to_string(all_data.size());
    
    if( SpecUtils::icontains(mime, "png") )
      filename += ".png";
    else if( SpecUtils::icontains(mime, "jpeg") )
      filename += ".jpeg";
    else if( SpecUtils::icontains(mime, "gif") )
      filename += ".gif";
    else if( SpecUtils::icontains(mime, "bmp") )
      filename += ".bmp";
    m_resource->suggestFileName( filename );
     
    m_image->show();
    m_error->hide();
    m_nav->setHidden( (all_data.size() < 2) );
    m_download->show();
    
    WString postxt = WString::tr("smmd-image-index")
                      .arg( static_cast<int>(m_current_index + 1) )
                      .arg( static_cast<int>(all_data.size()) );
    m_pos_txt->setText( postxt );
    
    
    bool have_add_info = false;
    if( SpecUtils::is_special( data->capture_start_time_ ) )
    {
      m_time->setText( "" );
      m_time->hide();
    }else
    {
      have_add_info = true;
      
      WString msg = WString::tr("smmd-time-label")
                     .arg( SpecUtils::to_common_string(data->capture_start_time_, true) );
      m_time->setText( msg );
      m_time->show();
    }
    
    if( data->remark_.empty() )
    {
      m_remark->hide();
      m_remark->setText( "" );
    }else
    {
      have_add_info = true;
      m_remark->show();
      WString msg = WString::tr("smmd-remark-label")
                     .arg( data->remark_ );
      m_remark->setText( msg );
    }
    
    if( data->descriptions_.empty() )
    {
      m_description->hide();
      m_description->setText( "" );
    }else
    {
      have_add_info = true;
      m_description->show();
      WString msg = WString::tr("smmd-desc-label")
                     .arg( data->descriptions_ );
      m_description->setText( msg );
    }
    
    m_add_info->setHidden( !have_add_info );
    
    // Just in case the image changes the dialogs size, trigger it to resize.
    wApp->doJavaScript( wApp->javaScriptClass() + ".TriggerResizeEvent();" );
  }//void setIndex( size_t index )
  
  
  void updateDisplay( std::shared_ptr<const SpecMeas> meas )
  {
    m_meas = meas;
    m_current_index = 0;
    setIndex( m_current_index );
  }//void updateDisplay( std::shared_ptr<const SpecMeas> meas )

 
};//class AnaResultDisplay
}//namespace


SimpleDialog *displayMultimedia( const std::shared_ptr<const SpecMeas> &spec )
{
  wApp->useStyleSheet( "InterSpec_resources/MultimediaDisplay.css" );
  
  
  auto dialog = SimpleDialog::make();
  //dialog->setModal( false ); //doesnt seem to have any effect
  dialog->addButton( WString::tr("Close") );
  
  WContainerWidget *contents = dialog->contents();
  const bool multiple_images = (spec && (spec->multimedia_data().size() > 1));
  const char *title_key = multiple_images ? "smmd-window-title-multiple" : "smmd-window-title-single";
  // I think its find that we may not have read MultimediaDisplay.xml yet - I dont think WString resolves the keys immediately
  WText *dialogTitle = contents->addNew<WText>( WString::tr(title_key) );
  dialogTitle->addStyleClass( "title MultimediaDialogTitle" );
  dialogTitle->setInline( false );

  MultimediaDisplay *display = contents->addNew<MultimediaDisplay>();
  display->updateDisplay( spec );
  
  return dialog;
}//SimpleDialog *displayMultimedia( const std::shared_ptr<const SpecMeas> &spec )
