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

#include "rapidxml/rapidxml.hpp"

#include <Wt/WText>
#include <Wt/Utils>
#include <Wt/WLink>
#include <Wt/WMenu>
#include <Wt/WImage>
#include <Wt/WAnchor>
#include <Wt/WString>
#include <Wt/Dbo/ptr>
#include <Wt/WTemplate>
#include <Wt/WCheckBox>
#include <Wt/WMenuItem>
#include <Wt/Json/Array>
#include <Wt/WAnimation>
#include <Wt/WFileUpload>
#include <Wt/WTabWidget>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/Json/Parser>
#include <Wt/Json/Object>
#include <Wt/WApplication>
#include <Wt/WEnvironment>
#include <Wt/WMediaPlayer>
#include <Wt/WRadioButton>
#include <Wt/WStandardItem>
#include <Wt/WStackedWidget>
#include <Wt/WContainerWidget>
#include <Wt/WStandardItemModel>
#include <Wt/WAbstractItemDelegate>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/UseInfoWindow.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/DbFileBrowser.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/RowStretchTreeView.h"


using namespace Wt;
using namespace std;


namespace
{
  class SampleSpecRow : public Wt::WContainerWidget
  {
  public:
    SampleSpecRow( WContainerWidget *parent = 0 )
    : WContainerWidget( parent )
    {
      m_text = new WText( this );
      m_load = new WPushButton( WString::tr("uiw-load-btn"), this );
      
      // We need this, or else by the time we get to loadSampleSelected(), no row will be
      this->touchStarted().preventPropagation();
      
      addStyleClass( "SampleSpecRow" );
      m_text->addStyleClass( "SampleSpecText" );
      m_load->addStyleClass( "SampleSpecLoadBtn" );
    }
    
    virtual ~SampleSpecRow(){};
    
    WText *m_text;
    WPushButton *m_load;
  };//class SampleSpecRow
  
  class SampleSpectrumDelegate : public Wt::WAbstractItemDelegate
  {
    UseInfoWindow *m_parent;
    
  public:
    SampleSpectrumDelegate( UseInfoWindow *parent )
    : Wt::WAbstractItemDelegate( parent ),
      m_parent( parent )
    {
    }
    
    virtual ~SampleSpectrumDelegate(){};
    
    virtual WWidget *update(WWidget *widget, const WModelIndex& index,
                            WFlags<ViewItemRenderFlag> flags)
    {
      SampleSpecRow *w = dynamic_cast<SampleSpecRow *>(widget);
      if( !w )
      {
        w = new SampleSpecRow();
        w->m_load->clicked().connect( m_parent, &UseInfoWindow::loadSampleSelected );
      }
      
      try
      {
        w->m_text->setText( asString(index.data()) );
        
        //It appears that when you select an item, this function does not get called - unexpected.
        //w->m_load->setHidden( !flags.testFlag(Wt::ViewItemRenderFlag::RenderSelected) );
        //if( flags.testFlag(Wt::ViewItemRenderFlag::RenderFocused) )
        //if( flags.testFlag(Wt::ViewItemRenderFlag::RenderInvalid) )
      }catch( std::exception &e )
      {
        cerr << "SampleSpectrumDelegate: caught exception converting to string: " << e.what() << endl;
        w->m_text->setText( "" );
      }
      
      return w;
    }//virtual WWidget *update(...)
  };//class SampleSpectrumDelegate
  
}//namespace


UseInfoWindow::UseInfoWindow( std::function<void(bool)> showAgainCallback,
                              InterSpec* viewer )
: AuxWindow( "", (AuxWindowProperties::IsModal
                         | AuxWindowProperties::EnableResize
                         | AuxWindowProperties::DisableCollapse
                         | AuxWindowProperties::SetCloseable) ),
  m_session(),
#if( USE_DB_TO_STORE_SPECTRA )
  m_snapshotModel( nullptr ),
#endif
  m_tableSample( nullptr ),
  m_messageModelSample( nullptr ),
  m_viewer( viewer ),
  m_menu( nullptr )
  //, m_videoTab( nullptr )
{
  m_viewer->useMessageResourceBundle( "UseInfoWindow" );
  
  rejectWhenEscapePressed();
  
  WStackedWidget *stack = new WStackedWidget();
  stack->addStyleClass( "UseInfoStack" );
  //stack->setOverflow( WContainerWidget::OverflowAuto );
  
  WAnimation animation(Wt::WAnimation::Fade, Wt::WAnimation::Linear, 200);
  stack->setTransitionAnimation( animation, true );
  
  m_menu = new WMenu( stack, Wt::Vertical );
  m_menu->addStyleClass( (m_viewer->isPhone() ? "VerticalNavMenuPhone HeavyNavMenuPhone SideMenuPhone" : "VerticalNavMenu HeavyNavMenu SideMenu") );
  
  WDialog::contents()->setOverflow( WContainerWidget::OverflowHidden );
  
  WGridLayout *layout = stretcher();
  layout->setContentsMargins( 9, 0, 9, 0 );
  layout->addWidget( m_menu, 0, 0, AlignLeft );
  const int nstack_rows = m_viewer->isMobile() ? 1 : 3;  //"License and Terms" and "More in depth information" links for non-mobile users
  layout->addWidget( stack, 0, 1, nstack_rows, 1 );
  layout->setRowStretch( 0, 1 );
  
  const string docroot = wApp->docRoot();
  const string bundle_file = SpecUtils::append_path(docroot, "InterSpec_resources/static_text/use_instructions" );
  m_resourceBundle.use(bundle_file,false);
  
  initVideoInfoMap();
  WContainerWidget* welcomeContainer = new WContainerWidget();
  //welcomeContainer->setOverflow(WContainerWidget::OverflowAuto);
  welcomeContainer->setOffsets(WLength(0,WLength::Pixel));
  welcomeContainer->setMargin(WLength(0,WLength::Pixel));
  
  Wt::WTabWidget *tabW = new Wt::WTabWidget();
  tabW->addStyleClass( "UseInfoTabsWidget" );
  tabW->contentsStack()->addStyleClass( "UseInfoTabsStack" );
  
  WContainerWidget *spectrumContainer = new WContainerWidget();
  WGridLayout *spectrumLayout = new WGridLayout();
  spectrumLayout->setContentsMargins( 3, 5, 4, 5);
  spectrumContainer->setLayout(spectrumLayout);
  spectrumContainer->setOverflow( WContainerWidget::OverflowHidden );
  
  WContainerWidget *samplesContainer = new WContainerWidget();
  WGridLayout *samplesLayout = new WGridLayout();
  samplesLayout->setContentsMargins( 0, 0, 0, 0 );
  samplesContainer->setLayout( samplesLayout );
  samplesContainer->setOverflow( WContainerWidget::OverflowAuto );
  
  WContainerWidget *importContainer = nullptr;
  WGridLayout *importLayout = nullptr;
  
  //The below tabs must use WTabWidget::PreLoading, or else we will get a
  //  JavaScript exception when showLoadingIndicator() is called (in the JS).
  //  This is super confusing, and I dont understand it.  This manifests
  //  particularly on Android native version of app, but I think it might
  //  elsewhere as well.  (wcjohns 20150124)
  tabW->addTab( spectrumContainer, WString::tr("uiw-tab-saved-states"), WTabWidget::PreLoading );
  tabW->addTab( samplesContainer, WString::tr("uiw-tab-example-spectra"), WTabWidget::PreLoading ); 
  
  if( m_viewer->isMobile() )
  {
    importContainer = new WContainerWidget();
    importLayout = new WGridLayout();
    importContainer->setLayout( importLayout );
    
    tabW->addTab( importContainer, WString::tr("uiw-tab-import-spectra"), WTabWidget::PreLoading );
  }//if( m_viewer->isMobile() )
  
  
  //----------Add workspace tab --------
  
  //We have to create a independant Dbo::Session for this class since the
  //  viewer->user.session() could be used in other threads, messing
  //  things up (Dbo::Session is not thread safe).
  m_session.reset( new DataBaseUtils::DbSession( *m_viewer->sql() ) );

  try
  {
#if( USE_DB_TO_STORE_SPECTRA )
    const Dbo::ptr<InterSpecUser> &user = m_viewer->user();
    SpecMeasManager *manager = m_viewer->fileManager();
    
    m_snapshotModel = new SnapshotBrowser( manager, m_viewer, nullptr );
    spectrumLayout->addWidget( m_snapshotModel, 0, 0 );
#endif
    
    const string docroot = wApp->docRoot();
    
    ///-----Sample -----
    m_messageModelSample = new Wt::WStandardItemModel( 0, 3, this );
    
    WStandardItem *parserType = new WStandardItem();
    parserType->setData(SpecUtils::ParserType::N42_2006); //SpeIaea
    vector<WStandardItem *> message( 3 );
    message[0] = new Wt::WStandardItem("Ba-133 (16k bin N42)");
    message[1] = new Wt::WStandardItem( SpecUtils::append_path(docroot, "example_spectra/ba133_source_640s_20100317.n42") );
    message[2] = parserType;
    m_messageModelSample->appendRow(message);
    
    
    parserType = new WStandardItem();
    parserType->setData(SpecUtils::ParserType::N42_2006); //SpeIaea
    message[0] = new Wt::WStandardItem("Passthrough (16k bin, 8 det, 1064 meas N42)");
    message[1] = new Wt::WStandardItem( SpecUtils::append_path(docroot, "example_spectra/passthrough.n42" ));
    message[2] = parserType;
    m_messageModelSample->appendRow(message);
  
    
    parserType = new WStandardItem();
    parserType->setData(SpecUtils::ParserType::N42_2006); //SpeIaea
    message[0] = new Wt::WStandardItem("Background (16k bin N42)");
    message[1] = new Wt::WStandardItem( SpecUtils::append_path(docroot, "example_spectra/background_20100317.n42") );
    message[2] = parserType;
    m_messageModelSample->appendRow(message);
    
    const pair<string, string> ref_spectra[7] = {
      {"Shielded Ir192", "Ir192_Shielded.txt"},
      {"Shielded Pu239", "Pu239_Shielded_100g.txt"},
      {"Tl201 and 202 in phantom", "Tl201wTl202_Phantom.txt"},
      {"2% Enriched Uranium", "U235_Unshielded_0200.txt"},
      {"90% Enriched Uranium", "U235_Unshielded_9000.txt"},
      {"Uranium ore","Uore_Unshielded.txt"},
      {"Shielded Np237","Np237_Shielded.txt"}
    };

    const pair<string, string> ref_detectors[] = {
      {"HPGe", "Common_Field_Nuclides/Fulcrum-40h"},
      {"NaI", "Common_Field_Nuclides/Verifinder-NaI"},
      {"CZT", "Common_Field_Nuclides/CZT_H3D_M400"}
    };

    const string static_data_dir = InterSpec::staticDataDirectory();
    const string ref_spectra_dir = SpecUtils::append_path( static_data_dir, "reference_spectra" );

    for( const pair<string, string> &descript_filename : ref_spectra )
    {
      const string &description = descript_filename.first;
      const string &filename = descript_filename.second;

      for( const pair<string, string> dettype_path : ref_detectors )
      {
        const string &det_type = dettype_path.first;
        const string path = SpecUtils::append_path( ref_spectra_dir, dettype_path.second );

        const string descrip = description + " (" + det_type + ")";
        const string fullpath = SpecUtils::append_path( path, filename );

        parserType = new WStandardItem();
        parserType->setData(SpecUtils::ParserType::Auto);
        message[0] = new Wt::WStandardItem(descrip);
        message[1] = new Wt::WStandardItem(fullpath);
        message[2] = parserType;
        m_messageModelSample->appendRow(message);
      }//for( const auto &[detector, path] : ref_detectors )
    }//for( const auto &[description, filename] : ref_spectra )


    //m_messageModelSample->setHeaderData(  0, Horizontal, WString("Example Spectra"), DisplayRole );
    m_tableSample = new RowStretchTreeView();
    m_tableSample->setRootIsDecorated( false ); //makes the tree look like a table! :)
    m_tableSample->setHeaderHeight( 0 );
    m_tableSample->addStyleClass( "DbSpecFileSelectTable" );
    m_tableSample->setModel( m_messageModelSample );
    m_tableSample->setOffsets( WLength(0, WLength::Pixel), Wt::Left | Wt::Top );
    m_tableSample->setColumnResizeEnabled(true);
    m_tableSample->setRowHeight( WLength(30,WLength::Pixel) );
    m_tableSample->setColumnAlignment(0, Wt::AlignLeft);
    //m_tableSample->setHeaderAlignment(0, Wt::AlignCenter);
    if( m_viewer->isPhone() )
      m_tableSample->setColumnWidth(0,240);
    else
      m_tableSample->setColumnWidth(0,340);
    
    m_tableSample->hideColumn(1);
    m_tableSample->hideColumn(2);
    
    //if( m_viewer->isPhone() )
    //  m_tableSample->setMaximumSize( WLength::Auto, 120 );
    
    m_tableSample->setAlternatingRowColors(true);
    m_tableSample->setSelectionMode( Wt::SingleSelection );
    m_tableSample->setEditTriggers( Wt::WAbstractItemView::NoEditTrigger );
    
    auto sampleDelegate = new SampleSpectrumDelegate( this );
    m_tableSample->setItemDelegate( sampleDelegate );
    
    samplesLayout->addWidget(m_tableSample,0,0);
    samplesLayout->setRowStretch(0,1);
    samplesLayout->setColumnStretch(0,1);
    
    m_tableSample->selectionChanged().connect( this, &UseInfoWindow::handleSampleSelectionChanged );
    m_tableSample->doubleClicked().connect( boost::bind( &UseInfoWindow::handleSampleDoubleClicked,
                                                        this, boost::placeholders::_1,
                                                        boost::placeholders::_2 ) );
    
    
    //Import tab
    if( importContainer )
    {
      string msg = WString::tr("uiw-tap-choose-file").toUTF8();
      WText *desc = new WText( msg );
      desc->setTextAlignment( Wt::AlignmentFlag::AlignCenter );
      importLayout->addWidget( desc, 0, 0, AlignCenter );
      
      WFileUpload *upload = new WFileUpload();
      importLayout->addWidget( upload, 1, 0, AlignCenter | AlignMiddle );
      
      upload->fileTooLarge().connect( std::bind( [desc,upload](){
        upload->hide();
        desc->setText( WString::tr("uiw-file-too-large") );
      }) );
      
      upload->uploaded().connect( std::bind( [this,upload](){
        const string filename = upload->spoolFileName();
        const string displayFileName = upload->clientFileName().toUTF8();
        m_viewer->userOpenFileFromFilesystem( filename, displayFileName );
      } ) );
                                            
      upload->changed().connect( upload, &WFileUpload::upload );
      
      
      
      msg = WString::tr("uiw-open-file-instructions").toUTF8();
      desc = new WText( msg );
      importLayout->addWidget( desc, 2, 0, AlignBottom );
    }//if( importContainer )
    
      
    const int sampleIndex = tabW->indexOf(samplesContainer);
    
#if( !USE_DB_TO_STORE_SPECTRA )
    tabW->setCurrentIndex( sampleIndex );
#else
    if( !m_snapshotModel->numSnaphots() )
      tabW->setCurrentIndex( sampleIndex );
#endif
  } //try
  catch( exception &e )
  {
    cerr << "Exception:" << e.what() << endl;
  }//catch
  
  
  //----------Add left menu --------
  SideMenuItem *item = new SideMenuItem( (m_viewer->isPhone() ? WString() : WString::tr("uiw-welcome")), welcomeContainer );
  if( m_viewer->isPhone() )
  {
    //item->setIcon( "InterSpec_resources/images/home_small.svg" );
    //item->addStyleClass( "WhiteIcon" );
    WAnchor *anchor = item->anchor();
    WImage *image = new WImage( "InterSpec_resources/images/home_small.svg" );
    image->addStyleClass( "WhiteIcon" );
    anchor->addWidget( image );
  }
  
  item->clicked().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
  item->mouseWentDown().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
  m_menu->addItem( item );
  
  
  // ---- Videos
  /*
   //Dont put the videos in for now - they arent that good - but leaving code
   //  commented out because the mechanism still works for in the future if
   //  I do improve the videos.
  {
    WContainerWidget* controlContainer = new WContainerWidget();
    controlContainer->setOverflow(WContainerWidget::OverflowHidden);
    controlContainer->setOffsets(WLength(0,WLength::Pixel));
    controlContainer->setMargin(WLength(0,WLength::Pixel));
    WGridLayout* controlLayout = new WGridLayout();
    controlLayout->setContentsMargins(0, 0, 0, 0);
    
    controlContainer->setLayout(controlLayout);
    
    m_videoTab = new Wt::WTabWidget();
      const WBorder border( WBorder::Solid, WBorder::Thin, Wt::gray);
      m_videoTab->contentsStack()->decorationStyle().setBorder( border,  Wt::Bottom | Wt::Right | Wt::Left );
      m_videoTab->setMargin(2);
      m_videoTab->setOffsets(2);

    controlLayout->addWidget(m_videoTab,1,0);
    controlLayout->setRowStretch(0,1);
    
    
    std::string controlInfo="";
    m_resourceBundle.resolveKey("videoinfo",controlInfo);
    
    try
    {
      Json::Object result;
      Json::parse( controlInfo, result );
      
      for( Wt::Json::Object::const_iterator iter = result.begin();
          iter != result.end(); ++iter )
      {
        const Json::Array &controls = iter->second;
        for( size_t i = controls.size(); i >0 ; i-- )
        {
          const Json::Object &info = controls[i-1];
          const Json::Value &title = info.get("title");
          
          if( title.isNull() )
          {
            cerr << "makeItem found a null key or file string, "
            "you may want to check use_instructions.xml" << endl;
            continue;
          }//if( binding.isNull() || file.isNull() )
          
          const WString titlestr = title;
          WContainerWidget* samplesContainer = new WContainerWidget();
          samplesContainer->setMargin(15,Wt::Top);
          samplesContainer->setOffsets(15,Wt::Top);
          samplesContainer->setStyleClass("centered");
          itemCreator(iter->first, samplesContainer, titlestr );
          WMenuItem * tab= m_videoTab->addTab(samplesContainer, titlestr, Wt::WTabWidget::LazyLoading);
          tab->clicked().connect( boost::bind( &UseInfoWindow::tab_select_item, this, tab) );
          tab->mouseWentDown().connect( boost::bind( &UseInfoWindow::tab_select_item, this, tab) );
          
        }//for( size_t i = 0; i < videos.size(); ++i )
      }//for( loop over sesults )
    }catch( std::exception &e )
    {
      cerr << "getVideoInfoMap() caught: " << e.what() << endl
      << "you may want to check use_instructions.xml" << endl;
    }//try / catch
   
    item = new SideMenuItem( m_viewer->isMobile()?"":"Tutorial videos", controlContainer );
    item->setId("TutorialVideos");
    item->clicked().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
    item->mouseWentDown().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
    //item->setIcon("InterSpec_resources/images/film.png");
    m_menu->addItem( item );
  } //videos
*/
  
  
  //----------Create Controls window--------
  
  //Select either mobile or desktop mouse/keyboard/gesture page
  if( m_viewer->isMobile() )
  {
    item = makeItem( m_viewer->isPhone() ? "" : "Use Info", "mobile-mouse-interactions" );
    //item->setIcon("InterSpec_resources/images/phone_tips.svg");
    //item->addStyleClass( "WhiteIcon" );
    
    if( m_viewer->isPhone() )
    {
      WAnchor *anchor = item->anchor();
      WImage *image = new WImage( "InterSpec_resources/images/phone_tips.svg" );
      image->addStyleClass( "WhiteIcon" );
      anchor->addWidget( image );
    }
  }else if( m_viewer->isDesktop() )
  {
    WContainerWidget* controlContainer = new WContainerWidget();
    controlContainer->setOverflow(WContainerWidget::OverflowHidden);
    WGridLayout* controlLayout = new WGridLayout();
    controlLayout->setContentsMargins(0, 0, 0, 0);
    controlContainer->setLayout(controlLayout);
    
    Wt::WTabWidget *controlTab = new Wt::WTabWidget();
      
    const WBorder border( WBorder::Solid, WBorder::Thin, Wt::gray);
    controlTab->contentsStack()
              ->decorationStyle()
              .setBorder( border,  Wt::Bottom | Wt::Right | Wt::Left );
    controlTab->setMargin( 2 );
    controlTab->setOffsets( 2 );

    //  tabW->setMinimumSize(300, 300);
    controlLayout->addWidget(controlTab,1,0);
    controlLayout->setRowStretch(0,1);
    
    std::string controlInfo = "";
    m_resourceBundle.resolveKey( "desktop-mouse-interactions-info", controlInfo );
    
    try
    {
      Json::Object result;
      Json::parse( controlInfo, result , true);
      
      for( Wt::Json::Object::const_iterator iter = result.begin();
          iter != result.end(); ++iter )
      {
        const Json::Object &info = iter->second;
        const Json::Value &title = info.get( "title" );
        const Json::Value &image = info.get( "img" );
        const Json::Value &desc_id = info.get( "desc_message_id" );
          
        if( title.isNull() && image.isNull() )
        {
          cerr << "makeItem found a null key or file string, "
                  "you may want to check use_instructions.xml" << endl;
          continue;
        }//if( binding.isNull() || file.isNull() )
          
        const WString titlestr = title.orIfNull( "" );
        const WString imagestr = image.orIfNull( "" );
        string descstr = "";
          
        if( !desc_id.isNull() )
          m_resourceBundle.resolveKey( desc_id, descstr );
        
        WContainerWidget *samplesContainer = new WContainerWidget();
        samplesContainer->addStyleClass("WelcomeControlsTabContents");
        
        if( !imagestr.empty() )
        {
          WImage* img = new WImage( imagestr.toUTF8() );
          img->setMaximumSize(WLength(300), WLength::Auto);
          samplesContainer->addWidget(img);
          img->setStyleClass("img-centered");
        }//if( !imagestr.empty() )
        
        if( !descstr.empty() )
        {
          WText* txt = new WText( descstr, XHTMLUnsafeText );
          txt->setStyleClass("img-centered");
          samplesContainer->addWidget(txt);
        } //descstr.empty()
        
        controlTab->addTab(samplesContainer, titlestr, Wt::WTabWidget::LazyLoading);
      }//for( loop over results )
    }catch( std::exception &e )
    {
      cerr << "processing control info map caught: " << e.what() << endl
      << "you may want to check use_instructions.xml" << endl;
    }//try / catch
    
    
    item = new SideMenuItem( WString::tr("uiw-controls"), controlContainer );
    item->clicked().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
    item->mouseWentDown().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
    m_menu->addItem( item );
  } //isDesktop()

  
  const string tips_file = "InterSpec_resources/static_text/faqs/faqs";
  SideMenuItem *faqsItem = makeTextItem( "FAQs", tips_file );
  
  // ---- Cheat sheet PDF
/*
 //Cheat sheet is currently outdated, so dont display it.
  if( !m_viewer->isMobile() )
  {
    WContainerWidget* controlContainer = new WContainerWidget();
    controlContainer->setOverflow(WContainerWidget::OverflowHidden);
    controlContainer->setOffsets(WLength(0,WLength::Pixel));
    controlContainer->setMargin(WLength(0,WLength::Pixel));
    WGridLayout* controlLayout = new WGridLayout();
    controlLayout->setContentsMargins(0, 0, 0, 0);
    
    controlContainer->setLayout(controlLayout);
    WAnchor* anchor = new WAnchor();
    anchor->setLink(WLink("InterSpec_resources/static_text/InteractionCheatSheet.pdf"));
    //anchor->setImage(new WImage("InterSpec_resources/images/pdf_page.png"));
    anchor->setTarget( Wt::TargetNewWindow );
    anchor->setText(WString::tr("uiw-cheat-sheet"));

    auto cheatTxt = new WText(WString::tr("uiw-cheat-sheet-desc"));
    controlLayout->addWidget(cheatTxt,0,0);
    controlLayout->addWidget(anchor,1,0);
    controlLayout->setRowStretch(1, 1);
    
    item = new SideMenuItem( (m_viewer->isPhone()?"":"Cheat Sheet"), controlContainer );
    item->clicked().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
    item->mouseWentDown().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
    if( m_viewer->isMobile() )
    {
      WImage *image = new WImage( "InterSpec_resources/images/pdf_page.png" );
      image->addStyleClass( "WhiteIcon" );
      item->anchor()->addWidget( image );
    }
    m_menu->addItem( item );
  } //Cheat sheet PDF
*/
  
  string msg;

  if( m_viewer->isMobile() )
  {
    msg += "<div align=\"left\">";
#if( USE_DB_TO_STORE_SPECTRA )
    if( m_snapshotModel && !m_snapshotModel->numSnaphots() )
      msg += WString::tr("uiw-mobile-start-example").toUTF8();
    else
      msg += WString::tr("uiw-mobile-saved-spectra").toUTF8();
#else
    msg += WString::tr("uiw-mobile-select-example").toUTF8();
#endif
   
  }else
  {
    msg += WString::tr("uiw-desktop-welcome-title").toUTF8();
    msg += WString::tr("uiw-desktop-welcome-text").toUTF8();
    msg += WString::tr("uiw-desktop-drag-drop").toUTF8();
  }

  msg += "</div>";
  
  WText *text = new WText( msg );
  
  WString mobilePostcriptText = WString::tr("uiw-mobile-phone-postscript");
  if( m_viewer->isTablet() )
    mobilePostcriptText = WString::tr("uiw-mobile-tablet-postscript");
  
  WGridLayout *rightLayout = nullptr;
  if( m_viewer->isPhone() )
  {
    rightLayout = new WGridLayout();
    rightLayout->setContentsMargins( 0, 0, 0, 0 );
    rightLayout->setVerticalSpacing( 0 );
    rightLayout->setHorizontalSpacing( 0 );
    welcomeContainer->setLayout(rightLayout);
    
    text->setInline( false );
    text->addStyleClass( "MobileWelcomeText" );
    tabW->addStyleClass( "UseInfoTabsMobile" );
    //welcomeContainer->addWidget( text );
    //welcomeContainer->addWidget( tabW );
    rightLayout->addWidget( text, 0, 0 );
    rightLayout->addWidget( tabW, 1, 0 );
    rightLayout->setRowStretch( 1, 1 );
    
    text = new WText( mobilePostcriptText, Wt::XHTMLUnsafeText );
    text->setInline( false );
    text->addStyleClass( "MobileWelcomePoscriptText" );
    //welcomeContainer->addWidget( text );
    rightLayout->addWidget( text, 2, 0 );
    
    faqsItem->addStyleClass( "PhoneSideMenuTextItem" );
    
    /*
     //"Spectra attached to emails can be copied into <code>InterSpec</code> from within your mail app.<br />"
     //"Once you have analyzed a spectrum, you can save a snapshot of your work, and it will show up in the <em>Saved States</em> tab."
     */
  }else
  {
    rightLayout = new WGridLayout();
    welcomeContainer->setLayout(rightLayout);
    
    rightLayout->addWidget( text, 0, 0 );
    rightLayout->addWidget( tabW, 1, 0 );
    rightLayout->setRowStretch( 1, 1 );
    rightLayout->setContentsMargins( 9, 1, 9, 1 );
    
     if( m_viewer->isMobile() )
     {
       text = new WText( mobilePostcriptText, Wt::XHTMLUnsafeText );
       text->addStyleClass( "MobileWelcomePoscriptText" );
       rightLayout->addWidget( text, 2, 0 );
     }//if( m_viewer->isMobile() )
  }


  m_menu->select( 0 );
  
  if( !m_viewer->isMobile() )
    layout->setRowStretch( 0, 1 );
  layout->setColumnStretch( 1, 1 );

  
  WContainerWidget *bottom = footer();

  if( m_viewer->isMobile() )
  {
    AuxWindow::addHelpInFooter(bottom, "");
  }else
  {
    // Note that clicking on either link below deletes this window, and creates a new window, which
    //  will be captured as two seperate undo steps, even though it is a single user action.
    //  Not worth the effort to fix this, atm.
    
    //Add link for license and terms
    WAnchor *more = new WAnchor( WLink(), WString::tr("uiw-license-and-terms") );
    
    auto showLicenceAndTerms = [](){
      // We'll keep from capturing closing the window and opening the new window as two
      //  separate undo/redo events.
      {
        UndoRedoManager::BlockUndoRedoInserts blocker;
        InterSpec *viewer = InterSpec::instance();
        viewer->deleteWelcomeDialog(false);
        viewer->showLicenseAndDisclaimersWindow();
      }
      
      UndoRedoManager *undoRedo = UndoRedoManager::instance();
      if( undoRedo && undoRedo->canAddUndoRedoNow() )
      {
        auto undo = [](){
          InterSpec *viewer = InterSpec::instance();
          viewer->deleteLicenseAndDisclaimersWindow();
          viewer->showWelcomeDialogWorker( true );
        };
        auto redo = [](){
          InterSpec *viewer = InterSpec::instance();
          viewer->deleteWelcomeDialog(false);
          viewer->showLicenseAndDisclaimersWindow();
        };
        undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Show License and Terms" );
      }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
    };//showLicenceAndTerms
    
    more->clicked().connect( std::bind(showLicenceAndTerms) );
    more->addStyleClass("InfoLink");
    layout->addWidget( more, 1, 0 );
    
    //regular locations
    more = new WAnchor( WLink(), WString::tr("uiw-more-info") );
    
    auto showMoreInDepth = [](){
      // We'll keep from capturing closing the window and opening the new window as two
      //  seperate undo/redo events.  And instead make them a single step below
      {
        UndoRedoManager::BlockUndoRedoInserts blocker;
        InterSpec *viewer = InterSpec::instance();
        viewer->deleteWelcomeDialog( false );
        viewer->showHelpWindow( "getting-started" );
      }
      
      UndoRedoManager *undoRedo = UndoRedoManager::instance();
      if( undoRedo && undoRedo->canAddUndoRedoNow() )
      {
        auto undo = [](){
          InterSpec *viewer = InterSpec::instance();
          viewer->closeHelpWindow();
          viewer->showWelcomeDialogWorker( true );
        };
        auto redo = [](){
          InterSpec *viewer = InterSpec::instance();
          viewer->deleteWelcomeDialog( false );
          viewer->showHelpWindow( "getting-started" );
        };
        undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Show help window." );
      }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
    };//auto showMoreInDepth
    
    more->clicked().connect( std::bind(showMoreInDepth) );
    more->addStyleClass("InfoLink");
    layout->addWidget( more, 2, 0 );
  }
  
  
  layout->setContentsMargins( 0, 9, 9, 0 );
    
    
  if( showAgainCallback )
  {
    //WGridLayout *rightLayout = nullptr;
    bool isVertical = false;
    if( m_viewer->isPhone() )
    {
      int w = m_viewer->renderedWidth();
      int h = m_viewer->renderedHeight();
      if( w < 100 )
      {
        w = wApp->environment().screenWidth();
        h = wApp->environment().screenHeight();
      }
      isVertical = ((w > 100) && (w < h));
    }//if( isPhone )
    
    
    WCheckBox *showAgainCb = nullptr;
    assert( rightLayout );
    if( isVertical && rightLayout )
    {
      showAgainCb = new WCheckBox( WString::tr("uiw-show-at-start-cb") );
      showAgainCb->addStyleClass( "CbNoLineBreak" );
      rightLayout->addWidget( showAgainCb, rightLayout->rowCount(), 0 );
    }else
    {
      // Put CB into header
      showAgainCb = new WCheckBox( WString::tr("uiw-show-at-start-cb"), bottom );
    }
    showAgainCb->setFloatSide( Left );
    
    try
    {
      const bool showAtStartup = UserPreferences::preferenceValue<bool>( "ShowSplashScreen", viewer );
      showAgainCb->setChecked( showAtStartup );
    }catch(...)
    {
      // probably wont ever get here, but JIC
    }
    

    showAgainCb->checked().connect( boost::bind( showAgainCallback, true ) );
    showAgainCb->unChecked().connect( boost::bind( showAgainCallback, false ) );
    if( m_viewer->isMobile() )
    {
      showAgainCb->setFloatSide(Right);
      showAgainCb->setMargin(10,Right);
      showAgainCb->setMargin(3,Top);
    }//isPhone
  }//if( !force )
  
  WPushButton *ok = addCloseButtonToFooter( WString::tr("Close") );
  ok->clicked().connect( boost::bind( &AuxWindow::hide, this ) );
  
  if( viewer
     //&& (viewer->renderedWidth() > 715.0) //At initial app load we actually dont yet know client browser size; should maybe do check clientside.
      //&& (viewer->renderedHeight() > 512.0)
      && !viewer->isMobile() )
  {
    if( viewer->renderedWidth() > 100 && viewer->renderedHeight() > 100 )
    {
      //We do know the window size, so lets use it.  This is very rough - (not tested as of 20190706)
      if( viewer->renderedWidth() > 725 )
        resizeWindow( std::max(700,viewer->renderedWidth()/2), std::max(0.75*viewer->renderedHeight(),512.0) );
      else
        resizeWindow( viewer->renderedWidth() - 30, viewer->renderedHeight() - 40 );
    }else
    {
      //We dont know the window size
      resizeScaledWindow( 0.5, 0.8 );
    }
    centerWindowHeavyHanded();
  }else
  {
    //On android, sometimes on startup the screen is a little offset and not
    //  totally sized right, so we'll be a bit heavy handed here.
    const string resizejs = "var res = function(){ " + resizeScaledWindowJs(1.0,1.0) + " };";
    const string reposjs = "var repo = function(){ " + repositionWindowJs(0,0) + " };";
    
    const string js = resizejs + reposjs + "var fcn = function(){ try{ res(); repo(); }catch(e){} };"
    "fcn(); setTimeout(fcn,500); setTimeout(fcn,1500); setTimeout(fcn,5000);";
    
    doJavaScript( js );
  }//if( viewer )
  
}//UseInfoWindow::UseInfoWindow( std::function<void(bool)> showAgainCallback,InterSpec* viewer ):



UseInfoWindow::~UseInfoWindow()
{
}//~UseInfoWindow()



/**
 Initializes m_videoInfos
 **/
void UseInfoWindow::initVideoInfoMap()
{
  std::string videoInfo="";
  m_resourceBundle.resolveKey("videoinfo",videoInfo);
  
  try
  {
    Json::Object result;
    Json::parse( videoInfo, result );
    
    for( Wt::Json::Object::const_iterator iter = result.begin();
        iter != result.end(); ++iter )
    {
      const Json::Array &videos = iter->second;
      for( size_t i = 0; i < videos.size(); ++i )
      {
        const Json::Object &info = videos[i];
        const Json::Value &key = info.get("key");
        const Json::Value &file1 = info.get("fileMP4");
        const Json::Value &file2 = info.get("fileOGV");
        const Json::Value &title = info.get("title");
        
        if( key.isNull() || file1.isNull() )
        {
          cerr << "makeItem found a null key or file string, "
          "you may want to check use_instructions.xml" << endl;
          continue;
        }//if( binding.isNull() || file.isNull() )
        
        const WString bindingstr = key;
        const WString filestrMP4 = file1;
        
        VideoInformation vidinfo;
        vidinfo.key = bindingstr.toUTF8();
        vidinfo.fileMP4 = filestrMP4.toUTF8();
        
#if( !defined(IOS) )
        if (!file2.isNull())
        {
          const WString filestrOGV = file2;
          vidinfo.fileOGV = filestrOGV.toUTF8();
        }//if (!file2.isNull())
#endif
        
        if( !title.isNull() )
        vidinfo.title = title;
        m_videoInfos[iter->first].push_back( vidinfo );
      }//for( size_t i = 0; i < videos.size(); ++i )
    }//for( loop over sesults )
  }catch( std::exception &e )
  {
    cerr << "getVideoInfoMap() caught: " << e.what() << endl
    << "you may want to check use_instructions.xml" << endl;
  }//try / catch
  
  return;
}//void UseInfoWindow::initVideoInfoMap()



void UseInfoWindow::itemCreator( const string &resource, Wt::WContainerWidget *parent, WString title)
{
  std::string resourceContent="";
  m_resourceBundle.resolveKey(resource,resourceContent);
  
  WTemplate* templ = new WTemplate( parent );
  templ->setTemplateText(WString(resourceContent, UTF8),XHTMLUnsafeText);
  const std::vector<VideoInformation> videoinfo = m_videoInfos[resource];
  for( size_t i = 0; i < videoinfo.size(); ++i )
  {
    Wt::WMediaPlayer *player = new Wt::WMediaPlayer( WMediaPlayer::Video );

    //resize video size, so fits on phone
    if (m_viewer->isPhone())
      player->setVideoSize(470, 110);
    
    //Note: Try to play OGV first...ordering seems to matter!
    if (videoinfo[i].fileOGV.length()>0)     //if have OGV, add that too
    {
      player->addSource( WMediaPlayer::OGV, WLink(videoinfo[i].fileOGV) );
    }// if (videoinfo[i].fileOGV.length()>0)
    
    // Last resort, go to MP4 (does not play in Firefox)
    // http://jplayer.org/latest/developer-guide/#jPlayer-media-encoding
    if (videoinfo[i].fileMP4.length()>0)
      player->addSource( WMediaPlayer::M4V, WLink(videoinfo[i].fileMP4) );

//    player->setTitle( videoinfo[i].title );
    templ->bindWidget( videoinfo[i].key, player );
    m_players[title.toUTF8()]=player; //save this list of players so we can stop them
    
  }//for( size_t i = 0; i < videoinfo.size(); ++i )
}//void UseInfoWindow::itemCreator( const string &resource, Wt::WContainerWidget *parent, WString title)



void UseInfoWindow::right_select_item(  WMenuItem *item )
{
  m_menu->select( item );
  item->triggered().emit( item ); //doenst look like this is emmitted either
                                  //when body of SideMenuItem is clicked
                                  //stop all players
  tab_select_item(item);
}//void UseInfoWindow::select_item(  SideMenuItem *item )


void UseInfoWindow::tab_select_item( WMenuItem *item )
{
  if( !item )
    return;
  
  std::map <string, Wt::WMediaPlayer*>::iterator iter;
  
  for (iter = m_players.begin(); iter != m_players.end(); ++iter)
  {
    if( item == NULL )
      (iter->second)->stop();
    
    string tab = item->text().toUTF8(); //which tab we clicked
    
    /*
    if (item->id().compare("TutorialVideos")==0)
    {
      //if we clicked on Tutorial videos side button, then we should check which
      //video tab was previously selected
      tab = m_videoTab->tabText(m_videoTab->currentIndex()).toUTF8();
    } //if (item->text().toUTF8().compare("Tutorial videos")==0)
    */
    
    //Either match clicked tab, or just check which is selected
    if (iter->first.compare( tab)==0)
    {
      (iter->second)->play();
    } //if (iter->first.compare(item->text().toUTF8())==0)
    else
    {
      (iter->second)->stop();
    }//else
  }// for (iter = m_players.begin(); iter != m_players.end(); ++iter)
}//void UseInfoWindow::tab_select_item(WMenuItem *item)


SideMenuItem * UseInfoWindow::makeItem( const WString &title, const string &resource)
{
  //Using a lambda to create the callback results in videos playing multiple times...
  //std::function<void(WContainerWidget *)> f = [this,title,resource](WContainerWidget *w){
    //itemCreator( resource, w, title );
  //};
  
  std::function<void(WContainerWidget *)> f = boost::bind( &UseInfoWindow::itemCreator, this,
                                                          resource, boost::placeholders::_1,
                                                          title );
  
  WWidget *w = deferCreate( f );
  w->addStyleClass( "UseInfoItem" );
  
  SideMenuItem *item = new SideMenuItem( title, w );
  item->clicked().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
  item->mouseWentDown().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
  
  m_menu->addItem( item );
  
  return item;
}//SideMenuItem * UseInfoWindow::makeItem( const WString &title, const string &resource)


SideMenuItem *UseInfoWindow::makeTextItem( const Wt::WString &title, const std::string &resource )
{
  std::function<void(WContainerWidget *)> f = boost::bind( &UseInfoWindow::textItemCreator, this,
                                                          resource, boost::placeholders::_1 );
  
  WWidget *w = deferCreate( f );
  w->addStyleClass( "UseInfoItem" );
  
  SideMenuItem *item = new SideMenuItem( title, w );
  item->clicked().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
  item->mouseWentDown().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
  
  m_menu->addItem( item );
  
  return item;
}//SideMenuItem * makeTextItem( const Wt::WString &title, const std::string &resource );


void UseInfoWindow::textItemCreator( const std::string &resource, Wt::WContainerWidget *parent )
{
  assert( parent );
  if( !parent )
    return;
  
  try
  {
    const string docroot = wApp->docRoot();
    const string resource_base = SpecUtils::append_path(docroot, resource);
    const string resource_name = SpecUtils::filename(resource);
    
    WMessageResourceBundle bundle;
    bundle.use( resource_base, true );
    
    vector<pair<string,string>> topic_keys;
    
    // We want to load the topics in the order they are in the file, not alphabetical like
    //  `bundle.keys(WMessageResourceBundle::Scope::Default)` would return.
    {// Begin scope to manually parse out XML file
      string resource_file = AppUtils::find_localized_xml_file( resource_base );
      
      vector<char> data;
      SpecUtils::load_file_data( resource_file.c_str(), data );
      
      rapidxml::xml_document<char> doc;
      doc.parse<rapidxml::parse_trim_whitespace>( data.data() );
      
      const rapidxml::xml_node<char> * const msgs_node = doc.first_node("messages");
      if( !msgs_node )
        throw runtime_error( "No 'messages' node." );
      
      for( rapidxml::xml_node<char> *node = msgs_node->first_node("message");
          node; node = node->next_sibling("message") )
      {
        const rapidxml::xml_attribute<char> * const id = node->first_attribute("id");
        const rapidxml::xml_attribute<char> * const title = node->first_attribute("title");
        assert( id && id->value_size() );
        
        if( id && id->value_size() && title && title->value_size() )
        {
          string key{id->value(), id->value() + id->value_size()};
          string desc{title->value(), title->value() + title->value_size()};
          
#if( PERFORM_DEVELOPER_CHECKS )
          for( const auto &i : topic_keys )
          {
            if( i.first == key )
              log_developer_error( __func__, ("FAQ topic ID " + key + "is non-unique").c_str() );
            assert( i.first != key );
          }//for( const auto &i : topic_keys )
#endif
          
          if( key != "intro" )
            topic_keys.emplace_back( std::move(key), std::move(desc) );
        }//if( we have key and title strings )
      }//for( loop over nodes )
    }// End scope to manually parse out XML file
    
    if( topic_keys.empty() )
      throw runtime_error( "Error getting 'message' nodes." );
    
    const bool has_css = SpecUtils::is_file(resource_base + ".css");
    assert( has_css );
    if( has_css )
      wApp->useStyleSheet( resource + ".css" ); //resource_base may be absolute filesystem path; resource is
    
    // Add some JS to scroll to the contents
    wApp->declareJavaScriptFunction( resource_name + "_scroll_to",
    "function(id){"
      //"const sec_el = document.getElementById(id);"
      "const sec_el = document.querySelector('[id^=' + id + ']');"  //match ID starting with specified ID
      "if(!sec_el){console.error('No element with id=', id); return;}"
      "const parent_el = document.getElementById('" + parent->id() + "');"
      "if( !parent_el ){console.error('No element with id=" + parent->id() + "'); return;}"
      "parent_el.scrollTop = sec_el.offsetTop;"
    "}");
    
    // Create intro
    WContainerWidget *into_div = new WContainerWidget( parent );
    into_div->addStyleClass( resource_name + "-intro" );
    
    std::string intro_txt;
    bundle.resolveKey( "intro", intro_txt );
    WText *intro = new WText( WString::fromUTF8(intro_txt), into_div );
    intro->setInline( false );
    intro->addStyleClass( resource_name + "-intro-txt" );
    
    // Create a table of contents
    WContainerWidget *into_toc = new WContainerWidget( into_div );
    into_toc->setList( true );
    into_div->addStyleClass( resource_name + "-toc" );
        
    // Add the contents
    WContainerWidget *subjects = new WContainerWidget( parent );
    subjects->addStyleClass( resource_name + "-subjects" );
    
    for( const pair<string,string> &key_title : topic_keys )
    {
      WContainerWidget *w = new WContainerWidget( subjects );
      w->addStyleClass( resource_name + "-subject" );
      
      WContainerWidget *title = new WContainerWidget( w );
      title->addStyleClass( resource_name + "-subject-title" );
      WPushButton *up_btn = new WPushButton( title );
      HelpSystem::attachToolTipOn( {up_btn}, WString::tr("uiw-back-to-top"), true );
      up_btn->setIcon( WLink("InterSpec_resources/images/minimal_go_up_arrow.svg") );
      up_btn->clicked().connect( "function(){ document.getElementById('" + parent->id() + "').scrollTop = '0px'; }" );
      WText *title_txt = new WText( key_title.second, title );
        
      std::string content;
      bundle.resolveKey( key_title.first, content );
      
      //cout << "key: " << key_title.first << " --> content: " << content << endl;
      
      WText *sub = new WText( content, w );
      sub->addStyleClass( resource_name + "-subject-content" );
      
      w->setObjectName( key_title.first );
      
      // Add to the TOC
      WContainerWidget *toc_item = new WContainerWidget( into_toc );
      const string this_js_call = "function(){"
        + wApp->javaScriptClass() + "." + resource_name + "_scroll_to('" + w->id() + "');"
      "}";
      toc_item->clicked().connect( this_js_call );
      WText *link = new WText( key_title.second, toc_item );
      link->addStyleClass( resource_name + "-toc-item LinkBtn" );
    }//for( const pair<string,string> &key_title : topic_keys )
  }catch( std::exception &e )
  {
    new WText( WString::tr("uiw-error-loading-faq").arg(e.what()), parent );
  }//try catch
}//void textItemCreator( const std::string &resource, Wt::WContainerWidget *parent )


void UseInfoWindow::handleSampleSelectionChanged()
{
}//void handleSampleSelectionChanged()


void UseInfoWindow::handleSampleDoubleClicked( WModelIndex index, WMouseEvent event )
{
  if( event.button() != WMouseEvent::LeftButton || !index.isValid() )
    return;
  
  //[for macOS build at least] There seems to be a wierd issue that if you
  //  double click on an item with out it being highlighted, then the spectrum
  //  loads, but doesnt show until you try to interact with it... weird!
  
  WModelIndexSet selected;
  selected.insert( index );
  m_tableSample->setSelectedIndexes( selected );
  loadSampleSelected();
}

void UseInfoWindow::handleAppUrl( std::string query_str )
{
  const map<string,string> parts = AppUtils::query_str_key_values( query_str );
  
  const auto topic_pos = parts.find("TOPIC");
  
  if( topic_pos != end(parts) )
  {
    for( WMenuItem *item : m_menu->items() )
    {
      //item->pathComponent() will be "faqs", and is probably a better thing to use
      if( SpecUtils::iequals_ascii(item->text().toUTF8(), topic_pos->second) )
      {
        right_select_item( item );
        
        if( SpecUtils::iequals_ascii("faqs", topic_pos->second) )
        {
          const auto subtopic_pos = parts.find("SUBTOPIC");
          if( subtopic_pos != end(parts) )
          {
            string topic_id = WString::fromUTF8(subtopic_pos->second).jsStringLiteral();
            doJavaScript( "setTimeout( function(){try{"
                         + wApp->javaScriptClass() + "." + "faqs_scroll_to(" + topic_id
                         + ");}catch(e){}}, 500 );" );
          }
        }//if( a FAQs topic )
        
        // TODO: add selecting tab on "Controls" and "Welcome" tabs.

        break;
      }
    }//for( loop over menu items )
  }//if( topic_pos != end(parts) )
}//void handleAppUrl( std::string query_str )


void UseInfoWindow::showFaqTab()
{
  for( WMenuItem *item : m_menu->items() )
  {
    //item->pathComponent() will be "faqs", and is probably a better thing to use
    if( item->text() != "FAQs" )
      continue;
    right_select_item( item );
    return;
  }//for( loop over menu items )
}//void showFaq( std::string &topic_id = "" );


void UseInfoWindow::loadSample( const Wt::WModelIndex index )
{
  if( !index.isValid() )
    return;
  
  WString val = boost::any_cast<WString> (m_messageModelSample->data(index.row(), 1, DisplayRole));
  SpecUtils::ParserType parserType = boost::any_cast<SpecUtils::ParserType> (m_messageModelSample->data(index.row(), 2, UserRole));
  
  SpecMeasManager* manager = m_viewer->fileManager();
  
  const string filepath = val.toUTF8();
  
  cout << "Will try to load '" << filepath << "'" << endl;
  
  manager->loadFromFileSystem( filepath, SpecUtils::SpectrumType::Foreground, parserType );
}//void UseInfoWindow::loadSample( const Wt::WModelIndex index )


void UseInfoWindow::loadSampleSelected()
{
  WModelIndexSet indices = m_tableSample->selectedIndexes();
  if( indices.empty() )
  {
    cout << "UseInfoWindow::loadSampleSelected(): No selected indices." << endl;
    return;
  }
  
  loadSample( *indices.begin() );
}//void loadSampleSelected( )


SideMenuItem::SideMenuItem( const WString &text, WWidget *contents )
: WMenuItem( text, contents, WMenuItem::LazyLoading )
{
} //SideMenuItem::SideMenuItem( const WString &text, WWidget *contents ): WMenuItem( text, contents, WMenuItem::LazyLoading )



SignalBase &SideMenuItem::activateSignal()
{
  return triggered();
} //SignalBase &SideMenuItem::activateSignal()
