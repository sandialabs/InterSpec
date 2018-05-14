/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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

#include <Wt/WText>
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
#include <Wt/WTabWidget>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/Json/Parser>
#include <Wt/Json/Object>
#include <Wt/WMediaPlayer>
#include <Wt/WRadioButton>
#include <Wt/WStandardItem>
#include <Wt/WStackedWidget>
#include <Wt/WContainerWidget>
#include <Wt/WStandardItemModel>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/UseInfoWindow.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"


using namespace Wt;
using namespace std;


UseInfoWindow::UseInfoWindow( std::function<void(bool)> showAgainCallback,
                              InterSpec* viewer ):
  AuxWindow("Welcome",true),
  m_session(),
  m_tableSample(NULL),
  m_loadSampleButton(NULL),
  m_messageModelSample(NULL),
  m_viewer(viewer),
  m_menu(NULL)
{
  rejectWhenEscapePressed();
  
  double width = 0.0, height = 0.0;
  bool hideStateTab = false, hideSpecTab = false;
  
  if( viewer )
  {
    width = 0.5*viewer->renderedWidth();
    height = 0.8*viewer->renderedHeight();
    if( height < 512.0 )
      height = 1.0*std::min( viewer->renderedHeight(), 512 );
    height = std::min( height, 1024.0 );  //1024 not actually tested, could maybye bee 800
  }//if( viewer )

  if( width < 7150.0 || height < 512.0 )
  {
    setMinimumSize(715,512);
    resize( WLength(50, WLength::FontEm), WLength(80,WLength::Percentage));
  }else
  {
    resize( WLength(width), WLength(height) );
  }
  
  centerWindow();
  setResizable( true );
  disableCollapse();
  
  WStackedWidget *stack = new WStackedWidget();
  stack->addStyleClass( "UseInfoStack" );
  stack->setOverflow( WContainerWidget::OverflowAuto);
  
  WAnimation animation(Wt::WAnimation::Fade, Wt::WAnimation::Linear, 200);
  stack->setTransitionAnimation( animation, true );
  
  m_menu = new WMenu( stack, Wt::Vertical );
  m_menu->addStyleClass( (m_viewer->isMobile() ? "SideMenuPhone" : "SideMenu") );
  
  WDialog::contents()->setOverflow( WContainerWidget::OverflowHidden );
  
  WGridLayout *layout = stretcher();
  layout->addWidget( m_menu, 0, 0, AlignLeft );
  layout->addWidget( stack, 0, 1, 2, 1 );
  layout->setRowStretch( 0, 1 );
  
  m_resourceBundle.use("InterSpec_resources/static_text/use_instructions",false);
  
  initVideoInfoMap();
  WContainerWidget* welcomeContainer = new WContainerWidget();
  welcomeContainer->setOverflow(WContainerWidget::OverflowAuto);
  welcomeContainer->setOffsets(WLength(0,WLength::Pixel));
  welcomeContainer->setMargin(WLength(0,WLength::Pixel));
  WGridLayout* rightLayout = new WGridLayout();
  welcomeContainer->setLayout(rightLayout);
  
  Wt::WTabWidget *tabW = new Wt::WTabWidget();
  const WBorder border( WBorder::Solid, WBorder::Thin, Wt::gray);
    tabW->contentsStack()->decorationStyle().setBorder( border,  Wt::Bottom | Wt::Right | Wt::Left );
    tabW->setMargin(2);
    tabW->setOffsets(2);
  rightLayout->addWidget(tabW,1,0);
  
  WContainerWidget* spectrumContainer = new WContainerWidget();
  WGridLayout* spectrumLayout = new WGridLayout();
  spectrumLayout->setContentsMargins( 3, 5, 4, 5);
  spectrumContainer->setLayout(spectrumLayout);
  spectrumContainer->setOverflow(WContainerWidget::OverflowHidden);
  spectrumContainer->setMargin(0);
  spectrumContainer->setOffsets(0);
  
  WContainerWidget* samplesContainer = new WContainerWidget();
  WGridLayout* samplesLayout = new WGridLayout();
  samplesLayout->setContentsMargins( 3, 5, 4, 0 );
  samplesContainer->setLayout(samplesLayout);
  samplesContainer->setOverflow(WContainerWidget::OverflowAuto);


  
  //The bellow tabs must use WTabWidget::PreLoading, or else we will get a
  //  JavaScript exception when showLoadingIndicator() is called (in the JS).
  //  This is super confusing, and I dont understand it.  This manifests
  //  particularly on Android native version of app, but I think it might
  //  elsewhere as well.  (wcjohns 20150124)
  WMenuItem *spectraItem = tabW->addTab(spectrumContainer, "Saved Snapshots", WTabWidget::PreLoading);
  spectraItem->setIcon("InterSpec_resources/images/book.png");
  WMenuItem *samplesItem = tabW->addTab(samplesContainer, "Samples", WTabWidget::PreLoading);
  samplesItem->setIcon("InterSpec_resources/images/images.png");
  
  //----------Add workspace tab --------
  
  //We have to create a independant Dbo::Session for this class since the
  //  viewer->user.session() could be used in other threads, messing
  //  things up (Dbo::Session is not thread safe).
  m_session.reset( new DataBaseUtils::DbSession( *m_viewer->sql() ) );

  try
  {
#if( USE_DB_TO_STORE_SPECTRA )
    Dbo::ptr<InterSpecUser> user = m_viewer->m_user;
    SpecMeasManager* manager = m_viewer->fileManager();
    WContainerWidget *bottom = new WContainerWidget();
    m_snapshotModel = new SnapshotFactory( manager, m_viewer, "", (SpectrumType)0,
                                           std::shared_ptr<SpectraFileHeader>(),
                                           spectrumLayout, 0, bottom );
    spectrumLayout->addWidget( bottom, spectrumLayout->rowCount()+1 , 0, AlignRight );
#endif
    
    ///-----Sample -----
    m_messageModelSample = new Wt::WStandardItemModel(0,3, this);
    
    WStandardItem* parserType = new WStandardItem();
    parserType->setData(k2006Icd1Parser); //kIaeaParser
    vector<WStandardItem*> message(3);
    message[0] = new Wt::WStandardItem("Ba-133 (16k bin N42)");
    message[1] = new Wt::WStandardItem("example_spectra/ba133_source_640s_20100317.n42");
    message[2] = parserType;
    m_messageModelSample->appendRow(message);
    
    if( viewer->isMobile() )
    {
      parserType = new WStandardItem();
      parserType->setData(k2006Icd1Parser); //kIaeaParser
      message[0] = new Wt::WStandardItem("Passthrough (16k bin 1064 meas N42)");
      message[1] = new Wt::WStandardItem("example_spectra/passthrough.n42");
      message[2] = parserType;
      m_messageModelSample->appendRow(message);
    }//if( viewer->isMobile() )
    else {
      parserType = new WStandardItem();
      parserType->setData(k2006Icd1Parser); //kIaeaParser
      message[0] = new Wt::WStandardItem("Passthrough (16k bin ICD1, 8 det., 133 samples)");
      message[1] = new Wt::WStandardItem("example_spectra/passthrough.n42");
      message[2] = parserType;
      m_messageModelSample->appendRow(message);
    }//else
    
    parserType = new WStandardItem();
    parserType->setData(k2006Icd1Parser); //kIaeaParser
    message[0] = new Wt::WStandardItem("Background (16k bin N42)");
    message[1] = new Wt::WStandardItem("example_spectra/background_20100317.n42");
    message[2] = parserType;
    m_messageModelSample->appendRow(message);
    
    if( viewer->isMobile() )
    {
      parserType = new WStandardItem();
      parserType->setData(kIaeaParser); //k2006Icd1Parser
      message[0] = new Wt::WStandardItem("Ba-133 (Low Res, No Calib)");
      message[1] = new Wt::WStandardItem("example_spectra/Ba133LowResNoCalib.spe");
      message[2] = parserType;
      m_messageModelSample->appendRow(message);
      
      parserType = new WStandardItem();
      parserType->setData(kIaeaParser); //k2006Icd1Parser
      message[0] = new Wt::WStandardItem("Co-60 (Low Res, No Calib)");
      message[1] = new Wt::WStandardItem("example_spectra/Co60LowResNoCalib.spe");
      message[2] = parserType;
      m_messageModelSample->appendRow(message);
      
      parserType = new WStandardItem();
      parserType->setData(kIaeaParser); //k2006Icd1Parser
      message[0] = new Wt::WStandardItem("Cs-137 (Low Res, No Calib)");
      message[1] = new Wt::WStandardItem("example_spectra/Cs137LowResNoCalib.spe");
      message[2] = parserType;
      m_messageModelSample->appendRow(message);
      
      parserType = new WStandardItem();
      parserType->setData(kIaeaParser); //k2006Icd1Parser
      message[0] = new Wt::WStandardItem("Th-232 (Low Res, No Calib)");
      message[1] = new Wt::WStandardItem("example_spectra/Th232LowResNoCalib.spe");
      message[2] = parserType;
      m_messageModelSample->appendRow(message);
    } // if( viewer->isMobile() )
    
    m_messageModelSample->setHeaderData(  0, Horizontal, WString("Spectra"), DisplayRole );
    m_tableSample = new RowStretchTreeView();
    m_tableSample->setRootIsDecorated	(	false); //makes the tree look like a table! :)
    m_tableSample->addStyleClass( "DbSpecFileSelectTable" );
    m_tableSample->setModel(m_messageModelSample);
    m_tableSample->setOffsets( WLength(0, WLength::Pixel), Wt::Left | Wt::Top );
    m_tableSample->setColumnResizeEnabled(true);
    m_tableSample->setColumnAlignment(0, Wt::AlignLeft);
    m_tableSample->setHeaderAlignment(0, Wt::AlignCenter);
    m_tableSample->setColumnWidth(0,340);
    m_tableSample->hideColumn(1);
    m_tableSample->hideColumn(2);
    
    m_tableSample->setAlternatingRowColors(true);
    m_tableSample->setSelectionMode(Wt::SingleSelection);
    m_tableSample->setEditTriggers(Wt::WAbstractItemView::NoEditTrigger);
    
    samplesLayout->addWidget(m_tableSample,0,0);
    samplesLayout->setRowStretch(0,1);
    samplesLayout->setColumnStretch(0,1);
    
    WContainerWidget* buttons = new WContainerWidget();
    m_loadSampleButton = new WPushButton( "Load", buttons );
    m_loadSampleButton->setStyleClass("DatabaseGoIcon");
    m_loadSampleButton->setFloatSide(Right);
    m_loadSampleButton->clicked().connect( this, &UseInfoWindow::loadSampleSelected );
    
    m_loadSampleButton->disable();
    m_tableSample->selectionChanged().connect( this, &UseInfoWindow::handleSampleSelectionChanged );
    m_tableSample->doubleClicked().connect( boost::bind(&UseInfoWindow::loadSample, this, _1 ) );
    
    samplesLayout->addWidget( buttons, samplesLayout->rowCount()+1 , 0, AlignRight );
      
      
    const int sampleIndex = tabW->indexOf(samplesContainer);
    
#if( !USE_DB_TO_STORE_SPECTRA )
    tabW->setCurrentIndex( sampleIndex );
#else
    if( m_snapshotModel->size() == 0 )
      tabW->setCurrentIndex( sampleIndex );
#endif
  } //try
  catch( exception &e )
  {
    cerr << "Exception:" << e.what() << endl;
  }//catch
  
  
  //----------Add left menu --------
  
  SideMenuItem *item = new SideMenuItem( m_viewer->isMobile()?"":"Welcome", welcomeContainer );
  item->clicked().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
  item->mouseWentDown().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
  item->setIcon("InterSpec_resources/images/user.png");
  m_menu->addItem( item );
  
  
  // ---- Videos
//QT5 does not play back videos correctly
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
    
/*
    item = new SideMenuItem( m_viewer->isMobile()?"":"Tutorial videos", controlContainer );
    item->setId("TutorialVideos");
    item->clicked().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
    item->mouseWentDown().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
    item->setIcon("InterSpec_resources/images/film.png");
    m_menu->addItem( item );
*/
  } //videos
    
  //----------Create Controls window--------
  
  //Select either mobile or desktop mouse/keyboard/gesture page
  if( m_viewer->isMobile() )
  {
    item = makeItem( "", "mobile-mouse-interactions" );
    item->setIcon("InterSpec_resources/images/iphone.png");
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
          
        if( title.isNull() || image.isNull() )
        {
          cerr << "makeItem found a null key or file string, "
                  "you may want to check use_instructions.xml" << endl;
          continue;
        }//if( binding.isNull() || file.isNull() )
          
        const WString titlestr = title;
        const WString imagestr = image;
        string descstr = "";
          
        if( !desc_id.isNull() )
          m_resourceBundle.resolveKey( desc_id, descstr );
        
        WContainerWidget* samplesContainer = new WContainerWidget();
        samplesContainer->setOffsets(15,Wt::Top);
        samplesContainer->setMargin(15,Wt::Top);
        samplesContainer->setOverflow(WContainerWidget::OverflowAuto);
        
        WImage* img = new WImage( imagestr.toUTF8() );
        img->setMaximumSize(WLength(300), WLength::Auto);
        samplesContainer->addWidget(img);
        img->setStyleClass("img-centered");
        if (!descstr.empty())
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
    
    
    item = new SideMenuItem( "Controls", controlContainer );
    item->clicked().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
    item->mouseWentDown().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
    item->setIcon("InterSpec_resources/images/mouse.png");
    m_menu->addItem( item );
  } //isDesktop()

  
  // ---- Cheat sheet PDF
  if( m_viewer->isSupportFile() )
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
    anchor->setImage(new WImage("InterSpec_resources/images/page_white_acrobat.png"));
    anchor->setTarget( Wt::TargetNewWindow );
    anchor->setText("Cheat Sheet (PDF)");

    controlLayout->addWidget(new WText("Export and print out this cheat sheet (Adobe Acrobat PDF format) to quickly reference actions, shortcut keys and troubleshooting."),0,0);
    controlLayout->addWidget(anchor,1,0);
    controlLayout->setRowStretch(1, 1);
    
    item = new SideMenuItem( m_viewer->isMobile()?"":"Cheat Sheet", controlContainer );
    item->clicked().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
    item->mouseWentDown().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
    item->setIcon("InterSpec_resources/images/page_white_acrobat.png");
    m_menu->addItem( item );
  } //Cheat sheet PDF
  
  char msg[512];
  size_t pos = 0;
  const char firstline[] = "<h1 align=\"center\">Welcome To <em>InterSpec</em>"
                           "</h1><div align=\"left\">";
  
  memcpy( msg+pos, firstline, sizeof(firstline)-1 );
  pos += sizeof(firstline)-1;
  
  if( m_viewer->isMobile() )
  {
    const char line[] = "Tutorials are on the left, while sample and previously"
                        " opened spectra and snapshots are bellow.";
    memcpy( msg+pos, line, sizeof(line)-1 );
    pos += sizeof(line)-1;
  }else if( hideStateTab && hideSpecTab )
  {
    const char line[] = "You can learn the basics of using <em>InterSpec</em> "
                        "by clicking through the tutorials, controls, and cheat"
                        " sheet on the left.";
    memcpy( msg+pos, line, sizeof(line)-1 );
    pos += sizeof(line)-1;
  }else
  {
    const char line[] = "You can use the references on the left to become "
                        "more familiar with <em>InterSpec</em>, or you can "
                        "pick up from a previous session or spectrum bellow.";
    memcpy( msg+pos, line, sizeof(line)-1 );
    pos += sizeof(line)-1;
  }

  if( !m_viewer->isMobile() )
  {
    const char lastline[] = " Alternatively, you can also drag and drop your "
                            "own spectrum file onto <em>InterSpec</em>.";
    memcpy( msg+pos, lastline, sizeof(lastline) );
    pos += sizeof(lastline)-1;
  }
  
  const char divend[] = "</div>";
  memcpy( msg+pos, divend, sizeof(divend) );  //copies '\0' as well
  pos += sizeof(divend)-1;
  
  assert( pos < 512 );  //strlen(msg) is either 294 or 318 currently
  
  WText *text = new WText( msg, XHTMLUnsafeText );
  rightLayout->addWidget( text, 0, 0 );
  rightLayout->setRowStretch( 1, 1 );
  rightLayout->setContentsMargins( 9, 1, 9, 1 );


  m_menu->select( 0 );
  
  if (!m_viewer->isMobile())
      layout->setRowStretch( 0, 1 );
  layout->setColumnStretch( 1, 1 );
  
  WContainerWidget *bottom = footer();

//  footer()->resize( WLength::Auto, WLength(50) );
  
    
   
    if (m_viewer->isMobile())
    {
        //phone position
        AuxWindow::addHelpInFooter(bottom, "");
    }
    else
    {
        //regular locations
        WAnchor *more = new WAnchor( WLink(), "More in depth information" );
        
        more->clicked().connect( boost::bind( &AuxWindow::hide, this ) );
        more->clicked().connect( boost::bind( &HelpSystem::createHelpWindow, string("setting-up") ) );
        more->addStyleClass("InfoLink");
        layout->addWidget( more, 1, 0 );
    }
    layout->setContentsMargins( 0, 9, 9, 0 );
    
    
  if( showAgainCallback )
  {
    WCheckBox *cb = new WCheckBox( "Always show at startup", bottom );
    cb->setFloatSide( Left );
    cb->setChecked(viewer->m_user->preferenceValue<bool>( "ShowSplashScreen" ));

    cb->checked().connect( boost::bind( showAgainCallback, true ) );
    cb->unChecked().connect( boost::bind( showAgainCallback, false ) );
    if (m_viewer->isMobile())
    {
        cb->setFloatSide(Right);
        cb->setMargin(10,Right);
//        cb->setMargin(WLength(20,WLength::Pixel),Wt::Left);
    } //isPhone
  }//if( !force )

  
        
//  {
//    WContainerWidget *dummy = new WContainerWidget();
//    SideMenuItem *item = new SideMenuItem( "More In Depth Info", dummy );
// //    item->setIcon("InterSpec_resources/images/page_white_acrobat.png");
//    item->clicked().connect( boost::bind( &AuxWindow::hide, this ) );
//    item->clicked().connect( boost::bind( &HelpSystem::createHelpWindow, string("setting-up") ) );
//    m_menu->addItem( item );
//  }
  
  WPushButton *ok = addCloseButtonToFooter();
  ok->clicked().connect( boost::bind( &AuxWindow::hide, this ) );
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
    
    if (item->id().compare("TutorialVideos")==0)
    {
      //if we clicked on Tutorial videos side button, then we should check which
      //video tab was previously selected
      tab = m_videoTab->tabText(m_videoTab->currentIndex()).toUTF8();
    } //if (item->text().toUTF8().compare("Tutorial videos")==0)
    
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
  
  std::function<void(WContainerWidget *)> f = boost::bind( &UseInfoWindow::itemCreator, this, resource, _1, title );
  
  WWidget *w = deferCreate( f );
  w->addStyleClass( "UseInfoItem" );
  
  SideMenuItem *item = new SideMenuItem( title, w );
  item->clicked().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
  item->mouseWentDown().connect( boost::bind( &UseInfoWindow::right_select_item, this, item) );
  
  m_menu->addItem( item );
  
  return item;
}//SideMenuItem * UseInfoWindow::makeItem( const WString &title, const string &resource)


void UseInfoWindow::handleSampleSelectionChanged()
{
  if( m_tableSample->selectedIndexes().size() )
    m_loadSampleButton->enable();
  else
    m_loadSampleButton->disable();
}//void handleSampleSelectionChanged()

void UseInfoWindow::loadSample( const Wt::WModelIndex index )
{
  if( !index.isValid() )
    return;
  
  WString val = boost::any_cast<WString> (m_messageModelSample->data(index.row(), 1, DisplayRole));
  ParserType parserType = boost::any_cast<ParserType> (m_messageModelSample->data(index.row(), 2, UserRole));
  
  SpecMeasManager* manager = m_viewer->fileManager();
  
  manager->loadFromFileSystem(val.toUTF8(), kForeground, parserType);
}//void UseInfoWindow::loadSample( const Wt::WModelIndex index )

void UseInfoWindow::loadSampleSelected( )
{
  WModelIndexSet indices = m_tableSample->selectedIndexes();
  if( !indices.size() )
  {
    m_loadSampleButton->disable();
    return;
  }//if( !indices.size() )
  
  loadSample( *indices.begin() );
//  if( manager->loadFromFileSystem(val.toUTF8(), kForeground, parserType) )
//    AuxWindow::deleteAuxWindow(this); //closes help window
} //void UseInfoWindow::loadSampleSelected( )

SideMenuItem::SideMenuItem( const WString &text, WWidget *contents )
: WMenuItem( text, contents, WMenuItem::LazyLoading )
{
} //SideMenuItem::SideMenuItem( const WString &text, WWidget *contents ): WMenuItem( text, contents, WMenuItem::LazyLoading )



SignalBase &SideMenuItem::activateSignal()
{
  return triggered();
} //SignalBase &SideMenuItem::activateSignal()
