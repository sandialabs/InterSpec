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
#include <fstream>

#include <Wt/WText>
#include <Wt/WTree>
#include <Wt/WTimer>
#include <Wt/WImage>
#include <Wt/WString>
#include <Wt/WTemplate>
#include <Wt/WTreeNode>
#include <Wt/WLineEdit>
#include <Wt/Json/Value>
#include <Wt/Json/Array>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/Json/Object>
#include <Wt/Json/Parser>
#include <Wt/WContainerWidget>
#include <Wt/WMessageResourceBundle>

#include "SpecUtils/DateTime.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/UndoRedoManager.h"

#if( USE_OSX_NATIVE_MENU )
#include "InterSpec/PopupDiv.h"
#include "target/osx/NativeMenu.h"
#endif


using namespace std;
using namespace Wt;


namespace
{
  /** A class used to track when the widget that a qtip (tool tip created by
      attachToolTipOn(..) ) is attached to, gets deleted.  When the original
      widget gets deleted, we should delete the qtip from the client-side DOM.
      Therefore, everytime we attach a qtip, we will create an instance of
      RmHelpFromDom as the child of the widget getting the tooltip attached,
      then from the destructor of RmHelpFromDom, we will issue the JS to cleanup
      the DOM.  This all requires the 'id' of the qtip is set to the 'id' of
      this class (the qtip library will append a 'qtip-' to the id, so there
      will be no DOM conflicts).
   */
  class RmHelpFromDom : public Wt::WObject
  {
  public:
    RmHelpFromDom( WObject *parent )
    : WObject( parent )
    {
    }
    
    virtual ~RmHelpFromDom()
    {
      WApplication *app = wApp;
      if( !app || !app->domRoot() )
        return;
      
      const string js = "try{ $('#qtip-" + id() + "').qtip('destroy', true);}catch(error){}";  //I dont think try/catch is necessary, but JIC
      app->doJavaScript( js );
    }
  };//RmHelpFromDom
  
  /** Returns the contents of the localized, specified file, as a string.
   
   It doesn't directly open the file you specify, but uses the same logic as wt-3.7.1/src/Wt/WMessageResources.C
   to find the version of the specified file closest to the current locale.
   
   Throws exception, with descriptive message, if couldn't open file.
   */
  string getInternationalizedFileContents( const string &orig_filepath )
  {
    string filepath = orig_filepath;
    string extension = SpecUtils::file_extension(filepath);
    if( !extension.empty() )
      filepath = filepath.substr(0,filepath.size() - extension.size());
    else
      extension = ".xml";
    
    string locale = WLocale::currentLocale().name();
    
    for( ; ; )
    {
      const string trial_path = filepath + (locale.length() > 0 ? "_" : "") + locale + extension;
      
#ifdef _WIN32
      const std::wstring wfilepath = SpecUtils::convert_from_utf8_to_utf16(trial_path);
      ifstream infile( wfilepath.c_str() );
#else
      ifstream infile( trial_path.c_str() );
#endif
      
      if( !infile.is_open() )
      {
        if( locale.empty() )
          break;
        
        /* try a lesser specified variant */
        std::string::size_type l = locale.rfind('-');
        if( l != std::string::npos )
          locale.erase(l);
        else
          locale = "";
        
        continue;
      }//if( !infile.is_open() )

      infile.seekg(0, std::ios::end);
      const size_t size = infile.tellg();
      infile.seekg(0);
      
      if( size < 128*1024 )
      {
        string contents( size, ' ' );
        infile.read( &(contents[0]), size );
        return contents;
      }
      
      throw runtime_error( "'" + filepath + "' is too large a file - not loading contents." );
    }//for( ; ; )
    
    throw runtime_error( "Could not open expected file: '" + orig_filepath + "'" );
    
    return "";
  };//getInternationalizedFileContents(...)
}//namespace


namespace HelpSystem
{
  void createHelpWindow( const std::string &preselect )
  {
    InterSpec::instance()->showHelpWindow( preselect );
  }//void createHelpWindow( const std::string &preselect )
  
  HelpWindow::HelpWindow(std::string preselect)
   : AuxWindow( WString::tr("window-title-help"),
               (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
                 | AuxWindowProperties::IsHelpWindow
                 | AuxWindowProperties::EnableResize
                 | AuxWindowProperties::SetCloseable) ),
     m_tree ( new Wt::WTree() ),
     m_searchText( new Wt::WLineEdit() )
  {
    UndoRedoManager::BlockUndoRedoInserts undo_blocker;
    
    wApp->useStyleSheet( "InterSpec_resources/HelpWindow.css" );
    InterSpec::instance()->useMessageResourceBundle( "HelpSystem" );
    
    m_tree->setSelectionMode(Wt::SingleSelection);
    
    //To Do add in special
    //"help-xml"
    //"help-xml-desktop"
    //"help-xml-mobile"
    //"help-xml-phone"
    //"help-xml-tablet"
    
    const string docroot = wApp->docRoot();
    const std::string help_json = SpecUtils::append_path(docroot,"InterSpec_resources/static_text/help.json");
    
    try
    {
      m_helpLookupTable = getInternationalizedFileContents(help_json);
    }catch( std::exception & )
    {
      passMessage( WString::tr("hw-err-opening-json").arg(help_json), 2 );
    }
    
    initialize();
    
    m_helpWindowContent = new WContainerWidget();
    m_helpWindowContent->setPadding(WLength(10,WLength::Pixel));
    
    m_tree->setMargin(WLength(0,WLength::Pixel));
    m_tree->addStyleClass("helpTree"); //both need to scroll
    m_helpWindowContent->addStyleClass("helpTree"); //both need to scroll
    m_tree->itemSelectionChanged().connect(boost::bind( &HelpWindow::selectHelpToShow, this ) );
    
    setTopic( preselect );
    
    //contents()->setMargin(0);
    //contents()->setPadding(0);
    
    contents()->keyPressed().connect( this, &HelpWindow::handleArrowPress );
    WGridLayout *layout = stretcher();
    m_searchText->setStyleClass("searchBox");
    m_searchText->setPlaceholderText( WString::tr("hw-search-empty-text") );
    m_searchText->setMinimumSize(Wt::WLength::Auto, WLength(1.5,Wt::WLength::FontEm));
    m_searchText->setMargin(WLength(5,WLength::Pixel));
//    m_searchText->keyPressed().connect( this, &HelpWindow::initialize );
    m_searchText->keyWentDown().connect( this, &HelpWindow::handleKeyPressInSearch );
//    m_searchText->keyWentUp().connect( this, &HelpWindow::initialize );
    
//    m_searchText->keyWentUp	().connect(std::bind([=] (const Wt::WKeyEvent& e) {
// //      if (e.key()==Wt::Key_Space || e.key()==Wt::Key_Enter || e.key()==Wt::Key_Delete|| e.key()==Wt::Key_Backspace)
//        initialize();
//    }, std::placeholders::_1));

//    m_searchText->keyWentUp().connect(std::bind([=]() {
//      if (m_searchText->text().empty())
//        initialize();
//    }));
    layout->addWidget( m_searchText, 0, 0 );
    layout->addWidget( m_tree, 1, 0 );
    layout->addWidget( m_helpWindowContent, 0, 1, 2, 1 );
    layout->setContentsMargins(0,0,0,0);
    layout->setRowStretch( 1, 1 );
    layout->setColumnResizable( 0, true, WLength(25,WLength::FontEx) );
    
    WContainerWidget *bottom = footer();
  //  bottom->setStyleClass("modal-footer");
    //bottom->resize( WLength::Auto, 50 );
    
    InterSpecApp *app = dynamic_cast<InterSpecApp *>(wApp);
    if( app && app->viewer() )
    {
      if( app->viewer()->isMobile() )
      {
        WPushButton *ancor = new WPushButton( bottom );
        ancor->setText( WString::tr("hw-welcome-link") );
        ancor->setStyleClass( "LinkBtn" );
        ancor->setFloatSide( Wt::Right );
        ancor->clicked().connect( boost::bind( &InterSpec::showWelcomeDialog, app->viewer(), true ) );
      }
      else
      {
          WAnchor *anchor = new WAnchor( WLink(), WString::tr("hw-tutorials-link"), bottom );
          anchor->setFloatSide( Wt::Left );
          anchor->addStyleClass( "InfoLink" );
          anchor->clicked().connect( boost::bind( &AuxWindow::hide, this ) );
          anchor->clicked().connect( boost::bind( &InterSpec::showWelcomeDialog, app->viewer(), true ) );
      }
    }//if( app && app->viewer() )
    
    
    WPushButton *ok = addCloseButtonToFooter();
    ok->clicked().connect( boost::bind( &AuxWindow::hide, this ) );
    
    if( !preselect.empty() || (app && app->isMobile()) )  //Keep keyboard from popping up
      ok->setFocus();
    
    rejectWhenEscapePressed();
    
    show();
    setMinimumSize(500,400);
    resizeScaledWindow(0.85, 0.85);
    resizeToFitOnScreen();
    centerWindow();
    //centerWindowHeavyHanded();  //why not?
  } //void HelpSystem::showHelpDialog()
  
  
  HelpWindow::~HelpWindow()
  {
  }
  
  void HelpWindow::setTopic( const std::string &preselect )
  {
    if( preselect.empty() )
      return;
    
    const auto pos = m_treeLookup.find(preselect);
    if( pos != end(m_treeLookup) )
      m_tree->select( pos->second ); // This will trigger the call to #HelpWindow::selectHelpToShow
    else
      passMessage( "The help instructions does not have an entry for " + preselect, 2 );
  }//void setTopic( const std::string &preselect )
  
  
  const std::string &HelpWindow::currentTopic() const
  {
    return m_displayedTopicName;
  }//const std::string &currentTopic() const;
  
  
  WTreeNode *first_child( WTreeNode *node )
  {
    const vector<WTreeNode *> &kids = node->childNodes();
    
    for( size_t i = 0; i < kids.size(); ++i )
      if( kids[i]->isVisible() )
        return first_child( kids[i] );
    return node;
  }
  
  void HelpWindow::handleArrowPress( const Wt::WKeyEvent e )
  {
    //This function not yet working!
    const Wt::Key key = e.key();
    
    const set<WTreeNode *> &selected = m_tree->selectedNodes();
    if( selected.empty() )
    {
      m_tree->select( first_child(m_tree->treeRoot()) );
      return;
    }
    
    WTreeNode *current = *selected.begin();
    WTreeNode *parent = current->parentNode();
    
    if( !parent || (parent==m_tree->treeRoot()) )
      return;
    
    const vector<WTreeNode *> &kids = current->childNodes();
    const vector<WTreeNode *> &siblings = parent->childNodes();
    const size_t pos = std::find(siblings.begin(), siblings.end(), current)
                       - siblings.begin();
    
    switch( key )
    {
      case Wt::Key_Left:
      {
        m_tree->select( parent );
        break;
      }//case Wt::Key_Left:
        
      case Wt::Key_Up:
      {
        WTreeNode *toselect = parent;
        
        for( int index = static_cast<int>(pos) - 1; !toselect && index >= 0; --index )
        {
          if( siblings[index]->isVisible() )
            toselect = siblings[index];
        }
        
        m_tree->select( toselect );
      }//case Wt::Key_Up:
        
      case Wt::Key_Right:
      {
        for( size_t i = 0; i < kids.size(); ++i )
        {
          if( kids[i]->isVisible() )
          {
            m_tree->select( kids[i] );
            break;
          }
        }
      }//case Wt::Key_Right:
        
      case Wt::Key_Down:
      {
        WTreeNode *firstkid = first_child( current );
        
        WTreeNode *toselect = (firstkid==current ? (WTreeNode *)0 : firstkid);
        
        for( size_t index = pos+1; !toselect && index < siblings.size(); ++index )
          if( siblings[index]->isVisible() )
            toselect = siblings[index];
      
        if( !toselect && parent->parentNode() )
        {
          const vector<WTreeNode *> &uncles = parent->parentNode()->childNodes();
          const size_t pos = std::find(uncles.begin(), uncles.end(), parent)
                             - uncles.begin();
          for( size_t i = pos; !toselect && i < uncles.size(); ++i )
          {
            if( uncles[i]->isVisible() )
              toselect = uncles[i];
          }
        }
        
        if( toselect )
          m_tree->select( toselect );
        
        break;
      }//case Wt::Key_Down:
        
      default:
        break;
    }//switch( key )
  }//void HelpWindow::handleArrowPress( const Wt::Key key )
  
  void HelpWindow::handleKeyPressInSearch( const Wt::WKeyEvent e )
  {
    switch( e.key() )
    {
      case Wt::Key_Left: case Wt::Key_Up: case Wt::Key_Right: case Wt::Key_Down:
        handleArrowPress( e );
        break;
      
      case Wt::Key_Home: case Wt::Key_Alt: case Wt::Key_Control:
      case Wt::Key_Shift: case Wt::Key_Tab: case Wt::Key_unknown:
      case Wt::Key_PageUp: case Wt::Key_PageDown:
        return;
        
      default:
        initialize();
        break;
    }//switch( e.key() )
  }//void handleKeyPressInSearch( const Wt::WKeyEvent e )
  
  void HelpWindow::initialize()
  {
    Wt::WTreeNode *node = new Wt::WTreeNode("Help");
    m_tree->setTreeRoot( node );
    node->label()->setTextFormat( Wt::PlainText );
    node->setLoadPolicy( Wt::WTreeNode::NextLevelLoading );
    m_tree->select( node ); //select first node
    
    try
    {
      Json::Value jsonResult;
      Json::parse( m_helpLookupTable, jsonResult );
      Json::Array res = jsonResult;
      populateTree( jsonResult, node );
      
      //Lets show the first result
      const vector<WTreeNode *> &kids = node->childNodes();
      if( !m_searchText->text().empty() && !kids.empty() )
        m_tree->select( first_child( kids[0] ) );
    }catch( std::exception &e )
    {
      cerr << "showHelpDialog() caught: " << e.what() << endl
      << "you may want to check help.json" << endl;
    }//try / catch
    
    node->setNodeVisible( false ); //hide the root note, not necessary
    
    m_tree->refresh();
  }//void initialize()
  
  
  
  void HelpWindow::selectHelpToShow()
  {
    m_tree->treeRoot()->setNodeVisible(false); //this is needed to keep the root hidden, otherwise it shows on click
    
    if( m_tree->selectedNodes().empty() )
      return;
    
    WTreeNode * node = *(m_tree->selectedNodes().begin());
    
    m_helpWindowContent->clear();
    
    string topic_name;
    Wt::WTemplate *content = nullptr;
    
    for( const auto &nameToNode : m_treeLookup )
    {
      if( nameToNode.second == node )
      {
        topic_name = nameToNode.first;
        if( topic_name == m_displayedTopicName )
        {
          cout << "Help content for '" << topic_name << "' should already be displayed." << endl;
          return;
        }
        
        content = getContentToDisplay( nameToNode.first );
        break;
      }//if( nameToNode.second == node )
    }//for( const auto &nameToNode : m_treeLookup )
    
    if( !content )
      content = new WTemplate( WString::tr("hw-err-opening-xml") );
    
    m_helpWindowContent->addWidget( content );
    
    //Note: this keeps the new pages scrolled to the top.
    doJavaScript("document.body.scrollTop = document.documentElement.scrollTop = 0;");
    
    const string prevTopic = m_displayedTopicName;
    m_displayedTopicName = topic_name;
    
    
    UndoRedoManager *undoRedo = UndoRedoManager::instance();
    if( undoRedo && undoRedo->canAddUndoRedoNow() && !prevTopic.empty() )
    {
      auto undo = [prevTopic](){
        InterSpec::instance()->showHelpWindow(prevTopic);
      };
      auto redo = [topic_name](){
        InterSpec::instance()->showHelpWindow(topic_name);
      };
      undoRedo->addUndoRedoStep( undo, redo, "Change help topic." );
    }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
  }//void selectHelpToShow()
  
  
  std::string HelpWindow::getHelpContents( const std::string &tag ) const
  {
    auto iter = m_treeLookup.find( tag );
    if( iter == end(m_treeLookup) )
      return "";
    
    WTreeNode *node = iter->second;
    
    std::string helpStr, filepath;
    auto pos = m_contentLookup.find( tag );
    if( pos != end(m_contentLookup) )
      filepath = pos->second;
    SpecUtils::trim( filepath );
    
    if( filepath.empty() )
    {
      const vector<WTreeNode *> &kids = node->childNodes();
      for( const WTreeNode *kid : kids )
      {
        for( const auto &nameToNode : m_treeLookup )
        {
          if( nameToNode.second == kid )
          {
            if( helpStr.empty() )
            {
              helpStr = getHelpContents( nameToNode.first );
            }else
            {
              helpStr += "<div style=\"margin-top: 25px;\">" + getHelpContents( nameToNode.first ) + "</div>";
            }
            
            break;
          }//if( we found the node that has appropriate key )
        }//for( const auto &nameToNode : m_treeLookup )
      }//for( const WTreeNode *kid : kids )
      
      if( helpStr.empty() )
        helpStr = "No help item associated with key '" + tag + "'";
    }else
    {
      //We have to grab the contents of this file, and place in helpStr...
      //All paths relative to the curent help.json file, which is in
      // "InterSpec_resources/static_text"
      
      //ToDo: should remove leading and trailing quotes I guess...
      
      const string docroot = wApp->docRoot();
      filepath = SpecUtils::append_path( "InterSpec_resources/static_text", filepath );
      filepath = SpecUtils::append_path(docroot, filepath);
      
      try
      {
        helpStr = getInternationalizedFileContents(filepath);
      }catch( std::exception &e )
      {
        helpStr = e.what();
      }
    }//if( filepath.empty() )
    
    return helpStr;
  }//std::string getHelpContents( const std::string &tag ) const
  
  
  Wt::WTemplate *HelpWindow::getContentToDisplay( const std::string &tag ) const
  {
    const string helpStr = getHelpContents( tag );
    
    WTemplate *t = new WTemplate();
    t->setTemplateText( WString(helpStr, UTF8), XHTMLUnsafeText );
    t->setInternalPathEncoding( true );
    
    return t;
  }//WTemplate *getContent( WTreeNode *node )
  
  void HelpWindow::populateTree(Json::Array &res, WTreeNode* parent)
  {
    InterSpecApp *app = dynamic_cast<InterSpecApp *>(wApp);
    
    const string searchTxt = m_searchText->text().toUTF8();
    const bool searching = (searchTxt.length() > 0);
    
    //go through each children element
    for( size_t i = 0; i < res.size(); ++i )
    { //res[i] is std::vector<Json::Value>
      const Json::Object &info = res[i];
      
      const Json::Value &ident = info.get("id");
      const Json::Value &title = info.get("title");
      const Json::Value &keywords = info.get("keywords");
      const Json::Value &children = info.get("children");
      
      if( ident.isNull() || title.isNull())
      {
        cerr << "found a null key or file string, "
        "you may want to check help.json" << endl;
        continue;
      }//if( id.isNull() || title.isNull())
      
      const WString identifier = ident;
      const WString titlestr = title;
    
      //This is a hack to keep help topics we arent using from showing - eh,
      //  not to clean, but whatever for now
      const string idstr = identifier.toUTF8();
      if( !USE_TERMINAL_WIDGET && SpecUtils::icontains(idstr,"terminal-dialog") )
      {
        continue;
      }else if( !USE_SPECRUM_FILE_QUERY_WIDGET && SpecUtils::icontains(idstr,"spectrum-file-query") )
      {
        continue;
      }else if( !USE_REMOTE_RID && SpecUtils::icontains(idstr,"external-rid") )
      {
        continue;
      }
      //Other potential topics to filter (if we add help content for them):
      //  USE_SEARCH_MODE_3D_CHART, USE_GOOGLE_MAP
        

      WTreeNode *insertNode = new Wt::WTreeNode( titlestr.toUTF8() );
      
      if( searching )
        insertNode->setHidden(true); //by default hidden, and just show nodes that need to be visible
      
      parent->addChildNode(insertNode);
      parent->expand();
      
      if( !keywords.isNull() && searching )
      {
        const WString keywordstr = keywords;
        
        typedef const boost::iterator_range<std::string::const_iterator> StringRange;
        std::string str1(keywordstr.toUTF8());
        std::string str2(searchTxt);
        SpecUtils::trim(str2);
        
        vector< std::string > searchWords;
        SpecUtils::split( searchWords, str2, " ," ); //split search terms
        
        bool childHide=false;
        for( size_t i = 0; i < searchWords.size(); ++i )
        {
          const std::string &keyword = searchWords[i];
          if (str1.length()==0 || (str1.length()>0 && !boost::ifind_first(StringRange(str1.begin(), str1.end()), StringRange(keyword.begin(), keyword.end()) )))
          {
            childHide=true;
            break;
          } //if found keyword in keywords
        } //for( std::string keyword : searchWords )
        
        if (!childHide)
          setPathVisible(insertNode);
      }//if (!keywords.isNull() && searching)

      const string identstr = identifier.toUTF8();
      m_treeLookup[identstr] = insertNode; //for future matching
      
      if( app )
      {
        WString xml_file;
        const Json::Value &xml = info.get("help-xml");
        const Json::Value &xml_desktop = info.get("help-xml-desktop");
        const Json::Value &xml_mobile = info.get("help-xml-mobile");
        const Json::Value &xml_phone = info.get("help-xml-phone");
        const Json::Value &xml_tablet = info.get("help-xml-tablet");
        
        if( app->isPhone() )
        {
          if( !xml_phone.isNull() )
            xml_file = xml_phone;
          else if( !xml_mobile.isNull() )
            xml_file = xml_mobile;
          else if( !xml.isNull() )
            xml_file = xml;
        }else if( app->isTablet() || app->isMobile() )
        {
          if( !xml_tablet.isNull() )
            xml_file = xml_tablet;
          else if( !xml_mobile.isNull() )
            xml_file = xml_mobile;
          else if( !xml.isNull() )
            xml_file = xml;
        }else  //desktop
        {
          if( !xml_desktop.isNull() )
            xml_file = xml_desktop;
          else if( !xml.isNull() )
            xml_file = xml;
        }
        
        m_contentLookup[identstr] = xml_file.toUTF8();
      }//if( app )
      
      
      if( !children.isNull() )
      {
        //has children
        Json::Array addNodes = children;
        populateTree(addNodes,insertNode);
      }//if (!children.isNull())
    }//for( size_t i = 0; i < res.size(); ++i )
  }//void populateTree(Json::Array &res, WTreeNode* parent, std::map <string, Wt::WTreeNode*> &treeLookup)
  
  /**
   Recursively set entire path from root to this node to be visible
   **/
  void HelpWindow::setPathVisible(Wt::WTreeNode* parent)
  {
    parent->setHidden(false);
    if (parent->parentNode()!=NULL)
      setPathVisible(parent->parentNode());
  }//  void HelpWindow::setPathVisible(WTreeNode* parent)
  
  

  /*
   Attach a tooltip onto a widget.  Depending on preference, it might show immediately.
   
   */
  void attachToolTipOn( Wt::WWebWidget* widget, const Wt::WString &text,
                        const bool enableShowing,
                        const ToolTipPosition position,
                        const ToolTipPrefOverride forceShowing )
  {
    attachToolTipOn( {widget}, text, enableShowing, position, forceShowing );
  }
  
  void attachToolTipOn( std::initializer_list<Wt::WWebWidget*> widgets,
                       const Wt::WString &text,
                        const bool enableShowing,
                        const ToolTipPosition position,
                        const ToolTipPrefOverride forceShowing )
  {
    assert( widgets.size() );
    if( !widgets.size() )
      return;
    
    //if is gesture controlled, we do not want to add tooltips to objects
    InterSpecApp *app = dynamic_cast<InterSpecApp *>( wApp );
    if (app && app->isMobile())
       return;
    
    Wt::WWebWidget *first_widget = *begin(widgets);
    
#if( USE_OSX_NATIVE_MENU )
    for( Wt::WWebWidget *w : widgets )
    {
      PopupDivMenuItem *menuitem = dynamic_cast<PopupDivMenuItem *>( w );
      if( menuitem && menuitem->getNsMenuItem() )
      {
        addOsxMenuItemToolTip( menuitem->getNsMenuItem(), text.toUTF8().c_str() );
        return;
      }//if( menuitem )
    }//for( Wt::WWebWidget *w : widgets )
#endif
    
    //get the InterSpecApp so we can access the value set (from the preferences) whether
    //we want to display the tooltips immediately for beginners.  Pro users also hide
    //the tooltip when they start typing.
    
    const bool overrideShow = (forceShowing == ToolTipPrefOverride::AlwaysShow);
    const bool instantShow = (forceShowing == ToolTipPrefOverride::InstantAlways);
    
    //Create popup notifications
    Wt::WStringStream strm;
    
    //need to escape the ' in the text message, and also remove new-lines from the string
    std::string val = text.toUTF8();
    boost::replace_all(val, "'", "\\'");
    boost::replace_all( val, "\r", "" );
    boost::replace_all( val, "\n", "" );
    
    string pos = "";
    
    switch( position )
    {
      case ToolTipPosition::Top:    pos = "my: 'bottom center', at: 'top center', "; break;
      case ToolTipPosition::Bottom: pos = "my: 'top center', at: 'bottom center', "; break;
      case ToolTipPosition::Right:  pos = "my: 'left center', at: 'right center', "; break;
      case ToolTipPosition::Left:   pos = "my: 'right center', at: 'left center', "; break;
    }//switch( position )
    
    string selector;
    for( const auto w : widgets )
      selector += (selector.empty() ? "#" : ",#") + w->id();
    
    RmHelpFromDom *remover = new RmHelpFromDom( first_widget );
    
    //Note: need to pre-render, as toggling requires tooltip already rendered.
    strm << "$('"<< selector <<"').qtip({ "
        "id: '" << remover->id() << "',"
        "prerender: true, "
        "content:  { text: '" << val << "'}, "
        "position: { " << pos <<
                    "viewport: $(window), "
                    "adjust: { "
                        "method: 'flipinvert flipinvert', "
                        "x:5} "
                    "},"
        "show:  {  event: '" << (enableShowing ? (string("mouseenter") + (overrideShow ? " focus" : "")) : string(""))
                << "', delay: " << string(instantShow ? "0" : "500") << " },"
        "hide:  {  fixed: true, event: 'mouseleave focusout" << string(instantShow ? "" : " keypress click") << "'  },"
        "style: { classes: 'qtip-rounded qtip-shadow" << string(overrideShow ? "" : " canDisableTt") << "',"
                  "tip: {"
                        "corner: true, "
                        "width: 18, "
                        "height: 12}"
                "}"
        "});";
    
    first_widget->doJavaScript( strm.str() );
  }//void attachToolTipOn( Wt::WWebWidget* widget, const std::string &text, bool force)

} //namespace HelpSystem
