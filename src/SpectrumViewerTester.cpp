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

#include <map>
#include <mutex>
#include <deque>
#include <vector>
#include <string>
#include <memory>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sys/stat.h>
#include <condition_variable>

#include <boost/ref.hpp>

#include <Wt/WText>
#include <Wt/WMenu>
#include <Wt/WImage>
#include <Wt/WTable>
#include <Wt/Dbo/Dbo>
#include <Wt/WServer>
#include <Wt/WPainter>
#include <Wt/WSvgImage>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WStringStream>
#include <Wt/WCssStyleSheet>
#include <Wt/WStackedWidget>
#include <Wt/WContainerWidget>
#include <Wt/Dbo/backend/Sqlite3>
#include <Wt/Test/WTestEnvironment>

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/PeakEdit.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/IsotopeId.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/InterSpec.h"
#include "SpecUtils/Filesystem.h"
#include "InterSpec/DecayWindow.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/UseInfoWindow.h"
#include "InterSpec/OneOverR2Calc.h"
#include "InterSpec/PeakFitChi2Fcn.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/PeakInfoDisplay.h"
#include "InterSpec/SpecFileSummary.h"
#include "InterSpec/GammaCountDialog.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/UnitsConverterTool.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/LocalTimeDelegate.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/CompactFileManager.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeSearchByEnergy.h"
#if( !ANDROID && !IOS )
#include "InterSpec/FileDragUploadResource.h"
#endif
#include "InterSpec/ShieldingSourceDisplay.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"
#include "InterSpec/GammaXsGui.h"

#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"


//#include <Wt/Dbo/WtSqlTraits>

#if( ALLOW_URL_TO_FILESYSTEM_MAP )
#include "InterSpec/DbToFilesystemLink.h"
#endif

#if( INCLUDE_ANALYSIS_TEST_SUITE )
#include "InterSpec/SpectrumViewerTester.h"
#endif

using namespace Wt;
using namespace std;

double SpectrumViewerTester::sm_mean_max_nsigma_diff      = 0.1;
double SpectrumViewerTester::sm_amplitude_max_frac_diff   = 0.05;
double SpectrumViewerTester::sm_width_max_frac_diff       = 0.05;
const double SpectrumViewerTester::sm_lower_roi_max_nsigma_diff = 0.1;
const double SpectrumViewerTester::sm_upper_roi_max_nsigma_diff = 0.1;
const double SpectrumViewerTester::sm_multipeakfit_nominal_weight = 0.5;
const double SpectrumViewerTester::sm_multipeakfit_range_max_vary_frac = 0.05;
const int    SpectrumViewerTester::sm_num_manual_click = 50;
const double SpectrumViewerTester::sm_nfwhm_manually_click = 0.75;
const double SpectrumViewerTester::sm_min_manually_click_corr_frac = 0.9;


typedef boost::function<void(WContainerWidget *)> TestDisplayCreateFunction;

namespace
{
  void select_item( WMenu *menu, SideMenuItem *item )
  {
    //boost::bind fails to compile for WMenu::select since it accepts a int or ptr
    menu->select( item );
    item->triggered().emit( item ); //doenst look like this is emmitted either
    //when body of SideMenuItem is clicked
  }//void select_item( WMenu *menu, SideMenuItem *item )
  
  std::pair<SideMenuItem *, DeferredWidget< TestDisplayCreateFunction > *>
                makeItem( SpectrumViewerTester::TestType type,
                          SpectrumViewerTester *tester,
                          WMenu *m )
  {
    TestDisplayCreateFunction f
          = boost::bind( &SpectrumViewerTester::doScorring, tester, type, _1 );
    
    DeferredWidget< TestDisplayCreateFunction > *w = deferCreate( f );
    w->addStyleClass( "UseInfoItem" );
    const char *title = SpectrumViewerTester::tostr( type );
    SideMenuItem *item = new SideMenuItem( title, w );
    item->clicked().connect( boost::bind( &select_item, m, item ) );
    m->addItem( item );
  
    std::pair<SideMenuItem *, DeferredWidget< TestDisplayCreateFunction > *> p;
    p.first = item;
    p.second = w;
    return p;
  }//makeItem(...)
  
}//namespace

bool SpectrumViewerTester::compatibleMeans( const PeakDef *original,
                                            const PeakDef *found )
{
  const double meanDiff = fabs(original->mean() - found->mean());
  const bool compatible = (meanDiff < sm_mean_max_nsigma_diff*original->sigma());
  //if( !compatible )
  //  cerr << "Not compatible mean: original->mean()=" << original->mean()
  //       << ", found->mean()=" << found->mean() << ", sigma=" << original->sigma() << endl;

  return compatible;
}

bool SpectrumViewerTester::compatibleAmplitudes( const PeakDef *original,
                                                 const PeakDef *found )
{
  const double ampDiff = fabs(original->amplitude() - found->amplitude());
  return ampDiff < sm_amplitude_max_frac_diff*original->amplitude();
}

bool SpectrumViewerTester::compatibleLowerRoi( const PeakDef *original,
                                               const PeakDef *found )
{
  const double lowRoiDiff = fabs(original->lowerX() - found->lowerX());
  return lowRoiDiff < sm_lower_roi_max_nsigma_diff*original->sigma();
}

bool SpectrumViewerTester::compatibleUpperRoi( const PeakDef *original,
                                               const PeakDef *found )
{
  const double highRoiDiff = fabs(original->upperX() - found->upperX());
  return highRoiDiff < sm_upper_roi_max_nsigma_diff*original->sigma();
}

bool SpectrumViewerTester::compatibleWidths( const PeakDef *orig,
                                             const PeakDef *found )
{
  const double widthDiff = fabs(orig->sigma() - found->sigma());
  return widthDiff < sm_width_max_frac_diff*orig->sigma();
}


bool SpectrumViewerTester::sameNuclideId( const PeakDef *orig,
                                          const PeakDef *found )
{
  return (orig->parentNuclide()==found->parentNuclide()
          && orig->reaction()==found->reaction()
          && orig->xrayElement()==found->xrayElement()
          && orig->nuclearTransition()==found->nuclearTransition()
          && orig->decayParticleIndex()==found->decayParticleIndex()
          && orig->sourceGammaType()==found->sourceGammaType()
          && fabs(orig->xrayEnergy()-found->xrayEnergy())<0.001
          && fabs(orig->reactionEnergy()-found->reactionEnergy())<0.001 );
}//sameNuclideId

#if( BUILD_AS_OFFLINE_ANALYSIS_TEST_SUITE )

void SpectrumViewerTester::writeHtmlHeader( std::ostream &strm,
                                            const string &title )
{
  strm << "<html>\n"
  << "\t<head>\n"
  << "\t\t<title>" << title << "</title>\n";
  
  //Try to directly include the CSS into the HTML file, but if not put a link to
  //  it
  ifstream incss( "../InterSpec_resources/SpectrumViewerTester.css" );
  if( !incss.is_open() )
  {
    cerr << "\n\n\n\nLinking to style sheet instead of including it!\n\n" << endl;
    strm << "\t\t<link rel=\"stylesheet\" type=\"text/css\""
            " href=\"../InterSpec_resources/SpectrumViewerTester.css\">\n";
  }else
  {
    strm  << "\t\t<style>\n";
    incss.seekg( 0, ios::end );
    const istream::pos_type eof_pos = incss.tellg();
    incss.seekg( 0, ios::beg );
    const size_t size = 0 + eof_pos;
    std::unique_ptr<char[]> buffer( new char [size+1] );
    incss.read( buffer.get(), size );
    strm.write( buffer.get(), size );
    
    strm << "\n\t\t</style>\n";
  }//if( incss.is_open() )

  strm << "\t</head>\n";
}//void writeHtmlHeader( std::ostream &strm )


void SpectrumViewerTester::writeStateStaticHtmlHeader( std::ostream &strm,
                                                       InterSpec *viewer,
                                                    const UserState *dbstate )
{
  const string username = dbstate->user ? dbstate->user->userName() : string("Unknown");
  const string now = WDateTime::currentDateTime().toString( DATE_TIME_FORMAT_STR ).toUTF8();
  
  writeHtmlHeader( strm, dbstate->name.toUTF8() );
  
  strm << "\t<center>\n"
       << "\t\t<h1>Test Results For " << dbstate->name.toUTF8() << "</h1>\n"
       << "\t\t<h2>Test Performed at " << now << "</h2>\n"
       << "\t\t<table class=\"StateTestStateInfoTable\">\n"
       << "\t\t\t<tr><th>Test Name</th><td>" << dbstate->name.toUTF8() << "</td></tr>\n"
       << "\t\t\t<tr><th>Test Creator</th><td>" << username << "</td></tr>\n"
       << "\t\t\t<tr><th>Creation Time</th><td>"
          << dbstate->creationTime.toString( DATE_TIME_FORMAT_STR ) << "</td></tr>\n"
       << "\t\t\t<tr><th>Last Modified</th><td>"
          << dbstate->serializeTime.toString( DATE_TIME_FORMAT_STR ) << "</td></tr>\n"
       << "\t\t\t<tr><th>Description</th><td>"
          << dbstate->description.toUTF8() << "</td></tr>\n"
       << "\t\t\t<tr><th>User State ID</th><td>"
          << viewer->m_currentStateID << "</td></tr>\n"
       << "\t\t\t<tr><th>Foreground ID</th><td>"
          << dbstate->foregroundId << "</td></tr>\n"
       << "\t\t\t<tr><th>Background ID</th><td>"
          << dbstate->backgroundId << "</td></tr>\n"
       << "\t\t\t<tr><th>1D Model ID</th><td>"
         << dbstate->shieldSourceModelId << "</td></tr>\n"
       << "\t\t</table>\n"
       << "\t</center>\n";
}//void writeStateStaticHtmlHeader( std::ostream &strm )


void SpectrumViewerTester::doOfflineTesting()
{
  vector<int> userStateIds;
  
  {//begin getting userStates
    Wt::Test::WTestEnvironment env( Wt::Application );
    std::unique_ptr<InterSpecApp> app( new InterSpecApp( env ) );
    InterSpec *viewer = app->m_viewer;
    std::shared_ptr<DataBaseUtils::DbSession> session = viewer->sql();

//    const string &dbLocation = DataBaseUtils::preferenceDatabaseFile();
//    Wt::Dbo::SqlConnection &database = *viewer->m_database;
//    session.setConnection( database );
//    mapDbClasses( &session );
//    try { session.createTables(); }catch(...){}
  
    char query[256];
    snprintf( query, sizeof(query), "StateType=%i", int(UserState::kForTest) );
  
    {//begin interaction with DB
      DataBaseUtils::DbTransaction transaction( *session );
      Dbo::collection<Dbo::ptr<UserState> > states = session->session()->find<UserState>().where(query);
  
      for( Dbo::collection<Dbo::ptr<UserState> >::iterator iter = states.begin();
          iter != states.end(); ++iter )
        userStateIds.push_back( iter->id() );
      transaction.commit();
    }//end interaction with DB
  }//end getting userStates
  
  const WDateTime now = WDateTime::currentDateTime();
  
  char todaystr[64], basedir[128];
  snprintf( todaystr, sizeof(todaystr), "%.4i%.2i%.2i",
            now.date().year(), now.date().month(), now.date().day() );
  snprintf( basedir, sizeof(basedir), "automated_test_%s", todaystr );
  

  if( !SpecUtils::is_directory(basedir) )
    SpecUtils::create_directory( basedir );
  
  int nStatesTested = 0;
  vector<string> medResRultFiles, highResResultFiles;
  map<string,string> htmlFileToTestName, fileToScoreSummary, fileToHttpRef;
  Score overallscore, overalmedres, overallhighres;
  Score allscore[NumTestType], medres[NumTestType], highres[NumTestType];
  
  for( size_t i = 0; i < userStateIds.size(); ++i )
  {
    const int dbstateid = userStateIds[i];
    Wt::Test::WTestEnvironment env( Wt::Application );
    std::unique_ptr<InterSpecApp> app( new InterSpecApp( env ) );
    InterSpec *viewer = app->m_viewer;
    
    Dbo::ptr<UserState> dbstate;
    std::shared_ptr<DataBaseUtils::DbSession> sql = viewer->sql();
    
    string username;
    
    {//begin code block to interact with database
      
      DataBaseUtils::DbTransaction transaction( *sql );
      
      dbstate = sql->session()->find<UserState>().where("id=?").bind( dbstateid );
      assert( dbstate );
      if( dbstate->user )
        username = dbstate->user->userName();
      viewer->loadStateFromDb( dbstate );
      transaction.commit();
    }//end code block to interact with database
    
    if( !viewer->m_dataMeasurement )
    {
      cerr << "State has no forground!" << endl;
      continue;
    }//if( !viewer->m_dataMeasurement )
    
    ++nStatesTested;
    const int width = 1024, height = 768;
    std::unique_ptr<SpectrumViewerTester> tester;
    tester.reset( new SpectrumViewerTester( viewer, width, height ) );
    
    char outfilename[128], fileName[256];
    snprintf( outfilename, sizeof(outfilename), "test_%i.html", int(i) );
    snprintf( fileName, sizeof(fileName), "%s/%s", basedir, outfilename );
    
    ofstream outfile( fileName );
    
    const bool isHigRes =(viewer->m_dataMeasurement->num_gamma_channels()>2048);
    
    Score sumscore;
    Score scores[NumTestType];
    for( TestType t = TestType(0); t < NumTestType; t = TestType(t+1) )
    {
      scores[t] = tester->doTest( t );
      const double w = scoreWeight(t);
      sumscore.addScore( scores[t], w );
      allscore[t].addScore( scores[t], 1.0 );
      
      if( isHigRes )
        highres[t].addScore( scores[t], 1.0 );
      else
        medres[t].addScore( scores[t], 1.0 );
    }//for( TestType t = TestType(0); t < NumTestType; t = TestType(t+1) )

    overallscore.addScore( sumscore, 1.0 );
    if( isHigRes )
    {
      overallhighres.addScore( sumscore, 1.0 );
      highResResultFiles.push_back( outfilename );
    }else
    {
      overalmedres.addScore( sumscore, 1.0 );
      medResRultFiles.push_back( outfilename );
    }
    
    
    
    writeStateStaticHtmlHeader( outfile, viewer, dbstate.get() );
    
    {//begin codeblock to write a summarry score
      const double corrPercent = (100.0*sumscore.m_ncorrect)/sumscore.m_ntest;
      stringstream scoreSum;
      scoreSum << sumscore.m_ncorrect << " of " << sumscore.m_ntest
               << " tests correct (" << std::setprecision(1) << std::fixed
               << corrPercent << "%); " << sumscore.m_nwrong
               << " issues";
      
      outfile << "\t<div class=\"OverallStaticScoreBox\">Overall "
              << scoreSum.str() << " found</div>\n";
      
      htmlFileToTestName[outfilename] = dbstate->name.toUTF8();
      
      fileToScoreSummary[outfilename] = scoreSum.str();
      
      const string dbidstr = std::to_string(dbstate.id());
      string url = "http://localhost:8080?teststate=" + dbidstr;
//      if( !wApp || wApp->environment().isTest() )
//      {
//#if( !BUILD_AS_OFFLINE_ANALYSIS_TEST_SUITE )
//        if( wApp && WServer::instance() )
//          url = wApp->environment().hostName() + "?teststate=" + dbidstr;
//#endif
        outfile << "\t\t<div class=\"ToTestFromResultLink\">Link to test state: "
                << "<a href=\"" << url << "\" class=\"NewTabLink\" target=\"_blank\">"
                << url << "</a></div>\n";
//      }//if( !wApp || wApp->environment().isTest() )
      
      fileToHttpRef[outfilename] = url;
    }//end codeblock to write a summarry score
    
    
    //Make a list of links for easy access to each section of the report
    outfile << "\t<div class=\"StaticLinkTOC\">"
            << "\t<div class=\"StaticLinkTOCTitle\">Tests</div>\n"
            << "\t\t<ul class=\"StaticLinkToSubResultsList\">\n"
            << "\t\t\t<li><a href=\"#SpectrumOverview\">Spectrum Overview</li>\n";
    
    for( TestType t = TestType(0); t < NumTestType; t = TestType(t+1) )
    {
      outfile << "\t\t\t<li><a href=\"#Test" << t << "\"><div>"
              << tostr(t) << " ("
              << scores[t].m_ncorrect << " of " << scores[t].m_ntest
              << " correct, " << scores[t].m_nwrong << " issues)</div>"
              << "</a></li>\n";
    }//for( TestType t = TestType(0); t < NumTestType; t = TestType(t+1) )
    
    outfile << "\t\t</ul>\n"
            << "\t</div>";
    

    {//begin codeblock to write spectrum overview image
      const int w = static_cast<int>(1.8 * tester->m_picWidth);
      const int h = static_cast<int>(1.2 * tester->m_picHeight);
      std::shared_ptr<WSvgImage> overviewSvg
                                 = tester->renderChartToSvg( -1.0, -1.0, w, h );
      outfile << "\t<a name=\"SpectrumOverview\"></a>\n"
              << "\t<div class=\"StaticSpectrumOverview\""
                 " style=\"width:"<< w << "px;\">\n"
              << "\t\t<div class=\"StaticSpectrumOverviewText\">Spectrum Overview</div>\n"
              << "\t\t<img>\n";
      overviewSvg->write( outfile );
      outfile << "\t\t</img>\n"
              << "\t</div>\n";
    }//end codeblock to write spectrum overview image
    
    const string getel = "document.getElementById('";
    const string dispnone = "').style.display='none'; ";
    const string dispshow = "').style.display=''; ";
    
    outfile << "<div class=\"StateTestStaticResults\""
               " style=\"width:" << (width-140) << "\">";
    for( TestType t = TestType(0); t < NumTestType; t = TestType(t+1) )
    {
      const string testnumstr = std::to_string(int(t));
      const string contentId = "TestContent" + testnumstr;
      const string hideId = "HideButton" + testnumstr;
      const string showId = "showButton" + testnumstr;
      
      ScoreDisplay display( scores[t] );
      
      outfile << "\t<a name=\"Test" << testnumstr << "\"></a>\n"
              << "<div class=\"ScoreDisplayStatic\"><div class=\"ScoreDisplayTitleStatic\">"
              << "<button onclick=\""
                 << getel << contentId << dispnone
                 << getel << hideId << dispnone
                 << getel << showId << dispshow
              << "\" id=\"" << hideId << "\" class=\"ScoreDisplayNoteContentShow\"></button>\n"
              << "<button onclick=\""
                 << getel << contentId << dispshow
                 << getel << hideId << dispshow
                 << getel << showId << dispnone
              << "\" id=\"" << showId << "\" style=\"display:none;\""
                 " class=\"ScoreDisplayNoteContentHide\"></button>\n"
             << "\t\t<span>" << tostr(t) << "</span></div>\n"
             << "\t<div id=\"" << contentId << "\" "
                "class=\"ScoreDisplayNoteContentStatic\">\n";
      
      display.renderScore( outfile );
      
      outfile << "\t</div>\n";
    }//for( TestType t = TestType(0); t < NumTestType; t = TestType(t+1) )
  
    outfile << "</div>\n"
            << "</html>";
  }//for( size_t i = 0; i < userStates.size(); ++i )
  
  
//We now have the following information to construct the main page from, that
//  will also link to the pages for each subpage
  const string indexfilename = basedir + string("/index.html");
  ofstream indexfile( indexfilename.c_str() );
  
  char title[256];
  snprintf( title, sizeof(title), "InterSpec Tests, %s" , todaystr );
  
  writeHtmlHeader( indexfile, title );
  
  indexfile << "\t<div class=\"ScoreOverviewTitle\">InterSpec Test Results from "
            << todaystr << "</div>\n"
            << "\t<div class=\"ScoreOverviewStaticSummarry\">\n"
            << "\t<div class=\"ScoreOverviewStaticTitle\">" << htmlFileToTestName.size()
            << " test states have been examined.</div>\n"
            << "\t<div class=\"StaticSummarryScore\">\n"
            << "\t\t<div><b>All Spectra:</b> "
            << std::setprecision(1) << std::fixed
            << overallscore.m_ncorrect << " of "
            << overallscore.m_ntest << " tests correct ("
            << ((100.0*overallscore.m_ncorrect)/overallscore.m_ntest) << "%), "
            << overallscore.m_nwrong << " issues found</div>\n"
            << "\t\t<div><b>Medium Resolution Spectra:</b> "
            << overalmedres.m_ncorrect << " of "
            << overalmedres.m_ntest << " tests correct ("
            << ((100.0*overalmedres.m_ncorrect)/overalmedres.m_ntest) << "%), "
            << overalmedres.m_nwrong << " issues found</div>\n"
            << "\t\t<div><b>High Resolution Spectra:</b> "
            << overallhighres.m_ncorrect << " of "
            << overallhighres.m_ntest << " tests correct ("
            << ((100.0*overallhighres.m_ncorrect)/overallhighres.m_ntest) << "%), "
            << overallhighres.m_nwrong << " issues found</div>\n"
            << "\t</div>\n"
            << "\t<table class=\"StaticSummarryTable\">\n"
            << "\t\t<tr><th>Test Type</th><th>All Spectra</th>"
               "<th>Medium Resolution</th><th>High Resolution</th></tr>\n";
  for( TestType t = TestType(0); t < NumTestType; t = TestType(t+1) )
  {
    indexfile << "\t\t<tr><th>" << tostr(t) << "</th>"
              << "<td>"
              << allscore[t].m_ncorrect << " of "
              << allscore[t].m_ntest << " correct ("
              << ((100.0*allscore[t].m_ncorrect)/allscore[t].m_ntest) << "%), "
              << allscore[t].m_nwrong << " issues"
              << "</td>"
              << "<td>"
              << medres[t].m_ncorrect << " of "
              << medres[t].m_ntest << " correct ("
              << ((100.0*medres[t].m_ncorrect)/medres[t].m_ntest) << "%), "
              << medres[t].m_nwrong << " issues"
              << "</td>"
              << "<td>"
              << highres[t].m_ncorrect << " of "
              << highres[t].m_ntest << " correct ("
              << ((100.0*highres[t].m_ncorrect)/highres[t].m_ntest) << "%), "
              << highres[t].m_nwrong << " issues"
              << "</td></tr>\n";
  }//for( TestType t = TestType(0); t < NumTestType; t = TestType(t+1) )
  
  indexfile << "\t</table>\n"
            << "\t<div class=\"StaticIndividualTestTitle\">Links to individual test</div>\n"
            << "\t<table class=\"StaticLinkToIndividualTestsTable\">\n"
            << "\t\t<tr><th colspan=\"3\" align=\"center\">Medium Resolution Individual Tests</th></tr>\n"
            << "\t\t<tr><th>Test Name/Link</th><th>Score Summary</th><th>Localhost Link</th></tr>\n";
  for( const string &name : medResRultFiles )
  {
    indexfile << "\t\t<tr><td><a href=\"" << name << "\""
              << " class=\"NewTabLink\" target=\"_blank\">" << htmlFileToTestName[name] << "</a></td>"
              << "<td>" << fileToScoreSummary[name] << "</td>"
              << "<td><a href=\"" << fileToHttpRef[name] << "\" class=\"NewTabLink\" target=\"_blank\">"
              << fileToHttpRef[name] << "</a></td>\n";
  }//for( const string &name : medResRultFiles )

  indexfile << "\t</table>\n"
            << "\t<table class=\"StaticLinkToIndividualTestsTable\">\n"
            << "\t\t<tr><th colspan=\"3\" align=\"center\">High Resolution Individual Tests</th></tr>\n"
            << "\t\t<tr><th>Test Name/Link</th><th>Score Summary</th><th>Localhost Link</th></tr>\n";
  for( const string &name : highResResultFiles )
  {
    indexfile << "\t\t<tr><td><a href=\"" << name << "\""
              << " class=\"NewTabLink\" target=\"_blank\">" << htmlFileToTestName[name] << "</a></td>"
              << "<td>" << fileToScoreSummary[name] << "</td>"
              << "<td><a href=\"" << fileToHttpRef[name] << "\" class=\"NewTabLink\" target=\"_blank\">"
              << fileToHttpRef[name] << "</a></td>\n";
  }//for( const string &name : medResRultFiles )
  
  indexfile << "\t</table>\n"
            << "</html>\n";
  
}//void SpectrumViewerTester::doOfflineTesting()
#endif //BUILD_AS_OFFLINE_ANALYSIS_TEST_SUITE


SpectrumViewerTester::Score::Score()
  : m_ntest(0.0),
    m_nwrong(0.0),
    m_ncorrect(0.0),
    m_dbid(-1),
    m_test( SpectrumViewerTester::NumTestType )
{
}//Score constructor


void SpectrumViewerTester::Score::addScore( const Score &rhs, double w )
{
  m_ntest    += w*rhs.m_ntest;
  m_nwrong   += w*rhs.m_nwrong;
  m_ncorrect += w*rhs.m_ncorrect;
}//void addScore( const Score &rhs, double w )


const char *SpectrumViewerTester::tostr( TestType type )
{
  switch( type )
  {
    case SpectrumViewerTester::AutomatedPeakSearch:
      return "Auto Peak Search";
      
    case SpectrumViewerTester::ManualPeakClicking:
      return "Manual Peak Clicking";
      
    case SpectrumViewerTester::MultiplePeakInRoiFit:
      return "Multiple Peak Fit";
      
    case SpectrumViewerTester::SourceShieldingFit:
      return "Shield/Source Fit";
      
    case SpectrumViewerTester::NumTestType:
      break;
  }//switch( type )
  
  cerr << "SpectrumViewerTester::tostr type=" << type << endl;
  throw runtime_error( "SpectrumViewerTester::tostr: invalid TestType" );
  return "";
}//const char *tostr( TestType type )


SpectrumViewerTesterWindow::SpectrumViewerTesterWindow( InterSpec *viewer )
  : AuxWindow( "InterSpec Tester" ),
    m_tester( 0 )
{
  addStyleClass( "SpectrumViewerTesterWindow" );
  
  double width  = 0.75*viewer->renderedWidth();
  double height = 0.8*viewer->renderedHeight();
  
  width  = std::max( 600.0, width );
  height = std::max( 450.0, height );
  
  
  m_tester = new SpectrumViewerTester( viewer,
                                       static_cast<int>(width),
                                       static_cast<int>(height) );
  m_tester->setSizeHint( static_cast<int>(width), static_cast<int>(height) );
  
  m_tester->m_done.connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );
  finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );
  
  WGridLayout *layout = stretcher();
  layout->addWidget( m_tester, 0, 0 );
  layout->setContentsMargins( 0, 0, 0, 0 );
  
  show();
  
  disableCollapse();
  
  resize( width, height );
  
  setResizable( true );
//  resizeToFitOnScreen();
  centerWindow();
}//SpectrumViewerTesterWindow

SpectrumViewerTester::ScoreNote::ScoreNote()
  : m_text(""),
    m_title(""),
    m_testImage(),
    m_originalImage()
{}


SpectrumViewerTester::ScoreNote::~ScoreNote()
{
//  if( m_testImage && m_testImage->parent() )
//  {
//    m_testImage.release();
//  }
}//~ScoreNote()


SpectrumViewerTester::ScoreDisplay::ScoreDisplay( const Score &score,
                                                Wt::WContainerWidget *parent )
: WContainerWidget( parent ),
  m_score( score ),
  m_notes( 0 )
{
  addStyleClass( "ScoreDisplay" );
  WText *title = new WText( tostr(m_score.m_test), this );
  title->addStyleClass( "ScoreDisplayTitle" );
  title->setInline( false );
  
  m_notes = new WContainerWidget( this );
  m_notes->addStyleClass( "ScoreDisplayNotes" );
  
  stringstream html;
  renderScore( html );
  m_notes->clear();
  new WText( html.str(), XHTMLUnsafeText, m_notes );
//  cerr << "\n\n\nHtmlText" << endl;
//  htmlText( cerr );
//  cerr << "EndHtmlText\n\n\n" << endl;
}//ScoreDisplay constructor


void SpectrumViewerTester::ScoreDisplay::setScore( const Score &score )
{
  m_score = score;
  
  stringstream html;
  renderScore( html );
  m_notes->clear();
  new WText( html.str(), XHTMLText, m_notes );
}//class ScoreDisplay


void SpectrumViewerTester::ScoreDisplay::renderScore( std::ostream &strm )
{
  const double corrPercent = (100.0*m_score.m_ncorrect)/m_score.m_ntest;
  
  strm << "<div class=\"ScoreDisplayNotes\" id=\"" << id() << "score\">\n"
  << "\t<div class=\"ScoreDisplayNoteTitle\">\n"
  << "\t\t<div class=\"ScoreSummarryTable\">\n"
  << "\t\t\t<div>" << m_score.m_ncorrect << " of "
                   << m_score.m_ntest << " tests correct ("
                   << std::setprecision(1) << std::fixed
                   << corrPercent << "%)</div>\n"
  << "\t\t\t<div>" << m_score.m_nwrong << " issues found</div>\n"
  << "\t\t</div>\n"
  << "\t</div>\n";
  
  int noten = 0;
  for( const ScoreNote &note : m_score.m_notes )
  {
    //Well put the entire contents, minus the title in a div with unique ID
    //  number, so that later on we can implement a method to collapse/expand
    //  the results section
    ++noten;
    const string contentId = id() + "N" + std::to_string(noten);
    const string hideId = id() + "H" + std::to_string(noten);
    const string showId = id() + "S" + std::to_string(noten);
    
    strm << "\t<div class=\"ScoreDisplayNoteItem\">\n"
         << "\t\t<div class=\"ScoreDisplayNoteTitle\">"
         << "<button onclick=\""
              "document.getElementById('" << contentId << "').style.display='none';"
              "document.getElementById('" << hideId << "').style.display='none';"
              "document.getElementById('" << showId << "').style.display='';\" "
              "id=\"" << hideId << "\" class=\"ScoreDisplayNoteContentShow\"></button>\n"
          << "<button onclick=\""
              "document.getElementById('" << contentId << "').style.display='';"
              "document.getElementById('" << hideId << "').style.display='';"
              "document.getElementById('" << showId << "').style.display='none';\" "
              "id=\"" << showId << "\" style=\"display:none;\""
              " class=\"ScoreDisplayNoteContentHide\"></button>\n"
         << "<span>" << note.m_title << "</span></div>\n"
         << "\t\t<div id=\"" << contentId << "\">\n"
         << "\t\t<table class=\"ScoreDisplayNoteTable\">\n"
         << "\t\t\t<tbody>\n"
         << "\t\t\t\t<tr>\n";

    if( note.m_testImage.get() )
    {
      strm << "\t\t\t\t\t<td class=\"ScoreDisplayNoteTestImage\">\n";
      strm << "\t\t\t\t\t\t<div class=\"ScoreImageTitle\">Automated Result</div>\n";
      strm << "\t\t\t\t\t\t<img>\n";
      note.m_testImage->write( strm );
      strm << "\t\t\t\t\t\t</img>\n";
      strm << "\t\t\t\t\t</td>\n";
      
      if( !note.m_originalImage.get() )
        strm << "\t\t\t\t\t<td>No Analysts Image Available</td>\n";
    }//if( note.m_testImage )
    
    if( note.m_originalImage.get() )
    {
      strm << "\t\t\t\t\t<td class=\"ScoreDisplayNoteOriginalImage\">\n";
      strm << "\t\t\t\t\t\t<div class=\"ScoreImageTitle\">Analysts Result</div>\n";
      strm << "\t\t\t\t\t\t<img>\n";
      note.m_originalImage->write( strm );
      strm << "\t\t\t\t\t\t</img>\n";
      strm << "\t\t\t\t\t</td>\n";
      
      if( !note.m_testImage.get() )
        strm << "\t\t\t\t\t<td>Test Image Not Applicable</td>\n";
    }//if( m_originalImage )
    
    strm << "\t\t\t\t</tr>\n"
         << "\t\t\t</tbody>\n"
         << "\t\t</table>\n"
         << "\t\t<div class=\"ScoreDisplayNoteText\">"
         << note.m_text << "</div>\n"
         << "\t\t</div>\n"
         << "\t</div>\n";
  }//for( const string &note : m_score.m_notes )
  
  strm << "</div>\n";
}//void renderScore()


const SpectrumViewerTester::Score &SpectrumViewerTester::ScoreDisplay::score() const
{
  return m_score;
}


void SpectrumViewerTester::updateOverview()
{
  if( !m_overview )
    return;
  
  char onstack[256];
  
  Score score;
  score.m_test = NumTestType;
  m_overview->clear();
  
  WText *txt = new WText( "Score Overview", m_overview );
  txt->setInline( false );
  txt->addStyleClass( "ScoreOverviewTitle" );

  int ncomplete = 0;
  for( TestType t = TestType(0); t < NumTestType; t = TestType(t+1) )
  {
    if( m_scoreDisplays[t] )
    {
      ncomplete += 1;
      const double w = scoreWeight(t);
      score.addScore( m_scores[t], w );
    }//if( m_scoreDisplays[t] )
  }//for( loop over TestType )

  if( ncomplete != int(NumTestType) )
  {
    snprintf( onstack, sizeof(onstack),
             "Have completed %i of %i tests, click buttons on left to run tests",
             ncomplete, int(NumTestType) );
    txt = new WText( onstack, m_overview );
    txt->setInline( false );
    txt->addStyleClass( "ScoreOverviewTestsComplete" );
  }//if( ncomplete != int(NumTestType) )
  
  if( ncomplete )
  {
    const double corrPercent = (100.0*score.m_ncorrect)/score.m_ntest;
    stringstream strm;
    strm << "\t\t<div class=\"OverViewScoreSummarryTable\">\n"
         << "\t\t\t<div>" << score.m_ncorrect << " of "
         << score.m_ntest << " tests correct ("
         << std::setprecision(1) << std::fixed
         << corrPercent << "%)</div>\n"
         << "\t\t\t<div>" << score.m_nwrong << " issues found</div>\n"
         << "\t\t</div>\n";
    /*txt =*/ new WText( strm.str(), XHTMLText, m_overview );
  }//if( ncomplete )
  
//  if( ncomplete != int(NumTestType) )
//  {
//    WPushButton *b = new WPushButton( "Run All Tests", m_overview );
//    b->addStyleClass( "RunAllTestsButton" );
//    b->clicked().connect( this, &SpectrumViewerTester::scoreAllTests );
//  }//if( ncomplete != int(NumTestType) )
}//void updateOverview()


double SpectrumViewerTester::scoreWeight( SpectrumViewerTester::TestType type )
{
  switch( type )
  {
    case SpectrumViewerTester::AutomatedPeakSearch:
    case SpectrumViewerTester::ManualPeakClicking:
    case SpectrumViewerTester::MultiplePeakInRoiFit:
    case SpectrumViewerTester::SourceShieldingFit:
      return 1.0;
    case SpectrumViewerTester::NumTestType:
      break;
  }//switch( type )
  
  return 1.0;
}//double scoreWeight( TestType type );


void SpectrumViewerTester::scoreAllTests()
{
  for( TestType t = TestType(0); t < NumTestType; t = TestType(t+1) )
  {
    m_menuItems[t]->contents()->load();
  }//for( loop over TestType )
}//void scoreAllTests()


SpectrumViewerTester::Score SpectrumViewerTester::doTest( SpectrumViewerTester::TestType type )
{
  Score score;
  
  switch( type )
  {
    case SpectrumViewerTester::AutomatedPeakSearch:
      score = testAutomatedPeakSearch();
      break;
      
    case SpectrumViewerTester::ManualPeakClicking:
      score = testManualPeakClicking();
      break;
  
    case SpectrumViewerTester::MultiplePeakInRoiFit:
      score = testMultiplePeakFit();
      break;
      
    case SpectrumViewerTester::SourceShieldingFit:
      score = testShieldSourceFit();
      break;
      
    case SpectrumViewerTester::NumTestType:
      throw runtime_error( "SpectrumViewerTester::doTest: invalid TestType" );
      break;
  }//switch( type )
  #if( USE_DB_TO_STORE_SPECTRA )
  score.m_dbid = m_viewer->m_currentStateID;
  #endif
  return score;
}//Score doTest( TestType type )


void SpectrumViewerTester::doScorring( SpectrumViewerTester::TestType type,
                                        Wt::WContainerWidget *parent )
{
  m_scores[type] = doTest(type);
  m_scoreDisplays[type] = new ScoreDisplay( m_scores[type], parent );
  if( m_menuItems[type] )
  {
    m_menuItems[type]->enable();
    m_menuItems[type]->removeStyleClass( "TestNotDone" );
  }//if( m_menuItems[type] )
  
  updateOverview();
}//doScorring(...)


SpectrumViewerTester::SpectrumViewerTester( InterSpec *viewer,
                                            int width, int height,
                                            WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_viewer( viewer ),
    m_overview( 0 ),
    m_sideMenu( 0 ),
    m_picWidth( 500 ),
    m_picHeight( 375 )
{
  setSizeHint( width, height );
  
  for( TestType t = TestType(0); t < NumTestType; t = TestType(t+1) )
  {
    m_scoreDisplays[t] = 0;
    m_menuItems[t] = 0;
  }//for( loop over TestType )
  
  if( wApp )
    wApp->useStyleSheet( "InterSpec_resources/SpectrumViewerTester.css" );
  
  WStackedWidget *stack = new WStackedWidget();
  stack->addStyleClass( "UseInfoStack" );
  WAnimation animation(Wt::WAnimation::Fade, Wt::WAnimation::Linear, 200);
  stack->setTransitionAnimation( animation, true );
  
  m_sideMenu = new WMenu( stack, Wt::Vertical );
  m_sideMenu->addStyleClass( "VerticalNavMenu HeavyNavMenu SideMenu" );
  
  WGridLayout *layout = new WGridLayout();
  setLayout( layout );
  layout->addWidget( m_sideMenu, 0, 0, AlignLeft );
  layout->addWidget( stack, 0, 1 );
  layout->setContentsMargins( 0, 0, 0, 0 );
  
  
  m_overview = new WContainerWidget();
  m_overview->addStyleClass( "ScoreOverview" );
  m_overview->addStyleClass( "UseInfoItem" );
  SideMenuItem *overviewItem = new SideMenuItem( "Overview", m_overview );
  overviewItem->clicked().connect( boost::bind( &select_item, m_sideMenu, overviewItem ) );
  
  m_sideMenu->addItem( overviewItem );
  
  for( TestType type = TestType(0); type < NumTestType; type= TestType(type+1) )
  {
    std::pair<SideMenuItem *, DeferredWidget< TestDisplayCreateFunction > *> iw;
    iw = makeItem( type, this, m_sideMenu );
    m_menuItems[type] = iw.first;
//    iw.second
    m_menuItems[type]->addStyleClass( "TestNotDone" );
  }//for( loop over TestType )
  
  updateOverview();
  m_sideMenu->select( 0 );
  layout->setRowStretch( 0, 1 );
  layout->setColumnStretch( 1, 1 );
}//SpectrumViewerTester constuctor


SpectrumViewerTester::~SpectrumViewerTester()
{
}//~SpectrumViewerTester()



bool SpectrumViewerTester::peaksAreSimilar( const PeakDef *l, const PeakDef *r )
{
  return compatibleMeans(l,r)
         && compatibleAmplitudes(l,r)
         && compatibleLowerRoi(l,r)
         && compatibleUpperRoi(l,r)
         && compatibleWidths(l,r);
}//bool peaksAreSimilar(...)


void SpectrumViewerTester::setSizeHint( int width, int height )
{
  m_picWidth = (width - 150 - 60) / 2;
  
  double hToW = 480.0/640.0;
  m_picHeight = static_cast<int>(hToW*m_picWidth);
}//void setSizeHint( int width, int height )


std::shared_ptr<Wt::WSvgImage> SpectrumViewerTester::renderChartToSvg( double lowx,
                                                       double upperx,
                                                       int width, int height )
{
  auto hist = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  
  shared_ptr<const deque< PeakModel::PeakShrdPtr > > peaks = m_viewer->peakModel()->peaks();
  std::vector<std::shared_ptr<const ReferenceLineInfo>> reflines;
  auto theme = m_viewer->getColorTheme();
  const bool compact = false;
  
  return PeakSearchGuiUtils::renderChartToSvg( hist, peaks, reflines, lowx, upperx, width,
                                               height, theme, compact );
}//renderChartToSvg(...)


string SpectrumViewerTester::makePeakSummarryTable( const PeakDef &peak,
                                        const std::shared_ptr<const SpecUtils::Measurement> &data )
{
  stringstream strm;
  
  const string startrow = "<tr><td class=\"PeakSummarryTableCell\"><b>";
  const string starttd = "</b></td><td class=\"PeakSummarryTableCell\">";
  
  strm << std::setprecision(3) << std::fixed <<
  "<table class=\"PeakSummarryTable\">\n"
    << startrow << "Mean" << starttd << peak.mean() << " keV</td></tr>\n"
    << startrow << "FWHM" << starttd << 2.35482*peak.sigma() << " keV</td></tr>\n"
    << startrow << "Area" << starttd << peak.amplitude() << "</td></tr>\n"
    << startrow << "Type" << starttd;
  
  switch( peak.type() )
  {
    case PeakDef::GaussianDefined:
      strm << "Gaussian</td></tr>\n";
      break;
    case PeakDef::DataDefined:
      strm << "Data Defined</td></tr>\n";
      break;
  }//switch( peak.type() )
  
  
  strm
    << startrow << "Lower ROI" << starttd << peak.lowerX() << " keV</td></tr>\n"
    << startrow << "Upper ROI" << starttd << peak.upperX() << " keV</td></tr>\n"
    << startrow << "Continuum" << starttd;

  switch( peak.continuum()->type() )
  {
    case PeakContinuum::NoOffset:     strm << "None</td></tr>\n";           break;
    case PeakContinuum::Constant:     strm << "Constant</td></tr>\n";       break;
    case PeakContinuum::Linear:       strm << "Linear</td></tr>\n";         break;
    case PeakContinuum::Quadratic:    strm << "Quadratic</td></tr>\n";      break;
    case PeakContinuum::Cubic:        strm << "Cubic</td></tr>\n";          break;
    case PeakContinuum::FlatStep:     strm << "Flat Step</td></tr>\n";      break;
    case PeakContinuum::LinearStep:   strm << "Linear Step</td></tr>\n";    break;
    case PeakContinuum::BiLinearStep: strm << "Bi-linear Step</td></tr>\n"; break;
    case PeakContinuum::External:     strm << "Global</td></tr>\n";         break;
  }//switch( peak.continuum()->type() )
  
  double offset_area = 0.0;
  
  try
  {
    peak.continuum()->offset_integral(peak.lowerX(), peak.upperX(), data);
  }catch(...)
  {
    //can only get here for PeakContinuum::FlatStep/LinearStep/BiLinearStep with invalid data
  }
  
  strm << startrow << "Cont. Area" << starttd << offset_area << "</td></tr>\n";
  if( peak.chi2Defined() )
    strm << startrow << "Chi2/DOF" << starttd << peak.chi2dof() << "</td></tr>\n";
  strm << "</table>";
  
  return strm.str();
}//string SpectrumViewerTester::makePeakSummarryTable( const PeakDef &peak )


struct PeakLessThan : binary_function <PeakDef,PeakDef,bool>
{
  bool operator() (const PeakDef& x, const PeakDef& y) const
  {
    return x.mean()<y.mean();
  }
};


SpectrumViewerTester::Score SpectrumViewerTester::testManualPeakClicking()
{
  Score score;
  score.m_test = ManualPeakClicking;
  PeakModel *peakModel = m_viewer->m_peakModel;
  
  const auto data = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  assert( data );
  
  if( !peakModel->peaks() || peakModel->peaks()->empty() )
    return score;
  
  vector<PeakDef> truthpeaks;
  for( const PeakModel::PeakShrdPtr &peak : *peakModel->peaks() )
    truthpeaks.push_back( *peak );
  
  peakModel->setPeaks( truthpeaks );
  const size_t numpeaks = peakModel->npeaks();
  
  for( size_t peakn = 0; peakn < numpeaks; ++peakn )
  {
    PeakModel::PeakShrdPtr peak, nextPeak, prevPeak;
    try
    {
      peak = peakModel->peakPtr( peakn );
      if( peakn > 0 )
        prevPeak = peakModel->peakPtr( peakn-1 );
      if( peakn < (numpeaks-1) )
        nextPeak = peakModel->peakPtr( peakn+1 );
    }catch( std::exception &e )
    {
      cerr << "SpectrumViewerTester::testManualPeakClicking() caught: "
           << e.what() << endl;
      continue;
    }//try / catch
    
    int nsuccess = 0, nfailures = 0, neffectOtherPeaks = 0, nnofit = 0, nnotsimilar = 0;
    
    int meanNotCompat = 0, ampNotCompat = 0,
        lowerRoiNotCompat = 0, upperRoiNotCompat = 0,
        widthsNotCompat = 0;
    
    // TODO: make sure the ROI is at least ~45px (a guess) wide so for tiny peaks the
    //       uncertainty applied to where the user clicked isnt too huge.
    const float wantedRoiWidthPx = 45;
    const float spectrumLowestEnergy = data->gamma_energy_min();
    const float spectrumHighestEnergy = data->gamma_energy_max();
    const double chartWidthPx = std::max( m_viewer->m_spectrum->chartWidthInPixels(), 640.0 );
    const double roiWidth = peak->upperX() - peak->lowerX();
    const double targetKevPerPx = roiWidth / wantedRoiWidthPx;
    const double targetEnergyRange = chartWidthPx * targetKevPerPx;
    const double targetDisplayLowerEnergy = peak->mean() - 0.5*targetEnergyRange;
    const double targetDisplayUpperEnergy = peak->mean() + 0.5*targetEnergyRange;
    cout << "targetDisplayLowerEnergy=" << targetDisplayLowerEnergy << ", targetDisplayUpperEnergy=" << targetDisplayUpperEnergy << ", peak->mean()=" << peak->mean() << ", chartWidthPx=" << chartWidthPx << ", roiWidth=" << roiWidth << endl;
    
    m_viewer->setDisplayedEnergyRange( targetDisplayLowerEnergy, targetDisplayUpperEnergy );
    
    for( int attempt = 0; attempt < sm_num_manual_click; ++attempt )
    {
      peakModel->setPeaks( truthpeaks );
      peakModel->removePeak( peakn );
      
      const double frac = double(attempt) / sm_num_manual_click;
      const double x = peak->mean() + (2.0*frac-1.0)*sm_nfwhm_manually_click*peak->fwhm();
      
      const size_t nprepeaks = peakModel->npeaks();
      m_viewer->searchForSinglePeak( x );
      
      const size_t npostpeaks = peakModel->npeaks();
      if( nprepeaks == npostpeaks )
      {
        ++nnofit;
        ++nfailures;
      }else
      {
        PeakModel::PeakShrdPtr nearestLook, nextLook, prevLook;

        nearestLook = peakModel->nearestPeak( peak->mean() );
        if( nextPeak )
          nextLook = peakModel->nearestPeak( nextPeak->mean() );
        if( prevPeak )
          prevLook = peakModel->nearestPeak( prevPeak->mean() );
        
        if( nearestLook == nextLook || nearestLook == prevLook
           || !nearestLook
           || (nextPeak && !nextLook)
           || (prevPeak && !prevLook) )
        {
          ++nfailures;
        }else if( !compatibleMeans( &(*peak), &(*nearestLook) )
                  || !compatibleAmplitudes( &(*peak), &(*nearestLook) )
                  || !compatibleWidths( &(*peak), &(*nearestLook) ) )
        {
          ++nfailures;
          ++nnotsimilar;
          
          meanNotCompat     += !compatibleMeans( &(*peak), &(*nearestLook) );
          ampNotCompat      += !compatibleAmplitudes( &(*peak), &(*nearestLook) );
          lowerRoiNotCompat += !compatibleLowerRoi( &(*peak), &(*nearestLook) );
          upperRoiNotCompat += !compatibleUpperRoi( &(*peak), &(*nearestLook) );
          widthsNotCompat   += !compatibleWidths( &(*peak), &(*nearestLook) );
          
          //Could add in the average percentage off here
          
        }else if( (nextPeak && !peaksAreSimilar( &(*nextPeak), &(*nextLook) ))
                  || (prevPeak && !peaksAreSimilar( &(*prevPeak), &(*prevLook) )) )
        {
          ++neffectOtherPeaks;
        }else
        {
          ++nsuccess;
        }
      }//if( nprepeaks == npostpeaks ) / else
    }//for( int attempt = 0; attempt < sm_num_manual_click; ++attempt )
    
    
    const double successfrac = double(nsuccess) / double(sm_num_manual_click);
    const double failfrac = double(nfailures) / double(sm_num_manual_click);
    const double nOtherEffect = double(neffectOtherPeaks) / double(sm_num_manual_click);
    const double noFitFrac = double(nnofit) / double(sm_num_manual_click);
    const double notsimilarFrac = double(nnotsimilar) / double(sm_num_manual_click);
    
    
    score.m_ntest  += 1.0;
    score.m_nwrong += (failfrac + nOtherEffect);
    score.m_ncorrect += successfrac;
    
    if( successfrac < sm_min_manually_click_corr_frac )
    {
      ScoreNote note;
    
      char onstack[256];
      if( neffectOtherPeaks )
      {
        snprintf( onstack, sizeof(onstack),
                  "Peak at %.1f keV fit in tolerances %i of %i tries, "
                  "effected other peaks %i times",
                  peak->mean(),
                  nsuccess, sm_num_manual_click, neffectOtherPeaks );
      }else
      {
        snprintf( onstack, sizeof(onstack),
                 "Peak at %.1f keV fit in tolerances %i of %i tries",
                 peak->mean(),
                 nsuccess, sm_num_manual_click );
      }//if( neffectOtherPeaks ) / else
      
      note.m_title = onstack;
      note.m_text = makePeakSummarryTable( *peak, data );
      
      double lowx = peak->mean() - sm_nfwhm_manually_click*peak->fwhm();
      double highx = peak->mean() + sm_nfwhm_manually_click*peak->fwhm();
      
      stringstream msg;
      msg << "<div>Simulated clicking " << sm_num_manual_click
          << " times ranging from " << std::setprecision(1) << std::fixed
          << lowx << " to " << highx << " keV</div>\n"
          << "<div>Fit within tolerances " << nsuccess << " times; "
          << 100.0*successfrac << "% of the time</div>\n";
      if( nnotsimilar )
      {
         msg << "<div>Fit, but not in tolerances " << nnotsimilar << " times; "
             << 100.0*notsimilarFrac << "% of the time"
              << "<div style=\"margin-left: 20px;\">Caused by mean being off: "
              << 100.0*double(meanNotCompat)/nnotsimilar << "% of these times</div>"
              << "<div style=\"margin-left: 20px;\">Caused by width being off: "
              << 100.0*double(widthsNotCompat)/nnotsimilar << "% of these times</div>"
              << "<div style=\"margin-left: 20px;\">Caused by amplitude being off: "
              << 100.0*double(ampNotCompat)/nnotsimilar << "% of these times</div>"
//              << "<div style=\"margin-left: 20px;\">Caused by ROI start being off: "
//              << 100.0*double(lowerRoiNotCompat)/nnotsimilar << "% of these times</div>"
//              << "<div style=\"margin-left: 20px;\">Caused by ROI end off: "
//              << 100.0*double(upperRoiNotCompat)/nnotsimilar << "% of these times</div>"
            <<"</div>\n";
      }//if( nnotsimilar )
      
      if( nnofit )
        msg << "<div>Failed to fit an additional peak " << nnofit << " times; "
            << 100.0*noFitFrac << "% of the time</div>\n";
      
      if( neffectOtherPeaks )
        msg << "<div>Effected the neighboring peaks " << neffectOtherPeaks << " times; "
            << 100.0*nOtherEffect << "% of the time</div>\n";
      
      msg << "<div>Failed (including above) " << nfailures << " times; "
          << 100.0*failfrac << "% of the time</div>\n";
      
      note.m_text += msg.str();
      
      const double diff = std::max( 5.0, highx - lowx );
      lowx -= 2.0*diff;
      highx += 2.0*diff;
    
//      peakModel->setPeaks( truthpeaks );
//      peakModel->removePeak( peakn );
//      note.m_testImage = renderChartToSvg( lowx, highx, m_picWidth, m_picHeight );
    
      peakModel->setPeaks( truthpeaks );
      note.m_originalImage = renderChartToSvg( lowx, highx, m_picWidth, m_picHeight );
    
      score.m_notes.push_back( note );
    }//if( fail )
  }//for( size_t peakn = 0; peakn < numpeaks; ++peakn )
  
  peakModel->setPeaks( truthpeaks );
  
  return score;
}//Score testManualPeakClicking()


SpectrumViewerTester::Score SpectrumViewerTester::testAutomatedPeakSearch()
{
  Score score;
  score.m_test = AutomatedPeakSearch;
  PeakModel *peakModel = m_viewer->m_peakModel;
  
  if( !peakModel->peaks() || peakModel->peaks()->empty() )
    return score;
  
  vector<PeakDef> truthpeaks, testpeaks;
  for( const PeakModel::PeakShrdPtr &peak : *peakModel->peaks() )
    truthpeaks.push_back( *peak );

  peakModel->setPeaks( std::vector<PeakDef>() );
  
  shared_ptr<const SpecUtils::Measurement> data = m_viewer->m_spectrum->data();
  try
  {
    //This will take a minute, and app will look paused, but whatever.
    auto resultpeaks = ExperimentalAutomatedPeakSearch::search_for_peaks( data, nullptr, nullptr, false );
    
    for( const auto &p : resultpeaks )
      testpeaks.push_back( *p );
  }catch( std::exception &e )
  {
    string msg = "InterSpec::testAutomatedPeakSearch(): caught exception: '";
    msg += e.what();
    msg += "'";
    cerr << msg << endl;
  }//try / catch

  
  for( const PeakModel::PeakShrdPtr &peak : *peakModel->peaks() )
    testpeaks.push_back( *peak );
  
  //Lets go through and match the test peaks to the truth peaks, by mean.
  // -Each truth peak will have a list of 0 or more test peaks associated w/ it.
  //  'truthToTest'
  // -Test peaks not assigned to truth peak will be put in seperate list.
  //  'testWithNoTruth'
  // -Truth peaks with no test peaks placed in seperate list.
  //  'truthWithNoTest'
  // -After assigning test to truth peaks, will then go through and remove
  //  pairs where there is a 1:1 correspondence between truth and test, and
  //  all aspects of the peaks are within tolerances.
  //  'truthToTest' and 'correctPeaks'.
  // -For ROIs with multiple truth peaks, obtain peaks with means in that range,
  //  and check to see if there is an equal number of test peaks that share a
  //  ROI; eliminate peaks for both these ROI for next step.
  // -For ROIs with mutliple test peaks, report this
  //
  //The reporting then will go as:
  // -List truth peaks with no associated test peaks
  //  'truthWithNoTest'
  // -List test peaks with no cooresponding truth peak.
  //  'testWithNoTruth'
  // -List truth peaks with mutliple test peaks
  //  'truthToTest'
  // -List truth peaks with 1 test peak, but test peak outside of tolerances
  //  'truthToTest'
  // -List truth ROIs with incorrect test peaks
  //  'truthContinuumToTest'
  // -List test ROIs that dont have cooresponding truth ROIs
  //  'testRoiWithNoTruthRoi'

  
  typedef map< PeakDef, vector<PeakDef>, PeakLessThan > TruthToTestMap;
  TruthToTestMap truthToTest;
  vector<PeakDef> testWithNoTruth, truthWithNoTest, correctPeaks;
  
  for( const PeakDef &peak : truthpeaks )
    truthToTest[peak];
  for( size_t i = 0; i < testpeaks.size(); ++i )
  {
    const PeakDef &test = testpeaks[i];
    double mindiff = 999999.9;
    size_t nearestInd = truthpeaks.size();
    for( size_t j = 0; j < truthpeaks.size(); ++j )
    {
      const PeakDef &truth = truthpeaks[j];
      const double diff = fabs(test.mean() - truth.mean()) / truth.sigma();
      if( diff < mindiff )
      {
        mindiff = diff;
        nearestInd = j;
      }
    }//for( size_t j = 0; j < truthpeaks.size(); ++j )
    
    if( mindiff < sm_mean_max_nsigma_diff )
      truthToTest[ truthpeaks[nearestInd] ].push_back( test );
    else
      testWithNoTruth.push_back( test );
  }//for( size_t i = 0; i < testpeaks.size(); ++i )
  
  for( const TruthToTestMap::value_type &vt : truthToTest )
  {
    const PeakDef &truth = vt.first;
    const vector<PeakDef> &tests = vt.second;
    
    if( tests.empty() )
    {
      truthWithNoTest.push_back( truth );
    }else if( tests.size() == 1 )
    {
      const PeakDef &test = tests[0];
      if( peaksAreSimilar( &truth, &test ) )
        correctPeaks.push_back( truth );
      //else - peak does not match
    }else
    {
      //nothing to do here - we have multiple test peaks for this truth peak
    }
  }//for( const TruthToTestMap::value_type &vt : truthToTest )
  
  typedef map<std::shared_ptr<const PeakContinuum>, vector<PeakDef> > ContinuumToPeaks;
  ContinuumToPeaks truthContinuumToPeaks, testContinuumToPeaks;
  
  for( const PeakDef &peak : truthpeaks )
    truthContinuumToPeaks[peak.continuum()].push_back( peak );
  for( const PeakDef &peak : testpeaks )
    testContinuumToPeaks[peak.continuum()].push_back( peak );

  //Remove the (truth) peaks from 'truthContinuumToPeaks' for the truth peaks
  //  that dont have multiple peaks for a ROI
  vector<std::shared_ptr<const PeakContinuum> > toRemove;
  for( const ContinuumToPeaks::value_type &t : truthContinuumToPeaks )
  {
    if( t.second.size() < 2 )
      toRemove.push_back( t.first );
  }//for( each truthContinuumToPeaks )
  for( size_t i = 0; i < toRemove.size(); ++i )
    truthContinuumToPeaks.erase( toRemove[i] );
  toRemove.clear();
  
  //Remove the (truth) peaks from 'testContinuumToPeaks' for the truth peaks
  //  that dont have multiple peaks for a ROI
  for( const ContinuumToPeaks::value_type &t : testContinuumToPeaks )
  {
    if( t.second.size() < 2 )
      toRemove.push_back( t.first );
  }//for( each testContinuumToPeaks )
  for( size_t i = 0; i < toRemove.size(); ++i )
    testContinuumToPeaks.erase( toRemove[i] );
  toRemove.clear();
  
  //Remove the 1:1 Truth:Test peaks from 'truthToTest' (no reporting necassary)
  for( const PeakDef &truth : correctPeaks )
    truthToTest.erase( truth );
  
  //Remove the truth peaks that have no test peaks from 'truthToTest' (will use
  //  'truthWithNoTest' to report these)
  for( const PeakDef &truth : truthWithNoTest )
    truthToTest.erase( truth );
  
  //Now go through and try to match up ROIs between the truth and test
  typedef map< std::shared_ptr<const PeakContinuum>, vector<std::shared_ptr<const PeakContinuum> > > ContinuumToCandMap;
  ContinuumToCandMap truthContinuumToTest;
  set<std::shared_ptr<const PeakContinuum> > testRoiWithNoTruthRoi;
  
  for( ContinuumToPeaks::value_type &testToPeaks : testContinuumToPeaks )
    testRoiWithNoTruthRoi.insert( testToPeaks.first );
  
  for( ContinuumToPeaks::value_type &truthToPeaks : truthContinuumToPeaks )
  {
    std::shared_ptr<const PeakContinuum> truth = truthToPeaks.first;
    if( !truth->energyRangeDefined() )
    {
      cerr << "\n\n\n\nTruth value ROI with multiple peaks doesnt have energy"
           << " range defined!\n\n" << endl;
      continue;
    }//if( !truth->energyRangeDefined() )
    
    const double lowerx = truth->lowerEnergy();
    const double upperx = truth->upperEnergy();
    
    for( ContinuumToPeaks::value_type &testToPeaks : testContinuumToPeaks )
    {
      std::shared_ptr<const PeakContinuum> test = testToPeaks.first;
      if( !test->energyRangeDefined() )
        continue;
      
      //Check if ranges overlap, if so, its a candidate continuum
      if( max(lowerx,test->lowerEnergy()) <= min(upperx,test->upperEnergy()) )
      {
        truthContinuumToTest[truth].push_back( test );
        if( testRoiWithNoTruthRoi.count(test) )
          testRoiWithNoTruthRoi.erase(test);
      }//if( ranges overlap )
    }//for( ContinuumToPeaks::value_type &truthToPeaks : testContinuumToPeaks )
  }//for( ContinuumToPeaks::value_type &truthToPeaks : truthContinuumToPeaks )
  
  vector<std::shared_ptr<const PeakContinuum> > correctNumTestPeaksInRoi;
  for( ContinuumToCandMap::value_type &vt : truthContinuumToTest )
  {
    std::shared_ptr<const PeakContinuum> truth = vt.first;
    vector<std::shared_ptr<const PeakContinuum> > &testConts = vt.second;
    
    for( size_t i = 0; i < testConts.size(); ++i )
    {
      const vector<PeakDef> &truths = truthContinuumToPeaks[truth];
      const vector<PeakDef> &tests = testContinuumToPeaks[testConts[i]];
    
      if( truths.size() == tests.size() )
        correctNumTestPeaksInRoi.push_back( truth );
    }//for( size_t i = 0; i < testConts.size(); ++i )
  }//for( ContinuumToCandMap::value_type &vt : truthContinuumToTest )
 
  for( size_t i = 0; i < correctNumTestPeaksInRoi.size(); ++i )
    truthContinuumToTest.erase( correctNumTestPeaksInRoi[i] );
  

  score.m_ntest  += static_cast<double>( truthpeaks.size() );
  score.m_nwrong += static_cast<double>( truthWithNoTest.size() );
  score.m_nwrong += static_cast<double>( testWithNoTruth.size() );
  score.m_nwrong += static_cast<double>( truthToTest.size() );
  score.m_nwrong += static_cast<double>( truthContinuumToTest.size() );
  score.m_nwrong += static_cast<double>( testRoiWithNoTruthRoi.size() );
  score.m_ncorrect += static_cast<double>( truthpeaks.size()
                                            - truthWithNoTest.size() );


  //-List truth peaks with no associated test peaks
  for( const PeakDef &peak : truthWithNoTest )
  {
    ScoreNote note;
    
    char onstack[256];
    snprintf( onstack, sizeof(onstack),
             "Peak not found at %.1f keV", peak.mean() );
    note.m_title = onstack;
    
    
    WStringStream peakTableStrm;
    
    note.m_text = makePeakSummarryTable( peak, data );
    
    const double lowx = 2.0*peak.lowerX() - peak.upperX();
    const double upperx = 2.0*peak.upperX() - peak.lowerX();
    
    peakModel->setPeaks( testpeaks );
    note.m_testImage = renderChartToSvg( lowx, upperx, m_picWidth, m_picHeight );
    
    peakModel->setPeaks( truthpeaks );
    note.m_originalImage = renderChartToSvg( lowx, upperx, m_picWidth, m_picHeight );
    
    score.m_notes.push_back( note );
  }//for( ContinuumToPeaks::value_type &vt : continuumToPeaks )
  
  //-List test peaks with no cooresponding truth peak.
  for( const PeakDef &test : testWithNoTruth )
  {
    ScoreNote note;
    
    char onstack[256];
    snprintf( onstack, sizeof(onstack),
             "Found unexpected peak at %.1f keV", test.mean() );
    note.m_title = onstack;
    note.m_text = makePeakSummarryTable( test, data );
    
    const double lowx = 2.0*test.lowerX() - test.upperX();
    const double upperx = 2.0*test.upperX() - test.lowerX();
    
    peakModel->setPeaks( testpeaks );
    note.m_testImage = renderChartToSvg( lowx, upperx, m_picWidth, m_picHeight );
    
    peakModel->setPeaks( truthpeaks );
    note.m_originalImage = renderChartToSvg( lowx, upperx, m_picWidth, m_picHeight );
    
    score.m_notes.push_back( note );
  }//for( ContinuumToPeaks::value_type &vt : continuumToPeaks )
  
  
  // -List truth peaks with mutliple test peaks
  for( const TruthToTestMap::value_type &vt : truthToTest )
  {
    const PeakDef &truth = vt.first;
    const vector<PeakDef> &tests = vt.second;
    if( tests.size() < 2 )
      continue;
    
    ScoreNote note;
    
    char onstack[256];
    snprintf( onstack, sizeof(onstack),
             "Found %i peaks where one was expected at %.1f keV",
             int(tests.size()), truth.mean() );
    note.m_title = onstack;
//    note.m_text = makePeakSummarryTable( truth, data );
    
    double lowxer = 99999.9, upperx = -999999.9;
    for( const PeakDef &p : tests )
    {
      lowxer = std::min(lowxer,p.lowerX());
      upperx = std::max(upperx,p.upperX());
    }//for( const PeakDef &p : tests )
    
    const double dx = upperx - lowxer;
    lowxer = lowxer - 1.0*dx;
    upperx = upperx + 1.0*dx;
    
    peakModel->setPeaks( testpeaks );
    note.m_testImage = renderChartToSvg( lowxer, upperx, m_picWidth, m_picHeight );
//    
    peakModel->setPeaks( truthpeaks );
    note.m_originalImage = renderChartToSvg( lowxer, upperx, m_picWidth, m_picHeight );
    
    score.m_notes.push_back( note );
  }//for( ContinuumToPeaks::value_type &vt : continuumToPeaks )

  // -List truth peaks with 1 test peak, but test peak outside of tolerances
  for( const TruthToTestMap::value_type &vt : truthToTest )
  {
    const PeakDef &truth = vt.first;
    const vector<PeakDef> &tests = vt.second;
    if( tests.size() != 1 )
      continue;
    
    const PeakDef &test = tests[0];
    
    ScoreNote note;
    
    char onstack[1024];
    snprintf( onstack, sizeof(onstack),
             "Found peak outside of tolerance at %.1f keV", truth.mean() );
    note.m_title = onstack;
    
  
    if( !compatibleMeans( &truth, &test ) )
    {
      snprintf( onstack, sizeof(onstack),
               "<div class=\"FailedConditionNote\">Means not compatible:"
               " expected %.1f, got %.1f</div>",
                truth.mean(), test.mean() );
      note.m_text += onstack;
    }//if( !compatibleMeans( &truth, &test ) )
    
    if( !compatibleAmplitudes( &truth, &test ) )
    {
      snprintf( onstack, sizeof(onstack),
               "<div class=\"FailedConditionNote\">Amplitudes not compatible:"
               " expected %.1f, got %.1f</div>",
               truth.amplitude(), test.amplitude() );
      note.m_text += onstack;
    }//if( !compatibleAmplitudes( &truth, &test ) )
    
    if( !compatibleLowerRoi( &truth, &test ) )
    {
      snprintf( onstack, sizeof(onstack),
               "<div class=\"FailedConditionNote\">ROI lower edge not compatible:"
               " expected %.1f, got %.1f</div>",
               truth.lowerX(), test.lowerX() );
      note.m_text += onstack;
    }//if( !compatibleLowerRoi( &truth, &test ) )
    
    if( !compatibleUpperRoi( &truth, &test ) )
    {
      snprintf( onstack, sizeof(onstack),
               "<div class=\"FailedConditionNote\">ROI upper edge not compatible:"
               " expected %.1f, got %.1f</div>",
               truth.upperX(), test.upperX() );
      note.m_text += onstack;
    }//if( !compatibleUpperRoi( &truth, &test ) )
    
    if( !compatibleWidths( &truth, &test ) )
    {
      snprintf( onstack, sizeof(onstack),
               "<div class=\"FailedConditionNote\">Peak width not compatible:"
               " expected %.1f, got %.1f</div>",
               truth.sigma(), test.sigma() );
      note.m_text += onstack;
    }//if( !compatibleWidths( &truth, &test ) )
    
    
    double lowxer = 99999.9, upperx = -999999.9;
    for( const PeakDef &p : tests )
    {
      lowxer = std::min(lowxer,p.lowerX());
      upperx = std::max(upperx,p.upperX());
    }//for( const PeakDef &p : tests )
    
    const double dx = upperx - lowxer;
    lowxer = lowxer - 1.0*dx;
    upperx = upperx + 1.0*dx;
    
    peakModel->setPeaks( testpeaks );
    note.m_testImage = renderChartToSvg( lowxer, upperx, m_picWidth, m_picHeight );
    
    peakModel->setPeaks( truthpeaks );
    note.m_originalImage = renderChartToSvg( lowxer, upperx, m_picWidth, m_picHeight );
    
    score.m_notes.push_back( note );
  }//for( ContinuumToPeaks::value_type &vt : continuumToPeaks )
  

  // -List truth ROIs with incorrect test peaks
  for( const ContinuumToCandMap::value_type &vt : truthContinuumToTest )
  {
    std::shared_ptr<const PeakContinuum> truthCont = vt.first;
    const vector<std::shared_ptr<const PeakContinuum> > &testConts = vt.second;
    const vector<PeakDef> &truths = truthContinuumToPeaks[truthCont];
    
    
    for( size_t i = 0; i < testConts.size(); ++i )
    {
      const vector<PeakDef> &tests = testContinuumToPeaks[testConts[i]];
     
      if( tests.size() == truths.size() )
        continue;
    
      ScoreNote note;
    
      char onstack[256];
      snprintf( onstack, sizeof(onstack),
               "Found %i peaks in ROI where %i is expected",
               int(tests.size()), int(truths.size()) );
      note.m_title = onstack;
      //    note.m_text = makePeakSummarryTable( truth, data );
    
      double lowxer = 99999.9, upperx = -999999.9;
      for( const PeakDef &p : tests )
      {
        lowxer = std::min(lowxer,p.lowerX());
        upperx = std::max(upperx,p.upperX());
      }//for( const PeakDef &p : tests )
      
      for( const PeakDef &p : truths )
      {
        lowxer = std::min(lowxer,p.lowerX());
        upperx = std::max(upperx,p.upperX());
      }//for( const PeakDef &p : tests )
    
      const double dx = upperx - lowxer;
      lowxer = lowxer - 1.0*dx;
      upperx = upperx + 1.0*dx;
    
      peakModel->setPeaks( testpeaks );
      note.m_testImage = renderChartToSvg( lowxer, upperx, m_picWidth, m_picHeight );
    
      peakModel->setPeaks( truthpeaks );
      note.m_originalImage = renderChartToSvg( lowxer, upperx, m_picWidth, m_picHeight );
    
      score.m_notes.push_back( note );
    }//for( size_t i = 0; i < testConts.size(); ++i )
  }//for( ContinuumToPeaks::value_type &vt : continuumToPeaks )
  
  
  // -List test ROIs that dont have cooresponding truth ROIs
  typedef std::shared_ptr<const PeakContinuum> ContPtr;
  for( const ContPtr &continuum : testRoiWithNoTruthRoi )
  {
    ScoreNote note;
    const vector<PeakDef> &tests = testContinuumToPeaks[continuum];
    
    char onstack[256];
    snprintf( onstack, sizeof(onstack),
             "Found unexpected ROI shared by %i peaks.", int(tests.size()) );
    note.m_title = onstack;
    
    snprintf( onstack, sizeof(onstack),
             "ROI lower energy: %.f, upper energy: %.1f",
             continuum->lowerEnergy(), continuum->upperEnergy() );
    note.m_text = onstack;
    
    const double dx = continuum->upperEnergy() - continuum->lowerEnergy();
    const double lowx = continuum->lowerEnergy() - dx;
    const double upperx = continuum->upperEnergy() + dx;
    
    peakModel->setPeaks( testpeaks );
    note.m_testImage = renderChartToSvg( lowx, upperx, m_picWidth, m_picHeight );
    
    peakModel->setPeaks( truthpeaks );
    note.m_originalImage = renderChartToSvg( lowx, upperx, m_picWidth, m_picHeight );
    
    score.m_notes.push_back( note );
  }//for( const ContPtr &continuum : testRoiWithNoTruthRoi )

    
  peakModel->setPeaks( truthpeaks );  //not necassarry
  
  return score;
}//Score testAutomatedPeakSearch()

SpectrumViewerTester::Score SpectrumViewerTester::testMultiplePeakFit()
{
  Score answer;
  answer.m_test = MultiplePeakInRoiFit;
  Score nominalscore = testMultiplePeakFitRangeVaried( 0.0, 0.0 );
  
  answer.addScore( nominalscore, sm_multipeakfit_nominal_weight );
  answer.m_notes.insert( answer.m_notes.end(),
                         nominalscore.m_notes.begin(),
                         nominalscore.m_notes.end() );
  
  const int numtest = 10;
  const double testweight = (1.0-sm_multipeakfit_nominal_weight)/numtest;
  for( int i = 0; i < numtest; ++i )
  {
    const double upperDelta = 0.0;
    const double lowerDelta = (i - 0.5*numtest + 1.0)
                            * sm_multipeakfit_range_max_vary_frac / (numtest-2);
    
    Score thisscore = testMultiplePeakFitRangeVaried( lowerDelta, upperDelta );
    answer.addScore( nominalscore, testweight );
    
    if( (thisscore.m_nwrong - nominalscore.m_nwrong) > 0.1 )
    {
      stringstream msg;
      msg << "The analysts range caused " << nominalscore.m_nwrong << " issues,"
          << " but varying the lower energy of ROI by " << (100.0*lowerDelta)
          << "% of ROI range, caused there to be " << thisscore.m_nwrong
          << " issues, contributing a weight of " << testweight*thisscore.m_nwrong
          << " towards this states total score.\n";
      if( thisscore.m_notes.size() )
        msg << "\t<ul>\n";
      for( size_t i = 0; i < thisscore.m_notes.size(); ++i )
        msg << "\t\t<li>" << thisscore.m_notes[i].m_text << "</li>\n";
      if( thisscore.m_notes.size() )
        msg << "\t</ul>\n";
      
      char title[512];
      snprintf( title, sizeof(title),
               "Varying lower energy by %.1f%% of ROI span caused issue",
               (100.0*lowerDelta) );
      
      ScoreNote note;
      note.m_title = title;
      note.m_text = msg.str();
      answer.m_notes.push_back( note );
    }//if( fabs(nominalscore.m_nwrong - thisscore.m_nwrong) > 0.1 )
  }//for( test variying lower range )
  
  for( int i = 0; i < numtest; ++i )
  {
    const double lowerDelta = 0.0;
    const double upperDelta = (i - 0.5*numtest + 1.0)
                             * sm_multipeakfit_range_max_vary_frac / (numtest-2);
    
    Score thisscore = testMultiplePeakFitRangeVaried( lowerDelta, upperDelta );
    answer.addScore( nominalscore, testweight );
    
    if( (thisscore.m_nwrong - nominalscore.m_nwrong) > 0.1 )
    {
      stringstream msg;
      msg << "The analysts range caused " << nominalscore.m_nwrong << " issues,"
      << " but varying the upper energy of ROI by " << (100.0*lowerDelta)
      << "% of ROI range, caused there to be " << thisscore.m_nwrong
      << " issues, contributing a weight of " << testweight*thisscore.m_nwrong
      << " towards this states total score.";
      
      if( thisscore.m_notes.size() )
        msg << "\t<ol>\n";
      for( size_t i = 0; i < thisscore.m_notes.size(); ++i )
        msg << "\t\t<li>" << thisscore.m_notes[i].m_text << "</li>\n";
      if( thisscore.m_notes.size() )
        msg << "\t</ol>\n";
      
      char title[512];
      snprintf( title, sizeof(title),
               "Varying upper energy by %.1f%% of ROI span caused issue",
               (100.0*lowerDelta) );
      
      ScoreNote note;
      note.m_title = title;
      note.m_text = msg.str();
      answer.m_notes.push_back( note );
    }//if( fabs(nominalscore.m_nwrong - thisscore.m_nwrong) > 0.1 )
  }//for( test variying upper range )
  
  return answer;
}//Score testMultiplePeakFit()


SpectrumViewerTester::Score SpectrumViewerTester::testShieldSourceFit()
{
  Score answer;
  answer.m_test = TestType::SourceShieldingFit;
  
  shared_ptr<SpecMeas> spec = m_viewer->measurment(SpecUtils::SpectrumType::Foreground);
  if( !spec )
  {
    ScoreNote note;
    note.m_title = "No foreground defined";
    note.m_text = "Couldnt perform fit";
    answer.m_notes.push_back( note );
    
    return answer;
  }//if( !spec )
  
  auto fitModel = spec->shieldingSourceModel();
  if( !fitModel )
  {
    ScoreNote note;
    note.m_title = "No model is defined for fitting";
    note.m_text = "Cant test";
    answer.m_notes.push_back( note );
    
    return answer;
  }//if( !fitModel )
  
  AuxWindow *window = nullptr;
  
  try
  {
    pair<ShieldingSourceDisplay *,AuxWindow *> widgets = ShieldingSourceDisplay::createWindow( m_viewer );
    window = get<1>(widgets);
    ShieldingSourceDisplay *disp = get<0>(widgets);
    
    std::tuple<int,int,bool> nfitpars = disp->numTruthValuesForFitValues();
    
    if( !std::get<2>(nfitpars) )
    {
      ScoreNote note;
      note.m_title = "Truth values not entered";
      note.m_text = "Not all required truth values or tolerances were entered.";
      answer.m_notes.push_back( note );
      
      return answer;
    }//if( we can actually test things )
    
    ScoreNote note;
    note.m_text = "";
    note.m_title = "Fit Chi2 Graphic";
    note.m_originalImage = make_shared<Wt::WSvgImage>(m_picWidth, m_picHeight);
    disp->updateChi2ChartActual();
    disp->renderChi2Chart( *note.m_originalImage );
    
    // Set all the quantities being fit for to some default values, so we wont be starting off in
    //  unfairly close to the real answer.
    disp->setFitQuantitiesToDefaultValues();
    
    shared_ptr<ShieldingSourceDisplay::ModelFitResults> result = disp->doModelFit( false );
    
    // Uhg, sometimes users have to hit  a couple times right now - this should be fixed
    result = disp->doModelFit( false );
    
    note.m_testImage = make_shared<Wt::WSvgImage>(m_picWidth, m_picHeight);
    disp->renderChi2Chart( *note.m_testImage );
    answer.m_notes.push_back( note );
    
    switch( result->succesful )
    {
      case ShieldingSourceDisplay::ModelFitResults::FitStatus::TimedOut:
      case ShieldingSourceDisplay::ModelFitResults::FitStatus::UserCancelled:
      case ShieldingSourceDisplay::ModelFitResults::FitStatus::InvalidOther:
      {
        answer.m_ntest += std::get<1>(nfitpars);
        answer.m_nwrong += std::get<1>(nfitpars);
        
        ScoreNote note;
        note.m_title = "Fit of model failed";
        note.m_text = "Failed to fit model, so saying all " + std::to_string(std::get<1>(nfitpars))
                      + " fit parameters were wrong.";
        
        if( !result->errormsgs.empty() )
          note.m_text += "<br /><b>Error Messages:</b>";
        
        for( const auto &msg : result->errormsgs )
          note.m_text += "<div>" + msg + "</div>";
        answer.m_notes.push_back( note );
        
        return answer;
        break;
      }//invalid fit result
        
      case ShieldingSourceDisplay::ModelFitResults::FitStatus::InterMediate:
        assert( 0 );
        break;
        
      case ShieldingSourceDisplay::ModelFitResults::FitStatus::Final:
        break;
    }//switch( result->succesful )
    
    const tuple<bool,int,int,vector<string>> testres = disp->testCurrentFitAgainstTruth();
    assert( get<0>(testres) );
    const int numPass = get<1>(testres);
    const int numTested = get<2>(testres);
    const vector<string> &remarks = get<3>(testres);
    
    answer.m_ncorrect += numPass;
    answer.m_ntest += numTested;
    answer.m_nwrong += (numTested - numPass);
    
    ScoreNote valuesNote;
    valuesNote.m_title = "Fit values";
    for( const auto &msg : remarks )
      valuesNote.m_text += "<div>" + msg + "</div>";
    
    answer.m_notes.push_back( valuesNote );
  }catch( std::exception &e )
  {
    ScoreNote note;
    note.m_title = "Failed test with exception";
    note.m_text = "Exception message: " + string( e.what() );
    answer.m_notes.push_back( note );
  }//try / catch
  
  AuxWindow::deleteAuxWindow( window );
  
  return answer;
}//Score testShieldSourceFit()



SpectrumViewerTester::Score SpectrumViewerTester::testMultiplePeakFitRangeVaried(
                            double lowerXChangeFrac, double upperXChangeFrac )
{
  Score score;
  score.m_test = MultiplePeakInRoiFit;
  
  PeakModel *peakModel = m_viewer->m_peakModel;
  
  const auto data = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  
  if( !peakModel->peaks() || peakModel->peaks()->empty() )
    return score;
  
  char onstack[512];
  vector<PeakDef> truthpeaks;
  typedef map< std::shared_ptr<const PeakContinuum>, vector<PeakModel::PeakShrdPtr> > ContinuumToPeakMap;
  ContinuumToPeakMap continuumToPeaks;
  
  for( const PeakModel::PeakShrdPtr &peak : *peakModel->peaks() )
  {
    truthpeaks.push_back( *peak );
    continuumToPeaks[peak->continuum()].push_back( peak );
  }//for( const PeakModel::PeakShrdPtr &peak : *peakModel->peaks() )
  

  for( ContinuumToPeakMap::value_type &vt : continuumToPeaks )
  {
    vector<PeakModel::PeakShrdPtr> peaks = vt.second;
    const int npeaks = static_cast<int>( peaks.size() );
    
    std::sort( peaks.begin(), peaks.end(), &PeakDef::lessThanByMeanShrdPtr );
    
    if( npeaks < 2 )
      continue;

    peakModel->setPeaks( truthpeaks );
    
    for( int peakn = 0; peakn < npeaks; ++peakn )
    {
      PeakModel::PeakShrdPtr peak = peaks[peakn];
      
      const int nrows = peakModel->rowCount();
      for( int i = 0; i < nrows; ++i )
      {
        Wt::WModelIndex ind = peakModel->index( i, 0 );
        PeakModel::PeakShrdPtr thispeak = peakModel->peak( ind );
        if( thispeak && ((thispeak->mean()-peak->mean())<0.1) )
        {
          peakModel->removePeak( ind );
          break;
        }//if( thispeak && ((thispeak->mean()-peak->mean())<0.1) )
      }//for( int i = 0; i < nrows; ++i )
    }//for( int peakn = 0; peakn < npeaks; ++peakn )

    if( (truthpeaks.size()-npeaks) != peakModel->peaks()->size() )
      throw runtime_error( "SpectrumViewerTester::testMultiplePeakFit: "
                           "Failed to remove all of the input peaks, or "
                           "something" );
    
    vector<PeakModel::PeakShrdPtr> preFitPeaks;
    for( const PeakModel::PeakShrdPtr &peak : *peakModel->peaks() )
      preFitPeaks.push_back( peak );
    std::sort( preFitPeaks.begin(), preFitPeaks.end(), &PeakDef::lessThanByMeanShrdPtr );
    
    double lowerx = vt.first->lowerEnergy();
    double upperx = vt.first->upperEnergy();
    const double energydelta = upperx - lowerx;
    lowerx = lowerx + lowerXChangeFrac*energydelta;
    upperx = upperx + upperXChangeFrac*energydelta;
    
    D3SpectrumDisplayDiv *specchart = m_viewer->spectrum();
    specchart->performDragCreateRoiWork( lowerx, upperx, npeaks, true, -1.0, -1.0 );
    
    vector<PeakDef> foundpeaks, refitOrigPeaks, allPostFitPeaks;
    for( const PeakModel::PeakShrdPtr &peak : *peakModel->peaks() )
    {
      allPostFitPeaks.push_back( *peak );
      if( !std::count( preFitPeaks.begin(), preFitPeaks.end(), peak) )
        foundpeaks.push_back( *peak );
      else
        refitOrigPeaks.push_back( *peak );
    }//for( const PeakModel::PeakShrdPtr &peak : *peakModel->peaks() )
    
    std::sort( foundpeaks.begin(), foundpeaks.end(), &PeakDef::lessThanByMean );
    std::sort( refitOrigPeaks.begin(), refitOrigPeaks.end(), &PeakDef::lessThanByMean );
    
    score.m_ntest += 1.0;
    
    if( foundpeaks.size() != peaks.size() )
    {
      score.m_nwrong += 1.0;
      ScoreNote note;
      
      snprintf( onstack, sizeof(onstack),
                "Failed to fit %i Peaks between %.1f and %.1f keV",
                int(npeaks), lowerx, upperx );
      note.m_title = onstack;
      
      snprintf( onstack, sizeof(onstack),
               "Fit for %i peaks instead of %i between %.1f and %.1f keV. "
               "The final chi2/dof was %.3f, compared to the expected %.3f.",
                int(foundpeaks.size()),
                int(npeaks), lowerx, upperx,
                (foundpeaks.size() ? foundpeaks[0].chi2dof() : -1.0),
                 peaks[0]->chi2dof()
               );
      
      note.m_text = onstack;
      
      const double x0 = lowerx - (upperx-lowerx);
      const double x1 = lowerx + (upperx-lowerx);
      note.m_testImage = renderChartToSvg( x0, x1, m_picWidth, m_picHeight );
      
      peakModel->setPeaks( truthpeaks );
      note.m_originalImage = renderChartToSvg( x0, x1, m_picWidth, m_picHeight );
      
      score.m_notes.push_back( note );
      
      continue;
    }//if( foundpeaks.size() != peaks.size() )
    
    //Check chi2 of region
    ScoreNote note;
    
    //Check mean, area, and width of each peak
    bool allPeaksMatch = true;
    for( size_t i = 0; i < peaks.size(); ++i )
    {
      const PeakDef *orig = peaks[i].get();
      const PeakDef *found = &foundpeaks[i];
      
      string issue;
      if( !compatibleMeans( orig, found ) )
        issue += "Mean Incompatible";
      if( !compatibleWidths( orig, found ) )
        issue = issue + (!issue.empty()?", ": "") + string("Width Incompatible");
      if( !compatibleAmplitudes( orig, found ) )
        issue = issue + (!issue.empty()?", ": "") + string("Amplitude Incompatible");
      if( !sameNuclideId( orig, found ) )
        issue = issue + (!issue.empty()?", ": "") + string("Nuclide ID Changed");
      
      if( !issue.empty() )
      {
        allPeaksMatch = false;
        
        stringstream msg;
        msg << "\t<div style=\"display: inline-block;\">\n"
            << "\t\t<div>Peak " << i << " at " << orig->mean() << " keV has: "
            << issue << "</div>"
            << "\t\t<div style=\"display: inline-block; padding-top:10px;\">"
               "<div style=\"color:green;\"><b>Expected</b></div>\n\t\t"
            << makePeakSummarryTable(*orig, data) << "</div>"
            << "\t\t<div style=\"display: inline-block; padding-top:10px;\">"
               "<div style=\"color:red;\"><b>Found</b></div>\n\t\t"
            << makePeakSummarryTable(*found, data) << "</div>\n";
        msg << "\t</div>\n";
        
        note.m_text += msg.str();
      }//if( !issue.empty() )
    }//for( size_t i = 0; i < peaks.size(); ++i )
    
    if( !allPeaksMatch )
    {
      score.m_nwrong += 1.0;
      snprintf( onstack, sizeof(onstack),
               "Peaks fit between %.1f and %.1f keV are a bit off",
               lowerx, upperx );
      note.m_title = onstack;
      
      snprintf( onstack, sizeof(onstack),
               "\t<div>The final chi2/dof was %.3f, compared to the expected %.3f.</div>\n",
               (foundpeaks.size() ? foundpeaks[0].chi2dof() : -1.0),
                peaks[0]->chi2dof()
                );
      note.m_text += onstack;
      
      const double x0 = lowerx - (upperx-lowerx);
      const double x1 = lowerx + (upperx-lowerx);
      
      peakModel->setPeaks( allPostFitPeaks );
      note.m_testImage = renderChartToSvg( x0, x1, m_picWidth, m_picHeight );
      
      peakModel->setPeaks( truthpeaks );
      note.m_originalImage = renderChartToSvg( x0, x1, m_picWidth, m_picHeight );
      
      score.m_notes.push_back( note );
      
      continue;
    }//if( !allPeaksMatch )
  
    if( preFitPeaks.size() != refitOrigPeaks.size() )
    {
      score.m_nwrong += 1.0;
      snprintf( onstack, sizeof(onstack),
               "Fitting for %i peaks between %.1f and %.1f changed number of "
               "other peaks", npeaks, lowerx, upperx );
      note.m_title = onstack;
      
      stringstream msg;
      msg << "<div>There were " << preFitPeaks.size() << " outside of the"
          << " region being fit for, however after the fit there where "
          << refitOrigPeaks.size() << "</div>" << endl;
      
      double minx = 9999999999.9, maxx = -9999999999.9;
      for( size_t i = 0; i < preFitPeaks.size(); ++i )
      {
        bool found = false;
        PeakModel::PeakShrdPtr pre = preFitPeaks[i];
        for( size_t j = 0; j < refitOrigPeaks.size(); ++j )
        {
          PeakDef post = refitOrigPeaks[j];
          if( fabs(post.mean() - pre->mean() ) < 0.1 )
            found = true;
        }//for( size_t j = 0; j < refitOrigPeaks.size(); ++j )
        
        if( !found )
        {
          minx = std::min( minx, pre->lowerX() );
          maxx = std::max( maxx, pre->upperX() );
          msg << "<div>Original peak at " << pre->mean() << " keV disapeared.</div>";
        }//if( !found )
      }//for( size_t i = 0; i < preFitPeaks.size(); ++i )
      
      for( size_t j = 0; j < refitOrigPeaks.size(); ++j )
      {
        bool found = false;
        PeakDef post = refitOrigPeaks[j];
        
        for( size_t i = 0; i < preFitPeaks.size(); ++i )
        {
          PeakModel::PeakShrdPtr pre = preFitPeaks[i];
          if( fabs(post.mean() - pre->mean() ) < 0.1 )
            found = true;
        }//for( size_t j = 0; j < refitOrigPeaks.size(); ++j )
        
        if( !found )
        {
          minx = std::min( minx, post.lowerX() );
          maxx = std::max( maxx, post.upperX() );
          msg << "<div>New peak appeared at " << post.mean() << " keV.</div>";
        }//if( !found )
      }//for( size_t i = 0; i < preFitPeaks.size(); ++i )
      
      SpectrumDataModel *model = m_viewer->m_spectrum->m_model;
      auto axisH = model->histUsedForXAxis();

      minx -= 20.0;
      maxx += 20.0;
      
      if( axisH )
      {
        const double chartmin = axisH->gamma_energy_min();
        const double chartmax = axisH->gamma_energy_max();
        
        minx = std::max( minx, chartmin );
        maxx = std::min( maxx, chartmax );
      }//if( axisH )
      
      peakModel->setPeaks( allPostFitPeaks );
      note.m_testImage = renderChartToSvg( minx, maxx, m_picWidth, m_picHeight );
      
      peakModel->setPeaks( truthpeaks );
      note.m_originalImage = renderChartToSvg( minx, maxx, m_picWidth, m_picHeight );
      
      note.m_text += msg.str();
      score.m_notes.push_back( note );
      
      continue;
    }//if( preFitPeaks.size() != refitOrigPeaks.size() )
  
//    preFitPeaks.size() != refitOrigPeaks.size() )
    int npeaksAltered = 0;
    for( size_t i = 0; i < preFitPeaks.size(); ++i )
    {
      PeakModel::PeakShrdPtr pre = preFitPeaks[i];
      PeakDef post = refitOrigPeaks[i];
      
      const bool compatible = compatibleMeans(pre.get(), &post)
                              && compatibleAmplitudes(pre.get(), &post)
                              && compatibleWidths(pre.get(), &post);
      
      if( !compatible )
      {
        ++npeaksAltered;
        ScoreNote thisnote;
        
        snprintf( onstack, sizeof(onstack),
                 "Peak at %f keV altered when fitting for %i peaks between"
                 " %.1f and %.1f keV", pre->mean(), npeaks, lowerx, upperx );
        thisnote.m_title = onstack;
        
        stringstream msg;
        msg << "\t<div>\n"
            << "\t\t<div style=\"display: inline-block;\">\n"
               "\t\t\t<div style=\"color:green;\"><b>Expected</b></div>\n\t\t\t"
            << makePeakSummarryTable(*pre, data) << "\n\t\t</div>"
            << "\t\t<div style=\"display: inline-block;\">\n"
               "\t\t\t<div style=\"color:red;\"><b>Found</b></div>\n\t\t\t"
            << makePeakSummarryTable(post, data) << "\n\t\t</div>\n"
            << "\t</div>\n";
        thisnote.m_text += msg.str();
        
        const double lowerx = std::min( pre->lowerX(), post.lowerX() );
        const double upperx = std::max( pre->upperX(), post.upperX() );
        const double minx = lowerx - (upperx - lowerx);
        const double maxx = upperx + (upperx - lowerx);
        
        peakModel->setPeaks( allPostFitPeaks );
        thisnote.m_testImage = renderChartToSvg( minx, maxx, m_picWidth, m_picHeight );
        
        peakModel->setPeaks( truthpeaks );
        thisnote.m_originalImage = renderChartToSvg( minx, maxx, m_picWidth, m_picHeight );
        
        score.m_notes.push_back( thisnote );
        score.m_nwrong += 1.0;
      }//if( peaksAreSimilar( pre.get(), &post ) )
    }//for( size_t i = 0; i < preFitPeaks.size(); ++i )
    
    if( npeaksAltered )
    {
      continue;
    }//if( npeaksAltered )
    
    score.m_ncorrect += 1.0;
  }//for( ContinuumToPeakMap::value_type &vt : continuumToPeaks )

  
  return score;
}//Score testMultiplePeakFit()



