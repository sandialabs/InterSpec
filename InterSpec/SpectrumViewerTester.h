#ifndef SpectrumViewerTester_h
#define SpectrumViewerTester_h
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
#include <vector>
#include <ostream>

#include <Wt/WContainerWidget>

#include "InterSpec/AuxWindow.h"

class PeakDef;
struct UserState;
class SideMenuItem;
class InterSpec;
class SpectrumViewerTester;

namespace Wt
{
  class WText;
  class WMenu;
  class WSvgImage;
}//namespace Wt


class SpectrumViewerTesterWindow : public AuxWindow
{
public:
  SpectrumViewerTesterWindow( InterSpec *viewer );
  
  SpectrumViewerTester *m_tester;
};//class SpectrumViewerTesterWindow





class SpectrumViewerTester : public Wt::WContainerWidget
{
public:
  enum TestType
  {
    AutomatedPeakSearch,
    ManualPeakClicking,
    MultiplePeakInRoiFit,
    
    //PeakRefitTest
    
//    NuclideIdGuessing,
    //SourceShieldingFit
    //RecalibrationSpectrumTesting
    //RecalibrationPeakTesting
    NumTestType
  };//enum TestType

  static const char *tostr( TestType type );
  
  struct ScoreNote
  {
    ScoreNote();
    ~ScoreNote();
    std::string m_text;
    std::string m_title;
    
    //Note that WSvgImage was chosen over WRasterImage (png/jpeg) for 2 reasons:
    //  -WRasterImage would cause a crash often times on my mac
    //  -We can embed an SVG image into the HTML
    std::shared_ptr<Wt::WSvgImage> m_testImage;
    std::shared_ptr<Wt::WSvgImage> m_originalImage;
  };//struct ScoreNote
  
  struct Score
  {
    Score();
    double m_ntest;     //total number of things tested for
    double m_nwrong;    //number of things 'wrong'
    double m_ncorrect;
    int m_dbid;
    TestType m_test;
    std::vector<ScoreNote> m_notes;
    
    //addScore(...): add the numeric scores (m_ntest, m_nwrong, m_ncorrect) to
    //  *this, with the weight w. Does not add m_notes.
    void addScore( const Score &rhs, double w );
  };//struct Score
  
  class ScoreDisplay : public Wt::WContainerWidget
  {
  public:
    ScoreDisplay( const Score &score,
                  Wt::WContainerWidget *parent = 0 );
    void setScore( const Score &score );
    
    //renderScore(): renders m_score to (X)HTML text.
    void renderScore( std::ostream &strm );
    
    const Score &score() const;
    
    Score m_score;
    Wt::WContainerWidget *m_notes;
  };//class ScoreDisplay

  
public:
#if( BUILD_AS_OFFLINE_ANALYSIS_TEST_SUITE )
  static void doOfflineTesting();
  static void writeCssToHtmlHeader( std::ostream &strm );
  static void writeHtmlHeader( std::ostream &strm, const std::string &title );
  static void writeStateStaticHtmlHeader( std::ostream &strm,
                                          InterSpec *viewer,
                                          const UserState *userState );
#endif
  
  SpectrumViewerTester( InterSpec *viewer,
                        int widthHint, int heightHint,
                        Wt::WContainerWidget *parent = 0 );
  virtual ~SpectrumViewerTester();

  Score testAutomatedPeakSearch();
  Score testManualPeakClicking();
  Score testMultiplePeakFit();
  
  void setSizeHint( int width, int height );
  
  std::shared_ptr<Wt::WSvgImage> renderChartToSvg( double lowx, double upperx,
                                                     int width, int height );
  
  //scoreAllTests(): XXX - untested, doesnt actually work yet
  void scoreAllTests();
  
  void doScorring( TestType type, Wt::WContainerWidget *parent );
  Score doTest( TestType type );
  
  void updateOverview();
  
  static std::string makePeakSummarryTable( const PeakDef &peak );
  
  static bool peaksAreSimilar( const PeakDef *original, const PeakDef *found );
  static bool compatibleMeans( const PeakDef *orig, const PeakDef *found );
  static bool compatibleWidths( const PeakDef *orig, const PeakDef *found );
  static bool compatibleAmplitudes( const PeakDef *orig, const PeakDef *found );
  static bool compatibleLowerRoi( const PeakDef *orig, const PeakDef *found );
  static bool compatibleUpperRoi( const PeakDef *orig, const PeakDef *found );
  static bool sameNuclideId( const PeakDef *orig, const PeakDef *found );
  
protected:
  void renderScore( const Score &score );
  static double scoreWeight( TestType type );
  
  Score testMultiplePeakFitRangeVaried( double lowerXChangeFrac,
                                        double upperXChangeFrac );
  
protected:
  InterSpec *m_viewer;
  
  Wt::Signal<> m_done;
  
  Score m_scores[NumTestType];
  SideMenuItem *m_menuItems[NumTestType];
  ScoreDisplay *m_scoreDisplays[NumTestType];
  Wt::WContainerWidget *m_overview;
  Wt::WMenu *m_sideMenu;
  
  int m_picWidth;
  int m_picHeight;

public:
  static double sm_mean_max_nsigma_diff;
  static double sm_amplitude_max_frac_diff;
  static double sm_width_max_frac_diff;
  static const double sm_lower_roi_max_nsigma_diff;
  static const double sm_upper_roi_max_nsigma_diff;
  static const double sm_multipeakfit_nominal_weight;
  static const double sm_multipeakfit_range_max_vary_frac;

  static const int    sm_num_manual_click;
  static const double sm_nfwhm_manually_click;
  static const double sm_min_manually_click_corr_frac;
  
  friend class SpectrumViewerTesterWindow;
};//class SpectrumViewerTester



#endif //#ifndef SpectrumViewerTester_h
