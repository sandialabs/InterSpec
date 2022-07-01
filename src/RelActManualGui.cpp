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

#include <Wt/WMenu>
#include <Wt/WLabel>
#include <Wt/WPanel>
#include <Wt/WTable>
#include <Wt/WLabel>
#include <Wt/WComboBox>
#include <Wt/WLineEdit>
#include <Wt/WMenuItem>
#include <Wt/WResource>
#include <Wt/WTableRow>
#include <Wt/WTableCell>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/Http/Request>
#include <Wt/WApplication>
#include <Wt/Http/Response>
#include <Wt/WItemDelegate>
#include <Wt/WStackedWidget>
#include <Wt/WRegExpValidator>
#include <Wt/WContainerWidget>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/RelEffChart.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/RelActManualGui.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/IsotopeSelectionAids.h"

using namespace Wt;
using namespace std;

namespace
{
class RelActManualReportResource : public Wt::WResource
{
  Wt::WApplication *m_app;
  InterSpec *m_interspec;
  RelActManualGui *m_tool;
  
public:
  RelActManualReportResource( RelActManualGui *tool, InterSpec *viewer, WObject* parent = nullptr )
  : WResource( parent ), m_app( WApplication::instance() ), m_interspec( viewer ), m_tool( tool )
  {
    assert( m_app );
    assert( m_tool );
    assert( m_interspec );
  }
  
  virtual ~RelActManualReportResource()
  {
    beingDeleted();
  }
  
  virtual void handleRequest( const Wt::Http::Request &request, Wt::Http::Response &response )
  {
    assert( m_app );
    assert( m_interspec );
    
    try
    {
      WApplication::UpdateLock lock( m_app );
      
      if( !lock )
        throw std::runtime_error( "Error grabbing application lock to from RelActManualReportResource resource." );
      
      shared_ptr<SpecMeas> meas = m_interspec->measurment( SpecUtils::SpectrumType::Foreground );
      if( !meas )
        throw std::runtime_error( "Error getting spectrum file currently being shown." );
      
      auto foreground = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
      
      string filename = meas->filename();
      if( filename.empty() )
        filename = "rel_act";
      const string orig_extension = SpecUtils::file_extension(filename);
      if( orig_extension.size() && (orig_extension.size() < filename.size()) )
        filename = filename.substr(0,filename.size() - orig_extension.size());
      const string orig_file_name = filename;
      filename += ".html";
      
      //Remove bad filename characters
      const string notallowed = "\\/:?\"<>|*";
      for( auto it = begin(filename) ; it < end(filename) ; ++it )
      {
        if( notallowed.find(*it) != string::npos )
          *it = ' ';
      }
      
      suggestFileName( filename, WResource::Attachment );
      response.setMimeType( "application/octet-stream" );
      
      const shared_ptr<const RelActCalcManual::RelEffSolution> solution = m_tool->currentSolution();
      
      if( !solution )
      {
        response.out() << "<!DOCTYPE html>\n"
        "\t<head><meta charset=\"utf-8\"><title>No Rel. Activity Solution Available</title></head>"
        "\t<body>"
        "\t\tSorry - no solution currently available."
        "\t</body>"
        "</html>";
        
        return;
      }//if( !solution )
      
      
      PeakModel *peakModel = m_interspec->peakModel();
      deque<PeakModel::PeakShrdPtr> all_peaks;
      if( peakModel && peakModel->peaks() )
        all_peaks = *peakModel->peaks();
      
      vector<shared_ptr<const PeakDef>> display_peaks;
      
      for( const shared_ptr<const PeakDef> &p : all_peaks )
      {
        bool use_peak = false;
        for( const auto &r : solution->m_input_peak )
        {
          if( fabs(p->mean() - r.m_energy) < 1.0 )
            use_peak = true;
        }
        
        if( use_peak )
          display_peaks.push_back( p );
      }//for( const shared_ptr<const PeakDef> &p : all_peaks )
        
      
      string title;
      if( foreground && (foreground->title().length() > 0) )
        title = foreground->title();
      else
        title = orig_file_name;
      
      solution->print_html_report( response.out(), title, foreground, display_peaks );
    }catch( std::exception &e )
    {
      log("error") << "Error handling request for RelActManualReportResource: " << e.what();
      response.out() << "Error creating HTML                           file: " << e.what()
      << "\n\nPlease report to InterSpec@sandia.gov.";
      
      //passMessage( "Error getting spectrum file currently being shown",
      //              "", WarningWidget::WarningMsgHigh );
      
      response.setStatus(500);
      assert( 0 );
    }//try / catch
  }//void handleRequest(...)
};//class RelActManualReportResource


class ManRelEffNucDisp : public Wt::WPanel
{
public:
  const SandiaDecay::Nuclide * const m_nuc;
  double m_current_age;
  const bool m_age_is_settable;
  Wt::WLineEdit *m_age_edit;
  Wt::WTableRow *m_age_row;
  
public:
  ManRelEffNucDisp( const SandiaDecay::Nuclide * const nuc, double age, WContainerWidget *parent = nullptr )
  : WPanel( parent ),
   m_nuc( nuc ),
   m_current_age( (nuc && (age < 0.0)) ? PeakDef::defaultDecayTime(nuc) : age ),
   m_age_is_settable( nuc ? !PeakDef::ageFitNotAllowed(nuc) : false ),
   m_age_edit( nullptr ),
   m_age_row( nullptr )
  {
    assert( m_nuc );
    
    if( !m_nuc )
      throw runtime_error( "ManRelEffNucDisp: null nuc." );
    
    addStyleClass( "ManRelEffNucDisp" );
    
    setTitle( m_nuc->symbol );
    setCollapsible( true );
    setCollapsed( true );
    setAnimation( { WAnimation::AnimationEffect::SlideInFromTop,
      WAnimation::TimingFunction::Linear, 250 } );
    
    WTable *content = new WTable();
    content->addStyleClass( "NucInfoTable" );
    setCentralWidget( content );
    
    WTableCell *cell = content->elementAt(0, 0);
    WLabel *label = new WLabel( "Age", cell );
  
    m_age_row = content->rowAt(0);
    
    cell = content->elementAt(0, 1);
    
    string agestr;
    if( nuc->decaysToStableChildren() )
    {
      m_current_age = 0.0;
      agestr = "0y";
    }else
    {
      m_current_age = PeakDef::defaultDecayTime( nuc, &agestr );
    }//if( decay to stable only ) / else
    
    cout << "PhysicalUnits::sm_timeDurationHalfLiveOptionalRegex='" << PhysicalUnits::sm_timeDurationHalfLiveOptionalRegex << "'" << endl;
    
    if( m_age_is_settable )
    {
      m_age_edit = new WLineEdit( agestr, cell );
      
      WRegExpValidator *validator = new WRegExpValidator( PhysicalUnits::sm_timeDurationHalfLiveOptionalRegex, m_age_edit );
      validator->setFlags(Wt::MatchCaseInsensitive);
      m_age_edit->setValidator(validator);
      m_age_edit->setAutoComplete( false );
      label->setBuddy( m_age_edit );
      m_age_edit->addStyleClass( "AgeEdit" );
      
      m_age_edit->changed().connect( this, &ManRelEffNucDisp::handleAgeChange );
      m_age_edit->blurred().connect( this, &ManRelEffNucDisp::handleAgeChange );
      m_age_edit->enterPressed().connect( this, &ManRelEffNucDisp::handleAgeChange );
    }else
    {
      label = new WLabel( agestr, cell );
      label->addStyleClass( "FixedAge" );
    }//if( m_age_is_settable ) / else
    
    cell = content->elementAt(1, 0);
    //label = new WLabel( "Half Life", cell );
    label = new WLabel( "<span style=\"font-size: small;\">&lambda;<sub>&frac12;</sub></span>", cell );
    
    cell = content->elementAt(1, 1);
    WText *txt = new WText( PhysicalUnits::printToBestTimeUnits(nuc->halfLife), cell );
    
    
    cell = content->elementAt(2, 0);
    label = new WLabel( "Spec. Act.", cell );
    
    cell = content->elementAt(2, 1);
    const bool useCurrie = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    const double specificActivity = nuc->activityPerGram() / PhysicalUnits::gram;
    const string sa = PhysicalUnits::printToBestSpecificActivityUnits( specificActivity, 3, useCurrie );
    txt = new WText( sa, cell );
    
    // We could maybe list which gammas are currently being used
  }//ManRelEffNucDisp(...)
  
  void handleAgeChange()
  {
    assert( m_age_edit );
  }//void handleAgeChange()
  
  void setAgeHidden( const bool hidden )
  {
    m_age_row->setHidden( hidden );
  }
  
};//class ManRelEffNucDisp

}//namespace


RelActManualGui::RelActManualGui( InterSpec *viewer, Wt::WContainerWidget *parent )
: WContainerWidget( parent ),
  m_renderFlags( 0 ),
  m_currentSolution{},
  m_interspec( viewer ),
  m_layout( nullptr ),
  m_optionsColumn( nullptr ),
  m_relEffEqnForm( nullptr ),
  m_relEffEqnOrder( nullptr ),
  m_nucDataSrcHolder( nullptr ),
  m_nucDataSrc( nullptr ),
  m_matchTolerance( nullptr ),
  m_addUncertainty( nullptr ),
  m_downloadHtmlReport( nullptr ),
  m_peakTableColumn( nullptr ),
  m_peakModel( viewer ? viewer->peakModel() : nullptr ),
  m_nuclidesDisp( nullptr ),
  m_nucAge{},
  m_resultMenu( nullptr ),
  m_chart( nullptr ),
  m_results( nullptr )
{
  assert( m_interspec && m_peakModel );
  if( !m_interspec || !m_peakModel )
    throw runtime_error( "No interspec/peak model" );
  
  wApp->useStyleSheet( "InterSpec_resources/RelActManualGui.css" );
  
  addStyleClass( "EnergyCalTool RelActManualGui" );
  
  init();
}//RelActManualGui constructor
  
  

void RelActManualGui::init()
{
  m_renderFlags |= RenderActions::UpdateCalc;
  scheduleRender();
  
  if( m_layout )
    delete m_layout;
  
  m_layout = new WGridLayout( this );
  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_interspec );
  
  m_layout->setContentsMargins( 0, 0, 0, 0 );
  m_layout->setVerticalSpacing( 0 );
  m_layout->setHorizontalSpacing( 0 );
  
  //Create the options column
  m_optionsColumn = new WContainerWidget();
  m_optionsColumn->addStyleClass( "CalColumn MoreActionCol" );
  m_layout->addWidget( m_optionsColumn, 0, 0 );
  
  WGridLayout *collayout = new WGridLayout( m_optionsColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );
  
  WText *header = new WText( "Options" );
  header->addStyleClass( "ColHeader" );
  collayout->addWidget( header, 0, 0 );
  
  //We will put the apply-to list inside a div so we can style consistently with other rows
  // (a <ul> element doesnt accept same css as <div>, apparently).
  WContainerWidget *optionsDiv = new WContainerWidget();
  optionsDiv->addStyleClass( "CalColContent MoreActionsMenuContent" );
  collayout->addWidget( optionsDiv, 1, 0 );
  collayout->setRowStretch( 1, 1 );
  
  WContainerWidget *optionsList = new WContainerWidget( optionsDiv );
  optionsList->addStyleClass( "MoreActionsMenuList" );
  optionsList->setList( true );
  

  WContainerWidget *holder = new WContainerWidget( optionsList );
  WLabel *label = new WLabel( "Rel. Eff. Eqn Form", holder );
  label->setInline( false );
  
  m_relEffEqnForm = new WComboBox( optionsList );
  m_relEffEqnForm->activated().connect( this, &RelActManualGui::relEffEqnFormChanged );
  
  const char *tooltip = "The functional form to use for the relative efficiciency curve.<br />"
  "Options are:"
  "<table style=\"margin-left: 10px;\">"
  "<tr><th>Log(energy):</th>               <th>y = a + b*ln(x) + c*(ln(x))^2 + d*(ln(x))^3 + ...</th></tr>"
  "<tr><th>Log(rel. eff.):</th>            <th>y = exp( a + b*x + c/x + d/x^2 + e/x^3 + ... )</th></tr>"
  "<tr><th>Log(energy)Log(rel. eff.):</th> <th>y = exp( a  + b*(lnx) + c*(lnx)^2 + d*(lnx)^3 + ... )</th></tr>"
  "<tr><th>FRAM Empirical:</th>            <th>y = exp( a + b/x^2 + c*(lnx) + d*(lnx)^2 + e*(lnx)^3 )</th></tr>"
  "</table>";
  HelpSystem::attachToolTipOn( holder, tooltip, showToolTips );
  
  
  // Will assume FramEmpirical is the highest
  static_assert( static_cast<int>(RelActCalc::RelEffEqnForm::FramEmpirical)
                 > static_cast<int>(RelActCalc::RelEffEqnForm::LnXLnY),
                "RelEffEqnForm was changed!"
  );
  
  
  for( int i = 0; i <= static_cast<int>(RelActCalc::RelEffEqnForm::FramEmpirical); ++i )
  {
    const auto eqn_form = RelActCalc::RelEffEqnForm( i );
    
    const char *txt = "";
    switch( eqn_form )
    {
      case RelActCalc::RelEffEqnForm::LnX:
        //y = a + b*ln(x) + c*(ln(x))^2 + d*(ln(x))^3 + ...
        txt = "Log(x)";
        break;
        
      case RelActCalc::RelEffEqnForm::LnY:
        //y = exp( a + b*x + c/x + d/x^2 + e/x^3 + ... )
        txt = "Log(y)";
        break;
        
      case RelActCalc::RelEffEqnForm::LnXLnY:
        //y = exp( a  + b*(lnx) + c*(lnx)^2 + d*(lnx)^3 + ... )
        txt = "Log(x)Log(y)";
        break;
        
      case RelActCalc::RelEffEqnForm::FramEmpirical:
        //y = exp( a + b/x^2 + c*(lnx) + d*(lnx)^2 + e*(lnx)^3 )
        txt = "Empirical";
        break;
    }
    
    m_relEffEqnForm->addItem( txt );
  }//for( loop over RelEffEqnForm )
  
  m_relEffEqnForm->setCurrentIndex( static_cast<int>(RelActCalc::RelEffEqnForm::LnX) );
  
  // Blah Blah Blah - add other options here
  holder = new WContainerWidget( optionsList );
  tooltip = "The order (how many energy-dependent terms) relative efficiency equation to use.";
  HelpSystem::attachToolTipOn( holder, tooltip, showToolTips );
  
  label = new WLabel( "Rel. Eff. Eqn Order", holder );
  
  m_relEffEqnOrder = new WComboBox( holder );
  m_relEffEqnForm->activated().connect( this, &RelActManualGui::relEffEqnOrderChanged );
  
  m_relEffEqnOrder->addItem( "1" );
  m_relEffEqnOrder->addItem( "2" );
  m_relEffEqnOrder->addItem( "3" );
  m_relEffEqnOrder->addItem( "4" );
  m_relEffEqnOrder->addItem( "5" );
  m_relEffEqnOrder->addItem( "6" );
  m_relEffEqnOrder->setCurrentIndex( 2 );
  
  
  m_nucDataSrcHolder = new WContainerWidget( optionsList );
  holder = m_nucDataSrcHolder;
  
  tooltip = "The source for gamma branching ratios to use.";
  HelpSystem::attachToolTipOn( holder, tooltip, showToolTips );
  
  label = new WLabel( "B.R. src", holder );
  m_nucDataSrc = new WComboBox( holder );
  label->setBuddy( m_nucDataSrc );
  m_nucDataSrc->activated().connect( this, &RelActManualGui::nucDataSrcChanged );
  
  using RelActCalcManual::PeakCsvInput::NucDataSrc;
  for( NucDataSrc src = NucDataSrc(0); src < NucDataSrc::Undefined; src = NucDataSrc(static_cast<int>(src) + 1) )
  {
    const char *src_label = "";
    switch( src )
    {
      case NucDataSrc::Icrp107_U:         src_label = "ICRP 107";   break;
      case NucDataSrc::Lanl_U:            src_label = "LANL/FRAM";  break;
      case NucDataSrc::IcrpLanlGadras_U:  src_label = "SNL Combo";  break;
      case NucDataSrc::SandiaDecay:       src_label = "InterSpec";  break;
      case NucDataSrc::Undefined:         assert( 0 );              break;
    }//switch( src )
    
    m_nucDataSrc->addItem( WString::fromUTF8(src_label) );
  }//for( loop over sources )
  
  m_nucDataSrc->setCurrentIndex( static_cast<int>(NucDataSrc::SandiaDecay) );
  
  
  holder = new WContainerWidget( optionsList );
  label = new WLabel( "Match tol. (FWHM)", holder );

  tooltip = "The number of FWHM, from the peak mean, to include source gammas from as contributing"
  " to a peaks area.<br />"
  "For some photopeaks of some nuclides multiple gammas that are close in energy may contribute"
  " to creating a detected peak area.  This match tolerance specifies how many FWHM from the"
  " observed peak mean source gammas should be summed to determine the branching ratio to use."
  "<br />Specifying a value of zero will will cause only the gamma energy assigned to a peak to"
  " be used, even if there are very nearby other gammas.";
  HelpSystem::attachToolTipOn( holder, tooltip, showToolTips );
  
  m_matchTolerance = new NativeFloatSpinBox( holder );
  label->setBuddy( m_matchTolerance );
  m_matchTolerance->setSpinnerHidden();
  m_matchTolerance->setWidth( 30 );
  m_matchTolerance->setRange( 0, 5 );
  m_matchTolerance->setValue( 0.5 ); //Other places we use 1.25/2.355 = 0.530786
  
  
  holder = new WContainerWidget( optionsList );
  label = new WLabel( "Add. Uncert", holder );
  
  tooltip = "An additional uncertainty to add to the relative efficiency line, for each fit"
  " photopeak.<br />"
  "Small deviations in efficiency of detection for one or a few high statistics peaks can cause the"
  " efficiency curve to notably deviate from the other points if only statical uncertainties are"
  " used; when an additional uncertainty is added the relative efficiency will then do a better job"
  " of visibly going through all the data points, and from limited testing produce more accurate"
  " results.  You can think of this as adding a systematic uncertainty to each detected photopeak,"
  " that is uncorrelated between peaks.  From limited testing the value used is not hugely"
  " important, just as long as there is something.  You can also choose to use an unweighted fit,"
  " where each peak will contribute to the fit equally, no matter its statistical uncertainty.";
  
  HelpSystem::attachToolTipOn( holder, tooltip, showToolTips );
  
  m_addUncertainty = new WComboBox( holder );
  label->setBuddy( m_addUncertainty );
  m_addUncertainty->activated().connect( this, &RelActManualGui::addUncertChanged );
  
  for( AddUncert i = AddUncert(0); i < AddUncert::NumAddUncert; i = AddUncert(static_cast<int>(i) + 1) )
  {
    const char *uncert_txt = "";
    switch( i )
    {
      case AddUncert::Unweighted:         uncert_txt = "Unweighted"; break;
      case AddUncert::StatOnly:           uncert_txt = "Stat. Only"; break;
      case AddUncert::OnePercent:         uncert_txt = "1%";         break;
      case AddUncert::FivePercent:        uncert_txt = "5%";         break;
      case AddUncert::TenPercent:         uncert_txt = "10%";        break;
      case AddUncert::TwentyFivePercent:  uncert_txt = "25%";        break;
      case AddUncert::FiftyPercent:       uncert_txt = "50%";        break;
      case AddUncert::SeventyFivePercent: uncert_txt = "75%";        break;
      case AddUncert::OneHundredPercent:  uncert_txt = "100%";       break;
      case AddUncert::NumAddUncert:       assert(0);                 break;
    }//switch( i )
    
    m_addUncertainty->addItem( WString::fromUTF8(uncert_txt) );
  }//for( loop over AddUncert )
  
  m_addUncertainty->setCurrentIndex( static_cast<int>(AddUncert::FiftyPercent) );
  
  
  WContainerWidget *btndiv = new WContainerWidget();
  btndiv->addStyleClass( "BtmBtnDiv" );
  collayout->addWidget( btndiv, 2, 0 );
  
  auto helpBtn = new WContainerWidget( btndiv );
  helpBtn->addStyleClass( "Wt-icon ContentHelpBtn" );
  helpBtn->clicked().connect( boost::bind( &HelpSystem::createHelpWindow, "rel-act-manual" ) );
  
  m_htmlResource = new RelActManualReportResource( this, m_interspec, btndiv );
  
#if( BUILD_AS_OSX_APP )
  m_downloadHtmlReport = new WAnchor( WLink(m_htmlResource), btndiv );
  m_downloadHtmlReport->setTarget( AnchorTarget::TargetNewWindow );
  m_downloadHtmlReport->setStyleClass( "LinkBtn DownloadLink CALp" );
  m_downloadHtmlReport->setText( "CALp" );
#else
  m_downloadHtmlReport = new WPushButton( "HTML Report", btndiv );
  m_downloadHtmlReport->setIcon( "InterSpec_resources/images/download_small.svg" );
  m_downloadHtmlReport->setLinkTarget( Wt::TargetNewWindow );
  m_downloadHtmlReport->setStyleClass( "LinkBtn DownloadBtn CALp" );
  m_downloadHtmlReport->setLink( WLink(m_htmlResource) );
  
#if( ANDROID )
  // Using hacked saving to temporary file in Android, instead of via network download of file.
  m_downloadHtmlReport->clicked().connect( std::bind([this](){
    android_download_workaround( m_calpResource, "rel_eff.html");
  }) );
#endif //ANDROID
  
#endif
  
  WContainerWidget *nucCol = new WContainerWidget();
  nucCol->addStyleClass( "CalColumn RelActNucCol" );
  
  collayout = new WGridLayout( nucCol );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );

  header = new WText( "Nuclides" );
  header->addStyleClass( "ColHeader" );
  collayout->addWidget( header, 0, 0 );
  
  m_nuclidesDisp = new WContainerWidget();
  m_nuclidesDisp->addStyleClass( "CalColContent ApplyToMenuContent" );
  collayout->addWidget( m_nuclidesDisp, 1, 0 );
  collayout->setRowStretch( 1, 1 );
  m_layout->addWidget( nucCol, 0, 1 );
  
  // Create the "Cal Peaks" table
  m_peakTableColumn = new WContainerWidget();
  m_peakTableColumn->addStyleClass( "CalColumn PeakTableCol" );
  m_layout->addWidget( m_peakTableColumn, 0, 2 );
  //m_layout->setColumnStretch( 2, 1 );
  
  
  collayout = new WGridLayout( m_peakTableColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );
  
  header = new WText( "Peaks to Use" );
  header->addStyleClass( "ColHeader" );
  collayout->addWidget( header, 0, 0 );
  
  m_peakTable = new RowStretchTreeView();
  m_peakTable->addStyleClass( "CalColContent PeakTable" );
  collayout->addWidget( m_peakTable, 1, 0 );
  collayout->setRowStretch( 1, 1 );
  
  m_peakTable->setRootIsDecorated( false ); //makes the tree look like a table! :)
  m_peakTable->setModel( m_peakModel );
  const int numModelCol = m_peakModel->columnCount();
  for( int col = 0; col < numModelCol; ++col )
    m_peakTable->setColumnHidden( col, true );
  
  m_peakTable->setSortingEnabled( true );
  m_peakTable->setAlternatingRowColors( true );
  m_peakTable->setSelectable( true );
  m_peakTable->setSelectionMode( SingleSelection );
  m_peakTable->setEditTriggers( WAbstractItemView::SingleClicked
                               | WAbstractItemView::DoubleClicked );
  
  m_peakTable->setColumnHidden( PeakModel::kUseForManualRelEff, false );
  m_peakTable->setColumnHidden( PeakModel::kMean, false );
  m_peakTable->setColumnHidden( PeakModel::kIsotope, false );
  m_peakTable->setColumnHidden( PeakModel::kCps, false );
  m_peakTable->setColumnHidden( PeakModel::kPhotoPeakEnergy, false );
  
  // For a table width of 450 px, a total column width sum of 384 pixels seems to keep it so there
  //  is no horizontal scroll bar.
  m_peakTableColumn->setWidth( 450 );
  m_peakTable->setColumnWidth( PeakModel::kUseForManualRelEff, 76 );
  m_peakTable->setColumnWidth( PeakModel::kMean, 66);
  m_peakTable->setColumnWidth( PeakModel::kIsotope, 72 );
  m_peakTable->setColumnWidth( PeakModel::kCps, 76 );
  m_peakTable->setColumnWidth( PeakModel::kPhotoPeakEnergy, 94 );

  
  WItemDelegate *dblDelagate = new WItemDelegate( m_peakTable );
  dblDelagate->setTextFormat( "%.2f" );
  m_peakTable->setItemDelegateForColumn( PeakModel::kMean, dblDelagate );
  
  PhotopeakDelegate *nuclideDelegate = new PhotopeakDelegate( PhotopeakDelegate::NuclideDelegate, true, m_peakTable );
  m_peakTable->setItemDelegateForColumn( PeakModel::kIsotope, nuclideDelegate );
  
  PhotopeakDelegate *photopeakDelegate = new PhotopeakDelegate( PhotopeakDelegate::GammaEnergyDelegate, true, m_peakTable );
  m_peakTable->setItemDelegateForColumn( PeakModel::kPhotoPeakEnergy, photopeakDelegate );
  
  m_peakModel->dataChanged().connect( this, &RelActManualGui::handlePeaksChanged );
  m_peakModel->rowsRemoved().connect( this, &RelActManualGui::handlePeaksChanged );
  m_peakModel->rowsInserted().connect( this, &RelActManualGui::handlePeaksChanged );
  m_peakModel->layoutChanged().connect( this, &RelActManualGui::handlePeaksChanged );
  
  m_interspec->displayedSpectrumChanged().connect( this, &RelActManualGui::displayedSpectrumChanged );
  
  
  
  // Create the results column
  WContainerWidget *resCol = new WContainerWidget();
  resCol->addStyleClass( "CalColumn RelActResCol" );
  
  collayout = new WGridLayout( resCol );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  
  header = new WText( "Results" );
  header->addStyleClass( "ColHeader" );
  collayout->addWidget( header, 0, 0 );
  
  
  
  WContainerWidget *resultContent = new WContainerWidget();
  resultContent->addStyleClass( "CalColContent ApplyToMenuContent" );
  collayout->addWidget( resultContent, 1, 0 );
  collayout->setRowStretch( 1, 1 );
  m_layout->addWidget( resCol, 0, 3 );
  m_layout->setColumnStretch( 3, 1 );
  
  WGridLayout *resultsLayout = new WGridLayout( resultContent );
  
  
  WStackedWidget *stack = new WStackedWidget();
  stack->addStyleClass( "ResultStack" );
  stack->setOverflow( WContainerWidget::OverflowHidden, Wt::Orientation::Horizontal );
  stack->setOverflow( WContainerWidget::OverflowAuto, Wt::Orientation::Vertical );
  
  WAnimation animation(Wt::WAnimation::Fade, Wt::WAnimation::Linear, 200);
  stack->setTransitionAnimation( animation, true );
  
  m_resultMenu = new WMenu( stack, Wt::Vertical );
  m_resultMenu->addStyleClass( "ResultMenu LightNavMenu" );
  
  
  resultsLayout->addWidget( m_resultMenu, 0, 0 );
  resultsLayout->addWidget( stack, 0, 1 );
  resultsLayout->setColumnStretch( 1, 1 );
  
  m_results = new WContainerWidget();
  m_results->addStyleClass( "ResultsTxt" );
  WMenuItem *item = new WMenuItem( "Results", m_results );
  m_resultMenu->addItem( item );
  
  m_chart = new RelEffChart();
  item = new WMenuItem( "Chart", m_chart );
  m_resultMenu->addItem( item );
  
  displayedSpectrumChanged();
}//void init()



void RelActManualGui::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool renderFull = (flags & Wt::RenderFlag::RenderFull);
  //const bool renderUpdate = (flags & Wt::RenderFlag::RenderUpdate);
  
  WContainerWidget::render( flags );
  
  //if( renderFull )
  //  defineJavaScript();
  
  if( m_renderFlags.testFlag(RelActManualGui::RenderActions::UpdateCalc) )
    calculateSolution();
  
  m_renderFlags = 0;
}


void RelActManualGui::calculateSolution()
{
  m_currentSolution.reset();
  m_chart->setData( {}, {}, "" );
  m_results->clear();
  
  // TODO: should do the actual computation not on the GUI thread!
  try
  {
    using namespace RelActCalcManual;
    vector<GenericPeakInfo> peak_infos;
    
    if( !m_peakModel || !m_peakModel->peaks() )
      throw runtime_error( "No peaks avaiable" );
    
    double addUncert = -2.0;
    
    const AddUncert addUncertType = AddUncert( m_addUncertainty->currentIndex() );
    switch( addUncertType )
    {
        case AddUncert::Unweighted:         addUncert = -1.0; break;
        case AddUncert::StatOnly:           addUncert = 0.0;  break;
        case AddUncert::OnePercent:         addUncert = 0.01; break;
        case AddUncert::FivePercent:        addUncert = 0.05; break;
        case AddUncert::TenPercent:         addUncert = 0.1;  break;
        case AddUncert::TwentyFivePercent:  addUncert = 0.25; break;
        case AddUncert::FiftyPercent:       addUncert = 0.5;  break;
        case AddUncert::SeventyFivePercent: addUncert = 0.75; break;
        case AddUncert::OneHundredPercent:  addUncert = 1.0;  break;
        case AddUncert::NumAddUncert:       assert(0);        break;
    }//switch( i )
    
    if( addUncert < -1.0 )
      throw runtime_error( "Invalid add. uncert. selected." );
      
      
    const double matchToleranceFwhm = m_matchTolerance->value();
    
    vector<SandiaDecayNuc> nuclides_to_match_to;

    
    set<string> unique_isotopes;
    for( const PeakModel::PeakShrdPtr &p : *m_peakModel->peaks() )
    {
      if( p && p->parentNuclide() && p->useForManualRelEff() )
      {
        GenericPeakInfo peak;
        peak.m_energy = p->mean();
        peak.m_fwhm = p->fwhm();
        peak.m_counts = p->amplitude();
        peak.m_counts_uncert = p->amplitudeUncert();
        peak.m_base_rel_eff_uncert = addUncert;
        
        SandiaDecayNuc nuc;
        bool hadNuclide = false;
        for( const auto &existing : nuclides_to_match_to )
        {
          if( existing.nuclide == p->parentNuclide() )
          {
            nuc = existing;
            hadNuclide = true;
            break;
          }
        }//for( const auto &existing : nuclides_to_match_to )
        
        if( !hadNuclide )
        {
          nuc.nuclide = p->parentNuclide();
          nuc.age = -1.0;
          
          for( auto w : m_nuclidesDisp->children() )
          {
            ManRelEffNucDisp *rr = dynamic_cast<ManRelEffNucDisp *>(w);
            if( rr && (rr->m_nuc == nuc.nuclide) )
              nuc.age = rr->m_current_age;
          }//for( auto w : m_nuclidesDisp->children() )
          
          assert( nuc.age >= 0.0 );
          if( nuc.age < 0.0 )
            throw runtime_error( "Error finding age for " + nuc.nuclide->symbol );
          
          nuclides_to_match_to.push_back( nuc );
        }//if( pos != end(nuclides_to_match_to) ) / else
        
        
        if( matchToleranceFwhm == 0.0 )
        {
          GenericLineInfo info;
          info.m_yield = 0.0;
          info.m_isotope = p->parentNuclide()->symbol;
          
          const double decay_act_mult = SandiaDecay::MBq;
          SandiaDecay::NuclideMixture mix;
          mix.addAgedNuclideByActivity( nuc.nuclide, decay_act_mult, nuc.age );
          const vector<SandiaDecay::EnergyRatePair> gammas = mix.gammas( 0.0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy, true );
          for( const auto &g : gammas )
          {
            if( fabs(g.energy - p->gammaParticleEnergy()) < 0.01 )
              info.m_yield += (g.numPerSecond / decay_act_mult);
          }
          
          // Or could use:
          //vector<EnergyYield> decay_gammas( const SandiaDecay::Nuclide * const parent,
          //                                        const double age,
          //                                        const std::vector<double> &gammas_to_exclude )
          
          peak.m_source_gammas.push_back( info );
        }//if( matchToleranceFwhm == 0.0 )
        
        peak_infos.push_back( peak );
        
        unique_isotopes.insert( p->parentNuclide()->symbol );
      }//
    }//for( const PeakModel::PeakShrdPtr &p : *m_peakModel->peaks() )
    
    const vector<SandiaDecayNuc> orig_nuclides_to_match_to = nuclides_to_match_to;
    
    bool validSrcData = false;
    const auto srcData = RelActCalcManual::PeakCsvInput::NucDataSrc( m_nucDataSrc->currentIndex() );
    
    switch( srcData )
    {
      case RelActCalcManual::PeakCsvInput::NucDataSrc::SandiaDecay:
        validSrcData = true;
        break;
        
      case RelActCalcManual::PeakCsvInput::NucDataSrc::Icrp107_U:
      case RelActCalcManual::PeakCsvInput::NucDataSrc::Lanl_U:
      case RelActCalcManual::PeakCsvInput::NucDataSrc::IcrpLanlGadras_U:
      {
        validSrcData = true;
        //Remove Uranium nuclides from nuclides_to_match_to
        const auto nucs_copy = nuclides_to_match_to;
        nuclides_to_match_to.clear();
        for( const auto &n : nucs_copy )
        {
          if( n.nuclide && (n.nuclide->atomicNumber != 92) )
            nuclides_to_match_to.push_back( n );
        }
        break;
      }//case( not using SandiaDecay for uranium )
        
      case RelActCalcManual::PeakCsvInput::NucDataSrc::Undefined:
        validSrcData = false;
        break;
    }//switch( srcData )
    
    if( (matchToleranceFwhm != 0.0) && nuclides_to_match_to.size() )
      peak_infos = add_nuclides_to_peaks( peak_infos, nuclides_to_match_to, matchToleranceFwhm*2.634 );
    
    switch( srcData )
    {
      case RelActCalcManual::PeakCsvInput::NucDataSrc::SandiaDecay:
      case RelActCalcManual::PeakCsvInput::NucDataSrc::Undefined:
        break;
        
      case RelActCalcManual::PeakCsvInput::NucDataSrc::Icrp107_U:
      case RelActCalcManual::PeakCsvInput::NucDataSrc::Lanl_U:
      case RelActCalcManual::PeakCsvInput::NucDataSrc::IcrpLanlGadras_U:
      {
        // should re-factor fill_in_nuclide_info to use SandiaDecay for non-U isotopes
        //  And see if it makes sense to use this function instead of add_nuclides_to_peaks(...) everywhere.
        std::vector<RelActCalcManual::PeakCsvInput::NucAndAge> isotopes;
        
        for( const auto &n : orig_nuclides_to_match_to )
        {
          if( n.nuclide && (n.nuclide->atomicNumber != 92) )
            isotopes.emplace_back( n.nuclide->symbol, n.age );
        }
        
        if( isotopes.empty() )
          throw runtime_error( "A uranium specialized data source was selected, but there is no uranium in this problem." );
        
        RelActCalcManual::PeakCsvInput::NucMatchResults matched_res
            = RelActCalcManual::PeakCsvInput::fill_in_nuclide_info( peak_infos, srcData,
                                                            {}, isotopes, matchToleranceFwhm, {} );
        if( matched_res.peaks_matched.empty() )
          throw runtime_error( "No uranium peaks were matched to the specialized uranium nuclear data." );
        
        peak_infos = matched_res.peaks_matched;
        peak_infos.insert( end(peak_infos),
                          begin(matched_res.peaks_not_matched),
                          end(matched_res.peaks_not_matched) );
        
        std::sort( begin(peak_infos), end(peak_infos),
          []( const RelActCalcManual::GenericPeakInfo &lhs, const RelActCalcManual::GenericPeakInfo &rhs ){
            return lhs.m_energy < rhs.m_energy;
        });
      }//case not SandiaDecay
    }//switch( srcData )
    
    
    
    if( !validSrcData )
      throw runtime_error( "Invalid SourceData" );
    
    
    bool validSelect = false;
    const auto eqn_form = RelActCalc::RelEffEqnForm( m_relEffEqnForm->currentIndex() );
    switch( eqn_form )
    {
      case RelActCalc::RelEffEqnForm::LnX:
      case RelActCalc::RelEffEqnForm::LnY:
      case RelActCalc::RelEffEqnForm::LnXLnY:
      case RelActCalc::RelEffEqnForm::FramEmpirical:
        validSelect = true;
        break;
    }//switch( eqn_form )
    
    if( !validSelect )
      throw runtime_error( "Invalid RelEffEqnForm" );
    
    const int orderIndex = m_relEffEqnOrder->currentIndex();
    if( (orderIndex < 0) || (orderIndex > 7) )
      throw runtime_error( "Invalid RelEffEqnOrder" );
    const size_t eqn_order = static_cast<size_t>(orderIndex) + 1;
    
    
    const RelEffSolution solution = solve_relative_efficiency( peak_infos, eqn_form, eqn_order );
  
    m_currentSolution = make_shared<RelEffSolution>( solution );
    
    updateGuiWithResults();
  }catch( std::exception &e )
  {
    cout << "Error doing calc: " << e.what() << endl;
    
    //Set error to results TXT and show results TXT tab.
    
    string errormsg = "Error calculating relative activity/efficiency: " + string(e.what());
    WText *error = new WText( errormsg, m_results );
    error->addStyleClass( "CalcError" );
    m_resultMenu->select( 0 );
  }//try / catch
}//void calculateSolution()


void RelActManualGui::updateGuiWithResults()
{
  m_results->clear();
  
  if( !m_currentSolution )
  {
    WText *txt = new WText( "No results available." );
    txt->setInline( false );
    txt->addStyleClass( "NoCalcResults" );
    return;
  }//
  
  const RelActCalcManual::RelEffSolution &solution = *m_currentSolution;
  
  switch( solution.m_status )
  {
    case RelActCalcManual::ManualSolutionStatus::NotInitialized:
    case RelActCalcManual::ManualSolutionStatus::ErrorInitializing:
    case RelActCalcManual::ManualSolutionStatus::ErrorFindingSolution:
    case RelActCalcManual::ManualSolutionStatus::ErrorGettingSolution:
    {
      m_chart->setData( {}, {}, "" );
      break;
    }
      
    case RelActCalcManual::ManualSolutionStatus::Success:
      break;
  }//switch( solution.m_status )
  
  
  if( !solution.m_error_message.empty() )
  {
    WText *errtxt = new WText( "Error: " + solution.m_error_message , m_results );
    errtxt->setInline( false );
    errtxt->addStyleClass( "CalcError" );
  }
  
  
  for( string warning : solution.m_warnings )
  {
    WText *warntxt = new WText( "Warning: " + warning, m_results );
    warntxt->setInline( false );
    warntxt->addStyleClass( "CalcWarning" );
  }
  
  
  if( solution.m_status != RelActCalcManual::ManualSolutionStatus::Success )
    return;
  
  // We'll first update the chart
  const string relEffEqn = RelActCalc::rel_eff_eqn_js_function( solution.m_rel_eff_eqn_form, solution.m_rel_eff_eqn_coefficients );
  
  deque<PeakModel::PeakShrdPtr> displayed_peaks;
  if( m_peakModel && m_peakModel->peaks() )
    displayed_peaks = *m_peakModel->peaks();
  
  map<string, pair<double,string>> relActsColors;
  for( const auto &act : solution.m_rel_activities )
  {
    // Look through PeakModel for a peak of this nuclide, and use that color; inefficient, but good
    //  enough for now.
    string color;
    for( const auto &p : displayed_peaks )
    {
      if( p && p->parentNuclide()
         && (p->parentNuclide()->symbol == act.m_isotope)
         && !p->lineColor().isDefault() )
      {
        color = p->lineColor().cssText();
        break;
      }
    }//for( loop over peaks to look for colors )
    
    relActsColors[act.m_isotope] = make_pair( act.m_rel_activity, color );
  }//for( const auto &act : solution.m_rel_activities )
  
  m_chart->setData( solution.m_input_peak, relActsColors, relEffEqn );
  
  
  // Now update the text
  // TODO: refactor getting these tables from the solution...
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  char buffer[2048] = {'\0'};
  stringstream results_html;
  results_html << "<div>&chi;<sup>2</sup>=" << PhysicalUnits::printCompact( solution.m_chi2, 4)
  << " and there were " << solution.m_dof << " DOF (&chi;<sup>2</sup>/<sub>DOF</sub>="
  << PhysicalUnits::printCompact(solution.m_chi2/solution.m_dof, 4) << ")</div>\n";
  
  results_html << "<div class=\"releffeqn\">Rel. Eff. Eqn: y = "
  << RelActCalc::rel_eff_eqn_text( solution.m_rel_eff_eqn_form, solution.m_rel_eff_eqn_coefficients )
  << "</div>\n";
  
  solution.get_mass_fraction_table( results_html );
  solution.get_mass_ratio_table( results_html );
  
  new WText( results_html.str(), m_results );
}//void updateGuiWithResults( const RelActCalcManual::RelEffSolution &solution );


void RelActManualGui::relEffEqnFormChanged()
{
  m_renderFlags |= RenderActions::UpdateCalc;
  scheduleRender();
}//void relEffEqnFormChanged()


void RelActManualGui::relEffEqnOrderChanged()
{
  m_renderFlags |= RenderActions::UpdateCalc;
  scheduleRender();
}


void RelActManualGui::nucDataSrcChanged()
{
  m_renderFlags |= RenderActions::UpdateCalc;
  scheduleRender();
}


void RelActManualGui::addUncertChanged()
{
  m_renderFlags |= RenderActions::UpdateCalc;
  scheduleRender();
}


void RelActManualGui::handlePeaksChanged()
{
  bool showAge = false;
  const auto srcData = RelActCalcManual::PeakCsvInput::NucDataSrc( m_nucDataSrc->currentIndex() );
  switch( srcData )
  {
    case RelActCalcManual::PeakCsvInput::NucDataSrc::SandiaDecay:
      showAge = true;
      break;
      
    case RelActCalcManual::PeakCsvInput::NucDataSrc::Icrp107_U:
    case RelActCalcManual::PeakCsvInput::NucDataSrc::Lanl_U:
    case RelActCalcManual::PeakCsvInput::NucDataSrc::IcrpLanlGadras_U:
    case RelActCalcManual::PeakCsvInput::NucDataSrc::Undefined:
      showAge = false;
      break;
  }//switch( srcData )
  
  map<const SandiaDecay::Nuclide *,ManRelEffNucDisp *> current_nucs;
  
  for( auto w : m_nuclidesDisp->children() )
  {
    ManRelEffNucDisp *rr = dynamic_cast<ManRelEffNucDisp *>(w);
    assert( rr && rr->m_nuc && !current_nucs.count(rr->m_nuc) );
    if( rr && rr->m_nuc )
      current_nucs[rr->m_nuc] = rr;
  }//for( auto w : m_nuclidesDisp->children() )
  
  set<const SandiaDecay::Nuclide *> nucs_in_peaks;
  PeakModel *peakModel = m_interspec->peakModel();
  const auto peaks = peakModel ? peakModel->peaks() : nullptr;
  if( peaks )
  {
    for( const auto &p : *peaks )
    {
      if( p->parentNuclide() && p->useForManualRelEff() )
        nucs_in_peaks.insert( p->parentNuclide() );
    }//
  }//if( peaks )
  
  // Go through and add any new nuclide from peaks
  for( const SandiaDecay::Nuclide *nuc : nucs_in_peaks )
  {
    auto pos = current_nucs.find( nuc );
    if( pos != end(current_nucs) )
      continue;
    
    double age = -1;
    const auto agePos = m_nucAge.find(nuc->symbol);
    if( agePos != end(m_nucAge) )
      age = agePos->second;
    
    // We'll insert the display sorted by {atomicNumber, atomicMass}
    int insert_index = 0;
    const auto kids = m_nuclidesDisp->children();
    for( size_t index = 0; index < kids.size(); ++index )
    {
      ManRelEffNucDisp *prev = dynamic_cast<ManRelEffNucDisp *>( kids[index] );
      assert( prev );
      if( !prev || !prev->m_nuc )
        continue;
      
      if( prev->m_nuc->atomicNumber > nuc->atomicNumber )
      {
        insert_index = static_cast<int>( index );
        break;
      }else if( (prev->m_nuc->atomicNumber == nuc->atomicNumber)
               && (prev->m_nuc->atomicMass > nuc->atomicMass) )
      {
        insert_index = static_cast<int>( index );
        break;
      }else
      {
        insert_index = static_cast<int>( index + 1 );
      }
    }//for( loop over existing displays to find position )
    
    ManRelEffNucDisp *rr = new ManRelEffNucDisp( nuc, age );
    m_nuclidesDisp->insertWidget( insert_index, rr );
  }//for( loop over to add displays for new nuclides )
  
  
  // Loop over and remove displays for any nuclides we no longer need
  for( const auto nuc_widget : current_nucs )
  {
    if( nucs_in_peaks.count(nuc_widget.first) )
      continue;
    
    m_nucAge[nuc_widget.first->symbol] = nuc_widget.second->m_current_age;
    
    delete nuc_widget.second;
  }//for( loop over to remove any nuclides )
  
  // We may have deleted some of current_nucs, so lets clear it, just to make sure we dont access
  current_nucs.clear();
  
  
  for( auto w : m_nuclidesDisp->children() )
  {
    ManRelEffNucDisp *rr = dynamic_cast<ManRelEffNucDisp *>(w);
    if( rr )
    {
      const bool isU = (rr->m_nuc->atomicNumber == 92);
      rr->setAgeHidden( isU ? !showAge : false );
    }else
    {
      assert( 0 );
    }
  }//for( auto w : m_nuclidesDisp->children() )
  
  
  m_renderFlags |= RenderActions::UpdateCalc;
  scheduleRender();
}//void handlePeaksChanged()


void RelActManualGui::displayedSpectrumChanged()
{
  handlePeaksChanged();
  
  m_renderFlags |= RenderActions::UpdateCalc;
  scheduleRender();
}//void displayedSpectrumChanged()


shared_ptr<const RelActCalcManual::RelEffSolution> RelActManualGui::currentSolution()
{
  return m_currentSolution;
}
