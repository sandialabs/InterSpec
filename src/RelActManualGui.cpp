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
#include <Wt/WCheckBox>
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
#include "SpecUtils/RapidXmlUtils.hpp"

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

const int RelActManualGui::sm_xmlSerializationMajorVersion = 0;
const int RelActManualGui::sm_xmlSerializationMinorVersion = 1;

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
  m_backgroundSubtract( nullptr ),
  m_backgroundSubtractHolder( nullptr ),
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
  m_renderFlags |= RenderActions::UpdateNuclides;
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
  m_optionsColumn->addStyleClass( "ToolTabTitledColumn OptionsCol" );
  m_layout->addWidget( m_optionsColumn, 0, 0 );
  
  WGridLayout *collayout = new WGridLayout( m_optionsColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );
  
  WText *header = new WText( "Options" );
  header->addStyleClass( "ToolTabColumnTitle" );
  collayout->addWidget( header, 0, 0 );
  
  //We will put the apply-to list inside a div so we can style consistently with other rows
  // (a <ul> element doesnt accept same css as <div>, apparently).
  WContainerWidget *optionsDiv = new WContainerWidget();
  optionsDiv->addStyleClass( "ToolTabTitledColumnContent OptionsColContent" );
  collayout->addWidget( optionsDiv, 1, 0 );
  collayout->setRowStretch( 1, 1 );
  
  WTable *optionsList = new WTable( optionsDiv );
  optionsList->addStyleClass( "OptionsList" );
  WLabel *label = new WLabel( "Eqn Form", optionsList->elementAt(0, 0) );
  
  m_relEffEqnForm = new WComboBox( optionsList->elementAt(0, 1) );
  m_relEffEqnForm->activated().connect( this, &RelActManualGui::relEffEqnFormChanged );
  
  const char *tooltip = "The functional form to use for the relative efficiciency curve.<br />"
  "Options are:"
  "<table style=\"margin-left: 10px;\">"
  "<tr><th>Log(energy):</th>               <th>y = a + b*ln(x) + c*(ln(x))^2 + d*(ln(x))^3 + ...</th></tr>"
  "<tr><th>Log(rel. eff.):</th>            <th>y = exp( a + b*x + c/x + d/x^2 + e/x^3 + ... )</th></tr>"
  "<tr><th>Log(energy)Log(rel. eff.):</th> <th>y = exp( a  + b*(lnx) + c*(lnx)^2 + d*(lnx)^3 + ... )</th></tr>"
  "<tr><th>FRAM Empirical:</th>            <th>y = exp( a + b/x^2 + c*(lnx) + d*(lnx)^2 + e*(lnx)^3 )</th></tr>"
  "</table>";
  HelpSystem::attachToolTipOn( optionsList->elementAt(0,0), tooltip, showToolTips );
  HelpSystem::attachToolTipOn( optionsList->elementAt(0,1), tooltip, showToolTips );
  
  
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
  
  label = new WLabel( "Eqn Order", optionsList->elementAt(1, 0) );
  
  m_relEffEqnOrder = new WComboBox( optionsList->elementAt(1, 1) );
  m_relEffEqnOrder->activated().connect( this, &RelActManualGui::relEffEqnOrderChanged );
  
  m_relEffEqnOrder->addItem( "1" );
  m_relEffEqnOrder->addItem( "2" );
  m_relEffEqnOrder->addItem( "3" );
  m_relEffEqnOrder->addItem( "4" );
  m_relEffEqnOrder->addItem( "5" );
  m_relEffEqnOrder->addItem( "6" );
  m_relEffEqnOrder->setCurrentIndex( 2 );
  
  
  tooltip = "The order (how many energy-dependent terms) relative efficiency equation to use.";
  HelpSystem::attachToolTipOn( optionsList->elementAt(1, 0), tooltip, showToolTips );
  HelpSystem::attachToolTipOn( optionsList->elementAt(1, 1), tooltip, showToolTips );
  
  label = new WLabel( "Yield Info", optionsList->elementAt(2, 0) );
  m_nucDataSrc = new WComboBox( optionsList->elementAt(2, 1) );
  label->setBuddy( m_nucDataSrc );
  m_nucDataSrc->activated().connect( this, &RelActManualGui::nucDataSrcChanged );
  
  tooltip = "The nuclear data source for gamma branching ratios.";
  HelpSystem::attachToolTipOn( optionsList->elementAt(2, 0), tooltip, showToolTips );
  HelpSystem::attachToolTipOn( optionsList->elementAt(2, 1), tooltip, showToolTips );
  
  
  using RelActCalcManual::PeakCsvInput::NucDataSrc;
  for( NucDataSrc src = NucDataSrc(0); src < NucDataSrc::Undefined; src = NucDataSrc(static_cast<int>(src) + 1) )
  {
    const char *src_label = "";
    switch( src )
    {
      case NucDataSrc::Icrp107_U:         src_label = "ICRP 107";   break;
      case NucDataSrc::Lanl_U:            src_label = "FRAM";       break;
      case NucDataSrc::IcrpLanlGadras_U:  src_label = "Combo";      break;
      case NucDataSrc::SandiaDecay:       src_label = "InterSpec";  break;
      case NucDataSrc::Undefined:         assert( 0 );              break;
    }//switch( src )
    
    m_nucDataSrc->addItem( WString::fromUTF8(src_label) );
  }//for( loop over sources )
  
  m_nucDataSrc->setCurrentIndex( static_cast<int>(NucDataSrc::SandiaDecay) );
  
  m_nucDataSrcHolder = optionsList->rowAt(2);

  label = new WLabel( "Match tol.", optionsList->elementAt(3, 0) ); //(FWHM)
  m_matchTolerance = new NativeFloatSpinBox( optionsList->elementAt(3, 1) );
  label->setBuddy( m_matchTolerance );
  m_matchTolerance->setSpinnerHidden();
  m_matchTolerance->setWidth( 35 );
  m_matchTolerance->setRange( 0, 5 );
  m_matchTolerance->setValue( 0.5 ); //Other places we use 1.25/2.355 = 0.530786
  label = new WLabel( "&nbsp;FWHM", optionsList->elementAt(3, 1) );
  label->setBuddy( m_matchTolerance );
  m_matchTolerance->valueChanged().connect( this, &RelActManualGui::matchToleranceChanged );
  
  
  tooltip = "The number of FWHM, from the peak mean, to include source gammas from as contributing"
  " to a peaks area.<br />"
  "For some photopeaks of some nuclides multiple gammas that are close in energy may contribute"
  " to creating a detected peak area.  This match tolerance specifies how many FWHM from the"
  " observed peak mean source gammas should be summed to determine the branching ratio to use."
  "<br />Specifying a value of zero will will cause only the gamma energy assigned to a peak to"
  " be used, even if there are very nearby other gammas.";
  HelpSystem::attachToolTipOn( optionsList->elementAt(3, 0), tooltip, showToolTips );
  HelpSystem::attachToolTipOn( optionsList->elementAt(3, 1), tooltip, showToolTips );
  
  label = new WLabel( "Add. Uncert", optionsList->elementAt(4, 0) );
  m_addUncertainty = new WComboBox( optionsList->elementAt(4, 1) );
  label->setBuddy( m_addUncertainty );
  
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
  
  HelpSystem::attachToolTipOn( optionsList->elementAt(4, 0), tooltip, showToolTips );
  HelpSystem::attachToolTipOn( optionsList->elementAt(4, 1), tooltip, showToolTips );
  
  
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
  
  
  m_backgroundSubtract = new WCheckBox( "Background Subtract", optionsList->elementAt(5, 0) );
  m_backgroundSubtract->addStyleClass( "BackSub" );
  optionsList->elementAt(5, 0)->setColumnSpan( 2 );
  m_backgroundSubtractHolder = optionsList->rowAt(5);
  m_backgroundSubtract->checked().connect( this, &RelActManualGui::backgroundSubtractChanged );
  m_backgroundSubtract->unChecked().connect( this, &RelActManualGui::backgroundSubtractChanged );
  
  
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
  m_downloadHtmlReport->setText( "HTML Report" );
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
  nucCol->addStyleClass( "ToolTabTitledColumn RelActNucCol" );
  
  collayout = new WGridLayout( nucCol );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );

  header = new WText( "Nuclides" );
  header->addStyleClass( "ToolTabColumnTitle" );
  collayout->addWidget( header, 0, 0 );
  
  m_nuclidesDisp = new WContainerWidget();
  m_nuclidesDisp->addStyleClass( "ToolTabTitledColumnContent" );
  collayout->addWidget( m_nuclidesDisp, 1, 0 );
  collayout->setRowStretch( 1, 1 );
  m_layout->addWidget( nucCol, 0, 1 );
  
  // Create the "Cal Peaks" table
  m_peakTableColumn = new WContainerWidget();
  m_peakTableColumn->addStyleClass( "ToolTabTitledColumn PeakTableCol" );
  m_layout->addWidget( m_peakTableColumn, 0, 2 );
  //m_layout->setColumnStretch( 2, 1 );
  
  
  collayout = new WGridLayout( m_peakTableColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );
  
  header = new WText( "Peaks to Use" );
  header->addStyleClass( "ToolTabColumnTitle" );
  collayout->addWidget( header, 0, 0 );
  
  m_peakTable = new RowStretchTreeView();
  m_peakTable->addStyleClass( "ToolTabTitledColumnContent PeakTable" );
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
  resCol->addStyleClass( "ToolTabTitledColumn RelActResCol" );
  
  collayout = new WGridLayout( resCol );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  
  header = new WText( "Results" );
  header->addStyleClass( "ToolTabColumnTitle" );
  collayout->addWidget( header, 0, 0 );
  
  
  
  WContainerWidget *resultContent = new WContainerWidget();
  resultContent->addStyleClass( "ToolTabTitledColumnContent ResultColumnContent" );
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
  
  // When outside the link area is clicked, the item doesnt get selected, so we'll work around this.
  item->clicked().connect( std::bind([this,item](){
    m_resultMenu->select( item );
    item->triggered().emit( item );
  }) );
  
  m_chart = new RelEffChart();
  item = new WMenuItem( "Chart", m_chart );
  
  item->clicked().connect( std::bind([this,item](){
    m_resultMenu->select( item );
    item->triggered().emit( item );
  }) );
  
  m_resultMenu->addItem( item );
  
  displayedSpectrumChanged();
}//void init()



void RelActManualGui::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  //const bool renderFull = (flags & Wt::RenderFlag::RenderFull);
  //if( renderFull )
  //  defineJavaScript();
  
  
  if( m_renderFlags.testFlag(RelActManualGui::RenderActions::UpdateSpectrumOptions) )
  {
    updateSpectrumBasedOptions();
    m_renderFlags |= RelActManualGui::RenderActions::UpdateCalc; //just to make sure
  }
    
  
  if( m_renderFlags.testFlag(RelActManualGui::RenderActions::UpdateNuclides) )
  {
    updateNuclides();
    m_renderFlags |= RelActManualGui::RenderActions::UpdateCalc; //just to make sure
  }
  
  
  if( m_renderFlags.testFlag(RelActManualGui::RenderActions::UpdateCalc) )
  {
    calculateSolution();
  }
  
  m_renderFlags = 0;
  
  WContainerWidget::render( flags );
}


void RelActManualGui::calculateSolution()
{
  m_currentSolution.reset();
  m_chart->setData( vector<RelActCalcManual::GenericPeakInfo>{}, {}, "" );
  m_results->clear();
  
  // TODO: should do the actual computation not on the GUI thread!
  try
  {
    using namespace RelActCalcManual;
    vector<string> prep_warnings;
    vector<GenericPeakInfo> peak_infos;
    
    InterSpec *viewer = InterSpec::instance();
    assert( viewer );
    if( !viewer )
      throw runtime_error( "Not in GUI thread???" );
    
    if( !m_peakModel || !m_peakModel->peaks() )
      throw runtime_error( "No peaks available" );
    
    // PeakModel stores all PeakDefs as const, so I think it is fine to just make a copy of
    //  the deque, and then it should be thread-safe to access the PeakDef objects from off
    //  of the main thread.
    const deque<shared_ptr<const PeakDef>> peaks = *m_peakModel->peaks();
    
    const shared_ptr<const SpecUtils::Measurement> fore_spec = viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    const shared_ptr<const SpecUtils::Measurement> back_spec = viewer->displayedHistogram(SpecUtils::SpectrumType::Background);
    const std::shared_ptr<const SpecMeas> back_meas = viewer->measurment(SpecUtils::SpectrumType::Background);
    
    const double back_sub_nsigma_near = 1.0; // Fairly arbitrary.  TODO: have this be a user setable value?
    const float foreground_live_time = fore_spec ? fore_spec->live_time() : 1.0f;
    const float background_live_time = back_spec ? back_spec->live_time() : 1.0f;
    const bool background_sub = (!m_backgroundSubtractHolder->isHidden()
                                 && m_backgroundSubtract->isChecked());
    
    deque<shared_ptr<const PeakDef>> background_peaks;
    if( background_sub && back_spec && back_meas && (background_live_time > 0.0) )
    {
      const auto &displayed = viewer->displayedSamples(SpecUtils::SpectrumType::Background);
      shared_ptr<const deque<shared_ptr<const PeakDef>>> backpeaks = back_meas->peaks( displayed );
      if( backpeaks )
        background_peaks = *backpeaks;
    }//if( background_sub && back_spec && back_meas )
    
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
      
    const double match_tol_sigma = 2.35482 * m_matchTolerance->value();
    
    const RelActCalcManual::PeakCsvInput::NucDataSrc srcData = nucDataSrc();
    
    map<const SandiaDecay::Nuclide *,double> nuclide_ages;
    for( auto w : m_nuclidesDisp->children() )
    {
      ManRelEffNucDisp *rr = dynamic_cast<ManRelEffNucDisp *>(w);
      if( rr && rr->m_nuc )
        nuclide_ages[rr->m_nuc] = rr->m_current_age;
    }//for( auto w : m_nuclidesDisp->children() )
    
    
    
    // I *think* (but worth a second check) everything below here, until results are updated to the
    //  GUI could be done off the main thread, and then post'd to update the gui at the end
    
    set<pair<float,float>> energy_cal_match_warning_energies;
    
    vector<SandiaDecayNuc> nuclides_to_match_to;
    
    size_t num_peaks_back_sub = 0;
    double lowest_energy_peak = 3000;
    bool has_U_or_Pu = false;
    set<string> unique_isotopes;
    for( const PeakModel::PeakShrdPtr &p : peaks )
    {
      if( p && p->parentNuclide() && p->useForManualRelEff() )
      {
        GenericPeakInfo peak;
        peak.m_energy = p->mean();
        peak.m_fwhm = p->gausPeak() ? p->fwhm() : (2.35482 * 0.25 * p->roiWidth());
        peak.m_counts = p->amplitude();
        peak.m_counts_uncert = p->amplitudeUncert();
        peak.m_base_rel_eff_uncert = addUncert;
        
        lowest_energy_peak = std::min( lowest_energy_peak, peak.m_energy );
        
        if( background_sub )
        {
          const double sigma = p->gausPeak() ? p->sigma() : 0.25*p->roiWidth();
          const double scale = foreground_live_time / background_live_time;
          
          double back_counts = 0.0, back_uncert_2 = 0.0;
          for( const shared_ptr<const PeakDef> &back_peak : background_peaks )
          {
            // In principle the peak shouldnt need to be a gaussian peak - but this has yet to be
            //  tested
            //if( !back_peak->gausPeak() )
            //  continue;
            
            if( fabs(back_peak->mean() - p->mean()) < (back_sub_nsigma_near * sigma) )
            {
              back_counts += scale * back_peak->peakArea();
              back_uncert_2 += scale * scale * back_peak->peakAreaUncert() * back_peak->peakAreaUncert();
            }//if( fabs(backPeak.mean()-peak.mean()) < sigma )
            
            if( back_counts > 0.0 )
            {
              num_peaks_back_sub += 1;
              peak.m_counts -= back_counts;
              peak.m_counts_uncert = sqrt( peak.m_counts_uncert*peak.m_counts_uncert + back_uncert_2 );
            }
          }//for( const PeakDef &peak : backPeaks )
          
          if( peak.m_counts <= 0.0 )
          {
            stringstream msg;
            msg << "After background subtraction, peak at "
            << std::fixed << std::setprecision(2) << peak.m_energy << " keV had negative counts"
            << " so was not used.";
            prep_warnings.push_back( msg.str() );
            continue;
          }
        }//if( background_sub )
        
        // If we are using SandiaDecay as our nuclear data source, we will use the energy of the
        //  gamma line, and not the fit-energy.  I'm a little torn on this, as the behavior could
        //  be unexpected to the user, either way.  But if we dont do this, and the energy
        //  calibration is slightly off for HPGe detectors, we will miss the match
        if( srcData == RelActCalcManual::PeakCsvInput::NucDataSrc::SandiaDecay )
        {
          if( (fabs(peak.m_energy - p->gammaParticleEnergy()) > (match_tol_sigma*(peak.m_fwhm/2.35482)))
             && (match_tol_sigma > 0.0) )
          {
            energy_cal_match_warning_energies.emplace( peak.m_energy, p->gammaParticleEnergy() );
          }
          
          peak.m_energy = p->gammaParticleEnergy();
        }//if( using SandiaDecay nuc data )
        
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
          
          const auto age_pos = nuclide_ages.find(nuc.nuclide);
          assert( age_pos != end(nuclide_ages) );
          
          if( age_pos != end(nuclide_ages) )
            nuc.age = age_pos->second;
          
          assert( nuc.age >= 0.0 );
          if( nuc.age < 0.0 )
            throw runtime_error( "Error finding age for " + nuc.nuclide->symbol );
          
          nuclides_to_match_to.push_back( nuc );
        }//if( pos != end(nuclides_to_match_to) ) / else
        
/*
 //Above we set the peak energy to the gamma line energy, so I believe
 //  `RelActCalcManual::PeakCsvInput::fill_in_nuclide_info(...)`, should math things correctly.
        if( match_tol_sigma == 0.0 )
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
        }//if( match_tol_sigma == 0.0 )
*/
        
        peak_infos.push_back( peak );
        
        has_U_or_Pu |= (p->parentNuclide()->atomicNumber == 92);
        has_U_or_Pu |= (p->parentNuclide()->atomicNumber == 94);
        
        unique_isotopes.insert( p->parentNuclide()->symbol );
      }//
    }//for( const PeakModel::PeakShrdPtr &p : *m_peakModel->peaks() )
    
    
    if( background_sub && !num_peaks_back_sub )
    {
      prep_warnings.push_back( "Subtraction of background peaks was selected, but no background"
                               " peaks matched to the selected foreground peaks." );
    }//if( user wanted to background subtract peaks, but no peaks matched up )
    
    if( has_U_or_Pu && (lowest_energy_peak < 122) )
    {
      prep_warnings.push_back( "The relative efficiency curve does not account for x-ray"
                              " absorption edges - using peaks under 120 keV for U or Pu problems"
                              " is not recommended.");
    }else if( lowest_energy_peak < 90 )
    {
      prep_warnings.push_back( "The relative efficiency curve does not account for x-ray"
                              " absorption edges of any potential shielding - please ensure the"
                              " peaks used do not span across absorption edges of any shielding.");
    }
    
    
    if( energy_cal_match_warning_energies.size() && (match_tol_sigma > 0.0) )
    {
      stringstream msg;
      const bool multiple = (energy_cal_match_warning_energies.size() > 1);
      msg << "The assigned gamma" << (multiple ? "s of {" : " of ");
      bool first = true;
      for( const auto &pp : energy_cal_match_warning_energies )
      {
        msg << (first ? "" : ", ") << PhysicalUnits::printCompact(pp.second, 4);
        first = false;
      }
      msg << (multiple ? "}" : "") << " keV are outside of the match tolerance with the peak"
      << (multiple ? "s of {" : " of ");
      first = true;
      for( const auto &pp : energy_cal_match_warning_energies )
      {
        msg << (first ? "" : ", ") << PhysicalUnits::printCompact(pp.first, 4);
        first = false;
      }
      msg << (multiple ? "}" : "")
          << "; note that the assigned gamma energy is used to compensate for nearby gammas,"
          << " and not the fit peak mean.";
      
      prep_warnings.push_back( msg.str() );
    }//if( energy_cal_match_warning_energies.size() )
    
    // Check that if we are using a specialized nuc data source, we actually have uranium in the
    //  problem
    switch( srcData )
    {
      case RelActCalcManual::PeakCsvInput::NucDataSrc::SandiaDecay:
        break;
        
      case RelActCalcManual::PeakCsvInput::NucDataSrc::Icrp107_U:
      case RelActCalcManual::PeakCsvInput::NucDataSrc::Lanl_U:
      case RelActCalcManual::PeakCsvInput::NucDataSrc::IcrpLanlGadras_U:
      {
        if( match_tol_sigma <= 0.0 )
          throw runtime_error( "You must have a non-zero match tolerance if you are not using InterSpecs default nuclear data." );
        
        size_t num_u_isos = 0;
        for( const auto &n : nuclides_to_match_to )
        {
          if( n.nuclide && (n.nuclide->atomicNumber == 92) )
            num_u_isos += 1;
        }
        
        if( !num_u_isos )
          throw runtime_error( "A uranium specialized data source was selected, but there is no uranium in this problem." );
        
        break;
      }//case( not using SandiaDecay for uranium )
        
      case RelActCalcManual::PeakCsvInput::NucDataSrc::Undefined:
        break;
    }//switch( srcData )
    
    if( peak_infos.empty() )
      throw runtime_error( "No peaks selected." );
    
    assert( !nuclides_to_match_to.empty() );
    if( nuclides_to_match_to.empty() )
      throw runtime_error( "No nuclides selected." );
    
    {// Begin code to fill in nuclide information
      vector<RelActCalcManual::PeakCsvInput::NucAndAge> isotopes;
      for( const auto &n : nuclides_to_match_to )
        isotopes.emplace_back( n.nuclide->symbol, n.age );
      
      // Make sure the match tolerance is ever so slightly above zero, for practical purposes
      const double tol = std::max( match_tol_sigma, 0.0001 );
      
      RelActCalcManual::PeakCsvInput::NucMatchResults matched_res
        = RelActCalcManual::PeakCsvInput::fill_in_nuclide_info( peak_infos, srcData,
                                                           {}, isotopes, tol, {} );
      
      if( matched_res.unused_isotopes.size() )
      {
        stringstream msg;
        msg << "Failed to match nuclide" << ((matched_res.unused_isotopes.size() > 1) ? "s" : "")
        << " to any peaks: ";
        for( size_t i = 0; i < matched_res.unused_isotopes.size(); ++i )
          msg << (i ? ", " : "") << matched_res.unused_isotopes[i];
        
        throw runtime_error( msg.str() );
      }//if( matched_res.unused_isotopes.size() )
      
      if( (srcData == RelActCalcManual::PeakCsvInput::NucDataSrc::SandiaDecay)
         && matched_res.peaks_not_matched.size() )
      {
        throw runtime_error( "logic error: not all input peaks were matched to a nuclide - even though they should have been." );
      }//if( we somehow didnt match a peak the user wanted to be used )
    
      peak_infos = matched_res.peaks_matched;
      //peak_infos.insert( end(peak_infos),
      //                begin(matched_res.peaks_not_matched),
      //                end(matched_res.peaks_not_matched) );
    
      //std::sort( begin(peak_infos), end(peak_infos),
      //        []( const RelActCalcManual::GenericPeakInfo &lhs, const RelActCalcManual::GenericPeakInfo &rhs ){
      //  return lhs.m_energy < rhs.m_energy;
      //});
    }// End code to fill in nuclide information
    
  
    const RelActCalc::RelEffEqnForm eqn_form = relEffEqnForm();
    const size_t eqn_order = relEffEqnOrder();
    
    
    RelEffSolution solution = solve_relative_efficiency( peak_infos, eqn_form, eqn_order );
    
    solution.m_warnings.insert(begin(solution.m_warnings), begin(prep_warnings), end(prep_warnings));
    
    m_currentSolution = make_shared<RelEffSolution>( solution );
    
    updateGuiWithResults();
  }catch( std::exception &e )
  {
    cout << "Error doing calc: " << e.what() << endl;
    
    //Set error to results TXT and show results TXT tab.
    
    string errormsg = "Error calculating relative activity/efficiency: " + string(e.what());
    WText *error = new WText( errormsg, m_results );
    error->setInline( false );
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
      m_chart->setData( vector<RelActCalcManual::GenericPeakInfo>{}, {}, "" );
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
    m_resultMenu->select( 0 );
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
  
  results_html << "<div class=\"ToolAlphaWarning\">"
  "Errors are statistical only, and have not been validated.<br />"
  "The tool is is in an alpha-preview state only."
  "</div>\n";
  
  
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


void RelActManualGui::matchToleranceChanged()
{
  m_renderFlags |= RenderActions::UpdateCalc;
  scheduleRender();
}//void matchToleranceChanged();


void RelActManualGui::addUncertChanged()
{
  m_renderFlags |= RenderActions::UpdateCalc;
  scheduleRender();
}

void RelActManualGui::backgroundSubtractChanged()
{
  m_renderFlags |= RenderActions::UpdateCalc;
  scheduleRender();
}


void RelActManualGui::handlePeaksChanged()
{
  m_renderFlags |= RenderActions::UpdateCalc;
  m_renderFlags |= RenderActions::UpdateNuclides;
  scheduleRender();
}



void RelActManualGui::updateNuclides()
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
  
  bool has_uranium = false;
  for( auto w : m_nuclidesDisp->children() )
  {
    ManRelEffNucDisp *rr = dynamic_cast<ManRelEffNucDisp *>(w);
    if( rr )
    {
      const bool isU = (rr->m_nuc->atomicNumber == 92);
      has_uranium |= isU;
      rr->setAgeHidden( isU ? !showAge : false );
    }else
    {
      assert( 0 );
    }
  }//for( auto w : m_nuclidesDisp->children() )
  
  m_nucDataSrcHolder->setHidden( !has_uranium );
}//void updateNuclides()


void RelActManualGui::updateSpectrumBasedOptions()
{
  InterSpec *interspec = InterSpec::instance();
  
  assert( interspec );
  
  if( !interspec )
    return;
  
  shared_ptr<const SpecMeas> back = interspec->measurment(SpecUtils::SpectrumType::Background);
  
  if( !back )
  {
    m_backgroundSubtractHolder->setHidden( true );
    return;
  }//if( !back )
  
  const set<int> &displayed = interspec->displayedSamples(SpecUtils::SpectrumType::Background);
  const auto peaks = back->peaks( displayed );
  
  if( !peaks || peaks->empty() )
  {
    m_backgroundSubtractHolder->setHidden( true );
    return;
  }//if( !peaks || peaks->empty() )
  
  // We have peaks - although they may not overlap
  m_backgroundSubtractHolder->setHidden( false );
}//void updateSpectrumBasedOptions();



void RelActManualGui::displayedSpectrumChanged()
{
  m_renderFlags |= RenderActions::UpdateCalc;
  m_renderFlags |= RenderActions::UpdateNuclides;
  m_renderFlags |= RenderActions::UpdateSpectrumOptions;
  
  scheduleRender();
}//void displayedSpectrumChanged()


RelActCalc::RelEffEqnForm RelActManualGui::relEffEqnForm() const
{
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
  
  return eqn_form;
}//RelEffEqnForm relEffEqnForm() const


size_t RelActManualGui::relEffEqnOrder() const
{
  const int orderIndex = m_relEffEqnOrder->currentIndex();
  if( (orderIndex < 0) || (orderIndex > 7) )
    throw runtime_error( "Invalid RelEffEqnOrder" );
  
  return static_cast<size_t>(orderIndex) + 1;
}//size_t relEffEqnOrder() const


RelActCalcManual::PeakCsvInput::NucDataSrc RelActManualGui::nucDataSrc() const
{
  bool validSrcData = false;
  const auto srcData = RelActCalcManual::PeakCsvInput::NucDataSrc( m_nucDataSrc->currentIndex() );
  
  switch( srcData )
  {
    case RelActCalcManual::PeakCsvInput::NucDataSrc::SandiaDecay:
    case RelActCalcManual::PeakCsvInput::NucDataSrc::Icrp107_U:
    case RelActCalcManual::PeakCsvInput::NucDataSrc::Lanl_U:
    case RelActCalcManual::PeakCsvInput::NucDataSrc::IcrpLanlGadras_U:
      validSrcData = true;
      break;
      
    case RelActCalcManual::PeakCsvInput::NucDataSrc::Undefined:
      validSrcData = false;
      break;
  }//switch( srcData )
  
  if( !validSrcData )
    throw runtime_error( "Invalid SourceData" );
  
  // Now check to see if there is any uranium in the problem, and if not, return SandiaDecay;
  //  Note, could call `m_nucDataSrcHolder->isHidden()`, but it may not be up to date until the
  //  render stage
  
  if( srcData != RelActCalcManual::PeakCsvInput::NucDataSrc::SandiaDecay )
  {
    PeakModel *peakModel = m_interspec->peakModel();
    const auto peaks = peakModel ? peakModel->peaks() : nullptr;
    if( peaks )
    {
      for( const auto &p : *peaks )
      {
        const auto parent = (p && p->parentNuclide()) ? p->parentNuclide() : nullptr;
        if( parent && p->useForManualRelEff() && (parent->atomicNumber == 92) )
          return srcData;
      }//
    }//if( peaks )
    
    // IF we made it here, there was no Uranium isotope, so use SandiaDecay
    return RelActCalcManual::PeakCsvInput::NucDataSrc::SandiaDecay;
  }//if( srcData != RelActCalcManual::PeakCsvInput::NucDataSrc::SandiaDecay )
  
  return srcData;
}//NucDataSrc nucDataSrc() const


shared_ptr<const RelActCalcManual::RelEffSolution> RelActManualGui::currentSolution()
{
  return m_currentSolution;
}


rapidxml::xml_node<char> *RelActManualGui::serialize( rapidxml::xml_node<char> *parent_node )
{
  rapidxml::xml_document<char> *doc = parent_node->document();
  assert( doc );
  
  const char *name = "RelActManualGui";
  rapidxml::xml_node<> * const base_node = doc->allocate_node( rapidxml::node_element, name );
  parent_node->append_node( base_node );
  
  //If you change the available options or formatting or whatever, increment the
  //  version field of the XML!
  const string versionstr = std::to_string(RelActManualGui::sm_xmlSerializationMajorVersion)
                            + "." + std::to_string(RelActManualGui::sm_xmlSerializationMinorVersion);
  const char *value = doc->allocate_string( versionstr.c_str() );
  rapidxml::xml_attribute<> *attr = doc->allocate_attribute( "version", value );
  base_node->append_attribute( attr );
  
  const RelActCalc::RelEffEqnForm eqn_form = relEffEqnForm();
  value = RelActCalc::to_str( eqn_form );
  rapidxml::xml_node<> *node = doc->allocate_node( rapidxml::node_element, "RelEffEqnForm", value );
  base_node->append_node( node );
  
  const size_t eqn_order = relEffEqnOrder();
  value = doc->allocate_string( std::to_string(eqn_order).c_str() );
  node = doc->allocate_node( rapidxml::node_element, "RelEffEqnOrder", value );
  base_node->append_node( node );
  
  const RelActCalcManual::PeakCsvInput::NucDataSrc srcData = nucDataSrc();
  
  value = RelActCalcManual::PeakCsvInput::to_str( srcData );
  node = doc->allocate_node( rapidxml::node_element, "NucDataSrc", value );
  base_node->append_node( node );
  
  const float match_tolerance = m_matchTolerance->value();
  const string match_tolerance_str = PhysicalUnits::printCompact(match_tolerance, 7);
  value = doc->allocate_string( match_tolerance_str.c_str() );
  node = doc->allocate_node( rapidxml::node_element, "MatchTolerance", value );
  base_node->append_node( node );
  
  value = RelActManualGui::to_str( AddUncert(m_addUncertainty->currentIndex()) );
  node = doc->allocate_node( rapidxml::node_element, "AddUncertainty", value );
  base_node->append_node( node );
  
  
  value = m_resultMenu->currentIndex() ? "1" : "0";
  node = doc->allocate_node( rapidxml::node_element, "ResultTabShowing", value );
  base_node->append_node( node );
  
  
  map<string,double> nuc_age_cache = m_nucAge;
  
  for( auto w : m_nuclidesDisp->children() )
  {
    const ManRelEffNucDisp *rr = dynamic_cast<const ManRelEffNucDisp *>(w);
    if( rr && rr->m_nuc )
      nuc_age_cache[rr->m_nuc->symbol] = rr->m_current_age;
  }//for( auto w : m_nuclidesDisp->children() )
  
  
  rapidxml::xml_node<> *nuc_ages_node = doc->allocate_node( rapidxml::node_element, "NuclideAges" );
  base_node->append_node( nuc_ages_node );
  for( const auto &n : nuc_age_cache )
  {
    rapidxml::xml_node<> *nuc_node = doc->allocate_node( rapidxml::node_element, "Nuclide" );
    nuc_ages_node->append_node( nuc_node );
    
    value = doc->allocate_string( n.first.c_str() );
    rapidxml::xml_node<> *name_node = doc->allocate_node( rapidxml::node_element, "Name", value );
    nuc_node->append_node(name_node);
    
    value = doc->allocate_string( PhysicalUnits::printCompact(n.second,8).c_str() );
    rapidxml::xml_node<> *age_node = doc->allocate_node( rapidxml::node_element, "Age", value );
    nuc_node->append_node(age_node);
  }
  
  return base_node;
}//serialize(...)


void RelActManualGui::deSerialize( const rapidxml::xml_node<char> *base_node )
{
  if( SpecUtils::xml_name_str(base_node) != "RelActManualGui" )
    throw runtime_error( "RelActManualGui::deSerialize: invalid base node passed in: '"
                        + SpecUtils::xml_name_str(base_node) + "'" );
    
  int version;
  const rapidxml::xml_attribute<char> *attr = XML_FIRST_ATTRIB(base_node, "version");
  if( !attr || !attr->value() || !(stringstream(attr->value()) >> version) )
    throw runtime_error( "Deserializing requires a version" );
  
  if( version != sm_xmlSerializationMajorVersion )
    throw runtime_error( "Invalid version of RelActManualGui XML" );
  
  m_nuclidesDisp->clear();
  //m_nucAge.clear();
  
  const rapidxml::xml_node<char> *RelEffEqnForm_node = XML_FIRST_NODE(base_node, "RelEffEqnForm");
  const rapidxml::xml_node<char> *RelEffEqnOrder_node = XML_FIRST_NODE(base_node, "RelEffEqnOrder");
  const rapidxml::xml_node<char> *NucDataSrc_node = XML_FIRST_NODE(base_node, "NucDataSrc");
  const rapidxml::xml_node<char> *MatchTolerance_node = XML_FIRST_NODE(base_node, "MatchTolerance");
  const rapidxml::xml_node<char> *AddUncertainty_node = XML_FIRST_NODE(base_node, "AddUncertainty");
  const rapidxml::xml_node<char> *ResultTabShowing_node = XML_FIRST_NODE(base_node, "ResultTabShowing");
  const rapidxml::xml_node<char> *NuclideAges_node = XML_FIRST_NODE(base_node, "NuclideAges");
  
  
  if( !RelEffEqnForm_node || !RelEffEqnOrder_node || !NucDataSrc_node || !MatchTolerance_node
     || !AddUncertainty_node || !ResultTabShowing_node || !NuclideAges_node )
    throw runtime_error( "RelActManualGui::deSerialize: missing required node" );
  
  bool got_eqn_form = false;
  RelActCalc::RelEffEqnForm eqn_form = RelActCalc::RelEffEqnForm::LnX;
  const string rel_eff_eqn_form_str = SpecUtils::xml_value_str( RelEffEqnForm_node );
  for( RelActCalc::RelEffEqnForm form = RelActCalc::RelEffEqnForm(0);
      form <= RelActCalc::RelEffEqnForm::FramEmpirical;
      form = RelActCalc::RelEffEqnForm(static_cast<int>(form)+1) )
  {
    const char * const test_form = RelActCalc::to_str( form );
    if( test_form == rel_eff_eqn_form_str )
    {
      got_eqn_form = true;
      eqn_form = form;
      break;
    }
  }//for( loop over RelEffEqnForm )
  
  if( !got_eqn_form )
    throw runtime_error( "RelActManualGui::deSerialize: '" + rel_eff_eqn_form_str + "' is invalid RelEffEqnForm" );
  
  
  size_t eqn_order;
  const string rel_eff_order_str = SpecUtils::xml_value_str(RelEffEqnOrder_node);
  if( !(stringstream(rel_eff_order_str) >> eqn_order)
     || (eqn_order <= 0)
     || (eqn_order > m_relEffEqnOrder->count() ) )
  {
    throw runtime_error( "RelActManualGui::deSerialize: '" + rel_eff_order_str + "' is invalid RelEffEqnOrder" );
  }
  
  bool got_data_src = false;
  auto data_src = RelActCalcManual::PeakCsvInput::NucDataSrc::Icrp107_U;
  const string data_src_str = SpecUtils::xml_value_str(NucDataSrc_node);
  for( auto src = RelActCalcManual::PeakCsvInput::NucDataSrc(0);
      src < RelActCalcManual::PeakCsvInput::NucDataSrc::Undefined;
      src = RelActCalcManual::PeakCsvInput::NucDataSrc(static_cast<int>(src) + 1) )
  {
    const char *val = RelActCalcManual::PeakCsvInput::to_str( src );
    if( val == data_src_str )
    {
      got_data_src = true;
      data_src = src;
      break;
    }
  }
  
  if( !got_data_src )
    throw runtime_error( "RelActManualGui::deSerialize: '" + data_src_str + "' is invalid NucDataSrc" );
  
  float match_tolerance;
  const string match_tol_str = SpecUtils::xml_value_str(MatchTolerance_node);
  if( !(stringstream( match_tol_str ) >> match_tolerance)
     || (match_tolerance < 0.0f)
     || (match_tolerance > 5.0f)
     || IsNan(match_tolerance))
  {
    throw runtime_error( "RelActManualGui::deSerialize: '" + match_tol_str + "' is invalid match tolerance" );
  }
  
  bool got_add_uncert = false;
  AddUncert add_uncert = AddUncert::FiftyPercent;
  const string add_uncert_str = SpecUtils::xml_value_str(AddUncertainty_node);
  for( AddUncert uncert = AddUncert(0);
      uncert <= AddUncert::NumAddUncert;
      uncert = AddUncert(static_cast<int>(uncert) + 1) )
  {
    const char *test_str = RelActManualGui::to_str( uncert );
    if( test_str == add_uncert_str )
    {
      got_add_uncert = true;
      add_uncert = uncert;
      break;
    }
  }//for( loop over AddUncert )
  
  if( !got_add_uncert )
    cerr << "RelActManualGui::deSerialize: failed to decode '"
         << add_uncert_str << "' into a AddUncert enum" << endl;
  
  int tab_showing = 0;
  const string tab_show_str = SpecUtils::xml_value_str(ResultTabShowing_node);
  if( !(stringstream(tab_show_str) >> tab_showing) )
    cerr << "RelActManualGui::deSerialize: invalid ResultTabShowing node '"
         << tab_show_str << "'" << endl;
  
  if( (tab_showing != 0) && (tab_showing != 1) )
    tab_showing = 0;
  
  XML_FOREACH_DAUGHTER( nuc_node, NuclideAges_node, "Nuclide" )
  {
    rapidxml::xml_node<> *name_node = XML_FIRST_NODE(nuc_node, "Name");
    rapidxml::xml_node<> *age_node = XML_FIRST_NODE(nuc_node, "Age");
    
    const string nuc = SpecUtils::xml_value_str(name_node);
    const string age_str = SpecUtils::xml_value_str(age_node);
    if( nuc.empty() || age_str.empty() )
    {
      cerr << "RelActManualGui::deSerialize: invalid Nuclide Name or Age '"
      << nuc << "', '" << age_str << "'" << endl;
      continue;
    }
    
    double age;
    if( !(stringstream(age_str) >> age) )
    {
      cerr << "RelActManualGui::deSerialize: Age '" << age_str << "'" << endl;
      continue;
    }
    
    m_nucAge[nuc] = age;
  }//for( loop over <Nuclide> nodes )
  
  
  // Now need to set state of widgets
  m_relEffEqnForm->setCurrentIndex( static_cast<int>(eqn_form) );
  m_relEffEqnOrder->setCurrentIndex( eqn_order - 1 );
  m_nucDataSrc->setCurrentIndex( static_cast<int>(data_src) );
  m_matchTolerance->setValue( match_tolerance );
  m_addUncertainty->setCurrentIndex( static_cast<int>(add_uncert) );
  m_resultMenu->select( tab_showing );
  
  // Schedult
  m_renderFlags |= RenderActions::UpdateCalc;
  m_renderFlags |= RenderActions::UpdateNuclides;
  scheduleRender();
}//deSerialize(...)


const char *RelActManualGui::to_str( const RelActManualGui::AddUncert val )
{
  switch( val )
  {
    case AddUncert::Unweighted:         return "Unweighted";
    case AddUncert::StatOnly:           return "StatOnly";
    case AddUncert::OnePercent:         return "OnePercent";
    case AddUncert::FivePercent:        return "FivePercent";
    case AddUncert::TenPercent:         return "TenPercent";
    case AddUncert::TwentyFivePercent:  return "TwentyFivePercent";
    case AddUncert::FiftyPercent:       return "FiftyPercent";
    case AddUncert::SeventyFivePercent: return "SeventyFivePercent";
    case AddUncert::OneHundredPercent:  return "OneHundredPercent";
    case AddUncert::NumAddUncert:       return "NumAddUncert";
  }//
  
  return "InvalidAddUncert";
}//to_str( const AddUncert val )
