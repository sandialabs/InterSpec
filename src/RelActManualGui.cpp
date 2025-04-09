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

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include <boost/math/distributions/chi_squared.hpp>

#include <Wt/WMenu>
#include <Wt/WLabel>
#include <Wt/WPanel>
#include <Wt/WTable>
#include <Wt/WLabel>
#include <Wt/WServer>
#include <Wt/WComboBox>
#include <Wt/WGroupBox>
#include <Wt/WLineEdit>
#include <Wt/WMenuItem>
#include <Wt/WResource>
#include <Wt/WTableRow>
#include <Wt/WCheckBox>
#include <Wt/WIOService>
#include <Wt/WTableCell>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/Http/Request>
#include <Wt/WApplication>
#include <Wt/Http/Response>
#include <Wt/WItemDelegate>
#include <Wt/WStackedWidget>
#include <Wt/WContainerWidget>
#include <Wt/WRegExpValidator>
#include <Wt/WSuggestionPopup>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/RelEffChart.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/SwitchCheckbox.h"
#include "InterSpec/RelActManualGui.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/RelEffShieldWidget.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeSelectionAids.h"
#include "InterSpec/PhysicalUnitsLocalized.h"

using namespace Wt;
using namespace std;

#if( ANDROID )
// Defined in target/android/android.cpp
extern void android_download_workaround( Wt::WResource *resource, std::string description );
#endif

const int RelActManualGui::sm_xmlSerializationMajorVersion = 0;
const int RelActManualGui::sm_xmlSerializationMinorVersion = 1;

namespace
{
struct DoWorkOnDestruct
{
  std::function<void()> m_worker;
  DoWorkOnDestruct( std::function<void()> &&worker ) : m_worker( std::move(worker) ){}
  ~DoWorkOnDestruct(){ if(m_worker) m_worker(); }
};//struct DoWorkOnDestruct


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
    setTakesUpdateLock( true );
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
        "\t<head><meta charset=\"utf-8\"><title>" 
            << WString::tr("ramrr-no-solution-title").toUTF8()
            << "</title></head>"
        "\t<body>"
        "\t\t" << WString::tr("ramrr-no-solution-text").toUTF8()
        << "\t</body>"
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
        for( const auto &r : solution->m_input.peaks )
        {
          if( fabs(p->mean() - r.m_mean) < 1.0 )
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
      
      auto background = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Background);
      const double back_sf = m_interspec->displayScaleFactor(SpecUtils::SpectrumType::Background);
      solution->print_html_report( response.out(), title, foreground, display_peaks, background, back_sf );
    }catch( std::exception &e )
    {
      cerr << "Error handling request for RelActManualReportResource: " << e.what() << endl;
      log("error") << "Error handling request for RelActManualReportResource: " << e.what();
      response.out() << "Error creating HTML file: " << e.what()
      << "\n\nPlease report to InterSpec@sandia.gov.";
      
      passMessage( WString::tr("ramrr-err-report").arg(e.what()), WarningWidget::WarningMsgHigh );
      
      response.setStatus(500);
      assert( 0 );
    }//try / catch
  }//void handleRequest(...)
};//class RelActManualReportResource


class ManRelEffNucDisp : public Wt::WPanel
{
public:
  const SandiaDecay::Nuclide * const m_nuc;
  const ReactionGamma::Reaction * const m_reaction;
  double m_current_age;
  const bool m_age_is_settable;
  Wt::WTable *m_nucContentTable;
  Wt::WLineEdit *m_age_edit_tmp;
  Wt::WTableRow *m_age_row;
  Wt::WCheckBox *m_decay_during_meas;
  
  Wt::Signal<> m_updated;
  
  
public:
  ManRelEffNucDisp( const SandiaDecay::Nuclide * const nuc,
                   const ReactionGamma::Reaction * const reaction,
                   double age, const float meas_time,
                   WContainerWidget *parent = nullptr )
  : WPanel( parent ),
   m_nuc( nuc ),
   m_reaction( reaction ),
   m_current_age( (nuc && (age < 0.0)) ? PeakDef::defaultDecayTime(nuc) : age ),
   m_age_is_settable( nuc ? !PeakDef::ageFitNotAllowed(nuc) : false ),
   m_nucContentTable( nullptr ),
   m_age_edit_tmp( nullptr ),
   m_age_row( nullptr ),
   m_decay_during_meas( nullptr ),
   m_updated( this )
  {
    assert( m_nuc || m_reaction );
    
    if( !m_nuc && !m_reaction )
      throw runtime_error( "ManRelEffNucDisp: null nuclide and reaction" );
    
    addStyleClass( "ManRelEffNucDisp" );
    
    setTitle( m_nuc ? m_nuc->symbol : m_reaction->name() );
    setCollapsible( true );
    setCollapsed( true );
    //setAnimation( { WAnimation::AnimationEffect::SlideInFromTop,
    //  WAnimation::TimingFunction::Linear, 250 } );
    
    if( m_nuc )
    {
      m_nucContentTable = new WTable();
      m_nucContentTable->addStyleClass( "NucInfoTable" );
      setCentralWidget( m_nucContentTable );
      
      WTableCell *cell = m_nucContentTable->elementAt(0, 0);
      WLabel *label = new WLabel( WString::tr("Age"), cell );
      
      m_age_row = m_nucContentTable->rowAt(0);
      
      cell = m_nucContentTable->elementAt(0, 1);
      
      string agestr;
      if( !nuc || nuc->decaysToStableChildren() )
      {
        m_current_age = 0.0;
        agestr = "0y";
      }else
      {
        m_current_age = PeakDef::defaultDecayTime( nuc, &agestr );
      }//if( decay to stable only ) / else
      
      if( m_age_is_settable )
      {
        m_age_edit_tmp = new WLineEdit( agestr, cell );
        
        WRegExpValidator *validator = new WRegExpValidator( PhysicalUnitsLocalized::timeDurationHalfLiveOptionalRegex(), m_age_edit_tmp );
        validator->setFlags(Wt::MatchCaseInsensitive);
        m_age_edit_tmp->setValidator(validator);
        label->setBuddy( m_age_edit_tmp );
        m_age_edit_tmp->addStyleClass( "AgeEdit" );
        
        m_age_edit_tmp->setAutoComplete( false );
        m_age_edit_tmp->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
        m_age_edit_tmp->setAttributeValue( "autocorrect", "off" );
        m_age_edit_tmp->setAttributeValue( "spellcheck", "off" );
#endif
        
        m_age_edit_tmp->changed().connect( this, &ManRelEffNucDisp::handleAgeChange );
        m_age_edit_tmp->blurred().connect( this, &ManRelEffNucDisp::handleAgeChange );
        m_age_edit_tmp->enterPressed().connect( this, &ManRelEffNucDisp::handleAgeChange );
      }else
      {
        label = new WLabel( agestr, cell );
        label->addStyleClass( "FixedAge" );
      }//if( m_age_is_settable ) / else
      
      
      cell = m_nucContentTable->elementAt(1, 0);
      //label = new WLabel( "Half Life", cell );
      label = new WLabel( WString("<span style=\"font-size: small;\">{1}</span>").arg(WString::tr("T1/2")), cell );
      
      cell = m_nucContentTable->elementAt(1, 1);
      WText *txt = new WText( PhysicalUnitsLocalized::printToBestTimeUnits(nuc->halfLife), cell );
      
      cell = m_nucContentTable->elementAt(2, 0);
      label = new WLabel( WString::tr("mrend-spec-act"), cell );
      
      cell = m_nucContentTable->elementAt(2, 1);
      const bool useCurie = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
      const double specificActivity = nuc->activityPerGram() / PhysicalUnits::gram;
      const string sa = PhysicalUnits::printToBestSpecificActivityUnits( specificActivity, 3, useCurie );
      txt = new WText( sa, cell );
      
      if( meas_time > 0.005*m_nuc->halfLife ) //0.005 times HL is arbitrary
      {
        showDecayDuringMeasurementCb();
        setDecayDuringMeasurement( (meas_time > 0.005*m_nuc->halfLife) ); //0.5% of HL is arbitrary
      }
      
      // We could maybe list which gammas are currently being used
    }else
    {
      string target;
      if( reaction->targetElement )
        target = reaction->targetElement->name;
      else if( reaction->targetNuclide )
        target = reaction->targetNuclide->symbol;
      
      WString infostr;
      switch( m_reaction->type )
      {
        case AlphaNeutron:            infostr = WString::tr("mrend-X(a,n)").arg(target); break;
        case NeutronAlpha:            infostr = WString::tr("mrend-X(n,a)").arg(target); break;
        case AlphaProton:             infostr = WString::tr("mrend-X(a,p)").arg(target); break;
        case NeutronCapture:          infostr = WString::tr("mrend-X(n,g)").arg(target); break;
        case NeutronInelasticScatter: infostr = WString::tr("mrend-X(n,n)").arg(target); break;
        case AlphaInelasticScatter:   infostr = WString::tr("mrend-X(a,a)").arg(target); break;
        case AnnihilationReaction:    infostr = WString::tr("mrend-annih");              break;
        case NumReactionType:         infostr = WString::tr("mrend-unknown-rxctn");      break;
      }//switch( m_reaction->type )
      
      WText *content = new WText( infostr );
      content->setInline( false );
      content->addStyleClass( "NucInfoTable" );
      setCentralWidget( content );
    }//if( m_nuc ) / else ( m_reaction )
  }//ManRelEffNucDisp(...)
  
  
  Wt::Signal<> &updated()
  {
    return m_updated;
  }
  
  void handleDecayDuringMeasurementChanged()
  {
    m_updated.emit();
  }
  
  void showDecayDuringMeasurementCb()
  {
    if( m_decay_during_meas || !m_nuc || !m_nucContentTable )
      return;
    
    WTableCell *cell = m_nucContentTable->elementAt(3, 0);
    cell->setColumnSpan( 2 );
    m_decay_during_meas = new WCheckBox( WString::tr("mrend-cb-decay-during-meas"), cell );
    m_decay_during_meas->addStyleClass( "CbNoLineBreak" );
    m_decay_during_meas->checked().connect( this, &ManRelEffNucDisp::handleDecayDuringMeasurementChanged );
    m_decay_during_meas->unChecked().connect( this, &ManRelEffNucDisp::handleDecayDuringMeasurementChanged );
    
    const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );
    
    HelpSystem::attachToolTipOn( m_decay_during_meas, WString::tr("mrend-tt-decay-during-meas"), 
                                showToolTips );
  }//void showDecayDuringMeasurementCb()
  
  
  void setDecayDuringMeasurement( const bool correct )
  {
    if( !m_nuc || !m_nucContentTable )
      return;
    
    const bool prev = (m_decay_during_meas && m_decay_during_meas->isChecked());
    if( prev == correct )
      return;
    
    showDecayDuringMeasurementCb();
    if( m_decay_during_meas )
      m_decay_during_meas->setChecked( correct );
  }//void setDecayDuringMeasurement( const bool correct )
  
  
  bool decayDuringMeasurement() const
  {
    return (m_nuc && m_decay_during_meas && m_decay_during_meas->isChecked());
  }
  
  void handleAgeChange()
  {
    assert( m_age_edit_tmp || m_reaction );
   
    if( m_reaction )
      return;
    
    assert( m_nuc );
    if( !m_nuc || !m_age_edit_tmp )
      return;
    
    try
    {
      double age = 0;
      const string agestr = m_age_edit_tmp->text().toUTF8();
      age = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( agestr, m_nuc->halfLife );
      if( age < 0 )
        throw runtime_error( "Negative age not allowed." );
      
      m_current_age = age;
    }catch( std::exception & )
    {
      if( m_current_age >= 0.0 )
      {
        string agestr = PhysicalUnitsLocalized::printToBestTimeUnits( m_current_age, 5 );
        m_age_edit_tmp->setText( WString::fromUTF8(agestr) );
      }else
      {
        string agestr;
        m_current_age = PeakDef::defaultDecayTime( m_nuc, &agestr );
        m_age_edit_tmp->setText( WString::fromUTF8(agestr) );
      }
    }//try / catch
    
    m_updated.emit();
  }//void handleAgeChange()
  
  void setAgeHidden( const bool hidden )
  {
    assert( m_age_row || m_reaction );
    
    if( m_age_row )
      m_age_row->setHidden( hidden );
  }
  
  void setAge( const double age )
  {
    assert( m_nuc );
    if( !m_nuc || !m_age_edit_tmp  || (age < 0) || IsInf(age) || IsNan(age) )
      return;
    
    const string agestr = PhysicalUnitsLocalized::printToBestTimeUnits( age );
    m_age_edit_tmp->setText( WString::fromUTF8(agestr) );
    m_current_age = age;
    setCollapsed( false );
  }//void setAge( const double age )
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
  m_physModelUseHoerl( nullptr ),
  m_physModelUseHoerlHolder( nullptr ),
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
  m_nucColumnTitle( nullptr ),
  m_physicalModelShields( nullptr ),
  m_selfAttenShield( nullptr ),
  m_extAttenShields( nullptr ),
  m_nucAge{},
  m_nucDecayCorrect{},
  m_resultMenu( nullptr ),
  m_chart( nullptr ),
  m_results( nullptr ),
  m_defaultDetector( nullptr )
{
  assert( m_interspec && m_peakModel );
  if( !m_interspec || !m_peakModel )
    throw runtime_error( "No interspec/peak model" );
  
  wApp->useStyleSheet( "InterSpec_resources/RelActManualGui.css" );
  
  addStyleClass( "EnergyCalTool RelActManualGui" );
    
  // If the widget gets to less than about 1145px wide, then Wt layout will start shrinking
  //  the columns, even though they are fixed or minimum sized.  When this happens if we dont
  //  set the vertical overflow to hidden, a useless horizantal scrollbar will show up on the
  //  entire bottom of the widget (but it only will ever scroll for like 5 px - probably just
  //  padding somewhere), even though each column will also get scrollbars or be squeezed.
  setOverflow( Overflow::OverflowHidden, Wt::Orientation::Horizontal );
  
  init();
}//RelActManualGui constructor
  
  

void RelActManualGui::init()
{
  m_interspec->useMessageResourceBundle( "RelActManualGui" );
  
  m_renderFlags |= RenderActions::UpdateCalc;
  m_renderFlags |= RenderActions::UpdateNuclides;
  scheduleRender();
  
  if( m_layout )
    delete m_layout;
  
  m_layout = new WGridLayout( this );
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", m_interspec );
  
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
  
  WText *header = new WText( WString::tr("ramg-options-label") );
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
  int row = optionsList->rowCount();
  WLabel *label = new WLabel( WString::tr("ramg-eqn-form-label"), optionsList->elementAt(row, 0) );
  
  m_relEffEqnForm = new WComboBox( optionsList->elementAt(row, 1) );
  m_relEffEqnForm->activated().connect( boost::bind(&RelActManualGui::relEffEqnFormChanged, this, true) );
  
  HelpSystem::attachToolTipOn( {optionsList->elementAt(row,0), optionsList->elementAt(row,1)},
                              WString::tr("ramg-tt-eqn-form"), showToolTips );
  
  
  // Will assume FramPhysicalModel is the highest
  static_assert( static_cast<int>(RelActCalc::RelEffEqnForm::FramPhysicalModel)
                 > static_cast<int>(RelActCalc::RelEffEqnForm::LnXLnY),
                "RelEffEqnForm was changed!"
  );
  
  
  for( int i = 0; i <= static_cast<int>(RelActCalc::RelEffEqnForm::FramPhysicalModel); ++i )
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

      case RelActCalc::RelEffEqnForm::FramPhysicalModel:
        //y = [(1 - exp(-mu_0 * x_0))/(mu_0 * x_0) * [exp(-mu_1 * x_1) * ...] * [Det Eff] * [(0.001*E)^b * c^(1000/E)]
        txt = "Physical";
        break;
    }//switch( eqn_form )
    
    m_relEffEqnForm->addItem( txt );
  }//for( loop over RelEffEqnForm )
  
  m_relEffEqnForm->setCurrentIndex( static_cast<int>(RelActCalc::RelEffEqnForm::LnX) );
  
  row = optionsList->rowCount();
  label = new WLabel( WString::tr("ramg-eqn-order-label"), optionsList->elementAt(row, 0) );
  
  m_relEffEqnOrder = new WComboBox( optionsList->elementAt(row, 1) );
  m_relEffEqnOrder->activated().connect( this, &RelActManualGui::relEffEqnOrderChanged );
  
  m_relEffEqnOrder->addItem( "0" );
  m_relEffEqnOrder->addItem( "1" );
  m_relEffEqnOrder->addItem( "2" );
  m_relEffEqnOrder->addItem( "3" );
  m_relEffEqnOrder->addItem( "4" );
  m_relEffEqnOrder->addItem( "5" );
  m_relEffEqnOrder->addItem( "6" );
  m_relEffEqnOrder->setCurrentIndex( 3 );
  
  HelpSystem::attachToolTipOn( {optionsList->elementAt(row, 0),optionsList->elementAt(row, 1)},
                              WString::tr("ramg-tt-eqn-order"), showToolTips );
  
  m_relEffEqnOrderHolder = optionsList->rowAt(row);

  row = optionsList->rowCount();
  m_physModelUseHoerl = new WCheckBox( WString::tr("ramg-phys-model-use-hoerl"), optionsList->elementAt(row, 0) );
  m_physModelUseHoerl->setChecked( true );
  m_physModelUseHoerl->checked().connect( this, &RelActManualGui::physModelUseHoerlChanged );
  m_physModelUseHoerl->unChecked().connect( this, &RelActManualGui::physModelUseHoerlChanged );
  optionsList->elementAt(row, 0)->setColumnSpan( 2 );
  HelpSystem::attachToolTipOn( {optionsList->elementAt(row, 0)},
                              WString::tr("ramg-tt-phys-model-use-hoerl"), showToolTips );
  m_physModelUseHoerl->addStyleClass( "PhysModelUseHoerl" );
  m_physModelUseHoerlHolder = optionsList->rowAt(row);

  row = optionsList->rowCount();
  label = new WLabel( WString::tr("ramg-yield-info-label"), optionsList->elementAt(row, 0) );
  m_nucDataSrc = new WComboBox( optionsList->elementAt(row, 1) );
  label->setBuddy( m_nucDataSrc );
  m_nucDataSrc->activated().connect( this, &RelActManualGui::nucDataSrcChanged );
  
  HelpSystem::attachToolTipOn( {optionsList->elementAt(row, 0),optionsList->elementAt(row, 1)},
                              WString::tr("ramg-tt-data-src"), showToolTips );
  
  
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
  
  m_nucDataSrcHolder = optionsList->rowAt(row);

  row = optionsList->rowCount();
  label = new WLabel( WString::tr("ramg-match-tol-label"), optionsList->elementAt(row, 0) ); //(FWHM)
  m_matchTolerance = new NativeFloatSpinBox( optionsList->elementAt(row, 1) );
  label->setBuddy( m_matchTolerance );
  m_matchTolerance->setSpinnerHidden();
  m_matchTolerance->setWidth( 35 );
  m_matchTolerance->setRange( 0, 5 );
  m_matchTolerance->setValue( 0.5 ); //Other places we use 1.25/2.355 = 0.530786
  label = new WLabel( WString("&nbsp;{1}").arg(WString::tr("FWHM")), optionsList->elementAt(row, 1) );
  label->setBuddy( m_matchTolerance );
  m_matchTolerance->valueChanged().connect( this, &RelActManualGui::matchToleranceChanged );
  
  
  HelpSystem::attachToolTipOn( {optionsList->elementAt(row, 0),optionsList->elementAt(row, 1)},
                              WString::tr("ramg-tt-match-tol"), showToolTips );
  
  row = optionsList->rowCount();
  label = new WLabel( WString::tr("ramg-add-uncert-label"), optionsList->elementAt(row, 0) );
  m_addUncertainty = new WComboBox( optionsList->elementAt(row, 1) );
  label->setBuddy( m_addUncertainty );
  
  HelpSystem::attachToolTipOn( {optionsList->elementAt(row, 0),optionsList->elementAt(row, 1)},
                              WString::tr("ramg-tt-add-uncert"), showToolTips );
  
  m_addUncertainty->activated().connect( this, &RelActManualGui::addUncertChanged );
  
  for( AddUncert i = AddUncert(0); i < AddUncert::NumAddUncert; i = AddUncert(static_cast<int>(i) + 1) )
  {
    WString uncert_txt;
    switch( i )
    {
      case AddUncert::Unweighted:         uncert_txt = WString::tr("ramg-unweighted-label"); break;
      case AddUncert::StatOnly:           uncert_txt = WString::tr("ramg-stat-only-label");  break;
      case AddUncert::OnePercent:         uncert_txt = WString::fromUTF8("1%");              break;
      case AddUncert::FivePercent:        uncert_txt = WString::fromUTF8("5%");              break;
      case AddUncert::TenPercent:         uncert_txt = WString::fromUTF8("10%");             break;
      case AddUncert::TwentyFivePercent:  uncert_txt = WString::fromUTF8("25%");             break;
      case AddUncert::FiftyPercent:       uncert_txt = WString::fromUTF8("50%");             break;
      case AddUncert::SeventyFivePercent: uncert_txt = WString::fromUTF8("75%");             break;
      case AddUncert::OneHundredPercent:  uncert_txt = WString::fromUTF8("100%");            break;
      case AddUncert::NumAddUncert:       assert(0);                 break;
    }//switch( i )
    
    m_addUncertainty->addItem( uncert_txt );
  }//for( loop over AddUncert )
  
  m_addUncertainty->setCurrentIndex( static_cast<int>(AddUncert::StatOnly) );
  
  row = optionsList->rowCount();
  m_backgroundSubtract = new WCheckBox( WString::tr("ramg-back-sub-cb"), optionsList->elementAt(row, 0) );
  m_backgroundSubtract->addStyleClass( "BackSub CbNoLineBreak" );
  optionsList->elementAt(row, 0)->setColumnSpan( 2 );
  m_backgroundSubtractHolder = optionsList->rowAt(row);

  m_backgroundSubtract->checked().connect( this, &RelActManualGui::backgroundSubtractChanged );
  m_backgroundSubtract->unChecked().connect( this, &RelActManualGui::backgroundSubtractChanged );

  
  WContainerWidget *btndiv = new WContainerWidget();
  btndiv->addStyleClass( "BtmBtnDiv" );
  collayout->addWidget( btndiv, 2, 0 );
  
  auto helpBtn = new WContainerWidget( btndiv );
  helpBtn->addStyleClass( "Wt-icon ContentHelpBtn" );
  helpBtn->clicked().connect( boost::bind( &HelpSystem::createHelpWindow, "rel-act-manual" ) );
  
  m_htmlResource = new RelActManualReportResource( this, m_interspec, btndiv );
  
#if( BUILD_AS_OSX_APP || IOS )
  m_downloadHtmlReport = new WAnchor( WLink(m_htmlResource), btndiv );
  m_downloadHtmlReport->setTarget( AnchorTarget::TargetNewWindow );
  m_downloadHtmlReport->setStyleClass( "LinkBtn DownloadLink CALp" );
#else
  m_downloadHtmlReport = new WPushButton( btndiv );
  m_downloadHtmlReport->setIcon( "InterSpec_resources/images/download_small.svg" );
  m_downloadHtmlReport->setLink( WLink( m_htmlResource ) );
  m_downloadHtmlReport->setLinkTarget( Wt::TargetNewWindow );
  m_downloadHtmlReport->setStyleClass( "LinkBtn DownloadBtn CALp" );
 
#if( ANDROID )
  // Using hacked saving to temporary file in Android, instead of via network download of file.
  m_downloadHtmlReport->clicked().connect( std::bind([this](){
    android_download_workaround( m_htmlResource, "rel_eff.html");
  }) );
#endif //ANDROID
#endif
  
  m_downloadHtmlReport->setText( WString::tr("ramg-html-export-label") );

  WContainerWidget *nucCol = new WContainerWidget();
  nucCol->addStyleClass( "ToolTabTitledColumn RelActNucCol" );
  
  collayout = new WGridLayout( nucCol );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );

  header = new WText( WString::tr("ramg-nucs-label") );
  header->addStyleClass( "ToolTabColumnTitle" );
  collayout->addWidget( header, 0, 0 );
  m_nucColumnTitle = header;
  
  m_nuclidesDisp = new WContainerWidget();
  m_nuclidesDisp->addStyleClass( "ToolTabTitledColumnContent" );
  collayout->addWidget( m_nuclidesDisp, 1, 0 );
  collayout->setRowStretch( 1, 1 );
  m_layout->addWidget( nucCol, 0, 1 );


  m_physicalModelShields = new WContainerWidget( m_nuclidesDisp );
  m_physicalModelShields->addStyleClass( "PhysicalModelShields" );
  m_selfAttenShield = nullptr; //We wont create it until we need it
  m_extAttenShields = new WContainerWidget( m_physicalModelShields );
  m_physicalModelShields->hide();


  // Create the "Peaks to Use" table
  m_peakTableColumn = new WContainerWidget();
  m_peakTableColumn->addStyleClass( "ToolTabTitledColumn PeakTableCol" );
  m_layout->addWidget( m_peakTableColumn, 0, 2 );
  //m_layout->setColumnStretch( 2, 1 );
  
  
  collayout = new WGridLayout( m_peakTableColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );
  
  header = new WText( WString::tr("ramg-peaks-to-use-label") );
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
  
  header = new WText( WString::tr("ramg-results-label") );
  header->addStyleClass( "ToolTabColumnTitle" );
  collayout->addWidget( header, 0, 0 );
  
  
  
  WContainerWidget *resultContent = new WContainerWidget();
  resultContent->addStyleClass( "ToolTabTitledColumnContent ResultColumnContent" );
  resultContent->setMinimumSize(350, WLength::Auto );
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
  WMenuItem *item = new WMenuItem( WString::tr("ramg-mi-results"), m_results );
  m_resultMenu->addItem( item );
  
  // When outside the link area is clicked, the item doesnt get selected, so we'll work around this.
  item->clicked().connect( std::bind([this,item](){
    m_resultMenu->select( item );
    item->triggered().emit( item );
  }) );
  // We wont have an undo/redo for changing results tabs, because sometimes we programmatically
  //  set to a specific tab, like if there is an error, and this can create an infinite cycle
  //  that will block previous undo/redo's
  //item->triggered().connect( this, &RelActManualGui::resultTabChanged );
  
  m_chart = new RelEffChart();
  item = new WMenuItem( WString::tr("ramg-mi-chart"), m_chart );
  m_resultMenu->addItem( item );
  item->clicked().connect( std::bind([this,item](){
    m_resultMenu->select( item );
    item->triggered().emit( item );
  }) );
  //item->triggered().connect( this, &RelActManualGui::resultTabChanged );
  
  displayedSpectrumChanged();
}//void init()


shared_ptr<const RelActManualGui::GuiState> RelActManualGui::getGuiState() const
{
  try
  {
    shared_ptr<GuiState> state = make_shared<GuiState>();
    
    state->m_relEffEqnFormIndex = relEffEqnForm();
    state->m_relEffEqnOrderIndex = static_cast<int>( relEffEqnOrder() );
    state->m_physModelUseHoerl = m_physModelUseHoerl->isChecked();
    state->m_nucDataSrcIndex = nucDataSrc();
    state->m_matchToleranceValue = m_matchTolerance->value();
    state->m_addUncertIndex = RelActManualGui::AddUncert(m_addUncertainty->currentIndex());
    state->m_resultTab = m_resultMenu->currentIndex();
    state->m_backgroundSubtract = (!m_backgroundSubtractHolder->isHidden()
                                   && m_backgroundSubtract->isChecked());
    state->nucAge = m_nucAge;
    state->nucDecayCorrect = m_nucDecayCorrect;
    
    for( auto w : m_nuclidesDisp->children() )
    {
      const ManRelEffNucDisp *rr = dynamic_cast<const ManRelEffNucDisp *>(w);
      if( rr && rr->m_nuc )
      {
        //state->m_nucAgesAndDecayCorrect.emplace_back( rr->m_nuc->symbol, rr->m_current_age,
        //                                                rr->decayDuringMeasurement() );
        
        state->nucAge[rr->m_nuc->symbol] = rr->m_current_age;
        
        if( rr->m_decay_during_meas )
          state->nucDecayCorrect[rr->m_nuc->symbol] = rr->m_decay_during_meas->isChecked();
        else
          state->nucDecayCorrect.erase( rr->m_nuc->symbol );
      }
    }//for( auto w : m_nuclidesDisp->children() )
    
    // We could store the state of the Physical model shields, or not
    if( state->m_relEffEqnFormIndex == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      if( m_selfAttenShield && m_selfAttenShield->nonEmpty() )
        state->m_selfAttenShield = m_selfAttenShield->state();
      
      // We'll store the state of the Physical model shields in the XML, even if they are not used right now
      auto ext_kids = m_extAttenShields ? m_extAttenShields->children() : vector<WWidget *>();
      for( auto w : ext_kids )
      {
        RelEffShieldWidget *rr = dynamic_cast<RelEffShieldWidget *>(w);
        if( rr && rr->nonEmpty() )
          state->m_externalShields.push_back( rr->state() );
      }//for( auto w : m_extAttenShields->children() )
    }//if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    
    return state;
  }catch( std::exception &e )
  {
    cerr << "getGuiState(): Unexpected exception: " << e.what() << endl;
    assert( 0 );
  }
  
  return nullptr;
}//std::shared_ptr<const GuiState> getGuiState() const


void RelActManualGui::setGuiState( const RelActManualGui::GuiState &state )
{
  assert( !m_renderFlags.testFlag(RenderActions::AddUndoRedoStep) );
  
  bool updateCalc = false;
  const auto eqn_form = RelActCalc::RelEffEqnForm( state.m_relEffEqnFormIndex );

  if( m_relEffEqnForm->currentIndex() != static_cast<int>(state.m_relEffEqnFormIndex) )
  {
    updateCalc = true;
    m_relEffEqnForm->setCurrentIndex( static_cast<int>(state.m_relEffEqnFormIndex) );
    relEffEqnFormChanged( false );
  }
  
  assert( !m_renderFlags.testFlag(RenderActions::AddUndoRedoStep) );
  
  if( m_relEffEqnOrder->currentIndex() != state.m_relEffEqnOrderIndex )
  {
    updateCalc = true;
    m_relEffEqnOrder->setCurrentIndex( state.m_relEffEqnOrderIndex );
  }
  
  assert( !m_renderFlags.testFlag(RenderActions::AddUndoRedoStep) );
  
  if( m_physModelUseHoerl->isChecked() != state.m_physModelUseHoerl )
  {
    updateCalc |= (eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel);
    m_physModelUseHoerl->setChecked( state.m_physModelUseHoerl );
  }

  assert( !m_renderFlags.testFlag(RenderActions::AddUndoRedoStep) );
  
  if( m_nucDataSrc->currentIndex() != static_cast<int>(state.m_nucDataSrcIndex) )
  {
    updateCalc = true;
    m_nucDataSrc->setCurrentIndex( static_cast<int>(state.m_nucDataSrcIndex) );
  }

  assert( !m_renderFlags.testFlag(RenderActions::AddUndoRedoStep) );
  
  if( m_matchTolerance->value() != state.m_matchToleranceValue )
  {
    updateCalc = true;
    m_matchTolerance->setValue( state.m_matchToleranceValue );
  }
  
  assert( !m_renderFlags.testFlag(RenderActions::AddUndoRedoStep) );
  
  if( m_addUncertainty->currentIndex() != static_cast<int>(state.m_addUncertIndex) )
  {
    updateCalc = true;
    m_addUncertainty->setCurrentIndex( static_cast<int>(state.m_addUncertIndex) );
  }
  
  assert( !m_renderFlags.testFlag(RenderActions::AddUndoRedoStep) );
  
  if( m_backgroundSubtract->isChecked() != state.m_backgroundSubtract )
  {
    updateCalc = true;
    m_backgroundSubtract->setChecked( state.m_backgroundSubtract );
  }
  
  assert( !m_renderFlags.testFlag(RenderActions::AddUndoRedoStep) );
  
  if( m_resultMenu->currentIndex() != state.m_resultTab )
    m_resultMenu->select( state.m_resultTab );
  
  assert( !m_renderFlags.testFlag(RenderActions::AddUndoRedoStep) );
  
  m_nucAge = state.nucAge;
  m_nucDecayCorrect = state.nucDecayCorrect;
  
  for( auto w : m_nuclidesDisp->children() )
  {
    ManRelEffNucDisp *rr = dynamic_cast<ManRelEffNucDisp *>(w);
    if( !rr || !rr->m_nuc )
      continue;
    
    const auto age_pos = state.nucAge.find( rr->m_nuc->symbol );
    if( (age_pos != end(state.nucAge)) && (age_pos->second != rr->m_current_age) )
    {
      rr->setAge( age_pos->second );
      m_renderFlags |= RenderActions::UpdateNuclides;
      updateCalc = true;
    }//
    
    const bool decayCorr = rr->decayDuringMeasurement();
    const auto correct_pos = state.nucDecayCorrect.find( rr->m_nuc->symbol );
    if( (correct_pos != end(state.nucDecayCorrect)) && (decayCorr != correct_pos->second) )
    {
      rr->setDecayDuringMeasurement( correct_pos->second );
      m_renderFlags |= RenderActions::UpdateNuclides;
      updateCalc = true;
    }
  }//for( auto w : m_nuclidesDisp->children() )
  
  assert( !m_renderFlags.testFlag(RenderActions::AddUndoRedoStep) );
  
  if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    updateCalc = true;
    initPhysicalModelAttenShieldWidgets();
    assert( m_selfAttenShield );
    
    if( state.m_selfAttenShield )
    {
      if( m_selfAttenShield )
        m_selfAttenShield->setState( *state.m_selfAttenShield );
    }else
    {
      if( m_selfAttenShield )
        m_selfAttenShield->resetState();
    }//if( state.m_selfAttenShield )

    assert( !m_renderFlags.testFlag(RenderActions::AddUndoRedoStep) );
    
    assert( m_extAttenShields );
    
    vector<RelEffShieldWidget *> existing_external;
    for( const auto &s : m_extAttenShields->children() )
    {
      RelEffShieldWidget *rr = dynamic_cast<RelEffShieldWidget *>(s);
      if( rr )
        existing_external.push_back( rr );
    }
    
    while( !existing_external.empty()
          && (existing_external.size() > state.m_externalShields.size()) )
    {
      RelEffShieldWidget *rr = existing_external.back();
      existing_external.resize( existing_external.size() - 1 );
    }
    assert( (existing_external.size() == 1) || (existing_external.size() <= state.m_externalShields.size()) );
    
    for( size_t i = 0; i < state.m_externalShields.size(); ++i )
    {
      const unique_ptr<RelEffShieldState> &shield_state = state.m_externalShields[i];
      assert( shield_state );
      
      if( i >= existing_external.size() )
      {
        RelEffShieldWidget *rr = new RelEffShieldWidget( RelEffShieldWidget::ShieldType::ExternalAtten, m_extAttenShields );
        rr->changed().connect( this, &RelActManualGui::handlePhysicalModelShieldChanged );
        existing_external.push_back( rr );
      }
      
      assert( i < existing_external.size() );
      
      if( shield_state )
        existing_external.at(i)->setState( *shield_state );
    }//for( size_t i = 0; i < state.m_externalShields.size(); ++i )
    
    if( state.m_externalShields.empty() && !existing_external.empty() )
    {
      assert( existing_external.size() == 1 );
      existing_external.front()->resetState();
    }
  
    assert( !m_renderFlags.testFlag(RenderActions::AddUndoRedoStep) );
    
    m_physModelUseHoerlHolder->setHidden( false );
    m_relEffEqnOrderHolder->setHidden( true );
  }else
  {
    m_physModelUseHoerlHolder->setHidden( true );
    m_relEffEqnOrderHolder->setHidden( false );
  }//if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel ) / else

  if( updateCalc )
  {
    m_renderFlags |= RenderActions::UpdateCalc;
    scheduleRender();
  }
  
  assert( !m_renderFlags.testFlag(RenderActions::AddUndoRedoStep) );
}//void setGuiState( const GuiState &state )


void RelActManualGui::addUndoRedoStep( const shared_ptr<const RelActManualGui::GuiState> &state )
{
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( !undoRedo || !undoRedo->canAddUndoRedoNow() || !state || !m_currentGuiState )
    return;
  
  if( (*state) == (*m_currentGuiState) )
    return;
  
  const shared_ptr<const GuiState> prev_state = m_currentGuiState;
  
  auto undo = [prev_state](){
    InterSpec *viewer = InterSpec::instance();
    RelActManualGui *gui = viewer->createRelActManualWidget();
    if( gui && prev_state )
      gui->setGuiState( *prev_state );
  };//undo
    
  auto redo = [state](){
    InterSpec *viewer = InterSpec::instance();
    RelActManualGui *gui = viewer->createRelActManualWidget();
    if( gui && state )
      gui->setGuiState( *state );
  };//redo
    
  undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Isotopics from peaks change." );
}//void addUndoRedoStep( const std::shared_ptr<const GuiState> &state )


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
  
  
  const shared_ptr<const GuiState> current_gui_state = getGuiState();
  if( current_gui_state && m_renderFlags.testFlag(RelActManualGui::RenderActions::AddUndoRedoStep ) )
  {
    addUndoRedoStep( current_gui_state );
  }
  if( current_gui_state )
    m_currentGuiState = current_gui_state;
  
  
  m_renderFlags = 0;
  
  WContainerWidget::render( flags );
}//render( Wt::WFlags<Wt::RenderFlag> )


RelActManualGui::RelActCalcRawInput RelActManualGui::get_raw_info_for_calc_input()
{
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  if( !viewer )
    throw runtime_error( "Not in GUI thread???" );

  assert( m_peakModel );
  if( !m_peakModel )
    throw runtime_error( "No peak model" );
  
  RelActCalcRawInput setup_input;
  
  setup_input.state = getGuiState();
  assert( setup_input.state );
  if( !setup_input.state )
    throw runtime_error( "RelActManualGui::get_raw_info_for_calc_input: GUI is in invalid state???" );
  
  // PeakModel stores all PeakDefs as const, so I think it is fine to just make a copy of
  //  the deque, and then it should be thread-safe to access the PeakDef objects from off
  //  of the main thread.
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = m_peakModel->peaks();
  if( !peaks )
    throw runtime_error( "No foreground peaks not set/" );
  
  setup_input.peaks = *peaks;
  
  setup_input.fore_spec = viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  setup_input.back_spec = viewer->displayedHistogram(SpecUtils::SpectrumType::Background);
  
  const float background_live_time = setup_input.back_spec ? setup_input.back_spec->live_time() : 1.0f;
  
  const std::shared_ptr<const SpecMeas> back_meas = viewer->measurment(SpecUtils::SpectrumType::Background);
  
  
  if( setup_input.state->m_backgroundSubtract
     && setup_input.back_spec
     && back_meas
     && (background_live_time > 0.0) )
  {
    const auto &displayed = viewer->displayedSamples(SpecUtils::SpectrumType::Background);
    shared_ptr<const deque<shared_ptr<const PeakDef>>> backpeaks = back_meas->peaks( displayed );
    if( backpeaks )
      setup_input.background_peaks = *backpeaks;
  }//if( background_sub && back_spec && back_meas )
  
  const auto fore = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  if( fore )
    setup_input.detector = fore->detector();
  
  if( setup_input.state->m_relEffEqnFormIndex == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    if( !setup_input.detector )
      setup_input.detector = m_defaultDetector;

    if( !setup_input.detector )
    {
      const string datadir = InterSpec::staticDataDirectory();
      const string basename = SpecUtils::append_path( datadir, "GenericGadrasDetectors" );
      const string detdir = SpecUtils::append_path( basename, "HPGe 40%" );
      
      try
      {
        shared_ptr<DetectorPeakResponse> det = make_shared<DetectorPeakResponse>();
        det->fromGadrasDirectory( detdir );
        det->setFwhmCoefficients( {}, DetectorPeakResponse::ResolutionFnctForm::kNumResolutionFnctForm );
        setup_input.detector = det;
        m_defaultDetector = det;
      }catch( std::exception &e )
      {
        cerr << "Failed to load default Gadras detector: " << e.what() << endl;
      }
    }//if( !phys_model_detector )
  }//if( FramPhysicalModel )
  
  return setup_input;
}//RelActCalcRawInput get_raw_info_for_calc_input() const



void RelActManualGui::prepare_calc_input( const RelActCalcRawInput &setup_input,
                                         MaterialDB *materialDB,
                                         RelActCalcManual::RelEffInput &setup_output )
{
  const shared_ptr<const GuiState> &state = setup_input.state;
  if( !state )
    throw runtime_error( "Invalid GUI state" );
  
  using namespace RelActCalcManual;
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  
  const shared_ptr<const SpecUtils::Measurement> &fore_spec = setup_input.fore_spec;
  const shared_ptr<const SpecUtils::Measurement> &back_spec = setup_input.back_spec;
  
  
  const deque<shared_ptr<const PeakDef>> &peaks = setup_input.peaks;
  const deque<shared_ptr<const PeakDef>> &background_peaks = setup_input.background_peaks;
  
  vector<GenericPeakInfo> peak_infos;
  
  const double back_sub_nsigma_near = 1.0; // Fairly arbitrary.  TODO: have this be a user settable value?
  const float foreground_live_time = fore_spec ? fore_spec->live_time() : 1.0f;
  const float background_live_time = back_spec ? back_spec->live_time() : 1.0f;
  const bool background_sub = state->m_backgroundSubtract;
  
  double addUncert = -2.0;
  
  const AddUncert addUncertType = state->m_addUncertIndex;
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
    
  const double match_tol_sigma = 2.35482 * state->m_matchToleranceValue;
  
  const RelActCalcManual::PeakCsvInput::NucDataSrc srcData = state->m_nucDataSrcIndex;
  
  map<const SandiaDecay::Nuclide *,double> nuclide_ages;
  for( const auto &kv : state->nucAge )
  {
    const SandiaDecay::Nuclide * const nuc = db->nuclide(kv.first);
    if( nuc )
      nuclide_ages[nuc] = kv.second;
  }
  
  map<const SandiaDecay::Nuclide *,bool> decay_correct_during_meas;
  for( const auto &kv : state->nucDecayCorrect )
  {
    const SandiaDecay::Nuclide * const nuc = db->nuclide(kv.first);
    if( nuc )
      decay_correct_during_meas[nuc] = kv.second;
  }
  
  bool has_reaction = false;
  set<pair<float,float>> energy_cal_match_warning_energies;
  
  vector<SandiaDecayNuc> nuclides_to_match_to;
  
  size_t num_peaks_back_sub = 0;
  double lowest_energy_peak = 3000;
  bool has_U_or_Pu = false;
  set<string> unique_isotopes;
  set<int> src_atomic_numbers;
  for( const PeakModel::PeakShrdPtr &p : peaks )
  {
    if( p && (p->parentNuclide() || p->reaction()) && p->useForManualRelEff() )
    {
      GenericPeakInfo peak;
      peak.m_mean = peak.m_energy = p->mean();
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
        }//for( const PeakDef &peak : backPeaks )
        
        if( back_counts > 0.0 )
        {
          num_peaks_back_sub += 1;
          peak.m_counts -= back_counts;
          peak.m_counts_uncert = sqrt( peak.m_counts_uncert*peak.m_counts_uncert + back_uncert_2 );
        }
        
        if( peak.m_counts <= 0.0 )
        {
          char buffer[32];
          snprintf( buffer, sizeof(buffer), "%.2f", peak.m_mean );
          
          setup_output.prep_warnings.push_back( WString::tr("ramg-back-sub-neg").arg( buffer ).toUTF8() );
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
          energy_cal_match_warning_energies.emplace( static_cast<float>(peak.m_energy), p->gammaParticleEnergy() );
        }
        
        peak.m_energy = p->gammaParticleEnergy();
      }//if( using SandiaDecay nuc data )
      
      SandiaDecayNuc nuc;
      bool hadNuclide = false;
      for( const auto &existing : nuclides_to_match_to )
      {
        if( (p->parentNuclide() && (existing.nuclide == p->parentNuclide()))
           || (p->reaction() && (existing.reaction == p->reaction())) )
        {
          nuc = existing;
          hadNuclide = true;
          break;
        }
      }//for( const auto &existing : nuclides_to_match_to )
      
      if( !hadNuclide )
      {
        nuc.age = -1.0;
        if( p->parentNuclide() )
        {
          nuc.nuclide = p->parentNuclide();
          const auto age_pos = nuclide_ages.find(nuc.nuclide);
          assert( age_pos != end(nuclide_ages) );
        
          if( age_pos != end(nuclide_ages) )
            nuc.age = age_pos->second;
        
          assert( nuc.age >= 0.0 );
          if( nuc.age < 0.0 )
            throw runtime_error( "Error finding age for " + nuc.nuclide->symbol );
          
          const auto decay_corr_pos = decay_correct_during_meas.find(nuc.nuclide);
          
          // This next assert will trigger sometime; if a nuclide display doesnt have the decay
          //  during option (like long half-life nuclides), then that nuclide entry will be
          //  removed from `state->nucDecayCorrect` (because we dont have user preference info about
          //  it - or something)
          //assert( decay_corr_pos != end(decay_correct_during_meas) );
          
          if( decay_corr_pos != end(decay_correct_during_meas) )
            nuc.correct_for_decay_during_meas = decay_corr_pos->second;
        }else
        {
          has_reaction = true;
          nuc.reaction = p->reaction();
        }
        
        assert( nuc.nuclide || nuc.reaction );
        
        nuclides_to_match_to.push_back( nuc );
      }//if( !hadNuclide )
      
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
      
      if( p->parentNuclide() )
      {
        has_U_or_Pu |= (p->parentNuclide()->atomicNumber == 92);
        has_U_or_Pu |= (p->parentNuclide()->atomicNumber == 94);
        
        unique_isotopes.insert( p->parentNuclide()->symbol );
        src_atomic_numbers.insert( p->parentNuclide()->atomicNumber );
      }else
      {
        assert( p->reaction() );
        unique_isotopes.insert( p->reaction()->name() );
      }
    }//
  }//for( const PeakModel::PeakShrdPtr &p : *m_peakModel->peaks() )
  
  if( has_reaction )
  {
    setup_output.prep_warnings.push_back( WString::tr("ramg-warn-reaction").toUTF8() );
  }//if( user is using reactions )
  
  if( background_sub && !num_peaks_back_sub )
  {
    setup_output.prep_warnings.push_back( WString::tr("ramg-warn-no-bkg-sub-used").toUTF8() );
  }//if( user wanted to background subtract peaks, but no peaks matched up )
  
  
  if( energy_cal_match_warning_energies.size() && (match_tol_sigma > 0.0) )
  {
    const bool multiple = (energy_cal_match_warning_energies.size() > 1);
    
    string nuc_energies, peak_energies;
    for( const auto &pp : energy_cal_match_warning_energies )
    {
      nuc_energies += string(nuc_energies.empty() ? "" : ", ") + SpecUtils::printCompact(pp.second, 4);
      peak_energies += string(peak_energies.empty() ? "" : ", ") + SpecUtils::printCompact(pp.first, 4);
    }
    
    if( multiple )
    {
      nuc_energies = "{" + nuc_energies + "}";
      peak_energies = "{" + peak_energies + "}";
    }
    
    WString msg = WString::tr("ramg-warn-match-outside-tol")
                    .arg( multiple ? "s" : "" )
                    .arg( nuc_energies )
                    .arg( peak_energies );

    setup_output.prep_warnings.push_back( msg.toUTF8() );
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
    {
      assert( n.nuclide || n.reaction );
      
      if( n.nuclide )
        isotopes.emplace_back( n.nuclide->symbol, n.age, n.correct_for_decay_during_meas );
      else
        isotopes.emplace_back( n.reaction->name(), -1.0, false );
    }//for( const auto &n : nuclides_to_match_to )
    
    const float meas_time = fore_spec ? fore_spec->real_time() : -1.0f;
    
    
    // Make sure the match tolerance is ever so slightly above zero, for practical purposes
    const double tol = std::max( match_tol_sigma, 0.0001 );
    
    const RelActCalcManual::PeakCsvInput::NucMatchResults matched_res
      = RelActCalcManual::PeakCsvInput::fill_in_nuclide_info( peak_infos, srcData,
                                                         {}, isotopes, tol, {}, meas_time );
    
    // Add decay correction factors to `matched_res`, and then copy over to RelActCalcManual::RelEffSolution, and then put in the results HTML table
    
    if( matched_res.unused_isotopes.size() )
    {
      string unused_nucs;
      for( size_t i = 0; i < matched_res.unused_isotopes.size(); ++i )
        unused_nucs += string(i ? ", " : "") + matched_res.unused_isotopes[i];
      
      WString msg = WString::tr("ramg-warn-failed-match")
                      .arg( (matched_res.unused_isotopes.size() > 1) ? "s" : "" )
                      .arg( unused_nucs );
      
      throw runtime_error( msg.toUTF8() );
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
    
    setup_output.peaks_before_decay_correction = matched_res.not_decay_corrected_peaks;
  }// End code to fill in nuclide information
  

  // We'll do the actual calculation off of the main thread; in order to make sure the widget
  //  still exists at the end of computations, we'll use WApplication::bind(), in combination
  //  with shared_ptrs's to make sure everything is okay
  const RelActCalc::RelEffEqnForm eqn_form = state->m_relEffEqnFormIndex;
  const size_t eqn_order = static_cast<size_t>( std::max(0,state->m_relEffEqnOrderIndex) );
  
  
  setup_output.peaks = peak_infos;
  setup_output.eqn_form = eqn_form;
  setup_output.eqn_order = eqn_order;
  // We will only use Ceres to fit equation parameters when we have to; using matrix math (i.e.
  //  Eigen) looks to be at least about twice as fast.
  setup_output.use_ceres_to_fit_eqn = (eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel);

  if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    setup_output.use_ceres_to_fit_eqn = true;
    setup_output.phys_model_use_hoerl = setup_input.state->m_physModelUseHoerl;
    setup_output.phys_model_detector = setup_input.detector;

    if( state->m_selfAttenShield )
    {
      setup_output.phys_model_self_atten = state->m_selfAttenShield->fitInput(materialDB);
      assert( setup_output.phys_model_self_atten );
    }

    for( const unique_ptr<RelEffShieldState> &w : state->m_externalShields )
    {
      shared_ptr<const RelActCalc::PhysicalModelShieldInput> rr = w ? w->fitInput(materialDB) : nullptr;
      assert( rr );
      if( rr )
        setup_output.phys_model_external_attens.push_back( rr );
    }//for( auto w : m_extAttenShields->children() )
  }//if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )


  // We will warn about the x-ray absorption edges if the user is not using the Physical Model
  //  and the lowest energy peak is below 122 keV for U or Pu, or below 90 keV for other nuclides.
  //  And even if using the Physical Model, we will warn unless a shielding is same atomic number
  //  as one of the nuclides (this doesnt catch all cases, but is maybe an okay tradef off of not
  //  giving too many false-positives to the user).
  bool warn_xray = has_U_or_Pu ? (lowest_energy_peak < 122) : (lowest_energy_peak < 90);

  if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    const auto shares_atomic_number = [&src_atomic_numbers]( const RelActCalc::PhysicalModelShieldInput &shield ){
      const auto mat = shield.material;

      if( mat )
      {
         for( const auto &nuc_frac : mat->nuclides )
         {
           if( src_atomic_numbers.count(nuc_frac.first->atomicNumber) )
             return true;
         }

         for( const auto &elem_frac : mat->elements )
         {
           if( src_atomic_numbers.count(elem_frac.first->atomicNumber) )
             return true;
         }
      }else
      {
        if( shield.fit_atomic_number )
          return true;

        const int lower_an = static_cast<int>(std::floor(shield.atomic_number));
        const int upper_an = static_cast<int>(std::ceil(shield.atomic_number));
        
        if( src_atomic_numbers.count(lower_an) || src_atomic_numbers.count(upper_an) )
          return true;
      }//if( mat ) / else

      return false;
    };//shares_atomic_number
    
    if( warn_xray && setup_output.phys_model_self_atten )
      warn_xray = (warn_xray && !shares_atomic_number( *setup_output.phys_model_self_atten ));

    for( const auto &s : setup_output.phys_model_external_attens )
    {
      if( !warn_xray )
        break;
      warn_xray = (warn_xray && !shares_atomic_number( *s ));
    }
  }//if( solution.m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  
  if( has_U_or_Pu && warn_xray )
  {
    setup_output.prep_warnings.push_back( WString::tr("ramg-warn-rel-eff-u/pu-xray").toUTF8() );
  }else if( warn_xray )
  {
    setup_output.prep_warnings.push_back( WString::tr("ramg-warn-rel-eff-other-xray").toUTF8() );
  }
}//RelActCalcSetupOutput prepare_calc_input( const RelActCalcRawInput &input )


void RelActManualGui::calculateSolution()
{
  m_currentSolution.reset();
  m_chart->setData( vector<RelActCalcManual::GenericPeakInfo>{}, {}, "", {}, "" );
  m_results->clear();
  
  try
  {
    const RelActCalcRawInput raw_info = get_raw_info_for_calc_input();

    const string sessionId = wApp->sessionId();
    auto solution = make_shared<RelActCalcManual::RelEffSolution>();
    auto updater = wApp->bind( boost::bind( &RelActManualGui::updateGuiWithResults, this, solution ) );
    
    auto errmsg = make_shared<WString>();
    auto err_updater = wApp->bind( boost::bind( &RelActManualGui::updateGuiWithError, this, errmsg ) );
    
    InterSpec *interspec = InterSpec::instance();
    std::shared_ptr<MaterialDB> materialdb = interspec ? interspec->materialDataBaseShared() : nullptr;
    
    // We could almost call `prepare_calc_input(...)` off the GUI thread, and I think the only
    //  impact would be warning messages wouldnt be localized.
    RelActCalcManual::RelEffInput setup_output;
    prepare_calc_input( raw_info, materialdb.get(), setup_output );
    

    WServer::instance()->ioService().boost::asio::io_service::post( std::bind(
      [setup_output, materialdb, sessionId, solution, updater, errmsg, err_updater](){
        try
        {
          *solution = solve_relative_efficiency( setup_output );

          WServer::instance()->post( sessionId, updater );
        }catch( std::exception &e )
        {
          const string exception_msg = e.what();
          
          auto updater = [errmsg, exception_msg, err_updater](){
            *errmsg = WString::tr("ramg-err-performing-calc").arg(exception_msg);
            err_updater();
          };
          
          WServer::instance()->post( sessionId, updater );
        }//try / catch
    } ) );
  }catch( std::exception &e )
  {
    cout << "Error setting up RelActManualGui calc: " << e.what() << endl;
    
    auto msg = std::make_shared<Wt::WString>(WString::tr("ramg-err-setting-up-calc").arg( e.what() ));
    
    updateGuiWithError( msg );
  }//try / catch
}//void calculateSolution()


void RelActManualGui::updateGuiWithError( std::shared_ptr<Wt::WString> error_msg )
{
  //Set error to results TXT and show results TXT tab.
  m_currentSolution = nullptr;
  m_results->clear();
  
  assert( error_msg );
  if( !error_msg )
    return;
  
  WText *error = new WText( *error_msg, m_results );
  error->setInline( false );
  error->addStyleClass( "CalcError" );
  m_resultMenu->select( 0 );
  
  WApplication *app = WApplication::instance();
  assert( app );
  if( app )
    app->triggerUpdate();
}//void updateGuiWithError( std::string error_msg )


void RelActManualGui::updateGuiWithResults( shared_ptr<RelActCalcManual::RelEffSolution> answer )
{
  //Lets make sure we trigger an update, no matter which path we take to leave this function
  DoWorkOnDestruct triggerUpdate( [](){
    WApplication *app = WApplication::instance();
    assert( app );
    if( app )
      app->triggerUpdate();
  });//triggerUpdate
  
  m_currentSolution = answer;
  
  m_results->clear();
  
  if( !m_currentSolution )
  {
    WText *txt = new WText( WString::tr("ramg-no-results-available") );
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
      m_chart->setData( vector<RelActCalcManual::GenericPeakInfo>{}, {}, "", {}, "" );
      break;
    }
      
    case RelActCalcManual::ManualSolutionStatus::Success:
      break;
  }//switch( solution.m_status )
  
  
  if( !solution.m_error_message.empty() )
  {
    WText *errtxt = new WText( WString::tr("ramg-result-error-msg").arg(solution.m_error_message), m_results );
    errtxt->setInline( false );
    errtxt->addStyleClass( "CalcError" );
    m_resultMenu->select( 0 );
  }
  
  
  for( string warning : solution.m_warnings )
  {
    WText *warntxt = new WText( WString::tr("ramg-result-warn-msg").arg(warning), m_results );
    warntxt->setInline( false );
    warntxt->addStyleClass( "CalcWarning" );
  }
  
  
  if( solution.m_status != RelActCalcManual::ManualSolutionStatus::Success )
    return;
  
  // We'll first update the chart
  string relEffEqn = solution.rel_eff_eqn_js_function();
  string relEffEqnUncert = solution.rel_eff_eqn_js_uncert_fcn();
  if( solution.m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    // Update shield widgets
    const auto &input = solution.m_input;
    assert( input.phys_model_external_attens.size() == solution.m_phys_model_external_atten_shields.size() );

    const auto update_shield = [&]( RelEffShieldWidget *w,
                      const unique_ptr<RelActCalcManual::RelEffSolution::PhysModelShieldFit> &fit,
                      const shared_ptr<const RelActCalc::PhysicalModelShieldInput> &input ){
      assert( w );
      if( !w )
        return;
      
      assert( input );
      
      if( !fit )
      {
        w->resetState();
      }else if( fit->m_material )
      {
        w->setMaterialSelected( true );
        if( w->materialNameTxt().toUTF8() != fit->m_material->name )
          w->setMaterial( fit->m_material->name );
        
        const double thickness = fit->m_areal_density / fit->m_material->density;
        
        if( !input || (input->fit_areal_density) )
        {
          w->setThickness( thickness );
        }else
        {
          // We'll just double check that thickness was actually not changed.
          assert( fabs(fit->m_areal_density - input->areal_density)
                 <= 0.00001*std::max(fit->m_areal_density,input->areal_density) );
          
          const double dx = w->thickness() - thickness;
          assert( fabs(dx) <= 0.00001*std::max(w->thickness(),thickness) );
          if( fabs(dx) > 0.00001*std::max(w->thickness(),thickness) )
            w->setThickness( thickness );
        }//if( input->fit_areal_density ) / else
        

        if( input )
          w->setFitThickness( input->fit_areal_density );
      }else
      {
        w->setMaterialSelected( false );
        
        if( !input || input->fit_atomic_number )
        {
          w->setAtomicNumber( fit->m_atomic_number );
        }else if( input )
        {
          const double prev_an = w->atomicNumber();
          const double new_an = fit->m_atomic_number;
          const double delta_an = fabs(new_an - prev_an);
          assert( delta_an <= 0.0001*std::max(prev_an, new_an) );
          if( delta_an > 0.0001*std::max(prev_an, new_an) )
            w->setAtomicNumber( new_an );
        }//if( fit_atomic_number ) / else
        
        const double new_ad = fit->m_areal_density / PhysicalUnits::g_per_cm2;
        if( !input || input->fit_areal_density )
        {
          w->setArealDensity( new_ad );
        }else if( input )
        {
          const double prev_ad = w->arealDensity();
          const double delta_ad = fabs(new_ad - prev_ad);
          assert( delta_ad <= 0.0001*std::max(prev_ad, new_ad) );
          if( delta_ad > 0.0001*std::max(prev_ad, new_ad) )
            w->setAtomicNumber( new_ad );
        }//if( fit_areal_density ) / else
        
        if( input )
        {
          w->setFitAtomicNumber( input->fit_atomic_number );
          w->setFitArealDensity( input->fit_areal_density );
        }
      }//if( fit->m_material ) / else
    };//update_shield lambda

    if( !m_selfAttenShield )
      initPhysicalModelAttenShieldWidgets();
    
    if( m_selfAttenShield )
      update_shield( m_selfAttenShield, solution.m_phys_model_self_atten_shield, input.phys_model_self_atten );
    
    // We shouldnt need to update the number of external attenuation shield widgets,
    //  but we'll go ahead and implement it anyway, incase something odd changed while 
    //  computation was happening.
    assert( m_extAttenShields );
    vector<RelEffShieldWidget *> ext_shields;
    vector<WWidget *> ext_kids = m_extAttenShields ? m_extAttenShields->children() : vector<WWidget *>{};
    for( auto w : ext_kids )
    {
      const auto shield = dynamic_cast<RelEffShieldWidget *>(w);
      if( shield && ((ext_shields.size() < solution.m_phys_model_external_atten_shields.size()) || shield->nonEmpty()) )
      {
        ext_shields.push_back( shield );
      }else if( shield && !shield->nonEmpty() )
      {
        assert( !shield->nonEmpty() );
        delete shield;
      }
    }//for( auto w : ext_kids )

    assert( ext_shields.size() == solution.m_phys_model_external_atten_shields.size() ); //Again, should be the right size
    // We shouldnt need to add any more shields, but we'll go ahead and implement doing this
    while( ext_shields.size() < solution.m_phys_model_external_atten_shields.size() )
    {
      auto shield = new RelEffShieldWidget( RelEffShieldWidget::ShieldType::ExternalAtten, m_extAttenShields );
      shield->changed().connect( this, &RelActManualGui::handlePhysicalModelShieldChanged );
      ext_shields.push_back( shield );
    }

    assert( ext_shields.size() == input.phys_model_external_attens.size() );
    assert( ext_shields.size() == solution.m_phys_model_external_atten_shields.size() );
    
    for( size_t i = 0; i < std::max(ext_shields.size(), solution.m_phys_model_external_atten_shields.size()); ++i )
    {
      if( i >= solution.m_phys_model_external_atten_shields.size() )
      {
        ext_shields[i]->resetState();
        continue;
      }
      const auto &fit_val = solution.m_phys_model_external_atten_shields[i];
      const shared_ptr<const RelActCalc::PhysicalModelShieldInput> in_shield = input.phys_model_external_attens[i];

      assert( fit_val );
      if( i < ext_shields.size() )
        update_shield( ext_shields[i], fit_val, in_shield );
    }//for( loop over external shields )

    ext_kids = m_extAttenShields ? m_extAttenShields->children() : vector<WWidget *>{};
    if( ext_kids.empty() )
    {
      auto shield = new RelEffShieldWidget( RelEffShieldWidget::ShieldType::ExternalAtten, m_extAttenShields );
      shield->changed().connect( this, &RelActManualGui::handlePhysicalModelShieldChanged );
    }
  }//if( !FramPhysicalModel ) / else
  
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
  
  WString pval_str;
  try
  {
    boost::math::chi_squared chi_squared_dist( solution.m_dof );
    const double prob = boost::math::cdf( chi_squared_dist, solution.m_chi2 );
    if( prob > 0.99 )
      pval_str = WString::tr("ramg-1-pval").arg( SpecUtils::printCompact(1.0 - prob, 3) );
    else
      pval_str = WString::tr("ramg-pval").arg( SpecUtils::printCompact(prob, 3) );
  }catch( std::exception &e )
  {
    pval_str = "";
  }
  
  WString chi2_title = WString::tr("ramg-chart-info-title");
  chi2_title.arg( SpecUtils::printCompact(solution.m_chi2, 3) )
            .arg( static_cast<int>(solution.m_dof) )
            .arg( pval_str );
  
  // If we have U or Pu, we'll give the enrichment, or if we have two nuclides we'll
  //  give their ratio.  If we have U and Pu, we wont give enrichment.
  set<string> isotopes;
  for( const auto &relact : solution.m_rel_activities )
    isotopes.insert( relact.m_isotope );
  if( (isotopes.count("U235") && isotopes.count("U238") && !isotopes.count("Pu239"))
     || (isotopes.count("Pu239") && isotopes.count("Pu240") && !isotopes.count("U235")) )
  {
    string enrich;
    const string iso = isotopes.count("U235") ? "U235" : "Pu239";
    
    try
    {
      const double nominal = solution.mass_fraction(iso);
      const double plus = solution.mass_fraction( iso, 1.0 );
      const double minus = solution.mass_fraction( iso, -1.0 );
      const double error = 0.5*( fabs(plus - nominal) + fabs(nominal - minus) );
      enrich = ", " + PhysicalUnits::printValueWithUncertainty(100.0*nominal, 100.0*error, 4)
               + "% "+ iso;
    }catch( std::exception & )
    {
      // We dont have the covariance matrix required to get mass fraction errors
      //  We'll try leaving out the errors
      try 
      {
        const double nominal = solution.mass_fraction(iso);
        enrich = ", " + SpecUtils::printCompact(100.0*nominal, 4) + "% "+ iso;
      }catch( std::exception & )
      {
        // I really dont think we should ever get here
        assert( 0 );
      }
    }//try / catch
    
    chi2_title.arg( enrich );
  }else if( isotopes.size() == 2 )
  {
    const vector<RelActCalcManual::IsotopeRelativeActivity> &rel_acts = solution.m_rel_activities;
    const int num_index = (rel_acts[0].m_rel_activity > rel_acts[1].m_rel_activity) ? 1 : 0;
    const int denom_index = (num_index ? 0 : 1);
    const string num_nuc = rel_acts[num_index].m_isotope;
    const string den_nuc = rel_acts[denom_index].m_isotope;
    const double ratio = solution.activity_ratio(num_nuc, den_nuc);
    const double uncert = solution.activity_ratio_uncert(num_nuc, den_nuc);
    
    string ratio_txt = ", act(" + num_nuc + "/" + den_nuc + ")="
                   + PhysicalUnits::printValueWithUncertainty(ratio, uncert, 4);
    chi2_title.arg( ratio_txt );
  }else
  {
    chi2_title.arg( "" );
  }
  
  m_chart->setData( solution.m_input.peaks, relActsColors, relEffEqn, chi2_title, relEffEqnUncert );
  
  
  // Now update the text
  // TODO: refactor getting these tables from the solution...
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  char buffer[2048] = {'\0'};
  stringstream results_html;
  
  solution.get_mass_fraction_table( results_html );
  results_html << "<br />";
  solution.get_mass_ratio_table( results_html );
  results_html << "<br />";
  
  results_html << "<div class=\"releffeqn\">Rel. Eff.: y = ";
  results_html << solution.rel_eff_eqn_txt(true);

  results_html << "</div>\n";
  
  results_html << "<br /> <div>&chi;<sup>2</sup>=" << SpecUtils::printCompact( solution.m_chi2, 4)
  << " " << WString::tr("ramg-and-there-were").toUTF8() << " " << solution.m_dof
  << " DOF (&chi;<sup>2</sup>/<sub>" << WString::tr("ramg-dof").toUTF8() << "</sub>="
  << SpecUtils::printCompact(solution.m_chi2/solution.m_dof, 4) << ")</div>\n";
  
  
  results_html << "<div class=\"ToolAlphaWarning\">";
  
  switch( AddUncert(m_addUncertainty->currentIndex()) )
  {
    case AddUncert::Unweighted:
      results_html << WString::tr("ramg-fit-unweighted-txt").toUTF8();
      break;
      
    case AddUncert::StatOnly:
      results_html << WString::tr("ramg-fit-stat-only");
      break;
      
    case AddUncert::OnePercent:
    case AddUncert::FivePercent:
    case AddUncert::TenPercent:
    case AddUncert::TwentyFivePercent:
    case AddUncert::FiftyPercent:
    case AddUncert::SeventyFivePercent:
    case AddUncert::OneHundredPercent:
      results_html << WString::tr("ramg-fit-uncert-increased").toUTF8();
      break;
      
    case AddUncert::NumAddUncert:
      assert( 0 );
      break;
  }//switch( add_uncert_type )
  
  results_html << "</div>\n";
  
  
  new WText( results_html.str(), m_results );
}//void updateGuiWithResults( const RelActCalcManual::RelEffSolution &solution );


void RelActManualGui::relEffEqnFormChanged( const bool user_action )
{
  const auto eqn_form = RelActCalc::RelEffEqnForm( m_relEffEqnForm->currentIndex() );
  switch( eqn_form )
  {
    case RelActCalc::RelEffEqnForm::LnX:
    case RelActCalc::RelEffEqnForm::LnY:
    case RelActCalc::RelEffEqnForm::LnXLnY:
    case RelActCalc::RelEffEqnForm::FramEmpirical:
      m_physicalModelShields->hide();
      m_nucColumnTitle->setText( WString::tr("ramg-nucs-label") );
      m_relEffEqnOrderHolder->show();
      m_physModelUseHoerlHolder->hide();
      break;
    
    case RelActCalc::RelEffEqnForm::FramPhysicalModel:
      m_physicalModelShields->show();
      m_nucColumnTitle->setText( WString::tr("ramg-nucs-shield-label") );
      m_relEffEqnOrderHolder->hide();
      m_physModelUseHoerlHolder->show();
      initPhysicalModelAttenShieldWidgets();
      break;
  }//switch( eqn_form )

  
  m_renderFlags |= RenderActions::UpdateCalc;
  if( user_action )
    m_renderFlags |= RenderActions::AddUndoRedoStep;
  scheduleRender();
}//void relEffEqnFormChanged()


void RelActManualGui::relEffEqnOrderChanged()
{
  m_renderFlags |= RenderActions::UpdateCalc;
  m_renderFlags |= RenderActions::AddUndoRedoStep;
  scheduleRender();
}


void RelActManualGui::physModelUseHoerlChanged()
{
  m_renderFlags |= RenderActions::UpdateCalc;
  m_renderFlags |= RenderActions::AddUndoRedoStep;
  scheduleRender();
}

void RelActManualGui::nucDataSrcChanged()
{
  m_renderFlags |= RenderActions::UpdateCalc;
  m_renderFlags |= RenderActions::AddUndoRedoStep;
  scheduleRender();
}

void RelActManualGui::handlePhysicalModelShieldChanged()
{
  m_renderFlags |= RenderActions::UpdateCalc;
  m_renderFlags |= RenderActions::AddUndoRedoStep;
  scheduleRender();
}//void handlePhysicalModelShieldChanged()


void RelActManualGui::matchToleranceChanged()
{
  m_renderFlags |= RenderActions::UpdateCalc;
  m_renderFlags |= RenderActions::AddUndoRedoStep;
  scheduleRender();
}//void matchToleranceChanged();


void RelActManualGui::addUncertChanged()
{
  m_renderFlags |= RenderActions::UpdateCalc;
  m_renderFlags |= RenderActions::AddUndoRedoStep;
  scheduleRender();
}


void RelActManualGui::backgroundSubtractChanged()
{
  m_renderFlags |= RenderActions::UpdateCalc;
  m_renderFlags |= RenderActions::AddUndoRedoStep;
  scheduleRender();
}

void RelActManualGui::initPhysicalModelAttenShieldWidgets()
{
  assert( m_physicalModelShields );

  if( m_selfAttenShield )
  {
    m_selfAttenShield->resetState();
  }else
  {      
    m_selfAttenShield = new RelEffShieldWidget( RelEffShieldWidget::ShieldType::SelfAtten );
    m_selfAttenShield->changed().connect( this, &RelActManualGui::handlePhysicalModelShieldChanged );
    m_physicalModelShields->insertWidget( 0, m_selfAttenShield );
  }
  
  assert( m_extAttenShields );
  if( m_extAttenShields )
  {
    m_extAttenShields->clear();
    auto shield = new RelEffShieldWidget( RelEffShieldWidget::ShieldType::ExternalAtten, m_extAttenShields );
    shield->changed().connect( this, &RelActManualGui::handlePhysicalModelShieldChanged );
  }//if( m_extAttenShields )
}//void initPhysicalModelAttenShieldWidgets()

//void RelActManualGui::resultTabChanged()
//{
//  m_renderFlags |= RenderActions::AddUndoRedoStep;
//  scheduleRender();
//}


void RelActManualGui::handlePeaksChanged()
{
  m_renderFlags |= RenderActions::UpdateCalc;
  m_renderFlags |= RenderActions::UpdateNuclides;
  m_renderFlags |= RenderActions::AddUndoRedoStep;
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
  map<const ReactionGamma::Reaction *,ManRelEffNucDisp *> current_rctns;
  
  for( auto w : m_nuclidesDisp->children() )
  {
    ManRelEffNucDisp *rr = dynamic_cast<ManRelEffNucDisp *>(w);
    if( !rr )
      continue; // Shielding inputs for Physical Model
    assert( rr && (rr->m_nuc || rr->m_reaction) && !current_nucs.count(rr->m_nuc) && !current_rctns.count(rr->m_reaction) );
    if( rr && rr->m_nuc )
      current_nucs[rr->m_nuc] = rr;
    else if( rr && rr->m_reaction )
      current_rctns[rr->m_reaction] = rr;
  }//for( auto w : m_nuclidesDisp->children() )
  
  shared_ptr<const SpecUtils::Measurement> foreground = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  const float meas_time = foreground ? foreground->real_time() : -1.0f;
  
  set<const SandiaDecay::Nuclide *> nucs_in_peaks;
  set<const ReactionGamma::Reaction *> reactions_in_peaks;
  PeakModel *peakModel = m_interspec->peakModel();
  const auto peaks = peakModel ? peakModel->peaks() : nullptr;
  if( peaks )
  {
    for( const auto &p : *peaks )
    {
      if( p->useForManualRelEff() )
      {
        if( p->parentNuclide() )
          nucs_in_peaks.insert( p->parentNuclide() );
        else if( p->reaction() )
          reactions_in_peaks.insert( p->reaction() );
      }//if( p->useForManualRelEff() )
    }//for( const auto &p : *peaks )
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
      assert( prev || kids[index]->hasStyleClass("PhysicalModelShields") );
      if( !prev || !prev->m_nuc )
      {
        if( kids[index] && kids[index]->hasStyleClass("PhysicalModelShields") )
          insert_index = static_cast<int>( index + 1 );  //Make sure shields are always before nuclides
        continue;
      }
      
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
    
    ManRelEffNucDisp *rr = new ManRelEffNucDisp( nuc, nullptr, age, meas_time );
    rr->updated().connect( this, &RelActManualGui::handlePeaksChanged );
    
    const auto decay_corr_pos = m_nucDecayCorrect.find(nuc->symbol);
    if( decay_corr_pos != end(m_nucDecayCorrect) )
    {
      if( rr->m_decay_during_meas )
        rr->setDecayDuringMeasurement( decay_corr_pos->second );
    }//if( we cached if we should decay correct this nuclide )
    
    m_nuclidesDisp->insertWidget( insert_index, rr );
  }//for( loop over to add displays for new nuclides )
  
  
  // Go through and add any new nuclide from peaks
  for( const ReactionGamma::Reaction *reaction : reactions_in_peaks )
  {
    auto pos = current_rctns.find( reaction );
    if( pos != end(current_rctns) )
      continue;
    
    const SandiaDecay::Element *current_el = reaction->targetElement;
    if( !current_el && reaction->productNuclide )
    {
      const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
      current_el = db->element( reaction->productNuclide->atomicNumber );
    }
    
    // We'll insert the display sorted by {atomicNumber, alphabetical}
    int insert_index = 0;
    const auto kids = m_nuclidesDisp->children();
    for( size_t index = 0; index < kids.size(); ++index )
    {
      ManRelEffNucDisp *prev = dynamic_cast<ManRelEffNucDisp *>( kids[index] );
      assert( prev );
      if( !prev || !prev->m_reaction )
      {
        if( kids[index] && kids[index]->hasStyleClass("PhysicalModelShields") )
          insert_index = static_cast<int>( index + 1 );  //Make sure shields are always before nuclides
        continue;
      }
      
      const SandiaDecay::Element *prev_el = prev->m_reaction->targetElement;
      if( !prev_el && prev->m_reaction->productNuclide )
      {
        const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
        prev_el = db->element( prev->m_reaction->productNuclide->atomicNumber );
      }
        
      if( prev_el && current_el )
      {
        if( (prev_el->atomicNumber > current_el->atomicNumber)
           || ((prev_el->atomicNumber == current_el->atomicNumber)
               && (prev->m_reaction->name() > reaction->name())) )
        {
          insert_index = static_cast<int>( index );
          break;
        }
      }else if( prev->m_reaction->name() > reaction->name() )
      {
        insert_index = static_cast<int>( index );
        break;
      }
      
      insert_index = static_cast<int>( index + 1 );
    }//for( loop over existing displays to find position )
    
    ManRelEffNucDisp *rr = new ManRelEffNucDisp( nullptr, reaction, 0.0, meas_time );
    //ManRelEffNucDisp::updated() is only emitted for age changes; not relevant for reactions,
    //  but we'll connect up JIC for the future.
    rr->updated().connect( this, &RelActManualGui::handlePeaksChanged );
    
    m_nuclidesDisp->insertWidget( insert_index, rr );
  }//for( loop over to add displays for new nuclides )
  
  
  // Loop over and remove displays for any nuclides we no longer need
  for( const auto nuc_widget : current_nucs )
  {
    if( nucs_in_peaks.count(nuc_widget.first) )
      continue;
    
    m_nucAge[nuc_widget.first->symbol] = nuc_widget.second->m_current_age;
    
    if( nuc_widget.second->m_decay_during_meas )
    {
      m_nucDecayCorrect[nuc_widget.first->symbol] = nuc_widget.second->decayDuringMeasurement();
    }else
    {
      auto decay_corr_pos = m_nucDecayCorrect.find(nuc_widget.first->symbol);
      if( decay_corr_pos != end(m_nucDecayCorrect) )
        m_nucDecayCorrect.erase( decay_corr_pos );
    }
    
    delete nuc_widget.second;
  }//for( loop over to remove any nuclides )
  
  // Loop over and remove displays for any reaction we no longer need
  for( const auto nuc_widget : current_rctns )
  {
    if( !reactions_in_peaks.count(nuc_widget.first) )
      delete nuc_widget.second;
  }//for( loop over to remove any nuclides )
  
  // We may have deleted some of current_nucs or current_rctns, so lets clear it,
  //  just to make sure we dont access
  current_nucs.clear();
  current_rctns.clear();
  
  bool has_uranium = false;
  for( auto w : m_nuclidesDisp->children() )
  {
    ManRelEffNucDisp *rr = dynamic_cast<ManRelEffNucDisp *>(w);
    assert( rr || w->hasStyleClass("PhysicalModelShields") );
    if( rr && rr->m_nuc )
    {
      const bool isU = (rr->m_nuc->atomicNumber == 92);
      has_uranium |= isU;
      rr->setAgeHidden( isU ? !showAge : false );
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
  
  
  // Lets give a friendly reminder of "Add. Uncert" wrt if Uranium is present
  auto spec = m_interspec ? m_interspec->measurment(SpecUtils::SpectrumType::Foreground) : nullptr;
  if( spec )
  {
    const bool stat_only = (m_addUncertainty->currentIndex() == static_cast<int>(AddUncert::StatOnly));
    const set<int> &displayed = m_interspec->displayedSamples(SpecUtils::SpectrumType::Foreground);
    const auto peaks = spec->peaks( displayed );
    
    if( peaks )
    {
      bool have_u = false;
      for( size_t i = 0; !have_u && (i < peaks->size()); ++i )
        have_u = ((*peaks)[i]->parentNuclide() && ((*peaks)[i]->parentNuclide()->atomicNumber == 92) );
      
      WString msg;
      if( have_u )
      {
        if( stat_only )
          msg = WString::tr("ramg-consider-add-uncert-u");
        else
          msg = WString::tr("ramg-you-using-add-uncert-u");
      }else if( !have_u && !stat_only )
      {
        msg = WString::tr("ramg-you-using-add-uncert-non-u");
      }
      
      if( !msg.empty() )
        m_interspec->logMessage( msg, 2 );
    }//if( !peaks || peaks->empty() )
  }//if( spec )
  
  
  m_addUncertainty->setCurrentIndex( static_cast<int>(AddUncert::StatOnly) );
  
  
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
    case RelActCalc::RelEffEqnForm::FramPhysicalModel:
      validSelect = true;
      break;
  }//switch( eqn_form )
  
  if( !validSelect )
    throw runtime_error( "Invalid RelEffEqnForm" );
  
  return eqn_form;
}//RelEffEqnForm relEffEqnForm() const


size_t RelActManualGui::relEffEqnOrder() const
{
  RelActCalc::RelEffEqnForm eqgnForm = relEffEqnForm();
  if( eqgnForm == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    return 0;

  const int orderIndex = m_relEffEqnOrder->currentIndex();
  if( (orderIndex < 0) || (orderIndex > 7) )
    throw runtime_error( "Invalid RelEffEqnOrder" );
  
  return static_cast<size_t>(orderIndex);
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


::rapidxml::xml_node<char> *RelActManualGui::serialize( ::rapidxml::xml_node<char> *parent_node )
{
  auto state = getGuiState();
  if( !state )
    return nullptr;
  
  return state->serialize( parent_node );
}//serialize(...)


void RelActManualGui::deSerialize( const ::rapidxml::xml_node<char> *base_node )
{
  // Clear out any existing nuclides
  GuiState state;
  state.deSerialize( base_node );
  
  setGuiState( state );
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


::rapidxml::xml_node<char> *RelActManualGui::GuiState::serialize( ::rapidxml::xml_node<char> *parent_node ) const
{
  ::rapidxml::xml_document<char> *doc = parent_node->document();
  assert( doc );
  
  const GuiState &state = *this;
  
  const char *name = "RelActManualGui";
  ::rapidxml::xml_node<char> * const base_node = doc->allocate_node( ::rapidxml::node_element, name );
  parent_node->append_node( base_node );
  
  //If you change the available options or formatting or whatever, increment the
  //  version field of the XML!
  const string versionstr = std::to_string(RelActManualGui::sm_xmlSerializationMajorVersion)
                            + "." + std::to_string(RelActManualGui::sm_xmlSerializationMinorVersion);
  const char *value = doc->allocate_string( versionstr.c_str() );
  ::rapidxml::xml_attribute<> *attr = doc->allocate_attribute( "version", value );
  base_node->append_attribute( attr );
  
  value = RelActCalc::to_str( state.m_relEffEqnFormIndex );
  ::rapidxml::xml_node<char> *node = doc->allocate_node( ::rapidxml::node_element, "RelEffEqnForm", value );
  base_node->append_node( node );
  
  value = doc->allocate_string( std::to_string(state.m_relEffEqnOrderIndex).c_str() );
  node = doc->allocate_node( ::rapidxml::node_element, "RelEffEqnOrder", value );
  base_node->append_node( node );
  
  value = state.m_physModelUseHoerl ? "true" : "false";
  node = doc->allocate_node( ::rapidxml::node_element, "PhysModelUseHoerl", value );
  base_node->append_node( node );

  value = RelActCalcManual::PeakCsvInput::to_str( state.m_nucDataSrcIndex );
  node = doc->allocate_node( ::rapidxml::node_element, "NucDataSrc", value );
  base_node->append_node( node );
  
  const string match_tolerance_str = SpecUtils::printCompact(state.m_matchToleranceValue, 7);
  value = doc->allocate_string( match_tolerance_str.c_str() );
  node = doc->allocate_node( ::rapidxml::node_element, "MatchTolerance", value );
  base_node->append_node( node );
  
  value = RelActManualGui::to_str( state.m_addUncertIndex );
  node = doc->allocate_node( ::rapidxml::node_element, "AddUncertainty", value );
  base_node->append_node( node );
  
  value = state.m_resultTab ? "1" : "0";
  node = doc->allocate_node( ::rapidxml::node_element, "ResultTabShowing", value );
  base_node->append_node( node );
  
  value = state.m_backgroundSubtract ? "1" : "0";
  node = doc->allocate_node( ::rapidxml::node_element, "BackgroundSubtract", value );
  base_node->append_node( node );
  
  ::rapidxml::xml_node<char> *nuc_ages_node = doc->allocate_node( ::rapidxml::node_element, "NuclideAges" );
  base_node->append_node( nuc_ages_node );
  for( const auto &n : state.nucAge )
  {
    ::rapidxml::xml_node<char> *nuc_node = doc->allocate_node( ::rapidxml::node_element, "Nuclide" );
    nuc_ages_node->append_node( nuc_node );
    
    value = doc->allocate_string( n.first.c_str() );
    ::rapidxml::xml_node<char> *name_node = doc->allocate_node( ::rapidxml::node_element, "Name", value );
    nuc_node->append_node(name_node);
    
    value = doc->allocate_string( SpecUtils::printCompact(n.second,8).c_str() );
    ::rapidxml::xml_node<char> *age_node = doc->allocate_node( ::rapidxml::node_element, "Age", value );
    nuc_node->append_node(age_node);
    
    const auto decay_corr_pos = state.nucDecayCorrect.find(n.first);
    if( decay_corr_pos != end(state.nucDecayCorrect) )
    {
      value = decay_corr_pos->second ? "true" : "false";
      ::rapidxml::xml_node<char> *corr_node = doc->allocate_node( ::rapidxml::node_element, "DecayDuringMeasurement", value );
      nuc_node->append_node(corr_node);
    }
  }//for( const auto &n : nuc_age_cache )
  
  if( state.m_selfAttenShield || !state.m_externalShields.empty() )
  {
    ::rapidxml::xml_node<char> *shields_node = doc->allocate_node( ::rapidxml::node_element, "PhysicalModelShields" );
    base_node->append_node( shields_node );
    
    if( state.m_selfAttenShield )
    {
      ::rapidxml::xml_node<char> *self_node = doc->allocate_node( ::rapidxml::node_element, "SelfAttenShield" );
      shields_node->append_node( self_node );
      state.m_selfAttenShield->toXml( self_node );
    }//if( state.self_atten_shield )
    
    if( !state.m_externalShields.empty() )
    {
      ::rapidxml::xml_node<char> *ext_node = doc->allocate_node( ::rapidxml::node_element, "ExternalAttenShields" );
      shields_node->append_node( ext_node );
      for( const unique_ptr<RelEffShieldState> &s : state.m_externalShields )
        s->toXml( ext_node );
    }//if( !ext_atten_shields.empty() )
  }//if( state.self_atten_shield || !state.ext_atten_shields.empty() )

  return base_node;
}//::rapidxml::xml_node<char> *GuiState::serialize( ::rapidxml::xml_node<char> *parent_node )


void RelActManualGui::GuiState::deSerialize( const ::rapidxml::xml_node<char> *base_node )
{
  if( SpecUtils::xml_name_str(base_node) != "RelActManualGui" )
    throw runtime_error( "RelActManualGui::deSerialize: invalid base node passed in: '"
                        + SpecUtils::xml_name_str(base_node) + "'" );
    
  int version;
  const ::rapidxml::xml_attribute<char> *attr = XML_FIRST_ATTRIB(base_node, "version");
  if( !attr || !attr->value() || !(stringstream(attr->value()) >> version) )
    throw runtime_error( "Deserializing requires a version" );
  
  if( version != RelActManualGui::sm_xmlSerializationMajorVersion )
    throw runtime_error( "Invalid version of RelActManualGui XML" );
  
  const ::rapidxml::xml_node<char> *RelEffEqnForm_node = XML_FIRST_NODE(base_node, "RelEffEqnForm");
  const ::rapidxml::xml_node<char> *RelEffEqnOrder_node = XML_FIRST_NODE(base_node, "RelEffEqnOrder");
  const ::rapidxml::xml_node<char> *PhysModelUseHoerl_node = XML_FIRST_NODE(base_node, "PhysModelUseHoerl");
  const ::rapidxml::xml_node<char> *NucDataSrc_node = XML_FIRST_NODE(base_node, "NucDataSrc");
  const ::rapidxml::xml_node<char> *MatchTolerance_node = XML_FIRST_NODE(base_node, "MatchTolerance");
  const ::rapidxml::xml_node<char> *AddUncertainty_node = XML_FIRST_NODE(base_node, "AddUncertainty");
  const ::rapidxml::xml_node<char> *ResultTabShowing_node = XML_FIRST_NODE(base_node, "ResultTabShowing");
  const ::rapidxml::xml_node<char> *BackgroundSubtract_node = XML_FIRST_NODE(base_node, "BackgroundSubtract");
  const ::rapidxml::xml_node<char> *NuclideAges_node = XML_FIRST_NODE(base_node, "NuclideAges");
  const ::rapidxml::xml_node<char> *PhysicalModelShields_node = XML_FIRST_NODE(base_node, "PhysicalModelShields");
  
  if( !RelEffEqnForm_node || !NucDataSrc_node || !MatchTolerance_node
     || !AddUncertainty_node || !ResultTabShowing_node || !NuclideAges_node )
    throw runtime_error( "RelActManualGui::deSerialize: missing required node" );
  
  GuiState &state = *this;
  
  bool got_eqn_form = false;
  state.m_relEffEqnFormIndex = RelActCalc::RelEffEqnForm::LnX;
  const string rel_eff_eqn_form_str = SpecUtils::xml_value_str( RelEffEqnForm_node );
  for( RelActCalc::RelEffEqnForm form = RelActCalc::RelEffEqnForm(0);
      form <= RelActCalc::RelEffEqnForm::FramPhysicalModel;
      form = RelActCalc::RelEffEqnForm(static_cast<int>(form)+1) )
  {
    const char * const test_form = RelActCalc::to_str( form );
    if( test_form == rel_eff_eqn_form_str )
    {
      got_eqn_form = true;
      state.m_relEffEqnFormIndex = form;
      break;
    }
  }//for( loop over RelEffEqnForm )
  
  if( !got_eqn_form )
    throw runtime_error( "RelActManualGui::deSerialize: '" + rel_eff_eqn_form_str + "' is invalid RelEffEqnForm" );
  
  state.m_relEffEqnOrderIndex = 0; //We'll let this node be optional for physical mode.
  if( RelEffEqnOrder_node || (state.m_relEffEqnFormIndex != RelActCalc::RelEffEqnForm::FramPhysicalModel) )
  {
    const string rel_eff_order_str = SpecUtils::xml_value_str(RelEffEqnOrder_node);
    if( !(stringstream(rel_eff_order_str) >> state.m_relEffEqnOrderIndex)
       || (state.m_relEffEqnOrderIndex < 0) )
    {
      throw runtime_error( "RelActManualGui::deSerialize: '" + rel_eff_order_str + "' is invalid RelEffEqnOrder" );
    }
  }//if( we expect, or are seeing equation order node ).

  state.m_physModelUseHoerl = true;
  if( PhysModelUseHoerl_node )
  {
    const string use_hoerl_str = SpecUtils::xml_value_str( PhysModelUseHoerl_node );
    state.m_physModelUseHoerl = ((use_hoerl_str == "true") || (use_hoerl_str == "1"));
  }

  
  bool got_data_src = false;
  state.m_nucDataSrcIndex = RelActCalcManual::PeakCsvInput::NucDataSrc::Icrp107_U;
  const string data_src_str = SpecUtils::xml_value_str(NucDataSrc_node);
  for( auto src = RelActCalcManual::PeakCsvInput::NucDataSrc(0);
      src < RelActCalcManual::PeakCsvInput::NucDataSrc::Undefined;
      src = RelActCalcManual::PeakCsvInput::NucDataSrc(static_cast<int>(src) + 1) )
  {
    const char *val = RelActCalcManual::PeakCsvInput::to_str( src );
    if( val == data_src_str )
    {
      got_data_src = true;
      state.m_nucDataSrcIndex = src;
      break;
    }
  }
  
  if( !got_data_src )
    throw runtime_error( "RelActManualGui::deSerialize: '" + data_src_str + "' is invalid NucDataSrc" );
  
  const string match_tol_str = SpecUtils::xml_value_str(MatchTolerance_node);
  if( !(stringstream( match_tol_str ) >> state.m_matchToleranceValue)
     || (state.m_matchToleranceValue < 0.0f)
     || (state.m_matchToleranceValue > 5.0f)
     || IsNan(state.m_matchToleranceValue))
  {
    throw runtime_error( "RelActManualGui::deSerialize: '" + match_tol_str + "' is invalid match tolerance" );
  }
  
  bool got_add_uncert = false;
  state.m_addUncertIndex = RelActManualGui::AddUncert::FiftyPercent;
  const string add_uncert_str = SpecUtils::xml_value_str(AddUncertainty_node);
  for( RelActManualGui::AddUncert uncert = RelActManualGui::AddUncert(0);
      uncert <= RelActManualGui::AddUncert::NumAddUncert;
      uncert = RelActManualGui::AddUncert(static_cast<int>(uncert) + 1) )
  {
    const char *test_str = RelActManualGui::to_str( uncert );
    if( test_str == add_uncert_str )
    {
      got_add_uncert = true;
      state.m_addUncertIndex = uncert;
      break;
    }
  }//for( loop over AddUncert )
  
  if( !got_add_uncert )
    cerr << "RelActManualGui::deSerialize: failed to decode '"
         << add_uncert_str << "' into a AddUncert enum" << endl;
  
  state.m_resultTab = 0;
  const string tab_show_str = SpecUtils::xml_value_str(ResultTabShowing_node);
  if( !(stringstream(tab_show_str) >> state.m_resultTab) )
    cerr << "RelActManualGui::deSerialize: invalid ResultTabShowing node '"
         << tab_show_str << "'" << endl;
  
  if( (state.m_resultTab != 0) && (state.m_resultTab != 1) )
    state.m_resultTab = 0;
  
  state.m_backgroundSubtract = false;
  if( BackgroundSubtract_node )
  {
    const string back_sub_str = SpecUtils::xml_value_str(BackgroundSubtract_node);
    state.m_backgroundSubtract = ((back_sub_str == "1") || SpecUtils::iequals_ascii(back_sub_str, "true"));
  }
  
  
  XML_FOREACH_CHILD( nuc_node, NuclideAges_node, "Nuclide" )
  {
    ::rapidxml::xml_node<char> *name_node = XML_FIRST_NODE(nuc_node, "Name");
    ::rapidxml::xml_node<char> *age_node = XML_FIRST_NODE(nuc_node, "Age");
    ::rapidxml::xml_node<char> *decay_corr_node = XML_FIRST_NODE(nuc_node, "DecayDuringMeasurement");
    
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
    
    state.nucAge[nuc] = age;
    
    if( decay_corr_node )
    {
      const string decay_corr = SpecUtils::xml_value_str(decay_corr_node);
      state.nucDecayCorrect[nuc] = ((decay_corr == "true") || (decay_corr == "1"));
    }
  }//for( loop over <Nuclide> nodes )
    
  if( PhysicalModelShields_node )
  {
    const ::rapidxml::xml_node<char> *self_node = XML_FIRST_NODE(PhysicalModelShields_node, "SelfAttenShield");
    const ::rapidxml::xml_node<char> *self_shield_node = self_node ? XML_FIRST_NODE(self_node, "RelEffShield") : nullptr;
    if( self_shield_node )
    {
      auto shield_state = make_unique<RelEffShieldState>();
      shield_state->fromXml( self_shield_node );
      state.m_selfAttenShield = std::move(shield_state);
    }
    
    const ::rapidxml::xml_node<char> *ext_node = XML_FIRST_NODE(PhysicalModelShields_node, "ExternalAttenShields");
    if( ext_node )
    {
      XML_FOREACH_CHILD( shield_node, ext_node, "RelEffShield" )
      {
        auto shield_state = make_unique<RelEffShieldState>();
        shield_state->fromXml( shield_node );
        state.m_externalShields.push_back( std::move(shield_state) );
      }
    }//if( ext_node )
  }//if( PhysicalModelShields_node )
}//void GuiState::deSerialize( const ::rapidxml::xml_node<char> *base_node )


bool RelActManualGui::GuiState::operator==( const RelActManualGui::GuiState &rhs ) const
{
  return (m_relEffEqnFormIndex == rhs.m_relEffEqnFormIndex)
    && (m_relEffEqnOrderIndex == rhs.m_relEffEqnOrderIndex)
    && (m_physModelUseHoerl == rhs.m_physModelUseHoerl)
    && (m_nucDataSrcIndex == rhs.m_nucDataSrcIndex)
    && (m_matchToleranceValue == rhs.m_matchToleranceValue)
    && (m_addUncertIndex == rhs.m_addUncertIndex)
    && (m_backgroundSubtract == rhs.m_backgroundSubtract)
    && (m_resultTab == rhs.m_resultTab)
    //&& (m_nucAgesAndDecayCorrect == rhs.m_nucAgesAndDecayCorrect)
    && (nucAge == rhs.nucAge)
    && (nucDecayCorrect == rhs.nucDecayCorrect)
  ;
}//RelActManualGui::GuiState::GuiState::operator==
