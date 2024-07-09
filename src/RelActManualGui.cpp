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

#include <Wt/WMenu>
#include <Wt/WLabel>
#include <Wt/WPanel>
#include <Wt/WTable>
#include <Wt/WLabel>
#include <Wt/WServer>
#include <Wt/WComboBox>
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
#include <Wt/WRegExpValidator>
#include <Wt/WContainerWidget>

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
#include "InterSpec/RelActCalc.h"
#include "InterSpec/RelEffChart.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/RelActManualGui.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/IsotopeSelectionAids.h"
#include "InterSpec/PhysicalUnitsLocalized.h"

using namespace Wt;
using namespace std;

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
      const bool useCurie = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
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
    m_decay_during_meas->checked().connect( this, &ManRelEffNucDisp::handleDecayDuringMeasurementChanged );
    m_decay_during_meas->unChecked().connect( this, &ManRelEffNucDisp::handleDecayDuringMeasurementChanged );
    
    const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );
    
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
  m_nucDecayCorrect{},
  m_resultMenu( nullptr ),
  m_chart( nullptr ),
  m_results( nullptr )
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
  WLabel *label = new WLabel( WString::tr("ramg-eqn-form-label"), optionsList->elementAt(0, 0) );
  
  m_relEffEqnForm = new WComboBox( optionsList->elementAt(0, 1) );
  m_relEffEqnForm->activated().connect( this, &RelActManualGui::relEffEqnFormChanged );
  
  HelpSystem::attachToolTipOn( {optionsList->elementAt(0,0), optionsList->elementAt(0,1)},
                              WString::tr("ramg-tt-eqn-form"), showToolTips );
  
  
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
    }//switch( eqn_form )
    
    m_relEffEqnForm->addItem( txt );
  }//for( loop over RelEffEqnForm )
  
  m_relEffEqnForm->setCurrentIndex( static_cast<int>(RelActCalc::RelEffEqnForm::LnX) );
  
  label = new WLabel( WString::tr("ramg-eqn-order-label"), optionsList->elementAt(1, 0) );
  
  m_relEffEqnOrder = new WComboBox( optionsList->elementAt(1, 1) );
  m_relEffEqnOrder->activated().connect( this, &RelActManualGui::relEffEqnOrderChanged );
  
  m_relEffEqnOrder->addItem( "0" );
  m_relEffEqnOrder->addItem( "1" );
  m_relEffEqnOrder->addItem( "2" );
  m_relEffEqnOrder->addItem( "3" );
  m_relEffEqnOrder->addItem( "4" );
  m_relEffEqnOrder->addItem( "5" );
  m_relEffEqnOrder->addItem( "6" );
  m_relEffEqnOrder->setCurrentIndex( 3 );
  
  HelpSystem::attachToolTipOn( {optionsList->elementAt(1, 0),optionsList->elementAt(1, 1)},
                              WString::tr("ramg-tt-eqn-order"), showToolTips );
  
  label = new WLabel( WString::tr("ramg-yield-info-label"), optionsList->elementAt(2, 0) );
  m_nucDataSrc = new WComboBox( optionsList->elementAt(2, 1) );
  label->setBuddy( m_nucDataSrc );
  m_nucDataSrc->activated().connect( this, &RelActManualGui::nucDataSrcChanged );
  
  HelpSystem::attachToolTipOn( {optionsList->elementAt(2, 0),optionsList->elementAt(2, 1)},
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
  
  m_nucDataSrcHolder = optionsList->rowAt(2);

  label = new WLabel( WString::tr("ramg-match-tol-label"), optionsList->elementAt(3, 0) ); //(FWHM)
  m_matchTolerance = new NativeFloatSpinBox( optionsList->elementAt(3, 1) );
  label->setBuddy( m_matchTolerance );
  m_matchTolerance->setSpinnerHidden();
  m_matchTolerance->setWidth( 35 );
  m_matchTolerance->setRange( 0, 5 );
  m_matchTolerance->setValue( 0.5 ); //Other places we use 1.25/2.355 = 0.530786
  label = new WLabel( WString("&nbsp;{1}").arg(WString::tr("FWHM")), optionsList->elementAt(3, 1) );
  label->setBuddy( m_matchTolerance );
  m_matchTolerance->valueChanged().connect( this, &RelActManualGui::matchToleranceChanged );
  
  
  HelpSystem::attachToolTipOn( {optionsList->elementAt(3, 0),optionsList->elementAt(3, 1)},
                              WString::tr("ramg-tt-match-tol"), showToolTips );
  
  label = new WLabel( WString::tr("ramg-add-uncert-label"), optionsList->elementAt(4, 0) );
  m_addUncertainty = new WComboBox( optionsList->elementAt(4, 1) );
  label->setBuddy( m_addUncertainty );
  
  HelpSystem::attachToolTipOn( {optionsList->elementAt(4, 0),optionsList->elementAt(4, 1)},
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
  
  
  m_backgroundSubtract = new WCheckBox( WString::tr("ramg-back-sub-cb"), optionsList->elementAt(5, 0) );
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
    android_download_workaround( m_calpResource, "rel_eff.html");
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
  shared_ptr<GuiState> state = make_shared<GuiState>();
  
  state->m_relEffEqnFormIndex = m_relEffEqnForm->currentIndex();
  state->m_relEffEqnOrderIndex = m_relEffEqnOrder->currentIndex();
  state->m_nucDataSrcIndex = m_nucDataSrc->currentIndex();
  state->m_matchToleranceValue = m_matchTolerance->value();
  state->m_addUncertIndex = m_addUncertainty->currentIndex();
  state->m_backgroundSubtract = m_backgroundSubtract->isChecked();
  state->m_resultTab = m_resultMenu->currentIndex();
  
  for( auto w : m_nuclidesDisp->children() )
  {
    ManRelEffNucDisp *rr = dynamic_cast<ManRelEffNucDisp *>(w);
    if( rr && rr->m_nuc )
      state->m_nucAgesAndDecayCorrect.emplace_back( rr->m_nuc->symbol, rr->m_current_age,
                                                    rr->decayDuringMeasurement() );
  }//for( auto w : m_nuclidesDisp->children() )
  
  return state;
}//std::shared_ptr<const GuiState> getGuiState() const


void RelActManualGui::setGuiState( const GuiState &state )
{
  bool updateCalc = false;
  if( m_relEffEqnForm->currentIndex() != state.m_relEffEqnFormIndex )
  {
    updateCalc = true;
    m_relEffEqnForm->setCurrentIndex( state.m_relEffEqnFormIndex );
  }
  
  if( m_relEffEqnOrder->currentIndex() != state.m_relEffEqnOrderIndex )
  {
    updateCalc = true;
    m_relEffEqnOrder->setCurrentIndex( state.m_relEffEqnOrderIndex );
  }
  
  if( m_nucDataSrc->currentIndex() != state.m_nucDataSrcIndex )
  {
    updateCalc = true;
    m_nucDataSrc->setCurrentIndex( state.m_nucDataSrcIndex );
  }

  if( m_matchTolerance->value() != state.m_matchToleranceValue )
  {
    updateCalc = true;
    m_matchTolerance->setValue( state.m_matchToleranceValue );
  }
  
  if( m_addUncertainty->currentIndex() != state.m_addUncertIndex )
  {
    updateCalc = true;
    m_addUncertainty->setCurrentIndex( state.m_addUncertIndex );
  }
  
  if( m_backgroundSubtract->isChecked() != state.m_backgroundSubtract )
  {
    updateCalc = true;
    m_backgroundSubtract->setChecked( state.m_backgroundSubtract );
  }
  
  if( m_resultMenu->currentIndex() != state.m_resultTab )
    m_resultMenu->select( state.m_resultTab );
  
  
  for( auto w : m_nuclidesDisp->children() )
  {
    ManRelEffNucDisp *rr = dynamic_cast<ManRelEffNucDisp *>(w);
    if( !rr || !rr->m_nuc )
      continue;
    
    for( const auto &i : state.m_nucAgesAndDecayCorrect )
    {
      if( std::get<0>(i) == rr->m_nuc->symbol )
      {
        if( std::get<1>(i) != rr->m_current_age )
        {
          rr->setAge( std::get<1>(i) );
          m_nucAge[std::get<0>(i)] = std::get<1>(i);
          m_renderFlags |= RenderActions::UpdateNuclides;
          updateCalc = true;
        }//if( i.second != rr->m_current_age )
        
        const bool decayCorr = rr->decayDuringMeasurement();
        if( std::get<2>(i) != decayCorr )
        {
          rr->setDecayDuringMeasurement( std::get<2>(i) );
          m_nucDecayCorrect[std::get<0>(i)] = std::get<2>(i);
          m_renderFlags |= RenderActions::UpdateNuclides;
          updateCalc = true;
        }
        
        break;
      }//if( i.first == rr->m_nuc->symbol )
    }//for( loop over previous ages )
  }//for( auto w : m_nuclidesDisp->children() )
  
  if( updateCalc )
  {
    m_renderFlags |= RenderActions::UpdateCalc;
    scheduleRender();
  }
  
  assert( !m_renderFlags.testFlag(RenderActions::AddUndoRedoStep) );
}//void setGuiState( const GuiState &state )


bool RelActManualGui::GuiState::operator==( const RelActManualGui::GuiState &rhs ) const
{
  return (m_relEffEqnFormIndex == rhs.m_relEffEqnFormIndex)
    && (m_relEffEqnOrderIndex == rhs.m_relEffEqnOrderIndex)
    && (m_nucDataSrcIndex == rhs.m_nucDataSrcIndex)
    && (m_matchToleranceValue == rhs.m_matchToleranceValue)
    && (m_addUncertIndex == rhs.m_addUncertIndex)
    && (m_backgroundSubtract == rhs.m_backgroundSubtract)
    && (m_resultTab == rhs.m_resultTab)
    && (m_nucAgesAndDecayCorrect == rhs.m_nucAgesAndDecayCorrect);
}//RelActManualGui::GuiState::operator==


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
  if( m_renderFlags.testFlag(RelActManualGui::RenderActions::AddUndoRedoStep ) )
  {
    addUndoRedoStep( current_gui_state );
  }
  m_currentGuiState = current_gui_state;
  
  
  m_renderFlags = 0;
  
  WContainerWidget::render( flags );
}//render( Wt::WFlags<Wt::RenderFlag> )


void RelActManualGui::calculateSolution()
{
  m_currentSolution.reset();
  m_chart->setData( vector<RelActCalcManual::GenericPeakInfo>{}, {}, "" );
  m_results->clear();
  
  // TODO: should do the actual computation not on the GUI thread!
  try
  {
    using namespace RelActCalcManual;
    bool has_reaction = false;
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
    map<const SandiaDecay::Nuclide *,bool> decay_correct_during_meas;
    for( auto w : m_nuclidesDisp->children() )
    {
      ManRelEffNucDisp *rr = dynamic_cast<ManRelEffNucDisp *>(w);
      if( rr && rr->m_nuc )
      {
        nuclide_ages[rr->m_nuc] = rr->m_current_age;
        decay_correct_during_meas[rr->m_nuc] = rr->decayDuringMeasurement();
      }
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
      if( p && (p->parentNuclide() || p->reaction()) && p->useForManualRelEff() )
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
            snprintf( buffer, sizeof(buffer), "%.2f", peak.m_energy );
            
            prep_warnings.push_back( WString::tr("ramg-back-sub-neg").arg( buffer ).toUTF8() );
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
            assert( decay_corr_pos != end(decay_correct_during_meas) );
            
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
        }else
        {
          assert( p->reaction() );
          unique_isotopes.insert( p->reaction()->name() );
        }
      }//
    }//for( const PeakModel::PeakShrdPtr &p : *m_peakModel->peaks() )
    
    if( has_reaction )
    {
      prep_warnings.push_back( WString::tr("ramg-warn-reaction").toUTF8() );
    }//if( user is using reactions )
    
    if( background_sub && !num_peaks_back_sub )
    {
      prep_warnings.push_back( WString::tr("ramg-warn-no-bkg-sub-used").toUTF8() );
    }//if( user wanted to background subtract peaks, but no peaks matched up )
    
    if( has_U_or_Pu && (lowest_energy_peak < 122) )
    {
      prep_warnings.push_back( WString::tr("ramg-warn-rel-eff-u/pu-xray").toUTF8() );
    }else if( lowest_energy_peak < 90 )
    {
      prep_warnings.push_back( WString::tr("ramg-warn-rel-eff-other-xray").toUTF8() );
    }
    
    
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

      prep_warnings.push_back( msg.toUTF8() );
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
    
    // Note peaks that have been corrected for decay during measurement; only peaks with corrections
    //  will get put in here - and will include info for only gamma lines that have been corrected.
    vector<RelActCalcManual::GenericPeakInfo> pre_decay_correction_info;
    
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
      
      pre_decay_correction_info = matched_res.not_decay_corrected_peaks;
    }// End code to fill in nuclide information
    

    // We'll do the actual calculation off of the main thread; in order to make sure the widget
    //  still exists at the end of computations, we'll use WApplication::bind(), in combination
    //  with shared_ptrs's to make sure everything is okay
    const RelActCalc::RelEffEqnForm eqn_form = relEffEqnForm();
    const size_t eqn_order = relEffEqnOrder();
    
    const string sessionId = wApp->sessionId();
    auto solution = make_shared<RelActCalcManual::RelEffSolution>();
    auto updater = wApp->bind( boost::bind( &RelActManualGui::updateGuiWithResults, this, solution ) );
    
    auto errmsg = make_shared<WString>();
    auto err_updater = wApp->bind( boost::bind( &RelActManualGui::updateGuiWithError, this, errmsg ) );
    
    WServer::instance()->ioService().boost::asio::io_service::post( std::bind(
      [peak_infos, pre_decay_correction_info, eqn_form, eqn_order, sessionId, solution, updater, prep_warnings, errmsg, err_updater](){
        try
        {
          *solution = solve_relative_efficiency( peak_infos, eqn_form, eqn_order );
          solution->m_warnings.insert(begin(solution->m_warnings), begin(prep_warnings), end(prep_warnings));
          solution->m_input_peaks_before_decay_corr = pre_decay_correction_info;
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
      m_chart->setData( vector<RelActCalcManual::GenericPeakInfo>{}, {}, "" );
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
  
  solution.get_mass_fraction_table( results_html );
  results_html << "<br />";
  solution.get_mass_ratio_table( results_html );
  results_html << "<br />";
  
  results_html << "<div class=\"releffeqn\">Rel. Eff.: y = "
  << RelActCalc::rel_eff_eqn_text( solution.m_rel_eff_eqn_form, solution.m_rel_eff_eqn_coefficients )
  << "</div>\n";
  
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


void RelActManualGui::relEffEqnFormChanged()
{
  m_renderFlags |= RenderActions::UpdateCalc;
  m_renderFlags |= RenderActions::AddUndoRedoStep;
  scheduleRender();
}//void relEffEqnFormChanged()


void RelActManualGui::relEffEqnOrderChanged()
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
        continue;
      
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
    assert( rr );
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
  ::rapidxml::xml_document<char> *doc = parent_node->document();
  assert( doc );
  
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
  
  const RelActCalc::RelEffEqnForm eqn_form = relEffEqnForm();
  value = RelActCalc::to_str( eqn_form );
  ::rapidxml::xml_node<char> *node = doc->allocate_node( ::rapidxml::node_element, "RelEffEqnForm", value );
  base_node->append_node( node );
  
  const size_t eqn_order = relEffEqnOrder();
  value = doc->allocate_string( std::to_string(eqn_order).c_str() );
  node = doc->allocate_node( ::rapidxml::node_element, "RelEffEqnOrder", value );
  base_node->append_node( node );
  
  const RelActCalcManual::PeakCsvInput::NucDataSrc srcData = nucDataSrc();
  
  value = RelActCalcManual::PeakCsvInput::to_str( srcData );
  node = doc->allocate_node( ::rapidxml::node_element, "NucDataSrc", value );
  base_node->append_node( node );
  
  const float match_tolerance = m_matchTolerance->value();
  const string match_tolerance_str = SpecUtils::printCompact(match_tolerance, 7);
  value = doc->allocate_string( match_tolerance_str.c_str() );
  node = doc->allocate_node( ::rapidxml::node_element, "MatchTolerance", value );
  base_node->append_node( node );
  
  value = RelActManualGui::to_str( AddUncert(m_addUncertainty->currentIndex()) );
  node = doc->allocate_node( ::rapidxml::node_element, "AddUncertainty", value );
  base_node->append_node( node );
  
  
  value = m_resultMenu->currentIndex() ? "1" : "0";
  node = doc->allocate_node( ::rapidxml::node_element, "ResultTabShowing", value );
  base_node->append_node( node );
  
  
  map<string,double> nuc_age_cache = m_nucAge;
  map<string,bool> nuc_decay_correct = m_nucDecayCorrect;
  
  for( auto w : m_nuclidesDisp->children() )
  {
    const ManRelEffNucDisp *rr = dynamic_cast<const ManRelEffNucDisp *>(w);
    if( rr && rr->m_nuc )
    {
      nuc_age_cache[rr->m_nuc->symbol] = rr->m_current_age;
      
      if( rr->m_decay_during_meas )
        nuc_decay_correct[rr->m_nuc->symbol] = rr->m_decay_during_meas->isChecked();
      else
        nuc_decay_correct.erase( rr->m_nuc->symbol );
    }
  }//for( auto w : m_nuclidesDisp->children() )
  
  
  ::rapidxml::xml_node<char> *nuc_ages_node = doc->allocate_node( ::rapidxml::node_element, "NuclideAges" );
  base_node->append_node( nuc_ages_node );
  for( const auto &n : nuc_age_cache )
  {
    ::rapidxml::xml_node<char> *nuc_node = doc->allocate_node( ::rapidxml::node_element, "Nuclide" );
    nuc_ages_node->append_node( nuc_node );
    
    value = doc->allocate_string( n.first.c_str() );
    ::rapidxml::xml_node<char> *name_node = doc->allocate_node( ::rapidxml::node_element, "Name", value );
    nuc_node->append_node(name_node);
    
    value = doc->allocate_string( SpecUtils::printCompact(n.second,8).c_str() );
    ::rapidxml::xml_node<char> *age_node = doc->allocate_node( ::rapidxml::node_element, "Age", value );
    nuc_node->append_node(age_node);
    
    const auto decay_corr_pos = nuc_decay_correct.find(n.first);
    if( decay_corr_pos != end(nuc_decay_correct) )
    {
      value = decay_corr_pos->second ? "true" : "false";
      ::rapidxml::xml_node<char> *corr_node = doc->allocate_node( ::rapidxml::node_element, "DecayDuringMeasurement", value );
      nuc_node->append_node(corr_node);
    }
  }//for( const auto &n : nuc_age_cache )
  
  return base_node;
}//serialize(...)


void RelActManualGui::deSerialize( const ::rapidxml::xml_node<char> *base_node )
{
  if( SpecUtils::xml_name_str(base_node) != "RelActManualGui" )
    throw runtime_error( "RelActManualGui::deSerialize: invalid base node passed in: '"
                        + SpecUtils::xml_name_str(base_node) + "'" );
    
  int version;
  const ::rapidxml::xml_attribute<char> *attr = XML_FIRST_ATTRIB(base_node, "version");
  if( !attr || !attr->value() || !(stringstream(attr->value()) >> version) )
    throw runtime_error( "Deserializing requires a version" );
  
  if( version != sm_xmlSerializationMajorVersion )
    throw runtime_error( "Invalid version of RelActManualGui XML" );
  
  m_nuclidesDisp->clear();
  //m_nucAge.clear();
  
  const ::rapidxml::xml_node<char> *RelEffEqnForm_node = XML_FIRST_NODE(base_node, "RelEffEqnForm");
  const ::rapidxml::xml_node<char> *RelEffEqnOrder_node = XML_FIRST_NODE(base_node, "RelEffEqnOrder");
  const ::rapidxml::xml_node<char> *NucDataSrc_node = XML_FIRST_NODE(base_node, "NucDataSrc");
  const ::rapidxml::xml_node<char> *MatchTolerance_node = XML_FIRST_NODE(base_node, "MatchTolerance");
  const ::rapidxml::xml_node<char> *AddUncertainty_node = XML_FIRST_NODE(base_node, "AddUncertainty");
  const ::rapidxml::xml_node<char> *ResultTabShowing_node = XML_FIRST_NODE(base_node, "ResultTabShowing");
  const ::rapidxml::xml_node<char> *NuclideAges_node = XML_FIRST_NODE(base_node, "NuclideAges");
  
  
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
  
  
  int eqn_order = -1;
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
    
    m_nucAge[nuc] = age;
    
    if( decay_corr_node )
    {
      const string decay_corr = SpecUtils::xml_value_str(decay_corr_node);
      m_nucDecayCorrect[nuc] = ((decay_corr == "true") || (decay_corr == "1"));
    }
  }//for( loop over <Nuclide> nodes )
  
  
  // Now need to set state of widgets
  m_relEffEqnForm->setCurrentIndex( static_cast<int>(eqn_form) );
  m_relEffEqnOrder->setCurrentIndex( eqn_order );
  m_nucDataSrc->setCurrentIndex( static_cast<int>(data_src) );
  m_matchTolerance->setValue( match_tolerance );
  m_addUncertainty->setCurrentIndex( static_cast<int>(add_uncert) );
  m_resultMenu->select( tab_showing );
  
  // Schedule calc/render
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
