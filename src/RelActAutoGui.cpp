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

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"
#include "rapidxml/rapidxml_utils.hpp"

#include <Wt/WMenu>
#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WPoint>
#include <Wt/WServer>
#include <Wt/WCheckBox>
#include <Wt/WComboBox>
#include <Wt/WGroupBox>
#include <Wt/WMenuItem>
#include <Wt/WResource>
#include <Wt/WIOService>
#include <Wt/WFileUpload>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WInPlaceEdit>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WStackedWidget>
#include <Wt/WContainerWidget>
#include <Wt/WRegExpValidator>
#include <Wt/WSuggestionPopup>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/PeakFit.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/ColorSelect.h"
#include "InterSpec/RelEffChart.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/EnergyCalTool.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/RelActAutoGui.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/SwitchCheckbox.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/RelActTxtResults.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/RelEffShieldWidget.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeNameFilterModel.h"
#include "InterSpec/PhysicalUnitsLocalized.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"
#include "InterSpec/RelActAutoGuiRelEffOptions.h"


using namespace Wt;
using namespace std;

#if( ANDROID )
// Defined in target/android/android.cpp
extern void android_download_workaround( Wt::WResource *resource, std::string description );
#endif

namespace
{
  //DeleteOnClosePopupMenu - same class as from D3SpectrumDisplayDiv... should refactor
  class DeleteOnClosePopupMenu : public PopupDivMenu
  {
    bool m_deleteWhenHidden;
  public:
    DeleteOnClosePopupMenu( WPushButton *p, const PopupDivMenu::MenuType t )
      : PopupDivMenu( p, t ), m_deleteWhenHidden( false ) {}
    virtual ~DeleteOnClosePopupMenu(){}
    void markForDelete(){ m_deleteWhenHidden = true; }
    virtual void setHidden( bool hidden, const WAnimation &a = WAnimation() )
    {
      PopupDivMenu::setHidden( hidden, a );
      if( hidden && m_deleteWhenHidden )
        delete this;
    }
  };//class PeakRangePopupMenu


  class RelActAutoReportResource : public Wt::WResource
  {
    Wt::WApplication *m_app;
    RelActAutoGui *m_tool;
    std::shared_ptr<const RelActCalcAuto::RelActAutoSolution> m_solution;
    
  public:
    RelActAutoReportResource( RelActAutoGui *tool, WObject* parent = nullptr )
    : WResource( parent ), m_app( WApplication::instance() ), m_tool( tool ), m_solution( nullptr )
    {
      assert( m_app );
      assert( m_tool );
    }
  
    virtual ~RelActAutoReportResource()
    {
      beingDeleted();
    }
  
    void updateSolution( const shared_ptr<const RelActCalcAuto::RelActAutoSolution> &solution )
    {
      m_solution = solution;
    }
    
    virtual void handleRequest( const Wt::Http::Request &request, Wt::Http::Response &response )
    {
      assert( m_app );
      
      try
      {
        WApplication::UpdateLock lock( m_app );
        
        if( !lock )
          throw std::runtime_error( "Error grabbing application lock to from RelActAutoReportResource resource." );
        
        //const shared_ptr<const RelActCalcAuto::RelActAutoSolution> solution = m_tool->getCurrentSolution();
        
        if( !m_solution )
        {
          response.out() << "<!DOCTYPE html>\n"
          "\t<head><meta charset=\"utf-8\"><title>No <em>Isotopics by nuclide</em> solution available</title></head>"
          "\t<body>"
          "\t\tSorry - no solution currently available."
          "\t</body>"
          "</html>";
          
          return;
        }//if( !m_solution )
        
        string filename;
        InterSpec *viewer = InterSpec::instance();
        shared_ptr<SpecMeas> meas = viewer ? viewer->measurment(SpecUtils::SpectrumType::Foreground) : nullptr;
        
        if( meas )
          filename = meas->filename();
        
        if( filename.empty() )
          filename = "rel_act_from_nuc";
        const string orig_extension = SpecUtils::file_extension(filename);
        if( orig_extension.size() && (orig_extension.size() < filename.size()) )
          filename = filename.substr(0,filename.size() - orig_extension.size());
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
              
        m_solution->print_html_report( response.out() );
      }catch( std::exception &e )
      {
        log("error") << "Error handling request for RelActAutoReportResource: " << e.what();
        response.out() << "Error creating HTML file: " << e.what()
        << "\n\nPlease report to InterSpec@sandia.gov.";
      }//try / catch
    }//void handleRequest(...)
  };//class RelActAutoReportResource


  class RelActAutoParamsResource : public Wt::WResource
  {
    Wt::WApplication *m_app;
    RelActAutoGui *m_tool;
    
  public:
    RelActAutoParamsResource( RelActAutoGui *tool, WObject* parent = nullptr )
    : WResource( parent ), m_app( WApplication::instance() ), m_tool( tool )
    {
      assert( m_app );
      assert( m_tool );
    }
    
    virtual ~RelActAutoParamsResource()
    {
      beingDeleted();
    }
    
    virtual void handleRequest( const Wt::Http::Request &request, Wt::Http::Response &response )
    {
      assert( m_app );
      
      try
      {
        WApplication::UpdateLock lock( m_app );
        
        if( !lock )
          throw std::runtime_error( "Error grabbing application lock to from RelActAutoReportResource resource." );
        
        //const shared_ptr<const RelActCalcAuto::RelActAutoSolution> solution = m_tool->getCurrentSolution();
        string filename = "isotopics_by_nuclides";
        InterSpec *viewer = InterSpec::instance();
        shared_ptr<SpecMeas> meas = viewer ? viewer->measurment(SpecUtils::SpectrumType::Foreground) : nullptr;
        
        if( meas && !meas->filename().empty() )
          filename += "_" + meas->filename();
        
        const string orig_extension = SpecUtils::file_extension(filename);
        if( orig_extension.size() && (orig_extension.size() < filename.size()) )
          filename = filename.substr(0,filename.size() - orig_extension.size());
        filename += "_releff.xml";
        
        suggestFileName( filename, WResource::Attachment );
        response.setMimeType( "application/xml" );
        
        std::unique_ptr<rapidxml::xml_document<char>> xml = m_tool->guiStateToXml();
        
        if( !xml )
        {
          response.out() << "Error getting XML state.\n";
          return;
        }
        
        rapidxml::print( response.out(), *xml, 0 );
      }catch( std::exception &e )
      {
        log("error") << "Error handling request for RelActAutoReportResource: " << e.what();
        response.out() << "Error creating XML parameter file: " << e.what()
        << "\n\nPlease report to InterSpec@sandia.gov.";
      }//try / catch
    }//void handleRequest(...)
  };//class RelActAutoReportResource


  class RelActAutoEnergyRange : public WContainerWidget
  {
    RelActAutoGui *m_gui;
    Wt::Signal<> m_updated;
    Wt::Signal<> m_remove_energy_range;
    
    NativeFloatSpinBox *m_lower_energy;
    NativeFloatSpinBox *m_upper_energy;
    WComboBox *m_continuum_type;
    WCheckBox *m_force_full_range;
    WPushButton *m_to_individual_rois;
    
    /// Used to track the Highlight region this energy region corresponds to in D3SpectrumDisplayDiv
    size_t m_highlight_region_id;
    
  public:
    RelActAutoEnergyRange( RelActAutoGui *gui, WContainerWidget *parent = nullptr )
    : WContainerWidget( parent ),
      m_gui( gui ),
      m_updated( this ),
      m_remove_energy_range( this ),
      m_lower_energy( nullptr ),
      m_upper_energy( nullptr ),
      m_continuum_type( nullptr ),
      m_force_full_range( nullptr ),
      m_to_individual_rois( nullptr ),
      m_highlight_region_id( 0 )
    {
      addStyleClass( "RelActAutoEnergyRange" );
      
      wApp->useStyleSheet( "InterSpec_resources/GridLayoutHelpers.css" );
      
      const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );
      
      WLabel *label = new WLabel( "Lower Energy", this );
      label->addStyleClass( "GridFirstCol GridFirstRow" );
      
      m_lower_energy = new NativeFloatSpinBox( this );
      m_lower_energy->addStyleClass( "GridSecondCol GridFirstRow" );
      m_lower_energy->setSpinnerHidden( true );
      label->setBuddy( m_lower_energy );
      
      label = new WLabel( "keV", this );
      label->addStyleClass( "GridThirdCol GridFirstRow" );
      
      label = new WLabel( "Upper Energy", this );
      label->addStyleClass( "GridFirstCol GridSecondRow" );
      
      m_upper_energy = new NativeFloatSpinBox( this );
      m_upper_energy->addStyleClass( "GridSecondCol GridSecondRow" );
      m_upper_energy->setSpinnerHidden( true );
      label->setBuddy( m_upper_energy );
      
      label = new WLabel( "keV", this );
      label->addStyleClass( "GridThirdCol GridSecondRow" );
      
      m_lower_energy->valueChanged().connect( this, &RelActAutoEnergyRange::handleEnergyChange );
      m_upper_energy->valueChanged().connect( this, &RelActAutoEnergyRange::handleEnergyChange );
      
      label = new WLabel( "Continuum Type:", this );
      label->addStyleClass( "GridFourthCol GridFirstRow" );
      m_continuum_type = new WComboBox( this );
      m_continuum_type->addStyleClass( "GridFifthCol GridFirstRow" );
      label->setBuddy( m_continuum_type );
      
      // We wont allow "External" here
      for( int i = 0; i < static_cast<int>(PeakContinuum::OffsetType::External); ++i )
      {
        const char *key = PeakContinuum::offset_type_label_tr( PeakContinuum::OffsetType(i) );
        m_continuum_type->addItem( WString::tr(key) );
      }//for( loop over PeakContinuum::OffsetType )
      
      m_continuum_type->setCurrentIndex( static_cast<int>(PeakContinuum::OffsetType::Linear) );
      m_continuum_type->changed().connect( this, &RelActAutoEnergyRange::handleContinuumTypeChange );

      m_force_full_range = new WCheckBox( "Force full-range", this );
      m_force_full_range->addStyleClass( "GridFourthCol GridSecondRow GridSpanTwoCol" );
      m_force_full_range->checked().connect( this, &RelActAutoEnergyRange::handleForceFullRangeChange );
      m_force_full_range->unChecked().connect( this, &RelActAutoEnergyRange::handleForceFullRangeChange );
      
      WPushButton *removeEnergyRange = new WPushButton( this );
      removeEnergyRange->setStyleClass( "DeleteEnergyRangeOrNuc GridSixthCol GridFirstRow Wt-icon" );
      removeEnergyRange->setIcon( "InterSpec_resources/images/minus_min_black.svg" );
      removeEnergyRange->clicked().connect( this, &RelActAutoEnergyRange::handleRemoveSelf );
      
      m_to_individual_rois = new WPushButton( this );
      m_to_individual_rois->setStyleClass( "ToIndividualRois GridSixthCol GridSecondRow Wt-icon" );
      m_to_individual_rois->setIcon( "InterSpec_resources/images/expand_list.svg" );
      m_to_individual_rois->clicked().connect( boost::bind( &RelActAutoGui::handleConvertEnergyRangeToIndividuals,
                                                           m_gui, static_cast<WWidget *>(this) ) );
      m_to_individual_rois->setHidden( m_force_full_range->isChecked() );
      
      const char *tooltip = "Converts this energy range into individual ROIs, based on which gammas should be"
      " grouped together vs apart (e.g., leaves gammas with overlapping peaks in a single ROI, while if two gammas"
      " are far apart, they will be split into seperate ROIs.";
      HelpSystem::attachToolTipOn( m_to_individual_rois, tooltip, showToolTips );
    }//RelActAutoEnergyRange constructor
    
    
    bool isEmpty() const
    {
      return ((fabs(m_lower_energy->value() - m_upper_energy->value() ) < 1.0)
              || (m_upper_energy->value() <= 0.0));
    }
    
    
    void handleRemoveSelf()
    {
      m_remove_energy_range.emit();
    }//void handleRemoveSelf()
    
    
    void handleContinuumTypeChange()
    {
      m_updated.emit();
    }
    
    void handleForceFullRangeChange()
    {
      m_to_individual_rois->setHidden( m_to_individual_rois->isDisabled() || m_force_full_range->isChecked() );
      
      m_updated.emit();
    }
    
    
    void enableSplitToIndividualRanges( const bool enable )
    {
      m_to_individual_rois->setHidden( !enable || m_force_full_range->isChecked() );
      m_to_individual_rois->setEnabled( enable );
    }
    
    
    void handleEnergyChange()
    {
      float lower = m_lower_energy->value();
      float upper = m_upper_energy->value();
      if( lower > upper )
      {
        m_lower_energy->setValue( upper );
        m_upper_energy->setValue( lower );
        
        std::swap( lower, upper );
      }//if( lower > upper )
      
      m_updated.emit();
    }//void handleEnergyChange()
    
    
    void setEnergyRange( float lower, float upper )
    {
      if( lower > upper )
        std::swap( lower, upper );
      
      m_lower_energy->setValue( lower );
      m_upper_energy->setValue( upper );
      
      m_updated.emit();
    }//void setEnergyRange( float lower, float upper )
    
    
    bool forceFullRange() const
    {
      return m_force_full_range->isChecked();
    }
    
    
    void setForceFullRange( const bool force_full )
    {
      if( force_full == m_force_full_range->isChecked() )
        return;
      
      m_force_full_range->setChecked( force_full );
      m_updated.emit();
    }
    
    
    void setContinuumType( const PeakContinuum::OffsetType type )
    {
      const int type_index = static_cast<int>( type );
      if( (type_index < 0)
         || (type_index >= static_cast<int>(PeakContinuum::OffsetType::External)) )
      {
        assert( 0 );
        return;
      }
      
      if( type_index == m_continuum_type->currentIndex() )
        return;
      
      m_continuum_type->setCurrentIndex( type_index );
      m_updated.emit();
    }//void setContinuumType( PeakContinuum::OffsetType type )
    
    
    void setHighlightRegionId( const size_t chart_id )
    {
      m_highlight_region_id = chart_id;
    }
    
    
    size_t highlightRegionId() const
    {
      return m_highlight_region_id;
    }
  
    float lowerEnergy() const
    {
      return m_lower_energy->value();
    }
    
    
    float upperEnergy() const
    {
      return m_upper_energy->value();
    }
    
    
    void setFromRoiRange( const RelActCalcAuto::RoiRange &roi )
    {
      if( roi.continuum_type == PeakContinuum::OffsetType::External )
        throw runtime_error( "setFromRoiRange: External continuum type not supported" );
      
      if( roi.allow_expand_for_peak_width )
        throw runtime_error( "setFromRoiRange: allow_expand_for_peak_width set to true not supported" );
      
      m_lower_energy->setValue( roi.lower_energy );
      m_upper_energy->setValue( roi.upper_energy );
      m_continuum_type->setCurrentIndex( static_cast<int>(roi.continuum_type) );
      m_force_full_range->setChecked( roi.force_full_range );
      
      enableSplitToIndividualRanges( !roi.force_full_range );
    }//setFromRoiRange(...)
    
    
    RelActCalcAuto::RoiRange toRoiRange() const
    {
      RelActCalcAuto::RoiRange roi;
      
      roi.lower_energy = m_lower_energy->value();
      roi.upper_energy = m_upper_energy->value();
      roi.continuum_type = PeakContinuum::OffsetType( m_continuum_type->currentIndex() );
      roi.force_full_range = m_force_full_range->isChecked();
      roi.allow_expand_for_peak_width = false; //not currently supported to the GUI, or tested in the calculation code
      
      return roi;
    }//RelActCalcAuto::RoiRange toRoiRange() const
    
    
    Wt::Signal<> &updated()
    {
      return m_updated;
    }
    
    
    Wt::Signal<> &remove()
    {
      return m_remove_energy_range;
    }
  };//class RelActAutoEnergyRange


  class RelActAutoNuclide : public WContainerWidget
  {
    RelActAutoGui *m_gui;
    
    WLineEdit *m_nuclide_edit;
    WLabel *m_age_label;
    WLineEdit *m_age_edit;
    WCheckBox *m_fit_age;
    ColorSelect *m_color_select;
    
    Wt::Signal<> m_updated;
    Wt::Signal<> m_remove;
    
  public:
    RelActAutoNuclide( RelActAutoGui *gui, WContainerWidget *parent = nullptr )
    : WContainerWidget( parent ),
    m_gui( gui ),
    m_nuclide_edit( nullptr ),
    m_age_label( nullptr ),
    m_age_edit( nullptr ),
    m_fit_age( nullptr ),
    m_color_select( nullptr ),
    m_updated( this ),
    m_remove( this )
    {
      addStyleClass( "RelActAutoNuclide" );
      
      const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );
      
      WLabel *label = new WLabel( "Nuclide:", this );
      m_nuclide_edit = new WLineEdit( "", this );
      
      m_nuclide_edit->setAutoComplete( false );
      m_nuclide_edit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
      m_nuclide_edit->setAttributeValue( "autocorrect", "off" );
      m_nuclide_edit->setAttributeValue( "spellcheck", "off" );
#endif
      label->setBuddy( m_nuclide_edit );
      
      m_nuclide_edit->changed().connect( this, &RelActAutoNuclide::handleIsotopeChange );
      
      // If we are typing in this box, we want to let app-hotkeys propagate up, but not arrow keys and
      //  stuff
      const string jsAppKeyDownFcn = wApp->javaScriptClass() + ".appKeyDown";
      const string keyDownJs = "function(s1,e1){"
      "if(e1 && e1.ctrlKey && e1.key && " + jsAppKeyDownFcn + ")"
      + jsAppKeyDownFcn + "(e1);"
      "}";
      m_nuclide_edit->keyWentDown().connect( keyDownJs );
      
      const char *tooltip = "ex. <b>U235</b>, <b>235 Uranium</b>"
      ", <b>U-235m</b> (meta stable state)"
      ", <b>Cs137</b>, Pb, Fe(n,n), etc.";
      HelpSystem::attachToolTipOn( m_nuclide_edit, tooltip, showToolTips );
      
      string replacerJs, matcherJs;
      IsotopeNameFilterModel::replacerJs( replacerJs );
      IsotopeNameFilterModel::nuclideNameMatcherJs( matcherJs );
      IsotopeNameFilterModel *isoSuggestModel = new IsotopeNameFilterModel( this );
      isoSuggestModel->excludeXrays( false );
      isoSuggestModel->excludeEscapes( false );
      isoSuggestModel->excludeReactions( false );
      
      WSuggestionPopup *nuclideSuggest = new WSuggestionPopup( matcherJs, replacerJs, this );
#if( WT_VERSION < 0x3070000 ) //I'm not sure what version of Wt "wtNoReparent" went away.
      nuclideSuggest->setJavaScriptMember("wtNoReparent", "true");
#endif
      nuclideSuggest->setMaximumSize( WLength::Auto, WLength(15, WLength::FontEm) );
      nuclideSuggest->setWidth( WLength(70, Wt::WLength::Unit::Pixel) );
      
      IsotopeNameFilterModel::setQuickTypeFixHackjs( nuclideSuggest );
      
      isoSuggestModel->filter( "" );
      nuclideSuggest->setFilterLength( -1 );
      nuclideSuggest->setModel( isoSuggestModel );
      nuclideSuggest->filterModel().connect( isoSuggestModel, &IsotopeNameFilterModel::filter );
      nuclideSuggest->forEdit( m_nuclide_edit, WSuggestionPopup::Editing );  // | WSuggestionPopup::DropDownIcon
      
      m_age_label = new WLabel( "Age:", this );
      m_age_edit = new WLineEdit( "", this );
      m_age_label->setBuddy( m_age_edit );
      
      WRegExpValidator *validator = new WRegExpValidator( PhysicalUnitsLocalized::timeDurationHalfLiveOptionalRegex(), this );
      validator->setFlags(Wt::MatchCaseInsensitive);
      m_age_edit->setValidator(validator);
      m_age_edit->setAutoComplete( false );
      m_age_edit->setAttributeValue( "ondragstart", "return false" );
      m_age_edit->changed().connect( this, &RelActAutoNuclide::handleAgeChange );
      
      m_fit_age = new WCheckBox( "Fit Age", this );
      m_fit_age->setWordWrap( false );
      m_fit_age->checked().connect( this, &RelActAutoNuclide::handleFitAgeChange );
      m_fit_age->unChecked().connect( this, &RelActAutoNuclide::handleFitAgeChange );
      
      WContainerWidget *spacer = new WContainerWidget( this );
      spacer->addStyleClass( "RelActAutoNuclideSpacer" );
      
      
      m_color_select = new ColorSelect( ColorSelect::PrefferNative, this );
      m_color_select->setColor( WColor("#FF6633") );
      
      m_color_select->cssColorChanged().connect( this, &RelActAutoNuclide::handleColorChange );
      
      
      WPushButton *removeEnergyRange = new WPushButton( this );
      removeEnergyRange->setStyleClass( "DeleteEnergyRangeOrNuc Wt-icon" );
      removeEnergyRange->setIcon( "InterSpec_resources/images/minus_min_black.svg" );
      removeEnergyRange->clicked().connect( this, &RelActAutoNuclide::handleRemoveSelf );
    }//RelActAutoNuclide
    
    
    void handleIsotopeChange()
    {
      const auto nuc_input = nuclide();

      const auto hide_age_stuff = [this,&nuc_input](){
        m_age_label->hide();
        m_age_edit->hide();
        m_fit_age->setUnChecked();
        m_fit_age->hide();

        const string nucstr = m_nuclide_edit->text().toUTF8();
        if( !nucstr.empty() && std::holds_alternative<std::monostate>(nuc_input) )
          passMessage( nucstr + " is not a valid nuclide, x-ray, or reaction.", WarningWidget::WarningMsgHigh );
      };//hide_age_stuff

      if( std::holds_alternative<std::monostate>(nuc_input) )
      {
        hide_age_stuff();
        m_updated.emit();
        return;
      }

      std::string src_name;
      const SandiaDecay::Nuclide *nuc = nullptr;
      const SandiaDecay::Element *el = nullptr;
      const ReactionGamma::Reaction * reaction = nullptr;

      if( std::holds_alternative<const SandiaDecay::Nuclide *>(nuc_input) )
      {
        nuc = std::get<const SandiaDecay::Nuclide *>(nuc_input);

        assert( nuc );
        if( !nuc )
          throw runtime_error( "No valid nuclide" );

        src_name = nuc->symbol;
      
        if( IsInf(nuc->halfLife) )
        {
          const string nucstr = m_nuclide_edit->text().toUTF8();
          passMessage( nucstr + " is a stable nuclide.", WarningWidget::WarningMsgHigh );
        
          m_nuclide_edit->setText( "" );
          m_nuclide_edit->validate();
          m_fit_age->setUnChecked();
          m_fit_age->hide();
        
          m_updated.emit();
          return;
        }//if( IsInf(nuc->halfLife) )
      
        const bool age_is_fittable = !PeakDef::ageFitNotAllowed( nuc );
      
        if( nuc->decaysToStableChildren() )
        {
          m_age_edit->setText( "0y" );
        }else
        {
          string agestr;
          PeakDef::defaultDecayTime( nuc, &agestr );
          m_age_edit->setText( agestr );
        }
      
        m_age_label->setHidden( !age_is_fittable );
        m_age_edit->setHidden( !age_is_fittable );
        m_fit_age->setHidden( !age_is_fittable );
        m_fit_age->setChecked( false );
      }
      
      if( std::holds_alternative<const SandiaDecay::Element *>(nuc_input) )
      {
        el = std::get<const SandiaDecay::Element *>(nuc_input);
        assert( el );
        src_name = el->symbol;
        hide_age_stuff();
      }

      if( std::holds_alternative<const ReactionGamma::Reaction *>(nuc_input) )
      {
        reaction = std::get<const ReactionGamma::Reaction *>(nuc_input);
        assert( reaction );
        src_name = reaction->name();
        hide_age_stuff();
      }
      
      bool haveFoundColor = false;
      
      // Check user is showing reference lines, displayed peaks, and previous user-selected colors
      if( !haveFoundColor )
      {
        const ReferencePhotopeakDisplay *refdisp = InterSpec::instance()->referenceLinesWidget();
        if( refdisp )
        {
          const Wt::WColor c = refdisp->suggestColorForSource( src_name );
          if( !c.isDefault() )
          {
            haveFoundColor = true;
            m_color_select->setColor( c );
          }
        }//if( refdisp )
      }//if( !haveFoundColor )
      
      // Check if the theme explicitly specifies a color for this nuclide
      if( !haveFoundColor )
      {
        shared_ptr<const ColorTheme> theme = InterSpec::instance()->getColorTheme();
        if( theme )
        {
          const map<string,WColor> &defcolors = theme->referenceLineColorForSources;
          const auto pos = defcolors.find( src_name );
          if( pos != end(defcolors) )
            m_color_select->setColor( pos->second );
        }//if( theme )
      }//if( !haveFoundColor )
      
      if( !haveFoundColor )
      {
        // TODO: create a cache of previous user-selected colors
      }//if( !haveFoundColor )
      
      m_updated.emit();
    }//void handleIsotopeChange()
    
    
    void setFitAgeVisible( bool visible, bool do_fit )
    {
      if( visible )
      {
        m_fit_age->setHidden( false );
        m_fit_age->setChecked( do_fit );
      }else
      {
        m_fit_age->setHidden( true );
        m_fit_age->setChecked( false );
      }
    }//setFitAgeVisible(...)
    
    
    void setAgeDisabled( bool disabled )
    {
      m_fit_age->setDisabled( disabled );
      m_age_edit->setDisabled( disabled );
    }//void setAgeDisabled( bool disabled )
    
    
    void setAge( const string &age )
    {
      m_age_edit->setValueText( WString::fromUTF8(age) );
    }
    
    void setNuclideEditFocus()
    {
      m_nuclide_edit->setFocus();
    }
    
    void handleAgeChange()
    {
      m_updated.emit();
    }//void handleAgeChange()
    
    
    void handleFitAgeChange()
    {
      
      m_updated.emit();
    }
    
    
    void handleColorChange()
    {
      m_updated.emit();
    }
    
    /** Will return std::monostate if invalid text entered. */
    std::variant<std::monostate, const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *> nuclide() const
    {
      const string nucstr = m_nuclide_edit->text().toUTF8();
      const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
      assert( db );
      if( !db )
        throw runtime_error( "Couldnt load decay DB" );
      
      const SandiaDecay::Nuclide *nuc = db->nuclide( nucstr );
      if( nuc )
        return nuc;
      
      const SandiaDecay::Element *el = db->element( nucstr );
      if( el )
        return el;
      
      const ReactionGamma * const reaction_db = ReactionGammaServer::database();
      assert( reaction_db );
      if( !reaction_db )
        throw runtime_error( "Couldnt load reaction DB" );
      
      try
      {
        const ReactionGamma::Reaction *reaction = nullptr;
        vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
        reaction_db->gammas( nucstr, possible_rctns );
        
        // TODO: we are currently taking the first reaction; however, in principle there could be multiple - however, `ReactionGamma` doesnt have an interface to just return a reaction by name, I guess because
        for( size_t i = 0; !reaction && (i < possible_rctns.size()); ++i )
          reaction = possible_rctns[i].reaction;
        
        if( reaction )
          return reaction;
      }catch( std::exception & )
      {
        //ReactionGamma::gammas(...) throws if not a valid reaction
      }//try / catch
      
      return std::monostate{};
    }
    
    WColor color() const
    {
      return m_color_select->color();
    }
    
    void setColor( const WColor &color )
    {
      m_color_select->setColor( color );
    }
    
    double age() const
    {
      const auto nuc_input = nuclide();
      if( !std::holds_alternative<const SandiaDecay::Nuclide *>(nuc_input) )
        return 0.0;
      
      const SandiaDecay::Nuclide * const nuc = std::get<const SandiaDecay::Nuclide *>(nuc_input);
      if( !nuc )
        return 0.0;
    
      if( m_age_edit->isHidden() || m_age_edit->text().empty() )
        return PeakDef::defaultDecayTime( nuc );
      
      double age = 0.0;
      try
      {
        const string agestr = m_age_edit->text().toUTF8();
        age = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( agestr, nuc->halfLife );
      }catch( std::exception & )
      {
        age = PeakDef::defaultDecayTime( nuc );
      }//try / catch
      
      return age;
    }//double age() const
    
    
    RelActCalcAuto::NucInputInfo toNucInputInfo() const
    {
      const auto nuc_input = nuclide();
      
   
      RelActCalcAuto::NucInputInfo nuc_info;

      if( std::holds_alternative<const SandiaDecay::Nuclide *>(nuc_input) )
      {
        nuc_info.nuclide = std::get<const SandiaDecay::Nuclide *>(nuc_input);
        nuc_info.age = age(); // Must not be negative.
        nuc_info.fit_age = m_fit_age->isChecked();
      }else if( std::holds_alternative<const SandiaDecay::Element *>(nuc_input) )
      {
        nuc_info.element = std::get<const SandiaDecay::Element *>(nuc_input);
      }else if( std::holds_alternative<const ReactionGamma::Reaction *>(nuc_input) )
      {
        nuc_info.reaction = std::get<const ReactionGamma::Reaction *>(nuc_input);
      }else
      {
        throw runtime_error( "No valid nuclide" );
      }
      
      // TODO: Not implemented: vector<double> gammas_to_exclude;
      //nuc_info.gammas_to_exclude = ;
      
      nuc_info.peak_color_css = m_color_select->color().cssText();

      return nuc_info;
    }//RelActCalcAuto::RoiRange toRoiRange() const
    
    
    void fromNucInputInfo( const RelActCalcAuto::NucInputInfo &info )
    {
      if( !info.nuclide && !info.element && !info.reaction )
      {
        m_nuclide_edit->setText( "" );
        if( !info.peak_color_css.empty() )
          m_color_select->setColor( WColor(info.peak_color_css) );
        
        m_age_label->hide();
        m_age_edit->hide();
        m_age_edit->setText( "0s" );
        m_fit_age->setUnChecked();
        m_fit_age->hide();
        
        m_updated.emit();
        
        return;
      }//if( !info.nuclide && !info.element && !info.reaction )
      
      if( info.nuclide )
        m_nuclide_edit->setText( WString::fromUTF8(info.nuclide->symbol) );
      else if( info.element )
        m_nuclide_edit->setText( WString::fromUTF8(info.element->symbol) );
      else if( info.reaction )
        m_nuclide_edit->setText( WString::fromUTF8(info.reaction->name()) );
      
      if( !info.peak_color_css.empty() )
        m_color_select->setColor( WColor(info.peak_color_css) );
      
      if( info.nuclide )
      {
        const SandiaDecay::Nuclide * const nuc = info.nuclide;
        
        const bool age_is_fittable = !PeakDef::ageFitNotAllowed(nuc);
        m_age_label->setHidden( !age_is_fittable );
        m_age_edit->setHidden( !age_is_fittable );
        m_fit_age->setHidden( !age_is_fittable );
        m_fit_age->setChecked( age_is_fittable && info.fit_age );
      
        if( !age_is_fittable || (info.age < 0.0) )
        {
          string agestr = "0s";
          if( nuc )
            PeakDef::defaultDecayTime( nuc, &agestr );
          m_age_edit->setText( WString::fromUTF8(agestr) );
        }else
        {
          const string agestr = PhysicalUnitsLocalized::printToBestTimeUnits(info.age);
          m_age_edit->setText( WString::fromUTF8(agestr) );
        }

         // TODO: blah blah blah - implement the below
         //std::optional<double> info.fit_age_min;  
         //std::optional<double> info.fit_age_max;
      }else 
      {
        m_age_label->hide();
        m_age_edit->hide();
        m_fit_age->setUnChecked();
        m_fit_age->hide();
      }//if( info.nuclide )
      
      // Not currently supported: info.gammas_to_exclude -> vector<double>;

      // TODO: blah blah blah - implement the below
      //std::optional<double> info.min_rel_act;
      //std::optional<double> info.max_rel_act;
      //std::optional<double> info.starting_rel_act;
      //std::vector<double> info.gammas_to_exclude;

      handleIsotopeChange();
    }//void fromNucInputInfo( const RelActCalcAuto::NucInputInfo &info )
    
    
    void handleRemoveSelf()
    {
      m_remove.emit();
    }
    
    
    Wt::Signal<> &updated()
    {
      return m_updated;
    }
    
    
    Wt::Signal<> &remove()
    {
      return m_remove;
    }

    void addActRatioConstraint( const RelActCalcAuto::RelEffCurveInput::ActRatioConstraint &constraint )
    {
      // TODO: blah blah blah - implement this
    }

    void addMassFractionConstraint( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &constraint )
    {
      // TODO: blah blah blah - implement this
    }

    void setIsInCurves( const std::set<size_t> &curves_with_nuc, size_t num_rel_eff_curves )
    {
      // TODO: blah blah blah - implement this
    }
  };//class RelActAutoNuclide

  
  
  
  class RelActFreePeak : public WContainerWidget
  {
    RelActAutoGui *m_gui;
    NativeFloatSpinBox *m_energy;
    WCheckBox *m_fwhm_constrained;
    WCheckBox *m_apply_energy_cal;
    WText *m_invalid;
    
    Wt::Signal<> m_updated;
    Wt::Signal<> m_remove;
    
  public:
    RelActFreePeak( RelActAutoGui *gui, WContainerWidget *parent = nullptr )
    : WContainerWidget( parent ),
    m_gui( gui ),
    m_energy( nullptr ),
    m_fwhm_constrained( nullptr ),
    m_apply_energy_cal( nullptr ),
    m_updated( this ),
    m_remove( this )
    {
      addStyleClass( "RelActFreePeak" );
      
      const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );
      
      WLabel *label = new WLabel( "Energy", this );
      label->addStyleClass( "GridFirstCol GridFirstRow" );
      
      m_energy = new NativeFloatSpinBox( this );
      label->setBuddy( m_energy );
      m_energy->valueChanged().connect( this, &RelActFreePeak::handleEnergyChange );
      m_energy->addStyleClass( "GridSecondCol GridFirstRow" );
      
      // Things are a little cramped if we include the keV
      //label = new WLabel( "keV", this );
      //label->addStyleClass( "GridThirdCol GridFirstRow" );
      
      WPushButton *removeFreePeak = new WPushButton( this );
      removeFreePeak->setStyleClass( "DeleteEnergyRangeOrNuc Wt-icon" );
      removeFreePeak->setIcon( "InterSpec_resources/images/minus_min_black.svg" );
      removeFreePeak->clicked().connect( this, &RelActFreePeak::handleRemoveSelf );
      removeFreePeak->addStyleClass( "GridThirdCol GridFirstRow" );
      
      m_fwhm_constrained = new WCheckBox( "Constrain FWHM", this );
      m_fwhm_constrained->setChecked( true );
      m_fwhm_constrained->checked().connect( this, &RelActFreePeak::handleFwhmConstrainChanged );
      m_fwhm_constrained->unChecked().connect( this, &RelActFreePeak::handleFwhmConstrainChanged );
      m_fwhm_constrained->addStyleClass( "FreePeakConstrain GridFirstCol GridSecondRow GridSpanThreeCol" );
      
      const char *tooltip = "When checked, this peak will be constrained to the FWHM functional form gamma and x-ray"
      "peaks are normally constrained to.<br />"
      "When un-checked, the FWHM for this peak will be fit from the data, independent of all other peak widths.<br />"
      "Un-checking this option is useful for annihilation and reaction photopeaks.";
      HelpSystem::attachToolTipOn( m_fwhm_constrained, tooltip, showToolTips );
      
      
      m_apply_energy_cal = new WCheckBox( "True Energy", this );
      m_apply_energy_cal->setStyleClass( "CbNoLineBreak" );
      m_apply_energy_cal->setChecked( true );
      m_apply_energy_cal->checked().connect( this, &RelActFreePeak::handleApplyEnergyCalChanged );
      m_apply_energy_cal->unChecked().connect( this, &RelActFreePeak::handleApplyEnergyCalChanged );
      m_apply_energy_cal->addStyleClass( "FreePeakConstrain GridFirstCol GridThirdRow GridSpanThreeCol" );
      tooltip = "Check this option if the peak is for a gamma of known energy.<br />"
      "Un-check this option if this peak is an observed peak in the spectrum with unknown true energy.";
      HelpSystem::attachToolTipOn( m_apply_energy_cal, tooltip, showToolTips );
      
      
      m_invalid = new WText( "Not in a ROI.", this );
      m_invalid->addStyleClass( "InvalidFreePeakEnergy GridFirstCol GridFourthRow GridSpanThreeCol" );
      m_invalid->hide();
    }//RelActFreePeak constructor
    
    
    float energy() const
    {
      return m_energy->value();
    }
    
    
    void setEnergy( const float energy )
    {
      m_energy->setValue( energy );
      m_updated.emit();
    }
    
    
    bool fwhmConstrained() const
    {
      return m_fwhm_constrained->isChecked();
    }
    
    
    void setFwhmConstrained( const bool constrained )
    {
      m_fwhm_constrained->setChecked( constrained );
      m_updated.emit();
    }
    
    
    bool applyEnergyCal() const
    {
      return (m_apply_energy_cal->isVisible() && m_apply_energy_cal->isChecked());
    }
    
    
    void setApplyEnergyCal( const bool apply )
    {
      m_apply_energy_cal->setChecked( apply );
    }
    
    
    void setInvalidEnergy( const bool invalid )
    {
      if( invalid != m_invalid->isVisible() )
        m_invalid->setHidden( !invalid );
    }
    
    
    void handleRemoveSelf()
    {
      m_remove.emit();
    }
    
    
    void handleEnergyChange()
    {
      m_updated.emit();
    }//
    
    
    void handleFwhmConstrainChanged()
    {
      m_updated.emit();
    }
    
    
    void handleApplyEnergyCalChanged()
    {
      m_updated.emit();
    }
    
    
    void setApplyEnergyCalVisible( const bool visible )
    {
      if( m_apply_energy_cal->isVisible() != visible )
        m_apply_energy_cal->setHidden( !visible );
    }
    
    
    Wt::Signal<> &updated()
    {
      return m_updated;
    }
    
    
    Wt::Signal<> &remove()
    {
      return m_remove;
    }
  };//RelActFreePeak
}//namespace


std::pair<RelActAutoGui *,AuxWindow *> RelActAutoGui::createWindow( InterSpec *viewer  )
{
  assert( viewer );
  
  AuxWindow *window = nullptr;
  RelActAutoGui *disp = nullptr;
  
  try
  {
    disp = new RelActAutoGui( viewer );
    
    window = new AuxWindow( "Relative Act. Isotopics", 
                           Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::SetCloseable)
                           | AuxWindowProperties::EnableResize );
    // We have to set minimum size before calling setResizable, or else Wt's Resizable.js functions
    //  will be called first, which will then default to using the initial size as minimum allowable
    window->setMinimumSize( 800, 480 );
    window->setResizable( true );
    window->contents()->setOffsets(WLength(0,WLength::Pixel));
    
    disp->setHeight( WLength(100, WLength::Percentage) );
    disp->setWidth( WLength(100, WLength::Percentage) );
    
    window->contents()->addWidget( disp );
    
    //window->stretcher()->addWidget( disp, 0, 0 );
    //window->stretcher()->setContentsMargins(0,0,0,0);
    //    window->footer()->resize(WLength::Auto, WLength(50.0));
    
    AuxWindow::addHelpInFooter( window->footer(), "rel-act-dialog" );
    
    disp->addDownloadAndUploadLinks( window->footer() );
    
    WPushButton *closeButton = window->addCloseButtonToFooter();
    closeButton->clicked().connect(window, &AuxWindow::hide);
    
    //window->rejectWhenEscapePressed();
    
    // TODO: Similar to activity shielding fit, should store the current widget state in the SpecMeas
    
    double windowWidth = viewer->renderedWidth();
    double windowHeight = viewer->renderedHeight();
    
    if( (windowHeight > 110) && (windowWidth > 110) )
    {
      if( !viewer->isPhone() )
      {
        // A size of 1050px by 555px is about the smallest that renders everything nicely.
        if( (windowWidth > (1050.0/0.8)) && (windowHeight > (555.0/0.8)) )
        {
          windowWidth = 0.8*windowWidth;
          windowHeight = 0.8*windowHeight;
          window->resizeWindow( windowWidth, windowHeight );
        }else
        {
          windowWidth = 0.9*windowWidth;
          windowHeight = 0.9*windowHeight;
          window->resizeWindow( windowWidth, windowHeight );
        }
      }//if( !viewer->isPhone() )
    }else if( !viewer->isPhone() )
    {
      //When loading an application state that is showing this window, we may
      //  not know the window size (e.g., windowWidth==windowHeight==0), so
      //  instead skip giving the initial size hint, and instead size things
      //  client side (maybe we should just do this always?)
      window->resizeScaledWindow( 0.90, 0.90 );
    }
    
    window->centerWindow();
    window->finished().connect( viewer, &InterSpec::closeShieldingSourceFit );
    
    window->WDialog::setHidden(false);
    window->show();
    window->centerWindow();
  }catch( std::exception &e )
  {
    passMessage( "Error creating Relative Act. Isotopics tool: " + string(e.what()),
                WarningWidget::WarningMsgHigh );
    
    if( disp )
      delete disp;
    disp = nullptr;
    
    if( window )
      AuxWindow::deleteAuxWindow( window );
    window = nullptr;
  }//try / catch
  
  return make_pair( disp, window );
}//createWindow( InterSpec *viewer  )



const char *RelActAutoGui::to_str( const RelActAutoGui::AddUncert val )
{
  switch( val )
  {
    case RelActAutoGui::AddUncert::StatOnly:           return "StatOnly";
    case RelActAutoGui::AddUncert::OneHundrethPercent: return "OneHundrethPercent";
    case RelActAutoGui::AddUncert::OneTenthPercent:    return "OneTenthPercent";
    case RelActAutoGui::AddUncert::OnePercent:         return "OnePercent";
    case RelActAutoGui::AddUncert::FivePercent:        return "FivePercent";
    case RelActAutoGui::AddUncert::TenPercent:         return "TenPercent";
    case RelActAutoGui::AddUncert::TwentyFivePercent:  return "TwentyFivePercent";
    case RelActAutoGui::AddUncert::FiftyPercent:       return "FiftyPercent";
    case RelActAutoGui::AddUncert::SeventyFivePercent: return "SeventyFivePercent";
    case RelActAutoGui::AddUncert::OneHundredPercent:  return "OneHundredPercent";
    case RelActAutoGui::AddUncert::NumAddUncert:       return "NumAddUncert";
  }//
  
  return "InvalidAddUncert";
}//to_str( const AddUncert val )



RelActAutoGui::RelActAutoGui( InterSpec *viewer, Wt::WContainerWidget *parent )
: WContainerWidget( parent ),
  m_render_flags( 0 ),
  m_default_par_sets_dir( "" ),
  m_user_par_sets_dir( "" ),
  m_interspec( viewer ),
  m_foreground( nullptr ),
  m_background( nullptr ),
  m_background_sf( 0.0 ),
  m_status_indicator( nullptr ),
  m_spectrum( nullptr ),
  m_peak_model( nullptr ),
  m_rel_eff_chart( nullptr ),
  m_txt_results( nullptr ),
  m_upper_menu( nullptr ),
  m_presets( nullptr ),
  m_loading_preset( false ),
  m_current_preset_index( -1 ),
  m_error_msg( nullptr ),
  m_fit_chi2_msg( nullptr ),
  m_rel_eff_opts_menu( nullptr ),
  m_rel_eff_opts_stack( nullptr ),
  m_fwhm_eqn_form( nullptr ),
  m_fwhm_estimation_method( nullptr ),
  m_fit_energy_cal( nullptr ),
  m_background_subtract( nullptr ),
  m_same_z_age( nullptr ),
  m_skew_type( nullptr ),
  m_add_uncert( nullptr ),
  m_more_options_menu( nullptr ),
  m_apply_energy_cal_item( nullptr ),
  m_show_ref_lines_item( nullptr ),
  m_hide_ref_lines_item( nullptr ),
  m_set_peaks_foreground( nullptr ),
  m_photopeak_widget(),
  m_clear_energy_ranges( nullptr ),
  m_show_free_peak( nullptr ),
  m_free_peaks_container( nullptr ),
//  m_u_pu_data_source( nullptr ),
  m_nuclides( nullptr ),
  m_energy_ranges( nullptr ),
  m_free_peaks( nullptr ),
  m_is_calculating( false ),
  m_cancel_calc{},
  m_solution{},
  m_calc_started( this ),
  m_calc_successful( this ),
  m_calc_failed( this ),
  m_solution_updated( this ),
  m_html_download_rsc( new RelActAutoReportResource( this, this ) ),
  m_xml_download_rsc( new RelActAutoParamsResource( this, this ) )
{
  assert( m_interspec );
  if( !m_interspec )
    throw runtime_error( "RelActAutoGui: requires pointer to InterSpec" );
  
  new UndoRedoManager::BlockGuiUndoRedo( this );
    
  wApp->useStyleSheet( "InterSpec_resources/RelActAutoGui.css" );
  wApp->useStyleSheet( "InterSpec_resources/GridLayoutHelpers.css" );
  
  addStyleClass( "RelActAutoGui" );
  
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", m_interspec );
  
  //WText *alpha_warning = new WText( "This tool is under active development - this is an early preview", this );
  //alpha_warning->addStyleClass( "RelActCalcAutoAlphaBuildWarning" );
    
  WContainerWidget *upper_div = new WContainerWidget( this );
  upper_div->addStyleClass( "RelActAutoUpperArea" );
  
  WStackedWidget *upper_stack = new WStackedWidget();
  upper_stack->addStyleClass( "UpperStack" );
  WAnimation animation(Wt::WAnimation::Fade, Wt::WAnimation::Linear, 200);
  upper_stack->setTransitionAnimation( animation, true );
  
  m_upper_menu = new WMenu( upper_stack, Wt::Vertical, upper_div );
  m_upper_menu->addStyleClass( "UpperMenu LightNavMenu" );
  upper_div->addWidget( upper_stack );

  
  m_spectrum = new D3SpectrumDisplayDiv();
  m_spectrum->setCompactAxis( true );
  m_spectrum->disableLegend();
  m_interspec->colorThemeChanged().connect( boost::bind( &D3SpectrumDisplayDiv::applyColorTheme, m_spectrum, boost::placeholders::_1 ) );
  m_spectrum->applyColorTheme( m_interspec->getColorTheme() );
  
  const bool logypref = UserPreferences::preferenceValue<bool>( "LogY", m_interspec );
  m_spectrum->setYAxisLog( logypref );
  
  auto set_log_y = wApp->bind( boost::bind( &D3SpectrumDisplayDiv::setYAxisLog, m_spectrum, true ) );
  auto set_lin_y = wApp->bind( boost::bind( &D3SpectrumDisplayDiv::setYAxisLog, m_spectrum, false ) );
  std::function<void (boost::any)> logy_fcn = [=](boost::any value){
    if( boost::any_cast<bool>(value) )
      set_log_y();
    else
      set_lin_y();
  };
  
  m_interspec->preferences()->addCallbackWhenChanged( "LogY", m_spectrum, 
                                                     &D3SpectrumDisplayDiv::setYAxisLog );
  
  m_peak_model = new PeakModel( m_spectrum );
  m_peak_model->setNoSpecMeasBacking();
  
  m_spectrum->setPeakModel( m_peak_model );
  
  m_spectrum->existingRoiEdgeDragUpdate().connect( this, &RelActAutoGui::handleRoiDrag );
  m_spectrum->dragCreateRoiUpdate().connect( this, &RelActAutoGui::handleCreateRoiDrag );
  m_spectrum->rightClicked().connect( this, &RelActAutoGui::handleRightClick );
  m_spectrum->shiftKeyDragged().connect( this, &RelActAutoGui::handleShiftDrag );
  m_spectrum->doubleLeftClick().connect( this, &RelActAutoGui::handleDoubleLeftClick );
  
  m_rel_eff_chart = new RelEffChart();
  m_txt_results = new RelActTxtResults();
  
  WMenuItem *item = new WMenuItem( "Spec.", m_spectrum );
  m_upper_menu->addItem( item );
  
  // When outside the link area is clicked, the item doesnt get selected, so we'll work around this.
  item->clicked().connect( std::bind([this,item](){
    m_upper_menu->select( item );
    item->triggered().emit( item );
  }) );
  
  item = new WMenuItem( "Rel. Eff.", m_rel_eff_chart );
  m_upper_menu->addItem( item );
  
  item->clicked().connect( std::bind([this,item](){
    m_upper_menu->select( item );
    item->triggered().emit( item );
  }) );
  
  
  item = new WMenuItem( "Result", m_txt_results );
  m_upper_menu->addItem( item );
  
  item->clicked().connect( std::bind([this,item](){
    m_upper_menu->select( item );
    item->triggered().emit( item );
  }) );
  
  m_upper_menu->select( static_cast<int>(0) );
  
  m_interspec->spectrumScaleFactorChanged().connect( boost::bind(
                 &RelActAutoGui::handleDisplayedSpectrumChange, this, boost::placeholders::_1) );
  m_interspec->displayedSpectrumChanged().connect( boost::bind(
                 &RelActAutoGui::handleDisplayedSpectrumChange, this, boost::placeholders::_1) );
  
  m_interspec->detectorChanged().connect( this, &RelActAutoGui::handleDetectorChange );
  m_interspec->detectorModified().connect( this, &RelActAutoGui::handleDetectorChange );

  m_default_par_sets_dir = SpecUtils::append_path( InterSpec::staticDataDirectory(), "rel_act" );
  const vector<string> default_par_sets = SpecUtils::recursive_ls( m_default_par_sets_dir, ".xml" );
  
  vector<string> user_par_sets;
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER )
  try
  {
    m_user_par_sets_dir = SpecUtils::append_path( InterSpec::writableDataDirectory(), "rel_act" );
    user_par_sets = SpecUtils::recursive_ls( m_user_par_sets_dir, ".xml" );
  }catch( std::exception & )
  {
    cerr << "RelActAutoGui: Writable data directory not set.\n";
  }
#endif
  
  WContainerWidget *presetDiv = new WContainerWidget( this );
  presetDiv->addStyleClass( "PresetsRow" );
  WLabel *label = new WLabel( "Presets", presetDiv );
  m_presets = new WComboBox( presetDiv );
  label->setBuddy( m_presets );
  
  m_presets->addItem( "Blank" );
  m_preset_paths.push_back( "" );
  
  for( const string &filename : default_par_sets )
  {
    string dispname = SpecUtils::filename(filename);
    if( dispname.size() > 4 )
      dispname = dispname.substr(0, dispname.size() - 4);
    
    m_presets->addItem( WString::fromUTF8(dispname) );
    m_preset_paths.push_back( filename );
  }//for( string filename : default_par_sets )
  
  for( const string &filename : user_par_sets )
  {
    string dispname = SpecUtils::filename(filename);
    if( dispname.size() > 4 )
      dispname = dispname.substr(0, filename.size() - 4);
    
    m_presets->addItem( WString::fromUTF8("User: " + dispname) );
    m_preset_paths.push_back( filename );
  }//for( string filename : default_par_sets )
  
  // m_presets->addItem( WString::fromUTF8("Custom") );
  
  m_presets->setCurrentIndex( 0 );
  m_current_preset_index = 0;
  m_presets->changed().connect( this, &RelActAutoGui::handlePresetChange );
  
  WContainerWidget *spacer = new WContainerWidget( presetDiv );
  spacer->addStyleClass( "RelActAutoSpacer" );

  m_error_msg = new WText( "Not Calculated.", presetDiv );
  m_error_msg->addStyleClass( "RelActAutoErrMsg" );

  m_fit_chi2_msg = new WText( "", presetDiv );
  m_fit_chi2_msg->addStyleClass( "RelActAutoChi2Msg" );

  m_status_indicator = new WText( "Calculating...", presetDiv );
  m_status_indicator->addStyleClass( "RelActAutoStatusMsg" );
  m_status_indicator->hide();
  
  // We'll take care of the options that apply to all types of Rel Eff curves now.
  WGroupBox *generalOptionsDiv = new WGroupBox( "Peak Fitting Options", this );
  generalOptionsDiv->addStyleClass( "RelActAutoGeneralOptionsRow" );

  m_fit_energy_cal = new WCheckBox( "Fit Energy Cal.", generalOptionsDiv );
  m_fit_energy_cal->checked().connect( this, &RelActAutoGui::handleFitEnergyCalChanged );
  m_fit_energy_cal->unChecked().connect( this, &RelActAutoGui::handleFitEnergyCalChanged );
  
  m_background_subtract = new WCheckBox( "Back. Sub.", generalOptionsDiv );
  m_background_subtract->checked().connect( this, &RelActAutoGui::handleBackgroundSubtractChanged );
  m_background_subtract->unChecked().connect( this, &RelActAutoGui::handleBackgroundSubtractChanged );
  
  WContainerWidget *fwhmEstDiv = new WContainerWidget( generalOptionsDiv );
  fwhmEstDiv->addStyleClass( "RelActAutoFwhmEstDiv" );
  label = new WLabel( "FWHM Est.", fwhmEstDiv );
  m_fwhm_estimation_method = new WComboBox( fwhmEstDiv );
  label->setBuddy( m_fwhm_estimation_method );

  for( int i = 0; i <= static_cast<int>(RelActCalcAuto::FwhmEstimationMethod::FixedToDetectorEfficiency); ++i )
  {
    const char *name = "";
    switch( RelActCalcAuto::FwhmEstimationMethod(i) )
    {
      case RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum: 
        name = "Start from Det/Peaks"; 
        break;

      case RelActCalcAuto::FwhmEstimationMethod::StartingFromAllPeaksInSpectrum: 
        name = "Start from Peaks"; 
        break;
      
      case RelActCalcAuto::FwhmEstimationMethod::FixedToAllPeaksInSpectrum: 
        name = "Fixed to Peaks"; 
        break;

      case RelActCalcAuto::FwhmEstimationMethod::StartingFromDetectorEfficiency: 
        name = "Start from Det. Eff."; 
        break;

      case RelActCalcAuto::FwhmEstimationMethod::FixedToDetectorEfficiency: 
        name = "Fixed to Det. Eff."; 
        break;
    }//switch( RelActCalcAuto::FwhmEstimationMethod(i) )
    
    m_fwhm_estimation_method->addItem( name );
  }//for( loop over RelActCalcAuto::FwhmEstimationMethod )

  WContainerWidget *fwhmFormDiv = new WContainerWidget( generalOptionsDiv );
  fwhmFormDiv->addStyleClass( "RelActAutoFwhmFormDiv" );
  label = new WLabel( "FWHM Form", fwhmFormDiv );
  
  m_fwhm_eqn_form = new WComboBox( fwhmFormDiv );
  label->setBuddy( m_fwhm_eqn_form );

  for( int i = 0; i <= static_cast<int>(RelActCalcAuto::FwhmForm::Polynomial_6); ++i )
  {
    const char *name = "";
    switch( RelActCalcAuto::FwhmForm(i) )
    {
      case RelActCalcAuto::FwhmForm::Gadras:        name = "Gadras"; break;
      case RelActCalcAuto::FwhmForm::SqrtEnergyPlusInverse:  name = "sqrt(A0 + A1*E + A2/E)"; break;
      case RelActCalcAuto::FwhmForm::ConstantPlusSqrtEnergy: name = "A0 + A1*sqrt(E)"; break;
      case RelActCalcAuto::FwhmForm::Polynomial_2:  name = "sqrt(A0 + A1*E)"; break;
      case RelActCalcAuto::FwhmForm::Polynomial_3:  name = "sqrt(A0 + A1*E + A2*E*E)"; break;
      case RelActCalcAuto::FwhmForm::Polynomial_4:  name = "sqrt(A0 + A1*E^1...A3*E^3)"; break;
      case RelActCalcAuto::FwhmForm::Polynomial_5:  name = "sqrt(A0 + A1*E^1...A4*E^4)"; break;
      case RelActCalcAuto::FwhmForm::Polynomial_6:  name = "sqrt(A0 + A1*E^1...A5*E^5)"; break;
      case RelActCalcAuto::FwhmForm::NotApplicable: name = "Use Det. Eff."; break;
    }//switch( RelActCalcAuto::FwhmForm(i) )
    
    m_fwhm_eqn_form->addItem( name );
  }//for( loop over RelActCalcAuto::FwhmForm )
  
  const char *tooltip = "The equation type used to model peak FWHM as a function of energy.";
  HelpSystem::attachToolTipOn( fwhmFormDiv, tooltip, showToolTips );
  
  // TODO: need to set m_fwhm_eqn_form based on energy ranges selected
  m_fwhm_eqn_form->setCurrentIndex( static_cast<int>(RelActCalcAuto::FwhmForm::SqrtEnergyPlusInverse) );
  m_fwhm_estimation_method->setCurrentIndex( static_cast<int>(RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum) );
  m_fwhm_eqn_form->changed().connect( this, &RelActAutoGui::handleFwhmFormChanged );
  m_fwhm_estimation_method->changed().connect( this, &RelActAutoGui::handleFwhmEstimationMethodChanged );
  
  WContainerWidget *skewDiv = new WContainerWidget( generalOptionsDiv );
  skewDiv->addStyleClass( "RelActAutoSkewDiv" );
  label = new WLabel( "Peak Skew", skewDiv );
  m_skew_type = new WComboBox( skewDiv );
  label->setBuddy( m_skew_type );
  m_skew_type->activated().connect( this, &RelActAutoGui::handleSkewTypeChanged );
  tooltip = "The type of skew to apply to the peaks; skew parameters will be fit.";
  HelpSystem::attachToolTipOn( {label,m_skew_type}, tooltip, showToolTips );
  for( auto st = PeakDef::SkewType(0); st <= PeakDef::SkewType::DoubleSidedCrystalBall;
      st = PeakDef::SkewType(st + 1) )
  {
    const char *label = PeakDef::to_label(st);
    m_skew_type->addItem( label );
  }//for( loop over SkewTypes )
    
  m_skew_type->setCurrentIndex( 0 );
  
  WContainerWidget *addUncertDiv = new WContainerWidget( generalOptionsDiv );
  addUncertDiv->addStyleClass( "RelActAutoAddUncertDiv" );
  label = new WLabel( "Add. Uncert", addUncertDiv );
  m_add_uncert = new WComboBox( addUncertDiv );
  label->setBuddy( m_add_uncert );
  m_add_uncert->activated().connect( this, &RelActAutoGui::handleAdditionalUncertChanged );
    
  for( RelActAutoGui::AddUncert i = RelActAutoGui::AddUncert(0);
      i < RelActAutoGui::AddUncert::NumAddUncert;
      i = RelActAutoGui::AddUncert(static_cast<int>(i) + 1) )
  {
    WString uncert_txt;
    switch( i )
    {
      case AddUncert::StatOnly:           uncert_txt = WString::fromUTF8("None");  break;
      case AddUncert::OneHundrethPercent: uncert_txt = WString::fromUTF8("0.01%"); break;
      case AddUncert::OneTenthPercent:    uncert_txt = WString::fromUTF8("0.1%");  break;
      case AddUncert::OnePercent:         uncert_txt = WString::fromUTF8("1%");    break;
      case AddUncert::FivePercent:        uncert_txt = WString::fromUTF8("5%");    break;
      case AddUncert::TenPercent:         uncert_txt = WString::fromUTF8("10%");   break;
      case AddUncert::TwentyFivePercent:  uncert_txt = WString::fromUTF8("25%");   break;
      case AddUncert::FiftyPercent:       uncert_txt = WString::fromUTF8("50%");   break;
      case AddUncert::SeventyFivePercent: uncert_txt = WString::fromUTF8("75%");   break;
      case AddUncert::OneHundredPercent:  uncert_txt = WString::fromUTF8("100%");  break;
      case AddUncert::NumAddUncert:       assert(0);                               break;
    }//switch( i )
       
    m_add_uncert->addItem( uncert_txt );
  }//for( loop over AddUncert )
     
  m_add_uncert->setCurrentIndex( static_cast<int>(RelActAutoGui::AddUncert::StatOnly) );
    
  
  WGroupBox *optionsDiv = new WGroupBox( "Relative Efficiency Curve Options", this );
  optionsDiv->addStyleClass( "RelActAutoOptions" );

  m_rel_eff_opts_stack = new WStackedWidget();
  m_rel_eff_opts_stack->addStyleClass( "RelEffCurveOptsStack" );
  //Do not set a transition animation, it will cause all elements of the stack to be hidden,
  //  and totally stuck hidden, when we remove an element from the WMenu/WStackedWidget.
  //m_rel_eff_opts_stack->setTransitionAnimation( animation, true );

  m_rel_eff_opts_menu = new WMenu( m_rel_eff_opts_stack, optionsDiv );
  m_rel_eff_opts_menu->addStyleClass( "RelEffCurveOptsMenu LightNavMenu" );
  optionsDiv->addWidget( m_rel_eff_opts_stack );

  handleAddRelEffCurve();

/*
  label = new WLabel( "Yield Info", optionsDiv );
  label->addStyleClass( "GridSeventhCol GridFirstRow" );
  
  m_u_pu_data_source = new WComboBox( optionsDiv );
  label->setBuddy( m_u_pu_data_source );
  m_u_pu_data_source->activated().connect( this, &RelActAutoGui::handleNucDataSrcChanged );
  m_u_pu_data_source->addStyleClass( "GridEighthCol GridFirstRow" );
  
  tooltip = "The nuclear data source for gamma branching ratios of uranium and plutonium.";
  HelpSystem::attachToolTipOn( {label, m_u_pu_data_source}, tooltip, showToolTips );
  
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
    
    m_u_pu_data_source->addItem( WString::fromUTF8(src_label) );
  }//for( loop over sources )
  
  m_u_pu_data_source->setCurrentIndex( static_cast<int>(NucDataSrc::SandiaDecay) );
  m_u_pu_data_source->changed().connect( this, &RelActAutoGui::handleDataSrcChanged );
 */
  
    
  WPushButton *more_btn = new WPushButton( this );
  more_btn->setIcon( "InterSpec_resources/images/more_menu_icon.svg" );
  more_btn->addStyleClass( "MoreMenuIcon Wt-icon" );
  
  m_more_options_menu = new PopupDivMenu( more_btn, PopupDivMenu::TransientMenu );
  const bool is_phone = false; //isPhone();
  if( is_phone )
    m_more_options_menu->addPhoneBackItem( nullptr );
  
  m_apply_energy_cal_item = m_more_options_menu->addMenuItem( "Apply Energy Cal" );
  m_apply_energy_cal_item->triggered().connect( this, &RelActAutoGui::startApplyFitEnergyCalToSpecFile );
  
  m_show_ref_lines_item = m_more_options_menu->addMenuItem( "Show Ref. Gamma Lines" );
  m_show_ref_lines_item->triggered().connect( boost::bind( &RelActAutoGui::handleShowRefLines, this, true ) );
  
  m_hide_ref_lines_item = m_more_options_menu->addMenuItem( "Hide Ref. Gamma Lines" );
  m_hide_ref_lines_item->triggered().connect( boost::bind( &RelActAutoGui::handleShowRefLines, this, false ) );
  m_hide_ref_lines_item->setHidden( true );
  m_hide_ref_lines_item->setDisabled( true );
  
  m_set_peaks_foreground = m_more_options_menu->addMenuItem( "Set Peaks to foreground" );
  m_set_peaks_foreground->triggered().connect( boost::bind( &RelActAutoGui::setPeaksToForeground, this ) );
  m_set_peaks_foreground->setDisabled( true );
  
  WContainerWidget *bottomArea = new WContainerWidget( this );
  bottomArea->addStyleClass( "EnergiesAndNuclidesHolder" );
  
  //WContainerWidget *nuclidesHolder = new WContainerWidget( bottomArea );
  WGroupBox *nuclidesHolder = new WGroupBox( "Nuclides", bottomArea );
  nuclidesHolder->addStyleClass( "NuclidesHolder" );
  
  //WContainerWidget *energiesHolder = new WContainerWidget( bottomArea );
  WGroupBox *energiesHolder = new WGroupBox( "Energy Ranges", bottomArea );
  energiesHolder->addStyleClass( "EnergiesHolder" );
  
  //m_free_peaks_container = new WContainerWidget( bottomArea );
  m_free_peaks_container = new WGroupBox( "Free Peaks", bottomArea );
  m_free_peaks_container->addStyleClass( "FreePeaksHolder" );
  m_free_peaks_container->hide();
  
  //WText *nuc_header = new WText( "Nuclides", nuclidesHolder );
  //nuc_header->addStyleClass( "EnergyNucHeader" );
  
  m_nuclides = new WContainerWidget( nuclidesHolder );
  m_nuclides->addStyleClass( "EnergyNucContent" );
  
  WContainerWidget *nuc_footer = new WContainerWidget( nuclidesHolder );
  nuc_footer->addStyleClass( "EnergyNucFooter" );
  
  WPushButton *add_nuc_icon = new WPushButton( nuc_footer );
  add_nuc_icon->setStyleClass( "AddEnergyRangeOrNuc Wt-icon" );
  add_nuc_icon->setIcon("InterSpec_resources/images/plus_min_black.svg");
  add_nuc_icon->clicked().connect( this, &RelActAutoGui::handleAddNuclide );
  tooltip = "Add a nuclide.";
  HelpSystem::attachToolTipOn( add_nuc_icon, tooltip, showToolTips );

  spacer = new WContainerWidget( nuc_footer );
  spacer->addStyleClass( "RelActAutoSpacer" );

  // same_z_age is something that _could_ be a per-relative-efficiency option, 
  //  but it's a little cleaner and maybe clearer to have it here, near the nuclides
  //  (and we are currently forcing nuclides to be same age between RelEff curves...)
  m_same_z_age = new WCheckBox( "Same El. Same Age", nuc_footer );
  m_same_z_age->addStyleClass( "SameZAgeCb CbNoLineBreak" );
  m_same_z_age->checked().connect( this, &RelActAutoGui::handleSameAgeChanged );
  m_same_z_age->unChecked().connect( this, &RelActAutoGui::handleSameAgeChanged );
  
  //WText *energy_header = new WText( "Energy Ranges", energiesHolder );
  //energy_header->addStyleClass( "EnergyNucHeader" );
  
  m_energy_ranges = new WContainerWidget( energiesHolder );
  m_energy_ranges->addStyleClass( "EnergyNucContent" );
  
  WContainerWidget *energies_footer = new WContainerWidget( energiesHolder );
  energies_footer->addStyleClass( "EnergyNucFooter" );
  
  WPushButton *add_energy_icon = new WPushButton( energies_footer );
  add_energy_icon->setStyleClass( "AddEnergyRangeOrNuc Wt-icon" );
  add_energy_icon->setIcon("InterSpec_resources/images/plus_min_black.svg");
  add_energy_icon->clicked().connect( this, &RelActAutoGui::handleAddEnergy );
  
  tooltip = "Add an energy range.";
  HelpSystem::attachToolTipOn( add_energy_icon, tooltip, showToolTips );
  
  spacer = new WContainerWidget( energies_footer );
  spacer->addStyleClass( "RelActAutoSpacer" );

  m_show_free_peak = new WPushButton( "add free peaks", energies_footer );
  m_show_free_peak->addStyleClass( "ShowFreePeaks LightButton" );
  m_show_free_peak->clicked().connect( this, &RelActAutoGui::handleShowFreePeaks );
  
  m_clear_energy_ranges = new WPushButton( "clear all ranges", energies_footer );
  m_clear_energy_ranges->addStyleClass( "ClearEnergyRanges LightButton" );
  m_clear_energy_ranges->clicked().connect( this, &RelActAutoGui::handleClearAllEnergyRanges );
  tooltip = "Removes all energy ranges.";
  HelpSystem::attachToolTipOn( m_clear_energy_ranges, tooltip, showToolTips );
  m_clear_energy_ranges->hide();
  
  //WText *free_peaks_header = new WText( "Free Peaks", m_free_peaks_container );
  //free_peaks_header->addStyleClass( "EnergyNucHeader" );
  m_free_peaks = new WContainerWidget( m_free_peaks_container );
  m_free_peaks->addStyleClass( "EnergyNucContent" );
  
  WContainerWidget *free_peaks_footer = new WContainerWidget( m_free_peaks_container );
  free_peaks_footer->addStyleClass( "EnergyNucFooter" );
  
  WPushButton *add_free_peak_icon = new WPushButton( free_peaks_footer );
  add_free_peak_icon->setStyleClass( "AddEnergyRangeOrNuc Wt-icon" );
  add_free_peak_icon->setIcon("InterSpec_resources/images/plus_min_black.svg");
  add_free_peak_icon->clicked().connect( boost::bind( &RelActAutoGui::handleAddFreePeak, this, 0.0, true, true ) );
  
  WPushButton *hide_free_peak = new WPushButton( "Close", free_peaks_footer );
  hide_free_peak->addStyleClass( "HideFreePeaks LightButton" );
  hide_free_peak->clicked().connect( this, &RelActAutoGui::handleHideFreePeaks );
  
  tooltip = "Remove all free peaks, and hide the free peaks input.";
  HelpSystem::attachToolTipOn( hide_free_peak, tooltip, showToolTips );
  
  
  tooltip = "&quot;Free peaks&quot; are peaks with a specific energy, that are not associated with any nuclide,"
  " and their amplitude is fit to the best value for the data.  The peaks may optionally be released from the"
  " functional FWHM constraint as well.<br />"
  "Free peaks are useful to handle peaks from reactions, or from unidentified nuclides.";
  HelpSystem::attachToolTipOn( {m_free_peaks_container,m_show_free_peak}, tooltip, showToolTips );
  
    
  auto html_rsc = dynamic_cast<RelActAutoReportResource *>( m_html_download_rsc );
  m_solution_updated.connect( boost::bind( &RelActAutoReportResource::updateSolution,
                                          html_rsc, boost::placeholders::_1 ) );
  
  m_render_flags |= RenderActions::UpdateSpectra;
  m_render_flags |= RenderActions::UpdateCalculations;
}//RelActAutoGui constructor
  

RelActAutoGui::~RelActAutoGui()
{
  // We need to manually manage any WPopupMenu's we create.
  if( m_more_options_menu && wApp && wApp->domRoot() )
    delete m_more_options_menu;
  m_more_options_menu = nullptr;
  
  if( m_cancel_calc )
    m_cancel_calc->store( true );
}//~RelActAutoGui();


void RelActAutoGui::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  if( m_render_flags.testFlag(RenderActions::UpdateSpectra) )
  {
    updateDuringRenderForSpectrumChange();
    m_render_flags |= RenderActions::UpdateCalculations;
  }
  
  if( m_render_flags.testFlag(RenderActions::UpdateFreePeaks)
     || m_render_flags.testFlag(RenderActions::UpdateEnergyRanges)
     || m_render_flags.testFlag(RenderActions::UpdateFitEnergyCal) )
  {
    updateDuringRenderForFreePeakChange();
    m_render_flags |= RenderActions::UpdateCalculations;
  }
  
  if( m_render_flags.testFlag(RenderActions::UpdateNuclidesPresent) )
  {
    updateDuringRenderForNuclideChange();
    m_render_flags |= RenderActions::UpdateCalculations;
  }
  
  if( m_render_flags.testFlag(RenderActions::UpdateNuclidesPresent)
     || m_render_flags.testFlag(RenderActions::UpdateRefGammaLines) )
  {
    updateDuringRenderForRefGammaLineChange();
  }
  
  if( m_render_flags.testFlag(RenderActions::UpdateEnergyRanges) )
  {
    updateDuringRenderForEnergyRangeChange();
    m_render_flags |= RenderActions::UpdateCalculations;
  }
  
  if( m_render_flags.testFlag(RenderActions::ChartToDefaultRange) )
  {
    updateSpectrumToDefaultEnergyRange();
  }
  
  if( m_render_flags.testFlag(RenderActions::UpdateCalculations) )
  {
    startUpdatingCalculation();
  }
  
  m_render_flags = 0;
  m_loading_preset = false;
  
  WContainerWidget::render( flags );
}//void render( Wt::WFlags<Wt::RenderFlag> flags )


void RelActAutoGui::handleDisplayedSpectrumChange( SpecUtils::SpectrumType type )
{
  switch( type )
  {
    case SpecUtils::SpectrumType::Foreground:
    case SpecUtils::SpectrumType::Background:
      break;
      
    case SpecUtils::SpectrumType::SecondForeground:
      return;
  }//switch( type )
  
  const auto fore = m_interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  const auto back = m_interspec->displayedHistogram( SpecUtils::SpectrumType::Background );
  const double backsf = m_interspec->displayScaleFactor( SpecUtils::SpectrumType::Background );
  
  const bool enableBackSub = (back && fore);
  if( enableBackSub != m_background_subtract->isEnabled() )
    m_background_subtract->setEnabled( enableBackSub );
  
  if( (fore == m_foreground) && (back == m_background) && (!back || (backsf == m_background_sf)) )
    return;
  
  
  if( fore != m_foreground )
  {
    m_cached_drf.reset();
    m_cached_all_peaks.clear();
    
    m_render_flags |= RenderActions::ChartToDefaultRange;
  }//if( fore != m_foreground )
  
  
  m_render_flags |= RenderActions::UpdateSpectra;
  scheduleRender();
}//void handleDisplayedSpectrumChange( SpecUtils::SpectrumType )


void RelActAutoGui::checkIfInUserConfigOrCreateOne( const bool force_create )
{
  if( m_loading_preset )
    return;
  
  const int index = m_presets->currentIndex();
  if( !force_create && (m_current_preset_index >= static_cast<int>(m_preset_paths.size())) )
  {
    // We are in a user-modified state, go ahead and return
    return;
  }
  
  string name;
  if( force_create )
    name = "Custom";
  else if( m_current_preset_index == 0 )
    name = "User Created";
  else
    name = "Modified " + m_presets->itemText(m_current_preset_index).toUTF8();
  
  m_presets->addItem( WString::fromUTF8(name) );
  
  m_current_preset_index = m_presets->count() - 1;
  m_presets->setCurrentIndex( m_current_preset_index );
}//void checkIfInUserConfigOrCreateOne()


RelActCalcAuto::Options RelActAutoGui::getCalcOptions() const
{
  RelActCalcAuto::Options options;
  
  
  options.fit_energy_cal = m_fit_energy_cal->isChecked();
  options.fwhm_form = RelActCalcAuto::FwhmForm( std::max(0,m_fwhm_eqn_form->currentIndex()) );
  options.fwhm_estimation_method = RelActCalcAuto::FwhmEstimationMethod( std::max(0,m_fwhm_estimation_method->currentIndex()) );
  if( options.fwhm_estimation_method == RelActCalcAuto::FwhmEstimationMethod::FixedToDetectorEfficiency )
    options.fwhm_form = RelActCalcAuto::FwhmForm::NotApplicable;
  
  const shared_ptr<const SpecUtils::Measurement> fore
                           = m_interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  const shared_ptr<const SpecMeas> meas
                                   = m_interspec->measurment( SpecUtils::SpectrumType::Foreground );
  if( fore && !fore->title().empty() )
    options.spectrum_title = fore->title();
  else if( meas && !meas->filename().empty() )
    options.spectrum_title = meas->filename();
  
  options.skew_type = PeakDef::SkewType::NoSkew;
  const int skew_index = m_skew_type->currentIndex();
  if( (skew_index >= 0) && (skew_index <= PeakDef::SkewType::DoubleSidedCrystalBall) )
    options.skew_type = PeakDef::SkewType( m_skew_type->currentIndex() );
  
  options.additional_br_uncert = -1.0;
  const auto add_uncert = RelActAutoGui::AddUncert(m_add_uncert->currentIndex());
  switch( add_uncert )
  {
    case AddUncert::StatOnly:           options.additional_br_uncert = 0.00;   break;
    case AddUncert::OneHundrethPercent: options.additional_br_uncert = 0.0001; break;
    case AddUncert::OneTenthPercent:    options.additional_br_uncert = 0.001;  break;
    case AddUncert::OnePercent:         options.additional_br_uncert = 0.01;   break;
    case AddUncert::FivePercent:        options.additional_br_uncert = 0.05;   break;
    case AddUncert::TenPercent:         options.additional_br_uncert = 0.10;   break;
    case AddUncert::TwentyFivePercent:  options.additional_br_uncert = 0.25;   break;
    case AddUncert::FiftyPercent:       options.additional_br_uncert = 0.50;   break;
    case AddUncert::SeventyFivePercent: options.additional_br_uncert = 0.75;   break;
    case AddUncert::OneHundredPercent:  options.additional_br_uncert = 1.00;   break;
    case AddUncert::NumAddUncert:       assert( 0 );                           break;
  }//switch( add_uncert )
  assert( options.additional_br_uncert >= 0.0 );
  options.additional_br_uncert = std::max( options.additional_br_uncert, 0.0 );
  
  options.floating_peaks = getFloatingPeaks();
  options.rois = getRoiRanges();

  const int num_rel_eff_curves = m_rel_eff_opts_menu->count();
  assert( num_rel_eff_curves > 0 );
  for( int rel_eff_curve_index = 0; rel_eff_curve_index < num_rel_eff_curves; ++rel_eff_curve_index )
  {
    const RelActAutoGuiRelEffOptions *opts = getRelEffCurveOptions( rel_eff_curve_index );
    if( !opts )
      throw runtime_error( "Failed to get RelActAutoGuiRelEffOptions" );

    RelActCalcAuto::RelEffCurveInput rel_eff_curve;
    rel_eff_curve.name = opts->name().toUTF8();
    rel_eff_curve.nuclides = getNucInputInfo( rel_eff_curve_index );
    rel_eff_curve.act_ratio_constraints = getActRatioConstraints( rel_eff_curve_index );
    rel_eff_curve.mass_fraction_constraints = getMassFractionConstraints( rel_eff_curve_index );
    rel_eff_curve.nucs_of_el_same_age = m_same_z_age->isChecked();
    rel_eff_curve.rel_eff_eqn_type = opts->rel_eff_eqn_form();
    rel_eff_curve.rel_eff_eqn_order = opts->rel_eff_eqn_order();
    rel_eff_curve.phys_model_self_atten = opts->phys_model_self_atten();
    rel_eff_curve.phys_model_external_atten = opts->phys_model_external_atten();
    rel_eff_curve.phys_model_use_hoerl = opts->phys_model_use_hoerl();
    rel_eff_curve.pu242_correlation_method = opts->pu242_correlation_method();

    options.rel_eff_curves.push_back( rel_eff_curve );
  }//for( int rel_eff_curve_index = 0; rel_eff_curve_index < num_rel_eff_curves; ++rel_eff_curve_index )

  return options;
}//RelActCalcAuto::Options getCalcOptions() const


vector<RelActCalcAuto::NucInputInfo> RelActAutoGui::getNucInputInfo( const int rel_eff_curve_index ) const
{
  vector<RelActCalcAuto::NucInputInfo> answer;
  
  const vector<WWidget *> &kids = m_nuclides->children();
  for( WWidget *w : kids )
  {
    const RelActAutoNuclide *nuc = dynamic_cast<const RelActAutoNuclide *>( w );
    assert( nuc );
    if( nuc && !std::holds_alternative<std::monostate>(nuc->nuclide()) )
      answer.push_back( nuc->toNucInputInfo() );
  }//for( WWidget *w : kids )
  
  return answer;
}//RelActCalcAuto::NucInputInfo getNucInputInfo() const


vector<RelActCalcAuto::RelEffCurveInput::ActRatioConstraint> RelActAutoGui::getActRatioConstraints( const int rel_eff_curve_index ) const
{
  vector<RelActCalcAuto::RelEffCurveInput::ActRatioConstraint> answer;
  
  // blah blah blah - need to implement this
  return answer;
}//vector<RelActCalcAuto::RelEffCurveInput::ActRatioConstraint> getActRatioConstraints() const


vector<RelActCalcAuto::RelEffCurveInput::MassFractionConstraint> RelActAutoGui::getMassFractionConstraints( const int rel_eff_curve_index ) const
{
  vector<RelActCalcAuto::RelEffCurveInput::MassFractionConstraint> answer;
  
  // blah blah blah - need to implement this

  return answer;
}//vector<RelActCalcAuto::MassFractionConstraint> getMassFractionConstraints() const


vector<RelActCalcAuto::RoiRange> RelActAutoGui::getRoiRanges() const
{
  vector<RelActCalcAuto::RoiRange> answer;
  const vector<WWidget *> &kids = m_energy_ranges->children();
  for( WWidget *w : kids )
  {
    const RelActAutoEnergyRange *roi = dynamic_cast<const RelActAutoEnergyRange *>( w );
    assert( roi );
    if( roi && !roi->isEmpty() )
      answer.push_back( roi->toRoiRange() );
  }//for( WWidget *w : kids )
  
  return answer;
}//RelActCalcAuto::RoiRange getRoiRanges() const


vector<RelActCalcAuto::FloatingPeak> RelActAutoGui::getFloatingPeaks() const
{
  // We will only return peaks within defined ROIs.
  vector<pair<float,float>> rois;
  const vector<WWidget *> &roi_widgets = m_energy_ranges->children();
  for( WWidget *w : roi_widgets )
  {
    const RelActAutoEnergyRange *roi = dynamic_cast<const RelActAutoEnergyRange *>( w );
    assert( roi );
    if( roi && !roi->isEmpty() )
      rois.emplace_back( roi->lowerEnergy(), roi->upperEnergy() );
  }//for( WWidget *w : kids )
  
  
  vector<RelActCalcAuto::FloatingPeak> answer;
  const std::vector<WWidget *> &kids = m_free_peaks->children();
  for( WWidget *w : kids )
  {
    RelActFreePeak *free_peak = dynamic_cast<RelActFreePeak *>(w);
    assert( free_peak );
    if( !free_peak )
      continue;
    
    const float energy = free_peak->energy();
    bool in_roi = false;
    for( const auto &roi : rois )
      in_roi = (in_roi || ((energy >= roi.first) && (energy <= roi.second)));
    
    if( !in_roi )
      continue;
    
    RelActCalcAuto::FloatingPeak peak;
    peak.energy = energy;
    peak.release_fwhm = !free_peak->fwhmConstrained();
    peak.apply_energy_cal_correction = free_peak->applyEnergyCal();
    
    answer.push_back( peak );
  }//for( loop over RelActFreePeak widgets )
  
  return answer;
}//RelActCalcAuto::FloatingPeak getFloatingPeaks() const


shared_ptr<const RelActCalcAuto::RelActAutoSolution> RelActAutoGui::getCurrentSolution() const
{
  return m_solution;
}


void RelActAutoGui::handleRoiDrag( double new_roi_lower_energy,
                   double new_roi_upper_energy,
                   double new_roi_lower_px,
                   double new_roi_upper_px,
                   const double original_roi_lower_energy,
                   const bool is_final_range )
{
  //cout << "RelActAutoGui::handleRoiDrag: original_roi_lower_energy=" << original_roi_lower_energy
  //<< ", new_roi_lower_energy=" << new_roi_lower_energy << ", new_roi_upper_energy=" << new_roi_upper_energy
  //<< ", is_final_range=" << is_final_range << endl;
  
  double min_de = 999999.9;
  RelActAutoEnergyRange *range = nullptr;
  
  const vector<WWidget *> &kids = m_energy_ranges->children();
  for( WWidget *w : kids )
  {
    RelActAutoEnergyRange *roi = dynamic_cast<RelActAutoEnergyRange *>( w );
    assert( roi );
    if( !roi || roi->isEmpty() )
      continue;
    
    RelActCalcAuto::RoiRange roi_range = roi->toRoiRange();
    
    const double de = fabs(roi_range.lower_energy - original_roi_lower_energy);
    if( de < min_de )
    {
      min_de = de;
      range = roi;
    }
  }//for( WWidget *w : kids )

  if( !range || (min_de > 2.5) )  // Sometimes the ROI might say its original lower energy is like 603 keV, but the RelActAutoEnergyRange might say 604.2.
  {
    cerr << "Unexpectedly couldnt find ROI in getRoiRanges()!" << endl;
    cout << "\t\toriginal_roi_lower_energy=" << original_roi_lower_energy
    << ", new_roi_lower_energy=" << new_roi_lower_energy << ", new_roi_upper_energy=" << new_roi_upper_energy
    << ", is_final_range=" << is_final_range << endl;
    for( WWidget *w : kids )
    {
      RelActAutoEnergyRange *roi = dynamic_cast<RelActAutoEnergyRange *>( w );
      assert( roi );
      if( !roi || roi->isEmpty() )
        continue;
      
      RelActCalcAuto::RoiRange roi_range = roi->toRoiRange();
      cout << "\t\tRange: " << roi_range.lower_energy << ", " << roi_range.upper_energy << endl;
    }
    cout << endl << endl;
    
    return;
  }//if( failed to find continuum )
  
  if( is_final_range && (new_roi_upper_px < new_roi_lower_px) )
  {
    handleRemoveEnergy( range );
    return;
  }//if( the user
  
  // We will only set RelActAutoEnergyRange energies when its final, otherwise the lower energy wont
  //  match on later updates
  if( is_final_range )
  {
    range->setEnergyRange( new_roi_lower_energy, new_roi_upper_energy );
    
    handleEnergyRangeChange();
  }else
  {
    shared_ptr<const deque<PeakModel::PeakShrdPtr>> origpeaks = m_peak_model->peaks();
    if( !origpeaks )
      return;
    
    double minDe = 999999.9;
    std::shared_ptr<const PeakContinuum> continuum;
    for( auto p : *origpeaks )
    {
      const double de = fabs( p->continuum()->lowerEnergy() - original_roi_lower_energy );
      if( de < min_de )
      {
        min_de = de;
        continuum = p->continuum();
      }
    }//for( auto p : *origpeaks )
    
    if( !continuum || min_de > 1.0 )  //0.001 would probably be fine instead of 1.0
    {
      m_spectrum->updateRoiBeingDragged( {} );
      return;
    }//if( failed to find continuum )
    
    
    auto new_continuum = std::make_shared<PeakContinuum>( *continuum );
    
    // re-use the c++ ROI value that we arent dragging, to avoid rounding or whatever
    const bool dragginUpperEnergy = (fabs(new_roi_lower_energy - continuum->lowerEnergy())
                                      < fabs(new_roi_upper_energy - continuum->upperEnergy()));
    
    if( dragginUpperEnergy )
      new_roi_lower_energy = continuum->lowerEnergy();
    else
      new_roi_upper_energy = continuum->upperEnergy();
    
    new_continuum->setRange( new_roi_lower_energy, new_roi_upper_energy );
    
    vector< shared_ptr<const PeakDef>> new_roi_initial_peaks;
    for( auto p : *origpeaks )
    {
      if( p->continuum() == continuum )
      {
        auto newpeak = make_shared<PeakDef>(*p);
        newpeak->setContinuum( new_continuum );
        new_roi_initial_peaks.push_back( newpeak );
      }
    }//for( auto p : *origpeaks )
    
    m_spectrum->updateRoiBeingDragged( new_roi_initial_peaks );
  }//if( is_final_range )
}//void handleRoiDrag(...)


void RelActAutoGui::handleCreateRoiDrag( const double lower_energy,
                         const double upper_energy,
                         const int num_peaks_to_force,
                         const bool is_final_range )
{
  //cout << "RelActAutoGui::handleCreateRoiDrag: lower_energy=" << lower_energy
  //<< ", upper_energy=" << upper_energy << ", num_peaks_to_force=" << num_peaks_to_force
  //<< ", is_final_range=" << is_final_range << endl;
  
  if( !m_foreground )
  {
    m_spectrum->updateRoiBeingDragged( {} );
    return;
  }
  
  const double roi_lower = std::min( lower_energy, upper_energy );
  const double roi_upper = std::max( lower_energy, upper_energy );
  
  if( is_final_range )
  {
    m_spectrum->updateRoiBeingDragged( {} );
    
    RelActAutoEnergyRange *energy_range = new RelActAutoEnergyRange( this, m_energy_ranges );
    
    energy_range->updated().connect( this, &RelActAutoGui::handleEnergyRangeChange );
    
    energy_range->remove().connect( boost::bind( &RelActAutoGui::handleRemoveEnergy,
                                                this, static_cast<WWidget *>(energy_range) ) );
    
    energy_range->setEnergyRange( roi_lower, roi_upper );
    energy_range->setForceFullRange( true );
    
    handleEnergyRangeChange();
    
    return;
  }//if( is_final_range )
  
  
  try
  {
    // Make a single peak with zero amplitude to at least give feedback (the continuum range line)
    //  to the user about where they have dragged.
    const double mean = 0.5*(roi_lower + roi_upper);
    const double sigma = 0.5*fabs( roi_upper - roi_lower );
    const double amplitude = 0.0;
    
    auto peak = make_shared<PeakDef>(mean, sigma, amplitude );
    
    peak->continuum()->setRange( roi_lower, roi_upper );
    peak->continuum()->calc_linear_continuum_eqn( m_foreground, mean, roi_lower, roi_upper, 3, 3 );
    
    m_spectrum->updateRoiBeingDragged( vector< shared_ptr<const PeakDef>>{peak} );
  }catch( std::exception &e )
  {
    m_spectrum->updateRoiBeingDragged( {} );
    cerr << "RelActAutoGui::handleCreateRoiDrag caught exception: " << e.what() << endl;
    return;
  }//try / catch
}//void handleCreateRoiDrag(...)


void RelActAutoGui::handleShiftDrag( const double lower_energy, const double upper_energy )
{
  const vector<WWidget *> kids = m_energy_ranges->children();
  for( WWidget *w : kids )
  {
    RelActAutoEnergyRange *roi = dynamic_cast<RelActAutoEnergyRange *>( w );
    assert( roi );
    if( !roi || roi->isEmpty() )
      continue;
    
    const double roi_lower = roi->lowerEnergy();
    const double roi_upper = roi->upperEnergy();
    
    // If the ranges intersect, deal with it
    if( (upper_energy >= roi_lower) && (lower_energy <= roi_upper) )
      handleRemovePartOfEnergyRange( w, lower_energy, upper_energy );
  }//for( WWidget *w : kids )
}//void handleShiftDrag( const double lower_energy, const double upper_energy )


void RelActAutoGui::handleDoubleLeftClick( const double energy, const double /* counts */ )
{
  try
  {
    char buffer[256] = { '\0' };
    
    // Check if click was in a ROI, and if so ignore it
    const vector<RelActCalcAuto::RoiRange> orig_rois = getRoiRanges();
    for( const RelActCalcAuto::RoiRange &roi : orig_rois )
    {
      if( (energy > roi.lower_energy) && (energy < roi.upper_energy) )
      {
        snprintf( buffer, sizeof(buffer),
                 "%.1f keV is already in the energy range [%.1f, %.1f] keV - no action will be taken.",
                 energy, roi.lower_energy, roi.upper_energy );
        
        passMessage( buffer, WarningWidget::WarningMsgMedium );
        return;
      }
    }//for( const RelActCalcAuto::RoiRange &roi : orig_rois )
    
    // If we're here, the double-click was not in an existing ROI.
    if( !m_foreground )
      return;
    
    const double xmin = m_spectrum->xAxisMinimum();
    const double xmax = m_spectrum->xAxisMaximum();
    
    const double specWidthPx = m_spectrum->chartWidthInPixels();
    const double pixPerKeV = (xmax > xmin && xmax > 0.0 && specWidthPx > 10.0) ? std::max(0.001,(specWidthPx/(xmax - xmin))): 0.001;
    
    // We'll prefer the DRF from m_solution, to what the user has picked
    shared_ptr<const DetectorPeakResponse> det = m_solution ? m_solution->m_drf : nullptr;
    
    if( !det || !det->hasResolutionInfo()  )
    {
      shared_ptr<const SpecMeas> meas = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
      det = meas ? meas->detector() : nullptr;
    }//if( solution didnt have DRF )
    
    const auto found_peaks = searchForPeakFromUser( energy, pixPerKeV, m_foreground, {}, det );
    
    // If we didnt fit a peak, and we dont
    double lower_energy = energy - 10;
    double upper_energy = energy + 10;
    shared_ptr<const PeakDef> peak = found_peaks.first.empty() ? nullptr : found_peaks.first.front();
    if( peak )
    {
      // Found a peak
      lower_energy = peak->lowerX();
      upper_energy = peak->upperX();
    }else if( det && det->hasResolutionInfo() )
    {
      // No peak found, but we have a DRF with FWHM info
      const float fwhm = std::max( 1.0f, det->peakResolutionFWHM(energy) );
      
      // We'll make a ROI that is +-4 FWHM (arbitrary choice)
      lower_energy = energy - 4*fwhm;
      upper_energy = energy + 4*fwhm;
    }else
    {
      // If we're here, were desperate - lets pluck some values out of the air
      //  (We could check on auto-searched peaks and estimate from those, but
      //   really, at this point, its not worth the effort as things are probably
      //   low quality)
      const bool highres = PeakFitUtils::is_high_res(m_foreground);
      if( highres )
      {
        lower_energy = energy - 5;
        upper_energy = energy + 5;
      }else
      {
        lower_energy = 0.8*energy;
        upper_energy = 1.2*energy;
      }
    }//if( peak ) / else if( det ) / else
    
    RelActCalcAuto::RoiRange new_roi;
    new_roi.lower_energy = lower_energy;
    new_roi.upper_energy = upper_energy;
    new_roi.continuum_type = PeakContinuum::OffsetType::Linear;
    new_roi.force_full_range = true;
    new_roi.allow_expand_for_peak_width = false;
    
    
    // Now check if we overlap with a ROI, or perhaps we need to expand an existing ROI
    
    const vector<WWidget *> prev_roi_widgets = m_energy_ranges->children();
    
    RelActAutoEnergyRange *new_roi_w = new RelActAutoEnergyRange( this, m_energy_ranges );
    new_roi_w->updated().connect( this, &RelActAutoGui::handleEnergyRangeChange );
    new_roi_w->remove().connect( boost::bind( &RelActAutoGui::handleRemoveEnergy, this, static_cast<WWidget *>(new_roi_w) ) );
    new_roi_w->setFromRoiRange( new_roi );
    
    vector<RelActAutoEnergyRange *> overlapping_rois;
    for( WWidget *w : prev_roi_widgets )
    {
      RelActAutoEnergyRange *roi = dynamic_cast<RelActAutoEnergyRange *>( w );
      assert( roi );
      if( !roi )
        continue;
      if( (new_roi.upper_energy > roi->lowerEnergy()) && (new_roi.lower_energy < roi->upperEnergy()) )
        overlapping_rois.push_back( roi );
    }
    
    bool combined_an_roi = false;
    for( RelActAutoEnergyRange *roi : overlapping_rois )
    {
      WWidget *w = handleCombineRoi( roi, new_roi_w );
      RelActAutoEnergyRange *combined_roi = dynamic_cast<RelActAutoEnergyRange *>( w );
      
      assert( combined_roi );
      if( !combined_roi )
        break;
      
      combined_an_roi = true;
      new_roi_w = combined_roi;
    }
    
    lower_energy = new_roi_w->lowerEnergy();
    upper_energy = new_roi_w->upperEnergy();
    
    if( combined_an_roi )
    {
      snprintf( buffer, sizeof(buffer),
               "Extended existing energy range to [%.1f, %.1f] keV.",
               lower_energy, upper_energy );
    }else
    {
      snprintf( buffer, sizeof(buffer),
               "Added a new new energy range from %.1f to %.1f keV.",
               lower_energy, upper_energy );
    }//if( combined_an_roi ) / else
        
    passMessage( buffer, WarningWidget::WarningMsgLow );
    
    checkIfInUserConfigOrCreateOne( false );
    m_render_flags |= RenderActions::UpdateEnergyRanges;
    m_render_flags |= RenderActions::UpdateCalculations;
    scheduleRender();
  }catch( std::exception &e )
  {
    passMessage( "handleDoubleLeftClick error: " + string(e.what()), WarningWidget::WarningMsgHigh );
  }//try / catch
}//void handleDoubleLeftClick( const double energy, const double counts )


void RelActAutoGui::handleRightClick( const double energy, const double counts,
                      const int page_x_px, const int page_y_px )
{
  vector<RelActAutoEnergyRange *> ranges;
  const vector<WWidget *> &kids = m_energy_ranges->children();
  for( WWidget *w : kids )
  {
    RelActAutoEnergyRange *roi = dynamic_cast<RelActAutoEnergyRange *>( w );
    assert( roi );
    if( roi && !roi->isEmpty() )
       ranges.push_back( roi );
  }//for( WWidget *w : kids )
  
  std::sort( begin(ranges), end(ranges),
             []( const RelActAutoEnergyRange *lhs, const RelActAutoEnergyRange *rhs) -> bool{
    return lhs->lowerEnergy() < rhs->lowerEnergy();
  } );
  
  
  RelActAutoEnergyRange *range = nullptr, *range_to_left = nullptr, *range_to_right = nullptr;
  for( size_t i = 0; i < ranges.size(); ++i )
  {
    RelActAutoEnergyRange *roi = ranges[i];
    RelActCalcAuto::RoiRange roi_range = roi->toRoiRange();
  
    if( (energy >= roi_range.lower_energy) && (energy < roi_range.upper_energy) )
    {
      range = roi;
      if( i > 0 )
        range_to_left = ranges[i-1];
      if( (i + 1) < ranges.size() )
        range_to_right = ranges[i+1];
      
      break;
    }//
  }//for( RelActAutoEnergyRange *roi : ranges )
  
  if( !range )
  {
    cerr << "Right-click is not in an ROI" << endl;
    return;
  }//if( failed to find continuum )
  
  
  DeleteOnClosePopupMenu *menu = new DeleteOnClosePopupMenu( nullptr, PopupDivMenu::TransientMenu );
  menu->aboutToHide().connect( menu, &DeleteOnClosePopupMenu::markForDelete );
  menu->setPositionScheme( Wt::Absolute );
  
  PopupDivMenuItem *item = nullptr;
  
  const bool is_phone = false; //isPhone();
  if( is_phone )
    item = menu->addPhoneBackItem( nullptr );
  
  item = menu->addMenuItem( "ROI options:" );
  item->disable();
  item->setSelectable( false );
  menu->addSeparator();
  
  const auto roi = range->toRoiRange();
  
  PopupDivMenu *continuum_menu = menu->addPopupMenuItem( "Set Continuum Type" );
  for( auto type = PeakContinuum::OffsetType(0);
      type < PeakContinuum::External; type = PeakContinuum::OffsetType(type+1) )
  {
    WMenuItem *item = continuum_menu->addItem( WString::tr(PeakContinuum::offset_type_label_tr(type)) );
    item->triggered().connect( boost::bind( &RelActAutoEnergyRange::setContinuumType, range, type ) );
    if( type == roi.continuum_type )
      item->setDisabled( true );
  }//for( loop over PeakContinuum::OffsetTypes )
  
  item = menu->addMenuItem( "Remove ROI" );
  item->triggered().connect( boost::bind( &RelActAutoGui::handleRemoveEnergy, this, static_cast<WWidget *>(range) ) );
  
  char buffer[128] = { '\0' };
  snprintf( buffer, sizeof(buffer), "Split ROI at %.1f keV", energy );
  item = menu->addMenuItem( buffer );
  item->triggered().connect( boost::bind( &RelActAutoGui::handleSplitEnergyRange, this, static_cast<WWidget *>(range), energy ) );
  
  const char *item_label = "";
  if( roi.force_full_range )
    item_label = "Don't force full-range";
  else
    item_label = "Force full-range";
  item = menu->addMenuItem( item_label );
  item->triggered().connect( boost::bind( &RelActAutoGui::handleToggleForceFullRange, this, static_cast<WWidget *>(range) ) );
  
  // TODO: we could be a little more intelligent about when offering to combine ROIs
  if( range_to_left )
  {
    item = menu->addMenuItem( "Combine with ROI to left" );
    item->triggered().connect( boost::bind( &RelActAutoGui::handleCombineRoi, this,
                                           static_cast<WWidget *>(range_to_left),
                                           static_cast<WWidget *>(range) ) );
  }//if( range_to_left )
  
  if( range_to_right )
  {
    item = menu->addMenuItem( "Combine with ROI to right" );
    item->triggered().connect( boost::bind( &RelActAutoGui::handleCombineRoi, this,
                                           static_cast<WWidget *>(range),
                                           static_cast<WWidget *>(range_to_right) ) );
  }//if( range_to_right )
  
  // TODO: Add floating peak item
  snprintf( buffer, sizeof(buffer), "Add free peak at %.1f keV", energy );
  item = menu->addMenuItem( buffer );
  item->triggered().connect( boost::bind( &RelActAutoGui::handleAddFreePeak, this, energy, true, true ) );
  
  
  if( is_phone )
  {
    menu->addStyleClass( " Wt-popupmenu Wt-outset" );
    menu->showMobile();
  }else
  {
    menu->addStyleClass( " Wt-popupmenu Wt-outset RelActAutoGuiContextMenu" );
    menu->popup( WPoint(page_x_px - 30, page_y_px - 30) );
  }
}//void handleRightClick(...)


void RelActAutoGui::setCalcOptionsGui( const RelActCalcAuto::Options &options )
{
  assert( options.rel_eff_curves.size() >= 1 );
  if( options.rel_eff_curves.empty() )
    throw runtime_error( "RelActAutoGui::setCalcOptionsGui: for dev, must have at least one rel-eff curve." );
  
  m_fit_energy_cal->setChecked( options.fit_energy_cal );
  
  m_fwhm_estimation_method->setCurrentIndex( static_cast<int>(options.fwhm_estimation_method) );
  const bool fixed_to_det_eff = (options.fwhm_estimation_method == RelActCalcAuto::FwhmEstimationMethod::FixedToDetectorEfficiency);
  
  m_fwhm_eqn_form->setHidden( fixed_to_det_eff );
  if( m_fwhm_eqn_form->label() )
    m_fwhm_eqn_form->label()->setHidden( fixed_to_det_eff );
  if( !fixed_to_det_eff && (options.fwhm_form != RelActCalcAuto::FwhmForm::NotApplicable) )
    m_fwhm_eqn_form->setCurrentIndex( static_cast<int>(options.fwhm_form) );
  
  m_skew_type->setCurrentIndex( static_cast<int>(options.skew_type) );
  
  // We'll just round add-uncert to the nearest-ish value we allow in the GUI
  RelActAutoGui::AddUncert add_uncert = AddUncert::NumAddUncert;
  if( options.additional_br_uncert <= 0.00005 )
    add_uncert = AddUncert::StatOnly;
  else if( options.additional_br_uncert <= 0.0005 )
    add_uncert = AddUncert::OneHundrethPercent;
  else if( options.additional_br_uncert <= 0.005 )
    add_uncert = AddUncert::OneTenthPercent;
  else if( options.additional_br_uncert <= 0.025 )
    add_uncert = AddUncert::OnePercent;
  else if( options.additional_br_uncert <= 0.075 )
    add_uncert = AddUncert::FivePercent;
  else if( options.additional_br_uncert <= 0.175 )
    add_uncert = AddUncert::TenPercent;
  else if( options.additional_br_uncert <= 0.375 )
    add_uncert = AddUncert::TwentyFivePercent;
  else if( options.additional_br_uncert <= 0.625 )
    add_uncert = AddUncert::FiftyPercent;
  else if( options.additional_br_uncert <= 0.875 )
    add_uncert = AddUncert::SeventyFivePercent;
  else
    add_uncert = AddUncert::OneHundredPercent;
  
  assert( add_uncert != AddUncert::NumAddUncert );
  if( add_uncert == AddUncert::NumAddUncert )
    add_uncert = AddUncert::StatOnly;
  m_add_uncert->setCurrentIndex( static_cast<int>(add_uncert) );
  
  m_nuclides->clear();

  // First, remove any extra Rel Eff curve GUIs
  const size_t num_rel_eff_curves = options.rel_eff_curves.size();
  while( m_rel_eff_opts_menu->count() > num_rel_eff_curves )
  {
    RelActAutoGuiRelEffOptions *opts = getRelEffCurveOptions( m_rel_eff_opts_menu->count() - 1 );
    assert( opts );
    handleDelRelEffCurve( opts );
  }
  
  // Add any new Rel Eff curve GUIs we need
  while( m_rel_eff_opts_menu->count() < num_rel_eff_curves )
  {
    handleAddRelEffCurve();
  }

  assert( m_rel_eff_opts_menu->count() == static_cast<int>(num_rel_eff_curves) );

  // Now, set the values for each Rel Eff curve
  for( size_t i = 0; i < num_rel_eff_curves; ++i )
  {
    const RelActCalcAuto::RelEffCurveInput &rel_eff = options.rel_eff_curves[i];
    RelActAutoGuiRelEffOptions *opts = getRelEffCurveOptions( i );
    assert( opts );
    opts->setRelEffCurveInput( rel_eff );
  }

  // We'll set the "Same El. Same Age" checkbox based on whether any of the
  //  Rel Eff curves had this option enabled.
  bool nucs_of_el_same_age = false;
  for( size_t i = 0; i < num_rel_eff_curves; ++i )
    nucs_of_el_same_age |= options.rel_eff_curves[i].nucs_of_el_same_age;
  m_same_z_age->setChecked( nucs_of_el_same_age );

  
  // Update the nuclide widgets
  //  If there are multiple Rel Eff curves, they may share identical nuclides, 
  //  in which case we will only show the nuclide once in the GUI.
  //  So we'll start off by making a copy of the nuclides, for each Rel Eff curve,
  //  and then we'll remove duplicates.
  vector<vector<RelActCalcAuto::NucInputInfo>> all_nuclides( num_rel_eff_curves );
  for( size_t i = 0; i < num_rel_eff_curves; ++i )
    all_nuclides[i] = options.rel_eff_curves[i].nuclides;

  // Now, remove duplicates
  for( size_t curve_index = 0; curve_index < num_rel_eff_curves; ++curve_index )
  {
    const vector<RelActCalcAuto::NucInputInfo> &nuclides = all_nuclides[curve_index];
    for( size_t nuc_index = 0; nuc_index < nuclides.size(); ++nuc_index )
    {
      const RelActCalcAuto::NucInputInfo &nuc = nuclides[nuc_index];
      
      bool is_in_all_curves = true;
      set<size_t> curves_with_nuc;
      curves_with_nuc.insert( curve_index );  // We know this curve has the nuc
      for( size_t other_curve_index = curve_index+1; other_curve_index < num_rel_eff_curves; ++other_curve_index )
      {
        const vector<RelActCalcAuto::NucInputInfo> &other_nuclides = all_nuclides[other_curve_index];
        const auto pos = std::find( begin(other_nuclides), end(other_nuclides), nuc );
        if( pos == end(other_nuclides) )
          is_in_all_curves = false;
        else  
         curves_with_nuc.insert( other_curve_index );
      }//for( size_t other_curve_index = curve_index+1; other_curve_index < num_rel_eff_curves; ++other_curve_index )
      
      // If the exact nuclide is in other curves, we'll delete it from all other ones, so we will only add it to the GUI once
      for( size_t other_curve_index : curves_with_nuc )
      {
        vector<RelActCalcAuto::NucInputInfo> &other_nuclides = all_nuclides[other_curve_index];
        auto pos = std::find( begin(other_nuclides), end(other_nuclides), nuc );
        assert( pos != end(other_nuclides) );
        if( pos != end(other_nuclides) )
          other_nuclides.erase( pos );
      }
      
      RelActAutoNuclide *nuc_widget = new RelActAutoNuclide( this, m_nuclides );
      nuc_widget->updated().connect( this, &RelActAutoGui::handleNuclidesChanged );
      nuc_widget->remove().connect( boost::bind( &RelActAutoGui::handleRemoveNuclide,
                                              this, static_cast<WWidget *>(nuc_widget) ) );
      nuc_widget->fromNucInputInfo( nuc );

      if( nuc.nuclide )
      {
        assert( curve_index < options.rel_eff_curves.size() );
        const RelActCalcAuto::RelEffCurveInput &rel_eff = options.rel_eff_curves[curve_index];
        for( const RelActCalcAuto::RelEffCurveInput::ActRatioConstraint &constraint : rel_eff.act_ratio_constraints )
        {
          if( constraint.constrained_nuclide == nuc.nuclide )
            nuc_widget->addActRatioConstraint( constraint );
        }//for( const auto &constraint : nuc.act_ratio_constraints )

        for( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &constraint : rel_eff.mass_fraction_constraints )
        {
          if( constraint.nuclide == nuc.nuclide )
            nuc_widget->addMassFractionConstraint( constraint );
        }//for( const auto &constraint : nuc.mass_fraction_constraints )  
      }//if( nuc.nuclide )

      nuc_widget->setIsInCurves( curves_with_nuc, num_rel_eff_curves );
    }//for( const size_t nuc_index = 0; nuc_index < nuclides.size(); ++nuc_index )
  }//for( loop over curve_index )

  
  m_energy_ranges->clear();
  for( const RelActCalcAuto::RoiRange &roi : options.rois )
  {
    RelActAutoEnergyRange *energy_range = new RelActAutoEnergyRange( this, m_energy_ranges );
    
    energy_range->updated().connect( this, &RelActAutoGui::handleEnergyRangeChange );
    
    energy_range->remove().connect( boost::bind( &RelActAutoGui::handleRemoveEnergy,
                                                this, static_cast<WWidget *>(energy_range) ) );
    
    energy_range->setFromRoiRange( roi );
  }//for( const RelActCalcAuto::RoiRange &roi : rois )
  
  
  // Free Peaks
  m_free_peaks->clear();
  if( options.floating_peaks.empty() && !m_free_peaks_container->isHidden() )
    handleHideFreePeaks();
  
  if( !options.floating_peaks.empty() && m_free_peaks_container->isHidden() )
    handleShowFreePeaks();
  
  for( const RelActCalcAuto::FloatingPeak &peak : options.floating_peaks )
    handleAddFreePeak( peak.energy, !peak.release_fwhm, peak.apply_energy_cal_correction );
  
  // options.spectrum_title
  m_render_flags |= RenderActions::UpdateNuclidesPresent;  //To trigger calling `updateDuringRenderForNuclideChange()`
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void setCalcOptionsGui( const RelActCalcAuto::Options &options )


void RelActAutoGui::handleToggleForceFullRange( Wt::WWidget *w )
{
  if( !w )
    return;
  
  const std::vector<WWidget *> &kids = m_energy_ranges->children();
  const auto pos = std::find( begin(kids), end(kids), w );
  if( pos == end(kids) )
  {
    cerr << "Failed to find a RelActAutoEnergyRange in m_energy_ranges!" << endl;
    assert( 0 );
    return;
  }
  
  assert( dynamic_cast<RelActAutoEnergyRange *>(w) );
  
  RelActAutoEnergyRange *roi = dynamic_cast<RelActAutoEnergyRange *>(w);
  assert( roi );
  if( !roi )
    return;
  
  roi->setForceFullRange( !roi->forceFullRange() );
}//void handleToggleForceFullRange( Wt::WWidget *w )


Wt::WWidget *RelActAutoGui::handleCombineRoi( Wt::WWidget *left_roi, Wt::WWidget *right_roi )
{
  if( !left_roi || !right_roi || (left_roi == right_roi) )
  {
    assert( 0 );
    return nullptr;
  }
  
  const std::vector<WWidget *> &kids = m_energy_ranges->children();
  const auto left_pos = std::find( begin(kids), end(kids), left_roi );
  const auto right_pos = std::find( begin(kids), end(kids), right_roi );
  if( (left_pos == end(kids)) || (right_pos == end(kids)) )
  {
    cerr << "Failed to find left or right RelActAutoEnergyRange in m_energy_ranges!" << endl;
    assert( 0 );
    return nullptr;
  }
  
  RelActAutoEnergyRange *left_range = dynamic_cast<RelActAutoEnergyRange *>(left_roi);
  RelActAutoEnergyRange *right_range = dynamic_cast<RelActAutoEnergyRange *>(right_roi);
  
  assert( left_range && right_range );
  if( !left_range || !right_range )
    return nullptr;
  
  const RelActCalcAuto::RoiRange lroi = left_range->toRoiRange();
  const RelActCalcAuto::RoiRange rroi = right_range->toRoiRange();
  
  RelActCalcAuto::RoiRange new_roi = lroi;
  new_roi.lower_energy = std::min( lroi.lower_energy, rroi.lower_energy );
  new_roi.upper_energy = std::max( lroi.upper_energy, rroi.upper_energy );
  if( rroi.force_full_range )
    new_roi.force_full_range = true;
  if( rroi.allow_expand_for_peak_width )
    new_roi.allow_expand_for_peak_width = true;
  new_roi.continuum_type = std::max( lroi.continuum_type, rroi.continuum_type );
  
  delete right_range;
  right_roi = nullptr;
  right_range = nullptr;
  
  left_range->setFromRoiRange( new_roi );
  
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateEnergyRanges;
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
  
  return left_range;
}//void handleCombineRoi( Wt::WWidget *left_roi, Wt::WWidget *right_roi );


rapidxml::xml_node<char> *RelActAutoGui::serialize( rapidxml::xml_node<char> *parent_node ) const
{
  RelActCalcAuto::RelActAutoGuiState state;
  state.options = getCalcOptions();
  state.background_subtract = (m_background_subtract->isEnabled() && m_background_subtract->isChecked());
  state.show_ref_lines = m_hide_ref_lines_item->isEnabled();
  state.lower_display_energy = m_spectrum->xAxisMinimum();
  state.upper_display_energy = m_spectrum->xAxisMaximum();
  
  return state.serialize( parent_node );
}//rapidxml::xml_node<char> *RelActAutoGui::serialize( rapidxml::xml_node<char> *parent )


void RelActAutoGui::deSerialize( const rapidxml::xml_node<char> *base_node )
{
  MaterialDB *materialDb = m_interspec->materialDataBase();
  
  RelActCalcAuto::RelActAutoGuiState state;
  state.deSerialize( base_node, materialDb );
  
  m_background_subtract->setChecked( state.background_subtract );
  
  m_show_ref_lines_item->setHidden( state.show_ref_lines );
  m_hide_ref_lines_item->setHidden( !state.show_ref_lines );
  m_show_ref_lines_item->setDisabled( state.show_ref_lines );
  m_hide_ref_lines_item->setDisabled( !state.show_ref_lines );
    
  if( state.lower_display_energy < state.upper_display_energy )
  {
    // Note: this next line only works when creating this widget and then loading its state
    //       because #updateDuringRenderForSpectrumChange checks if the widget has rendered
    //       yet, and if not, and it looks like a custom range has been set, then it wont
    //       reset the range.
    
    m_spectrum->setXAxisRange( state.lower_display_energy, state.upper_display_energy );
  }
    
  m_loading_preset = true;
  
  setCalcOptionsGui( state.options );
  
  m_solution.reset();
  m_peak_model->setPeaks( vector<PeakDef>{} );
  
  m_solution_updated.emit( m_solution );
  m_calc_failed.emit();
  
  m_render_flags |= RenderActions::UpdateNuclidesPresent;
  m_render_flags |= RenderActions::UpdateEnergyRanges;
  m_render_flags |= RenderActions::UpdateCalculations;
  m_render_flags |= RenderActions::UpdateFreePeaks;
  m_render_flags |= RenderActions::UpdateRefGammaLines;
  if( state.lower_display_energy >= state.upper_display_energy )
    m_render_flags |= RenderActions::ChartToDefaultRange;
  
  scheduleRender();
}//void deSerialize( const rapidxml::xml_node<char> *base_node )


std::unique_ptr<rapidxml::xml_document<char>> RelActAutoGui::guiStateToXml() const
{
  std::unique_ptr<rapidxml::xml_document<char>> doc( new rapidxml::xml_document<char>() );
  
  serialize( doc.get() );
  
  return std::move( doc );
}//std::unique_ptr<rapidxml::xml_document<char>> guiStateToXml() const


void RelActAutoGui::setGuiStateFromXml( const rapidxml::xml_document<char> *doc )
{
  if( !doc )
    throw runtime_error( "RelActAutoGui::setGuiStateFromXml: nullptr passed in." );
  
  const rapidxml::xml_node<char> *base_node = doc->first_node( "RelActCalcAuto" );
  if( !base_node )
    throw runtime_error( "RelActAutoGui::setGuiStateFromXml: couldnt find <RelActCalcAuto> node." );
  
  deSerialize( base_node );
}//void setGuiStateFromXml( const rapidxml::xml_node<char> *node );


void RelActAutoGui::handlePresetChange()
{
  // Right now this function just handles setting GUI to default values, but eventually it will
  //  check if the user is in a modified parameter set, and if so save it to memory, so it can go
  //  back to it later
  
  const int index = m_presets->currentIndex();
  if( index == m_current_preset_index )
    return;
  
  m_render_flags |= RenderActions::UpdateNuclidesPresent;
  m_render_flags |= RenderActions::UpdateEnergyRanges;
  m_render_flags |= RenderActions::UpdateCalculations;
  m_render_flags |= RenderActions::ChartToDefaultRange;
  scheduleRender();
  
  if( m_current_preset_index >= static_cast<int>(m_preset_paths.size()) )
    m_previous_presets[m_current_preset_index] = guiStateToXml();
  
  m_current_preset_index = index;
  
  m_nuclides->clear();
  m_energy_ranges->clear();
  
  while( m_rel_eff_opts_menu->count() > 1 )
  {
    RelActAutoGuiRelEffOptions *opts = getRelEffCurveOptions( m_rel_eff_opts_menu->count() - 1 );
    assert( opts );
    handleDelRelEffCurve( opts );
  } 

  if( index <= 0 )
  {
    // Clear everything out!
    return;
  }
  
  if( index >= m_preset_paths.size() )
  {
    // TODO: let users download config, or clone them
    
    const auto iter = m_previous_presets.find(index);
    if( iter == std::end(m_previous_presets) )
    {
      passMessage( "Expected state information for '" + m_presets->currentText().toUTF8()
                   + "' is not available - this is not expected - sorry!", WarningWidget::WarningMsgHigh );
      
      return;
    }//if( iter == std::end(m_previous_presets) )
    
    if( !iter->second )
    {
      passMessage( "State information was not previously able to be saved for '"
                  + m_presets->currentText().toUTF8() + "' - this is not expected - sorry!",
                  WarningWidget::WarningMsgHigh );
      
      return;
    }//if( !iter->second )
    
    try
    {
      setGuiStateFromXml( iter->second.get() );
    }catch( std::exception &e )
    {
      passMessage( "Error de-serializing tool state: " + string(e.what()),
                  WarningWidget::WarningMsgHigh );
    }
    
    return;
  }//if( index >= m_preset_paths.size() )
  
  assert( index < m_preset_paths.size() && (index > 0) );
  if( index >= m_preset_paths.size() || (index <= 0) )
    throw runtime_error( "RelActAutoGui::handlePresetChange: invalid selection index " );

  unique_ptr<rapidxml::file<char>> input_file; // define file out here to keep in scope for catch
  
  try
  {
    const string xml_path = m_preset_paths[index];
    input_file.reset( new rapidxml::file<char>( xml_path.c_str() ) );
    rapidxml::xml_document<char> doc;
    doc.parse<rapidxml::parse_trim_whitespace>( input_file->data() );
   
    setGuiStateFromXml( &doc );
  }catch( rapidxml::parse_error &e )
  {
    string msg = "Error parsing preset XML: " + string(e.what());
    const char * const position = e.where<char>();
    if( position && *position )
    {
      const char *end_pos = position;
      for( size_t i = 0; (*end_pos) && (i < 80); ++i )
        end_pos += 1;
      msg += "<br />&nbsp;&nbsp;At: " + std::string(position, end_pos);
    }//if( position )
    
    passMessage( msg, WarningWidget::WarningMsgHigh );
  }catch( std::exception &e )
  {
    passMessage( "Error loading preset: " + string(e.what()), WarningWidget::WarningMsgHigh );
  }//try / cat to read the XML
}//void RelActAutoGui::handlePresetChange()


void RelActAutoGui::handleRelEffEqnFormChanged()
{
  for( int i = 0; i < m_rel_eff_opts_stack->count(); ++i )
  {
    RelActAutoGuiRelEffOptions *opts = getRelEffCurveOptions( i );
    assert( opts );
    opts->showAndHideOptionsForEqnType();
  }
  
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleRelEffEqnFormChanged();


void RelActAutoGui::handleRelEffEqnOrderChanged()
{
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleRelEffEqnOrderChanged();


void RelActAutoGui::handleFwhmFormChanged()
{
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleFwhmFormChanged()


void RelActAutoGui::handleFwhmEstimationMethodChanged()
{
  const RelActCalcAuto::FwhmEstimationMethod index 
           = static_cast<RelActCalcAuto::FwhmEstimationMethod>( m_fwhm_estimation_method->currentIndex() );
  
  const bool fixed_to_det_eff = (index == RelActCalcAuto::FwhmEstimationMethod::FixedToDetectorEfficiency);
  if( m_fwhm_eqn_form->label() )
    m_fwhm_eqn_form->label()->setHidden( fixed_to_det_eff );
  m_fwhm_eqn_form->setHidden( fixed_to_det_eff );
  
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleFwhmEstimationMethodChanged()


void RelActAutoGui::handleFitEnergyCalChanged()
{
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateFitEnergyCal;
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleFitEnergyCalChanged();


void RelActAutoGui::handleBackgroundSubtractChanged()
{
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}


void RelActAutoGui::handleSameAgeChanged()
{
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}

void RelActAutoGui::handlePuByCorrelationChanged()
{
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}

void RelActAutoGui::handleSkewTypeChanged()
{
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}

void RelActAutoGui::handleNucDataSrcChanged()
{
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleNucDataSrcChanged()


void RelActAutoGui::handleNuclidesChanged()
{
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateNuclidesPresent;
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleNuclidesChanged()


void RelActAutoGui::handleNuclidesInfoEdited()
{
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleNuclidesInfoEdited()


void RelActAutoGui::handleEnergyRangeChange()
{
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateEnergyRanges;
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleEnergyRangeChange()


void RelActAutoGui::handleFreePeakChange()
{
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateFreePeaks;
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleFreePeakChange()


void RelActAutoGui::handleAdditionalUncertChanged()
{
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleAdditionalUncertChanged()


void RelActAutoGui::setOptionsForNoSolution()
{
  //if( !m_spectrum->chartTitle().empty() )
  //  m_spectrum->setChartTitle( WString() );
  m_fit_chi2_msg->setText( "" );
  m_fit_chi2_msg->hide();

  makeZeroAmplitudeRoisToChart();
  m_set_peaks_foreground->setDisabled( true );

  m_rel_eff_chart->setData( 0.0, {}, {}, "", "" );

  for( WWidget *w : m_energy_ranges->children() )
  {
    RelActAutoEnergyRange *roi = dynamic_cast<RelActAutoEnergyRange *>( w );
    assert( roi );
    if( roi )
      roi->enableSplitToIndividualRanges( false );
  }//for( WWidget *w : kids )
  
}//void setOptionsForNoSolution()


void RelActAutoGui::setOptionsForValidSolution()
{
  assert( m_solution && (m_solution->m_status == RelActCalcAuto::RelActAutoSolution::Status::Success) );
  if( !m_solution || (m_solution->m_status != RelActCalcAuto::RelActAutoSolution::Status::Success) )
    return;
    
  // Check if we should allow setting energy calibration from fit solution
  const bool fit_energy_cal = (m_solution->m_fit_energy_cal[0] || m_solution->m_fit_energy_cal[1]);
  m_apply_energy_cal_item->setDisabled( !fit_energy_cal );
  m_set_peaks_foreground->setDisabled( m_solution->m_fit_peaks_in_spectrums_cal.empty() );
  
  for( WWidget *w : m_energy_ranges->children() )
  {
    RelActAutoEnergyRange *roi = dynamic_cast<RelActAutoEnergyRange *>( w );
    assert( roi );
    if( !roi )
      continue;
    
    const float lower_energy = roi->lowerEnergy();
    const float upper_energy = roi->upperEnergy();
    
    // Only enable splitting the energy range if it will split into more than one sub-range
    size_t num_sub_ranges = 0;
    for( const RelActCalcAuto::RoiRange &range : m_solution->m_final_roi_ranges )
    {
      const double mid_energy = 0.5*(range.lower_energy + range.upper_energy);
      num_sub_ranges += ((mid_energy > lower_energy) && (mid_energy < upper_energy));
    }
    
    roi->enableSplitToIndividualRanges( (num_sub_ranges > 1) );
  }//for( WWidget *w : kids )
  
}//void setOptionsForValidSolution()


void RelActAutoGui::makeZeroAmplitudeRoisToChart()
{
  m_peak_model->setPeaks( vector<PeakDef>{} );
  const vector<RelActCalcAuto::RoiRange> rois = getRoiRanges();
  
  if( !m_foreground )
    return;
  
  vector<shared_ptr<const PeakDef>> peaks;
  for( const auto &roi : rois )
  {
    try
    {
      const double mean = 0.5*(roi.lower_energy + roi.upper_energy);
      const double sigma = 0.5*fabs( roi.upper_energy - roi.lower_energy );
      const double amplitude = 0.0;
      
      auto peak = make_shared<PeakDef>(mean, sigma, amplitude );
      
      peak->continuum()->setRange( roi.lower_energy, roi.upper_energy );
      peak->continuum()->calc_linear_continuum_eqn( m_foreground, mean, roi.lower_energy, roi.upper_energy, 3, 3 );
      
      peaks.push_back( peak );
    }catch( std::exception &e )
    {
      m_spectrum->updateRoiBeingDragged( {} );
      cerr << "RelActAutoGui::makeZeroAmplitudeRoisToChart caught exception: " << e.what() << endl;
      return;
    }//try / catch
  }//for( const auto &roi : rois )
  
  m_peak_model->setPeaks( peaks );
}//void RelActAutoGui::makeZeroAmplitudeRoisToChart()


void RelActAutoGui::handleAddNuclide()
{
  const vector<WWidget *> prev_kids = m_nuclides->children();
  
  RelActAutoNuclide *nuc_widget = new RelActAutoNuclide( this, m_nuclides );
  nuc_widget->updated().connect( this, &RelActAutoGui::handleNuclidesChanged );
  nuc_widget->remove().connect( boost::bind( &RelActAutoGui::handleRemoveNuclide,
                                              this, static_cast<WWidget *>(nuc_widget) ) );
  nuc_widget->setNuclideEditFocus();
  
  shared_ptr<const ColorTheme> theme = InterSpec::instance()->getColorTheme();
  if( theme )
  {
    vector<WColor> colors = theme->referenceLineColor;
    
    for( WWidget *w : prev_kids )
    {
      RelActAutoNuclide *nuc = dynamic_cast<RelActAutoNuclide *>( w );
      assert( nuc );
      if( !nuc )
        continue;
      
      const WColor color = nuc->color();
      const auto pos = std::find( begin(colors), end(colors), color );
      if( pos != end(colors) )
        colors.erase( pos );
    }
    
    if( colors.size() )
      nuc_widget->setColor( colors.front() );
  }//if( theme )
  
  
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateNuclidesPresent;
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleAddNuclide()


void RelActAutoGui::handleAddEnergy()
{
  const int nprev = m_energy_ranges->count();
  
  RelActAutoEnergyRange *energy_range = new RelActAutoEnergyRange( this, m_energy_ranges );
  
  energy_range->updated().connect( this, &RelActAutoGui::handleEnergyRangeChange );
  
  energy_range->remove().connect( boost::bind( &RelActAutoGui::handleRemoveEnergy,
                                      this, static_cast<WWidget *>(energy_range) ) );
  
  if( nprev == 0 )
  {
    const auto cal = m_foreground ? m_foreground->energy_calibration() : nullptr;
    const float upper_energy = (cal && cal->valid()) ? cal->upper_energy() : 3000.0f;
    energy_range->setEnergyRange( 125.0f, upper_energy );
  }else
  {
    energy_range->setForceFullRange( true );
  }//if( this is the first energy range ) / else
  
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateEnergyRanges;
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleAddEnergy()


void RelActAutoGui::handleClearAllEnergyRanges()
{
  SimpleDialog *dialog = new SimpleDialog( "Clear energy ranges?", "&nbsp;" );
  WPushButton *yes = dialog->addButton( "Yes" );
  dialog->addButton( "No" );
  yes->clicked().connect( this, &RelActAutoGui::removeAllEnergyRanges );
}//void handleClearAllEnergyRanges()


void RelActAutoGui::removeAllEnergyRanges()
{
  m_energy_ranges->clear();
  
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateEnergyRanges;
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void removeAllEnergyRanges()


void RelActAutoGui::handleShowFreePeaks()
{
  m_show_free_peak->hide();
  m_free_peaks_container->show();
}//void handleShowFreePeaks()


void RelActAutoGui::handleHideFreePeaks()
{
  const int nfree_peaks = m_free_peaks->count();
  m_free_peaks->clear();
  m_show_free_peak->show();
  m_free_peaks_container->hide();
  
  if( nfree_peaks )
  {
    checkIfInUserConfigOrCreateOne( false );
    m_render_flags |= RenderActions::UpdateCalculations;
    m_render_flags |= RenderActions::UpdateFreePeaks;
    scheduleRender();
  }
}//void handleHideFreePeaks()


void RelActAutoGui::handleRemoveFreePeak( Wt::WWidget *w )
{
  if( !w )
    return;
  
  const std::vector<WWidget *> &kids = m_free_peaks->children();
  const auto pos = std::find( begin(kids), end(kids), w );
  if( pos == end(kids) )
  {
    cerr << "Failed to find a RelActFreePeak in m_free_peaks!" << endl;
    assert( 0 );
    return;
  }
  
  assert( dynamic_cast<RelActFreePeak *>(w) );
  
  delete w;
  
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateEnergyRanges;
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleRemoveFreePeak( Wt::WWidget *w )


void RelActAutoGui::handleRemoveEnergy( WWidget *w )
{
  if( !w )
    return;
  
  const std::vector<WWidget *> &kids = m_energy_ranges->children();
  const auto pos = std::find( begin(kids), end(kids), w );
  if( pos == end(kids) )
  {
    cerr << "Failed to find a RelActAutoEnergyRange in m_energy_ranges!" << endl;
    assert( 0 );
    return;
  }
  
  assert( dynamic_cast<RelActAutoEnergyRange *>(w) );
  
  delete w;
  
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateEnergyRanges;
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleRemoveEnergy( Wt::WContainerWidget *w )


void RelActAutoGui::handleSplitEnergyRange( Wt::WWidget *w, const double energy )
{
  handleRemovePartOfEnergyRange( w, energy, energy );
}


void RelActAutoGui::handleConvertEnergyRangeToIndividuals( Wt::WWidget *w )
{
  RelActAutoEnergyRange *energy_range = dynamic_cast<RelActAutoEnergyRange *>(w);
  assert( energy_range );
  
  const shared_ptr<const RelActCalcAuto::RelActAutoSolution> solution = m_solution;
  if( !solution || (solution->m_status != RelActCalcAuto::RelActAutoSolution::Status::Success) )
  {
    // TODO: just hide/disable the button untill we have a valid solution
    SimpleDialog *dialog = new SimpleDialog( "Can't perform this action.",
                                            "Sorry, a valid solution is needed before an energy range can be split." );
    dialog->addButton( "Continue" );
    
    return;
  }//if( !solution || (solution->m_status != RelActCalcAuto::RelActAutoSolution::Status::Success) )
  
  if( !energy_range || energy_range->isEmpty() )
  {
    SimpleDialog *dialog = new SimpleDialog( "Can't perform this action.",
                                            "Sorry, energy range is currently not valid." );
    dialog->addButton( "Continue" );
    return;
  }
  
  const float lower_energy = energy_range->lowerEnergy();
  const float upper_energy = energy_range->upperEnergy();
  
  vector<RelActCalcAuto::RoiRange> to_ranges;
  for( const RelActCalcAuto::RoiRange &range : solution->m_final_roi_ranges )
  {
    // If the center of `range` falls between `lower_energy` and `upper_energy`, we'll
    //  assume its a match.  This is strictly true, as `range.allow_expand_for_peak_width`
    //  could be true, and/or another ROI can slightly overlap the original one we are
    //  interested in.
    //  TODO: improve the robustness of the matching between the initial ROI, and auto-split ROIs
    
    const double mid_energy = 0.5*(range.lower_energy + range.upper_energy);
    //range.continuum_type = PeakContinuum::OffsetType::;
    //range.force_full_range = false;
    //range.allow_expand_for_peak_width = false;
    
    if( (mid_energy > lower_energy) && (mid_energy < upper_energy) )
      to_ranges.push_back( range );
  }//for( loop over m_final_roi_ranges )
  
  
  // We'll sort the ranges into reverse energy order so when we insert them at a fixed index,
  //  they will be in increasing energy order.
  std::sort( begin(to_ranges), end(to_ranges), []( const RelActCalcAuto::RoiRange &lhs, const RelActCalcAuto::RoiRange &rhs ) -> bool {
    return (lhs.lower_energy + lhs.upper_energy) > (rhs.lower_energy + rhs.upper_energy);
  } );
  
  
  char buffer[512] = { '\0' };
  if( to_ranges.empty() )
  {
    snprintf( buffer, sizeof(buffer),
             "The energy range %.1f keV to %.1f keV did not contain any significant gamma contributions.",
             lower_energy, upper_energy);
    SimpleDialog *dialog = new SimpleDialog( "Can't perform this action.", buffer );
    dialog->addButton( "Continue" );
    return;
  }//if( to_ranges.empty() )
  

  
  snprintf( buffer, sizeof(buffer),
           "This action will divide the energy range from %.1f keV to %.1f keV into %i seperate energy ranges.<br />"
           "<p>Would you like to continue?</p>",
           lower_energy, upper_energy, static_cast<int>(to_ranges.size()) );
  
  SimpleDialog *dialog = new SimpleDialog( "Divide Energy Range Up?", buffer );
  WPushButton *yes_button = dialog->addButton( "Yes" );
  dialog->addButton( "No" );
  
  
  const auto on_yes = [this,w,to_ranges](){
    
    const std::vector<WWidget *> &kids = m_energy_ranges->children();
    const auto pos = std::find( begin(kids), end(kids), w );
    if( pos == end(kids) )
    {
      SimpleDialog *dialog = new SimpleDialog( "Error", "There was an unexpected error finding the original"
                                              " energy range - sorry, cant complete operation." );
      dialog->addButton( "Continue" );
      return;
    }//
    
    const int orig_w_index = static_cast<int>( pos - begin(kids) );
    
    delete w;
    
    for( const RelActCalcAuto::RoiRange &range : to_ranges )
    {
      RelActAutoEnergyRange *roi = new RelActAutoEnergyRange( this );
      roi->updated().connect( this, &RelActAutoGui::handleEnergyRangeChange );
      roi->remove().connect( boost::bind( &RelActAutoGui::handleRemoveEnergy, this, static_cast<WWidget *>(roi) ) );
      
      roi->setFromRoiRange( range );
      roi->setForceFullRange( true );
      
      m_energy_ranges->insertWidget( orig_w_index, roi );
    }//for( const RelActCalcAuto::RoiRange &range : to_ranges )
    
    checkIfInUserConfigOrCreateOne( false );
    m_render_flags |= RenderActions::UpdateEnergyRanges;
    m_render_flags |= RenderActions::UpdateCalculations;
    scheduleRender();
  };//on_yes lamda
  
  
  yes_button->clicked().connect( std::bind(on_yes) );
}//void RelActAutoGui::handleConvertEnergyRangeToIndividuals( Wt::WWidget *w )


void RelActAutoGui::handleAddFreePeak( const double energy, const bool constrain_fwhm, const bool apply_cal )
{
  if( m_free_peaks_container->isHidden() )
    handleShowFreePeaks();
  
  RelActFreePeak *peak = new RelActFreePeak( this, m_free_peaks );
  peak->updated().connect( this, &RelActAutoGui::handleFreePeakChange );
  peak->remove().connect( boost::bind( &RelActAutoGui::handleRemoveFreePeak, this, static_cast<WWidget *>(peak) ) );
  peak->setEnergy( static_cast<float>(energy) );
  peak->setFwhmConstrained( constrain_fwhm );
  peak->setApplyEnergyCal( apply_cal );
  
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= UpdateFreePeaks;
  scheduleRender();
}//void handleAddFreePeak( const double energy, const bool constrain_fwhm )


void RelActAutoGui::handleRemovePartOfEnergyRange( Wt::WWidget *w,
                                                  double lower_energy,
                                                  double upper_energy )
{
  if( !w )
    return;
 
  if( upper_energy < lower_energy )
    std::swap( lower_energy, upper_energy );
  
  const std::vector<WWidget *> &kids = m_energy_ranges->children();
  const auto pos = std::find( begin(kids), end(kids), w );
  if( pos == end(kids) )
  {
    cerr << "Failed to find a RelActAutoEnergyRange in m_energy_ranges (handleRemovePartOfEnergyRange)!" << endl;
    assert( 0 );
    return;
  }
  
  const int orig_w_index = pos - begin( kids );
  RelActAutoEnergyRange *range = dynamic_cast<RelActAutoEnergyRange *>( w );
  assert( range );
  if( !range )
    return;
  
  RelActCalcAuto::RoiRange roi = range->toRoiRange();
  
  if( (upper_energy < roi.lower_energy) || (lower_energy > roi.upper_energy) )
  {
    assert( 0 );
    return;
  }
  
  delete w;
  
  // Check if we want the whole energy range removed
  if( (lower_energy <= roi.lower_energy) && (upper_energy >= roi.upper_energy) )
  {
    // TODO: remove peaks from ROI from PeakModel
    handleEnergyRangeChange();
    return;
  }
  
  const bool is_in_middle = ((upper_energy < roi.upper_energy) && (lower_energy > roi.lower_energy));
  const bool is_left = ((upper_energy > roi.lower_energy) && (lower_energy <= roi.lower_energy));
  const bool is_right = ((lower_energy > roi.lower_energy) && (upper_energy >= roi.upper_energy));
  
  assert( is_in_middle || is_left || is_right );
  assert( (int(is_in_middle) + int(is_left) + int(is_right)) == 1 );
  
  if( is_in_middle )
  {
    RelActAutoEnergyRange *left_range = new RelActAutoEnergyRange( this );
    left_range->updated().connect( this, &RelActAutoGui::handleEnergyRangeChange );
    left_range->remove().connect( boost::bind( &RelActAutoGui::handleRemoveEnergy,
                                              this, static_cast<WWidget *>(left_range) ) );
    
    RelActAutoEnergyRange *right_range = new RelActAutoEnergyRange( this );
    right_range->updated().connect( this, &RelActAutoGui::handleEnergyRangeChange );
    right_range->remove().connect( boost::bind( &RelActAutoGui::handleRemoveEnergy,
                                               this, static_cast<WWidget *>(right_range) ) );
    
    left_range->setFromRoiRange( roi );
    right_range->setFromRoiRange( roi );
    left_range->setEnergyRange( roi.lower_energy, lower_energy );
    right_range->setEnergyRange( upper_energy, roi.upper_energy );
    
    m_energy_ranges->insertWidget( orig_w_index, right_range );
    m_energy_ranges->insertWidget( orig_w_index, left_range );
    
    // TODO: we could split PeakModels ROI peaks here and set them to provide instant feedback during computation
  }else if( is_left )
  {
    RelActAutoEnergyRange *right_range = new RelActAutoEnergyRange( this );
    right_range->updated().connect( this, &RelActAutoGui::handleEnergyRangeChange );
    right_range->remove().connect( boost::bind( &RelActAutoGui::handleRemoveEnergy,
                                               this, static_cast<WWidget *>(right_range) ) );
    right_range->setEnergyRange( upper_energy, roi.upper_energy );
    right_range->setFromRoiRange( roi );
    m_energy_ranges->insertWidget( orig_w_index, right_range );
    
    // TODO: we could update PeakModels peaks/range here and set them to provide instant feedback during computation
  }else if( is_right )
  {
    RelActAutoEnergyRange *left_range = new RelActAutoEnergyRange( this );
    left_range->updated().connect( this, &RelActAutoGui::handleEnergyRangeChange );
    left_range->remove().connect( boost::bind( &RelActAutoGui::handleRemoveEnergy,
                                              this, static_cast<WWidget *>(left_range) ) );
    left_range->setFromRoiRange( roi );
    left_range->setEnergyRange( roi.lower_energy, lower_energy );
    m_energy_ranges->insertWidget( orig_w_index, left_range );
    
    // TODO: we could update PeakModels peaks/range here and set them to provide instant feedback during computation
  }
  
  handleEnergyRangeChange();
}//void handleSplitEnergyRange( Wt::WWidget *energy_range, const double energy )


void RelActAutoGui::handleRemoveNuclide( Wt::WWidget *w )
{
  if( !w )
    return;
  
  assert( dynamic_cast<RelActAutoNuclide *>(w) );
  
  const std::vector<WWidget *> &kids = m_nuclides->children();
  const auto pos = std::find( begin(kids), end(kids), w );
  if( pos == end(kids) )
  {
    cerr << "Failed to find a RelActAutoNuclide in m_nuclides!" << endl;
    assert( 0 );
    return;
  }
  
  delete w;
  
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateNuclidesPresent;
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleRemoveNuclide( Wt::WWidget *w )


void RelActAutoGui::startApplyFitEnergyCalToSpecFile()
{
  const bool fit_offset = (m_solution && m_solution->m_fit_energy_cal[0]);
  const bool fit_gain = (m_solution && m_solution->m_fit_energy_cal[1]);
  if( !fit_offset && !fit_gain )
    return;
  
  string msg = "This will";
  
  char buffer[128] = { '\0' };
  
  bool printed_some = false;
  if( fit_offset )
  {
    double offset = -(m_solution->m_energy_cal_adjustments[0]/RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0)
                      * RelActCalcAuto::RelActAutoSolution::sm_energy_offset_range_keV;
    snprintf( buffer, sizeof(buffer), " add an offset of %.2f keV", offset );
    msg += buffer;
    printed_some = true;
  }
  
  if( fit_gain )
  {
    double gain = -(m_solution->m_energy_cal_adjustments[1]/RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0)
                      * RelActCalcAuto::RelActAutoSolution::sm_energy_gain_range_keV;
    snprintf( buffer, sizeof(buffer), "%s increase gain by %.5f", (printed_some ? " and" : ""), gain );
    msg += buffer;
    printed_some = true;
  }//if( fit_gain )
  
  if( m_solution && (m_solution->m_fit_energy_cal.size() > 2) && m_solution->m_fit_energy_cal[2] )
  {
    double quad = -(m_solution->m_energy_cal_adjustments[2]/RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0)
                      * RelActCalcAuto::RelActAutoSolution::sm_energy_quad_range_keV;
    snprintf( buffer, sizeof(buffer), "%s increase quadratic by %.5f", (printed_some ? " and" : ""), quad );
    msg += buffer;
    printed_some = true;
  }//if( fit_gain )
  
  
  msg += " for the primary foreground";
  const bool has_back = !!m_interspec->displayedHistogram(SpecUtils::SpectrumType::Background);
  if( has_back )
    msg += " and background";
  msg += has_back ? " files" : " file";
  msg += ".<br />Would you like to do this?";
  
  
  SimpleDialog *dialog = new SimpleDialog( "Apply fit energy calibration?", msg );
  WPushButton *yes = dialog->addButton( "Yes" );
  dialog->addButton( "No" );
  yes->clicked().connect( this, &RelActAutoGui::applyFitEnergyCalToSpecFile );
}//void startApplyFitEnergyCalToSpecFile();


void RelActAutoGui::applyFitEnergyCalToSpecFile()
{
  // We will apply to currently displayed spectra; its too complicated to
  //  give the user all the options of what to apply it to, like the
  //  energy calibration tool
  bool ownEnergyCal = false;
  EnergyCalTool *tool = nullptr;
  try
  {
    assert( m_solution->m_foreground );
    if( !m_solution || !m_solution->m_foreground )
      throw runtime_error( "Solution foreground not set???" );
    
    const auto orig_cal = m_solution->m_foreground->energy_calibration();
    assert( orig_cal && orig_cal->valid() );
    if( !orig_cal || !orig_cal->valid() )
      throw runtime_error( "Solution foreground energy calibration invalid???" );
    
    const size_t nchannel = orig_cal->num_channels();
    
    shared_ptr<const SpecUtils::EnergyCalibration> new_cal = m_solution->m_spectrum
                                      ? m_solution->m_spectrum->energy_calibration() : nullptr; // or m_solution->get_adjusted_energy_cal()
    
    assert( new_cal && new_cal->valid() );
    if( !new_cal || !new_cal->valid() ) //shouldnt ever happen
      throw runtime_error( "Updated energy calibration is invalid - not applying" );
    
    tool = m_interspec->energyCalTool();
    if( !tool )
    {
      ownEnergyCal = true;
      tool = new EnergyCalTool( m_interspec, m_interspec->peakModel(), nullptr );
    }
    
    MeasToApplyCoefChangeTo fore, back;
    fore.meas = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
    fore.sample_numbers = m_interspec->displayedSamples(SpecUtils::SpectrumType::Foreground);
    const auto fore_dets = m_interspec->detectorsToDisplay(SpecUtils::SpectrumType::Foreground);
    fore.detectors = set<string>( begin(fore_dets), end(fore_dets) );
    
    back.meas = m_interspec->measurment(SpecUtils::SpectrumType::Background);
    back.sample_numbers = m_interspec->displayedSamples(SpecUtils::SpectrumType::Background);
    const auto back_dets = m_interspec->detectorsToDisplay(SpecUtils::SpectrumType::Background);
    back.detectors = set<string>( begin(back_dets), end(back_dets) );
    
    // Should we apply it to the secondary spectrum?
    
    vector<MeasToApplyCoefChangeTo> change_meas;
    if( fore.meas )
      change_meas.push_back( fore );
    if( back.meas )
      change_meas.push_back( back );
    
    if( change_meas.empty() ) // never expect to happe
      throw runtime_error( "Somehow invalid foreground or background SpecMeas?" );
    
    tool->applyCalChange( orig_cal, new_cal, change_meas, false );
    
    string msg = "Have updated energy calibration for displayed foreground";
    if( back.meas )
      msg += " and background.";
    
    passMessage( msg, WarningWidget::WarningMsgInfo );
    m_fit_energy_cal->setChecked( false );
  }catch( std::exception &e )
  {
    passMessage( "Error applying energy calibration: " + string(e.what()), WarningWidget::WarningMsgHigh );
  }// try / catch
  
  if( ownEnergyCal && tool )
  {
    delete tool;
    tool = nullptr;
  }
  
  m_render_flags |= RenderActions::UpdateSpectra;
  m_render_flags |= RenderActions::UpdateCalculations;
  m_render_flags |= RenderActions::ChartToDefaultRange;
  scheduleRender();
}//void applyFitEnergyCalToSpecFile()


void RelActAutoGui::handleShowRefLines( const bool show )
{
  if( show == m_hide_ref_lines_item->isEnabled() )
    return;
  
  m_show_ref_lines_item->setHidden( show );
  m_hide_ref_lines_item->setHidden( !show );
  
  // We will use the enabled/disabled state to track which option
  //  the user wants.  The isHidden() state is affected by the popup
  //  menu status, so we cant use it.
  m_show_ref_lines_item->setDisabled( show );
  m_hide_ref_lines_item->setDisabled( !show );
  
  m_render_flags |= RenderActions::UpdateRefGammaLines;
  scheduleRender();
}//void RelActAutoGui::handleShowRefLines()


void RelActAutoGui::setPeaksToForeground()
{
  assert( m_solution && !m_solution->m_fit_peaks_in_spectrums_cal.empty() );
  if( !m_solution || m_solution->m_fit_peaks_in_spectrums_cal.empty() )
  {
    SimpleDialog *dialog = new SimpleDialog( "Can't Continue", "No peaks in current solution" );
    dialog->addButton( "Close" );
    return;
  }//if( no solution peaks )
  
  PeakModel *peak_model = m_interspec->peakModel();
  assert( peak_model );
  if( !peak_model )
    return;
  
  SimpleDialog *dialog = new SimpleDialog( "Use peaks with foreground?", "" );
  dialog->addStyleClass( "SetToPeaksDialog" );
  WText *message = new WText( "Peaks will not have uncertainties unless you choose"
                             " to re-fit them, in which case the re-fit peaks may"
                             " differ from the current peaks since they will no"
                             " longer be constrained by the Relative Efficiency"
                             " curve or spectrum wide FWHM response.", dialog->contents() );
  message->addStyleClass( "content" );
  message->setInline( false );
  
  SwitchCheckbox *replace_or_add = nullptr;
  const vector<PeakDef> previous_peaks = peak_model->peakVec();
  if( !previous_peaks.empty() )
  {
    WContainerWidget *holder = new WContainerWidget( dialog->contents() );
    holder->addStyleClass( "AddOrReplaceSwitchRow" );
    
    replace_or_add = new SwitchCheckbox( "Add peaks", "Replace peaks", holder );
    replace_or_add->setChecked( true ); //Make "Replace peaks" the default answer
  }//if( we have peaks )
  
  WContainerWidget *refit_holder = new WContainerWidget( dialog->contents() );
  refit_holder->addStyleClass( "AddOrReplaceRefitRow" );
  WCheckBox *refit_peaks = new WCheckBox( "Refit Peaks", refit_holder );
  
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );
  const char *tooltip = 
  "When checked, peaks will be refit without the constraints of the relative efficiency curve,"
  " expected branching ratios, or FWHM constraints from other ROIs - allowing statistical"
  " uncertainties to be assigned to the peaks amplitude.<br/>"
  "Peaks that are within 0.53 FWHM (1.25 sigma) of each other will be merged together.<br/>"
  "  Peak mean, FWHM, and continuums will be refit, with usually peak means limited to be changed"
  " by no more than 0.11 times the peak FWHM, but if this fails,"
  " then the limit may be increased to 0.21 times the peak FWHM. <br/>"
  "Fit peak amplitudes may also be limits in how much can be changed from the relative efficiency"
  " peak amplitude, so you may need to manually refit some peaks again.<br/>";
  HelpSystem::attachToolTipOn( refit_holder, tooltip, showToolTips );
  
  
  dialog->addButton( "No" );
  WPushButton *yes = dialog->addButton( "Yes" );
  
  
  const vector<PeakDef> solution_peaks = m_solution->m_fit_peaks_in_spectrums_cal;
  std::shared_ptr<const DetectorPeakResponse> ana_drf = m_solution->m_drf;
  
  if( m_solution->m_options.fit_energy_cal )
  {
    // The fit peaks have already been adjusted for energy calibration, so I dont
    //  think we need to update them here
  }//if( m_solution->m_options.fit_energy_cal )
  
  yes->clicked().connect( std::bind([solution_peaks, replace_or_add, refit_peaks, previous_peaks, ana_drf](){
    const bool replace_peaks = (!replace_or_add || replace_or_add->isChecked());
    
    InterSpec *interpsec = InterSpec::instance();
    assert( interpsec );
    if( !interpsec )
      return;
    
    shared_ptr<const SpecUtils::Measurement> foreground = interpsec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    
    vector<PeakDef> final_peaks;
    if( foreground && refit_peaks->isChecked() )
    {
      map< shared_ptr<const PeakContinuum>,vector<shared_ptr<const PeakDef>> > rois;
     
      for( const PeakDef &peak : solution_peaks )
        rois[peak.continuum()].push_back( make_shared<PeakDef>(peak) );
      
      vector< vector<shared_ptr<const PeakDef>> > fit_peaks( rois.size() );
      
      SpecUtilsAsync::ThreadPool pool;
      size_t roi_num = 0;
      for( const auto &cont_peaks : rois )
      {
        const vector<shared_ptr<const PeakDef>> *peaks = &(cont_peaks.second);
        
        pool.post( [&fit_peaks, roi_num, foreground, peaks, ana_drf](){
          
          // If two peaks are near each other, we wont be able to resolve them in the fit,
          //  so just get rid of the smaller amplitude peak
          vector<shared_ptr<const PeakDef>> peaks_to_filter = *peaks;
          std::sort( begin(peaks_to_filter), end(peaks_to_filter),
                []( const shared_ptr<const PeakDef> &lhs, const shared_ptr<const PeakDef> &rhs ) -> bool {
            if( (lhs->type() != PeakDef::GaussianDefined) 
               || (rhs->type() != PeakDef::GaussianDefined) )
            {
              return (lhs->type() < rhs->type());
            }
            return lhs->amplitude() > rhs->amplitude();
          } ); //sort(...)
          
          
          vector<shared_ptr<const PeakDef>> peaks_to_refit;
          for( const auto &to_add : peaks_to_filter )
          {
            if( to_add->type() != PeakDef::GaussianDefined )
            {
              peaks_to_refit.push_back( to_add );
              continue;
            }
            
            bool keep = true;
            for( const auto &already_added : peaks_to_refit )
            {
              if( already_added->type() != PeakDef::GaussianDefined )
                continue;
              
              // Using the default value of ShieldingSourceFitOptions::photopeak_cluster_sigma,
              //  1.25, to decide if we should keep this peak or not
              if( fabs(to_add->mean() - already_added->mean()) < 1.25*already_added->sigma() )
              {
                keep = false;
                break;
              }
            }//for( const auto &already_added : peaks_to_refit )
            
            if( keep )
              peaks_to_refit.push_back( to_add );
          }//for( const auto &to_add : peaks_to_filter )
          
          std::sort( begin(peaks_to_refit), end(peaks_to_refit), &PeakDef::lessThanByMeanShrdPtr );
          
          
          const double meanSigmaVary = 0.25; //arbitrary
          fit_peaks[roi_num] = refitPeaksThatShareROI( foreground, ana_drf, peaks_to_refit, meanSigmaVary );
          
          if( fit_peaks[roi_num].size() != peaks_to_refit.size() )
          {
            cout << "refitPeaksThatShareROI gave " << fit_peaks[roi_num].size() << " peaks, while"
            << " we wanted " << peaks->size() << ", will try fitPeaksInRange(...)" << endl;
            vector<PeakDef> input_peaks, fixed_peaks;
            for( const auto &p : peaks_to_refit )
              input_peaks.push_back( *p );
            
            const double lx = input_peaks.front().lowerX();
            const double ux = input_peaks.front().upperX();
            
            const double ncausality = 10;
            const double stat_threshold = 0.5;
            const double hypothesis_threshold = -1.0;
            
            const vector<PeakDef> retry_peak = fitPeaksInRange( lx, ux, ncausality, stat_threshold,
                                                          hypothesis_threshold, input_peaks,
                                                          foreground, fixed_peaks, true );
            
            if( (retry_peak.size() == peaks_to_refit.size())
               || (fit_peaks[roi_num].empty() && !retry_peak.empty()) )
            {
              // This is *usually* the case
              fit_peaks[roi_num].clear();
              for( const auto &p : retry_peak )
                fit_peaks[roi_num].push_back( make_shared<PeakDef>(p) );
            }else if( !fit_peaks[roi_num].empty() )
            {
              // Maybe a peak became insignificant or something, just go with it
            }else
            {
              cerr << "fitPeaksInRange(...) also failed us, giving " << retry_peak.size()
              << " peaks when we wanted " << peaks_to_refit.size()
              << ", will just use Rel. Eff. peaks." << endl;
              fit_peaks[roi_num] = *peaks;
            }//if( retry_peak.size() == peaks->size() ) / else
          }//if( fit_peaks[roi_num].size() != peaks->size() )
        } );
        
        ++roi_num;
      }//for( auto &cont_peaks : rois )

      pool.join();
      
      for( const auto &pvec : fit_peaks )
      {
        for( const auto &p : pvec )
          final_peaks.push_back( *p );
      }
      std::sort( begin(final_peaks), end(final_peaks), &PeakDef::lessThanByMean );
    }else
    {
      final_peaks = solution_peaks;
    }
    
    PeakModel *peak_model = interpsec->peakModel();
    assert( peak_model );
    if( !peak_model )
      return;
    
    if( replace_peaks )
      peak_model->setPeaks( final_peaks );
    else
      peak_model->addPeaks( final_peaks );
    
    // Potential TODO: we could do something more sophisticated when adding peaks - kinda merge
    //  them like we do when we do a peak search, see the `PeakSelectorWindow` class constructor
    //  for how its done there
  }) );
  
}//void setPeaksToForeground()


void RelActAutoGui::handleRelEffModelOptionsChanged()
{
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleRelEffModelOptionsChanged()



void RelActAutoGui::handleDetectorChange()
{
  // We could sometimes get away with not updating calculations if we arent using a physical model,
  //  but I think having the detectors FWHM may slightly impact the auto-search peaks (not 100% sure),
  //  so we'll just always refresh calculations on detector changes.

  m_cached_drf = nullptr;
  m_cached_all_peaks.clear();

  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void RelActAutoGui::handleDetectorChange()


void RelActAutoGui::handleAddRelEffCurve()
{ 
  set<string> existing_curve_names;
  for( int index = 0; index < m_rel_eff_opts_menu->count(); ++index )
  {
    WMenuItem *this_item = m_rel_eff_opts_menu->itemAt( index );
    existing_curve_names.insert( this_item->text().toUTF8() );
  }

  string name;
  for( int index = 0; index < 10; ++index )
  {
    name = "Curve " + std::to_string(index);
    if( existing_curve_names.count(name) == 0 )
      break;
  }

  RelActAutoGuiRelEffOptions *rel_eff_curve = new RelActAutoGuiRelEffOptions( this, name, nullptr );  
  WMenuItem *item = new WMenuItem( name, rel_eff_curve, WMenuItem::LoadPolicy::PreLoading );
  m_rel_eff_opts_menu->addItem( item );
  assert( item->contents() == static_cast<WWidget *>(rel_eff_curve) );
  
  // When outside the link area is clicked, the item doesnt get selected, so we'll work around this.
  item->clicked().connect( std::bind([this,item](){
    m_rel_eff_opts_menu->select( item );
    item->triggered().emit( item );
  }) );

  rel_eff_curve->addRelEffCurve().connect( this, &RelActAutoGui::handleAddRelEffCurve );
  rel_eff_curve->delRelEffCurve().connect( boost::bind( &RelActAutoGui::handleDelRelEffCurve, this, boost::placeholders::_1 ) );
  rel_eff_curve->nameChanged().connect( this, &RelActAutoGui::handleRelEffCurveNameChanged );
  rel_eff_curve->optionsChanged().connect( this, &RelActAutoGui::handleRelEffModelOptionsChanged );

  const bool single_curve = (m_rel_eff_opts_menu->count() == 1);

  m_rel_eff_opts_menu->setHidden( single_curve );

  for( int index = 0; index < m_rel_eff_opts_menu->count(); ++index )
  {
    WMenuItem *this_item = m_rel_eff_opts_menu->itemAt( index );
    RelActAutoGuiRelEffOptions *this_curve = dynamic_cast<RelActAutoGuiRelEffOptions *>( this_item->contents() );
    assert( this_curve );
    if( this_curve )
      this_curve->setIsOnlyOneRelEffCurve( single_curve );
  }

  m_rel_eff_opts_menu->select( item );

  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleAddRelEffCurve( RelActAutoGuiRelEffOptions *curve )


void RelActAutoGui::handleDelRelEffCurve( RelActAutoGuiRelEffOptions *curve )
{
  assert( curve );
  if( !curve )
    return;

  assert( m_rel_eff_opts_menu->count() > 1 ); 
  if( m_rel_eff_opts_menu->count() <= 1 )
    return;

  WMenuItem *item = nullptr;
  for( int index = 0; !item && (index < m_rel_eff_opts_menu->count()); ++index )
  {
    WMenuItem *this_item = m_rel_eff_opts_menu->itemAt( index );
    if( this_item->contents() == curve )
      item = this_item;
  }
  
  assert( item );
  if( !item )
    return;

  m_rel_eff_opts_menu->removeItem( item );
  delete curve;
  delete item;

  
  const bool single_curve = (m_rel_eff_opts_menu->count() == 1);
  m_rel_eff_opts_menu->setHidden( single_curve );
  for( int index = 0; index < m_rel_eff_opts_menu->count(); ++index )
  {
    WMenuItem *this_item = m_rel_eff_opts_menu->itemAt( index );
    RelActAutoGuiRelEffOptions *this_curve = dynamic_cast<RelActAutoGuiRelEffOptions *>( this_item->contents() );
    assert( this_curve );
    if( this_curve )
      this_curve->setIsOnlyOneRelEffCurve( single_curve );
  }

  m_rel_eff_opts_menu->select( -1 );
  m_rel_eff_opts_menu->select( m_rel_eff_opts_menu->count() - 1 );
  
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleDelRelEffCurve( RelActAutoGuiRelEffOptions *curve )


void RelActAutoGui::handleRelEffCurveNameChanged( RelActAutoGuiRelEffOptions *curve, const Wt::WString &name )
{
  assert( curve );
  if( !curve )
    return;
  
  WMenuItem *item = nullptr;
  for( int index = 0; !item && (index < m_rel_eff_opts_menu->count()); ++index )
  {
    WMenuItem *this_item = m_rel_eff_opts_menu->itemAt( index );
    if( this_item->contents() == curve )
      item = this_item;
  }
  
  assert( item );
  if( !item )
    return;

  item->setText( name );

  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handleRelEffCurveNameChanged( RelActAutoGuiRelEffOptions *curve, const Wt::WString &name )


Wt::Signal<> &RelActAutoGui::calculationStarted()
{
  return m_calc_started;
}

Wt::Signal<> &RelActAutoGui::calculationSuccessful()
{
  return m_calc_successful;
}


Wt::Signal<> &RelActAutoGui::calculationFailed()
{
  return m_calc_failed;
}


Signal<shared_ptr<const RelActCalcAuto::RelActAutoSolution> > &RelActAutoGui::solutionUpdated()
{
  return m_solution_updated;
}


void RelActAutoGui::addDownloadAndUploadLinks( Wt::WContainerWidget *parent )
{
  if( !parent )
    return;
  
#if( BUILD_AS_OSX_APP || IOS )
  WAnchor *btn = new WAnchor( WLink(m_html_download_rsc), parent );
  btn->setTarget( AnchorTarget::TargetNewWindow );
  btn->setStyleClass( "LinkBtn DownloadLink RelActDownload" );
#else
  WPushButton *btn = new WPushButton( parent );
  btn->setIcon( "InterSpec_resources/images/download_small.svg" );
  btn->setLink( WLink( m_html_download_rsc ) );
  btn->setLinkTarget( Wt::TargetNewWindow );
  btn->setStyleClass( "LinkBtn DownloadBtn RelActDownload" );

#if( ANDROID )
  // Using hacked saving to temporary file in Android, instead of via network download of file.
  btn->clicked().connect( std::bind([this](){
    android_download_workaround( m_html_download_rsc, "isotopics_by_nuclide.html");
  }) );
#endif //ANDROID
#endif

  btn->setText( "HTML Report" );
  
  m_calc_started.connect( btn, &WWidget::disable );
  m_calc_failed.connect( btn, &WWidget::disable );
  m_calc_successful.connect( btn, &WWidget::enable );
  
#if( BUILD_AS_OSX_APP || IOS )
  btn = new WAnchor( WLink(m_xml_download_rsc), parent );
  btn->setTarget( AnchorTarget::TargetNewWindow );
  btn->setStyleClass( "LinkBtn DownloadLink RelActDownload" );
  btn->setText( "XML Config" );
#else
  btn = new WPushButton( "XML Config", parent );
  btn->setIcon( "InterSpec_resources/images/download_small.svg" );
  btn->setLinkTarget( Wt::TargetNewWindow );
  btn->setStyleClass( "LinkBtn DownloadBtn RelActDownload" );
  btn->setLink( WLink(m_xml_download_rsc) );
  
#if( ANDROID )
  // Using hacked saving to temporary file in Android, instead of via network download of file.
  btn->clicked().connect( std::bind([this](){
    android_download_workaround( m_xml_download_rsc, "isotopics_by_nuclide_config.html");
  }) );
#endif //ANDROID
  
#endif

  // TODO: add XML upload...
  WPushButton *uploadbtn = new WPushButton( parent );
  uploadbtn->setIcon( "InterSpec_resources/images/upload_small.svg" );
  uploadbtn->setStyleClass( "LinkBtn UploadBtn RelActDownload" );
  uploadbtn->clicked().connect( this, &RelActAutoGui::handleRequestToUploadXmlConfig );
}//void addDownloadAndUploadLinks( Wt::WContainerWidet *parent )


void RelActAutoGui::handleRequestToUploadXmlConfig()
{
  SimpleDialog *dialog = new SimpleDialog();
  WPushButton *closeButton = dialog->addButton( "Cancel" );
  WGridLayout *stretcher = new WGridLayout();
  stretcher->setContentsMargins( 0, 0, 0, 0 );
  dialog->contents()->setLayout( stretcher );
  dialog->contents()->setOverflow( WContainerWidget::Overflow::OverflowVisible,
                                  Wt::Horizontal | Wt::Vertical );
  WText *title = new WText( "Import XML config file" );
  title->addStyleClass( "title" );
  stretcher->addWidget( title, 0, 0 );
  
  WText *t = new WText( "<p>Select the <em>Isotopics by nuclide</em> XML file to use</p>" );
  stretcher->addWidget( t, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
  t->setTextAlignment( Wt::AlignCenter );
  
  
  WFileUpload *upload = new WFileUpload();
  upload->fileTooLarge().connect( std::bind( [=](){
    dialog->contents()->clear();
    dialog->footer()->clear();
    
    WPushButton *closeButton = dialog->addButton( "Close" );
    WGridLayout *stretcher = new WGridLayout();
    stretcher->setContentsMargins( 0, 0, 0, 0 );
    dialog->contents()->setLayout( stretcher );
    WText *title = new WText( "File to large to upload" );
    title->addStyleClass( "title" );
    stretcher->addWidget( title, 0, 0 );
  }) );
  
  upload->changed().connect( upload, &WFileUpload::upload );
  upload->uploaded().connect( std::bind( [this,dialog,upload](){
    
    try
    {
      const string xml_path = upload->spoolFileName();
      rapidxml::file<char> input_file( xml_path.c_str() );;
      rapidxml::xml_document<char> doc;
      doc.parse<rapidxml::parse_trim_whitespace>( input_file.data() );
      
      setGuiStateFromXml( &doc );
    }catch( rapidxml::parse_error &e )
    {
      string msg = "Error parsing config XML: " + string(e.what());
      const char * const position = e.where<char>();
      if( position && *position )
      {
        const char *end_pos = position;
        for( size_t i = 0; (*end_pos) && (i < 80); ++i )
          end_pos += 1;
        msg += "<br />&nbsp;&nbsp;At: " + std::string(position, end_pos);
      }//if( position )
      
      passMessage( msg, WarningWidget::WarningMsgHigh );
    }catch( std::exception &e )
    {
      passMessage( "Error loading <em>Isotopics by nuclide</em> XML config file: "
                  + string(e.what()), WarningWidget::WarningMsgHigh );
    }//try / cat to read the XML
    
    dialog->accept();
    
    //wApp->doJavaScript( "$('.Wt-dialogcover').hide();" ); // JIC
    //dialog->done( Wt::WDialog::DialogCode::Accepted );
  } ) );
  
  stretcher->addWidget( upload, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
  
  InterSpec *interspec = InterSpec::instance();
  if( interspec && !interspec->isPhone() )
  {
    t = new WText( "<p style=\"font-size: small;\">Note: you can also drag-n-drop the XML config files onto InterSpec<br /></p>" );
    stretcher->addWidget( t, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
    t->setTextAlignment( Wt::AlignCenter );
  }
  
  /*
   //In case we want to use AuxWindow instead of SimpleDialog
   AuxWindow *window = new AuxWindow( "Import CALp file",
   (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
   | AuxWindowProperties::PhoneNotFullScreen
   | AuxWindowProperties::DisableCollapse
   | AuxWindowProperties::SetCloseable) );
   
   //...
   
   window->rejectWhenEscapePressed();
   window->show();
   window->resizeToFitOnScreen();
   window->centerWindow();
   
   WPushButton *close = window->addCloseButtonToFooter( "Cancel" );
   close->clicked().connect( boost::bind( &AuxWindow::hide, window ) );
   
   window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
   
   // TODO: add link to relevant section of documentation
   //AuxWindow::addHelpInFooter( window->footer(), "energy-cal-CALp" );
   */
}//void handleRequestToUploadCALp();


RelActAutoGuiRelEffOptions *RelActAutoGui::getRelEffCurveOptions( const int index )
{
  assert( m_rel_eff_opts_menu );
  assert( index >= 0 );
  assert( index < m_rel_eff_opts_menu->count() );
  if( (index < 0) || (index >= m_rel_eff_opts_menu->count()) )
    return nullptr;
    
  Wt::WMenuItem * const item = m_rel_eff_opts_menu->itemAt( index );
  assert( item );
  if( !item )
    return nullptr;
  
  WWidget * const contents = item->contents();
  assert( contents );
  if( !contents )
    return nullptr;
  
  RelActAutoGuiRelEffOptions *ptr = dynamic_cast<RelActAutoGuiRelEffOptions *>( contents );
  if( !ptr )
    throw runtime_error( "Failed to cast to RelActAutoGuiRelEffOptions" );

  return ptr;
}//RelActAutoGuiRelEffOptions *RelActAutoGui::getRelEffCurveOptions( const size_t index )


const RelActAutoGuiRelEffOptions *RelActAutoGui::getRelEffCurveOptions( const int index ) const
{
  assert( m_rel_eff_opts_menu );
  assert( index >= 0 );
  assert( index < m_rel_eff_opts_menu->count() );
  if( (index < 0) || (static_cast<int>(index) >= m_rel_eff_opts_menu->count()) )
    return nullptr;
  
  Wt::WMenuItem * const item = m_rel_eff_opts_menu->itemAt( index );
  assert( item );
  if( !item )
    return nullptr;
  
  const WWidget * const contents = item->contents();
  assert( contents );
  if( !contents )
    return nullptr;
  
  const RelActAutoGuiRelEffOptions * const ptr = dynamic_cast<const RelActAutoGuiRelEffOptions *>( contents );
  if( !ptr )
    throw runtime_error( "Failed to cast to RelActAutoGuiRelEffOptions" );

  return ptr;
}//RelActAutoGuiRelEffOptions *RelActAutoGui::getRelEffCurveOptions( const size_t index )


void RelActAutoGui::updateDuringRenderForSpectrumChange()
{
  const shared_ptr<const SpecUtils::Measurement> fore
                   = m_interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  const shared_ptr<const SpecUtils::Measurement> back
                   = m_interspec->displayedHistogram( SpecUtils::SpectrumType::Background );
  
  bool foreground_same = (fore == m_foreground);
  
  // If we are creating this widget, and deSerializing, we have already set the energy range we
  //  want to display, so dont overwrite that, if it looks like it has been set.
  if( !foreground_same
     && !m_spectrum->data()
     && fore
     && !m_spectrum->isRendered()
     && (m_spectrum->xAxisMaximum() > (5.0 + m_spectrum->xAxisMinimum()))
     && (m_spectrum->xAxisMinimum() >= fore->gamma_energy_min())
     && (m_spectrum->xAxisMaximum() <= fore->gamma_energy_max())
     )
  {
    foreground_same = true;
  }
  
  m_foreground = fore;
  m_background = back;
  if( m_background )
    m_background_sf = m_interspec->displayScaleFactor( SpecUtils::SpectrumType::Background );
  else
    m_background_sf = 1.0;
  
  m_spectrum->setData( m_foreground, foreground_same );
  m_spectrum->setBackground( m_background );
  m_spectrum->setDisplayScaleFactor( m_background_sf, SpecUtils::SpectrumType::Background );
  
  m_spectrum->removeAllDecorativeHighlightRegions();
  //m_spectrum->setChartTitle( WString() );
  m_fit_chi2_msg->setText( "" );
  m_fit_chi2_msg->hide();

  m_background_subtract->setHidden( !m_background );
}//void updateDuringRenderForSpectrumChange()


void RelActAutoGui::updateSpectrumToDefaultEnergyRange()
{
  if( !m_foreground || (m_foreground->gamma_energy_max() < 1.0f) )
  {
    m_spectrum->setXAxisRange( 0, 3000 );
    return;
  }
  
  const double spec_min = m_foreground->gamma_energy_min();
  const double spec_max = m_foreground->gamma_energy_max();
  
  const vector<RelActCalcAuto::RoiRange> rois = getRoiRanges();
  if( rois.empty() )
  {
    m_spectrum->setXAxisRange( spec_min, spec_max );
    return;
  }
  
  double min_energy = 3000, max_energy = 0;
  for( const RelActCalcAuto::RoiRange &roi : rois )
  {
    min_energy = std::min( min_energy, roi.lower_energy );
    max_energy = std::max( max_energy, roi.upper_energy );
  }
  
  //Make sure max energy is greater than min energy by at least ~10 keV; this is arbitrary, but we
  //  dont ant the spectrum hyper-zoomed-in to like a single channel
  if( (max_energy > min_energy) && ((max_energy - min_energy) > 10.0) )
  {
    const double range = max_energy - min_energy;
    min_energy -= 0.1*range;
    max_energy += 0.1*range;
    
    min_energy = std::max( min_energy, spec_min );
    max_energy = std::min( max_energy, spec_max );
    
    m_spectrum->setXAxisRange( min_energy, max_energy );
  }else
  {
    m_spectrum->setXAxisRange( spec_min, spec_max );
  }
}//void updateSpectrumToDefaultEnergyRange()


void RelActAutoGui::updateDuringRenderForNuclideChange()
{ 
  bool has_multiple_nucs_of_z = false;
  map<short,int> z_to_num_isotopes;
  
  const int num_rel_eff_curves = m_rel_eff_opts_menu->count();
  for( int rel_eff_curve_index = 0; rel_eff_curve_index < num_rel_eff_curves; ++rel_eff_curve_index )
  {
    set<string> nuc_names;

    RelActAutoGuiRelEffOptions *opts = getRelEffCurveOptions( rel_eff_curve_index );
    assert( opts );
    if( !opts )
      continue;

    const vector<RelActCalcAuto::NucInputInfo> nuclides = getNucInputInfo( rel_eff_curve_index );
    
    opts->updatePuCorrelationOptions( nuclides );

    for( const RelActCalcAuto::NucInputInfo &nuc : nuclides )
    {
      if( !nuc.nuclide )
        continue;
    
      const short z = nuc.nuclide->atomicNumber;
      if( !z_to_num_isotopes.count(z) )
        z_to_num_isotopes[z] = 0;
    
      int &num_this_z = z_to_num_isotopes[z];
      num_this_z += 1;
    
      has_multiple_nucs_of_z = has_multiple_nucs_of_z || (num_this_z > 1);
    }//for( const RelActCalcAuto::NucInputInfo &nuc : nuclides )
  }//for( int rel_eff_curve_index = 0; rel_eff_curve_index < num_rel_eff_curves; ++rel_eff_curve_index )

  if( m_same_z_age->isVisible() != has_multiple_nucs_of_z )
    m_same_z_age->setHidden( !has_multiple_nucs_of_z );
}//void updateDuringRenderForNuclideChange()


void RelActAutoGui::updateDuringRenderForRefGammaLineChange()
{
  // Determine if we should show/hide the reference lines.
  const bool show_ref_lines = m_hide_ref_lines_item->isEnabled();
  vector<RelActCalcAuto::NucInputInfo> nuclides;

  const int num_rel_eff_curves = m_rel_eff_opts_menu->count();
  for( int rel_eff_curve_index = 0; rel_eff_curve_index < num_rel_eff_curves; ++rel_eff_curve_index )
  {
    const vector<RelActCalcAuto::NucInputInfo> nucs = getNucInputInfo( rel_eff_curve_index );
    for( const RelActCalcAuto::NucInputInfo &nuc : nucs )
    {
      const auto pos = std::find_if( begin(nuclides), end(nuclides),
        [&nuc]( const RelActCalcAuto::NucInputInfo &nuc2 ){
          return (nuc.nuclide == nuc2.nuclide) && (nuc.element == nuc2.element) && (nuc.reaction == nuc2.reaction);
      } );

      if( pos == end(nuclides) )
        nuclides.push_back( nuc );
    }//for( const RelActCalcAuto::NucInputInfo &nuc : nucs )
  }//for( int rel_eff_curve_index = 0; rel_eff_curve_index < num_rel_eff_curves; ++rel_eff_curve_index )
  
  if( !show_ref_lines || nuclides.empty() )
  {
    if( m_photopeak_widget )
      m_photopeak_widget->clearAllLines();
  }else
  {
    if( !m_photopeak_widget )
    {
      MaterialDB *materialDb = m_interspec->materialDataBase();
      Wt::WSuggestionPopup *materialSuggest = nullptr;
      WContainerWidget *parent = nullptr;
      m_photopeak_widget.reset( new ReferencePhotopeakDisplay( m_spectrum, materialDb, materialSuggest,
                                                              m_interspec, parent ) );
      
      // We need to set peaks getting assigned ref. line color, or else our call to
      //  #setColorsForSpecificSources will be useless.
      m_photopeak_widget->setPeaksGetAssignedRefLineColor( true );
    }//if( !m_photopeak_widget )
    
    assert( m_photopeak_widget );
    
    // First, clear out all the old lines
    m_photopeak_widget->clearAllLines();
    
    map<std::string,Wt::WColor> nuclide_colors;
    for( const RelActCalcAuto::NucInputInfo &nuc : nuclides )
    {
      if( (nuc.nuclide || nuc.element || nuc.reaction) && !nuc.peak_color_css.empty() )
        nuclide_colors[nuc.name()] = WColor( nuc.peak_color_css );
    }
    
    // If the user has a ColorTheme preference to assign a specific color to a nuclide, that
    //  will over-ride the colors we are about to set, but this is a minor detail to ignore at
    //  the moment.
    m_photopeak_widget->setColorsForSpecificSources( nuclide_colors );
    
    for( size_t i = 0; i < nuclides.size(); ++i )
    {
      const RelActCalcAuto::NucInputInfo &nuc = nuclides[i];
      
      if( i )
        m_photopeak_widget->persistCurentLines();
      
      if( nuc.nuclide )
      {
        m_photopeak_widget->setIsotope( nuc.nuclide, nuc.age );
      }else if( nuc.element )
      {
        m_photopeak_widget->setElement( nuc.element );
      }else if( nuc.reaction )
      {
        m_photopeak_widget->setReaction( nuc.reaction );
      }
    }//
  }//if( !show_ref_lines ) / else
}//void updateDuringRenderForRefGammaLineChange()


void RelActAutoGui::updateDuringRenderForFreePeakChange()
{
  // Check that free peaks are within ROIs, and if not, mark them as such

  if( m_free_peaks_container->isHidden() )
    return;
  
  vector<pair<float,float>> rois;
  const vector<WWidget *> &roi_widgets = m_energy_ranges->children();
  for( WWidget *w : roi_widgets )
  {
    const RelActAutoEnergyRange *roi = dynamic_cast<const RelActAutoEnergyRange *>( w );
    assert( roi );
    if( roi && !roi->isEmpty() )
      rois.emplace_back( roi->lowerEnergy(), roi->upperEnergy() );
  }//for( WWidget *w : kids )
  
  
  const bool fit_energy_cal = m_fit_energy_cal->isChecked();
  
  const std::vector<WWidget *> &kids = m_free_peaks->children();
  for( WWidget *w : kids )
  {
    RelActFreePeak *free_peak = dynamic_cast<RelActFreePeak *>(w);
    assert( free_peak );
    if( !free_peak )
      continue;
    
    const float energy = free_peak->energy();
    bool in_roi = false;
    for( const auto &roi : rois )
      in_roi = (in_roi || ((energy >= roi.first) && (energy <= roi.second)));
    
    free_peak->setInvalidEnergy( !in_roi );
    free_peak->setApplyEnergyCalVisible( fit_energy_cal );
  }//for( loop over RelActFreePeak widgets )
}//void updateDuringRenderForFreePeakChange()


void RelActAutoGui::updateDuringRenderForEnergyRangeChange()
{
  // Check if we need to show/hide/edit:
  // - m_presets
  // - m_rel_eff_eqn_order
  // - m_fwhm_eqn_form
  //
  // Need to make sure energy ranges arent overlapping tooo much
  
  const bool show_clear_ranges = (m_energy_ranges->count() > 0);
  if( m_clear_energy_ranges->isVisible() != show_clear_ranges )
    m_clear_energy_ranges->setHidden( !show_clear_ranges );
}//void updateDuringRenderForEnergyRangeChange()


void RelActAutoGui::startUpdatingCalculation()
{
  m_error_msg->setText("");
  m_error_msg->hide();
  m_fit_chi2_msg->setText("");
  m_fit_chi2_msg->hide();
  m_status_indicator->hide();
  
  // Disable being able to apply energy calibration fit from solution, until we get a new solution
  if( !m_apply_energy_cal_item->isDisabled() )
    m_apply_energy_cal_item->disable();
  
  shared_ptr<const SpecUtils::Measurement> foreground = m_foreground;
  shared_ptr<const SpecUtils::Measurement> background = m_background;
  RelActCalcAuto::Options options;
  
  try
  {
    if( !foreground )
      throw runtime_error( "No foreground spectrum is displayed." );
    
    if( !m_solution || !m_peak_model->rowCount() )
      makeZeroAmplitudeRoisToChart();
    
    if( m_background_subtract->isChecked() )
      background = m_interspec->displayedHistogram( SpecUtils::SpectrumType::Background );
    
    options = getCalcOptions();
    
    if( options.rel_eff_curves.empty() )
      throw runtime_error( "No relative efficiency curves defined." );

    for( const auto &rel_eff_curve : options.rel_eff_curves )
    {
      if( rel_eff_curve.nuclides.empty() )
        throw runtime_error( "No nuclides defined for relative efficiency curve." );
    }

    if( options.rois.empty() )
      throw runtime_error( "No energy ranges defined." );
  }catch( std::exception &e )
  {
    m_is_calculating = false;
    m_error_msg->setText( e.what() );
    m_error_msg->show();

    setOptionsForNoSolution();
    
    return;
  }//try / catch
  
  m_status_indicator->setText( "Calculating..." );
  m_status_indicator->show();
  
  if( m_cancel_calc )
    m_cancel_calc->store( true );
  
  const string sessionid = wApp->sessionId();
  m_is_calculating = true;
  m_cancel_calc = make_shared<std::atomic_bool>();
  m_cancel_calc->store( false );
  shared_ptr<atomic_bool> cancel_calc = m_cancel_calc;
  
  shared_ptr<const DetectorPeakResponse> cached_drf = m_cached_drf;
  if( !cached_drf )
  {
    auto m = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
    if( m && m->detector() && m->detector()->isValid() )
      cached_drf = m->detector();
  }
  
  WApplication *app = WApplication::instance();
  const string sessionId = app->sessionId();
  
  vector<shared_ptr<const PeakDef>> cached_all_peaks = m_cached_all_peaks;
  
  auto solution = make_shared<RelActCalcAuto::RelActAutoSolution>();
  auto error_msg = make_shared<string>();
  
  auto gui_update_callback = wApp->bind( boost::bind( &RelActAutoGui::updateFromCalc, this, solution, cancel_calc) );
  auto error_callback = wApp->bind( boost::bind( &RelActAutoGui::handleCalcException, this, error_msg, cancel_calc) );
  
  
  auto worker = [=](){
    try
    {
      RelActCalcAuto::RelActAutoSolution answer
        = RelActCalcAuto::solve( options, foreground, background, cached_drf, cached_all_peaks, cancel_calc );
      
      WServer::instance()->post( sessionId, [=](){
        WApplication *app = WApplication::instance();
        
        if( app )
        {
          *solution = answer;
          gui_update_callback();
          app->triggerUpdate();
        }else
        {
          cerr << "Failed to get WApplication::instance() for worker in RelActAutoGui::startUpdatingCalculation" << endl;
          assert( 0 );
        }//if( lock ) / else
      } );
    }catch( std::exception &e )
    {
      const string msg = e.what();
      cout << "Caught exception: " << msg << endl;
      cout << endl;
      
      WServer::instance()->post( sessionId, [=](){
        WApplication *app = WApplication::instance();
        
        if( app )
        {
          *error_msg = msg;
          error_callback();
          app->triggerUpdate();
        }else
        {
          cerr << "Failed to get WApplication::UpdateLock for worker in RelActAutoGui::startUpdatingCalculation" << endl;
          assert( 0 );
        }//if( lock ) / else
      } );
    }//try / catch
  };//auto worker
  
  m_calc_started.emit();
  
  WServer::instance()->ioService().boost::asio::io_service::post( worker );
}//void startUpdatingCalculation()


void RelActAutoGui::updateFromCalc( std::shared_ptr<RelActCalcAuto::RelActAutoSolution> answer,
                                    std::shared_ptr<std::atomic_bool> cancel_flag )
{
  assert( answer );
  if( !answer )
    return;
  
  // If we started a new calculation between when this one was started, and right now, dont
  //  do anything to the GUI state.
  if( cancel_flag != m_cancel_calc )
    return;
  
  if( cancel_flag && (*cancel_flag) )
    return;
  
  m_is_calculating = false;
  m_status_indicator->hide();
  
  switch( answer->m_status )
  {
    case RelActCalcAuto::RelActAutoSolution::Status::Success:
      break;
      
    case RelActCalcAuto::RelActAutoSolution::Status::NotInitiated:
    case RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem:
    case RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem:
    case RelActCalcAuto::RelActAutoSolution::Status::UserCanceled:
    {
      string msg = "Calculation didn't complete.";
      if( !answer->m_error_message.empty() )
        msg += ": " + answer->m_error_message;
      
      m_error_msg->setText( msg );
      m_error_msg->show();

      m_txt_results->setNoResults();

      setOptionsForNoSolution();
      
      return;
    }//if( calculation wasnt successful )
  }//switch( answer->m_status )
  
  if( answer->m_drf )
    m_cached_drf = answer->m_drf;
  
  if( !answer->m_spectrum_peaks.empty() )
    m_cached_all_peaks = answer->m_spectrum_peaks;
  
  m_solution = answer;
  
  m_txt_results->updateResults( *answer );
  
  m_peak_model->setPeaks( answer->m_fit_peaks_in_spectrums_cal );
  
  cout << "\n\n\nCalc finished: \n";
  answer->print_summary( std::cout );
  cout << "\n\n\n";
  
  const string rel_eff_eqn_js = answer->rel_eff_eqn_js_function(0);
  
  const double live_time = answer->m_foreground ? answer->m_foreground->live_time() : 1.0f;

  
  WString chi2_title("χ²/dof = {1}/{2}{3}");
  chi2_title.arg( SpecUtils::printCompact(answer->m_chi2, 3) )
            .arg( static_cast<int>(answer->m_dof) );

  // If we have U or Pu, we'll give the enrichment, or if we have two nuclides we'll
  //  give their ratio
  set<const SandiaDecay::Nuclide *> isotopes;
  for( const auto &relact : answer->m_rel_activities[0] )
  {
    if( relact.nuclide )
      isotopes.insert( relact.nuclide );
  }
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  assert( db );
  const SandiaDecay::Nuclide * const u235 = db->nuclide( "U235" );
  const SandiaDecay::Nuclide * const u238 = db->nuclide( "U238" );
  const SandiaDecay::Nuclide * const pu239 = db->nuclide( "Pu239" );
  const SandiaDecay::Nuclide * const pu240 = db->nuclide( "Pu240" );
  assert( u235 && u238 && pu239 && pu240 );
  
  if( (isotopes.count(u235) && isotopes.count(u238) && !isotopes.count(pu239))
     || (isotopes.count(pu239) && isotopes.count(pu240) && !isotopes.count(u235)) )
  {
    const SandiaDecay::Nuclide * const iso = isotopes.count(u235) ? u235 : pu239;
    const double nominal = answer->mass_enrichment_fraction( iso, 0 );
    // TODO: add errors
    string enrich = ", " + SpecUtils::printCompact(100.0*nominal, 4)
                        // + "±" + SpecUtils::printCompact(100.0*error, 4)
                         + "% " + iso->symbol;
    chi2_title.arg( enrich );
  }else if( isotopes.size() == 2 )
  {
    const vector<RelActCalcAuto::NuclideRelAct> &rel_acts = answer->m_rel_activities.at(0);
    const RelActCalcAuto::NuclideRelAct *num_rel_act = nullptr, *denom_rel_act = nullptr;
    for( size_t i = 0; i < rel_acts.size(); ++i )
    {
      if( rel_acts[i].nuclide )
      {
        if( !num_rel_act )
          num_rel_act = &(rel_acts[i]);
        else
          denom_rel_act = &(rel_acts[i]);
      }
    }//for( size_t i = 0; i < rel_acts.size(); ++i )
    
    assert( num_rel_act && denom_rel_act );
    
    if( num_rel_act && denom_rel_act )
    {
      if( num_rel_act->rel_activity > denom_rel_act->rel_activity )
        std::swap( num_rel_act, denom_rel_act );
      
      const SandiaDecay::Nuclide * const num_nuc = num_rel_act->nuclide;
      const SandiaDecay::Nuclide * const den_nuc = denom_rel_act->nuclide;
      assert( num_nuc && den_nuc );
      
      const double ratio = answer->activity_ratio( num_nuc, den_nuc, 0 );
      // TODO: add errors
      string ratio_txt = ", act(" + num_nuc->symbol + "/" + den_nuc->symbol + ")="
      + SpecUtils::printCompact(ratio, 4);
      
      chi2_title.arg( ratio_txt );
    }else
    {
      chi2_title.arg( "" );
    }
  }else
  {
    chi2_title.arg( "" );
  }

  string rel_eff_eqn_txt;
  try
  {
    rel_eff_eqn_txt = "y = " + answer->rel_eff_txt( false, 0);
  }catch( std::exception & )
  {
    // Oh well
    assert( 0 );
  }

  m_rel_eff_chart->setData( live_time, answer->m_fit_peaks, answer->m_rel_activities.at(0),
                           rel_eff_eqn_js, rel_eff_eqn_txt );

  m_fit_chi2_msg->setText( chi2_title );
  m_fit_chi2_msg->show();

  bool any_nucs_updated = false;
  for( const RelActCalcAuto::NuclideRelAct &fit_nuc : m_solution->m_rel_activities.at(0) )
  {
    if( !fit_nuc.age_was_fit || (fit_nuc.age < 0.0) )
      continue;
    
    any_nucs_updated = true;
    for( WWidget *w : m_nuclides->children() )
    {
      RelActAutoNuclide *nuc = dynamic_cast<RelActAutoNuclide *>( w );
      assert( nuc );
      
      if( !nuc || std::holds_alternative<std::monostate>(nuc->nuclide()) )
        continue;
      
      const std::variant<std::monostate, const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *>
           src_info = nuc->nuclide();
      const SandiaDecay::Nuclide *nuclide = nullptr;
      const SandiaDecay::Element *element = nullptr;
      const ReactionGamma::Reaction * reaction = nullptr;
      
      std::visit( [&]( auto &&arg ) {
        using T = std::decay_t<decltype(arg)>;
        if constexpr( std::is_same_v<T,std::monostate> )
        {
          
        }if constexpr( std::is_same_v<T, const SandiaDecay::Nuclide *> )
        {
          nuclide = std::get<const SandiaDecay::Nuclide *>(src_info);
        }else if constexpr( std::is_same_v<T, const SandiaDecay::Element *> )
        {
          element = std::get<const SandiaDecay::Element *>(src_info);
        }else if constexpr( std::is_same_v<T, const ReactionGamma::Reaction *> )
        {
          reaction = std::get<const ReactionGamma::Reaction *>(src_info);
        }else
        {
          assert( 0 );
          //static_assert( false, "Non-exhaustive visitor.");
        }
      }, src_info );

      if( (fit_nuc.nuclide == nuclide)
          && (fit_nuc.element == element)
          && (fit_nuc.reaction == reaction) )
      {
        const string agestr = PhysicalUnitsLocalized::printToBestTimeUnits( fit_nuc.age, 3 );
        nuc->setAge( agestr );
        break;
      }//if( this is the widget for this nuclide )
    }//for( WWidget *w : kids )
  }//for( const RelActCalcAuto::NuclideRelAct &fit_nuc : m_solution->m_rel_activities )
  
  
  setOptionsForValidSolution();
  
  // Check if we need to update Physical model shieldings
  assert( m_solution->m_options.rel_eff_curves.size() >= 1 );
  const size_t num_rel_eff_curves = std::max( size_t(1), m_solution->m_options.rel_eff_curves.size() );
  assert( m_rel_eff_opts_menu->count() == static_cast<int>(num_rel_eff_curves) );

  // Remove any extra rel eff curve options
  while( m_rel_eff_opts_menu->count() > num_rel_eff_curves )
  {
    RelActAutoGuiRelEffOptions *opts = getRelEffCurveOptions( m_rel_eff_opts_menu->count() - 1 );
    assert( opts );
    handleDelRelEffCurve( opts );
  }

  // Add any missing rel eff curve options
  while( m_rel_eff_opts_menu->count() < num_rel_eff_curves )
  {
    handleAddRelEffCurve();
  }

  for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff_curves; ++rel_eff_index )
  {
    if( rel_eff_index >= m_solution->m_options.rel_eff_curves.size() )
    {
      continue; // I guess do nothing...
    }

    RelActAutoGuiRelEffOptions *opts = getRelEffCurveOptions( rel_eff_index );
    assert( opts );
    if( !opts )
      throw std::runtime_error( "Failed to get RelActAutoGuiRelEffOptions" );


    const RelActCalcAuto::RelEffCurveInput &rel_eff = m_solution->m_options.rel_eff_curves[rel_eff_index];
    
    std::optional<RelActCalcAuto::RelActAutoSolution::PhysicalModelFitInfo> phys_info;
    if( rel_eff_index < m_solution->m_phys_model_results.size() )
      phys_info = m_solution->m_phys_model_results[rel_eff_index];

    assert( (rel_eff.rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel) || phys_info.has_value() );
    if( (rel_eff.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel) && phys_info.has_value() )
    {
      const RelActCalcAuto::RelActAutoSolution::PhysicalModelFitInfo &phys_info_ref = *phys_info;
      opts->update_shield_widgets( phys_info_ref.self_atten, phys_info_ref.ext_shields );
   
      //phys_info_ref.hoerl_b/hoerl_c will have (double) values iff the correction function was used
      assert( phys_info_ref.hoerl_b.has_value() == phys_info_ref.hoerl_c.has_value() );
      assert( phys_info_ref.hoerl_b.has_value() == opts->phys_model_use_hoerl() );
      if( phys_info_ref.hoerl_b.has_value() != opts->phys_model_use_hoerl() )
        opts->setPhysModelUseHoerl( phys_info_ref.hoerl_b.has_value() );
    }//if( Physical model fit info is available )
  }//for( size_t rel_eff_index = 0; rel_eff_index < num_rel_eff_curves; ++rel_eff_index )

  
  m_solution_updated.emit( m_solution );
  if( m_solution && (m_solution->m_status == RelActCalcAuto::RelActAutoSolution::Status::Success) )
    m_calc_successful.emit();
  else
    m_calc_failed.emit();
  
  if( any_nucs_updated )
  {
    m_render_flags |= RenderActions::UpdateRefGammaLines;
    scheduleRender();
  }
}//void updateFromCalc( std::shared_ptr<RelActCalcAuto::RelActAutoSolution> answer )


void RelActAutoGui::handleCalcException( std::shared_ptr<std::string> message,
                                        std::shared_ptr<std::atomic_bool> cancel_flag )
{
  assert( message );
  if( !message )
    return;
  
  // If we started a new calculation between when this one was started, and right now, dont
  //  do anything to the GUI state.
  if( cancel_flag != m_cancel_calc )
    return;
  
  m_is_calculating = false;
  m_status_indicator->hide();
  
  string msg = "Calculation error: ";
  msg += *message;
  
  m_error_msg->setText( msg );
  m_error_msg->show();

  m_fit_chi2_msg->setText( "" ); // should already be hidden, but JIC
  m_fit_chi2_msg->hide();

  m_solution.reset();
  
  setOptionsForNoSolution();
  
  m_solution_updated.emit( m_solution );
  m_calc_failed.emit();
}//void handleCalcException( std::shared_ptr<std::string> message )
