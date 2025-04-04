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
#include <Wt/WMenuItem>
#include <Wt/WResource>
#include <Wt/WIOService>
#include <Wt/WFileUpload>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
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
      label->setBuddy( m_lower_energy );
      
      label = new WLabel( "keV", this );
      label->addStyleClass( "GridThirdCol GridFirstRow" );
      
      label = new WLabel( "Upper Energy", this );
      label->addStyleClass( "GridFirstCol GridSecondRow" );
      
      m_upper_energy = new NativeFloatSpinBox( this );
      m_upper_energy->addStyleClass( "GridSecondCol GridSecondRow" );
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
      ", <b>Cs137</b>, etc.";
      HelpSystem::attachToolTipOn( m_nuclide_edit, tooltip, showToolTips );
      
      string replacerJs, matcherJs;
      IsotopeNameFilterModel::replacerJs( replacerJs );
      IsotopeNameFilterModel::nuclideNameMatcherJs( matcherJs );
      IsotopeNameFilterModel *isoSuggestModel = new IsotopeNameFilterModel( this );
      isoSuggestModel->excludeXrays( true );
      isoSuggestModel->excludeEscapes( true );
      isoSuggestModel->excludeReactions( true );
      
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
      const SandiaDecay::Nuclide *nuc = nuclide();
      if( !nuc )
      {
        const string nucstr = m_nuclide_edit->text().toUTF8();
        
        m_fit_age->setUnChecked();
        m_fit_age->hide();
        
        if( !nucstr.empty() )
          passMessage( nucstr + " is not a valid nuclide.", WarningWidget::WarningMsgHigh );
        
        m_updated.emit();
        return;
      }//if( !nuc )
      
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
      
      bool haveFoundColor = false;
      
      // Check if user is showing reference lines for this nuclide, and if so, use that color
      if( !haveFoundColor )
      {
        vector<ReferenceLineInfo> reflines;
        const ReferencePhotopeakDisplay *refdisp = InterSpec::instance()->referenceLinesWidget();
        if( refdisp )
          reflines = refdisp->showingNuclides();
        
        for( const auto &refline : reflines )
        {
          if( (refline.m_nuclide == nuc) && !refline.m_input.m_color.isDefault() )
          {
            haveFoundColor = true;
            m_color_select->setColor( refline.m_input.m_color );
            break;
          }//if( we found a match )
        }//for( const auto &refline : reflines )
        
        // TODO: look in ReferencePhotopeakDisplay cache of previous colors
      }//if( !haveFoundColor )
      
      // Check for user-fit peaks
      if( !haveFoundColor )
      {
        PeakModel *peakmodel = InterSpec::instance()->peakModel();
        if( peakmodel && peakmodel->peaks() )
        {
          for( auto peakshrdptr : *peakmodel->peaks() )
          {
            if( !peakshrdptr )
              continue;
            
            if( (peakshrdptr->parentNuclide() == nuc) && !peakshrdptr->lineColor().isDefault() )
            {
              haveFoundColor = true;
              m_color_select->setColor( peakshrdptr->lineColor() );
              break;
            }//if( we found a matching nuclide )
          }//for( loop over user peaks )
        }//if( PeakModel is valid  - should always be )
      }//if( !haveFoundColor )
      
      // Finally check
      if( !haveFoundColor )
      {
        shared_ptr<const ColorTheme> theme = InterSpec::instance()->getColorTheme();
        if( theme )
        {
          const map<string,WColor> &defcolors = theme->referenceLineColorForSources;
          const auto pos = defcolors.find( nuc->symbol );
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
    
    /** Will return nullptr if invalid text entered. */
    const SandiaDecay::Nuclide *nuclide() const
    {
      const string nucstr = m_nuclide_edit->text().toUTF8();
      const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
      
      return db->nuclide( nucstr );
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
      const SandiaDecay::Nuclide * const nuc = nuclide();
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
      const SandiaDecay::Nuclide * const nuc = nuclide();
      if( !nuc )
        throw runtime_error( "No valid nuclide" );
      
      RelActCalcAuto::NucInputInfo nuc_info;

      nuc_info.nuclide = nuc;
      nuc_info.age = age(); // Must not be negative.
      nuc_info.fit_age = m_fit_age->isChecked();
      
      // TODO: Not implemented: vector<double> gammas_to_exclude;
      //nuc_info.gammas_to_exclude = ;
      
      nuc_info.peak_color_css = m_color_select->color().cssText();

      return nuc_info;
    }//RelActCalcAuto::RoiRange toRoiRange() const
    
    
    void fromNucInputInfo( const RelActCalcAuto::NucInputInfo &info )
    {
      if( !info.nuclide )
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
      }//if( !info.nuclide )
      
      m_nuclide_edit->setText( WString::fromUTF8(info.nuclide->symbol) );
      handleIsotopeChange();
      
      if( !info.peak_color_css.empty() )
        m_color_select->setColor( WColor(info.peak_color_css) );
      
      const SandiaDecay::Nuclide * const nuc = nuclide();
      
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
      // Not currently supported: info.gammas_to_exclude -> vector<double>;
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
  m_rel_eff_eqn_form( nullptr ),
  m_rel_eff_eqn_order_label( nullptr ),
  m_rel_eff_eqn_order( nullptr ),
  m_fwhm_eqn_form( nullptr ),
  m_fwhm_estimation_method( nullptr ),
  m_fit_energy_cal( nullptr ),
  m_background_subtract( nullptr ),
  m_same_z_age( nullptr ),
  m_pu_corr_method( nullptr ),
  m_skew_type( nullptr ),
  m_add_uncert( nullptr ),
  m_phys_model_opts( nullptr ),
  m_phys_model_shields( nullptr ),
  m_phys_model_self_atten( nullptr ),
  m_phys_ext_attens( nullptr ),
  m_phys_model_use_hoerl( nullptr ),
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
  WLabel *label = new WLabel( "Option Presets", presetDiv );
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
  
  m_error_msg = new WText( "Not Calculated.", presetDiv );
  m_error_msg->addStyleClass( "RelActAutoErrMsg" );

  m_fit_chi2_msg = new WText( "", presetDiv );
  m_fit_chi2_msg->addStyleClass( "RelActAutoChi2Msg" );

  m_status_indicator = new WText( "Calculating...", presetDiv );
  m_status_indicator->addStyleClass( "RelActAutoStatusMsg" );
  m_status_indicator->hide();
  
  
  WContainerWidget *optionsDiv = new WContainerWidget( this );
  optionsDiv->addStyleClass( "RelActAutoOptions" );
  
  label = new WLabel( "Eqn Type", optionsDiv );
  label->addStyleClass( "GridFirstCol GridFirstRow" );
  
  m_rel_eff_eqn_form = new WComboBox( optionsDiv );
  m_rel_eff_eqn_form->addStyleClass( "GridSecondCol GridFirstRow" );
  label->setBuddy( m_rel_eff_eqn_form );
  m_rel_eff_eqn_form->activated().connect( this, &RelActAutoGui::handleRelEffEqnFormChanged );
  
  const char *tooltip = "The functional form to use for the relative efficiciency curve.<br />"
  "Options are:"
  "<table style=\"margin-left: 10px;\">"
  "<tr><th>Log(energy):</th>               <th>y = a + b*ln(x) + c*(ln(x))^2 + d*(ln(x))^3 + ...</th></tr>"
  "<tr><th>Log(rel. eff.):</th>            <th>y = exp( a + b*x + c/x + d/x^2 + e/x^3 + ... )</th></tr>"
  "<tr><th>Log(energy)Log(rel. eff.):</th> <th>y = exp( a  + b*(lnx) + c*(lnx)^2 + d*(lnx)^3 + ... )</th></tr>"
  "<tr><th>FRAM Empirical:</th>            <th>y = exp( a + b/x^2 + c*(lnx) + d*(lnx)^2 + e*(lnx)^3 )</th></tr>"
  "</table>";
  HelpSystem::attachToolTipOn( {label,m_rel_eff_eqn_form}, tooltip, showToolTips );
  
  // Will assume FramEmpirical is the highest
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
        txt = "Physical";
        break;
    }
    
    m_rel_eff_eqn_form->addItem( txt );
  }//for( loop over RelEffEqnForm )
  
  m_rel_eff_eqn_form->setCurrentIndex( static_cast<int>(RelActCalc::RelEffEqnForm::LnX) );
  
  m_rel_eff_eqn_order_label = new WLabel( "Eqn Order", optionsDiv );
  m_rel_eff_eqn_order_label->addStyleClass( "GridThirdCol GridFirstRow" );
  
  m_rel_eff_eqn_order = new WComboBox( optionsDiv );
  m_rel_eff_eqn_order->addStyleClass( "GridFourthCol GridFirstRow" );
  m_rel_eff_eqn_order_label->setBuddy( m_rel_eff_eqn_order );
  m_rel_eff_eqn_order->activated().connect( this, &RelActAutoGui::handleRelEffEqnOrderChanged );
  
  m_rel_eff_eqn_order->addItem( "0" );
  m_rel_eff_eqn_order->addItem( "1" );
  m_rel_eff_eqn_order->addItem( "2" );
  m_rel_eff_eqn_order->addItem( "3" );
  m_rel_eff_eqn_order->addItem( "4" );
  m_rel_eff_eqn_order->addItem( "5" );
  m_rel_eff_eqn_order->addItem( "6" );
  m_rel_eff_eqn_order->setCurrentIndex( 3 );
  
  
  tooltip = "The order (how many energy-dependent terms) relative efficiency equation to use.";
  HelpSystem::attachToolTipOn( {label, m_rel_eff_eqn_order}, tooltip, showToolTips );
  

  label = new WLabel( "FWHM Est.", optionsDiv );
  label->addStyleClass( "GridFirstCol GridSecondRow" );
  m_fwhm_estimation_method = new WComboBox( optionsDiv );
  m_fwhm_estimation_method->addStyleClass( "GridSecondCol GridSecondRow" );
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

  label = new WLabel( "FWHM Form", optionsDiv );
  label->addStyleClass( "GridThirdCol GridSecondRow" );
  
  m_fwhm_eqn_form = new WComboBox( optionsDiv );
  m_fwhm_eqn_form->addStyleClass( "GridFourthCol GridSecondRow" );
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
  
  tooltip = "The equation type used to model peak FWHM as a function of energy.";
  HelpSystem::attachToolTipOn( {label, m_fwhm_eqn_form}, tooltip, showToolTips );
  
  // TODO: need to set m_fwhm_eqn_form based on energy ranges selected
  m_fwhm_eqn_form->setCurrentIndex( static_cast<int>(RelActCalcAuto::FwhmForm::SqrtEnergyPlusInverse) );
  m_fwhm_estimation_method->setCurrentIndex( static_cast<int>(RelActCalcAuto::FwhmEstimationMethod::StartFromDetEffOrPeaksInSpectrum) );
  m_fwhm_eqn_form->changed().connect( this, &RelActAutoGui::handleFwhmFormChanged );
  m_fwhm_estimation_method->changed().connect( this, &RelActAutoGui::handleFwhmEstimationMethodChanged );
  

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
  
  m_fit_energy_cal = new WCheckBox( "Fit Energy Cal.", optionsDiv );
  m_fit_energy_cal->addStyleClass( "GridFirstCol GridThirdRow GridSpanTwoCol" );
  m_fit_energy_cal->checked().connect( this, &RelActAutoGui::handleFitEnergyCalChanged );
  m_fit_energy_cal->unChecked().connect( this, &RelActAutoGui::handleFitEnergyCalChanged );
  
  
  m_background_subtract = new WCheckBox( "Back. Sub.", optionsDiv );
  m_background_subtract->addStyleClass( "GridThirdCol GridThirdRow GridSpanTwoCol" );
  m_background_subtract->checked().connect( this, &RelActAutoGui::handleBackgroundSubtractChanged );
  m_background_subtract->unChecked().connect( this, &RelActAutoGui::handleBackgroundSubtractChanged );
  
  
  m_same_z_age = new WCheckBox( "Same El. Same Age", optionsDiv );
  m_same_z_age->addStyleClass( "GridFifthCol GridThirdRow GridSpanTwoCol" );
  m_same_z_age->checked().connect( this, &RelActAutoGui::handleSameAgeChanged );
  m_same_z_age->unChecked().connect( this, &RelActAutoGui::handleSameAgeChanged );
  
  label = new WLabel( "Pu242 corr", optionsDiv );
  label->addStyleClass( "GridSeventhCol GridSecondRow" );
  m_pu_corr_method = new WComboBox( optionsDiv );
  m_pu_corr_method->addStyleClass( "GridEighthCol GridSecondRow" );
  label->setBuddy( m_pu_corr_method );
  m_pu_corr_method->activated().connect( this, &RelActAutoGui::handlePuByCorrelationChanged );
  tooltip = "Pu-242 is often not directly observable in gamma spectra.  However, to"
  " correct for this isotope when calculating enrichment, the expected contributions of this"
  " isotope can be inferred from the other Pu isotopes."
  "  This form allows you to select the correction method.";
  HelpSystem::attachToolTipOn( {label,m_pu_corr_method}, tooltip, showToolTips );
  m_pu_corr_method->hide();
  label->hide();
  
    
  label = new WLabel( "Peak Skew", optionsDiv );
  label->addStyleClass( "GridFifthCol GridSecondRow" );
  m_skew_type = new WComboBox( optionsDiv );
  m_skew_type->addStyleClass( "GridSixthCol GridSecondRow" );
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
  
    
  label = new WLabel( "Add. Uncert", optionsDiv );
  label->addStyleClass( "GridFifthCol GridFirstRow" );
  m_add_uncert = new WComboBox( optionsDiv );
  m_add_uncert->addStyleClass( "GridSixthCol GridFirstRow" );
  label->setBuddy( m_add_uncert );
    m_add_uncert->activated().connect( this, &RelActAutoGui::handleAdditionalUncertChanged );
    
  for( RelActAutoGui::AddUncert i = RelActAutoGui::AddUncert(0);
      i < RelActAutoGui::AddUncert::NumAddUncert;
      i = RelActAutoGui::AddUncert(static_cast<int>(i) + 1) )
  {
    WString uncert_txt;
    switch( i )
    {
      case AddUncert::StatOnly:           uncert_txt = WString::fromUTF8("None"); break;
      case AddUncert::OnePercent:         uncert_txt = WString::fromUTF8("1%");   break;
      case AddUncert::FivePercent:        uncert_txt = WString::fromUTF8("5%");   break;
      case AddUncert::TenPercent:         uncert_txt = WString::fromUTF8("10%");  break;
      case AddUncert::TwentyFivePercent:  uncert_txt = WString::fromUTF8("25%");  break;
      case AddUncert::FiftyPercent:       uncert_txt = WString::fromUTF8("50%");  break;
      case AddUncert::SeventyFivePercent: uncert_txt = WString::fromUTF8("75%");  break;
      case AddUncert::OneHundredPercent:  uncert_txt = WString::fromUTF8("100%"); break;
      case AddUncert::NumAddUncert:       assert(0);                              break;
    }//switch( i )
       
    m_add_uncert->addItem( uncert_txt );
  }//for( loop over AddUncert )
     
  m_add_uncert->setCurrentIndex( static_cast<int>(RelActAutoGui::AddUncert::StatOnly) );
    
    
  m_phys_model_opts = new WContainerWidget( optionsDiv );
  m_phys_model_opts->addStyleClass( "PhysicalModelOpts GridEleventhCol GridFirstRow GridSpanThreeRows" );
  
  m_phys_model_shields = new WContainerWidget( m_phys_model_opts );
  m_phys_model_shields->addStyleClass( "PhysicalModelShields" );
    
  m_phys_ext_attens = new WContainerWidget( m_phys_model_shields );
  
  // We will pu the Use 
  WContainerWidget *phys_opt_row = new WContainerWidget( m_phys_model_opts );
  phys_opt_row->setStyleClass( "PhysicalModelOptRow" );
  
  m_phys_model_use_hoerl = new WCheckBox( "Use Corr. Fcn.", phys_opt_row );
  m_phys_model_use_hoerl->addStyleClass( "UseCorrFcnCb" );
  m_phys_model_use_hoerl->setChecked( true );
  m_phys_model_use_hoerl->setWordWrap( false );
  m_phys_model_use_hoerl->checked().connect( this, &RelActAutoGui::handlePhysModelUseHoerlChange );
  m_phys_model_use_hoerl->unChecked().connect( this, &RelActAutoGui::handlePhysModelUseHoerlChange );

  SpecMeasManager *spec_manager = m_interspec->fileManager();
  SpectraFileModel *spec_model = spec_manager ? spec_manager->model() : nullptr;
  DetectorDisplay *det_disp = new DetectorDisplay( m_interspec, spec_model, phys_opt_row );
  m_interspec->detectorChanged().connect( this, &RelActAutoGui::handleDetectorChange );
  m_interspec->detectorModified().connect( this, &RelActAutoGui::handleDetectorChange );

  m_phys_model_opts->hide();
    
    
  WPushButton *more_btn = new WPushButton( optionsDiv );
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
    
    
  // TODO: add item to account for Rel Eff, and Rel Acts
  
  WContainerWidget *bottomArea = new WContainerWidget( this );
  bottomArea->addStyleClass( "EnergiesAndNuclidesHolder" );
  
  WContainerWidget *nuclidesHolder = new WContainerWidget( bottomArea );
  nuclidesHolder->addStyleClass( "NuclidesHolder" );
  
  WContainerWidget *energiesHolder = new WContainerWidget( bottomArea );
  energiesHolder->addStyleClass( "EnergiesHolder" );
  
  m_free_peaks_container = new WContainerWidget( bottomArea );
  m_free_peaks_container->addStyleClass( "FreePeaksHolder" );
  m_free_peaks_container->hide();
  
  WText *nuc_header = new WText( "Nuclides", nuclidesHolder );
  nuc_header->addStyleClass( "EnergyNucHeader" );
  
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
  
  WText *energy_header = new WText( "Energy Ranges", energiesHolder );
  energy_header->addStyleClass( "EnergyNucHeader" );
  
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
  
  m_show_free_peak = new WPushButton( "add free peaks", energies_footer );
  m_show_free_peak->addStyleClass( "ShowFreePeaks LightButton" );
  m_show_free_peak->clicked().connect( this, &RelActAutoGui::handleShowFreePeaks );
  
  m_clear_energy_ranges = new WPushButton( "clear all ranges", energies_footer );
  m_clear_energy_ranges->addStyleClass( "ClearEnergyRanges LightButton" );
  m_clear_energy_ranges->clicked().connect( this, &RelActAutoGui::handleClearAllEnergyRanges );
  tooltip = "Removes all energy ranges.";
  HelpSystem::attachToolTipOn( m_clear_energy_ranges, tooltip, showToolTips );
  m_clear_energy_ranges->hide();
  
  WText *free_peaks_header = new WText( "Free Peaks", m_free_peaks_container );
  free_peaks_header->addStyleClass( "EnergyNucHeader" );
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
    case AddUncert::StatOnly:           options.additional_br_uncert = 0.00;  break;
    case AddUncert::OnePercent:         options.additional_br_uncert = 0.01;  break;
    case AddUncert::FivePercent:        options.additional_br_uncert = 0.05;  break;
    case AddUncert::TenPercent:         options.additional_br_uncert = 0.10;  break;
    case AddUncert::TwentyFivePercent:  options.additional_br_uncert = 0.25;  break;
    case AddUncert::FiftyPercent:       options.additional_br_uncert = 0.50;  break;
    case AddUncert::SeventyFivePercent: options.additional_br_uncert = 0.75;  break;
    case AddUncert::OneHundredPercent:  options.additional_br_uncert = 1.00;  break;
    case AddUncert::NumAddUncert:       assert( 0 );                          break;
  }//switch( add_uncert )
  assert( options.additional_br_uncert >= 0.0 );
  options.additional_br_uncert = std::max( options.additional_br_uncert, 0.0 );
  
  options.floating_peaks = getFloatingPeaks();
  options.rois = getRoiRanges();


  RelActCalcAuto::RelEffCurveInput rel_eff_curve;
  rel_eff_curve.nuclides = getNucInputInfo();
  rel_eff_curve.nucs_of_el_same_age = m_same_z_age->isChecked();
  rel_eff_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm( std::max( 0, m_rel_eff_eqn_form->currentIndex() ) );
  rel_eff_curve.rel_eff_eqn_order = 0;
  if( rel_eff_curve.rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    rel_eff_curve.rel_eff_eqn_order = static_cast<size_t>( std::max( 0, m_rel_eff_eqn_order->currentIndex() ) );
  
  if( rel_eff_curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    if( m_phys_model_self_atten && m_phys_model_self_atten->nonEmpty() )
      rel_eff_curve.phys_model_self_atten = m_phys_model_self_atten->fitInput();
    
    for( WWidget *w : m_phys_ext_attens->children() )
    {
      RelEffShieldWidget *sw = dynamic_cast<RelEffShieldWidget *>( w );
      if( sw && sw->nonEmpty() )
      {
        auto input = sw->fitInput();
        if( input )
          rel_eff_curve.phys_model_external_atten.push_back( input );
      }
    }//for( WWidget *w : m_phys_ext_attens->children() )
    
    rel_eff_curve.phys_model_use_hoerl = m_phys_model_use_hoerl->isChecked();
  }//if( options.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )

  if( m_pu_corr_method->isVisible() )
  {
    const string currtxt = m_pu_corr_method->currentText().toUTF8();
    for( int i = 0; i <= static_cast<int>(RelActCalc::PuCorrMethod::NotApplicable); ++i )
    {
      const auto method = RelActCalc::PuCorrMethod(i);
      const string &desc = RelActCalc::to_description( method );
      if( desc == currtxt )
      {
        rel_eff_curve.pu242_correlation_method = method;
        break;
      }
    }//for( loop over RelActCalc::PuCorrMethod )
  }//if( m_pu_corr_method->isVisible() )

  options.rel_eff_curves = { rel_eff_curve };

  return options;
}//RelActCalcAuto::Options getCalcOptions() const


vector<RelActCalcAuto::NucInputInfo> RelActAutoGui::getNucInputInfo() const
{
  vector<RelActCalcAuto::NucInputInfo> answer;
  
  const vector<WWidget *> &kids = m_nuclides->children();
  for( WWidget *w : kids )
  {
    const RelActAutoNuclide *nuc = dynamic_cast<const RelActAutoNuclide *>( w );
    assert( nuc );
    if( nuc && nuc->nuclide() )
      answer.push_back( nuc->toNucInputInfo() );
  }//for( WWidget *w : kids )
  
  return answer;
}//RelActCalcAuto::NucInputInfo getNucInputInfo() const


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
  assert( options.rel_eff_curves.size() == 1 );
  if( options.rel_eff_curves.size() != 1 )
    throw runtime_error( "RelActAutoGui::setCalcOptionsGui: for dev, must have exactly one rel-eff curve." );
  
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
  if( options.additional_br_uncert <= 0.005 )
    add_uncert = AddUncert::StatOnly;
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
  
  const RelActCalcAuto::RelEffCurveInput &rel_eff = options.rel_eff_curves[0];
  m_same_z_age->setChecked( rel_eff.nucs_of_el_same_age );
  m_rel_eff_eqn_form->setCurrentIndex( static_cast<int>(rel_eff.rel_eff_eqn_type) );
  if( rel_eff.rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    m_rel_eff_eqn_order->setCurrentIndex( static_cast<int>(rel_eff.rel_eff_eqn_order) );
  
  const shared_ptr<const RelActCalc::PhysicalModelShieldInput> &self_atten = rel_eff.phys_model_self_atten;
  const vector<shared_ptr<const RelActCalc::PhysicalModelShieldInput>> &ext_attens = rel_eff.phys_model_external_atten;
  
  if( (rel_eff.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel)
     || self_atten || !ext_attens.empty() )
  {
    initPhysModelShields();  //Sets shielding widgets to default state
    
    assert( m_phys_model_self_atten );
    if( !m_phys_model_self_atten || !m_phys_ext_attens )
      throw std::logic_error( "!m_phys_model_self_atten || !m_phys_ext_attens ?!?" );
    assert( m_phys_ext_attens->count() == 1 );
    
    
    if( self_atten )
    {
      RelEffShieldState state;
      state.setStateFromFitInput( *self_atten );
      m_phys_model_self_atten->setState( state );
    }
    
    // Fill out the external shielding widgets
    vector<RelEffShieldWidget *> ext_shields;
    for( WWidget *w : m_phys_ext_attens->children() )
    {
      RelEffShieldWidget *sw = dynamic_cast<RelEffShieldWidget *>( w );
      assert( sw );
      if( sw )
        ext_shields.push_back( sw );
    }
    
    assert( ext_shields.size() == 1 );
    
    // Remove any extra shielding widgets
    while( (ext_shields.size() > 1) && (ext_shields.size() > ext_attens.size()) )
    {
      delete ext_shields.back();
      ext_shields.pop_back();
    }
    
    if( (ext_shields.size()) == 1 && ext_attens.empty() )
      ext_shields[0]->resetState();
    
    // Set state, reusing as many shielding widgets as we can
    const size_t num_ext_reused = std::min( ext_shields.size(), ext_attens.size() );
    for( size_t i = 0; i < num_ext_reused; ++i )
    {
      RelEffShieldState state;
      state.setStateFromFitInput( *ext_attens[i] );
      ext_shields[i]->setState( state );
    }
    
    // Add any new shielding widgets we need
    for( size_t i = num_ext_reused; i < ext_attens.size(); ++i )
    {
      RelEffShieldWidget *sw = new RelEffShieldWidget( RelEffShieldWidget::ShieldType::ExternalAtten, m_phys_ext_attens );
      sw->changed().connect( this, &RelActAutoGui::handlePhysModelShieldChange );
      ext_shields.push_back( sw );
      
      RelEffShieldState state;
      state.setStateFromFitInput( *ext_attens[i] );
      sw->setState( state );
    }
  }else
  {
    // Normally we will keep Physical model shielding around, incase the user is
    //  messing around, but if we are setting state from XML, we'll hard reset them
    if( m_phys_model_self_atten )
    {
      delete m_phys_model_self_atten;
      m_phys_model_self_atten = nullptr;
    }
    while( !m_phys_ext_attens->children().empty() )
      delete m_phys_ext_attens->children().front();
  }//if( Physical Model ) / else

  showAndHideOptionsForEqnType();
  
  for( const RelActCalcAuto::NucInputInfo &nuc : rel_eff.nuclides )
  {
    RelActAutoNuclide *nuc_widget = new RelActAutoNuclide( this, m_nuclides );
    nuc_widget->updated().connect( this, &RelActAutoGui::handleNuclidesChanged );
    nuc_widget->remove().connect( boost::bind( &RelActAutoGui::handleRemoveNuclide,
                                              this, static_cast<WWidget *>(nuc_widget) ) );
    nuc_widget->fromNucInputInfo( nuc );
  }//for( const RelActCalcAuto::NucInputInfo &nuc : nuclides )
  
  
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
  
  
  // We will be a bit cheap here, and set one entry in pu correlation selector, and then
  //  rely on `updateDuringRenderForNuclideChange()` to input all the actual options
  m_pu_corr_method->setHidden( rel_eff.pu242_correlation_method == RelActCalc::PuCorrMethod::NotApplicable );
  m_pu_corr_method->clear();
  m_pu_corr_method->addItem( RelActCalc::to_description(rel_eff.pu242_correlation_method) );
  m_pu_corr_method->setCurrentIndex( 0 );
  if( m_pu_corr_method->label() )
    m_pu_corr_method->label()->setHidden( m_pu_corr_method->isHidden() );
  
  
  
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
  
  if( index <= 0 )
  {
    // Clear everything out!
    return;
  }
  
  if( index >= m_preset_paths.size() )
  {
    // TODO: let users download config, or clone them
    
    m_nuclides->clear();
    m_energy_ranges->clear();
    
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
  showAndHideOptionsForEqnType();
  
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


void RelActAutoGui::initPhysModelShields()
{
  if( m_phys_model_self_atten )
  {
    m_phys_model_self_atten->resetState();
  }else
  {
    m_phys_model_self_atten = new RelEffShieldWidget( RelEffShieldWidget::ShieldType::SelfAtten );
    m_phys_model_shields->insertWidget( 0, m_phys_model_self_atten );
    m_phys_model_self_atten->changed().connect( this, &RelActAutoGui::handlePhysModelShieldChange );
  }
  
  vector<RelEffShieldWidget *> starting_ext_shields;
  const vector<WWidget *> &ext_atten_widgets = m_phys_ext_attens->children();
  for( WWidget *w : ext_atten_widgets )
  {
    RelEffShieldWidget *sw = dynamic_cast<RelEffShieldWidget *>( w );
    assert( sw );
    if( sw )
      starting_ext_shields.push_back( sw );
  }
  
  if( starting_ext_shields.empty() )
  {
    RelEffShieldWidget *sw = new RelEffShieldWidget( RelEffShieldWidget::ShieldType::ExternalAtten,
                                                    m_phys_ext_attens );
    sw->changed().connect( this, &RelActAutoGui::handlePhysModelShieldChange );
  }else
  {
    starting_ext_shields[0]->resetState();
    for( size_t i = 1; i < starting_ext_shields.size(); ++i )
      delete starting_ext_shields[i];
  }//if( starting_ext_shields.empty() )
}//void initPhysModelShields()


void RelActAutoGui::handlePhysModelUseHoerlChange()
{
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handlePhysModelUseHoerlChange()


void RelActAutoGui::handlePhysModelShieldChange()
{
  checkIfInUserConfigOrCreateOne( false );
  m_render_flags |= RenderActions::UpdateCalculations;
  scheduleRender();
}//void handlePhysModelShieldChange()


void RelActAutoGui::showAndHideOptionsForEqnType()
{
  const RelActCalc::RelEffEqnForm eqn_type
                  = RelActCalc::RelEffEqnForm( std::max( 0, m_rel_eff_eqn_form->currentIndex() ) );
  
  const bool is_physical = (eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel);
  
  m_rel_eff_eqn_order->setDisabled( is_physical );
  m_rel_eff_eqn_order_label->setDisabled( is_physical );
  m_phys_model_opts->setHidden( !is_physical );
  if( is_physical && !m_phys_model_self_atten )
    initPhysModelShields();
}//void showAndHideOptionsForEqnType()


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
  // Check if we need to show/hide
  //  - m_same_z_age
  //  - m_u_pu_by_correlation - and also update its text
  //  - m_u_pu_data_source
  //
  // Need to check nuclides arent duplicated
  // Make sure "Fit Age" checkbox is in a correct state.  I.e.
  //  - make sure if each nuclide *could* have ages fit, that there is peaks from at least two progeny in the used energy ranges
  //  - make sure "Same El Same Age" is checked, then all nuclide ages for an element have same age and same m_fit_age value
  
  
  // For a first go, we'll show #m_same_z_age if we hadve more than one nuclide for a given Z.
  //   However, we really need to be more complex about this, and also hide "Fit Age" and disable age input if this option is selected, on all but the first entry of this Z
  bool has_multiple_nucs_of_z = false;
  map<short,int> z_to_num_isotopes;
  
  const vector<RelActCalcAuto::NucInputInfo> nuclides = getNucInputInfo();
  set<string> nuc_names;
  for( const RelActCalcAuto::NucInputInfo &nuc : nuclides )
  {
    if( !nuc.nuclide )
      continue;
    
    nuc_names.insert( nuc.nuclide->symbol );
    
    const short z = nuc.nuclide->atomicNumber;
    if( !z_to_num_isotopes.count(z) )
      z_to_num_isotopes[z] = 0;
    
    int &num_this_z = z_to_num_isotopes[z];
    num_this_z += 1;
    
    has_multiple_nucs_of_z = has_multiple_nucs_of_z || (num_this_z > 1);
  }//for( const RelActCalcAuto::NucInputInfo &nuc : nuclides )
  
  if( m_same_z_age->isVisible() != has_multiple_nucs_of_z )
    m_same_z_age->setHidden( !has_multiple_nucs_of_z );
  
  
  // We need Pu238, Pu239, and Pu240 for Bignan95_PWR and Bignan95_BWR.
  // We need Pu239 and at least one other isotope for ByPu239Only.
  const bool can_bignan = (nuc_names.count("Pu238") && nuc_names.count("Pu239")
                           && nuc_names.count("Pu240"));
  const bool can_239only = (nuc_names.count("Pu239")
                             && (nuc_names.count("Pu238") || nuc_names.count("Pu240")
                                  || nuc_names.count("Pu241") || nuc_names.count("Am241")) );
  
  const bool show_pu_corr = (can_bignan || can_239only);
  m_pu_corr_method->setHidden( !show_pu_corr );
  WLabel *combo_label = m_pu_corr_method->label();
  if( combo_label )
    combo_label->setHidden( !show_pu_corr );
  

  if( show_pu_corr )
  {
    const int nentries = m_pu_corr_method->count();
    const int nentries_needed = 1 + (can_239only ? 1 : 0) + (can_bignan ? 2 : 0);
    if( nentries != nentries_needed )
    {
      string current_txt;
      if( !m_pu_corr_method->count() )
        current_txt = m_pu_corr_method->currentText().toUTF8();
      
      m_pu_corr_method->clear();
      
      const string &none = RelActCalc::to_description( RelActCalc::PuCorrMethod::NotApplicable );
      const string &byPu239Only = RelActCalc::to_description( RelActCalc::PuCorrMethod::ByPu239Only );
      const string &bignanBwr = RelActCalc::to_description( RelActCalc::PuCorrMethod::Bignan95_BWR );
      const string &bignanPwr = RelActCalc::to_description( RelActCalc::PuCorrMethod::Bignan95_PWR );
      
      m_pu_corr_method->addItem( WString::fromUTF8(none) );
      int next_index = 0;
      if( can_239only )
      {
        if( current_txt == byPu239Only )
          next_index = m_pu_corr_method->count();
        m_pu_corr_method->addItem( WString::fromUTF8(byPu239Only) );
      }//if( can_239only )
      
      if( can_bignan )
      {
        if( current_txt == bignanPwr )
          next_index = m_pu_corr_method->count();
        m_pu_corr_method->addItem( WString::fromUTF8(bignanPwr) );
        
        if( current_txt == bignanBwr )
          next_index = m_pu_corr_method->count();
        m_pu_corr_method->addItem( WString::fromUTF8(bignanBwr) );
      }//if( can_bignan )
      
      m_pu_corr_method->setCurrentIndex( next_index );
    }//if( nentries != nentries_needed )
  }//if( show_pu_corr )
}//void updateDuringRenderForNuclideChange()


void RelActAutoGui::updateDuringRenderForRefGammaLineChange()
{
  // Determine if we should show/hide the reference lines.
  const bool show_ref_lines = m_hide_ref_lines_item->isEnabled();
  const vector<RelActCalcAuto::NucInputInfo> nuclides = getNucInputInfo();
  
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
      if( nuc.nuclide && !nuc.peak_color_css.empty() )
        nuclide_colors[nuc.nuclide->symbol] = WColor( nuc.peak_color_css );
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
      m_photopeak_widget->setIsotope( nuc.nuclide, nuc.age );
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
    
    if( options.rel_eff_curves.size() != 1 )
      throw runtime_error( "No relative efficiency curves defined - currently you must have exactly one curve." );

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
    if( m && m->detector() && m->detector()->isValid() && m->detector()->hasResolutionInfo() )
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
    const int num_index = (rel_acts[0].rel_activity > rel_acts[1].rel_activity) ? 1 : 0;
    const int denom_index = (num_index ? 0 : 1);
    const SandiaDecay::Nuclide * const num_nuc = rel_acts[num_index].nuclide;
    const SandiaDecay::Nuclide * const den_nuc = rel_acts[denom_index].nuclide;
    assert( num_nuc && den_nuc );
    if( num_nuc && den_nuc )
    {
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
      if( !nuc || !nuc->nuclide() )
        continue;
      
      if( fit_nuc.nuclide == nuc->nuclide() )
      {
        const string agestr = PhysicalUnitsLocalized::printToBestTimeUnits( fit_nuc.age, 3 );
        nuc->setAge( agestr );
        break;
      }//if( this is the widget for this nuclide )
    }//for( WWidget *w : kids )
  }//for( const RelActCalcAuto::NuclideRelAct &fit_nuc : m_solution->m_rel_activities )
  
  
  setOptionsForValidSolution();
  
  // Check if we need to update Physical model shieldings
  assert( m_solution->m_options.rel_eff_curves.size() == 1 );
  if( m_solution->m_options.rel_eff_curves.size() != 1 )
    throw runtime_error( "updateFromCalc: only a single rel eff curve is supported" );
  const RelActCalcAuto::RelEffCurveInput &rel_eff = m_solution->m_options.rel_eff_curves[0];
    
    
  if( (rel_eff.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel)
     && m_solution->m_phys_model_results.at(0).has_value() )
  {
    const RelActCalcAuto::RelActAutoSolution::PhysicalModelFitInfo &phys_info = *m_solution->m_phys_model_results.at(0);
    
    assert( m_phys_model_self_atten );
    if( !m_phys_model_self_atten )
      initPhysModelShields();
    assert( m_phys_model_self_atten );
    
    const auto update_shield_widget = []( const RelActCalcAuto::RelActAutoSolution::PhysicalModelFitInfo::ShieldInfo &shield,
                                      RelEffShieldWidget *w ) {
      if( !w )
        return;
      
      const Material * const widget_mat = w->material();
      assert( (!!shield.material) == (!!widget_mat) );
      assert( !widget_mat || !shield.material || (widget_mat->name == shield.material->name) );
      
      if( (!!shield.material) != (!!widget_mat) )
        w->setMaterialSelected( !!shield.material );
      
      if( shield.material && (!widget_mat || (widget_mat->name != shield.material->name)) )
        w->setMaterial( shield.material->name );
      
      if( shield.atomic_number )
      {
        assert( shield.atomic_number_was_fit == w->fitAtomicNumber() );
        assert( shield.atomic_number_was_fit || ((*shield.atomic_number) == w->atomicNumber()));
        if( (*shield.atomic_number) != w->atomicNumber() )
          w->setAtomicNumber( *shield.atomic_number );
        
        if( shield.atomic_number_was_fit != w->fitAtomicNumber() )
          w->setFitAtomicNumber( shield.atomic_number_was_fit );
      }//if( shield.atomic_number )
      
      
      if( shield.areal_density )
      {
        if( shield.material )
        {
          assert( shield.areal_density_was_fit == w->fitThickness() );
          const double thick = shield.areal_density / shield.material->density;
          
          if( !shield.areal_density_was_fit && (thick != w->thickness()) )
          {
            cerr << "Found fixed shielding thickness has changed for '" << shield.material->name << "'; was "
            << PhysicalUnits::printToBestLengthUnits(w->thickness()) << " but returned answer is "
            << PhysicalUnits::printToBestLengthUnits(thick) << endl;
            cerr << endl;
          }
          assert( shield.areal_density_was_fit || (thick == w->thickness()));
          
          w->setThickness( thick );
          if( shield.areal_density_was_fit != w->fitThickness() )
            w->setFitThickness( shield.areal_density_was_fit );
        }else
        {
          const double ad_g_cm2 = shield.areal_density / PhysicalUnits::g_per_cm2;
          
          assert( shield.areal_density_was_fit == w->fitArealDensity() );
          assert( shield.areal_density_was_fit || (ad_g_cm2 == w->arealDensity()));
          
          if( ad_g_cm2 != w->arealDensity() )
            w->setArealDensity( ad_g_cm2 );
          
          if( shield.areal_density_was_fit != w->fitArealDensity() )
            w->setFitArealDensity( shield.areal_density_was_fit );
        }//if( shield.material ) / else
      }//if( shield.areal_density )
    };//set_shield_widget lambda
   
    if( phys_info.self_atten.has_value() )
      update_shield_widget( *phys_info.self_atten, m_phys_model_self_atten );
    else if( m_phys_model_self_atten )
      m_phys_model_self_atten->resetState();
    
    vector<RelEffShieldWidget *> ext_shields;
    for( WWidget *w : m_phys_ext_attens->children() )
    {
      RelEffShieldWidget *sw = dynamic_cast<RelEffShieldWidget *>( w );
      if( sw )
        ext_shields.push_back( sw );
    }//for( WWidget *w : m_phys_ext_attens->children() )
    
    assert( ext_shields.size() >= phys_info.ext_shields.size() );
    while( ext_shields.size() < phys_info.ext_shields.size() )
    {
      RelEffShieldWidget *sw = new RelEffShieldWidget( RelEffShieldWidget::ShieldType::ExternalAtten, m_phys_ext_attens );
      sw->changed().connect( this, &RelActAutoGui::handlePhysModelShieldChange );
      ext_shields.push_back( sw );
    }
    
    while( (ext_shields.size() > 1) && (ext_shields.size() > phys_info.ext_shields.size()) )
    {
      delete ext_shields.back();
      ext_shields.pop_back();
    }
    
    if( (ext_shields.size() == 1) && phys_info.ext_shields.empty() )
      ext_shields[0]->resetState();
    
    assert( ext_shields.size() >= phys_info.ext_shields.size() );
    for( size_t i = 0; i < phys_info.ext_shields.size(); ++i )
      update_shield_widget( phys_info.ext_shields[i], ext_shields[i] );
   
    assert( phys_info.hoerl_b.has_value() == phys_info.hoerl_c.has_value() );
    assert( !m_phys_model_use_hoerl || (phys_info.hoerl_b.has_value() == m_phys_model_use_hoerl->isChecked()) );
    if( m_phys_model_use_hoerl && (m_phys_model_use_hoerl->isChecked() != phys_info.hoerl_b.has_value()) )
      m_phys_model_use_hoerl->setChecked( phys_info.hoerl_b.has_value() );
  }//if( m_solution->m_options.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  
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
