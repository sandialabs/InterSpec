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

#include <Wt/Utils>
#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WTable>
#include <Wt/WLineEdit>
#include <Wt/WTextArea>
#include <Wt/WApplication>
#include <Wt/WContainerWidget>

#include "InterSpec/SpecMeas.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/InterSpec.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/ShowRiidInstrumentsAna.h"

using namespace std;
using namespace Wt;

namespace
{

// Eventually may enable allowing to edit values, like when shown from the SpecFileSummary widget
//  Some work was already done to allow this, but putting behind this preprocessor variable for
//  the moment to keep things clearer
#define ENABLE_EDIT_RESULTS 0

//Class to display a DetectorAnalysis; just a first go at it
//  Some sort of model and table, or something might be a better implementation.
class AnaResultDisplay : public WContainerWidget
{
  WText *m_summary;
  std::shared_ptr<const SpecMeas> m_meas;
  
  WTable *m_table;
  
#if( ENABLE_EDIT_RESULTS )
  bool m_modifiable;
#endif
  
  enum AnaResultEditableFields
  {
    AlgorithmName,
    AlgorithmCreator,
    AlgorithmVersion,
    AlgorithmDescription,
    AlgorithmResultDescription,
    AlgorithmRemarks
  };//enum AnaResultEditableFields
  
#if( ENABLE_EDIT_RESULTS )
  WLineEdit *m_algorithm_name;
  WLineEdit *m_algorithm_version;
  WLineEdit *m_algorithm_creator;
  WLineEdit *m_algorithm_description;
  WLineEdit *m_algorithm_result_description;
  WTextArea *m_algorithm_remarks;
#else
  WText *m_algorithm_name;
  WText *m_algorithm_version;
  WText *m_algorithm_creator;
  WText *m_algorithm_description;
  WText *m_algorithm_result_description;
  WText *m_algorithm_remarks;
#endif
  
  
  void handleFieldUpdate( AnaResultEditableFields field )
  {
    switch( field )
    {
      case AlgorithmName:
      case AlgorithmVersion:
      case AlgorithmCreator:
      case AlgorithmDescription:
      case AlgorithmResultDescription:
      case AlgorithmRemarks:
        break;
    }//switch( field )
    
    InterSpec::instance()->logMessage( "Editing Analysis Results Not implemented - sorry", "", 0 );
  }//void handleFieldUpdate( AnaResultEditableFields field )
  
#if( ENABLE_EDIT_RESULTS )
  // For when we want to enable editing
  template <class T>
  void addField( T *&edit, WTable *table, const WString &labelstr,
                int row, int col, int rowspan = 1, int colspan = 1 )
  {
    WLabel *label = new WLabel(labelstr);
    label->setStyleClass("AlgoInfoLabel");
    WTableCell *cell = table->elementAt(row, col);
    cell->setRowSpan( rowspan );
    cell->addWidget( label );
    edit = new T();
    cell = table->elementAt(row, col+1);
    cell->setRowSpan( rowspan );
    cell->setColumnSpan( colspan );
    cell->addWidget( edit );
    label->setBuddy( edit );
    edit->disable();
    edit->setStyleClass( "AlgoInfoEdit" );
  }//WLineEdit *addField(...)
#endif
  
  void addField( WText *&edit, WTable *table, const WString &labelstr,
                int row, int col, int rowspan = 1, int colspan = 1 )
  {
    WLabel *label = new WLabel(labelstr);
    label->setStyleClass( "AlgoInfoLabel" );
    WTableCell *cell = table->elementAt(row, col);
    cell->setRowSpan( rowspan );
    cell->addWidget( label );
    
    cell = table->elementAt(row, col+1);
    cell->setRowSpan( rowspan );
    cell->setColumnSpan( colspan );
    edit = new WText();
    cell->addWidget( edit );
    edit->setStyleClass( "AlgoInfoTxt" );
  }
  
public:
  AnaResultDisplay( WContainerWidget *parent = nullptr )
  : WContainerWidget( parent ),
  m_summary( nullptr ),
  m_meas{ nullptr },
  m_table( nullptr ),
#if( ENABLE_EDIT_RESULTS )
  m_modifiable( false ),
#endif
  m_algorithm_name( nullptr ),
  m_algorithm_version( nullptr ),
  m_algorithm_creator( nullptr ),
  m_algorithm_description( nullptr ),
  m_algorithm_result_description( nullptr ),
  m_algorithm_remarks( nullptr )
  {
    wApp->useStyleSheet( "InterSpec_resources/ShowRiidInstrumentsAna.css" );
    
    addStyleClass( "AnaResultDisplay" );
    
    WText *algoInfo = new WText( "Algorithm Information:", this );
    algoInfo->addStyleClass( "AlgoInfoTitle" );
    algoInfo->setInline( false );
    
    m_table = new WTable( this );
    m_table->addStyleClass( "AlgoInfoTbl" );
    m_table->setHeaderCount( 1, Wt::Orientation::Vertical );
    
    addField( m_algorithm_name, m_table, "Name", AlgorithmName, 0 );
    addField( m_algorithm_creator, m_table, "Creator", AlgorithmCreator, 0 );
    addField( m_algorithm_version, m_table, "Version", AlgorithmVersion, 0 );
    addField( m_algorithm_description, m_table, "Desc.", AlgorithmDescription, 0 );
    addField( m_algorithm_result_description, m_table, "Desc.", AlgorithmResultDescription, 0 );
    addField( m_algorithm_remarks, m_table, "Remarks", AlgorithmRemarks, 0 );
    
#if( ENABLE_EDIT_RESULTS )
    m_algorithm_name->changed().connect( boost::bind( &AnaResultDisplay::handleFieldUpdate, this, AlgorithmName) );
    m_algorithm_version->changed().connect( boost::bind( &AnaResultDisplay::handleFieldUpdate, this,AlgorithmVersion ) );
    m_algorithm_creator->changed().connect( boost::bind( &AnaResultDisplay::handleFieldUpdate, this, AlgorithmCreator ) );
    m_algorithm_description->changed().connect( boost::bind( &AnaResultDisplay::handleFieldUpdate, this, AlgorithmDescription ) );
    m_algorithm_result_description->changed().connect( boost::bind( &AnaResultDisplay::handleFieldUpdate, this, AlgorithmResultDescription ) );
    m_algorithm_remarks->changed().connect( boost::bind( &AnaResultDisplay::handleFieldUpdate, this, AlgorithmRemarks ) );
#endif
  }//constructor
  
  void updateDisplay( std::shared_ptr<const SpecMeas> meas )
  {
    m_meas = meas;
    std::shared_ptr<const SpecUtils::DetectorAnalysis> ana;
    if( meas )
      ana = meas->detectors_analysis();
    
    if( m_summary )
      delete m_summary;
    m_summary = nullptr;
    
    if( !ana )
    {
      m_algorithm_name->setText( "" );
      m_algorithm_version->setText( "" );
      m_algorithm_creator->setText( "" );
      m_algorithm_description->setText( "" );
      m_algorithm_result_description->setText( "" );
      m_algorithm_remarks->setText( "" );
      return;
    }
    
    m_algorithm_name->setText( ana->algorithm_name_ );
    m_table->rowAt(AlgorithmName)->setHidden( m_algorithm_name->text().empty() );
    
    string algovrsn;
    for( const auto &nv : ana->algorithm_component_versions_ )
      algovrsn += (algovrsn.empty() ? "" : ", ") + nv.first + ": " + nv.second;
    
    m_algorithm_version->setText( algovrsn );
    m_table->rowAt(AlgorithmVersion)->setHidden( algovrsn.empty() );
    
    m_algorithm_creator->setText( ana->algorithm_creator_ );
    m_table->rowAt(AlgorithmCreator)->setHidden( m_algorithm_creator->text().empty() );
    
    m_algorithm_description->setText( ana->algorithm_description_ );
    m_table->rowAt(AlgorithmDescription)->setHidden( m_algorithm_description->text().empty() );
    
    m_algorithm_result_description->setText( ana->algorithm_result_description_ );
    m_table->rowAt(AlgorithmResultDescription)->setHidden( m_algorithm_result_description->text().empty() );
    
#if( ENABLE_EDIT_RESULTS )
    string remarktxt;
    for( string r : ana->remarks_ )
    {
      SpecUtils::ireplace_all( r, "&#009;", "  " );
      remarktxt += (remarktxt.size() ? "\n" : "") + r;
    }
    m_algorithm_remarks->setText( remarktxt );
#else
    string remarktxt;
    for( string r : ana->remarks_ )
    {
      if( r.empty() )
        continue;
      
      // Lets make lines line "Top Line 1: Am241&#009;662.45&#009;400.81" semi-reasonable, and have
      //  displayed all on the same line.
      //  TODO: there are probably a lot better ways to do this, like first check if there are any
      //    lines like this, and if so, make a html table.
      //  TODO: there are probably a lot of other special cases we can handle.
      SpecUtils::ireplace_all( r, "&#009;", "---mytab---" );
      r = Wt::Utils::htmlEncode( r, Wt::Utils::HtmlEncodingFlag::EncodeNewLines );
      SpecUtils::ireplace_all( r, "---mytab---", "&nbsp;&nbsp;&nbsp;" );
      SpecUtils::ireplace_all( r, ": ", ":&nbsp;" );
      SpecUtils::ireplace_all( r, "Top Line ", "Top&nbsp;Line&nbsp;" );
      
      remarktxt += "<div>" + r + "</div>";
    }//for( string r : ana->remarks_ )
    
    m_algorithm_remarks->setText( WString::fromUTF8(remarktxt) );
#endif
    
    m_table->rowAt(AlgorithmRemarks)->setHidden( m_algorithm_remarks->text().empty() );
    
    //Now Display the Nuclide results - hacking for now - should make dedicated
    //  widget
    string anastr;
    
    for( const SpecUtils::DetectorAnalysisResult &res : ana->results_ )
    {
      string result;
      if( res.nuclide_.size() )
        result = "<tr><th>Nuclide</th><td>" + res.nuclide_ + "</td></tr>";
      if( res.nuclide_type_.size() )
        result += "<tr><th>Category</th><td>" + res.nuclide_type_ + "</td></tr>";
      if( res.id_confidence_.size() )
        result += "<tr><th>Confidence</th><td>" + res.id_confidence_ + "</td></tr>";
      if( res.detector_.size() )
        result += "<tr><th>Detector</th><td>" + res.detector_ + "</td></tr>";
      
      if( res.dose_rate_ > 0.0 )
        result += "<tr><th>Dose</th><td>"
        + PhysicalUnits::printToBestEquivalentDoseRateUnits( 1.0E-6 * res.dose_rate_*PhysicalUnits::sievert/PhysicalUnits::hour ) + "</td></tr>";
      if( res.distance_ > 0.0 )
        result += "<tr><th>Distance</th><td>" + PhysicalUnits::printToBestLengthUnits(0.1*res.distance_) + "</td></tr>";
      if( res.activity_ > 0.0 )
      {
        const bool useBq = InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
        result += "<tr><th>Activity</th><td>" + PhysicalUnits::printToBestActivityUnits( res.activity_, 2, !useBq, 1.0 ) + "</td></tr>";
      }
      if( res.remark_.size() > 0 )
        result += "<tr><th>Remark</th><td>" + res.remark_ + "</td></tr>";
      
      //  float real_time_;           //in units of seconds (eg: 1.0 = 1 s)
      //  boost::posix_time::ptime start_time_;
      
      if( !result.empty() )
        anastr += "<table class=\"RiidNuclideTable\">" + result + "</table>";
    }//for( const SpecUtils::DetectorAnalysisResult &res : ana->results_ )
    
    if( ana->results_.empty() )
      anastr = "<div class=\"RiidNoResultsTxt\">No nuclide identifications provided.</div>";
    
    if( anastr.size() > 2 )
    {
      anastr = "<div class=\"RiidNuclideResults\"><div class=\"AlgoInfoTitle\">Algorithm Results:</div>" + anastr + "</div>";
      m_summary = new WText( anastr );
      m_summary->setInline( false );
      this->addWidget( m_summary );
    }else
    {
      m_summary = new WText( "No algorithm results available." );
      m_summary->addStyleClass( "AlgoInfoTitle" );
      m_summary->setInline( false );
    }
  }//void updateDisplay( std::shared_ptr<const SpecMeas> meas )

  /*
  void allowModifiable( bool allow )
  {
    //unimplemented
    m_modifiable = allow;
    
    m_algorithm_name->setEnabled( allow );
    m_algorithm_version->setEnabled( allow );
    m_algorithm_creator->setEnabled( allow );
    m_algorithm_description->setEnabled( allow );
    m_algorithm_result_description->setEnabled( allow );
    m_algorithm_remarks->setEnabled( allow );
    
    m_table->rowAt(AlgorithmName)->setHidden( !allow && m_algorithm_name->text().empty() );
    m_table->rowAt(AlgorithmVersion)->setHidden( !allow && m_algorithm_version->text().empty() );
    m_table->rowAt(AlgorithmCreator)->setHidden( !allow && m_algorithm_creator->text().empty() );
    m_table->rowAt(AlgorithmDescription)->setHidden( !allow && m_algorithm_description->text().empty() );
    m_table->rowAt(AlgorithmResultDescription)->setHidden( !allow && m_algorithm_result_description->text().empty() );
  }//void allowModifiable( bool allow )
  */
};//class AnaResultDisplay
}//namespace


std::string riidAnaSummary( const std::shared_ptr<const SpecMeas> &spec )
{
  string summary;
  
  shared_ptr<const SpecUtils::DetectorAnalysis> ana = spec ? spec->detectors_analysis() : nullptr;
  if( !ana )
    return summary;

  
  for( const SpecUtils::DetectorAnalysisResult &res : ana->results_ )
  {
    if( !res.isEmpty() && !res.nuclide_.empty() )
      summary += (summary.empty() ? "" : ", ") + res.nuclide_;
  }
  
  if( summary.length() > 64 )
    summary = summary.substr(0,61) + "...";
  
  if( summary.empty() && ana->results_.empty() )
    summary = "no nuclides identified.";
  
  return Wt::Utils::htmlEncode( WString::fromUTF8(summary),0).toUTF8();
}//riidAnaSummary(...)


void showRiidInstrumentsAna( const std::shared_ptr<const SpecMeas> &spec )
{
  auto dialog = new SimpleDialog();
  //dialog->setModal( false ); //doesnt seem to have any effect
  dialog->addButton( "Close" );
  
  WContainerWidget *contents = dialog->contents();
  WText *dialogTitle = new WText( "The Detectors RIID Results", contents );
  dialogTitle->addStyleClass( "title RiidDialogTitle" );
  dialogTitle->setInline( false );
  
  AnaResultDisplay *display = new AnaResultDisplay( contents );
  display->updateDisplay( spec );
}//void showRiidInstrumentsAna( const std::shared_ptr<const SpecMeas> &spec )
