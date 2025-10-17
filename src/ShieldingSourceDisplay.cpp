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

#include <memory>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <iostream>

#include <boost/scope_exit.hpp>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

//Roots Minuit2 includes
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"


#include <Wt/WText>
#include <Wt/WTime>
#include <Wt/WImage>
#include <Wt/WTable>
#include <Wt/WLabel>
#include <Wt/WAnchor>
#include <Wt/WBorder>
#include <Wt/WServer>
#include <Wt/WPainter>
#include <Wt/WGroupBox>
#include <Wt/WTextArea>
#include <Wt/WRectArea>
#include <Wt/WResource>
#include <Wt/WSvgImage>
#include <Wt/WTableCell>
#include <Wt/WTabWidget>
#include <Wt/WIOService>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WJavaScript>
#include <Wt/WFileUpload>
#include <Wt/WEnvironment>
#include <Wt/Http/Response>
#include <Wt/WSelectionBox>
#include <Wt/WStandardItem>
#include <Wt/WItemDelegate>
#include <Wt/WDoubleSpinBox>
#include <Wt/Dbo/QueryModel>
#include <Wt/WMemoryResource>
#include <Wt/WSuggestionPopup>
#include <Wt/WRegExpValidator>
#include <Wt/WDoubleValidator>
#include <Wt/WStandardItemModel>
#if( WT_VERSION >= 0x3030600 )
#include <Wt/Chart/WStandardChartProxyModel>
#endif
#include <Wt/WCssDecorationStyle>
#include <Wt/Chart/WChartPalette>
#include <Wt/Chart/WCartesianChart>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/SpecUtilsAsync.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/SwitchCheckbox.h"
#include "InterSpec/ShieldingSelect.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/IsotopeSelectionAids.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/PhysicalUnitsLocalized.h"
#include "InterSpec/ShieldingSourceDisplay.h"

using namespace Wt;
using namespace std;


#if( ANDROID )
// Defined in target/android/android.cpp
extern void android_download_workaround( Wt::WResource *resource, std::string description );
#endif


typedef std::shared_ptr<const PeakDef> PeakShrdPtr;

const int ShieldingSourceDisplay::sm_xmlSerializationMajorVersion = 0;

/** Change log:
 - 20211031, Version 1: added "SourceType" node under the <Nuclide> to say whether its a point, intrinsic, or trace source.
                  added "Geometry" node under the "ShieldingSourceFit" node
 - 20230712, no version change: split out to SourceFitDef::sm_xmlSerializationMinorVersion.
 */
const int ShieldingSourceDisplay::sm_xmlSerializationMinorVersion = 1;


using GammaInteractionCalc::GeometryType;
using GammaInteractionCalc::TraceActivityType;

namespace
{
  const std::string ns_no_uncert_info_txt = "Perform model fit to update and get uncertainties.";
  
  /** If a distance WLineEdit has a number, but no distance units, will add a " cm" to the text value. */
  void make_sure_distance_units_present( Wt::WLineEdit *edit )
  {
    if( !edit )
      return;
      
    string diststr = edit->text().toUTF8();
    SpecUtils::trim( diststr );
    
    if( diststr.empty() )
      return;
    
    if( diststr.find_first_not_of( " \t0123456789.eE+-\n" ) == string::npos )
    {
      diststr += " cm";
      edit->setText( diststr );
    }
  }//void make_sure_distance_units_present( Wt::WLineEdit *edit )
}//namespace


WT_DECLARE_WT_MEMBER
(ShowChi2Info, Wt::JavaScriptFunction, "ShowChi2Info",
function(id,info)
{
  $("<div id=\"" + id + "inf\" class=\"ChartMouseOverInfo\"></div>").html(info).appendTo($('#'+id));
}
);


namespace
{
  class StringDownloadResource : public Wt::WResource
  {
    ShieldingSourceDisplay *m_display;
    Wt::WApplication *m_app;
    
  public:
    StringDownloadResource( ShieldingSourceDisplay *parent )
      : WResource( parent ), m_display( parent ), m_app( WApplication::instance() )
    {
      assert( m_app );
    }
    
    virtual ~StringDownloadResource()
    {
      beingDeleted();
    }
    
    virtual void handleRequest( const Wt::Http::Request &request,
                               Wt::Http::Response &response )
    {
      WApplication::UpdateLock lock( m_app );
      
      if( !lock )
      {
        log("error") << "Failed to WApplication::UpdateLock in StringDownloadResource.";
        response.out() << "Error grabbing application lock to form StringDownloadResource resource; please report to InterSpec@sandia.gov.";
        response.setStatus(500);
        assert( 0 );
        return;
      }//if( !lock )
      
      suggestFileName( "shielding_source_fit_model.xml", WResource::Attachment );
      response.setMimeType( "application/xml" );
      if( !m_display )
        return;
      rapidxml::xml_document<char> doc;
      m_display->serialize( &doc );
      rapidxml::print( response.out(), doc, 0 );
    }
  };//class PeakCsvResource : public Wt::WResource
  
  
  /** Struct that saves the ShieldingSourceDisplay state to XML, when this struct is first constructed, and then again
   when the struct destructs; it then uses these two XML states to create a undo/redo point.
   If blocks all other undo/redo step insertions until this object is destructed.
   
   You are looking at at least about 15 kb memory for each undo/redo step, just for the XML; for small things, it may be more
   efficient to not use this mechanism.
   */
  struct ShieldSourceChange
  {
    ShieldingSourceDisplay *m_display;
    string m_description;
    shared_ptr<const string> m_pre_doc;
    
    // Make sure we dont insert multiple undo/redo steps, so using a BlockUndoRedoInserts.
    unique_ptr<UndoRedoManager::BlockUndoRedoInserts> m_block;
    
    ShieldSourceChange( ShieldingSourceDisplay *display, const string &descrip )
      : m_display( display ),
        m_description( descrip ),
        m_pre_doc( nullptr ),
        m_block( nullptr )
    {
      UndoRedoManager *undoRedo = UndoRedoManager::instance();
      if( undoRedo && undoRedo->canAddUndoRedoNow() )
      {
        m_block = make_unique<UndoRedoManager::BlockUndoRedoInserts>();
        shared_ptr<string> doc_xml = make_shared<string>();
        doSerialization( *doc_xml );
        m_pre_doc = doc_xml;
      }
    }//ShieldSourceChange
    
    
    void doSerialization( std::string &xmldoc )
    {
      try
      {
        rapidxml::xml_document<char> doc;
        m_display->serialize( &doc );
        rapidxml::print( std::back_inserter(xmldoc), doc, rapidxml::print_no_indenting );
      }catch( std::exception &e )
      {
        passMessage( "Error adding undo step for " + m_description + ": "
                    + string(e.what()), WarningWidget::WarningMsgHigh );
        xmldoc.clear(); //jic
      }
    }//void doSerialization( std::string &xmldoc )
    
    
    static void doDeSerialization( ShieldingSourceDisplay *display,
                                   const shared_ptr<const string> xml_doc )
    {
      if( !display || !xml_doc || xml_doc->empty() )
        throw logic_error( "Error with input logic" );
      
      display->cancelModelFitWithNoUpdate();
      
      rapidxml::xml_document<char> new_doc;
      const int flags = rapidxml::parse_normalize_whitespace | rapidxml::parse_trim_whitespace;
      string pre_doc_cpy = *xml_doc;
      new_doc.parse<flags>( &(pre_doc_cpy[0]) );
      display->deSerialize( new_doc.first_node() );
    }//static void doDeSerialization( ShieldingSourceDisplay *display, const shared_ptr<const string> xml_doc )
    
    
    ~ShieldSourceChange()
    {
      UndoRedoManager *undoRedo = UndoRedoManager::instance();
      if( !m_display || !m_pre_doc || m_pre_doc->empty() || !undoRedo )
        return;
      
      shared_ptr<string> post_doc_xml = make_shared<string>();
      try
      {
        rapidxml::xml_document<char> doc;
        m_display->serialize( &doc );
        rapidxml::print( std::back_inserter(*post_doc_xml), doc, rapidxml::print_no_indenting );
      }catch( std::exception &e )
      {
        passMessage( "Error adding undo/redo step for " + m_description + ": "
                    + string(e.what()), WarningWidget::WarningMsgHigh );
        return;
      }
      
      const string descrip = std::move( m_description );
      const shared_ptr<const string> post_doc = post_doc_xml;
      const shared_ptr<const string> pre_doc = m_pre_doc;
      
      if( (*post_doc) == (*pre_doc) )
      {
        Wt::log("debug") << "ShieldSourceChange: no changes found; not inserting undo/redo step";
        return;
      }
      
      // Now need to remove the block to inserting undo/redo steps
      m_block.reset();
      
      //cout << "Adding undo/redo step will take approx "
      //     << (m_description.size() + pre_doc->size() + post_doc->size())
      //     << " bytes of memory" << endl;
      
      auto undo = [pre_doc, descrip](){
        ShieldingSourceDisplay *display = InterSpec::instance()->shieldingSourceFit();
        if( !display || !pre_doc || pre_doc->empty() )
          return;
        
        try
        {
          doDeSerialization( display, pre_doc );
        }catch( std::exception &e )
        {
          passMessage( "Error undoing " + descrip + ": " + string(e.what()),
                      WarningWidget::WarningMsgHigh );
        }
      };
      
      auto redo = [post_doc, descrip](){
        ShieldingSourceDisplay *display = InterSpec::instance()->shieldingSourceFit();
        if( !display || !post_doc || post_doc->empty() )
          return;
        
        try
        {
          doDeSerialization( display, post_doc );
        }catch( std::exception &e )
        {
          passMessage( "Error redoing" + descrip + ": " + string(e.what()),
                      WarningWidget::WarningMsgHigh );
        }
      };
      
      undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), descrip );
    }//~ShieldSourceChange()
    
    ShieldSourceChange( const ShieldSourceChange & ) = delete; // non construction-copyable
    ShieldSourceChange &operator=( const ShieldSourceChange & ) = delete; // non copyable
  };//struct PeakModelChange
}//namespace



SourceFitModel::SourceFitModel( PeakModel *peakModel,
                                const bool sameAgeIsotopes,
                                Wt::WObject *parent )
  : WAbstractItemModel( parent ),
    m_sortOrder( Wt::AscendingOrder ),
    m_sortColumn( kIsotope ),
    m_peakModel( peakModel ),
    m_sameAgeForIsotopes( sameAgeIsotopes ),
    m_det_type( DetectorPeakResponse::EffGeometryType::FarField )
{
  auto interspec = InterSpec::instance();
  if( !interspec )
  {
    m_displayCuries = true;
  }else
  {
    interspec->useMessageResourceBundle( "ShieldingSourceDisplay" ); //jic
    m_displayCuries = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", interspec );
    interspec->preferences()->addCallbackWhenChanged( "DisplayBecquerel",
                                          this, &SourceFitModel::displayUnitsChanged );
  }//if( !interspec ) / else
  
  peakModel->rowsAboutToBeRemoved().connect( this, &SourceFitModel::peakModelRowsRemovedCallback );
  peakModel->rowsInserted().connect( this, &SourceFitModel::peakModelRowsInsertedCallback );
  peakModel->dataChanged().connect( this, &SourceFitModel::peakModelDataChangedCallback );
  peakModel->modelReset().connect( this, &SourceFitModel::peakModelResetCallback );
  repopulateIsotopes();
}//SourceFitModel(...)


SourceFitModel::~SourceFitModel()
{
  //nothing to do here
}//~SourceFitModel()


void SourceFitModel::displayUnitsChanged( bool useBq )
{
  //cout << "in SourceFitModel::displayUnitsChanged" << endl;
  try
  {
    if( useBq == m_displayCuries )
    {
      m_displayCuries = !useBq;
      const int nrow = rowCount();
      if( nrow > 0 )
        dataChanged().emit( createIndex(0,0,nullptr), createIndex(nrow-1,kNumColumns-1,nullptr) );
    }
  }catch( std::exception &e )
  {
    //Shouldnt ever happen, but print a little something out JIC
    cerr << "SourceFitModel::displayUnitsChanged: Failed to convert boost any: " << e.what() << endl;
  }
  
  //cout << "m_displayCuries is now: " << m_displayCuries << endl;
}//void SourceFitModel::displayUnitsChanged( boost::any value )


const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &SourceFitModel::underlyingData() const
{
  return m_nuclides;
}


void SourceFitModel::setUnderlyingData( const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &input_data )
{
  auto blocker = make_shared<UndoRedoManager::BlockUndoRedoInserts>();
  
  const vector<ShieldingSourceFitCalc::IsoFitStruct> orig_rows = m_nuclides;
  
  vector<ShieldingSourceFitCalc::IsoFitStruct> data;
  shared_ptr<const deque<PeakModel::PeakShrdPtr>> peaks = m_peakModel->peaks();
  const size_t npeaks = peaks ? peaks->size() : size_t(0);
  for( ShieldingSourceFitCalc::IsoFitStruct in : input_data )
  {
    bool hasNuc = false;
    set<const SandiaDecay::Nuclide *> progeny;
    for( size_t i = 0; i < npeaks; ++i )
    {
      const auto &p = (*peaks)[i];
      if( p->useForShieldingSourceFit() && (in.nuclide == p->parentNuclide()) )
      {
        hasNuc = true;
        const SandiaDecay::Transition * const trans = p->nuclearTransition();
        if( trans && trans->parent )
          progeny.insert( trans->parent );
      }//
    }//for( size_t i = 0; i < npeaks; ++i )
    
    if( hasNuc )
    {
      in.numProgenyPeaksSelected = progeny.size();
      data.push_back( in );
    }
  }//for( const ShieldingSourceFitCalc::IsoFitStruct &in : input_data )
  
  std::sort( begin(data), end(data),
             boost::bind( &SourceFitModel::compare, boost::placeholders::_1,
                         boost::placeholders::_2, m_sortColumn, m_sortOrder ) );
  
  
  if( data.empty() && m_nuclides.empty() )
    return;
  
  const bool changingNumRows = (data.size() != m_nuclides.size());
  
  if( changingNumRows && m_nuclides.size() )
  {
    const int lastrow = static_cast<int>(m_nuclides.size() - 1);
    beginRemoveRows( WModelIndex(), 0, lastrow );
    m_nuclides.clear();
    endRemoveRows();
  }//if( changingNumRows && m_nuclides.size() )
  
  if( data.empty() )
    return;
  
  if( changingNumRows )
  {
    beginInsertRows( WModelIndex(), 0, static_cast<int>(data.size() - 1) );
    m_nuclides = data;
    endInsertRows();
  }else
  {
    m_nuclides = data;
    const WModelIndex topLeft = index( 0, 0 );
    const WModelIndex bottomRight = index( static_cast<int>(data.size() - 1), columnCount() - 1 );
    dataChanged().emit( topLeft, bottomRight );
  }
  
  blocker.reset();
  
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    auto undo = [orig_rows](){
      InterSpec *viewer = InterSpec::instance();
      ShieldingSourceDisplay *srcfit = viewer ? viewer->shieldingSourceFit() : nullptr;
      SourceFitModel *srcmodel = srcfit ? srcfit->sourceFitModel() : nullptr;
      if( srcmodel )
        srcmodel->setUnderlyingData( orig_rows );
    };
    
    auto redo = [data](){
      InterSpec *viewer = InterSpec::instance();
      ShieldingSourceDisplay *srcfit = viewer ? viewer->shieldingSourceFit() : nullptr;
      SourceFitModel *srcmodel = srcfit ? srcfit->sourceFitModel() : nullptr;
      if( srcmodel )
        srcmodel->setUnderlyingData( data );
    };
    
    undoRedo->addUndoRedoStep( undo, redo, "Update Shielding/Source nuclide info" );
  }//if( undoRedo && !undoRedo->isInUndoOrRedo() )
}//void setUnderlyingData( std::vector<ShieldingSourceFitCalc::IsoFitStruct> &data )


void SourceFitModel::setDetectorType( const DetectorPeakResponse::EffGeometryType det_type )
{
  if( m_det_type == det_type )
    return;
  
  m_det_type = det_type;
  headerDataChanged().emit( Orientation::Horizontal,
                           static_cast<int>(kActivity), static_cast<int>(kActivity) );
  
  if( !m_nuclides.empty() )
  {
    int col = static_cast<int>(kActivity);
    const WModelIndex act_first = index( 0, col );
    const WModelIndex act_second = index( static_cast<int>(m_nuclides.size()-1), col );
    dataChanged().emit(act_first, act_second);
    
    col = static_cast<int>(kActivityUncertainty);
    const WModelIndex uncert_first = index( 0, col );
    const WModelIndex uncert_second = index( static_cast<int>(m_nuclides.size()-1), col );
    dataChanged().emit(uncert_first, uncert_second);
    
    switch( det_type )
    {
      case DetectorPeakResponse::EffGeometryType::FarField:
        break;
        
      case DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct:
      case DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2:
      case DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2:
      case DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram:
      {
        for( const ShieldingSourceFitCalc::IsoFitStruct &iso : m_nuclides )
        {
          switch( iso.sourceType )
          {
            case ShieldingSourceFitCalc::ModelSourceType::Point:
              break;
              
            case ShieldingSourceFitCalc::ModelSourceType::Intrinsic:
            case ShieldingSourceFitCalc::ModelSourceType::Trace:
              setSourceType( iso.nuclide, ShieldingSourceFitCalc::ModelSourceType::Point );
              break;
          }//switch( iso.sourceType )
        }//for( loop over m_nuclides )
        
        break;
      }//case Fixed Geometry
    }//switch( det_type )
  }//if( !m_nuclides.empty() )
}//void setDetectorType( const int detector_EffGeometryType )


DetectorPeakResponse::EffGeometryType SourceFitModel::detType() const
{
  return m_det_type;
}

int SourceFitModel::numNuclides() const
{
  return static_cast<int>( m_nuclides.size() );
}//int numNuclides() const


int SourceFitModel::row( const SandiaDecay::Nuclide *nuclide ) const
{
  if( !nuclide )
    return -1;

  const int nrow = static_cast<int>( m_nuclides.size() );
  for( int i = 0; i < nrow; ++i )
    if( m_nuclides[i].nuclide == nuclide )
      return i;
  return -1;
}//int row( const SandiaDecay::Nuclide *nuclide ) const;


int SourceFitModel::nuclideIndex( const SandiaDecay::Nuclide *nuc ) const
{
  for( size_t i = 0; i < m_nuclides.size(); ++i )
    if( m_nuclides[i].nuclide == nuc )
      return static_cast<int>(i);
  throw std::runtime_error( "SourceFitModel::nuclideIndex with invalid nuclide (" + (nuc ? nuc->symbol : string("null")) + ")" );
  return -1;
}//int nuclideNumber( const SandiaDecay::Nuclide *nuc ) const


const SandiaDecay::Nuclide *SourceFitModel::nuclide( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].nuclide;
}//const SandiaDecay::Nuclide *nuclide( int nuc ) const


double SourceFitModel::activity( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].activity;
}//double activity( int nuc ) const

double SourceFitModel::activityUncert( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].activityUncertainty;
}//double activityUncert( int nuc ) const


bool SourceFitModel::fitActivity( int nuc ) const
{
  if( (nuc < 0) || (nuc >= static_cast<int>(m_nuclides.size())) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  
  switch( m_nuclides[nuc].sourceType )
  {
    case ShieldingSourceFitCalc::ModelSourceType::Point:
    case ShieldingSourceFitCalc::ModelSourceType::Trace:
      return m_nuclides[nuc].fitActivity;
      
    case ShieldingSourceFitCalc::ModelSourceType::Intrinsic:
      return false;
  }//switch( m_nuclides[nuc].sourceType )
  
  assert( 0 );
  return false;
}//bool fitActivity( int nuc ) const


double SourceFitModel::age( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].age;
}//double age( int nuc ) const


double SourceFitModel::ageUncert( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].ageUncertainty;
}//double ageUncert( int nuc ) const


bool SourceFitModel::fitAge( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  
  const auto &nucibj = m_nuclides[nuc];
  
  return ((nucibj.numProgenyPeaksSelected > 1 || nucibj.ageDefiningNuc) && nucibj.fitAge);
}//bool fitAge( int nuc ) const

#if( INCLUDE_ANALYSIS_TEST_SUITE )
boost::optional<double> SourceFitModel::truthActivity( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].truthActivity;
}

boost::optional<double> SourceFitModel::truthActivityTolerance( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].truthActivityTolerance;
}

boost::optional<double> SourceFitModel::truthAge( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].truthAge;
}

boost::optional<double> SourceFitModel::truthAgeTolerance( int nuc ) const
{
  if( (nuc < 0) || (nuc >= static_cast<int>(m_nuclides.size())) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].truthAgeTolerance;
}
#endif  //#if( INCLUDE_ANALYSIS_TEST_SUITE )


const SandiaDecay::Nuclide *SourceFitModel::ageDefiningNuclide(
                                  const SandiaDecay::Nuclide *dependantNuc ) const
{
  const int nuc = nuclideIndex( dependantNuc );
  if( (nuc < 0) || (nuc >= static_cast<int>(m_nuclides.size())) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  
  if( m_nuclides[nuc].ageDefiningNuc )
    return m_nuclides[nuc].ageDefiningNuc;
  return m_nuclides[nuc].nuclide;
}//bool ageDefiningNuclide( int nuc ) const


ShieldingSourceFitCalc::ModelSourceType SourceFitModel::sourceType( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  
  return m_nuclides[nuc].sourceType;
}//ModelSourceType sourceType( int nuc ) const



bool SourceFitModel::isVolumetricSource( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  
  switch( m_nuclides[nuc].sourceType )
  {
    case ShieldingSourceFitCalc::ModelSourceType::Point:
      return false;
      
    case ShieldingSourceFitCalc::ModelSourceType::Intrinsic:
    case ShieldingSourceFitCalc::ModelSourceType::Trace:
      return true;
  }//switch( m_nuclides[nuc].sourceType )
  
  assert( 0 );
  return false;
}//bool isVolumetricSource(...)


void SourceFitModel::setSharedAgeNuclide( const SandiaDecay::Nuclide *dependantNuc,
                                           const SandiaDecay::Nuclide *definingNuc )
{
  const int row = nuclideIndex( dependantNuc );
  if( row < 0 )
    return;
  
  if( PeakDef::ageFitNotAllowed( dependantNuc ) )
    definingNuc = nullptr;
  
  if( definingNuc )
  {
    const int definingRow = nuclideIndex( definingNuc );
    if( definingRow < 0 )
      throw runtime_error( "SourceFitModel::setSharedAgeNuclide: defining"
                           " nuclide must also be in the model" );
  }//if( definingNuc )
 
  if( definingNuc == dependantNuc )
    definingNuc = nullptr;
  
  ShieldingSourceFitCalc::IsoFitStruct &iso = m_nuclides[row];
  if( iso.ageDefiningNuc == definingNuc )
    return;
  
  if( definingNuc && dependantNuc && (definingNuc->atomicNumber != dependantNuc->atomicNumber) )
    throw runtime_error( "SourceFitModel::setSharedAgeNuclide: dependant and"
                         " defining nuclides must have same atomic number" );
  
  iso.ageDefiningNuc = definingNuc;
  dataChanged().emit( createIndex(row,kFitAge,(void *)0),
                      createIndex(row,kAgeUncertainty,(void *)0) );
}//void makeAgeFitable( const SandiaDecay::Nuclide *nuc, bool fit )


void SourceFitModel::setSourceType( const SandiaDecay::Nuclide *nuc, ShieldingSourceFitCalc::ModelSourceType type )
{
  const int nrows = static_cast<int>( m_nuclides.size() );

  for( int row = 0; row < nrows; ++row )
  {
    ShieldingSourceFitCalc::IsoFitStruct &iso = m_nuclides[row];

    if( iso.nuclide && (iso.nuclide==nuc) )
    {
      if( iso.sourceType != type )
      {
        iso.sourceType = type;
        dataChanged().emit( createIndex(row,0,nullptr), createIndex(row,columnCount()-1,nullptr) );
        return;
      }
    }//if( this is the nuclide we want )
  }//for( const ShieldingSourceFitCalc::IsoFitStruct &iso : m_nuclides )
}//void setSourceType( const SandiaDecay::Nuclide *nuc, ShieldingSourceFitCalc::ModelSourceType type )


void SourceFitModel::setUseSameAgeForIsotopes( bool useSame )
{
  if( m_sameAgeForIsotopes == useSame )
    return;

  m_sameAgeForIsotopes = useSame;
  
  if( !m_sameAgeForIsotopes )
  {
    for( size_t i = 0; i < m_nuclides.size(); ++i )
      setSharedAgeNuclide( m_nuclides[i].nuclide, nullptr );
    return;
  }//if( !m_sameAgeForIsotopes )
  
  //construct map from element number to isotopes being fit for
  //if there is more than one isotope, choose the youngest age to assign to
  //all of the isotopes.  Do the assignment and update the chart
  typedef map< int, vector<const SandiaDecay::Nuclide *> > ElToNucMap_t;
  ElToNucMap_t elToNucMap;
  
  for( const ShieldingSourceFitCalc::IsoFitStruct &n : m_nuclides )
  {
    if( n.nuclide )
      elToNucMap[n.nuclide->atomicNumber].push_back( n.nuclide );
  }
  
  for( const ElToNucMap_t::value_type &vt : elToNucMap )
  {
    const vector<const SandiaDecay::Nuclide *> &nucs = vt.second;
    if( nucs.size() < 2 )
    {
      if( !nucs.empty() )
        setSharedAgeNuclide( nucs[0], nullptr );
      continue;
    }
    
    bool fitAgeWanted = false;
    size_t minIndex = 0, numSharingAge = 0;
    double minAge = std::numeric_limits<double>::infinity();
    for( size_t i = 0; i < nucs.size(); ++i )
    {
      const SandiaDecay::Nuclide *nuc = nucs[i];
      if( PeakDef::ageFitNotAllowed( nuc ) )
      {
        setSharedAgeNuclide( nuc, nullptr ); // jic, I guess, but not really needed, probably
        continue;
      }
        
      const int ind = nuclideIndex( nuc );
      const double thisage = age( ind );
      if( thisage < minAge )
      {
        minIndex = i;
        minAge = thisage;
      }//if( age < minAge )
      
      numSharingAge += 1;
      fitAgeWanted = (fitAgeWanted || fitAge(ind));
    }//for( const SandiaDecay::Nuclide *nuc : nucs )
    
    for( size_t i = 0; i < nucs.size(); ++i )
    {
      const SandiaDecay::Nuclide *nuc = nucs[i];
      if( PeakDef::ageFitNotAllowed( nuc ) )
        continue;
        
      if( numSharingAge < 2 )
      {
        setSharedAgeNuclide( nuc, nullptr );
      }else if( i == minIndex )
      {
        WModelIndex ind = index( nuc, SourceFitModel::kAge );
        setData( ind, fitAgeWanted );
        ind = index( nuc, SourceFitModel::kAgeUncertainty );
        setData( ind, boost::any() );
        setSharedAgeNuclide( nuc, NULL );
      }else
      {
        setSharedAgeNuclide( nuc, nucs[minIndex] );
      }//if( i == 0 )
    }//for( const SandiaDecay::Nuclide *nuc : nucs )
  }//for( const ElToNucMap_t::value_type &vt : elToNucMap )
}//void setUseSameAgeForIsotopes( bool useSame )


void SourceFitModel::insertPeak( const PeakShrdPtr peak )
{
  if( !peak )
  {  //probably wont ever happen
    cerr << "SourceFitModel::insertPeak(...)\n\tNot a valid row being added" << endl;
    return;
  }//if( !peak )

  const SandiaDecay::Nuclide *nuclide = peak->parentNuclide();
  if( !nuclide )
    return;
  
  if( !peak->useForShieldingSourceFit() )
    return;

  for( const ShieldingSourceFitCalc::IsoFitStruct &iso : m_nuclides )
    if( iso.nuclide == peak->parentNuclide() )
      return;

  ShieldingSourceFitCalc::IsoFitStruct newIso;
  newIso.activity = 0.001*PhysicalUnits::curie;
  newIso.fitActivity = true;
  newIso.nuclide = peak->parentNuclide();
  newIso.age = PeakDef::defaultDecayTime( newIso.nuclide );
  newIso.ageDefiningNuc = nullptr;
  newIso.fitAge = false;
  newIso.ageIsFittable = !PeakDef::ageFitNotAllowed( newIso.nuclide );
  
  std::map<const SandiaDecay::Nuclide *, ShieldingSourceFitCalc::IsoFitStruct>::iterator oldval
                                    = m_previousResults.find( newIso.nuclide );
  if( oldval != m_previousResults.end() )
  {
    const ShieldingSourceFitCalc::IsoFitStruct &oldIso = oldval->second;
    newIso.activity = oldIso.activity;
    newIso.fitActivity = oldIso.fitActivity;
    newIso.age = oldIso.age;
    newIso.fitAge = oldIso.fitAge;
    newIso.ageIsFittable = oldIso.ageIsFittable;
    
    m_previousResults.erase( oldval );
  }//if( m_previousResults.count(newIso.nuclide) )
  
  if( m_previousResults.size() > 100 )
    m_previousResults.clear();
  
  if( m_sameAgeForIsotopes && newIso.ageIsFittable )
  {
    vector<size_t> thisElementIndexs;
    for( size_t i = 0; i < m_nuclides.size(); ++i )
    {
      const ShieldingSourceFitCalc::IsoFitStruct &iso = m_nuclides[i];
      if( iso.ageIsFittable && (iso.nuclide->atomicNumber == newIso.nuclide->atomicNumber) )
      {
        thisElementIndexs.push_back( i );
        if( !iso.ageDefiningNuc )
          newIso.ageDefiningNuc = iso.nuclide;
      }
    }//for( const ShieldingSourceFitCalc::IsoFitStruct &iso : m_nuclides )
    
    if( !newIso.ageDefiningNuc && !thisElementIndexs.empty() )
    {
      const ShieldingSourceFitCalc::IsoFitStruct &previso = m_nuclides[thisElementIndexs[0]];
      newIso.ageDefiningNuc = previso.nuclide;
      //There is a slight hickup with emitting a datachanged and then inserting
      //  a row; if there wasnt this hickup, it would be nice to use the below
//      if( newIso.age < previso.age )
//        setSharedAgeNuclide( previso.nuclide, newIso.nuclide );
//      else
//        newIso.ageDefiningNuc = previso.nuclide;
    }//if( !newIso.ageDefiningNuc && !thisElementIndexs.empty() )
  }//if( m_sameAgeForIsotopes && newIso.ageIsFittable )
  
  
  
  const int npeaks = m_peakModel->rowCount();
  set<const SandiaDecay::Nuclide *> progeny;
  for( int peakn = 0; peakn < npeaks; ++peakn )
  {
    const PeakShrdPtr peak = m_peakModel->peak( m_peakModel->index(peakn,0) );
    if( peak && (peak->parentNuclide() == newIso.nuclide) && peak->useForShieldingSourceFit() )
    {
      const SandiaDecay::Transition * const trans = peak->nuclearTransition();
      if( trans && trans->parent )
        progeny.insert( trans->parent );
    }
  }//for( loop over peaks)
  
  newIso.numProgenyPeaksSelected = progeny.size();
  
  std::vector<ShieldingSourceFitCalc::IsoFitStruct>::iterator pos;
  pos = std::lower_bound( m_nuclides.begin(), m_nuclides.end(), newIso,
                          boost::bind( &SourceFitModel::compare, boost::placeholders::_1,
                                      boost::placeholders::_2, m_sortColumn, m_sortOrder ) );

  const int row = static_cast<int>( pos - m_nuclides.begin() );
  beginInsertRows( WModelIndex(), row, row );
  m_nuclides.insert( pos, newIso );
  endInsertRows();
}//void SourceFitModel::insertPeak( const PeakShrdPtr peak )


void SourceFitModel::peakModelRowsInsertedCallback( Wt::WModelIndex /*parent*/,
                                                     int firstRow, int lastRow )
{
//  if( firstRow != lastRow )
//  {
//    WStringStream msg;
//    msg << "SourceFitModel::peakModelRowsInsertedCallback(...)"
//            " first row should equal last row: "
//         << firstRow << " != " << lastRow;
//    cerr << msg.str() << endl;
//    throw std::runtime_error( msg.str() );
//  }//if( firstRow != lastRow )

  for( int row = firstRow; row <= lastRow; ++row )
  {
    const PeakShrdPtr peak = m_peakModel->peak( m_peakModel->index(row,0) );
    if( peak && peak->parentNuclide() )
      insertPeak( peak );
   }//for( int row = firstRow; row <= lastRow; ++row )
}//void peakModelRowsInsertedCallback(...)


void SourceFitModel::peakModelRowsRemovedCallback( Wt::WModelIndex /*index*/,
                                                    int firstRow, int lastRow )
{
  if( m_nuclides.empty() )
    return;

  vector<PeakShrdPtr> livingPeaks;
  const int nrow = m_peakModel->rowCount();
  
//  for( int row = firstRow; row <= lastRow; ++row )
  for( int row = 0; row < nrow; ++row )
  {
    if( row >= firstRow && row <= lastRow )
      continue;

    const PeakShrdPtr peak = m_peakModel->peak( m_peakModel->index(row,0) );

    if( !peak )
    {  //probably wont ever happen
      cerr << "SourceFitModel::peakModelRowsRemovedCallback(...)\n\tNot a valid row being removed" << endl;
      continue;
    }//if( !peak )

    livingPeaks.push_back( peak );
  }//for( int row = firstRow; row <= lastRow; ++row )
  
  vector<int> removedRows;
  vector<const SandiaDecay::Nuclide *> removed;
  
  for( int row = firstRow; row <= lastRow; ++row )
  {
    const PeakShrdPtr peak = m_peakModel->peak( m_peakModel->index(row,0) );

    if( !peak )
      continue;

    const SandiaDecay::Nuclide *nuc = peak->parentNuclide();

    if( !nuc )
      continue;

    bool shouldRemove = true;
    for( PeakShrdPtr p : livingPeaks )
    {
      if( (p != peak) && (p->parentNuclide() == nuc) )
      {
        shouldRemove = false;
        break;
      }//if( (p != peak) && (p->parentNuclide() == nuc) )
    }//for( PeakShrdPtr p : livingPeaks )

    if( !shouldRemove )
      continue;

    if( std::find(removed.begin(), removed.end(), nuc) != removed.end() )
      continue;

    size_t index = m_nuclides.size();
    for( size_t i = 0; i < index; ++i )
      if( m_nuclides[i].nuclide == peak->parentNuclide() )
        index = i;
    
    if( index == m_nuclides.size() )
    {
      cerr << "SourceFitModel::peakModelRowsRemovedCallback(...)\n\tSerious logic error, fix this ish" << endl;
      cerr << "m_nuclides.size()=" << m_nuclides.size()
           << ", m_peakModel->npeaks()=" << m_peakModel->npeaks()
           << " nuc->symbol=" << nuc->symbol << endl;
      continue;
    }//if( index == m_nuclides.size() )

    removed.push_back( nuc );
    
    m_previousResults[nuc] = m_nuclides[index];
    
    const int iindex = static_cast<int>(index);
    removedRows.push_back( iindex );
    beginRemoveRows( WModelIndex(), iindex, iindex );
    m_nuclides.erase( m_nuclides.begin() + iindex );
    endRemoveRows();
  }//for( int row = firstRow; row <= lastRow; ++row )
  
  //Now we have to go back through and make sure the removed nuclides werent
  //  actually the age defining for an existing nuclide
  for( const SandiaDecay::Nuclide *nuc : removed )
  {
    const SandiaDecay::Nuclide *newDefining = nullptr;
    double minage = std::numeric_limits<double>::infinity();
    vector<const SandiaDecay::Nuclide *> nucstochange;
    
    for( ShieldingSourceFitCalc::IsoFitStruct &ifs : m_nuclides )
    {
      if( ifs.ageDefiningNuc == nuc )
      {
        nucstochange.push_back( ifs.nuclide );
        if( ifs.age < minage )
        {
          minage = ifs.age;
          newDefining = ifs.nuclide;
        }
      }//if( ifs.ageDefiningNuc == nuc )
    }//for( ShieldingSourceFitCalc::IsoFitStruct &ifs : m_nuclides )
    
    for( const SandiaDecay::Nuclide *nuc : nucstochange )
      setSharedAgeNuclide( nuc, newDefining );
  }//for( const SandiaDecay::Nuclide *nuc : removed )
  
//  for( const int iindex : removedRows )
//  {
//    beginRemoveRows( WModelIndex(), iindex, iindex );
//    m_nuclides.erase( m_nuclides.begin() + iindex );
//    endRemoveRows();
//  }//for( const int iindex : removeRows )
}//void peakModelRowsRemovedCallback( Wt::WModelIndex index, int firstRow, int lastRow )


bool is_within( int a, int start, int end )
{
  return (a>=start && a<=end);
}


void SourceFitModel::peakModelDataChangedCallback( Wt::WModelIndex topLeft,
                                                    Wt::WModelIndex bottomRight )
{
  const PeakModel::Columns colstart = PeakModel::Columns( topLeft.column() );
  const PeakModel::Columns colend = PeakModel::Columns( bottomRight.column() );
  const bool changedIso = is_within( PeakModel::kIsotope, colstart, colend);
  const bool changedFit = is_within( PeakModel::kUseForShieldingSourceFit, colstart, colend);
  
  const bool dirty = (changedIso || changedFit);
  
  if( !dirty )
  {
    const bool changedPhotoPeak
                    = is_within( PeakModel::kPhotoPeakEnergy, colstart, colend);
    const bool changedAmp = is_within( PeakModel::kAmplitude, colstart, colend);
    
    //Not dirty as far as this SourceFitModel is concerned, but as a convience
    //  to things hooked to this model, lets see if we chould redraw the chi2
    //  chart, and if do, emit a dataChanged() signal to indicate this.
    if( changedAmp || changedPhotoPeak )
    {
      try
      {
        cerr << "topLeft.isValid()=" << topLeft.isValid() << endl;
        const PeakModel::PeakShrdPtr peak = m_peakModel->peak(topLeft);
        cerr << "!!peak=" << !!peak << endl;
        if( !!peak )
          cerr << "peak->useForShieldingSourceFit()=" << peak->useForShieldingSourceFit() << endl;
        
        const bool beingUsed = (!!peak && peak->useForShieldingSourceFit());
        if( beingUsed )
          dataChanged().emit( WModelIndex(), WModelIndex() );
      }catch(...)
      {
        cerr << "Unexpected exception" << endl;
      }
    }//if( peak amplitude changed )
    
    return;
  }//if( !dirty )

  set<const SandiaDecay::Nuclide *> prenucs, postnucs;
  
  vector<const SandiaDecay::Nuclide *> preisotopes, postisotopes;
  for( const ShieldingSourceFitCalc::IsoFitStruct &ifs : m_nuclides )
  {
    prenucs.insert( ifs.nuclide );
    preisotopes.push_back( ifs.nuclide );
  }
  
  const size_t npreisotopes = m_nuclides.size();
  
  repopulateIsotopes();
  
  for( const ShieldingSourceFitCalc::IsoFitStruct &ifs : m_nuclides )
  {
    postnucs.insert( ifs.nuclide );
    postisotopes.push_back( ifs.nuclide );
  }
  
  
  if( m_sameAgeForIsotopes && (prenucs != postnucs) )
  {
    vector<const SandiaDecay::Nuclide *> removednucs, addednucs;
    for( const SandiaDecay::Nuclide *nuc : preisotopes )
    {
      if( !std::count(postisotopes.begin(),postisotopes.end(),nuc) )
        removednucs.push_back( nuc );
    }
    
    for( const SandiaDecay::Nuclide *nuc : postisotopes )
    {
      if( !std::count(preisotopes.begin(),preisotopes.end(),nuc) )
        addednucs.push_back( nuc );
    }
    
    for( const SandiaDecay::Nuclide *nuc : removednucs )
    {
      bool removedADefining = false;
      const SandiaDecay::Nuclide *defining = NULL;
      double minage = std::numeric_limits<double>::infinity();
      
      for( ShieldingSourceFitCalc::IsoFitStruct &ifs : m_nuclides )
      {
        if( ifs.ageIsFittable && (ifs.nuclide->atomicNumber == nuc->atomicNumber) )
        {
          if( ifs.age < minage )
          {
            minage = ifs.age;
            defining = ifs.nuclide;
          }
          removedADefining = (removedADefining || (ifs.ageDefiningNuc==nuc));
        }//if( ifs.nuclide->atomicNumber == nuc->atomicNumber )
      }//for( ShieldingSourceFitCalc::IsoFitStruct &ifs : m_nuclides )
      
      if( removedADefining )
      {
        for( ShieldingSourceFitCalc::IsoFitStruct &ifs : m_nuclides )
        {
          if( ifs.ageIsFittable && (ifs.nuclide->atomicNumber == nuc->atomicNumber) )
            setSharedAgeNuclide( ifs.nuclide, defining );
        }//for( ShieldingSourceFitCalc::IsoFitStruct &ifs : m_nuclides )
      }//if( removedADefining )
    }//for( const SandiaDecay::Nuclide *nuc : removednucs )
    
    for( const SandiaDecay::Nuclide *nuc : addednucs )
    {
      if( PeakDef::ageFitNotAllowed( nuc ) )
        continue;
      
      const SandiaDecay::Nuclide *defining = NULL;
      double minage = std::numeric_limits<double>::infinity();
      
      for( ShieldingSourceFitCalc::IsoFitStruct &ifs : m_nuclides )
      {
        if( ifs.ageIsFittable
           && (ifs.age < minage)
           && (ifs.nuclide->atomicNumber == nuc->atomicNumber) )
        {
          minage = ifs.age;
          defining = ifs.nuclide;
        }
      }//for( ShieldingSourceFitCalc::IsoFitStruct &ifs : m_nuclides )
      
      if( !defining )
        defining = nuc;
      
      for( ShieldingSourceFitCalc::IsoFitStruct &ifs : m_nuclides )
      {
        if( ifs.ageIsFittable && (ifs.nuclide->atomicNumber == nuc->atomicNumber) )
          setSharedAgeNuclide( ifs.nuclide, defining );
      }
    }//for( const SandiaDecay::Nuclide *nuc : addednucs )
  }//if( m_sameAgeForIsotopes && (preisotopes != postisotopes) )
  
  
  if( npreisotopes == m_nuclides.size() )
  { //No signals got emmited to update the chi2Chart, so lets check if we should
    try
    {
      const PeakModel::PeakShrdPtr &peak = m_peakModel->peak(topLeft);
      if( changedFit || (peak && peak->useForShieldingSourceFit()) )
        dataChanged().emit( index(0,0), index(rowCount()-1,columnCount()-1) );
      
    }catch(...){}
  }//if( npreisotopes == m_nuclides.size() )
  
}//void peakModelDataChangedCallback(...);


void SourceFitModel::peakModelResetCallback()
{
  if( m_nuclides.empty() )
    return;

  for( const ShieldingSourceFitCalc::IsoFitStruct &iss : m_nuclides )
    m_previousResults[iss.nuclide] = iss;
  
  beginRemoveRows( WModelIndex(), 0, static_cast<int>(m_nuclides.size())-1 );
  m_nuclides.clear();
  endRemoveRows();
}//void peakModelResetCallback();


void SourceFitModel::repopulateIsotopes()
{
  //go through and see if we have to delete any peaks
  vector<int> indexs_to_remove;
  const int norigoiso = static_cast<int>( m_nuclides.size() );

  for( int ison = 0; ison < norigoiso; ++ison )
  {
    ShieldingSourceFitCalc::IsoFitStruct &isof = m_nuclides[ison];

    size_t numSourcePeaks = 0; // I think we can remove this and instead use only progeny.size(), but havent tested
    set<const SandiaDecay::Nuclide *> progeny;
    const int npeaks = m_peakModel->rowCount();
    for( int peakn = 0; peakn < npeaks; ++peakn )
    {
      const PeakShrdPtr peak = m_peakModel->peak( m_peakModel->index(peakn,0) );

      if( (peak->parentNuclide() == isof.nuclide) && peak->useForShieldingSourceFit() )
      {
        const SandiaDecay::Transition * const trans = peak->nuclearTransition();
        if( trans && trans->parent )
          progeny.insert( trans->parent );
        numSourcePeaks += 1;
      }
    }//for( loop over peaks)

    isof.numProgenyPeaksSelected = progeny.size();
    
    if( numSourcePeaks == 0 )
    {
      m_previousResults[isof.nuclide] = isof;
      indexs_to_remove.push_back( ison );
    }
  }//for( const ShieldingSourceFitCalc::IsoFitStruct &isof : m_nuclides )

//  cerr << "indexs_to_remove.size()=" << indexs_to_remove.size() << endl;

  for( vector<int>::reverse_iterator iter = indexs_to_remove.rbegin();
       iter != indexs_to_remove.rend(); ++iter )
  {
    const int iindex = *iter;
    beginRemoveRows( WModelIndex(), iindex, iindex );
    m_nuclides.erase( m_nuclides.begin() + iindex );
    endRemoveRows();
  }//if( found != m_peakModel->rowCount() )


  //go through and see if we have to add any peaks
  for( int peakn = 0; peakn < m_peakModel->rowCount(); ++peakn )
  {
    const PeakShrdPtr peak = m_peakModel->peak( m_peakModel->index(peakn,0) );
    if( !peak || !peak->useForShieldingSourceFit() )
       continue;

    bool found = false;
    for( const ShieldingSourceFitCalc::IsoFitStruct &isof : m_nuclides )
    {
      found = (peak->parentNuclide() == isof.nuclide);
      if( found )
        break;
    }//for( const ShieldingSourceFitCalc::IsoFitStruct &isof : m_nuclides )

    if( !found )
      insertPeak( peak );
  }//for( loop over peaks)
    
  // Double check that things get refreshed.
  //reset();
  layoutAboutToBeChanged().emit();
  layoutChanged().emit();
}//void SourceFitModel::repopulateIsotopes()



int SourceFitModel::columnCount( const Wt::WModelIndex & ) const
{
  return static_cast<int>( kNumColumns );
}//int columnCount( const Wt::WModelIndex & ) const


int SourceFitModel::rowCount( const Wt::WModelIndex & parent) const
{
  if (parent.isValid())
  {
    return 0;
  } //parent.isValid()
  return static_cast<int>( m_nuclides.size() );
}//int rowCount( const Wt::WModelIndex & ) const


Wt::WFlags<Wt::ItemFlag> SourceFitModel::flags( const Wt::WModelIndex &index ) const
{
  const int row = index.row();
  const int nrow = static_cast<int>( m_nuclides.size() );

  if( row>=nrow || row<0 )
    return WFlags<ItemFlag>();

  const ShieldingSourceFitCalc::IsoFitStruct &iso = m_nuclides[row];

  switch( index.column() )
  {
    case kIsotope:
      break;
    case kActivity:
      switch( iso.sourceType )
      {
        case ShieldingSourceFitCalc::ModelSourceType::Point:
          return WFlags<ItemFlag>(Wt::ItemIsEditable);
          
        case ShieldingSourceFitCalc::ModelSourceType::Intrinsic:
        case ShieldingSourceFitCalc::ModelSourceType::Trace:
          return WFlags<ItemFlag>();
      }//switch( iso.sourceType )
      break;
      
    case kFitActivity:
      switch( iso.sourceType )
      {
        case ShieldingSourceFitCalc::ModelSourceType::Point:
          return WFlags<ItemFlag>(ItemIsUserCheckable);
          
        case ShieldingSourceFitCalc::ModelSourceType::Intrinsic:
        case ShieldingSourceFitCalc::ModelSourceType::Trace:
          return WFlags<ItemFlag>();
      }//switch( iso.sourceType )
      break;
      
    case kAge:
      if( !iso.ageIsFittable && !iso.ageDefiningNuc )
        return WFlags<ItemFlag>();
      
      if( iso.ageDefiningNuc )
      {
        for( const ShieldingSourceFitCalc::IsoFitStruct &isodef : m_nuclides )
        {
          if( isodef.nuclide == iso.ageDefiningNuc )
            return isodef.ageIsFittable ? WFlags<ItemFlag>(Wt::ItemIsEditable) : WFlags<ItemFlag>();
        }
        return WFlags<ItemFlag>();
      }//if( iso.ageDefiningNuc )
      
      return WFlags<ItemFlag>(Wt::ItemIsEditable);
    case kFitAge:
//      if( iso.shieldingIsSource )
//        return WFlags<ItemFlag>();
      return WFlags<ItemFlag>(ItemIsUserCheckable);
    case kIsotopeMass:
      return WFlags<ItemFlag>();
    case kActivityUncertainty:
    case kAgeUncertainty:
      return WFlags<ItemFlag>();
      
#if( INCLUDE_ANALYSIS_TEST_SUITE )
    case kTruthActivity:
    case kTruthActivityTolerance:
    case kTruthAge:
    case kTruthAgeTolerance:
      return WFlags<ItemFlag>(Wt::ItemIsEditable);
#endif
      
    case kNumColumns:
      break;
  }//switch( section )

  return WFlags<ItemFlag>();
}//Wt::WFlags<Wt::ItemFlag> flags( const Wt::WModelIndex &index ) const




boost::any SourceFitModel::headerData( int section, Orientation orientation, int role ) const
{
  //When orientation is Horizontal, section is a column number,
  //  when orientation is Vertical, section is a row (peak) number.

  if( role == LevelRole )
    return 0;

  if( (orientation != Horizontal)
      || ((role != DisplayRole)
           &&  (role != Wt::ToolTipRole))  )
    return WAbstractItemModel::headerData( section, orientation, role );

  if( role == Wt::ToolTipRole )
  {
    const char *tooltip_key = nullptr;
    switch( section )
    {
      case kIsotope:             tooltip_key = "sfm-tt-header-isotope";    break;
      case kActivity:            tooltip_key = "sfm-tt-header-activity";   break;
      case kAge:                 tooltip_key = "sfm-tt-header-age";        break;
      case kFitActivity:         tooltip_key = "sfm-tt-header-fit-act";    break;
      case kFitAge:              tooltip_key = "sfm-tt-header-fit-age";    break;
      case kIsotopeMass:         tooltip_key = "sfm-tt-header-iso-mass";   break;
      case kActivityUncertainty: tooltip_key = "sfm-tt-header-act-uncert"; break;
      case kAgeUncertainty:      tooltip_key = "sfm-tt-header-age-uncert"; break;

#if( INCLUDE_ANALYSIS_TEST_SUITE )
      case kTruthActivity:
      case kTruthActivityTolerance:
      case kTruthAge:
      case kTruthAgeTolerance:
        break;
#endif
        
      case kNumColumns:
        break;
    }//switch( section )
    if( !tooltip_key )
      return boost::any();
    
    return boost::any( WString::tr(tooltip_key) );
  }//if( role == Wt::ToolTipRole )

  //If we're here, role==DisplayRole
  switch( section )
  {
    case kIsotope:     return boost::any( WString::tr("Nuclide") );
    case kActivity:
      switch( m_det_type )
      {
        case DetectorPeakResponse::EffGeometryType::FarField:
        case DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct:
          return boost::any( WString::tr("Activity") );
          
        case DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2:
        case DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2:
        case DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram:
          return boost::any( WString("{1}{2}")
                            .arg( WString::tr("Activity"))
                            .arg( DetectorPeakResponse::det_eff_geom_type_postfix(m_det_type) ) );
      }//switch( m_det_type )
      assert( 0 );
      break;
      
    case kAge:         return boost::any( WString::tr("Age") );
    case kFitActivity: return boost::any( WString::tr("sfm-header-fit-act") );
    case kFitAge:      return boost::any( WString::tr("sfm-header-fit-age") );
    case kIsotopeMass: return boost::any( WString::tr("Mass") );
    case kActivityUncertainty: return boost::any( WString::tr("sfm-header-act-uncert") );
    case kAgeUncertainty:      return boost::any( WString::tr("sfm-header-age-uncert") );
#if( INCLUDE_ANALYSIS_TEST_SUITE )
    case kTruthActivity:          return boost::any( WString("Truth Act.") );
    case kTruthActivityTolerance: return boost::any( WString("Truth Act. Tol.") );
    case kTruthAge:               return boost::any( WString("Truth Age") );
    case kTruthAgeTolerance:      return boost::any( WString("Truth Age Tol.") );
#endif
    case kNumColumns:  return boost::any();
  }//switch( section )

  return boost::any();
}//headerData(...)


Wt::WModelIndex SourceFitModel::parent( const Wt::WModelIndex &index ) const
{
  return WModelIndex();
}//Wt::WModelIndex parent(...) const


boost::any SourceFitModel::data( const Wt::WModelIndex &index, int role ) const
{
  //should consider implementing ToolTipRole
  if( (role != Wt::DisplayRole) && (role != Wt::EditRole) && (role != Wt::ToolTipRole)
     && (role != (Wt::ItemDataRole::UserRole + 10))
    && !((role==Wt::CheckStateRole) && ((index.column()==kFitActivity) || (index.column()==kFitAge))) )
  {
    return boost::any();
  }

  
  const int row = index.row();
  const int column = index.column();
  const int nrows = static_cast<int>( m_nuclides.size() );
  
  if( row<0 || column<0 || column>=kNumColumns || row>=nrows )
    return boost::any();

  const bool extra_precision = (role == (Wt::ItemDataRole::UserRole + 10));
  const ShieldingSourceFitCalc::IsoFitStruct &isof = m_nuclides[row];

  if( role == Wt::ToolTipRole )
  {
    if( column == kIsotope )
    {
      WString msg = WString("{1}={2}")
                      .arg( WString::tr("T1/2") )
                      .arg( PhysicalUnitsLocalized::printToBestTimeUnits( isof.nuclide->halfLife, 2, PhysicalUnits::second ) );
      // TODO: should put in dominant photopeaks or other information here
      return boost::any( msg );
    }else if( column == kAge )
    {
      if( !isof.ageIsFittable )
        return boost::any( WString::tr("sfm-tt-aging-not-allowed") );
    }//if / else to determine column

    return boost::any();
  }//if( role == Wt::ToolTipRole )

  switch( column )
  {
    case kIsotope:
      return boost::any( WString(isof.nuclide->symbol) );
      
    case kActivity:
    {
      const double act = isof.activity;
      string ans = PhysicalUnits::printToBestActivityUnits( act, (extra_precision ? 8 : 3), m_displayCuries );
      ans += DetectorPeakResponse::det_eff_geom_type_postfix(m_det_type);
      
      // We'll require the uncertainty to be non-zero to show it - 1 micro-bq is an arbitrary cutoff to
      //  consider anything below it zero.
      const double uncert = isof.activityUncertainty;
      
      bool shouldHaveUncert = ( (uncert > 1.0E-6*PhysicalUnits::bq)
                               || ((uncert > 0.0) && (uncert > 1.0E-6*fabs(act)) ) );
      
      switch( isof.sourceType )
      {
        case ShieldingSourceFitCalc::ModelSourceType::Intrinsic:
          // TODO: check if dimensions are being fit - but I *think* uncert will only be >0 if this is the case - so we are good?
          break;
          
        case ShieldingSourceFitCalc::ModelSourceType::Trace:
          // TODO: check if activity is being fit, and if not, then if dimensions are being fit - but I *think* uncert will only be >0 if one of these is the case anyway - so we are good?
          break;
          
        case ShieldingSourceFitCalc::ModelSourceType::Point:
          // Dont show uncertainty unless we fit for it - although in principle it should be zero, right
          if( !fitActivity(row) )
            shouldHaveUncert = false;
          break;
      }//switch( iso.sourceType )
      
      if( shouldHaveUncert )
      {
        ans += " \xC2\xB1 " + PhysicalUnits::printToBestActivityUnits( uncert, (extra_precision ? 5 : 2), m_displayCuries );
        ans += DetectorPeakResponse::det_eff_geom_type_postfix(m_det_type);
      }//if( shouldHaveUncert )
      
      return boost::any( WString(ans) );
    }//case kActivity:

    case kFitActivity:
    {
      switch( isof.sourceType )
      {
        case ShieldingSourceFitCalc::ModelSourceType::Point:
          return boost::any( isof.fitActivity );
          
        case ShieldingSourceFitCalc::ModelSourceType::Intrinsic:
        case ShieldingSourceFitCalc::ModelSourceType::Trace:
          return boost::any();
      }//switch( iso.sourceType )
    }//case kFitActivity:

    case kAge:
    {
//      if( isof.shieldingIsSource )
//        return boost::any();
      double age = 0.0, uncert = 0.0;
      const SandiaDecay::Nuclide *nuc = nullptr;
      if( isof.ageDefiningNuc && (isof.ageDefiningNuc != isof.nuclide) )
      {
        const int ind = nuclideIndex( isof.ageDefiningNuc );
        if( ind >= 0 )
        {
          nuc = isof.ageDefiningNuc;
          age = m_nuclides[ind].age;
          if( m_nuclides[ind].ageIsFittable && m_nuclides[ind].fitAge )
            uncert = m_nuclides[ind].ageUncertainty;
        }else
        {
          nuc = isof.nuclide;
          age = isof.age;
          if( isof.ageIsFittable && isof.fitAge )
            uncert = isof.ageUncertainty;
          cerr << "SourceFitModel::data: ran into error when retrieving"
               << " age for a nuclide with a defining age isotope that isnt in"
               << " the model; charging on!" << endl;
        }//if( ind >= 0 )
      }else
      {
        nuc = isof.nuclide;
        age = isof.age;
        if( isof.ageIsFittable && isof.fitAge )
          uncert = isof.ageUncertainty;
      }//if( isof.ageDefiningNuc && (isof.ageDefiningNuc!=isof.nuclide) )
      
      if( !isof.ageIsFittable )
        return boost::any( WString::tr("NA") );
      
      string ans = PhysicalUnitsLocalized::printToBestTimeUnits( age, (extra_precision ? 8 : 2) );
      if( uncert > 0.0 )
        ans += " \xC2\xB1 " + PhysicalUnitsLocalized::printToBestTimeUnits( uncert, (extra_precision ? 8 : 1) );
      
      return boost::any( WString(ans) );
    }//case kAge:

    case kFitAge:
    {
//      if( isof.shieldingIsSource )
//        return boost::any();
      if( isof.ageDefiningNuc && (isof.ageDefiningNuc != isof.nuclide) )
        return boost::any();
      
      // Make sure there is more than two peaks being fitted for this nuclide to enable fitting age
      //  (I guess you could fix activity, and shielding, and select a progeny peak, and fit for
      //   age based on that peak growing in, but this probably isnt realistically ever done, but if
      //   you did want to do it, you could round-about calculate it)
      if( (isof.numProgenyPeaksSelected <= 1) && !isof.ageDefiningNuc )
        return boost::any();
      
      if( !isof.ageIsFittable )
        return boost::any();
      
      return boost::any( isof.fitAge );
    }//case kFitAge:

    case kIsotopeMass:
    {
      const double act = isof.activity;
      const double mass_grams = act / isof.nuclide->activityPerGram();

      if( IsInf(mass_grams) || IsNan(mass_grams) )
        return boost::any();

      return boost::any( WString(PhysicalUnits::printToBestMassUnits(mass_grams,(extra_precision ? 8 : 3),1.0)) );
    }//case kIsotopeMass:

    case kActivityUncertainty:
    {
      if( isof.activityUncertainty < 0.0 )
        return boost::any();
      
      double act = isof.activityUncertainty;
      string ans = PhysicalUnits::printToBestActivityUnits( act, (extra_precision ? 8 : 2), m_displayCuries );
      ans += DetectorPeakResponse::det_eff_geom_type_postfix(m_det_type);
      
      return boost::any( WString(ans) );
    }//case kActivityUncertainty:

    case kAgeUncertainty:
    {
      if( (!isof.ageIsFittable) || isof.ageUncertainty < 0.0 )
        return boost::any();
      const string ans = PhysicalUnitsLocalized::printToBestTimeUnits( isof.ageUncertainty, (extra_precision ? 8 : 2) );
      return boost::any( WString(ans) );
    }//case kAgeUncertainty:

#if( INCLUDE_ANALYSIS_TEST_SUITE )
    case kTruthActivity:
    {
      if( !isof.truthActivity )
        return boost::any();
      
      string ans = PhysicalUnits::printToBestActivityUnits( *isof.truthActivity, (extra_precision ? 8 : 4), m_displayCuries );
      ans += DetectorPeakResponse::det_eff_geom_type_postfix(m_det_type);
      
      return boost::any( WString(ans) );
    }
      
    case kTruthActivityTolerance:
    {
      if( !isof.truthActivityTolerance )
        return boost::any();
      string ans = PhysicalUnits::printToBestActivityUnits( *isof.truthActivityTolerance, (extra_precision ? 8 : 4), m_displayCuries );
      ans += DetectorPeakResponse::det_eff_geom_type_postfix(m_det_type);
      
      return boost::any( WString(ans) );
    }
      
    case kTruthAge:
    {
      if( !isof.truthAge )
        return boost::any();
      const string ans = PhysicalUnitsLocalized::printToBestTimeUnits( *isof.truthAge, (extra_precision ? 8 : 4) );
      return boost::any( WString(ans) );
    }
      
    case kTruthAgeTolerance:
    {
      if( !isof.truthAgeTolerance )
        return boost::any();
      const string ans = PhysicalUnitsLocalized::printToBestTimeUnits( *isof.truthAgeTolerance, (extra_precision ? 8 : 4) );
      return boost::any( WString(ans) );
    }
#endif  //#if( INCLUDE_ANALYSIS_TEST_SUITE )
      
    case kNumColumns:
      return boost::any();
  }//switch( column )

  return boost::any();
}//boost::any data(...) const


Wt::WModelIndex SourceFitModel::index( int row, int column,
                                        const Wt::WModelIndex &/*parent*/ ) const
{
  return WAbstractItemModel::createIndex( row, column, (void *)0 );
}//Wt::WModelIndex index(...) const



WModelIndex SourceFitModel::index( const SandiaDecay::Nuclide *nuc,
                                       SourceFitModel::Columns column ) const
{
  const int nrow = static_cast<int>( m_nuclides.size() );
  for( int row = 0; row < nrow; ++row )
    if( m_nuclides[row].nuclide == nuc )
      return index( row, column );
  return WModelIndex();
}//WModelIndex index(...) const


WModelIndex SourceFitModel::index( const string &symbol,
                                   SourceFitModel::Columns column ) const
{
  const int nrow = static_cast<int>( m_nuclides.size() );
  for( int row = 0; row < nrow; ++row )
    if( m_nuclides[row].nuclide && (m_nuclides[row].nuclide->symbol == symbol) )
      return index( row, column );
  return WModelIndex();
}//WModelIndex index(...) const


bool SourceFitModel::setData( const Wt::WModelIndex &index, const boost::any &value, int role )
{
  try
  {
    if( !index.isValid() )
      return false;

    if( (role != Wt::EditRole) && (role != Wt::CheckStateRole) )
      return false;

    const int row = index.row();
    const int column = index.column();
    const int nrows = static_cast<int>( m_nuclides.size() );

#if( INCLUDE_ANALYSIS_TEST_SUITE )
    if( value.empty() )
    {
      switch( SourceFitModel::Columns(column) )
      {
        case kTruthActivity:
        case kTruthActivityTolerance:
        case kTruthAge:
        case kTruthAgeTolerance:
        case kActivityUncertainty:
        case kAgeUncertainty:
          break;
        default:
          return false;
      }//switch( column )
    }//if( value.empty() )
#else
    if( value.empty() && (column != kAgeUncertainty) && (column != kActivityUncertainty) )
      return false;
#endif  //#if( INCLUDE_ANALYSIS_TEST_SUITE )
    
    
    if( (row < 0) || (column < 0) || (column >= kNumColumns) || (row >= nrows) )
      return false;
    
    if( (role == Wt::CheckStateRole) && (column != kFitActivity) && (column != kFitAge) )
      return false;
    
    // In principle we should only let (role == Wt::CheckStateRole) allow to change this value,
    //  but it doesnt seem any harm in also allowing EditRole (the default value for the fcn call)
    //  affect this value, as long as a bool is passed in.
    bool boolean_val = false;
    if( (column == kFitActivity) || (column == kFitAge) )
    {
      if( value.type() != boost::typeindex::type_id<bool>().type_info() )
        return false;
      
      boolean_val = boost::any_cast<bool>( value );
    }//if( (column==kFitActivity) || (column==kFitAge) )

    const WString txt_val = asString( value );
    string utf_str = txt_val.toUTF8();

    if( (column != kFitActivity) && (column != kFitAge) && (column != kActivityUncertainty) && txt_val.empty() )
      return false;

    // filter out the +- and everything after.
    if( (column == kActivity) || (column == kIsotopeMass) || (column == kAge) )
    {
      const auto pos = utf_str.find( "\xC2\xB1" );
      if( pos != string::npos )
        utf_str = utf_str.substr(0, pos);
    }//if( we might have the +- )
    
    ShieldingSourceFitCalc::IsoFitStruct &iso = m_nuclides[row];
    
    
    // To facilitate undo/redo we could grab value for the row/column we are changing, and then
    //  set things back, but if we do this we will lose uncertainty, and maybe other stuff; we
    //  could work around this, but instead for the moment we will just copy all the data
#define SOURCE_FIT_MODEL_FULL_COPY_UNDO_REDO 1
#if( SOURCE_FIT_MODEL_FULL_COPY_UNDO_REDO )
    const auto prev_data = make_shared<const vector<ShieldingSourceFitCalc::IsoFitStruct>>( m_nuclides );
#else
    const boost::any prev_value = SourceFitModel::data( index, Wt::ItemDataRole::UserRole + 10 );
    const boost::any new_value = value;
#endif
    

    //If were here, all is fine
    switch( column )
    {
      case kIsotope:
        return false;
        
      case kActivity:
      {
        iso.activity = PhysicalUnits::stringToActivity( utf_str );
        
        if( iso.activityUncertainty >= 0.0 ) // For activity we will emit the whole row changed below
          iso.activityUncertainty = -1.0;
      
        break;
      }//case kActivity:

      case kFitActivity:
        iso.fitActivity = boolean_val;
      break;

      case kAge:
      {
        if( !iso.ageIsFittable )
        {
          break;
        }else
        {
          const double hl = (iso.nuclide ? iso.nuclide->halfLife : -1.0);
          iso.age = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( utf_str, hl );

          if( iso.ageUncertainty >= 0.0 )
          {
            iso.ageUncertainty = -1.0;
            WModelIndex uncertIndex = createIndex(row, kAgeUncertainty,(void *)0);
            dataChanged().emit( uncertIndex, uncertIndex );
          }//if( iso.ageUncertainty >= 0.0 )
          
          if( m_sameAgeForIsotopes )
          {
            const SandiaDecay::Nuclide *ageNuc = iso.ageDefiningNuc ? iso.ageDefiningNuc : iso.nuclide;
            assert( ageNuc );
            
            for( size_t i = 0; i < m_nuclides.size(); ++i )
            {
              const int thisrow = static_cast<int>(i);
              
              if( (m_nuclides[i].ageDefiningNuc == ageNuc) || (m_nuclides[i].nuclide == ageNuc) )
              {
                m_nuclides[i].age = iso.age;
                
                WModelIndex ind = createIndex( thisrow, kAge, (void *)0 );
                if( thisrow != row )  // We'll emit for 'row' later on
                  dataChanged().emit( ind, ind );
                
                if( iso.ageUncertainty >= 0.0 )
                {
                  ind = createIndex( thisrow, kAgeUncertainty, (void *)0 );
                  dataChanged().emit( ind, ind );
                }//if( iso.ageUncertainty >= 0.0 )
              }//if( nuc.ageDefiningNuc == iso.nuclide )
            }//for( ShieldingSourceFitCalc::IsoFitStruct *nuc : m_nuclides )
          }//if( m_sameAgeForIsotopes )
        }//if( decays to stable children / else )
        break;
      }//case kAge:

      case kFitAge:
        if( !iso.ageIsFittable )
        {
          if( iso.fitAge )
            iso.fitAge = false;
          else
            return false;
        }else
        {
          if(  boolean_val == iso.fitAge )  //dont change this if we dont have to
            return false;
          iso.fitAge = boolean_val;
        }//if( !ageIsFittable ) / else
      break;

      case kIsotopeMass:
        cerr << "SourceFitModel::setData(...)\n\tYou shouldnt be trying to set kIsotopeMass"
             << endl;
      break;

      case kActivityUncertainty:
        iso.activityUncertainty = -1.0;
        if( !value.empty() )
        {
          iso.activityUncertainty = PhysicalUnits::stringToActivity( utf_str );
        }//if( !value.empty() )
      break;

      case kAgeUncertainty:
      {
        iso.ageUncertainty = -1.0;
        if( iso.ageIsFittable && !value.empty() )
        {
          const double hl = (iso.nuclide ? iso.nuclide->halfLife : -1.0);
          iso.ageUncertainty = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( utf_str, hl );
        }//if( decays to stable children / else )
        
        if( m_sameAgeForIsotopes )
        {
          for( size_t i = 0; i < m_nuclides.size(); ++i )
          {
            if( m_nuclides[i].ageDefiningNuc == iso.nuclide )
            {
              WModelIndex ind = createIndex( static_cast<int>(i), kAgeUncertainty, (void *)0);
              dataChanged().emit( ind, ind );
            }//if( nuc.ageDefiningNuc == iso.nuclide )
          }//for( ShieldingSourceFitCalc::IsoFitStruct *nuc : m_nuclides )
        }//if( m_sameAgeForIsotopes )
        
        break;
      }//case kAgeUncertainty:

#if( INCLUDE_ANALYSIS_TEST_SUITE )
      case kTruthActivity:
        if( value.empty() )
          iso.truthActivity.reset();
        else
          iso.truthActivity = PhysicalUnits::stringToActivity( utf_str );
        break;
        
      case kTruthActivityTolerance:
        if( value.empty() )
          iso.truthActivityTolerance.reset();
        else
          iso.truthActivityTolerance = PhysicalUnits::stringToActivity( utf_str );
        break;
        
      case kTruthAge:
      {
        const double hl = (iso.nuclide ? iso.nuclide->halfLife : -1.0);
        if( value.empty() )
          iso.truthAge.reset();
        else
          iso.truthAge = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( utf_str, hl );
        break;
      }
         
      case kTruthAgeTolerance:
      {
        const double hl = (iso.nuclide ? iso.nuclide->halfLife : -1.0);
        if( value.empty() )
          iso.truthAgeTolerance.reset();
        else
          iso.truthAgeTolerance = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( utf_str, hl );
        break;
      }
#endif
        
      case kNumColumns:
        return false;
    }//switch( column )

    // We will emit that all columns have been updated, since setting activity uncert will mean
    //  we have to set activity text again (since it might have a +- entry), and similar
    {
      UndoRedoManager::BlockUndoRedoInserts block;
      dataChanged().emit( createIndex(row,0,nullptr), createIndex(row,kNumColumns-1,nullptr) );
    }
    //
    // We could be more selective and do something like:
    //if( column == kActivity )
    //  dataChanged().emit( createIndex(row,0,nullptr), createIndex(row,kNumColumns-1,nullptr) );
    //else
    //  dataChanged().emit( index, index );
    
    UndoRedoManager *undoRedo = UndoRedoManager::instance();
    if( undoRedo && undoRedo->canAddUndoRedoNow() )
    {
#if( SOURCE_FIT_MODEL_FULL_COPY_UNDO_REDO )
      const auto current_data = make_shared<const vector<ShieldingSourceFitCalc::IsoFitStruct>>( m_nuclides );
      auto undo_redo = [prev_data, current_data]( const bool is_undo ){
        InterSpec *viewer = InterSpec::instance();
        ShieldingSourceDisplay *srcfit = viewer ? viewer->shieldingSourceFit() : nullptr;
        SourceFitModel *srcmodel = srcfit ? srcfit->sourceFitModel() : nullptr;
        if( !srcmodel )
          return;
        
        srcmodel->layoutAboutToBeChanged().emit();
        srcmodel->m_nuclides = *(is_undo ? prev_data : current_data);
        std::sort( begin(srcmodel->m_nuclides), end(srcmodel->m_nuclides),
                   boost::bind( &SourceFitModel::compare,
                              boost::placeholders::_1, boost::placeholders::_2,
                               srcmodel->m_sortColumn, srcmodel->m_sortOrder ) );
        srcmodel->layoutChanged().emit();
      };
      undoRedo->addUndoRedoStep( [=](){ undo_redo(true); }, [=](){ undo_redo(false); },
                                "Set Shielding/Source nuclide info" );
#else
      auto undo = [prev_value,index,role](){
        InterSpec *viewer = InterSpec::instance();
        ShieldingSourceDisplay *srcfit = viewer ? viewer->shieldingSourceFit() : nullptr;
        SourceFitModel *srcmodel = srcfit ? srcfit->sourceFitModel() : nullptr;
        if( srcmodel )
          srcmodel->setData( index, prev_value, role );
      };
      
      auto redo = [new_value,index,role](){
        InterSpec *viewer = InterSpec::instance();
        ShieldingSourceDisplay *srcfit = viewer ? viewer->shieldingSourceFit() : nullptr;
        SourceFitModel *srcmodel = srcfit ? srcfit->sourceFitModel() : nullptr;
        if( srcmodel )
          srcmodel->setData( index, new_value, role );
      };
      
      undoRedo->addUndoRedoStep( undo, redo, "Set Shielding/Source nuclide info" );
  #endif
    }//if( undoRedo && !undoRedo->isInUndoOrRedo() )
    
  }catch( exception &e )
  {
    cerr << "SourceFitModel::setData(...)\n\tWarning: exception caught; what="
         << e.what() << endl;
    return false;
  }//try/catch

  return true;
}//bool setData(...)


bool SourceFitModel::compare( const ShieldingSourceFitCalc::IsoFitStruct &lhs_input,
                             const ShieldingSourceFitCalc::IsoFitStruct &rhs_input,
                             Columns sortColumn, Wt::SortOrder order )
{
#if( INCLUDE_ANALYSIS_TEST_SUITE )
  auto optionalLess = []( const boost::optional<double> &olhs, const boost::optional<double> &orhs) -> bool{
    if( (!olhs) != (!orhs) )
      return !olhs;
    
    if( !olhs )
      return false;
    
    return ((*olhs) < (*orhs));
  };
#endif
  
  const bool ascend = (order == Wt::AscendingOrder);
  const ShieldingSourceFitCalc::IsoFitStruct &lhs = ascend ? lhs_input : rhs_input;
  const ShieldingSourceFitCalc::IsoFitStruct &rhs = ascend ? rhs_input : lhs_input;
  
  switch( sortColumn )
  {
    case kIsotope:     return (lhs.nuclide->symbol < rhs.nuclide->symbol);
    case kActivity:    return (lhs.activity < rhs.activity);
    case kFitActivity: return (lhs.fitActivity < rhs.fitActivity);
    case kAge:         return (lhs.age < rhs.age);
    case kFitAge:      return (lhs.fitAge < rhs.fitAge);
    case kIsotopeMass: return ((lhs.activity/lhs.nuclide->activityPerGram())
                               < (rhs.activity/rhs.nuclide->activityPerGram()) );
    case kActivityUncertainty: return (lhs.activityUncertainty < rhs.activityUncertainty);
    case kAgeUncertainty:      return (lhs.ageUncertainty < rhs.ageUncertainty);
      
#if( INCLUDE_ANALYSIS_TEST_SUITE )
    case kTruthActivity:          return optionalLess( lhs.truthActivity, rhs.truthActivity );
    case kTruthActivityTolerance: return optionalLess( lhs.truthActivityTolerance, rhs.truthActivityTolerance );
    case kTruthAge:               return optionalLess( lhs.truthAge, rhs.truthAge );
    case kTruthAgeTolerance:      return optionalLess( lhs.truthAgeTolerance, rhs.truthAgeTolerance );
#endif
      
    case kNumColumns:  return false;
  }//switch( sortColumn )
  assert( 0 );
  return false;
}//bool compare(...);

void SourceFitModel::sort( int column, Wt::SortOrder order )
{
  layoutAboutToBeChanged().emit();
  m_sortOrder = order;
  m_sortColumn = Columns( column );
  std::sort( m_nuclides.begin(), m_nuclides.end(),
             boost::bind( &SourceFitModel::compare, boost::placeholders::_1,
                         boost::placeholders::_2, m_sortColumn, m_sortOrder ) );
  layoutChanged().emit();
}//void sort(...)


ShieldingSourceDisplay::Chi2Graphic::Chi2Graphic( Wt::WContainerWidget *parent )
  : Wt::Chart::WCartesianChart( parent ),
    m_nFitForPar( 0 ),
    m_showChi( true ),
    m_textPenColor( Wt::black )
{
  setPreferredMethod( WPaintedWidget::HtmlCanvas );
  LOAD_JAVASCRIPT( wApp, "shieldingSourceDisplay.cpp", "Chi2Graphic", wtjsShowChi2Info);
  
  auto interpsec = InterSpec::instance();
  if( interpsec )
    interpsec->useMessageResourceBundle( "ShieldingSourceDisplay" ); //JIC
  
  // We will use calcAndSetAxisRanges() to set the axis ranges, so turn off auto range.
  axis(Chart::YAxis).setAutoLimits( 0 );
  axis(Chart::YAxis).setRange( -1.0, 1.0 );
  
  axis(Chart::XAxis).setAutoLimits( 0 );
  axis(Chart::XAxis).setRange( 0, 3000.0 );
}

ShieldingSourceDisplay::Chi2Graphic::~Chi2Graphic()
{
}


void ShieldingSourceDisplay::Chi2Graphic::setColorsFromTheme( std::shared_ptr<const ColorTheme> theme )
{
  if( !theme )
    return;
  
  //if( !theme->foregroundLine.isDefault() )
  //  m_chartEnergyLineColor = theme->foregroundLine;
  
  if( !theme->spectrumChartText.isDefault() )
  {
    WPen txtpen(theme->spectrumChartText);
    this->setTextPen( txtpen );
    this->axis(Chart::XAxis).setTextPen( txtpen );
    this->axis(Chart::YAxis).setTextPen( txtpen );
    this->setTextPenColor( theme->spectrumChartText );
  }
  
  if( theme->spectrumChartBackground.isDefault() )
    this->setBackground( Wt::NoBrush );
  else
    this->setBackground( WBrush(theme->spectrumChartBackground) );
  
  if( (theme->spectrumChartMargins.isDefault() && !theme->spectrumChartBackground.isDefault()) )
  {
    //theme->spectrumChartBackground
  }else if( !theme->spectrumChartMargins.isDefault() )
  {
    //theme->spectrumChartMargins
  }
  
  if( !theme->spectrumAxisLines.isDefault() )
  {
    WPen defpen = this->axis(Chart::XAxis).pen();
    defpen.setColor( theme->spectrumAxisLines );
    this->axis(Chart::XAxis).setPen( defpen );
    this->axis(Chart::YAxis).setPen( defpen );
  }
}//ShieldingSourceDisplay::Chi2Graphic::setColorsFromTheme( theme )


void ShieldingSourceDisplay::Chi2Graphic::setNumFitForParams( unsigned int npar )
{
  m_nFitForPar = static_cast<int>( npar );
}

void ShieldingSourceDisplay::Chi2Graphic::setTextPenColor( const Wt::WColor &color )
{
  m_textPenColor = color;
}

void ShieldingSourceDisplay::Chi2Graphic::setShowChiOnChart( const bool show_chi )
{
  m_showChi = show_chi;
  removeSeries(1);
  removeSeries(2);
  setXSeriesColumn( 0 );
  
  Chart::WDataSeries series( (show_chi ? 1 : 2), Chart::PointSeries );
  series.setMarkerSize( 10.0 );
  
  addSeries( series );
}


void ShieldingSourceDisplay::Chi2Graphic::calcAndSetAxisRanges()
{
#if( WT_VERSION < 0x3030600 )
  WAbstractItemModel *theModel = model();
#else
  Wt::Chart::WAbstractChartModel *theModel = model();
  assert( theModel );
#endif

  if( !theModel )
    return;
  
  double ymin = DBL_MAX, ymax = -DBL_MAX, xmin = DBL_MAX, xmax = -DBL_MAX;
  
  const int nrow = theModel->rowCount();
  const int ycol = m_showChi ? 1 : 2;
  for( int row = 0; row < nrow; ++row )
  {
    try
    {
#if( WT_VERSION >= 0x3030600 )
      const double thischi = theModel->data(row,ycol);
      const double energy = theModel->data(row,0);
#else
      WModelIndex index = theModel->index(row,ycol);
      const double thischi = boost::any_cast<double>( theModel->data(index) );
      index = theModel->index(row,0);
      const double energy = boost::any_cast<double>( theModel->data(index) );
#endif
      xmin = std::min( xmin, energy );
      xmax = std::max( xmax, energy );
      
      ymin = std::min( ymin, thischi );
      ymax = std::max( ymax, thischi );
    }catch(...)
    {
    }
  }//for( int row = 0; row < nrow; ++row )
  
  if( !nrow || (ymin == DBL_MAX) || (ymax == -DBL_MAX) || (xmin == DBL_MAX) || (xmax == -DBL_MAX) )
  {
    // There are no points to display.
    xmin = 0.0;
    xmax = 3000.0;
    ymin = (m_showChi ? -1.0 : 0.0);
    ymax = (m_showChi ?  1.0 : 2.0);
  }else if( (nrow == 1) || (ymin == ymax) )
  {
    // There is one point to display.
    xmin -= 10.0;
    xmax += 10.0;
    
    //  For chi, we should be really close to 0.0, and for multiple should be close to 1.0
    ymin -= 0.25*ymin;
    ymax += 0.25*ymax;
    
    if( m_showChi )
    {
      ymin = std::min( ymin, -1.0 );
      ymax = std::max( ymax, 1.0 );
    }else
    {
      ymin = std::min( ymin, 0.5 );
      ymax = std::max( ymax, 1.5 );
    }
  }else
  {
    // There are many points to display
    const double xrange = std::max( xmax - xmin, 10.0 );
    xmin = ((xmin > 0.1*xrange) ? (xmin - 0.1*xrange) : 0.0);
    xmax += 0.1*xrange;
    
    
    const double yrange = std::max( ymax - ymin, 0.1 );
    ymin -= 0.25*yrange;
    ymax += 0.25*yrange;
    
    if( m_showChi )
    {
      ymin = std::min( ymin, -1.0 );
      ymax = std::max( ymax, 1.0 );
    }else
    {
      ymin = std::min( ymin, 0.5 );
      ymax = std::max( ymax, 1.5 );
    }
  }//if( no rows ) / else / else

  // Round up/down to the nearest 10 keV
  xmin = 10.0 * std::floor( 0.1*xmin );
  xmax = 10.0 * std::ceil( 0.1*xmax );
  
  // Round up/down to the nearest 0.1
  ymin = 0.1 * std::floor( 10.0*ymin );
  ymax = 0.1 * std::ceil( 10.0*ymax );
  
  //cout << "Setting xrange=" << xmin << ", " << xmax << "], yrange=[" << ymin << ", " << ymax << "]" << endl;
  
  const auto xLocation = m_showChi ? Chart::AxisValue::ZeroValue : Chart::AxisValue::MinimumValue;
  axis(Chart::XAxis).setLocation( xLocation );
  
  axis(Chart::XAxis).setRange( xmin, xmax );
  axis(Chart::YAxis).setRange( ymin, ymax );
}//void calcAndSetAxisRanges();



void ShieldingSourceDisplay::Chi2Graphic::calcAndSetAxisPadding( double yHeightPx )
{
  double ymin = axis(Chart::YAxis).minimum();
  double ymax = axis(Chart::YAxis).maximum();
  
  if( fabs(ymax - ymin) < 0.1 )
  {
    ymin = 0.0;
    ymax = 100.0;
  }//

//Calculate number of pixels we need to pad, for x axis to be at 45 pixels from
//  the bottom of the chart; if less than 10px, pad at least 10px, or at most
//  45px.
  const int topPadding = plotAreaPadding(Top);
  const double fracY = -ymin / (ymax - ymin);
  double pxToXAxis = (fracY >= 1.0) ? 0.0 : (45.0 - fracY*(yHeightPx - topPadding) ) / (1.0-fracY);
  if( IsInf(pxToXAxis) || IsNan(pxToXAxis) ) //JIC, but shouldnt be needed
    pxToXAxis = 10;
  pxToXAxis = std::floor( pxToXAxis + 0.5 );
  pxToXAxis = std::max( pxToXAxis, 10.0 );
  pxToXAxis = std::min( pxToXAxis, 45.0 );
  const int bottomPadding = static_cast<int>(pxToXAxis);
  
  if( bottomPadding != plotAreaPadding(Bottom) )
    setPlotAreaPadding( bottomPadding, Bottom );
  
  yHeightPx -= (topPadding + bottomPadding);
  
  {//Begin codeblock to determin min/max label values, following logic
   //  determined from round*125() and other places in WAxis.C
    const int numLabels = yHeightPx / 25.0;
    const double range = ymax - ymin;
    const double rangePerLabel = range / numLabels;

    double renderinterval = 1.0;
    double n = std::pow(10, std::floor(std::log10(rangePerLabel)));
    double msd = rangePerLabel / n;
    
    if (msd < 1.5)
      renderinterval = n;
    else if (msd < 3.3)
      renderinterval = 2*n;
    else if (msd < 7)
      renderinterval = 5*n;
    else
      renderinterval = 10*n;
    
    ymin += std::numeric_limits<double>::epsilon();
    ymax -= std::numeric_limits<double>::epsilon();
    ymin = renderinterval * std::floor( ymin / renderinterval);
    ymax = renderinterval * std::ceil( ymax / renderinterval);
  }//End codeblock to determin min/max label values

  
  const int oldleft = plotAreaPadding(Wt::Left);
  const WString minlabel = axis(Chart::Y1Axis).label(ymin);
  const WString maxlabel = axis(Chart::Y1Axis).label(ymax);
  const size_t maxnchars = minlabel.narrow().length();
  const size_t minnchars = maxlabel.narrow().length();
  const size_t nchars = std::max( minnchars, maxnchars );
  const int left = 32 + 6*std::max(static_cast<int>(nchars),2);
  
  if( left != oldleft )
    setPlotAreaPadding( left, Wt::Left );
}//void calcAndSetAxisPadding()


void ShieldingSourceDisplay::Chi2Graphic::paintEvent( WPaintDevice *device )
{
  calcAndSetAxisRanges();
  if( device )
    calcAndSetAxisPadding( device->height().toPixels() );
  Wt::Chart::WCartesianChart::paintEvent( device );
}//void paintEvent( Wt::WPaintDevice *paintDevice )


void ShieldingSourceDisplay::Chi2Graphic::paint( Wt::WPainter &painter,
                      const Wt::WRectF &rectangle ) const
{
  WCartesianChart::paint( painter, rectangle );

  //cout << "plot padding: [" << plotAreaPadding(Top) << ", " << plotAreaPadding(Right) << ", " << plotAreaPadding(Bottom) << ", " << plotAreaPadding(Left) << "]" << endl;
  
  //I think removing of the areas() is already done by
  //  WCartesianChart::paintEvent(...), but jic
  while( !areas().empty() )
    delete areas().front();
  
#if( WT_VERSION < 0x3030600 )
  WStandardItemModel *chi2Model = dynamic_cast<WStandardItemModel *>( model() );
#else
  auto proxyChi2Model = dynamic_cast<Wt::Chart::WStandardChartProxyModel *>( model() );
  assert( proxyChi2Model );
  assert( proxyChi2Model->sourceModel() );
  WStandardItemModel *chi2Model = dynamic_cast<WStandardItemModel *>( proxyChi2Model->sourceModel() );
#endif
    
  if( !chi2Model )  //prob never happen, but JIC
    return;
  
  const WPointF br = painter.window().bottomRight();
  const int width = br.x();

  const int nrow = chi2Model->rowCount();
  double chi2 = 0.0;
  
  for( int row = 0; row < nrow; ++row )
  {
    try
    {
      WModelIndex index = chi2Model->index(row,1);
      const double thischi = boost::any_cast<double>( chi2Model->data(index) );
      chi2 += thischi*thischi;
      

      WColor color;
      try
      {
        boost::any color_any = chi2Model->data(index, Wt::MarkerPenColorRole );
        color = boost::any_cast<WColor>(color_any);
        if( color.isDefault() )
          throw runtime_error("");
      }catch(...)
      {
        //I dont think we will ever get here, but JIC I guess.
        Wt::Chart::WChartPalette *pal = palette();
        if( pal )
          color = pal->brush(0).color();
        else
          color = WColor( Wt::darkRed );
      }//try / catch, get the color
      
      index = chi2Model->index(row,2);
      const double thisscale = boost::any_cast<double>( chi2Model->data(index) );
      
      index = chi2Model->index(row,0);
      double energy = boost::any_cast<double>( chi2Model->data(index) );
      
      
      const double yval = m_showChi ? thischi : thisscale;
      const WPointF pos = mapToDevice( energy, yval );
      
      index = chi2Model->index(row,3);
      const string nucname = Wt::asString( chi2Model->data(index) ).toUTF8();
      energy = ((100.0*energy+0.5)/100.0);
      
      char mouseoverjs[256];
      if( nucname.empty() || nucname.size() > 8 )
      {
        snprintf( mouseoverjs, sizeof(mouseoverjs),
                  "function(){Wt.WT.ShowChi2Info('%s','%.2f keV, &sigma;<sub>%i</sub>=%.2f, %.2fx model');}",
                  id().c_str(), energy, row, thischi, thisscale );
      }else
      {
        snprintf( mouseoverjs, sizeof(mouseoverjs),
                 "function(){Wt.WT.ShowChi2Info('%s','%s, %.2f keV, &sigma;<sub>%i</sub>=%.2f, %.2fx model');}",
                 id().c_str(), nucname.c_str(), energy, row, thischi, thisscale );
      }
      
      //We could potentially include more info here, but it will be kinda hard,
      //  so lets not worry about it now.
      WRectArea *area = new WRectArea( pos.x()-10.0, pos.y()-10.0, 20.0, 20.0 );
      const string mouseoutjs = "function(){$('#" + id() + "inf').remove();}";
      area->mouseWentOver().connect( string(mouseoverjs) );
      area->mouseWentOut().connect( mouseoutjs );
      //area->clicked()
      
      //I hate doing const_casts, but the Wt source code itself does this
      //  for adding image areas in Wt 3.3.2
      const_cast<Chi2Graphic*>(this)->addArea( area );
      
      if( !m_showChi )
      {
        index = chi2Model->index(row,4);
        const double scale_uncert = boost::any_cast<double>( chi2Model->data(index) );
        
        const WPointF upper_uncert = mapToDevice( energy, yval + scale_uncert );
        const WPointF lower_uncert = mapToDevice( energy, yval - scale_uncert );
        const WPen oldPen = painter.pen();
        painter.setPen( WPen(color) );
        painter.drawLine( upper_uncert, lower_uncert );
        painter.setPen( oldPen );
      }//if( !m_showChi )
    }catch( std::exception &e )
    {
      cerr << "Caught exception drawing Chi2 chart: " << e.what() << endl;
    }
  }//for( int row = 0; row < nrow; ++row )
  
  if( nrow > 0 && !IsNan(sqrt(chi2)) )
  {
    char buffer[64];
    //snprintf( buffer, sizeof(buffer), "&lt;dev&gt;=%.2g&sigma;", (sqrt(chi2)/nrow) );  //displays as literal tesxt, and not the symbols on android
    //snprintf( buffer, sizeof(buffer), "\x3c\xCF\x87\x3E=%.2g\xCF\x83", (sqrt(chi2)/nrow) ); //<χ>=13.2σ
    snprintf( buffer, sizeof(buffer), "\x3c\x64\x65\x76\x3E=%.2g\xCF\x83", (sqrt(chi2)/nrow) ); //<dev>=13.2σ
  
    WString text = WString::fromUTF8(buffer);
    const size_t msglen = SpecUtils::utf8_str_len( buffer, strlen(buffer) );
  
    const double rightPadding = static_cast<double>( plotAreaPadding(Right) );
    const double charwidth = 16;
    const double charheight = 15;
    const double x = width - charwidth*msglen - rightPadding - 5;
    const double y = plotAreaPadding(Top) + 5;
    const double twidth = charwidth*msglen;
    const double theight = charheight + 2;
    
    WPen oldPen = painter.pen();
    painter.setPen( WPen(m_textPenColor) );
    painter.drawText( x, y, twidth, theight, AlignRight, TextSingleLine, text );
    painter.setPen( oldPen );
  }//if( nrow > 0 && !IsNan(sqrt(chi2)) )

  const double yAxisMin = axis(Chart::YAxis).minimum();
  const double yAxisMax = axis(Chart::YAxis).maximum();
  
  if( !m_showChi && (yAxisMin < 1.0) && (yAxisMax > 1.0) )
  {
    const double xAxisMin = axis(Chart::XAxis).minimum();
    const double xAxisMax = axis(Chart::XAxis).maximum();
    
    WPointF left = mapToDevice( xAxisMin, 1.0 );
    WPointF right = mapToDevice( xAxisMax, 1.0 );
    
    //cout << "xAxisMin=" << xAxisMin << ", xAxisMax=" << xAxisMax << ", yAxisMin=" << yAxisMin << ", yAxisMax=" << yAxisMax << endl;
    //cout << "left={" << left.x() << "," << left.y() << "}, right={" << right.x() << "," << right.y() << "}" << endl;
    
    left.setY( left.y() + 0.5 );
    right.setY( right.y() + 0.5 );
    
    WPen pen( GlobalColor::lightGray );
    pen.setWidth( 1 );
    pen.setStyle( Wt::PenStyle::DashLine );
    const WPen oldPen = painter.pen();
    painter.setPen( pen );
    painter.drawLine( left, right );
    painter.setPen( oldPen );
  }//if( m_showChi )
  
  //Draw the y-axis label
  painter.rotate( -90 );
  painter.setPen( WPen() );
  painter.setBrush( WBrush() );
  
  
  
  const char *yaxistitle_key = m_showChi ? "x2g-yaxis-title-chi" : "x2g-yaxis-title-mult";
  const char * const yaxistooltip_key = m_showChi ? "x2g-tt-chi-axis" : "x2g-tt-scale-axis";
  
  WPen oldPen = painter.pen();
  painter.setPen( WPen(m_textPenColor) );
  painter.drawText( -0.45*painter.viewPort().height(), 0.0, 0.2, 0.1,
                     AlignCenter, WString::tr(yaxistitle_key) );
  painter.setPen( oldPen );
  
  const double yAxisWiddth = plotAreaPadding(Wt::Left);
  const double chartHeight = painter.viewPort().height();
  
  WRectArea *yAxisArea = new WRectArea( 0.0, 0.0, yAxisWiddth, chartHeight );
  yAxisArea->setToolTip( WString::tr(yaxistooltip_key) );
  const_cast<Chi2Graphic*>(this)->addArea( yAxisArea );
  
  painter.restore();
}//Chi2Graphic::paint(


pair<ShieldingSourceDisplay *,AuxWindow *> ShieldingSourceDisplay::createWindow( InterSpec *viewer )
{
  assert( viewer );
  
  // Its fine to not call `viewer->useMessageResourceBundle( "ShieldingSourceDisplay" );` yet
  //  - just needs to be done by render time
  
  AuxWindow *window = nullptr;
  ShieldingSourceDisplay *disp = nullptr;
  
  try
  {
    MaterialDB *matdb = viewer->materialDataBase();
    PeakModel *peakModel = viewer->peakModel();
    WSuggestionPopup *shieldSuggest = viewer->shieldingSuggester();
    
    disp = new ShieldingSourceDisplay( peakModel, viewer, shieldSuggest, matdb );
    window = new AuxWindow( WString::tr("window-title-act-shield-fit"),
                           (AuxWindowProperties::SetCloseable | AuxWindowProperties::EnableResize) );
    // We have to set minimum size before calling setResizable, or else Wt's Resizable.js functions
    //  will be called first, which will then default to using the initial size as minimum allowable
    if( !viewer->isPhone() )
    {
      window->setMinimumSize( 800, 480 );
      window->setResizable( true );
    }
    
    window->contents()->setOffsets(WLength(0,WLength::Pixel));
    window->stretcher()->addWidget( disp, 0, 0 );
    window->stretcher()->setContentsMargins(0,0,0,0);
  //    window->footer()->resize(WLength::Auto, WLength(50.0));
      
    WPushButton *closeButton = window->addCloseButtonToFooter();
    closeButton->clicked().connect(window, &AuxWindow::hide);
      
    AuxWindow::addHelpInFooter( window->footer(), "activity-shielding-dialog" );
    
    window->rejectWhenEscapePressed();
      
    //Should take lock on m_dataMeasurement->mutex_
    
    shared_ptr<SpecMeas> meas = viewer->measurment(SpecUtils::SpectrumType::Foreground);
    rapidxml::xml_document<char> *shield_source = nullptr;
    if( meas )
      shield_source = meas->shieldingSourceModel();
      
    if( shield_source && shield_source->first_node() )
    {
      //string msg = "Will try to deserailize: \n";
      //rapidxml::print( std::back_inserter(msg), *shield_source->first_node(), 0 );
      //cout << msg << endl << endl;
      try
      {
        disp->deSerialize( shield_source->first_node() );
      }catch( std::exception &e )
      {
        string xmlstring;
        rapidxml::print(std::back_inserter(xmlstring), *shield_source, 0);
        stringstream debugmsg;
        debugmsg << "Error loading Shielding/Source model: "
                    "\n\tError Message: " << e.what()
                 << "\n\tModel XML: " << xmlstring;
  #if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, debugmsg.str().c_str() );
  #else
        cerr << debugmsg.str() << endl;
  #endif
        passMessage( "There was an error loading the previous shielding/source model for this file - model state is suspect!",
                      WarningWidget::WarningMsgHigh );
      }
    }//if( shield_source )
    
  //    m_shieldingSourceFitWindow->resizeScaledWindow( 0.75, 0.75 );
      
    double windowWidth = viewer->renderedWidth();
    double windowHeight = viewer->renderedHeight();
    
  //    double footerheight = m_shieldingSourceFitWindow->footer()->height().value();
  //    m_shieldingSourceFitWindow->setMinimumSize( WLength(200), WLength(windowHeight) );
      
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

      //Give the m_shieldingSourceFitWindow a hint about what size it will be
      //  rendered at so it can decide what widgets should be rendered - acounting
      //  for borders and stuff (roughly)
      disp->initialSizeHint( windowWidth - 12, windowHeight - 28 - 50 );
    }else if( !viewer->isPhone() )
    {
      //When loading an application state that is showing this window, we may
      //  not know the window size (e.g., windowWidth==windowHeight==0), so
      //  instead skip giving the initial size hint, and instead size things
      //  client side (maybe we should just do this always?)
      window->resizeScaledWindow( 0.90, 0.90 );
    }
        
#if( INCLUDE_ANALYSIS_TEST_SUITE )
    WPushButton *setTruth = new WPushButton( "Set Truth Values" );
    window->footer()->insertWidget( 0, setTruth );
    setTruth->clicked().connect( disp, &ShieldingSourceDisplay::showInputTruthValuesWindow );
#endif
    
    window->centerWindow();
  //   m_shieldingSourceFitWindow->contents()->  setHeight(WLength(windowHeight));
    window->finished().connect( viewer, &InterSpec::closeShieldingSourceFit );
      
    window->WDialog::setHidden(false);
    window->show();
    window->centerWindow();
  }catch( std::exception &e )
  {
    passMessage( WString::tr("ssd-error-create-tool").tr(e.what()), WarningWidget::WarningMsgHigh );
    
    if( disp )
      delete disp;
    disp = nullptr;
    
    if( window )
      AuxWindow::deleteAuxWindow( window );
    window = nullptr;
  }//try / catch
  
  return make_pair( disp, window );
}//createWindow( InterSpec *viewer  )




ShieldingSourceDisplay::ShieldingSourceDisplay( PeakModel *peakModel,
                                                InterSpec *specViewer,
                                                WSuggestionPopup *matSuggest,
                                                MaterialDB *materialDB,
                                                WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_chi2ChartNeedsUpdating( true ),
    m_width( 0 ),
    m_height( 0 ),
    m_nResizeSinceHint( 0 ),
    m_modifiedThisForeground( false ),
    m_peakModel( peakModel ),
    m_specViewer( specViewer ),
    m_sourceModel( nullptr ),
    m_peakView( nullptr ),
    m_sourceView( nullptr ),
    m_detectorDisplay( nullptr ),
    m_distanceLabel( nullptr ),
    m_prevDistStr(),
    m_distanceEdit( nullptr ),
    m_addMaterialShielding( nullptr ),
    m_addGenericShielding( nullptr ),
    m_layout( nullptr ),
    m_addItemMenu( nullptr ),
#if( USE_DB_TO_STORE_SPECTRA )
    m_saveAsNewModelInDb( nullptr ),
#endif
    m_materialSuggest( matSuggest ),
    m_shieldingSelects( nullptr ),
    m_geometryLabel( nullptr ),
    m_prevGeometry( GammaInteractionCalc::GeometryType::Spherical ),
    m_geometrySelect( nullptr ),
    m_fixedGeometryTxt( nullptr ),
    m_showChi2Text( nullptr ),
    m_chi2Model( nullptr ),
    m_chi2Graphic( nullptr ),
    m_multiIsoPerPeak( nullptr ),
    m_clusterWidth( nullptr ),
    m_attenForAir( nullptr ),
    m_backgroundPeakSub( nullptr ),
    m_sameIsotopesAge( nullptr ),
    m_decayCorrect( nullptr ),
    m_showChiOnChart( nullptr ),
    m_optionsDiv( nullptr ),
    m_photopeak_cluster_sigma( 1.25 ),
    m_multithread_computation( true ),
    m_showLog( nullptr ),
    m_logDiv( nullptr ),
    m_calcLog{},
    m_peakCalcLogInfo{},
    m_modelUploadWindow( nullptr ),
#if( USE_DB_TO_STORE_SPECTRA )
    m_modelDbBrowseWindow( nullptr ),
    m_modelDbSaveWindow( nullptr ),
#endif
    m_materialDB( materialDB ),
    m_fitModelButton( nullptr ),
    m_fitProgressTxt( nullptr ),
    m_cancelfitModelButton( nullptr )
{
  wApp->useStyleSheet( "InterSpec_resources/ShieldingSourceDisplay.css" );
  
  addStyleClass( "ShieldingSourceDisplay" );
  
  assert( m_specViewer );
  m_specViewer->useMessageResourceBundle( "ShieldingSourceDisplay" );
      
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", m_specViewer );
  
  setLayoutSizeAware( true );
  const bool isotopesHaveSameAge = true;
  m_sourceModel = new SourceFitModel( m_peakModel, isotopesHaveSameAge, this );
  m_peakView = new RowStretchTreeView();
  m_peakView->setRootIsDecorated	(	false); //makes the tree look like a table! :)
  
  m_peakView->setAlternatingRowColors( true );
  m_peakView->setEditTriggers( WAbstractItemView::SingleClicked | WAbstractItemView::DoubleClicked );
  m_peakView->setModel( m_peakModel );
  m_peakView->addStyleClass( "PeakView" );

  for( PeakModel::Columns col = PeakModel::Columns(0);
       col < PeakModel::kNumColumns;
       col = PeakModel::Columns(col+1) )
  {
    switch( col )
    {
      case PeakModel::kMean: case PeakModel::kUseForShieldingSourceFit:
      case PeakModel::kPhotoPeakEnergy: case PeakModel::kIsotope:
        m_peakView->setColumnHidden( col, false );
      break;
      default:
        m_peakView->setColumnHidden( col, true );
      break;
    }//switch( col )
  }//for( loop over peak columns )

  WItemDelegate *dblDelagate = new WItemDelegate( this );
  dblDelagate->setTextFormat( "%.2f" );
  m_peakView->setItemDelegateForColumn( PeakModel::kMean, dblDelagate );

  m_peakView->setColumnWidth( PeakModel::kMean, WLength(7.5,WLength::FontEx) );
  m_peakView->setColumnWidth( PeakModel::kPhotoPeakEnergy, WLength(10,WLength::FontEx) );
  m_peakView->setColumnWidth( PeakModel::kIsotope, WLength(9,WLength::FontEx) );
  m_peakView->setColumnWidth( PeakModel::kUseForShieldingSourceFit, WLength(7,WLength::FontEx) );


  PhotopeakDelegate *nuclideDelegate = new PhotopeakDelegate( PhotopeakDelegate::NuclideDelegate, true, m_peakView );
  m_peakView->setItemDelegateForColumn( PeakModel::kIsotope, nuclideDelegate );


  m_sourceView = new RowStretchTreeView();
  m_sourceView->setRootIsDecorated	(	false); //makes the tree look like a table! :)
  

  m_sourceView->setModel( m_sourceModel );
  m_sourceView->setSortingEnabled( true );
  m_sourceView->setAlternatingRowColors( true );
  m_sourceView->addStyleClass( "SourceView" );

  for( SourceFitModel::Columns col = SourceFitModel::Columns(0);
       col < SourceFitModel::kNumColumns;
       col = SourceFitModel::Columns(col+1) )
  {
    switch( col )
    {
      case SourceFitModel::kActivity:
        //need make custom delegate
      break;
      case SourceFitModel::kAge:
        //need to make custom delegate
      break;

      case SourceFitModel::kActivityUncertainty:
      case SourceFitModel::kAgeUncertainty:
        m_sourceView->setColumnHidden( static_cast<int>( col ), true );
        break;
        
      case SourceFitModel::kFitActivity:
      case SourceFitModel::kFitAge:
      case SourceFitModel::kIsotope:
      case SourceFitModel::kIsotopeMass:
#if( INCLUDE_ANALYSIS_TEST_SUITE )
      case SourceFitModel::kTruthActivity: case SourceFitModel::kTruthActivityTolerance:
      case SourceFitModel::kTruthAge: case SourceFitModel::kTruthAgeTolerance:
#endif
      case SourceFitModel::kNumColumns:
      break;
    }//case( col )
  }//for( loop over SourceFitModel columns )


  m_sourceView->setColumnWidth( SourceFitModel::kActivity, WLength(150,WLength::Pixel) );
  m_sourceView->setColumnWidth( SourceFitModel::kAge, WLength(9,WLength::FontEx) );
  m_sourceView->setColumnWidth( SourceFitModel::kFitAge, WLength(10,WLength::FontEx) );
  m_sourceView->setColumnWidth( SourceFitModel::kFitActivity, WLength(10,WLength::FontEx) );
  m_sourceView->setColumnWidth( SourceFitModel::kIsotope, WLength(9,WLength::FontEx) );
  m_sourceView->setColumnWidth( SourceFitModel::kIsotopeMass, WLength(9,WLength::FontEx) );

  m_sourceView->setColumnWidth( SourceFitModel::kActivityUncertainty, WLength(10,WLength::FontEx) );
  m_sourceView->setColumnWidth( SourceFitModel::kAgeUncertainty, WLength(10,WLength::FontEx) );

#if( INCLUDE_ANALYSIS_TEST_SUITE )
  m_sourceView->setColumnWidth( SourceFitModel::kTruthActivity, WLength(13,WLength::FontEx) );
  m_sourceView->setColumnWidth( SourceFitModel::kTruthActivityTolerance, WLength(14,WLength::FontEx) );
  m_sourceView->setColumnWidth( SourceFitModel::kTruthAge, WLength(13,WLength::FontEx) );
  m_sourceView->setColumnWidth( SourceFitModel::kTruthAgeTolerance, WLength(14,WLength::FontEx) );
#endif
  
  
  m_detectorDisplay = new DetectorDisplay( m_specViewer, m_specViewer->fileManager()->model() );
  m_detectorDisplay->setInline( true );

  Wt::WPushButton *addItemMenubutton = new WPushButton();
  addItemMenubutton->setStyleClass( "RoundMenuIcon InvertInDark" );
  addItemMenubutton->clicked().preventPropagation();
  m_addItemMenu = new PopupDivMenu( addItemMenubutton, PopupDivMenu::TransientMenu );

  //this validates floating point numbers followed by a distance unit
  WRegExpValidator *distValidator = new WRegExpValidator( PhysicalUnits::sm_distanceUnitOptionalRegex, this );
  distValidator->setFlags( Wt::MatchCaseInsensitive );
  m_distanceLabel = new WLabel( WString("{1}:").arg(WString::tr("Distance")) );
  
  m_distanceEdit = new WLineEdit( "100 cm" );
  
  m_distanceEdit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_distanceEdit->setAttributeValue( "autocorrect", "off" );
  m_distanceEdit->setAttributeValue( "spellcheck", "off" );
#endif
  m_distanceEdit->setTextSize( 5 );

  m_distanceLabel->setBuddy( m_distanceEdit );
  m_distanceEdit->setValidator( distValidator );
  HelpSystem::attachToolTipOn( m_distanceEdit, WString::tr("ssd-tt-distance"), showToolTips );

  m_geometryLabel = new WLabel( WString("{1}:").arg(WString::tr("Geometry")) );
  m_geometrySelect = new WComboBox();
  
  // We want the current index of m_geometrySelect to correspond to the GeometryType enum value
  for( GeometryType type = GeometryType(0);
      type != GeometryType::NumGeometryType;
      type = GeometryType( static_cast<int>(type) + 1) )
  {
    const char *lbl_key = "";
    switch( GeometryType(type) )
    {
      case GeometryType::Spherical:      lbl_key = "ssd-geom-point";    break;
      case GeometryType::CylinderEndOn:  lbl_key = "ssd-geom-cyl-end";  break;
      case GeometryType::CylinderSideOn: lbl_key = "ssd-geom-cyl-side"; break;
      case GeometryType::Rectangular:    lbl_key = "ssd-rectangular";   break;
      case GeometryType::NumGeometryType: assert( 0 ); break;
      default: assert( 0 ); throw runtime_error(""); break;
    }//switch( GeometryType(type) )
    
    m_geometrySelect->addItem( WString::tr(lbl_key) );
  }//for( int type = 0; type < 10; ++type )
  
  m_prevGeometry = GeometryType::Spherical;
  m_geometrySelect->setCurrentIndex( static_cast<int>(GeometryType::Spherical) );
  m_geometrySelect->changed().connect( this, &ShieldingSourceDisplay::handleGeometryTypeChange );

  m_fixedGeometryTxt = new WText( WString::tr("ssd-fixed-geom-notification") );
  m_fixedGeometryTxt->hide();

  m_shieldingSelects = new WContainerWidget();
  m_shieldingSelects->setStyleClass( "ShieldingSelectContainer" );

  HelpSystem::attachToolTipOn( m_shieldingSelects, WString::tr("ssd-tt-shieldings"), showToolTips );

  WLabel *addShieldingLabel = new WLabel( WString::tr("ssd-add-shield-label") );
  m_addMaterialShielding = new WPushButton( WString::tr("Material") );
  HelpSystem::attachToolTipOn( m_addMaterialShielding, WString::tr("ssd-tt-add-shield"),
                              showToolTips, HelpSystem::ToolTipPosition::Top  );
  m_addMaterialShielding->setIcon( "InterSpec_resources/images/shield_white.png" );
  m_addMaterialShielding->clicked().connect( this,
                                      &ShieldingSourceDisplay::doAddShielding );
  
  m_addGenericShielding = new WPushButton( WString::tr("ssd-generic-btn") );
  HelpSystem::attachToolTipOn( m_addGenericShielding, WString::tr("ssd-tt-generic"),
                              showToolTips , HelpSystem::ToolTipPosition::Top );
  m_addGenericShielding->setIcon( "InterSpec_resources/images/atom_white.png" );
  m_addGenericShielding->clicked().connect( this,
                                     &ShieldingSourceDisplay::addGenericShielding );
  
  
  m_fitModelButton = new WPushButton( WString::tr("ssd-perform-fit-btn") );
  m_fitModelButton->clicked().connect( boost::bind(&ShieldingSourceDisplay::doModelFit, this, true) );

  m_fitProgressTxt = new WText();
  m_fitProgressTxt->hide();
  
  m_cancelfitModelButton = new WPushButton( WString::tr("ssd-cancel-fit") );
  m_cancelfitModelButton->clicked().connect( boost::bind( &ShieldingSourceDisplay::cancelModelFit, this ) );
  m_cancelfitModelButton->hide();
  
  m_showLog = m_addItemMenu->addMenuItem( WString::tr("ssd-mi-calc-log") );
  m_showLog->triggered().connect( this, &ShieldingSourceDisplay::showCalcLog );
  m_showLog->disable();

  PopupDivMenuItem *item = NULL;
//  PopupDivMenuItem *item = m_addItemMenu->addMenuItem( "Test Serialization" );
//  item->triggered().connect( this, &ShieldingSourceDisplay::testSerialization );

  item = m_addItemMenu->addMenuItem( WString::tr("ssd-mi-import-model") );
  item->triggered().connect( this, &ShieldingSourceDisplay::startModelUpload );
  
  StringDownloadResource *xmlResource = new StringDownloadResource( this );
  item = m_addItemMenu->addMenuItem( WString::tr("ssd-mi-export-model") );
  item->setLink( WLink( xmlResource ) );
  item->setLinkTarget(Wt::TargetNewWindow);
  
#if( ANDROID )
  // Using hacked saving to temporary file in Android, instead of via network download of file.
  item->clicked().connect( std::bind([xmlResource](){
    android_download_workaround(xmlResource, "fit_model.xml");
  }) );
#endif //ANDROID
  
  
#if( USE_DB_TO_STORE_SPECTRA )
  item = m_addItemMenu->addMenuItem( WString::tr("ssd-mi-from-db") );
  item->triggered().connect( this,
                          &ShieldingSourceDisplay::startBrowseDatabaseModels );
  
  item = m_addItemMenu->addMenuItem( WString::tr("ssd-mi-save-to-db") );
  item->triggered().connect(
                boost::bind( &ShieldingSourceDisplay::startSaveModelToDatabase,
                             this, false) );
  
  m_saveAsNewModelInDb = m_addItemMenu->addMenuItem( WString::tr("ssd-mi-clone-db-entry") );
  m_saveAsNewModelInDb->triggered().connect( this,
                            &ShieldingSourceDisplay::saveCloneModelToDatabase );
  m_saveAsNewModelInDb->disable();
#endif //#if( USE_DB_TO_STORE_SPECTRA )
  
  m_showChi2Text = new WText( WString::tr("ssd-to-small-for-chart"), XHTMLText );
  m_showChi2Text->setInline( false );
  m_showChi2Text->hide();

  m_chi2Graphic = new Chi2Graphic();
    
  m_chi2Model = new WStandardItemModel( 0, 6, parent );

#if( WT_VERSION < 0x3030600 )
  m_chi2Graphic->setModel( m_chi2Model );
#else
  auto proxyModel = new Wt::Chart::WStandardChartProxyModel( m_chi2Model, this );
  m_chi2Graphic->setModel( m_chi2Model );
#endif

  m_chi2Graphic->setPlotAreaPadding( 12, Right );
  m_chi2Graphic->setPlotAreaPadding(  2, Top );
  
  // Left and bottom paddings will be set by Chi2Graphic::calcAndSetAxisPadding(...)
  //m_chi2Graphic->setPlotAreaPadding( 50, Left );
  //m_chi2Graphic->setPlotAreaPadding( 40, Bottom );
  //m_chi2Graphic->setAutoLayoutEnabled();
  
  m_chi2Graphic->addSeries( Chart::WDataSeries(1, Chart::PointSeries) );
  m_chi2Graphic->setXSeriesColumn( 0 );
  m_chi2Graphic->axis(Chart::XAxis).setTitle( WString::tr("Energy (keV)") );
  
  WFont font( WFont::Default );
  font.setSize( WFont::Small );
  m_chi2Graphic->axis(Chart::YAxis).setTitleFont(font);
  m_chi2Graphic->axis(Chart::XAxis).setTitleFont(font);
  
  WFont labelFont( WFont::Default );
  labelFont.setSize( WFont::Medium );
  m_chi2Graphic->axis(Chart::YAxis).setLabelFont(font);
  m_chi2Graphic->axis(Chart::XAxis).setLabelFont(font);
  
  m_chi2Graphic->axis(Chart::XAxis).setLocation( Chart::AxisValue::ZeroValue );
  m_chi2Graphic->axis(Chart::YAxis).setLocation( Chart::AxisValue::MinimumValue );

  m_chi2Graphic->axis(Chart::XAxis).setScale( Chart::LinearScale );
  m_chi2Graphic->axis(Chart::YAxis).setScale( Chart::LinearScale );

  m_chi2Graphic->setType( Chart::ScatterPlot );
//  m_chi2Graphic->setMinimumSize( WLength(200), WLength(175) );
  
  
  //We should check the color theme for colors
  m_specViewer->colorThemeChanged().connect( m_chi2Graphic, &Chi2Graphic::setColorsFromTheme );
  m_chi2Graphic->setColorsFromTheme( m_specViewer->getColorTheme() );
  
  
  //The next line is kinda inefficient because if all that changed was
  //  fit activity or fit age, then we dont really need to update the chi2 chart
  m_sourceModel->dataChanged().connect( this, &ShieldingSourceDisplay::updateChi2Chart );

  //XXX -
  //  When the below is connected, the order the isotpes are added to
  //  m_sourceModel is somewhat off, so there is an issue updating the
  //  m_chi2Graphic because of a mismatch, in what isotopes are actually in
  //  m_sourceModel
  m_sourceModel->rowsInserted().connect( this, &ShieldingSourceDisplay::updateChi2Chart );
  
  //XXX - this next line causes a runtime_exception to be raised under some
  //      circumstances, causing the Chi2 chart to not be updated - leaving for
  //      now, but I'm not really sure as to the sequence of events that causes
  //      this to happen (model gets out of sync with the chi2 functions)
  m_sourceModel->rowsRemoved().connect( this, &ShieldingSourceDisplay::updateChi2Chart );

  m_sourceModel->rowsInserted().connect( this, &ShieldingSourceDisplay::addSourceIsotopesToShieldings );
  m_sourceModel->rowsAboutToBeRemoved().connect( this, &ShieldingSourceDisplay::removeSourceIsotopesFromShieldings );

  m_sourceModel->layoutChanged().connect( this, &ShieldingSourceDisplay::updateChi2Chart );
  
  m_distanceEdit->changed().connect( this, &ShieldingSourceDisplay::handleUserDistanceChange );
  m_distanceEdit->enterPressed().connect( this, &ShieldingSourceDisplay::handleUserDistanceChange );
  
  m_specViewer->detectorChanged().connect( boost::bind( &ShieldingSourceDisplay::handleDetectorChanged, this, boost::placeholders::_1 ) );
  m_specViewer->detectorModified().connect( boost::bind( &ShieldingSourceDisplay::handleDetectorChanged, this, boost::placeholders::_1 ) );
  
  m_showChiOnChart = new SwitchCheckbox( "Mult.", "&chi;" );
  m_showChiOnChart->setChecked();
  HelpSystem::attachToolTipOn( m_showChiOnChart, WString::tr("ssd-tt-chi2-switch"),
                                  showToolTips, HelpSystem::ToolTipPosition::Right );
  //m_showChiOnChart->setToolTip( WString::tr("ssd-tt-chi2-switch") );
  m_showChiOnChart->checked().connect( this, &ShieldingSourceDisplay::showGraphicTypeChanged );
  m_showChiOnChart->unChecked().connect( this, &ShieldingSourceDisplay::showGraphicTypeChanged );
  
  WContainerWidget *allPeaksDiv = new WContainerWidget();
  WCheckBox *allpeaks = new WCheckBox( WString::tr("ssd-cb-all-peaks"), allPeaksDiv );
  allpeaks->addStyleClass( "CbNoLineBreak AllPeaksCb" );
  allpeaks->setTristate( true );
  allpeaks->changed().connect( boost::bind( &ShieldingSourceDisplay::toggleUseAll, this, allpeaks ) );
  m_peakModel->dataChanged().connect( boost::bind( &ShieldingSourceDisplay::updateAllPeaksCheckBox, this, allpeaks ) );
  updateAllPeaksCheckBox( allpeaks ); //initialize
  
  
  m_optionsDiv = new WGroupBox( WString::tr("ssd-options-title") );
  m_optionsDiv->addStyleClass( "FitOptions" );
      
  //The ToolTip of WCheckBoxes is a bit finicky, and only works over the
  //  checkbox itself, so lets make it work over the label to, via lineDiv
  WContainerWidget *lineDiv = new WContainerWidget( m_optionsDiv );
  lineDiv->addStyleClass( "FitOptionsRow" );
  m_multiIsoPerPeak = new WCheckBox( WString::tr("ssd-multi-iso-per-peak"), lineDiv );
  m_multiIsoPerPeak->addStyleClass( "CbNoLineBreak" );
  HelpSystem::attachToolTipOn( lineDiv, WString::tr("ssd-tt-multi-iso-per-peak"),
                                      showToolTips, HelpSystem::ToolTipPosition::Right );
  m_multiIsoPerPeak->setChecked();
  m_multiIsoPerPeak->checked().connect( this, &ShieldingSourceDisplay::multiNucsPerPeakChanged );
  m_multiIsoPerPeak->unChecked().connect( this, &ShieldingSourceDisplay::multiNucsPerPeakChanged );
  
  lineDiv = new WContainerWidget( m_optionsDiv );
  lineDiv->addStyleClass( "NumInputOptLine FitOptionsRow" );
  WLabel *clusterWidthLabel = new WLabel( WString::tr("ssd-cluster-width-label"), lineDiv );
  lineDiv->addWidget( clusterWidthLabel );
  m_clusterWidth = new NativeFloatSpinBox( lineDiv );
  m_clusterWidth->addStyleClass( "CbNoLineBreak" );
  m_clusterWidth->setRange( 0.0f, 10.0f );
  m_clusterWidth->setWidth( 50 );
  m_clusterWidth->setSpinnerHidden();
  clusterWidthLabel->setBuddy( m_clusterWidth );
  HelpSystem::attachToolTipOn( lineDiv, WString::tr("ssd-tt-cluster-width"),
                                showToolTips, HelpSystem::ToolTipPosition::Right );
  m_clusterWidth->setValue( m_photopeak_cluster_sigma );
  m_clusterWidth->valueChanged().connect( this, &ShieldingSourceDisplay::clusterWidthChanged );
  m_clusterWidth->setDisabled( !m_multiIsoPerPeak->isChecked() );
  if( m_clusterWidth->label() )
    m_clusterWidth->label()->setDisabled( !m_multiIsoPerPeak->isChecked() );

  lineDiv = new WContainerWidget( m_optionsDiv );
  lineDiv->addStyleClass( "FitOptionsRow" );
  m_attenForAir = new WCheckBox( WString::tr("ssd-cb-atten-for-air"), lineDiv );
  m_attenForAir->addStyleClass( "CbNoLineBreak" );
  //lineDiv->setToolTip( WString::tr("ssd-tt-atten-for-air") );
  HelpSystem::attachToolTipOn( lineDiv, WString::tr("ssd-tt-atten-for-air"),
                              showToolTips, HelpSystem::ToolTipPosition::Right );
  m_attenForAir->setChecked();
  m_attenForAir->checked().connect( this, &ShieldingSourceDisplay::attenuateForAirChanged );
  m_attenForAir->unChecked().connect( this, &ShieldingSourceDisplay::attenuateForAirChanged );
  
  
  lineDiv = new WContainerWidget( m_optionsDiv );
  lineDiv->addStyleClass( "FitOptionsRow" );
  m_backgroundPeakSub = new WCheckBox( WString::tr("ssd-cb-sub-back-peaks"), lineDiv );
  m_backgroundPeakSub->addStyleClass( "CbNoLineBreak" );
  lineDiv->setToolTip( WString::tr("ssd-tt-sub-back-peaks") );
  HelpSystem::attachToolTipOn( lineDiv, WString::tr("ssd-tt-sub-back-peaks"),
                                  showToolTips, HelpSystem::ToolTipPosition::Right );
  m_backgroundPeakSub->checked().connect( this, &ShieldingSourceDisplay::backgroundPeakSubChanged );
  m_backgroundPeakSub->unChecked().connect( this, &ShieldingSourceDisplay::backgroundPeakSubChanged );
  
  
  lineDiv = new WContainerWidget( m_optionsDiv );
  lineDiv->addStyleClass( "FitOptionsRow" );
  m_sameIsotopesAge = new WCheckBox( WString::tr("ssd-cb-same-el-same-age"), lineDiv );
  m_sameIsotopesAge->addStyleClass( "CbNoLineBreak" );
  //lineDiv->setToolTip( WString::tr("ssd-tt-same-el-same-age") );
  HelpSystem::attachToolTipOn( lineDiv, WString::tr("ssd-tt-same-el-same-age"),
                                showToolTips, HelpSystem::ToolTipPosition::Right );
  bool account_for_decay_during_meas = false;
  m_sameIsotopesAge->setChecked( account_for_decay_during_meas );
  m_sameIsotopesAge->checked().connect( this, &ShieldingSourceDisplay::sameIsotopesAgeChanged );
  m_sameIsotopesAge->unChecked().connect( this, &ShieldingSourceDisplay::sameIsotopesAgeChanged );

      
  lineDiv = new WContainerWidget( m_optionsDiv );
  lineDiv->addStyleClass( "FitOptionsRow" );
  m_decayCorrect = new WCheckBox( WString::tr("ssd-cb-corr-for-decay"), lineDiv );
  m_decayCorrect->addStyleClass( "CbNoLineBreak" );
  //lineDiv->setToolTip( WString::tr("ssd-tt-corr-for-decay") );
  HelpSystem::attachToolTipOn( lineDiv, WString::tr("ssd-tt-corr-for-decay"),
                                showToolTips, HelpSystem::ToolTipPosition::Right );
  
  // We'll set decay correction on by default, only if the measurement time is at least
  //  0.5% (arbitrary) of any of the nuclide half-lives.
  //  This isnt actually super great, because it wont catch if nuclides are added later,
  //  or something.
  bool decayCorrectActivity = false;
  const auto foreground = m_specViewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  const float real_time = foreground ? foreground->real_time() : 0.0f;
  if( real_time > 0.0f )
  {
    for( int nuc_index = 0; nuc_index < m_sourceModel->numNuclides(); ++nuc_index )
    {
      const SandiaDecay::Nuclide * const nuc = m_sourceModel->nuclide(nuc_index);
      decayCorrectActivity |= (nuc && (real_time > 0.005*nuc->halfLife));
    }//for( loop over nuclides )
  }//if( we have a foreground with non-zero real time )
      
  m_decayCorrect->setChecked( decayCorrectActivity );
  m_decayCorrect->checked().connect( this, &ShieldingSourceDisplay::decayCorrectChanged );
  m_decayCorrect->unChecked().connect( this, &ShieldingSourceDisplay::decayCorrectChanged );
      
  
  WContainerWidget *detectorDiv = new WContainerWidget();
  WGridLayout *detectorLayout = new WGridLayout();
  detectorDiv->setLayout( detectorLayout );
  
  
  WContainerWidget* smallerContainer = new WContainerWidget();
  WGridLayout *smallLayout = new WGridLayout();
  smallerContainer->setLayout(smallLayout);
  
  smallLayout->addWidget( m_distanceLabel,         0, 0, AlignRight | AlignMiddle );
  smallLayout->addWidget( m_distanceEdit,          0, 1, 1, 2);
  smallLayout->addWidget( m_geometryLabel,         1, 0, AlignRight | AlignMiddle );
  smallLayout->addWidget( m_geometrySelect,        1, 1, 1, 2);
  smallLayout->addWidget( m_fixedGeometryTxt,      2, 0, 1, 3, AlignCenter );
  smallLayout->addWidget( addShieldingLabel,       3, 0, AlignRight | AlignMiddle );
  smallLayout->addWidget( m_addMaterialShielding,  3, 1);
  smallLayout->addWidget( m_addGenericShielding,   3, 2);
  smallLayout->setContentsMargins( 0, 5, 0, 5 );
  smallerContainer->setPadding(0);
  
  //m_geometryLabel->setText( "Shield Geometry" );
  HelpSystem::attachToolTipOn( m_geometrySelect, WString::tr("ssd-tt-geometry"), showToolTips );
  
  //---------------
  
  WContainerWidget *peakDiv = new WContainerWidget();
  
  WGridLayout *peakGrid = new Wt::WGridLayout();
  peakGrid->setContentsMargins( 0, 0, 0, 0 );
  peakGrid->setVerticalSpacing( 0 );
  peakGrid->setHorizontalSpacing( 0 );
  peakGrid->setRowStretch( 0, 1 );
  peakGrid->setColumnStretch( 0, 1 );
  peakDiv->setLayout( peakGrid );
  peakGrid->addWidget( m_peakView, 0, 0 );
  peakGrid->addWidget( allPeaksDiv, 1, 0 );
      
  
  m_layout = new WGridLayout();
  m_layout->setContentsMargins( 0, 0, 0, 0 );
  setLayout( m_layout );
  
  WContainerWidget *peaksDiv = new WContainerWidget();
  WGridLayout *peaksLayout = new WGridLayout();
  peaksDiv->setLayout(peaksLayout);
  peaksLayout->addWidget( peakDiv, 0, 0 );
  peaksLayout->setRowStretch( 0, 1 );
  peaksLayout->setVerticalSpacing( 0 );
  peaksLayout->setHorizontalSpacing( 0 );
  peaksLayout->setContentsMargins( 0, 0, 0, 0 );
      
  WContainerWidget *sourceDiv = new WContainerWidget();
  WGridLayout *sourceGrid = new Wt::WGridLayout();
  sourceGrid->setRowStretch(0, 1);
  sourceGrid->setColumnStretch(0, 1);
  sourceDiv->setLayout(sourceGrid);
  sourceGrid->addWidget(m_sourceView,0,0);
  sourceGrid->setContentsMargins( 0, 0, 0, 0 );
  sourceGrid->setVerticalSpacing( 0 );
  sourceGrid->setHorizontalSpacing( 0 );
      
  
  if( m_specViewer->isPhone() )
  {
    //phone layout
    detectorLayout->addWidget( smallerContainer,           0, 0);
    detectorLayout->addWidget( m_shieldingSelects,         1, 0);
    
    detectorLayout->setRowStretch( 1, 1 );
    detectorLayout->setColumnStretch( 0, 1 );
    detectorLayout->setHorizontalSpacing( 0 );
    detectorLayout->setVerticalSpacing( 0 );
    detectorLayout->setContentsMargins( 1, 1, 1, 1 );
    
    WTabWidget *tab = new WTabWidget();
    tab->setMargin( 0 );
    m_layout->addWidget( tab, 0, 0 );
    m_layout->setColumnStretch( 0, 1 );
    m_layout->setRowStretch( 0, 1 );
    
    tab->addTab(peaksDiv, WString::tr("ssd-phone-tab-source-peaks"), Wt::WTabWidget::PreLoading);
    tab->addTab(sourceDiv, WString::tr("ssd-phone-tab-source-isotopes"), Wt::WTabWidget::PreLoading);
    tab->addTab(detectorDiv, WString::tr("ssd-phone-tab-shielding"), Wt::WTabWidget::PreLoading);
    
    WContainerWidget *chartDiv = new WContainerWidget();
    chartDiv->setOffsets(0);
    chartDiv->setMargin(0);
    chartDiv->setPadding(5);
    WGridLayout *chartLayout = new WGridLayout();
    chartDiv->setLayout(chartLayout);
    chartLayout->setContentsMargins(0, 0, 0, 0);
    
    chartLayout->addWidget( m_detectorDisplay,      0, 0, AlignLeft );
    chartLayout->addWidget( addItemMenubutton,      0, 1, AlignRight);
    chartLayout->addWidget( m_chi2Graphic,          1, 0, 1, 2 );
    m_showChiOnChart->setWidth( 130 );
    chartLayout->addWidget( m_showChiOnChart,       2, 1, AlignRight );
    chartLayout->addWidget( m_optionsDiv,           3, 0, 1, 2 );
    chartLayout->addWidget( m_fitModelButton,       4, 0, 1, 2, AlignCenter );
    chartLayout->addWidget( m_fitProgressTxt,       5, 0, 1, 2, AlignCenter );
    chartLayout->addWidget( m_cancelfitModelButton, 6, 0, 1, 2, AlignCenter );
    
    chartLayout->setRowStretch(1, 1);
    tab->addTab(chartDiv,"Fit", Wt::WTabWidget::PreLoading);
  }else
  {
    //regular layout
    
    // We'll put the detector and menu icon in a flexbox layout, which is less hassle than
    //  the Wt layout
    WContainerWidget *toprow = new WContainerWidget();
    toprow->addStyleClass( "DetAndMenu" );
    toprow->addWidget( m_detectorDisplay );
    toprow->addWidget( addItemMenubutton );
    
    
    detectorLayout->addWidget( toprow,                     0, 0 );
    detectorLayout->addWidget( smallerContainer,           1, 0 );
    detectorLayout->addWidget( m_shieldingSelects,         2, 0 );
    detectorLayout->addWidget( m_fitModelButton,           3, 0, AlignCenter );
    detectorLayout->addWidget( m_fitProgressTxt,           4, 0 );
    detectorLayout->addWidget( m_cancelfitModelButton,     5, 0, AlignCenter );
    detectorLayout->addWidget( m_showChi2Text,             6, 0 );
    
    detectorLayout->setRowStretch( 2, 1 );
    detectorLayout->setHorizontalSpacing( 0 );
    detectorLayout->setVerticalSpacing( 0 );
    detectorLayout->setContentsMargins( 1, 1, 1, 1 );
    
    WGridLayout *tablesLayout = new WGridLayout();
    
    tablesLayout->addWidget( peaksDiv,   0, 0 );
    tablesLayout->addWidget( sourceDiv, 0, 1 );
    tablesLayout->setColumnResizable( 0, true, WLength(340,WLength::Pixel) ); //335px seems to be the limit where the peak table will get horizontal scroll-bars
    tablesLayout->setHorizontalSpacing( 5 );
    tablesLayout->setVerticalSpacing( 0 );
    tablesLayout->setContentsMargins( 0, 0, 0, 0 );
    
    
    // We will put the chart in a div that will also hold the Rel/Chi switch; its a bit of a hack
    //  to get the chart type switch near the chart; it would probably be best to have the switch be
    //  at the top of the chart, but that doesnt work so well because the chart <img> will be over
    //  the switch if we want the switch to be over the image...  probably something better to do
    //  here
    WContainerWidget *chartHolder = new WContainerWidget();
    WGridLayout *chartLayout = new WGridLayout( chartHolder );
    chartLayout->setVerticalSpacing( 0 );
    chartLayout->setHorizontalSpacing( 0 );
    chartLayout->setContentsMargins( 0, 0, 0, 0 );
    chartLayout->addWidget( m_chi2Graphic, 0, 0 );
    
    //Put the switch in a <div> (which will be 0x0 px), so we can position the switch using absolute
    WContainerWidget *switchHolder = new WContainerWidget();
    switchHolder->addWidget( m_showChiOnChart );
    m_showChiOnChart->setAttributeValue( "style", "position: absolute; bottom: 0px; right: 20px" );
    switchHolder->setHeight( 0 );
    chartLayout->addWidget( switchHolder, 1, 0, AlignRight );
    chartLayout->setRowStretch( 0, 1 );
    
    WGridLayout *leftLayout = new WGridLayout();
    leftLayout->addWidget( chartHolder,    0, 0 );
    leftLayout->addLayout( tablesLayout,   1, 0 );
    leftLayout->addWidget( m_optionsDiv,   2, 0 );
    leftLayout->setRowResizable( 0, true, WLength(40.0,WLength::Percentage) );
    leftLayout->setRowStretch( 1, 1 );
    leftLayout->setHorizontalSpacing( 5 );
    leftLayout->setVerticalSpacing( 5 );
    leftLayout->setContentsMargins( 0, 0, 0, 0 );
    
    WContainerWidget *leftDiv = new WContainerWidget();
    leftDiv->setLayout(leftLayout);
    leftDiv->setOverflow(WContainerWidget::OverflowHidden);
    
    detectorDiv->setOverflow( WContainerWidget::OverflowHidden );
    detectorDiv->setWidth( 290 );
    
    m_layout->addWidget( leftDiv, 0, 0);
    m_layout->addWidget( detectorDiv, 0, 1 );
    m_layout->setColumnStretch( 0, 1 );
    m_layout->setHorizontalSpacing( 0 );
    m_layout->setVerticalSpacing( 0 );
    
    setOverflow( WContainerWidget::OverflowVisible );
    setOffsets( WLength(0,WLength::Pixel) );
  } //regular layout
  
  handleDetectorChanged( m_detectorDisplay->detector() ); // Will also call updateChi2Chart()
}//ShieldingSourceDisplay constructor


//When the button is triggered, update model
void ShieldingSourceDisplay::toggleUseAll( Wt::WCheckBox *button )
{
  UndoRedoManager::PeakModelChange peak_undo_creator;
  
  const bool useForFit = button->isChecked();
  const size_t npeaks = m_peakModel->npeaks();
  
  for( size_t i = 0; i < npeaks; ++i )
  {
    try
    {
      const PeakModel::PeakShrdPtr peak = m_peakModel->peakPtr( i );
      WModelIndex index = m_peakModel->indexOfPeak( peak );
      
      if( !!peak && index.isValid() && peak->parentNuclide() )
      {
        index = m_peakModel->index( index.row(), PeakModel::kUseForShieldingSourceFit );
        m_peakModel->setData( index, useForFit );
      }
    }catch( std::exception & )
    {
      //shouldnt ever happen, wont worry about
    }
  }//for( size_t i = 0; i < npeaks; ++i )
}//void ShieldingSourceDisplay::toggleUseAll(Wt::WCheckBox* button)


//When model is changed, update checkbox tristate
void ShieldingSourceDisplay::updateAllPeaksCheckBox( WCheckBox *but)
{
    bool allon=true;
    bool alloff=true;
    
    for( int peakn = 0; peakn < m_peakModel->rowCount(); ++peakn )
    {
        WModelIndex index = m_peakModel->index( peakn,
                                               PeakModel::kUseForShieldingSourceFit );
        const PeakModel::PeakShrdPtr &peak = m_peakModel->peak( index);
        if( peak->useForShieldingSourceFit() && alloff)
        {
            alloff=false;
        } //if( peak->useForShieldingSourceFit() && alloff)
        else if( !peak->useForShieldingSourceFit() && allon)
        {
            allon=false;
        } //if( !peak->useForShieldingSourceFit() && allon)
    }//for
    
    if (alloff && !allon)
        but->setCheckState(Wt::Unchecked);
    else if (allon && !alloff)
        but->setCheckState(Wt::Checked);
    else
        but->setCheckState(Wt::PartiallyChecked);
} //updateAllPeaksCheckBox( WCheckBox *but)


void ShieldingSourceDisplay::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool renderFull = (flags & Wt::RenderFlag::RenderFull);
  
  if( m_chi2ChartNeedsUpdating )
  {
    updateChi2ChartActual( nullptr );
    m_chi2ChartNeedsUpdating = false;
  }
  
  WContainerWidget::render( flags );
}//void render( Wt::WFlags<Wt::RenderFlag> flags )


SourceFitModel *ShieldingSourceDisplay::sourceFitModel()
{
  return m_sourceModel;
}


ShieldingSourceFitCalc::ShieldingSourceFitOptions ShieldingSourceDisplay::fitOptions() const
{
  const shared_ptr<const DetectorPeakResponse> det = m_detectorDisplay->detector();
  const bool fixed_geom = (det && det->isFixedGeometry());
  assert( (m_photopeak_cluster_sigma - m_photopeak_cluster_sigma) < 1.0E-6 );
  
  ShieldingSourceFitCalc::ShieldingSourceFitOptions options;
  options.multiple_nucs_contribute_to_peaks = m_multiIsoPerPeak->isChecked();
  options.attenuate_for_air = (!fixed_geom && m_attenForAir->isChecked());
  options.account_for_decay_during_meas = m_decayCorrect->isChecked();
  options.multithread_self_atten = m_multithread_computation;
  options.photopeak_cluster_sigma = m_photopeak_cluster_sigma;
  options.background_peak_subtract = m_backgroundPeakSub->isChecked();
  options.same_age_isotopes = m_sameIsotopesAge->isChecked();
  
  return options;
}//ShieldingSourceFitOptions fitOptions() const



ShieldingSourceDisplay::~ShieldingSourceDisplay() noexcept(true)
{
  {//begin make sure calculation is cancelled
    std::lock_guard<std::mutex> lock( m_currentFitFcnMutex );
    if( m_currentFitFcn )
    {
      try
      {
        m_currentFitFcn->cancelFit();
      }catch( ... )
      {
        cerr << "Caught exception call m_currentFitFcn->cancelFit(), which probably shouldnt happen" << endl;
      }
    }
  }//end make sure calculation is cancelled
  
  if( m_addItemMenu )
  {
    delete m_addItemMenu;
    m_addItemMenu = NULL;
  }//if( m_addItemMenu )
  
  closeModelUploadWindow();
#if( USE_DB_TO_STORE_SPECTRA )
  closeBrowseDatabaseModelsWindow();
  closeSaveModelToDatabaseWindow();
#endif
}//ShieldingSourceDisplay destructor constructor


pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters>
                                                      ShieldingSourceDisplay::shieldingFitnessFcn()
{
  //make sure fitting for at least one nuclide:
  std::vector<ShieldingSourceFitCalc::ShieldingInfo> initial_shieldings;
  
  for( WWidget *widget : m_shieldingSelects->children() )
  {
    ShieldingSelect *select = dynamic_cast<ShieldingSelect *>( widget );
    if( select )
    {
      ShieldingSourceFitCalc::ShieldingInfo info = select->toShieldingInfo();
      
      if( !info.m_isGenericMaterial && !info.m_material )
        info.m_material = make_shared<Material>( *m_materialDB->material( "void" ) );
      
      initial_shieldings.push_back( info );
    }//if( select )
  }//for( WWidget *widget : m_shieldingSelects->children() )
    
  const GeometryType geom = geometry();
  const string distanceStr = m_distanceEdit->text().toUTF8();
  const double distance = PhysicalUnits::stringToDistance( distanceStr );
  shared_ptr<SpecMeas> meas_file = m_specViewer->measurment(SpecUtils::SpectrumType::Foreground);
  shared_ptr<const SpecUtils::Measurement> foreground = m_specViewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  shared_ptr<const SpecUtils::Measurement> background = m_specViewer->displayedHistogram(SpecUtils::SpectrumType::Background);
  shared_ptr<const DetectorPeakResponse> detector = meas_file ? meas_file->detector() : nullptr;
  shared_ptr<const deque< PeakModel::PeakShrdPtr > > peaks = m_peakModel->peaks();
  if( !peaks )
    throw runtime_error( "No peaks." );
    
  shared_ptr<const deque<shared_ptr<const PeakDef>>> background_peaks;
    
  if( m_backgroundPeakSub->isChecked() )
  {
    shared_ptr<const SpecMeas> back = m_specViewer->measurment(SpecUtils::SpectrumType::Background);
      
    if( back && background )
    {
      const auto &displayed = m_specViewer->displayedSamples(SpecUtils::SpectrumType::Background);
      background_peaks = back->peaks( displayed );
    }//if( back )
      
    if( !background_peaks || background_peaks->empty() )
    {
      m_backgroundPeakSub->setUnChecked();
      passMessage( WString::tr("ssd-warn-no-back-peaks"), WarningWidget::WarningMsgInfo );
    }//if( !peaks || peaks->empty() )
  }//if( m_backgroundPeakSub->isChecked() )
    
  vector<ShieldingSourceFitCalc::SourceFitDef> src_definitions;
  // Could we instead just directly use m_sourceModel->underlyingData()?
  //  Maybe with trace sources getting slightly modified to be quantity being fit, and not total activity
  for( int ison = 0; ison < m_sourceModel->numNuclides(); ++ison )
  {
    const SandiaDecay::Nuclide * const nuclide = m_sourceModel->nuclide(ison);
    assert( nuclide );
    if( !nuclide )
      throw runtime_error( "Invalid source file nuclide" );
      
    ShieldingSourceFitCalc::SourceFitDef srcdef;
    srcdef.nuclide = nuclide;
    srcdef.age = m_sourceModel->age( ison );
    srcdef.fitAge = m_sourceModel->fitAge( ison );
    srcdef.activity = m_sourceModel->activity( ison );
    srcdef.ageDefiningNuc = m_sourceModel->ageDefiningNuclide( nuclide );
    srcdef.sourceType = m_sourceModel->sourceType(ison);
      
#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
    for( const ShieldingSourceFitCalc::IsoFitStruct &underdata : m_sourceModel->underlyingData() )
    {
      if( underdata.nuclide == nuclide )
      {
        srcdef.truthActivity = underdata.truthActivity;
        srcdef.truthActivityTolerance = underdata.truthActivityTolerance;
        srcdef.truthAge = underdata.truthAge;
        srcdef.truthAgeTolerance = underdata.truthAgeTolerance;
        break;
      }
    }//for( const ShieldingSourceFitCalc::IsoFitStruct &underdata : m_sourceModel->underlyingData() )
#endif
    
    
    switch( srcdef.sourceType )
    {
      case ShieldingSourceFitCalc::ModelSourceType::Point:
        srcdef.fitActivity = m_sourceModel->fitActivity( ison );
        break;
          
      case ShieldingSourceFitCalc::ModelSourceType::Intrinsic:
        srcdef.fitActivity = false;
        break;
        
      case ShieldingSourceFitCalc::ModelSourceType::Trace:
      {
        // Go through shieldings and get display activity, so we can fit for that.
        int numShieldingsTraceSrcFor = 0;
        for( const ShieldingSourceFitCalc::ShieldingInfo &shield : initial_shieldings )
        {
          const ShieldingSourceFitCalc::TraceSourceInfo *trace_info = nullptr;
          for( size_t i = 0; !trace_info && (i < shield.m_traceSources.size()); ++i )
            trace_info = (shield.m_traceSources[i].m_nuclide == srcdef.nuclide) ? &(shield.m_traceSources[i]) : nullptr;
          
          if( trace_info )
          {
            numShieldingsTraceSrcFor += 1;
            srcdef.activity = trace_info->m_activity;
            srcdef.fitActivity = trace_info->m_fitActivity;
            
            // Even though it doesnt really matter, lets try to keep the model in sync with trace
            //  widget, so we'll toss in a development check for it
            if( srcdef.fitActivity != m_sourceModel->fitActivity(static_cast<int>(ison)) )
            {
              cerr << "\n\n\n\nTemporarily disabling assert 'fitAct=" << srcdef.fitActivity << "'- reaenable\n\n\n" << endl;
  //           assert( fitAct == m_sourceModel->fitActivity(nucn) );
            }
          }//if( this shielding has the nuclide as a trace source )
        }//for( WWidget *w : m_shieldingSelects->children() )
          
        assert( numShieldingsTraceSrcFor == 1 );
          
        if( numShieldingsTraceSrcFor != 1 )
          throw runtime_error( "Unexpected inconsistent state - couldnt find trace source widget for " + nuclide->symbol );
        break;
      }//case ShieldingSourceFitCalc::ModelSourceType::Trace:
    }//switch( m_sourceModel->sourceType(ison) )
      
    // If we are fitting activity or age, we can possibly get into a situation where the
    //  values have become NaN - if this is the case, lets put in a number, to hopefully
    //  help get out of this badness
    if( srcdef.fitAge && (IsInf(srcdef.age) || IsNan(srcdef.age)) )
      srcdef.age = PeakDef::defaultDecayTime( nuclide, nullptr );
    
    if( srcdef.fitActivity && (IsInf(srcdef.activity) || IsNan(srcdef.activity)) )
      srcdef.activity = 1.0E-6 * PhysicalUnits::curie;
    
    src_definitions.push_back( srcdef );
  }//for( const SandiaDecay::Nuclide *nuc : nuclides )
    
  const ShieldingSourceFitCalc::ShieldingSourceFitOptions options = fitOptions();
    
  return GammaInteractionCalc::ShieldingSourceChi2Fcn::create( distance, geom,
                                  initial_shieldings, src_definitions,
                                  detector, foreground, background, *peaks,
                                  background_peaks, options );
}//pair<ShieldingSourceChi2Fcn,ROOT::Minuit2::MnUserParameters> shieldingFitnessFcn()

  
#if( INCLUDE_ANALYSIS_TEST_SUITE )
void ShieldingSourceDisplay::showInputTruthValuesWindow()
{
  //The error handling and display quality is minimal since this is a developer-only tool, and will
  //  not be used by end-users.
  //Also, if you change the model any while this window is open - bad things will happen.
  
  AuxWindow *window = new AuxWindow( "Input Truth Values",
                                     (AuxWindowProperties::IsModal | AuxWindowProperties::TabletNotFullScreen) );
  
  WContainerWidget *contents = window->contents();
  
  try
  {
    WTable *table = new WTable( contents );
    table->addStyleClass( "TruthValueTable" );
    table->setHeaderCount( 1 );
    new WLabel( "Quantity", table->elementAt(0, 0) );
    new WLabel( "Value", table->elementAt(0, 1) );
    new WLabel( "Tolerance", table->elementAt(0, 2) );
    
    const int nnuc = m_sourceModel->numNuclides();
    for( int i = 0; i < nnuc; ++i )
    {
      const SandiaDecay::Nuclide *nuc = m_sourceModel->nuclide( i );
      assert( nuc );
      
      const SandiaDecay::Nuclide *ageNuc = m_sourceModel->ageDefiningNuclide( nuc );
      const ShieldingSourceFitCalc::ModelSourceType sourceType = m_sourceModel->sourceType( i );
      const bool selfAttNuc = (sourceType == ShieldingSourceFitCalc::ModelSourceType::Intrinsic);
      
//      if( selfAttNuc )
//        throw runtime_error( "Model is not candidate for truth-level info<br />"
//                             "Self-attuating sources not implemented yet" );
//      if( ageNuc && (ageNuc != nuc) )
//        throw runtime_error( "Model is not candidate for truth-level info<br />"
//                             "Shared-age nuclides not allowed" );
      
      // For self-attenuating shieldings, we'll just test the shielding thickness
      // For nuclides whose age is controlled by another nuclide, we dont need to test age.
      if( selfAttNuc || (ageNuc && (ageNuc != nuc)) )
        continue;
      
      const bool fitAct = m_sourceModel->fitActivity(i);
      const bool fitAge = m_sourceModel->fitAge(i);
      
      auto setFieldValue = [nuc,this]( WLineEdit *valuefld, const SourceFitModel::Columns type ){
        const int srcrow = m_sourceModel->nuclideIndex(nuc);
        auto index = m_sourceModel->index( srcrow, type );
        valuefld->setText( asString( m_sourceModel->data(index) ) );
      };//setFieldValue(...)
      
      auto fieldUpdate = [this,nuc]( WLineEdit *valuefld, const SourceFitModel::Columns type ){
        const string valtxt = valuefld->text().toUTF8();
        
        const int srcrow = m_sourceModel->nuclideIndex(nuc);
        WModelIndex index = m_sourceModel->index( srcrow, type );
        
        try
        {
          m_sourceModel->setData(index, valtxt);
        }catch(...)
        {
          passMessage( "'" + valtxt + "' is not a valid entry", WarningWidget::WarningMsgHigh );
          valuefld->setText( asString( m_sourceModel->data(index) ) );
        }
      };//fieldUpdate
      
      if( fitAct )
      {
        const int row = table->rowCount();
        WLabel *label = new WLabel( nuc->symbol + " Activity", table->elementAt(row, 0) );
        WLineEdit *value = new WLineEdit( table->elementAt(row, 1) );
        label->setBuddy( value );
        
        value->setAutoComplete( false );
        value->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
        value->setAttributeValue( "autocorrect", "off" );
        value->setAttributeValue( "spellcheck", "off" );
#endif
        setFieldValue( value, SourceFitModel::Columns::kTruthActivity );
        
        auto valueUpdate = [fieldUpdate,value](){
          fieldUpdate( value, SourceFitModel::Columns::kTruthActivity );
        };
        
        value->changed().connect( std::bind(valueUpdate) );
        value->enterPressed().connect( std::bind(valueUpdate) );
        
        
        WLineEdit *tolerance = new WLineEdit( table->elementAt(row, 2) );
        
        tolerance->setAutoComplete( false );
        tolerance->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
        tolerance->setAttributeValue( "autocorrect", "off" );
        tolerance->setAttributeValue( "spellcheck", "off" );
#endif
        setFieldValue( tolerance, SourceFitModel::Columns::kTruthActivityTolerance );
        auto toleranceUpdate = [fieldUpdate,tolerance](){
          fieldUpdate( tolerance, SourceFitModel::Columns::kTruthActivityTolerance );
        };
        
        tolerance->changed().connect( std::bind(toleranceUpdate) );
        tolerance->enterPressed().connect( std::bind(toleranceUpdate) );
      }//if( fitAct )
      
      
      if( fitAge )
      {
        const int row = table->rowCount();
        new WLabel( nuc->symbol + " Age", table->elementAt(row, 0) );
        WLineEdit *value = new WLineEdit( table->elementAt(row, 1) );
        
        value->setAutoComplete( false );
        value->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
        value->setAttributeValue( "autocorrect", "off" );
        value->setAttributeValue( "spellcheck", "off" );
#endif
        setFieldValue( value, SourceFitModel::Columns::kTruthAge );
        
        auto valueUpdate = [fieldUpdate,value](){
          fieldUpdate( value, SourceFitModel::Columns::kTruthAge );
        };
        
        value->changed().connect( std::bind(valueUpdate) );
        value->enterPressed().connect( std::bind(valueUpdate) );
        
        
        WLineEdit *tolerance = new WLineEdit( table->elementAt(row, 2) );
        
        tolerance->setAutoComplete( false );
        tolerance->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
        tolerance->setAttributeValue( "autocorrect", "off" );
        tolerance->setAttributeValue( "spellcheck", "off" );
#endif
        setFieldValue( tolerance, SourceFitModel::Columns::kTruthAgeTolerance );
        auto toleranceUpdate = [fieldUpdate,tolerance](){
          fieldUpdate( tolerance, SourceFitModel::Columns::kTruthAgeTolerance );
        };
        
        tolerance->changed().connect( std::bind(toleranceUpdate) );
        tolerance->enterPressed().connect( std::bind(toleranceUpdate) );
      }//if( fitAge )
    }//for( int i = 0; i < nnuc; ++i )
    
    
    for( WWidget *widget : m_shieldingSelects->children() )
    {
      ShieldingSelect *select = dynamic_cast<ShieldingSelect *>(widget);
      if( !select )
        continue;
      
      if( select->isGenericMaterial() )
      {
        if( select->fitArealDensity() )
        {
          const int row = table->rowCount();
          WLabel *label = new WLabel( "Areal Density", table->elementAt(row, 0) );
          WLineEdit *value = new WLineEdit( table->elementAt(row, 1) );
          
          value->setAutoComplete( false );
          value->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
          value->setAttributeValue( "autocorrect", "off" );
          value->setAttributeValue( "spellcheck", "off" );
#endif
          WDoubleValidator *dblValidator = new WDoubleValidator( 0, 500, value );
          value->setValidator( dblValidator );
          value->addStyleClass( "numberValidator"); //used to detect mobile keyboard
          label->setBuddy( value );
          
          if( select->truthAD )
            value->setText( std::to_string(*select->truthAD) );
          
          auto updateVal = [select,value](){
            double answer = 0;
            if( (stringstream(value->text().toUTF8()) >> answer) )
              select->truthAD = answer;
            else if( select->truthAD )
              value->setText( std::to_string(*select->truthAD) );
            else
              value->setText( "" );
          };//updateVal(...)
          
          value->changed().connect( std::bind(updateVal) );
          value->enterPressed().connect( std::bind(updateVal) );
          
          WLineEdit *tolerance = new WLineEdit( table->elementAt(row, 2) );
          
          tolerance->setAutoComplete( false );
          tolerance->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
          tolerance->setAttributeValue( "autocorrect", "off" );
          tolerance->setAttributeValue( "spellcheck", "off" );
#endif
          dblValidator = new WDoubleValidator( 0, 100, tolerance );
          tolerance->setValidator( dblValidator );
          tolerance->addStyleClass( "numberValidator"); //used to detect mobile keyboard
          
          auto updateTolerance = [select,tolerance](){
            double answer = 0;
            if( (stringstream(tolerance->text().toUTF8()) >> answer) )
              select->truthADTolerance = answer;
            else if( select->truthADTolerance )
              tolerance->setText( std::to_string(*select->truthADTolerance) );
            else
              tolerance->setText( "" );
          };//updateVal(...)
          
          if( select->truthADTolerance )
            tolerance->setText( std::to_string(*select->truthADTolerance) );
          
          tolerance->changed().connect( std::bind(updateTolerance) );
          tolerance->enterPressed().connect( std::bind(updateTolerance) );
        }//if( fit AD )
        
        if( select->fitAtomicNumber() )
        {
          const int row = table->rowCount();
          WLabel *label = new WLabel( "Atomic Number", table->elementAt(row, 0) );
          WLineEdit *value = new WLineEdit( table->elementAt(row, 1) );
          
          value->setAutoComplete( false );
          value->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
          value->setAttributeValue( "autocorrect", "off" );
          value->setAttributeValue( "spellcheck", "off" );
#endif
          WDoubleValidator *dblValidator = new WDoubleValidator( MassAttenuation::sm_min_xs_atomic_number,
                                                                 MassAttenuation::sm_max_xs_atomic_number, value );
          value->setValidator( dblValidator );
          value->addStyleClass( "numberValidator"); //used to detect mobile keyboard
          label->setBuddy( value );
          
          if( select->truthAN )
            value->setText( std::to_string(*select->truthAN) );
          
          auto updateVal = [select,value](){
            double answer = 0;
            if( (stringstream(value->text().toUTF8()) >> answer) )
              select->truthAN = answer;
            else if( select->truthAN )
              value->setText( std::to_string(*select->truthAN) );
            else
              value->setText( "" );
          };//updateVal(...)
          
          value->changed().connect( std::bind(updateVal) );
          value->enterPressed().connect( std::bind(updateVal) );
          
          WLineEdit *tolerance = new WLineEdit( table->elementAt(row, 2) );
          
          tolerance->setAutoComplete( false );
          tolerance->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
          tolerance->setAttributeValue( "autocorrect", "off" );
          tolerance->setAttributeValue( "spellcheck", "off" );
#endif
          dblValidator = new WDoubleValidator( 0, 100, tolerance );
          tolerance->setValidator( dblValidator );
          tolerance->addStyleClass( "numberValidator"); //used to detect mobile keyboard
          
          auto updateTolerance = [select,tolerance](){
            double answer = 0;
            const string txt = tolerance->text().toUTF8();
            if( txt.empty() )
            {
              select->truthANTolerance.reset();
            }else if( (stringstream(txt) >> answer) )
            {
              select->truthANTolerance = answer;
            }else if( select->truthANTolerance )
            {
              tolerance->setText( std::to_string(*select->truthANTolerance) );
            }else
            {
              select->truthANTolerance.reset();
              tolerance->setText( "" );
            }
          };//updateVal(...)
          
          if( select->truthANTolerance )
            tolerance->setText( std::to_string(*select->truthANTolerance) );
          
          tolerance->changed().connect( std::bind(updateTolerance) );
          tolerance->enterPressed().connect( std::bind(updateTolerance) );
        }//if( fit AN )
      }else
      {
        shared_ptr<const Material> mat = select->material();
        if( !mat )
          throw runtime_error( "There is a non-generic material that is blank" );
        
//        vector<const SandiaDecay::Nuclide *> srcnucs = select->selfAttenNuclides();
//        if( !srcnucs.empty() )
//          throw runtime_error( "A shieldings used as sources is not yet implemented" );
        
        const auto geom = geometry();
        
        auto setupLength = [this, select, table, geom]( const int dim ){
          assert( (dim >= 0) && (dim <= 2) );
          
          const char *labeltxt = "Thickness";
          
          switch( geom )
          {
            case GammaInteractionCalc::GeometryType::Spherical:
              assert( dim == 0 );
              labeltxt = "Thickness";
              if( (dim != 0) || !select->fitThickness() )
                return;
              break;
              
            case GammaInteractionCalc::GeometryType::CylinderEndOn:
            case GammaInteractionCalc::GeometryType::CylinderSideOn:
              assert( (dim == 0) || (dim == 1) );
              if( (dim != 0) && (dim != 1) )
                return;
              
              if( (dim == 0) && !select->fitCylindricalRadiusThickness() )
                return;
              
              if( (dim == 1) && !select->fitCylindricalLengthThickness() )
                return;
              
              labeltxt = (dim == 0) ? "Rad. Thick." : "Cyl. Len.";
              break;
              
            case GammaInteractionCalc::GeometryType::Rectangular:
              assert( (dim == 0) || (dim == 1) || (dim == 2) );
              if( (dim != 0) && (dim != 1) && (dim != 2) )
                return;
              
              if( (dim == 0) && !select->fitRectangularWidthThickness() )
                return;
              
              if( (dim == 1) && !select->fitRectangularHeightThickness() )
                return;
              
              if( (dim == 2) && !select->fitRectangularDepthThickness() )
                return;
              
              labeltxt = (dim == 0) ? "Width" : ((dim == 1) ? "Height" : "Depth");
              break;
              
            case GammaInteractionCalc::GeometryType::NumGeometryType:
              assert( 0 );
              return;
              break;
          }//switch( geom )
      
          
          boost::optional<double> *thicknessVal = nullptr;
          boost::optional<double> *toleranceVal = nullptr;
          
          
          switch( dim )
          {
            case 0:
              thicknessVal = &(select->truthThickness);
              toleranceVal = &(select->truthThicknessTolerance);
              break;
              
            case 1:
              thicknessVal = &(select->truthThicknessD2);
              toleranceVal = &(select->truthThicknessD2Tolerance);
              break;
              
            case 2:
              thicknessVal = &(select->truthThicknessD3);
              toleranceVal = &(select->truthThicknessD3Tolerance);
              break;
              
            default:
              throw std::logic_error( "Invalid truth dim" );
          }//switch( dim )
          
          assert( thicknessVal && toleranceVal );
          
          const int row = table->rowCount();
          WLabel *label = new WLabel( labeltxt, table->elementAt(row, 0) );
          WLineEdit *value = new WLineEdit( table->elementAt(row, 1) );
          
          value->setAutoComplete( false );
          value->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
          value->setAttributeValue( "autocorrect", "off" );
          value->setAttributeValue( "spellcheck", "off" );
#endif
          WRegExpValidator *validator = new WRegExpValidator( PhysicalUnits::sm_distanceUncertaintyUnitsOptionalRegex, value );
          validator->setFlags( Wt::MatchCaseInsensitive );
          value->setValidator( validator );
          label->setBuddy( value );
          
          if( *thicknessVal )
            value->setText( PhysicalUnits::printToBestLengthUnits( **thicknessVal, 4 ) );
          
          auto updateVal = [value,thicknessVal](){
            const string txt = value->text().toUTF8();
            if( txt.empty() )
            {
              thicknessVal->reset();
              return;
            }
            
            try
            {
              (*thicknessVal) = PhysicalUnits::stringToDistance( txt );
            }catch( ... )
            {
              if( *thicknessVal )
                value->setText( PhysicalUnits::printToBestLengthUnits( **thicknessVal, 4) );
              else
                value->setText( "" );
            }//try / catch
          };//updateVal(...)
          
          value->changed().connect( std::bind(updateVal) );
          value->enterPressed().connect( std::bind(updateVal) );
          
          WLineEdit *tolerance = new WLineEdit( table->elementAt(row, 2) );
          
          tolerance->setAutoComplete( false );
          tolerance->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
          tolerance->setAttributeValue( "autocorrect", "off" );
          tolerance->setAttributeValue( "spellcheck", "off" );
#endif
          
          validator = new WRegExpValidator( PhysicalUnits::sm_distanceUncertaintyUnitsOptionalRegex, tolerance );
          validator->setFlags( Wt::MatchCaseInsensitive );
          tolerance->setValidator( validator );
          
          auto updateTolerance = [tolerance,toleranceVal](){
            const string txt = tolerance->text().toUTF8();
            if( txt.empty() )
            {
              toleranceVal->reset();
              return;
            }
            
            try
            {
              (*toleranceVal) = PhysicalUnits::stringToDistance( txt );
            }catch( ... )
            {
              if( *toleranceVal )
                tolerance->setText( PhysicalUnits::printToBestLengthUnits( **toleranceVal,4) );
              else
                tolerance->setText( "" );
            }//try / catch
          };//updateVal(...)
          
          if( *toleranceVal )
            tolerance->setText( PhysicalUnits::printToBestLengthUnits( **toleranceVal ) );
          
          tolerance->changed().connect( std::bind(updateTolerance) );
          tolerance->enterPressed().connect( std::bind(updateTolerance) );
        };//setupLength lamda
        
        switch( geom )
        {
          case GammaInteractionCalc::GeometryType::Spherical:
            setupLength(0);
            break;
            
          case GammaInteractionCalc::GeometryType::CylinderEndOn:
          case GammaInteractionCalc::GeometryType::CylinderSideOn:
            setupLength(0);
            setupLength(1);
            break;
            
          case GammaInteractionCalc::GeometryType::Rectangular:
            setupLength(0);
            setupLength(1);
            setupLength(2);
            break;
            
          case GammaInteractionCalc::GeometryType::NumGeometryType:
            assert( 0 );
            break;
        }//switch( geometry() )
        
        
        if( select->fitForAnyMassFractions() )
        {
          const map<const SandiaDecay::Element *,vector<ShieldingSelect::NucMasFrac>> el_currentMassFractions
                                                        = select->sourceNuclideMassFractions();
          
          // Get rid of any truth mass-fractions, that are no longer self-attenuating sources
          set<const SandiaDecay::Nuclide *> nonexistent_nucs;
          for( const auto &el_nucfracs : el_currentMassFractions )
          {
            const vector<ShieldingSelect::NucMasFrac> &currentMassFractions = el_nucfracs.second;
            
            for( const auto &prev_el_to_nuc_fracs : select->truthFitMassFractions )
            {
              const map<const SandiaDecay::Nuclide *,std::pair<double,double>> &prevTruthMassFrac 
                                                                  = prev_el_to_nuc_fracs.second;
              
              for( const auto &prev_nuc_frac : prevTruthMassFrac )
              {
                bool is_current = false;
                for( size_t i = 0; !is_current && (i < currentMassFractions.size()); ++i )
                  is_current = (get<0>(currentMassFractions[i]) == prev_nuc_frac.first);
                if( !is_current
                   || !prev_nuc_frac.first
                   || (prev_nuc_frac.second.first < 0.0)
                   || (prev_nuc_frac.second.first > 1.0)
                   || (prev_nuc_frac.second.second < 0.0)
                   || (prev_nuc_frac.second.second > 1.0) )
                {
                  nonexistent_nucs.insert( prev_nuc_frac.first );
                }
              }//for( const auto &prev_nuc_frac : select->truthFitMassFractions )
            }//for( const auto &prev_el_to_nuc_fracs : select->truthFitMassFractions )
          }//for( const auto &el_nucfracs : el_currentMassFractions )
          
          for( const auto nuc : nonexistent_nucs )
          {
            for( auto &el_vals : select->truthFitMassFractions )
            {
              el_vals.second.erase( nuc ); //a little wasteful to cal for every element, but whatever
            }
          }
          
          for( const auto &el_nucfracs : el_currentMassFractions )
          {
            const SandiaDecay::Element *el = el_nucfracs.first;
            const vector<ShieldingSelect::NucMasFrac> &currentMassFractions = el_nucfracs.second;
            
            for( const tuple<const SandiaDecay::Nuclide *,double,bool> &mass_frac : currentMassFractions )
            {
              const SandiaDecay::Nuclide * const nuclide = get<0>(mass_frac);
              
              const int row = table->rowCount();
              WLabel *label = new WLabel( (nuclide ? nuclide->symbol : "Other") + " mass frac.", table->elementAt(row, 0) );
              
              NativeFloatSpinBox *value = new NativeFloatSpinBox( table->elementAt(row, 1) );
              value->setSpinnerHidden( true );
              value->setRange( 0.0, 1.0 );
              value->setText( "" );
              label->setBuddy( value );
              
              NativeFloatSpinBox *tolerance = new NativeFloatSpinBox( table->elementAt(row, 2) );
              tolerance->setSpinnerHidden( true );
              tolerance->setRange( 0.0, 1.0 );
              tolerance->setText( "" );
              
              map<const SandiaDecay::Nuclide *,pair<double,double>> &truthNucs = select->truthFitMassFractions[el];
              
              const auto pos = truthNucs.find(nuclide);
              if( pos != end(truthNucs) )
              {
                const double truthval = pos->second.first;
                const double truthtol = pos->second.second;
                value->setValue( truthval );
                tolerance->setValue( truthtol );
              }
              
              auto updateValAndTol = [select,value,tolerance,nuclide,&truthNucs](){
                if( value->text().empty() || tolerance->text().empty() )
                {
                  truthNucs.erase(nuclide);
                  return;
                }
                
                const double truthval = value->value();
                const double truthtol = tolerance->value();
                truthNucs[nuclide] = make_pair(truthval, truthtol);
              };//updateValAndTol(...)
              
              value->valueChanged().connect( std::bind(updateValAndTol) );
              tolerance->valueChanged().connect( std::bind(updateValAndTol) );
            }//for( const pair<const SandiaDecay::Nuclide *,double> &mass_frac : currentMassFractions )
          }//for( const auto &el_nucfracs : el_currentMassFractions )
        }//if( select->fitForAnyMassFractions() )
      }//if( generic material ) / else
    }//for( WWidget *widget : m_shieldingSelects->children() )
  }catch( std::exception &e )
  {
    contents->clear();
    WText *txt = new WText( e.what() , contents );
    txt->setInline( false );
  }//try / catch

  
  WPushButton *button = window->addCloseButtonToFooter("Okay");
  button->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
    
  window->centerWindow();
  window->disableCollapse();
  window->rejectWhenEscapePressed();
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  window->resizeToFitOnScreen();
  window->show();
}//showInputTruthValuesWindow()


void ShieldingSourceDisplay::setFitQuantitiesToDefaultValues()
{
  const int nnuc = m_sourceModel->numNuclides();
  for( int i = 0; i < nnuc; ++i )
  {
    const SandiaDecay::Nuclide *nuc = m_sourceModel->nuclide( i );
    assert( nuc );
    
    const SandiaDecay::Nuclide *ageNuc = m_sourceModel->ageDefiningNuclide( nuc );
    const ShieldingSourceFitCalc::ModelSourceType sourceType = m_sourceModel->sourceType( i );
    const bool selfAttNuc = (sourceType == ShieldingSourceFitCalc::ModelSourceType::Intrinsic);
    
    // For self-attenuating shieldings, we'll just test the shielding thickness
    // For nuclides whose age is controlled by another nuclide, we dont need to test age.
    if( selfAttNuc || (ageNuc && (ageNuc != nuc)) )
      continue;
    
    if( m_sourceModel->fitActivity(i) )
    {
      WModelIndex index = m_sourceModel->index( i, SourceFitModel::kActivity );
      const bool useCi = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", m_specViewer );
      if( useCi )
        m_sourceModel->setData( index, "1 mCi" );
      else
        m_sourceModel->setData( index, "37 MBq" );
    }//if( fit activity )
    
    if( m_sourceModel->fitAge(i) )
    {
      string agestr;
      PeakDef::defaultDecayTime( nuc, &agestr );
      WModelIndex index = m_sourceModel->index( i, SourceFitModel::kAge );
      m_sourceModel->setData( index, agestr );
    }//if( fit age )
  }//for( int i = 0; i < nnuc; ++i )
  
  
  for( WWidget *widget : m_shieldingSelects->children() )
  {
    ShieldingSelect *select = dynamic_cast<ShieldingSelect *>(widget);
    if( !select )
      continue;
    
    if( select->isGenericMaterial() )
    {
      if( select->fitArealDensity() )
        select->arealDensityEdit()->setText( "0" );
      if( select->fitAtomicNumber() )
        select->atomicNumberEdit()->setText( "26" );
    }else
    {
      shared_ptr<const Material> mat = select->material();
      if( !mat || select->fitForAnyMassFractions() )
        continue;
      
      switch( geometry() )
      {
        case GammaInteractionCalc::GeometryType::Spherical:
          if( select->fitThickness() )
            select->setSphericalThickness( 1.0*PhysicalUnits::cm );
          break;
          
        case GammaInteractionCalc::GeometryType::CylinderEndOn:
        case GammaInteractionCalc::GeometryType::CylinderSideOn:
          if( select->fitCylindricalRadiusThickness() )
            select->setCylindricalRadiusThickness( 0.5*PhysicalUnits::cm );
          if( select->fitCylindricalLengthThickness() )
            select->setCylindricalLengthThickness( 0.5*PhysicalUnits::cm );
          break;
          
        case GammaInteractionCalc::GeometryType::Rectangular:
          if( select->fitRectangularWidthThickness() )
            select->setRectangularWidthThickness( 0.5*PhysicalUnits::cm );
          if( select->fitRectangularHeightThickness() )
            select->setRectangularHeightThickness( 0.5*PhysicalUnits::cm );
          if( select->fitRectangularDepthThickness() )
            select->setRectangularDepthThickness( 0.5*PhysicalUnits::cm );
          break;
          
        case GammaInteractionCalc::GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }//switch( geometry() )
    }//if( generic material ) / else
  }//for( WWidget *widget : m_shieldingSelects->children() )
}//void setFitQuantitiesToDefaultValues()



std::tuple<int,int,bool> ShieldingSourceDisplay::numTruthValuesForFitValues()
{
  bool isValid = true;
  int nFitQuantities = 0, nQuantitiesCan = 0;
  
  const int nnuc = m_sourceModel->numNuclides();
  for( int i = 0; i < nnuc; ++i )
  {
    const SandiaDecay::Nuclide *nuc = m_sourceModel->nuclide( i );
    assert( nuc );
    
    const SandiaDecay::Nuclide *ageNuc = m_sourceModel->ageDefiningNuclide( nuc );
    const ShieldingSourceFitCalc::ModelSourceType sourceType = m_sourceModel->sourceType( i );
    const bool selfAttNuc = (sourceType == ShieldingSourceFitCalc::ModelSourceType::Intrinsic);
    
    // For self-attenuating shieldings, we'll just test the shielding thickness
    // For nuclides whose age is controlled by another nuclide, we dont need to test age.
    if( selfAttNuc || (ageNuc && (ageNuc != nuc)) )
      continue;
    
    if( m_sourceModel->fitActivity(i) )
    {
      const boost::optional<double> activity = m_sourceModel->truthActivity(i);
      const boost::optional<double> tolerance = m_sourceModel->truthActivityTolerance(i);
      nFitQuantities += 1;
      nQuantitiesCan += (activity && tolerance);
      if( !(activity && tolerance) )
      {
        auto actindex = m_sourceModel->index( i, SourceFitModel::kTruthActivity );
        auto tolindex = m_sourceModel->index( i, SourceFitModel::kTruthActivityTolerance );
        
        cerr << "Dont have: (activity && tolerance): (" << !!activity << " && " << !!tolerance << ") -> via data -> ("
        << asString(actindex.data()) << ", " << asString(tolindex.data()) << ")" << endl;
      }
    }//if( fit activity )
    
    if( m_sourceModel->fitAge(i) )
    {
      const boost::optional<double> age = m_sourceModel->truthAge(i);
      const boost::optional<double> tolerance = m_sourceModel->truthAgeTolerance(i);
      nFitQuantities += 1;
      nQuantitiesCan += (age && tolerance);
      if( !(age && tolerance) )
        cerr << "Dont have: (age && tolerance)" << endl;
    }//if( fit age )
  }//for( int i = 0; i < nnuc; ++i )
  
  
  for( WWidget *widget : m_shieldingSelects->children() )
  {
    ShieldingSelect *select = dynamic_cast<ShieldingSelect *>(widget);
    if( !select )
      continue;
    
    if( select->isGenericMaterial() )
    {
      if( select->fitArealDensity() )
      {
        if( select->truthAD && select->truthADTolerance )
          nQuantitiesCan += 1;
        nFitQuantities += 1;
      }//if( fit AD )
      
      if( select->fitAtomicNumber() )
      {
        if( select->truthAN && select->truthANTolerance )
          nQuantitiesCan += 1;
        nFitQuantities += 1;
      }//if( fit AN )
    }else
    {
      shared_ptr<const Material> mat = select->material();
      if( !mat )
      {
        cerr << "Dont have: Couldnt get material" << endl;
        isValid = false;
        continue;
      }
      
      switch( geometry() )
      {
        case GammaInteractionCalc::GeometryType::Spherical:
          if( select->fitThickness() )
          {
            if( select->truthThickness && select->truthThicknessTolerance )
              nQuantitiesCan += 1;
            nFitQuantities += 1;
          }//if( fit thickness )
          break;
          
        case GammaInteractionCalc::GeometryType::CylinderEndOn:
        case GammaInteractionCalc::GeometryType::CylinderSideOn:
          if( select->fitCylindricalRadiusThickness() )
          {
            if( select->truthThickness && select->truthThicknessTolerance )
              nQuantitiesCan += 1;
            nFitQuantities += 1;
          }//if( fit thickness )
          
          if( select->fitCylindricalLengthThickness() )
          {
            if( select->truthThicknessD2 && select->truthThicknessD2Tolerance )
              nQuantitiesCan += 1;
            nFitQuantities += 1;
          }//if( fit thickness )
          break;
          
          
        case GammaInteractionCalc::GeometryType::Rectangular:
          if( select->fitRectangularWidthThickness() )
          {
            if( select->truthThickness && select->truthThicknessTolerance )
              nQuantitiesCan += 1;
            nFitQuantities += 1;
          }
          
          if( select->fitRectangularHeightThickness() )
          {
            if( select->truthThicknessD2 && select->truthThicknessD2Tolerance )
              nQuantitiesCan += 1;
            nFitQuantities += 1;
          }
          
          if( select->fitRectangularDepthThickness() )
          {
            if( select->truthThicknessD3 && select->truthThicknessD3Tolerance )
              nQuantitiesCan += 1;
            nFitQuantities += 1;
          }
          break;
          
        case GammaInteractionCalc::GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }//switch( geometry() )
      
      
      map<const SandiaDecay::Element *,vector<ShieldingSelect::NucMasFrac>> currentElAndMassFractions
                                                            = select->sourceNuclideMassFractions();
        
      for( const auto &el_nucs : currentElAndMassFractions )
      {
        const SandiaDecay::Element * const el = el_nucs.first;
        const vector<ShieldingSelect::NucMasFrac> &currentMassFractions = el_nucs.second;
        for( const auto &current_nuc_frac : currentMassFractions )
        {
          nFitQuantities += 1;
          
          const auto is0To1 = []( const double v ) -> bool { return (v >= 0.0) && (v <= 1.0); };
          if( select->truthFitMassFractions.count(el) )
          {
            const auto &nucs = select->truthFitMassFractions[el];
            auto pos = nucs.find( get<0>(current_nuc_frac) );
            nQuantitiesCan += ((pos != end(nucs))
                               && is0To1(pos->second.first) && is0To1(pos->second.second) );
          }
        }//for( const auto &prev_nuc_frac : select->truthFitMassFractions )
      }//for( const auto &el_nucs : currentMassFractions )
    }//if( generic material ) / else
  }//for( WWidget *widget : m_shieldingSelects->children() )
  
  if( nQuantitiesCan != nFitQuantities )
  {
    cerr << "Dont have: nQuantitiesCan != nFitQuantities (" << nQuantitiesCan << " != " << nFitQuantities << ")" << endl;
    isValid = false;
  }
  
  if( !nQuantitiesCan )
  {
    cerr << "Dont have: !nQuantitiesCan: " << nQuantitiesCan << endl;
    isValid = false;
  }
  
  return std::tuple<int,int,bool>( nQuantitiesCan, nFitQuantities, isValid );
}//bool haveTruthValuesForAllFitValues()


void ShieldingSourceDisplay::renderChi2Chart( Wt::WSvgImage &image )
{
  WPainter p( &image );
  m_chi2Graphic->paint( p );
  p.end();
}//void renderChi2Chart( Wt::WSvgImage &image );


tuple<bool,int,int,vector<string>> ShieldingSourceDisplay::testCurrentFitAgainstTruth()
{
  bool successful = true;
  int numCorrect = 0, numTested = 0;
  vector<string> textInfoLines;

  try
  {
    const int nnuc = m_sourceModel->numNuclides();
    for( int i = 0; i < nnuc; ++i )
    {
      const SandiaDecay::Nuclide *nuc = m_sourceModel->nuclide( i );
      assert( nuc );
      
      const SandiaDecay::Nuclide *ageNuc = m_sourceModel->ageDefiningNuclide( nuc );
      const ShieldingSourceFitCalc::ModelSourceType sourceType = m_sourceModel->sourceType( i );
      const bool selfAttNuc = (sourceType == ShieldingSourceFitCalc::ModelSourceType::Intrinsic);
      
      // For self-attenuating shieldings, we'll just test the shielding thickness
      // For nuclides whose age is controlled by another nuclide, we dont need to test age.
      if( selfAttNuc || (ageNuc && (ageNuc != nuc)) )
        continue;
      
      if( m_sourceModel->fitActivity(i) )
      {
        boost::optional<double> truthAct = m_sourceModel->truthActivity(i);
        boost::optional<double> tolerance = m_sourceModel->truthActivityTolerance(i);
        
        if( !truthAct || !tolerance )
        {
          successful = false;
          textInfoLines.push_back( "Did not have truth value for " + nuc->symbol + " activity." );
          continue;
        }
        
        const double fitAct = m_sourceModel->activity(i);
        const bool closeEnough = (fabs(*truthAct - fitAct) < *tolerance);
        
        numTested += 1;
        numCorrect += closeEnough;
        
        textInfoLines.push_back( "For " + nuc->symbol + " fit activity "
                                + PhysicalUnits::printToBestActivityUnits(fitAct)
                                + DetectorPeakResponse::det_eff_geom_type_postfix(m_sourceModel->detType()) + " with the"
                                " truth value of "
                                + PhysicalUnits::printToBestActivityUnits(*truthAct)
                                + DetectorPeakResponse::det_eff_geom_type_postfix(m_sourceModel->detType())
                                + " and tolerance "
                                + PhysicalUnits::printToBestActivityUnits(*tolerance)
                                + DetectorPeakResponse::det_eff_geom_type_postfix(m_sourceModel->detType())
                                + (closeEnough ? " - within tolerance." : " - out of tolerance." )
                                );
      }//if( fit activity )
      
      if( m_sourceModel->fitAge(i) )
      {
        const boost::optional<double> truthAge = m_sourceModel->truthAge(i);
        const boost::optional<double> tolerance = m_sourceModel->truthAgeTolerance(i);
        
        if( !truthAge || !tolerance )
        {
          successful = false;
          textInfoLines.push_back( "Did not have truth value for " + nuc->symbol + " age." );
          continue;
        }
        
        const double fitAge = m_sourceModel->age(i);
        const bool closeEnough = (fabs(*truthAge - fitAge) < *tolerance);
        
        numTested += 1;
        numCorrect += closeEnough;
        
        textInfoLines.push_back( "For " + nuc->symbol + " fit age "
                                + PhysicalUnits::printToBestTimeUnits(fitAge) + " with the"
                                " truth value of "
                                + PhysicalUnits::printToBestTimeUnits(*truthAge)
                                + " and tolerance "
                                + PhysicalUnits::printToBestTimeUnits(*tolerance)
                                + (closeEnough ? " - within tolerance." : " - out of tolerance." )
                                );
      }//if( fit age )
    }//for( int i = 0; i < nnuc; ++i )
    
    
    for( WWidget *widget : m_shieldingSelects->children() )
    {
      ShieldingSelect *select = dynamic_cast<ShieldingSelect *>(widget);
      if( !select )
        continue;
      
      if( select->isGenericMaterial() )
      {
        if( select->fitArealDensity() )
        {
          if( !select->truthAD || !select->truthADTolerance )
          {
            successful = false;
            textInfoLines.push_back( "Did not have truth AD for generic shielding" );
            continue;
          }
          
          const double fitAD = select->arealDensity();
          const bool closeEnough = (fabs(*select->truthAD - fitAD) < *select->truthADTolerance);
          
          numTested += 1;
          numCorrect += closeEnough;
          
          textInfoLines.push_back( "For Generic Shielding fit AN " + std::to_string(fitAD)
                                  + " with the truth value of " + std::to_string(*select->truthAD)
                                  + " and tolerance " + std::to_string(*select->truthADTolerance)
                                  + (closeEnough ? " - within tolerance." : " - out of tolerance." )
                                  );
        }//if( fit AD )
        
        if( select->fitAtomicNumber() )
        {
          if( !select->truthAN || !select->truthANTolerance )
          {
            successful = false;
            textInfoLines.push_back( "Did not have truth AN for generic shielding" );
            continue;
          }
          
          const double fitAN = select->atomicNumber();
          const bool closeEnough = (fabs(*select->truthAN - fitAN) < *select->truthANTolerance);
          
          numTested += 1;
          numCorrect += closeEnough;
          
          textInfoLines.push_back( "For Generic Shielding fit AN " + std::to_string(fitAN)
                                  + " with the truth value of " + std::to_string(*select->truthAN)
                                  + " and tolerance " + std::to_string(*select->truthANTolerance)
                                  + (closeEnough ? " - within tolerance." : " - out of tolerance." )
                                  );
        }//if( fit AN )
      }else
      {
        shared_ptr<const Material> mat = select->material();
        if( !mat )
        {
          successful = false;
          textInfoLines.push_back( "There was an invalid material." );
          continue;
        }
        
        const auto geom = geometry();
        auto checkDimension = [select, geom, mat, &successful, &textInfoLines, &numTested, &numCorrect]( const int dim ){
          
          string labeltxt;
          bool fitdim = false;
          double fitValue = 0.0;
          boost::optional<double> *thicknessVal = nullptr;
          boost::optional<double> *toleranceVal = nullptr;
          
          switch( geom )
          {
            case GammaInteractionCalc::GeometryType::Spherical:
              labeltxt = "thickness";
              fitdim = select->fitThickness();
              fitValue = select->thickness();
              thicknessVal = &(select->truthThickness);
              toleranceVal = &(select->truthThicknessTolerance);
              break;
              
            case GammaInteractionCalc::GeometryType::CylinderEndOn:
            case GammaInteractionCalc::GeometryType::CylinderSideOn:
              if( dim == 0 )
              {
                labeltxt = "cyl. radius";
                fitdim = select->fitCylindricalRadiusThickness();
                fitValue = select->cylindricalRadiusThickness();
                thicknessVal = &(select->truthThickness);
                toleranceVal = &(select->truthThicknessTolerance);
              }else if( dim == 1 )
              {
                labeltxt = "cyl. length";
                fitdim = select->fitCylindricalLengthThickness();
                fitValue = select->cylindricalLengthThickness();
                thicknessVal = &(select->truthThicknessD2);
                toleranceVal = &(select->truthThicknessD2Tolerance);
              }else
              {
                assert( 0 );
                throw std::logic_error( "invalid dim" );
              }
              break;
              
            case GammaInteractionCalc::GeometryType::Rectangular:
              if( dim == 0 )
              {
                labeltxt = "rect. width";
                fitdim = select->fitRectangularWidthThickness();
                fitValue = select->rectangularWidthThickness();
                thicknessVal = &(select->truthThickness);
                toleranceVal = &(select->truthThicknessTolerance);
              }else if( dim == 1 )
              {
                labeltxt = "rect. height";
                fitdim = select->fitRectangularHeightThickness();
                fitValue = select->rectangularHeightThickness();
                thicknessVal = &(select->truthThicknessD2);
                toleranceVal = &(select->truthThicknessD2Tolerance);
              }else if( dim == 2 )
              {
                labeltxt = "rect. depth";
                fitdim = select->fitRectangularDepthThickness();
                fitValue = select->rectangularDepthThickness();
                thicknessVal = &(select->truthThicknessD3);
                toleranceVal = &(select->truthThicknessD3Tolerance);
              }else
              {
                assert( 0 );
                throw std::logic_error( "invalid dim" );
              }
              break;
              
            case GammaInteractionCalc::GeometryType::NumGeometryType:
              assert( 0 );
              break;
          }//switch( geometry() )
          
          assert( thicknessVal && toleranceVal );
          
          if( !fitdim || !thicknessVal || !toleranceVal )
            return;
          
          if( !(*thicknessVal) || !(*toleranceVal) )
          {
            successful = false;
            textInfoLines.push_back( "Missing truth " + labeltxt + " for shielding '" + mat->name + "'" );
            return;
          }
          
          const bool closeEnough = (fabs((**thicknessVal) - fitValue) < (**toleranceVal));
          
          numTested += 1;
          numCorrect += closeEnough;
          
          textInfoLines.push_back( "For shielding '" + mat->name + "' fit " + labeltxt + " "
                                  + PhysicalUnits::printToBestLengthUnits(fitValue,4)
                                  + " with the truth value of "
                                  + PhysicalUnits::printToBestLengthUnits(**thicknessVal,4)
                                  + " and tolerance "
                                  + PhysicalUnits::printToBestLengthUnits(**toleranceVal)
                                  + (closeEnough ? " - within tolerance." : " - out of tolerance." )
                                  );
        };//auto checkDimension
        
        
        switch( geometry() )
        {
          case GammaInteractionCalc::GeometryType::Spherical:
            checkDimension(0);
            break;
            
          case GammaInteractionCalc::GeometryType::CylinderEndOn:
          case GammaInteractionCalc::GeometryType::CylinderSideOn:
            checkDimension(0);
            checkDimension(1);
            break;
            
          case GammaInteractionCalc::GeometryType::Rectangular:
            checkDimension(0);
            checkDimension(1);
            checkDimension(2);
            break;
            
          case GammaInteractionCalc::GeometryType::NumGeometryType:
            break;
        }//switch( geometry() )
        
        
        {//Begin check mass-fractions
          auto checkMassFrac = [select, mat, &successful, &textInfoLines, &numTested, &numCorrect](
                                    const SandiaDecay::Element * const el,
                                    const SandiaDecay::Nuclide * const nuc, const double value ){
            assert( el );
            const auto el_pos = select->truthFitMassFractions.find(el);
            if( el_pos == end(select->truthFitMassFractions) )
            {
              successful = false;
              textInfoLines.push_back( "Missing truth mass-fraction for element " + (el ? el->symbol : string("nullptr")) );
              return;
            }
            
            const string nuc_name = (nuc ? nuc->symbol : string("other nucs"));
            const map<const SandiaDecay::Nuclide *,pair<double,double>> &nucs = el_pos->second;
            const auto truth_pos = nucs.find(nuc);
            if( truth_pos == end(nucs) )
            {
              successful = false;
              textInfoLines.push_back( "Missing truth mass-fraction for " + nuc_name );
              return;
            }
                                  
            const double truthval = truth_pos->second.first;
            const double truthtol = truth_pos->second.second;
            if( (truthval < 0.0) || (truthval > 1.0) || (truthtol < 0.0) || (truthtol > 1.0) )
            {
              successful = false;
              textInfoLines.push_back( "Invalid truth mass-fraction for " + nuc_name );
              return;
            }
            
            numTested += 1;
            const bool closeEnough = ((value >= (truthval - truthtol)) && (value <= (truthval + truthtol)));
            numCorrect += closeEnough;
            
            textInfoLines.push_back( "For shielding '" + mat->name + "' fit " + nuc_name
                                     + " to have mass fraction " + SpecUtils::printCompact(value, 5)
                                     + " with the truth value of " + SpecUtils::printCompact(truthval,5)
                                     + " and tolerance " + SpecUtils::printCompact(truthtol,5)
                                     + (closeEnough ? " - within tolerance." : " - out of tolerance." ) );
          };//checkMassFrac( ... )
          
          const map<const SandiaDecay::Element *,vector<ShieldingSelect::NucMasFrac>> gui_mass_fracs
                                        = select->sourceNuclideMassFractions();
          for( const auto &el_nucs : gui_mass_fracs )
          {
            const SandiaDecay::Element * const el = el_nucs.first;
            const vector<ShieldingSelect::NucMasFrac> &nucs = el_nucs.second;
            for( const auto &nuc_frac : nucs )
            {
              const SandiaDecay::Nuclide * const nuc = std::get<0>(nuc_frac);
              const double frac = std::get<1>(nuc_frac);
              const bool fit = std::get<2>(nuc_frac);
              if( fit )
                checkMassFrac( el, nuc, frac );
            }//for( const auto &prev_nuc_frac : select->truthFitMassFractions )
          }//for( const auto &el_nucs : gui_mass_fracs )
        }//End check mass-fractions
        
      }//if( generic material ) / else
    }//for( WWidget *widget : m_shieldingSelects->children() )
    
    successful = (successful && numTested);
  }catch( std::exception &e )
  {
    successful = false;
    textInfoLines.push_back( "Caught exception during testing: " + string(e.what()) );
  }//try / catch
  
  return tuple<bool,int,int,vector<string>>( successful, numCorrect, numTested, textInfoLines );
}//std::tuple<bool,int,int,std::vector<std::string>> testCurrentFitAgainstTruth();



#endif //INCLUDE_ANALYSIS_TEST_SUITE



void ShieldingSourceDisplay::multiNucsPerPeakChanged()
{
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && !undoRedo->isInUndoOrRedo() )
  {
    auto undo_redo = [](){
      ShieldingSourceDisplay *display = InterSpec::instance()->shieldingSourceFit();
      if( display )
      {
        display->m_multiIsoPerPeak->setChecked( !display->m_multiIsoPerPeak->isChecked() );
        display->multiNucsPerPeakChanged();
      }
    };
    
    undoRedo->addUndoRedoStep( undo_redo, undo_redo, "Multiple nuclides per peak changed" );
  }//if( undoRedo )
  
  m_clusterWidth->setDisabled( !m_multiIsoPerPeak->isChecked() );
  if( m_clusterWidth->label() )
    m_clusterWidth->label()->setDisabled( !m_multiIsoPerPeak->isChecked() );
  
  updateChi2Chart();
}//void multiNucsPerPeakChanged()


void ShieldingSourceDisplay::clusterWidthChanged()
{
  const double prev_width = m_photopeak_cluster_sigma;

  const double new_width = std::min( 10.0, std::max( 0.0, static_cast<double>(m_clusterWidth->value()) ) );
  if( new_width == m_clusterWidth->value() )
    m_clusterWidth->setValue( new_width );

  m_photopeak_cluster_sigma = new_width;

  
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && !undoRedo->isInUndoOrRedo() )
  {
    auto undo_redo = [prev_width, new_width]( const bool undo ){
      ShieldingSourceDisplay *display = InterSpec::instance()->shieldingSourceFit();
      if( display )
      {
        display->m_photopeak_cluster_sigma = undo ? prev_width : new_width;
        display->m_clusterWidth->setValue( display->m_photopeak_cluster_sigma );
        display->updateChi2Chart();
      }
    };
    
    undoRedo->addUndoRedoStep( [undo_redo](){undo_redo(true);},
                              [undo_redo](){undo_redo(false);},
                              "Cluster width changed" );
  }//if( undoRedo )

  updateChi2Chart();
}//void clusterWidthChanged()

void ShieldingSourceDisplay::attenuateForAirChanged()
{
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && !undoRedo->isInUndoOrRedo() )
  {
    auto undo_redo = [](){
      ShieldingSourceDisplay *display = InterSpec::instance()->shieldingSourceFit();
      if( display )
      {
        display->m_attenForAir->setChecked( !display->m_attenForAir->isChecked() );
        display->attenuateForAirChanged();
      }
    };
    
    undoRedo->addUndoRedoStep( undo_redo, undo_redo, "Attenuate for air changed" );
  }//if( undoRedo )
  
  updateChi2Chart();
}


void ShieldingSourceDisplay::backgroundPeakSubChanged()
{
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && !undoRedo->isInUndoOrRedo() )
  {
    auto undo_redo = [](){
      ShieldingSourceDisplay *display = InterSpec::instance()->shieldingSourceFit();
      if( display )
      {
        display->m_backgroundPeakSub->setChecked( !display->m_backgroundPeakSub->isChecked() );
        display->backgroundPeakSubChanged();
      }
    };
    
    undoRedo->addUndoRedoStep( undo_redo, undo_redo, "Background peak subtraction changed." );
  }//if( undoRedo )
  
  
  if( m_backgroundPeakSub->isChecked() )
  {
    std::shared_ptr<const SpecMeas> back = m_specViewer->measurment(SpecUtils::SpectrumType::Background);
    
    if( !back )
    {
      m_backgroundPeakSub->setUnChecked();
      passMessage( WString::tr("ssd-warn-no-back-spec"), WarningWidget::WarningMsgHigh );
      return;
    }//if( !back )
    
    std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaks;
    
    const set<int> &displayed = m_specViewer->displayedSamples(SpecUtils::SpectrumType::Background);
    peaks = back->peaks( displayed );
    
    if( !peaks || peaks->empty() )
    {
      m_backgroundPeakSub->setUnChecked();
      passMessage( WString::tr("ssd-warn-no-back-peaks-toggle"), WarningWidget::WarningMsgHigh );
      return;
    }//if( !peaks || peaks->empty() )
  }//if( m_backgroundPeakSub->isChecked() )
  
  updateChi2Chart();
}//void backgroundPeakSubChanged()


void ShieldingSourceDisplay::sameIsotopesAgeChanged()
{
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && !undoRedo->isInUndoOrRedo() )
  {
    auto undo_redo = [](){
      ShieldingSourceDisplay *display = InterSpec::instance()->shieldingSourceFit();
      if( display )
      {
        display->m_sameIsotopesAge->setChecked( !display->m_sameIsotopesAge->isChecked() );
        display->sameIsotopesAgeChanged();
      }
    };
    
    undoRedo->addUndoRedoStep( undo_redo, undo_redo, "Isotope age grouping changed." );
  }//if( undoRedo )
  
  
  m_sourceModel->setUseSameAgeForIsotopes( m_sameIsotopesAge->isChecked() );
  updateChi2Chart();
}//void sameIsotopesAgeChanged()


void ShieldingSourceDisplay::decayCorrectChanged()
{
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && !undoRedo->isInUndoOrRedo() )
  {
    auto undo_redo = [](){
      ShieldingSourceDisplay *display = InterSpec::instance()->shieldingSourceFit();
      if( display )
      {
        display->m_decayCorrect->setChecked( !display->m_decayCorrect->isChecked() );
        display->decayCorrectChanged();
      }
    };
    
    undoRedo->addUndoRedoStep( undo_redo, undo_redo, "Decay correct activity changed." );
  }//if( undoRedo )
  
  updateChi2Chart();
}//void decayCorrectChanged()


void ShieldingSourceDisplay::showGraphicTypeChanged()
{
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && !undoRedo->isInUndoOrRedo() )
  {
    auto undo_redo = [](){
      ShieldingSourceDisplay *display = InterSpec::instance()->shieldingSourceFit();
      if( display )
      {
        display->m_showChiOnChart->setChecked( !display->m_showChiOnChart->isChecked() );
        display->showGraphicTypeChanged();
      }
    };
    
    undoRedo->addUndoRedoStep( undo_redo, undo_redo, "Isotope age grouping changed." );
  }//if( undoRedo )
  
  
  const bool chi = m_showChiOnChart->isChecked();
  m_chi2Graphic->setShowChiOnChart( chi );
  updateChi2Chart();
}

void ShieldingSourceDisplay::handleUserDistanceChange()
{
  string distanceStr = m_distanceEdit->text().toUTF8();
  
  //Default to cm if no distance is given
  SpecUtils::trim( distanceStr );
  if( distanceStr.find_first_not_of( " \t0123456789.eE+-\n" ) == string::npos )
  {
    distanceStr += " cm";
    m_distanceEdit->setText( distanceStr );
  }
  
  try
  {
    const double distance = PhysicalUnits::stringToDistance( distanceStr );
    if( distance <= 0.0 )
      throw runtime_error( "Must have a non-zero distance" );
    
    UndoRedoManager *undoRedo = UndoRedoManager::instance();
    if( undoRedo && !undoRedo->isInUndoOrRedo() )
    {
      const string prevValue = m_prevDistStr;
      auto undo_redo = [prevValue, distanceStr](){
        UndoRedoManager *undoRedo = UndoRedoManager::instance();
        ShieldingSourceDisplay *display = InterSpec::instance()->shieldingSourceFit();
        if( display && undoRedo )
        {
          const string &value = undoRedo->isInUndo() ? prevValue : distanceStr;
          display->m_distanceEdit->setText( WString::fromUTF8(value) );
          display->handleUserDistanceChange();
        }
      };
      
      undoRedo->addUndoRedoStep( undo_redo, undo_redo, "Change source distance." );
    }//if( undoRedo )
    
    m_prevDistStr = distanceStr;
  }catch(...)
  {
    m_distanceEdit->setText( m_prevDistStr );
  }
  
  updateChi2Chart();
}//void ShieldingSourceDisplay::handleUserDistanceChange()


GeometryType ShieldingSourceDisplay::geometry() const
{
  int currentIndex = m_geometrySelect->currentIndex();
  if( currentIndex < 0 || currentIndex >= static_cast<int>(GeometryType::NumGeometryType) )
    currentIndex = static_cast<int>(GeometryType::Spherical);
  
  return GeometryType(currentIndex);
}//GeometryType geometry() const;


void ShieldingSourceDisplay::handleGeometryTypeChange()
{
  const GeometryType type = geometry();
  
  //cout << "ShieldingSourceDisplay::handleGeometryTypeChange(): Changing to "
  //     << GammaInteractionCalc::to_str(type) << endl;
  
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  unique_ptr<ShieldSourceChange> state_undo_creator;
  
  const bool is_same_geometry = (type == m_prevGeometry);
  
  if( !is_same_geometry && undoRedo && !undoRedo->isInUndoOrRedo() )
  {
    m_geometrySelect->setCurrentIndex( static_cast<int>(m_prevGeometry) );
    
    state_undo_creator = make_unique<ShieldSourceChange>( this, "Change geometry" );
    
    m_geometrySelect->setCurrentIndex( static_cast<int>(type) );
  }//if( (type != m_prevGeometry) && undoRedo && !undoRedo->isInUndoOrRedo() )
  
  m_prevGeometry = type;
  
  for( WWidget *widget : m_shieldingSelects->children() )
  {
    ShieldingSelect *select = dynamic_cast<ShieldingSelect *>(widget);
    assert( select );
    
    if( select )
      select->setGeometry( type );
  }//for( WWidget *widget : m_shieldingSelects->children() )
  
  if( !is_same_geometry )
  {
    try
    {
      checkDistanceAndThicknessConsistent();
    }catch( std::exception &e )
    {
      passMessage( e.what(), WarningWidget::WarningMsgMedium );
    }//try / catch
  }//if( !is_same_geometry )
  
  handleShieldingChange();
  
  updateChi2Chart(); //I think this should have already been called, but JIC since its a cheap call
}//void handleGeometryTypeChange()


void ShieldingSourceDisplay::checkDistanceAndThicknessConsistent()
{
  const char * const contained_err_msg_key = "ssd-update-dim-to-be-min";
  const char * const scaled_err_msg_key = "ssd-scaled-shield-to-be-less-dist";
 
  
  // TODO: We should probably check that we arent trying to fit multiple AN of generic shieldings,
  //       or the AD of two generic shieldings that have similar AN, or many other potentially
  //       degenerate cases.
  
  const GeometryType type = geometry();
  
  // Lets grab the shielding's we care about
  vector<ShieldingSelect *> shieldings;
  shieldings.reserve( m_shieldingSelects->children().size() );
  for( WWidget *widget : m_shieldingSelects->children() )
  {
    ShieldingSelect *select = dynamic_cast<ShieldingSelect *>(widget);
    assert( select );
    if( select && !select->isGenericMaterial() )
      shieldings.push_back( select );
  }//for( WWidget *widget : m_shieldingSelects->children() )
  
  
  // Dont let geometry dimensions be fit for, that we know we dont have the power to fit for
  //  e.g., rectangular width/height outside of last source shell
  int outer_source_shell = -1;
  for( int i = 0; i < static_cast<size_t>(shieldings.size()); ++i )
  {
    if( !shieldings[i]->traceSourceNuclides().empty()
       || !shieldings[i]->selfAttenNuclides().empty() )
      outer_source_shell = i;
  }//for( int i = 0; i < static_cast<size_t>(shieldings.size()); ++i )
  
  
  for( int i = 0; i < static_cast<size_t>(shieldings.size()); ++i )
  {
    const bool enable = (i <= outer_source_shell);
    switch( type )
    {
      case GammaInteractionCalc::GeometryType::Spherical:
        break;
        
      case GammaInteractionCalc::GeometryType::CylinderEndOn:
        shieldings[i]->setFitCylindricalRadiusEnabled( enable );
        break;
        
      case GammaInteractionCalc::GeometryType::CylinderSideOn:
        shieldings[i]->setFitCylindricalLengthEnabled( enable );
        break;
        
      case GammaInteractionCalc::GeometryType::Rectangular:
        shieldings[i]->setFitRectangularWidthEnabled( enable );
        shieldings[i]->setFitRectangularHeightEnabled( enable );
        break;
        
      case GammaInteractionCalc::GeometryType::NumGeometryType:
        assert( 0 );
        break;
    }//switch( type )
  }//for( int i = outer_source_shell + 1; i < static_cast<size_t>(shieldings.size()); ++i )
  
  
  
  // Now check shielding doesnt go past detector, and if it does, shrink things down.
  double shieldrad = 0.0, distance = 0.0;
  
  const double tolerance = PhysicalUnits::um;
  
  int shield_num = 0;
  bool updated_a_dim = false;
  for( ShieldingSelect *select : shieldings )
  {
    assert( type == select->geometry() );
    if( type != select->geometry() )
      throw runtime_error( "A shieldings geometry didnt match expected." );
    
    for( WLineEdit *edit : select->distanceEdits() )
      make_sure_distance_units_present( edit );
    
    // Check to make sure this shielding is larger than all shielding it contains
    switch( type )
    {
      case GeometryType::Spherical:
      {
        double thick = select->thickness();
        if( thick < tolerance )
        {
          updated_a_dim = true;
          thick = tolerance;
          select->setSphericalThickness( thick );
        }
        
        break;
      }//case GeometryType::Spherical:
      
      case GeometryType::CylinderEndOn:
      case GeometryType::CylinderSideOn:
      {
        const double min_delta = (shield_num ? 0.0 : tolerance);
        
        double rad = select->cylindricalRadiusThickness();
        double len = select->cylindricalLengthThickness();
        
        // If this is the first shielding, make sure it is at least 1 um rad/width.
        //  After that, it makes sense to maybe have one of the dimensions be the same length
        //  (e.x., in a hollow pipe, the inner void will be same length as metal tube).
        if( rad < min_delta )
        {
          rad = min_delta;
          updated_a_dim = true;
          select->setCylindricalRadiusThickness( min_delta );
        }
        
        if( len < min_delta )
        {
          len = min_delta;
          updated_a_dim = true;
          select->setCylindricalLengthThickness( min_delta );
        }
        
        // Dont let the dimensions both be zero, even after first shielding... I think the code
        //  would be fine if they are both zero, but this just doesnt quite seem right...
        //  TODO: need to think about implications of (not) letting, both dimensions be zero
        if( (rad == 0.0) && (len == 0.0) )
        {
          // I think we could leave
          len = rad = tolerance;
          updated_a_dim = true;
          select->setCylindricalLengthThickness( tolerance );
          select->setCylindricalLengthThickness( tolerance );
        }
        
        break;
      }//case GeometryType::CylinderSideOn:
        
      case GeometryType::Rectangular:
      {
        const double min_delta = (shield_num ? 0.0 : tolerance);
        
        double width = select->rectangularWidthThickness();
        double height = select->rectangularHeightThickness();
        double depth = select->rectangularDepthThickness();
        
        if( width < min_delta )
        {
          width = min_delta;
          updated_a_dim = true;
          select->setRectangularWidthThickness( min_delta );
        }
        
        if( height < min_delta )
        {
          height = min_delta;
          updated_a_dim = true;
          select->setRectangularHeightThickness( min_delta );
        }
        
        if( depth < min_delta )
        {
          depth = min_delta;
          updated_a_dim = true;
          select->setRectangularDepthThickness( min_delta );
        }
        
        // Dont let the dimensions both be zero, even after first shielding.
        //  See notes for cylindrical case.
        if( (width == 0.0) && (height == 0.0) && (depth == 0.0) )
        {
          width = height = depth = tolerance;
          updated_a_dim = true;
          select->setRectangularWidthThickness( tolerance );
          select->setRectangularHeightThickness( tolerance );
          select->setRectangularDepthThickness( tolerance );
        }
        
        break;
      }//case GeometryType::CylinderSideOn:
        
      case GeometryType::NumGeometryType:
        assert( 0 );
        break;
    }//switch( type )
    
    
    // Get the extent of the shielding, in the direction of the detector
    switch( type )
    {
      case GeometryType::Spherical:
        shieldrad += select->thickness();
        break;
        
      case GeometryType::CylinderEndOn:
        shieldrad += select->cylindricalLengthThickness();
        break;
        
      case GeometryType::CylinderSideOn:
        shieldrad += select->cylindricalRadiusThickness();
        break;
        
      case GeometryType::Rectangular:
        shieldrad += select->rectangularDepthThickness();
        break;
        
      case GeometryType::NumGeometryType:
        assert( 0 );
        break;
    }//switch( type )
    
    shield_num += 1;
  }//for( WWidget *widget : m_shieldingSelects->children() )
  
  const string distanceStr = m_distanceEdit->text().toUTF8();
  distance = PhysicalUnits::stringToDistance( distanceStr );
  if( distance <= 0.0 )
    throw runtime_error( "Distance must be greater than zero" );
  
  if( shieldrad <= distance )
  {
    if( updated_a_dim )
    {
      handleShieldingChange();
      
      throw runtime_error( WString::tr(contained_err_msg_key).toUTF8() );
    }
  }else
  {
    const double scale = 0.95*distance/shieldrad;
    for( ShieldingSelect *select : shieldings )
    {
      switch( type )
      {
        case GeometryType::Spherical:
          select->setSphericalThickness( scale * select->thickness() );
          break;
          
        case GeometryType::CylinderEndOn:
          select->setCylindricalLengthThickness( scale * select->cylindricalLengthThickness() );
          break;
          
        case GeometryType::CylinderSideOn:
          select->setCylindricalRadiusThickness( scale * select->cylindricalRadiusThickness() );
          break;
          
        case GeometryType::Rectangular:
          select->setRectangularDepthThickness( scale * select->rectangularDepthThickness() );
          break;
          
        case GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }//switch( type )
    }//for( WWidget *widget : m_shieldingSelects->children() )
    
    WString msg;
    if( updated_a_dim )
    {
      msg = WString("{1}<br />{2}")
        .arg(WString::tr(contained_err_msg_key))
        .arg( WString::tr(scaled_err_msg_key) );
    }else
    {
      msg = WString::tr(scaled_err_msg_key);
    }
    
    handleShieldingChange();
    
    throw runtime_error( msg.toUTF8() );
  }//if( shieldrad < distance )
}//void checkDistanceAndThicknessConsistent()


void ShieldingSourceDisplay::checkForMultipleGenericMaterials()
{
  int num_fit_for_ad = 0;
  for( WWidget *widget : m_shieldingSelects->children() )
  {
    ShieldingSelect *thisSelect = dynamic_cast<ShieldingSelect *>(widget);
    if( thisSelect && thisSelect->isGenericMaterial()
       && thisSelect->fitAtomicNumber() )
    {
      num_fit_for_ad += 1;
    }
  }//for( WWidget *widget : m_shieldingSelects->children() )
  
  if( num_fit_for_ad > 1 )
    throw runtime_error( WString::tr("ssd-err-fit-only-one-an").toUTF8() );
}//void checkForMultipleGenericMaterials()


void ShieldingSourceDisplay::checkAndWarnZeroMassFraction()
{
  //ToDo: will needs to implement, the point of this function will be to allow
  //  making sure the user hasnt selected a shielding to be a self attenuating
  //  source, but only allowed the source isotopes to be a mass fraction of zero.
  //  Should through std::exception with a use-friendly descriptive message on
  //  error.
  
  cerr << "ShieldingSourceDisplay::checkAndWarnZeroMassFraction() unimplemented!" << endl;
  
  
  
}//void checkAndWarnZeroMassFraction()


void ShieldingSourceDisplay::handleShieldingChange()
{
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  
  // A radius inside a trace or self-attenuating source could have changed, potentially changing
  //  the trace/intrinsic activity, so we'll go through and update every shielding widget, on every
  //  change to any of them
  for( WWidget *widget : m_shieldingSelects->children() )
  {
    ShieldingSelect *select = dynamic_cast<ShieldingSelect *>(widget);
    if( !select || select->isGenericMaterial() )
      continue;
    
    const vector<const SandiaDecay::Nuclide *> trace_srcs = select->traceSourceNuclides();
    const vector<const SandiaDecay::Nuclide *> self_atten_srcs = select->selfAttenNuclides();
    
    for( const SandiaDecay::Nuclide *nuc : trace_srcs )
      updateActivityOfShieldingIsotope( select, nuc );
    
    for( const SandiaDecay::Nuclide *nuc : self_atten_srcs )
      updateActivityOfShieldingIsotope( select, nuc );
  }//for( WWidget *widget : m_shieldingSelects->children() )
  
  updateChi2Chart();
}//void handleShieldingChange()


void ShieldingSourceDisplay::handleDetectorChanged( std::shared_ptr<DetectorPeakResponse> det )
{
  // `m_detectorDisplay` probably hasnt been updated yet, because of the order in signal/slot,
  //  so we'll update it now, even though everything wont be updated until render, so it
  //  would probably be fine anyway
  m_detectorDisplay->setDetector( det );
  
  unique_ptr<ShieldSourceChange> state_undo_creator;
  
  const bool fixed_geom = (det && det->isFixedGeometry());
  const DetectorPeakResponse::EffGeometryType det_type = det ? det->geometryType()
                                                 : DetectorPeakResponse::EffGeometryType::FarField;
  m_sourceModel->setDetectorType( det_type );
  
  if( fixed_geom /* && !m_distanceLabel->isHidden() */ )
  {
    // If there are any self-attenuating, or trace sources defined, or geometry is non-spherical,
    //  we should add an undo step, to go along with the change of detector undo.
    //  But this is a little tenuous, as it relies on the undo step of this function being inserted
    //  before the undo of the changed detector, and even then it will appear as an extra step to
    //  the user.
    
    bool any_volume_src = false;
    for( int nuc_num = 0; nuc_num < m_sourceModel->numNuclides(); ++nuc_num )
      any_volume_src |= m_sourceModel->isVolumetricSource( nuc_num );
    
    UndoRedoManager *undoRedo = UndoRedoManager::instance();
    if( undoRedo && undoRedo->canAddUndoRedoNow()
       && (any_volume_src || (geometry() != GeometryType::Spherical)) )
    {
      state_undo_creator = make_unique<ShieldSourceChange>( this, "Change to fixed-geometry DRF" );
    }
    
    m_attenForAir->hide();
    m_distanceLabel->hide();
    m_distanceEdit->hide();
    m_geometryLabel->hide();
    m_geometrySelect->hide();
    m_fixedGeometryTxt->show();
    
    m_geometrySelect->setCurrentIndex( static_cast<int>(GeometryType::Spherical) );
    handleGeometryTypeChange();
    
    for( WWidget *widget : m_shieldingSelects->children() )
    {
      ShieldingSelect *select = dynamic_cast<ShieldingSelect *>(widget);
      if( !select )
        continue;
      
      // If we were fitting a nuclide via trace source, or self-attenuation (fitting dimensions),
      //  then we will need to tell the nuclide model to fit for the activity there as well.
      set<const SandiaDecay::Nuclide *> nucs_fitting_act;
      
      if( !select->isGenericMaterial() )
      {
        for( const SandiaDecay::Nuclide *nuc : select->traceSourceNuclides() )
        {
          if( select->fitTraceSourceActivity(nuc) )
            nucs_fitting_act.insert( nuc );
        }
        
        bool fitting_dims = false;
        switch( select->geometry() )
        {
          case GammaInteractionCalc::GeometryType::Spherical:
            fitting_dims |= select->fitThickness();
            break;
            
          case GammaInteractionCalc::GeometryType::CylinderEndOn:
          case GammaInteractionCalc::GeometryType::CylinderSideOn:
            fitting_dims |= select->fitCylindricalRadiusThickness();
            fitting_dims |= select->fitCylindricalLengthThickness();
            break;
            
          case GammaInteractionCalc::GeometryType::Rectangular:
            fitting_dims |= select->fitRectangularWidthThickness();
            fitting_dims |= select->fitRectangularHeightThickness();
            fitting_dims |= select->fitRectangularDepthThickness();
            break;
            
          case GammaInteractionCalc::GeometryType::NumGeometryType:
            break;
        }//switch( select->geometry() )
        
        if( fitting_dims )
        {
          const vector<const SandiaDecay::Nuclide *> self_atten_nucs = select->selfAttenNuclides();
          nucs_fitting_act.insert( begin(self_atten_nucs), end(self_atten_nucs) );
        }//if( fitting_dims )
      }//if( !select->isGenericMaterial() )
      
      // This next call removed all intrinsic and trace sources
      select->setFixedGeometry( true );
      
      assert( select->traceSourceNuclides().empty() );
      assert( select->selfAttenNuclides().empty() );
      
      for( const SandiaDecay::Nuclide *nuc : nucs_fitting_act )
      {
        WModelIndex index = m_sourceModel->index( nuc, SourceFitModel::Columns::kFitActivity );
        m_sourceModel->setData( index, boost::any(true), Wt::CheckStateRole );
      }
    }//for( WWidget *widget : m_shieldingSelects->children() )
    
    
    const int nnuc = m_sourceModel->numNuclides();
    for( int nucn = 0; nucn < nnuc; ++nucn )
    {
      const SandiaDecay::Nuclide *nuc = m_sourceModel->nuclide(nucn);
      if( nuc )
        m_sourceModel->setSourceType( nuc, ShieldingSourceFitCalc::ModelSourceType::Point );
    }
  }//if( DRF is fixed geometry, but layout is for non-fixed geometry )
  
  if( !fixed_geom /* && m_distanceLabel->isHidden() */ )
  {
    m_attenForAir->show();
    m_distanceLabel->show();
    m_distanceEdit->show();
    m_geometryLabel->show();
    m_geometrySelect->show();
    m_fixedGeometryTxt->hide();
    
    for( WWidget *widget : m_shieldingSelects->children() )
    {
      ShieldingSelect *thisSelect = dynamic_cast<ShieldingSelect *>(widget);
      if( thisSelect )
        thisSelect->setFixedGeometry( false );
    }//for( WWidget *widget : m_shieldingSelects->children() )
  }//if( DRF is non-fixed geometry, but layout is for fixed geometry )
  
  updateChi2Chart();
}//void handleDetectorChanged()

void ShieldingSourceDisplay::updateChi2Chart()
{
  m_chi2ChartNeedsUpdating = true;
  scheduleRender(); //trigger re-render
}//void updateChi2Chart()


void ShieldingSourceDisplay::updateChi2ChartActual( std::shared_ptr<const ShieldingSourceFitCalc::ModelFitResults> results )
{
  try
  {
    checkDistanceAndThicknessConsistent();
  }catch( exception &e )
  {
    passMessage( e.what(), WarningWidget::WarningMsgHigh );
  }//try / catch
  
  try
  {
    if( m_logDiv )
    {
      m_logDiv->contents()->clear();
      m_logDiv->hide();
    }//if( m_logDiv )
    
    unsigned int ndof = 1;
    vector<GammaInteractionCalc::PeakResultPlotInfo> chis;
    
    if( results && results->peak_comparisons && results->peak_calc_details
       && (results->successful == ShieldingSourceFitCalc::ModelFitResults::FitStatus::Final) )
    {
      ndof = results->numDOF;
      chis = *results->peak_comparisons;
      m_calcLog = results->peak_calc_log;
      m_peakCalcLogInfo.reset( new vector<GammaInteractionCalc::PeakDetail>( *results->peak_calc_details ) );
    }else
    {
      m_calcLog.clear();
      m_peakCalcLogInfo.reset();
     
      auto fcnAndPars = shieldingFitnessFcn();
      
      std::shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn> &chi2Fcn = fcnAndPars.first;
      ROOT::Minuit2::MnUserParameters &inputPrams = fcnAndPars.second;
      
      ndof = inputPrams.VariableParameters();
      
      const vector<double> params = inputPrams.Params();
      const vector<double> errors = inputPrams.Errors();
      GammaInteractionCalc::ShieldingSourceChi2Fcn::NucMixtureCache mixcache;
      
      vector<GammaInteractionCalc::PeakDetail> calcLog;
      chis = chi2Fcn->energy_chi_contributions( params, errors, mixcache, &m_calcLog, &calcLog );
      
      m_peakCalcLogInfo.reset( new vector<GammaInteractionCalc::PeakDetail>( calcLog ) );
    }
    
    m_showLog->setDisabled( m_calcLog.empty() );

    vector<GammaInteractionCalc::PeakResultPlotInfo> keeper_points;

    for( size_t row = 0; row < chis.size(); ++row )
    {
      const double energy = chis[row].energy;
      const double chi = chis[row].numSigmaOff;
      const double scale = chis[row].observedOverExpected;
      const WColor &color = chis[row].peakColor;
      const double scale_uncert = chis[row].observedOverExpectedUncert;

      if( fabs(chi) < 1.0E5 && !IsInf(chi) && !IsNan(chi)
          && !IsInf(energy) && !IsNan(energy) )
        keeper_points.push_back( chis[row] );
    }//for( size_t row = 0; row < chis.size(); ++row )

    m_chi2Graphic->setNumFitForParams( ndof );
    if( !m_calcLog.empty() )
    {
      char buffer[64];
      snprintf( buffer, sizeof(buffer), "There %s %i parameter%s fit for",
                (ndof>1 ? "were" : "was"), int(ndof), (ndof>1 ? "s" : "") );
      m_calcLog.push_back( "&nbsp;" );
      m_calcLog.push_back( buffer );
    }//if( !m_calcLog.empty() )
    
    const int nrow = static_cast<int>( keeper_points.size() );
    const int nStartRows = m_chi2Model->rowCount();
    if( nStartRows < nrow )
      m_chi2Model->insertRows( nStartRows, nrow - nStartRows );

    if( nStartRows > nrow )
      m_chi2Model->removeRows( nrow, nStartRows - nrow );

    std::shared_ptr<const deque< PeakModel::PeakShrdPtr > > peaks = m_peakModel->peaks();
    
    for( int row = 0; row < nrow; ++row  )
    {
      const GammaInteractionCalc::PeakResultPlotInfo &p = keeper_points[row];
      
      const double &energy = p.energy;
      const double &chi = p.numSigmaOff;
      const double &scale = p.observedOverExpected;
      WColor color = p.peakColor;
      const double &scale_uncert = p.observedOverExpectedUncert;
      
      if( IsNan(energy) || IsInf(chi) )
      {
        passMessage( WString::tr("ssd-warn-invalid-chi2-for-energy").arg( round(100*energy)/100 ),
                    WarningWidget::WarningMsgHigh );
        continue;
      }//if( IsNan(p.second) || IsInf(p.second) )
      
      m_chi2Model->setData( row, 0, boost::any(energy) );
      m_chi2Model->setData( row, 1, boost::any(chi) );
      m_chi2Model->setData( row, 2, boost::any(scale) );
      
      if( color.isDefault() )
        color = m_specViewer->getColorTheme()->defaultPeakLine;
      color.setRgb( color.red(), color.green(), color.blue(), 255 );
      
      m_chi2Model->setData( row, 1, boost::any(color), Wt::MarkerPenColorRole );
      m_chi2Model->setData( row, 1, boost::any(color), Wt::MarkerBrushColorRole );
      m_chi2Model->setData( row, 2, boost::any(color), Wt::MarkerPenColorRole );
      m_chi2Model->setData( row, 2, boost::any(color), Wt::MarkerBrushColorRole );
      
      //If we wanted to include the nuclide in the model, we would have to loop
      //  over photopeaks in m_peakModel to try and match things up
      WString nuclidename;
      if( !!peaks )
      {
        for( const PeakModel::PeakShrdPtr &peak : *peaks )
        {
          if( peak->useForShieldingSourceFit()
              && peak->parentNuclide() && peak->decayParticle() )
          {
            if( fabs(energy - peak->decayParticle()->energy) < 0.001 )
            {
              nuclidename = peak->parentNuclide()->symbol;
              break;
            }
          }
        }//for( const PeakModel::PeakShrdPtr &p : *peaks )
        
        m_chi2Model->setData( row, 3, boost::any(nuclidename) );
        m_chi2Model->setData( row, 4, boost::any(scale_uncert) );
      }//if( !!peaks )
    }//for( int row = 0; row < nrow; ++row  )
  }catch( std::exception &e )
  {
    //One reason we may have made it here is if there are no peaks selected
    //Another is if we have removed using a peak from the peak, and it was an
    //  age defining, and we havent updated the dependant ages yet (but will in this
    //  event loop).
    if( m_chi2Model->rowCount() )
      m_chi2Model->removeRows( 0, m_chi2Model->rowCount() );
    cerr << "ShieldingSourceDisplay::updateChi2ChartActual()\n\tCaught:" << e.what() << endl;
  }
  
  
  const std::vector<WAbstractArea *> oldareas = m_chi2Graphic->areas();
  for( WAbstractArea *a : oldareas )
  {
    m_chi2Graphic->removeArea( a );
    delete a;
  }

//  cerr << "m_chi2Model->rowCount()=" << m_chi2Model->rowCount()
//       << "m_chi2Model->columnCount()=" << m_chi2Model->columnCount() << endl;

  m_calcLog.push_back( ns_no_uncert_info_txt );
}//void ShieldingSourceDisplay::updateChi2ChartActual()


void ShieldingSourceDisplay::showCalcLog()
{
  if( !m_logDiv )
  {
    m_logDiv = new AuxWindow( WString::tr("ssd-calc-log-window-title") );
    m_logDiv->contents()->addStyleClass( "CalculationLog" );
    m_logDiv->disableCollapse();
    m_logDiv->rejectWhenEscapePressed();
    //set min size so setResizable call before setResizable so Wt/Resizable.js wont cause the initial
    //  size to be the min-size
    m_logDiv->setMinimumSize( 640, 480 );
    m_logDiv->setResizable( true );
  }//if( !m_logDiv )
  
  m_logDiv->contents()->clear();
  vector<uint8_t> totaldata;
  for( const string &str : m_calcLog )
  {
    auto line = new WText( str, m_logDiv->contents() );
    line->setInline( false );
    totaldata.insert( end(totaldata), begin(str), end(str) );
    totaldata.push_back( static_cast<uint8_t>('\r') );
    totaldata.push_back( static_cast<uint8_t>('\n') );
  }
  
  m_logDiv->show();
  
  // Add a link to download this log file
  m_logDiv->footer()->clear();
  
  auto downloadResource = new WMemoryResource( "text/plain", m_logDiv->footer() );
  downloadResource->setData( totaldata );
  const int offset = wApp->environment().timeZoneOffset();
  const auto nowTime = WDateTime::currentDateTime().addSecs( 60 * offset );
  string filename = "act_shield_fit_" + nowTime.toString( "yyyyMMdd_hhmmss" ).toUTF8() + ".txt";
  downloadResource->suggestFileName( filename, WResource::DispositionType::Attachment );

#if( BUILD_AS_OSX_APP || IOS )
  WAnchor *logDownload = new WAnchor( WLink( downloadResource ), m_logDiv->footer() );
  logDownload->setStyleClass( "LinkBtn" );
  logDownload->setTarget( AnchorTarget::TargetNewWindow );
#else
  WPushButton *logDownload = new WPushButton( m_logDiv->footer() );
  logDownload->setIcon( "InterSpec_resources/images/download_small.svg" );
  logDownload->setLink( WLink( downloadResource ) );
  logDownload->setLinkTarget( Wt::TargetNewWindow );  //Note: we need to set new window after setLink, or else this wont actually get set
  logDownload->setStyleClass( "LinkBtn DownloadBtn" );
#endif
  
  logDownload->setText( WString::tr("ssd-calc-log-export-txt") );
  logDownload->setFloatSide( Wt::Side::Left );
    
#if( ANDROID )
  // Using hacked saving to temporary file in Android, instead of via network download of file.
  logDownload->clicked().connect( std::bind([downloadResource](){
    android_download_workaround(downloadResource, "fit_log.txt");
  }) );
#endif //ANDROID
  
  
  WPushButton *close = m_logDiv->addCloseButtonToFooter();
  close->clicked().connect( boost::bind( &AuxWindow::hide, m_logDiv ) );
  m_logDiv->finished().connect( this, &ShieldingSourceDisplay::closeCalcLogWindow );
  
  const int wwidth = m_specViewer->renderedWidth();
  const int wheight = m_specViewer->renderedHeight();
  m_logDiv->setMaximumSize( 0.8*wwidth, 0.8*wheight );
  m_logDiv->resizeToFitOnScreen();
  m_logDiv->centerWindow();
  
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    auto undo = [](){
      ShieldingSourceDisplay *shieldSourceFit = InterSpec::instance()->shieldingSourceFit();
      if( shieldSourceFit )
        shieldSourceFit->closeCalcLogWindow();
    };
    
    auto redo = [](){
      ShieldingSourceDisplay *shieldSourceFit = InterSpec::instance()->shieldingSourceFit();
      if( shieldSourceFit )
        shieldSourceFit->showCalcLog();
    };
    
    undoRedo->addUndoRedoStep( undo, redo, "Show calculation log." );
  }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
}//void showCalcLog()


void ShieldingSourceDisplay::closeCalcLogWindow()
{
  assert( m_logDiv );
  if( !m_logDiv )
    return;
  
  AuxWindow::deleteAuxWindow( m_logDiv );
  m_logDiv = nullptr;
  
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    auto undo = [](){
      ShieldingSourceDisplay *shieldSourceFit = InterSpec::instance()->shieldingSourceFit();
      if( shieldSourceFit )
        shieldSourceFit->showCalcLog();
    };
    
    auto redo = [](){
      ShieldingSourceDisplay *shieldSourceFit = InterSpec::instance()->shieldingSourceFit();
      if( shieldSourceFit )
        shieldSourceFit->closeCalcLogWindow();
    };
    
    undoRedo->addUndoRedoStep( undo, redo, "Close calculation log." );
  }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
}//void closeCalcLogWindow()


const ShieldingSelect *ShieldingSourceDisplay::innerShielding( const ShieldingSelect * const select ) const
{
  const ShieldingSelect *previous = nullptr;
  for( WWidget *widget : m_shieldingSelects->children() )
  {
    ShieldingSelect *thisSelect = dynamic_cast<ShieldingSelect *>(widget);
    if( !thisSelect || thisSelect->isGenericMaterial() )
      continue;

    if( thisSelect == select )
      return previous;

    previous = thisSelect;
  }//for( WWidget *widget : m_shieldingSelects->children() )

  throw runtime_error( "innerShielding(ShieldingSelect *): invalid select passed in" );
  
  return nullptr;
}//const ShieldingSelect *innerShielding( const ShieldingSelect * const select ) const


void ShieldingSourceDisplay::closeModelUploadWindow()
{
  AuxWindow::deleteAuxWindow( m_modelUploadWindow );
  m_modelUploadWindow = nullptr;
}//void closeModelUploadWindow()


void ShieldingSourceDisplay::finishModelUpload( WFileUpload *upload )
{
  ShieldSourceChange change( this, "Upload Activity/Shielding fit model." );
  
  shared_ptr<const string> pre_doc;
  if( change.m_pre_doc && !change.m_pre_doc->empty() )
  {
    pre_doc = change.m_pre_doc;
  }else
  {
    shared_ptr<string> xmldoc = make_shared<string>();
    change.doSerialization( *xmldoc );
    pre_doc = xmldoc;
  }
  
  try
  {
    const std::string filename = upload->spoolFileName();
    
    std::vector<char> data;
    SpecUtils::load_file_data( filename.c_str(), data );
    
    rapidxml::xml_document<char> new_doc;
    const int flags = rapidxml::parse_normalize_whitespace
                      | rapidxml::parse_trim_whitespace;
    new_doc.parse<flags>( &data.front() );
    deSerialize( new_doc.first_node() );
    
    m_modifiedThisForeground = true;
  }catch( std::exception &e )
  {
    passMessage( WString::tr("ssd-err-opening-import-model").arg(e.what()), 
                WarningWidget::WarningMsgHigh );
    
    try
    {
      ShieldSourceChange::doDeSerialization( this, pre_doc );
    }catch( std::exception & )
    {
      passMessage( "Even worse, there was an error trying to recover",
                   WarningWidget::WarningMsgHigh );
    }//try / catch
    
    // We wont insert an undo/redo step here, but this means if the user does an undo, then
    //  there will be a blank step (i.e., will try to close the upload window when there is
    //  none to close)
    // lets keep the ShieldSourceChange from inserting a undo/redo step.
    change.m_pre_doc.reset();
  }//try / catch
  
  closeModelUploadWindow();
}//void finishModelUpload(...)


void ShieldingSourceDisplay::modelUploadError( const ::int64_t size_tried )
{
  stringstream msg;
  const int max_size = static_cast<int>( wApp->maximumRequestSize() );
  msg << "Error uploading Source Shielding Fit Model.  Tried to upload "
      << size_tried << " (max size " << max_size << ")";
  passMessage( msg.str(), WarningWidget::WarningMsgHigh );
  
  closeModelUploadWindow();
}//void modelUploadError( const ::int64_t size_tried );


void ShieldingSourceDisplay::startModelUpload()
{
  if( m_modelUploadWindow )
    return;
  
  m_modelUploadWindow = new AuxWindow( WString::tr("ssd-import-model-window-title"),
                      (AuxWindowProperties::IsModal | AuxWindowProperties::TabletNotFullScreen) );
  
  WContainerWidget *contents = m_modelUploadWindow->contents();
  WFileUpload *upload = new WFileUpload( contents );
  upload->setInline( false );
  
  upload->uploaded().connect( boost::bind( &ShieldingSourceDisplay::finishModelUpload, this, upload ) );
  upload->fileTooLarge().connect( boost::bind( &ShieldingSourceDisplay::modelUploadError, this,
                                              boost::placeholders::_1 ) );
  upload->changed().connect( upload, &WFileUpload::upload );
  
  
  WPushButton *button = m_modelUploadWindow->addCloseButtonToFooter( WString::tr("Cancel") );
  button->clicked().connect( m_modelUploadWindow, &AuxWindow::hide );
  
  m_modelUploadWindow->centerWindow();
  m_modelUploadWindow->disableCollapse();
  
  m_modelUploadWindow->rejectWhenEscapePressed();
  m_modelUploadWindow->finished().connect( this, &ShieldingSourceDisplay::closeModelUploadWindow );

  m_modelUploadWindow->show();
  
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    auto undo = [](){
      ShieldingSourceDisplay *display = InterSpec::instance()->shieldingSourceFit();
      if( display )
        display->closeModelUploadWindow();
    };
    
    auto redo = [](){
      ShieldingSourceDisplay *display = InterSpec::instance()->shieldingSourceFit();
      if( display )
        display->startModelUpload();
    };
    
    undoRedo->addUndoRedoStep( undo, redo, "Open model upload window." );
  }//if( can insert undo/redo step )
}//void startModelUpload()


#if( USE_DB_TO_STORE_SPECTRA )
typedef Dbo::QueryModel< Dbo::ptr<ShieldingSourceModel> > QueryModel_t;

void updateDescription( WSelectionBox *select,
                        WSelectionBox *other_select,
                        WTextArea *summary,
                        WPushButton *button,
                        InterSpec *specViewer )
{
  button->disable();
  summary->setValueText( "" );
  if( other_select )
    other_select->clearSelection();
  
  const int row = select->currentIndex();
  if( row < 0 )
    return;
  
  QueryModel_t *querymodel = dynamic_cast<QueryModel_t *>( select->model() );
  if( !querymodel )
    throw runtime_error( "updateDescription(...): invalid input" );
  
  Dbo::ptr<ShieldingSourceModel> shieldmodel = querymodel->resultRow( row );
  if( !shieldmodel )
  {
    cerr << "updateDescription(): Unable to load selected model - sorry :(\n";
    return;
  }//if( !shieldmodel )
  
  WString descrip = shieldmodel->description;
  
  //XXX - I dont understand why the date isnt always valid, unless there is a
  //  bug in Wt with serializing dates to the database somehow.
  if( shieldmodel->creationTime.isValid() )
  {
    descrip += (WString("\nCreated ")
               + shieldmodel->creationTime.toString(DATE_TIME_FORMAT_STR));
    if( abs(shieldmodel->creationTime.secsTo( shieldmodel->serializeTime )) > 30 )
      descrip += (WString(", Saved ")
                 + shieldmodel->serializeTime.toString(DATE_TIME_FORMAT_STR));
  }//if( shieldmodel->creationTime.isValid() )
  
  
  {//begin interact with database
    std::shared_ptr<DataBaseUtils::DbSession> sql = specViewer->sql();
    DataBaseUtils::DbTransaction transaction( *sql );
    const Dbo::collection<Dbo::ptr<UserFileInDb> > &files
                                                  = shieldmodel->filesUsedWith;
    if( files.size() )
    {
      descrip += "\nSaved while working with spectra:";
      for( Dbo::collection<Dbo::ptr<UserFileInDb> >::const_iterator iter = files.begin();
           iter != files.end(); ++iter )
      {
        descrip += "\n    " + (*iter)->filename;
//                + " (" + (*iter)->serializeTime.toString(DATE_TIME_FORMAT_STR)
//                + ")";
      }//for( loop over files )
    }//if( files.size() )
    transaction.commit();
  }//end interact with database
  
  summary->setValueText( descrip );
  
  button->enable();
}//updateDescription(...)


void ShieldingSourceDisplay::removeModelFromDb( WSelectionBox *selec1,
                                                WSelectionBox *selec2 )
{
  WSelectionBox *selec = selec1;
  int row = selec ? selec->currentIndex() : -1;
  if( row < 0 )
  {
    selec = selec2;
    row = selec ? selec->currentIndex() : -1;
  }

  if( !selec || row < 0 )
    return;

  QueryModel_t *querymodel = dynamic_cast<QueryModel_t *>( selec->model() );  
  if( !querymodel )
    throw runtime_error( "removeModelFromDb(...): invalid input" );
  
  Dbo::ptr<ShieldingSourceModel> shieldmodel = querymodel->resultRow( row );
  if( !shieldmodel || !shieldmodel.session() )
  {
    cerr << "Error removeModelFromDb: shieldmodel=" << shieldmodel << " and "
         << " session " << shieldmodel.session() << endl;
    return;
  }

  {
    std::shared_ptr<DataBaseUtils::DbSession> sql = m_specViewer->sql();
    DataBaseUtils::DbTransaction transaction( *sql );
    shieldmodel.remove();
    transaction.commit();
    querymodel->reload();
  }
  
  
  if( selec1 )
    querymodel = dynamic_cast<QueryModel_t *>( selec1->model() );
  if( querymodel )
    querymodel->reload();
  if( selec2 )
    querymodel = dynamic_cast<QueryModel_t *>( selec2->model() );
  if( querymodel )
    querymodel->reload();
}//void ShieldingSourceDisplay::removeModelFromDb( WSelectionBox *selec )

bool ShieldingSourceDisplay::loadModelFromDb( Dbo::ptr<ShieldingSourceModel> shieldmodel )
{
  if( !shieldmodel )
    return false;
  
  rapidxml::xml_document<char> original_doc;
  try
  {
    serialize( &original_doc );
    rapidxml::xml_document<char> new_doc;
    const int flags = rapidxml::parse_normalize_whitespace
                      | rapidxml::parse_trim_whitespace;
    string data = shieldmodel->xmlData;
    if( data.size() )
      new_doc.parse<flags>( &(data[0]) );
    deSerialize( new_doc.first_node() );
  }catch( std::exception &e )
  {
    passMessage( WString::tr("ssd-err-opening-db-model").arg(e.what()), 
                WarningWidget::WarningMsgHigh );
    
    try
    {
      deSerialize( &original_doc );
    }catch( std::exception & )
    {
      passMessage( "Even worse, there was an error trying to recover",
                  WarningWidget::WarningMsgHigh );
    }//try / catch
    
    return false;
  }//try / catch
  
  return true;
}//void loadModelFromDb( Wt::Dbo::ptr<ShieldingSourceModel> entry )

void ShieldingSourceDisplay::finishLoadModelFromDatabase( WSelectionBox *first_selct,
                                                          WSelectionBox *other_select )
{
  ShieldSourceChange change( this, "Load model from database." );
  
  WSelectionBox *selec = first_selct;
  int row = selec ? selec->currentIndex() : -1;
  if( row < 0 )
  {
    selec = other_select;
    row = selec ? selec->currentIndex() : -1;
  }//if( row < 0 )
  
  if( !selec )
    throw runtime_error( "finishLoadModelFromDatabase(...): invalid slection box" );
  
  QueryModel_t *querymodel = dynamic_cast<QueryModel_t *>( selec->model() );
  
  if( row < 0 )
    return;
  if( !querymodel )
    throw runtime_error( "finishLoadModelFromDatabase(...): invalid input" );
    
  Dbo::ptr<ShieldingSourceModel> shieldmodel = querymodel->resultRow( row );
  if( !shieldmodel )
  {
    passMessage( WString::tr("ssd-err-loading-model"), WarningWidget::WarningMsgHigh );
    return;
  }//if( !shieldmodel )
  
  if( !loadModelFromDb( shieldmodel ) )
    change.m_pre_doc.reset(); //keep ShieldSourceChange from inserting undo/redo step
  
  closeBrowseDatabaseModelsWindow();
}//void finishLoadModelFromDatabase()


void ShieldingSourceDisplay::closeBrowseDatabaseModelsWindow()
{
  AuxWindow::deleteAuxWindow( m_modelDbBrowseWindow );
  m_modelDbBrowseWindow = nullptr;
}//void closeBrowseDatabaseModelsWindow();


void ShieldingSourceDisplay::startBrowseDatabaseModels()
{
  if( !m_specViewer || !m_specViewer->user() )
    throw runtime_error( "startBrowseDatabaseModels(): invalid user" );
  
  if( m_modelDbBrowseWindow )
    return;
  
  WTextArea *summary = NULL;
  WPushButton *accept = NULL, *cancel = NULL, *del = NULL;
  m_modelDbBrowseWindow = new AuxWindow( WString::tr("ssd-prev-saved-window-title"),
              (AuxWindowProperties::IsModal | AuxWindowProperties::TabletNotFullScreen) );
  m_modelDbBrowseWindow->finished().connect( this, &ShieldingSourceDisplay::closeBrowseDatabaseModelsWindow );
  
  try
  {
    WContainerWidget *contents = m_modelDbBrowseWindow->contents();
    summary = new WTextArea();
//    summary->setColumns( 30 );
    summary->setWidth( 316 );
    summary->setHeight( 50 );
    summary->disable();
    accept = new WPushButton( WString::tr("Load") );
    accept->disable();
  
    cancel = new WPushButton( WString::tr("Cancel") );
    cancel->clicked().connect( m_modelDbBrowseWindow, &AuxWindow::hide );
  
    del = new WPushButton( WString::tr("Delete") );
    del->setIcon( "InterSpec_resources/images/minus_min_white.png" );
    del->disable();

    Dbo::ptr<UserFileInDb> dbmeas;
    try
    {
      dbmeas = m_specViewer->measurementFromDb( SpecUtils::SpectrumType::Foreground, false );
    }catch( std::exception & )
    {
    }
    
    size_t nfileprev[2];
    WSelectionBox *selections[2] = { (WSelectionBox *)0, (WSelectionBox *)0 };
    
    {//begin codeblock for database interaction
      std::shared_ptr<DataBaseUtils::DbSession> sql = m_specViewer->sql();
      DataBaseUtils::DbTransaction transaction( *sql );
      nfileprev[0] = dbmeas ? dbmeas->modelsUsedWith.size() : 0;
      nfileprev[1] = m_specViewer->user()->shieldSrcModels().size();
      transaction.commit();
    }//end codeblock for database interaction
    
    
    for( size_t i = 0; i < 2; ++i )
    {
      if( !nfileprev[i] )
        continue;
      
     	WSelectionBox *selection = new WSelectionBox();
      selections[i] = selection;
      QueryModel_t *model = new QueryModel_t( selection );
      if( i == 0 )
        model->setQuery( dbmeas->modelsUsedWith.find() );
      else
        model->setQuery( m_specViewer->user()->shieldSrcModels().find() );
      model->addColumn( "Name" );
      selection->setModel( model );
      selection->setModelColumn( 0 );
      selection->setHeight( 100 );
      selection->setWidth( 320 );
      selection->setMaximumSize( 320, 100 );
      selection->setSelectionMode( Wt::SingleSelection );
      
      const char *msg_key = ((i==0) ? "ssd-models-prev-for-fore"
                                : "ssd-models-prev-saved");
      WText *title = new WText( WString::tr(msg_key), contents );
      title->setAttributeValue( "style", "font-weight:bold;" );
      title->setMargin( (i?8:12), Wt::Top );
      title->setInline( false );
      selection->setInline( false );
      selection->setMargin( 3, Wt::Top );
      selection->setMargin( 9, Wt::Left );
      contents->addWidget( selection );
    }//for( size_t i = 0; i < 2; ++i )
  
    if( selections[0] )
    {
      selections[0]->activated().connect( boost::bind( &updateDescription,
                selections[0], selections[1], summary, accept, m_specViewer ) );
      selections[0]->doubleClicked().connect(
            boost::bind( &ShieldingSourceDisplay::finishLoadModelFromDatabase,
                         this, selections[0], selections[1] ) );
      accept->clicked().connect(
              boost::bind( &ShieldingSourceDisplay::finishLoadModelFromDatabase,
                           this, selections[0], selections[1] ) );
    }//if( selections[0] )
    
    if( selections[1] )
    {
      selections[1]->activated().connect( boost::bind( &updateDescription,
                selections[1], selections[0], summary, accept, m_specViewer ) );
      selections[1]->doubleClicked().connect(
             boost::bind( &ShieldingSourceDisplay::finishLoadModelFromDatabase,
                          this, selections[1], selections[0] ) );
      accept->clicked().connect(
             boost::bind( &ShieldingSourceDisplay::finishLoadModelFromDatabase,
                          this, selections[1], selections[0] ) );
    }//if( selections[1] )
    
    del->clicked().connect(
                        boost::bind( &ShieldingSourceDisplay::removeModelFromDb,
                                       this, selections[0], selections[1] ) );
    
    if( nfileprev[0] || nfileprev[1] )
    {
      WText *title = new WText( WString::tr("ssd-model-desc") );
      title->setAttributeValue( "style", "font-weight:bold;margin-top:12px;" );

      WCheckBox *cb = new WCheckBox( WString::tr("ssd-cb-allow-del") );
      
      m_modelDbBrowseWindow->footer()->addWidget( cb );
      
      if( !m_specViewer->isMobile() )
      {
        cb->setFloatSide(Left);
        del->setFloatSide(Left);
      }
        
      m_modelDbBrowseWindow->footer()->addWidget( del );
        
      cb->checked().connect( del, &WPushButton::enable );
      cb->unChecked().connect( del, &WPushButton::disable );
      m_modelDbBrowseWindow->footer()->addWidget( cancel );

      m_modelDbBrowseWindow->footer()->addWidget( accept );

      title->setInline( false );
      summary->setInline( false );
      summary->setMargin( 3, Wt::Top );
      summary->setMargin( 9, Wt::Left );
      
      contents->addWidget( title );
      contents->addWidget( summary );
    }else
    {
      WText *info = new WText( WString::tr("ssd-no-model-in-db") );
      info->setInline( false );
      contents->addWidget( info );
      contents->addWidget( cancel );
      
      if( accept )
      {
        delete accept;
        accept = NULL;
      }//if( accept )
      if( summary )
      {
        delete summary;
        summary = NULL;
      }
    }//if( nfileprev[0] || nfileprev[1] )
    
    m_modelDbBrowseWindow->rejectWhenEscapePressed();
    
    m_modelDbBrowseWindow->setWidth( 350 );
    m_modelDbBrowseWindow->disableCollapse();
    m_modelDbBrowseWindow->centerWindow();
    m_modelDbBrowseWindow->show();
    
    UndoRedoManager *undoRedo = UndoRedoManager::instance();
    if( undoRedo && undoRedo->canAddUndoRedoNow() )
    {
      auto undo = [](){
        ShieldingSourceDisplay *display = InterSpec::instance()->shieldingSourceFit();
        if( display )
          display->closeBrowseDatabaseModelsWindow();
      };
      
      auto redo = [](){
        ShieldingSourceDisplay *display = InterSpec::instance()->shieldingSourceFit();
        if( display )
          display->startBrowseDatabaseModels();
      };
      
      undoRedo->addUndoRedoStep( undo, redo, "Show Activity/Shielding fit model browser." );
    }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
  }catch( std::exception &e )
  {
    if( accept )
      delete accept;
    if( cancel )
      delete cancel;
    if( summary )
     delete summary;
    if( del )
      delete del;
    
    AuxWindow::deleteAuxWindow( m_modelDbBrowseWindow );
    m_modelDbBrowseWindow = nullptr;
    
    passMessage( WString::tr("ssd-err-creating-db-browser"), WarningWidget::WarningMsgHigh );
    cerr << "\n\nShieldingSourceDisplay::startBrowseDatabaseModels() caught: "
         << e.what() << endl << endl;
  }//try / catch
}//void ShieldingSourceDisplay::startBrowseDatabaseModels()


Wt::Dbo::ptr<ShieldingSourceModel> ShieldingSourceDisplay::modelInDb()
{
  saveModelIfAlreadyInDatabase();
  if( m_modelInDb )
    return m_modelInDb;
  
  finishSaveModelToDatabase( defaultModelName(), defaultModelDescription() );
  
  return m_modelInDb;
}//Wt::Dbo::ptr<ShieldingSourceModel> modelInDb()
#endif //#if( USE_DB_TO_STORE_SPECTRA )

std::string ShieldingSourceDisplay::defaultModelName() const
{
  std::shared_ptr<const SpecMeas> meas = m_specViewer->measurment( SpecUtils::SpectrumType::Foreground );
  string name;
  if( meas && !meas->filename().empty() )
  {
    name = "1D model for " + meas->filename() + " ";
    name +=  WDateTime::currentDateTime().toString( "yyyyMMdd" ).toUTF8();
    if( name.length() > 255 )
      name = name.substr( 0, 255 );
  }//if( meas && !meas->filename().empty() )
  return name;
}//defaultModelName()


std::string ShieldingSourceDisplay::defaultModelDescription() const
{
  string descrip;
  for( int row = 0; row < m_sourceModel->rowCount(); ++row )
  {
    const SandiaDecay::Nuclide *nuc = m_sourceModel->nuclide( row );
    if( nuc )
    {
      if( row )
        descrip += ", ";
      descrip += nuc->symbol;
    }//if( nuc )
  }//for( int row = 0; row < m_sourceModel->rowCount(); ++row )
  
  descrip += "; shielded by ";
  
  int nshield = 0;
  for( WWidget *w : m_shieldingSelects->children() )
  {
    ShieldingSelect *select = dynamic_cast<ShieldingSelect *>( w );
    if( select && select->material() )
    {
      if( nshield )
        descrip += " and ";
      descrip += "'" + select->material()->name + "'";
      ++nshield;
    }//if( select )
  }//for( WWidget *w : m_shieldingSelects->children() )
  
  if( descrip.size() > 511 )
    descrip = descrip.substr( 0, 511 );
  
  return descrip;
}//std::string ShieldingSourceDisplay::defaultModelDescription() 

#if( USE_DB_TO_STORE_SPECTRA )
void ShieldingSourceDisplay::closeSaveModelToDatabaseWindow()
{
  AuxWindow::deleteAuxWindow( m_modelDbSaveWindow );
  m_modelDbSaveWindow = nullptr;
}//void closeSaveModelToDatabaseWindow()


void ShieldingSourceDisplay::startSaveModelToDatabase( bool prompt )
{
  if( m_modelInDb && !prompt )
  {
    finishGuiSaveModelToDatabase( nullptr, nullptr );
    return;
  }//if( m_modelInDb && !prompt )
  
  if( m_modelDbSaveWindow )
    return;
  
  m_modelDbSaveWindow = new AuxWindow( WString::tr("ssd-save-model-to-db-window-title"),
                  (AuxWindowProperties::IsModal
                   | AuxWindowProperties::TabletNotFullScreen
                   | AuxWindowProperties::DisableCollapse) );
  WContainerWidget *contents = m_modelDbSaveWindow->contents();
  m_modelDbSaveWindow->centerWindow();
  
  m_modelDbSaveWindow->rejectWhenEscapePressed();
  m_modelDbSaveWindow->finished().connect( this, &ShieldingSourceDisplay::closeSaveModelToDatabaseWindow );
  m_modelDbSaveWindow->show();
 
  WLabel *label = new WLabel( WString::tr("ssd-save-to-db-name"), contents );
  label->setInline( false );

  WLineEdit *nameEdit = new WLineEdit( contents );
  nameEdit->setAttributeValue( "ondragstart", "return false" );
  nameEdit->setInline( false );
  if( m_modelInDb )
    nameEdit->setValueText( m_modelInDb->name );

  if( nameEdit->valueText().empty() )
    nameEdit->setValueText( defaultModelName() );
  
  label = new WLabel( WString::tr("ssd-save-to-db-desc"), contents );
  label->setInline( false );
  
  WLineEdit *descEdit = new WLineEdit( contents );
  descEdit->setAttributeValue( "ondragstart", "return false" );
  descEdit->setInline( false );
  if( m_modelInDb )
    descEdit->setValueText( m_modelInDb->description );
  
  if( descEdit->valueText().empty() )
    descEdit->setValueText( defaultModelDescription() );
  
  nameEdit->enterPressed().connect( boost::bind( &WFormWidget::setFocus, descEdit, true ) );
  
  descEdit->setTextSize( 32 );
  nameEdit->setTextSize( 32 );

 

  WPushButton *button = new WPushButton( WString::tr("Save"), m_modelDbSaveWindow->footer() );
  button->setIcon( "InterSpec_resources/images/disk2.png" );
  
  button->clicked().connect(
              boost::bind( &ShieldingSourceDisplay::finishGuiSaveModelToDatabase,
                           this, nameEdit, descEdit ) );
  descEdit->enterPressed().connect( boost::bind( &WFormWidget::setFocus, button, true ) );
  
  
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    auto undo = [](){
      ShieldingSourceDisplay *shieldSourceFit = InterSpec::instance()->shieldingSourceFit();
      if( shieldSourceFit )
        shieldSourceFit->closeSaveModelToDatabaseWindow();
    };
    
    auto redo = [prompt](){
      ShieldingSourceDisplay *shieldSourceFit = InterSpec::instance()->shieldingSourceFit();
      if( shieldSourceFit )
        shieldSourceFit->startSaveModelToDatabase(prompt);
    };
    
    undoRedo->addUndoRedoStep( undo, redo, "Show model save to database dialog." );
  }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
  
}//void startSaveModelToDatabase()


void ShieldingSourceDisplay::finishGuiSaveModelToDatabase( WLineEdit *name_edit,
                                                        WLineEdit *desc_edit )
{
  if( !m_modelDbSaveWindow )
    return;
  
  if( name_edit && name_edit->valueText().empty() )
  {
    WText *txt = new WText( WString::tr("ssd-must-enter-name"), m_modelDbSaveWindow->contents() );
    txt->setInline( false );
    txt->setAttributeValue( "style", "color:red;" );
    return;
  }//if( name_edit && name_edit->valueText().empty() )
  
  //Check that the name is unique
  //  I'm not sure if I want to actually enforce this, or at a minimum should
  //  give a way to overide this
//  if( !m_modelInDb )
//  {
//    DbTransaction transaction( viewer );
//    const size_t nexisting = m_specViewer->m_user->m_shieldSrcModels.find()
//                            .where( "Name = ?" ).bind( name_edit->valueText() )
//                            .resultList().size();
//    transaction.commit();
//    
//    if( nexisting )
//    {
//      vector<WWidget *> children = window->contents()->children();
//      WText *txt = children.size() ? dynamic_cast<WText *>(children.back()) : 0;
//      if( !txt )
//      {
//        txt = new WText( window->contents() );
//        txt->setInline( false );
//        txt->setAttributeValue( "style", "color:red;" );
//      }//if( !txt )
//      txt->setText( "You already have a model with this name" );
//      return;
//    }//if( nexisting )
//  }//if( check if a unique name )
  
  WString name, description;
  if( name_edit )
    name = name_edit->valueText();
  if( desc_edit )
    description = desc_edit->valueText();
  
  finishSaveModelToDatabase( name, description );
  
  closeSaveModelToDatabaseWindow();
}//finishGuiSaveModelToDatabase(...)


bool ShieldingSourceDisplay::finishSaveModelToDatabase( const Wt::WString &name,
                                                       const Wt::WString &desc )
{
  ShieldingSourceModel *model = NULL;
  
  std::shared_ptr<DataBaseUtils::DbSession> sql = m_specViewer->sql();
  DataBaseUtils::DbTransaction transaction( *sql );

  try
  {
    if( !m_modelInDb )
    {
      model = new ShieldingSourceModel();
      model->user = m_specViewer->user();
      model->serializeTime = WDateTime::currentDateTime();
      m_modelInDb.reset( new ShieldingSourceModel() );
      m_modelInDb = sql->session()->add( model );
    }//if( m_modelInDb ) / else
    
    model = m_modelInDb.modify();
    model->serializeTime = WDateTime::currentDateTime();
    
    if( !name.empty() )
      model->name = name;
    if( !desc.empty() )
      model->description = desc;
    
    const string utfname = model->name.toUTF8();
    const string utfdesc = model->description.toUTF8();
    if( utfname.size() > 255 )
      model->name.fromUTF8( utfname.substr( 0, 255 ) );
    if( utfdesc.size() > 511 )
      model->description.fromUTF8( utfdesc.substr( 0, 511 ) );
    
    rapidxml::xml_document<char> doc;
    serialize( &doc );
    model->xmlData.clear();
    rapidxml::print(std::back_inserter(model->xmlData), doc, 0);

    Dbo::ptr<UserFileInDb> dbmeas;
    try
    {
      dbmeas = m_specViewer->measurementFromDb( SpecUtils::SpectrumType::Foreground, false );
    }catch( std::exception & )
    {
    }
    
    if( dbmeas )
      model->filesUsedWith.insert( dbmeas );
    
    transaction.commit();
    m_saveAsNewModelInDb->enable();
  }catch( std::exception & )
  {
    m_modelInDb.reset();
    m_saveAsNewModelInDb->disable();
    transaction.rollback();
    return false;
  }//try / catch
  
  return true;
}//void finishSaveModelToDatabase( AxuWindow *window, Wt::WLineEdit *edit )


void ShieldingSourceDisplay::saveCloneModelToDatabase()
{
  if( !m_modelInDb )
    return;
  
  std::shared_ptr<DataBaseUtils::DbSession> sql = m_specViewer->sql();
  DataBaseUtils::DbTransaction transaction( *sql );

  try
  {
    m_modelInDb.reread();
    ShieldingSourceModel *model = new ShieldingSourceModel();
    model->shallowEquals( *m_modelInDb );
    
    model->user = m_specViewer->user();
    model->serializeTime = WDateTime::currentDateTime();
    model->name = model->name + " Clone";
    
    const string utfname = model->name.toUTF8();
    if( utfname.size() > 255 )
      model->name.fromUTF8( utfname.substr(0,255) );
    
    m_modelInDb = sql->session()->add( model );
    
    Dbo::ptr<UserFileInDb> dbmeas;
    try
    {
      dbmeas = m_specViewer->measurementFromDb( SpecUtils::SpectrumType::Foreground, false );
    }catch( std::exception & )
    {
    }
    
    if( dbmeas )
      m_modelInDb.modify()->filesUsedWith.insert( dbmeas );
    transaction.commit();
  }catch( std::exception &e )
  {
    m_modelInDb.reset();
    m_saveAsNewModelInDb->disable();
    cerr << "\n\nException caught in ShieldingSourceDisplay::saveCloneModelToDatabase(): "
         << e.what() << endl;
    transaction.rollback();
    passMessage( WString::tr("ssd-err-saving-to-db"), WarningWidget::WarningMsgHigh );
    return;
  }//try / catch
  
  cerr << "\n\nFinishing in ShieldingSourceDisplay::saveCloneModelToDatabase()" << endl;
  startSaveModelToDatabase( true );
}//void saveCloneModelToDatabase()


void ShieldingSourceDisplay::saveModelIfAlreadyInDatabase()
{
  if( !m_modelInDb )
    return;
  
  std::shared_ptr<DataBaseUtils::DbSession> sql = m_specViewer->sql();
  DataBaseUtils::DbTransaction transaction( *sql );
  
  try
  {
//    m_modelInDb.reread();
    m_modelInDb.modify()->serializeTime = WDateTime::currentDateTime();
    rapidxml::xml_document<char> doc;
    serialize( &doc );
    m_modelInDb.modify()->xmlData.clear();
    rapidxml::print(std::back_inserter(m_modelInDb.modify()->xmlData), doc, 0);
    transaction.commit();
    cerr << "\n\nFinishing in ShieldingSourceDisplay::saveModelIfAlreadyInDatabase()" << endl;
  }catch( std::exception &e )
  {
    cerr << "\n\nException caught in ShieldingSourceDisplay::saveModelIfAlreadyInDatabase(): "
    << e.what() << endl;
    transaction.rollback();
  }//try / catch
}//void saveModelIfAlreadyInDatabase()
#endif //#if( USE_DB_TO_STORE_SPECTRA )


void ShieldingSourceDisplay::testSerialization()
{
  rapidxml::xml_document<char> doc;
  
  cerr << "Staring serialization" << endl;
  serialize( &doc );
  std::string s;
  rapidxml::print(std::back_inserter(s), doc, 0);
  cerr << "Serialized:\n" << doc << endl;
  cerr << "\n\nStarting de-serialization" << endl;
  deSerialize( &doc );
  cerr << "Done deserializing\n\n" << endl;
}//testSerialization()





void ShieldingSourceDisplay::deSerializePeaksToUse(
                                    const rapidxml::xml_node<char> *peaks_node )
{
  const rapidxml::xml_attribute<char> *attr;
  const rapidxml::xml_node<char> *symbol_node, *energy_node, *peak_node;
  
  for( peak_node = peaks_node->first_node( "Peak", 4 );
       peak_node; peak_node = peak_node->next_sibling( "Peak", 4 ) )
  {
    attr = peak_node->first_attribute( "Use", 3 );
    symbol_node = peak_node->first_node( "Nuclide", 7 );
    energy_node = peak_node->first_node( "Energy", 6 );
    
    bool use = true;
    double energy;
    if( attr && attr->value() && !(stringstream(attr->value())>>use) )
      throw runtime_error( "Invalid 'Use' attribute in Peak XML element" );
    if( !use )
      continue;
    if( !symbol_node || !energy_node
       || !symbol_node->value() || !energy_node->value()
       || !(stringstream(energy_node->value()) >> energy) )
      throw runtime_error( "Invalid or missing node for Peak XML element" );
    
    //    PeakModel::PeakShrdPtr nearestPeak = m_peakModel->nearestPeak( energy );
    
    double nearestE = DBL_MAX;
    WModelIndex nearest_index;
    PeakModel::PeakShrdPtr nearest_peak;
    for( int peakn = 0; peakn < m_peakModel->rowCount(); ++peakn )
    {
      WModelIndex index = m_peakModel->index( peakn, PeakModel::kUseForShieldingSourceFit );
      PeakModel::PeakShrdPtr peak = m_peakModel->peak( index );
      const SandiaDecay::Nuclide *nuclide = peak->parentNuclide();
      
      try
      {
        const float gamenergy = peak->gammaParticleEnergy();
        const double dE = fabs(gamenergy - energy);
        if( nuclide && (dE < nearestE) )
        {
          nearestE = dE;
          nearest_peak = peak;
          nearest_index = index;
        }//if( this is nearest candidate peak )
      }catch( std::exception & )
      {
        
      }//try / catch
    }//for( int peakn = 0; peakn < m_peakModel->rowCount(); ++peakn )
    
    if( nearest_peak && (nearestE < 0.1) )
      m_peakModel->setData( nearest_index, true );
  }//for( node = peaks_node->first_node(); node; node = node->next_sibling() )
}//void deSerializePeaksToUse( rapidxml::xml_node<char> *peaks_node )


void ShieldingSourceDisplay::deSerializeSourcesToFitFor( const rapidxml::xml_node<char> *sources )
{
  vector<ShieldingSourceFitCalc::IsoFitStruct> &model_isos = m_sourceModel->m_nuclides;
  vector<bool> model_row_modded( model_isos.size(), false );
  
  for( const rapidxml::xml_node<> *src_node = sources->first_node( "Nuclide", 7 );
       src_node; src_node = src_node->next_sibling( "Nuclide", 7 ) )
  {
    ShieldingSourceFitCalc::IsoFitStruct row;
    row.deSerialize( src_node );
    
    for( size_t i = 0; i < model_isos.size(); ++i )
    {
      ShieldingSourceFitCalc::IsoFitStruct &cand = model_isos[i];
      if( model_row_modded[i] || (cand.nuclide != row.nuclide) )
        continue;
      row.numProgenyPeaksSelected = cand.numProgenyPeaksSelected;
      cand = row;
      model_row_modded[i] = true;
    }//for( size_t i = 0; i < m_sourceModel->m_nuclides.size(); ++i )
  }//for( loop over source isotopes )
  
  
  //We should have modded all rows that are currently in m_sourceModel, but
  //  we're not garunteed all rows that where in the XML are also in m_sourceModel
  //Note 20200929: this check below is a little over-optimistic, so I commented it out; if you
  //  stumble on this code in like a year, and its still commentd out, just delete it, because:
  //  If m_sourceModel had "Tl201" because that is the peaks in our current spectrum (and that is
  //  the peak that is marked to be used in a fit), but the XML model we're loading had "I131" in
  //  it (because thats what we fit for in the previous record of the spectrum file), then
  //  model_row_modded for this model will still be false.
  //for( const bool modded : model_row_modded )
  //{
  //  if( !modded )
  //    throw runtime_error( "Inconsistent state of source model and serialized from XML" );
  //}//for( const bool modded : model_row_modded )
  
  if( m_sourceModel->rowCount() )
  {
    WModelIndex topLeft = m_sourceModel->index( 0, 0 );
    WModelIndex bottomRight = m_sourceModel->index( m_sourceModel->rowCount()-1,
                                              m_sourceModel->columnCount()- 1 );
    m_sourceModel->dataChanged().emit( topLeft, bottomRight );
  }//if( m_sourceModel->rowCount() )
}//void deSerializeSourcesToFitFor( rapidxml::xml_node<char> *sources );


void ShieldingSourceDisplay::deSerializeShieldings( const rapidxml::xml_node<char> *shiledings )
{
  const rapidxml::xml_node<> *shield_node;
  
  const shared_ptr<const DetectorPeakResponse> det = m_detectorDisplay->detector();
  
  // If our Detector Peak Efficiency is for fixed-geometry, we wont allow intrinsic
  //  or trace sources
  const bool is_fixed_geom = det && det->isFixedGeometry();
  
  for( shield_node = shiledings->first_node( "Shielding", 9 );
      shield_node; shield_node = shield_node->next_sibling( "Shielding", 9 ) )
  {
    ShieldingSelect *select = addShielding( nullptr, false );
    select->deSerialize( shield_node, is_fixed_geom );
    
#if( PERFORM_DEVELOPER_CHECKS )
    for( const SandiaDecay::Nuclide *nuc : select->traceSourceNuclides() )
    {
      const int nucIndex = m_sourceModel->nuclideIndex(nuc);
      const bool guiFit = select->fitTraceSourceActivity(nuc);
      const bool modelFit = m_sourceModel->fitActivity( nucIndex );
      assert( guiFit == modelFit );
      
      const double modelActivity = m_sourceModel->activity(nucIndex);
      const double guiActivity = select->traceSourceTotalActivity(nuc);
      assert( fabs(modelActivity - guiActivity) < 0.0001*std::max(modelActivity, guiActivity) );
    }//for( const SandiaDecay::Nuclide *nuc : select->traceSourceNuclides() )
#endif
  }//for( loop over shieldings )
  
}//void deSerializeShieldings( const rapidxml::xml_node<char> *shiledings )


void ShieldingSourceDisplay::reset()
{
  m_modifiedThisForeground = false;
  m_prevDistStr = "100 cm";
  m_distanceEdit->setValueText( m_prevDistStr );
  const vector<WWidget *> shieldings = m_shieldingSelects->children();
  for( WWidget *child : shieldings )
  {
    const ShieldingSelect *select = dynamic_cast<ShieldingSelect *>( child );
    if( select )
      delete select;
  }//for( WWebWidget *child : shieldings )
  
  for( int peakn = 0; peakn < m_peakModel->rowCount(); ++peakn )
  {
    WModelIndex index = m_peakModel->index( peakn,
                                         PeakModel::kUseForShieldingSourceFit );
    const PeakModel::PeakShrdPtr &peak = m_peakModel->peak( index );
    if( peak->useForShieldingSourceFit() )
      m_peakModel->setData( index, false );
  }//for( int peakn = 0; peakn < m_peakModel->rowCount(); ++peakn )
  
  //We shouldnt actually need the next line
//  m_sourceModel->repopulateIsotopes();
  handleDetectorChanged( m_detectorDisplay->detector() );
}//void reset()


void ShieldingSourceDisplay::newForegroundSet()
{
  m_modifiedThisForeground = false;
}//void newForegroundSet()


bool ShieldingSourceDisplay::userChangedDuringCurrentForeground() const
{
  return m_modifiedThisForeground;
}//bool userChangedDuringCurrentForeground() const;


void ShieldingSourceDisplay::deSerialize( const rapidxml::xml_node<char> *base_node )
{
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  
  const rapidxml::xml_node<char> *geom_node, *muti_iso_node, *atten_air_node, *back_sub_node,
                                 *peaks_node, *isotope_nodes, *shieldings_node,
                                 *dist_node, *same_age_node, *decay_corr_node, *chart_disp_node;
  const rapidxml::xml_attribute<char> *attr;
  
  if( !base_node )
    throw runtime_error( "No ShieldingSourceFit node" );
  
  if( !rapidxml::internal::compare( base_node->name(), base_node->name_size(), "ShieldingSourceFit", 18, true) )
    throw runtime_error( "ShieldingSourceDisplay::deSerialize: invalid node name passed in: '"
                        + std::string(base_node->name(),base_node->name()+base_node->name_size()) + "'" );
  
  geom_node       = base_node->first_node( "Geometry", 8 );
  muti_iso_node   = base_node->first_node( "MultipleIsotopesPerPeak", 23 );
  atten_air_node  = base_node->first_node( "AttenuateForAir", 15 );
  back_sub_node   = base_node->first_node( "BackgroundPeakSubtraction", 25 );
  same_age_node   = base_node->first_node( "SameAgeIsotopes", 15 );
  decay_corr_node = base_node->first_node( "DecayCorrect", 12 );
  chart_disp_node = base_node->first_node( "ShowChiOnChart", 14 );
  peaks_node      = base_node->first_node( "Peaks", 5 );
  isotope_nodes   = base_node->first_node( "Nuclides", 8 );
  shieldings_node = base_node->first_node( "Shieldings", 10 );
  dist_node       = base_node->first_node( "Distance", 8 );
  
  if( !peaks_node || !isotope_nodes || !shieldings_node )
    throw runtime_error( "Missing necessary XML node" );
  
  int version;
  bool show_chi_on_chart = true;
  
  attr = base_node->first_attribute( "version", 7 );
  if( !attr || !attr->value() || !(stringstream(attr->value())>>version) )
    throw runtime_error( "Deserializing requires a version" );
  
  if( version != sm_xmlSerializationMajorVersion )
    throw runtime_error( "Invalid version of ShieldingSourceDisplay XML" );
  
  // Note that version is either "0" (implied minor version 0), or "0.1" or something; we could read
  //  in the minor version to sm_xmlSerializationMinorVersion and use it, but theres really no need.
  GeometryType geom_type = GeometryType::Spherical;
  if( geom_node && geom_node->value() )
  {
    bool found = false;
    const string val( geom_node->value(), geom_node->value() + geom_node->value_size() );
    
    for( GeometryType i = GeometryType(0);
        !found && (i != GeometryType::NumGeometryType);
        i = GeometryType(static_cast<int>(i) + 1) )
    {
      found = (val == GammaInteractionCalc::to_str(i));
      if( found )
        geom_type = i;
    }//for( loop over GeometryType )
    
    if( !found )
      throw runtime_error( "Invalid geometry specified in XML: " + val );
  }//if( geom_node )
  
  
  ShieldingSourceFitCalc::ShieldingSourceFitOptions options;
  options.deSerialize( base_node );
  
  options.background_peak_subtract = (options.background_peak_subtract
                                 && m_specViewer->measurment(SpecUtils::SpectrumType::Background));
  
  m_multithread_computation = options.multithread_self_atten;
  m_photopeak_cluster_sigma = options.photopeak_cluster_sigma;
  
  if( chart_disp_node && chart_disp_node->value() )
  {
    if( !(stringstream(chart_disp_node->value()) >> show_chi_on_chart) )
      throw runtime_error( "Invalid ShowChiOnChart node" );
  }
  
  if( !dist_node || !dist_node->value() )
    throw runtime_error( "Invalid or missing Distance node" );
  
  //clear out the GUI
  reset();
  
  
  m_prevGeometry = geom_type;
  m_geometrySelect->setCurrentIndex( static_cast<int>(geom_type) );
  
  m_multiIsoPerPeak->setChecked( options.multiple_nucs_contribute_to_peaks );
  m_clusterWidth->setValue( options.photopeak_cluster_sigma );
  m_clusterWidth->setDisabled( !options.multiple_nucs_contribute_to_peaks );
  if( m_clusterWidth->label() )
    m_clusterWidth->label()->setDisabled( !m_multiIsoPerPeak->isChecked() );
  m_photopeak_cluster_sigma = options.photopeak_cluster_sigma;
  m_attenForAir->setChecked( options.attenuate_for_air );
  m_backgroundPeakSub->setChecked( options.background_peak_subtract );
  m_sameIsotopesAge->setChecked( options.same_age_isotopes );
  m_decayCorrect->setChecked( options.account_for_decay_during_meas );
  m_showChiOnChart->setChecked( show_chi_on_chart );
  m_chi2Graphic->setShowChiOnChart( show_chi_on_chart );
  m_distanceEdit->setValueText( WString::fromUTF8(dist_node->value()) );
  m_prevDistStr = dist_node->value();
  
  deSerializePeaksToUse( peaks_node );
  m_sourceModel->repopulateIsotopes();
  deSerializeSourcesToFitFor( isotope_nodes );
  deSerializeShieldings( shieldings_node );
  
  m_modifiedThisForeground = true;

  // Make sure if DRF is fixed/not-fixed, that the GUI state will be consistent with this
  //  (it may not have been the same when serialized)
  handleDetectorChanged( m_detectorDisplay->detector() );
  
  updateChi2Chart();
}//void deSerialize( rapidxml::xml_document<char> &doc )


::rapidxml::xml_node<char> *ShieldingSourceDisplay::serialize( rapidxml::xml_node<char> *parent_node )
{
  rapidxml::xml_document<char> *doc = parent_node->document();
  
  const char *name, *value;
  rapidxml::xml_node<> *base_node, *node, *peaks_node, *isotope_nodes, *shieldings_node;
  rapidxml::xml_attribute<> *attr;
  
  name = "ShieldingSourceFit";
  base_node = doc->allocate_node( rapidxml::node_element, name );
  parent_node->append_node( base_node );
  
  //If you change the available options or formatting or whatever, increment the
  //  version field of the XML!
  const string versionstr = std::to_string(ShieldingSourceDisplay::sm_xmlSerializationMajorVersion)
                          + "." + std::to_string(ShieldingSourceDisplay::sm_xmlSerializationMinorVersion);
  value = doc->allocate_string( versionstr.c_str() );
  attr = doc->allocate_attribute( "version", value );
  base_node->append_attribute( attr );
  
  // TODO: maybe add detector names, and sample numbers being used, for both the foeground
  //       and background - this will make it easier to extract the model out of the N42 file.
  //       Or rather, could associate the shield/source model with foreground sample
  //       numbers/detector, like peaks (but then would still need background info).
  
  value = GammaInteractionCalc::to_str( geometry() );
  node = doc->allocate_node( rapidxml::node_element, "Geometry", value );
  base_node->append_node( node );
  
  const ShieldingSourceFitCalc::ShieldingSourceFitOptions options = fitOptions();
  options.serialize( base_node );
  
  name = "ShowChiOnChart";
  value = m_showChiOnChart->isChecked() ? "1" : "0";
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "Distance";
  value = doc->allocate_string( m_distanceEdit->valueText().toUTF8().c_str() );
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );

  name = "Shieldings";
  shieldings_node = doc->allocate_node( rapidxml::node_element, name );
  base_node->append_node( shieldings_node );

  
  const vector<WWidget *> shieldings = m_shieldingSelects->children();
  for( WWidget *child : shieldings )
  {
    const ShieldingSelect *select = dynamic_cast<ShieldingSelect *>( child );
    if( select )
      select->serialize( shieldings_node );
  }//for( WWebWidget *child : shieldings )
  
  peaks_node = doc->allocate_node( rapidxml::node_element, "Peaks" );
  base_node->append_node( peaks_node );
  
  for( int peakn = 0; peakn < m_peakModel->rowCount(); ++peakn )
  {
    const PeakDef &peak = m_peakModel->peak( peakn );
    const SandiaDecay::Nuclide *nuclide = peak.parentNuclide();
    const SandiaDecay::RadParticle *particle = peak.decayParticle();
    const bool annhilation = peak.sourceGammaType()==PeakDef::AnnihilationGamma;
    
    if( peak.useForShieldingSourceFit() && nuclide && (particle || annhilation) )
    {
      rapidxml::xml_node<> *peak_node, *nuc_node, *energy_node;
      peak_node = doc->allocate_node( rapidxml::node_element, "Peak" );
      peaks_node->append_node( peak_node );
      
      //meh, lets be explicit about using this peak, although were only writing
      //  peaks were using for the fit
      attr = doc->allocate_attribute( "Use", "1" );
      peak_node->append_attribute( attr );
      
      value = doc->allocate_string( nuclide->symbol.c_str() );
      nuc_node = doc->allocate_node( rapidxml::node_element, "Nuclide", value );
      peak_node->append_node( nuc_node );
      
      const float energy = peak.gammaParticleEnergy();
      value = doc->allocate_string( std::to_string(energy).c_str() );
      energy_node = doc->allocate_node( rapidxml::node_element, "Energy", value );
      peak_node->append_node( energy_node );
    }//if( peak.useForShieldingSourceFit() )
  }//for( int peakn = 0; peakn < m_peakModel->rowCount(); ++peakn )
  
  isotope_nodes = doc->allocate_node( rapidxml::node_element, "Nuclides" );
  base_node->append_node( isotope_nodes );
  
  const vector<ShieldingSourceFitCalc::IsoFitStruct> &srcData = m_sourceModel->underlyingData();
  assert( srcData.size() == m_sourceModel->rowCount() );
  for( size_t nuc = 0; nuc < srcData.size(); ++nuc )
    srcData[nuc].serialize( isotope_nodes );
  
  
  try
  {
    auto fcnAndPars = shieldingFitnessFcn();
    
    std::shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn> &chi2Fcn = fcnAndPars.first;
    ROOT::Minuit2::MnUserParameters &inputPrams = fcnAndPars.second;
    
    const unsigned int ndof = inputPrams.VariableParameters();
    const vector<double> params = inputPrams.Params();
    const vector<double> errors = inputPrams.Errors();
    GammaInteractionCalc::ShieldingSourceChi2Fcn::NucMixtureCache mixcache;
    const vector<GammaInteractionCalc::PeakResultPlotInfo> chis
                  = chi2Fcn->energy_chi_contributions( params, errors, mixcache, nullptr, nullptr );
    
    if( chis.size() )
    {
      char buffer[64] = { 0 };
      
      //We will write this information - for the record, JIC, but we will not
      //  read it back in later, as we will re-calculate it from the actual
      //  data present when this model gets deserialized.
      rapidxml::xml_node<> *chi2_node = doc->allocate_node( rapidxml::node_element, "Chi2Elements" );
      base_node->append_node( chi2_node );
      
      double chi2 = 0.0;
      rapidxml::xml_node<> *node = nullptr;
      rapidxml::xml_node<> *eval_node = doc->allocate_node( rapidxml::node_element, "EvaluatedEnergies" );
      
      for( const GammaInteractionCalc::PeakResultPlotInfo &p : chis )
      {
        rapidxml::xml_node<> *point_node = doc->allocate_node( rapidxml::node_element, "EvalPoint" );
        eval_node->append_node( point_node );
        
        const double energy = p.energy;
        const double chi = p.numSigmaOff;
        const double scale = p.observedOverExpected;
        const WColor &color = p.peakColor;
        const double scale_uncert = p.observedOverExpectedUncert;
        
        if( !IsInf(chi) && !IsNan(chi) )
          chi2 += chi*chi;
        
        node = doc->allocate_node( rapidxml::node_element, "Energy" );
        snprintf( buffer, sizeof(buffer), "%f", energy );
        value = doc->allocate_string( buffer );
        node->value( value );
        point_node->append_node( node );
        
        node = doc->allocate_node( rapidxml::node_element, "Chi" );
        snprintf( buffer, sizeof(buffer), "%f", chi );
        value = doc->allocate_string( buffer );
        node->value( value );
        point_node->append_node( node );
        
        node = doc->allocate_node( rapidxml::node_element, "Scale" );
        snprintf( buffer, sizeof(buffer), "%f", scale );
        value = doc->allocate_string( buffer );
        node->value( value );
        point_node->append_node( node );
        
        node = doc->allocate_node( rapidxml::node_element, "Color" );
        snprintf( buffer, sizeof(buffer), "%s", (color.isDefault() ? "" : color.cssText().c_str()) );
        value = doc->allocate_string( buffer );
        node->value( value );
        point_node->append_node( node );
        
        node = doc->allocate_node( rapidxml::node_element, "ScaleUncert" );
        snprintf( buffer, sizeof(buffer), "%f", scale_uncert );
        value = doc->allocate_string( buffer );
        node->value( value );
        point_node->append_node( node );
      }//for( const auto &p : chis )
      
      node = doc->allocate_node( rapidxml::node_element, "Chi2" );
      snprintf( buffer, sizeof(buffer), "%f", chi2 );
      value = doc->allocate_string( buffer );
      node->value( value );
      chi2_node->append_node( node );
      
      node = doc->allocate_node( rapidxml::node_element, "NumParamFit" );
      snprintf( buffer, sizeof(buffer), "%u", ndof );
      value = doc->allocate_string( buffer );
      node->value( value );
      chi2_node->append_node( node );
      
      chi2_node->append_node( eval_node );
    }//if( chis.size() )
  }catch( std::exception &e )
  {
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, ("Failed to get chi2 info during serialization - caught exception: " + string(e.what())).c_str() );
#endif
  }
  
  return parent_node;
}//::rapidxml::xml_node<char> * serialize()



void ShieldingSourceDisplay::layoutSizeChanged( int width, int height )
{
    //If phone, we just force to show since we don't care about the screen size
    if (m_specViewer->isPhone())
      return;
    
  m_width = width;
  m_height = height;

  if( m_nResizeSinceHint > 3 )
  {
    if( m_height < 320 && !m_chi2Graphic->isHidden() )
    {
      m_chi2Graphic->hide();
//      m_layout->setRowStretch( 1, 0 );
      m_showChi2Text->show();
    }else if( m_chi2Graphic->isHidden() && m_height >= 320 )
    {
      m_chi2Graphic->show();
      m_showChi2Text->hide();
//      m_layout->setRowStretch( 1, 6 );
//      leftLayout->setRowResizable( 0, true );
//      leftLayout->setRowResizable( 1, true );
    }//if / else
  }//if( m_nResizeSinceHint > 3 )

  ++m_nResizeSinceHint;
}//void layoutSizeChanged( int width, int height )


void ShieldingSourceDisplay::initialSizeHint( int width, int height )
{
  m_nResizeSinceHint = 5;
  layoutSizeChanged( width, height );
  m_nResizeSinceHint = 0;
}//void initialSizeHint( int width, int height )


void ShieldingSourceDisplay::doAddShielding()
{
  addShielding( nullptr, true );
}

void ShieldingSourceDisplay::doAddShieldingBefore( ShieldingSelect *select )
{
  addShielding( select, true );
}


void ShieldingSourceDisplay::doAddShieldingAfter( ShieldingSelect *select )
{
  const vector<WWidget *> selects = m_shieldingSelects->children();
  
  int index = -1;
  for( size_t i = 0; i < selects.size(); ++i )
  {
    ShieldingSelect *s = dynamic_cast<ShieldingSelect *>( selects[i] );
      
    if( s == select )
      index = static_cast<int>( i );
    
    if( s && index>=0 && (i < (selects.size()-1)) )
    {
      ShieldingSelect *next = dynamic_cast<ShieldingSelect *>( selects[i+1] );
      addShielding( next, true );
      return;
    }
  }//for( size_t i = 0; i < selects.size(); ++i )
  
  if( index < 0 )
    cerr << "ShieldingSourceDisplay::doAddShieldingAfter(): Warning, failed to"
    << " find passed in select" << endl;
  
  //if index is >= 0, then user requested to add a shielding after the last
  //  shielding
  addShielding( nullptr, true );
}//void doAddShieldingAfter( ShieldingSelect *select )


size_t ShieldingSourceDisplay::numberShieldings() const
{
  // A simple count of children widgets should be good enough, but maybe to just be thorough,
  //  we'll do a dynamic cast to make sure each child is a ShieldingSelect
  //return m_shieldingSelects->children().size();
  
  size_t nwidgets = 0;
  for( WWidget *w : m_shieldingSelects->children() )
    nwidgets += (dynamic_cast<ShieldingSelect *>(w) != nullptr);
    
  assert( nwidgets == m_shieldingSelects->children().size() );
  
  return nwidgets;
}//size_t numberShieldings() const;


void ShieldingSourceDisplay::addGenericShielding()
{
  ShieldSourceChange state_undo_creator( this, "Add generic shielding" );
  
  ShieldingSelect *temp = addShielding( nullptr, false );
  if( !temp->isGenericMaterial() )
    temp->handleToggleGeneric();
}//void addGenericShielding()



ShieldingSelect *ShieldingSourceDisplay::addShielding( ShieldingSelect *before,
                                                       const bool addUndoRedo )
{
  // Handle undo/redo if the user clicked a button to add a shielding (all paths to here
  //  will have `addUndoRedo == true` in this case, and all paths to here
  //  not from a user explicitly clicking a button, will have this as false).
  unique_ptr<ShieldSourceChange> state_undo_creator;
  if( addUndoRedo )
    state_undo_creator = make_unique<ShieldSourceChange>( this, "Add Shielding" );
  
  m_modifiedThisForeground = true;
  
  ShieldingSelect *select = new ShieldingSelect( m_materialDB, m_sourceModel, m_materialSuggest, this );

  if( before && m_shieldingSelects->indexOf(before) >= 0 )
    m_shieldingSelects->insertBefore( select, before );
  else
    m_shieldingSelects->addWidget( select );
  
  select->setGeometry( geometry() );
  
  const auto det = m_detectorDisplay->detector();
  select->setFixedGeometry( det && det->isFixedGeometry() );
  
  
  select->addShieldingBefore().connect( boost::bind( &ShieldingSourceDisplay::doAddShieldingBefore,
                                                    this, boost::placeholders::_1 ) );
  select->addShieldingAfter().connect( boost::bind( &ShieldingSourceDisplay::doAddShieldingAfter,
                                                   this, boost::placeholders::_1 ) );
  
  //connect up signals of select and such
  select->m_arealDensityEdit->valueChanged().connect( this, &ShieldingSourceDisplay::updateChi2Chart );
  select->m_atomicNumberEdit->valueChanged().connect( this, &ShieldingSourceDisplay::updateChi2Chart );

  select->materialChanged().connect( this, &ShieldingSourceDisplay::handleShieldingChange );
  select->materialModified().connect( this, &ShieldingSourceDisplay::handleShieldingChange );
  
  select->remove().connect( this, &ShieldingSourceDisplay::removeShielding );
  select->materialModified().connect( this, &ShieldingSourceDisplay::materialModifiedCallback );
  select->materialChanged().connect( this, &ShieldingSourceDisplay::materialChangedCallback );

  select->addingIsotopeAsSource().connect( boost::bind(
      &ShieldingSourceDisplay::isotopeIsBecomingVolumetricSourceCallback, this,
      select, boost::placeholders::_1, boost::placeholders::_2 ) );
  select->removingIsotopeAsSource().connect( boost::bind(
      &ShieldingSourceDisplay::isotopeRemovedAsVolumetricSourceCallback, this,
      select, boost::placeholders::_1, boost::placeholders::_2 ) );
  select->activityFromVolumeNeedUpdating().connect( boost::bind(
      &ShieldingSourceDisplay::updateActivityOfShieldingIsotope, this,
      boost::placeholders::_1, boost::placeholders::_2 ) );

  
  Signal<ShieldingSelect *, shared_ptr<const string>, shared_ptr<const string>> &
      undoSignal = select->userChangedStateSignal();
  undoSignal.connect( boost::bind(
              &ShieldingSourceDisplay::handleShieldingUndoRedoPoint, this,
              boost::placeholders::_1, boost::placeholders::_2, boost::placeholders::_3 ) );
  
  
  handleShieldingChange();
  select->setTraceSourceBtnStatus();
  
  updateChi2Chart();
  
  return select;
}//void addShielding()



void ShieldingSourceDisplay::addSourceIsotopesToShieldings( Wt::WModelIndex, int firstRow, int lastRow )
{
//  m_modifiedThisForeground = true;
  
  firstRow = max( 0, firstRow );
  lastRow = min( lastRow, m_sourceModel->rowCount() );

  vector<const SandiaDecay::Nuclide *> isotopes;
  for( int row = firstRow; row <= lastRow; ++row )
  {
    const SandiaDecay::Nuclide *nuc = m_sourceModel->nuclide( row );
    if( nuc )
      isotopes.push_back( nuc );
  }//for( int row = firstRow; row <= lastRow; ++row )

  if( isotopes.empty() )
    return;

  for( WWidget *widget : m_shieldingSelects->children() )
  {
    ShieldingSelect *thisSelect = dynamic_cast<ShieldingSelect *>(widget);
    if( !thisSelect || thisSelect->isGenericMaterial() )
      continue;
    
    for( const SandiaDecay::Nuclide *iso : isotopes )
      thisSelect->modelNuclideAdded( iso );
  }//for( WWidget *widget : m_shieldingSelects->children() )
}//void addSourceIsotopesToShieldings( Wt::WModelIndex, int firstRow, int lastRow )


void ShieldingSourceDisplay::removeSourceIsotopesFromShieldings( Wt::WModelIndex, int firstRow, int lastRow )
{
//  m_modifiedThisForeground = true;
  
  const int nrow = m_sourceModel->rowCount();
  if( firstRow<0 || firstRow>lastRow || lastRow<0 || lastRow>=nrow )
  {
    stringstream msg;
    msg << "ShieldingSourceDisplay::removeSourceIsotopesFromShieldings(...)\n\tfunction called with invalid argument "
        << "- this is a major logic error I think";
    throw std::runtime_error( msg.str() );
  }//if( invalid row )

  vector<const SandiaDecay::Nuclide *> isotopes;
  for( int row = firstRow; row <= lastRow; ++row )
  {
    const SandiaDecay::Nuclide *nuc = m_sourceModel->nuclide( row );
    if( nuc )
      isotopes.push_back( nuc );
  }//for( int row = firstRow; row <= lastRow; ++row )


  for( WWidget *widget : m_shieldingSelects->children() )
  {
    ShieldingSelect *thisSelect = dynamic_cast<ShieldingSelect *>(widget);
    if( !thisSelect || thisSelect->isGenericMaterial() )
      continue;

    for( const SandiaDecay::Nuclide *iso : isotopes )
      thisSelect->sourceRemovedFromModel( iso );
  }//for( WWidget *widget : m_shieldingSelects->children() )
}//void removeSourceIsotopesFromShieldings( Wt::WModelIndex, int firstRow, int lastRow )


void ShieldingSourceDisplay::materialModifiedCallback( ShieldingSelect *select )
{
  if( !select )
  {
    cerr << "ShieldingSourceDisplay::materialModifiedCallback(...)\n\tShouldnt be here!" << endl;
    return;
  }//if( !select )

  
  m_modifiedThisForeground = true;
  
  //I meant to do some more work here, but dont recall what

  updateChi2Chart();
}//void materialModifiedCallback( ShieldingSelect *select )



void ShieldingSourceDisplay::materialChangedCallback( ShieldingSelect *select )
{
  if( !select )
  {
    cerr << "ShieldingSourceDisplay::materialChangedCallback(...)\n\tShouldnt be here!" << endl;
    return;
  }//if( !select )
  
  //The select has already removed any isotopes as sources which arent in the
  //  current material, however, we have to add in all other candidate isotopes

  m_modifiedThisForeground = true;
  
  const int nrow = m_sourceModel->rowCount();
  for( int row = 0; row < nrow; ++row )
  {
    const SandiaDecay::Nuclide *nuc = m_sourceModel->nuclide( row );
    if( nuc )
      select->modelNuclideAdded( nuc );
  }//for( int row = 0; row < nrow; ++row )

  updateChi2Chart();
}//void materialChangedCallback( ShieldingSelect *select )



void ShieldingSourceDisplay::updateActivityOfShieldingIsotope( ShieldingSelect *select,
                                       const SandiaDecay::Nuclide *nuc )
{
  // This function gets called when shielding changes for either self-attenuating or trace sources.
  if( !nuc ) //This can happen if fitting for "other" non-source nuclides fraction.
    return;
  
  shared_ptr<const Material> material = select ? select->material() : nullptr;
  
  if( !select || select->isGenericMaterial() || !material || !nuc )
  {
    assert(0);
    cerr << "updateActivityOfShieldingIsotope()\n\tShould not be here!" << endl;
    return;
  }//if( !select || select->isGenericMaterial() )

  const shared_ptr<const DetectorPeakResponse> det = m_detectorDisplay->detector();
  if( det && det->isFixedGeometry() )
  {
    // This happens during deSerialization, if we serialized with a non-fixed geometry DRF
    //  but then are deSerializing with a fixed geometry DRF
    return;
  }
  
  char buffer[64] = { '\0' };
  const int row = m_sourceModel->row( nuc );
  const WModelIndex index = m_sourceModel->index( row, SourceFitModel::kActivity );
  
  if( select->isTraceSourceForNuclide(nuc) )
  {
    // JIC there were any geometry changes, force the ShieldingSelect to update its total activity
    const double activity = select->updateTotalTraceSourceActivityForGeometryChange(nuc);
    
    // Setting the activity to the model will clear the uncertatnty
    //   If user is typing in activty to the widget, we probably want this.
    //   But if this path is some autoamted update, then this would be bad (and hopefully we're avoiding)
    snprintf( buffer, sizeof(buffer), "%1.8E ci", activity/PhysicalUnits::curie );
    m_sourceModel->setData( index, string(buffer) );
    
    return;
  }//if( select->isTraceSourceForNuclide(nuc) )
  
  const double density = material->density;
  const double volume = select->shieldingVolume();
  
  bool foundNucInMaterial = false;
  double weight = density * volume;
  
  const map<const SandiaDecay::Element *,vector<ShieldingSelect::NucMasFrac>> srcNucFracs
                                               = select->sourceNuclideMassFractions();
  for( const auto &el_nucs : srcNucFracs )
  {
    for( const ShieldingSelect::NucMasFrac &nuc_frac : el_nucs.second )
    {
      if( get<0>(nuc_frac) == nuc )
      {
        foundNucInMaterial = true;
        weight *= get<1>(nuc_frac);
        break;
      }
    }//for( const auto &nuc_frac : srcNucFracs )
  }//for( const auto &el_nucs : srcNucFracs )
  
  if( !foundNucInMaterial )
  {
    const string msg = "ShieldingSourceDisplay::updateActivityOfShieldingIsotope(...)\n"
                       "\tShould not be here! nuc=" + nuc->symbol
                       + " and material=" + material->name;
    throw std::runtime_error( msg );
  }//if( !nuclide )

  const double weight_in_grams = weight / PhysicalUnits::gram;
  const double activity_per_gram = nuc->activityPerGram();
  const double activity = activity_per_gram * weight_in_grams / PhysicalUnits::curie;

//  cerr << "updateActivityOfShieldingIsotope: weight_in_grams=" << weight_in_grams
//  << ", activity_per_gram=" << activity_per_gram
//  << ", activity=" << activity << " ci" << endl;
  
  if( row < 0 )
  {
    const char *msg = "ShieldingSourceDisplay::updateActivityOfShieldingIsotope(...)\n\tShould not be here!";
    throw std::runtime_error( msg );
  }//if( row < 0 )

  
  if( IsNan(activity) || IsInf(activity) )
  {
    string msg = "An invalid activity was calculated for " + nuc->symbol
                 + ", other results may be invalid";
    passMessage( msg, WarningWidget::WarningMsgHigh );
    return;
  }//if( IsNan(p.second) || IsInf(p.second) )
  
  
  snprintf( buffer, sizeof(buffer), "%1.8E ci", activity );
  m_sourceModel->setData( index, string(buffer) );
}//void updateActivityOfShieldingIsotope(..)


void ShieldingSourceDisplay::isotopeIsBecomingVolumetricSourceCallback(
                                                ShieldingSelect *caller,
                                                const SandiaDecay::Nuclide *nuc,
                                                const ShieldingSourceFitCalc::ModelSourceType type )
{
  assert( nuc );
  assert( caller );
  if( !caller )
    return;
  
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  
  const vector<WWidget *> &children = m_shieldingSelects->children();
  
  if( !caller->material() )
  {
    
  }//if( !caller->material() )
  
  //Make sure no other selects have this isotope selected
  for( WWidget *widget : children )
  {
    ShieldingSelect *select = dynamic_cast<ShieldingSelect *>( widget );
    if( !select || (select==caller) )
      continue;

    select->uncheckSourceIsotopeCheckBox( nuc );
  }//for( WWidget *widget : children )

  //Set appropriate flags in the SourceFitModel so activity wont be editiable
  m_sourceModel->setSourceType( nuc, type );
  
  switch( type )
  {
    case ShieldingSourceFitCalc::ModelSourceType::Point:
      break;
      
    case ShieldingSourceFitCalc::ModelSourceType::Intrinsic:
      break;
      
    case ShieldingSourceFitCalc::ModelSourceType::Trace:
    {
      // Sync the model up so fitting activity is same as trace-source widget
      const bool fitAct = caller->fitTraceSourceActivity(nuc);
      WModelIndex fitActIndex = m_sourceModel->index( nuc, SourceFitModel::kFitActivity );
      m_sourceModel->setData( fitActIndex, fitAct, Wt::CheckStateRole );
      break;
    }
  }//switch( type )
  
  // Reset the age to something reasonable.
  //  TODO: do we really want to do this?  Maybe not.
  string agestr = "";
  if( nuc )  //nuc should always be non-NULL, but just incase
  {
    PeakDef::defaultDecayTime( nuc, &agestr );
    WModelIndex index = m_sourceModel->index( nuc, SourceFitModel::kAge );
    m_sourceModel->setData( index, agestr );
  }//if( nuc )

  WModelIndex fitAgeIndex = m_sourceModel->index( nuc, SourceFitModel::kFitAge );
  m_sourceModel->setData( fitAgeIndex, false );

  
  // Update if other shieldings can have add trace source menu item enabled
  for( WWidget *widget : children )
  {
    ShieldingSelect *select = dynamic_cast<ShieldingSelect *>( widget );
    if( select )
      select->setTraceSourceBtnStatus();
  }//for( WWidget *widget : children )
  
  
  //Update activity displayed
  updateActivityOfShieldingIsotope( caller, nuc );

  //Update the Chi2
  updateChi2Chart();
}//isotopeIsBecomingVolumetricSourceCallback(...)


void ShieldingSourceDisplay::isotopeRemovedAsVolumetricSourceCallback(
                                            ShieldingSelect *select,
                                            const SandiaDecay::Nuclide *nuc,
                                            const ShieldingSourceFitCalc::ModelSourceType type )
{
  //Set appropriate flags in the SourceFitModel so activity will be editiable
  
  assert( !nuc || (m_sourceModel->sourceType(m_sourceModel->nuclideIndex(nuc)) == type) );

  m_sourceModel->setSourceType( nuc, ShieldingSourceFitCalc::ModelSourceType::Point );
  
  const vector<WWidget *> &children = m_shieldingSelects->children();
  for( WWidget *widget : children )
  {
    ShieldingSelect *select = dynamic_cast<ShieldingSelect *>( widget );
    if( select )
      select->setTraceSourceBtnStatus();
  }//for( WWidget *widget : children )
  
  //Update the Chi2
  updateChi2Chart();
}//isotopeRemovedAsVolumetricSourceCallback(...)


void ShieldingSourceDisplay::handleShieldingUndoRedoPoint( const ShieldingSelect * const select,
                                  const shared_ptr<const string> &prev_state,
                                  const shared_ptr<const string> &current_state )
{
  if( !select )
    return;
  
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( !undoRedo || !undoRedo->canAddUndoRedoNow() )
    return;
  
  //Find index of the select
  const vector<WWidget *> &kids = m_shieldingSelects->children();
  size_t index = 0;
  for( ; index < kids.size(); ++index )
  {
    if( select == dynamic_cast<ShieldingSelect *>( kids[index] ) )
      break;
  }
  
  if( index >= kids.size() )
  {
    Wt::log("error") << "ShieldingSourceDisplay::handleShieldingUndoRedoPoint: couldnt find passed in ShieldingSelect.";
    return;
  }
  
  const shared_ptr<const DetectorPeakResponse> det = m_detectorDisplay->detector();
  const bool is_fixed_geom = det && det->isFixedGeometry();
  
  auto undo_redo = [prev_state, current_state, index, is_fixed_geom]( const bool is_undo ){
    
    // TODO: We arent totally getting things right for trace and self-attenuating sources.
    // TODO: I think we are getting multiple undo/redo steps for a single user action for some things related to trace and self-attenuating sources.
    
    try
    {
      const auto xml_state = is_undo ? prev_state : current_state;
      
      ShieldingSourceDisplay *display = InterSpec::instance()->shieldingSourceFit();
      if( !xml_state || xml_state->empty() || !display || !display->m_shieldingSelects )
      {
        Wt::log("error") << "ShieldingSourceDisplay undo/redo - invalid step???";
        return;
      }
      
      const vector<WWidget *> &kids = display->m_shieldingSelects->children();
      
      ShieldingSelect * const select = (index < kids.size())
      ? dynamic_cast<ShieldingSelect *>( kids[index] ) : nullptr;
      
      if( !select )
      {
        Wt::log("error") << "ShieldingSourceDisplay undo/redo - not enough shieldings???";
        return;
      }
      
      const vector<const SandiaDecay::Nuclide *> pre_trace_nucs = select->traceSourceNuclides();
      const vector<const SandiaDecay::Nuclide *> pre_self_att_nucs = select->selfAttenNuclides();
      
      string xml_state_cpy = *xml_state;
      
      rapidxml::xml_document<char> doc;
      const int flags = rapidxml::parse_normalize_whitespace | rapidxml::parse_trim_whitespace;
      doc.parse<flags>( &(xml_state_cpy[0]) );
      select->deSerialize( doc.first_node(), is_fixed_geom );
      
      const vector<const SandiaDecay::Nuclide *> post_trace_nucs = select->traceSourceNuclides();
      const vector<const SandiaDecay::Nuclide *> post_self_att_nucs = select->selfAttenNuclides();
      
      assert( !is_fixed_geom || post_trace_nucs.empty() );
      assert( !is_fixed_geom || post_self_att_nucs.empty() );
      
      for( auto *pre_nuc : pre_trace_nucs )
      {
        if( std::find(begin(post_trace_nucs), end(post_trace_nucs), pre_nuc) == end(post_trace_nucs) )
        {
          // Somehow maybe things have already been in display->m_sourceModel...
          const auto nuc_index = display->m_sourceModel->nuclideIndex(pre_nuc);
          const ShieldingSourceFitCalc::ModelSourceType type = display->m_sourceModel->sourceType( nuc_index );
          if( type != ShieldingSourceFitCalc::ModelSourceType::Point )
            display->isotopeRemovedAsVolumetricSourceCallback( select, pre_nuc, ShieldingSourceFitCalc::ModelSourceType::Trace );
        }
      }//for( auto *pre_nuc : pre_trace_nucs )
      
      for( auto *post_nuc : post_trace_nucs )
      {
        if( std::find(begin(pre_trace_nucs), end(pre_trace_nucs), post_nuc) == end(pre_trace_nucs) )
          display->isotopeIsBecomingVolumetricSourceCallback( select, post_nuc, ShieldingSourceFitCalc::ModelSourceType::Trace );
      }
      
      for( auto *pre_nuc : pre_self_att_nucs )
      {
        if( std::find(begin(post_self_att_nucs), end(post_self_att_nucs), pre_nuc) == end(post_self_att_nucs) )
        {
          // Somehow maybe things have already been in display->m_sourceModel...
          const auto nuc_index = display->m_sourceModel->nuclideIndex(pre_nuc);
          const ShieldingSourceFitCalc::ModelSourceType type = display->m_sourceModel->sourceType( nuc_index );
          if( type != ShieldingSourceFitCalc::ModelSourceType::Point )
            display->isotopeRemovedAsVolumetricSourceCallback( select, pre_nuc, ShieldingSourceFitCalc::ModelSourceType::Intrinsic );
        }
      }//for( auto *pre_nuc : pre_self_att_nucs )
      
      for( auto *post_nuc : post_self_att_nucs )
      {
        if( std::find(begin(pre_self_att_nucs), end(pre_self_att_nucs), post_nuc) == end(pre_self_att_nucs) )
          display->isotopeIsBecomingVolumetricSourceCallback( select, post_nuc, ShieldingSourceFitCalc::ModelSourceType::Intrinsic );
      }

      display->handleShieldingChange();
    }catch( std::exception &e )
    {
      Wt::log("error") << "Error executing ShieldingSelect undo/redo: " << e.what();
    }//try / catch
  };//auto undo_redo
  
  auto undo = [undo_redo](){ undo_redo(true); };
  auto redo = [undo_redo](){ undo_redo(false); };
  
  undoRedo->addUndoRedoStep( undo, redo, "Update shielding" );
}//void handleShieldingUndoRedoPoint(...)


void ShieldingSourceDisplay::removeShielding( ShieldingSelect *select )
{
  if( !select )
    return;

  ShieldSourceChange state_undo_creator( this, "Remove Shielding" );
  
  m_modifiedThisForeground = true;
  
  bool foundShielding = false;
  const vector<WWidget *> &children = m_shieldingSelects->children();
  for( WWidget *widget : children )
  {
    ShieldingSelect *shielding = dynamic_cast<ShieldingSelect *>( widget );
    if( shielding == select )
    {
      delete shielding;
      handleShieldingChange();
      foundShielding = true;
      break;
    }//if( shielding == select )
  }//for( WWidget *widget : children )

  // Now go through and update if we can add a trace source to the other shieldings
  for( WWidget *widget : children )
  {
    ShieldingSelect *shielding = dynamic_cast<ShieldingSelect *>( widget );
    if( shielding )
      shielding->setTraceSourceBtnStatus();
  }//for( WWidget *widget : children )
  
  assert( foundShielding );
  if( !foundShielding )
    cerr << "\n\nCouldnt finding select to delete" << endl;
}//void removeShielding( ShieldingSelect *select )








void ShieldingSourceDisplay::setWidgetStateForFitStarting()
{
  m_fitModelButton->hide();
  m_fitProgressTxt->show();
  m_fitProgressTxt->setText("");
  m_cancelfitModelButton->show();
  
  m_peakView->disable();
  m_sourceView->disable();
  m_optionsDiv->disable();
  m_addItemMenu->disable();
  m_distanceEdit->disable();
  m_detectorDisplay->disable();
  m_shieldingSelects->disable();
  m_addGenericShielding->disable();
  m_addMaterialShielding->disable();
#if( USE_DB_TO_STORE_SPECTRA )
  m_saveAsNewModelInDb->disable();
#endif
}//void setWidgetStateForFitStarting()


void ShieldingSourceDisplay::setWidgetStateForFitBeingDone()
{
  m_fitModelButton->show();
  m_fitProgressTxt->setText("");
  m_fitProgressTxt->hide();
  m_cancelfitModelButton->hide();
  
  m_peakView->enable();
  m_sourceView->enable();
  m_optionsDiv->enable();
  m_addItemMenu->enable();
  m_distanceEdit->enable();
  m_detectorDisplay->enable();
  m_shieldingSelects->enable();
  m_addGenericShielding->enable();
  m_addMaterialShielding->enable();
#if( USE_DB_TO_STORE_SPECTRA )
  m_saveAsNewModelInDb->enable();
#endif
}//void setWidgetStateForFitBeingDone();



void ShieldingSourceDisplay::cancelModelFit()
{
  std::lock_guard<std::mutex> lock( m_currentFitFcnMutex );
  if( m_currentFitFcn )
    m_currentFitFcn->cancelFit();
}//void cancelModelFit()


void ShieldingSourceDisplay::cancelModelFitWithNoUpdate()
{
  std::lock_guard<std::mutex> lock( m_currentFitFcnMutex );
  if( m_currentFitFcn )
  {
    m_currentFitFcn->cancelFit();
    m_currentFitFcn = nullptr;
    
    setWidgetStateForFitBeingDone();
  }//if( m_currentFitFcn )
}//void cancelModelFitWithNoUpdate()


void ShieldingSourceDisplay::updateGuiWithModelFitProgress( std::shared_ptr<ShieldingSourceFitCalc::ModelFitProgress> progress )
{
  ShieldingSourceFitCalc::ModelFitProgress status;
  
  {
    std::lock_guard<std::mutex> lock( progress->m_mutex );
    status.chi2 = progress->chi2;
    status.numFcnCalls = progress->numFcnCalls;
    status.elapsedTime = progress->elapsedTime;
    status.parameters = progress->parameters;
  }

  status.chi2 = floor(10.0*status.chi2 + 0.5) / 10.0;
  status.elapsedTime = floor(10.0*status.elapsedTime + 0.5) / 10.0;
  
  char buffer[128];
  snprintf( buffer, sizeof(buffer),
            "%s &chi;&sup2;: %.1f, %i fcn calls in %.1fs",
            WString::tr("ssd-Best").toUTF8().c_str(),
            status.chi2, static_cast<int>(status.numFcnCalls), status.elapsedTime );
  
  m_fitProgressTxt->setText( buffer );
  
  wApp->triggerUpdate();
}//void updateGuiWithModelFitProgress( std::shared_ptr<ModelFitProgress> progress )


void ShieldingSourceDisplay::updateGuiWithModelFitResults( std::shared_ptr<ShieldingSourceFitCalc::ModelFitResults> results )
{
  WApplication *app = wApp;
  WApplication::UpdateLock applock( app );
  
  if( !applock )
  {
    // Shouldnt ever get here
    cerr << "Failed to get application lock!" << endl;
    return;
  }
  
  // Make sure we trigger a app update
  BOOST_SCOPE_EXIT(app){
    if( app )
      app->triggerUpdate();
  } BOOST_SCOPE_EXIT_END
 
  ShieldSourceChange state_undo_creator( this, "Fit activity/shielding" );
  
  assert( results );
  
  setWidgetStateForFitBeingDone();
  
  std::lock( results->m_mutex, m_currentFitFcnMutex );
  std::lock_guard<std::mutex> result_lock( results->m_mutex, std::adopt_lock );
  std::lock_guard<std::mutex> lock( m_currentFitFcnMutex, std::adopt_lock );
  
  const ShieldingSourceFitCalc::ModelFitResults::FitStatus status = results->successful;
  const vector<ShieldingSourceFitCalc::ShieldingInfo> &initial_shieldings = results->initial_shieldings;
  const vector<ShieldingSourceFitCalc::FitShieldingInfo> &final_shieldings = results->final_shieldings;
  const vector<double> &paramValues = results->paramValues;
  const vector<double> &paramErrors = results->paramErrors;
  const vector<string> &errormsgs = results->errormsgs;
  
  assert( (status != ShieldingSourceFitCalc::ModelFitResults::FitStatus::Final)
         || (initial_shieldings.size() == final_shieldings.size()) );
  
  if( !m_currentFitFcn )
  {
    passMessage( "Programming Logic Error - received model fit results at an invalid time.", WarningWidget::WarningMsgHigh );
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, "Programming Logic Error - received model fit results at an invalid time." );
#endif
    return;
  }//if( !m_currentFitFcn )
  
  vector<ShieldingSelect *> gui_shieldings;
  for( WWidget *widget : m_shieldingSelects->children() )
  {
    ShieldingSelect *select = dynamic_cast<ShieldingSelect *>(widget);
    if( select )
      gui_shieldings.push_back( select );
  }//for( WWidget *widget : m_shieldingSelects->children() )
  
  
  if( initial_shieldings.size() != gui_shieldings.size() )
  {
    passMessage( "Programming Logic Error - number of shieldings have changed from "
                + std::to_string(initial_shieldings.size())
                + " to "
                + std::to_string(gui_shieldings.size())
                , WarningWidget::WarningMsgHigh );
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, "Programming Logic Error - number of shieldings changed - fit results when model was no longer valid." );
#endif
    m_currentFitFcn.reset();
    return;
  }//if( initial_shieldings.size() != gui_shieldings.size() )
  
#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
  for( size_t i = 0; i < initial_shieldings.size(); ++i )
  {
    const ShieldingSourceFitCalc::ShieldingInfo now_info = gui_shieldings[i]->toShieldingInfo();
    const ShieldingSourceFitCalc::ShieldingInfo orig_info = initial_shieldings[i];
    
    try
    {
      ShieldingSourceFitCalc::ShieldingInfo::equalEnough( now_info, orig_info );
    }catch( std::exception &e )
    {
      cerr << "shieldings changed: " << e.what() << endl;
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, "Programming Logic Error - shieldings changed - fit results when model was no longer valid." );
#endif
      assert( 0 );
    }//try / catch
  }//for( size_t i = 0; i < initial_shieldings.size(); ++i )
#endif
  
  
  switch( status )
  {
    case ShieldingSourceFitCalc::ModelFitResults::FitStatus::TimedOut:
      // TODO: check if best fit Chi2 is better than current, and if so, use new values, if available.
      //       Note that the chi2 on the chart is sqrt(results->chi2) / number_peaks
      
    case ShieldingSourceFitCalc::ModelFitResults::FitStatus::UserCancelled:
    case ShieldingSourceFitCalc::ModelFitResults::FitStatus::InvalidOther:
    {
      string msg = "<b>" + WString::tr("ssd-fit-failed").toUTF8() + "</b>.";
      for( auto &s : errormsgs )
        msg += "<div>&nbsp;&nbsp;" + s + "</div>";
      
      passMessage( msg, WarningWidget::WarningMsgHigh );
      
      m_currentFitFcn.reset();
      return;
    }//
     
    case ShieldingSourceFitCalc::ModelFitResults::FitStatus::InterMediate:
    {
      passMessage( "Intermediate Fit status not handled yet.", WarningWidget::WarningMsgHigh );
      return;
    }//case ModelFitResults::FitStatus::InterMediate:
      
    case ShieldingSourceFitCalc::ModelFitResults::FitStatus::Final:
      break;
  }//switch( status )
  
  
  for( auto &s : errormsgs )
    passMessage( s + "<br />" + WString::tr("ssd-using-fit-anyway").toUTF8(), WarningWidget::WarningMsgHigh );
  
  try
  {
    const size_t nshieldings = gui_shieldings.size();
    assert( gui_shieldings.size() == final_shieldings.size() );
    if( gui_shieldings.size() != final_shieldings.size() )
      throw logic_error( "Number of shieldings changed during fitting - should not happen." );
    
    // First we'll update mass-fractions of self-attenuating sources, if we were fitting any of them
    const vector<const Material *> massfracFitMaterials
                                           = m_currentFitFcn->materialsFittingMassFracsFor();
    
    for( size_t shielding_index = 0; shielding_index < nshieldings; ++shielding_index )
    {
      ShieldingSelect *select = gui_shieldings[shielding_index];
      assert( select );
      
      if( select->isGenericMaterial() )
        continue;
      
      shared_ptr<const Material> usrmaterial = select->material();
      assert( usrmaterial );
      if( !usrmaterial )
        continue;
      
      const bool calcFitMassFrac = std::count(begin(massfracFitMaterials), end(massfracFitMaterials), usrmaterial.get());
      if( calcFitMassFrac != select->fitForAnyMassFractions() )
      {
        throw logic_error( "GUI fit mass fraction for material '" + usrmaterial->name
                                 + "' doesn't match calculation fit mass fraction status." );
      }
      
      if( !calcFitMassFrac )
        continue;
      
      
      const map<const SandiaDecay::Element *,vector<const SandiaDecay::Nuclide *>> fit_el_to_nucs
                                    = m_currentFitFcn->nuclideFittingMassFracFor( shielding_index );
                                  
      const vector<const SandiaDecay::Element *> fitOtherEls
                                = m_currentFitFcn->elementsFittingOtherFracFor( shielding_index );
      
      vector<const SandiaDecay::Nuclide *> guiNucsFittingFrac;
      vector<const SandiaDecay::Element *> guiElFittingOtherFrac;
      
      const map<const SandiaDecay::Element *,vector<ShieldingSelect::NucMasFrac>> all_gui_self_atten
                                                             = select->sourceNuclideMassFractions();
      for( const auto &el_nucs : all_gui_self_atten )
      {
        const SandiaDecay::Element * const this_el = el_nucs.first;
        assert( this_el );
        if( !this_el )
          continue;
        
        for( const ShieldingSelect::NucMasFrac &nuc_mass_frac : el_nucs.second )
        {
          const SandiaDecay::Nuclide * const this_nuc = get<0>(nuc_mass_frac);
          const double this_frac = get<1>(nuc_mass_frac);
          const bool this_fit = get<2>(nuc_mass_frac);
          
          if( this_fit )
          {
            if( this_nuc )
              guiNucsFittingFrac.push_back( this_nuc );
            else
              guiElFittingOtherFrac.push_back( this_el );
          }//if( this_fit )
        }//for( const ShieldingSelect::NucMasFrac &nuc_mass_frac : all_gui_self_atten )
      }//for( const auto &el_nucs : all_gui_self_atten )
      
      vector<const SandiaDecay::Nuclide *> fitnucs; //Will not include "other" non-srce components
      // Since not all self-attenuating nuclides will have mass-fractions fit for, we will remove
      //  the ones we are fitting for from this next variable, to get down to what we want.
      vector<const SandiaDecay::Nuclide *> selfAttenNotFitNucs
                                      = m_currentFitFcn->selfAttenuatingNuclides( shielding_index );
      for( const auto &el_nucs : fit_el_to_nucs )
      {
        for( const SandiaDecay::Nuclide * nuc : el_nucs.second )
        {
          if( nuc )
            fitnucs.push_back( nuc );
          const auto new_end = std::remove(begin(selfAttenNotFitNucs), end(selfAttenNotFitNucs), nuc);
          selfAttenNotFitNucs.erase( new_end, end(selfAttenNotFitNucs) );
        }//for( const SandiaDecay::Nuclide * nuc : el_nucs.second )
      }//for( const auto &el_nucs : fitnucs )
      
      
      //We'll do some sanity checks to make sure the model we fit is consistent with the GUI - which
      //  it should always be
      if( fitnucs.size() != guiNucsFittingFrac.size() )  //never expect this to happen
        throw logic_error( "Number of calc self-atten nuclides does not equal num GUI self-atten nucs." );
        
      if( fitOtherEls.size() != guiElFittingOtherFrac.size() ) //never expect this to happen
        throw logic_error( "Number of non-source element being fit does not equal the GUI." );
      
      for( const SandiaDecay::Nuclide *nuc : fitnucs )
      {
        assert( nuc );
        if( !nuc )
          throw logic_error( "Null nuclide in fit result." );
        
        const auto pos = std::find( begin(guiNucsFittingFrac), end(guiNucsFittingFrac), nuc );
        if( pos == end(guiNucsFittingFrac) )
          throw logic_error( "The fit and GUI nuclides are not equal." ); //never expect this to happen
      }//for( const SandiaDecay::Nuclide *nuc : fitnucs )
      
      for( const SandiaDecay::Element *el : fitOtherEls )
      {
        assert( el );
        if( !el )
          throw logic_error( "Null element in fit result." );
        
        const auto pos = std::find( begin(guiElFittingOtherFrac), end(guiElFittingOtherFrac), el );
        if( pos == end(guiElFittingOtherFrac) )
          throw logic_error( "The fit and GUI other element components are not equal." ); //never expect this to happen
        
        //Make sure elements have at least one nuclide matching them.
        size_t num_nucs_of_el = 0;
        for( const SandiaDecay::Nuclide *nuc : fitnucs )
          num_nucs_of_el += (nuc->atomicNumber == el->atomicNumber);
        
        if( num_nucs_of_el < 1 )
          throw logic_error( "Element " + el->symbol + " fit non-source component, but didnt have"
                              " any source nuclides of this element." );
      }//for( const SandiaDecay::Element *el : fitOtherEls )

      
      map<short int,double> prefitsums, postfitsums; //atomic number to sum
      
      // The double is a fundamental type, so they would be initialized to 0.0 the first time
      //  they are accessed, but this makes me queasy, so we'll explicitly do it
      for( const SandiaDecay::Nuclide *nuc : fitnucs )
      {
        prefitsums[nuc->atomicNumber] = 0.0;
        postfitsums[nuc->atomicNumber] = 0.0;
      }

      const map<const SandiaDecay::Element *,vector<tuple<const SandiaDecay::Nuclide *,double,bool>>>
          guiMassFracs = select->sourceNuclideMassFractions();
      
      for( const auto &el_nucs : guiMassFracs )
      {
        const SandiaDecay::Element * const el = el_nucs.first;
        assert( el );
        if( !el )
          throw std::logic_error( "nullptr Element???" );
        
        for( const tuple<const SandiaDecay::Nuclide *,double,bool> &nuc_frac_fit : el_nucs.second )
        {
          if( get<2>(nuc_frac_fit) )
            prefitsums[el->atomicNumber] += get<1>(nuc_frac_fit);
        }//for( const auto &nuc_frac_fit : guiMassFracs )
      }//for( const auto &el_nucs : guiMassFracs )
      
      const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
      map<const SandiaDecay::Element *,vector<ShieldingSelect::MassFracInfo>> fitMassFracs;
      
      for( const auto &el_nucs : fit_el_to_nucs )
      {
        const SandiaDecay::Element * const el = el_nucs.first;
        assert( el );
        
        for( const SandiaDecay::Nuclide * nuc : el_nucs.second )
        {
          double frac, uncert;
          m_currentFitFcn->massFractionOfElement( frac, uncert, shielding_index, nuc, el, paramValues, paramErrors );
          
          postfitsums[el->atomicNumber] += frac;
          const bool isFit = nuc ? m_currentFitFcn->isVariableMassFraction( shielding_index, nuc )
                                 : m_currentFitFcn->isVariableOtherMassFraction( shielding_index, el );
          
          ShieldingSelect::MassFracInfo info;
          info.m_nuclide = nuc;
          info.m_fraction = frac;
          info.m_frac_uncert = uncert;
          info.m_fit_mass_frac = isFit;
          info.m_use_as_source = true;
          
          fitMassFracs[el].push_back( info );
        }//for( const SandiaDecay::Nuclide * nuc : el_nucs.second )
      }//for( const auto &el_nucs : fit_el_to_nucs )
      
      for( const SandiaDecay::Element *el : fitOtherEls )
      {
        double frac, uncert;
        const SandiaDecay::Nuclide * const nuc = nullptr;
        
        m_currentFitFcn->massFractionOfElement( frac, uncert, shielding_index, nuc, el,
                                      paramValues, paramErrors );
      
        postfitsums[el->atomicNumber] += frac;
        
        ShieldingSelect::MassFracInfo info;
        info.m_nuclide = nuc;
        info.m_fraction = frac;
        info.m_frac_uncert = uncert;
        info.m_fit_mass_frac = true;
        info.m_use_as_source = true;
        
        fitMassFracs[el].push_back( info );
      }//for( const SandiaDecay::Element *el : fitOtherEls )
      
      for( const SandiaDecay::Nuclide *nuc : selfAttenNotFitNucs )
      {
        const SandiaDecay::Element * const el = db->element( nuc->atomicNumber );
        assert( el );
        
        const double frac = m_currentFitFcn->massFractionOfElement( shielding_index, nuc, paramValues );
        
        ShieldingSelect::MassFracInfo info;
        info.m_nuclide = nuc;
        info.m_fraction = frac;
        info.m_frac_uncert = 0.0;
        info.m_fit_mass_frac = false;
        info.m_use_as_source = true;
        
        fitMassFracs[el].push_back( info );
      }//for( const SandiaDecay::Nuclide *nuc : selfAttenNotFitNucs )
      
      
      select->setMassFractions( fitMassFracs );
      
      assert( postfitsums.size() == prefitsums.size() );
      
      for( const auto &el_sum : prefitsums )
      {
        short int atomic_number = el_sum.first;
        assert( postfitsums.count(atomic_number) );
        if( !postfitsums.count(atomic_number) )
          throw logic_error( "postfitsums is missing an entry prefitsum has." );
        
        const double prefitsum = el_sum.second;
        const double sumfracs = postfitsums[atomic_number];
        
        const double frac_diff = fabs(sumfracs - prefitsum);
        assert( fabs(frac_diff) < 1.0E-5*std::max(sumfracs,prefitsum) );
        
        if( (frac_diff > 1.0E-12) && ((frac_diff/std::max(sumfracs,prefitsum)) > 1.0E-5) ) //limits chosen arbitrarily
          throw logic_error( "Mass fraction for of self-atten src atomic number " + std::to_string(atomic_number)
                            + " should be " + to_string(prefitsum) + " but calculation yielded "
                            + to_string(sumfracs) );
      }//for( const auto &el_sum : prefitsums )
    }//for( int i = 0; i < nshieldings; ++i )
    
    
    // Next we'll update shielding thicknesses, as this could affect setting trace source values
    for( int i = 0; i < nshieldings; ++i )
    {
      ShieldingSelect *select = gui_shieldings[i];
      if( select->isGenericMaterial() )
      {
        const double adUnits = PhysicalUnits::gram / PhysicalUnits::cm2;
        const double an = m_currentFitFcn->atomicNumber( i, paramValues );
        const double ad = m_currentFitFcn->arealDensity( i, paramValues ) / adUnits;
        
        select->m_atomicNumberEdit->setValue( static_cast<float>(an) );
        select->m_arealDensityEdit->setValue( static_cast<float>(ad) );
      }else
      {
        
        assert( geometry() == m_currentFitFcn->geometry() );
        if( geometry() != m_currentFitFcn->geometry() )
          throw runtime_error( "Geometry somehow changed during fitting process" );
        
        auto setEditTxt = []( WLineEdit *edit, const double val, const double err ){
          WString txt;
          if( err > 0.0 )
            txt = PhysicalUnits::printToBestLengthUnits( val, err );
          else
            txt = PhysicalUnits::printToBestLengthUnits( val );
          
          edit->setText( txt );
        };//setEditTxt lambda
        
        
        switch( geometry() )
        {
          case GeometryType::Spherical:
          {
            const double thickness = final_shieldings[i].m_dimensions[0];
            const double thicknessErr = final_shieldings[i].m_dimensionUncerts[0];
            assert( thickness == m_currentFitFcn->sphericalThickness(i, paramValues) );
            assert( std::max(0.0,thicknessErr) == std::max(0.0,m_currentFitFcn->sphericalThickness(i,paramErrors)) );
            
            //const double thickness = m_currentFitFcn->sphericalThickness( i, paramValues );
            //const double thicknessErr = m_currentFitFcn->sphericalThickness( i, paramErrors );
            
            setEditTxt( select->m_thicknessEdit, thickness, thicknessErr );
            
            break;
          }//case GeometryType::Spherical:
            
          case GeometryType::CylinderEndOn:
          case GeometryType::CylinderSideOn:
          {
            const double dr = final_shieldings[i].m_dimensions[0];
            const double drErr = final_shieldings[i].m_dimensionUncerts[0];
            assert( dr == m_currentFitFcn->cylindricalRadiusThickness( i, paramValues ) );
            
            //const double dr = m_currentFitFcn->cylindricalRadiusThickness( i, paramValues );
            //const double drErr = m_currentFitFcn->cylindricalRadiusThickness( i, paramErrors );
            
            setEditTxt( select->m_cylRadiusEdit, dr, drErr );
            
            const double dz = final_shieldings[i].m_dimensions[1];
            const double dzErr = final_shieldings[i].m_dimensionUncerts[1];
            assert( dz == m_currentFitFcn->cylindricalLengthThickness( i, paramValues ) );
            
            //const double dz = m_currentFitFcn->cylindricalLengthThickness( i, paramValues );
            //const double dzErr = m_currentFitFcn->cylindricalLengthThickness( i, paramErrors );
            setEditTxt( select->m_cylLengthEdit, dz, dzErr );
            
            break;
          }//case GeometryType::CylinderEndOn or CylinderSideOn:
            
            
          case GeometryType::Rectangular:
          {
            const double dw = final_shieldings[i].m_dimensions[0];
            const double dwErr = final_shieldings[i].m_dimensionUncerts[0];
            assert( dw == m_currentFitFcn->rectangularWidthThickness( i, paramValues ) );
            
            //const double dw = m_currentFitFcn->rectangularWidthThickness( i, paramValues );
            //const double dwErr = m_currentFitFcn->rectangularWidthThickness( i, paramErrors );
            setEditTxt( select->m_rectWidthEdit, dw, dwErr );
            
            const double dh = final_shieldings[i].m_dimensions[1];
            const double dhErr = final_shieldings[i].m_dimensionUncerts[1];
            assert( dh == m_currentFitFcn->rectangularHeightThickness( i, paramValues ) );
            
            //const double dh = m_currentFitFcn->rectangularHeightThickness( i, paramValues );
            //const double dhErr = m_currentFitFcn->rectangularHeightThickness( i, paramErrors );
            setEditTxt( select->m_rectHeightEdit, dh, dhErr );
            
            const double dd = final_shieldings[i].m_dimensions[2];
            const double ddErr = final_shieldings[i].m_dimensionUncerts[2];
            assert( dd == m_currentFitFcn->rectangularDepthThickness( i, paramValues ) );
            
            //const double dd = m_currentFitFcn->rectangularDepthThickness( i, paramValues );
            //const double ddErr = m_currentFitFcn->rectangularDepthThickness( i, paramErrors );
            setEditTxt( select->m_rectDepthEdit, dd, ddErr );
            
            break;
          }//case GeometryType::Rectangular:
            
          case GeometryType::NumGeometryType:
            assert( 0 );
            break;
        }//switch( geometry() )
        
        
        //        const double volume = (PhysicalUnits::pi*4.0/3.0) * ( pow(radius + thickness,3) - pow(radius,3) );
        //        const double density = select->material()->density;
        //        cerr << "From Geometry, material " << select->material()->name << " has "
        //             << "mass " << (volume*density)/PhysicalUnits::gram << " g" << endl;
        //        radius += thickness;
      }//if( genericMaterial ) / else
    }//for( size_t i = 0; i < shieldings.size(); ++i )
    

    
    // Finally we'll set activities and ages
    assert( results->fit_src_info.size() == m_sourceModel->underlyingData().size() );
    m_sourceModel->setUnderlyingData( results->fit_src_info );
    
    const size_t nnucs = m_currentFitFcn->numNuclides();
    
    // Go through and update trace source activities displayed in the shieldings
    for( const ShieldingSourceFitCalc::IsoFitStruct &src : results->fit_src_info )
    {
      if( src.sourceType != ShieldingSourceFitCalc::ModelSourceType::Trace )
        continue;
      
      for( WWidget *widget : m_shieldingSelects->children() )
      {
        ShieldingSelect *select = dynamic_cast<ShieldingSelect *>(widget);
        if( select && select->isTraceSourceForNuclide(src.nuclide) )
        {
          select->setTraceSourceTotalActivity( src.nuclide, src.activity, src.activityUncertainty, false );
          break;
        }
      }//for( WWidget *widget : m_shieldingSelects->children() )
    }//for( const ShieldingSourceFitCalc::IsoFitStruct &src : results->fit_src_info )
    
    
    /*
    //Go through and set the ages and activities fit for
    for( int nucn = 0; nucn < static_cast<int>(nnucs); ++nucn )
    {
      // TODO: 20230712: do we need to set ShieldingSelect TraceSource activities? Its possible they got set in the setUnderlyingData(...) call
      //
      const SandiaDecay::Nuclide * const nuc = m_sourceModel->nuclide(nucn);
      const ShieldingSourceFitCalc::ModelSourceType sourceType = m_sourceModel->sourceType(nucn);
      
      //Even if we didnt explicitly fit for the activity parameter, we may still have effectively
      //  fit for the activity if we fit for any shielding dimensions, so rather than trying to
      //  detect when we may have effectively fit for activity, we'll just set for anything more
      //  then the most obvios "we didnt fit for it" sceneriou.
      bool effectivelyFitActivity = true;
      switch( sourceType )
      {
        case ShieldingSourceFitCalc::ModelSourceType::Point:
          effectivelyFitActivity = m_sourceModel->fitActivity(nucn);
          break;
          
        case ShieldingSourceFitCalc::ModelSourceType::Intrinsic:
          break;
          
        case ShieldingSourceFitCalc::ModelSourceType::Trace:
        {
          success = false;
          for( WWidget *widget : m_shieldingSelects->children() )
          {
            ShieldingSelect *select = dynamic_cast<ShieldingSelect *>(widget);
            
            if( select && select->isTraceSourceForNuclide(nuc) )
            {
              success = true;
              select->setTraceSourceTotalActivity( nuc, total_activity );
            }
          }//for( WWidget *widget : m_shieldingSelects->children() )
          
          break;
        }
      }//switch( sourceType )
          
      if( effectivelyFitActivity )
      {
        switch( sourceType )
        {
          case ShieldingSourceFitCalc::ModelSourceType::Point:
          case ShieldingSourceFitCalc::ModelSourceType::Intrinsic:
            break;
            
          case ShieldingSourceFitCalc::ModelSourceType::Trace:
          {
            
            break;
          }//case ShieldingSourceFitCalc::ModelSourceType::Trace:
        }//switch( m_sourceModel->sourceType(ison) )
        
        if( !success )
        {
          const string msg = "An invalid activity was calculated for " + nuc->symbol
          + ", other results may be invalid to";
          passMessage( msg, WarningWidget::WarningMsgHigh );
        }//if( we failed to set activity value. )
        
        
        try
        {
          const double activityUncert = m_currentFitFcn->totalActivityUncertainty( nuc, paramValues, paramErrors );
          
          if( activityUncert < FLT_EPSILON )
          {
            m_sourceModel->setData( actUncertIndex, boost::any() );
          }else
          {
            char actUncertStr[64] = { '\0' };
            snprintf( actUncertStr, sizeof(actUncertStr), "%f bq", (activityUncert/PhysicalUnits::becquerel) );
            
            success = m_sourceModel->setData( actUncertIndex, WString(actUncertStr) );
            if( !success )
            {
              const string msg = "Calculated activity uncertainty for " + nuc->symbol
              + " is invalid, other results may be invalid to.";
              passMessage( msg, WarningWidget::WarningMsgHigh );
            }//if( IsNan(p.second) || IsInf(p.second) )
          }//if( activity uncert < 0 ) / else
        }catch( std::exception &e )
        {
          const string msg = "Unexpected error calculating activity uncertainty for " + nuc->symbol
          + ": " + string(e.what());
          passMessage( msg, WarningWidget::WarningMsgHigh );
        }
      }else if( m_sourceModel->activityUncert(ison) >= 0.0 )
      {
        m_sourceModel->setData( actUncertIndex, boost::any() );
      }//if( effectivelyFitActivity )
    }//for( int ison = 0; ison < niso; ++ison )
    */
    
    
    updateChi2ChartActual( results );
    m_chi2ChartNeedsUpdating = false;
    updateCalcLogWithFitResults( m_currentFitFcn, results, m_calcLog );
  }catch( std::exception &e )
  {
    passMessage( "Programming issue - caught exception: " + string(e.what())
                 + "<br />Application state may be suspect!", WarningWidget::WarningMsgHigh );
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, ("Programming Issue - caught exception: " + string(e.what())).c_str() );
#endif
  }//try / catch
  
  m_currentFitFcn.reset();
}//void updateGuiWithModelFitResults( std::vector<double> paramValues, paramErrors )


std::shared_ptr<ShieldingSourceFitCalc::ModelFitResults> ShieldingSourceDisplay::doModelFit( const bool fitInBackground )
{
  try
  {
    checkAndWarnZeroMassFraction();
  }catch( std::exception &e )
  {
    passMessage( WString::tr("ssd-err-fit-not-performed").arg(e.what()), WarningWidget::WarningMsgHigh );
    return nullptr;
  }
    
  try
  {
    checkForMultipleGenericMaterials();
  }catch( exception &e )
  {
    passMessage( WString::tr("ssd-err-fit-not-performed").arg(e.what()), WarningWidget::WarningMsgHigh );
    return nullptr;
  }//try / catch
  
  m_modifiedThisForeground = true;
  
  try
  {
    checkDistanceAndThicknessConsistent();
  }catch( exception &e )
  {
    passMessage( WString::tr("ssd-err-before-fit").arg(e.what()), WarningWidget::WarningMsgHigh );
  }//try / catch
  
  shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn> chi2Fcn;
  vector<string> errormsgs;
  vector<ShieldingSelect *> shieldings;
  
  //make sure fitting for at least one nuclide:
  
  auto inputPrams = make_shared<ROOT::Minuit2::MnUserParameters>();
  std::vector<ShieldingSourceFitCalc::ShieldingInfo> initial_shieldings;
  
  try
  {
    for( WWidget *widget : m_shieldingSelects->children() )
    {
      ShieldingSelect *select = dynamic_cast<ShieldingSelect *>( widget );
      if( select )
        initial_shieldings.push_back( select->toShieldingInfo() );
    }//for( WWidget *widget : m_shieldingSelects->children() )
    
    auto fcnAndPars = shieldingFitnessFcn();
    
    chi2Fcn = fcnAndPars.first;
    *inputPrams = fcnAndPars.second;
  }catch( std::exception &e )
  {
    
    passMessage( WString::tr("ssd-err-couldnt-make-chi2-fcn").arg(e.what()),
                WarningWidget::WarningMsgHigh );
    return nullptr;
  }//try / catch
  
  setWidgetStateForFitStarting();

  //Need to disable "All Peaks", Detector, Distance, and "Material", and "Generic"
  
  {
    std::lock_guard<std::mutex> lock( m_currentFitFcnMutex );
    m_currentFitFcn = chi2Fcn;
  }
  
  auto results = make_shared<ShieldingSourceFitCalc::ModelFitResults>();
  results->successful = ShieldingSourceFitCalc::ModelFitResults::FitStatus::InvalidOther;
  
  auto progress = std::make_shared<ShieldingSourceFitCalc::ModelFitProgress>();
  boost::function<void()> progress_updater = wApp->bind( boost::bind( &ShieldingSourceDisplay::updateGuiWithModelFitProgress, this, progress ) );
  
  //Wrap the GUI update with WApplication::bind in case this
  //  ShieldingSourceDisplay widget gets deleted before the computation is over
  boost::function<void()> gui_updater = wApp->bind( boost::bind( &ShieldingSourceDisplay::updateGuiWithModelFitResults, this, results ) );
  
  const string sessionid = wApp->sessionId();
  if( fitInBackground )
  {
    Wt::WServer *server = Wt::WServer::instance();
    server->ioService().boost::asio::io_service::post( boost::bind( &ShieldingSourceFitCalc::fit_model,
                            sessionid, chi2Fcn, inputPrams, progress, progress_updater, results, gui_updater ) );
  }else
  {
    ShieldingSourceFitCalc::fit_model( sessionid, chi2Fcn, inputPrams, progress, progress_updater, results, gui_updater );
  }
  
  return results;
}//void doModelFit()





void ShieldingSourceDisplay::updateCalcLogWithFitResults(
                                  shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn> chi2Fcn,
                                    std::shared_ptr<ShieldingSourceFitCalc::ModelFitResults> results,
                                                         vector<string> &calcLog )
{
  // This function is not internationalized - the plan is to eventually totally overhaul this
  //  how results are logged, so it isnt worth internationalizing this function now (and also,
  //  probably only the primary developer uses this calc log).
  assert( chi2Fcn );
  assert( results );
  const std::vector<double> &params = results->paramValues;
  const std::vector<double> &errors = results->paramErrors;
  
  if( calcLog.size() && calcLog.back() == ns_no_uncert_info_txt )
    calcLog.erase( calcLog.end()-1, calcLog.end() );
  
  const shared_ptr<const DetectorPeakResponse> &det = chi2Fcn->detector();
  
  const DetectorPeakResponse::EffGeometryType detType = (det && det->isValid())
                                                  ? det->geometryType()
                                                  : DetectorPeakResponse::EffGeometryType::FarField;
  
  try
  {
    for( size_t shielding_index = 0; shielding_index < chi2Fcn->numMaterials(); ++shielding_index )
    {
      if( !chi2Fcn->hasVariableMassFraction(shielding_index) )
        continue;
      
      const Material *mat = chi2Fcn->material(shielding_index);
      assert( mat );
      if( !mat )
        continue;
      
      stringstream msg;
      msg << "Shielding material " << mat->name << " fit mass fractions for isotopes:";
      
      const map<const SandiaDecay::Element *,vector<const SandiaDecay::Nuclide *>> el_to_nucs
                                            = chi2Fcn->nuclideFittingMassFracFor( shielding_index );
      for( const auto &el_nuc : el_to_nucs )
      {
        const SandiaDecay::Element * const el = el_nuc.first;
        
        for( const SandiaDecay::Nuclide *n : el_nuc.second )
        {
          double frac, uncert;
          chi2Fcn->massFractionOfElement( frac, uncert, shielding_index, n, el, params, errors );
          
          msg << " " << n->symbol << "(massfrac=" << frac << "+-" << uncert << "),";
        }//for( size_t shielding_index = 0; shielding_index < chi2Fcn->numMaterials(); ++shielding_index )
      }
      
      calcLog.push_back( msg.str() );
    }//for( const Material *mat : chi2Fcn->materialsFittingMassFracsFor() )
    
  {//begin add chi2 line
    stringstream msg;
    msg << "It took " << results->num_fcn_calls
        << " solution trials to reach chi2=" << results->chi2
        << " with an estimated distance to minumum of " << results->edm;
    calcLog.push_back( msg.str() );
  }//end add chi2 line
    
  //Need to list fit parameters and uncertainties here
  const size_t nnuc = chi2Fcn->numNuclides();
  for( size_t nucn = 0; nucn < nnuc; ++nucn )
  {
    const SandiaDecay::Nuclide *nuc = chi2Fcn->nuclide( nucn );
    if( nuc )
    {
      const bool useCi = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
      const double act = chi2Fcn->activity( nuc, params );
      const string actStr = PhysicalUnits::printToBestActivityUnits( act, 2, useCi );
      
      const double actUncert = chi2Fcn->activityUncertainty( nuc, params, errors );
      const string actUncertStr = PhysicalUnits::printToBestActivityUnits( actUncert, 2, useCi );
      
      const double mass = (act / nuc->activityPerGram()) * PhysicalUnits::gram;
      const std::string massStr = PhysicalUnits::printToBestMassUnits( mass, 2, PhysicalUnits::gram );
      
      const double age = chi2Fcn->age( nuc, params );
      const double ageUncert = chi2Fcn->age( nuc, errors );
      const string ageStr = PhysicalUnitsLocalized::printToBestTimeUnits( age, 2 );
      const string ageUncertStr = PhysicalUnitsLocalized::printToBestTimeUnits( ageUncert, 2 );
  
      string act_postfix = DetectorPeakResponse::det_eff_geom_type_postfix(detType), trace_total = "";
      if( chi2Fcn->isTraceSource(nuc) )
      {
        const double total_act = chi2Fcn->totalActivity(nuc,params);
        trace_total = "Total activity "
                      + PhysicalUnits::printToBestActivityUnits( total_act, 2, useCi )
                      + DetectorPeakResponse::det_eff_geom_type_postfix(detType) + ", ";
        
        switch( chi2Fcn->traceSourceActivityType(nuc) )
        {
          case TraceActivityType::TotalActivity:
            trace_total = "";
            break;
            
          case TraceActivityType::ActivityPerCm3:
            act_postfix = " per cm3";
            break;
          
          case TraceActivityType::ExponentialDistribution:
            act_postfix = " per m2, with relaxation length "
                          + PhysicalUnits::printToBestLengthUnits(chi2Fcn->relaxationLength(nuc) );
            break;
            
          case GammaInteractionCalc::TraceActivityType::ActivityPerGram:
            act_postfix = " per gram shielding";
            break;
            
          case GammaInteractionCalc::TraceActivityType::NumTraceActivityType:
            assert(0);
            break;
        }//switch( chi2Fcn->traceSourceActivityType(nuc) )
      }//if( chi2Fcn->isTraceSource(nuc) )
      
      stringstream msg;
      msg << nuc->symbol << " fit activity " << actStr << act_postfix
          << " (" << trace_total << nuc->symbol << " mass: " << massStr
          << act_postfix << ") with uncertainty " << actUncertStr << act_postfix
          << " (" << floor(0.5 + 10000*actUncert/act)/100.0 << "%)";
      
      if( ageUncert <= DBL_EPSILON )
        msg << " at assumed age " << ageStr;
      else
        msg << " with age " << ageStr << "+- " << ageUncertStr;
      
      if( chi2Fcn->isSelfAttenSource(nuc) )
        msg << ", a self attenuating source";
      else if( chi2Fcn->isTraceSource(nuc) )
        msg << ", a trace source";
  
      msg << ".";
      
      calcLog.push_back( msg.str() );
    }//if( nuc )
  }//for( size_t nucn = 0; nucn < nnuc; ++nucn )
  
    calcLog.push_back( "Geometry: " + string(GammaInteractionCalc::to_str(chi2Fcn->geometry())) );
    
  const int nmat = static_cast<int>( chi2Fcn->numMaterials() );
  for( int matn = 0; matn < nmat; ++matn )
  {
    stringstream msg;
    const Material *mat = chi2Fcn->material( matn );
    if( !mat )
    {
      const double adUnits = PhysicalUnits::gram / PhysicalUnits::cm2;
      const double ad = chi2Fcn->arealDensity( matn, params ) / adUnits;
      const double adUncert = chi2Fcn->arealDensity( matn, errors ) / adUnits;
      
      const double an = chi2Fcn->atomicNumber( matn, params );
      const double anUncert = chi2Fcn->atomicNumber( matn, errors );
      
      msg << std::setprecision(3) << "Shielding " << matn+1 << " has "
      << "AtomicNumber=" << an;
      if( anUncert > DBL_EPSILON )
        msg << " (+-" << anUncert << ")";
      msg << " and ArealDensity=" << ad;
      if( adUncert > DBL_EPSILON )
        msg << " (+-" << adUncert << ")";
      msg << " g/cm2";
      
    }else //if( !mat )
    {
      //assert( geometry() == chi2Fcn->geometry() );
      
      const double density = mat->density * PhysicalUnits::cm3 / PhysicalUnits::gram;
      msg << mat->name << " has density " << std::setprecision(3) << density << "g/cm3 ";
      
      switch( chi2Fcn->geometry() )
      {
        case GeometryType::Spherical:
        {
          const double thickness = chi2Fcn->sphericalThickness( matn, params );
          const double thicknessUncert = chi2Fcn->sphericalThickness( matn, errors );
          
          if( thicknessUncert > DBL_EPSILON )
          {
            msg << "and fit thickness "
                << PhysicalUnits::printToBestLengthUnits(thickness,thicknessUncert) << ".";
          }else
          {
            msg << "and has fixed thickness " << PhysicalUnits::printToBestLengthUnits(thickness)
                 << ".";
          }
          
          break;
        }//case GeometryType::Spherical:
          
        case GeometryType::CylinderEndOn:
        case GeometryType::CylinderSideOn:
        {
          msg << "and dimensions [";
          const double radThickness = chi2Fcn->cylindricalRadiusThickness( matn, params );
          const double radThicknessUncert = chi2Fcn->cylindricalRadiusThickness( matn, errors );
          const double lenThickness = chi2Fcn->cylindricalLengthThickness( matn, errors );
          const double halfLenUncert = chi2Fcn->cylindricalLengthThickness( matn, errors );
          
          const char *radLabel = matn ? "radial-thickness=" : "radius=";
          const char *lengthLabel = matn ? "length-thickness=" : "half-length=";
          
          if( radThicknessUncert > DBL_EPSILON )
            msg << radLabel << PhysicalUnits::printToBestLengthUnits(radThickness,radThicknessUncert);
          else
            msg << radLabel << PhysicalUnits::printToBestLengthUnits(radThickness) << " (fixed), ";
          
          if( lenThickness > DBL_EPSILON )
            msg << lengthLabel << PhysicalUnits::printToBestLengthUnits(radThickness,radThicknessUncert);
          else
            msg << lengthLabel << PhysicalUnits::printToBestLengthUnits(radThickness) << " (fixed)";
          
          msg << "]";
          
          break;
        }//case GeometryType::CylinderEndOn and CylinderSideOn:
                    
        case GeometryType::Rectangular:
        {
          const double widthThickness = chi2Fcn->rectangularWidthThickness( matn, params );
          const double widthThicknessUncert = chi2Fcn->rectangularWidthThickness( matn, errors );
          const double heightThickness = chi2Fcn->rectangularHeightThickness( matn, params );
          const double heightThicknessUncert = chi2Fcn->rectangularHeightThickness( matn, errors );
          const double depthThickness = chi2Fcn->rectangularDepthThickness( matn, params );
          const double depthThicknessUncert = chi2Fcn->rectangularDepthThickness( matn, errors );
          
          const char *widthLabel  = matn ? "width-thickness="  : "half-width=";
          const char *heightLabel = matn ? "height-thickness=" : "half-height=";
          const char *depthLabel  = matn ? "depth-thickness="  : "half-depth=";
          
          msg << widthLabel;
          if( widthThicknessUncert > DBL_EPSILON )
            msg << PhysicalUnits::printToBestLengthUnits(widthThickness,widthThicknessUncert) << ", ";
          else
            msg << PhysicalUnits::printToBestLengthUnits(widthThickness) << " (fixed), ";
          
          msg << heightLabel;
          if( heightThicknessUncert > DBL_EPSILON )
            msg << PhysicalUnits::printToBestLengthUnits(heightThickness,heightThicknessUncert) << ", ";
          else
            msg << PhysicalUnits::printToBestLengthUnits(heightThickness) << " (fixed), ";
          
          msg << depthLabel;
          if( depthThicknessUncert > DBL_EPSILON )
            msg << PhysicalUnits::printToBestLengthUnits(depthThickness,depthThicknessUncert);
          else
            msg << PhysicalUnits::printToBestLengthUnits(depthThickness) << " (fixed)";
          
          break;
        }//case GeometryType::Rectangular:
          
        case GammaInteractionCalc::GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }//switch( geometry() )
    }//if( !mat ) / else
    
    calcLog.push_back( msg.str() );
  }//for( size_t matn = 0; matn < nmat; ++matn )
  }catch( std::exception & )
  {
    calcLog.push_back( "There was an error and log may not be complete." );
  }
}//updateCalcLogWithFitResults(...)









