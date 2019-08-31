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

#include <vector>
#include <algorithm>

#include <boost/tuple/tuple.hpp>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WCheckBox>
#include <Wt/WLineEdit>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WStringStream>
#include <Wt/WDoubleSpinBox>
#include <Wt/WContainerWidget>
#include <Wt/WSuggestionPopup>
#include <Wt/WRegExpValidator>

#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/ColorSelect.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DetectorEdit.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/SpectrumChart.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/SpecMeasManager.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/CanvasForDragging.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeSelectionAids.h"
#include "InterSpec/ShieldingSourceDisplay.h"
#include "InterSpec/IsotopeNameFilterModel.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"

using namespace std;
using namespace Wt;

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

const int ReferencePhotopeakDisplay::sm_xmlSerializationVersion = 0;

const int DecayParticleModel::RowData::XRayDecayMode = 1000;
const int DecayParticleModel::RowData::ReactionToGammaMode = 1001;
const int DecayParticleModel::RowData::NormGammaDecayMode = 1002;

namespace
{
  //See: https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
  //  for a pallete of distint colors
  //Or http://artshacker.com/wp-content/uploads/2014/12/Kellys-22-colour-chart.jpg
  const static vector<Wt::WColor> ns_def_line_colors{
    {"#0000FF"}, {"#006600"}, {"#0099FF"}, {"#9933FF"},
    {"#FF66FF"}, {"#CC3333"}, {"#FF6633"}, {"#FFFF99"},
    {"#CCFFCC"}, {"#0000CC"}, {"#666666"}, {"#003333"}
  };
}//namespace

bool DecayParticleModel::less_than( const DecayParticleModel::RowData &lhs,
                                    const DecayParticleModel::RowData &rhs,
                                    const DecayParticleModel::Column c,
                                    const SortOrder order )
{

  bool less = false;

  switch( c )
  {
    case kEnergy:         less = (lhs.energy < rhs.energy);           break;
    case kBranchingRatio: less = (lhs.branchRatio < rhs.branchRatio); break;
    case kResponsibleNuc:
      if( !lhs.responsibleNuc || !rhs.responsibleNuc )
        less = false;
      else
        less = SandiaDecay::Nuclide::lessThan( lhs.responsibleNuc, rhs.responsibleNuc );
    break;
    case kDecayMode:      less = (lhs.decayMode < rhs.decayMode);     break;
    case kParticleType:   less = (lhs.particle < rhs.particle);       break;
    case kNumColumn:      return false;
  }//switch( r )

  if( order == Wt::AscendingOrder )
    return less;
  return !less;
}//less_than(...)


DecayParticleModel::DecayParticleModel( Wt::WObject *parent )
  : WAbstractItemModel( parent ),
    m_sortColumn( DecayParticleModel::kEnergy ),
    m_sortOrder( Wt::AscendingOrder )
{
}

DecayParticleModel::~DecayParticleModel()
{
}


WFlags<ItemFlag> DecayParticleModel::flags( const WModelIndex &p ) const
{
  if( p.isValid() && (p.column()==kDecayMode || p.column()==kParticleType) )
    return ItemIsXHTMLText;
  return WFlags<ItemFlag>();
}


boost::any DecayParticleModel::data( const WModelIndex &index, int role ) const
{
  using namespace SandiaDecay;

  if( role != DisplayRole )
    return boost::any();

  const int row = index.row();
  const int nrow = static_cast<int>( m_data.size() );
  if( row < 0 || row >= nrow )
    return boost::any();

  const int column = index.column();
  if( column < 0 || column >= kNumColumn )
    return boost::any();

  const RowData &dataRow = m_data[row];
  
  char buffer[64];
  
  switch( column )
  {
    case kEnergy:
      snprintf( buffer, sizeof(buffer), "%.4f", dataRow.energy );
      return WString( buffer );
      
    case kBranchingRatio:
      snprintf( buffer, sizeof(buffer), "%.4g", dataRow.branchRatio );
      return WString( buffer );

    case kResponsibleNuc:
      if( dataRow.responsibleNuc )
        return WString( dataRow.responsibleNuc->symbol );
      return boost::any();

    case kDecayMode:
    {
      switch( dataRow.decayMode )
      {
#ifndef WT_NO_STD_WSTRING
        case SandiaDecay::AlphaDecay:               return WString( L"\x03B1" );
        case SandiaDecay::BetaDecay:                return WString( L"\x03B2<sup>-</sup>" );
        case SandiaDecay::BetaPlusDecay:            return WString( L"\x03B2<sup>+</sup>" );
        case SandiaDecay::DoubleBetaDecay:          return WString( L"double \x03B2" );
#else
        case SandiaDecay::AlphaDecay:               return WString( "&alpha;" );
        case SandiaDecay::BetaDecay:                return WString( "&beta;<sup>-</sup>" );
        case SandiaDecay::BetaPlusDecay:            return WString( "&beta;<sup>+</sup>" );
        case SandiaDecay::DoubleBetaDecay:          return WString( "double &beta;" );
#endif
        case SandiaDecay::IsometricTransitionDecay: return WString( "Iso" );
        case SandiaDecay::ElectronCaptureDecay:     return WString( "e.c." );
        case SandiaDecay::ProtonDecay:              return WString( "proton" );
        case SandiaDecay::SpontaneousFissionDecay:  return WString( "s.f." );
        case SandiaDecay::Carbon14Decay:            return WString( "C14" );
        case RowData::XRayDecayMode:                return WString( "xray" );
        case RowData::ReactionToGammaMode:          return WString( "Reaction" );
        case RowData::NormGammaDecayMode:           return WString( "NORM" );
      }//switch( dataRow.decayMode )

      return boost::any();
    }//case kDecayMode:

    case kParticleType:
    {
      switch( dataRow.particle )
      {
#ifndef WT_NO_STD_WSTRING
        case BetaParticle:   return WString( L"\x03B2<sup>-</sup>" );
        case GammaParticle:  return WString( L"\x0263" );
        case AlphaParticle:          return WString( L"\x03B1" );
#else
        case BetaParticle:   return WString( "&beta<sup>-</sup>" );
        case GammaParticle:  return WString( "&gamma;" );
        case AlphaParticle:  return WString( "&alpha;" );
#endif
        case PositronParticle:        return WString( "e<sup>+</sup>" );
        case CaptureElectronParticle: return WString( "ec" );
       case XrayParticle:            return WString( "xray" );
      }//switch( dataRow.particle )
      return boost::any();
    }//case kParticleType:
  }//switch( column )

  return boost::any();
}//boost::any data( const Wt::WModelIndex &index, int role  )


boost::any DecayParticleModel::headerData( int column,
                                           Orientation orientation,
                                           int role ) const
{
  if( role == LevelRole )
    return 0;

  if( (orientation != Horizontal)
      || ((role != DisplayRole) && (role != ToolTipRole)) )
    return WAbstractItemModel::headerData( column, orientation, role );

  //If we are here, we want the column title
  if( role == DisplayRole )
  {
    switch( column )
    {
      case kEnergy:         return WString( "Energy (keV)" );
      case kBranchingRatio: return WString( "B.R." ); //"Phot/Decay" &gamma;/Decay  \u03B3/Decay  &#947;/Decay
      case kResponsibleNuc: return WString( "Parent" );
      case kDecayMode:      return WString( "Mode" );
      case kParticleType:   return WString( "Particle" );
      case kNumColumn:      return boost::any();
    }//switch( column )
  }else if( role == ToolTipRole )
  {
    switch( column )
    {
      case kEnergy:
        return WString( "Energy of the particle produced" );
      case kBranchingRatio:
        return WString( "Intensity of the particle, relative to the highest"
                        " intensity of that particle. For gammas, this is after"
                        " the optional shielding and detector effects are"
                        " applied." );
      case kResponsibleNuc:
        return WString( "Actual nuclide which decayed to give this particle" );
      case kDecayMode:
        return WString( "Decay mode of the parent nuclide, which produced this"
                        " particle" );
      case kParticleType:
        return WString( "The type of particle produced" );
      case kNumColumn:
      return boost::any();
    }//switch( column )
  }//if( role == DisplayRole ) / else

  return boost::any();
}//headerData(...)


int DecayParticleModel::columnCount( const Wt::WModelIndex &parent ) const
{
  if( parent.isValid() )
    return 0;

  return kNumColumn;
}


int DecayParticleModel::rowCount( const Wt::WModelIndex &parent ) const
{
  if( parent.isValid() )
    return 0;
  return static_cast<int>( m_data.size() );
}


WModelIndex DecayParticleModel::parent( const WModelIndex & ) const
{
  return WModelIndex();
}


WModelIndex DecayParticleModel::index( int row, int column,
                                       const WModelIndex &parent ) const
{
  if( parent.isValid() )
    return WModelIndex();

  if( row >= rowCount() || column >= columnCount() )
    return WModelIndex();

  return createIndex( row, column, (void *)NULL );
}//WModelIndex index( int row, int column, const WModelIndex & ) const


void DecayParticleModel::sort( int column, Wt::SortOrder order )
{
  m_sortOrder = order;
  m_sortColumn = Column(column);

  if( m_data.empty() )
    return;
  
  vector<RowData> data = m_data;

  boost::function<bool(const RowData&,const RowData&)> comparer;
  comparer = boost::bind(&DecayParticleModel::less_than, _1, _2,
                                                  m_sortColumn, m_sortOrder);

  stable_sort( data.begin(), data.end(), comparer );
  m_data.swap( data );
  dataChanged().emit( index(0,0), index(rowCount()-1,columnCount()-1) );
//  clear();
//  if( data.size() )
//  {
//    beginInsertRows( WModelIndex(), 0, static_cast<int>(data.size())-1 );
//    m_data = data;
//    endInsertRows();
//  }//if( data.size() )
}//sort(...)

void DecayParticleModel::clear()
{
  const int lastRow = rowCount();
  if( lastRow <= 0 )
    return;
  
  beginRemoveRows( WModelIndex(), 0, lastRow -1);
  m_data.clear();
  endRemoveRows();
}//clear()



void DecayParticleModel::setRowData( const std::vector<RowData> &newData )
{
  if( m_data.size() )
    clear();
  beginInsertRows( WModelIndex(), 0, static_cast<int>(newData.size()) -1);
  m_data = newData;
  endInsertRows();

  sort( m_sortColumn, m_sortOrder );
  layoutChanged().emit();
}//setRowData(...)



ReferencePhotopeakDisplay::ReferencePhotopeakDisplay(
#if ( USE_SPECTRUM_CHART_D3 )
                                            D3SpectrumDisplayDiv *chart,
#else
                                            SpectrumDisplayDiv *chart,
#endif
                                            MaterialDB *materialDB,
                                            WSuggestionPopup *materialSuggest,
                                            InterSpec *specViewer,
                                            WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_chart( chart ),
    m_spectrumViewer( specViewer ),
    m_nuclideEdit( NULL ),
    m_nuclideSuggest( NULL ),
    m_ageEdit( NULL ),
    m_lowerBrCuttoff( NULL ),
    m_promptLinesOnly( NULL ),
    m_halflife( NULL ),
    m_persistLines( NULL ),
    m_clearLines( NULL ),
    //m_fitPeaks( NULL ),
    m_showGammas( NULL ),
    m_showXrays( NULL ),
    m_showAlphas( NULL ),
    m_showBetas( NULL ),
    m_detectorDisplay( NULL ),
    m_materialDB( materialDB ),
    m_materialSuggest( materialSuggest ),
    m_shieldingSelect( NULL ),
    m_particleView( NULL ),
    m_particleModel( NULL ),
    m_currentlyShowingNuclide(),
    m_colorSelect( nullptr ),
    m_userHasPickedColor( false ),
    m_peaksGetAssignedRefLineColor( false ),
    m_lineColors{ ns_def_line_colors }
{
  const char *tooltip = nullptr;
  
  m_currentlyShowingNuclide.reset();
  
  if( !chart )
    throw runtime_error( "ReferencePhotopeakDisplay: a valid chart"
                         " must be passed in" );

  const bool showToolTipInstantly
      = InterSpecUser::preferenceValue<bool>( "ShowTooltips", specViewer );
  
  //The inputDiv/Layout is the left side of the widget that holds all the
  //  nuclide input,age, color picker, DRF, etc
  WContainerWidget *inputDiv = new WContainerWidget();
  WGridLayout *inputLayout = new WGridLayout();
  inputLayout->setContentsMargins( 0, 0, 0, 0 );
  inputDiv->setLayout( inputLayout );
    
    
  const bool isPhone = m_spectrumViewer->isPhone();
    
  if( isPhone )
  {
    addStyleClass( "ReferencePhotopeakDisplayMobile" );
    //m_layout->setHorizontalSpacing( 5 );
  }else
  {
    addStyleClass( "ReferencePhotopeakDisplay" );
  }//if( specViewer->isPhone() ) / else
  
  
  
  const WLength labelWidth(3.5,WLength::FontEm), fieldWidth(4,WLength::FontEm);
  const WLength optionWidth(5.25,WLength::FontEm), buttonWidth(5.25,WLength::FontEm);
  
  //added overflow to prevent scroll bars
  setOverflow( Wt::WContainerWidget::OverflowHidden );
  WLabel *nucInputLabel = new WLabel( "Nuclide:" );
  nucInputLabel->setMinimumSize( labelWidth, WLength::Auto );
  m_nuclideEdit = new WLineEdit( "" );
  m_nuclideEdit->setMargin( 1 );
//  m_nuclideEdit->setMinimumSize( WLength(10,WLength::FontEx), WLength::Auto );
  m_nuclideEdit->setMinimumSize( fieldWidth, WLength::Auto );
  m_nuclideEdit->setAutoComplete( false );
  
//  m_nuclideEdit->changed().connect( boost::bind( &ReferencePhotopeakDisplay::handleIsotopeChange, this, false ) );
//  m_nuclideEdit->blurred().connect( boost::bind( &ReferencePhotopeakDisplay::handleIsotopeChange, this, false ) );
//  m_nuclideEdit->enterPressed().connect( boost::bind( &ReferencePhotopeakDisplay::handleIsotopeChange, this, false ) );
  m_nuclideEdit->changed().connect( boost::bind( &ReferencePhotopeakDisplay::handleIsotopeChange, this, false ) );
  
//  m_nuclideEdit->selected().connect( boost::bind( &ReferencePhotopeakDisplay::handleIsotopeChange, this, false ) );
  
  inputLayout->addWidget( nucInputLabel, 0, 0, AlignMiddle );
  inputLayout->addWidget( m_nuclideEdit, 0, 1 );
  
  Wt::JSlot *hotKeySlot = m_spectrumViewer->hotkeyJsSlot();
  if( hotKeySlot )
    m_nuclideEdit->keyWentDown().connect( "function(s1,e1)"+ hotKeySlot->execJs( "s1", "e1" ) );
  
  tooltip = "ex. <b>U235</b>, <b>235 Uranium</b>, <b>U</b> (x-rays only)"
            ", <b>Uranium</b> (x-rays), <b>U-235m</b> (meta stable state)"
            ", <b>Cs137</b>, <b>background</b>, <b>H(n,g)</b>, etc.";
  HelpSystem::attachToolTipOn( m_nuclideEdit, tooltip, showToolTipInstantly );
  
  string replacerJs, matcherJs;
  IsotopeNameFilterModel::replacerJs( replacerJs );
  IsotopeNameFilterModel::nuclideNameMatcherJs( matcherJs );
  IsotopeNameFilterModel *isoSuggestModel = new IsotopeNameFilterModel( this );
  isoSuggestModel->addCustomSuggestPossibility( "background" );
  m_nuclideSuggest = new WSuggestionPopup( matcherJs, replacerJs, this );
  m_nuclideSuggest->setJavaScriptMember("wtNoReparent", "true");
  m_nuclideSuggest->setMaximumSize( WLength::Auto, WLength(15, WLength::FontEm) );

  //WSuggestionPopup refilter popup can fail when the user types in relly fast,
  //  so here we will do a shitty hack to catch this.
  string js = INLINE_JAVASCRIPT(
  var addTryCatch = function( elid )
  {
    var dofix = function(elid)
    {
      var el = Wt.WT.getElement(elid);
      var self = el ? jQuery.data(el, 'obj') : null;
      if( !self )
      {
        //Apparently not immediately available even though m_nuclideSuggest
        //  should be in the DOM by the time this JS gets executed.
        setTimeout( function(){dofix(elid);}, 100 );
        return;
      }
      
      var oldfcn = self.refilter;
      self.refilter = function(value){ try{ oldfcn(value); }catch(e){ console.log('My refilter caught: ' + e ); } };
    };
    dofix(elid);
  };
  );
  
  m_nuclideSuggest->doJavaScript( js + " addTryCatch('" + m_nuclideSuggest->id() + "');" );
  
  
  isoSuggestModel->filter( "" );
  m_nuclideSuggest->setFilterLength( -1 );
  m_nuclideSuggest->setModel( isoSuggestModel );
  m_nuclideSuggest->filterModel().connect( isoSuggestModel, &IsotopeNameFilterModel::filter );
  m_nuclideSuggest->forEdit( m_nuclideEdit, WSuggestionPopup::Editing );  // | WSuggestionPopup::DropDownIcon


  WLabel *ageInputLabel = new WLabel( "Age:" );
  m_ageEdit = new WLineEdit( "" );
  WRegExpValidator *validator = new WRegExpValidator( PhysicalUnits::sm_timeDurationHalfLiveOptionalRegex, this );
  validator->setFlags(Wt::MatchCaseInsensitive);
  m_ageEdit->setValidator(validator);
  m_ageEdit->setAutoComplete( false );

  
  inputLayout->addWidget( ageInputLabel, 1, 0, AlignMiddle );
  inputLayout->addWidget( m_ageEdit, 1, 1 );
  m_ageEdit->changed().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  m_ageEdit->blurred().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  m_ageEdit->enterPressed().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  
  tooltip = "<div>Age can be specified using a combination of time units, "
            "similar to '<b>5.3y 8d 22m</b>' or in half lives like "
            "'<b>2.5 HL</b>'.</div>"
            "<div>"
            "Acceptible time units: <b>year</b>, <b>yr</b>, <b>y</b>, <b>day</b>, <b>d</b>, <b>hrs</b>, <b>hour</b>, <b>h</b>, <b>minute</b>, "
            "<b>min</b>, <b>m</b>, <b>second</b>, <b>s</b>, <b>ms</b>, <b>microseconds</b>, <b>us</b>, <b>nanoseconds</b>, <b>ns</b>, or "
            "you can specify time period by <b>hh:mm:ss</b>. Half life units can be "
            "specified using <b>hl</b>, <b>halflife</b>, <b>halflives</b>, <b>half-life</b>, <b>half-lives</b>, "
            "<b>half lives</b>, or <b>half life</b>."
            "</div>"
            "<div>"
           "Half life units or time periods can not be mixed with "
           "other units. When multiple time periods are "
            "specified, they are summed, e.x. '1y6months 3m' is interpreted as "
            "18 months and 3 minutes"
            "</div>";
  //NAZ I think we can shorten this tooltip by taking out: "<div>Age can be specified using a combination of time units, " "similar to '<b>5.3y 8d 22m</b>' or in half lives like " "'<b>2.5 HL</b>'.</div>" portion. It is pretty much covered in the paragraph under it - just a suggestion :)

  HelpSystem::attachToolTipOn( m_ageEdit, tooltip, showToolTipInstantly );
  
  m_promptLinesOnly = new WCheckBox( "Prompt Only" );  //ɣ
  m_promptLinesOnly->setMargin( 5, Wt::Left );

  tooltip = "Gammas from only the original nuclide, and the daughters until one"
            " of them has a longer half-life than the original nuclide; the"
            " decay chain is in equilirium till that point.";
  HelpSystem::attachToolTipOn( m_promptLinesOnly, tooltip, showToolTipInstantly );
//m_promptLinesOnly->setHiddenKeepsGeometry( true );  // causes display to function badly; why?
  m_promptLinesOnly->checked().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  m_promptLinesOnly->unChecked().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  m_promptLinesOnly->hide();
  
  //m_layout->addWidget( m_promptLinesOnly, 3, 2, AlignMiddle );
  
  WContainerWidget *hlRow = new WContainerWidget();
  inputLayout->addWidget( hlRow, 2, 0, 1, 3 );
  
  m_halflife = new WText( hlRow );
  //Hack to make the text roughly vertically align in middle due to increased color input size...
  m_halflife->setAttributeValue( "style", "display: inline-block; margin-top: 5px" );
  //m_layout->addWidget( m_halflife, 2, 0, 1, 2, AlignMiddle | AlignCenter );
  hlRow->addWidget( m_promptLinesOnly );
  //should add a slider here

  //should add prompt and 'bare' only lines options
//  label = new WLabel( "Lowest I:" );
  //WLabel *minAmpLabel = new WLabel(   "Min Amp:" );

  //tooltip = "The minimum relative gamma amplitude to display; the most intense"
  //          " gamma ray will have value 1.0.  Amplitude is calculated after"
  //          " the optional shielding and detector effects are applied.";
  //HelpSystem::attachToolTipOn( minAmpLabel, tooltip, showToolTipInstantly );
  //m_lowerBrCuttoff = new WDoubleSpinBox();
  //HelpSystem::attachToolTipOn( m_lowerBrCuttoff, tooltip, showToolTipInstantly );
  //m_lowerBrCuttoff->setValue( 0.0 );
  //m_lowerBrCuttoff->setSingleStep( 0.01 );
  //m_lowerBrCuttoff->setRange( 0.0, 1.0 );
  //m_lowerBrCuttoff->valueChanged().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );

  //m_layout->addWidget( minAmpLabel, 3, 0, AlignMiddle );
  //m_layout->addWidget( m_lowerBrCuttoff, 3, 1 );
  
  SpectraFileModel *specFileModel = specViewer->fileManager()->model();
  m_detectorDisplay = new DetectorDisplay( specViewer, specFileModel );
  
  specViewer->detectorChanged().connect( boost::bind( &ReferencePhotopeakDisplay::handleIsotopeChange, this, true ) );
  specViewer->detectorModified().connect( boost::bind( &ReferencePhotopeakDisplay::handleIsotopeChange, this, true ) );
  
  inputLayout->addWidget( m_detectorDisplay, 3, 0,  1, 3 );

  m_shieldingSelect = new ShieldingSelect( m_materialDB, NULL, m_materialSuggest, false );
  m_shieldingSelect->materialEdit()->setEmptyText( "<shielding material>" );
  m_shieldingSelect->materialChanged().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  m_shieldingSelect->materialModified().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  inputLayout->addWidget( m_shieldingSelect, 4, 0, 1, 3 );


  //m_persistLines = new WPushButton( "Persist" );
  m_persistLines = new WPushButton( "Add Another" );
  tooltip = "Keep the currently displayed lines and add a new nuclide/source to display.";
  HelpSystem::attachToolTipOn( m_persistLines, tooltip, showToolTipInstantly );

  m_persistLines->clicked().connect( this, &ReferencePhotopeakDisplay::persistCurentLines );
  inputLayout->addWidget( m_persistLines, 0, 2 );
  m_persistLines->disable();

  //Well make the "Clear All" button a little bit wider for phones so that
  //  "Gammas" will be on same line as its check box, a little bit hacky
  if( specViewer->isPhone() )
    m_clearLines = new WPushButton( "Remove" );
  else
    m_clearLines = new WPushButton( "Clear" );
  m_clearLines->disable();
  
  if( specViewer->isMobile() )
  {
    //Get rid of the software keyboard on mobile devices.  It would be nice to
    //  simulate a screen tap on android devices as well, to get rid of system
    //  navigation UI
    m_nuclideEdit->enterPressed().connect( boost::bind( &WPushButton::setFocus, m_clearLines, true ) );
  }
  
  tooltip = "Clears all persisted lines, as well as the current non-persisted"
            " lines.";
  HelpSystem::attachToolTipOn( m_clearLines, tooltip, showToolTipInstantly );
  m_clearLines->clicked().connect( this, &ReferencePhotopeakDisplay::clearAllLines );
  inputLayout->addWidget( m_clearLines, 1, 2 );

  //m_fitPeaks = new WPushButton( "Fit Peaks" );
  //tooltip = "Fits dominant peaks for primary nuclide";
  //HelpSystem::attachToolTipOn( m_fitPeaks, tooltip, showToolTipInstantly );
  //m_fitPeaks->clicked().connect( this, &ReferencePhotopeakDisplay::fitPeaks );
  //inputLayout->addWidget( m_fitPeaks, 2, 2 );
  //m_fitPeaks->disable();
  
  if( m_lineColors.empty() )
    m_lineColors = ns_def_line_colors;
    
  m_colorSelect = new ColorSelect(ColorSelect::PrefferNative);
  m_colorSelect->setColor( m_lineColors[0] );
  
  if( ColorSelect::willUseNativeColorPicker() )
  {
    m_colorSelect->setFloatSide( Wt::Right );
    hlRow->addWidget( m_colorSelect );
  }else
  {
    WContainerWidget *w = new WContainerWidget();
    w->addWidget( m_colorSelect );
    hlRow->addWidget( w );
    w->setFloatSide( Wt::Right );
  }
  
  m_colorSelect->cssColorChanged().connect( boost::bind( &ReferencePhotopeakDisplay::userColorSelectCallback, this, _1 ) );
  m_currentlyShowingNuclide.lineColor = m_lineColors[0];
  
  WContainerWidget *whatToShow = new WContainerWidget();
  m_showGammas = new WCheckBox( "gammas", whatToShow );
  m_showXrays = new WCheckBox( "x-rays", whatToShow );
  m_showAlphas = new WCheckBox( "alphas", whatToShow );
  m_showBetas = new WCheckBox( "betas", whatToShow );
  HelpSystem::attachToolTipOn( whatToShow, "If checked, selection will be shown.  Gammas and "
                          "x-rays are shown in the table and on the chart, "
                          "alphas and betas only in the table.",
                          showToolTipInstantly );
  whatToShow->hide(); //(20190408): Ehh, not sure these are actually useful, btu not sure, so keeping around, just hidden for a while, and then we can delete them
  
  m_showGammas->setChecked();
  m_showXrays->setChecked();
  m_showGammas->checked().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  m_showGammas->unChecked().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  m_showXrays->checked().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  m_showXrays->unChecked().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  m_showAlphas->checked().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  m_showAlphas->unChecked().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  m_showBetas->checked().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  m_showBetas->unChecked().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  
  inputLayout->addWidget( whatToShow, 5, 0, 1, 3 );

  
  m_showAlphas->checked().connect( std::bind([](){
    passMessage( "Alphas are only be shown in the table, not on the spectrum.",
                "", WarningWidget::WarningMsgInfo );
  }) );
  m_showBetas->checked().connect( std::bind([](){
    passMessage( "Betas are only be shown in the table, not on the spectrum.",
                "", WarningWidget::WarningMsgInfo );
  }) );
  
  
  m_particleView = new RowStretchTreeView();
  
  m_particleView->setRootIsDecorated(	false ); //makes the tree look like a table! :)
  
  m_particleView->addStyleClass( "ParticleViewTable" );
  
  m_particleModel = new DecayParticleModel( this );
  m_particleView->setModel( m_particleModel );
  m_particleView->setAlternatingRowColors( true );
  m_particleView->setSortingEnabled( true );
  m_particleView->setColumnWidth( DecayParticleModel::kEnergy,
                                  WLength(6.5,WLength::FontEm) );
  m_particleView->setColumnWidth( DecayParticleModel::kBranchingRatio,
                                  WLength(6,WLength::FontEm) );
  m_particleView->setColumnWidth( DecayParticleModel::kResponsibleNuc,
                                 WLength(5,WLength::FontEm) );
  m_particleView->setColumnWidth( DecayParticleModel::kDecayMode,
                                 WLength(4,WLength::FontEm) );
  m_particleView->setColumnWidth( DecayParticleModel::kParticleType,
                                 WLength(5,WLength::FontEm) );
  
  if( !isPhone )
    setOverflow( WContainerWidget::OverflowAuto, Wt::Vertical );

  inputLayout->addWidget( new WContainerWidget(), 6, 0, 1, 3 );
  inputLayout->setRowStretch( 6, 1 );
    
  WGridLayout *overallLayout = new WGridLayout();
  overallLayout->setContentsMargins( 0, 0, 0, 0 );
  setLayout( overallLayout );
    
  overallLayout->addWidget( inputDiv, 0, 0 );
  overallLayout->addWidget( m_particleView, 0, 1 );
  overallLayout->setColumnStretch( 1, 1 );
}//ReferencePhotopeakDisplay constructor


ReferencePhotopeakDisplay::~ReferencePhotopeakDisplay()
{
  //I think the DOM root should take care of deleting m_nuclideSuggest
  if( m_nuclideSuggest )
    delete m_nuclideSuggest;
}//~ReferencePhotopeakDisplay()


const ReferenceLineInfo &ReferencePhotopeakDisplay::currentlyShowingNuclide() const
{
  return m_currentlyShowingNuclide;
}


const std::vector<ReferenceLineInfo> &
                                 ReferencePhotopeakDisplay::persistedNuclides() const
{
  return m_persisted;
}

std::vector<ReferenceLineInfo> ReferencePhotopeakDisplay::showingNuclides() const
{
  std::vector<ReferenceLineInfo> answer;
  
  if( !m_currentlyShowingNuclide.labelTxt.empty() )
    answer.push_back( m_currentlyShowingNuclide );
  answer.insert( answer.end(), m_persisted.begin(), m_persisted.end() );
  
  return answer;
}//std::vector<ReferenceLineInfo> showingNuclides() const;

void ReferencePhotopeakDisplay::setFocusToIsotopeEdit()
{
  m_nuclideEdit->setFocus();
  
#if( WT_VERSION >= 0x3030400 )
  InterSpecApp *app = dynamic_cast<InterSpecApp *>( wApp );
  const bool isMobile = (app && app->isMobile());
  const int nchar = static_cast<int>(m_nuclideEdit->text().narrow().length());
  if( nchar && !isMobile )
    m_nuclideEdit->setSelection( 0, nchar );
#endif
}//void setFocusToIsotopeEdit()


void ReferencePhotopeakDisplay::handleIsotopeChange( const bool useCurrentAge )
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const string isotopeLabel = m_nuclideEdit->text().toUTF8();
  const SandiaDecay::Nuclide *nuc = db->nuclide( isotopeLabel );
  const string agestr = m_ageEdit->text().toUTF8();
  const SandiaDecay::Element *el = NULL;

  if( nuc )
  {
    bool prompt = (nuc->canObtainPromptEquilibrium()
                   && m_promptLinesOnly->isChecked());
    if( (prompt || nuc->decaysToStableChildren()) && m_ageEdit->isEnabled() )
      m_ageEdit->disable();
    else if( !m_ageEdit->isEnabled() )
      m_ageEdit->enable();
  }//if( nuc )
  
  try
  {
  if( nuc && !IsInf(nuc->halfLife) && !nuc->decaysToChildren.empty() && !useCurrentAge )
  {
    
    if( nuc->decaysToStableChildren() )
    {
      m_ageEdit->setText( "0y" );
    }else if( nuc->canObtainPromptEquilibrium() && m_promptLinesOnly->isChecked() )
    {
      WString hlstr = PhysicalUnits::printToBestTimeUnits(
                                          5.0*nuc->promptEquilibriumHalfLife(),
                                          2, SandiaDecay::second );
      m_ageEdit->setText( hlstr );
    }else if( m_currentlyShowingNuclide.nuclide != nuc )
    {
//      const double age = PeakDef::defaultDecayTime( nuc );
//      WString agestr = PhysicalUnits::printToBestTimeUnits( age, 2,
//                                                          SandiaDecay::second );
      
      string agestr;
      PeakDef::defaultDecayTime( nuc, &agestr );
      m_ageEdit->setText( agestr );
    }else
    {
      const double hl = (nuc ? nuc->halfLife : -1.0);
      double age = PhysicalUnits::stringToTimeDurationPossibleHalfLife( agestr, hl );
      if( age > 100.0*nuc->halfLife || age < 0.0 )
        throw std::runtime_error( "" );
    }//if( nuc->decaysToStableChildren() ) / else
  }else if( nuc && useCurrentAge )
  {
    double age = PhysicalUnits::stringToTimeDurationPossibleHalfLife( agestr, nuc->halfLife );
    if( age > 100.0*nuc->halfLife || age < 0.0 )
      throw std::runtime_error( "" );
  }else if( nuc )
  {
    if( IsInf(nuc->halfLife) )
    {
      passMessage( isotopeLabel + " is stable",
                   "", WarningWidget::WarningMsgHigh );
    }else
    {
      passMessage( isotopeLabel + " is missing decay data",
                   "", WarningWidget::WarningMsgHigh );
    }
    
    m_nuclideEdit->setText( "" );
  }else
  {
    if( isotopeLabel.find_first_of( "0123456789" ) == string::npos )
      el = db->element( isotopeLabel );

    if( !el && !isotopeLabel.empty() )
    {
      std::vector<ReactionGamma::ReactionPhotopeak> reactions;
      
      try
      {
        const ReactionGamma *rctnDb = ReactionGammaServer::database();
        if( rctnDb )
          rctnDb->gammas( isotopeLabel, reactions );
      }catch(...)
      {}
      
      if( isotopeLabel.size() && reactions.empty()
          && !UtilityFunctions::icontains( isotopeLabel, "background") )
      {
        passMessage( isotopeLabel +
                     " is not a valid isotope, element, or reaction", "",
                     WarningWidget::WarningMsgHigh );
      }
    }
  }//if( nuc && !IsInf(nuc->halfLife) ) / else
  }catch(...)
  {
    if( nuc )
    {
      string agestr;
      PeakDef::defaultDecayTime( nuc, &agestr );
      passMessage( "Changed age to a more reasonable value for " + nuc->symbol,
                   " from '" + agestr + "' to '" + agestr + "'",
                   WarningWidget::WarningMsgLow );
      m_ageEdit->setText( agestr );
    }else
    {
      m_ageEdit->setText( "0y" );
    }//if( nuc ) / else
  }//try / catch
  
  updateDisplayChange();
  
  if( nuc || el )
    m_displayingNuclide.emit();
  else if( m_persisted.empty() )
    m_nuclidesCleared.emit();
}//void handleIsotopeChange();


std::map<std::string,std::vector<Wt::WColor>> ReferencePhotopeakDisplay::currentlyUsedPeakColors()
{
  std::map<string,vector<WColor>> answer;
  
  shared_ptr<const deque<PeakModel::PeakShrdPtr>> peaks;
  PeakModel *peakModel = m_spectrumViewer->peakModel();
  if( peakModel )
    peaks = peakModel->peaks();
  
  if( peaks )
  {
    for( const auto &p : *peaks )
    {
      string src;
      const WColor &color = p->lineColor();
      if( color.isDefault() )
        continue;
      if( p->parentNuclide() )
        src = p->parentNuclide()->symbol;
      else if( p->xrayElement() )
        src = p->xrayElement()->symbol;
      else if( p->reaction() )
        src = p->reaction()->name();
      else
        continue;
      vector<WColor> &colors = answer[src];
      if( std::find(begin(colors), end(colors), color) == end(colors) )
        colors.push_back( color );
    }//for( const auto &p : *peaks )
  }//if( peaks )

  for( const auto &p : m_persisted )
  {
    if( p.lineColor.isDefault() )
      continue;
    
    vector<WColor> &colors = answer[p.labelTxt];
    if( std::find(begin(colors), end(colors), p.lineColor) == end(colors) )
      colors.push_back( p.lineColor );
  }//for( const auto &p : m_persisted )
  
  return answer;
}//currentlyUsedPeakColors()


void ReferencePhotopeakDisplay::updateDisplayChange()
{
  bool show = true;
  show = (show && (!m_lowerBrCuttoff || m_lowerBrCuttoff->validate()==WValidator::Valid));

  const string isotxt = m_nuclideEdit->text().toUTF8();
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide *nuc = db->nuclide( isotxt );

  const bool isSameSrc = (m_currentlyShowingNuclide.labelTxt == isotxt);
  
  const SandiaDecay::Element *el = NULL;
  if( isotxt.find_first_of( "0123456789" ) == string::npos )
    el = db->element( isotxt );

  string reactions;
  vector<ReactionGamma::ReactionPhotopeak> rctnGammas;
  const string::size_type open_paren = isotxt.find( "(" );
  if( open_paren != string::npos )
  {
    try
    {
      const string::size_type close_paren = isotxt.find( ")", open_paren );
      if( close_paren == string::npos )
        throw runtime_error( "No closing paren" );

      const ReactionGamma *rctnDb = ReactionGammaServer::database();
      if( rctnDb )
      {
        reactions = rctnDb->gammas( isotxt, rctnGammas );
        nuc = NULL;
        el = NULL;
        //XXX - should use regex below to properly escape Fe(n,n')
        UtilityFunctions::ireplace_all( reactions, "'", "" );
//        UtilityFunctions::replace_all( reactions, "'", "\'" );
      }
    }catch( std::exception &e )
    {
      cerr << "ReferencePhotopeakDisplay::updateDisplayChange(): " <<e.what()<< endl;
    }//try / catch
  }//if( isotxt.find("(") != string::npos )
  
  bool isBackground = UtilityFunctions::icontains( isotxt, "background" );
  
  const bool showGammaChecked = (!m_showGammas || m_showGammas->isChecked());
  const bool showXrayChecked = (!m_showXrays || m_showXrays->isChecked());
  
  show = (show && (nuc || (el && showXrayChecked) || isBackground
                   || (!rctnGammas.empty() && showGammaChecked)) );

  double age = -1.0;
  bool canHavePromptEquil = false;

  if( nuc )
  {
    canHavePromptEquil = nuc->canObtainPromptEquilibrium();
    if( canHavePromptEquil == m_promptLinesOnly->isHidden() )
      m_promptLinesOnly->setHidden( !canHavePromptEquil );

    if( !canHavePromptEquil )
      m_promptLinesOnly->setUnChecked();

    if( canHavePromptEquil && m_promptLinesOnly->isChecked() )
    {
      //XXX - this next line doesnt have any effect, since the edit will be
      //      disabled anyway, or something...
      m_ageEdit->setText( "" );
      if( m_ageEdit->isEnabled() )
        m_ageEdit->disable();
    }else if( m_ageEdit->isDisabled() )
    {
      m_ageEdit->enable();
    }
    
    WString hlstr = PhysicalUnits::printToBestTimeUnits( nuc->halfLife,
                                                      2, SandiaDecay::second );
    //hlstr = L" \x03BB=" + hlstr;
    hlstr = "<span style=\"font-size: small;\">&lambda;<sub>&frac12;</sub>=" + hlstr + "</span>";
    m_halflife->setText( hlstr );
  }else
  {
    m_halflife->setText( "" );
  }//if( nuc ) / else

  if( nuc )
  {
    try
    {
      const string agestr = m_ageEdit->text().toUTF8();
      const double hl = (nuc ? nuc->halfLife : -1.0);
      age = PhysicalUnits::stringToTimeDurationPossibleHalfLife( agestr, hl );
      
      if( age > 100.0*nuc->halfLife || age < 0.0 )
      {
        string agestr;
        age = PeakDef::defaultDecayTime( nuc, &agestr );
        passMessage("Changed age to a more reasonable value for " + nuc->symbol,
                     " from '" + agestr + "' to " + agestr,
                     WarningWidget::WarningMsgLow );
        m_ageEdit->setText( agestr );
      }
    }catch(...)
    {
      if( m_ageEdit->text().toUTF8() == "" )
      {
        m_ageEdit->setText( "5 HL" );
        age = 5.0 * nuc->halfLife;
      }else
      {
        show = false;
      }
    }//try /catch to get the age
  }//if( !el )
  
  
  //m_fitPeaks->setDisabled( !show );
  
  //const string oldLabel = m_currentlyShowingNuclide.labelTxt;
  const bool userPickedOld = m_userHasPickedColor;
  
  //Default to using the old color, unless we find a reason to change it.
  WColor color = m_currentlyShowingNuclide.lineColor;
  
  //Mar
  m_currentlyShowingNuclide.reset();
  //cout << "Ref line widget has isSameSrc=" << isSameSrc << ", show=" << show << endl;
  
  //What we need to do here is search through
  if( isSameSrc && show && !color.isDefault() )
  {
    //dont need to do anything to the colors here I think
  }else
  {
    if( !isSameSrc )
      m_userHasPickedColor = false;
    
    if( m_peaksGetAssignedRefLineColor )
    {
#ifndef _MSC_VER
#warning "Need to test string comparison for sources always work.  E.g., need case-insensitive, etc."
#endif    
//cout << "Peak will get assigned color from ref line" << endl;
  
      const string src = m_nuclideEdit->text().toUTF8();
      const map<string,vector<WColor>> usedColors = currentlyUsedPeakColors();
      
      auto hasBeenUsed = [&usedColors](const WColor &color)->bool{
        for( const auto &s : usedColors )
          for( const auto &c : s.second )
            if( c == color )
              return true;
        return false;
      };
      
      const auto usedIter = usedColors.find(isotxt);
      const auto previter = m_previouslyPickedSourceColors.find(isotxt);
      const auto specificiter = m_specificSourcelineColors.find(isotxt);
      
      if( userPickedOld && !hasBeenUsed(color)
          && (specificiter == end(m_specificSourcelineColors)) )
      {
        //If the user picked a color, but didnt fit any poeaks, or persist, the
        //  lines, and the color theme doesnt call out this current source
        //  explicitly, then use the previous color.
        //  This maybe seems a little more intuitive from the users perspective.
        
        //We also need to propogate m_userHasPickedColor==true forward for when
        //  user keeps entering different sources.
        m_userHasPickedColor = true;
      }else if( usedIter != end(usedColors) && !usedIter->second.empty() )
      {
        color = usedIter->second[0];
        //cout << "Source " << isotxt << " has color " << color.cssText() << " in the spectrum already" << endl;
      }else if( previter != end(m_previouslyPickedSourceColors)
               && !hasBeenUsed(previter->second) )
      {
        color = previter->second;
        //cout << "Source " << isotxt << " has previously had color " << color.cssText() << " picked" << endl;
      }else if( specificiter != end(m_specificSourcelineColors) )
      {
        color = specificiter->second;
        //cout << "Source " << isotxt << " has a specific color, " << color.cssText() << ", in the ColorTheme" << endl;
      }else
      {
        auto colorcopy = m_lineColors;
        for( const auto &p : usedColors )
        {
          for( const auto &c : p.second )
          {
            auto pos = std::find( begin(colorcopy), end(colorcopy), c );
            if( pos != end(colorcopy) )
              colorcopy.erase(pos);
          }
        }//for( loop over colors used for peaks already )
        
        for( const auto &p : m_persisted )
        {
          auto pos = std::find( begin(colorcopy), end(colorcopy), p.lineColor );
          if( pos != end(colorcopy) )
            colorcopy.erase(pos);
        }//for( const auto &p : m_persisted )
        
        if( colorcopy.empty() )
          color = m_lineColors[m_persisted.size() % m_lineColors.size()];
        else
          color = colorcopy[0];
        //cout << "Source " << isotxt << " will select color, " << color.cssText()
        //     << ", from default list (len=" << colorcopy.size() << " of "
        //     << m_lineColors.size() << ")" << endl;
      }//if( src has been seen ) / else (user picked color previously)
    }else
    {
      color = m_lineColors[m_persisted.size() % m_lineColors.size()];
    }
  }//if( isSameSrc ) / else

  m_currentlyShowingNuclide.lineColor = color;
  m_colorSelect->setColor( color );
  
  m_currentlyShowingNuclide.displayLines = show;
  if( show )
  {
    m_currentlyShowingNuclide.nuclide         = nuc;
    m_currentlyShowingNuclide.element         = el;
    m_currentlyShowingNuclide.reactionGammas  = rctnGammas;
    m_currentlyShowingNuclide.reactionsTxt    = reactions;
    m_currentlyShowingNuclide.isBackground    = isBackground;
    m_currentlyShowingNuclide.isReaction      = !rctnGammas.empty();
    m_currentlyShowingNuclide.age             = age;
    m_currentlyShowingNuclide.labelTxt        = m_nuclideEdit->text().toUTF8();
    m_currentlyShowingNuclide.lowerBrCuttoff  = (m_lowerBrCuttoff ? m_lowerBrCuttoff->value() : 0.0);
    m_currentlyShowingNuclide.showGammas      = (!m_showGammas || m_showGammas->isChecked());
    m_currentlyShowingNuclide.showXrays       = (!m_showXrays || m_showXrays->isChecked());
    m_currentlyShowingNuclide.showAlphas      = m_showAlphas->isChecked();
    m_currentlyShowingNuclide.showBetas       = m_showBetas->isChecked();
    m_currentlyShowingNuclide.showLines       = true;
    m_currentlyShowingNuclide.promptLinesOnly
                       = (canHavePromptEquil && m_promptLinesOnly->isChecked());
    
    m_currentlyShowingNuclide.shieldingName = "";
    m_currentlyShowingNuclide.shieldingThickness = 0.0;
    
    try
    {
      if( m_shieldingSelect->isGenericMaterial() )
      {
        const double an = m_shieldingSelect->atomicNumber();
        const double ad = m_shieldingSelect->arealDensity()
                          / (PhysicalUnits::gram/PhysicalUnits::cm2);
        char buffer[128];
        snprintf( buffer, sizeof(buffer), "%.2f, %.2f g/cm2", an, ad );
        m_currentlyShowingNuclide.shieldingName = buffer;
      }else
      {
        std::shared_ptr<Material> m = m_shieldingSelect->material();
        const double thickness = m_shieldingSelect->thickness();
        if( !!m )
        {
          m_currentlyShowingNuclide.shieldingName = m->name;
          m_currentlyShowingNuclide.shieldingThickness = thickness;
        }
      }//if( m_shieldingSelect->isGenericMaterial() )
    }catch( std::exception &e )
    {
      cerr << "Failed to get shielding defintion for reference gamma lines: "
           << e.what() << endl;
    }//try/ catch
  
    m_currentlyShowingNuclide.detectorName = "";
    std::shared_ptr<DetectorPeakResponse> det = m_detectorDisplay->detector();
    if( !!det )
      m_currentlyShowingNuclide.detectorName = det->name();
    
    
    if( isBackground )
    {
      for( const BackgroundLine &bl : BackgroundLines )
      {
        if( std::get<3>(bl)!=BackgroundXRay || (!m_showXrays || m_showXrays->isChecked()) )
          m_currentlyShowingNuclide.backgroundLines.push_back( &bl );
      }//for( const BackgroundLine &bl : BackgroundLines )
    }//if( isBackground )
  }else
  {
    //m_currentlyShowingNuclide is already reset.
  }

  if( !show )
  {
    if( m_persistLines->isEnabled() )
      m_persistLines->disable();
    m_clearLines->setDisabled( m_persisted.empty() );
    
#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )
    ReferenceLineInfo emptylines;
    m_chart->setReferncePhotoPeakLines( emptylines );
#endif
    
    m_particleModel->clear();
    
#if( !USE_SPECTRUM_CHART_D3 )
    CanvasForDragging *canvas = m_chart->overlayCanvas();
    if( canvas )
    {
      const string canvasjs = "$('#c" + canvas->id() + "')";
      string js = canvasjs + ".data('gammalines',null);"
                  "Wt.WT.DrawGammaLines('c" + canvas->id() + "',true);";
      doJavaScript( js );
    }//canvas
#endif
    
    if( age < 0.0 || !nuc )
    {
      m_particleModel->clear();
      return;
    }//if( age < 0.0 || !nuc )
  }//if( !show )

  vector<DecayParticleModel::RowData> inforows;

  const double brCutoff = (m_lowerBrCuttoff ? m_lowerBrCuttoff->value() : 0.0);

//  bool islogy = m_chart->yAxisIsLog();
//  double chartMaxSf = (islogy ? log(2.5) : 1.0/1.1);

  SandiaDecay::NuclideMixture mixture;

  if( nuc && canHavePromptEquil && m_promptLinesOnly->isChecked() )
  {
    age = 0.0;
    mixture.addNuclideInPromptEquilibrium( nuc, 1.0E-3 * SandiaDecay::curie );
  }else if( nuc )
  {
    mixture.addNuclideByActivity( nuc, 1.0E-3 * SandiaDecay::curie );
  }//if( we want promt only ) / else

  const vector<SandiaDecay::NuclideActivityPair> activities
                                                     = mixture.activity( age );
  const double parent_activity = nuc ? mixture.activity( age, nuc ) : 0.0;

  vector<double> energies, branchratios;
  vector<SandiaDecay::ProductType> particle_type;
  vector<const SandiaDecay::Transition *> transistions;
  vector<SandiaDecay::ProductType> types;
  vector<const ReactionGamma::Reaction *> reactionPeaks;
  std::vector<const BackgroundLine *> backgroundLines;

  if( !m_showGammas || m_showGammas->isChecked() )
  {
    types.push_back( SandiaDecay::GammaParticle );
    types.push_back( SandiaDecay::PositronParticle );
  }
  
  if( m_showAlphas->isChecked() )
    types.push_back( SandiaDecay::AlphaParticle );
  if( m_showBetas->isChecked() )
    types.push_back( SandiaDecay::BetaParticle );
  if( (!m_showXrays || m_showXrays->isChecked()) )
    types.push_back( SandiaDecay::XrayParticle );
  
  DecayParticleModel::RowData positronrow;
  positronrow.energy      = static_cast<float>( 510.9989 * PhysicalUnits::keV );
  positronrow.branchRatio = 0.0f;
  positronrow.particle    = SandiaDecay::GammaParticle; //SandiaDecay::positron
  positronrow.responsibleNuc = 0;
  std::set<const SandiaDecay::Nuclide *> positronparents;
  std::set<const SandiaDecay::Transition *> positrontrans;

  for( SandiaDecay::ProductType type : types )
  {
    for( size_t nucIndex = 0; nucIndex < activities.size(); ++nucIndex )
    {
      const SandiaDecay::Nuclide *nuclide = activities[nucIndex].nuclide;
      const double activity = activities[nucIndex].activity;

      const size_t n_decaysToChildren = nuclide->decaysToChildren.size();
      
      for( size_t decayIndex = 0; decayIndex < n_decaysToChildren; ++decayIndex )
      {
        const SandiaDecay::Transition *transition
                                       = nuclide->decaysToChildren[decayIndex];
        const size_t n_products = transition->products.size();

        for( size_t productNum = 0; productNum < n_products; ++productNum )
        {
          const SandiaDecay::RadParticle &particle
                                            = transition->products[productNum];
          if( type == SandiaDecay::PositronParticle && particle.type == SandiaDecay::PositronParticle )
          {
            const double br = activity * particle.intensity
                              * transition->branchRatio / parent_activity;
            positronrow.branchRatio += 2.0*br;
            positronrow.decayMode   = transition->mode;
            positronparents.insert( transition->parent );
            positrontrans.insert( transition );
          }else if( (particle.type == type) && (particle.type == SandiaDecay::XrayParticle) )
          {
            size_t index = 0;
            for( ; index < energies.size(); ++index )
            {
              if( fabs(energies[index] - particle.energy) < 1.0E-6 )
                break;
            }
            
            const double br = activity * particle.intensity
                                   * transition->branchRatio / parent_activity;
            
            if( index < energies.size() )
            {
              branchratios[index] += br;
              inforows[index].branchRatio += br;
            }else
            {
              transistions.push_back( NULL );
              energies.push_back( particle.energy );
              branchratios.push_back( br );
              particle_type.push_back( SandiaDecay::XrayParticle );
              reactionPeaks.push_back( NULL );
              backgroundLines.push_back( NULL );
            
              DecayParticleModel::RowData row;
              row.energy      = particle.energy;
              row.branchRatio = br;
              row.particle    = SandiaDecay::XrayParticle;
              row.decayMode   = DecayParticleModel::RowData::XRayDecayMode;
              row.responsibleNuc = nuc;
            
              inforows.push_back( row );
            }
          }else if( particle.type == type )
          {
            transistions.push_back( transition );
            energies.push_back( particle.energy );
            const double br = activity * particle.intensity
                                   * transition->branchRatio / parent_activity;
            branchratios.push_back( br );
            particle_type.push_back( type );
            
            reactionPeaks.push_back( NULL );
            backgroundLines.push_back( NULL );

            DecayParticleModel::RowData row;
            row.energy      = particle.energy;
            row.branchRatio = br;
            row.particle    = particle.type;
            row.decayMode   = transition->mode;
            row.responsibleNuc = transition->parent;
            inforows.push_back( row );
          }//if( particle.type == type )
        }//for( size_t productNum = 0; productNum < n_products; ++productNum )
      }//for( size_t decayIndex = 0; decayIndex < n_decaysToChildren; ++decayIndex )
    }//for( size_t nucIndex = 0; nucIndex < activities.size(); ++nucIndex )
  }//for( SandiaDecay::ProductType type : types )

  if( positronrow.branchRatio > 0.0 )
  {
    if( positronparents.size() == 1 )
      positronrow.responsibleNuc = *positronparents.begin();
    
    if( positrontrans.size() == 1 )
      transistions.push_back( *positrontrans.begin() );
    else
      transistions.push_back( NULL );
    energies.push_back( positronrow.energy );
    branchratios.push_back( positronrow.branchRatio );
    particle_type.push_back( SandiaDecay::GammaParticle );
    reactionPeaks.push_back( NULL );
    backgroundLines.push_back( NULL );
    
    inforows.push_back( positronrow );
  }//if( positronrow.branchRatio > 0.0 )
  
  const bool showXrays = ((el && !nuc) && (!m_showXrays || m_showXrays->isChecked()));
  
  if( showXrays )
  {
    const SandiaDecay::Element *element = el;
    if( !element )
      element = db->element( nuc->atomicNumber );

    for( const SandiaDecay::EnergyIntensityPair &eip : element->xrays )
    {
      transistions.push_back( NULL );
      energies.push_back( eip.energy );
      branchratios.push_back( eip.intensity );
      particle_type.push_back( SandiaDecay::XrayParticle );
      reactionPeaks.push_back( NULL );
      backgroundLines.push_back( NULL );

      DecayParticleModel::RowData row;
      row.energy      = eip.energy;
      row.branchRatio = eip.intensity;
      row.particle    = SandiaDecay::XrayParticle;
      row.decayMode   = DecayParticleModel::RowData::XRayDecayMode;
      row.responsibleNuc = nuc;

      inforows.push_back( row );
    }//for( const SandiaDecay::EnergyIntensityPair &eip : element->xrays )
  }//if( m_showXrays->isChecked() )

  
  for( const ReactionGamma::ReactionPhotopeak &eip : rctnGammas )
  {
    transistions.push_back( NULL );
    energies.push_back( eip.energy );
    branchratios.push_back( eip.abundance );
    particle_type.push_back( SandiaDecay::GammaParticle );
    reactionPeaks.push_back( eip.reaction );
    backgroundLines.push_back( NULL );

    DecayParticleModel::RowData row;
    row.energy      = eip.energy;
    row.branchRatio = eip.abundance;
    row.particle    = SandiaDecay::GammaParticle;
    row.decayMode   = DecayParticleModel::RowData::ReactionToGammaMode;
    row.responsibleNuc = nuc;

    inforows.push_back( row );
  }//for( const SandiaDecay::EnergyIntensityPair &eip : element->xrays )


  for( const BackgroundLine *bl : m_currentlyShowingNuclide.backgroundLines )
  {
    transistions.push_back( NULL );
    energies.push_back( std::get<0>(*bl) );
    branchratios.push_back( std::get<1>(*bl) );
    particle_type.push_back( SandiaDecay::GammaParticle );
    reactionPeaks.push_back( NULL );
    backgroundLines.push_back( bl );
    
    DecayParticleModel::RowData row;
    row.energy      = std::get<0>(*bl);
    row.branchRatio = std::get<1>(*bl);
    row.particle    = SandiaDecay::GammaParticle;
    
    switch( std::get<3>(*bl) )
    {
      case U238Series: case U235Series: case Th232Series: case Ra226Series:
      case K40Background: case OtherBackground:
        row.decayMode = DecayParticleModel::RowData::NormGammaDecayMode;
        break;
      case BackgroundXRay:
        row.decayMode   = DecayParticleModel::RowData::XRayDecayMode;
      break;
      case BackgroundReaction:
        row.decayMode   = DecayParticleModel::RowData::ReactionToGammaMode;
      break;
    }//switch( get<3>(*bl) )
    
    row.responsibleNuc = db->nuclide( std::get<2>(*bl) );
    
    inforows.push_back( row );
  }//for( m_currentlyShowingNuclide.backgroundLines )

//  double maxOrigbr = 0.0;
//  for( size_t i = 0; i < inforows.size(); ++i )
//    maxOrigbr = max( double(inforows[i].branchRatio), maxOrigbr );
//  for( size_t i = 0; i < inforows.size(); ++i )
//    inforows[i].branchRatio /= maxOrigbr;

  //Lets get rid of brannching ratios that are incredible close to zero
  const float abs_min_br = FLT_MIN; //1.0E-15
  vector<DecayParticleModel::RowData> inforowstouse;
  for( size_t i = 0; i < inforows.size(); ++i )
    if( inforows[i].branchRatio > abs_min_br && inforows[i].branchRatio >= brCutoff )
      inforowstouse.push_back( inforows[i] );

  m_particleModel->setRowData( inforowstouse );

  //fold in detector response
  std::shared_ptr<DetectorPeakResponse> det = m_detectorDisplay->detector();

  //Wider peaks mean not as large value of 'y' for the peaks
  if( det && det->isValid() )
  {
    for( size_t i = 0; i < branchratios.size(); ++i )
      if( (particle_type[i] == SandiaDecay::GammaParticle) || (particle_type[i] == SandiaDecay::XrayParticle) )
        branchratios[i] *= det->efficiency( energies[i], PhysicalUnits::m );
  }//if( detector )

  
  //fold in shielding here....
  try
  {
    std::shared_ptr<const Material> material;
    boost::function<double(float)> att_coef_fcn;

    if( m_shieldingSelect->isGenericMaterial() )
    {
      const float atomic_number = static_cast<float>(m_shieldingSelect->atomicNumber());
      const float areal_density = static_cast<float>(m_shieldingSelect->arealDensity());
      att_coef_fcn
          = boost::bind( &GammaInteractionCalc::transmition_coefficient_generic,
                         atomic_number, areal_density, _1 );
    }else
    {
      material = m_shieldingSelect->material();
      if( !!material )
      {
        const float thick = static_cast<float>(m_shieldingSelect->thickness());
        att_coef_fcn
          = boost::bind( &GammaInteractionCalc::transmition_coefficient_material,
                          material.get(), _1, thick );
      }//if( !!material )
    }//if( isGenericMaterial ) / else

    if( !att_coef_fcn.empty() )
    {
      for( size_t i = 0; i < branchratios.size(); ++i )
        if( (particle_type[i] == SandiaDecay::GammaParticle) || (particle_type[i] == SandiaDecay::XrayParticle) )
          branchratios[i] *= exp( -1.0 * att_coef_fcn( energies[i] ) );
    }//if( att_coef_fcn )
  }catch( MassAttenuation::ErrorLoadingDataException & )
  {
    throw runtime_error( "Failed to open gamma XS data file" );
  }catch( std::exception &e )
  {
    cerr << SRC_LOCATION << ": caught error " << e.what() << endl;
#if( PERFORM_DEVELOPER_CHECKS )
    char msg[512];
    snprintf( msg, sizeof(msg), "Error caclulating attenuation: %s", e.what() );
    log_developer_error( BOOST_CURRENT_FUNCTION, msg );
#endif
  }

  
  //Peak height is: area*(1/(sigma*sqrt(2*pi)))*exp( -(x-mean)^2 / (2*sigma^2) ),
  //  therefore peak height is proportianal to area/sigma, lets correct for this
  if( det && det->isValid() && det->hasResolutionInfo() )
  {
    const vector<double> origbr = branchratios;
    
    try
    {
      for( size_t i = 0; i < branchratios.size(); ++i )
      {
        const double sigma = det->peakResolutionSigma( energies[i] );
        if( sigma <= 0.0 )
          throw exception();
        branchratios[i] /= sigma;
      }//for( size_t i = 0; i < branchratios.size(); ++i )
    }catch(...)
    {
      branchratios = origbr;
      cerr << "Encountered a negative or zero peak width, not taking detector "
           << "resolution into account, sorry :(" << endl;
    }//try / catch
  }//if( detector->hasResolutionInfo() )

  //Some decays may not produce gammas, but do produce xrays (not verified) so
  //  we want to normalize gammas and xrays relative to the largest branching
  //  ratio gamma or xray.
  double max_gamma_xray = 0.0;

  map<SandiaDecay::ProductType,double> maxbrs;
  
  for( size_t i = 0; i < branchratios.size(); ++i )
  {
    const SandiaDecay::ProductType type = particle_type[i];
    
    if( !maxbrs.count(type) )
      maxbrs[type] = 0.0;
//    if( (transistions[i]
//         || inforows[i].decayMode==DecayParticleModel::RowData::ReactionToGammaMode
//         || backgroundLines[i]) )
      maxbrs[type] = max(maxbrs[type],branchratios[i]);

    if( type == SandiaDecay::GammaParticle || type == SandiaDecay::XrayParticle )
      max_gamma_xray = std::max( max_gamma_xray, branchratios[i] );
  }//for( size_t i = 0; i < branchratios.size(); ++i )
  
  for( size_t i = 0; i < branchratios.size(); ++i )
  {
    const double energy = energies[i];
    
    //If this is an xray caused by a decay, lets normalize its amplitude relative
    //  to the gamma amplitudes.  If we are displaying just the xrays of an
    //  element, than we will normalize them to go between zero and one.
    const bool is_gamma = (particle_type[i] == SandiaDecay::GammaParticle);
    const bool is_xray = (particle_type[i] == SandiaDecay::XrayParticle);
    const bool is_decay_xray_gamma = (nuc && (is_gamma || is_xray));
    const double br = branchratios[i] / (is_decay_xray_gamma ? max_gamma_xray : maxbrs[particle_type[i]]);
    
    const SandiaDecay::Transition *transition = transistions[i];
    const BackgroundLine *backLine = backgroundLines.at(i);
    
    if( (transition || backLine) && (br <= brCutoff || IsInf(br) || IsNan(br)) )
      continue;
    
    string particlestr, decaystr, elstr;
    if( !is_xray )
    {
      const SandiaDecay::ProductType parttype
                               = SandiaDecay::ProductType( particle_type[i] );
      particlestr = SandiaDecay::to_str( parttype );
      
      if( transition )
      {
        if( transition->parent )
          decaystr = transition->parent->symbol;
        if( transition->child )
          decaystr += " to " + transition->child->symbol;
        decaystr += string(" via ") + SandiaDecay::to_str(transition->mode);
      }else if( reactionPeaks.at(i) )
      {
        decaystr = reactionPeaks[i]->name();
      }else if( backLine )
      {
        const string &symbol = std::get<2>(*backLine);
        if( symbol.size() )
          decaystr += symbol + ", ";
        
        switch( std::get<3>(*backLine) )
        {
          case U238Series:    decaystr += "U238 series";         break;
          case U235Series:    decaystr += "U235 series";         break;
          case Th232Series:   decaystr += "Th232 series";        break;
          case Ra226Series:   decaystr += "U238 (Ra226) series"; break;
          case K40Background: decaystr += "Primordial";          break;
          case OtherBackground: case BackgroundXRay: case BackgroundReaction:
            decaystr += std::get<4>(*backLine);    break;
        }//switch( get<3>(backgroundLines[i]) )
      }//if( transition ) / else ...
    }else
    {
      const SandiaDecay::Element *element = el;
      if( !element )
        element = db->element( nuc->atomicNumber );
      
      particlestr = "xray";
      decaystr = "xray";
      if( element )
        elstr = element->name;
      else if( backLine )
        elstr = std::get<4>(*backLine);
    }//if( xray ) / else
    
    m_currentlyShowingNuclide.energies.push_back(     energy );
    m_currentlyShowingNuclide.intensities.push_back(  br );
    m_currentlyShowingNuclide.particlestrs.push_back( particlestr );
    m_currentlyShowingNuclide.decaystrs.push_back(    decaystr );
    m_currentlyShowingNuclide.elementstrs.push_back(  elstr );
  }//for( size_t i = 0; i < branchratios.size(); ++i )
  
  typedef std::map<SandiaDecay::ProductType,double>::const_iterator MaxBrIter;
  
  for( MaxBrIter iter = maxbrs.begin(); iter != maxbrs.end(); ++iter )
  {
    const char *typestr = SandiaDecay::to_str( SandiaDecay::ProductType(iter->first) );
    m_currentlyShowingNuclide.particle_sf[typestr] = iter->second;
    
    if( iter->first == SandiaDecay::GammaParticle || iter->first == SandiaDecay::XrayParticle )
      m_currentlyShowingNuclide.particle_sf[typestr] = max_gamma_xray;
  }
  
  
#if( PERFORM_DEVELOPER_CHECKS )
  for( const string &s : m_currentlyShowingNuclide.particlestrs )
  {
    if( !m_currentlyShowingNuclide.particle_sf.count(s) )
    {
      char msg[512];
      snprintf( msg, sizeof(msg), "Missing particlestrs (%s) in particle_sf.", s.c_str() );
//      throw runtime_error(msg);
      log_developer_error( BOOST_CURRENT_FUNCTION, msg );
    }
  }
#endif
  
  
  
  //Clientside javascript currently doesnt know about this garuntee that gamma
  //  lines will be sorted by energy.
  m_currentlyShowingNuclide.sortByEnergy();
  
  //Also, we could play some tricks to eliminate some of the gamma lines that
  //  are so small in amplitude, they would never impact the user
  
#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )
  if( show )
    m_chart->setReferncePhotoPeakLines( m_currentlyShowingNuclide );
  else
    m_chart->setReferncePhotoPeakLines( ReferenceLineInfo() );
#endif
  
  
  string json;
#if( !USE_SPECTRUM_CHART_D3 )
  CanvasForDragging *canvas = m_chart->overlayCanvas();
  if( canvas )
  {
    const string canvasjs = "$('#c" + canvas->id() + "')";
    json = canvasjs + ".data('gammalines',";
    m_currentlyShowingNuclide.toJson( json );
    json += ");";
    json += "Wt.WT.DrawGammaLines('c" + canvas->id() + "',true);";
  }
#endif

  if( show )
  {
    if( m_persistLines->isDisabled() )
      m_persistLines->enable();
    if( m_clearLines->isDisabled() )
      m_clearLines->enable();

    if( m_spectrumViewer && m_spectrumViewer->isPhone() )
      m_clearLines->setText( m_persisted.empty() ? "Remove" : "Remove All" );
    else
      m_clearLines->setText( m_persisted.empty() ? "Clear" : "Clear All" );
    
    
    if( json.size() )
      doJavaScript( json );
  }//if( show )
}//void updateDisplayChange()



void ReferencePhotopeakDisplay::persistCurentLines()
{
  if( m_currentlyShowingNuclide.energies.empty() )
    return;

#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )
  m_chart->persistCurrentReferncePhotoPeakLines();
#endif

#if( !USE_SPECTRUM_CHART_D3 )
  CanvasForDragging *canvas = m_chart->overlayCanvas();
  if( canvas )
  {
    const string canvasjs = "$('#c" + canvas->id() + "')";
  
    WStringStream js;
    js << "var a=" << canvasjs << ".data('gammalines');"
       << "var b=" << canvasjs << ".data('persistedLines');"
       << "if(a && !b)b=[];"
       << "if(a)b[a.parent]=a;"
       << canvasjs << ".data('persistedLines',b);"
       << canvasjs << ".data('gammalines',null);";
    
    doJavaScript( js.str() );
  }//if( canvas )
#endif
  
  vector<ReferenceLineInfo>::iterator pos;
  pos = std::find( m_persisted.begin(), m_persisted.end(),
                   m_currentlyShowingNuclide );

  if( pos == m_persisted.end() && m_currentlyShowingNuclide.displayLines )
    m_persisted.push_back( m_currentlyShowingNuclide );

  m_currentlyShowingNuclide.reset();
  m_userHasPickedColor = false;
  
  m_nuclideEdit->setText( "" );
  m_persistLines->disable();
  m_clearLines->enable();
  if( m_spectrumViewer && m_spectrumViewer->isPhone() )
    m_clearLines->setText( "Remove All" );
  else
    m_clearLines->setText( "Clear All" );
  
  updateDisplayChange();
}//void persistCurentLines()


void ReferencePhotopeakDisplay::setColors( const std::vector<Wt::WColor> &referenceLineColor )
{
  m_lineColors.clear();
  m_lineColors.reserve( referenceLineColor.size() );
  for( const auto &i : referenceLineColor )
    if( !i.isDefault() )
      m_lineColors.push_back( i );
  
  if( m_lineColors.empty() )
    m_lineColors = ns_def_line_colors;
  
  //should update currently displayed line colors here, but there are bigger fish to fry ATM
}//void setColors( const std::vector<Wt::WColor> &referenceLineColor )


void ReferencePhotopeakDisplay::setColorsForSpecificSources( const std::map<std::string,Wt::WColor> &referenceLineColorForSources )
{
  m_specificSourcelineColors.clear();
  
  for( const auto &i : referenceLineColorForSources )
  {
    if( !i.second.isDefault() )
    m_specificSourcelineColors[i.first] = i.second;
  }
  
  //should update currently displayed line colors here, but there are bigger fish to fry ATM
}//void setColorsForSpecificSources( const std::map<std::string,Wt::WColor> &referenceLineColorForSources )


void ReferencePhotopeakDisplay::setPeaksGetAssignedRefLineColor( const bool theydo )
{
  m_peaksGetAssignedRefLineColor = theydo;
}


Wt::Signal<> &ReferencePhotopeakDisplay::displayingNuclide()
{
  return m_displayingNuclide;
}


Wt::Signal<> &ReferencePhotopeakDisplay::nuclidesCleared()
{
  return m_nuclidesCleared;
}

void ReferencePhotopeakDisplay::serialize( std::string &xml_data  ) const
{
  rapidxml::xml_document<char> doc;
  serialize( &doc );
  xml_data.clear();
  rapidxml::print(std::back_inserter(xml_data), doc, 0);
}//void serialize( std::string &xml_data  ) const


void ReferencePhotopeakDisplay::serialize(
                                  rapidxml::xml_node<char> *parent_node ) const
{
  rapidxml::xml_document<char> *doc = parent_node->document();
  
  const char *name, *value;
  rapidxml::xml_node<char> *base_node, *node, *element;
  rapidxml::xml_attribute<char> *attr;
  
  name = "ReferencePhotopeakDisplay";
  base_node = doc->allocate_node( rapidxml::node_element, name );
  parent_node->append_node( base_node );
  
  //If you change the available options or formatting or whatever, increment the
  //  version field of the XML!
  value = doc->allocate_string( std::to_string(sm_xmlSerializationVersion).c_str() );
  attr = doc->allocate_attribute( "version", value );
  base_node->append_attribute( attr );

  if( m_currentlyShowingNuclide.empty() )
  {
    node = doc->allocate_node( rapidxml::node_element, "CurrentGuiState" );
    base_node->append_node( node );

    name = "Nuclide";
    value = doc->allocate_string( m_nuclideEdit->text().toUTF8().c_str() );
    element = doc->allocate_node( rapidxml::node_element, name, value );
    node->append_node( element );
    
    name = "Age";
    value = doc->allocate_string( m_ageEdit->text().toUTF8().c_str() );
    element = doc->allocate_node( rapidxml::node_element, name, value );
    node->append_node( element );

    name = "LowestBranchRatio";
    value = doc->allocate_string( (m_lowerBrCuttoff ? m_lowerBrCuttoff->text().toUTF8().c_str() : "0.0") );
    element = doc->allocate_node( rapidxml::node_element, name, value );
    node->append_node( element );
    
    name = "PromptLinesOnly";
    value = (m_promptLinesOnly->isChecked() ? "1" : "0");
    element = doc->allocate_node( rapidxml::node_element, name, value );
    node->append_node( element );
    
    name = "ShowGammas";
    value = ((!m_showGammas || m_showGammas->isChecked()) ? "1" : "0");
    element = doc->allocate_node( rapidxml::node_element, name, value );
    node->append_node( element );

    name = "ShowXrays";
    value = ((!m_showXrays || m_showXrays->isChecked()) ? "1" : "0");
    element = doc->allocate_node( rapidxml::node_element, name, value );
    node->append_node( element );
    
    name = "ShowAlphas";
    value = (m_showAlphas->isChecked() ? "1" : "0");
    element = doc->allocate_node( rapidxml::node_element, name, value );
    node->append_node( element );
    
    name = "ShowBetas";
    value = (m_showBetas->isChecked() ? "1" : "0");
    element = doc->allocate_node( rapidxml::node_element, name, value );
    node->append_node( element );
  }else
  {
    node = doc->allocate_node( rapidxml::node_element, "CurrentLines" );
    base_node->append_node( node );
    m_currentlyShowingNuclide.serialize( node );
  }//if( !m_currentlyShowingNuclide.empty() )
  
  if( !m_persisted.empty() )
  {
    node = doc->allocate_node( rapidxml::node_element, "PersistedLines" );
    base_node->append_node( node );
    for( const ReferenceLineInfo &n : m_persisted )
      n.serialize( node );
  }//if( !m_persisted.empty() )
  
  WColor color;
  if( m_currentlyShowingNuclide.empty() )
    color = m_colorSelect->color();
  else
    color = m_currentlyShowingNuclide.lineColor;
  
  if( !color.isDefault() )
  {
    name = "CurrentLineColor";
    value = doc->allocate_string( color.cssText(false).c_str() );
    element = doc->allocate_node( rapidxml::node_element, name, value );
    attr = doc->allocate_attribute( "UserSelected", (m_userHasPickedColor ? "1" : "0") );
    element->append_attribute( attr );
    base_node->append_node( element );
  }//if( !color.isDefault() )

  //I guess we dont actually need to serialize m_lineColors,
  //  m_peaksGetAssignedRefLineColor, or m_specificSourcelineColors.size as
  //  these should be captured in the ColorTheme.
  /*
  rapidxml::xml_node<char> *colorsToUse = doc->allocate_node( rapidxml::node_element, "ColorsToUse" );
  attr = doc->allocate_attribute( "PeaksGetAssignedRefLineColor", (m_peaksGetAssignedRefLineColor ? "1" : "0") );
  colorsToUse->append_attribute( attr );
  base_node->append_node( colorsToUse );
  
  for( const auto &c : m_lineColors )
  {
    value = doc->allocate_string( c.cssText(false).c_str() );
    element = doc->allocate_node( rapidxml::node_element, "Color", value );
    colorsToUse->append_node( element );
  }
  
  if( m_specificSourcelineColors.size() )
  {
    rapidxml::xml_node<char> *specificSrcColor = doc->allocate_node( rapidxml::node_element, "SpecificSrcColors" );
    base_node->append_node( specificSrcColor );
    
    for( const auto &p : m_specificSourcelineColors )
    {
      auto src = doc->allocate_string( p.first.c_str() );
      auto color = doc->allocate_string( p.second.cssText(false).c_str() );
      element = doc->allocate_node( rapidxml::node_element, "SrcColor" );
      attr = doc->allocate_attribute( "src", src );
      element->append_attribute( attr );
      attr = doc->allocate_attribute( "color", color );
      element->append_attribute( attr );
      specificSrcColor->append_node( element );
    }
  }//if( m_specificSourcelineColors.size() )
  
  if( m_previouslyPickedSourceColors.size() )
  {
    rapidxml::xml_node<char> *prevSrcColor = doc->allocate_node( rapidxml::node_element, "PreviousSrcColors" );
    base_node->append_node( prevSrcColor );
    
    for( const auto &p : m_previouslyPickedSourceColors )
    {
      auto src = doc->allocate_string( p.first.c_str() );
      auto color = doc->allocate_string( p.second.cssText(false).c_str() );
      element = doc->allocate_node( rapidxml::node_element, "SrcColor" );
      attr = doc->allocate_attribute( "src", src );
      element->append_attribute( attr );
      attr = doc->allocate_attribute( "color", color );
      element->append_attribute( attr );
      prevSrcColor->append_node( element );
    }
  }//if( m_previouslyPickedSourceColors.size() )
  */
  
  m_shieldingSelect->serialize( base_node );
}//void serialize( rapidxml::xml_document<char> &doc )


void ReferencePhotopeakDisplay::deSerialize( std::string &xml_data  )
{
  clearAllLines();
  
  try
  {
    rapidxml::xml_document<char> doc;
    const int flags = rapidxml::parse_normalize_whitespace
                      | rapidxml::parse_trim_whitespace;

    //cout << "xml_data=" << xml_data << endl;
    
    if( xml_data.size() )
      doc.parse<flags>( &(xml_data[0]) );
    
    rapidxml::xml_attribute<char> *attr;
    rapidxml::xml_node<char> *base_node, *node, *showing_node,
                             *gui_node, *persisted_node;
    
    base_node = doc.first_node( "ReferencePhotopeakDisplay", 25 );
    if( !base_node )
      base_node = doc.first_node( "PhotopeakLineDisplay", 20 );
    
    if( !base_node )
      throw runtime_error( "Couldnt get base node, ReferencePhotopeakDisplay or PhotopeakLineDisplay" );
    
    int version = 0;
    attr = base_node->first_attribute( "version", 7 );
    if( !attr || !attr->value()
        || !(stringstream(attr->value()) >> version)
        || (version != sm_xmlSerializationVersion) )
      throw runtime_error( "Mising or invalid ReferencePhotopeakDisplay version" );
    
    gui_node = base_node->first_node( "CurrentGuiState", 15 );
    showing_node = base_node->first_node( "CurrentLines", 12 );
    
    if( (gui_node && showing_node) || !(gui_node || showing_node) )
      throw runtime_error( "Inconsistent saving of photpeak lines state" );
    
    if( showing_node && showing_node->first_node("DisplayedSource",15) )
      gui_node = showing_node->first_node("DisplayedSource",15);
    
    if( gui_node )
    {
      const SandiaDecay::Nuclide *nuc = nullptr;
      node = gui_node->first_node( "Nuclide", 7 );
      if( node && node->value() )
      {
        const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
        if( db )
          nuc = db->nuclide( node->value() );
        m_nuclideEdit->setText( node->value() );
      }
      
      node = gui_node->first_node( "Age", 3 );
      if( node && node->value() && nuc )
        m_ageEdit->setText( node->value() );
      else
        m_ageEdit->setText( "" );
      
      node = gui_node->first_node( "LowestBranchRatio", 17 );
      if( node && node->value() && m_lowerBrCuttoff )
        m_lowerBrCuttoff->setText( node->value() );
      
      node = gui_node->first_node( "PromptLinesOnly", 15 );
      if( node && node->value() && strlen(node->value()))
        m_promptLinesOnly->setChecked( (node->value()[0] == '1') );
      
      node = gui_node->first_node( "ShowGammas", 10 );
      if( node && node->value() && strlen(node->value()) && m_showGammas )
        m_showGammas->setChecked( (node->value()[0] == '1') );
      
      node = gui_node->first_node( "ShowXrays", 9 );
      if( node && node->value() && strlen(node->value()) && m_showXrays )
        m_showXrays->setChecked( (node->value()[0] == '1') );
      
      node = gui_node->first_node( "ShowAlphas", 10 );
      if( node && node->value() && strlen(node->value()))
        m_showAlphas->setChecked( (node->value()[0] == '1') );
      
      node = gui_node->first_node( "ShowBetas", 9 );
      if( node && node->value() && strlen(node->value()))
        m_showBetas->setChecked( (node->value()[0] == '1') );
    }//if( gui_node )
    
    if( showing_node )
    {
      m_currentlyShowingNuclide.reset();
      node = showing_node->first_node( "DisplayedSource", 15 );
      if( node )
        m_currentlyShowingNuclide.deSerialize( node );
    }//if( showing_node )
    
    m_persisted.clear();
    persisted_node = base_node->first_node( "PersistedLines", 14 );
    if( persisted_node )
    {
      for( node = persisted_node->first_node( "DisplayedSource", 15 );
           node; node = node->next_sibling( "DisplayedSource", 15 ) )
      {
        ReferenceLineInfo nuc;
        nuc.deSerialize( node );
        m_persisted.push_back( nuc );
      }//for( loop over DisplayedSource nodes )
    }//if( persisted_node )
    
    node = base_node->first_node( "Shielding", 9 );
    if( node )
      m_shieldingSelect->deSerialize( node );
//    else
//      m_shieldingSelect-reset();
    
    auto currentColor = base_node->first_node( "CurrentLineColor", 16 );
    if( currentColor && currentColor->value_size() )
    {
      const string value( currentColor->value(), currentColor->value() + currentColor->value_size() );
      try{ m_colorSelect->setColor( WColor(value) ); }catch(...){}
      
      auto pickedAttrib = currentColor->first_attribute( "UserSelected", 12 );
      if( pickedAttrib && pickedAttrib->value_size() )
        m_userHasPickedColor = (pickedAttrib->value()[0] == '1');
    }else if( !m_currentlyShowingNuclide.lineColor.isDefault() )
    {
      m_colorSelect->setColor( m_currentlyShowingNuclide.lineColor );
    }
    
    //The quantities m_lineColors, m_peaksGetAssignedRefLineColor,
    //  m_previouslyPickedSourceColors, and m_specificSourcelineColors should
    //  get set by the color theme.
    
    refreshLinesDisplayedToGui( 100 );
//    updateDisplayChange();
    if( m_currentlyShowingNuclide.empty() && m_persisted.empty() )
      m_nuclidesCleared.emit();
    else
      m_displayingNuclide.emit();
  }catch( std::exception &e )
  {
    cerr << "ReferencePhotopeakDisplay::deSerialize() caught: " << e.what() << endl;
    stringstream msg;
    msg << "Error opening displayed photopeaks from database for display: " << e.what();
    passMessage( msg.str(), "", WarningWidget::WarningMsgHigh );
  }//try / catch
}//void deSerialize( std::string &xml_data  )


std::string ReferencePhotopeakDisplay::jsonReferenceLinesArray()
{
  string answer = "[";
  if( m_currentlyShowingNuclide.energies.size()
      && m_currentlyShowingNuclide.displayLines )
  {
    m_currentlyShowingNuclide.toJson(answer);
  }
  
  for( size_t i = 0; i < m_persisted.size(); ++i )
  {
    const ReferenceLineInfo &ref = m_persisted[i];
    if (m_currentlyShowingNuclide.energies.empty()
        || ref.parentLabel() != m_currentlyShowingNuclide.parentLabel())
    {
      answer += ",";
      ref.toJson(answer);
    }
  }
  answer += "]";
  return answer;
}//jsonReferenceLines()


std::map<std::string,std::string> ReferencePhotopeakDisplay::jsonReferenceLinesMap()
{
  std::map<std::string,std::string> answer;
  
  if( m_currentlyShowingNuclide.parentLabel() != "" )
    m_currentlyShowingNuclide.toJson( answer[m_currentlyShowingNuclide.parentLabel()] );
    
  for( size_t i = 0; i < m_persisted.size(); ++i )
  {
    const ReferenceLineInfo &ref = m_persisted[i];
    if( ref.parentLabel() != m_currentlyShowingNuclide.parentLabel() )
      ref.toJson( answer[ref.parentLabel()] );
  }
  
  return answer;
}//std::map<std::string,std::string> jsonReferenceLinesMap();


void ReferencePhotopeakDisplay::refreshLinesDisplayedToGui( int millisecdelay )
{
  //Note that the millisecond delay _may_ be vestigulal (untested), from when
  //  this->doJavaScript() was being called instead of wApp->doJavaScript().
  
#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )
  m_chart->clearAllReferncePhotoPeakLines();
  for( size_t i = 0; i < m_persisted.size(); ++i )
  {
    m_chart->setReferncePhotoPeakLines( m_persisted[i] );
    m_chart->persistCurrentReferncePhotoPeakLines();
  }
  m_chart->setReferncePhotoPeakLines( m_currentlyShowingNuclide );
#endif

#if( USE_SPECTRUM_CHART_D3 )
  return;
#else
  CanvasForDragging *canvas = m_chart->overlayCanvas();
  if( !canvas )
    return;
  
  const string canvasjs = "$('#c" + canvas->id() + "')";
  
  WStringStream js;
  
  //Arrrg, we need a setTimeout(...) here since when are calling this function
  //  from deSerialize(...), the call to clearAllLines() will wipe out the
  //  effect of this function {this is just a guess, may be something else}
  if( millisecdelay > 0 )
    js << "setTimeout(function(){";
  if( m_currentlyShowingNuclide.energies.size()
      && m_currentlyShowingNuclide.displayLines )
  {
    string json = canvasjs + ".data('gammalines',";
    m_currentlyShowingNuclide.toJson( json );
    json += ");";
    js << json;
    if( m_persistLines->isDisabled() )
       m_persistLines->enable();
    if( m_clearLines->isDisabled() )
      m_clearLines->enable();
  }else
  {
    m_persistLines->disable();
    m_clearLines->setDisabled( m_persisted.empty() );
  }//if( m_currentlyShowingNuclide.energies.size() )
  
  if( m_spectrumViewer && m_spectrumViewer->isPhone() )
    m_clearLines->setText( m_persisted.empty() ? "Remove" : "Remove All" );
  else
    m_clearLines->setText( m_persisted.empty() ? "Clear" : "Clear All" );
  
  for( size_t i = 0; i < m_persisted.size(); ++i )
  {
    if( !i )
      js << "var b=[];";
    
    string json;
    m_persisted[i].toJson( json );
    js << "b['" << m_persisted[i].parentLabel() << "']=" << json << ";";
  }//for( const ReferenceLineInfo *disp : m_persisted )

  if( m_persisted.size() )
    js << canvasjs << ".data('persistedLines',b);";
  else
    js << canvasjs << ".data('persistedLines',null);";
  
  js << "Wt.WT.DrawGammaLines('c" + canvas->id() + "',true);";
  
  if( millisecdelay > 0 )
    js << "}, " << millisecdelay << " );";
  
//  cerr << "Sending: " << js.str() << endl;
  
  wApp->doJavaScript( js.str() );
#endif  //#if( USE_SPECTRUM_CHART_D3 ) / else
}//void refreshLinesDisplayedToGui()


void ReferencePhotopeakDisplay::setIsotope( const SandiaDecay::Nuclide *nuc,
                                       double age )
{
  if( nuc )
  {
    m_nuclideEdit->setText( nuc->symbol );
    if( age >= 0.0 )
      m_ageEdit->setText( PhysicalUnits::printToBestTimeUnits(age) );
  }else
  {
    m_nuclideEdit->setText( "" );
  }//if( nuc ) / else
  
  handleIsotopeChange( true );
}//void setIsotope(...)


void ReferencePhotopeakDisplay::setNuclideAndAge( const string &name,
                                             const string &age )
{
  m_nuclideEdit->setText( name );
  m_ageEdit->setText( age );
  handleIsotopeChange( !age.empty() );
}//void setNuclideAndAge( const std::string &name, const std::string &age )


void ReferencePhotopeakDisplay::setShieldingMaterialAndThickness( const string &name,
                                                     const string &thickness )
{
  m_shieldingSelect->setMaterialNameAndThickness( name, thickness );
  updateDisplayChange();
}//void setShieldingMaterialAndThickness(...)


void ReferencePhotopeakDisplay::setElement( const SandiaDecay::Element *el )
{
  if( el )
  {
    m_nuclideEdit->setText( el->symbol );
  }else
  {
    m_nuclideEdit->setText( "" );
  }//if( nuc ) / else
  
  handleIsotopeChange( false );
}//void setIsotope(...)


void ReferencePhotopeakDisplay::setReaction( const ReactionGamma::Reaction *rctn )
{
  if( rctn )
  {
    m_nuclideEdit->setText( rctn->name() );
  }else
  {
    m_nuclideEdit->setText( "" );
  }//if( nuc ) / else
  
  handleIsotopeChange( false );
}//void setReaction(...)


void ReferencePhotopeakDisplay::fitPeaks()
{
  passMessage( "ReferencePhotopeakDisplay::fitPeaks() not implemented yet,"
               " but it will be really cool once it is", "", 2 );
  
  //-Get nuclide being displayed, including persisted, then do a whole lot of work...
  //-Get the existing peaks, and (temporarily fix their means).
  //-If the existing detector response (m_detectorDisplay->detector()) does not
  //   have resolution information, do a pre-fit of the highest amplitude
  //   candiates (above 100 keV, highest amp candidate in each 50 keV span
  //   maybye), as use this as a starting width response.  If this fales, than
  //   just use a generic NaI, HPGe, or LaBr resonse.
  //-Use this detector response to do the 'skyline' type filtering of potential
  //   peaks.  E.g. loop over each gamma line, and use the same guessing
  //   algorithm that asigning a nuclide line to a peak does, to see if that
  //   gamma line would be selected, if so, fit for that peak, if not, dont
  //   (btw, Ba133 window identifiaction method fails for NaI Ba133, so check
  //    that out and improve it)
  //-Now try to figure out which peaks should be fit in the same ROIs, and
  //   co-fit them, while indiviually fitting the others.  The existing peaks
  //   will have to be be tossed into the fit as well, but with fixed means.
  //   This will have to be re-evaluted after the initial fit.
  //-To actually assign the isotopes to the peaks, should consider the case
  //   of a single peak for multipe isotopes being shown: at first just identify
  //   candate nuclides, then disregard shielding and check the amplitude
  //   relative to the nearest candidate peak of each isotope (at an assumed
  //   age of course), and assign it to the peak that is closest.  Could also
  //   consider trying to make sure to assign at least some peaks to each
  //   nuclide.
  //-
  
  
//  m_spectrumViewer
  
}//void fitPeaks()

void ReferencePhotopeakDisplay::clearAllLines()
{
  //m_fitPeaks->disable();
  m_persisted.clear();
  m_currentlyShowingNuclide.reset();
  m_particleModel->clear();
  
  m_currentlyShowingNuclide.lineColor = m_lineColors[0];
  m_colorSelect->setColor( m_lineColors[0] );
  
  m_nuclideEdit->setText( "" );
  m_persistLines->disable();
  m_clearLines->disable();
  
  if( m_spectrumViewer && m_spectrumViewer->isPhone() )
    m_clearLines->setText( "Clear" );
  else
    m_clearLines->setText( "Remove" );

#if( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE )
  m_chart->clearAllReferncePhotoPeakLines();
#endif

#if( !USE_SPECTRUM_CHART_D3 )
  CanvasForDragging *canvas = m_chart->overlayCanvas();
  if( canvas )
  {
    const string canvasjs = "$('#c" + canvas->id() + "')";

    WStringStream js;
    js << canvasjs << ".data('gammalines',null);"
       << canvasjs << ".data('persistedLines',null);"
       << "Wt.WT.DrawGammaLines('c" << canvas->id() << "',true);";

    //We capp wApp->doJavaScript(...) incase this object will be deleted before
    //  the javaScript is sent to the client (in which case it wont be sent if
    //  calling WWidget::doJavaScript(...) )
    wApp->doJavaScript( js.str() );
  }//if( canvas )
#endif //USE_SPECTRUM_CHART_D3
  
  m_nuclidesCleared.emit();
}//void clearAllLines()




void ReferencePhotopeakDisplay::userColorSelectCallback( const WColor &color )
{
  if( color.isDefault() )
  {
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( BOOST_CURRENT_FUNCTION, "Was passed an invalid color." );
#endif
    m_userHasPickedColor = false;
    m_colorSelect->setColor( m_currentlyShowingNuclide.lineColor );
  }else
  {
    m_userHasPickedColor = true;
    m_currentlyShowingNuclide.lineColor = color;
    m_previouslyPickedSourceColors[m_currentlyShowingNuclide.labelTxt] = color;
    refreshLinesDisplayedToGui( 0 );
  }
}//void userColorSelectCallback( const std::string &color )



