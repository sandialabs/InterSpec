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
#include <string>
#include <vector>
#include <memory>
#include <iostream>

#include <Wt/WMenu>
#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WAnchor>
#include <Wt/WCheckBox>
#include <Wt/WGridLayout>
#include <Wt/WApplication>
#include <Wt/WItemDelegate>
#include <Wt/WStackedWidget>
#include <Wt/WContainerWidget>

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/EnergyCalTool.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/WarningWidget.h"
#include "SpecUtils/EnergyCalibration.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/IsotopeSelectionAids.h"

using namespace std;
using namespace Wt;

namespace
{
template<class T> struct index_compare_assend
{
  index_compare_assend(const T arr) : arr(arr) {} //pass the actual values you want sorted into here
  bool operator()(const size_t a, const size_t b) const
  {
    return arr[a] < arr[b];
  }
  const T arr;
};//struct index_compare
}

namespace EnergyCalImp
{
  class DeviationPairDisplay;
  class DevPair : public Wt::WContainerWidget
  {
  protected:
    DevPair( Wt::WContainerWidget *parent = 0 );
    void setDevPair( const std::pair<float,float> &d );
    std::pair<float,float> devPair() const;
    void visuallyIndicateChanged();
    //Wt::WDoubleSpinBox *m_energy, *m_offset;
    NativeFloatSpinBox *m_energy, *m_offset;
    Wt::WContainerWidget *m_delete;
    friend class DeviationPairDisplay;
  };//class DevPair
    
  class DeviationPairDisplay : public Wt::WContainerWidget
  {
  public:
    DeviationPairDisplay( Wt::WContainerWidget *parent = 0 );
    void setDeviationPairs( std::vector< std::pair<float,float> > d );
    std::vector< std::pair<float,float> > deviationPairs() const;
    void removeDevPair( DevPair *devpair );
    DevPair *newDevPair( const bool emitChangedNow );
    Wt::Signal<> &changed(); //emits when dev pair is added, deleted, or changed
    
    void setInvalidValues();
    void setValidValues();
    
  protected:
    void sortDisplayOrder( const bool indicateVisually );
    
    void emitChanged();
      
    Wt::Signal<> m_changed;
    Wt::WContainerWidget *m_pairs;
  };//class DeviationPairDisplay




DevPair::DevPair( Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    //m_energy( new WDoubleSpinBox() ),
    //m_offset( new WDoubleSpinBox() ),
    m_energy( new NativeFloatSpinBox() ),
    m_offset( new NativeFloatSpinBox() ),
    m_delete( new WContainerWidget() )
{
  WGridLayout* layout = new WGridLayout();
  layout->setContentsMargins(0, 0, 0, 0);
  layout->setVerticalSpacing( 0 );
  
  setLayout(layout);
  setStyleClass( "DevPair" );
  
  layout->addWidget(m_energy,0,0);
  layout->addWidget(m_offset,0,1);
  layout->addWidget(m_delete,0,2);
  layout->setColumnStretch(0, 1);
  layout->setColumnStretch(1, 1);
  layout->setColumnStretch(2, 0);
          
  m_energy->setStyleClass( "DevEnergy" );
  m_offset->setStyleClass( "DevOffset" );
  m_energy->setPlaceholderText( "Energy" );
  m_offset->setPlaceholderText( "Offset" );
  m_energy->setValueText( "" );
  m_offset->setValueText( "" );
          
  m_delete->addStyleClass( "Wt-icon DeleteDevPair" );
}//DevPair constructor


void DevPair::setDevPair( const std::pair<float,float> &d )
{
  m_energy->setValue( d.first );
  m_offset->setValue( d.second );
  
  auto printval = []( float val ) -> std::string {
    char buffer[64];
    const float fraction = val - std::floor(val);
    if( fraction == 0.0 )
      snprintf( buffer, sizeof(buffer), "%.0f", val );
    else if( fraction == 0.1 )
      snprintf( buffer, sizeof(buffer), "%.1f", val );
    else
      snprintf( buffer, sizeof(buffer), "%.2f", val );
    return buffer;
  };
 
  m_energy->setText( printval(d.first) );
  m_offset->setText( printval(d.second) );
}


std::pair<float,float> DevPair::devPair() const
{
  const float energy = m_energy->value();
  const float offset = m_offset->value();
  return std::make_pair( energy, offset );
}


void DevPair::visuallyIndicateChanged()
{
  doJavaScript( "$('#" + id() + "').fadeIn(100).fadeOut(100).fadeIn(100).fadeOut(100).fadeIn(100);" );
}//void visuallyIndicateChanged();


DeviationPairDisplay::DeviationPairDisplay( Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_pairs( NULL )
{
  addStyleClass( "DevPairDisplay" );
  WLabel *title = new WLabel( "Deviation Pairs", this );
  title->setStyleClass( "Wt-itemview Wt-header Wt-label DevPairTitle" );
  title->setInline( false );

  m_pairs = new WContainerWidget(this);
  m_pairs->setStyleClass( "DevPairsContainer" );
          
  auto footer = new WContainerWidget(this);
  footer->setStyleClass( "DevPairsFooter" );
  
  auto addBtn = new WContainerWidget( footer );
  addBtn->addStyleClass( "Wt-icon AddDevPair" );
  addBtn->clicked().connect( boost::bind(&DeviationPairDisplay::newDevPair, this, true) );
  addBtn->setToolTip( "Add another deviation pair" );
}//DeviationPairDisplay constructor


void DeviationPairDisplay::setDeviationPairs( vector< pair<float,float> > d )
{
  m_pairs->clear();
  std::sort( d.begin(), d.end() );
  
  for( size_t i = 0; i < d.size(); ++i )
  {
    DevPair *dev = newDevPair( false );
    dev->setDevPair( d[i] );
  }//for( size_t i = 0; i < d.size(); ++i )
  
  sortDisplayOrder(false);
}//setDeviationPairs(...)


void DeviationPairDisplay::sortDisplayOrder( const bool indicateVisually )
{
  vector<size_t> sort_indices;
  vector<DevPair *> displays;
  
  vector<float> offsets;
  const vector<WWidget *> childs = m_pairs->children();
  for( WWidget *t : childs )
  {
    DevPair *p = dynamic_cast<DevPair *>( t );
    if( p )
    {
      const size_t index = sort_indices.size();
      offsets.push_back( p->devPair().first );
      displays.push_back( p );
      sort_indices.push_back( index );
    }
  }//for( WWidget *t : childs )
  
  std::stable_sort( sort_indices.begin(), sort_indices.end(),
                    index_compare_assend<vector<float>&>(offsets) );

  bool order_changed = false;
  for( size_t i = 0; i < sort_indices.size(); ++i )
    order_changed |= (sort_indices[i] != i );
  
  if( !order_changed )
    return;
          
  for( size_t i = 0; i < displays.size(); ++i )
    m_pairs->removeWidget( displays[i] );
    
  for( size_t i = 0; i < displays.size(); ++i )
  {
    DevPair *p = displays[ sort_indices[i] ];
    m_pairs->addWidget( p );
    if( indicateVisually && (i != sort_indices[i]) )
      p->visuallyIndicateChanged();
  }
}//void sortDisplayOrder()


void DeviationPairDisplay::emitChanged()
{
  m_changed.emit();
}


vector< pair<float,float> > DeviationPairDisplay::deviationPairs() const
{
  vector< pair<float,float> > answer;
  const vector<WWidget *> childs = m_pairs->children();
  for( const WWidget *t : childs )
  {
    const DevPair *p = dynamic_cast<const DevPair *>( t );
    if( p )
      answer.push_back( p->devPair() );
  }//for( WWidget *t : childs )
  
  std::sort( answer.begin(), answer.end() );
  
  return answer;
}//deviationPairs()


void DeviationPairDisplay::removeDevPair( DevPair *devpair )
{
  if( !devpair )
    return;
  
  const vector<WWidget *> childs = m_pairs->children();
  for( WWidget *t : childs )
  {
    if( devpair == dynamic_cast<DevPair *>( t ) ) //dynamic_cast prob not necessary
    {
      delete devpair;
      sortDisplayOrder(false);
      emitChanged();
      return;
    }
  }//for( WWidget *t : childs )
}//removeDevPair(...)


void DeviationPairDisplay::setInvalidValues()
{
  addStyleClass( "InvalidDevPairs" );
}


void DeviationPairDisplay::setValidValues()
{
  removeStyleClass( "InvalidDevPairs" );
}


DevPair *DeviationPairDisplay::newDevPair( const bool emitChangedNow )
{
  DevPair *dev = new DevPair( m_pairs );
  dev->m_delete->clicked().connect( boost::bind( &DeviationPairDisplay::removeDevPair, this, dev ) );
  //dev->m_energy->valueChanged().connect( this, &DeviationPairDisplay::emitChanged );
  //dev->m_offset->valueChanged().connect( this, &DeviationPairDisplay::emitChanged );
  dev->m_energy->valueChanged().connect( this, &DeviationPairDisplay::emitChanged );
  dev->m_offset->valueChanged().connect( this, &DeviationPairDisplay::emitChanged );
  dev->m_energy->blurred().connect( boost::bind(&DeviationPairDisplay::sortDisplayOrder, this, true) );
  
  if( emitChangedNow )
    emitChanged();
  return dev;
}//newDevPair()

Wt::Signal<> &DeviationPairDisplay::changed()
{
  return m_changed;
}


class CoefDisplay : public WContainerWidget
{
public:
  WLabel *m_label;
  WCheckBox *m_fit;
  NativeFloatSpinBox *m_value;
  
  CoefDisplay( const size_t order, WContainerWidget *parent = nullptr )
  : WContainerWidget( parent ),
    m_label( nullptr ),
    m_fit( nullptr ),
    m_value( nullptr )
  {
    addStyleClass( "CoefDisplay" );
    
    switch( order )
    {
      case 0:  m_label = new WLabel( "Offset", this );    break;
      case 1:  m_label = new WLabel( "Linear", this );    break;
      case 2:  m_label = new WLabel( "Quadratic", this ); break;
      case 3:  m_label = new WLabel( "Quartic", this );   break;
      default: m_label = new WLabel( std::to_string(order) + "'th order", this ); break;
    }//switch( order )
    
    m_label->addStyleClass( "CoefLabel" );
    
    m_value = new NativeFloatSpinBox( this );
    m_value->addStyleClass( "CoefInput" );
    
    m_fit = new WCheckBox( "Fit", this );
    m_fit->addStyleClass( "CoefFit" );
  }//CoefDisplay
};//class CoefDisplay

class CalDisplay : public WContainerWidget
{
  static const size_t sm_min_coef_display ;
  
  EnergyCalTool *m_tool;
  WText *m_type;
  WContainerWidget *m_coefficients;
  DeviationPairDisplay *m_devPairs;
  shared_ptr<const SpecUtils::EnergyCalibration> m_cal;
  
public:
  CalDisplay( EnergyCalTool *tool, WContainerWidget *parent = nullptr )
  : WContainerWidget( parent ),
   m_tool( tool ),
   m_type( nullptr ),
   m_coefficients( nullptr ),
   m_devPairs( nullptr )
  {
    addStyleClass( "CalDisplay" );
    
    WGridLayout *layout = new WGridLayout( this );
    layout->setContentsMargins( 0, 0, 0, 0 );
    layout->setVerticalSpacing( 0 );
    layout->setHorizontalSpacing( 0 );
    
    WContainerWidget *coefDiv = new WContainerWidget();
    coefDiv->addStyleClass( "CoefCol" );
    layout->addWidget( coefDiv, 0, 0 );
    m_type = new WText();
    m_type->setInline( false );
    m_type->addStyleClass( "CalType" );
    m_coefficients = new WContainerWidget( coefDiv );
    m_coefficients->addStyleClass( "CoefContent" );
    
    m_devPairs = new DeviationPairDisplay();
    layout->addWidget( m_devPairs, 0, 1 );
    //m_devPairs->changed().connect(<#WObject *target#>, <#WObject::Method method#>)
  }//CalDisplay( constructor )
  
  void updateToGui( const shared_ptr<const SpecUtils::EnergyCalibration> &cal )
  {
    m_cal = cal;
    
    if( !m_cal )
    {
      m_type->setText( "No Calibration" );
      m_coefficients->clear();
      m_devPairs->setDeviationPairs( {} );
      return;
    }//if( !m_cal )
    
    const char *typetxt = "";
    switch( m_cal->type() )
    {
      case SpecUtils::EnergyCalType::LowerChannelEdge:    typetxt = "Lower Channel Energy"; break;
      case SpecUtils::EnergyCalType::InvalidEquationType: typetxt = "Not Defined";          break;
      case SpecUtils::EnergyCalType::Polynomial:          typetxt = "Polynomial";           break;
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
                                                          typetxt = "Default Polynomial";   break;
      case SpecUtils::EnergyCalType::FullRangeFraction:   typetxt = "Full Range Fraction";  break;
    }//switch( m_cal->type() )
    
    m_type->setText( typetxt );
    
    switch( m_cal->type() )
    {
      case SpecUtils::EnergyCalType::LowerChannelEdge:
      case SpecUtils::EnergyCalType::InvalidEquationType:
        m_coefficients->clear();
        m_devPairs->setDeviationPairs( {} );
        return;
              
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::FullRangeFraction:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        break;
    }//switch( m_cal->type() )
    
    const auto &devpairs = m_cal->deviation_pairs();
    const vector<float> &coeffs = m_cal->coefficients();
    
    m_devPairs->setDeviationPairs( devpairs );
    m_devPairs->changed().connect( std::bind( &EnergyCalTool::userChangedDeviationPair, m_tool, this) );
    
    
    const size_t num_coef_disp = std::max( coeffs.size(), sm_min_coef_display );
    
    size_t coefnum = 0;
    const auto existing = m_coefficients->children();
    for( size_t i = 0; i < existing.size(); ++i )
    {
      auto ww = dynamic_cast<CoefDisplay *>( existing[i] );
      assert( ww || !existing[i] );
      
      if( ww )
      {
        if( coefnum >= num_coef_disp )
        {
          m_coefficients->removeWidget( existing[i] );
          delete existing[i];
        }else
        {
          const float value = (coefnum < coeffs.size()) ? coeffs[coefnum] : 0.0f;
          ww->m_value->setValue( value );
        }
        
        ++coefnum;
      }//if( ww )
    }//for( size_t i = 0; i < existing.size(); ++i )
    
    for( ; coefnum < num_coef_disp; ++coefnum )
    {
      CoefDisplay *disp = new CoefDisplay( coefnum, m_coefficients );
      const float value = (coefnum < coeffs.size()) ? coeffs[coefnum] : 0.0f;
      disp->m_value->setValue( value );
      disp->m_fit->setChecked( (coefnum < 2) );
      
      disp->m_fit->checked().connect( m_tool, &EnergyCalTool::fitCoefficientCBChanged );
      disp->m_fit->unChecked().connect( m_tool, &EnergyCalTool::fitCoefficientCBChanged );
      
      disp->m_value->valueChanged().connect( boost::bind(&EnergyCalTool::userChangedCoefficient, m_tool,coefnum, this) );
    }
  }//updateToGui(...)
  
  
};//class CalDisplay

const size_t CalDisplay::sm_min_coef_display = 4;


}//namespace


EnergyCalTool::EnergyCalTool( InterSpec *viewer, PeakModel *peakModel, WContainerWidget *parent )
: WContainerWidget( parent ),
  m_interspec( viewer ),
  m_peakModel( peakModel ),
  m_peakTable( nullptr ),
  m_specTypeMenu( nullptr ),
  m_specTypeMenuStack( nullptr ),
  m_detectorMenu{ nullptr },
  m_calInfoDisplayStack( nullptr ),
  m_noCalTxt( nullptr ),
  m_moreActionsColumn( nullptr ),
  m_applyToColumn( nullptr ),
  m_detColumn( nullptr ),
  m_calColumn( nullptr ),
  m_peakTableColumn( nullptr ),
  m_layout( nullptr ),
  m_applyToCbs{ nullptr },
  m_moreActions{ nullptr }
{
  wApp->useStyleSheet( "InterSpec_resources/EnergyCalTool.css" );
  
  addStyleClass( "EnergyCalTool" );
  
  m_layout = new WGridLayout( this );
  m_layout->setContentsMargins( 0, 0, 0, 0 );
  m_layout->setVerticalSpacing( 0 );
  m_layout->setHorizontalSpacing( 0 );
  
  m_noCalTxt = new WText( "No spectrum loaded" );
  m_noCalTxt->addStyleClass( "NoCalContentTxt" );
  m_layout->addWidget( m_noCalTxt, 0, 0, AlignmentFlag::AlignCenter | AlignmentFlag::AlignMiddle );
  
  //Create the more actions column...
  m_moreActionsColumn = new WContainerWidget();
  m_moreActionsColumn->addStyleClass( "CalColumn MoreActionCol" );
  m_layout->addWidget( m_moreActionsColumn, 0, 1 );
  
  WGridLayout *collayout = new WGridLayout( m_moreActionsColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );
  
  WText *header = new WText( "More Actions" );
  header->addStyleClass( "ColHeader" );
  collayout->addWidget( header, 0, 0 );
  
  //We will put the apply-to list inside a div so we can style consistently with other rows
  // (a <ul> element doesnt accept same css as <div>, apparently).
  WContainerWidget *moreActionsDiv = new WContainerWidget();
  moreActionsDiv->addStyleClass( "CalColContent MoreActionsMenuContent" );
  collayout->addWidget( moreActionsDiv, 1, 0 );
  
  WContainerWidget *moreActionsList = new WContainerWidget( moreActionsDiv );
  moreActionsList->addStyleClass( "MoreActionsMenuList" );
  moreActionsList->setList( true );
  
  
  for( MoreActionsIndex index = static_cast<MoreActionsIndex>(0);
      index < MoreActionsIndex::NumMoreActionsIndex;
      index = MoreActionsIndex(index + 1) )
  {
    const char *label = "";
    switch( index )
    {
      case Linearize:       label = "Linearize..."; break;
      case Truncate:        label = "Truncate Energy..."; break;
      case CombineChannels: label = "Combine Channels..."; break;
      case ConvertToFrf:    label = "To FRF..."; break;
      case ConvertToPoly:   label = "To Polynomial..."; break;
      case NumMoreActionsIndex:
        assert(0);
        break;
    }//switch( index )
    
    WContainerWidget *holder = new WContainerWidget( moreActionsList );
    m_moreActions[index] = new WAnchor( WLink(), label, holder );
    m_moreActions[index]->clicked().connect( boost::bind(&EnergyCalTool::moreActionBtnClicked, this, index) );
  }//for( loop over more actions )
  
  
  
  // Create the "Apply To" column that determines what to apply changes to
  m_applyToColumn = new WContainerWidget();
  m_applyToColumn->addStyleClass( "CalColumn ApplyToCol" );
  m_layout->addWidget( m_applyToColumn, 0, 2 );
  
  collayout = new WGridLayout( m_applyToColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );
  
  header = new WText( "Apply Changes To" );
  header->addStyleClass( "ColHeader" );
  collayout->addWidget( header, 0, 0 );
  
  //We will put the apply-to list inside a div so we can style consistently with other rows
  // (a <ul> element doesnt accept same css as <div>, apparently).
  WContainerWidget *applyToDiv = new WContainerWidget();
  applyToDiv->addStyleClass( "CalColContent ApplyToMenuContent" );
  collayout->addWidget( applyToDiv, 1, 0 );
  
  WContainerWidget *applyToList = new WContainerWidget( applyToDiv );
  applyToList->addStyleClass( "ApplyToMenuList" );
  applyToList->setList( true );
  
  
  for( ApplyToCbIndex index = static_cast<ApplyToCbIndex>(0);
      index < ApplyToCbIndex::NumApplyToCbIndex;
      index = ApplyToCbIndex(index+1) )
  {
    const char *label = "";
    switch( index )
    {
      case ApplyToCbIndex::ApplyToForeground:         label = "Foreground";          break;
      case ApplyToCbIndex::ApplyToBackground:         label = "Background";          break;
      case ApplyToCbIndex::ApplyToSecondary:          label = "Secondary";           break;
      case ApplyToCbIndex::ApplyToDisplayedDetectors: label = "Displayed Detectors"; break;
      case ApplyToCbIndex::ApplyToAllDetectors:       label = "All Detectors";       break;
      case ApplyToCbIndex::ApplyToDisplayedSamples:   label = "Displayed Samples";   break;
      case ApplyToCbIndex::ApplyToAllSamples:         label = "All Samples";         break;
      case ApplyToCbIndex::NumApplyToCbIndex:
        assert( 0 );
        break;
    }//switch( index )
    
    WContainerWidget *item = new WContainerWidget( applyToList );
    item->addStyleClass( "ApplyToItem" );
    auto cb = new WCheckBox( label , item );
    cb->setWordWrap( false );
    cb->addStyleClass( "ApplyToItem" );
    cb->setInline( false );
    
    
    switch( index )
    {
      case ApplyToCbIndex::ApplyToForeground:
      case ApplyToCbIndex::ApplyToBackground:
      case ApplyToCbIndex::ApplyToSecondary:
      case ApplyToCbIndex::ApplyToAllDetectors:
      case ApplyToCbIndex::ApplyToAllSamples:
        cb->setChecked( true );
        break;
        
      case ApplyToCbIndex::ApplyToDisplayedDetectors:
      case ApplyToCbIndex::ApplyToDisplayedSamples:
      case ApplyToCbIndex::NumApplyToCbIndex:
        cb->setChecked( false );
        break;
    }//switch( index )
    
    cb->checked().connect( boost::bind( &EnergyCalTool::applyToCbChanged, this, index ) );
    cb->unChecked().connect( boost::bind( &EnergyCalTool::applyToCbChanged, this, index ) );
    
    m_applyToCbs[index] = cb;
  }//for( loop over ApplyToCbIndex )
  
  
  // Create the "Detector" column that determines which coefficients to show
  m_detColumn = new WContainerWidget();
  m_detColumn->addStyleClass( "CalColumn DetCol" );
  m_layout->addWidget( m_detColumn, 0, 3 );
  
  collayout = new WGridLayout( m_detColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 2, 1 );
  
  header = new WText( "Detector" );
  header->addStyleClass( "ColHeader" );
  collayout->addWidget( header, 0, 0 );
  
  WAnimation animation(Wt::WAnimation::Fade, Wt::WAnimation::Linear, 200);
  
  
  //m_specTypeMenuStack = new WStackedWidget();
  //m_specTypeMenuStack->addStyleClass( "CalColContent CalSpecStack" );
  //m_specTypeMenuStack->setTransitionAnimation( animation );
  
  //m_specTypeMenu = new WMenu( m_specTypeMenuStack );
  //collayout->addWidget( m_specTypeMenu, 1, 0 );
  //m_specTypeMenu->addStyleClass( "CalSpecMenu" );
  //m_specTypeMenu->itemSelected().connect( this, &EnergyCalTool::specTypeToDisplayForChanged );
  
  //collayout->addWidget( m_specTypeMenuStack, 2, 0 );
  
  //m_calInfoDisplayStack = new WStackedWidget();
  //m_calInfoDisplayStack->addStyleClass( "CalColContent CalStack" );
  //m_calInfoDisplayStack->setTransitionAnimation( animation );
  
  /*
  for( int i = 0; i < 3; ++i )
  {
    WContainerWidget *detMenuDiv = new WContainerWidget();
    detMenuDiv->addStyleClass( "CalColContent DetMenuDiv" );
    
    const char *label = "";
    switch( i )
    {
      case 0: label = "For."; break;
      case 1: label = "Back"; break;
      case 2: label = "Sec."; break;
    }
    WMenuItem *item = m_specTypeMenu->addItem( label, detMenuDiv, WMenuItem::LoadPolicy::PreLoading );
    //Fix issue, for Wt 3.3.4 at least, if user doesnt click exactly on the <a> element
    item->clicked().connect( boost::bind(&WMenuItem::select, item) );
    
    m_detectorMenu[i] = new WMenu( m_calInfoDisplayStack, detMenuDiv );
    m_detectorMenu[i]->addStyleClass( "VerticalMenu DetCalMenu" );
  }//for( int i = 0; i < 3; ++i )
   */
  
  
  // Create the "Coefficients" column that show the polynomial/FRF coefficents.
  m_calColumn = new WContainerWidget();
  m_calColumn->addStyleClass( "CalColumn CoefColumn" );
  m_layout->addWidget( m_calColumn, 0, 4 );
  
  
  collayout = new WGridLayout( m_calColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );
  
  header = new WText( "Calibration Coefficents" );
  header->addStyleClass( "ColHeader" );
  
  collayout->addWidget( header, 0, 0 );
  //collayout->addWidget( m_calInfoDisplayStack, 1, 0 );

  
  // Create the "Cal Peaks" table
  m_peakTableColumn = new WContainerWidget();
  m_peakTableColumn->addStyleClass( "CalColumn PeakTableCol" );
  m_layout->addWidget( m_peakTableColumn, 0, 5 );
  m_layout->setColumnStretch( 5, 1 );
  
  collayout = new WGridLayout( m_peakTableColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );
  
  header = new WText( "Calibration Peaks" );
  header->addStyleClass( "ColHeader" );
  collayout->addWidget( header, 0, 0 );
  
  WContainerWidget *peakTabelHolder = new WContainerWidget();
  peakTabelHolder->addStyleClass( "CalColContent" );
  collayout->addWidget( peakTabelHolder, 1, 0 );
  
  m_peakTable = new RowStretchTreeView();
  m_peakTable->addStyleClass( "PeakTable" );
  //collayout->addWidget( m_peakTable, 1, 0 );
  peakTabelHolder->addWidget( m_peakTable );
  
  
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
  
  m_peakTable->setColumnHidden( PeakModel::kUseForCalibration, false );
  m_peakTable->setColumnHidden( PeakModel::kMean, false );
  m_peakTable->setColumnHidden( PeakModel::kIsotope, false );
  m_peakTable->setColumnHidden( PeakModel::kPhotoPeakEnergy, false );
  m_peakTable->setColumnHidden( PeakModel::kDifference, false );
  
  
  m_peakTable->setColumnWidth( PeakModel::kUseForCalibration, WLength(3.7, WLength::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kMean, WLength(4.5, WLength::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kIsotope, WLength(4.5, WLength::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kPhotoPeakEnergy, WLength(7.25, WLength::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kDifference, WLength(6, WLength::FontEm) );
  
  
  
  WItemDelegate *dblDelagate = new WItemDelegate( m_peakTable );
  dblDelagate->setTextFormat( "%.2f" );
  m_peakTable->setItemDelegateForColumn( PeakModel::kMean, dblDelagate );
  
  PhotopeakDelegate *nuclideDelegate = new PhotopeakDelegate( PhotopeakDelegate::NuclideDelegate, true, m_peakTable );
  m_peakTable->setItemDelegateForColumn( PeakModel::kIsotope, nuclideDelegate );
  
  PhotopeakDelegate *photopeakDelegate = new PhotopeakDelegate( PhotopeakDelegate::GammaEnergyDelegate, true, m_peakTable );
  m_peakTable->setItemDelegateForColumn( PeakModel::kPhotoPeakEnergy, photopeakDelegate );
  
  
  m_interspec->displayedSpectrumChanged().connect( boost::bind( &EnergyCalTool::displayedSpectrumChanged, this, _1, _2, _3, _4 ) );
  refreshGuiFromFiles();
}//EnergyCalTool

EnergyCalTool::~EnergyCalTool()
{
}
  


void EnergyCalTool::fitCoefficientCBChanged()
{
  cerr << "fitCoefficientCBChanged()" << endl;
}//void fitCoefficientCBChanged()


void EnergyCalTool::userChangedCoefficient( const size_t coefnum, EnergyCalImp::CalDisplay *display )
{
  cerr << "userChangedCoefficient( " << coefnum << ", " << " )" << endl;
}//userChangedCoefficient(...)


void EnergyCalTool::userChangedDeviationPair( EnergyCalImp::CalDisplay *display )
{
  cerr << "userChangedDeviationPair( " << ", " << " )" << endl;
}//void userChangedDeviationPair( CalDisplay *display )


void EnergyCalTool::displayedSpectrumChanged( const SpecUtils::SpectrumType type,
                                              const std::shared_ptr<SpecMeas> &meas,
                                              const std::set<int> &samples,
                                              const std::vector<std::string> &detectors )
{
  static_assert( static_cast<int>(SpecUtils::SpectrumType::Foreground) == 0, "" );
  static_assert( static_cast<int>(SpecUtils::SpectrumType::SecondForeground) == 1, "" );
  static_assert( static_cast<int>(SpecUtils::SpectrumType::Background) == 2, "" );

  const int index = static_cast<int>( type );
  assert( index >= 0 && index < 3 );
  
  if( meas != m_currentSpecMeas[index] )
  {
    //whole new file
    cout << "EnergyCalTool::displayedSpectrumChanged: new file" << endl;
    
    //We want to cache original energy calibration, if we havent already
  }else if( samples != m_currentSampleNumbers[index] )
  {
    // Just changed what was displayed
    cout << "EnergyCalTool::displayedSpectrumChanged: changed sample numbers" << endl;
  }else
  {
    //no change...
    cout << "EnergyCalTool::displayedSpectrumChanged: same file and sample numbers" << endl;
  }

  m_currentSpecMeas[index] = meas;
  m_currentSampleNumbers[index] = samples;
  
  refreshGuiFromFiles();
}//void displayedSpectrumChanged(...)


void EnergyCalTool::setShowNoCalInfo( const bool nocal )
{
  m_noCalTxt->setHidden( !nocal );
  m_moreActionsColumn->setHidden( nocal );
  m_applyToColumn->setHidden( nocal );
  m_detColumn->setHidden( nocal );
  m_calColumn->setHidden( nocal );
  m_peakTableColumn->setHidden( nocal );
}//void setShowNoCalInfo( const bool nocal )


void EnergyCalTool::specTypeToDisplayForChanged()
{
  const int selectedType = m_specTypeMenu->currentIndex();
  if( selectedType < 0 || selectedType > 2 )
  {
    assert( 0 );  //shouldnt ever happen, right?
    return;
  }
  
  WMenu *detMenu = m_detectorMenu[selectedType];
  WMenuItem *detitem = detMenu->currentItem();
  if( !detitem && detMenu->count() )
    detMenu->select( 0 );
  if( detitem )
    detMenu->select( detitem );
}//void specTypeToDisplayForChanged();


void EnergyCalTool::displayedSpecChangedCallback( const SpecUtils::SpectrumType,
                                                  const std::shared_ptr<SpecMeas>,
                                                  const std::set<int>,
                                                  const std::vector<std::string> )
{
  /// \TODO: set the various m_applyToCbs if it is a new spectrum being shown.
  /// \TODO: if this is the first time seeing a SpecMeas, cache all of its energy calibration
  ///        information
  
  // \TODO: we could maybe save a little time by inspecting what was changed, but the added
  //        complexity probably isnt worth it, so we'll skip this.
  refreshGuiFromFiles();
}//void displayedSpecChangedCallback(...)


void EnergyCalTool::refreshGuiFromFiles()
{
  m_renderFlags |= EnergyCalToolRenderFlags::FullGuiUpdate;
  scheduleRender();
}//void refreshGuiFromFiles()


void EnergyCalTool::doRefreshFromFiles()
{
  string prevdet[3];
  int previousSpecInd = m_specTypeMenu ? m_specTypeMenu->currentIndex() : 0;
  
  const SpecUtils::SpectrumType spectypes[3] = {
    SpecUtils::SpectrumType::Foreground,
    SpecUtils::SpectrumType::Background,
    SpecUtils::SpectrumType::SecondForeground
  };
  
  shared_ptr<const SpecMeas> specfiles[3];
  for( int i = 0; i < 3; ++i )
    specfiles[i] = m_interspec->measurment( spectypes[i] );
  
  set<string> specdetnames[3]; //Just the names of gamma detectors with at least 4 channels
  
/*
  //Delete calibration contents and menu items - we will add in all the current ones below.
  //  - this appears to not work well - the stacks need replacing...
  for( int i = 0; i < 3; ++i)
  {
    Wt::WMenu *menu = m_detectorMenu[i];
    
    if( menu->currentItem() )
      prevdet[i] = menu->currentItem()->text().toUTF8();
    
    for( WMenuItem *item : menu->items() )
    {
      WWidget *content = item->contents();
      assert( content );
      menu->removeItem( item );
      
      delete item;
      delete content;
    }//
  }//for( Wt::WMenu *menu : m_detectorMenu )
*/

  
  //If we try to re-use the menu and stacks, for some reason the calibration coefficents wont show
  //  up if we alter which  background/secondary spectra are showing... not sure why, but for the
  //  moment we'll just re-create the menus and stacks... not great, but works, for the moment.
  {
    for( int i = 0; i < 3; ++i )
    {
      if( m_detectorMenu[i] )
      {
        if( m_detectorMenu[i]->currentItem() )
          prevdet[i] = m_detectorMenu[i]->currentItem()->text().toUTF8();
        delete m_detectorMenu[i];
      }
      m_detectorMenu[i] = nullptr;
    }
    delete m_specTypeMenuStack;
    delete m_specTypeMenu;
    m_specTypeMenuStack = nullptr;
    m_specTypeMenu = nullptr;
    
    WAnimation animation(Wt::WAnimation::Fade, Wt::WAnimation::Linear, 200);
    
    m_specTypeMenuStack = new WStackedWidget();
    m_specTypeMenuStack->addStyleClass( "CalColContent CalSpecStack" );
    m_specTypeMenuStack->setTransitionAnimation( animation );
    
    auto detlayout = dynamic_cast<WGridLayout *>( m_detColumn->layout() );
    auto callayout = dynamic_cast<WGridLayout *>( m_calColumn->layout() );
    assert( detlayout );
    assert( callayout );
    
    m_specTypeMenu = new WMenu( m_specTypeMenuStack );
    m_specTypeMenu->addStyleClass( "CalSpecMenu" );
    m_specTypeMenu->itemSelected().connect( this, &EnergyCalTool::specTypeToDisplayForChanged );
    detlayout->addWidget( m_specTypeMenu, 1, 0 );
    detlayout->addWidget( m_specTypeMenuStack, 2, 0 );
    
    if( m_calInfoDisplayStack )
      delete m_calInfoDisplayStack;
    m_calInfoDisplayStack = new WStackedWidget();
    m_calInfoDisplayStack->addStyleClass( "CalColContent CalStack" );
    m_calInfoDisplayStack->setTransitionAnimation( animation );
    callayout->addWidget( m_calInfoDisplayStack, 1, 0 );
    
    const char * const labels[3] = {"For.","Back","Sec."};
    
    /// \TODO: only create these menus when actually needed, so we wont need to
    for( int i = 0; i < 3; ++i )
    {
      if( !specfiles[i] )
        continue;
      
      WContainerWidget *detMenuDiv = new WContainerWidget();  //this holds the WMenu for this SpecFile
      detMenuDiv->addStyleClass( "CalColContent DetMenuDiv" );
      
      WMenuItem *item = m_specTypeMenu->addItem( labels[i], detMenuDiv, WMenuItem::LoadPolicy::PreLoading );
      //Fix issue, for Wt 3.3.4 at least, if user doesnt click exactly on the <a> element
      item->clicked().connect( boost::bind(&WMenuItem::select, item) );
      
      m_detectorMenu[i] = new WMenu( m_calInfoDisplayStack, detMenuDiv );
      m_detectorMenu[i]->addStyleClass( "VerticalMenu DetCalMenu" );
    }//for( int i = 0; i < 3; ++i )
  }
  
  
  if( !specfiles[0] )
  {
    setShowNoCalInfo( true );
    return;
  }
  
  if( previousSpecInd < 0 ||  previousSpecInd > 2
     || (previousSpecInd == 1 && !specfiles[1])
     || (previousSpecInd == 2 && !specfiles[2])  )
  {
    previousSpecInd = 0;
  }
  
  int nFilesWithCalInfo = 0;
  bool selectedDetToShowCalFor = false;
  bool hasFRFCal = false, hasPolyCal = false, hasLowerChanCal = false;
  
  
  for( int i = 0; i < 3; ++i )
  {
    Wt::WMenu *detMenu = m_detectorMenu[i];
    if( !detMenu )
      continue;
    
    WMenuItem *specItem = m_specTypeMenu->itemAt(i);
    assert( specItem );
    
    const SpecUtils::SpectrumType type = spectypes[i];
    shared_ptr<const SpecMeas> meas = specfiles[i];
    assert( meas );
    
    
    //We want the names of just the detectors that have gamma calibration information, of the
    //  currently displayed samples.
    //  We will assume the first Measurement for a given named detector will have gamma data if
    //  any of the Measurements from that detector will.
    //  \TODO: evaluate if that is true.
    set<string> detectors, nongammadets;
    const vector<string> &detnames = meas->gamma_detector_names();
    
    const set<int> &samples = m_interspec->displayedSamples(type);
    const vector<string> displayedDets = m_interspec->detectorsToDisplay(type);
    
    for( const int sample : samples )
    {
      for( const string &name : detnames )
      {
        if( detectors.count(name) || nongammadets.count(name) )
          continue;
        
        auto m = meas->measurement( sample, name );
        if( m && (m->num_gamma_channels() > 4) )
          detectors.insert( name );
        else if( m && !detectors.count(name) )
          nongammadets.insert(name);
      }//for( const string &name : detnames )
      
      if( (detectors.size() + nongammadets.size()) == detnames.size() )
        break;
    }//for( const int sample : samples )
    
    specdetnames[i] = detectors;
    
    if( detectors.empty() )
    {
      if( m_specTypeMenu->currentIndex() == i )
        m_specTypeMenu->select(0);
      specItem->setHidden( true );
      continue;
    }
    
    nFilesWithCalInfo += 1;
    

    assert( specItem );
    specItem->setHidden( false );
    
    for( const string &detname : detectors )
    {
      for( const int sample : samples )
      {
        auto m = meas->measurement( sample, detname );
        if( !m || (m->num_gamma_channels() <= 4) )
          continue;
        
        auto calcontent = new EnergyCalImp::CalDisplay( this );
        WMenuItem *item = detMenu->addItem( WString::fromUTF8(detname), calcontent, WMenuItem::LoadPolicy::PreLoading );
        
        //Fix issue, for Wt 3.3.4 at least, if user doesnt click exactly on the <a> element
        item->clicked().connect( boost::bind(&WMenuItem::select, item) );
        
        shared_ptr<const SpecUtils::EnergyCalibration> energycal = m->energy_calibration();
        calcontent->updateToGui( energycal );
        
        if( energycal )
        {
          switch( energycal->type() )
          {
            case SpecUtils::EnergyCalType::Polynomial:
            case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
              hasPolyCal = true;
              break;
              
            case SpecUtils::EnergyCalType::FullRangeFraction:
              hasFRFCal = true;
              break;
              
            case SpecUtils::EnergyCalType::LowerChannelEdge:
              hasLowerChanCal = true;
              break;
            
            case SpecUtils::EnergyCalType::InvalidEquationType:
              break;
          }//switch( energycal->type() )
        }//if( energycal )
        
        
        if( (detname == prevdet[i]) && (i == previousSpecInd) )
        {
          m_specTypeMenu->select( previousSpecInd );
          detMenu->select( item );
          selectedDetToShowCalFor = true;
        }
        
        break;
      }
    }//for( const string &detname : detectors )
  }//for( int i = 0; i < 3; ++i )
  
  if( !selectedDetToShowCalFor && m_detectorMenu[0] && m_detectorMenu[0]->count() )
  {
    m_specTypeMenu->select( 0 );
    m_detectorMenu[0]->select( 0 );
  }
  
  setShowNoCalInfo( !nFilesWithCalInfo );
  
  
  //DOnt show spectype menu if we dont need to
  const bool showSpecType = ( (nFilesWithCalInfo < 2)
                              || ( (!specfiles[1] || (specfiles[0]==specfiles[1]))
                                    && (!specfiles[2] || (specfiles[0]==specfiles[2]))) );
  m_specTypeMenu->setHidden( showSpecType );
  
  bool anyApplyToCbShown = false;
  for( ApplyToCbIndex index = static_cast<ApplyToCbIndex>(0);
      index < ApplyToCbIndex::NumApplyToCbIndex;
      index = ApplyToCbIndex(index+1) )
  {
    Wt::WCheckBox *cb = m_applyToCbs[index];
    const auto cbparent = cb->parent();
    assert( cbparent );
    assert( dynamic_cast<WContainerWidget *>(cbparent) );
    
    bool hideRow = false;
    switch( index )
    {
      case ApplyToCbIndex::ApplyToForeground:
        hideRow = (!specfiles[0] || (!specfiles[1] && !specfiles[2]));
        break;
        
      case ApplyToCbIndex::ApplyToBackground:
        // We'll ignore case where foreground and background is same {SpecFile, Samples, Detectors}
        hideRow = !specfiles[1];
        break;
        
      case ApplyToCbIndex::ApplyToSecondary:
        // We'll ignore case where secondary and background is same {SpecFile, Samples, Detectors}
        hideRow = !specfiles[2];
        break;
        
      case ApplyToCbIndex::ApplyToDisplayedDetectors:
      {
        bool displayingAll = true;
        for( int i = 0; displayingAll && (i < 3); ++i )
        {
          const auto meas = specfiles[i];
          if( !meas )
            continue;
          
          const auto type = spectypes[i];
          const vector<string> displayed = m_interspec->detectorsToDisplay( type );
          for( const auto &name : meas->gamma_detector_names() )
          {
            if( std::find( begin(displayed), end(displayed), name ) == end(displayed) )
              displayingAll = false;
          }//for( const auto &name : meas->gamma_detector_names() )
        }//for( loop over the types of spectrum files )
        
        hideRow = displayingAll;
        if( displayingAll )
          cb->setChecked( false );
        
        break;
      }//case ApplyToCbIndex::ApplyToDisplayedDetectors:
        
      case ApplyToCbIndex::ApplyToAllDetectors:
      {
        static_assert( ApplyToCbIndex::ApplyToDisplayedDetectors
                       < ApplyToCbIndex::ApplyToAllDetectors, "" );
        
        const WCheckBox *dispDetsCb = m_applyToCbs[ApplyToCbIndex::ApplyToDisplayedDetectors];
        assert( dispDetsCb->parent() );
        
        const bool dispDetsHid = dispDetsCb->parent()->isHidden();
        
        hideRow = dispDetsHid;
        if( dispDetsHid )
          cb->setChecked( true );
        
        if( dispDetsCb->isChecked() )
          cb->setChecked( false );
        
        break;
      }//case ApplyToAllDetectors:
        
      case ApplyToCbIndex::ApplyToDisplayedSamples:
      {
        // We want to check that for each displayed unique SpecMeas we are displaying all samples.
        /// \TODO: avoid allocating all the set<int>'s below, and then the looping threw every value
        ///        Can probably easily use a hyristic to avoid most of the time, or change how
        ///        things are tracked to avoid probably all the time (at the cost of adding
        ///        complexity to the code)
        map<shared_ptr<const SpecMeas>,set<int>> undisplayed;
        for( int i = 0; i < 3; ++i )
        {
          const auto meas = specfiles[i];
          if( !meas )
            continue;
          
          if( !undisplayed.count(meas) )
            undisplayed[meas] = meas->sample_numbers();
          
          set<int> &undispsamples = undisplayed[meas];
          for( const int sample : m_interspec->displayedSamples(spectypes[i]) )
            undispsamples.erase( sample );
        }//for( loop over the types of spectrum files )
        
        bool displayingAll = true;
        for( const auto &p : undisplayed )
          displayingAll = (displayingAll && p.second.empty());
        
        hideRow = displayingAll;
        if( displayingAll )
          cb->setChecked( false );
        
        break;
      }//case ApplyToDisplayedSamples:
        
      case ApplyToCbIndex::ApplyToAllSamples:
      {
        static_assert( ApplyToCbIndex::ApplyToDisplayedSamples
                         < ApplyToCbIndex::ApplyToAllSamples, "" );
        
        const WCheckBox *dispSamplesCb = m_applyToCbs[ApplyToCbIndex::ApplyToDisplayedSamples];
        assert( dispSamplesCb->parent() );
        const bool dispSamplesHid = dispSamplesCb->parent()->isHidden();
        
        hideRow = dispSamplesHid;
        if( dispSamplesHid )
          cb->setChecked( true );
        
        if( dispSamplesCb->isChecked() )
          cb->setChecked( false );
        
        break;
      }//case ApplyToCbIndex::ApplyToAllSamples:
        
        
      case ApplyToCbIndex::NumApplyToCbIndex:
        assert( 0 );
        break;
    }//switch( index )
    
    cbparent->setHidden( hideRow );
    
    if( !hideRow )
      anyApplyToCbShown = true;
  }//for( loop over ApplyToCbIndex )
  
  m_applyToColumn->setHidden( !anyApplyToCbShown );
  
  bool hideDetCol = true;
  if( specfiles[0] && specdetnames[0].size() > 1 )
    hideDetCol = false;
  if( specfiles[1] && (specfiles[0] != specfiles[1]) )
    hideDetCol = false;
  if( specfiles[2] && (specfiles[0] != specfiles[2]) )
    hideDetCol = false;
  
  m_detColumn->setHidden( hideDetCol );
  
  for( MoreActionsIndex index = MoreActionsIndex(0);
      index < MoreActionsIndex::NumMoreActionsIndex;
      index = MoreActionsIndex(index + 1) )
  {
    Wt::WAnchor *anchor = m_moreActions[index];
    assert( anchor );
    auto aparent = anchor->parent();
    assert( dynamic_cast<WContainerWidget *>(aparent) );
    
    switch( index )
    {
      case MoreActionsIndex::Linearize:
      case MoreActionsIndex::Truncate:
      case MoreActionsIndex::CombineChannels:
        aparent->setHidden( !specfiles[0] );
        break;
        
      case MoreActionsIndex::ConvertToFrf:
        aparent->setHidden( hasPolyCal );
        break;
        
      case MoreActionsIndex::ConvertToPoly:
        aparent->setHidden( hasFRFCal || hasLowerChanCal );
        break;
        
      case MoreActionsIndex::NumMoreActionsIndex:
        break;
    }//switch( index )
  }//for( loop over
  
  
  //const int currentwidget = m_detectorMenu[0]->contentsStack()->currentIndex();
  //cout << "currentwidget=" << currentwidget << endl;
}//void doRefreshFromFiles()



void EnergyCalTool::moreActionBtnClicked( const EnergyCalTool::MoreActionsIndex index )
{
  cerr << "moreActionBtnClicked: " << index << endl;
}//void moreActionBtnClicked( const MoreActionsIndex index )


void EnergyCalTool::render( Wt::WFlags<Wt::RenderFlag> flags)
{
  //flags.testFlag(RenderFlag::RenderFull) will only be true on initial rending of widget, and
  //  after that only the RenderFlag::RenderUpdate flag will be set
  
  if( flags.testFlag(Wt::RenderFlag::RenderFull)
      || m_renderFlags.testFlag(EnergyCalToolRenderFlags::FullGuiUpdate) )
  {
    doRefreshFromFiles();
    m_renderFlags.clear( EnergyCalToolRenderFlags::FullGuiUpdate );
  }
  
  WContainerWidget::render(flags);
}//void render( Wt::WFlags<Wt::RenderFlag> flags)


void EnergyCalTool::applyToCbChanged( const EnergyCalTool::ApplyToCbIndex index )
{
  assert( index >= 0 && index <= EnergyCalTool::NumApplyToCbIndex );
  WCheckBox *cb = m_applyToCbs[index];
  const bool isChecked = cb->isChecked();
  
  switch( index )
  {
    case EnergyCalTool::ApplyToForeground:
    case EnergyCalTool::ApplyToBackground:
    case EnergyCalTool::ApplyToSecondary:
      break;
    
    case EnergyCalTool::ApplyToDisplayedDetectors:
      m_applyToCbs[EnergyCalTool::ApplyToAllDetectors]->setChecked( !isChecked );
      break;
      
    case EnergyCalTool::ApplyToAllDetectors:
      m_applyToCbs[EnergyCalTool::ApplyToDisplayedDetectors]->setChecked( !isChecked );
      break;
      
    case EnergyCalTool::ApplyToDisplayedSamples:
      m_applyToCbs[EnergyCalTool::ApplyToAllSamples]->setChecked( !isChecked );
      break;
      
    case EnergyCalTool::ApplyToAllSamples:
      m_applyToCbs[EnergyCalTool::ApplyToDisplayedSamples]->setChecked( !isChecked );
      break;
      
    case EnergyCalTool::NumApplyToCbIndex:
      break;
  }//switch( index )
}//void applyToCbChanged( const ApplyToCbIndex index )
