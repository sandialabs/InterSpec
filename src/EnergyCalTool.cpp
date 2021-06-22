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
#include <ctime>
#include <memory>
#include <iostream>

#include <Wt/WMenu>
#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WAnchor>
#include <Wt/WCheckBox>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WItemDelegate>
#include <Wt/WStackedWidget>
#include <Wt/WContainerWidget>

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/EnergyCalTool.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/SpectraFileModel.h"
#include "SpecUtils/EnergyCalibration.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/EnergyCalGraphical.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/EnergyCalAddActions.h"
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
    void setInvalidValues();
    void setValidValues();
    
    //void setMsg( const string &msg );
    
    enum class UserFieldChanged : int
    {
      AddedDeviationPair,
      RemovedDeviationPair,
      EnergyChanged,
      OffsetChanged
    };//enum UserFieldChanged
    
    /** Signal that gets emited when dev pair is added, deleted, or changed.
     Int argument is of type UserFieldChanged endum
     */
    Wt::Signal<int> &changed();
    
  protected:
    void sortDisplayOrder( const bool indicateVisually );
    
    void emitChanged( const UserFieldChanged whatChanged );
      
    Wt::Signal<int> m_changed;
    Wt::WContainerWidget *m_pairs;
    //Wt::WText *m_msg;
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
    //, m_msg( nullptr )
{
  addStyleClass( "DevPairDisplay" );
  WLabel *title = new WLabel( "Deviation Pairs", this );
  title->setStyleClass( "Wt-itemview Wt-header Wt-label DevPairTitle" );
  title->setInline( false );

  m_pairs = new WContainerWidget(this);
  m_pairs->setStyleClass( "DevPairsContainer" );
          
  //m_msg = new WText( "&nbsp;", this );
  //m_msg->addStyleClass( "DevPairMsg" );
  //m_msg->setHidden( true );
  
  auto footer = new WContainerWidget(this);
  footer->setStyleClass( "DevPairsFooter" );
  
  auto addBtn = new WContainerWidget( footer );
  addBtn->addStyleClass( "Wt-icon AddDevPair" );
  addBtn->clicked().connect( boost::bind(&DeviationPairDisplay::newDevPair, this, true) );
  addBtn->setToolTip( "Add another deviation pair" );
}//DeviationPairDisplay constructor

/*
void DeviationPairDisplay::setMsg( const string &msg )
{
  if( msg == m_msg->text().toUTF8() )
    return;
  
  const bool wasHidden = m_msg->isHidden();
  const bool shouldBeHidden = msg.empty();
  if( wasHidden != shouldBeHidden )
  {
    //WAnimation anim (WAnimation::AnimationEffect::Fade | WAnimation::AnimationEffect::SlideInFromBottom, WAnimation::Linear, 200);
    m_msg->setHidden( shouldBeHidden, anim );
  }
  m_msg->setText( WString::fromUTF8(msg) );
}
 */

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


void DeviationPairDisplay::emitChanged( const UserFieldChanged whatChanged )
{
  m_changed.emit( static_cast<int>(whatChanged) );
}


vector< pair<float,float> > DeviationPairDisplay::deviationPairs() const
{
  vector< pair<float,float> > answer;
  const vector<WWidget *> childs = m_pairs->children();
  for( const WWidget *t : childs )
  {
    const DevPair *p = dynamic_cast<const DevPair *>( t );
    
    if( !p )
      continue;
    
    // We will only insert a deviation pair if both fields have text entered.
    if( p->m_energy->text().empty() || p->m_offset->text().empty() )
      continue;
    
    answer.push_back( p->devPair() );
  }//for( WWidget *t : childs )
  
  std::sort( answer.begin(), answer.end() );
  
  // The EnergyCalibration code implicitly will add in a {0,0} deviation pair if its needed, so
  //  might as well just add it in now so it isnt implict to the user
  //if( (answer.size() == 1) && (answer[0].first > 0.01f) )
  //  answer.insert( begin(answer), {0.0f,0.0f} );
  
  // Remove duplicates, usually {0,0}
  //const double epsilon = 1.0;
  //for( size_t i = 1; i < answer.size(); ++i )
  //{
  //  const auto &prev = answer[i-1];
  //  if( (fabs(answer[i].first - prev.first) < epsilon)
  //      && (fabs(answer[i].second - prev.second) < epsilon) )
  //    answer.erase( begin(answer) + i );
  //}
  
  return answer;
}//deviationPairs()


void DeviationPairDisplay::removeDevPair( DevPair *devpair )
{
  if( !devpair )
    return;
  
  const vector<WWidget *> childs = m_pairs->children();
  for( WWidget *t : childs )
  {
    auto tt = dynamic_cast<DevPair *>( t ); //dynamic_cast prob not necessary
    if( devpair == tt )
    {
      delete devpair;
      sortDisplayOrder(false);
      emitChanged( UserFieldChanged::RemovedDeviationPair );
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
  dev->m_energy->valueChanged().connect( boost::bind( &DeviationPairDisplay::emitChanged, this, UserFieldChanged::EnergyChanged ) );
  dev->m_offset->valueChanged().connect( boost::bind( &DeviationPairDisplay::emitChanged, this, UserFieldChanged::OffsetChanged ) );
  dev->m_energy->blurred().connect( boost::bind(&DeviationPairDisplay::sortDisplayOrder, this, true) );
  
  if( emitChangedNow )
    emitChanged( UserFieldChanged::AddedDeviationPair );
  
  return dev;
}//newDevPair()

Wt::Signal<int> &DeviationPairDisplay::changed()
{
  return m_changed;
}


class CoefDisplay : public WContainerWidget
{
public:
  const size_t m_order;
  WLabel *m_label;
  WCheckBox *m_fit;
  NativeFloatSpinBox *m_value;
  
  CoefDisplay( const size_t order, WContainerWidget *parent = nullptr )
  : WContainerWidget( parent ),
    m_order( order ),
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
    m_value->setSpinnerHidden( true );
    
    m_fit = new WCheckBox( "Fit", this );
    m_fit->addStyleClass( "CoefFit" );
  }//CoefDisplay
};//class CoefDisplay

class CalDisplay : public WContainerWidget
{
  static const size_t sm_min_coef_display ;
  
  EnergyCalTool *m_tool;
  const SpecUtils::SpectrumType m_cal_type;
  const std::string m_det_name;
  
  WText *m_type;
  WText *m_convertMsg;
  WContainerWidget *m_coefficients;
  DeviationPairDisplay *m_devPairs;
  
  //I cant decide if I like hiding deviation pair widget when there is none, or not
#define HIDE_EMPTY_DEV_PAIRS 0
#if( HIDE_EMPTY_DEV_PAIRS )
  WPushButton *m_addPairs;
#endif
  
  shared_ptr<const SpecUtils::EnergyCalibration> m_cal;
  
public:
  CalDisplay( EnergyCalTool *tool,
             const SpecUtils::SpectrumType type,
             const std::string &detname,
             const bool isWideLayout,
             WContainerWidget *parent = nullptr )
  : WContainerWidget( parent ),
   m_tool( tool ),
   m_cal_type( type ),
   m_det_name( detname ),
   m_type( nullptr ),
   m_convertMsg( nullptr ),
   m_coefficients( nullptr ),
   m_devPairs( nullptr )
#if( HIDE_EMPTY_DEV_PAIRS )
   , m_addPairs( nullptr )
#endif
  {
    addStyleClass( "CalDisplay" );
    
    WGridLayout *layout = new WGridLayout( this );
    layout->setContentsMargins( 0, 0, 0, 0 );
    layout->setVerticalSpacing( 0 );
    layout->setHorizontalSpacing( 0 );
    
    WContainerWidget *coefDiv = new WContainerWidget();
    coefDiv->addStyleClass( "CoefCol" );
    layout->addWidget( coefDiv, 0, 0 );
    
    m_type = new WText( "&nbsp;", coefDiv );
    m_type->setInline( false );
    m_type->addStyleClass( "Wt-itemview Wt-header Wt-label CalType" );
    
    m_coefficients = new WContainerWidget( coefDiv );
    m_coefficients->addStyleClass( "CoefContent" );
    
    m_devPairs = new DeviationPairDisplay();
    
#if( HIDE_EMPTY_DEV_PAIRS )
    //For files with multiple detectors, the "Add dev. pairs" buttons doesnt show up right
    // for the detectors not currently showing - I guess should toggle dev pairs for all detectors.
    m_devPairs->setHidden( true );
    m_addPairs = new WPushButton( "Add dev. pairs" );
    m_addPairs->addStyleClass( "LinkBtn" );
    //m_addPairs->setIcon( "InterSpec_resources/images/plus_min_white.svg" );
    m_addPairs->setHidden( true );
    m_addPairs->clicked().connect( this, &CalDisplay::showDevPairs );
#endif
    
    if( isWideLayout )
    {
      layout->addWidget( m_devPairs, 0, 1 );
#if( HIDE_EMPTY_DEV_PAIRS )
      layout->addWidget( m_addPairs, 1, 0, AlignmentFlag::AlignRight );
      layout->setRowStretch( 0, 1 );
#endif
    }else
    {
#if( HIDE_EMPTY_DEV_PAIRS )
      layout->addWidget( m_addPairs, 1, 0, AlignmentFlag::AlignCenter );
      layout->addWidget( m_devPairs, 2, 0 );
#else
      layout->addWidget( m_devPairs, 1, 0 );
#endif
      m_devPairs->setHeight( 100 );
    }
    
    m_devPairs->changed().connect( boost::bind( &EnergyCalTool::userChangedDeviationPair, m_tool, this, _1 ) );
  }//CalDisplay( constructor )
  
  SpecUtils::SpectrumType spectrumType() const { return m_cal_type; }
  const std::string &detectorName() const { return m_det_name; }
  
    
  //void setDeviationPairMsg( const std::string &msg )
  //{
  //  m_devPairs->setMsg( msg );
  //}//void setDeviationPairMsg( const std::string &msg )
  
  
  shared_ptr<const SpecUtils::EnergyCalibration> lastSetCalibration()
  {
    return m_cal;
  }
  
  /// @param fitfor The order coefficients that should be set checked to fit for
  void setFitFor( const set<size_t> &fitfor )
  {
    for( auto w : m_coefficients->children() )
    {
      auto ww = dynamic_cast<const CoefDisplay *>( w );
      assert( ww );
      if( ww )
        ww->m_fit->setChecked( fitfor.count(ww->m_order) );
    }//for( auto w : m_coefficients->children() )
  }//void setFitFor( const set<size_t> &fitfor )
  
#if( HIDE_EMPTY_DEV_PAIRS )
  void showDevPairs()
  {
    m_addPairs->setHidden( true );
    m_devPairs->setHidden( false, WAnimation(WAnimation::Fade, WAnimation::Linear, 200) );
  }
#endif
  
  
  /// @returns The order coefficents that are checked to be fit for
  set<size_t> fitForCoefficents() const
  {
    set<size_t> coeffs;
    
    for( auto w : m_coefficients->children() )
    {
      auto ww = dynamic_cast<const CoefDisplay *>( w );
      assert( ww );
      if( !ww )
        continue;
      
      if( ww->m_fit->isChecked() )
        coeffs.insert( ww->m_order );
    }//for( auto w : m_coefficients->children() )
    
    return coeffs;
  }//vector<bool> fitForCoefficents() const
  
  
  vector<float> displayedCoefficents()
  {
    vector<float> coeffs;
  
    const auto existing = m_coefficients->children();
    for( size_t i = 0; i < existing.size(); ++i )
    {
      if( !existing[i] )
        continue;  //shouldnt happen, but JIC
      
      auto ww = dynamic_cast<const CoefDisplay *>( existing[i] );
      assert( ww );
      if( !ww )
        continue;
    
      //NativeFloatSpinBox avoids round-off errors, and we will rely on this here.
      const float dispvalue = ww->m_value->value();
      coeffs.push_back( dispvalue );
    }//for( size_t i = 0; i < existing.size(); ++i )
    
    //remove trailing zeros
    while( !coeffs.empty() && coeffs.back()==0.0f )
      coeffs.resize( coeffs.size() - 1 );
    
    return coeffs;
  }//vector<float> displayedCoefficents() const
  
  
  std::vector< std::pair<float,float> > displayedDeviationPairs() const
  {
    return m_devPairs->deviationPairs();
  }
  
  void setDeviationPairsInvalid()
  {
    m_devPairs->setInvalidValues();
  }
  
  void setDeviationPairsValid()
  {
    m_devPairs->setValidValues();
  }
  
  void updateToGui( const shared_ptr<const SpecUtils::EnergyCalibration> &cal )
  {
#if( HIDE_EMPTY_DEV_PAIRS )
    const bool hadCal = !!m_cal;
#endif
    
    m_cal = cal;
    
    if( !m_cal )
    {
      m_type->setText( "No Calibration" );
      m_coefficients->clear();
      m_devPairs->setDeviationPairs( {} );
#if( HIDE_EMPTY_DEV_PAIRS )
      m_devPairs->setHidden( true );
      m_addPairs->setHidden( true );
#endif
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
        if( !m_devPairs->isHidden() )
          m_devPairs->hide();

        if( !m_convertMsg )
        {
          if( auto p = dynamic_cast<WContainerWidget *>( m_type->parent() ) )
          {
            m_convertMsg = new WText( "Please convert to Polynomial calibration to edit", p );
            m_convertMsg->addStyleClass( "ConvertToPolyMsg" );
            m_convertMsg->setInline( false );
          }
        }//if( !m_convertMsg )
        
#if( HIDE_EMPTY_DEV_PAIRS )
        m_addPairs->setHidden( true );
        m_devPairs->setHidden( true );
#endif
        return;
              
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::FullRangeFraction:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
#if( !HIDE_EMPTY_DEV_PAIRS )
        if( m_devPairs->isHidden() )
          m_devPairs->show();
#endif
        if( m_convertMsg )
        {
          delete m_convertMsg;
          m_convertMsg = nullptr;
        }
        break;
    }//switch( m_cal->type() )
    
    const auto &devpairs = m_cal->deviation_pairs();
    const vector<float> &coeffs = m_cal->coefficients();
    
#if( HIDE_EMPTY_DEV_PAIRS )
    // Once deviation pairs are showing, we will leave them showing, even if there are no deviation
    //  pairs; this is because if you click to add deviation pairs, dont add any, then adjust the
    //  the gain, the deviation pair display getting hidden, causing a layout update, is really
    //  jarring.
    const auto anim = hadCal ? WAnimation(WAnimation::Fade, WAnimation::Linear, 200) : WAnimation{};
    if( m_devPairs->isHidden() && !devpairs.empty() )
      m_devPairs->setHidden( false, anim );
    m_addPairs->setHidden( !m_devPairs->isHidden(), anim );
#endif
    
    m_devPairs->setDeviationPairs( devpairs );
    //m_devPairs->changed().connect( std::bind( &EnergyCalTool::userChangedDeviationPair, m_tool, this) );
    
    const size_t num_coef_disp = std::max( coeffs.size(), sm_min_coef_display );
    vector<CoefDisplay *> coef_disps( num_coef_disp, nullptr );
    
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
          assert( coefnum < coef_disps.size() );
          coef_disps[coefnum] = ww;
          const float value = (coefnum < coeffs.size()) ? coeffs[coefnum] : 0.0f;
          ww->m_value->setValue( value );
        }
        
        ++coefnum;
      }//if( ww )
    }//for( size_t i = 0; i < existing.size(); ++i )
    
    for( ; coefnum < num_coef_disp; ++coefnum )
    {
      CoefDisplay *disp = new CoefDisplay( coefnum, m_coefficients );
      coef_disps[coefnum] = disp;
      const float value = (coefnum < coeffs.size()) ? coeffs[coefnum] : 0.0f;
      disp->m_value->setValue( value );
      disp->m_fit->setChecked( (coefnum < 2) );
      
      disp->m_fit->checked().connect( m_tool, &EnergyCalTool::updateFitButtonStatus );
      disp->m_fit->unChecked().connect( m_tool, &EnergyCalTool::updateFitButtonStatus );
      
      /* Note: if the user uses the up.down arrows in a NativeFloatSpinBox to change values, things
               get all messed up (new values get set via c++ messing  up current values, or the
               valueChanged() callback gets called like 10 times per second, causing changes faster
               than everything can keep up, and just generally poor working), so for now I have
               disabled these spinners via #NativeFloatSpinBox::setSpinnerHidden()
       */
      disp->m_value->valueChanged().connect( boost::bind(&EnergyCalTool::userChangedCoefficient, m_tool, coefnum, this) );
    }
    
    
    //Set the step size to move the upper range of energy by about 1 keV per step
    // Set up the little tick/spin/whatever boxes
    /*
     //The NativeFloatSpinBox::setSpinnerHidden() call is currently removing the spin-box up/down
     //  arrow, so we wont set the step-size, as on Firefox if we fit for a value, then it will turn
     //  red if the new value doesnt hit on the step size
    for( size_t i = 0; i < coef_disps.size(); ++i )
    {
      CoefDisplay *disp = coef_disps[i];
      assert( disp );
      
      float stepsize = 1.0f;
      if( m_cal && m_cal->num_channels() > 4 )
      {
        const size_t nchannel = m_cal->num_channels();
        switch( m_cal->type() )
        {
          case SpecUtils::EnergyCalType::Polynomial:
          case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
            stepsize = 1.0f / std::pow(nchannel,i);
            break;
            
          case SpecUtils::EnergyCalType::FullRangeFraction:
            stepsize = 1.0;
            break;
            
          case SpecUtils::EnergyCalType::InvalidEquationType:
          case SpecUtils::EnergyCalType::LowerChannelEdge:
            stepsize = 0.0f;
            break;
        }//switch( m_coeffEquationType )
      }//if( valid calibration )
      
      disp->m_value->setSingleStep( stepsize );
    }//for( int i = 0; i < sm_numCoefs; ++i )
     */
    
  }//updateToGui(...)
  
  
};//class CalDisplay

const size_t CalDisplay::sm_min_coef_display = 4;


}//namespace


EnergyCalTool::EnergyCalTool( InterSpec *viewer, PeakModel *peakModel, WContainerWidget *parent )
: WContainerWidget( parent ),
  m_interspec( viewer ),
  m_peakModel( peakModel ),
  m_tallLayoutContent( nullptr ),
  m_peakTable( nullptr ),
  m_specTypeMenu( nullptr ),
  m_specTypeMenuStack( nullptr ),
  m_detectorMenu{ nullptr },
  m_calInfoDisplayStack( nullptr ),
  m_noCalTxt( nullptr ),
  m_moreActionsColumn( nullptr ),
  m_applyToColumn( nullptr ),
  m_detColumn( nullptr ),
  m_detColLayout( nullptr ),
  m_calColumn( nullptr ),
  m_peakTableColumn( nullptr ),
  m_layout( nullptr ),
  m_applyToCbs{ nullptr },
  m_moreActions{ nullptr },
  m_fitCalBtn( nullptr ),
  m_lastGraphicalRecal( 0 ),
  m_lastGraphicalRecalType( EnergyCalGraphicalConfirm::NumRecalTypes ),
  m_lastGraphicalRecalEnergy( -999.0f ),
  m_graphicalRecal( nullptr )
{
  wApp->useStyleSheet( "InterSpec_resources/EnergyCalTool.css" );
  
  addStyleClass( "EnergyCalTool" );
  
  initWidgets( EnergyCalTool::LayoutType::Wide );
}


void EnergyCalTool::initWidgets( EnergyCalTool::LayoutType layoutType )
{
  const bool wide = (layoutType == LayoutType::Wide);
  
  if( (wide && !m_tallLayoutContent && m_layout) || (!wide && m_tallLayoutContent) )
    return;
  
  if( m_graphicalRecal )
  {
    AuxWindow::deleteAuxWindow( m_graphicalRecal );
    m_graphicalRecal = nullptr;
  }//if( m_graphicalRecal )
  
  
  
  if( wide )
  {
    removeStyleClass( "TallEnergyCal" );
    if( m_tallLayoutContent )
      delete m_tallLayoutContent;
    m_tallLayoutContent = nullptr;
    
    m_layout = new WGridLayout( this );
  }else
  {
    addStyleClass( "TallEnergyCal" );
    if( m_layout )
      delete m_layout;
    
    m_tallLayoutContent = new WContainerWidget( this );
    m_layout = new WGridLayout( m_tallLayoutContent );
  }//if( wide ) / else
  
  
  // \TODO: null out all the other memener variables
  m_peakTable = nullptr;
  m_specTypeMenu = nullptr;
  m_specTypeMenuStack = nullptr;
  for( auto &m : m_detectorMenu )
    m = nullptr;
  m_calInfoDisplayStack = nullptr;
  m_noCalTxt = nullptr;
  m_moreActionsColumn = nullptr;
  m_applyToColumn = nullptr;
  m_detColumn = nullptr;
  m_detColLayout = nullptr;
  m_calColumn = nullptr;
  m_peakTableColumn = nullptr;
  for( auto &m : m_applyToCbs )
    m = nullptr;
  for( auto &m : m_moreActions )
    m = nullptr;
  m_fitCalBtn = nullptr;
  
  
  
  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_interspec );
  
  m_layout->setContentsMargins( 0, 0, 0, 0 );
  m_layout->setVerticalSpacing( 0 );
  m_layout->setHorizontalSpacing( 0 );
  
  m_noCalTxt = new WText( "No spectrum loaded" );
  m_noCalTxt->addStyleClass( "NoCalContentTxt" );
  if( wide )
    m_layout->addWidget( m_noCalTxt, 0, 0, AlignmentFlag::AlignCenter | AlignmentFlag::AlignMiddle );
  else
    m_layout->addWidget( m_noCalTxt, 0, 0, 1, 2, AlignmentFlag::AlignCenter | AlignmentFlag::AlignMiddle );
  
  //Create the more actions column...
  m_moreActionsColumn = new WContainerWidget();
  m_moreActionsColumn->addStyleClass( "CalColumn MoreActionCol" );
  if( wide )
    m_layout->addWidget( m_moreActionsColumn, 0, 1 );
  else
    m_layout->addWidget( m_moreActionsColumn, 2, 0 );
  
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
  collayout->setRowStretch( 1, 1 );
  
  WContainerWidget *moreActionsList = new WContainerWidget( moreActionsDiv );
  moreActionsList->addStyleClass( "MoreActionsMenuList" );
  moreActionsList->setList( true );
  
  
  for( MoreActionsIndex index = static_cast<MoreActionsIndex>(0);
      index < MoreActionsIndex::NumMoreActionsIndex;
      index = MoreActionsIndex(static_cast<int>(index) + 1) )
  {
    const char *label = "", *tooltip = nullptr;
    switch( index )
    {
      case MoreActionsIndex::Linearize:
        label = "Linearize...";
        tooltip = "Linearizes spectra so that each energy channel has the same width.";
        break;
        
      case MoreActionsIndex::Truncate:
        label = "Truncate Energy...";
        tooltip = "Truncates the energy range of the spectrum by discarding data channels.";
        break;
        
      case MoreActionsIndex::CombineChannels:
        label = "Combine Channels...";
        tooltip = "Combines energy channels together so spectrum will have less data channels";
        break;
        
      case MoreActionsIndex::ConvertToFrf:
        label = "To FRF...";
        tooltip = "Converts the energy calibration type to Full Range Fraction";
        break;
        
      case MoreActionsIndex::ConvertToPoly:
        label = "To Polynomial...";
        tooltip = "Converts the energy calibration type to Polynomial.";
        break;
        
      case MoreActionsIndex::MultipleFilesCal:
        label = "Multi File Cal...";
        tooltip = "Allows using multiple spectra with peaks fit, from one or more spectrum files,"
                  " to perform an energy calibration.  Can especially be useful for"
                  " lower-resolution detectors where each spectrum has a different test source";
        break;
        
      case MoreActionsIndex::NumMoreActionsIndex:
        assert(0);
        break;
    }//switch( index )
    
    
    
    WContainerWidget *holder = new WContainerWidget( moreActionsList );
    m_moreActions[static_cast<int>(index)] = new WAnchor( WLink(), label, holder );
    m_moreActions[static_cast<int>(index)]->clicked().connect( boost::bind(&EnergyCalTool::moreActionBtnClicked, this, index) );
    
    if( tooltip )
      HelpSystem::attachToolTipOn( holder, tooltip, showToolTips );
  }//for( loop over more actions )
  
  WContainerWidget *btndiv = new WContainerWidget();
  btndiv->addStyleClass( "BtmBtnDiv" );
  collayout->addWidget( btndiv, 2, 0 );
  
  auto helpBtn = new WContainerWidget( btndiv );
  helpBtn->addStyleClass( "Wt-icon ContentHelpBtn" );
  helpBtn->clicked().connect( boost::bind( &HelpSystem::createHelpWindow, "energy-calibration" ) );
  
  m_fitCalBtn = new WPushButton( "Fit Coeffs", btndiv );
  m_fitCalBtn->addStyleClass( "FitCoefBtn" );
  m_fitCalBtn->clicked().connect( this, &EnergyCalTool::fitCoefficients );
  m_fitCalBtn->setDisabled( true );
  
  
  HelpSystem::attachToolTipOn( m_fitCalBtn, "Uses the expected energy of photopeaks "
  "associated with the fit peaks to fit for the energy coefficients.  This button is disabled if no"
  " coefficients are selected to fit for, or less peaks than coefficents are selected.",
                              showToolTips );
  
  
  
  // Create the "Apply To" column that determines what to apply changes to
  m_applyToColumn = new WContainerWidget();
  m_applyToColumn->addStyleClass( "CalColumn ApplyToCol" );
  if( wide )
    m_layout->addWidget( m_applyToColumn, 0, 2 );
  else
    m_layout->addWidget( m_applyToColumn, 2, 1 );
  
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
  
  WAnimation animation(Wt::WAnimation::Fade, Wt::WAnimation::Linear, 200);

  
  // Create the "Coefficients" column that show the polynomial/FRF coefficents.
  m_calColumn = new WContainerWidget();
  m_calColumn->addStyleClass( "CalColumn CoefColumn" );
  if( wide )
    m_layout->addWidget( m_calColumn, 0, 3 );
  else
    m_layout->addWidget( m_calColumn, 1, 0, 1, 2 );
  
  
  collayout = new WGridLayout( m_calColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  collayout->setRowStretch( 1, 1 );
  
  if( wide )
    collayout->setColumnStretch( 1, 1 );
  
  header = new WText( "Calibration Coefficents" );
  header->addStyleClass( "ColHeader" );
  
  collayout->addWidget( header, 0, 0, 1, 2 );
  //collayout->addWidget( m_calInfoDisplayStack, 1, 0 );
  
  
  // Create the "Detector" column that determines which coefficients to show
  m_detColumn = new WContainerWidget();
  m_detColumn->addStyleClass( "DetCol" );
  collayout->addWidget( m_detColumn, 1, 0 );
  
  m_detColLayout = new WGridLayout( m_detColumn );
  m_detColLayout->setContentsMargins( 0, 0, 0, 0 );
  m_detColLayout->setVerticalSpacing( 0 );
  m_detColLayout->setHorizontalSpacing( 0 );
  
  auto detheader = new WText( "Detector" );
  detheader->setInline( false );
  detheader->addStyleClass( "DetHdr Wt-itemview Wt-header Wt-label" );
  //detheader->resize( WLength::Auto, WLength(20,WLength::Unit::Pixel) );
  //collayout->addWidget( detheader, 0, 0 );
  m_detColLayout->addWidget( detheader, 0, 0  );
  m_detColLayout->setRowStretch( 2, 1 );
  
  
  // Create the "Cal Peaks" table
  m_peakTableColumn = new WContainerWidget();
  m_peakTableColumn->addStyleClass( "CalColumn PeakTableCol" );
  if( wide )
  {
    m_layout->addWidget( m_peakTableColumn, 0, 4 );
    m_layout->setColumnStretch( 4, 1 );
  }else
  {
    m_layout->addWidget( m_peakTableColumn, 3, 0, 1, 2 );
  }

  collayout = new WGridLayout( m_peakTableColumn );
  collayout->setContentsMargins( 0, 0, 0, 0 );
  collayout->setVerticalSpacing( 0 );
  collayout->setHorizontalSpacing( 0 );
  if( wide )
    collayout->setRowStretch( 1, 1 );
  
  header = new WText( "Calibration Peaks" );
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
  
  m_peakTable->setColumnHidden( PeakModel::kUseForCalibration, false );
  m_peakTable->setColumnHidden( PeakModel::kMean, false );
  m_peakTable->setColumnHidden( PeakModel::kIsotope, false );
  m_peakTable->setColumnHidden( PeakModel::kPhotoPeakEnergy, false );
  m_peakTable->setColumnHidden( PeakModel::kDifference, false );
  
  
  m_peakTable->setColumnWidth( PeakModel::kUseForCalibration, WLength(3.7, WLength::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kMean, WLength(4.5, WLength::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kIsotope, WLength(4.5, WLength::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kPhotoPeakEnergy, WLength(6.25, WLength::FontEm) );
  m_peakTable->setColumnWidth( PeakModel::kDifference, WLength(5, WLength::FontEm) );
  
  
  
  WItemDelegate *dblDelagate = new WItemDelegate( m_peakTable );
  dblDelagate->setTextFormat( "%.2f" );
  m_peakTable->setItemDelegateForColumn( PeakModel::kMean, dblDelagate );
  
  PhotopeakDelegate *nuclideDelegate = new PhotopeakDelegate( PhotopeakDelegate::NuclideDelegate, true, m_peakTable );
  m_peakTable->setItemDelegateForColumn( PeakModel::kIsotope, nuclideDelegate );
  
  PhotopeakDelegate *photopeakDelegate = new PhotopeakDelegate( PhotopeakDelegate::GammaEnergyDelegate, true, m_peakTable );
  m_peakTable->setItemDelegateForColumn( PeakModel::kPhotoPeakEnergy, photopeakDelegate );
  
  m_peakModel->dataChanged().connect( this, &EnergyCalTool::updateFitButtonStatus );
  m_peakModel->rowsRemoved().connect( this, &EnergyCalTool::updateFitButtonStatus );
  m_peakModel->rowsInserted().connect( this, &EnergyCalTool::updateFitButtonStatus );
  m_peakModel->layoutChanged().connect( this, &EnergyCalTool::updateFitButtonStatus );
  
  m_interspec->displayedSpectrumChanged().connect( boost::bind( &EnergyCalTool::displayedSpectrumChanged, this, _1, _2, _3, _4 ) );
  
  m_renderFlags |= EnergyCalToolRenderFlags::FullGuiUpdate;
  scheduleRender();
}//void initWidgets( EnergyCalTool::LayoutType layout )


void EnergyCalTool::setWideLayout()
{
  initWidgets( EnergyCalTool::LayoutType::Wide );
}//void setWideLayout()


void EnergyCalTool::setTallLayout()
{
  initWidgets( EnergyCalTool::LayoutType::Tall );
}//void setTallLayout()


EnergyCalTool::~EnergyCalTool()
{
}
  

set<string> EnergyCalTool::gammaDetectorsForDisplayedSamples( const SpecUtils::SpectrumType type )
{
  //We want the names of just the detectors that have gamma calibration information, of the
  //  currently displayed samples.
  //  We will assume the first Measurement for a given named detector will have gamma data if
  //  any of the Measurements from that detector will.
  //  \TODO: evaluate if that is true.
  
  auto meas = m_interspec->measurment( type );
  
  if( !meas )
    return {};
  
  const vector<string> &detnames = meas->gamma_detector_names();
  const vector<string> displayedDets = m_interspec->detectorsToDisplay(type);
  const set<int> &samples = m_interspec->displayedSamples(type);
  
  set<string> detectors, nongammadets;
  
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
  
  return detectors;
}//set<string> gammaDetectorsForDisplayedSamples( const SpecUtils::SpectrumType type )


vector<MeasToApplyCoefChangeTo> EnergyCalTool::measurementsToApplyCoeffChangeTo()
{
  std::vector<MeasToApplyCoefChangeTo> answer;
  
  //Lets loop over spectrum types (Foreground, Background, Secondary), and decide if we should apply
  //  changes to that file, and if so, decide which sample numbers/detectors
  const SpecUtils::SpectrumType spectypes[3] = {
    SpecUtils::SpectrumType::Foreground,
    SpecUtils::SpectrumType::Background,
    SpecUtils::SpectrumType::SecondForeground
  };
  
  for( const auto spectype : spectypes )
  {
    auto meas = m_interspec->measurment( spectype );
    if( !meas )
      continue;
    
    switch( spectype )
    {
      case SpecUtils::SpectrumType::Foreground:
        if( !m_applyToCbs[ApplyToCbIndex::ApplyToForeground]->isChecked() )
          continue;
        break;
        
      case SpecUtils::SpectrumType::SecondForeground:
        if( !m_applyToCbs[ApplyToCbIndex::ApplyToSecondary]->isChecked() )
          continue;
        break;
        
      case SpecUtils::SpectrumType::Background:
        if( !m_applyToCbs[ApplyToCbIndex::ApplyToBackground]->isChecked() )
          continue;
        break;
    }//switch( spectype )
    
    //If we're here, we should apply changes to this spectrum type.
    
    //It could be that the background SpecFile is the same as the foreground, so check if we
    //  already have an entry for this file in answer, and if so, use it.  We dont want duplicate
    //  entries for SpecFiles since this could cause us to maybe move peaks multiple times or
    //  something
    //  \TODO: if the backgeround and foreground use different detectors (no overlap) or different
    //         sample numbers (no overlap) should return multiple entries for the SpecFile to handle
    //         the edge-case correctly
    MeasToApplyCoefChangeTo *changes = nullptr;
    for( size_t i = 0; !changes && i < answer.size(); ++i )
      changes = (answer[i].meas == meas) ? &(answer[i]) : changes;
    
    if( !changes )
    {
      //We havent seen this SpecFile yet, create an entry in answer for it
      answer.emplace_back();
      changes = &(answer.back()); //C++14 returns a reference to the emplaced object.
      changes->meas = meas;
    }
    
    assert( changes );
    
    //m_applyToCbs[ApplyToCbIndex::ApplyToDisplayedDetectors]->parent()->isHidden();
    const set<string> detectors = gammaDetectorsForDisplayedSamples(spectype);
    const vector<string> displayed_dets = m_interspec->detectorsToDisplay(spectype);
    
    bool displayingAllDets = true;
    for( const auto &det : detectors )
    {
      if( std::find(begin(displayed_dets), end(displayed_dets), det) == end(displayed_dets) )
        displayingAllDets = false;
    }
    
    if( displayingAllDets )
    {
      //We will insert detector names since different SpecType from the same file could have
      //  different detector names available.
      // \TODO: if InterSpec class is upgraded to select detector by SpecType, we will have to
      //        upgrade this part of the code.
      for( const auto &det : detectors )
        changes->detectors.insert( det );
    }else
    {
      const auto applyToAllCb = m_applyToCbs[ApplyToCbIndex::ApplyToAllDetectors];
      const auto applyToDisplayedCb = m_applyToCbs[ApplyToCbIndex::ApplyToDisplayedDetectors];
      
      const bool toAll = (applyToAllCb->parent()
                          && !applyToAllCb->parent()->isHidden()
                          && applyToAllCb->isChecked());
      const bool toDisplayed = (applyToDisplayedCb->parent()
                                && !applyToDisplayedCb->parent()->isHidden()
                                && applyToDisplayedCb->isChecked());
      if( toAll == toDisplayed )
      {
        cerr << "EnergyCalTool::measurementsToApplyCoeffChangeTo:"
                " got (toAll == toDisplayed) which shouldnt have happended" << endl;
      }
        
      if( toAll || (toAll == toDisplayed) )
      {
        for( const auto &det : detectors )
          changes->detectors.insert( det );
      }else
      {
        for( const auto &dispdet : displayed_dets )
        {
          if( detectors.count(dispdet) )
            changes->detectors.insert( dispdet );
        }
      }//if( apply to all detectors ) / else ( only displayed detectors )
    }//if( displayingAllDets ) / else
    
    
    const auto toAllSamplesCb = m_applyToCbs[ApplyToCbIndex::ApplyToAllSamples];
    const bool onlyDispSamples = (toAllSamplesCb->parent()
                                  && !toAllSamplesCb->parent()->isHidden()
                                  && !toAllSamplesCb->isChecked());
    
    if( onlyDispSamples )
    {
      const set<int> &displayed_samples = m_interspec->displayedSamples(spectype);
      for( const int sample : displayed_samples )
        changes->sample_numbers.insert( sample );
    }else
    {
      changes->sample_numbers = meas->sample_numbers();
    }
  }//for( const auto spectype : spectypes )
  
  
  return answer;
}//std::vector<MeasToApplyCoefChangeTo> measurementsToApplyCoeffChangeTo()


void EnergyCalTool::applyCalChange( std::shared_ptr<const SpecUtils::EnergyCalibration> disp_prev_cal,
                                    std::shared_ptr<const SpecUtils::EnergyCalibration> new_disp_cal,
                                    const bool isOffsetOnly )
{
  using namespace SpecUtils;
  
  assert( disp_prev_cal && disp_prev_cal->valid() );
  assert( new_disp_cal && new_disp_cal->valid() );
  
  const auto forgrnd = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  const auto backgrnd = m_interspec->measurment(SpecUtils::SpectrumType::Background);
  const auto secgrnd = m_interspec->measurment(SpecUtils::SpectrumType::SecondForeground);
  
  const set<int> &foresamples = m_interspec->displayedSamples( SpectrumType::Foreground );
  
  
  // Create a cache of modified calibration both to save time/memory, but also keep it so previous
  //  samples that share a energy calibration will continue to do so (if possible based on what user
  //  wanted calibration applied to).  Also, we wont set any new calibrations until we know all
  //  updated calibrations and peaks are valid
  //Note: we could take this oppritunity to share calibration across SpecFile objects by not just
  //      comparing pointers, but also the actual EnergyCalibration object.  But for now we'll
  //      skip this to avoid trouble, and it isnt clear that it would actually be overall beneficial
  map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> updated_cals;
  
  // We will store updated peaks and not set any of them until we know all the energy calibrations
  //  and peak shifts were sucessfully done.
  map<shared_ptr<deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_peaks;
  
  const vector<MeasToApplyCoefChangeTo> changemeas = measurementsToApplyCoeffChangeTo();
  
  //We will loop over the changes to apply twice.  Once to calculate new calibrations, and make sure
  //  they are valid, then a second time to actually set them.  If a new calibration is invalid,
  //  an exception will be thrown so we will catch that.
  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    assert( change.meas );
    
    string dbgmsg = "For '" + change.meas->filename() + "' will apply changes to Detectors: {";
    for( auto iter = begin(change.detectors); iter != end(change.detectors); ++iter )
      dbgmsg += (iter==begin(change.detectors) ? "" : ",") + (*iter);
    dbgmsg += "} and Samples: {";
    for( auto iter = begin(change.sample_numbers); iter != end(change.sample_numbers); ++iter )
      dbgmsg += (iter==begin(change.sample_numbers) ? "" : ",") + std::to_string(*iter);
    dbgmsg += "}";
    cout << dbgmsg << endl;
    
    try
    {
      for( const int sample : change.sample_numbers )
      {
        for( const string &detname : change.detectors )
        {
          auto m = change.meas->measurement( sample, detname );
          if( !m || m->num_gamma_channels() <= 4 )
            continue;
          
          const auto meas_old_cal = m->energy_calibration();
          assert( meas_old_cal );
          
          if( !meas_old_cal || !meas_old_cal->valid() )
            continue;
          
          //If we have already computed the new calibration for a EnergyCalibration object, lets not
          //  re-due it.
          if( updated_cals.count(meas_old_cal) )
            continue;
          
          shared_ptr<const SpecUtils::EnergyCalibration> new_meas_cal;
          if( meas_old_cal == disp_prev_cal )
          {
            new_meas_cal = new_disp_cal;
          }else if( isOffsetOnly )
          {
            const vector<float> &new_disp_coefs = new_disp_cal->coefficients();
            const vector<float> &prev_disp_coefs = disp_prev_cal->coefficients();
            vector<float> new_coefs = meas_old_cal->coefficients();
            new_coefs[0] += (new_disp_coefs[0] - prev_disp_coefs[0]);
            
            auto cal = make_shared<EnergyCalibration>();
            switch( meas_old_cal->type() )
            {
              case SpecUtils::EnergyCalType::Polynomial:
              case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
                cal->set_polynomial( meas_old_cal->num_channels(), new_coefs, meas_old_cal->deviation_pairs() );
                break;
                
              case SpecUtils::EnergyCalType::FullRangeFraction:
                cal->set_full_range_fraction( meas_old_cal->num_channels(), new_coefs, meas_old_cal->deviation_pairs() );
                break;
                
              case SpecUtils::EnergyCalType::LowerChannelEdge:
                cal->set_lower_channel_energy( meas_old_cal->num_channels(), new_coefs ); //eh, whatever
                break;
                
              case SpecUtils::EnergyCalType::InvalidEquationType:
                assert( 0 );
                break;
            }//switch( meas_old_cal->type() )
            
            new_meas_cal = cal;
          }else
          {
            new_meas_cal = EnergyCal::propogate_energy_cal_change( disp_prev_cal, new_disp_cal, meas_old_cal );
          }
          assert( new_meas_cal && new_meas_cal->valid() );
          updated_cals[meas_old_cal] = new_meas_cal;
        }//for( const string &detname : change.detectors )
      }//for( loop over sample numbers )
    }catch( std::exception &e )
    {
      string msg = "Calibration change made a energy calibration become invalid";
      if( (backgrnd && (backgrnd != change.meas)) || (secgrnd && (secgrnd != change.meas)) )
      {
        if( change.meas == forgrnd )
          msg += " for the foreground";
        else if( change.meas == backgrnd )
          msg += " for the background";
        else if( change.meas == secgrnd )
          msg += " for the secondary spectrum";
      }//if( it is necassry to say which spectrum had the error )
      
      msg += ".  Error: ";
      msg += e.what();
      
      throw runtime_error( msg );
    }//try catch
    
    
    // Now go through and translate the peaks, but we wont actually update them to the SpecMeas
    //  until we know we can update all the peaks
    const set<set<int>> peaksamples = change.meas->sampleNumsWithPeaks();
    const vector<string> detnamesv( begin(change.detectors), end(change.detectors) );
    
    for( const set<int> &samples : peaksamples )
    {
      //If there is any overlap between 'samples' and 'change.sample_numbers', then apply the change
      //  Note: this isnt correct, but I cant think of a better solution at the moment.
      auto oldpeaks = change.meas->peaks(samples);
      auto oldcal = change.meas->suggested_sum_energy_calibration( samples, detnamesv );
      
      if( !oldpeaks || oldpeaks->empty()
         || !oldcal || (oldcal->type() == EnergyCalType::InvalidEquationType) )
      {
        if( !oldpeaks || !oldcal || !oldcal->valid() )
          cerr << "Failed to get peaks or oldcal!" << endl; //just for development
        continue;
      }
      
      auto newcal = updated_cals[oldcal];
      if( !newcal || !newcal->valid() )
      {
        cerr << "Failed to get newcal for peaks shift!" << endl; //just for development, shouldnt happen I dont think
        continue;
      }
      
      if( oldcal == newcal )
      {
        cerr << __func__ <<  ": oldcal == newcal - skipping shifting peak" << endl;
        continue;
      }
      
      try
      {
        auto newpeaks = EnergyCal::translatePeaksForCalibrationChange( *oldpeaks, oldcal, newcal );
        updated_peaks[oldpeaks] = newpeaks;
      }catch( std::exception &e )
      {
        string msg = "There was an issue translating peaks for this energy change;"
        " not applying change.  Error: " + string(e.what());
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, msg.c_str() );
#endif
        
        throw runtime_error( msg );
      }//try / catch
      
    }//for( const set<int> &samples : peaksampels )
  }//for( const MeasToApplyCoefChangeTo &change : changemeas )
  
  if( updated_cals.find(disp_prev_cal) == end(updated_cals) )
  {
    //Shouldnt ever happen; check is for development
    string msg = "There was an internal error updating energy calibration - energy cal"
    " associated with GUI wasn't updated - energy calibration state is suspect";
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, msg.c_str() );
#endif
    
    m_interspec->logMessage( msg, "", 3 );
    
#if( PERFORM_DEVELOPER_CHECKS )
    assert( 0 );
#endif
  }//if( updated_cals.find(disp_prev_cal) == end(updated_cals) )
  
  
  
  // Now go through and actually set the energy calibrations; they should all be valid and computed,
  //  as should all the shifted peaks.
  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    assert( change.meas );
    for( const int sample : change.sample_numbers )
    {
      for( const string &detname : change.detectors )
      {
        auto m = change.meas->measurement( sample, detname );
        if( !m || m->num_gamma_channels() <= 4 )
          continue;
        
        const auto measoldcal = m->energy_calibration();
        assert( measoldcal );
        
        auto iter = updated_cals.find( measoldcal );
        if( iter == end(updated_cals) )
        {
          //Shouldnt ever happen
          string msg = "There was an internal error updating energy calibration - precomputed"
          " calibration couldnt be found - energy calibation will not be fully updated";
#if( PERFORM_DEVELOPER_CHECKS )
          log_developer_error( __func__, msg.c_str() );
#endif
          
          m_interspec->logMessage( msg, "", 3 );
          assert( 0 );
          continue;
        }//if( we havent already computed a new energy cal )
        
        assert( iter->second );
        assert( iter->second->num_channels() == m->num_gamma_channels() );
        
        change.meas->set_energy_calibration( iter->second, m );
      }//for( loop over detector names )
    }//for( loop over sample numbers )
    
    
    //Now actually set the updated peaks
    const set<set<int>> peaksamples = change.meas->sampleNumsWithPeaks();
    
    for( const set<int> &samples : peaksamples )
    {
      auto oldpeaks = change.meas->peaks(samples);
      if( !oldpeaks )
        continue;
      
      const auto pos = updated_peaks.find(oldpeaks);
      if( pos == end(updated_peaks) )
      {
        if( !oldpeaks->empty() )
          cerr << "Couldnt find an expected entry in updated_peaks" << endl;
        continue;
      }
      
      change.meas->setPeaks( pos->second, samples );
      if( m_peakModel && (change.meas == forgrnd) && (samples == foresamples) )
        m_peakModel->setPeakFromSpecMeas(forgrnd, foresamples);
    }//for( const set<int> &samples : peaksampels )
  }//for( loop over SpecFiles for change )
  
  
  // \TODO: Set an undu point for each measurement, etc.
  
  m_interspec->refreshDisplayedCharts();
  doRefreshFromFiles();
}//applyCalChange(...)


void EnergyCalTool::addDeviationPair( const std::pair<float,float> &new_pair )
{
  using namespace SpecUtils;
  
  const auto forgrnd = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  const auto backgrnd = m_interspec->measurment(SpecUtils::SpectrumType::Background);
  const auto secgrnd = m_interspec->measurment(SpecUtils::SpectrumType::SecondForeground);
  
  const set<int> &foresamples = m_interspec->displayedSamples( SpectrumType::Foreground );
  
  // We will calculate all new energy calibrations, to make sure we actually can, befor setting them
  map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> updated_cals;
  
  // We will also pre-calculate updated peaks
  map<shared_ptr<deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_peaks;
  
  const vector<MeasToApplyCoefChangeTo> changemeas = measurementsToApplyCoeffChangeTo();
  
  // Do first loop to calculate new calibrations
  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    assert( change.meas );
    
    try
    {
      for( const int sample : change.sample_numbers )
      {
        for( const string &detname : change.detectors )
        {
          auto m = change.meas->measurement( sample, detname );
          if( !m || m->num_gamma_channels() <= 4 )
            continue;
          
          const auto old_cal = m->energy_calibration();
          assert( old_cal );
          
          if( !m || (m->num_gamma_channels() <= 4) || !old_cal || !old_cal->valid()
             || (old_cal->type() == EnergyCalType::LowerChannelEdge) )
            continue;
          
          //If we have already computed the new calibration for a EnergyCalibration object, lets not
          //  re-due it.
          if( updated_cals.count(old_cal) )
            continue;
          
          const size_t nchannel = old_cal->num_channels();
          const vector<float> &coefs = old_cal->coefficients();
          vector<pair<float,float>> dev_pairs = old_cal->deviation_pairs();
          if( dev_pairs.empty() )
            dev_pairs.push_back( {0.0f, 0.0f} );
          dev_pairs.push_back( new_pair );
          
          std::sort( begin(dev_pairs), end(dev_pairs) );
          
          shared_ptr<EnergyCalibration> new_cal = make_shared<EnergyCalibration>();
          
          switch( old_cal->type() )
          {
            case SpecUtils::EnergyCalType::Polynomial:
            case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
              new_cal->set_polynomial( nchannel, coefs, dev_pairs );
              break;
              
            case SpecUtils::EnergyCalType::FullRangeFraction:
              new_cal->set_full_range_fraction( nchannel, coefs, dev_pairs );
              break;
              
            case SpecUtils::EnergyCalType::LowerChannelEdge:
            case SpecUtils::EnergyCalType::InvalidEquationType:
              assert( 0 );
              break;
          }//switch( meas_old_cal->type() )
          
          assert( new_cal && new_cal->valid() );
          updated_cals[old_cal] = new_cal;
        }//for( const string &detname : change.detectors )
      }//for( loop over sample numbers )
    }catch( std::exception &e )
    {
      string msg = "Adding deviation pair made energy calibration become invalid";
      if( (backgrnd && (backgrnd != change.meas)) || (secgrnd && (secgrnd != change.meas)) )
      {
        if( change.meas == forgrnd )
          msg += " for the foreground";
        else if( change.meas == backgrnd )
          msg += " for the background";
        else if( change.meas == secgrnd )
          msg += " for the secondary spectrum";
      }//if( it is necassry to say which spectrum had the error )
      
      msg += ".  Error: ";
      msg += e.what();
      
      throw runtime_error( msg );
    }//try catch
    
    
    // Now go through and translate the peaks, but we wont actually update them to the SpecMeas
    //  until we know we can update all the peaks
    const set<set<int>> peaksamples = change.meas->sampleNumsWithPeaks();
    const vector<string> detnamesv( begin(change.detectors), end(change.detectors) );
    
    for( const set<int> &samples : peaksamples )
    {
      //If there is any overlap between 'samples' and 'change.sample_numbers', then apply the change
      //  Note: this isnt correct, but I cant think of a better solution at the moment.
      auto oldpeaks = change.meas->peaks(samples);
      auto oldcal = change.meas->suggested_sum_energy_calibration( samples, detnamesv );
      
      if( !oldpeaks || oldpeaks->empty()
         || !oldcal || (oldcal->type() == EnergyCalType::InvalidEquationType) )
      {
        if( !oldpeaks || !oldcal || !oldcal->valid() )
          cerr << "Failed to get peaks or oldcal!" << endl; //just for development
        continue;
      }
      
      auto newcal = updated_cals[oldcal];
      if( !newcal || !newcal->valid() )
      {
        cerr << "Failed to get newcal for peaks shift!" << endl; //just for development, shouldnt happen I dont think
        continue;
      }
      
      if( oldcal == newcal )
      {
        cerr << __func__ <<  ": oldcal == newcal - skipping shifting peak" << endl;
        continue;
      }
      
      try
      {
        auto newpeaks = EnergyCal::translatePeaksForCalibrationChange( *oldpeaks, oldcal, newcal );
        updated_peaks[oldpeaks] = newpeaks;
      }catch( std::exception &e )
      {
        string msg = "There was an issue translating peaks for this energy change;"
        " not applying change.  Error: " + string(e.what());
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, msg.c_str() );
#endif
        
        throw runtime_error( msg );
      }//try / catch
      
    }//for( const set<int> &samples : peaksampels )
  }//for( const MeasToApplyCoefChangeTo &change : changemeas )
  
  
  // Now go through and actually set the energy calibrations; they should all be valid and computed,
  //  as should all the shifted peaks.
  for( const MeasToApplyCoefChangeTo &change : changemeas )
  {
    assert( change.meas );
    for( const int sample : change.sample_numbers )
    {
      for( const string &detname : change.detectors )
      {
        auto m = change.meas->measurement( sample, detname );
        if( !m || m->num_gamma_channels() <= 4 )
          continue;
        
        const auto measoldcal = m->energy_calibration();
        assert( measoldcal );
        
        auto iter = updated_cals.find( measoldcal );
        if( iter == end(updated_cals) )
        {
          //Shouldnt ever happen
          string msg = "There was an internal error updating energy calibration - precomputed"
          " calibration couldnt be found - energy calibation will not be fully updated";
#if( PERFORM_DEVELOPER_CHECKS )
          log_developer_error( __func__, msg.c_str() );
#endif
          
          m_interspec->logMessage( msg, "", 3 );
          assert( 0 );
          continue;
        }//if( we havent already computed a new energy cal )
        
        assert( iter->second );
        assert( iter->second->num_channels() == m->num_gamma_channels() );
        
        change.meas->set_energy_calibration( iter->second, m );
      }//for( loop over detector names )
    }//for( loop over sample numbers )
    
    
    //Now actually set the updated peaks
    const set<set<int>> peaksamples = change.meas->sampleNumsWithPeaks();
    
    for( const set<int> &samples : peaksamples )
    {
      auto oldpeaks = change.meas->peaks(samples);
      if( !oldpeaks || oldpeaks->empty() )
        continue;
      
      const auto pos = updated_peaks.find(oldpeaks);
      if( pos == end(updated_peaks) )
      {
        if( !oldpeaks->empty() )
          cerr << "Couldnt find an expected entry in updated_peaks" << endl;
        continue;
      }
      
      change.meas->setPeaks( pos->second, samples );
      if( m_peakModel && (change.meas == forgrnd) && (samples == foresamples) )
        m_peakModel->setPeakFromSpecMeas(forgrnd, foresamples);
    }//for( const set<int> &samples : peaksampels )
  }//for( loop over SpecFiles for change )
  
  
  // \TODO: Set an undu point for each measurement, etc.
  
  m_interspec->refreshDisplayedCharts();
  refreshGuiFromFiles();
}//void addDeviationPair( const std::pair<float,float> &new_pair );



void EnergyCalTool::userChangedCoefficient( const size_t coefnum, EnergyCalImp::CalDisplay *display )
{
  cout << "EnergyCalTool::userChangedCoefficient" << endl;
  using namespace SpecUtils;
  assert( coefnum < 10 );  //If we ever allow lower channel energy adjustment this will need to be removed
  
  shared_ptr<const EnergyCalibration> disp_prev_cal = display->lastSetCalibration();
  if( !disp_prev_cal )
  {
    cerr << "unexpected error getting updaettd energy calibration coefficents" << endl;
    m_interspec->logMessage( "Unexpected error retrieving previous calibration coefficients - not applying changes", "", 2 );
    doRefreshFromFiles();
    return;
  }//if( !disp_prev_cal )
  
  
  {// Begin check to make sure the changed energy cal is actually checked for it to be applied to...
    bool willBeAppliedToDisplay = false;
    const std::string &detname = display->detectorName();
    const SpecUtils::SpectrumType type = display->spectrumType();
    std::shared_ptr<SpecMeas> cal_disp_meas = m_interspec->measurment(type);
    assert( cal_disp_meas );
    
    for( const MeasToApplyCoefChangeTo &delta : measurementsToApplyCoeffChangeTo() )
    {
      const shared_ptr<SpecMeas> &meas = delta.meas;
      if( meas != cal_disp_meas )
        continue;
      
      const set<string> &detectors = delta.detectors;
      if( !detectors.count(detname) )
        continue;
      
      // Actually, I think if we're here, we're probably good, but we'll check a little deeper to
      //   make sure check in EnergyCalTool::applyCalChange will be satisfied
      const set<int> &samples = delta.sample_numbers;
      for( const int sample : samples )
      {
        auto m = meas->measurement( sample, detname);
        willBeAppliedToDisplay = (m && (m->energy_calibration() == disp_prev_cal));
        if( willBeAppliedToDisplay )
          break;
      }//for( const int sample : samples )
    
      if( willBeAppliedToDisplay )
        break;
    }//for( loop over changes )
    
    if( !willBeAppliedToDisplay )
    {
      m_interspec->logMessage( "It looks like the energy calibration you changed isn't marked"
                               " for changes to be applied to; please correct that.", "", 2 );
      doRefreshFromFiles();
      return;
    }
  }// End check to make sure the changed energy cal is actually checked for it to be applied to...
  
  
  m_lastGraphicalRecal = 0;
  m_lastGraphicalRecalType = EnergyCalGraphicalConfirm::NumRecalTypes;
  m_lastGraphicalRecalEnergy = -999.0f;
  
  vector<float> dispcoefs = display->displayedCoefficents();
  if( dispcoefs.size() <= coefnum )
    dispcoefs.resize( coefnum+1, 0.0f );
  
  vector<float> prev_disp_coefs = disp_prev_cal->coefficients();
  if( prev_disp_coefs.size() <= coefnum )
    prev_disp_coefs.resize( coefnum+1, 0.0f );
  
  vector<float> new_disp_coefs = prev_disp_coefs;
  new_disp_coefs[coefnum] = dispcoefs[coefnum];
  
  const size_t dispnchannel = disp_prev_cal->num_channels();
  const auto &disp_dev_pairs = disp_prev_cal->deviation_pairs();
  
  shared_ptr<const EnergyCalibration> new_disp_cal;
  try
  {
    auto cal = make_shared<EnergyCalibration>();
    switch( disp_prev_cal->type() )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        cal->set_polynomial( dispnchannel, new_disp_coefs, disp_dev_pairs );
        break;
        
      case SpecUtils::EnergyCalType::FullRangeFraction:
        cal->set_full_range_fraction( dispnchannel, new_disp_coefs, disp_dev_pairs );
        break;
      
      case SpecUtils::EnergyCalType::LowerChannelEdge:
      case SpecUtils::EnergyCalType::InvalidEquationType:
        throw runtime_error( "Invalid calibration type changed?  Something is way wack." );
        break;
    }//switch( disp_prev_cal->type() )
    
    new_disp_cal = cal;
  }catch( std::exception &e )
  {
    display->updateToGui( disp_prev_cal );
    
    string msg = "Calibration change made energy calibration become invalid.  Error: ";
    msg += e.what();
    m_interspec->logMessage( msg, "", 2 );
    
    return;
  }//try / catch to create new_disp_cal
  
  assert( new_disp_cal && new_disp_cal->valid() );
  
  try
  {
    applyCalChange( disp_prev_cal, new_disp_cal, coefnum==0 );
  }catch( std::exception &e )
  {
    display->updateToGui( disp_prev_cal );
    m_interspec->logMessage( e.what(), "", 2 );
  }//try / catch
}//userChangedCoefficient(...)


void EnergyCalTool::userChangedDeviationPair( EnergyCalImp::CalDisplay *display, const int fieldTypeChanged )
{
  using namespace SpecUtils;
  
  const auto field = EnergyCalImp::DeviationPairDisplay::UserFieldChanged( fieldTypeChanged );
  
  switch( field )
  {
    case EnergyCalImp::DeviationPairDisplay::UserFieldChanged::AddedDeviationPair:
      return;  //hasnt been filled out yet; no need to do anything
      
    case EnergyCalImp::DeviationPairDisplay::UserFieldChanged::RemovedDeviationPair:
    case EnergyCalImp::DeviationPairDisplay::UserFieldChanged::EnergyChanged:
    case EnergyCalImp::DeviationPairDisplay::UserFieldChanged::OffsetChanged:
      break;
  };//enum UserFieldChanged

  
  m_lastGraphicalRecal = 0;
  m_lastGraphicalRecalType = EnergyCalGraphicalConfirm::NumRecalTypes;
  m_lastGraphicalRecalEnergy = -999.0f;
  
  const auto forgrnd = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  const set<int> &foresamples = m_interspec->displayedSamples( SpectrumType::Foreground );
  
  const shared_ptr<const EnergyCalibration> old_cal = display->lastSetCalibration();
  const vector<pair<float,float>> old_dev_pairs = old_cal
                                        ? old_cal->deviation_pairs() : vector<pair<float,float>>{};
  const vector<pair<float,float>> new_dev_pairs = display->displayedDeviationPairs();
  
  // After the user clicks to add a new deviation pair, and then fills in one of the fields, if the
  //  other field is not yet filled in, then the GUI wont insert that deviation pair; so for this
  //  case where the deviation pairs havent yet changed, lets skip doing anything yet
  if( new_dev_pairs.size() == old_dev_pairs.size() )
  {
    bool equal = true;
    for( size_t i = 0; equal && (i < new_dev_pairs.size()); ++i )
      equal = ( (fabs(new_dev_pairs[i].first - old_dev_pairs[i].first) < 1.0E-4)
                && (fabs(new_dev_pairs[i].second - old_dev_pairs[i].second) < 1.0E-4) );
    if( equal )
      return;
  }//if( new_dev_pairs.size() == old_dev_pairs.size() )
  
  
  const SpectrumType type = display->spectrumType();
  const std::string &detname = display->detectorName();
  
  auto specfile = m_interspec->measurment( type );
  if( !specfile )  //Shouldnt ever happen
  {
    display->updateToGui( old_cal );
    m_interspec->logMessage( "Internal error retrieveing correct measurement", "", 2 );
    return;
  }

  // We will store updated calibrations and not set any of them until we know they can all be
  //  sucessfully altered
  map<shared_ptr<const EnergyCalibration>,shared_ptr<const EnergyCalibration>> old_to_new_cal;
  
  try
  {
    for( auto &m : specfile->measurements() )
    {
      if( (m->detector_name() != detname) || (m->num_gamma_channels() < 5) )
        continue;
      
      const auto cal = m->energy_calibration();
      if( !cal || !cal->valid() || (cal->type() == EnergyCalType::LowerChannelEdge) )
        continue;
      
      if( old_to_new_cal.find(cal) != end(old_to_new_cal) )
        continue;
      
      const size_t nchannel = cal->num_channels();
      const vector<float> &coefficients = cal->coefficients();
      
      auto new_cal = make_shared<EnergyCalibration>();
      switch( cal->type() )
      {
        case EnergyCalType::Polynomial:
        case EnergyCalType::UnspecifiedUsingDefaultPolynomial:
          new_cal->set_polynomial( nchannel, coefficients, new_dev_pairs );
          break;
          
        case EnergyCalType::FullRangeFraction:
          new_cal->set_full_range_fraction( nchannel, coefficients, new_dev_pairs );
          break;
          
        case EnergyCalType::LowerChannelEdge:
        case EnergyCalType::InvalidEquationType:
          assert( 0 );
          break;
      }//switch( old_cal->type() )
      
      assert( new_cal->valid() );
      old_to_new_cal[cal] = new_cal;
    }//for( auto &m : specfile->measurements() )
  }catch( std::exception &e )
  {
    display->updateToGui( old_cal );
    //display->setDeviationPairsInvalid();
    
    string msg = "Changing the deviation pair made energy calibration become invalid.  Error: ";
    msg += e.what();
    m_interspec->logMessage( msg, "", 2 );
    
    return;
  }//try / catch
  
  // We will store updated peaks and not set any of them until we know all the energy calibrations
  //  and peak shifts were sucessfully done.
  map<shared_ptr<deque<shared_ptr<const PeakDef>>>,deque<shared_ptr<const PeakDef>>> updated_peaks;
  
  //const vector<string> &detnames = specfile->gamma_detector_names();
  const vector<string> &detnames = m_interspec->detectorsToDisplay( type );
  const set<set<int>> samplesWithPeaks = specfile->sampleNumsWithPeaks();
  for( const set<int> &samples : samplesWithPeaks )
  {
    try
    {
      auto dispcal = specfile->suggested_sum_energy_calibration( samples, detnames );
      
      const auto dispcaliter = old_to_new_cal.find(dispcal);
      if( dispcaliter == end(old_to_new_cal) )
        continue;
      
      bool detInSample = false;
      for( auto iter = begin(samples); !detInSample && (iter != end(samples)); ++iter )
      {
        auto m = specfile->measurement( *iter, detname );
        detInSample = !!m;
      }
      
      if( !detInSample )
        continue;
      
      // I *think* that applying the update to deviation pairs shouldnt compound when you then apply
      //  them to another detector since the dispaly energy cal should update, but could there be
      //  any edge-cases?
      auto oldpeaks = specfile->peaks(samples);
      if( !oldpeaks || oldpeaks->empty() )
        continue;
      
    
      auto newpeaks = EnergyCal::translatePeaksForCalibrationChange( *oldpeaks, dispcaliter->first,
                                                                     dispcaliter->second );
      updated_peaks[oldpeaks] = newpeaks;
    }catch( std::exception &e )
    {
      //display->updateToGui( old_cal );
      display->setDeviationPairsInvalid();
      
      string msg = "There was an issue translating peaks for this deviation pair change;"
                   " not applying change.  Error: " + string(e.what());
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, msg.c_str() );
#endif
      m_interspec->logMessage( msg, "", 2 );
      
      return;
    }//try / catch
  }//for( const set<int> &samples : samplesWithPeaks )
  
  display->setDeviationPairsValid();
  
  
  size_t num_updated = 0;
  for( auto &m : specfile->measurements() )
  {
    // I'm a little torn if we should update just the one energy calibration, or all occurances of
    //  the detectors energy calibration.
    //  Maybe we should update according to the checked GUI, but also restrict on name as well?
    //auto cal = m->energy_calibration();
    //if( cal == old_cal )
    if( (m->detector_name() != detname) || (m->num_gamma_channels() < 5) )
      continue;
    
    const auto cal = m->energy_calibration();
    if( !cal || !cal->valid() || (cal->type() == EnergyCalType::LowerChannelEdge) )
      continue;
    
    auto calpos = old_to_new_cal.find(cal);
    if( (calpos == end(old_to_new_cal)) || !calpos->second || !calpos->second->valid() )
    {
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, "Unexpectedly found invalid calibation in old_to_new_cal" );
#endif
      continue;
    }//if( sanity check that new calibration is valid - should always be )
    
    specfile->set_energy_calibration( calpos->second, m );
    ++num_updated;
  }//for( loop over measurements )
  
  if( num_updated == 0 )
  {
    display->updateToGui( old_cal );
    m_interspec->logMessage( "There was an error setting deviation pairs for this detector.", "", 2 );
    return;
  }
  
  
  for( const set<int> &samples : samplesWithPeaks )
  {
    auto oldpeaks = specfile->peaks(samples);
    if( !oldpeaks || oldpeaks->empty() )
      continue;
    
    const auto peakpos = updated_peaks.find(oldpeaks);
    if( peakpos == end(updated_peaks) )
    {
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, "Unexpectedly coudlnt find peaks in updated_peaks" );
#endif
      continue;
    }//if( sanity check that shouldnt ever happen )
      
    
    specfile->setPeaks( peakpos->second, samples );
    if( m_peakModel && (specfile == forgrnd) && (samples == foresamples) )
      m_peakModel->setPeakFromSpecMeas(forgrnd, foresamples);
  }//for( const set<int> &samples : samplesWithPeaks )
  
  const size_t ndets = specfile->gamma_detector_names().size();
  const size_t nsamples = specfile->sample_numbers().size();
  
  if( (ndets > 1) || (nsamples > 1) )
  {
    string msg;
    if( ndets > 1 )
      msg = "Deviation pairs applied only to the '" + detname + "' detector";
    if( nsamples > 1 )
      msg += string(msg.empty() ? "Deviation pairs applied" : ", but") + " to all sample numbers";
    
    int nfiles = 0;
    for( auto t : {0,1,2} )
      nfiles += (m_interspec->measurment(static_cast<SpectrumType>(t)) != specfile);
    
    if( nfiles && !msg.empty() )
    {
      switch( type )
      {
        case SpectrumType::Foreground:       msg += " of foreground"; break;
        case SpectrumType::SecondForeground: msg += " of secondary"; break;
        case SpectrumType::Background:       msg += " of background"; break;
      }//switch( type )
    }//if( nfiles )
  
    /// \TODO: keep from issuing this message for ever single change!
    if( !msg.empty() )
      m_interspec->logMessage( msg, "", 1 );
    
    //Calling setDeviationPairMsg(...) is useless because we currently completely replace 
    //display->setDeviationPairMsg( msg );
  }else
  {
    //display->setDeviationPairMsg( "" );
  }//if( more than one gamma detector and more than one sample number ) / else
  
  m_interspec->refreshDisplayedCharts();
  refreshGuiFromFiles();
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

  m_lastGraphicalRecal = 0;
  m_lastGraphicalRecalType = EnergyCalGraphicalConfirm::NumRecalTypes;
  m_lastGraphicalRecalEnergy = -999.0f;
  
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


void EnergyCalTool::setWasGraphicalRecal( int type, float energy )
{
  time( &m_lastGraphicalRecal );
  m_lastGraphicalRecalType = type;
  m_lastGraphicalRecalEnergy = energy;
}//void setWasGraphicalRecal( int type, double energy )


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
  
  updateFitButtonStatus();
}//void specTypeToDisplayForChanged();


bool EnergyCalTool::canDoEnergyFit()
{
  // Check that "Apply Changes To" for the "Foreground" is checked, otherwise fitting makes no sense
  if( !m_applyToCbs[ApplyToCbIndex::ApplyToForeground]->isChecked() )
    return false;
  
  // Check if there are any peaks currently showing.
  shared_ptr<const deque<PeakModel::PeakShrdPtr>> peaks = m_peakModel->peaks();
  if( !peaks )
    return false;
  
  size_t nPeaksToUse = 0;
  for( const PeakModel::PeakShrdPtr &p : *peaks )
    nPeaksToUse += (p && p->useForCalibration());
  
  if( nPeaksToUse < 1 )
    return false;
  
  if( !m_calInfoDisplayStack )
    return false;
  
  // We are actually going to fit the coefficients for the currently showing CalDisplay, so only
  //  consult the checkboxes on that one display.
  auto caldisp = dynamic_cast<EnergyCalImp::CalDisplay *>( m_calInfoDisplayStack->currentWidget() );
  if( !caldisp )
    return false;
  
  auto cal = caldisp->lastSetCalibration();
  switch( cal->type() )
  {
    case SpecUtils::EnergyCalType::LowerChannelEdge:
    case SpecUtils::EnergyCalType::InvalidEquationType:
      return false;
      
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::FullRangeFraction:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      break;
  }//switch( cal->type() )
  
  const set<size_t> ordersToFit = caldisp->fitForCoefficents();
  if( ordersToFit.empty() || ordersToFit.size() > nPeaksToUse )
    return false;
  
  return true;
}//bool canDoEnergyFit()


void EnergyCalTool::fitCoefficients()
{
  try
  {
    if( !canDoEnergyFit() )
    {
      m_interspec->logMessage( "Can not fit calibration coefficients right now.  You must select at"
                              " least as many peaks for fitting as coefficients you have selected"
                              " to fit for.", "", 2 );
      return;
    }//if( double check we can actually do the fit )
    
    
    // Check if there are any peaks currently showing.
    shared_ptr<const deque<PeakModel::PeakShrdPtr>> peaks = m_peakModel->peaks();
    if( !peaks || peaks->empty() )  //shouldnt ever happen.
      throw runtime_error( "No peaks available." );
    
    // We are actually going to fit the coefficients for the currently showing CalDisplay, so only
    //  consult the checkboxes on that one display.
    auto caldisp = dynamic_cast<EnergyCalImp::CalDisplay *>( m_calInfoDisplayStack->currentWidget() );
    if( !caldisp )  //shouldnt ever happen, but JIC
      throw runtime_error( "Unexpected error determining current calibration" );
    
    // The spectrum may not be displaying the detector we are currently seeing the calibration for
    //  lets make sure of this, and if so, switch to one we are displaying
    string detname = caldisp->detectorName();
    int previousSpecTypeInd = m_specTypeMenu ? m_specTypeMenu->currentIndex() : 0;
    const auto type = (previousSpecTypeInd == 1)
                          ? SpecUtils::SpectrumType::Background
                          : (previousSpecTypeInd==2
                            ? SpecUtils::SpectrumType::SecondForeground
                            : SpecUtils::SpectrumType::Foreground);
    const vector<string> displayed = m_interspec->detectorsToDisplay( type );
    
    if( std::find(begin(displayed), end(displayed), detname) == end(displayed) )
    {
      if( previousSpecTypeInd >= 0
         && previousSpecTypeInd < 3
         && m_detectorMenu[previousSpecTypeInd] )
      {
        for( WMenuItem *item : m_detectorMenu[previousSpecTypeInd]->items() )
        {
          const string thisdetname = item->text().toUTF8();
          if( std::find(begin(displayed), end(displayed), thisdetname) != end(displayed) )
          {
            item->select();
            cout << "Starting from " << caldisp << endl;
            caldisp = dynamic_cast<EnergyCalImp::CalDisplay *>( m_calInfoDisplayStack->currentWidget() );
            cout << "  we moved to " << caldisp << endl;
            if( !caldisp )  //shouldnt ever happen, but JIC
              throw runtime_error( "Unexpected error determining current calibration" );
            
            detname = caldisp->detectorName();
          }//if( we found a displayed detector )
        }//for( loop over menu items to find a displayed detector )
      }else
      {
        throw runtime_error( "Please select calibration of a displayed detector" );
      }
    }//if( std::find(begin(displayed), end(displayed), detname) == end(displayed) )
    
        
    auto original_cal = caldisp->lastSetCalibration();
    
    
    
    switch( original_cal->type() )
    {
      case SpecUtils::EnergyCalType::LowerChannelEdge:
      case SpecUtils::EnergyCalType::InvalidEquationType:
        throw runtime_error( "Unexpected calibration type from display" ); //shouldnt ever happen, but JIC
        
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::FullRangeFraction:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        break;
    }//switch( cal->type() )
    
    if( original_cal->num_channels() < 5 )
      throw runtime_error( "Not enough channels in the data." );
    
    
    const set<size_t> orders_to_fit = caldisp->fitForCoefficents();
    if( orders_to_fit.empty() )  //shouldnt ever happen
      throw runtime_error( "You must select at least one coefficient to fit for." );
    
    
    //TODO: meansFitError will currently contain only values of 1.0, eventually
    //      will contain the error of the fit mean for that peak
    vector<EnergyCal::RecalPeakInfo> peakInfos;
    
    for( const auto &peakptr : *peaks )
    {
      if( !peakptr )  //shouldnt be necassary, but JIC
        continue;
      const PeakDef &peak = *peakptr;
      
      if( !peak.useForCalibration() )
        continue;
      
      const double wantedEnergy = peak.gammaParticleEnergy();
      
      EnergyCal::RecalPeakInfo peakInfo;
      peakInfo.peakMean = peak.mean();
      peakInfo.peakMeanUncert = max( peak.meanUncert(), 0.25 );
      if( IsInf(peakInfo.peakMeanUncert) || IsNan(peakInfo.peakMeanUncert) )
        peakInfo.peakMeanUncert = 0.5;
      
      peakInfo.photopeakEnergy = wantedEnergy;
      peakInfo.peakMeanBinNumber = original_cal->channel_for_energy( peak.mean() );
      
      peakInfos.push_back( peakInfo );
    }//for( int col = 0; col < numModelCol; ++col )
    
    if( orders_to_fit.size() > peakInfos.size() )
      throw runtime_error( "You must use at least as many peaks associated with "
                          "nuclides, as the number of coefficients you want to "
                          "fit for." );
    
    auto answer = make_shared<SpecUtils::EnergyCalibration>();
    
    const size_t eqn_order = std::max( original_cal->coefficients().size(), (*orders_to_fit.rbegin()) + 1 );
    const size_t nchannel = original_cal->num_channels();
    const auto &devpairs = original_cal->deviation_pairs();
    
    double chi2 = -999;
    
    try
    {
      vector<bool> fitfor( eqn_order, false );
      
      for( auto order : orders_to_fit )
      {
        assert( order < fitfor.size() );
        fitfor[order] = true;
      }
      
      vector<float> coefficent_uncerts;
      vector<float> coefficents = original_cal->coefficients();
      if( coefficents.size() < eqn_order )
        coefficents.resize( eqn_order, 0.0f );
      
      switch ( original_cal->type() )
      {
        case SpecUtils::EnergyCalType::Polynomial:
        case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
          chi2 = EnergyCal::fit_energy_cal_poly( peakInfos, fitfor, nchannel, devpairs,
                                                coefficents, coefficent_uncerts );
          answer->set_polynomial( nchannel, coefficents, devpairs );
          break;
          
        case SpecUtils::EnergyCalType::FullRangeFraction:
          chi2 = EnergyCal::fit_energy_cal_frf( peakInfos, fitfor, nchannel, devpairs,
                                               coefficents, coefficent_uncerts );
          answer->set_full_range_fraction( nchannel, coefficents, devpairs );
          break;
          
        case SpecUtils::EnergyCalType::LowerChannelEdge:
        case SpecUtils::EnergyCalType::InvalidEquationType:
          throw runtime_error( "Didnt expect lower channel or invalid eqn type." );
          break;
      }//switch ( original_cal->type() )
      
      //Print some developer info to terminal
      stringstream msg;
      msg << "\nfit_energy_cal_poly gave chi2=" << chi2 << " with coefs={";
      for( size_t i = 0; i < coefficents.size(); ++i )
        msg << coefficents[i] << "+-" << coefficent_uncerts[i] << ", ";
      msg << "}\n";
      cout << msg.str() << endl;
      
    }catch( std::exception &e )
    {
      cerr << "fit_energy_cal_poly threw: " << e.what() << endl;
#if( PERFORM_DEVELOPER_CHECKS )
      char buffer[512] = { '\0' };
      snprintf( buffer, sizeof(buffer)-1, "fit_energy_cal_poly threw: %s", e.what() );
      log_developer_error( __func__, buffer );
#endif
    }//try / catch fit for coefficents using least linear squares
    
    
    if( !answer->valid() )
    {
      vector<bool> fitfor( eqn_order, false );
      
      for( auto order : orders_to_fit )
        fitfor[order] = true;
      
      vector<float> calib_coefs = original_cal->coefficients();
      if( calib_coefs.size() < eqn_order )
        calib_coefs.resize( eqn_order, 0.0f );
      
      std::string warning_msg;
      std::vector<float> coefs, coefs_uncert;
      chi2 = EnergyCal::fit_energy_cal_iterative( peakInfos, nchannel, original_cal->type(), fitfor,
                                          calib_coefs, devpairs, coefs, coefs_uncert, warning_msg );
      
      if( warning_msg.size() )
        m_interspec->logMessage( warning_msg, "", 3 );
      
      switch ( original_cal->type() )
      {
        case SpecUtils::EnergyCalType::Polynomial:
        case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
          answer->set_polynomial( nchannel, coefs, devpairs );
          break;
          
        case SpecUtils::EnergyCalType::FullRangeFraction:
          answer->set_full_range_fraction( nchannel, coefs, devpairs );
          break;
          
        case SpecUtils::EnergyCalType::LowerChannelEdge:
        case SpecUtils::EnergyCalType::InvalidEquationType:
          assert( 0 );
          break;
      }//switch ( original_cal->type() )
    }//if( !fit_coefs )
    
    if( !answer->valid() )
      throw runtime_error( "ErrorMsg Failed to fit energy calibration using both matrix and"
                           " iterative approaches" );
    
  
    if( !answer || !answer->valid() )
      return;
    
    applyCalChange( original_cal, answer, false );
    
    
    string msg = "Energy calibration has been updated from fitting.";
    
    //To show Chi2 in the message, uncomment out this next section
    /*
    if( peakInfos.size() > orders_to_fit.size() )
    {
      double dof = peakInfos.size() - 1;
      dof -= orders_to_fit.size();
      dof = (dof < 1) ? 1.0 : dof;
      
      char buffer[64];
      snprintf( buffer, sizeof(buffer), " &chi;&sup2;/dof=%.2g", (chi2/dof) );
      msg += buffer;
    }
    */
    
    m_interspec->logMessage( WString::fromUTF8(msg), "", 1 );
  }catch( std::exception &e )
  {
    string msg = "Failed calibration by fitting peak means. ";
    msg += e.what();
    
    cerr << "EnergyCalTool::fitCoefficients():\n\tCaught: " << msg << endl;
    m_interspec->logMessage( msg, "", 3 );
  }//try / catch
}//void fitCoefficients()


void EnergyCalTool::updateFitButtonStatus()
{
  const bool canFit = canDoEnergyFit();
  
  if( canFit != m_fitCalBtn->isEnabled() )
    m_fitCalBtn->setDisabled( !canFit );
}//void updateFitButtonStatus()


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


void EnergyCalTool::handleGraphicalRecalRequest( double xstart, double xfinish )
{
  try
  {
    auto foreground = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    shared_ptr<const SpecUtils::EnergyCalibration> energycal = foreground
                                                               ? foreground->energy_calibration()
                                                               : nullptr;
    if( !energycal || !energycal->valid()
        || (energycal->type() == SpecUtils::EnergyCalType::LowerChannelEdge) )
      return;
  
    if( m_graphicalRecal )
    {
      m_graphicalRecal->setEnergies( xstart, xfinish );
    }else
    {
      m_graphicalRecal = new EnergyCalGraphicalConfirm( xstart, xfinish, this,
                                    m_lastGraphicalRecal,
                                    EnergyCalGraphicalConfirm::RecalTypes(m_lastGraphicalRecalType),
                                    m_lastGraphicalRecalEnergy );
    
      m_graphicalRecal->finished().connect( this, &EnergyCalTool::deleteGraphicalRecalConfirmWindow );
    }
  }catch( std::runtime_error & )
  {
    m_interspec->logMessage( "Internal error doing graphical recal; sorry :(", "", 3 );
  }
}//void handleGraphicalRecalRequest( double xstart, double xfinish )


void EnergyCalTool::deleteGraphicalRecalConfirmWindow()
{
  if( m_graphicalRecal )
  {
    AuxWindow::deleteAuxWindow( m_graphicalRecal );
    m_graphicalRecal = nullptr;
  }//if( m_graphicalRecal )
  
  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_interspec );
  if( showToolTips )
  {
#if( USE_SPECTRUM_CHART_D3 )
    m_interspec->logMessage( "If you recalibrate again by ALT+CTRL+DRAG on"
                            " another portion of the spectrum, you will be given the"
                            " option of preserving the effects of this calibration"
                            " as well.", "", 1 );
#else
    m_interspec->logMessage( "If you recalibrate again by right-click dragging on"
                            " another portion of the spectrum, you will be given the"
                            " option of preserving the effects of this calibration"
                            " as well.", "", 1 );
#endif
  }//if( showToolTips )
}//void deleteGraphicalRecalConfirmWindow()


string EnergyCalTool::applyToSummaryTxt() const
{
  string answer;
  
  if( m_applyToColumn->isHidden() )
    return answer;
  
  for( ApplyToCbIndex index = static_cast<ApplyToCbIndex>(0);
      index < ApplyToCbIndex::NumApplyToCbIndex;
      index = ApplyToCbIndex(index+1) )
  {
    Wt::WCheckBox *cb = m_applyToCbs[index];
    const auto cbparent = cb->parent();
    assert( cbparent );
    assert( dynamic_cast<WContainerWidget *>(cbparent) );
    
    if( cbparent->isHidden() || !cb->isChecked() )
      continue;
    
    if( !answer.empty() )
      answer += ", ";
    
    switch( index )
    {
      case ApplyToCbIndex::ApplyToForeground:         answer += "foreground"; break;
      case ApplyToCbIndex::ApplyToBackground:         answer += "background"; break;
      case ApplyToCbIndex::ApplyToSecondary:          answer += "secondary";  break;
      case ApplyToCbIndex::ApplyToDisplayedDetectors: answer += "displayed detectors"; break;
      case ApplyToCbIndex::ApplyToAllDetectors:       answer += "all detectors"; break;
      case ApplyToCbIndex::ApplyToDisplayedSamples:   answer += "displayed samples"; break;
      case ApplyToCbIndex::ApplyToAllSamples:         answer += "all samples"; break;
      case ApplyToCbIndex::NumApplyToCbIndex:
        break;
    }//switch( index )
  }//for( loop over ApplyToCbIndex )
  
  return answer;
}//string applyToSummaryTxt() const


void EnergyCalTool::doRefreshFromFiles()
{
  //Labels for horizantal labels when you have mutliple spectra shown, and at least one of them has
  //  more than one detectors
  const char * const spec_type_labels[3] = {"For.","Back","Sec."};
  const char * const spec_type_labels_vert[3] = {"Foreground", "Background", "Secondary"};
  
  
  
  
  
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
  
  const bool isWide = (m_tallLayoutContent ? false : true);
  
  //If each spectrum file has only a single detector, then instead of having vertical menu display
  //  detector name, we will have "For.", "Back", "Sec."
  bool specTypeInForgrndMenu = true;
  
  //If we try to re-use the menu and stacks, for some reason the calibration coefficents wont show
  //  up if we alter which background/secondary spectra are showing... not sure why, but for the
  //  moment we'll just re-create the menus and stacks... not great, but works, for the moment.
  bool needStackRefresh = false;
  
  // Get the detector names, for the displayed sample numbers, for each spectrum type
  set<string> disp_det_names[3];
  
  int nFilesWithCalInfo = 0;
  
  // We want to preserve wich "Fit" check boxes are set.
  map< pair<SpecUtils::SpectrumType,string>, set<size_t> > set_fit_for_cbs;
  
  
  {//begin code-block to see if we need to refresh stack
    
    //TODO: need to check that isWide hasnt changed, and if it has set needStackRefresh=true
    
    //Check if we can put "For.", "Back", "Sec." on the vertical menu instead of detector names
    for( int i = 0; i < 3; ++i )
    {
      disp_det_names[i] = gammaDetectorsForDisplayedSamples( spectypes[i] );
      specTypeInForgrndMenu = (specTypeInForgrndMenu && (disp_det_names[i].size() <= 1));
      
      if( !disp_det_names[i].empty() )
        nFilesWithCalInfo += 1;
    }//for( int i = 0; i < 3; ++i )
    
    if( !m_specTypeMenu )
    {
      //cout << "needStackRefresh: m_specTypeMenu == nullptr" << endl;
      needStackRefresh = true;
    }
    
    for( int i = 0; !needStackRefresh && (i < 3); ++i )
    {
      WMenuItem *typeItem = m_specTypeMenu->itemAt(i);
      if( !typeItem )
      {
        //cout << "needStackRefresh: !typeItem" << endl;
        needStackRefresh = true;
        continue;
      }
      
      WMenu *detMenu = m_detectorMenu[(specTypeInForgrndMenu ? 0 : i)];
      
      if( !detMenu )
      {
        //cout << "needStackRefresh: !detMenu" << endl;
        needStackRefresh = true;
        continue;
      }
      
      if( !specTypeInForgrndMenu && ((!detMenu) != (!specfiles[i])) )
      {
        //cout << "needStackRefresh: (!detMenu != !specfiles[i])" << endl;
        needStackRefresh = true;
        continue;
      }
      
      if( specTypeInForgrndMenu )
      {
        if( specfiles[i] )
        {
          //Make sure one of the widgets have spec_type_labels_vert[i] in it
          needStackRefresh = true;
          for( int w = 0; needStackRefresh && (w < detMenu->count()); ++w )
          {
            auto item = detMenu->itemAt(w);
            needStackRefresh = !(item && (item->text() == WString(spec_type_labels_vert[i])));
          }
          
          //if( needStackRefresh )
          //  cout << "needStackRefresh: Did not have a menu entry for specfile[" << i << "]" << endl;
        }else
        {
          //Make sure none of the widgets have spec_type_labels_vert[i] in them
          for( int w = 0; !needStackRefresh && (w < detMenu->count()); ++w )
          {
            auto item = detMenu->itemAt(w);
            needStackRefresh = (!item || (item->text() == WString(spec_type_labels_vert[i])));
          }
          
          //if( needStackRefresh )
          //  cout << "needStackRefresh: For empty place " << i << " had a menu entry, but no spectrum file" << endl;
        }//if( specfiles[i] ) / else
      }else
      {
        assert( m_specTypeMenu );
        WMenuItem *item = m_specTypeMenu->itemAt(i);
        assert( item );
     
        if( !item )
        {
          //Shouldnt ever get here I think
          needStackRefresh = true;
          continue;
        }//if( !item )
        
        if( !specfiles[i] )
        {
          //If there is no spectrum file for index i, item should be hidden, and if not need a refresh
          needStackRefresh = !item->isHidden();
          //if( needStackRefresh )
          //  cout << "needStackRefresh: item is not hidden for missing spectrum " << i << endl;
          continue;
        }//if( !specfiles[i] )
        
        if( item->isHidden() )
        {
          //If item is hidden, but here we do have a spectrum file, we need a refresh
          //cout << "needStackRefresh: item is hidden for " << i << endl;
          needStackRefresh = true;
          continue;
        }//if( item->isHidden() )
        
        
        //Need to check all the detectors are the same names
        set<string> detsInMenu;
        for( int j = 0; j < detMenu->count(); ++j )
        {
          auto detitem = detMenu->itemAt(j);
          if( detitem )
            detsInMenu.insert( detitem->text().toUTF8() );
        }//for( int j = 0; j < detMenu->count(); ++j )
        
        needStackRefresh = (detsInMenu != disp_det_names[i]);
        //if( needStackRefresh )
        //  cout << "needStackRefresh: New det names dont equal old for type=" << i << endl;
      }//if( specTypeInForgrndMenu ) / else
    }//for( int i = 0; i < 3; ++i )

    
  }//end code-block to see if we need to refresh stack
  
  //Dont show spectype menu (the vertical "For.", "Back", "Sec." menu), if we dont need to
  const bool hideSpecType = ( specTypeInForgrndMenu
                             || (nFilesWithCalInfo < 2)
                             || ( (!specfiles[1] || (specfiles[0]==specfiles[1]))
                                 && (!specfiles[2] || (specfiles[0]==specfiles[2]))) );
  
  if( !m_specTypeMenu || (m_specTypeMenu->isHidden() != hideSpecType) )
  {
    //cout << "needStackRefresh: m_specTypeMenu->isHidden() != hideSpecType" << endl;
    needStackRefresh = true;
  }
  
  //cout << "needStackRefresh=" << needStackRefresh << endl;
  
  // TODO: instead of the above logic to catch when we need to refresh (\e.g., create all new)
  //       widgets, should combine it with the logic to create new widgets, but only delte or create
  
  if( needStackRefresh )
  {
    for( int i = 0; i < 3; ++i )
    {
      if( m_detectorMenu[i] )
      {
        for( WMenuItem *item : m_detectorMenu[i]->items() )
        {
          /// \TODO: the item text isnt necassarily the detector name - when specTypeInForgrndMenu
          ///        is true, this breaks down - need to fix this
          const string detname = item->text().toUTF8();
          auto display = dynamic_cast<EnergyCalImp::CalDisplay *>( item->contents() );
          if( display )
            set_fit_for_cbs[{spectypes[i],detname}] = display->fitForCoefficents();
          else
            cerr << "Unexpected widget type as a sub-menu!" << endl;
        }//for( loop over menu types )
        
        if( m_detectorMenu[i]->currentItem() )
          prevdet[i] = m_detectorMenu[i]->currentItem()->text().toUTF8();
        delete m_detectorMenu[i];
      }
      m_detectorMenu[i] = nullptr;
    }//for( int i = 0; i < 3; ++i )
    
    delete m_specTypeMenuStack;
    delete m_specTypeMenu;
    m_specTypeMenuStack = nullptr;
    m_specTypeMenu = nullptr;
    
    WAnimation animation(Wt::WAnimation::Fade, Wt::WAnimation::Linear, 200);
    
    m_specTypeMenuStack = new WStackedWidget();
    m_specTypeMenuStack->addStyleClass( "CalSpecStack" );
    m_specTypeMenuStack->setTransitionAnimation( animation );
    
    auto callayout = dynamic_cast<WGridLayout *>( m_calColumn->layout() );
    assert( callayout );
    
    m_specTypeMenu = new WMenu( m_specTypeMenuStack );
    m_specTypeMenu->addStyleClass( "CalSpecMenu" );
    m_specTypeMenu->itemSelected().connect( this, &EnergyCalTool::specTypeToDisplayForChanged );
    m_detColLayout->addWidget( m_specTypeMenu, 1, 0 );
    m_detColLayout->addWidget( m_specTypeMenuStack, 2, 0 );
    
    if( m_calInfoDisplayStack )
      delete m_calInfoDisplayStack;
    m_calInfoDisplayStack = new WStackedWidget();
    m_calInfoDisplayStack->addStyleClass( "CalColContent CalStack" );
    m_calInfoDisplayStack->setTransitionAnimation( animation );
    callayout->addWidget( m_calInfoDisplayStack, 1, 1 );
    
    /// \TODO: only create these menus when actually needed, so we wont need to
    for( int i = 0; i < 3; ++i )
    {
      if( !specfiles[i] )
      {
        //Add a dummy entry into the menu or else the 'm_specTypeMenu->itemAt(i)' call below will
        //  segfault or not necassarily give the wanted answer.
        WContainerWidget *detMenuDiv = new WContainerWidget();
        WMenuItem *item = m_specTypeMenu->addItem( "", detMenuDiv );
        item->setHidden( true );
        continue;
      }
      
      WContainerWidget *detMenuDiv = new WContainerWidget();  //this holds the WMenu for this SpecFile
      detMenuDiv->addStyleClass( "DetMenuDiv" );
      
      WMenuItem *item = m_specTypeMenu->addItem( spec_type_labels[i], detMenuDiv, WMenuItem::LoadPolicy::PreLoading );
      //Fix issue, for Wt 3.3.4 at least, if user doesnt click exactly on the <a> element
      item->clicked().connect( boost::bind(&WMenuItem::select, item) );
      
      m_detectorMenu[i] = new WMenu( m_calInfoDisplayStack, detMenuDiv );
      m_detectorMenu[i]->addStyleClass( "VerticalNavMenu HeavyNavMenu DetCalMenu" );
      
      m_detectorMenu[i]->itemSelected().connect( this, &EnergyCalTool::updateFitButtonStatus );
    }//for( int i = 0; i < 3; ++i )
  }//if( needStackRefresh )
  

  
  if( !specfiles[0] )
  {
    setShowNoCalInfo( true );
    m_fitCalBtn->disable();
    return;
  }
  
  if( previousSpecInd < 0 ||  previousSpecInd > 2
     || (previousSpecInd == 1 && !specfiles[1])
     || (previousSpecInd == 2 && !specfiles[2])  )
  {
    previousSpecInd = 0;
  }
  
  
  bool selectedDetToShowCalFor = false;
  bool hasFRFCal = false, hasPolyCal = false, hasLowerChanCal = false;
  
  
  for( int i = 0; i < 3; ++i )
  {
    if( !specfiles[i] )
      continue;
    
    Wt::WMenu *detMenu = m_detectorMenu[(specTypeInForgrndMenu ? 0 : i)];
    if( !detMenu )
      continue;
    
    WMenuItem *specItem = m_specTypeMenu->itemAt(i);
    assert( specItem );
    
    const SpecUtils::SpectrumType type = spectypes[i];
    shared_ptr<const SpecMeas> meas = specfiles[i];
    assert( meas );
    
    const set<int> &samples = m_interspec->displayedSamples(type);
    const set<string> &detectors = disp_det_names[i];
    
    specdetnames[i] = detectors;
    
    if( detectors.empty() )
    {
      if( m_specTypeMenu->currentIndex() == i )
        m_specTypeMenu->select(0);
      specItem->setHidden( true );
      continue;
    }
    

    assert( specItem );
    specItem->setHidden( false );
    
    for( const string &detname : detectors )
    {
      for( const int sample : samples )
      {
        auto m = meas->measurement( sample, detname );
        if( !m || (m->num_gamma_channels() <= 4) )
          continue;
        
        shared_ptr<const SpecUtils::EnergyCalibration> energycal = m->energy_calibration();
        
        string displayname = detname;
        
        if( specTypeInForgrndMenu )
        {
          /// \TODO: when specTypeInForgrndMenu is true, we may not match detector name because
          ///        we are getting the menu item text above and assuming the detector name.
          ///        should fix this.
          displayname = spec_type_labels_vert[i];
        }//if( specTypeInForgrndMenu )
        
        WMenuItem *item = nullptr;
        if( !needStackRefresh )
        {
          for( int j = 0; !item && (j < detMenu->count()); ++j )
          {
            auto jitem = detMenu->itemAt(j);
            if( jitem && (jitem->text() == WString(displayname)) )
              item = jitem;
          }
        }//if( !needStackRefresh )
        
        if( item )
        {
          WWidget *ww = item->contents();
          EnergyCalImp::CalDisplay *calcontent = dynamic_cast<EnergyCalImp::CalDisplay *>(ww);
          assert( calcontent );
          if( calcontent )
            calcontent->updateToGui( energycal );
          else
            item = nullptr;
        }//if( item )
        
        if( !item )
        {
          auto calcontent = new EnergyCalImp::CalDisplay( this, type, detname, isWide );
          item = detMenu->addItem( WString::fromUTF8(displayname), calcontent, WMenuItem::LoadPolicy::PreLoading );
          //Fix issue, for Wt 3.3.4 at least, if user doesnt click exactly on the <a> element
          item->clicked().connect( boost::bind(&WMenuItem::select, item) );
          calcontent->updateToGui( energycal );
          
          const auto fitfor_iter = set_fit_for_cbs.find( {type,displayname} );
          if( fitfor_iter != end(set_fit_for_cbs) )
            calcontent->setFitFor( fitfor_iter->second );
          
          if( (displayname == prevdet[i]) && (i == previousSpecInd) )
          {
            m_specTypeMenu->select( previousSpecInd );
            detMenu->select( item );
            selectedDetToShowCalFor = true;
          }
        }//if( !item )
        
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
      
        break;
      }
    }//for( const string &detname : detectors )
  }//for( int i = 0; i < 3; ++i )
  
  if( needStackRefresh && !selectedDetToShowCalFor && m_detectorMenu[0] && m_detectorMenu[0]->count() )
  {
    m_specTypeMenu->select( 0 );
    m_detectorMenu[0]->select( 0 );
  }
  
  setShowNoCalInfo( !nFilesWithCalInfo );
  m_specTypeMenu->setHidden( hideSpecType );
  
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
      index = MoreActionsIndex(static_cast<int>(index) + 1) )
  {
    Wt::WAnchor *anchor = m_moreActions[static_cast<int>(index)];
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
        aparent->setHidden( !hasPolyCal );
        break;
        
      case MoreActionsIndex::ConvertToPoly:
        aparent->setHidden( !(hasFRFCal || hasLowerChanCal) );
        break;
        
      case MoreActionsIndex::MultipleFilesCal:
      {
        SpecMeasManager *manager = m_interspec->fileManager();
        SpectraFileModel *fmodel = manager ? manager->model() : nullptr;
        const int nfiles = fmodel ? fmodel->rowCount() : 0;
        int nRecordsWithPeaks = 0;
        for( int row = 0; row < nfiles && (nRecordsWithPeaks < 2); ++row )
        {
          shared_ptr<SpectraFileHeader> header = fmodel->fileHeader( row );
          if( !header )
            continue;
          
          int nsamples = header->numSamples();
          shared_ptr<SpecMeas> meas = header->measurementIfInMemory();
          if( meas )
            nRecordsWithPeaks += meas->sampleNumsWithPeaks().size();
          else
            nRecordsWithPeaks += nsamples; //not worth readin file from disk, so we'll be hopeful
        }
        
        aparent->setHidden( nRecordsWithPeaks < 2 );
        break;
      }//case MoreActionsIndex::MultipleFilesCal:
        
      case MoreActionsIndex::NumMoreActionsIndex:
        break;
    }//switch( index )
  }//for( loop over
  
  // Update the "Fit Coeffs" button to be enabled/disabled
  updateFitButtonStatus();
  
  //const int currentwidget = m_detectorMenu[0]->contentsStack()->currentIndex();
  //cout << "currentwidget=" << currentwidget << endl;
}//void doRefreshFromFiles()



void EnergyCalTool::moreActionBtnClicked( const MoreActionsIndex index )
{
  const vector<MeasToApplyCoefChangeTo> measToChange = measurementsToApplyCoeffChangeTo();

  new EnergyCalAddActionsWindow( index, measToChange, this );
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
      updateFitButtonStatus();
      break;
      
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
