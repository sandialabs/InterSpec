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
#include <cmath>
#include <tuple>
#include <vector>
#include <utility>
#include <algorithm>

#include <Wt/WLabel>
#include <Wt/WCheckBox>
#include <Wt/WGridLayout>
#include <Wt/WContainerWidget>

#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/EnergyCalDevPairWidget.h"

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
}// namespace


namespace EnergyCalImp
{

DevPair::DevPair( const bool show_fit_offset, Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_energy( new NativeFloatSpinBox() ),
    m_offset( new NativeFloatSpinBox() ),
    m_fitOffset( nullptr ),
    m_delete( new WContainerWidget() )
{
  WGridLayout* layout = new WGridLayout();
  layout->setContentsMargins(0, 0, 0, 0);
  layout->setVerticalSpacing( 0 );

  setLayout(layout);
  setStyleClass( "DevPair" );

  int col = 0;
  layout->addWidget( m_energy, 0, col ); layout->setColumnStretch( col, 1 ); ++col;
  layout->addWidget( m_offset, 0, col ); layout->setColumnStretch( col, 1 ); ++col;

  if( show_fit_offset )
  {
    m_fitOffset = new WCheckBox( "Fit offset" );
    m_fitOffset->addStyleClass( "DevFitOffset CbNoLineBreak" );
    layout->addWidget( m_fitOffset, 0, col );
    layout->setColumnStretch( col, 0 );
    ++col;
  }//if( show_fit_offset )

  layout->addWidget( m_delete, 0, col );
  layout->setColumnStretch( col, 0 );

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
    else if( fabs(fraction - 0.1f) < 1.0E-4f )
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


bool DevPair::fitOffset() const
{
  return (m_fitOffset && m_fitOffset->isChecked());
}


void DevPair::setFitOffset( const bool fit )
{
  if( m_fitOffset )
    m_fitOffset->setChecked( fit );
}


void DevPair::visuallyIndicateChanged()
{
  doJavaScript( "$('#" + id() + "').fadeIn(100).fadeOut(100).fadeIn(100).fadeOut(100).fadeIn(100);" );
}//void visuallyIndicateChanged();


DeviationPairDisplay::DeviationPairDisplay( const bool show_fit_offsets, Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_pairs( NULL ),
    m_show_fit_offsets( show_fit_offsets )
{
  addStyleClass( "DevPairDisplay" );
  WLabel *title = new WLabel( WString::tr("ect-deviation-pairs"), this );
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
  std::sort( d.begin(), d.end() );

  // Update the existing DevPair rows in-place, and only create/delete rows as needed, so we dont
  //  disturb the widget the user may be interacting with (e.g., cursor position would be lost if
  //  we re-set the text of the field they are typing in).  Rows with either field blank (e.g.,
  //  the user clicked add, but hasnt filled things out yet) dont correspond to an actual
  //  deviation pair, so they are left untouched.
  //
  // The positional in-place update re-purposes a row to a different (energy,offset) pair, which
  //  would otherwise leave that rows "Fit offset" checkbox attached to whatever pair lands in that
  //  slot - wrong when the pair set changes (e.g., a fit prepends a {0,0} spline anchor).  So we
  //  snapshot the checkbox state keyed by energy, then re-apply it by energy afterwards.  Deviation
  //  pair energies pass through a fit unchanged, so an exact match is correct.
  std::map<float,bool> fit_offset_by_energy;

  vector<DevPair *> complete_rows;
  const vector<WWidget *> childs = m_pairs->children();
  for( WWidget *t : childs )
  {
    DevPair * const p = dynamic_cast<DevPair *>( t );
    if( p && !p->m_energy->text().empty() && !p->m_offset->text().empty() )
    {
      if( m_show_fit_offsets )
        fit_offset_by_energy[p->devPair().first] = p->fitOffset();
      complete_rows.push_back( p );
    }
  }//for( WWidget *t : childs )

  // Rows are kept sorted by energy (see #sortDisplayOrder), same as `d`, so we can just compare
  //  and assign positionally.
  for( size_t i = 0; i < std::min(complete_rows.size(), d.size()); ++i )
  {
    DevPair * const p = complete_rows[i];
    const pair<float,float> current = p->devPair();
    if( (current.first != d[i].first) || (current.second != d[i].second) )
      p->setDevPair( d[i] );
  }//for( rows to update in-place )

  for( size_t i = d.size(); i < complete_rows.size(); ++i )
    delete complete_rows[i];

  for( size_t i = complete_rows.size(); i < d.size(); ++i )
  {
    DevPair * const dev = newDevPair( false );
    dev->setDevPair( d[i] );
  }//for( rows to add )

  sortDisplayOrder(false);

  // Re-attach each rows "Fit offset" checkbox to its deviation pair (matched by energy); pairs not
  //  present before (e.g., a fit-inserted {0,0} anchor) default to un-checked.
  if( m_show_fit_offsets )
  {
    const vector<WWidget *> rows = m_pairs->children();
    for( WWidget *t : rows )
    {
      DevPair * const p = dynamic_cast<DevPair *>( t );
      if( !p || p->m_energy->text().empty() || p->m_offset->text().empty() )
        continue;
      const std::map<float,bool>::const_iterator pos = fit_offset_by_energy.find( p->devPair().first );
      p->setFitOffset( (pos != end(fit_offset_by_energy)) && pos->second );
    }//for( WWidget *t : rows )
  }//if( m_show_fit_offsets )
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

  return answer;
}//deviationPairs()


vector< tuple<float,float,bool> > DeviationPairDisplay::deviationPairsAndFit() const
{
  vector< tuple<float,float,bool> > answer;
  const vector<WWidget *> childs = m_pairs->children();
  for( const WWidget *t : childs )
  {
    const DevPair *p = dynamic_cast<const DevPair *>( t );
    if( !p || p->m_energy->text().empty() || p->m_offset->text().empty() )
      continue;

    const pair<float,float> dp = p->devPair();
    answer.emplace_back( dp.first, dp.second, p->fitOffset() );
  }//for( WWidget *t : childs )

  std::sort( begin(answer), end(answer),
             []( const tuple<float,float,bool> &a, const tuple<float,float,bool> &b ){
               return std::get<0>(a) < std::get<0>(b);
             } );

  return answer;
}//deviationPairsAndFit()


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


void DeviationPairDisplay::setPairsAreaMaxHeight( const int max_height_px )
{
  m_pairs->setMaximumSize( WLength::Auto, WLength(max_height_px, WLength::Pixel) );
  m_pairs->setOverflow( WContainerWidget::OverflowAuto, Wt::Vertical );
}//setPairsAreaMaxHeight(...)


DevPair *DeviationPairDisplay::newDevPair( const bool emitChangedNow )
{
  DevPair *dev = new DevPair( m_show_fit_offsets, m_pairs );
  dev->m_delete->clicked().connect( boost::bind( &DeviationPairDisplay::removeDevPair, this, dev ) );
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

}//namespace EnergyCalImp
