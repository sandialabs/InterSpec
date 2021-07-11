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
#include <memory>
#include <iostream>
#include <stdexcept>

// Disable streamsize <=> size_t warnings in boost
#pragma warning(disable:4244)

#include <boost/any.hpp>

#include <Wt/WSignal>
#include <Wt/WString>
#include <Wt/WModelIndex>
#include <Wt/Chart/WDataSeries>
#include <Wt/WAbstractItemModel>

#include "SpecUtils/SpecFile.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/SpectrumDataModel.h"

using namespace std;
using namespace Wt;

using SpecUtils::Measurement;

namespace
{
  const Wt::WColor ns_default_foreground_color( black );
  const Wt::WColor ns_default_background_color( cyan );
  const Wt::WColor ns_default_secondary_color( darkGreen );
}


SpectrumDataModel::SpectrumDataModel( Wt::WObject *parent )
  : Wt::WAbstractItemModel( parent ),
    m_rebinFactor( 1 ),
    m_dataLiveTime( -1.0 ),
    m_dataRealTime( -1.0 ),
    m_dataNeutronCounts( -1.0 ),
    m_backgroundLiveTime( -1.0 ),
    m_backgroundRealTime( -1.0 ),
    m_backgroundNeutronCounts( -1.0 ),
    m_secondDataLiveTime( -1.0 ),
    m_secondDataRealTime( -1.0 ),
    m_secondDataNeutronCounts( -1.0 ),
    m_backgroundSF( -1.0 ),
    m_secondSF( -1.0 ),
    m_secondDataOwnAxis(  true ),
    m_backgroundSubtract( false ),
    m_addHistIntegralToLegend( true ),
    m_dataSet( this ),
    m_seriesColors{ ns_default_foreground_color, ns_default_background_color, ns_default_secondary_color }
{
} // SpectrumDataModel constructor

SpectrumDataModel::~SpectrumDataModel()
{
  // no-op
} // SpectrumDataModel destructor


Wt::Signal<SpectrumDataModel::ColumnType> &SpectrumDataModel::dataSet()
{
  return m_dataSet;
}//Wt::Signal<ColumnType> &dataSet()

void SpectrumDataModel::setDataHistogram( std::shared_ptr<Measurement> hist )
{
  const float liveTime = hist ? hist->live_time() : 0.0f;
  const float realTime = hist ? hist->real_time() : 0.0f;
  const float neutronCounts = hist ? hist->neutron_counts_sum() : -1.0f;
  
  // Store the data.
  m_data = hist;
  m_dataLiveTime = liveTime;
  m_dataRealTime = realTime;
  m_dataNeutronCounts = neutronCounts;
  
  WModelIndex dataStart = index( 0, DATA_COLUMN );
  WModelIndex dataEnd   = index( rowCount() - 1, DATA_COLUMN );
  dataChanged().emit( dataStart, dataEnd );

  WModelIndex axisStart = index( 0, X_AXIS_COLUMN );
  WModelIndex axisEnd   = index( rowCount() - 1, X_AXIS_COLUMN );
  dataChanged().emit( axisStart, axisEnd );
  
  float newBackSF = m_dataLiveTime/m_backgroundLiveTime;
  if( IsInf(newBackSF) || IsNan(newBackSF) )
    newBackSF = 1.0f;
  if( !!m_background && fabs(m_backgroundSF-newBackSF)>0.0001 )
  {
    m_backgroundSF = newBackSF;
    WModelIndex backStart = index( 0, BACKGROUND_COLUMN );
    WModelIndex backEnd   = index( rowCount() - 1, BACKGROUND_COLUMN );
    dataChanged().emit( backStart, backEnd );
  }
  
  float newSecondSF = m_dataLiveTime/m_backgroundLiveTime;
  if( IsInf(newSecondSF) || IsNan(newSecondSF) )
    newSecondSF = 1.0f;
  if( !!m_secondData && fabs(m_secondSF-newSecondSF)>0.0001 )
  {
    m_secondSF = newSecondSF;
    m_secondSF = m_dataLiveTime/m_secondDataLiveTime;
    WModelIndex secondStart = index( 0, BACKGROUND_COLUMN );
    WModelIndex secondEnd   = index( rowCount() - 1, BACKGROUND_COLUMN );
    dataChanged().emit( secondStart, secondEnd );
  }
  
  m_dataSet.emit( DATA_COLUMN );
}//void setDataHistogram( std::shared_ptr<Measurement> hist )


void SpectrumDataModel::setSecondDataHistogram( std::shared_ptr<Measurement> hist, const bool ownAxis )
{
  const float liveTime = hist ? hist->live_time() : 0.0f;
  const float realTime = hist ? hist->real_time() : 0.0f;
  const float neutronCounts = hist ? hist->neutron_counts_sum() : -1.0f;
  
  m_secondDataLiveTime = liveTime;
  m_secondDataRealTime = realTime;
  m_secondDataNeutronCounts = neutronCounts;

  if( !m_secondData && !hist )
    return;

  const bool isAdding   = !m_secondData;
  const bool isDeleting = !hist;

  if( isAdding )
    beginInsertColumns( WModelIndex(), SECOND_DATA_COLUMN, SECOND_DATA_COLUMN );
  else if( isDeleting )
    beginRemoveColumns( WModelIndex(), SECOND_DATA_COLUMN, SECOND_DATA_COLUMN );

  // Store the data.
  m_secondData = hist;
  
  if( m_dataLiveTime > 0.0 )
    m_secondSF = m_dataLiveTime/m_secondDataLiveTime;
  else
    m_secondSF = 1.0f;

  if( IsInf(m_secondSF) || IsNan(m_secondSF) )
    m_secondSF = 1.0f;
  
  if( isAdding )
    endInsertColumns();
  else if( isDeleting )
    endRemoveColumns();

  if( !isAdding && m_secondData )
  {
    WModelIndex start = index( 0, 2 );
    WModelIndex end   = index( rowCount() - 1, 2 );

    dataChanged().emit( start, end );
  } // if( background data was changed )

  m_secondDataOwnAxis = ownAxis;
  
  m_dataSet.emit( SECOND_DATA_COLUMN );
}//void setSecondDataHistogram(...);


void SpectrumDataModel::setBackgroundHistogram( std::shared_ptr<Measurement> hist )
{
  const float liveTime = hist ? hist->live_time() : 0.0f;
  const float realTime = hist ? hist->real_time() : 0.0f;
  const float neutronCounts = hist ? hist->neutron_counts_sum() : -1.0f;
  
  m_backgroundLiveTime = liveTime;
  m_backgroundRealTime = realTime;
  m_backgroundNeutronCounts = neutronCounts;

  if( !m_background && !hist )
    return;

  const bool isAdding   = !m_background;
  const bool isDeleting = !hist;

  if( isAdding )
    beginInsertColumns( WModelIndex(), BACKGROUND_COLUMN, BACKGROUND_COLUMN );
  else if( isDeleting )
    beginRemoveColumns( WModelIndex(), BACKGROUND_COLUMN, BACKGROUND_COLUMN );

  // Store the data.
  m_background = hist;
  
  if( (m_dataLiveTime > DBL_EPSILON) && (m_backgroundLiveTime > DBL_EPSILON) )
    m_backgroundSF = m_dataLiveTime/m_backgroundLiveTime;
  else
    m_backgroundSF = 1.0f;

  if( IsInf(m_backgroundSF) || IsNan(m_backgroundSF) )
    m_backgroundSF = 1.0f;
  
  
  if( isAdding )
    endInsertColumns();
  else if( isDeleting )
    endRemoveColumns();

  if( !isAdding && m_background )
  {
    WModelIndex start = index( 0, BACKGROUND_COLUMN );
    WModelIndex end   = index( rowCount() - 1, BACKGROUND_COLUMN );

    dataChanged().emit( start, end );
  } // if( the background data was modified )
  
  m_dataSet.emit( BACKGROUND_COLUMN );
}//void setBackgroundHistogram( std::shared_ptr<Measurement> hist )



bool SpectrumDataModel::backgroundSubtract() const
{
  return m_backgroundSubtract;
}//bool SpectrumDataModel::backgroundSubtract() const


void SpectrumDataModel::setBackgroundSubtract( const bool subtract )
{
  // If it's already set to the value, don't do anything.
  if( m_backgroundSubtract == subtract ) return;

  // Otherwise, clear it,
  m_backgroundSubtract = subtract;
  // quit if it's not actually a valid field,
  if( !m_background ) return;
  // and get to removing it:
  const int lastRow = rowCount() - 1;

  // Flag changes in other columns. backgroundSubtract affects data and continuum.
  WModelIndex start = index(       0, DATA_COLUMN );
  WModelIndex end   = index( lastRow, DATA_COLUMN );
  dataChanged().emit( start, end );

  start = index(       0, BACKGROUND_COLUMN );
  end   = index( lastRow, BACKGROUND_COLUMN );
  dataChanged().emit( start, end );
} // void SpectrumDataModel::setBackgroundSubtract( const bool subtract )



// The following 8 functions will return empty std::shared_ptr<Measurement> variables if the
// corresponding data histogram was not set.
//std::shared_ptr<Measurement> SpectrumDataModel::getData()
//{
//  return m_data;
//}
//std::shared_ptr<Measurement> SpectrumDataModel::getSecondData() {  return m_secondData; }
//std::shared_ptr<Measurement> SpectrumDataModel::getBackground() {  return m_background; }

std::shared_ptr<const Measurement> SpectrumDataModel::getData() const       { return m_data;       }
std::shared_ptr<const Measurement> SpectrumDataModel::getSecondData() const { return m_secondData; }
std::shared_ptr<const Measurement> SpectrumDataModel::getBackground() const { return m_background; }

float SpectrumDataModel::dataRealTime() const
{
  return m_dataRealTime;
}


float SpectrumDataModel::dataLiveTime() const
{
  return m_dataLiveTime;
}


float SpectrumDataModel::dataNeutronCounts() const
{
  return m_dataNeutronCounts;
}


float SpectrumDataModel::secondDataScaledBy() const
{
  return m_secondSF;
//  if( m_secondDataLiveTime > 0.0 && m_dataLiveTime > 0.0 )
//    return m_secondSF*m_dataLiveTime/m_secondDataLiveTime;
//  return m_secondSF;
}//float secondDataScaledBy() const


float SpectrumDataModel::backgroundScaledBy() const
{
  return m_backgroundSF;
  
//  if( m_backgroundLiveTime > 0.0 && m_dataLiveTime > 0.0 )
//    return m_backgroundSF*m_dataLiveTime/m_backgroundLiveTime;
//  return m_backgroundSF;
}//float backgroundScaledBy() const


void SpectrumDataModel::setSecondDataScaleFactor( const float sf )
{
  if( !!m_secondData && m_secondSF!=sf )
  {
    m_secondSF = sf;
    WModelIndex secondStart = index( 0, BACKGROUND_COLUMN );
    WModelIndex secondEnd   = index( rowCount() - 1, BACKGROUND_COLUMN );
    dataChanged().emit( secondStart, secondEnd );
  }else
  {
    m_secondSF = sf;
  }
}//void setSecondDataScaleFactor( const float sf )


void SpectrumDataModel::setBackgroundDataScaleFactor( const float sf )
{
  if( !!m_background && m_backgroundSF!=sf )
  {
    m_backgroundSF = sf;
    WModelIndex backStart = index( 0, BACKGROUND_COLUMN );
    WModelIndex backEnd   = index( rowCount() - 1, BACKGROUND_COLUMN );
    dataChanged().emit( backStart, backEnd );
  }else
  {
    m_backgroundSF = sf;
  }
}//void setBackgroundDataScaleFactor( const float sf )


float SpectrumDataModel::secondDataRealTime() const
{
  return m_secondDataRealTime;
}

float SpectrumDataModel::secondDataLiveTime() const
{
  return m_secondDataLiveTime;
}

float SpectrumDataModel::secondDataNeutronCounts() const
{
  return m_secondDataNeutronCounts;
}

float SpectrumDataModel::backgroundRealTime() const
{
  return m_backgroundRealTime;
}

float SpectrumDataModel::backgroundLiveTime() const
{
  return m_backgroundLiveTime;
}

float SpectrumDataModel::backgroundNeutronCounts() const
{
  return m_backgroundNeutronCounts;
}


// The following either return their column or -1 if there is no data.
int SpectrumDataModel::dataColumn() const
{
  return m_data ? DATA_COLUMN : -1;
} // int SpectrumDataModel::dataColumn() const


int SpectrumDataModel::secondDataColumn() const
{
  return m_secondData ? SECOND_DATA_COLUMN : -1;
} // int SpectrumDataModel::secondDataColumn() const


int SpectrumDataModel::backgroundColumn() const
{
  return m_background ? BACKGROUND_COLUMN : -1;
} // int SpectrumDataModel::backgroundColumn() const



// Everything uses the row convention found at the top of the page.
int SpectrumDataModel::rowCount( const WModelIndex & parent) const
{
  if (parent.isValid()) {
    return 0;
  }
  const int rebin = ( m_rebinFactor > 1 ) ? m_rebinFactor : 1;
  // Note: m_rebinFactor >= 1 always, this is for testing.
  std::shared_ptr<const Measurement> xHist = histUsedForXAxis();
  if( !xHist ) return 0;

  return ( !xHist ) ? 0 : (static_cast<int>(xHist->num_gamma_channels()) / rebin);
} // int SpectrumDataModel::rowCount( const Wt::WModelIndex & ) const


int SpectrumDataModel::columnCount( const WModelIndex & ) const
{
  // There are now always 5 columns. This function is left in for flexibility in the future,
  // if the amount of columns ever changes.
  return 5;
} // int SpectrumDataModel::columnCount( const Wt::WModelIndex & ) const


double SpectrumDataModel::rowWidth( int row ) const
{
  std::shared_ptr<const Measurement> xHist = histUsedForXAxis();
  if( !xHist || row < 0 )
    return 0.0;

  const int rebin = ( m_rebinFactor > 1 ) ? m_rebinFactor : 1;
  // Note: m_rebinFactor >= 1 always, this is for testing.

  const size_t firstBin = static_cast<size_t>(row * rebin);
  const size_t lastBin  = firstBin + rebin - 1;

  const double xMin = xHist->gamma_channel_lower( firstBin );
  const double xMax = xHist->gamma_channel_upper( lastBin );

  return xMax - xMin;
} // double SpectrumDataModel::rowWidth( int row ) const


double SpectrumDataModel::rowCenter( int row ) const
{
  return rowLowEdge( row ) + 0.5 * rowWidth( row );
} // double SpectrumDataModel::rowCenter( int row ) const


double SpectrumDataModel::rowLowEdge( int row ) const
{
  std::shared_ptr<const Measurement> xHist = histUsedForXAxis();
  if( !xHist || row < 0 )
    return 0.0;

  const int rebin = ( m_rebinFactor > 1 ) ? m_rebinFactor : 1;
  // Note: m_rebinFactor >= 1 always, this is for testing.

  const size_t bin = static_cast<size_t>( row * rebin );

  return xHist->gamma_channel_lower( bin );
} // double SpectrumDataModel::rowLowEdge( int row ) const


int SpectrumDataModel::findRow( const double x ) const
{
  std::shared_ptr<const Measurement> xHist = histUsedForXAxis();
  if( !xHist )
    return -1;

  const size_t channel = xHist->find_gamma_channel( (float)x );
  
  const int rebin = ( m_rebinFactor > 1 ) ? m_rebinFactor : 1;
  // Note: m_rebinFactor >= 1 always, this is for testing.

  // Return either bin - 1 or the inverse of row * rebin + 1.
  return (rebin == 1) ? channel : (channel / rebin);
} // int findRow( const double x ) const


// { TODO this function can probably be made more efficient }
void SpectrumDataModel::yRangeInXRange( const double xMin, const double xMax,
                            double &yMin, double &yMax ) const
{
  const int lowRow  = findRow( xMin );
  const int highRow = findRow( xMax );
  const int columns = columnCount();

  yMax = -numeric_limits<double>::max();
  yMin =  numeric_limits<double>::max();

  for( int column = 1; column < columns; ++column )
  {
    // Make sure that column is active.
    if( !columnHasData( column ) ) continue;

    for( int row = lowRow; row <= highRow; ++row )
    {
      try
      {
        const double val = boost::any_cast<double>( data( index( row, column ), DisplayRole ) );
        yMax = max( yMax, val );
        yMin = min( yMin, val );
      } catch(...) { }
    } // for( looping over the rows )
  } // for( looping over the columns )

  // WCartesianChart breaks for very large numbers
  if( yMax == -numeric_limits<double>::max() ) yMax = 0.0;
  if( yMin ==  numeric_limits<double>::max() ) yMin = 0.0;
} // double SpectrumDataModel::yRangeInXRange( const double xMin, const double xMax,
  //                               double &yMin, double &yMax ) const


void SpectrumDataModel::setRebinFactor( const int factor )
{
  // The logic of this function has not been tested, and issues with
  // rowsAboutToBeInserted(), rowsAboutToBeRemoved(), endInsertRows(),
  // endRemoveRows(), and dataChanged() are not always obviously
  // manifest in the GUI charts. { TODO address this - (maybye fixed 20130811) }

  int newFactor = ( factor < 1 ) ? 1 : factor;

  const int oldNRows = rowCount();

  if( newFactor != m_rebinFactor )
  {
    int newNRows = 0;
    std::shared_ptr<const Measurement> xHist = histUsedForXAxis();
    if( xHist )
      newNRows = static_cast<int>(xHist->num_gamma_channels()) / newFactor;

    if( oldNRows > newNRows )
    {
      beginRemoveRows( WModelIndex(), newNRows, oldNRows - 1 );
      m_rebinFactor = newFactor;
      endRemoveRows();
    }else
    {
      beginInsertRows( WModelIndex(), oldNRows, newNRows - 1 );
      m_rebinFactor = newFactor;
      endInsertRows();
    } // if( oldNRows > newNRows ) else

    // Refresh every column!
    dataChanged().emit( index( 0, 0 ), index( newNRows - 1, columnCount() - 1 ) );
  } // if( there is a change )
} // void SpectrumDataModel::setRebinFactor( const int factor )


WModelIndex SpectrumDataModel::parent( const WModelIndex & ) const
{
  // Not visibly called, but needed by Wt.
  return WModelIndex();
} // WModelIndex SpectrumDataModel::parent( const WModelIndex &index ) const


double SpectrumDataModel::data( int row, int column ) const
{
  try
  {
    return boost::any_cast<double>( data( index( row, column ), DisplayRole ) );
  }
  catch(...)
  {
    // no-op
  }

  return 0.0;
} // double SpectrumDataModel::data( int row, int column ) const


boost::any SpectrumDataModel::displayBinValue( int row,
                                               SpectrumDataModel::ColumnType column ) const
{
  //Could optimize this function to be a bit more efficient since we no longer
  //  support using TH1Fs with the SpectrumDataModel class
  std::shared_ptr<const Measurement> xHist = histUsedForXAxis();

  if( !xHist )
    return boost::any();

  const size_t nchannel = xHist->num_gamma_channels();
  const int numRows = static_cast<int>(nchannel) / m_rebinFactor;
  if( (row < 0) || (row >= numRows) )
    return boost::any();

  const size_t newxAxisFirstBin = std::min( static_cast<size_t>(row * m_rebinFactor), nchannel-1);
  const size_t newxAxisLastBin = std::min( static_cast<size_t>(newxAxisFirstBin + m_rebinFactor - 1), nchannel-1 );
  
  std::shared_ptr<const Measurement> hist;

  switch( column )
  {
    case X_AXIS_COLUMN:
    {
      const float xMin = xHist->gamma_channel_lower( newxAxisFirstBin );
      const float xMax = xHist->gamma_channel_upper( newxAxisLastBin );
      return boost::any(0.5*(xMax+xMin));
    }
      
    case DATA_COLUMN:        hist = m_data;       break;
    case SECOND_DATA_COLUMN: hist = m_secondData; break;
    case BACKGROUND_COLUMN:  hist = m_background; break;
  }//switch( column )

  if( !hist )
    return boost::any();
  
  double integral = 0.0;
  const vector<float> &channel_contents = *(hist->gamma_channel_contents());
    
  if( hist->channel_energies() == xHist->channel_energies() )
  {
    for( size_t bin = newxAxisFirstBin; bin <= newxAxisLastBin; ++bin )
      integral += channel_contents[bin];
  }else
  {
    const float xmin = xHist->gamma_channel_lower( newxAxisFirstBin );
    const float xmax = xHist->gamma_channel_upper( newxAxisLastBin );
    
    const size_t firstChannel = hist->find_gamma_channel( xmin );
    const size_t lastChannel = hist->find_gamma_channel( xmax );
      
    for( size_t bin = firstChannel; bin <= lastChannel; ++bin )
    {
      const double value = channel_contents[bin];
    
      if( bin == firstChannel )
      {
        const double width = hist->gamma_channel_width( bin );
        const double upperx = ((bin==lastChannel)
                              ? xmax : hist->gamma_channel_upper( bin ));
        integral += value * (upperx-xmin) / width;
      }else if( bin == lastChannel )
      {
        const double width = hist->gamma_channel_width( bin );
        const double lowerx = hist->gamma_channel_lower( bin );
        integral += value * (xmax-lowerx) / width;
      }else
      {
        integral += value;
      }
    }//for( int bin = firstBin; bin <= lastBin; ++bin )
      
#if( PERFORM_DEVELOPER_CHECKS )
    const double otherintegral = hist->gamma_integral( xmin, xmax );
    const double denom = std::max( fabs(otherintegral), fabs(integral) );
    if( denom > 1.0E-6 && (fabs(otherintegral - integral)/denom)>1.0E-5 )
    {
      char buffer[512];
      snprintf( buffer, sizeof(buffer),
                "Gamma bin contents varied from using gamma_integral(...);"
                " %f (this method) vs %f (summing); xlow=%f, xhigh=%f" ,
                integral, otherintegral, xmin, xmax );
      log_developer_error( __func__, buffer );
    }//if( error is too large )
#endif
  }//if( hist->channel_energies() == xHist->channel_energies() ) / else
  
  
  switch( column )
  {
    case X_AXIS_COLUMN:      break;
    case DATA_COLUMN:        break;
    case SECOND_DATA_COLUMN: integral *= secondDataScaledBy(); break;
    case BACKGROUND_COLUMN:  integral *= backgroundScaledBy(); break;
  }//switch( column )

  return boost::any( integral );
}//displayBinValue(...)


boost::any SpectrumDataModel::data( const WModelIndex &index, int role ) const
{
  if( role != Wt::DisplayRole )
    return boost::any();

  const int row    = index.row();
  const int column = index.column();

  const int rebin = ( m_rebinFactor > 1 ) ? m_rebinFactor : 1;
  // Note: m_rebinFactor >= 1 always, this is for testing.

  std::shared_ptr<const Measurement> xHist = histUsedForXAxis();
  const int numRows    = xHist ? static_cast<int>(xHist->num_gamma_channels()) / rebin : 0;
  const int numColumns = columnCount();

  // Make sure it's within standard bounds.
  if( ( row    < 0 ) || ( row    >= numRows    ) )
    return boost::any();
  if( ( column < 0 ) || ( column >= numColumns ) )
    return boost::any();

  if( column == X_AXIS_COLUMN )
    return displayBinValue( row, X_AXIS_COLUMN );

  if( (column == DATA_COLUMN) && !!m_data )
  {
    if( !m_backgroundSubtract || !m_background )
      return displayBinValue( row, DATA_COLUMN );
    
    const double value = asNumber( displayBinValue( row, DATA_COLUMN ) )
                         - asNumber( displayBinValue( row, BACKGROUND_COLUMN ) );
    return boost::any( value );
  }else if( (column == SECOND_DATA_COLUMN) && !!m_secondData )
  {
    if( !m_backgroundSubtract || m_secondDataOwnAxis || !m_background )
      return displayBinValue( row, SECOND_DATA_COLUMN );
    
    const double value = asNumber( displayBinValue( row, SECOND_DATA_COLUMN ) )
                         - asNumber( displayBinValue( row, BACKGROUND_COLUMN ) );
    return boost::any( value );
  }else if( (column == BACKGROUND_COLUMN) && m_background )
  {
    return displayBinValue( row, BACKGROUND_COLUMN );
  }
  
  return boost::any();
} //boost::any SpectrumDataModel::data( const WModelIndex &index, int role ) const


WModelIndex SpectrumDataModel::index( int row, int column, const WModelIndex & ) const
{
  std::shared_ptr<const Measurement> xHist = histUsedForXAxis();

  const int numRows    = xHist ? static_cast<int>(xHist->num_gamma_channels()) : -1;

  // Check to see if it's within bounds for the bins
  if( ( row < 0 ) || ( row >= numRows ) )
    return WModelIndex();

  // Check to make sure that it's trying to retrieve valid data.
  if( !columnHasData( column ) )
    return WModelIndex();

  return WAbstractItemModel::createIndex( row, column, (void *)this );
} // WModelIndex SpectrumDataModel::index( int row, int column, const WModelIndex & ) const


bool SpectrumDataModel::columnHasData( int column ) const
{
  switch( column )
  {
  case X_AXIS_COLUMN:
    if( !histUsedForXAxis() ) return false;
    break;
  case DATA_COLUMN:
    if( !m_data )             return false;
    break;
  case SECOND_DATA_COLUMN:
    if( !m_secondData )       return false;
    break;
  case BACKGROUND_COLUMN:
    if( !m_background )       return false;
    break;
  default:
    return false;
  }

  // If none of those flagged it, we're good.
  return true;
}


boost::any SpectrumDataModel::headerData( int section, Orientation orientation, int role ) const
{
  // If orientation is horizontal, the section is a column (histogram) number.
  // If orientation is vertical,   the section is a row    (bin)       number.
  if( role == LevelRole ) return 0;

  if( ( orientation != Horizontal ) || ( role != DisplayRole ) )
    return WAbstractItemModel::headerData( section, orientation, role );

  if( section == X_AXIS_COLUMN )
    return WString( "Energy (keV)" );

  else if( section == DATA_COLUMN )
  {
    // If there is data and it's trying to add it to the legend
    if( m_data && m_addHistIntegralToLegend )
    {
      char buffer[32];
      snprintf( buffer, sizeof(buffer), " (%.2g counts)", m_data->gamma_count_sum() );
      return WString(  m_data->title() + buffer );
    }//if( m_data && m_addHistIntegralToLegend )
    
    // If the data exists, get its title
    else if( m_data )
      return WString( m_data->title() );
    // Otherwise, just print Data
    else
      return WString( "Data" );
  } // else if( section == DATA_COLUMN )

  else if( m_secondData && section == SECOND_DATA_COLUMN )
  {
    // Grab the title, if there is one
    string title = m_secondData->title();
    if( title == "" )
      title = "Second Data";

    // If there is data and it's trying to add it to the legend
    if( m_addHistIntegralToLegend )
    {
      char buffer[32];
      snprintf( buffer, sizeof(buffer), " (%.2g counts)", (secondDataScaledBy()*m_secondData->gamma_count_sum()) );
      return WString( title + buffer );
    }//if( m_addHistIntegralToLegend )

    // Otherwise, just give back the title.
    return WString( title );
  } // else if( m_secondData && section == SECOND_DATA_COLUMN )

  else if( m_background && section == BACKGROUND_COLUMN )
  {
    // If there is data and it's trying to add it to the legend
    if( m_addHistIntegralToLegend )
    {
      char buffer[32];
      snprintf( buffer, sizeof(buffer), " (%.2g counts)", (backgroundScaledBy()*m_background->gamma_count_sum()) );
      return WString( m_background->title() + buffer );
    }
    
    // Otherwise, just give back the title.
    return WString( m_background->title() );
  } // else if( m_background && section == BACKGROUND_COLUMN )
  
  // This should never happen, but just in case.
  return boost::any();
} // boost::any SpectrumDataModel::headerData( int section, Orientation orientation, int role ) const


void SpectrumDataModel::reset()
{
  // Reset all of the data fields
  setDataHistogram( nullptr );
  setSecondDataHistogram( nullptr, false );
  setBackgroundHistogram( nullptr );

  modelReset().emit();
} // voit SpectrumDataModel::reset()

bool SpectrumDataModel::secondDataOwnAxis() const
{
  return ( m_secondData && m_secondDataOwnAxis );
} // bool SpectrumDataModel::secondDataOwnAxis() const


void SpectrumDataModel::setForegroundSpectrumColor( const Wt::WColor &color )
{
  m_seriesColors[0] = color.isDefault() ? ns_default_foreground_color : color;
}
void SpectrumDataModel::setBackgroundSpectrumColor( const Wt::WColor &color )
{
  m_seriesColors[1] = color.isDefault() ? ns_default_background_color : color;
}
void SpectrumDataModel::setSecondarySpectrumColor( const Wt::WColor &color )
{
  m_seriesColors[2] = color.isDefault() ? ns_default_secondary_color : color;
}

vector< Chart::WDataSeries > SpectrumDataModel::suggestDataSeries() const
{
  const int binwidth = 1;
  
  vector< Chart::WDataSeries > answer;

  // Potentially add in all of the datafields
  if( m_background )
  {
    Chart::WDataSeries series( BACKGROUND_COLUMN, Chart::LineSeries );
    series.setStacked( false );
    
    // Force a black line and other line-based variables
    WPen dataPen( m_seriesColors[1] );
    dataPen.setWidth( WLength( binwidth ) );
    series.setPen( dataPen );
    series.setXSeriesColumn( 0 );
    answer.push_back( series );
  } // if( m_background )

  if( m_secondData )
  {
    Chart::WDataSeries series( SECOND_DATA_COLUMN, Chart::LineSeries );
    series.setStacked( false );
    
    // Force a darkGreen line and other line-based variables
    WPen dataPen( m_seriesColors[2] );
    dataPen.setWidth( WLength( binwidth ) );
    series.setPen( dataPen );
    series.setXSeriesColumn( 0 );

    // If it's on its own, note that.
    if( m_secondDataOwnAxis )
      series.bindToAxis( Chart::Y2Axis );

    answer.push_back( series );
  } // if( m_secondData )

  {
    Chart::WDataSeries series( DATA_COLUMN, Chart::LineSeries );
    series.setStacked( false );
    
    // Force a black line and other line-based variables
    WPen dataPen( m_seriesColors[0] );
    dataPen.setWidth( WLength( binwidth ) );
    series.setPen( dataPen );
    series.setXSeriesColumn( 0 );
    answer.push_back( series );
  } // nesting for m_data

  return answer;
} // vector< Chart::WDataSeries > SpectrumDataModel::suggestDataSeries() const


std::shared_ptr<const Measurement> SpectrumDataModel::histUsedForXAxis() const
{
  // Go through the priorities.
  if( m_data )
    return m_data;

/*
  if( !!m_background ) return m_background;
  // Second data is last because it's probably on another axis.
  if( !!m_secondData ) return m_secondData;
*/

  return nullptr;
}// std::shared_ptr<const Measurement> SpectrumDataModel::histUsedForXAxis()


void SpectrumDataModel::addIntegralOfHistogramToLegend( const bool doIt )
{
  if( doIt == m_addHistIntegralToLegend ) return;

  m_addHistIntegralToLegend = doIt;
  headerDataChanged().emit( Horizontal, 0, rowCount() );
} // void SpectrumDataModel::addIntegralOfHistogramToLegend( const bool doIt )

int SpectrumDataModel::rebinFactor() const
{
  return m_rebinFactor;
} // int SpectrumDataModel::rebinFactor() const



