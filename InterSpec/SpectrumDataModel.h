#ifndef SpectrumDataModel_h
#define SpectrumDataModel_h
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

#include <boost/any.hpp>

#include <Wt/WColor>
#include <Wt/WSignal>
#include <Wt/WModelIndex>
#include <Wt/Chart/WDataSeries>
#include <Wt/WAbstractItemModel>

class SpectrumChart;
namespace SpecUtils{ class Measurement; }


class SpectrumDataModel: public Wt::WAbstractItemModel
{
  /********\
  | This class is a holding area for data on the graph. It allows for
  | accessing, modifying, etc. the data points in order to shift data
  | and then draw it.
  |
  | When using any of the WAbstractItemModel or WModelIndex function calls,
  | this class has each column correspond to a piece of data:
  |   0th column: x-axis of the chart
  |   1st column: data histogram
  |   2nd column: second data histogram (optional)
  |   3rd column: background data
  |
  | The given x value for a bin is the center of the bin, not either edge.
  | Note that ROOT TH1s use the lower edge as the x value, and this is _not_
  | the case.
  |
  \********/

public:
  // Some convenient constants for addressing the columns.
  enum ColumnType
  {
    X_AXIS_COLUMN,
    DATA_COLUMN,
    SECOND_DATA_COLUMN,
    BACKGROUND_COLUMN
  };//enum ColumnType


public:
  SpectrumDataModel( Wt::WObject *parent = 0 );
  virtual ~SpectrumDataModel();

  // If one of the following functions that add or change a histogram would
  // encounter a binning inconsistency, it might throw an std::runtime_error
  
  //setDataHistogram(): also sets m_secondSF and m_backgroundSF to
  //  m_dataLiveTime/m_{Second|Background}LiveTime
  void setDataHistogram( std::shared_ptr<SpecUtils::Measurement> hist );
  
  //setSecondDataHistogram(): also sets m_secondSF to m_dataLiveTime/m_secondDataLiveTime
  void setSecondDataHistogram( std::shared_ptr<SpecUtils::Measurement> hist,
                                      const bool ownAxis );
  
  //setBackgroundHistogram(): also sets m_backgroundSF to
  //  m_dataLiveTime/m_backgroundLiveTime
  void setBackgroundHistogram( std::shared_ptr<SpecUtils::Measurement> hist );
  
  // Flags for whether the background are getting subtracted.
  void setBackgroundSubtract( const bool subtract = true );

  // Returns the values for whether backgrounds are being subtracted.
  bool backgroundSubtract() const;

  // The following functions will return an empty std::shared_ptr<SpecUtils::Measurement> if the data histogram
  // was not previously set.
  std::shared_ptr<SpecUtils::Measurement>      getData();
  std::shared_ptr<SpecUtils::Measurement>      getSecondData();
  std::shared_ptr<SpecUtils::Measurement>      getBackground();

  std::shared_ptr<const SpecUtils::Measurement> getData()       const;
  std::shared_ptr<const SpecUtils::Measurement> getSecondData() const;
  std::shared_ptr<const SpecUtils::Measurement> getBackground() const;

  float dataRealTime() const;
  float dataLiveTime() const;
  float dataNeutronCounts() const;

  float secondDataRealTime() const;
  float secondDataLiveTime() const;
  float secondDataNeutronCounts() const; //actual measurment, NOT scaled by secondDataScaledBy();

  float backgroundRealTime() const;
  float backgroundLiveTime() const;
  float backgroundNeutronCounts() const;  //actual measurment, NOT scaled by backgroundScaledBy();

  float secondDataScaledBy() const;  //returns m_secondSF
  float backgroundScaledBy() const;  //returns m_backgroundSF
  
  //Scale factors set in setSecondDataScaleFactor() and
  //  setBackgroundDataScaleFactor() multiply the live time scale factor.
  void setSecondDataScaleFactor( const float sf );
  void setBackgroundDataScaleFactor( const float sf );

  // Functions to grab column numbers. Returns -1 if there is nothing for that dataset.
  int dataColumn() const;
  int secondDataColumn() const;
  int backgroundColumn() const;

  // The following parcel of functions uses the column convention at the beginning
  // of the class. i.e., { x-axis, data, 2nd data, background }.
  virtual int rowCount(    const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  virtual int columnCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  virtual Wt::WModelIndex parent( const Wt::WModelIndex &index ) const;
  virtual boost::any data( const Wt::WModelIndex &index, int role = Wt::DisplayRole ) const;
  virtual Wt::WModelIndex index( int row, int column,
                                 const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  virtual boost::any headerData( int section, Wt::Orientation orientation = Wt::Horizontal,
                                 int role = Wt::DisplayRole ) const;
  virtual void reset();

  /** Changes the color of the series #suggestDataSeries() will return.
      You will need to set the new series into the chart for changes to take
      effect.
   */
  void setForegroundSpectrumColor( const Wt::WColor &color );
  void setBackgroundSpectrumColor( const Wt::WColor &color );
  void setSecondarySpectrumColor( const Wt::WColor &color );
  
  std::vector< Wt::Chart::WDataSeries > suggestDataSeries() const;
  bool secondDataOwnAxis() const;

  // Find which histogram is defining the x-axis by moving through priorities.
  std::shared_ptr<SpecUtils::Measurement>      histUsedForXAxis();
  std::shared_ptr<const SpecUtils::Measurement> histUsedForXAxis() const;

  // The functions that follow are convenience functions that take m_rebinFactor into
  // account and return quantities about the x-axis. These are similar in nature to
  // TH1::FindBin, TH1::GetBinCenter, TH1::GetBinWidth, and TH1::GetBinLowEdge, except
  // they use 0-based indices of 'row' rather than 1-based indices of 'bin.'
  double rowWidth(   int row ) const;
  double rowCenter(  int row ) const;
  double rowLowEdge( int row ) const;


  // This finds the row corresponding to a given x value.
  int findRow( const double x ) const;

  // Similarly, this finds the range of y values in a given x range.
  void yRangeInXRange( const double xMin, const double xMax,
                       double &yMin, double &yMax ) const;

  //Simplified convenience functions to grab data without futzing around with
  // either WModelIndex or boost::any
  virtual double data( int row, int column ) const;

  //displayBinValue(...): returns the data value for the display bin (e.g.
  //  taking m_rebinFactor into account) of the respective column.
  //  If the bin x-values of the respective std::shared_ptr<SpecUtils::Measurement> do not align with
  //  the edges of histUsedForXAxis(), then simple linear interpolation is used.
  //  Also, this function does not take into account background subracting.
  //  An empty boost::any() is returned if data is not avaliable
  virtual boost::any displayBinValue( int row, ColumnType column ) const;


  void addIntegralOfHistogramToLegend( const bool doIt = true );

  int  rebinFactor() const;
  void setRebinFactor( const int factor );
  
  // Checks for whether a given column is actually populated.
  bool columnHasData( int column ) const;

  
  //dataSet(): emitted when a foreground, backgorund, or secondary spectrum is
  //  changed.
  Wt::Signal<ColumnType> &dataSet();
  
protected:
  int        m_rebinFactor;
  std::shared_ptr<SpecUtils::Measurement> m_data;
  std::shared_ptr<SpecUtils::Measurement> m_secondData;
  std::shared_ptr<SpecUtils::Measurement> m_background;
  
  float m_dataLiveTime, m_dataRealTime, m_dataNeutronCounts;
  float m_backgroundLiveTime, m_backgroundRealTime, m_backgroundNeutronCounts;
  float m_secondDataLiveTime, m_secondDataRealTime, m_secondDataNeutronCounts;

  //m_backgroundSF and m_secondSF are the scale factors to be applied to
  //  the background and second spectrums, on top of live time normaliztion
  float m_backgroundSF, m_secondSF;
  
  // This both controls whether the second data owns the axis and whether it should
  // background subtract.
  bool m_secondDataOwnAxis;

  // This does not affect the number of columns. When true, boost::any() is returned for
  // calls to data(...) requesting background or continuum.
  bool m_backgroundSubtract;

  bool m_addHistIntegralToLegend;
  
  Wt::Signal<ColumnType> m_dataSet;
  
  //Foreground==0, Background==1, Secondary==2
  Wt::WColor m_seriesColors[3];
};//class SpectrumDataModel

#endif


