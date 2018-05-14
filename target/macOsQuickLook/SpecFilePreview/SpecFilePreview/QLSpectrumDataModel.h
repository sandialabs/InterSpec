#ifndef QLSpectrumDataModel_h
#define QLSpectrumDataModel_h

/* QLSpectrumDataModel is a version of InterSpecs SpectrumDataModel class that was
 branched 20171227 to create a version slightly more acceptable for QuickLook
 functionality.
 */


#include <vector>
#include <memory>

#include <boost/any.hpp>

#include <Wt/WSignal>
#include <Wt/WModelIndex>
#include <Wt/Chart/WDataSeries>
#include <Wt/WAbstractItemModel>

class QLSpectrumChart;

class QLSpectrumDataModel: public Wt::WAbstractItemModel
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
  QLSpectrumDataModel( Wt::WObject *parent = 0 );
  virtual ~QLSpectrumDataModel();

  // If one of the following functions that add or change a histogram would
  // encounter a binning inconsistency, it might throw an std::runtime_error
  
  //setDataHistogram(): also sets m_secondSF and m_backgroundSF to
  //  m_dataLiveTime/m_{Second|Background}LiveTime
  virtual void setDataHistogram( std::shared_ptr<Measurement> hist,
                                 float liveTime,
                                 float realTime,
                                 float neutronCounts );
  
  //setSecondDataHistogram(): also sets m_secondSF to
  //  m_dataLiveTime/m_secondDataLiveTime
  virtual void setSecondDataHistogram( std::shared_ptr<Measurement> hist,
                                       float liveTime,
                                       float realTime,
                                       float neutronCounts,
                                       bool ownAxis );
  
  //setBackgroundHistogram(): also sets m_backgroundSF to
  //  m_dataLiveTime/m_backgroundLiveTime
  virtual void setBackgroundHistogram( std::shared_ptr<Measurement> hist,
                                       float liveTime,
                                       float realTime,
                                       float neutronCounts );
  
  // Flags for whether the background are getting subtracted.
  void setBackgroundSubtract( const bool subtract = true );

  // Returns the values for whether backgrounds are being subtracted.
  bool backgroundSubtract() const;

  // The following functions will return an empty std::shared_ptr<Measurement> if the data histogram
  // was not previously set.
  virtual std::shared_ptr<Measurement>      getData();
  virtual std::shared_ptr<Measurement>      getSecondData();
  virtual std::shared_ptr<Measurement>      getBackground();

  virtual std::shared_ptr<const Measurement> getData()       const;
  virtual std::shared_ptr<const Measurement> getSecondData() const;
  virtual std::shared_ptr<const Measurement> getBackground() const;

  virtual float dataRealTime() const;
  virtual float dataLiveTime() const;
  virtual float dataNeutronCounts() const;

  virtual float secondDataRealTime() const;
  virtual float secondDataLiveTime() const;
  virtual float secondDataNeutronCounts() const; //actual measurment, NOT scaled by secondDataScaledBy();

  virtual float backgroundRealTime() const;
  virtual float backgroundLiveTime() const;
  virtual float backgroundNeutronCounts() const;  //actual measurment, NOT scaled by backgroundScaledBy();

  virtual float secondDataScaledBy() const;  //returns m_secondSF
  virtual float backgroundScaledBy() const;  //returns m_backgroundSF
  
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

  std::vector< Wt::Chart::WDataSeries > suggestDataSeries() const;
  virtual bool secondDataOwnAxis() const;

  // Find which histogram is defining the x-axis by moving through priorities.
  std::shared_ptr<Measurement>      histUsedForXAxis();
  std::shared_ptr<const Measurement> histUsedForXAxis() const;

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
  //  If the bin x-values of the respective std::shared_ptr<Measurement> do not align with
  //  the edges of histUsedForXAxis(), then simple linear interpolation is used.
  //  Also, this function does not take into account background subracting.
  //  An empty boost::any() is returned if data is not avaliable
  virtual boost::any displayBinValue( int row, ColumnType column ) const;


  void addIntegralOfHistogramToLegend( const bool doIt = true );

  virtual int  rebinFactor() const;
  virtual void setRebinFactor( const int factor );
  
  // Checks for whether a given column is actually populated.
  bool columnHasData( int column ) const;

  
  //dataSet(): emitted when a foreground, backgorund, or secondary spectrum is
  //  changed.
  Wt::Signal<ColumnType> &dataSet();
  
protected:
  int        m_rebinFactor;
  std::shared_ptr<Measurement> m_data;
  std::shared_ptr<Measurement> m_secondData;
  std::shared_ptr<Measurement> m_background;
  
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
};//class QLSpectrumDataModel

#endif


