#ifndef SearchMode3DModel_H
#define SearchMode3DModel_H
/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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

#include <boost/multi_array.hpp>

#include <Wt/WModelIndex>
#include <Wt/WAbstractTableModel>

//a forward declaration
class InterSpec;

class SearchMode3DDataModel : public Wt::WAbstractTableModel
{
  /*
   This class holds the data the SearchMode3DChart will plot; we will set this
   data using the spectrum the user is displaying.
   
   We have a slight problem in that we should use a scatter plot to represent
   our data, rather than a category plot, because a category plot makes each
   row/column the same size, however the time (row) or energy (column) range of
   may vary for our data.  We will therfore use a surface scatter plot to
   represent the data.  However, we want the chart to look like a bar chart,
   with flat bar tops, so we cant just tell Wt the counts for each time and
   energy value, because this would be the point in the center of the bin, which
   the mesh would then connect each of theses center values, and therefore be
   hard for analysts to intrepret how many counts are in that bin.
   We will therefore define counts at both the lower energy, and upper energy
   of each channel, as well as the lower time and upper time of each sample,
   meaning we will actually plot 4 points for each bin.  This is maybe a bit
   confusing, so we can talk about this.
   
   Possible TODO Items:
     -Make this model return counts per second, instead of counts.  This is 
      so count rates of time bins of differening intervals can be directly 
      compared.
     -Make it so if greater than a certain number of data points are defined (by
      having rowCount()*columnCount() be large), then only a single value per
      bin gets returned, instead of the 4 now.
     -Add coloring of the bins to the model.
     -Add mechanism to detect if an update actually needs to be performed when
      update(...) is called (e.g. if the foreground SpecMeas has actually 
      changed).
     -Improve how time samples are combined, so that if there is one large 
      background time sample, then dont combine that one with the other time 
      samples.  Or maybe dont include the background samples at all, but this 
      will require some further thought.
   */
  
public:
  //SearchMode3DDataModel constructor: parent passed in is for memory managment
  //  purposes only (when its constructor gets called, this
  //  SearchMode3DDataModel will be deleted).
  SearchMode3DDataModel( Wt::WObject *parent = 0 );
  
  //Number of rows in this model will be twice the number of time-samples in
  //  the data (the left and right, or lower and upper, values for each time
  //  sample) plus one more for the end time of the last sample, plus 1 more
  //  that is used to specify energy of the channels (row==0 values are the
  //  energy values).
  virtual int rowCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  
  //Number of columns in this model will be twice the number of energy
  //  channels in the data (the lower and upper energy values for each
  //  channel), plus another one for the upper energy of the last channel,
  //  plus 1 more (channel==0 values are the time values).
  //  m_energies contians the lower energies for each of the channels, plus
  //  the upper energy of the last channel.
  virtual int columnCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  
  //data(): returns the energy and time values for the axiss, as well as the
  //  counts in each bin.  If row==0, then the energy for the channel=column/2
  //  will be reurned.  If column==0, then the time for sample=row/2 will be
  //  returned.
  //  Otherwise the counts for sample=(row-1)/2, and channel=(column-1)/2 will
  //  be returned.
  //  If invalid row/column is specified, or role!=DisplayRole, then an empty
  //  boost::any() will be returned.  WModelIndex 'parent' should always be
  //  invalid (e.g. blank).
  virtual boost::any data( int row, int column,
                           int role = Wt::DisplayRole,
                           const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  
  //data(): equivalent to other data() functions.
  virtual boost::any data( const Wt::WModelIndex &index,
                           int role = Wt::DisplayRole ) const;
  
  //headerData(): not implemented.
  virtual boost::any headerData( int section,
                                 Wt::Orientation orientation = Wt::Horizontal,
                                 int role = Wt::DisplayRole) const;
  
  //minTime()/maxTime() return the min/max x-axis values of current data
  float minTime() const;
  float maxTime() const;

  //minEnergy()/maxEnergy() return the min/max y-axis values of current data
  float minEnergy() const;
  float maxEnergy() const;
  
  //minCounts()/maxCounts() return min/max z-axix values of current data
  float minCounts() const;
  float maxCounts() const;
  
  
  //Returns the min and max counts within a given time and energy range.
  //  20180131 - NOT TESTED WELL, and there apears to be a Wt bug for toggling between log and linear views
  std::pair<float,float> minMaxCounts( const float time_min, const float time_max,
                                       const float e_min, const float e_max ) const;
  
  
  
  //maxNumTimeSamples(): returns the maximum number of time samples the model
  //  is currently configured to allow.  If there are more time samples than
  //  this for the current data, then they are combined according to the
  //  comments for m_maxNumSamples.
  int maxNumTimeSamples() const;
  
  //setMaxNumTimeSamples(): sets the maximum number of time samples this model
  //  will display.  See comments for m_maxNumSamples.
  //  'num' must be 1 or greater, or an exception will be thrown.
  //SearchMode3DDataModel::update(InterSpec*) must be called before changes
  //  will take effect.
  //Default value is 60.
  void setMaxNumTimeSamples( const int num );
  
  //maxNumEnergyChannels(): returns the maximum number of energy channels the
  //  model is currently configured to allow.  If the data contains more energy
  //  channels than this number, they are combined as explained in the comments
  //  for m_maxNumChannels.
  int maxNumEnergyChannels() const;
  
  //setMaxNumEnergyChannels(): sets the maximum number of energy channels that
  //  will be displayed.  'num' must be 1 or greater, or an exception will be
  //  thrown. See comments for m_maxNumChannels.
  //SearchMode3DDataModel::update(InterSpec*) must be called before changes
  //  will take effect.
  //Default value is 128
  void setMaxNumEnergyChannels( const int num );
  
  //update(): updates the current data, based off of the foreground spectra
  //  of the passed in InterSpec by first removing all the data (and
  //  emmitting appropriate signals to this effect) and then adding all the new
  //  data (again with appropriate signals being emmitted).
  void update( InterSpec *viewer );
  
  
protected:
  //m_minCounts: holds the current minimum number of counts of any data bin
  float m_minCounts;
  
  //m_maxCounts: holds the current maximum number of counts of any data bin
  float m_maxCounts;
  
  //m_maxNumChannels: the maximum number of energy channels the chart should
  //  display.  If the data has more memory channels than this, then they are
  //  combined together (by progressive factors of two) until there is less than
  //  or equal to this many channels.
  //Defaults to 128 channels.
  int m_maxNumChannels;
  
  //m_maxNumSamples: the maximum number of time samples that the chart should
  //  display.  If the data has more data samples that this, then the channels
  //  will be combined to have as close to this number as possible; this may
  //  mean the last time row will have a differnt number of actual time samples
  //  than the rest.
  //Defaults to 60 samples.
  int m_maxNumSamples;
  
  //m_times: holds the time in seconds (starting from when first sample was
  //  started), of each sample start time.  Indexed as m_times[row/2].  Note
  //  that this vector will have one more element than the number of samples
  //  in the data (the end time of the last sample).
  //  You could also maybe store and return WDateTimes instead.
  std::vector<float> m_times;
  
  //m_energies: holds the lower energy of each channel.  Index as
  //  m_energies[column/2].  Note that this vector will have one more element
  //  than than the number of energy channels in the spectrum (the upper energy
  //  of the last channel).
  std::vector<float> m_energies;
  
  //m_counts: holds the number of counts detected for a given time sample
  //  and energy channel.  Indexed as m_counts[row/2][column/2],
  //  or equivalently m_counts[sample_number][energy_channel]
  boost::multi_array<float, 2> m_counts;
};//class SearchMode3DDataModel


#endif
