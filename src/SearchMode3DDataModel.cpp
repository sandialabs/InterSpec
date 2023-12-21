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

#include <set>
#include <string>
#include <vector>
#include <cfloat>

#include <Wt/WColor>
#include <Wt/WString>
#include <Wt/WModelIndex>
#include <Wt/WAbstractTableModel>

#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "SpecUtils/EnergyCalibration.h"
#include "InterSpec/SearchMode3DDataModel.h"

using namespace Wt;
using namespace std;


SearchMode3DDataModel::SearchMode3DDataModel( WObject *parent )
  : WAbstractTableModel( parent ),
  m_minCounts( 0.0f ),
  m_maxCounts( 0.0f ),
  m_maxNumChannels( 128 ),
  m_maxNumSamples( 60 )
{
}


int SearchMode3DDataModel::rowCount( const Wt::WModelIndex &parent ) const
{
  if( m_times.empty() )
    return 0;    
  return 2*static_cast<int>(m_times.size()) - 1;
}
  
  
int SearchMode3DDataModel::columnCount( const Wt::WModelIndex &parent ) const
{
  if( m_energies.empty() )
    return 0;
    
  return 2*static_cast<int>(m_energies.size()) - 1;
}
  

boost::any SearchMode3DDataModel::data( int row, int column, int role,
                                        const WModelIndex &parent ) const
{
  return data( createIndex(row, column, (void*)0), role );
}
  
  
boost::any SearchMode3DDataModel::data( const WModelIndex &index,
                                        int role ) const
{
  if( !index.isValid() )
    return boost::any();
    
  const int row = index.row();
  const int column = index.column();
  
  if( role == MarkerBrushColorRole )
    return boost::any();  //Kevin: you might want to customize the color here, or you could create your own version of WStandardColorMap
  else if( role != DisplayRole )
    return boost::any();
    
  //Lets check to make sure data requested is for valid data
  const int nrow = rowCount();
  const int ncolumn = columnCount();
    
  if( row < 0 || column < 0 || row >= nrow || column >= ncolumn )
    return boost::any();
    
  if( row == 0 )  //for row==0, we will return the energy for this column
    return (column==0 ? boost::any() : boost::any(m_energies[column/2]));
    
  if( column == 0 )  //for column == 0, we will return the time (since measurment started) for this row
    return boost::any(m_times[row/2]);
    
  const int sample = (row-1) / 2;
  const int channel = (column-1) / 2;
    
  if( channel == (m_energies.size()-1) || sample == (m_times.size()-1) )
    return boost::any( 0.0 );
  
  if( channel >= m_energies.size() || sample >= m_times.size() )
    return boost::any();
    
  return boost::any( m_counts[sample][channel] );
}//data()
  
  
boost::any SearchMode3DDataModel::headerData( int section,
                                              Orientation orientation,
                                              int role ) const
{
  return 0.0;  //un-implemented
}
  
float SearchMode3DDataModel::minCounts() const { return m_minCounts; }
float SearchMode3DDataModel::maxCounts() const { return m_maxCounts; }
float SearchMode3DDataModel::minEnergy() const { return (m_energies.empty() ? 0.0f : m_energies[0]); }
float SearchMode3DDataModel::maxEnergy() const { return (m_energies.empty() ? 1.0f : m_energies.back()); }
float SearchMode3DDataModel::minTime() const { return (m_times.empty() ? 0.0f : m_times[0]); }
float SearchMode3DDataModel::maxTime() const { return (m_times.empty() ? 1.0f : m_times.back()); }
int SearchMode3DDataModel::maxNumTimeSamples() const { return m_maxNumSamples; }
int SearchMode3DDataModel::maxNumEnergyChannels() const { return m_maxNumChannels; }


std::pair<float,float> SearchMode3DDataModel::minMaxCounts( const float time_min, const float time_max,
                                                            const float e_min, const float e_max ) const
{
  std::pair<float,float> answer( 0.0f, 1.0f );
  
  auto timeStart = lower_bound( begin(m_times), end(m_times), time_min );
  auto timeEnd = upper_bound( begin(m_times), end(m_times), time_max );
  
  auto eStart = lower_bound( begin(m_energies), end(m_energies), e_min );
  auto eEnd = upper_bound( begin(m_energies), end(m_energies), e_max );
  
  if( eStart == end(m_energies) )
    return answer;
  
  if( timeStart == end(m_times) )
    return answer;
  
  if( timeEnd != end(m_times) && ((*timeEnd) < (*timeStart)) )
    swap( timeEnd, timeStart );
  
  if( eEnd != end(m_energies) && ((*eEnd) < (*eStart)) )
    swap( eEnd, eStart );
  
  const size_t start_time_index = (timeStart - begin(m_times));
  const size_t end_time_index = 1 + (timeEnd - begin(m_times));
  
  const size_t start_energy_index = (eStart - begin(m_energies));
  const size_t end_energy_index = 1 + (eEnd - begin(m_energies));
  
  const size_t num_samples = (m_times.size() > 0 ? (m_times.size() - 1) : size_t(0));
  const size_t num_energies = (m_energies.size() > 0 ? (m_energies.size() - 1) : size_t(0));
  
  answer.first = FLT_MAX;
  answer.second = -FLT_MAX;
  for( size_t sample = start_time_index; sample < end_time_index && sample < num_samples; ++sample )
  {
    for( size_t channel = start_energy_index; channel < end_energy_index && channel < num_energies; ++channel )
    {
      const float value = m_counts[sample][channel];
      answer.first = std::min( answer.first, value );
      answer.second = std::max( answer.second, value );
    }
  }
  
  return answer;
}//minMaxCounts(...)


void SearchMode3DDataModel::setMaxNumTimeSamples( const int num )
{
  if( num < 1 )
    throw runtime_error( "SearchMode3DDataModel::setMaxNumTimeSamples(num): "
                         "num must be >= 1" );
  m_maxNumSamples = num;
}


void SearchMode3DDataModel::setMaxNumEnergyChannels( const int num )
{
  if( num < 1 )
    throw runtime_error( "SearchMode3DDataModel::setMaxNumEnergyChannels(num): "
                         "num must be >= 1" );
  m_maxNumChannels = num;
}


void SearchMode3DDataModel::update( InterSpec *viewer )
{
  //If the model currently contains any data, we will remove it, and notify any
  //  views that may be using the model.  Then if there is data to display
  //  we will add it to the model and notify the views of the newly available
  //  data.
  
  const int nrow = rowCount();
  const int ncol = columnCount();
  
  if( nrow )
  {
    beginRemoveRows( WModelIndex(), 0, nrow-1 );
    m_times.clear();
    endRemoveRows();
  }
    
  if( ncol )
  {
    beginRemoveColumns( WModelIndex(), 0, ncol-1 );
    m_energies.clear();
    endRemoveColumns();
  }
  
  if( nrow || ncol )
  {
    boost::multi_array<float, 2>::extent_gen extents;
    m_counts.resize( extents[0][0] );
  }
    
  m_minCounts = 0.0f;
  m_maxCounts = 1.0f;
    
    
  std::shared_ptr<const SpecMeas> meas = viewer->measurment( SpecUtils::SpectrumType::Foreground );
  const vector<string> det_to_use = viewer->detectorsToDisplay(SpecUtils::SpectrumType::Foreground);
  const set<int> sample_numbers = meas->sample_numbers();
  const vector<int> sample_numbers_vec( sample_numbers.begin(), sample_numbers.end() );
    
  //foreground_samples: samples the user has summed to display the spectrum of
  //    const set<int> foreground_samples = viewer->displayedSamples( SpecUtils::SpectrumType::Foreground );
  try
  {
    if( !meas || sample_numbers.empty() || det_to_use.empty() )
      throw runtime_error( "No data to display" );
    
    auto energy_cal = meas->suggested_sum_energy_calibration( sample_numbers, det_to_use );
    
    if( !energy_cal || energy_cal->num_channels() < 4 )
      throw runtime_error( "Not enough gamma channels to plot" );
    
    const size_t nsamples = sample_numbers.size();
    
    size_t nenergies = energy_cal->num_channels();
    size_t ncombine = 1;
    while( (ncombine < nenergies) && (nenergies / ncombine) > m_maxNumChannels )
    {
      ncombine += 1;
    }
      
    if( ncombine != 1 )
    {
      energy_cal = energy_cal_combine_channels( *energy_cal, ncombine );
      nenergies = energy_cal->num_channels();
    }
      
    m_minCounts = FLT_MAX;
    m_maxCounts = 0.0f;
    vector<float> newtimes;
    boost::multi_array<float, 2>::extent_gen extentgen;
    m_counts.resize( extentgen[nsamples][nenergies] );
      
    int samplen = 0;
    double cumulativeRealTime = 0.0;
    
    size_t sampleNumDelta = 1;
    //while loop inefficient, but whatever
    while( (sample_numbers_vec.size()/sampleNumDelta) >= m_maxNumSamples )
      ++sampleNumDelta;
    
    //Think of case where m_maxNumSamples==20, and sample_numbers_vec.size()==100,
    // then sampleNumDelta will be 6, but we want 5.
    if( sampleNumDelta > 1 && ((sample_numbers_vec.size()%sampleNumDelta)==0) )
      --sampleNumDelta;
    
    const size_t numSampleNums = sample_numbers_vec.size();
    for( size_t sampleNumIndex = 0; sampleNumIndex < numSampleNums; sampleNumIndex += sampleNumDelta )
    {
      set<int> thissamplenum;
        
      //Kevin: note that the last displayed time period may not have as many
      //  samples as the other time periods if (sample_numbers.size() % m_maxNumSamples) != 0
      for( int i = 0; (i < sampleNumDelta) && ((sampleNumIndex+i) < numSampleNums); ++i )
        thissamplenum.insert( sample_numbers_vec[sampleNumIndex + i] );
      
      float realtime = FLT_MAX;
      for( const string &detnamme : det_to_use )
      {
        for( const int samplenum : thissamplenum )
        {
          auto m = meas->measurement( samplenum, detnamme );
          if( m )
            realtime = std::min( realtime, m->real_time() );
        }//for( const int samplenum : thissamplenum )
      }//for( const string &detnamme : det_to_use )
        
      if( realtime > 1.0E+6f )
        realtime = 0.0f;
        
      auto summed = meas->sum_measurements( thissamplenum, det_to_use, energy_cal );
        
      if( !summed || !summed->gamma_channel_contents()
          || summed->gamma_channel_contents()->size() != nenergies )
        throw runtime_error( "Summing results have unexpected issues" );
        
      const vector<float> &counts = *summed->gamma_channel_contents();
      for( size_t i = 0; i < nenergies; ++i )
      {
        m_minCounts = std::min( m_minCounts, counts[i] );
        m_maxCounts = std::max( m_maxCounts, counts[i] );
        
        m_counts[samplen][i] = counts[i];
      }
      newtimes.push_back( cumulativeRealTime );
        
      ++samplen;
      cumulativeRealTime += realtime;
    }//for( const int samplenum : sample_numbers )
      
    newtimes.push_back( cumulativeRealTime );
    
    if( (energy_cal->num_channels() > 1) && (newtimes.size() > 1) )  //probably always true
    {
      const auto &new_energies = energy_cal->channel_energies();
      beginInsertColumns( WModelIndex(), 0, int(2*new_energies->size()+1) );
      beginInsertRows( WModelIndex(), 0, int(2*newtimes.size()+1) );
      
      m_times = newtimes;
      m_energies = *new_energies;
  
      endInsertRows();
      endInsertColumns();
    }//if( newenergies.size() && newtimes.size() )
    
    // The tabular Wt widgets are a bit sticky refreshing sometimes, so we'll also force the
    //  refresh here as well, even though I dont know if that issue is applicable.
    //reset();
    layoutAboutToBeChanged().emit();
    layoutChanged().emit();
  }catch( std::exception &e )
  {
      cerr << "SearchMode3DDataModel::update() caught: " << e.what() << endl;
  }
}//void update( InterSpec *viewer )


