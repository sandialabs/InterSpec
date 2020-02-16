#ifndef IsotopeSearchByEnergyModel_h
#define IsotopeSearchByEnergyModel_h
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

#include <Wt/WString>
#include <Wt/WModelIndex>
#include <Wt/WContainerWidget>
#include <Wt/WAbstractItemModel>

#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/ReactionGamma.h"

class SpecMeas;
class DetectorPeakResponse;
namespace SpecUtils{ class Measurement; }

namespace SandiaDecay
{
  struct Element;
  struct Nuclide;
  struct Transition;
  struct RadParticle;
}//namespace SandiaDecay



//ToDo: move IsotopeSearchByEnergyModel into its own header/source
class IsotopeSearchByEnergyModel : public Wt::WAbstractItemModel
{
public:
  enum Column
  {
    ParentIsotope, Energy, Distance, BranchRatio, ProfileDistance,
    SpecificIsotope, ParentHalfLife, AssumedAge,
    NumColumns
  };//enum Column
  
  enum RadSource
  {
    kGamma    = 0x1,
    kXRay     = 0x2,
    kReaction = 0x4
  };//enum RadSource
  
  struct IsotopeMatch
  {
    IsotopeMatch();
    
    double m_distance;    //sum of distance over all energies searched
    
    double m_age;         //age assumed for listing things
    double m_branchRatio; //branching ratio
    
    double m_profileDistance;
    
    //Only one of the following will be valid: m_nuclide, m_element, m_reaction
    const SandiaDecay::Nuclide     *m_nuclide;
    const SandiaDecay::Transition  *m_transition;
    const SandiaDecay::RadParticle *m_particle;
    PeakDef::SourceGammaType m_sourceGammaType;
    
    const SandiaDecay::Element *m_element;
    const SandiaDecay::EnergyIntensityPair *m_xray;
    
    const ReactionGamma::Reaction *m_reaction;
    ReactionGamma::EnergyAbundance m_reactionEnergy;
    
    Wt::WString m_displayData[NumColumns];
  };//struct IsotopeMatch
  
public:
  IsotopeSearchByEnergyModel( Wt::WObject *parent = 0 );
  virtual ~IsotopeSearchByEnergyModel();
  
  void clearResults();
  
  
  //SearchWorkingSpace: this is a place to store search results, and what led to
  //  them, before changing the data of the IsotopeSearchByEnergyModel.  This
  //  struct is necasary since the search results are computed in a background
  //  thread, and then the results are posted to the WApplications event loop,
  //  of which the call to (including alll input parameters) must be bound
  //  before posting the search to the background thread (for safety against
  //  calling updateSearchResults() after *this has been delted).
  struct SearchWorkingSpace
  {
    std::vector<double> energies;
    std::vector<double> windows;
    std::shared_ptr<const DetectorPeakResponse> detector_response_function;
    std::vector<std::shared_ptr<const PeakDef>> user_peaks;
    std::vector<std::shared_ptr<const PeakDef>> automated_search_peaks;
    std::shared_ptr<const SpecUtils::Measurement> displayed_measurement;
    std::shared_ptr<SpecMeas> foreground;
    std::set<int> foreground_samplenums;
    std::vector< std::vector<IsotopeMatch> > matches;
    Column sortColumn;
    Wt::SortOrder sortOrder;
    boost::function< void(void) > searchdoneCallback;
  };//struct SearchWorkingSpace
  
  void updateSearchResults( std::shared_ptr<SearchWorkingSpace> workingspace );
  static void setSearchEnergies( std::shared_ptr<SearchWorkingSpace> workingspace,
                                const double minbr,
                                const double minHalfLife,
                                Wt::WFlags<RadSource> radiation,
                                const std::string appid,
                                boost::function< void(void) > updatefcn );
  
  
  virtual int columnCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  
  virtual int rowCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  
  virtual Wt::WModelIndex parent( const Wt::WModelIndex &index ) const;
  virtual boost::any data( const Wt::WModelIndex &index,
                          int role = Wt::DisplayRole) const;
  
  //nuclide(...): returns the nuclide IF it is a nuclide defined row, otherwise
  //  NULL
  const SandiaDecay::Nuclide *nuclide( const Wt::WModelIndex &index ) const;
  
  //nuclide(...): returns the element IF it is a xray defined row, otherwise
  //  NULL
  const SandiaDecay::Element *xrayElement( const Wt::WModelIndex &index ) const;
  
  //nuclide(...): returns the reaction IF it is a reaction defined row,
  //  otherwise NULL
  const ReactionGamma::Reaction *reaction( const Wt::WModelIndex &index ) const;
  
  double assumedAge( const Wt::WModelIndex &index ) const;
  
  virtual boost::any headerData( int section,
                                Wt::Orientation orientation = Wt::Horizontal,
                                int role = Wt::DisplayRole ) const;
  
  virtual Wt::WModelIndex index( int row, int column,
                                const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  
  virtual void sort( int column, Wt::SortOrder order = Wt::AscendingOrder );
  
  virtual Wt::WFlags<Wt::ItemFlag> flags( const Wt::WModelIndex &index ) const;
  
  
  
  static void sortData( std::vector< std::vector<IsotopeMatch> > &input,
                       const std::vector<double> &energies,
                       int column, Wt::SortOrder order );
  
  typedef std::map<const SandiaDecay::Nuclide *, std::set<double> > NucToEnergiesMap;
  typedef std::vector< std::vector<IsotopeSearchByEnergyModel::IsotopeMatch> > SearchResults;
  
  //nuclidesWithAllEnergies(...): allows x-rays to be considered gamma rays
  static void nuclidesWithAllEnergies( const NucToEnergiesMap &candidates,
                                      const std::vector<double> &energies,
                                      const std::vector<double> &windows,
                                      const double minbr,
                                      const std::shared_ptr<const DetectorPeakResponse> detector_response_function,
                                      const std::shared_ptr<const SpecUtils::Measurement> displayed_measurement,
                                      const std::vector<std::shared_ptr<const PeakDef>> &user_peaks,
                                      const std::vector<std::shared_ptr<const PeakDef>> &automated_search_peaks,
                                      SearchResults &answer );
  
  //xraysWithAllEnergies(...): fairly inefficient
  static void xraysWithAllEnergies( const std::vector<double> &energies,
                                   const std::vector<double> &windows,
                                   const std::shared_ptr<const DetectorPeakResponse> detector_response_function,
                                   const std::shared_ptr<const SpecUtils::Measurement> displayed_measurement,
                                   const std::vector<std::shared_ptr<const PeakDef>> &user_peaks,
                                   const std::vector<std::shared_ptr<const PeakDef>> &automated_search_peaks,
                                   SearchResults &answer );
  
  //reactionsWithAllEnergies(...): not very well yet.  Does not take into
  //  account the minimum branching ratio desired
  static void reactionsWithAllEnergies( const std::vector<double> &energies,
                                       const std::vector<double> &windows,
                                       const std::shared_ptr<const DetectorPeakResponse> detector_response_function,
                                       const std::shared_ptr<const SpecUtils::Measurement> displayed_measurement,
                                       const std::vector<std::shared_ptr<const PeakDef>> &user_peaks,
                                       const std::vector<std::shared_ptr<const PeakDef>> &automated_search_peaks,
                                       SearchResults &answer );
  
  Column sortColumn() const;
  Wt::SortOrder sortOrder() const;
  
protected:
  Column m_sortColumn;
  Wt::SortOrder m_sortOrder;
  std::vector<double> m_windows;
  std::vector<double> m_energies;
  std::vector< std::vector<IsotopeMatch> > m_matches;
};//IsotopeSearchByEnergyModel


#endif  //ifndef IsotopeSearchByEnergyModel_h
