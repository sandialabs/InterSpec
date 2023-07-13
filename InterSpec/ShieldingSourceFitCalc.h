#ifndef ShieldingSourceFitCalc_h
#define ShieldingSourceFitCalc_h
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
#include <memory>
#include <vector>

#include <boost/function.hpp>

#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
#include <boost/optional.hpp>
#endif


//Forward declarations
struct Material;

class PeakDef;
class MaterialDB;
class PopupDivMenu;
class DetectorPeakResponse;

namespace SandiaDecay
{
  struct Nuclide;
  struct Element;
}//namespace SandiaDecay

namespace ROOT
{
  namespace Minuit2
  {
    class MnUserParameters;
  }//namespace Minuit2
}//namespace ROOT

namespace GammaInteractionCalc
{
  enum class GeometryType : int;
  enum class TraceActivityType : int;
  class ShieldingSourceChi2Fcn;
};


/** This namespace is structs that represent the data users input in the `ShieldingSelect` class, and the inputs of
 ShieldingSourceDisplay.h/.cpp
 */
namespace ShieldingSourceFitCalc
{
  /** Enum that classifies the type of source. */
  enum class ModelSourceType : int
  {
    /** A point source at the center of the shielding. */
    Point,
      
    /** A nuclide in the material itself is the source; e.g., a self-attenuating source like U, Pu, Th, etc. */
    Intrinsic,
      
    /** A trace source in a shielding.  Does not effect transport of gammas through the material, but is just a source term. */
    Trace
  };//enum class ModelSourceType
  
  
  struct SourceFitDef
  {
    const SandiaDecay::Nuclide *nuclide = nullptr;
    
    //activity: in units of PhysicalUnits
    double activity;
    bool fitActivity;
    
    //age: in units of PhysicalUnits::second
    double age;
    bool fitAge;
    
    //ageDefiningNuc: specifies if the age of nuclide should be tied to the age
    //  a different nuclide instead.  Will be NULL if this is not the case.
    const SandiaDecay::Nuclide *ageDefiningNuc;
        
    ModelSourceType sourceType;
    
  #if( INCLUDE_ANALYSIS_TEST_SUITE )
      boost::optional<double> truthActivity, truthActivityTolerance;
      boost::optional<double> truthAge, truthAgeTolerance;
  #endif
    
    SourceFitDef();
    
    virtual void deSerialize( const ::rapidxml::xml_node<char> *parent_node );
    virtual ::rapidxml::xml_node<char> *serialize( rapidxml::xml_node<char> *parent_node ) const;
    
    static const int sm_xmlSerializationMajorVersion;
    static const int sm_xmlSerializationMinorVersion;
  };//struct SourceFitDef
    
    
    
  struct IsoFitStruct : public SourceFitDef
  {
    //numProdigenyPeaksSelected: The number of different progeny selected to be included in the fit
    //  through all the peaks with the parent nuclide as assigned.
    size_t numProgenyPeaksSelected;
      
    //ageIsNotFittable: update this whenever you set the nuclide.  Intended to
    //  indicate nuclides where the spectrum doesnt change with time (ex Cs137,
    //  W187, etc).  Not rock solid yet (not set true as often as could be), but
    //  getting there.  See also PeakDef::ageFitNotAllowed(...).
    bool ageIsFittable;
      
    double activityUncertainty;
    double ageUncertainty;
      
    IsoFitStruct();
    
    virtual void deSerialize( const ::rapidxml::xml_node<char> *parent_node );
    virtual ::rapidxml::xml_node<char> *serialize( rapidxml::xml_node<char> *parent_node ) const;
  };//struct IsoFitStruct
    
  
  
  /** Struct holding information cooresponding to the `TraceSrcDisplay` class defined in ShieldingSelect.cpp;
   represents information about a trace-source in a shielding (e.g., a volumetric source distributed unifrmly in a shielding
   material, but does not effect the attenuation or density of that material.
   */
  struct TraceSourceInfo
  {
    GammaInteractionCalc::TraceActivityType m_type;
    bool m_fitActivity;
    const SandiaDecay::Nuclide *m_nuclide;
    double m_activity; //units are according to #m_type
    float m_relaxationDistance; //only applicable to #TraceActivityType::ExponentialDistribution
    
    TraceSourceInfo();
    void serialize( rapidxml::xml_node<char> *parent_node ) const;
    void deSerialize( const rapidxml::xml_node<char> *shield_node );
    
#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
    static void equalEnough( const TraceSourceInfo &lhs, const TraceSourceInfo &rhs );
#endif
  };//struct TraceSourceInfo

  
  /** Holds the information corresponding to the `ShieldingSelect` widget. */
  struct ShieldingInfo
  {
    /** The geometry of the shielding.
     If a generic material, then must be #GammaInteractionCalc::GeometryType::NumGeometryType.
     */
    GammaInteractionCalc::GeometryType m_geometry;
    
    /** Wether is a generic shielding (i.e., no physical extent, but specified by atomic number and areal density),
     or a physical shielding, that should have a valid material.
     */
    bool m_isGenericMaterial;
    
    /** A kinda vestigial variable that indicates if the user input that made this info was Shielding/Source fit tool,
     or somewhere else where fitting the material dimesnions or AN/AD wasnt intended.
     */
    bool m_forFitting;
    
    /** Material shielding is made out of; will be nullptr if generic material. */
    std::shared_ptr<const Material> m_material;
    
    /** Dimensions of this shielding; the meaning of the entries differs depending on the geometry,
     or if a generic material.
     
     Spherical: ['Thickness', n/a, n/a]
     Cylinder:  ['Radius','Length',n/a]
     Rectangle: ['Width','Height','Depth']
     Generic:   ['AtomicNumber','ArealDensity',n/a]
     */
    double m_dimensions[3];
    bool m_fitDimensions[3];
    
#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
    boost::optional<double> m_truthDimensions[3];
    boost::optional<double> m_truthDimensionsTolerances[3];
    std::map<const SandiaDecay::Nuclide *,std::pair<double,double>> m_truthFitMassFractions;
#endif
    
    // Self-atten source stuff
    bool m_fitMassFrac;
    
    /** Nuclide mass-fractions are only fit within the same element, and they are constrained to be the sum fraction
     of all the nuclides, for that element, that are being fit.
     
     That is, if you have {{I131, 0.1}, {I124,0.01},{Cs131,0.01},{Cs137,0.02}} (and 0.89 stable I127, and 0.97 Cs133), then
     if you fit mass fractions, I131+I124 will always sum to be 0.11 fraction of the Iodine, and Cs137+Cs131 will always
     sum to be 0.03 of the Cesium.
     */
    std::map<const SandiaDecay::Nuclide *,double> m_nuclideFractions;
    
    // Trace-source stuff
    std::vector<TraceSourceInfo> m_traceSources;
    
    
    ShieldingInfo();
    
    rapidxml::xml_node<char> *serialize( rapidxml::xml_node<char> *parent_node ) const;
    void deSerialize( const rapidxml::xml_node<char> *shield_node, MaterialDB *materialDb );
    
    /** Encodes current tool state to app-url format.  Returned string is just the query portion of URL;
     so will look something like "V=1&G=S&D1=1.2cm", and it will not be url-encoded.
     
     TODO: Currently does not encode self-attenuating source, trace-source information, or "truth" values.
     */
    std::string encodeStateToUrl() const;
    
    void handleAppUrl( std::string query_str, MaterialDB *materialDb );
    
#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
    static void equalEnough( const ShieldingInfo &lhs, const ShieldingInfo &rhs );
#endif
    
    static const int sm_xmlSerializationMajorVersion;
    static const int sm_xmlSerializationMinorVersion;
  };//struct ShieldingInfo
  
  
  /** Struct to capture uncertainties associated with fi*/
  struct FitShieldingInfo : ShieldingInfo
  {
    /** 1-sigma uncertainties of fit dimensions or AN/AD. */
    double m_dimensionUncerts[3];
    
    /** 1-sigma uncertainties for fit mass fractions. */
    std::map<const SandiaDecay::Nuclide *,double> m_nuclideFractionUncerts;
    
    /** 1-sigma uncertainties for fit trace source activities (duplicate information available in #IsoFitStruct) */
    std::map<const SandiaDecay::Nuclide *,double> m_traceSourceActivityUncerts;
    
    FitShieldingInfo();
  };//struct FitShieldingInfo : ShieldingInfo
  
  
  /** A struct to share progress of fitting a model; gives current best chi2 (and parameters giving that),
   as well as  elapsed time, and number of function calls.
   
   Use the mutex to access any member variables.
   
   You will usually create a shared pointer to this object, and pass it to the `fit_model` function, as
   capture the pointer for your `progress_fcn` passed to `fit_model`.
   */
  struct ModelFitProgress
  {
    std::mutex m_mutex;
    double chi2;
    double elapsedTime;
    size_t numFcnCalls;
    std::vector<double> parameters;
    
    ModelFitProgress();
  };//struct ModelFitProgress
  
  
  /** A struct to store the results of the model fit.
   
   */
  struct ModelFitResults
  {
    std::mutex m_mutex;
    
    enum class FitStatus{ UserCancelled, TimedOut, InvalidOther, InterMediate, Final };
    FitStatus successful;
    
    double edm;  //estimated distance to minimum.
    double chi2;
    int num_fcn_calls;
    std::vector<double> paramValues;
    std::vector<double> paramErrors;
    std::vector<std::string> errormsgs;
    
    std::vector<PeakDef> foreground_peaks;
    std::vector<PeakDef> background_peaks;
    std::vector<ShieldingSourceFitCalc::ShieldingInfo> initial_shieldings;
    std::vector<ShieldingSourceFitCalc::FitShieldingInfo> final_shieldings;
    
    std::vector<ShieldingSourceFitCalc::IsoFitStruct> fit_src_info;
  };//struct ModelFitResults
    
  
  /** Function that does the actual model fitting; does not need to be in the main GUI thread.
   
      \param wtsession The Wt session id of the current WApplication.
             If an empty string, then the `progress_fcn` and `finished_fcn`
             function will be called from the current thread; otherwise will
             post the calls of these functions to the Wt a
      \param chi2Fcn The initialized ShieldingSourceChi2Fcn object
      \param inputPrams The fit input parameters as filled out by #shieldingFitnessFcn
      \param progress Pointer to location to put the intermediate status of the
             fit (if desired)
      \param progress_fcn The function that will get called roughly every `sm_model_update_frequency_ms`
             milliseconds.  The `progress` object will be updated, then this function will be called, according to
             `wtsession`.  If using with the Wt GUI, you should wrap this function by WApplication::bind() in case this
             widget gets deleted (and hence why the pointer to ModelFitProgress
             must be passed separately, so you can have that in the WApplication
             bind call which must be done before calling this function).
      \param Pointer to location to put the results.  Should have
             #ModelFitResults::shieldings shielding already filled out
      \param finished_fcn Function to call to post to the WServer so the GUI can
             be updated in the main thread; should be wrapped by
             WApplication::bind() in case this widget gets deleted
             (the wrapped call is #updateGuiWithModelFitResults).
             If the fit is canceled via #ShieldingSourceChi2Fcn::cancelFitWithNoUpdate, then
             this function wont be called.
   */
  void fit_model( const std::string wtsession,
                  std::shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn> chi2Fcn,
                  std::shared_ptr<ROOT::Minuit2::MnUserParameters> inputPrams,
                  std::shared_ptr<ModelFitProgress> progress,
                  boost::function<void()> progress_fcn,
                  std::shared_ptr<ModelFitResults> results,
                  boost::function<void()> finished_fcn );
  
  /** The maximum time (in milliseconds) a model fit can take before the fit is
      aborted.  This generally will only ever be applicable to fits with
      self-attenuators, where there is a ton of peaks, or things go really
      haywire.
   
      Initialized to 120 seconds (e.g., 120*1000)
   */
#ifdef NDEBUG
  // Give up after two minutes for release builds
  const size_t sm_max_model_fit_time_ms = 120*1000;
#else
  // For debug builds we'll let it go 7 times longer, which is about the debug slow down
  const size_t sm_max_model_fit_time_ms = 7*120*1000;
#endif

  
  /** How often (in milliseconds) to update the GUI during a model fit.
      This generally will only ever be applicable to fits with self-attenuators.
      
      Initialized to 2000 (e.g., every two seconds)
   */
  const size_t sm_model_update_frequency_ms = 2000;
}//namespace ShieldingSourceFitCalc

#endif //ShieldingSourceFitCalc_h
