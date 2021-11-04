#ifndef GammaInteractionCalc_h
#define GammaInteractionCalc_h
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
#include <set>
#include <tuple>
#include <atomic>
#include <vector>
#include <utility>

#include <boost/asio/deadline_timer.hpp>

#include <Wt/WColor>

#include "Minuit2/FCNBase.h"


class PeakDef;
struct Material;
class DetectorPeakResponse;

namespace SandiaDecay
{
  struct Nuclide;
  class NuclideMixture;
}

namespace GammaInteractionCalc
{

/** We will let trace source activities be specified in a few different ways, since it makes sense from a users-perspective, and because
 we may want total activity to vary independently from physical volume and then this just makes things a little easier to track in the Minuit
 variables.
 */
enum class TraceActivityType : int
{
  TotalActivity,
  ActivityPerCm3,
  ActivityPerGram, //Needs to come last as wont be available if a "void" material
  //ActivityPPM,
  NumTraceActivityType
};//enum class TraceActivityType

/** Gives string representation to a TraceActivityType value. */
const char *to_str( TraceActivityType type );



//Returned in units of 1.0/[Length], so that
//  exp( -transmition_length_coefficient(...) * thickness)
//  gives you the probability a gamma of given energy will go through the
//  material of given thickness.
//Energy units should be in SandiaDecay/PhysicalUnits units.
double transmition_length_coefficient( const Material *material, float energy );
double transmition_coefficient_material( const Material *material, float energy,
                                float length );

//Returned in units of [Length]^2/[mass], so that
//  exp( -mass_attenuation_coef * areal_density )
//  gives you the probability a gamma of given energy will go through the
//  material with given atomic number and areal_density.
//  The quantity retuned by this function is commonly labeled μ
//Energy units should be in SandiaDecay/PhysicalUnits units.
double mass_attenuation_coef( float atomic_number, float energy );
double transmition_coefficient_generic( float atomic_number, float areal_density,
                                float energy );


void example_integration();
  
int Integrand( const int *ndim, const double xx[],
                      const int *ncomp, double ff[], void *userdata );


//point_to_line_dist(...): calculates the distance from 'point' to the line
//  passing through 'l0' and 'l1'
double point_to_line_dist( const double point[3],
                           const double l0[3], const double l1[3] );

//exit_point_of_sphere_z(...): Makes exit_point[3]=<x1,y1,z1> be the point of
//  intersection of a sphere or radius 'sphere_rad' centered at the origin, with
//  the line pointing from source_point[3]=<x0,y0,z0> towards
//  <0,0,observation_dist>.  Returns the distance between source_point[3]
//  and exit_point[3].
//  If postiveSolution==false, then the solution at a negative z will be given
//Note: CAN be used with a 2-dimensional integral which holds phi constant
//Note: it is safe to pass the same array as source_point[3] and exit_point[3]
double exit_point_of_sphere_z( const double source_point[3],
                               double exit_point[3],
                               double sphere_rad,
                               double observation_dist,
                               bool postiveSolution = true );

//distance(...): returns distance between two points specified in terms of x, y,
//  and z.
double distance( const double point_a[3], const double point_b[3] );


/*
 *The following functions are depreciated, but left commented out since I will
 *likely need to refer to them later - wcjohns 20121127
//exit_point_of_sphere_x(...): similar to exit_point_of_sphere_z(...), but
//  for the detector at <observation_dist,0,0>
//Makes x0,y0,z0 be the point of intersection of a sphere or radius
//  'sphere_rad' centered at the origin, with the line pointing from <x0,y0,z0>
//  towards <observation_dist,0,0>
//  Note: can not be used with a 2-dimensional integral which holds phi or theta
//        constant
void exit_point_of_sphere_x( double &x0, double &y0, double &z0,
                            double sphere_rad, double observation_dist );

double len_in_sphere_x( double r, double theta,
                        double phi, double sphere_rad, double observation_dist );
*/

struct DistributedSrcCalc
{
  //Right now this struct assumes sources are solid, in terms of the attenuation
  //  calculation
  DistributedSrcCalc();
  void eval( const double xx[], const int *ndimptr,
             double ff[], const int *ncompptr ) const;

  size_t m_sourceIndex;
  double m_detectorRadius;
  double m_observationDist;
  std::vector<std::pair<double,double> > m_sphereRadAndTransLenCoef;

  double energy;
  double integral;
  double srcVolumumetricActivity;

  const SandiaDecay::Nuclide *nuclide;
};//struct DistributedSrcCalc



class PointSourceShieldingChi2Fcn
    : public ROOT::Minuit2::FCNBase
{
//This class evaluated the chi2 of a given hypothesis, where it is assumed the
//  radioactive source is a point source located at the center of concentric
//  spherical shells consisting of non-radioactive materials, that may either
//  be defined via the 'Material' class, or be a generic material defined by
//  atomic number and areal density which is indicated by a NULL pointer.
//
//parameters:
//-Activity (in MBq) nuclide 0, (nuclides are sorted alphebaetically by name)
//-Age of nuclide 0
//-Activity (in MBq) nuclide 1
//-Age of nuclide 1 (if negative, must be negative value of one plus index of
//                   defining nuclide; e.g., a negative int.  Hacky, but whatever.
//                   To get age for this case, do `age = 1 + 2*((-1*parameter) - 1)` )
// ...
//if material 0 normal material (if Material* is non-NULL pointer)
//    -Material Thickness
//    -ignored parameter here
//else if generic material (if Material* is NULL)
//    -atomic number
//    -areal density  (in units of PhysicalUnits,
//                     e.g. to print out to user divide by g/cm2)
//if material 1 normal material (if Material* is non-NULL pointer)
//    -Material Thinckness
// ...
//
//could add another member variable that holds pointer to source isotopes to
//  fit for mass fractions of.  Then there would be M-1 additional parameters,
//  where M is the number of isotopes to fit for in a material.
//  Mass fraction is then ...
  
public:
  //Inorder to keep numbers roughly where Minuit2 can handle them, we
  //  have to work in units of 1.0E6 becquerel
  static const double sm_activityUnits;  //SandiaDecay::MBq


  /** A struct to hold information about the material shieldings are made out of, and their respective self-attenuating and trace sources.
   
   Thicknesses and activitiy levels are specified by the fitting parameters, so not tracked in this struct.
   */
  struct ShieldingInfo
  {
    /** The material this struct holds info for.
     If nullptr it represents generic shielding (e.g., AN/AD), and niether #self_atten_sources or #trace_sources may have entries.
     */
    const Material *material;
    
    /** The nuclides which act as self-attenuating sources within the material.  The Material object must contain these nuclides
     as components, as well as at least one peak being used have this nuclide as its assigned nuclide.
     */
    std::vector<const SandiaDecay::Nuclide *> self_atten_sources;
    
    /** The trace sources within this shielding.
     
     At least one peak being used have this nuclide as its assigned nuclide.
     */
    std::vector<std::pair<const SandiaDecay::Nuclide *,TraceActivityType>> trace_sources;
  };//struct ShieldingInfo
  
  
  PointSourceShieldingChi2Fcn(
                      double distance, double liveTime,
                      const std::vector<PeakDef> &peaks,
                      std::shared_ptr<const DetectorPeakResponse> detector,
                      const std::vector<ShieldingInfo> &materials,
                      bool allowMultipleNucsContribToPeaks );
  virtual ~PointSourceShieldingChi2Fcn();

  /** Causes exception to be thrown if DoEval() is called afterwards. */
  void cancelFit();
  
  /** Information tracked during test evaluation of the Chi2.  Set this object
   using #setGuiProgressUpdater so you can properly bind things to the GUI.
   */
  struct GuiProgressUpdateInfo
  {
    GuiProgressUpdateInfo( const size_t updateFreqMs,
                          std::function<void(size_t ncalls, double elapsed_time, double chi2, std::vector<double> pars)> updater );
    
    void fitting_starting();
    void completed_eval( const double chi2, const std::vector<double> &pars );
    
  private:
    /** A mutex that protects all member variables.
     Lock is taken in #fitting_starting and #completed_eval functions, and while m_gui_updater is being called.
     */
    std::mutex m_mutex;
    
    /** The function called periodically to update the best chi2 found. */
    const std::function<void(size_t, double, double, std::vector<double>)> m_gui_updater;
    
    const size_t m_update_frequency_ms;
    
    size_t m_num_fcn_calls;
    
    double m_fitStartTime;
    double m_currentTime;
    double m_lastGuiUpdateTime;
    
    double m_bestChi2;
    std::vector<double> m_bestParameters;
  };//struct GuiProgressUpdateInfo
  
  /** Call this function to have the Chi2 call the specified callback with the
      best solution found so far, at the specified intervals.  The callback will
      be posted to be executed in the specified Wt session using WServer.
   */
  void setGuiProgressUpdater( std::shared_ptr<GuiProgressUpdateInfo> updateInfo );
  
  /** Call this function *just* before starting to fit; it will set the zombie
      timer for the specified delay.  The zombie timer is how long the fit can
      go on before its considered a failure and is aborted.
      If the guiProgressUpdater() has bee set, then the best chi2 seen so far
      and, the timer, and parameters are reset.
   */
  void fittingIsStarting( const size_t deadlineMs );
  
  /** Call this function as soon as fitting is done - it cancels the zombie
      timer.
   */
  void fittingIsFinished();
  
  /** Sets whether to use a SpecUtilsAsync::ThreadPool when calculating contributions for each peak from self-attenuating and/or
   trace sources.
   If you are doing multiple parallel fits, you may want to disable multithread to better use the cpu.
   
   Default is true.
   
   \sa m_self_att_multithread
   */
  void setSelfAttMultiThread( const bool do_multithread );
  

/*
   Need to add method to extract mass fraction for isotopes fitting for the mass
   fraction
   still need to modify:
     --energy_chi_contributions(...)
     --Anywhere that references m_materials
*/
  void massFraction( double &massFrac, double &uncert,
                     const Material *material,
                     const SandiaDecay::Nuclide *nuc,
                     const std::vector<double> &pars,
                     const std::vector<double> &errors ) const;
  
  double massFraction( const Material *material,
                          const SandiaDecay::Nuclide *nuc,
                          const std::vector<double> &pars ) const;
  double massFractionUncert( const Material *material,
                             const SandiaDecay::Nuclide *nuc,
                             const std::vector<double> &pars,
                             const std::vector<double> &error ) const;
  
  bool isVariableMassFraction( const Material *material,
                               const SandiaDecay::Nuclide *nuc ) const;
  bool hasVariableMassFraction( const Material *material ) const;
  std::shared_ptr<Material> variedMassFracMaterial( const Material *material,
                                        const std::vector<double> &x ) const;
  void setNuclidesToFitMassFractionFor( const Material *material,
                   const std::vector<const SandiaDecay::Nuclide *> &nuclides );
  std::vector<const Material *> materialsFittingMassFracsFor() const;
  const std::vector<const SandiaDecay::Nuclide *> &nuclideFittingMassFracFor(
                                              const Material *material ) const;
  

  
  //setBackgroundPeaks(...): if you wish to correct for background counts, you
  //  can set that here.  The peaks you pass in should be the original
  //  background peaks for the background; similarly for the live time.  This
  //  function will scale peak areas/uncertainties for the live time
  void setBackgroundPeaks( const std::vector<PeakDef> &peaks, double liveTime );
  
  virtual double operator()( const std::vector<double> &x ) const;
  virtual double DoEval( const std::vector<double> &x ) const;

  //energy_chi_contributions(...): gives the chi2 contributions for each energy
  //  of peak, for the parameter values of x.
  //Does not take into account self attenuation.
  //'mixturecache' is to speed up multiple computations, and may be empty at
  //  first.
  typedef std::map< const SandiaDecay::Nuclide *, SandiaDecay::NuclideMixture> NucMixtureCache;
  
  //If 'info' is non-null then it will be filled with information about how much
  //  each nuclide/peak was attributed to each detected peak (currently not
  //  implemented)
  //Each retured enry is {energy,chi,scale,PeakColor,scale_uncert}.  Scale is obs/expected
  std::vector< std::tuple<double,double,double,Wt::WColor,double> > energy_chi_contributions(
                                  const std::vector<double> &x,
                                  NucMixtureCache &mixturecache,
                                  const bool allow_multiple_iso_contri,
                                  std::vector<std::string> *info = 0 ) const;

  PointSourceShieldingChi2Fcn&	operator=( const PointSourceShieldingChi2Fcn & );
  virtual double Up() const;

  size_t numExpectedFitParameters() const;

  //nuclide(...) throws std::runtime_error if an invalid number is passed in
  const SandiaDecay::Nuclide *nuclide( const int nucN ) const;

  //activity(...) and age(...) will throw runtime_exception if params is wrong
  //  size or invalid nuclide passed in
  double activity( const SandiaDecay::Nuclide *nuclide,
                   const std::vector<double> &params ) const;
  double activityUncertainty( const SandiaDecay::Nuclide *nuclide,
                  const std::vector<double> &params,
                  const std::vector<double> &errors ) const;
  double age( const SandiaDecay::Nuclide *nuclide,
                   const std::vector<double> &params ) const;
  size_t nuclideIndex( const SandiaDecay::Nuclide *nuclide ) const;

  /** Returns whether or not the nuclide is a self-attenuating source. */
  bool isSelfAttenSource( const SandiaDecay::Nuclide *nuc ) const;
  
  /** Returns whether or not the nuclide is a trace source. */
  bool isTraceSource( const SandiaDecay::Nuclide *nuc ) const;
  
  /** Returns the trace source activity type for a nuc.
   
   Throws exception if nuc is not a trace source.
   */
  GammaInteractionCalc::TraceActivityType traceSourceActivityType(
                                                          const SandiaDecay::Nuclide *nuc ) const;
  
  /** Returns whether or not the nuclide is a self-attenuating OR trace source. */
  bool isVolumetricSource( const SandiaDecay::Nuclide *nuc ) const;
  
  
  /** Returns the activity of the specified nuclide that is a self-attenuating source.
   
   Throws exception if nuclide is nullptr, or not a self-attenuating source.
   */
  double activityOfSelfAttenSource( const SandiaDecay::Nuclide *nuclide,
                                      const std::vector<double> &params ) const;
  
  /** Returns the total activity of the specified nuclide that is a trace source.
   
   Note: the #activity function returns the display activity (so either total, per cc, or per g) for the trace sources.
   
   Throws exception if nuclide is nullptr, or not a trace source.
   */
  double totalActivityOfTraceSource( const SandiaDecay::Nuclide *nuclide,
                                   const std::vector<double> &params ) const;
  
  size_t numNuclides() const;
  size_t numMaterials() const;
  
  bool isSpecificMaterial( int materialNum ) const;
  bool isGenericMaterial( int materialNum ) const;
  
  //material(): will throw exception if invalid materialNum, and will return
  //  NULL if a generic material.
  const Material *material( int materialNum ) const;
  
  //thickness(...): will throw std::runtime_exception if material is a generic
  //  material
  double thickness( int materialNum,
                    const std::vector<double> &params ) const;

  //arealDensity(...): will throw std::runtime_exception if material is a
  //  specific material
  double arealDensity( int materialNum,
                       const std::vector<double> &params ) const;

  //atomicNumber(...): will throw std::runtime_exception if material is a
  //  specific material
  double atomicNumber( int materialNum,
                       const std::vector<double> &params ) const;

  const std::vector<PeakDef> &peaks() const;

  static void selfShieldingIntegration( DistributedSrcCalc &calculator );

protected:
  
  void zombieCallback( const boost::system::error_code &ec );

  //observedPeakEnergyWidths(): get energy sorted pairs of peak means and widths
  static std::vector< std::pair<double,double> > observedPeakEnergyWidths(
                                           const std::vector<PeakDef> &peaks );

  //cluster_peak_activities(...): clusters the number of decays per second
  //  by energies.  If photopeakClusterSigma>0.0, then all gamma lines nearby
  //  a peaks mean will be considered to contribute to that peak.  'Nearby'
  //  is defined by 'energie_widths' - the energy and sigma of fit peaks.
  //  The same 'energie_widths' must be used with all subsequent calls using
  //  the same 'energy_count_map' (this is unchecked, so dont violate it).
  //  energie_widths should consist of the true photopeak energies, and not
  //  the detected energy.
  //  It is also assumed that the Nuclides have been added to the mixture
  //  with an activity of 1.0*sm_activityUnits (this will be divided out then
  //  multiplied by 'act'), at an initial age of 0.0.
  //  There must be exactly one parent nuclide in 'mixture' or an exception will
  //  be thrown.
  //  If energyToCluster is > 0.0, then only photopeaks within
  //  photopeakClusterSigma in mixture will be clustered and added to
  //  energy_count_map.  If energyToCluster <= 0.0, then all photopeaks in
  //  mixture will be clustered and added to energy_count_map.
  static void cluster_peak_activities( std::map<double,double> &energy_count_map,
                  const std::vector< std::pair<double,double> > &energie_widths,
                  SandiaDecay::NuclideMixture &mixture,
                  const double act,
                  const double thisAge,
                  const double photopeakClusterSigma,
                  const double energyToCluster,
                  std::vector<std::string> *info
              );


  //returns the chi computed from the expected verses observed counts; one
  //  chi2 for each peak energy.  Each returned entry is {energy,chi,scale,PeakColor,ScaleUncert},
  //  where scale is observed/expected
  static std::vector< std::tuple<double,double,double,Wt::WColor,double> > expected_observed_chis(
                              const std::vector<PeakDef> &peaks,
                              const std::vector<PeakDef> &backgroundPeaks,
                              const std::map<double,double> &energy_count_map,
                              std::vector<std::string> *info = 0 );


protected:
  std::atomic<int> m_cancel;  //0==KeepGoing, 1==UserCancelled, 2==CalcTimeout
  
  //Used to determine if m_guiUpdateInfo related stuff should be checked on
  std::atomic<bool> m_isFitting;
  
  std::shared_ptr<GuiProgressUpdateInfo> m_guiUpdateInfo;
  
  //Sometimes self attenuating fits can run-away and never end - protect against it
  std::mutex m_zombieCheckTimerMutex;
  std::shared_ptr<boost::asio::deadline_timer> m_zombieCheckTimer;
  
  double m_distance;
  double m_liveTime;

  //m_photopeakClusterSigma: if >=0.0 then not just the photpeak the fit peak
  //  is assigned to will be used, but also other photopeaks within
  //  m_photopeakClusterSigma of the fit peaks sigma will be used to calc
  //  expected as well.  The only place this variable is currently referenced is
  //  in energy_chi_contributions(...).  Note that if this value is large then
  //  there is the possibility of double counting the expected contributions
  //  from
  double m_photopeakClusterSigma;
  std::vector<PeakDef> m_peaks;
  std::vector<PeakDef> m_backgroundPeaks;
  std::shared_ptr<const DetectorPeakResponse> m_detector;
  std::vector<ShieldingInfo> m_materials;
  std::vector<const SandiaDecay::Nuclide *> m_nuclides; //sorted alphebetically and unique
  
  typedef std::map<const Material *,std::vector<const SandiaDecay::Nuclide *> > MaterialToNucsMap;
  //Nuclides will stored sorted alphebetically
  MaterialToNucsMap m_nuclidesToFitMassFractionFor;
  
  bool m_allowMultipleNucsContribToPeaks;
  
  /** Wether to use SpecUtilsAsync::ThreadPool to calculate self-attenuation peak values.
   
   TODO: Figure out if this should be an std::atomic
   
   Default is true.
   
   \sa setSelfAttMultiThread
   */
  bool m_self_att_multithread;
  
  
  //A cache of nuclide mixtures to
  mutable NucMixtureCache m_mixtureCache;
  static const size_t sm_maxMixtureCacheSize = 10000;
};//class PointSourceShieldingChi2Fcn

}//namespace GammaInteractionCalc

#endif  //GammaInteractionCalc_h





