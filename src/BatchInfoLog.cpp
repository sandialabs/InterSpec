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
#include <fstream>
#include <iostream>
#include <exception>


#include "external_libs/SpecUtils/3rdparty/inja/inja.hpp"

#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/D3SpectrumExport.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/BatchPeak.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/BatchInfoLog.h"
#include "InterSpec/BatchActivity.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShowRiidInstrumentsAna.h"

using namespace std;

namespace BatchInfoLog
{
  
  // Adds the basic direct info on a source (nuclide name, activity, age, etc), but does not
  //  Add which peaks it contributes to, or any information on gammas
void add_basic_src_details( const GammaInteractionCalc::SourceDetails &src,
                          const std::shared_ptr<const DetectorPeakResponse> &drf,
                          const bool useBq,
                          const std::vector<GammaInteractionCalc::ShieldingDetails> *shield_details,
                          nlohmann::basic_json<> &src_json )
{
  const DetectorPeakResponse::EffGeometryType eff_type = drf ? drf->geometryType()
                                              : DetectorPeakResponse::EffGeometryType::FarField;
  
  const string act_postfix = DetectorPeakResponse::det_eff_geom_type_postfix(eff_type);
  
  
  src_json["Nuclide"] = src.nuclide->symbol;
  src_json["Activity"] = PhysicalUnits::printToBestActivityUnits(src.activity,4,!useBq) + act_postfix;
  src_json["Activity_bq"] = src.activity / PhysicalUnits::bq;
  src_json["Activity_kBq"] = src.activity / PhysicalUnits::kBq;
  src_json["Activity_MBq"] = src.activity / PhysicalUnits::MBq;
  src_json["Activity_GBq"] = src.activity / PhysicalUnits::GBq;
  src_json["Activity_ci"] = src.activity / PhysicalUnits::ci;
  src_json["Activity_mCi"] = src.activity / PhysicalUnits::mCi;
  src_json["Activity_uCi"] = src.activity / PhysicalUnits::microCi;
  src_json["Activity_pCi"] = src.activity / PhysicalUnits::pCi;
  src_json["ActivityPostFix"] = act_postfix;
  
  bool activityIsFit = src.activityIsFit;
  if( src.isSelfAttenSource )
  {
    activityIsFit |= src.isSelfAttenVariableMassFrac;
   
    assert( shield_details );
    if( shield_details )
    {
      // TODO: check if mass fraction is being fit, or if a shielding dimension is being fit
      assert( src.selfAttenShieldIndex < shield_details->size() );
      const GammaInteractionCalc::ShieldingDetails &shield = shield_details->at(src.selfAttenShieldIndex);
      assert( shield.m_name == src.selfAttenShieldName );
      for( unsigned int i = 0; i < shield.m_num_dimensions; ++i )
        activityIsFit |= shield.m_fit_dimension[i];
    }//if( shield_details )
  }//if( src.isSelfAttenSource )
  
  src_json["ActivityIsFit"] = activityIsFit;
  
  
  if( activityIsFit )
  {
    src_json["ActivityUncert"] = PhysicalUnits::printToBestActivityUnits(src.activityUncertainty,4,!useBq) + act_postfix;
    
    src_json["ActivityUncert_bq"]  = src.activityUncertainty / PhysicalUnits::bq;
    src_json["ActivityUncert_kBq"] = src.activityUncertainty / PhysicalUnits::kBq;
    src_json["ActivityUncert_MBq"] = src.activityUncertainty / PhysicalUnits::MBq;
    src_json["ActivityUncert_GBq"] = src.activityUncertainty / PhysicalUnits::GBq;
    src_json["ActivityUncert_ci"]  = src.activityUncertainty / PhysicalUnits::ci;
    src_json["ActivityUncert_mCi"] = src.activityUncertainty / PhysicalUnits::mCi;
    src_json["ActivityUncert_uCi"] = src.activityUncertainty / PhysicalUnits::microCi;
    src_json["ActivityUncert_pCi"] = src.activityUncertainty / PhysicalUnits::pCi;
    
    const double act_uncert_percent = 100.0 * src.activityUncertainty / src.activity;
    src_json["ActivityUncertPercent"] = SpecUtils::printCompact( act_uncert_percent, 4 );
  }else
  {
    assert( src.activityUncertainty <= 0.0 );
  }
  
  src_json["NuclideMass"] = PhysicalUnits::printToBestMassUnits(src.nuclideMass, 4);
  src_json["Age"] = PhysicalUnits::printToBestTimeUnits(src.age, 4);
  src_json["AgeSeconds"] = src.age / PhysicalUnits::second;
  src_json["AgeDays"] = src.age / PhysicalUnits::day;
  src_json["AgeYears"] = src.age / PhysicalUnits::year;
  src_json["AgeIsFittable"] = src.ageIsFittable;
  src_json["AgeIsFit"] = src.ageIsFit;
  
  if( src.ageIsFit )
  {
    src_json["AgeUncert"] = PhysicalUnits::printToBestTimeUnits(src.ageUncertainty, 4);
    src_json["AgeUncertSeconds"] = src.ageUncertainty / PhysicalUnits::second;
    src_json["AgeUncertDays"] = src.ageUncertainty / PhysicalUnits::day;
    src_json["AgeUncertYears"] = src.ageUncertainty / PhysicalUnits::year;
  }else
  {
    assert( src.ageUncertainty <= 0.0 );
  }
  
  if( src.ageDefiningNuc )
    src_json["AgeDefiningNuclide"] = src.ageDefiningNuc->symbol;
  
  src_json["IsTraceSource"] = src.isTraceSource;
  
  if( src.isTraceSource )
  {
    string trace_src_postfix = "";
    
    switch( src.traceActivityType )
    {
      case GammaInteractionCalc::TraceActivityType::TotalActivity:
        trace_src_postfix = "";
        break;
        
      case GammaInteractionCalc::TraceActivityType::ActivityPerCm3:
        trace_src_postfix = "/cm^3";
        break;
      
      case GammaInteractionCalc::TraceActivityType::ExponentialDistribution:
        trace_src_postfix = "/m2 exp";
        break;
        
      case GammaInteractionCalc::TraceActivityType::ActivityPerGram:
        trace_src_postfix = "/g";
        break;
        
      case GammaInteractionCalc::TraceActivityType::NumTraceActivityType:
        break;
    }//switch( type )
    
    src_json["TraceActivityType"] = GammaInteractionCalc::to_str(src.traceActivityType);
    src_json["TraceDisplayActivity"] = PhysicalUnits::printToBestActivityUnits(src.traceSrcDisplayAct,4,!useBq) + trace_src_postfix;
    src_json["TraceDisplayActivity_bq"] = src.traceSrcDisplayAct / PhysicalUnits::bq;
    src_json["TraceDisplayActivity_kBq"] = src.traceSrcDisplayAct / PhysicalUnits::kBq;
    src_json["TraceDisplayActivity_MBq"] = src.traceSrcDisplayAct / PhysicalUnits::MBq;
    src_json["TraceDisplayActivity_GBq"] = src.traceSrcDisplayAct / PhysicalUnits::GBq;
    src_json["TraceDisplayActivity_ci"]  = src.traceSrcDisplayAct / PhysicalUnits::ci;
    src_json["TraceDisplayActivity_mCi"] = src.traceSrcDisplayAct / PhysicalUnits::mCi;
    src_json["TraceDisplayActivity_uCi"] = src.traceSrcDisplayAct / PhysicalUnits::microCi;
    src_json["TraceDisplayActivity_pCi"] = src.traceSrcDisplayAct / PhysicalUnits::pCi;
    
    src_json["TraceActivityPostFix"] = trace_src_postfix;
    
    if( src.activityIsFit )
    {
      src_json["TraceDisplayActivityUncert"] = PhysicalUnits::printToBestActivityUnits(src.traceSrcDisplayActUncertainty,4,!useBq) + trace_src_postfix;
      src_json["TraceDisplayActivityUncert_bq"] = src.traceSrcDisplayActUncertainty / PhysicalUnits::bq;
      src_json["TraceDisplayActivityUncert_kBq"] = src.traceSrcDisplayActUncertainty / PhysicalUnits::kBq;
      src_json["TraceDisplayActivityUncert_MBq"] = src.traceSrcDisplayActUncertainty / PhysicalUnits::MBq;
      src_json["TraceDisplayActivityUncert_GBq"] = src.traceSrcDisplayActUncertainty / PhysicalUnits::GBq;
      src_json["TraceDisplayActivityUncert_ci"]  = src.traceSrcDisplayActUncertainty / PhysicalUnits::ci;
      src_json["TraceDisplayActivityUncert_mCi"] = src.traceSrcDisplayActUncertainty / PhysicalUnits::mCi;
      src_json["TraceDisplayActivityUncert_uCi"] = src.traceSrcDisplayActUncertainty / PhysicalUnits::microCi;
      src_json["TraceDisplayActivityUncert_pCi"] = src.traceSrcDisplayActUncertainty / PhysicalUnits::pCi;
      
    }else
    {
      assert( src.traceSrcDisplayActUncertainty <= 0.0 );
    }
    
    if( src.traceActivityType == GammaInteractionCalc::TraceActivityType::ExponentialDistribution )
    {
      src_json["TraceRelaxationLength"] = PhysicalUnits::printToBestLengthUnits(src.traceRelaxationLength, 4);
      src_json["TraceRelaxationLength_mm"] = src.traceRelaxationLength / PhysicalUnits::mm;
      src_json["TraceRelaxationLength_cm"] = src.traceRelaxationLength / PhysicalUnits::cm;
      src_json["TraceRelaxationLength_m"] = src.traceRelaxationLength / PhysicalUnits::m;
      src_json["TraceRelaxationLength_inch"] = src.traceRelaxationLength / (2.54*PhysicalUnits::cm);
    }
  }//if( src.isTraceSource )
  
  src_json["IsSelfAttenSource"] = src.isSelfAttenSource;
  if( src.isSelfAttenSource )
  {
    src_json["SelfAttenShieldIndex"] = static_cast<int>( src.selfAttenShieldIndex );
    src_json["SelfAttenShieldName"] = src.selfAttenShieldName;
    src_json["SelfAttenIsVariableMassFrac"] = src.isSelfAttenVariableMassFrac;
    src_json["SelfAttenMassFrac"] = src.selfAttenMassFrac;
    if( src.isSelfAttenVariableMassFrac )
      src_json["SelfAttenMassFracUncert"] = src.selfAttenMassFracUncertainty;
  }
}//add_basic_src_details( SourceDetails &src, src_json )
  
  
  
  void add_act_shield_fit_options_to_json( const ShieldingSourceFitCalc::ShieldingSourceFitOptions &options,
                               const double distance,
                               const GammaInteractionCalc::GeometryType geometry,
                               const std::shared_ptr<const DetectorPeakResponse> &drf,
                               nlohmann::json &data )
  {
    const bool fixedGeom = (drf && drf->isFixedGeometry());
    
    data["FixedGeometryDetector"] = fixedGeom;
    
    auto &fit_setup = data["ActShieldFitSetup"];
    if( !fixedGeom )
    {
      fit_setup["Distance"] = PhysicalUnits::printToBestLengthUnits( distance, 3 );
      fit_setup["Distance_mm"] = distance / PhysicalUnits::mm;
      fit_setup["Distance_cm"] = distance / PhysicalUnits::cm;
      fit_setup["Distance_m"] = distance / PhysicalUnits::meter;
      fit_setup["Distance_km"] = distance / (1000.0*PhysicalUnits::meter);
      fit_setup["Distance_inch"] = distance / (2.54*PhysicalUnits::cm);
      fit_setup["Distance_feet"] = distance / (12.0*2.54*PhysicalUnits::cm);
      fit_setup["Geometry"] = GammaInteractionCalc::to_str( geometry );
    }else
    {
      string gem_desc;
      switch( drf->geometryType() )
      {
        case DetectorPeakResponse::EffGeometryType::FarField:
          assert( drf->geometryType() != DetectorPeakResponse::EffGeometryType::FarField );
          break;
          
        case DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct:
          gem_desc = "total activity";
          break;
          
        case DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2:
          gem_desc = "activity per square centimeter";
          break;
          
        case DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2:
          gem_desc = "activity per square meter";
          break;
          
        case DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram:
          gem_desc = "activity per gram";
          break;
      }//switch( drf->geometryType() )
      
      fit_setup["FixedGeometryType"] = gem_desc;
    }//if( !fixedGeom ) / else
    
    auto &fit_options = fit_setup["FitOptions"];
    fit_options["InterferenceCorrection"]   = options.multiple_nucs_contribute_to_peaks;
    fit_options["AttenuateForAir"]          = (options.attenuate_for_air && (!drf || !drf->isFixedGeometry()));
    fit_options["DecayDuringMeasurement"]   = options.account_for_decay_during_meas;
    fit_options["MultithreadSelfAttenCalc"] = options.multithread_self_atten;
    fit_options["PhotopeakClusterSigma"]    = options.photopeak_cluster_sigma;
    fit_options["BackgroundPeakSubtract"]   = options.background_peak_subtract;
    fit_options["ElementNuclidesSameAge"]   = options.same_age_isotopes;
  }//void add_act_shield_fit_options_to_json(...)
  
  
  /** Adds basic information about a peak (energy, fwhm, counts, etc), but not any information
     about gammas that contribute to it, etc
   */
  void add_basic_peak_info( const GammaInteractionCalc::PeakDetail &peak, nlohmann::basic_json<> &peak_json )
  {
    char buffer[64] = { '\0' };
    snprintf( buffer, sizeof(buffer), "%.2f", peak.energy );
    peak_json["Energy"] = string(buffer);
    peak_json["Energy_keV"] = peak.energy;
    
    snprintf( buffer, sizeof(buffer), "%.3f", peak.decayParticleEnergy );
    peak_json["DecayParticleEnergy"] = string(buffer);
    peak_json["DecayParticleEnergy_keV"] = peak.decayParticleEnergy;
    peak_json["AssignedNuclide"] = peak.assignedNuclide;
    peak_json["FWHM"] = peak.fwhm;
    
    
    peak_json["Counts"] = peak.counts;
    peak_json["CountsStr"] = SpecUtils::printCompact(peak.counts, 4);
    peak_json["CountsUncert"] = peak.countsUncert;
    peak_json["CountsUncertStr"] = SpecUtils::printCompact(peak.countsUncert, 4);
    peak_json["Cps"] = peak.cps;
    peak_json["CpsStr"] = SpecUtils::printCompact(peak.cps,4);
    peak_json["CpsUncert"] = peak.cpsUncert;
    peak_json["ShieldAttenuations"] = peak.m_attenuations;
    
    if( peak.backgroundCounts > 0.0 )
    {
      peak_json["BackgroundCounts"] = peak.backgroundCounts;
      peak_json["BackgroundCountsStr"] = SpecUtils::printCompact(peak.backgroundCounts, 4);
      
      if( peak.backgroundCountsUncert > 0.0 )
      {
        peak_json["BackgroundCountsUncert"] = peak.backgroundCountsUncert;
        peak_json["BackgroundCountsUncertStr"] = SpecUtils::printCompact(peak.backgroundCountsUncert, 4);
      }
    }//if( peak.backgroundCounts > 0.0 )
    
    peak_json["SignalCounts"] = peak.observedCounts;
    peak_json["SignalCountsStr"] = SpecUtils::printCompact(peak.observedCounts, 4);
    peak_json["SignalCountsUncert"] = peak.observedUncert;
    peak_json["SignalCountsUncertStr"] = SpecUtils::printCompact(peak.observedUncert, 4);
    
    peak_json["PredictedCounts"] = peak.expectedCounts;
    peak_json["PredictedNumSigmaOff"] = peak.numSigmaOff;
    
    peak_json["ObservedOverPredicted"] = peak.observedOverExpected;
    peak_json["ObservedOverPredictedUncert"] = peak.observedOverExpectedUncert;
    
    peak_json["DetectorSolidAngleFraction"] = peak.detSolidAngle;
    peak_json["DetectorIntrinsicEff"] = peak.detIntrinsicEff;
    peak_json["DetectorEff"] = peak.detEff;
  }//add_basic_peak_info( const PeakDetail &peak, nlohmann::basic_json<> &peak_json )
   
  
  void add_gamma_info_for_peak( const GammaInteractionCalc::PeakDetailSrc &ps,
                    const GammaInteractionCalc::SourceDetails * const src,
                    const std::shared_ptr<const DetectorPeakResponse> &drf,
                    const bool useBq,
                    const std::vector<GammaInteractionCalc::ShieldingDetails> * const shield_details,
                    nlohmann::basic_json<> &gamma_json )
  {
    assert( ps.nuclide );
    
    char buffer[64] = { '\0' };
    snprintf( buffer, sizeof(buffer), "%.3f", ps.energy );
    
    gamma_json["Nuclide"] = ps.nuclide ? ps.nuclide->symbol : string("null");
    gamma_json["Energy"] = string(buffer);
    gamma_json["Energy_keV"] = ps.energy;
    
    gamma_json["BranchingRatio"] = ps.br;
    gamma_json["BranchingRatioStr"] = SpecUtils::printCompact( ps.br, 5 );
    
    gamma_json["CpsContributedToPeak"] = ps.cps;
    gamma_json["CpsContributedToPeakStr"] = SpecUtils::printCompact( ps.cps, 5 );
    
    gamma_json["CountsContributedToPeak"] = ps.counts;
    gamma_json["CountsContributedToPeakStr"] = SpecUtils::printCompact( ps.counts, 5 );
    
    gamma_json["HasDecayCorrection"] = (ps.decayCorrection > 0.0);
    if( ps.decayCorrection > 0.0 )
    {
      gamma_json["DecayCorrection"] = ps.decayCorrection;
      gamma_json["DecayCorrectionStr"] = SpecUtils::printCompact(ps.decayCorrection, 4);
    }
    
    // The other information in `PeakDetailSrc` should be repeat of
    //  `GammaInteractionCalc::SourceDetails`.  We _could_ access all this information
    //  from the templating code, since we are in a loop over `SourceDetails`, but
    //  to make things easier on people, we'll just re-include it here.
    if( src )
      add_basic_src_details( *src, drf, useBq, shield_details, gamma_json );
  };//add_gamma_info_for_peak(...)
  
  
  void shield_src_fit_results_to_json( const ShieldingSourceFitCalc::ModelFitResults &results,
                                      const std::shared_ptr<const DetectorPeakResponse> &drf,
                                      const bool useBq,
                                      nlohmann::json &data )
  {
    const DetectorPeakResponse::EffGeometryType eff_type = drf ? drf->geometryType()
                                                : DetectorPeakResponse::EffGeometryType::FarField;
    
    const string act_postfix = DetectorPeakResponse::det_eff_geom_type_postfix(eff_type);
    
    add_act_shield_fit_options_to_json( results.options, results.distance, results.geometry, drf, data );
    
    // Now put in results
    data["FitChi2"] = results.chi2;
    data["EstimatedDistanceToMinimum"] = results.edm;
    data["NumberFcnCalls"] = results.num_fcn_calls;
    data["NumDof"] = results.numDOF;
    
    auto &fit_pars = data["RawFitParameter"];
    fit_pars["Values"] = results.paramValues;
    fit_pars["Errors"] = results.paramErrors;
    
    if( !results.errormsgs.empty() )
      data["ErrorMessages"] = results.errormsgs;
    
    int num_sources = 0;
    bool hasAnyTraceSrc = false, hasAnyVolumetricSrc = false, hasFitAnyAge = false;
    if( results.source_calc_details )
    {
      for( const GammaInteractionCalc::SourceDetails &src : *results.source_calc_details )
      {
        data["Sources"].push_back( {} );
        auto &src_json = data["Sources"].back();
       
        ++num_sources;
        hasAnyTraceSrc |= src.isTraceSource;
        hasAnyVolumetricSrc |= src.isSelfAttenSource;
        hasFitAnyAge |= src.ageIsFit;
        
        add_basic_src_details( src, drf, useBq, results.shield_calc_details.get(), src_json );
        
        // We wont put peaks into this struct, but instead when we make the JSON, we'll
        //  insert peaks from `PeakDetail` as `PeakDetailSrc::nuclide` match this nuclide.
        if( results.peak_calc_details )
        {
          for( const GammaInteractionCalc::PeakDetail &peak : *results.peak_calc_details )
          {
            bool src_contribs_to_peak = false;
            for( const GammaInteractionCalc::PeakDetailSrc &ps :  peak.m_sources )
            {
              src_contribs_to_peak = (ps.nuclide == src.nuclide);
              if( src_contribs_to_peak )
                break;
            }
            if( !src_contribs_to_peak )
              continue;
            
            src_json["PeaksThisNucContributesTo"].push_back( {} );
            nlohmann::basic_json<> &peak_json = src_json["PeaksThisNucContributesTo"].back();
            
            add_basic_peak_info( peak, peak_json );
            
            for( const GammaInteractionCalc::PeakDetailSrc &ps :  peak.m_sources )
            {
              if( ps.nuclide != src.nuclide )
                continue;
              
              peak_json["ThisNucsGammasForPeak"].push_back( {} );
              nlohmann::basic_json<> &gamma_json = peak_json["ThisNucsGammasForPeak"].back();
              
              add_gamma_info_for_peak( ps, &src, drf, useBq, results.shield_calc_details.get(), gamma_json );
            }
          }//for( loop over results.peak_calc_details )
        }//if( results.peak_calc_details )
      }//for( loop over results.source_calc_details )
    }//if( results.source_calc_details )
    
    data["NumSources"] = num_sources;
    
    data["HasTraceSource"] = hasAnyTraceSrc;
    data["HasSelfAttenSource"] = hasAnyVolumetricSrc;
    data["HasVolumetricSource"] = (hasAnyTraceSrc || hasAnyVolumetricSrc);
    data["AnySourceAgeFit"] = hasFitAnyAge;
    
    auto &shieldings_json = data["Shieldings"];
    shieldings_json["Geometry"] = GammaInteractionCalc::to_str(results.geometry);
    
    bool fitAnyShielding = false;
    if( !results.shield_calc_details )
    {
      shieldings_json["NumberShieldings"] = 0;
    }else
    {
      const vector<GammaInteractionCalc::ShieldingDetails> &shield_details
                                                            = *results.shield_calc_details;
      shieldings_json["NumberShieldings"] = static_cast<int>( shield_details.size() );
      switch( results.geometry )
      {
        case GammaInteractionCalc::GeometryType::Spherical:
          shieldings_json["DimensionMeanings"] = vector<string>{"Radius"};
          shieldings_json["NumDimensions"] = 1;
          break;
          
        case GammaInteractionCalc::GeometryType::CylinderEndOn:
        case GammaInteractionCalc::GeometryType::CylinderSideOn:
          shieldings_json["DimensionMeanings"] = vector<string>{"Radius", "Length"};
          shieldings_json["NumDimensions"] = 2;
          break;
          
        case GammaInteractionCalc::GeometryType::Rectangular:
          shieldings_json["DimensionMeanings"] = vector<string>{"Width", "Height", "Depth"};
          shieldings_json["NumDimensions"] = 3;
          break;
          
        case GammaInteractionCalc::GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }//switch( m_geometry )
      
      for( size_t shield_num = 0; shield_num < shield_details.size(); ++shield_num )
      {
        const GammaInteractionCalc::ShieldingDetails &shield = shield_details[shield_num];
        
        shieldings_json["Shields"].push_back({});
        auto &shield_json = shieldings_json["Shields"].back();
        
        shield_json["Name"] = shield.m_name;
        shield_json["ShieldingNumber"] = static_cast<int>( shield_num );
        shield_json["IsGeneric"] = shield.m_is_generic;
        if( !shield.m_is_generic )
        {
          shield_json["Formula"] = shield.m_chemical_formula;
          const double density = shield.m_density * PhysicalUnits::cm3 / PhysicalUnits::g;
          shield_json["Density_gPerCm3"] = density;
          shield_json["DensityStr"] = SpecUtils::printCompact(density, 5) + " g/cm3";
        }
        
        if( shield.m_an > 0 )
          shield_json["AN"] = SpecUtils::printCompact(shield.m_an, 4);
        if( shield.m_ad > 0 )
          shield_json["AD"] = SpecUtils::printCompact(shield.m_ad*PhysicalUnits::cm2/PhysicalUnits::g, 4);
        
        if( shield.m_is_generic )
        {
          const bool fitAn = shield.m_fit_dimension[0];
          const bool fitAD = shield.m_fit_dimension[1];
          shield_json["FitAN"] = fitAn;
          shield_json["FitAD"] = fitAD;
          fitAnyShielding |= fitAn;
          fitAnyShielding |= fitAD;
        }
        
        const vector<bool> fit_dim( shield.m_fit_dimension, shield.m_fit_dimension + shield.m_num_dimensions );
        shield_json["DimensionIsFit"] = fit_dim;
        
        shield_json["NumDimensions"] = static_cast<int>( shield.m_num_dimensions );
        shield_json["Geometry"] = GammaInteractionCalc::to_str(shield.m_geometry);
        
        shield_json["Thickness"] = PhysicalUnits::printToBestLengthUnits(shield.m_thickness,3);
        shield_json["Thickness_mm"] = shield.m_thickness / PhysicalUnits::mm;
        shield_json["Thickness_cm"] = shield.m_thickness / PhysicalUnits::cm;
        shield_json["Thickness_m"] = shield.m_thickness / PhysicalUnits::meter;
        shield_json["Thickness_inch"] = shield.m_thickness / (2.54*PhysicalUnits::cm);
        shield_json["Thickness_feet"] = shield.m_thickness / (12.0*2.54*PhysicalUnits::cm);
        shield_json["VolumeCm3"] = shield.m_volume / PhysicalUnits::cm3;
        shield_json["VolumeUncertCm3"] = shield.m_volume_uncert / PhysicalUnits::cm3;
        
        shield_json["InnerRadius"] = PhysicalUnits::printToBestLengthUnits(shield.m_inner_rad, 3);
        shield_json["OuterRadius"] = PhysicalUnits::printToBestLengthUnits(shield.m_inner_rad + shield.m_thickness, 3);
        
        const vector<double> inner_dims{ shield.m_inner_dimensions, shield.m_inner_dimensions + shield.m_num_dimensions };
        const vector<double> outer_dims{ shield.m_outer_dimensions, shield.m_outer_dimensions + shield.m_num_dimensions };
        vector<double> thicknesses( shield.m_num_dimensions );
        const vector<double> dim_uncerts( shield.m_dimension_uncert, shield.m_dimension_uncert + shield.m_num_dimensions );
        
        vector<double> inner_dims_mm = inner_dims;
        vector<double> outer_dims_mm = outer_dims;
        vector<double> thicknesses_mm = thicknesses;
        vector<double> dim_uncerts_mm = dim_uncerts;
        
        vector<double> inner_dims_cm = inner_dims, outer_dims_cm = outer_dims;
        vector<double> thicknesses_cm = thicknesses, dim_uncerts_cm = dim_uncerts;
        
        vector<double> inner_dims_m = inner_dims, outer_dims_m = outer_dims;
        vector<double> thicknesses_m = thicknesses, dim_uncerts_m = dim_uncerts;
        
        vector<double> inner_dims_inch = inner_dims, outer_dims_inch = outer_dims;
        vector<double> thicknesses_inch = thicknesses, dim_uncerts_inch = dim_uncerts;
        
        
        vector<string> inner_dims_strs( shield.m_num_dimensions );
        vector<string> outer_dims_strs( shield.m_num_dimensions );
        vector<string> thicknesses_strs( shield.m_num_dimensions );
        vector<string> dim_uncerts_strs( shield.m_num_dimensions );
        
        for( unsigned int dim = 0; dim < shield.m_num_dimensions; ++dim )
        {
          fitAnyShielding |= fit_dim[dim];
          thicknesses[dim] = outer_dims[dim] - inner_dims[dim];
          thicknesses_strs[dim] = PhysicalUnits::printToBestLengthUnits( thicknesses[dim], 5 );
          inner_dims_strs[dim] = PhysicalUnits::printToBestLengthUnits( inner_dims[dim], 5 );
          outer_dims_strs[dim] = PhysicalUnits::printToBestLengthUnits( outer_dims[dim], 5 );
          if( fit_dim[dim] )
            dim_uncerts_strs[dim] = PhysicalUnits::printToBestLengthUnits( dim_uncerts[dim], 5 );
          
          inner_dims_mm[dim] /= PhysicalUnits::mm;
          outer_dims_mm[dim] /= PhysicalUnits::mm;
          thicknesses_mm[dim] /= PhysicalUnits::mm;
          dim_uncerts_mm[dim] /= PhysicalUnits::mm;
          
          inner_dims_cm[dim] /= PhysicalUnits::cm;
          outer_dims_cm[dim] /= PhysicalUnits::cm;
          thicknesses_cm[dim] /= PhysicalUnits::cm;
          dim_uncerts_cm[dim] /= PhysicalUnits::cm;
          
          inner_dims_m[dim] /= PhysicalUnits::m;
          outer_dims_m[dim] /= PhysicalUnits::m;
          thicknesses_m[dim] /= PhysicalUnits::m;
          dim_uncerts_m[dim] /= PhysicalUnits::m;
          
          inner_dims_inch[dim] /= (2.54*PhysicalUnits::cm);
          outer_dims_inch[dim] /= (2.54*PhysicalUnits::cm);
          thicknesses_inch[dim] /= (2.54*PhysicalUnits::cm);
          dim_uncerts_inch[dim] /= (2.54*PhysicalUnits::cm);
        }//for( unsigned int dim = 0; dim < shield.m_num_dimensions; ++dim )
        
        shield_json["ThicknessesUncerts"] = dim_uncerts_strs;
        shield_json["InnerDims"]          = inner_dims_strs;
        shield_json["OuterDims"]          = outer_dims_strs;
        shield_json["Thicknesses"]        = thicknesses_strs;
        
        shield_json["InnerDims_mm"] = inner_dims_mm;
        shield_json["OuterDims_mm"] = outer_dims_mm;
        shield_json["Thicknesses_mm"] = thicknesses_mm;
        shield_json["ThicknessesUncerts_mm"] = dim_uncerts_mm;
        
        shield_json["InnerDims_cm"] = inner_dims_cm;
        shield_json["OuterDims_cm"] = outer_dims_cm;
        shield_json["Thicknesses_cm"] = thicknesses_cm;
        shield_json["ThicknessesUncerts_cm"] = dim_uncerts_cm;
        
        shield_json["InnerDims_m"] = inner_dims_m;
        shield_json["OuterDims_m"] = outer_dims_m;
        shield_json["Thicknesses_m"] = thicknesses_m;
        shield_json["ThicknessesUncerts_m"] = dim_uncerts_m;
        
        shield_json["InnerDims_inch"] = inner_dims_inch;
        shield_json["OuterDims_inch"] = outer_dims_inch;
        shield_json["Thicknesses_inch"] = thicknesses_inch;
        shield_json["ThicknessesUncerts_inch"] = dim_uncerts_inch;
        
        bool fittingAnyMassFrac = false;
        for( const GammaInteractionCalc::ShieldingDetails::SelfAttenComponent &comp : shield.m_mass_fractions )
        {
          fittingAnyMassFrac |= comp.m_is_fit;
          
          shield_json["SelfAttenSources"].push_back( {} );
          auto &self_atten_json = shield_json["SelfAttenSources"].back();
        
          assert( comp.m_nuclide );
          self_atten_json["Nuclide"] = comp.m_nuclide ? comp.m_nuclide->symbol : string("null");
          self_atten_json["IsFittingMassFraction"] = comp.m_is_fit;
          self_atten_json["MassFraction"] = comp.m_mass_frac;
          self_atten_json["MassFractionStr"] = SpecUtils::printCompact(comp.m_mass_frac, 5);
          if( comp.m_is_fit )
          {
            self_atten_json["MassFractionUncert"] = comp.m_mass_frac_uncert;
            self_atten_json["MassFractionUncertStr"] = SpecUtils::printCompact(comp.m_mass_frac_uncert, 5);
          }
          
          assert( results.source_calc_details );
          if( results.source_calc_details )
          {
            const vector<GammaInteractionCalc::SourceDetails> &srcs = *results.source_calc_details;
            auto pos = std::find_if( begin(srcs), end(srcs),
              [&comp]( const GammaInteractionCalc::SourceDetails &src ){
                return (src.nuclide == comp.m_nuclide);
            } );
            assert( pos != end(srcs) );
            if( pos != end(srcs) )
              add_basic_src_details( *pos, drf, useBq, results.shield_calc_details.get(), self_atten_json );
          }//if( results.source_calc_details )
        }//for( SelfAttenComponent & comp: shield.m_mass_fractions )
        
        shield_json["FitAnyMassFraction"] = fittingAnyMassFrac;
        
        
        for( const GammaInteractionCalc::ShieldingDetails::TraceSrcDetail &trace : shield.m_trace_sources )
        {
          shield_json["TraceSources"].push_back( {} );
          auto &trace_src_json = shield_json["TraceSources"].back();
          
          assert( trace.m_nuclide );
          trace_src_json["Nuclide"] = trace.m_nuclide ? trace.m_nuclide->symbol : string("null");
          trace_src_json["TraceSourceType"] = GammaInteractionCalc::to_str( trace.m_trace_type );
          if( trace.m_is_exp_dist )
          {
            trace_src_json["RelaxationLength"] = PhysicalUnits::printToBestLengthUnits( trace.m_relaxation_length, 4 );
            trace_src_json["RelaxationLength_mm"] = trace.m_relaxation_length / PhysicalUnits::mm;
            trace_src_json["RelaxationLength_cm"] = trace.m_relaxation_length / PhysicalUnits::cm;
            trace_src_json["RelaxationLength_m"] = trace.m_relaxation_length / PhysicalUnits::m;
            trace_src_json["RelaxationLength_inch"] = trace.m_relaxation_length / (2.54*PhysicalUnits::cm);
          }
          
          assert( results.source_calc_details );
          if( results.source_calc_details )
          {
            const vector<GammaInteractionCalc::SourceDetails> &srcs = *results.source_calc_details;
            auto pos = std::find_if( begin(srcs), end(srcs),
              [&trace]( const GammaInteractionCalc::SourceDetails &src ){
                return (src.nuclide == trace.m_nuclide);
            } );
            assert( pos != end(srcs) );
            if( pos != end(srcs) )
              add_basic_src_details( *pos, drf, useBq, results.shield_calc_details.get(), trace_src_json );
          }//if( results.source_calc_details )
        }//for( const TraceSrcDetail &src : shield.m_trace_sources )
        
      }//for( const GammaInteractionCalc::ShieldingDetails &shield : *results.shield_calc_details )
    }//if( results.shield_calc_details )
    
    data["AnyShieldingFit"] = fitAnyShielding;
    
    
    // Add Peak Details to JSON
    if( results.peak_calc_details )
    {
      auto &peaks_json = data["PeaksUsedForActivityFitting"];
      
      for( const GammaInteractionCalc::PeakDetail &peak : *results.peak_calc_details )
      {
        peaks_json["Peaks"].push_back( {} );
        nlohmann::basic_json<> &peak_json = peaks_json["Peaks"].back();
        add_basic_peak_info( peak, peak_json );
        
        for( const GammaInteractionCalc::PeakDetailSrc &pksrc : peak.m_sources )
        {
          assert( pksrc.nuclide && results.source_calc_details );
          if( !pksrc.nuclide || !results.source_calc_details )
            continue;
          
          const vector<GammaInteractionCalc::SourceDetails> &srcs = *results.source_calc_details;
          auto pos = std::find_if( begin(srcs), end(srcs),
            [&peak]( const GammaInteractionCalc::SourceDetails &src ){
              return src.nuclide && (src.nuclide->symbol == peak.assignedNuclide);
          } );
          assert( pos != end(srcs) );
          if( pos == end(srcs) )
            continue;
          
          
          peak_json["Sources"].push_back( {} );
          nlohmann::basic_json<> &src_json = peak_json["Sources"].back();
          
          add_basic_src_details( *pos, drf, useBq, results.shield_calc_details.get(), src_json );
          
          src_json["HasDecayCorrection"] = (pksrc.decayCorrection > 0.0);
          if( pksrc.decayCorrection > 0.0 )
          {
            peak_json["DecayCorrection"] = pksrc.decayCorrection;
            peak_json["DecayCorrectionStr"] = SpecUtils::printCompact(pksrc.decayCorrection, 4);
          }
          
          char buffer[64] = { '\0' };
          snprintf( buffer, sizeof(buffer), "%.2f", pksrc.energy );
          
          src_json["Energy"] = string(buffer);
          src_json["Energy_keV"] = pksrc.energy;
          
          src_json["Cps"] = pksrc.cps;
          src_json["CpsStr"] = SpecUtils::printCompact(pksrc.cps, 4);
          
          src_json["BranchingRatio"] = pksrc.br;
          src_json["BranchingRatioStr"] = SpecUtils::printCompact( pksrc.br, 5 );
          
          snprintf( buffer, sizeof(buffer), "%.2f", pksrc.counts );
          src_json["Counts"] = pksrc.counts;
          src_json["CountsStr"] = string(buffer);
        }//for( const GammaInteractionCalc::PeakDetailSrc &pksrc : peak.m_sources )
        
        // We could loop over the volumetric sources, and add info, but...
        /*
        struct VolumeSrc
        {
          bool trace;
          double integral;
          double srcVolumetricActivity;
          bool inSituExponential;
          double inSituRelaxationLength;
          double detIntrinsicEff;
          std::string sourceName;
        };//struct VolumeSrc
        
        std::vector<VolumeSrc> m_volumetric_srcs;
         */
      }//for( loop over results.peak_calc_details )
    }//if( results.peak_calc_details )
                                      
    auto add_detector_to_json = []( nlohmann::json &data, const shared_ptr<const DetectorPeakResponse> &drf ){
      if( !drf || !drf->isValid() )
        return;
                                        
      auto &drf_obj = data["Detector"];
      drf_obj["Name"] = drf->name();
      drf_obj["Description"] = drf->description();
      drf_obj["Diameter_mm"] = drf->detectorDiameter() / PhysicalUnits::mm;
      drf_obj["Diameter_cm"] = drf->detectorDiameter() / PhysicalUnits::cm;
      drf_obj["Diameter_m"] = drf->detectorDiameter() / PhysicalUnits::m;
      drf_obj["Diameter_inch"] = drf->detectorDiameter() / (2.54*PhysicalUnits::cm);
      drf_obj["Radius_mm"] = 0.5*drf->detectorDiameter() / PhysicalUnits::mm;
      drf_obj["Radius_cm"] = 0.5*drf->detectorDiameter() / PhysicalUnits::cm;
      drf_obj["Radius_m"] = 0.5*drf->detectorDiameter() / PhysicalUnits::m;
      drf_obj["Radius_inch"] = 0.5*drf->detectorDiameter() / (2.54*PhysicalUnits::cm);
      drf_obj["Diameter"] = PhysicalUnits::printToBestLengthUnits( drf->detectorDiameter(), 3 );
      drf_obj["Radius"] = PhysicalUnits::printToBestLengthUnits( 0.5*drf->detectorDiameter(), 3 );
      drf_obj["FixedGeometry"] = drf->isFixedGeometry();
    };
             
      if( drf )
        add_detector_to_json( data, drf );
    
    if( results.peak_comparisons )
    {
      // This info is already in the Peaks info, but we'll put in anyways
      auto &energy_obj = data["PeakToModelComparison"];
      
      for( const GammaInteractionCalc::PeakResultPlotInfo &p : *results.peak_comparisons )
      {
        energy_obj["UsedPeaks"].push_back( {} );
        nlohmann::basic_json<> &p_json = energy_obj["UsedPeaks"].back();
        
        char buffer[64] = { '\0' };
        snprintf( buffer, sizeof(buffer), "%.2f", p.energy );
        p_json["Energy"] = string(buffer);
        p_json["NumSigmaOff"] = SpecUtils::printCompact( p.numSigmaOff, 6 ); //`(observed_counts - expected_counts) / observed_uncertainty`
        p_json["ObservedOverExpected"] = SpecUtils::printCompact( p.observedOverExpected, 6 );
        p_json["ObservedOverExpectedUncert"] = SpecUtils::printCompact( p.observedOverExpectedUncert, 6 );
      }
    }//if( results.peak_comparisons )
    
  }//void shield_src_fit_results_to_json()
  
  
  void add_hist_to_json( nlohmann::json &data,
                       const bool is_background,
                       const shared_ptr<const SpecUtils::Measurement> &spec_ptr,
                       const shared_ptr<const SpecMeas> &spec_file,
                       const std::set<int> &sample_numbers,
                       const string &filename,
                       const shared_ptr<const BatchPeak::BatchPeakFitResult> &peak_fit )
  {
    if( !spec_ptr )
      return;
    
    const SpecUtils::Measurement &spec = *spec_ptr;
    
    const double lt = spec.live_time();
    const double rt = spec.real_time();
    const string lt_str = PhysicalUnits::printToBestTimeUnits(lt,3);
    const string rt_str = PhysicalUnits::printToBestTimeUnits(rt,3);
    const vector<int> sample_nums( begin(sample_numbers), end(sample_numbers) );
    
    D3SpectrumExport::D3SpectrumOptions spec_json_options;
    if( peak_fit )
    {
      const BatchPeak::BatchPeakFitResult &fit_res = *peak_fit;
      const deque<std::shared_ptr<const PeakDef>> &fore_peaks = fit_res.fit_peaks;
      const vector<shared_ptr<const PeakDef> > inpeaks( begin(fore_peaks), end(fore_peaks) );
      spec_json_options.peaks_json = PeakDef::peak_json( inpeaks, spec_ptr );
    }//if( fit_results.m_peak_fit_results )
    
    //spec_json_options.line_color = "rgb(0,0,0)"; //black
    spec_json_options.peak_color = "rgba(0,51,255,0.6)";
    spec_json_options.title = "";
    spec_json_options.display_scale_factor = 1.0;
    spec_json_options.spectrum_type = SpecUtils::SpectrumType::Foreground;
    
    //We will only have foreground or background on this spectrum, so even if we say
    //  background ID is 2, instead of a negative number, things should be fine
    const size_t specID = is_background ? 2 : 1;
    const int backID = is_background ? -1 : 2;
    
    stringstream spec_strm;
    D3SpectrumExport::write_spectrum_data_js( spec_strm, spec, spec_json_options, specID, backID );
    const string spectrum_json_str = spec_strm.str();
    
    nlohmann::json spectrum_json;
    try
    {
      if( !spectrum_json_str.empty() )
        spectrum_json = nlohmann::json::parse( spectrum_json_str );
    }catch( std::exception &e )
    {
      cerr << "Failed to parse spectrum JSON: " << e.what()
      << "\n\nJSON: " << spectrum_json_str << endl << endl;
      assert( 0 );
    }
    
    const char * const label = is_background ? "background" : "foreground";
    auto &spec_obj = data[label];
    
    spec_obj["LiveTime"] = lt_str;
    spec_obj["RealTime"] = rt_str;
    spec_obj["LiveTime_s"] = lt;
    spec_obj["RealTime_s"] = rt;
    spec_obj["DeadTime"] = PhysicalUnits::printToBestTimeUnits(rt - lt);
    spec_obj["DeadTime_s"] = (rt - lt)/PhysicalUnits::second;
    spec_obj["DeadTime_percent"] = 100.0*(rt - lt) / rt;
    spec_obj["StartTime"] = SpecUtils::to_extended_iso_string( spec.start_time() );
    spec_obj["StartTime_iso"] = SpecUtils::to_iso_string( spec.start_time() );
    spec_obj["StartTime_vax"] = SpecUtils::to_vax_string( spec.start_time() );
    spec_obj["StartTimeIsValid"] = !SpecUtils::is_special( spec.start_time() );
    spec_obj["LowerSpectrumEnergy"] = spec.gamma_channel_lower(0);
    spec_obj["UpperSpectrumEnergy"] = spec.gamma_channel_upper(spec.num_gamma_channels() - 1);
    spec_obj["NumberChannels"] = (int)spec.num_gamma_channels();
    spec_obj["GammaSum"] = spec.gamma_count_sum();
    spec_obj["GammaCps"] = (spec.live_time() > 0.0) ? (spec.gamma_count_sum() / spec.live_time()) : -1.0;
    spec_obj["SampleNumbers"] = sample_nums;
    //detector names?
    spec_obj["Filename"] = filename;
    spec_obj["spectrum"] = spectrum_json;
    spec_obj["HasGps"] = spec.has_gps_info();
    if( spec.has_gps_info() )
    {
      spec_obj["Longitude"] = spec.longitude();
      spec_obj["Latitude"] = spec.latitude();
    }
    spec_obj["HasNeutrons"] = spec.contained_neutron();
    if( spec.contained_neutron() )
    {
      spec_obj["NeutronCounts"] = spec.neutron_counts_sum();
      spec_obj["NeutronLiveTime"] = spec.neutron_live_time();
      spec_obj["NeutronCps"] = spec.neutron_counts_sum() / spec.neutron_live_time();
    }
    
    if( !spec.title().empty() )
      spec_obj["SpectrumTitle"] = spec.title();
    
    //The measured ambient radiation dose equivalent rate value, in microsieverts per hour (μSv/h).
    if( spec.dose_rate() >= 0.0 )
      spec_obj["DoseRate_uSvPerHour"] = spec.dose_rate();
    
    //The measured radiation exposure rate value, in milliroentgen per hour (mR/h).
    if( spec.exposure_rate() >= 0.0 )
      spec_obj["ExposureRate_mRPerHour"] = spec.exposure_rate();
    
    if( !spec.detector_type().empty() )
      spec_obj["DetectorTypeDesc"] = spec.detector_type();
    
    if( !spec.detector_name().empty() )
      spec_obj["DetectorName"] = spec.detector_name();
    
    //SourceType source_type() const;
    
    if( !spec.remarks().empty() )
      spec_obj["SpectrumRemarks"] = spec.remarks();
    
    if( !spec.parse_warnings().empty() )
      spec_obj["SpectrumParseWarnings"] = spec.parse_warnings();
    
    spec_obj["IsDerivedData"] = (spec.derived_data_properties() != 0);
    
    spec_obj["DerivedIoiSum"] = static_cast<bool>(spec.derived_data_properties()
                                                  & static_cast<uint32_t>(SpecUtils::Measurement::DerivedDataProperties::ItemOfInterestSum));
    spec_obj["DerivedUsedForRidAnalysis"] = static_cast<bool>(spec.derived_data_properties()
                                                              & static_cast<uint32_t>(SpecUtils::Measurement::DerivedDataProperties::UsedForAnalysis));
    spec_obj["DerivedProcessedFurther"] = static_cast<bool>(spec.derived_data_properties()
                                                            & static_cast<uint32_t>(SpecUtils::Measurement::DerivedDataProperties::ProcessedFurther));
    spec_obj["DerivedBackgroundSub"] = static_cast<bool>(spec.derived_data_properties()
                                                         & static_cast<uint32_t>(SpecUtils::Measurement::DerivedDataProperties::BackgroundSubtracted));
    spec_obj["DerivedIsBackground"] = static_cast<bool>(spec.derived_data_properties()
                                                        & static_cast<uint32_t>(SpecUtils::Measurement::DerivedDataProperties::IsBackground));
    
    auto cal = spec.energy_calibration();
    auto &cal_obj = spec_obj["EnergyCal"];
    
    cal_obj["NumChannels"] = cal ? static_cast<int>(cal->num_channels()) : 0;
    const SpecUtils::EnergyCalType cal_type = cal ? cal->type() : SpecUtils::EnergyCalType::InvalidEquationType;
    switch( cal_type )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        cal_obj["Type"] = "Polynomial";
        cal_obj["Coefficients"] = vector<double>{ begin(cal->coefficients()), end(cal->coefficients()) };
        break;
        
      case SpecUtils::EnergyCalType::FullRangeFraction:
        cal_obj["Type"] = "FullRangeFraction";
        cal_obj["Coefficients"] = vector<double>{ begin(cal->coefficients()), end(cal->coefficients()) };
        break;
      case SpecUtils::EnergyCalType::LowerChannelEdge:
        cal_obj["Type"] = "LowerChannelEdge";
        break;
      case SpecUtils::EnergyCalType::InvalidEquationType:
        cal_obj["Type"] = "Invalid";
        break;
    }//switch( cal->type() )
    
    if( cal && !cal->deviation_pairs().empty() )
    {
      vector<vector<double>> pairs;
      for( const pair<float,float> &p : cal->deviation_pairs() )
        pairs.push_back( { static_cast<double>(p.first), static_cast<double>(p.second) } );
      cal_obj["DeviationPairs"] = pairs;
    }//if( dev pairs )
    
    if( spec_file )
    {
      spec_obj["InstrumentModel"] = spec_file->instrument_model();
      spec_obj["SerialNumber"] = spec_file->instrument_id();
      spec_obj["Manufacturer"] = spec_file->manufacturer();
      spec_obj["InstrumentType"] = spec_file->instrument_type();
      spec_obj["DetectorType"] = SpecUtils::detectorTypeToString( spec_file->detector_type() );
      spec_obj["NumberRecordsInFile"] = static_cast<int>( spec_file->num_measurements() );
      spec_obj["RemarksInFile"] = spec_file->remarks();
      spec_obj["ParseWarningsForFile"] = spec_file->parse_warnings();
      
      spec_obj["HasInstrumentRid"] = !!spec_file->detectors_analysis();
      if( spec_file->detectors_analysis() )
        spec_obj["InstrumentRidSummary"] = riidAnaSummary( spec_file );
    }//if( spec_file )
  }//add_hist_to_json(...)
  
  
  void add_activity_fit_options_to_json( nlohmann::json &data,
                                      const BatchActivity::BatchActivityFitOptions &options )
  {
    add_peak_fit_options_to_json( data, options );
    
    auto &options_obj = data["PeakFitOptions"];
    
    const bool overiding_dist = options.distance_override.has_value();
    options_obj["IsSpecifyingDistance"] = overiding_dist;
    if( overiding_dist )
    {
      const double dist = options.distance_override.value();
      options_obj["SpecifiedDistance_m"] = dist / PhysicalUnits::m;
      options_obj["SpecifiedDistance_cm"] = dist / PhysicalUnits::cm;
      options_obj["SpecifiedDistance"] = PhysicalUnits::printToBestLengthUnits( dist, 6 );
    }
    
    options_obj["UseBq"] = options.use_bq;
    options_obj["IsSpecifyingDetector"] = !!options.drf_override;
    if( options.drf_override )
      options_obj["SpecifiedDetectorName"] = options.drf_override->name();
    
    options_obj["HardBackgroundSubtracted"] = !!options.hard_background_sub;
  }//add_activity_fit_options_to_json(...)
  
  
  void add_exe_info_to_json( nlohmann::json &data )
  {
    data["InterSpecCompileDate"] = string(__DATE__);
    data["InterSpecCompileDateIso"] = to_string( AppUtils::compile_date_as_int() );
    const auto now = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
    data["AnalysisTime"] = SpecUtils::to_iso_string( now );
    data["CurrentWorkingDirectory"] = SpecUtils::get_working_path();
  #if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
    try
    {
      const string exe_path = AppUtils::current_exe_path();
      data["InterSpecExecutablePath"] = exe_path;
    }catch( std::exception & )
    {
      assert( 0 );
    }
  #endif //if( not web or mobile )
  }//void add_exe_info_to_json( nlohmann::json &data )
  

  void add_peak_fit_options_to_json( nlohmann::json &data, const BatchPeak::BatchPeakFitOptions &options )
  {
    auto &options_obj = data["PeakFitOptions"];
    options_obj["RefitEnergyCal"] = options.refit_energy_cal;
    options_obj["UseExemplarEnergyCal"] = options.use_exemplar_energy_cal;
    options_obj["WriteN42WithResults"] = options.write_n42_with_results;
    options_obj["ShowNonFitPeaks"] = options.show_nonfit_peaks;
    options_obj["OutputDir"] = options.output_dir;
    options_obj["CreateCsvOutput"] = options.create_csv_output;
    options_obj["CreateJsonOutput"] = options.create_json_output;
    options_obj["OverwriteOutputFiles"] = options.overwrite_output_files;
    
    if( !options.background_subtract_file.empty() )
    {
      options_obj["BackgroundSubFile"] = options.background_subtract_file;
      if( !options.background_subtract_samples.empty() )
      {
        options_obj["BackgroundSubSamples"] = vector<int>{ begin(options.background_subtract_samples),
          end(options.background_subtract_samples)
        };
      }
      options_obj["UsedExistingBackgroundPeak"] = options.use_existing_background_peaks;
      options_obj["UseExemplarEnergyCalForBackground"] = options.use_exemplar_energy_cal_for_background;
    }//if( !options.background_subtract_file.empty() )
    
    options_obj["ReportTemplateIncludeDir"] = options.template_include_dir;
    options_obj["ReportTemplates"] = options.report_templates;
    options_obj["ReportSummaryTemplates"] = options.summary_report_templates;
    options_obj["PeakStatThreshold"] = options.peak_stat_threshold;
    options_obj["PeakShapeThreshold"] = options.peak_hypothesis_threshold;
  }//add_peak_fit_options_to_json(...)
}//namespace BatchInfoLog
