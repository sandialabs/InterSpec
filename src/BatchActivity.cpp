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

#include "rapidxml/rapidxml.hpp"

#include "Minuit2/MnUserParameters.h"

#include "external_libs/SpecUtils/3rdparty/inja/inja.hpp"

#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"
#include "SpecUtils/D3SpectrumExport.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/BatchPeak.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/BatchActivity.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"


using namespace std;

namespace BatchActivity
{

shared_ptr<DetectorPeakResponse> init_drf_from_name( std::string drf_file, std::string drf_name )
{
  if( drf_file.empty() && drf_name.empty() )
    return nullptr;
    
  //  Check if we are over-riding the DRF to use
  if( !drf_file.empty() )
  {
    // If user specified a directory, assume its a GADRAS DRF
    if( SpecUtils::is_directory(drf_file) )
    {
      try
      {
        auto drf = make_shared<DetectorPeakResponse>();
        drf->fromGadrasDirectory( drf_file );
        return drf;
      }catch( std::exception &e )
      {
        throw runtime_error( "Directory '" + drf_file + "' wasnt a GADRAS DRF: " + string(e.what()) );
      }
    }//if( SpecUtils::is_directory(drf_file) )
    
    
    // Try a URI file
    if( SpecUtils::is_file(drf_file) && (SpecUtils::file_size(drf_file) < 64*1024) )
    {
      try
      {
        std::vector<char> data;
        SpecUtils::load_file_data( drf_file.c_str(), data );
        string url_query( begin(data), end(data) );
        SpecUtils::trim( url_query );
        auto drf = DetectorPeakResponse::parseFromAppUrl( url_query );
        if( drf )
          return drf;
      }catch( std::exception &e )
      {
        // Was not a file with URI as its contents
      }
    }else
    {
      try
      {
        string url_query = drf_file;
        SpecUtils::trim( url_query );
        auto drf = DetectorPeakResponse::parseFromAppUrl( url_query );
        if( drf )
          return drf;
      }catch( std::exception &e )
      {
        // Was not a URI
      }
    }//if( was a small file ) / else
    
    
    
    // Try a CSV/TSV file that may have multiple DRF in it
    fstream input( drf_file.c_str(), ios::binary | ios::in );
    if( !input )
      throw runtime_error( "Couldnt open DRF file '" + drf_file + "'." );
      
      
    vector<shared_ptr<DetectorPeakResponse>> drfs;
    try
    {
      vector<string> credits;
        
      // This should also cover files that could be parsed by
      //  `DetectorPeakResponse::parseSingleCsvLineRelEffDrf(string &)`
        
      DetectorPeakResponse::parseMultipleRelEffDrfCsv( input, credits, drfs );
    }catch(...)
    {
    }//try / catch for parseMultipleRelEffDrfCsv
      
      
    if( drfs.size() == 1 )
    {
      if( drf_name.empty() )
      {
        return drfs[0];
      }else
      {
        if( (drfs[0]->name() != drf_name) || (drfs[0]->description() != drf_name) )
          throw runtime_error( "The DRF in '" + drf_file + "' is named '"
                              + drfs[0]->name() + "', but you specified a name of '"
                                + drf_name + "'" );
        return drfs[0];
      }//if( drf_name.empty() ) / else
    }else if( drfs.size() > 1 )
    {
      string names;
      for( size_t i = 0; i < drfs.size(); ++i )
      {
        names += string(names.empty() ? "" : ", ") + "'" + drfs[i]->name() + "'";
        if( !drf_name.empty() && (drfs[i]->name() == drf_name) )
          return drfs[i];  // TODO: check if multiple DRFs with same name
      }
        
      if( drf_name.empty() )
        throw runtime_error( "No DRF name specified; possible DRFs in '" + drf_file
                              + "', are: " + names );
    }//if( drfs.size() == 1 ) / else
    
    
    // Try a Rel Eff CSV file
    input.seekg( 0, ios::beg );
      
    try
    {
      auto det = DrfSelect::parseRelEffCsvFile( drf_file );
      if( det )
        return det;
    }catch( std::exception &e )
    {
      
    }//try catch
    
    
    // Try a ECC file
    try
    {
      input.seekg( 0, ios::beg );
        
      auto answer = DetectorPeakResponse::parseEccFile( input );
      if( std::get<0>(answer) )
        return std::get<0>(answer);
    }catch( std::exception &e )
    {
      // Not a .ECC file
    }
  }else if( !drf_name.empty() )  //if( !drf_file.empty() )
  {
    try
    {
      string manufacturer = "";
      string model = drf_name;
      auto drf = DrfSelect::initARelEffDetector( SpecUtils::DetectorType::Unknown,
                                                manufacturer, model );
      if( drf )
        return drf;
    }catch( std::exception &e )
    {
      
    }//try /catch to load a standard Detector by name
    
    
    // TODO: Look for previous DRFs in the DB, that the user has used
      
  }//if( !drf_file.empty() ) / else if( !drf_name.empty() )
  
  throw runtime_error( "Could not load specified detector efficiency function."
                        + (drf_file.empty() ? string("") : " Filename='" + drf_file + "'.")
                        + (drf_name.empty() ? string("") : " Name='" + drf_name + "'.") );
  
  return nullptr;
}//shared_ptr<DetectorPeakResponse> init_drf_from_name( std::string drf_file, std::string drf-name )
  

  
  // Adds the basic direct info on a source (nuclide name, activity, age, etc), but does not
  //  Add which peaks it contributes to, or any information on gammas
void add_basic_src_details( const GammaInteractionCalc::SourceDetails &src,
                           const std::shared_ptr<const DetectorPeakResponse> &drf,
                           const bool useBq, 
                           nlohmann::basic_json<> &src_json )
{
  const DetectorPeakResponse::EffGeometryType eff_type = drf ? drf->geometryType()
                                              : DetectorPeakResponse::EffGeometryType::FarField;
  
  const string act_postfix = DetectorPeakResponse::det_eff_geom_type_postfix(eff_type);
  
  
  src_json["Nuclide"] = src.nuclide->symbol;
  src_json["Activity"] = PhysicalUnits::printToBestActivityUnits(src.activity,4,!useBq); //I dont think we need `act_postfix` here since I *think* this is total activity.
  src_json["ActivityBq"] = src.activity / PhysicalUnits::bq;
  src_json["ActivityCi"] = src.activity / PhysicalUnits::ci;
  src_json["ActivityIsFit"] = src.activityIsFit;
  if( src.activityIsFit )
  {
    src_json["ActivityUncert"] = PhysicalUnits::printToBestActivityUnits(src.activityUncertainty,4,!useBq);//ditto on `act_postfix` not here
    src_json["ActivityUncertBq"] = src.activityUncertainty / PhysicalUnits::bq;
    src_json["ActivityUncertCi"] = src.activityUncertainty / PhysicalUnits::ci;
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
    src_json["TraceActivityType"] = GammaInteractionCalc::to_str(src.traceActivityType);
    src_json["TraceDisplayActivity"] = PhysicalUnits::printToBestActivityUnits(src.traceSrcDisplayAct,4,!useBq) + act_postfix;
    src_json["TraceDisplayActivityBq"] = src.traceSrcDisplayAct / PhysicalUnits::bq;
    src_json["TraceDisplayActivityCi"] = src.traceSrcDisplayAct / PhysicalUnits::ci;
    
    if( src.ageIsFit )
    {
      src_json["TraceDisplayActivityUncert"] = PhysicalUnits::printToBestActivityUnits(src.traceSrcDisplayActUncertainty,4,!useBq) + act_postfix;
      src_json["TraceDisplayActivityUncertBq"] = src.traceSrcDisplayActUncertainty / PhysicalUnits::bq;
      src_json["TraceDisplayActivityUncertCi"] = src.traceSrcDisplayActUncertainty / PhysicalUnits::ci;
    }else
    {
      assert( src.traceSrcDisplayActUncertainty <= 0.0 );
    }
    
    if( src.traceActivityType == GammaInteractionCalc::TraceActivityType::ExponentialDistribution )
    {
      src_json["TraceRelaxationLength"] = PhysicalUnits::printToBestLengthUnits(src.traceRelaxationLength, 4);
      src_json["TraceRelaxationLengthCm"] = src.traceRelaxationLength / PhysicalUnits::cm;
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
  
  
/** Adds basic information about a peak (energy, fwhm, counts, etc), but not any information
   about gammas that contribute to it, etc
 */
void add_basic_peak_info( const GammaInteractionCalc::PeakDetail &peak, nlohmann::basic_json<> &peak_json )
{
  char buffer[64] = { '\0' };
  snprintf( buffer, sizeof(buffer), "%.2f", peak.energy );
  peak_json["Energy"] = string(buffer);
  peak_json["EnergyKeV"] = peak.energy;
  
  snprintf( buffer, sizeof(buffer), "%.3f", peak.decayParticleEnergy );
  peak_json["DecayParticleEnergy"] = string(buffer);
  peak_json["AssignedNuclide"] = peak.assignedNuclide;
  peak_json["FWHM"] = SpecUtils::printCompact(peak.fwhm, 4);
  
  snprintf( buffer, sizeof(buffer), "%.2f", peak.counts );
  peak_json["Counts"] = string(buffer);
  peak_json["CountsStr"] = SpecUtils::printCompact(peak.counts, 4);
  snprintf( buffer, sizeof(buffer), "%.2f", peak.countsUncert );
  peak_json["CountsUncert"] = string(buffer);
  peak_json["CountsUncertStr"] = SpecUtils::printCompact(peak.countsUncert, 4);
  peak_json["Cps"] = SpecUtils::printCompact(peak.cps, 4);
  peak_json["CpsUncert"] = SpecUtils::printCompact(peak.cpsUncert, 4);
  peak_json["ShieldAttenuations"] = peak.m_attenuations;
  
  if( peak.backgroundCounts > 0.0 )
  {
    snprintf( buffer, sizeof(buffer), "%.2f", peak.backgroundCounts );
    peak_json["BackgroundCounts"] = string(buffer);
    peak_json["BackgroundCountsStr"] = SpecUtils::printCompact(peak.backgroundCounts, 4);
    
    if( peak.backgroundCountsUncert > 0.0 )
    {
      snprintf( buffer, sizeof(buffer), "%.2f", peak.backgroundCountsUncert );
      peak_json["BackgroundCountsUncert"] = string(buffer);
      peak_json["BackgroundCountsUncertStr"] = SpecUtils::printCompact(peak.backgroundCountsUncert, 4);
    }
  }//if( peak.backgroundCounts > 0.0 )
  
  snprintf( buffer, sizeof(buffer), "%.2f", peak.observedCounts );
  peak_json["SignalCounts"] = string(buffer);
  peak_json["SignalCountsStr"] = SpecUtils::printCompact(peak.observedCounts, 4);
  snprintf( buffer, sizeof(buffer), "%.2f", peak.observedUncert );
  peak_json["SignalCountsUncert"] = string(buffer);
  peak_json["SignalCountsUncertStr"] = SpecUtils::printCompact(peak.observedUncert, 4);
  
  snprintf( buffer, sizeof(buffer), "%.2f", peak.expectedCounts );
  peak_json["PredictedCounts"] = string(buffer);
  peak_json["PredictedNumSigmaOff"] = SpecUtils::printCompact(peak.numSigmaOff, 4);
  
  peak_json["ObservedOverPredicted"] = SpecUtils::printCompact(peak.observedOverExpected, 5);
  peak_json["ObservedOverPredictedUncert"] = SpecUtils::printCompact(peak.observedOverExpectedUncert, 5);
  
  peak_json["DetectorSolidAngleFraction"] = SpecUtils::printCompact(peak.detSolidAngle, 5);
  peak_json["DetectorIntrinsicEff"] = SpecUtils::printCompact(peak.detIntrinsicEff, 5);
  peak_json["DetectorEff"] = SpecUtils::printCompact(peak.detEff, 5);
}//add_basic_peak_info( const PeakDetail &peak, nlohmann::basic_json<> &peak_json )
  
  
void add_gamma_info_for_peak( const GammaInteractionCalc::PeakDetail::PeakSrc &ps, 
                             const GammaInteractionCalc::SourceDetails * const src,
                             const std::shared_ptr<const DetectorPeakResponse> &drf,
                             const bool useBq,
                             nlohmann::basic_json<> &gamma_json )
{
  assert( ps.nuclide );
  
  char buffer[64] = { '\0' };
  snprintf( buffer, sizeof(buffer), "%.3f", ps.energy );
  
  gamma_json["Nuclide"] = ps.nuclide ? ps.nuclide->symbol : string("null");
  gamma_json["Energy"] = string(buffer);
  gamma_json["EnergyKeV"] = ps.energy;
  
  gamma_json["BranchingRatio"] = SpecUtils::printCompact( ps.br, 5 );
  gamma_json["BranchingRatioFloat"] = ps.br;
  
  gamma_json["CpsContributedToPeak"] = SpecUtils::printCompact( ps.cps, 5 );
  gamma_json["CpsContributedToPeakFloat"] = ps.cps;
  
  gamma_json["CountsContributedToPeak"] = SpecUtils::printCompact( ps.counts, 5 );
  gamma_json["CountsContributedToPeakFloat"] = ps.counts;
  
  gamma_json["HasDecayCorrection"] = (ps.decayCorrection > 0.0);
  if( ps.decayCorrection > 0.0 )
  {
    gamma_json["DecayCorrection"] = SpecUtils::printCompact(ps.decayCorrection, 4);
    gamma_json["DecayCorrectionFloat"] = ps.decayCorrection;
  }
  
  // The other information in `PeakDetail::PeakSrc` should be repeat of
  //  `GammaInteractionCalc::SourceDetails`.  We _could_ access all this information
  //  from the templating code, since we are in a loop over `SourceDetails`, but
  //  to make things easier on people, we'll just re-include it here.
  if( src )
    add_basic_src_details( *src, drf, useBq,  gamma_json );
};//add_gamma_info_for_peak(...)
  
  
void add_fit_options_to_json( const ShieldingSourceFitCalc::ShieldingSourceFitOptions &options,
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
    fit_setup["Distance"] = distance / PhysicalUnits::cm;
    fit_setup["DistanceStr"] = PhysicalUnits::printToBestLengthUnits( distance, 3 );
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
}//void add_fit_options_to_json(...)
  
  
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

  
void shield_src_fit_results_to_json( const ShieldingSourceFitCalc::ModelFitResults &results,
                                    const std::shared_ptr<const DetectorPeakResponse> &drf,
                                    const bool useBq,
                                    nlohmann::json &data )
{
  // Some info about the compiled application
  add_exe_info_to_json( data );
  
  //Now add info about the analysis setup
  add_fit_options_to_json( results.options, results.distance, results.geometry, drf, data );
  
  
  const DetectorPeakResponse::EffGeometryType eff_type = drf ? drf->geometryType()
                                              : DetectorPeakResponse::EffGeometryType::FarField;
  
  const string act_postfix = DetectorPeakResponse::det_eff_geom_type_postfix(eff_type);
  
  
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
  
  
  if( results.source_calc_details )
  {
    for( const GammaInteractionCalc::SourceDetails &src : *results.source_calc_details )
    {
      data["Sources"].push_back( {} );
      auto &src_json = data["Sources"].back();
      
      add_basic_src_details( src, drf, useBq, src_json );
      
      // We wont put peaks into this struct, but instead when we make the JSON, we'll
      //  insert peaks from `PeakDetail` as `PeakDetail::PeakSrc::nuclide` match this nuclide.
      if( results.peak_calc_details )
      {
        for( const GammaInteractionCalc::PeakDetail &peak : *results.peak_calc_details )
        {
          bool src_contribs_to_peak = false;
          for( const GammaInteractionCalc::PeakDetail::PeakSrc &ps :  peak.m_sources )
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
          
          for( const GammaInteractionCalc::PeakDetail::PeakSrc &ps :  peak.m_sources )
          {
            if( ps.nuclide != src.nuclide )
              continue;
            
            peak_json["ThisNucsGammasForPeak"].push_back( {} );
            nlohmann::basic_json<> &gamma_json = peak_json["ThisNucsGammasForPeak"].back();
            
            add_gamma_info_for_peak( ps, &src, drf, useBq, gamma_json );
          }
        }//for( loop over results.peak_calc_details )
      }//if( results.peak_calc_details )
    }//for( loop over results.source_calc_details )
  }//if( results.source_calc_details )
  
  
  if( results.shield_calc_details )
  {
    const vector<GammaInteractionCalc::ShieldingDetails> &shield_details
                                                          = *results.shield_calc_details;
    auto &shieldings_json = data["Shieldings"];
    
    shieldings_json["Geometry"] = GammaInteractionCalc::to_str(results.geometry);
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
        shield_json["Density"] = SpecUtils::printCompact(density, 5);
        shield_json["DensityFloat"] = density;
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
      }
      
      const vector<bool> fit_dim( shield.m_fit_dimension, shield.m_fit_dimension + shield.m_num_dimensions );
      shield_json["DimensionIsFit"] = fit_dim;
      
      shield_json["NumDimensions"] = static_cast<int>( shield.m_num_dimensions );
      shield_json["Geometry"] = GammaInteractionCalc::to_str(shield.m_geometry);
      
      shield_json["Thickness"] = PhysicalUnits::printToBestLengthUnits(shield.m_thickness,3);
      shield_json["ThicknessCm"] = shield.m_thickness / PhysicalUnits::cm;
      shield_json["VolumeCm3"] = shield.m_volume / PhysicalUnits::cm3;
      shield_json["VolumeUncertCm3"] = shield.m_volume_uncert / PhysicalUnits::cm3;
      
      shield_json["InnerRadius"] = PhysicalUnits::printToBestLengthUnits(shield.m_inner_rad, 3);
      shield_json["OuterRadius"] = PhysicalUnits::printToBestLengthUnits(shield.m_inner_rad + shield.m_thickness, 3);
      
      vector<double> inner_dims{ shield.m_inner_dimensions, shield.m_inner_dimensions + shield.m_num_dimensions };
      vector<double> outer_dims{ shield.m_outer_dimensions, shield.m_outer_dimensions + shield.m_num_dimensions };
      vector<double> thicknesses( shield.m_num_dimensions );
      vector<double> dim_uncerts( shield.m_dimension_uncert, shield.m_dimension_uncert + shield.m_num_dimensions );
      
      vector<string> inner_dims_strs( shield.m_num_dimensions );
      vector<string> outer_dims_strs( shield.m_num_dimensions );
      vector<string> thicknesses_strs( shield.m_num_dimensions );
      vector<string> dim_uncerts_strs( shield.m_num_dimensions );
      
      for( unsigned int dim = 0; dim < shield.m_num_dimensions; ++dim )
      {
        thicknesses[dim] = outer_dims[dim] - inner_dims[dim];
        thicknesses_strs[dim] = PhysicalUnits::printToBestLengthUnits( thicknesses[dim], 5 );
        inner_dims_strs[dim] = PhysicalUnits::printToBestLengthUnits( inner_dims[dim], 5 );
        outer_dims_strs[dim] = PhysicalUnits::printToBestLengthUnits( outer_dims[dim], 5 );
        if( fit_dim[dim] )
          dim_uncerts_strs[dim] = PhysicalUnits::printToBestLengthUnits( dim_uncerts[dim], 5 );
        
        inner_dims[dim] /= PhysicalUnits::cm;
        outer_dims[dim] /= PhysicalUnits::cm;
        thicknesses[dim] /= PhysicalUnits::cm;
        dim_uncerts[dim] /= PhysicalUnits::cm;
      }//for( unsigned int dim = 0; dim < shield.m_num_dimensions; ++dim )
      
      shield_json["InnerDimsCm"] = inner_dims;
      shield_json["OuterDimsCm"] = outer_dims;
      shield_json["ThicknessesCm"] = thicknesses;
      shield_json["ThicknessesUncerts"] = dim_uncerts_strs;
      shield_json["ThicknessesUncertsCm"] = dim_uncerts;
      shield_json["InnerDims"]   = inner_dims_strs;
      shield_json["OuterDims"]   = outer_dims_strs;
      shield_json["Thicknesses"] = thicknesses_strs;
      
      
      for( const GammaInteractionCalc::ShieldingDetails::SelfAttenComponent &comp : shield.m_mass_fractions )
      {
        shield_json["SelfAttenSources"].push_back( {} );
        auto &self_atten_json = shield_json["SelfAttenSources"].back();
      
        assert( comp.m_nuclide );
        self_atten_json["Nuclide"] = comp.m_nuclide ? comp.m_nuclide->symbol : string("null");
        self_atten_json["IsFittingMassFraction"] = comp.m_is_fit;
        self_atten_json["MassFraction"] = SpecUtils::printCompact(comp.m_mass_frac, 5);
        self_atten_json["MassFractionFloat"] = comp.m_mass_frac;
        if( comp.m_is_fit )
        {
          self_atten_json["MassFractionUncert"] = SpecUtils::printCompact(comp.m_mass_frac_uncert, 5);
          self_atten_json["MassFractionUncertFloat"] = comp.m_mass_frac_uncert;
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
            add_basic_src_details( *pos, drf, useBq, self_atten_json );
        }//if( results.source_calc_details )
      }//for( SelfAttenComponent & comp: shield.m_mass_fractions )
      
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
          trace_src_json["RelaxationLengthCm"] = trace.m_relaxation_length / PhysicalUnits::cm;
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
            add_basic_src_details( *pos, drf, useBq, trace_src_json );
        }//if( results.source_calc_details )
      }//for( const TraceSrcDetail &src : shield.m_trace_sources )
      
    }//for( const GammaInteractionCalc::ShieldingDetails &shield : *results.shield_calc_details )
  }//if( results.shield_calc_details )
  
  
  // Add Peak Details to JSON
  if( results.peak_calc_details )
  {
    auto &peaks_json = data["PeaksUsedForActivityFitting"];
    
    for( const GammaInteractionCalc::PeakDetail &peak : *results.peak_calc_details )
    {
      peaks_json["Peaks"].push_back( {} );
      nlohmann::basic_json<> &peak_json = peaks_json["Peaks"].back();
      add_basic_peak_info( peak, peak_json );
      
      for( const GammaInteractionCalc::PeakDetail::PeakSrc &pksrc : peak.m_sources )
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
        
        add_basic_src_details( *pos, drf, useBq, src_json );
        
        peak_json["HasDecayCorrection"] = (pksrc.decayCorrection > 0.0);
        if( pksrc.decayCorrection > 0.0 )
        {
          peak_json["DecayCorrection"] = SpecUtils::printCompact(pksrc.decayCorrection, 4);
          peak_json["DecayCorrectionFloat"] = pksrc.decayCorrection;
        }
        
        char buffer[64] = { '\0' };
        snprintf( buffer, sizeof(buffer), "%.2f", pksrc.energy );
        
        peak_json["Energy"] = string(buffer);
        peak_json["EnergyKeV"] = pksrc.energy;
        
        peak_json["Cps"] = SpecUtils::printCompact(pksrc.cps, 4);
        peak_json["CpsFloat"] = pksrc.cps;
        
        peak_json["BranchingRatio"] = SpecUtils::printCompact( pksrc.br, 5 );
        peak_json["BranchingRatioFloat"] = pksrc.br;
        
        snprintf( buffer, sizeof(buffer), "%.2f", pksrc.counts );
        peak_json["Counts"] = string(buffer);
        peak_json["CountsFloat"] = pksrc.counts;
      }//for( const GammaInteractionCalc::PeakDetail::PeakSrc &pksrc : peak.m_sources )
      
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
    drf_obj["Diameter"] = static_cast<double>( drf->detectorDiameter() );
    drf_obj["Radius"] = 0.5*drf->detectorDiameter();
    drf_obj["DiameterStr"] = PhysicalUnits::printToBestLengthUnits( drf->detectorDiameter(), 3 );
    drf_obj["RadiusStr"] = PhysicalUnits::printToBestLengthUnits( 0.5*drf->detectorDiameter(), 3 );
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

  
void fit_activities_in_files( const std::string &exemplar_filename,
                          const std::set<int> &exemplar_sample_nums,
                          const std::vector<std::string> &files,
                          const BatchActivityFitOptions &options )
{
  vector<string> warnings;
  
  if( files.empty() )
    throw runtime_error( "No input files specified." );
  
  if( !options.output_dir.empty() && !SpecUtils::is_directory(options.output_dir) )
    throw runtime_error( "Output directory ('" + options.output_dir + "'), is not a directory." );
    
  if( options.write_n42_with_results && options.output_dir.empty() )
    throw runtime_error( "If you specify to write N42 files with the fit peaks, you must specify an output directory." );
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  if( !db )
  {
    const char *message = "Unable to load the nuclide decay database."
    " Run the executable from a working directory where the 'data' folder is "
    " located - the program expects to load data/sandia.decay.xml.";
    
    throw runtime_error( message );
  }//if( !db )
  
  
  // Load report templates, and setup inja::environment
  //InterSpec::setStaticDataDirectory( SpecUtils::append_path(datadir,"data") );
  string default_tmplt_dir = InterSpec::staticDataDirectory();
  assert( !default_tmplt_dir.empty() );
  if( default_tmplt_dir.empty() )
    default_tmplt_dir = "./data";
  default_tmplt_dir = SpecUtils::append_path( default_tmplt_dir, ".." ); //Also see: `WServer::instance()->appRoot()`, which isnt valid right now
  default_tmplt_dir = SpecUtils::append_path( default_tmplt_dir, "InterSpec_resources" );
  default_tmplt_dir = SpecUtils::append_path( default_tmplt_dir, "static_text" );
  default_tmplt_dir = SpecUtils::append_path( default_tmplt_dir, "ShieldSourceFitLog" );
  
  // inja assumes trailing path separator, for its template path
#if( defined(_WIN32) )
  const char path_sep = '\\';
#else
  const char path_sep = '/';
#endif
  
  if( !default_tmplt_dir.empty() && (default_tmplt_dir.back() != path_sep) )
    default_tmplt_dir += path_sep;

  
  string tmplt_dir;
  
  if( SpecUtils::iequals_ascii( options.template_include_dir, "default" ) )
  {
    tmplt_dir = default_tmplt_dir;
  }else if( !SpecUtils::iequals_ascii( options.template_include_dir, "none" ) )
  {
    tmplt_dir = options.template_include_dir;
    if( !tmplt_dir.empty() && (tmplt_dir.back() != path_sep) )
      tmplt_dir += path_sep;
  }
  
  if( !tmplt_dir.empty() && !SpecUtils::is_directory(tmplt_dir) )
    throw runtime_error( string("Template include directory, '") + tmplt_dir
                        + "', doesnt look to be a valid directory - not performing analysis." );
  
  
  inja::Environment env{ tmplt_dir };
  env.set_trim_blocks( true ); // remove the first newline after a block
  //env.set_lstrip_blocks( true ); //strip the spaces and tabs from the start of a line to a block
  
#if( BUILD_FOR_WEB_DEPLOYMENT )
    env.set_search_included_templates_in_files( false );
#else
  //To think about more: is there any security issues with allowing.
  //  It doesnt *look* like inja prevents using things like "../../../SomeSensitveFile.txt" as templates.
  env.set_search_included_templates_in_files( !tmplt_dir.empty() );
#endif
  
  try
  {
    // If we're using a custom include path, opening templates from the default template location
    //  is problematic, so we'll just create a new `inja::Environment` to open the default
    //  templates, and then add them to `env` - not really tested yet.
    const string def_tmplt = SpecUtils::append_path( default_tmplt_dir, "std_fit_log.tmplt.txt" );
    inja::Environment sub_env;
    inja::Template tmplt = sub_env.parse_template( def_tmplt );
    env.include_template( "default-act-fit-results", tmplt );
     
    //TODO: Work on other templates we can include
  }catch( std::exception &e )
  {
    cerr << "Error loading default template: " << e.what() << endl;
    warnings.push_back( "Error loading default template: " + string(e.what()) );
  }
  
  
  auto output_filename = [&options]( const string &filename, const string tmplt ) -> string {
    string outname = SpecUtils::filename( filename );
    const string file_ext = SpecUtils::file_extension(outname);
    if( !file_ext.empty() )
      outname = outname.substr(0, outname.size() - file_ext.size());
    
    string tmplt_name = tmplt;
    string tmplt_ext = SpecUtils::file_extension(tmplt_name);
    size_t pos = SpecUtils::ifind_substr_ascii(tmplt_name, "tmplt");
    if( pos == string::npos )
      pos = SpecUtils::ifind_substr_ascii(tmplt_name, "template");
    if( pos != string::npos )
      tmplt_name = tmplt_name.substr(0, pos);
    if( SpecUtils::iends_with(tmplt_name, "_") || SpecUtils::iends_with(tmplt_name, ".") || SpecUtils::iends_with(tmplt_name, "-") )
      tmplt_name = tmplt_name.substr(0, tmplt_name.size() - 1);
    
    if( tmplt_ext.empty()
       || SpecUtils::iequals_ascii(tmplt_ext, "tmplt" )
       || SpecUtils::iequals_ascii(tmplt_ext, "template" ) )
      tmplt_ext = SpecUtils::file_extension(tmplt_name);
    
    if( tmplt_ext.empty() )
      tmplt_ext = ".txt";
    
    outname += (outname.empty() ? "" : "_") + tmplt_name + tmplt_ext;
    return SpecUtils::append_path(options.output_dir, outname );
  };//output_filename(...)
  
  
  
  bool set_setup_info_to_summary_json = false;
  nlohmann::json summary_json;
  
  add_exe_info_to_json( summary_json );
  
  summary_json["ExemplarFile"] = exemplar_filename;
  if( !exemplar_sample_nums.empty() )
    summary_json["ExemplarSampleNumbers"] = vector<int>{begin(exemplar_sample_nums), end(exemplar_sample_nums)};
  summary_json["InputFiles"] = files;
  
  std::shared_ptr<const SpecMeas> cached_exemplar_n42;
  for( const string filename : files )
  {
    const BatchActivityFitResult fit_results
                 = fit_activities_in_file( exemplar_filename, exemplar_sample_nums,
                                     cached_exemplar_n42, filename, options );
    
    if( !cached_exemplar_n42 )
      cached_exemplar_n42 = fit_results.m_exemplar_file;
    warnings.insert(end(warnings), begin(fit_results.m_warnings), end(fit_results.m_warnings) );
    
    if( !set_setup_info_to_summary_json && fit_results.m_fit_results )
    {
      std::shared_ptr<const DetectorPeakResponse> drf = options.drf_override;
      if( !drf && cached_exemplar_n42 )
        drf = cached_exemplar_n42->detector();
      
      const double distance = fit_results.m_fit_results->distance;
      const GammaInteractionCalc::GeometryType geometry = fit_results.m_fit_results->geometry;
      const ShieldingSourceFitCalc::ShieldingSourceFitOptions &fit_opts = fit_results.m_fit_results->options;
      add_fit_options_to_json( fit_opts, distance, geometry, drf, summary_json );
      set_setup_info_to_summary_json = true;
    }//if( !set_setup_info_to_summary_json )
    
    
    assert( (fit_results.m_result_code != BatchActivity::BatchActivityFitResult::ResultCode::Success)
             || fit_results.m_fit_results );
    
    if( (fit_results.m_result_code == BatchActivity::BatchActivityFitResult::ResultCode::Success)
       && fit_results.m_fit_results )
    {
      
      const bool useBq = InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
      
      cout << "Success analyzing '" << filename << "'!" << endl;
      assert( fit_results.m_fit_results );
      
      // The goal is to create a template that prints out the exact same info as as the current GUI
      //  log file, and then upgrade from this to hit HTML with the charts and peak fits and stuff.
      nlohmann::json data;
      
      auto add_hist_to_json = []( nlohmann::json &data,
                         const bool is_background,
                         const shared_ptr<const SpecUtils::Measurement> &spec_ptr,
                         const std::set<int> &sample_numbers,
                         const string &filename,
                         const shared_ptr<const BatchPeak::BatchPeakFitResult> &peak_fit ){
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
        data[label] = {
          { "LiveTimeSeconds", lt },
          { "LiveTimeStr", lt_str },
          { "RealTimeSeconds", rt },
          { "RealTimeStr", rt_str },
          { "LowerSpectrumEnergy", spec.gamma_channel_lower(0) },
          { "UpperSpectrumEnergy", spec.gamma_channel_upper(spec.num_gamma_channels() - 1) },
          { "NumberChannels", (int)spec.num_gamma_channels() },
          { "GammaSum", spec.gamma_count_sum() },
          { "ContainedNeutrons", spec.contained_neutron() },
          { "NeutronSum", spec.neutron_counts_sum() },
          { "SampleNumbers", sample_nums },
          //detector names?
          { "Filename", filename },
          { "spectrum", spectrum_json }
        };
      };//add_hist_to_json(...)
      
      auto add_peak_fit_options_to_json = []( nlohmann::json &data, const BatchPeak::BatchPeakFitOptions &options ){
        auto &options_obj = data["PeakFitOptions"];
        options_obj["RefitEnergyCal"] = options.refit_energy_cal;
        options_obj["UseExemplarEnergyCal"] = options.use_exemplar_energy_cal;
        options_obj["WriteN42WithResults"] = options.write_n42_with_results;
        options_obj["ShowNonFitPeaks"] = options.show_nonfit_peaks;
        options_obj["OutputDir"] = options.output_dir;
        options_obj["CreateCsvOutput"] = options.create_csv_output;
        options_obj["OverwriteOutputFiles"] = options.overwrite_output_files;
        
        if( !options.background_subtract_file.empty() )
          options_obj["BackgroundSubFile"] = options.background_subtract_file;
        if( !options.background_subtract_samples.empty() )
        {
          options_obj["BackgroundSubSamples"] = vector<int>{ begin(options.background_subtract_samples),
            end(options.background_subtract_samples)
          };
        }
        
        options_obj["ReportTemplateIncludeDir"] = options.template_include_dir;
        options_obj["ReportTemplates"] = options.report_templates;
        options_obj["ReportSummaryTemplate"] = options.summary_report_template;
      };
      
      
      if( fit_results.m_foreground )
      {
        add_hist_to_json( data, false, fit_results.m_foreground,
                         fit_results.m_foreground_sample_numbers, fit_results.m_filename,
                         fit_results.m_peak_fit_results );
      }//if( fit_results.m_foreground )
      
      if( fit_results.m_background )
      {
        string filename = fit_results.m_background_file ? fit_results.m_background_file->filename() 
                          : string();
        add_hist_to_json( data, true, fit_results.m_background,
                         fit_results.m_background_sample_numbers, filename,
                         fit_results.m_background_peak_fit_results );
      }
      
      
      add_peak_fit_options_to_json( data, fit_results.m_options );
    
      shared_ptr<const DetectorPeakResponse> drf = fit_results.m_options.drf_override;
      if( !drf && fit_results.m_exemplar_file )
        drf = fit_results.m_exemplar_file->detector();
      const bool use_bq = fit_results.m_options.use_bq;
      
      shield_src_fit_results_to_json( *fit_results.m_fit_results, drf, use_bq, data );
      
      data["Filepath"] = filename;
      data["Filename"] = SpecUtils::filename( filename );
      data["ParentDir"] = SpecUtils::parent_path( filename );
      
      
      summary_json["Files"].push_back( data );
      
      //std::cout << std::setw(4) << data << std::endl;
      
      for( string tmplt : options.report_templates )
      {
        try
        {
          string rpt;
          
          if( SpecUtils::iequals_ascii(tmplt, "default" ) )
          {
            rpt = env.render("{% include \"default-act-fit-results\" %}", data);
          }else
          {
            const bool is_in_inc = SpecUtils::is_file( SpecUtils::append_path(tmplt_dir, tmplt) );
            
            inja::Template injatmplt;
            if( is_in_inc )
            {
              injatmplt = env.parse_template( tmplt );
            }else
            {
              bool is_file = SpecUtils::is_file( tmplt );
              if( !is_file )
              {
                const string tmplt_in_def_path = SpecUtils::append_path(default_tmplt_dir, tmplt);
                
                is_file = SpecUtils::is_file( tmplt_in_def_path );
                
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
                // TODO: consider using `AppUtils::locate_file(...)` to find the file
#endif
                
                if( !is_file )
                  throw runtime_error( "Could not find template file '" + tmplt + "'."
                                      " Please specify full path to file, or use the"
                                      " 'report-template-include-dir' option to specify"
                                      " directory where reports are located." );
                tmplt = tmplt_in_def_path;
              }
                
              inja::Environment sub_env;
              injatmplt = sub_env.parse_template( tmplt );
            }//
            
            rpt = env.render( injatmplt, data );
          }//if( default report format ) / else
          
          if( options.to_stdout )
            cout << "\n\n" << rpt << endl << endl;
          
          if( !options.output_dir.empty() )
          {
            const string out_file = output_filename( filename, tmplt );
            
            if( SpecUtils::is_file(out_file) && !options.overwrite_output_files )
            {
              warnings.push_back( "Not writing '" + out_file + "', as it would overwrite a file."
                                 " See the '--overwrite-output-files' option to force writing." );
            }else
            {
#ifdef _WIN32
              const std::wstring woutcsv = SpecUtils::convert_from_utf8_to_utf16(out_file);
              std::ofstream output( woutcsv.c_str(), ios::binary | ios::out );
#else
              std::ofstream output( out_file.c_str(), ios::binary | ios::out);
#endif
              if( !output )
                warnings.push_back( "Failed to open report output '" + out_file + "', for writing.");
              else
                output.write( rpt.c_str(), rpt.size() );
            }//if( is file ) / else write file
          }//if( !options.output_dir.empty() )
        }catch( inja::InjaError &e )
        {
          const string msg = "Error templating results (" + e.type + ": line "
          + std::to_string(e.location.line) + ", column " + std::to_string(e.location.column)
          + " of '" + tmplt + "'): " + e.message + ". While processing '" + filename + "'.";
          
          cerr << msg << endl;
          warnings.push_back( msg );
        }catch( std::exception &e )
        {
          cerr << "Error templating results: " << e.what() << endl;
          warnings.push_back( "Error templating results: " + string(e.what()) );
        }
      }//for( const string &tmplt : options.report_templates )
      
      
      
      /*
       TODO: need to have `fit_activities_in_file` add peaks and shielding model to a file,
             perhaps in `fit_results.m_peak_fit_results->measurement`, or maybe better yet,
             create whole new std::shared_ptr<SpecMeas> in BatchActivityFitResult
      assert( fit_results.measurement );
      
      if( options.write_n42_with_results && fit_results.measurement )
      {
        string outn42 = SpecUtils::append_path(options.output_dir, SpecUtils::filename(filename) );
        if( !SpecUtils::iequals_ascii(SpecUtils::file_extension(filename), ".n42") )
          outn42 += ".n42";
        
        if( SpecUtils::is_file(outn42) && !options.overwrite_output_files )
        {
          warnings.push_back( "Not writing '" + outn42 + "', as it would overwrite a file."
                             " See the '--overwrite-output-files' option to force writing." );
        }else
        {
          if( !fit_results.measurement->save2012N42File( outn42 ) )
            warnings.push_back( "Failed to write '" + outn42 + "'.");
          else
            cout << "Have written '" << outn42 << "' with peaks" << endl;
        }
      }//if( options.write_n42_with_results )
       */
      
      if( fit_results.m_peak_fit_results )
      {
        const BatchPeak::BatchPeakFitResult &peak_fit_res = *fit_results.m_peak_fit_results;
        deque<shared_ptr<const PeakDef>> fit_peaks = peak_fit_res.fit_peaks;
        if( options.show_nonfit_peaks )
        {
          for( const auto p : peak_fit_res.unfit_exemplar_peaks )
          {
            auto np = make_shared<PeakDef>(*p);
            np->setAmplitude( 0.0 );
            np->setAmplitudeUncert( 0.0 );
            np->setSigmaUncert( 0.0 );
            auto cont = make_shared<PeakContinuum>( *np->continuum() );
            const size_t num_cont_pars = PeakContinuum::num_parameters(cont->type());
            for( size_t i = 0; i < num_cont_pars; ++i )
            {
              cont->setPolynomialCoef( i, 0.0 );
              cont->setPolynomialUncert( i, -1.0 );
            }
            np->setContinuum( cont );
            np->set_coefficient( -1.0, PeakDef::Chi2DOF );
            fit_peaks.push_back(np);
          }
          std::sort( begin(fit_peaks), end(fit_peaks), &PeakDef::lessThanByMeanShrdPtr );
        }//if( make_nonfit_peaks_zero )
        
        
        if( !options.output_dir.empty() && options.create_csv_output )
        {
          const string leaf_name = SpecUtils::filename(filename);
          string outcsv = SpecUtils::append_path(options.output_dir, leaf_name) + "_peaks.CSV";
          
          if( SpecUtils::is_file(outcsv) && !options.overwrite_output_files )
          {
            warnings.push_back( "Not writing '" + outcsv + "', as it would overwrite a file."
                               " See the '--overwrite-output-files' option to force writing." );
          }else
          {
#ifdef _WIN32
            const std::wstring woutcsv = SpecUtils::convert_from_utf8_to_utf16(outcsv);
            std::ofstream output_csv( woutcsv.c_str(), ios::binary | ios::out );
#else
            std::ofstream output_csv( outcsv.c_str(), ios::binary | ios::out );
#endif
            
            if( !output_csv )
            {
              warnings.push_back( "Failed to open '" + outcsv + "', for writing.");
            }else
            {
              PeakModel::write_peak_csv( output_csv, leaf_name, fit_peaks, peak_fit_res.spectrum );
              cout << "Have written '" << outcsv << "'" << endl;
            }
          }//if( SpecUtils::is_file( outcsv ) ) / else
        }//if( !options.output_dir.empty() )
        
        if( options.to_stdout )
        {
          const string leaf_name = SpecUtils::filename(filename);
          cout << "peaks for '" << leaf_name << "':" << endl;
          PeakModel::write_peak_csv( cout, leaf_name, fit_peaks, peak_fit_res.spectrum );
          cout << endl;
        }
      }//if( fit_results.m_peak_fit_results )
    }else
    {
      cout << "Failure: " << fit_results.m_error_msg << endl;
      warnings.push_back( "Failed in analyzing '" + filename + "': " + fit_results.m_error_msg );
    }
  }//for( const string filename : files )
    
  if( !options.summary_report_template.empty() 
     && !SpecUtils::iequals_ascii( options.summary_report_template, "none" ) )
  {
    try
    {
      string summary_tmplt = options.summary_report_template;
      
      if( SpecUtils::iequals_ascii( summary_tmplt, "default" ) )
        summary_tmplt = "std_summary.tmplt.txt";
      
      inja::Template injatmplt;
      
      const string inc_dir_summary_tmplt = SpecUtils::append_path(tmplt_dir, summary_tmplt);
      bool is_file = SpecUtils::is_file( inc_dir_summary_tmplt );
      if( is_file )
      {
        injatmplt = env.parse_template( summary_tmplt );
      }else
      {
        string found_path = summary_tmplt;
        is_file = SpecUtils::is_file( found_path );
        if( !is_file )
        {
          found_path = SpecUtils::append_path(default_tmplt_dir, summary_tmplt);
          is_file = SpecUtils::is_file( found_path );
        }
        
        if( !is_file )
          throw runtime_error( "Failed to find file '" + summary_tmplt + "'."
                              " Please specify full path to file, or use the"
                              " 'report-template-include-dir' option to specify"
                              " directory where the report is located." );
        
        inja::Environment sub_env;
        injatmplt = sub_env.parse_template( found_path );
      }
      
      const string rpt = env.render( injatmplt, summary_json );
      
      
      if( options.to_stdout )
        cout << "\n\n" << rpt << endl << endl;
      
      if( !options.output_dir.empty() )
      {
        const string out_file = output_filename( "", summary_tmplt );
        
        if( SpecUtils::is_file(out_file) && !options.overwrite_output_files )
        {
          warnings.push_back( "Not writing '" + out_file + "', as it would overwrite a file."
                             " See the '--overwrite-output-files' option to force writing." );
        }else
        {
#ifdef _WIN32
          const std::wstring woutcsv = SpecUtils::convert_from_utf8_to_utf16(out_file);
          std::ofstream output( woutcsv.c_str(), ios::binary | ios::out );
#else
          std::ofstream output( out_file.c_str(), ios::binary | ios::out );
#endif
          if( !output )
            throw runtime_error( "Failed to open summary report output, '" + out_file + "'" );
          
          output.write( rpt.c_str(), rpt.size() );
        }//if( file exists and dont overwrite ) / else
      }
    }catch( inja::InjaError &e )
    {
      const string msg = "Error templating summary output (" + e.type + ": line "
      + std::to_string(e.location.line) + ", column " + std::to_string(e.location.column)
      + " of '" + options.summary_report_template + "'): " + e.message + ".";
      
      cerr << msg << endl;
      warnings.push_back( msg );
    }catch( std::exception &e )
    {
      warnings.push_back( "Error making summary output: " + string(e.what()) );
    }
  }//if( !options.summary_report_template.empty() )
  
  
  if( !warnings.empty() )
    cerr << endl << endl;
  for( const auto warn : warnings )
    cerr << warn << endl;
}//fit_activities_in_files(...)
  
  

BatchActivityFitResult fit_activities_in_file( const std::string &exemplar_filename,
                          std::set<int> exemplar_sample_nums,
                          std::shared_ptr<const SpecMeas> cached_exemplar_n42,
                          const std::string &filename,
                          const BatchActivityFitOptions &options )
{
  //  TODO: allow specifying, not just in the exemplar N42.  Also note InterSpec defines a URL encoding for model as well
  BatchActivityFitResult result;
  result.m_options = options;
  result.m_result_code = BatchActivityFitResult::ResultCode::UnknownStatus;
  result.m_filename = filename;
  result.m_exemplar_sample_numbers = exemplar_sample_nums;
  
  
  if( !cached_exemplar_n42 )
  {
    auto tmp = make_shared<SpecMeas>();
    if( !tmp->load_file( exemplar_filename, SpecUtils::ParserType::Auto ) )
    {
      result.m_error_msg = "Could not load exemplar '" + exemplar_filename + "'.";
      result.m_result_code = BatchActivityFitResult::ResultCode::CouldntOpenExemplar;
      return result;
    }//if( !meas->load_file( filename, ParserType::Auto ) )
    
    cached_exemplar_n42 = tmp;
  }//if( !cached_exemplar_n42 )
  
  result.m_exemplar_file = cached_exemplar_n42;
  
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
  {
    result.m_error_msg = "Could not initialize nuclide DecayDataBase.";
    result.m_result_code = BatchActivityFitResult::ResultCode::CouldntInitializeStaticResources;
    return result;
  }
  
  MaterialDB matdb;
  const string materialfile = SpecUtils::append_path( InterSpec::staticDataDirectory(), "MaterialDataBase.txt" );
  try
  {
    matdb.parseGadrasMaterialFile( materialfile, db, false );
  }catch( std::exception &e )
  {
    result.m_error_msg = "Could not initialize shielding material database.";
    result.m_result_code = BatchActivityFitResult::ResultCode::CouldntInitializeStaticResources;
    return result;
  }
  
  const Material * const iron = matdb.material("Fe (iron)");
  if( !iron )
  {
    result.m_error_msg = "Couldnt access materials from shielding database.";
    result.m_result_code = BatchActivityFitResult::ResultCode::CouldntInitializeStaticResources;
    return result;
  }
  
  if( !cached_exemplar_n42 )
  {
    result.m_error_msg = "No exemplar spectrum file provided";
    result.m_result_code = BatchActivityFitResult::ResultCode::NoExemplar;
    return result;
  }//if( !cached_exemplar_n42 )
  
  
  auto specfile = make_shared<SpecMeas>();
  if( !specfile->load_file( filename, SpecUtils::ParserType::Auto ) )
  {
    result.m_error_msg = "Could not load foreground '" + filename + "'.";
    result.m_result_code = BatchActivityFitResult::ResultCode::CouldntOpenInputFile;
    return result;
  }//if( !meas->load_file( filename, ParserType::Auto ) )
  result.m_foreground_file = specfile;
  
  
  shared_ptr<SpecMeas> backfile;
  if( !options.background_subtract_file.empty() )
  {
    if( options.background_subtract_file == filename )
    {
      backfile = specfile;
      result.m_background_file = specfile;
    }else if( options.background_subtract_file == exemplar_filename )
    {
      //backfile = cached_exemplar_n42;
      result.m_background_file = cached_exemplar_n42;
    }else
    {
      auto tmp = make_shared<SpecMeas>();
      if( !tmp->load_file(options.background_subtract_file, SpecUtils::ParserType::Auto) )
      {
        result.m_error_msg = "Could not load background '" + options.background_subtract_file + "'.";
        result.m_result_code = BatchActivityFitResult::ResultCode::CouldntOpenBackgroundFile;
        return result;
      }
      
      backfile = tmp;
      result.m_background_file = backfile;
    }//if( options.background_subtract_file == filename ) / else
  }//if( !options.background_subtract_file.empty() )
  
  
  // Find the sample number to use for either the foreground, or background measurement
  auto find_sample = []( shared_ptr<const SpecMeas> meas, const SpecUtils::SourceType wanted ) -> int {
    // This logic is probably repeated in a number of places throughout InterSpec now...
    assert( (wanted == SpecUtils::SourceType::Foreground)
           || (wanted == SpecUtils::SourceType::Background) );
    
    if( (wanted != SpecUtils::SourceType::Foreground)
           && (wanted != SpecUtils::SourceType::Background) )
    {
      throw std::logic_error( "Invalid src type" );
    }
    
    if( !meas )
      throw runtime_error( "No SpecMeas passed in." );
    
    const vector<string> &detectors = meas->detector_names();
    const set<int> &sample_nums = meas->sample_numbers();
    
    set<int> foreground_samples, background_samples;
    
    for( const int sample : sample_nums )
    {
      bool classified_sample_num = false;
      for( const string &det : detectors )
      {
        auto m = meas->measurement( sample, det );
        if( !m )
          continue;
        
        switch( m->source_type() )
        {
          case SpecUtils::SourceType::IntrinsicActivity:
          case SpecUtils::SourceType::Calibration:
            break;
            
          case SpecUtils::SourceType::Background:
            classified_sample_num = true;
            background_samples.insert(sample);
            break;
            
          case SpecUtils::SourceType::Foreground:
          case SpecUtils::SourceType::Unknown:
            classified_sample_num = true;
            foreground_samples.insert(sample);
            break;
        }//switch( m->source_type() )
        
        if( classified_sample_num )
          break;
      }//for( const string &det : detectors )
    }//for( const int sample : sample_nums )
    
    const set<int> &samples = (wanted == SpecUtils::SourceType::Foreground) ? foreground_samples 
                                                                            : background_samples;
    if( samples.size() != 1 )
      throw runtime_error( "Sample number to use could not be uniquely identified." );
    
    return *begin(samples);
  };//int find_sample(...)
  
  
  shared_ptr<const SpecUtils::Measurement> foreground;
  try
  {
    const int sample_num = find_sample( specfile, SpecUtils::SourceType::Foreground );
    
    try
    {
      vector<shared_ptr<const SpecUtils::Measurement>> meass = specfile->sample_measurements(sample_num);
      if( meass.size() == 1 )
        foreground = meass[0];
      else
        foreground = specfile->sum_measurements( {sample_num}, specfile->detector_names(), nullptr );
      if( !foreground )
        throw runtime_error( "Missing measurement." );
    }catch( std::exception &e )
    {
      result.m_error_msg = "Could getting foreground measurement: " + string(e.what());
      result.m_result_code = BatchActivityFitResult::ResultCode::ForegroundSampleNumberUnderSpecified;
      return result;
    }//Try / catch get
    
    result.m_foreground = foreground;
    result.m_foreground_sample_numbers = set<int>{ sample_num };
    assert( result.m_foreground_file );
  }catch( std::exception &e )
  {
    result.m_error_msg = "Could not determine foreground: " + string(e.what());
    result.m_result_code = BatchActivityFitResult::ResultCode::ForegroundSampleNumberUnderSpecified;
    return result;
  }// try / catch (find foreground)

  
  set<int> background_sample_nums;
  shared_ptr<const SpecUtils::Measurement> background;
  try
  {
    if( options.background_subtract_file == exemplar_filename )
    {
      if( options.background_subtract_samples.empty() )
      {
        const int sample_num = find_sample( cached_exemplar_n42, SpecUtils::SourceType::Background );
        background_sample_nums = set<int>{sample_num};
      }else
      {
        background_sample_nums = options.background_subtract_samples;
      }
      
      try
      {
        background = specfile->sum_measurements( background_sample_nums, cached_exemplar_n42->detector_names(), nullptr );
      }catch( std::exception &e )
      {
        result.m_error_msg = "Error summing background measurements from exemplar: " + string(e.what());
        result.m_result_code = BatchActivityFitResult::ResultCode::BackgroundSampleNumberUnderSpecified;
        return result;
      }//Try / catch get
    }else if( backfile )
    {
      if( !options.background_subtract_samples.empty() )
      {
        const SpecUtils::SourceType t = (backfile == specfile)
        ? SpecUtils::SourceType::Background
        : SpecUtils::SourceType::Foreground;
        
        const int sample_num = find_sample( backfile, t );
        
        background_sample_nums = set<int>{sample_num};
        
        try
        {
          vector<shared_ptr<const SpecUtils::Measurement>> meass = backfile->sample_measurements(sample_num);
          if( meass.size() == 1 )
            background = meass[0];
          else
            background = backfile->sum_measurements( background_sample_nums, backfile->detector_names(), nullptr );
          if( !background )
            throw runtime_error( "Missing measurement." );
        }catch( std::exception &e )
        {
          result.m_error_msg = "Could not get background measurement: " + string(e.what());
          result.m_result_code = BatchActivityFitResult::ResultCode::BackgroundSampleNumberUnderSpecified;
          return result;
        }//Try / catch get
      }else
      {
        background_sample_nums = options.background_subtract_samples;
        try
        {
          background = backfile->sum_measurements( background_sample_nums, backfile->detector_names(), nullptr );
          if( !background )
            throw runtime_error( "Missing measurement." );
        }catch( std::exception &e )
        {
          result.m_error_msg = "Could not sum background samples: " + string(e.what());
          result.m_result_code = BatchActivityFitResult::ResultCode::BackgroundSampleNumberUnderSpecified;
          return result;
        }//Try / catch get
      }//if( options.background_subtract_samples.empty() ) / else
    }//if( options.background_subtract_file == exemplar_filename ) / else
    
    result.m_background = background;
    result.m_background_sample_numbers = background_sample_nums;
    assert( (!!result.m_background_file) == (!!result.m_background) );
    assert( options.background_subtract_file.empty() == (!background) );
    assert( options.background_subtract_file.empty() == (!result.m_background_file) );
  }catch( std::exception &e )
  {
    result.m_error_msg = "Could not determine background: " + string(e.what());
    result.m_result_code = BatchActivityFitResult::ResultCode::BackgroundSampleNumberUnderSpecified;
    return result;
  }// try / catch (find foreground)
  
  
  // Now need to figure out the foreground and possibly background sample numbers in the exemplar
  BatchPeak::BatchPeakFitOptions peak_fit_options = options;
  peak_fit_options.background_subtract_file = "";
  

  const BatchPeak::BatchPeakFitResult foreground_peak_fit_result
              = BatchPeak::fit_peaks_in_file( exemplar_filename, exemplar_sample_nums,
                                            cached_exemplar_n42, filename, specfile,
                                            result.m_foreground_sample_numbers, peak_fit_options );
  
  result.m_peak_fit_results = make_shared<BatchPeak::BatchPeakFitResult>( foreground_peak_fit_result );
  
  if( !foreground_peak_fit_result.success )
  {
    result.m_error_msg = "Fitting of foreground peaks failed.";
    result.m_result_code = BatchActivityFitResult::ResultCode::ForegroundPeakFitFailed;
  }//if( !foreground_peaks.success )
  
  
  if( !options.background_subtract_file.empty() )
  {
    //if( options.background_subtract_file == exemplar_filename )
    //First see if peaks are already fit - if so, use them, otherwise try to fit them
    
    // TODO: we are fitting the foreground peaks to the background - we should allow specifying exemplar background peaks - either to refit, or to just use
    cerr << "TODO: we are fitting the foreground peaks to the background - we should allow specifying exemplar background peaks - either to refit, or to just use" << endl;
    const BatchPeak::BatchPeakFitResult background_peaks
                = BatchPeak::fit_peaks_in_file( exemplar_filename, 
                                               exemplar_sample_nums,
                                              cached_exemplar_n42, 
                                               options.background_subtract_file,
                                               backfile, 
                                               result.m_background_sample_numbers,
                                               peak_fit_options );
    
    
    
    result.m_background_peak_fit_results = make_shared<BatchPeak::BatchPeakFitResult>( background_peaks );
  }//if( backfile )
  
  
  const rapidxml::xml_document<char> *model_xml = cached_exemplar_n42->shieldingSourceModel();
  const rapidxml::xml_node<char> *base_node = model_xml ? model_xml->first_node() : nullptr;
  if( !base_node )
  {
    result.m_error_msg = "Exemplar file did not have a Shielding/Source fit model in it";
    result.m_result_code = BatchActivityFitResult::ResultCode::NoInputSrcShieldModel;
    return result;
  }
  
  
  deque<shared_ptr<const PeakDef>> foreground_peaks = foreground_peak_fit_result.fit_peaks;
  if( foreground_peaks.empty() )
  {
    result.m_error_msg = "No foreground peaks fit.";
    result.m_result_code = BatchActivityFitResult::ResultCode::NoFitForegroundPeaks;
    return result;
  }
  
  shared_ptr<const DetectorPeakResponse> detector = options.drf_override;
  if( !detector )
    detector = cached_exemplar_n42->detector();
  
  if( !detector )
  {
    result.m_error_msg = "No detector efficiency function specified.";
    result.m_result_code = BatchActivityFitResult::ResultCode::NoDetEffFnct;
    return result;
  }
  
  
  double distance = 1*PhysicalUnits::meter;
  if( !detector->isFixedGeometry() )
  {
    try
    {
      const rapidxml::xml_node<char> *dist_node = base_node->first_node( "Distance" );
      const string diststr = SpecUtils::xml_value_str( dist_node );
      distance = PhysicalUnits::stringToDistance( diststr );
    }catch( std::exception &e )
    {
      result.m_error_msg = "Failed to get distance.";
      result.m_result_code = BatchActivityFitResult::ResultCode::InvalidDistance;
      return result;
    }//try / catch to get distance
  }//if( !detector->isFixedGeometry() )
  
  GammaInteractionCalc::GeometryType geometry = GammaInteractionCalc::GeometryType::NumGeometryType;
  
  const rapidxml::xml_node<char> *geom_node = base_node->first_node( "Geometry" );
  const string geomstr = SpecUtils::xml_value_str( geom_node );
  
  for( GammaInteractionCalc::GeometryType type = GammaInteractionCalc::GeometryType(0);
      type != GammaInteractionCalc::GeometryType::NumGeometryType;
      type = GammaInteractionCalc::GeometryType(static_cast<int>(type) + 1) )
  {
    if( SpecUtils::iequals_ascii(geomstr, GammaInteractionCalc::to_str(type)) )
    {
      geometry = type;
      break;
    }
  }//for( loop over geometry types )
  
  if( geometry == GammaInteractionCalc::GeometryType::NumGeometryType )
  {
    result.m_error_msg = "Geometry type not specified.";
    result.m_result_code = BatchActivityFitResult::ResultCode::InvalidGeometry;
    return result;
  }//if( geometry == GammaInteractionCalc::GeometryType::NumGeometryType )
  
  
  vector<ShieldingSourceFitCalc::ShieldingInfo> shield_definitions;
  vector<ShieldingSourceFitCalc::SourceFitDef> src_definitions;
  ShieldingSourceFitCalc::ShieldingSourceFitOptions fit_options;
  try
  {
    fit_options.deSerialize( base_node );
  }catch( std::exception &e )
  {
    result.m_error_msg = "Exemplar file '" + exemplar_filename + "' Options invalid: " + string(e.what());
    result.m_result_code = BatchActivityFitResult::ResultCode::InvalidFitOptions;
    return result;
  }

  // Check that if options said to use background-subtraction, that we actually can.
  //  I'm a little mixed here - perhaps if no background is specified, we should just ignore this,
  //  or just create a warning???
  //  For now - lets just warn.
  shared_ptr<const deque<std::shared_ptr<const PeakDef>>> background_peaks;
  if( background && result.m_background_peak_fit_results && result.m_background_peak_fit_results->success )
  {
    background_peaks = make_shared<deque<shared_ptr<const PeakDef>>>(result.m_background_peak_fit_results->fit_peaks);
  }
  
  
  
  if( fit_options.background_peak_subtract && (!background_peaks || background_peaks->empty()) )
  {
    fit_options.background_peak_subtract = false;
    string msg = "Exemplar file '" + exemplar_filename + "' Options indicated"
    " background-subtract, but either no background was given, or no background peaks fit.";
    result.m_warnings.push_back( msg );
    //result.m_error_msg = msg;
    //result.m_result_code = BatchActivityFitResult::ResultCode::ExemplarUsedBackSubButNoBackground;
    //return result;
  }
                                           
  const rapidxml::xml_node<char> *shieldings_node = base_node->first_node( "Shieldings" );
  if( !shieldings_node )
  {
    result.m_error_msg = "No Shieldings node specified in Shield/Src fit model.";
    result.m_result_code = BatchActivityFitResult::ResultCode::NoShieldingsNode;
    return result;
  }
      
  XML_FOREACH_CHILD(shield_node, shieldings_node, "Shielding")
  {
    try
    {
      ShieldingSourceFitCalc::ShieldingInfo info;
      info.deSerialize( shield_node, &matdb );
      shield_definitions.push_back( std::move(info) );
    }catch( std::exception &e )
    {
      result.m_error_msg = "Invalid Shielding node: " + string(e.what());
      result.m_result_code = BatchActivityFitResult::ResultCode::ErrorParsingShielding;
      return result;
    }
  }//XML_FOREACH_CHILD(shield_node, shieldings_node, "Shielding")
      
      
  const rapidxml::xml_node<char> *srcs_node = base_node->first_node( "Nuclides" );
  if( !srcs_node )
  {
    result.m_error_msg = "Exemplar file '" + exemplar_filename + "' missing Nuclides (sources) node.";
    result.m_result_code = BatchActivityFitResult::ResultCode::MissingNuclidesNode;
    return result;
  }//if( !srcs_node )
      
  XML_FOREACH_CHILD( src_node, srcs_node, "Nuclide" )
  {
    try
    {
      ShieldingSourceFitCalc::SourceFitDef info;
      info.deSerialize( src_node );
      src_definitions.push_back( std::move(info) );
    }catch( std::exception &e )
    {
      result.m_error_msg = "Exemplar file '" + exemplar_filename + "' has invalid Nuclides node: " + string(e.what());
      result.m_result_code = BatchActivityFitResult::ResultCode::InvalidNuclideNode;
      return result;
    }// try / catch
  }//XML_FOREACH_CHILD( src_node, srcs_node, "Nuclide" )
      
  
  if( src_definitions.empty() )
  {
    result.m_error_msg = "Exemplar file '" + exemplar_filename + "' no sources defined.";
    result.m_result_code = BatchActivityFitResult::ResultCode::NoSourceNuclides;
    return result;
  }//if( src_definitions.empty() )
  
  // TODO: We may not have fit peaks for all `src_definitions`.  We should remove these sources, and add warnings about it

  // We have all the parts, lets do the computation:
  pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars =
  GammaInteractionCalc::ShieldingSourceChi2Fcn::create( distance, geometry,
                                                           shield_definitions, src_definitions, detector,
                                                           foreground, background, foreground_peaks, background_peaks, fit_options );
      
  auto inputPrams = make_shared<ROOT::Minuit2::MnUserParameters>();
  *inputPrams = fcn_pars.second;
      
  auto progress = make_shared<ShieldingSourceFitCalc::ModelFitProgress>();
  auto fit_results = make_shared<ShieldingSourceFitCalc::ModelFitResults>();
  result.m_fit_results = fit_results;
  
  auto progress_fcn = [progress](){
    // We dont really care too much probably right now
  };
      
  bool finished_fit_called = false;
  auto finished_fcn = [fit_results, &finished_fit_called](){
    finished_fit_called = true;
  };
      
  ShieldingSourceFitCalc::fit_model( "", fcn_pars.first, inputPrams, progress, progress_fcn, fit_results, finished_fcn );
      
  assert( finished_fit_called );
  
  result.m_result_code = BatchActivityFitResult::ResultCode::Success;
  
  if( fit_results->successful != ShieldingSourceFitCalc::ModelFitResults::FitStatus::Final )
  {
    result.m_error_msg = "Fit not successful.";
    result.m_result_code = BatchActivityFitResult::ResultCode::FitNotSuccessful;
    return result;
  }
  
      
  for( const ShieldingSourceFitCalc::SourceFitDef insrc : src_definitions )
  {
    assert( insrc.nuclide );
    bool found_in_output = false;
    for( const ShieldingSourceFitCalc::IsoFitStruct &fitsrc : fit_results->fit_src_info )
    {
      assert( fitsrc.nuclide );
      found_in_output = (fitsrc.nuclide == insrc.nuclide);
      if( found_in_output )
        break;
    }//for( const ShieldingSourceFitCalc::IsoFitStruct &fitsrc : fit_results->fit_src_info )
    
    if( !found_in_output )
    {
      result.m_warnings.push_back( "Fit results does not include " + insrc.nuclide->symbol );
      result.m_result_code = BatchActivityFitResult::ResultCode::DidNotFitAllSources;
    }//if( !found_in_output )
  }//for( const ShieldingSourceFitCalc::SourceFitDef insrc : src_definitions )
  
  // TODO: create an output file that has peaks, drf, and the fit model.  This will mean creating the XML that represents the fit model then turning it into a string
  
  cout << "Finished fitting for '" << filename << "'" << endl;
  
  return result;
}//fit_activities_in_file(...)
  
}//namespace BatchPeak

