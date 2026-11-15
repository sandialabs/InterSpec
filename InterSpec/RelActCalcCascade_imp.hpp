#ifndef RelActCalcCascade_imp_hpp
#define RelActCalcCascade_imp_hpp
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

#include <cmath>
#include <vector>
#include <cassert>
#include <optional>
#include <algorithm>
#include <stdexcept>
#include <type_traits>

#include "cascade/CascadeSummingProgram.h"

#include "InterSpec/MaterialDB.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/RelActCalc_imp.hpp"
#include "InterSpec/RelActCalcCascade.h"
#include "InterSpec/CascadeSummingCalc.h"

/** Templated implementations of RelActCalcCascade.h; `T` is `double` or a `ceres::Jet`. */
namespace RelActCalcCascade
{
  /** The value of `v`, without derivatives. */
  template<typename T>
  double scalar_value( const T &v )
  {
    if constexpr ( std::is_arithmetic_v<T> )
      return static_cast<double>( v );
    else
      return v.a;
  }//scalar_value(...)


  template<typename T>
  T on_axis_solid_angle( const T &distance, const double detector_radius )
  {
    using std::sqrt;  //Jets find theirs by argument dependent lookup

    // 0.5*(1 - d/s), s = sqrt(d^2 + r^2); written so it does not cancel at large distance.
    const T s = sqrt( distance*distance + detector_radius*detector_radius );
    return (0.5*detector_radius*detector_radius) / (s*(s + distance));
  }//on_axis_solid_angle(...)


  template<typename T>
  T on_axis_distance( const T &solid_angle, const double detector_radius )
  {
    using std::sqrt;

    // cos(theta) = 1 - 2g, with theta the half-angle of the disk; distance = r cos(theta)/sin(theta).
    return (detector_radius*(1.0 - 2.0*solid_angle)) / (2.0*sqrt( solid_angle*(1.0 - solid_angle) ));
  }//on_axis_distance(...)


  template<typename T>
  void SummingModel::set_attenuation( const std::optional<RelActCalc::PhysModelShield<T>> &self_atten,
                                      const std::vector<RelActCalc::PhysModelShield<T>> &external_attens,
                                      EvalInput<T> &input ) const
  {
    using std::exp;

    const size_t num_energies = m_energies.size();

    input.self_mu.clear();
    input.self_areal_density = T(0.0);
    input.self_atomic_number = T(0.0);

    if( self_atten.has_value() )
    {
      input.self_mu.resize( num_energies, T(0.0) );
      input.self_areal_density = self_atten->areal_density;

      if( self_atten->material )
      {
        const std::shared_ptr<const std::vector<double>> mu = material_mu( self_atten->material );
        for( size_t i = 0; i < num_energies; ++i )
          input.self_mu[i] = T( (*mu)[i] );
        input.self_atomic_number = T( static_cast<double>( self_atten->material->massWeightedAtomicNumber() ) );
      }else
      {
        for( size_t i = 0; i < num_energies; ++i )
          input.self_mu[i] = RelActCalc::get_atten_coef_for_an( self_atten->atomic_number, m_energies[i] );
        input.self_atomic_number = self_atten->atomic_number;
      }
    }//if( self_atten.has_value() )

    input.external_transmission.assign( num_energies, T(1.0) );
    input.external_areal_density = T(0.0);
    input.external_an_areal_density = T(0.0);
    for( const RelActCalc::PhysModelShield<T> &shield : external_attens )
    {
      input.external_areal_density += shield.areal_density;

      if( shield.material )
      {
        const std::shared_ptr<const std::vector<double>> mu = material_mu( shield.material );
        for( size_t i = 0; i < num_energies; ++i )
          input.external_transmission[i] *= exp( -(*mu)[i] * shield.areal_density );
        input.external_an_areal_density
                  += static_cast<double>( shield.material->massWeightedAtomicNumber() ) * shield.areal_density;
      }else
      {
        for( size_t i = 0; i < num_energies; ++i )
          input.external_transmission[i]
                    *= exp( -RelActCalc::get_atten_coef_for_an( shield.atomic_number, m_energies[i] ) * shield.areal_density );
        input.external_an_areal_density += shield.atomic_number * shield.areal_density;
      }
    }//for( loop over external attenuators )
  }//set_attenuation(...)


  template<typename T>
  T SummingModel::near_field_ratio( const std::vector<double> &ratios, const std::vector<double> &slopes,
                                    const T &solid_angle ) const
  {
    const std::vector<double> &grid = m_near_solid_angle;
    assert( (grid.size() >= 2) && (ratios.size() == grid.size()) && (slopes.size() == grid.size()) );

    // Toward far away (and for a solid angle that is not a number, which then comes back as one) carry on
    //  along the first slope; past touching, where no source can be, hold the last value.
    const double value = scalar_value( solid_angle );
    if( !(value > grid.front()) )
      return ratios.front() + slopes.front()*(solid_angle - grid.front());
    if( value >= grid.back() )
      return T( ratios.back() );

    const size_t upper = static_cast<size_t>( std::upper_bound( begin(grid), end(grid), value ) - begin(grid) );
    const size_t lower = upper - 1;
    const double width = grid[upper] - grid[lower];
    const T t = (solid_angle - grid[lower]) / width;
    const T t2 = t*t, t3 = t2*t;

    // Cubic Hermite: continuous, with a continuous first derivative.
    return (2.0*t3 - 3.0*t2 + 1.0)*ratios[lower] + (t3 - 2.0*t2 + t)*(width*slopes[lower])
           + (-2.0*t3 + 3.0*t2)*ratios[upper] + (t3 - t2)*(width*slopes[upper]);
  }//near_field_ratio(...)


  template<typename T>
  void SummingModel::point_efficiencies( const T &solid_angle, std::vector<T> &fep, std::vector<T> &total ) const
  {
    const size_t num_energies = m_energies.size();
    fep.assign( num_energies, T(0.0) );
    total.assign( num_energies, T(0.0) );
    for( size_t i = 0; i < num_energies; ++i )
    {
      fep[i] = solid_angle * m_detector.fep[i];
      total[i] = solid_angle * m_detector.total[i];

      // A geometry-aware response's own near-field efficiencies, as ratios to "solid angle x intrinsic".
      if( !m_near_solid_angle.empty() )
      {
        fep[i] *= near_field_ratio( m_near_fep_ratio[i], m_near_fep_slope[i], solid_angle );
        total[i] *= near_field_ratio( m_near_tot_ratio[i], m_near_tot_slope[i], solid_angle );
      }
    }//for( loop over energies )
  }//point_efficiencies(...)


  template<typename T>
  void SummingModel::evaluate( const EvalInput<T> &input, EvalResult<T> &result ) const
  {
    using std::exp;
    using std::sqrt;

    const size_t num_energies = m_energies.size();
    const bool self_atten = !input.self_mu.empty();

    assert( input.external_transmission.size() == num_energies );
    assert( !self_atten || (input.self_mu.size() == num_energies) );
    if( (input.external_transmission.size() != num_energies)
       || (self_atten && (input.self_mu.size() != num_energies)) )
      throw std::logic_error( "SummingModel::evaluate: attenuation not given at the model's energies." );

    // Without self-attenuation every point of the source is equivalent: one node.
    const DepthNodes &nodes = depth_nodes();
    const size_t num_nodes = self_atten ? nodes.fraction.size() : size_t(1);

    // The unattenuated absolute efficiencies; and the full-energy one behind the external attenuators, which
    //  is the same for every node.
    const T &front_solid_angle = input.front_solid_angle;
    std::vector<T> fep_bare, tot_bare;
    point_efficiencies( front_solid_angle, fep_bare, tot_bare );
    std::vector<T> fep_const( num_energies );
    for( size_t i = 0; i < num_energies; ++i )
      fep_const[i] = fep_bare[i] * input.external_transmission[i];

    // Photons that scatter in what is in front of a node, and still reach the detector, add to the total
    //  efficiency; the table is in g/cm2.
    const bool scatter = (m_scatter && m_scatter->valid());
    const double per_gcm2 = PhysicalUnits::cm2 / PhysicalUnits::g;
    std::vector<T> scattered;

    const size_t num_sources = m_sources.size();
    result.summing_out.resize( num_sources );
    result.sum_pairs.resize( num_sources );

    // Per-line accumulators of (full-energy weight x factor) and (full-energy weight).
    std::vector<std::vector<T>> line_num( num_sources ), line_den( num_sources );
    // Per-pair accumulators over the nodes, of the node-dependent part (the escape) only:
    //  weight*a*b*survival, weight*a, and weight*b.
    struct PairSums{ T both = T(0.0), a = T(0.0), b = T(0.0); };
    std::vector<std::vector<std::vector<PairSums>>> pair_sums( num_sources );

    for( size_t src = 0; src < num_sources; ++src )
    {
      const Source &source = m_sources[src];
      line_num[src].assign( source.input.line_energies.size(), T(0.0) );
      line_den[src].assign( source.input.line_energies.size(), T(0.0) );
      pair_sums[src].resize( source.sum_peaks.size() );
      for( size_t p = 0; p < source.sum_peaks.size(); ++p )
        pair_sums[src][p].assign( source.sum_peaks[p].pairs.size(), PairSums() );
    }

    // `escape`: the fraction of a node's photons that leave the source unscattered; `tot_node`: the absolute
    //  total efficiency for a photon from the node.
    std::vector<T> escape( num_energies, T(1.0) ), tot_node( num_energies ), program_tot, program_out, pair_survival;
    ceelo::SummingOutProgram::Workspace<T> workspace;

    for( size_t node = 0; node < num_nodes; ++node )
    {
      const double weight = self_atten ? nodes.weight[node] : 1.0;

      T depth( 0.0 );  //areal density of the source in front of the node
      if( self_atten )
      {
        depth = nodes.fraction[node] * input.self_areal_density;
        for( size_t i = 0; i < num_energies; ++i )
          escape[i] = exp( -input.self_mu[i] * depth );
      }//if( self_atten )

      scattered.clear();
      const T areal_density = depth + input.external_areal_density;
      if( scatter && (scalar_value( areal_density ) > 0.0) )
      {
        const T atomic_number = (depth*input.self_atomic_number + input.external_an_areal_density) / areal_density;
        m_scatter->evaluateAll( atomic_number, T(areal_density * per_gcm2), scattered );
        assert( scattered.size() == num_energies );
      }

      for( size_t i = 0; i < num_energies; ++i )
      {
        T reaching = input.external_transmission[i] * escape[i];
        if( !scattered.empty() )
          reaching += scattered[i];
        tot_node[i] = tot_bare[i] * reaching;
      }

      for( size_t src = 0; src < num_sources; ++src )
      {
        const Source &source = m_sources[src];
        const ceelo::SummingOutProgram &program = *source.program;

        program_tot.assign( source.program_energy.size(), T(0.0) );
        for( size_t i = 0; i < source.program_energy.size(); ++i )
        {
          const int index = source.program_energy[i];
          if( index >= 0 )
            program_tot[i] = tot_node[index];
        }

        program_out.assign( program.num_lines(), T(1.0) );
        pair_survival.assign( program.num_pairs(), T(1.0) );
        program.evaluate( program_tot.data(), program_out.data(), workspace,
                          pair_survival.empty() ? nullptr : pair_survival.data() );

        const size_t num_lines = source.input.line_energies.size();
        for( size_t line = 0; line < num_lines; ++line )
        {
          const T fep_weight = weight * escape[source.program_line_energy[line]];
          line_num[src][line] += fep_weight * program_out[line];
          line_den[src][line] += fep_weight;
        }

        for( size_t p = 0; p < source.sum_peaks.size(); ++p )
        {
          const std::vector<SumPair> &pairs = source.sum_peaks[p].pairs;
          for( size_t q = 0; q < pairs.size(); ++q )
          {
            const SumPair &pair = pairs[q];
            const T &node_a = escape[pair.energy_a];
            const T &node_b = escape[pair.energy_b];

            // Coincident survival: exact where the program has the pair; else from each gamma member's
            //  summing-out with its partner's share removed.
            T survival( 1.0 );
            T surv_a( 1.0 ), surv_b( 1.0 );
            const bool estimate = (pair.program_pair < 0);
            if( !estimate )
              survival = pair_survival[pair.program_pair];

            if( estimate && (pair.program_line_a >= 0) )
            {
              T partner = 1.0 - pair.prob_b_given_a * tot_node[pair.energy_b];
              if( partner < 0.05 )
                partner = T(0.05);
              surv_a = program_out[pair.program_line_a] / partner;
              if( surv_a > 1.0 )
                surv_a = T(1.0);
            }
            if( estimate && (pair.program_line_b >= 0) )
            {
              T partner = 1.0 - pair.prob_a_given_b * tot_node[pair.energy_a];
              if( partner < 0.05 )
                partner = T(0.05);
              surv_b = program_out[pair.program_line_b] / partner;
              if( surv_b > 1.0 )
                surv_b = T(1.0);
            }

            if( estimate && (pair.program_line_a >= 0) && (pair.program_line_b >= 0) )
            {
              const T product = surv_a * surv_b;
              survival = (product > 1.0E-12) ? T(sqrt( product )) : T(0.0);
            }else if( estimate && (pair.program_line_a >= 0) )
            {
              survival = surv_a;
            }else if( estimate && (pair.program_line_b >= 0) )
            {
              survival = surv_b;
            }

            PairSums &sums = pair_sums[src][p][q];
            sums.both += weight * node_a * node_b * survival;
            sums.a += weight * node_a;
            sums.b += weight * node_b;
          }//for( loop over pairs )
        }//for( loop over sum peaks )
      }//for( loop over sources )
    }//for( loop over nodes )

    for( size_t src = 0; src < num_sources; ++src )
    {
      const Source &source = m_sources[src];
      const size_t num_lines = source.input.line_energies.size();
      result.summing_out[src].assign( num_lines, T(1.0) );
      for( size_t line = 0; line < num_lines; ++line )
      {
        if( line_den[src][line] > 0.0 )
          result.summing_out[src][line] = line_num[src][line] / line_den[src][line];
      }

      result.sum_pairs[src].resize( source.sum_peaks.size() );
      for( size_t p = 0; p < source.sum_peaks.size(); ++p )
      {
        const size_t num_pairs = source.sum_peaks[p].pairs.size();
        result.sum_pairs[src][p].assign( num_pairs, SumPairValue<T>() );
        for( size_t q = 0; q < num_pairs; ++q )
        {
          const PairSums &sums = pair_sums[src][p][q];
          SumPairValue<T> &value = result.sum_pairs[src][p][q];
          // both/a is the node part of b's efficiency, averaged over where a's detected counts come from.
          const SumPair &pair = source.sum_peaks[p].pairs[q];
          if( sums.a > 0.0 )
            value.b_given_a = fep_const[pair.energy_b] * sums.both / sums.a;
          if( sums.b > 0.0 )
            value.a_given_b = fep_const[pair.energy_a] * sums.both / sums.b;
        }
      }//for( loop over sum peaks )
    }//for( loop over sources )
  }//SummingModel::evaluate(...)
}//namespace RelActCalcCascade

#endif //RelActCalcCascade_imp_hpp
