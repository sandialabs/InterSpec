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

// Checks peak files against InterSpec itself (linked to InterSpec's library; see run.sh, which sets
//  INTERSPEC_DATA_DIR to InterSpec's data directory):
//   interspec_check read <file>                   - loads with SpecMeas, prints the peaks InterSpec sees
//   interspec_check write-spe <in.n42> <out.spe> [samples] - writes the given (comma separated), or
//                                                   displayed, samples as SPE, with their peaks
//   interspec_check add-peaks <in.n42> <out.n42>   - adds reaction, element x-ray, S.E., D.E., and
//                                                   annihilation peaks to the displayed samples, and
//                                                   writes N42 (how tests/data/interspec_sources.n42 was made)
// Annihilation peaks are printed without their nuclear transition: InterSpec accepts them without
//  one, and the light page (which doesnt know the positron's decay) writes them that way.
#include <set>
#include <deque>
#include <memory>
#include <fstream>
#include <iostream>

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SandiaDecay/SandiaDecay.h"

using namespace std;

string src_str( const PeakDef &p )
{
  if( !p.hasSourceGammaAssigned() )
    return (p.parentNuclide() ? ("(nuc " + p.parentNuclide()->symbol + " w/o transition)") : string("-"));
  string s;
  if( p.parentNuclide() ) s = p.parentNuclide()->symbol;
  else if( p.xrayElement() ) s = p.xrayElement()->symbol + "-xray";
  else if( p.reaction() ) s = p.reaction()->name();
  const char *t = "";
  switch( p.sourceGammaType() )
  {
    case PeakDef::NormalGamma: break;
    case PeakDef::AnnihilationGamma: t = " annih"; break;
    case PeakDef::SingleEscapeGamma: t = " SE"; break;
    case PeakDef::DoubleEscapeGamma: t = " DE"; break;
    case PeakDef::XrayGamma: t = " xray"; break;
  }
  char buf[64];
  snprintf( buf, sizeof(buf), " %.2f", p.gammaParticleEnergy() );
  string trans;
  if( p.nuclearTransition() && (p.sourceGammaType() != PeakDef::AnnihilationGamma) )
    trans = string(" [") + p.nuclearTransition()->parent->symbol + "->" + (p.nuclearTransition()->child ? p.nuclearTransition()->child->symbol : string("?")) + "]";
  return s + buf + t + trans;
}

void print_peaks( const SpecMeas &meas )
{
  for( const string &w : meas.parse_warnings() )
    cout << "WARNING: " << w << endl;
  cout << "displayed:";
  for( int s : meas.displayedSampleNumbers() )
    cout << " " << s;
  cout << endl;
  for( const set<int> &samples : meas.sampleNumsWithPeaks() )
  {
    auto peaks = meas.peaks( samples );
    cout << "samples:";
    for( int s : samples )
      cout << " " << s;
    cout << " -> " << (peaks ? peaks->size() : 0) << " peaks" << endl;
    if( !peaks )
      continue;
    // Sorted by energy, as InterSpec shows them (the light page keeps them sorted)
    vector<shared_ptr<const PeakDef>> sorted( begin(*peaks), end(*peaks) );
    std::sort( begin(sorted), end(sorted), &PeakDef::lessThanByMeanShrdPtr );
    for( const auto &p : sorted )
      printf( "  %9.3f area=%.1f fwhm=%.3f %-34s cont=%s skew=%s cal=%d color=%s\n", p->mean(), p->peakArea(),
             p->fwhm(), src_str(*p).c_str(), PeakContinuum::offset_type_str(p->continuum()->type()),
             PeakDef::to_string(p->skewType()), int(p->useForEnergyCalibration()), p->lineColor().cssText().c_str() );
  }
}

int main( int argc, char **argv )
{
  const char * const data_dir = getenv( "INTERSPEC_DATA_DIR" );
  if( (argc < 3) || !data_dir )
  {
    cerr << "Usage: INTERSPEC_DATA_DIR=<InterSpec>/data interspec_check <command> <file> [...]" << endl;
    return 1;
  }
  InterSpec::setStaticDataDirectory( data_dir );
  const string cmd = argv[1];
  auto meas = make_shared<SpecMeas>();
  if( !meas->load_file( argv[2], SpecUtils::ParserType::Auto ) )
  {
    cerr << "load failed" << endl;
    return 1;
  }

  if( cmd == "read" )
  {
    print_peaks( *meas );
    return 0;
  }

  // Files written below record just the input file's name, not where it is on this computer
  meas->set_filename( SpecUtils::filename( meas->filename() ) );

  set<int> samples = meas->displayedSampleNumbers();
  if( samples.empty() )
    samples = meas->sample_numbers();
  if( (cmd == "write-spe") && (argc > 4) )
  {
    samples.clear();
    vector<int> nums;
    SpecUtils::split_to_ints( argv[4], strlen(argv[4]), nums );
    samples.insert( begin(nums), end(nums) );
  }

  if( cmd == "write-spe" )
  {
    ofstream out( argv[3], ios::binary );
    set<int> dets( begin(meas->detector_numbers()), end(meas->detector_numbers()) );
    return meas->write_iaea_spe( out, samples, dets ) ? 0 : 2;
  }

  if( cmd == "add-peaks" )
  {
    auto existing = meas->peaks( samples );
    deque<shared_ptr<const PeakDef>> peaks;
    if( existing )
      peaks = *existing;

    const auto add = [&peaks]( const double mean, const string &src_txt ) -> shared_ptr<PeakDef> {
      auto p = make_shared<PeakDef>( mean, 0.8, 1000.0 );
      p->continuum()->setRange( mean - 5, mean + 5 );
      p->continuum()->setType( PeakContinuum::Linear );
      p->continuum()->setParameters( mean, vector<double>{10.0, 0.0}, vector<double>{} );
      if( !src_txt.empty() )
      {
        const auto result = PeakModel::setNuclideXrayReaction( *p, src_txt, 4.0 );
        if( result == PeakModel::SetGammaSource::FailedSourceChange )
          cerr << "Failed to set source '" << src_txt << "'" << endl;
      }
      peaks.push_back( p );
      return p;
    };

    add( 846.77, "Fe(n,n) 846.77 keV" );
    add( 74.97, "Pb xray 74.97 keV" );
    add( 2103.53, "Th232 S.E. 2614.53 keV" );  //the photon energy, as InterSpec's text parsing takes
    add( 1592.53, "Th232 D.E. 2614.53 keV" );
    add( 2223.25, "H(n,g) 2223.25 keV" );

    // Annihilation from Na22's positrons
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    const SandiaDecay::Nuclide *na22 = db->nuclide( "Na22" );
    auto annih = add( 511.0, "" );
    for( const SandiaDecay::Transition *t : na22->decaysToChildren )
    {
      for( size_t i = 0; i < t->products.size(); ++i )
      {
        if( (t->products[i].type == SandiaDecay::PositronParticle) && !annih->hasSourceGammaAssigned() )
          annih->setNuclearTransition( na22, t, static_cast<int>(i), PeakDef::AnnihilationGamma );
      }
    }

    std::sort( begin(peaks), end(peaks), &PeakDef::lessThanByMeanShrdPtr );
    meas->setPeaks( peaks, samples );
    ofstream out( argv[3], ios::binary );
    return meas->write_2012_N42( out ) ? 0 : 2;
  }

  cerr << "Unknown command" << endl;
  return 1;
}
