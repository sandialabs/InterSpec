#include <cmath>
#include <vector>
#include <string>
#include <iomanip>
#include <float.h>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <assert.h>

#include "SandiaDecay.h"


//To compile on UNIX like systems:
//  g++ sandia_decay.cpp SandiaDecay.cpp -o sandia_decay.exe

using namespace std;

//Forward declarations
int scratch_work();
int exampleCoding();
int lookForIssues();
int testAllDecays();
int print_nuclides_with_gammas_in_energy_range();
int prinNuclideInfo( int argc, const char * argv[] );
int printNuclideDecayChain( const char *nuclideName );
void addElementWithNaturalAbundance( SandiaDecay::NuclideMixture &mixture,
                                    const SandiaDecay::Element *element,
                                    double parent_element_activity );
std::vector<const SandiaDecay::Nuclide *>
nuclidesWithGammaInRange( const float lowE, const float highE,
                         const std::vector<const SandiaDecay::Nuclide *> &candidates,
                         const bool allowaging );


int main( int argc, const char * argv[] )
{
  //If no arguments are passed in, the results of the exampleCoding()
  //  is printed out.
  //If a single argument is passed in, it is assumed to be an isotope, and
  //  information about that istotope is printed out
  //If the argument 'decay' followed by a isotope is entered in, then the decay
  //  chain of that isotope is printed out
  //If a single argument of either 'test' or 'issues' is used, then the relavent
  //  funtions are called

//  return print_nuclides_with_gammas_in_energy_range();

  if( argc == 1 )
    return exampleCoding();

  if( strcmp( argv[1], "test" ) == 0 )
    return testAllDecays();

  if( strcmp( argv[1], "issues" ) == 0 )
    return lookForIssues();

  if( (argc>2) && (strcmp( argv[1], "decay" )==0) )
    return printNuclideDecayChain( argv[2] );

  return prinNuclideInfo( argc, argv );
}//int main( int argc, const char * argv[] )





int exampleCoding()
{
  using namespace SandiaDecay;

  SandiaDecayDataBase database( "sandia.decay.xml" );

  //One way to retrieve Nuclide objects from the database (note the string
  //  can be formatted like U-238 92-U-238, uranium_238 238U, etc etc
  const Nuclide *u238 = database.nuclide( "U238" );
  const Nuclide *u235 = database.nuclide( "U235" );

  //a second way to get Nuclides from the database
  int atomicNum, massNum, iso;
  assert( u238 == database.nuclide( atomicNum = 92, massNum = 238, iso = 0 ) );


  //Print some information out about U238, using the overloaded operator<<
  cout << SandiaDecay::human_str_summary(*u238) << endl << endl;

  //Manually print out partial info on U238 to give idea of what can be accessed
  cout << u238->symbol << ": Atomic Number " << u238->atomicNumber
      << ", Mass Number " << u238->massNumber << ", ISO " << u238->isomerNumber
      << ", " << u238->atomicMass << " AMU, Half Life " << u238->halfLife
      << " seconds, has activity " << u238->activityPerGram()/curie
      << " currie/gram, with" << u238->decaysToChildren.size() << " decays, and "
      << u238->decaysFromParents.size() << " parents.";
  //Print out some information on the second transistion of U238
  const Transition *secondDecay = u238->decaysToChildren.at(1);
  cout << " The second transition of " << secondDecay->parent->symbol
       << " is to "
       << (secondDecay->child?secondDecay->child->symbol: string(" itself)"))
       << " through " << secondDecay->mode
       << " decay, with a branching ratio " << secondDecay->branchRatio
       << " and one of the particles emmited is a ";
  const RadParticle &firstProduct = secondDecay->products.at(0);
  cout << firstProduct.type << " with energy "
       << firstProduct.energy/SandiaDecay::keV
       << " keV and intensity " << firstProduct.intensity
       << endl << endl; //you can also access hinderance, forbiddenness, & logFT



  //Example of how to create a mixture of nuclides, decay them, and then
  //  retrieve what gammas there will be at a given time; you can also get
  //  the alphas and betas as well
  NuclideMixture mixture;
  const Nuclide *np235 = database.nuclide( "Np-235" );
  const Nuclide *co60 = database.nuclide( "cobalt-60" );
  mixture.addNuclideByActivity( np235, 0.001*SandiaDecay::curie );
  mixture.addNuclideByActivity( co60, 0.01*SandiaDecay::curie );
  const double time = 1.0*SandiaDecay::year;
  const vector<EnergyRatePair> gammas = mixture.gammas( time, SandiaDecay::NuclideMixture::OrderByEnergy, true );
  cout << "After " << time/SandiaDecay::year << " year of aging, the original "
       << "1 mCi of Np235 and 10 mCi of" << " Co60, have the gamma lines:\n";
  for( size_t i = 0; i < gammas.size(); ++i )
    cout << "\t" << gammas[i].energy / keV << " keV with intesnity "
         << gammas[i].numPerSecond << "/sec" << endl;



  //An example of a different way to decay nuclides, starting from a single
  //  nuclide - with this method of calculating, you would have to manually
  //  multiple out activity*branchingRatio*intentsity to find abundance of
  //  gammas/beta/alphas.
  vector<NuclideTimeEvolution> evolutions
                        = SandiaDecayDataBase::getTimeEvolution( u238, 1.0*Bq );
  cout << "Starting from " << u238->symbol << " we get:" << endl;
  for( size_t i = 0; i < evolutions.size(); ++i )
  {
    const NuclideTimeEvolution &evo = evolutions[i];
    const Nuclide *nuclide = evo.nuclide;
    const double time = 100000.0*SandiaDecay::month;
    const double activity = evo.activity( time );
    cout << "    " <<  nuclide->symbol << " : " <<  activity/SandiaDecay::Bq
         << " bq at " << time/SandiaDecay::year << " years" << endl;
  }//for( loop over decayed-to element, i )
  cout << endl << endl;



  //Example of how get the fraction of time one nuclde decays through or from
  //  another nuclide
  const Nuclide *hg206 = database.nuclide( "80 Hg 206" );
  cout << "The U238 decay chain proceeds through Hg206 "
       << 100.0 * u238->branchRatioToDecendant( hg206 )
       << " (" <<  100.0*hg206->branchRatioFromForebear( u238 ) << " )"
       << "% of the time" << endl << endl;


  //An example of how to get all nuclides that may have decayed into the nuclide
  //  of interest, in this case Pb206; Note pb206Forebearers is pre-sorted to be
  //  in roughly parent-->child order.  You can also call Nuclide::decendants()
  const Nuclide *pb206 = database.nuclide( "Pb206" );
  const vector<const Nuclide *> pb206Forebearers = pb206->forebearers();
  cout << "Pb206 can come from: \n\t";
  for( size_t desc = 0; desc < pb206Forebearers.size(); ++desc )
    cout << pb206Forebearers[desc]->symbol << ", ";
  cout << endl << endl;

  /*
  NuclideMixture::addElementUsingNaturalAbundance(...) not quite implemented yet
  //An example of creating an isotope mixture, using the natural abundance of
  //  isotopes (note, this method is beta)
  NuclideMixture naturalAbundanceMix;
  const Element *uranium = database.element( 92 );
  const double n_initial_atom = 10000.0;
  naturalAbundanceMix.addElementUsingNaturalAbundance( uranium, n_initial_atom );
  cout << "Natural uranium mixture looks like:" << endl
       << naturalAbundanceMix.info(0.0) << endl;
  */


  //A simple example of given an enrichment level (by mass), getting the
  //  relative amplitudes of the _activity_ for each of the elements
  NuclideMixture enrichMix;

  const Element *u = database.element( 92 );
  const std::vector<NuclideAbundancePair> &u_isotopes = u->isotopes;
  cout << "Natural " << u->name << " is composed of:" << endl;
  for( size_t i = 0; i < u_isotopes.size(); ++i )
      cout << "\t" << u_isotopes[i].abundance << " of " << u_isotopes[i].nuclide->symbol << endl;
  cout << endl;

  addElementWithNaturalAbundance( enrichMix, u, 1.0*curie );
  
//  const double total_mass = 1000.0;
//  enrichMix.addNuclideByAbundance( u235, total_mass*0.00711*u235->atomsPerGram() );
//  enrichMix.addNuclideByAbundance( u238, total_mass*0.99284*u238->atomsPerGram() );
//  cout << "ratio us activiuty(u235)/activity(u238)="
//       << enrichMix.activity(0, u235)/curie << "/" << enrichMix.activity(0, u238)/curie
//       << "=" << enrichMix.activity(0, u235)/enrichMix.activity(0, u238) << endl;

  cout << "Pure U238: " << u238->activityPerGram()/becquerel << " becquerel/g" << endl
       << "Pure U235: " << u235->activityPerGram()/becquerel << " becquerel/g" << endl;
  cout << "For natural uranium:" << endl;
  const double age = 0.0*year;
  for( size_t i = 0; i < u_isotopes.size(); ++i )
  {
    const Nuclide *nuc = u_isotopes[i].nuclide;
    cout << "\t" << nuc->symbol << " activity = "
         << enrichMix.activity( age, nuc )/curie << " curie" << endl;
  }//for( size_t i = 0; i < u_isotopes.size(); ++i )




  cout << "\tAnd weighs " << enrichMix.totalMassInGrams( age ) << " g" << endl;


  return EXIT_SUCCESS;
}//int exampleCoding()



int prinNuclideInfo( int argc, const char * argv[] )
{
  if( argc < 2 )
  {
    cerr << "Example Usage: ./" << argv[0] << " Co60" << endl;
    return EXIT_FAILURE;
  }//if( argc < 2 )

  using namespace SandiaDecay;
  SandiaDecayDataBase database( "sandia.decay.xml" );

  const Nuclide *nuclide = database.nuclide( argv[1] );

  if( !nuclide )
  {
    cerr << "Could not find isotope '" << argv[1] << "' in the database"
         << endl;
    return EXIT_FAILURE;
  }//if( !nuclide )

  cout << SandiaDecay::human_str_summary(*nuclide) << endl;

  return EXIT_SUCCESS;
}//int prinNuclideInfo(int argc, const char * argv[])



int printNuclideDecayChain( const char *nuclideName )
{
  using namespace SandiaDecay;

  SandiaDecayDataBase database( "sandia.decay.xml" );

  const Nuclide *nuclide = database.nuclide( nuclideName );

  if( !nuclide )
  {
    cerr << "Could not find isotope '" << nuclideName << "' in the database"
         << endl;
    return EXIT_FAILURE;
  }//if( !nuclide )


  cout << setw(6) << nuclide->symbol << " decays through (% of time):" << endl;
  const vector<const Nuclide *> decendants = nuclide->descendants();

  for( size_t i = 0; i < decendants.size(); ++i )
  {
    const double br = nuclide->branchRatioToDecendant( decendants[i] );
    cout << "    " << setw(6) << left << decendants[i]->symbol << " ("
         << fixed << setprecision(2) << setw(6) << right << 100.0*br
         << "%)" << endl;
  }//for( loop over decendants, i )


  return EXIT_SUCCESS;
}//int printNuclideDecayChain( const char *nuclideName );

typedef std::vector<const SandiaDecay::Nuclide *> NuclideVec;

int testAllDecays()
{
  //This function just makes sure that the decays for all nuclides can be
  //  evaluated; does not check to makesure answers are necassarily correct.
  try
  {
    using namespace SandiaDecay;
    SandiaDecayDataBase database( "sandia.decay.xml" );
   
    const NuclideVec &nuclidesVec = database.nuclides();

    size_t num_nuclieds_decayed = 0;
    NuclideVec::const_iterator iter;
    for( iter = nuclidesVec.begin(); iter != nuclidesVec.end(); ++iter )
    {
      const Nuclide *nuclide = *iter;
      const double halfLife = nuclide->halfLife;

      if( !(halfLife <= DBL_MAX && halfLife >= -DBL_MAX)  ) //isinf() dosnt seem to be available in Microsoft VC++
        continue;

      ++num_nuclieds_decayed;
      NuclideMixture mixture;
      mixture.addNuclideByActivity( nuclide, 1.0*SandiaDecay::curie );
      for( size_t i = 0; i < 100; ++i )
      {
        const double time = i*7.0*halfLife/100.0;

        const vector<NuclideActivityPair> activities = mixture.activity( time );
        const vector<EnergyRatePair> gammas = mixture.gammas( time, SandiaDecay::NuclideMixture::OrderByEnergy, true );
      }//for( test the decay for 7 halfLives )
    }//for( iterate over nuclides available, iter )

    cout << "All " << num_nuclieds_decayed << " radioactive nuclides have been "
         << " decayed and evaluated at 100 different times, for both daughter "
         << " products, and gammas produced" << endl;
  }catch( std::exception &e )
  {
    cerr << "Problem teting decays: " << e.what() << endl;
    return EXIT_FAILURE;
  }//try / catch

  return EXIT_SUCCESS;
}//int testAllDecays();



void addElementWithNaturalAbundance( SandiaDecay::NuclideMixture &mixture,
                                    const SandiaDecay::Element *element,
                                    double parent_element_activity )
{
  //XXX
  //  Below I am assuming sandia.decay.xml lists natural abundance by number
  //  of atoms... otherhise should use Nuclide::atomsPerGram()
  const vector<SandiaDecay::NuclideAbundancePair> &isotopes = element->isotopes;
  vector<double> iso_activity( isotopes.size() );
  
  double total_act_perM = 0.0;
  
  for( size_t i = 0; i < isotopes.size(); ++i )
  {
    const SandiaDecay::NuclideAbundancePair &iso = isotopes[i];
    iso_activity[i] =  iso.nuclide->numAtomsToActivity( iso.abundance*1.0E6 ); //1.0E6 chosen with no special insight
    total_act_perM += iso_activity[i];
  }//for( size_t i = 0; i < n_iso; ++i )
  
  const double multiple_needed = parent_element_activity / total_act_perM;
  
  for( size_t i = 0; i < isotopes.size(); ++i )
  {
    const SandiaDecay::NuclideAbundancePair &info = isotopes[i];
    mixture.addNuclideInSecularEquilibrium( info.nuclide, multiple_needed * iso_activity[i] );
  }//for( size_t i = 0; i < isotopes.size(); ++i )
}//void addElementWithNaturalAbundance(...)



int lookForIssues()
{
  //This is a scratchwork function used to look for issues in sandia.decay.xml
  //  or the SandiaDecay code.
  using namespace SandiaDecay;

  SandiaDecayDataBase database( "sandia.decay.xml" );

  const NuclideVec &nuclidesVec = database.nuclides();

  NuclideVec::const_iterator iter;
  for( iter = nuclidesVec.begin(); iter != nuclidesVec.end(); ++iter )
  {
    const Nuclide *nuclide = *iter;
    const double halfLife = nuclide->halfLife;

    if( !(halfLife <= DBL_MAX && halfLife >= -DBL_MAX) ) //isinf() dosnt seem to be available in Microsoft VC++
      continue;


    double totalBr = 0.0;
    for( size_t i = 0; i < nuclide->decaysToChildren.size(); ++i )
    {
      //      if( nuclide->decaysToChildren[i]->mode != kElectronCapture )
      //      if( nuclide->decaysToChildren[i]->mode != kSf )
      totalBr += nuclide->decaysToChildren[i]->branchRatio;
    }

    if( (fabs(1.0-totalBr)>0.001) && (totalBr!=0.0) && nuclide->halfLife>300.0 )
    {
      const Nuclide *meta1 = database.nuclide( nuclide->atomicNumber,
                                               nuclide->massNumber, 1 );
      const Nuclide *meta2 = database.nuclide( nuclide->atomicNumber,
                                               nuclide->massNumber, 2 );
      cout << nuclide->symbol << " has a BR=" << totalBr;

      if( meta1 ) cout << " has Meta-1";
      else cout << " doesnt have Meta-1";

      if( meta2 ) cout << " has Meta-2";
      else cout << " doesnt have Meta-2";

      cout << " halfLife=" << nuclide->halfLife;

      cout << endl;
    }

  }//for( iterate over nuclides available, iter )

  return EXIT_SUCCESS;
}//int lookForIssues();



int print_nuclides_with_gammas_in_energy_range()
{
  using namespace SandiaDecay;
  SandiaDecayDataBase database( "sandia.decay.xml" );

  const double lower_energy = 10.0*keV;
  const double upper_energy = 3000.0*keV;
  const bool allowaging = false;
  const vector<const Nuclide *> &all_nuclides = database.nuclides();

  vector<const Nuclide *> nuclides;
  nuclides = nuclidesWithGammaInRange( lower_energy, upper_energy,
                                       all_nuclides, allowaging );

  cout << "# Nuclides with gammas in the range " << lower_energy/keV << " keV"
       << " to " << upper_energy/keV << " keV." << "\r\n"
       << "# This file generated from sandia.decay.xml and the sandia_decay c++ interface"
       << " on 20121025." << "\r\n"
       << "# For each nuclide the gamma energy in keV, number per decay for that energy, and decay"
       << " mode is given." << "\r\n"
       << "# The sandia decay project is implemented by the Secondary Reachback "
       << "project " << "\r\n"
       << "# at Sandia National labs from the ENSDF dataset and GADRAS library." << "\r\n"
       << "# Results have not been thoroughly tested or validated." << "\r\n"
       << "# Errors and mistakes may be blamed on Will Johnson org. 08131" << "\r\n";


  for( size_t i = 0; i < nuclides.size(); ++i )
  {
    const Nuclide *nuc = nuclides[i];
    cout << "# ##############################################################################################" << "\r\n";
    cout << nuc->symbol << "\t" << nuc->halfLife << " seconds" << "\r\n";
    const vector<const Transition *> &transistions = nuc->decaysToChildren;

    for( size_t j = 0; j < transistions.size(); ++j )
    {
      const Transition *trans = transistions[j];
//      cout << "\tDecayMode=" << trans->mode << "\t" << "br=" << trans->branchRatio
//           << " DaughterNuclide=";
//     if( trans->child )
//       cout << trans->child->symbol << "\r\n";
//     else
//      cout << "FissionProducts" << "\r\n";
//    cout << "\tGammas for this decay mode:" << "\r\n";
      const vector<RadParticle> &products = trans->products;
      for( size_t k = 0; k < products.size(); ++k )
      {
        const RadParticle &part = products[k];
        if( part.type == SandiaDecay::GammaParticle )
        {
          cout << "\t" << (part.energy/keV)
               << "\t" << (trans->branchRatio * part.intensity)
               << "\t" << trans->mode
               << "\r\n";
        }
      }//for( size_t k = 0; k < products.size(); ++k )
    }//for( size_t j = 0; j < transistions.size(); ++j )
  }//for( size_t i = 0; i < nuclides.size(); ++i )

  cout << flush;

  return 1;
}//int print_nuclides_with_gammas_in_energy_range()

std::vector<const SandiaDecay::Nuclide *>
nuclidesWithGammaInRange( const float lowE, const float highE,
                          const std::vector<const SandiaDecay::Nuclide *> &candidates,
                          const bool allowaging )
{
  using namespace SandiaDecay;
  
  std::vector<const Nuclide *> answer;
  
  if( allowaging )
  {
    for( size_t i = 0; i < candidates.size(); ++i )
    {
      const Nuclide * const nuc = candidates[i];
      NuclideMixture mix;
      mix.addNuclideByActivity( nuc, 1.0E6 * becquerel );
      const vector<EnergyRatePair> gammas = mix.gammas( nuc->halfLife, NuclideMixture::OrderByAbundance, true );
      
      for( size_t i = 0; i < gammas.size(); ++i )
      {
        const EnergyRatePair &a = gammas[i];
        if( a.energy >= lowE && a.energy <= highE )
          if( find(answer.begin(), answer.end(), nuc) == answer.end() )
            answer.push_back( nuc );
      }//for( size_t i = 0; i < gammas.size(); ++i )
    }//for( size_t i = 0; i < candidates.size(); ++i )
  }else
  {
    for( size_t i = 0; i < candidates.size(); ++i )
    {
      const Nuclide * const nuc = candidates[i];
      const std::vector<const Transition *> &trans = nuc->decaysToChildren;
      for( size_t j = 0; j < trans.size(); ++j )
      {
        const Transition * const tran = trans[j];
        const std::vector<RadParticle> &particles = tran->products;
        for( size_t k = 0; k < particles.size(); ++k )
        {
          const RadParticle &particle = particles[k];
          if( particle.type != GammaParticle )
            continue;
          
          if( particle.energy >= lowE && particle.energy <= highE )
            if( find(answer.begin(),answer.end(),nuc) == answer.end() )
              answer.push_back( nuc );
        }//for( size_t k = 0; k < particles.size(); ++k )
      }//for( size_t j = 0; j < trans.size(); ++j )
    }//for( size_t i = 0; i < candidates.size(); ++i )
  }//if( allowaging ) / else
  
  return answer;
}//nuclidesWithGammaInRange(...)


#include <fstream>
int scratch_work()
{
  using namespace SandiaDecay;

  NuclideMixture mix, mix2, natMix, promptMix;
  SandiaDecayDataBase database( "sandia.decay.xml" );
  const Nuclide *nuc = database.nuclide( "U235" );

/*
  const SymbolToNuclideMap &nucmap = database.symbolToNuclideMap();

  for( SymbolToNuclideMap::const_iterator iter = nucmap.begin(); iter != nucmap.end(); ++iter )
  {
    const Nuclide *nuclide = iter->second;

    const std::vector<const Transition *> &transitions = nuclide->decaysToChildren;

    for( size_t i = 0; i < transitions.size(); ++i )
    {
      const Transition *trans = transitions[i];
      if( trans->mode == SandiaDecay::kBeta )
      {
        const std::vector<RadParticle> &products = trans->products;

        for( size_t j = 0; j < products.size(); ++j )
        {
          const RadParticle &particle = products[j];
          if( particle.energy > 2200  && (nuclide->halfLife > (4*3600)) )
          {
            cout << *nuclide << endl;
            i = transitions.size();
            break;
          }//if( particle.energy > 3000 )
        }//for( size_t j = 0; j < products.size(); ++j )
      }//if( trans->mode == SandiaDecay::kBeta )
    }//for( size_t i = 0; i < transitions.size(); ++i )
  }//for( SymbolToNuclideMap::iterator iter = nucmap.begin(); iter != nucmap.end(); ++iter )

  return 1;
*/

  promptMix.addNuclideInPromptEquilibrium( nuc, 1.0*becquerel );
//  promptMix.addNuclideByActivity( nuc, 1.0 );
  vector<EnergyRatePair> promptGammas = promptMix.gammas( 0.0*year, SandiaDecay::NuclideMixture::OrderByEnergy, true );
  for( size_t i = 0; i < promptGammas.size(); ++i )
  {
    if( promptGammas[i].numPerSecond > 0.00001 )
      cerr << promptGammas[i].energy << " keV --> " << 100.0*promptGammas[i].numPerSecond << "/s" << endl;
  }


  cerr << nuc->symbol << " secularEquilibriumHalfLife = " << nuc->secularEquilibriumHalfLife() << endl;
  cerr << nuc->symbol << " promptEquilibriumHalfLife = " << nuc->promptEquilibriumHalfLife() << endl;

  const Element *el = database.element( 92 );

  addElementWithNaturalAbundance( natMix, el, 1.0*becquerel );
  
  ofstream logfile( (nuc->symbol+"_secular_eq_log.txt").c_str() );
  nuc = database.nuclide( "U238" );
  mix.addNuclideInSecularEquilibrium( nuc, 1.0*becquerel );
  logfile << "By addNuclideInSecularEquilibrium(...):" << endl;
  logfile << "Total activity=" << mix.totalActivity(0.0)/becquerel << " becquerel" << endl;
  vector<EnergyRatePair> gammas = mix.gammas( 0.0*year, SandiaDecay::NuclideMixture::OrderByEnergy, true );
  for( size_t i = 0; i < gammas.size(); ++i )
    logfile << gammas[i].energy << " keV --> " << gammas[i].numPerSecond << "/s" << endl;

  mix2.addNuclideByActivity( nuc, pow(2.0,8.0) );
  logfile << endl << endl << "aging a long time yields:" << endl;
  logfile << "Total activity=" << mix2.totalActivity(8.0*nuc->halfLife)/becquerel << " becquerel" << endl;
  vector<EnergyRatePair> gammas2 = mix2.gammas( 8.0*nuc->halfLife, SandiaDecay::NuclideMixture::OrderByEnergy, true );
  for( size_t i = 0; i < gammas2.size(); ++i )
    logfile << gammas2[i].energy << " keV --> " << gammas2[i].numPerSecond << "/s" << endl;

  return EXIT_SUCCESS;
}//int scratch_work()

double det_response( double energy )
{
  static vector<double> energies, response;
  if( energies.empty() )
  {
    energies.push_back(10); response.push_back( 0.01*0.0000	    );
    energies.push_back(20	); response.push_back( 0.01*0.0000);
    energies.push_back(25	); response.push_back( 0.01*0.0402);
    energies.push_back(30	); response.push_back( 0.01*0.8865);
    energies.push_back(35	); response.push_back( 0.01*4.3059);
    energies.push_back(40	); response.push_back( 0.01*10.471);
    energies.push_back(45	); response.push_back( 0.01*17.919);
    energies.push_back(50	); response.push_back( 0.01*25.275);
    energies.push_back(55	); response.push_back( 0.01*31.811);
    energies.push_back(60	); response.push_back( 0.01*37.325);
    energies.push_back(70	); response.push_back( 0.01*45.530);
    energies.push_back(80	); response.push_back( 0.01*50.889);
    energies.push_back(90	); response.push_back( 0.01*54.204);
    energies.push_back(100	); response.push_back( 0.01*56.158);
    energies.push_back(120	); response.push_back( 0.01*57.227);
    energies.push_back(140	); response.push_back( 0.01*55.590);
    energies.push_back(170	); response.push_back( 0.01*51.082);
    energies.push_back(200	); response.push_back( 0.01*45.385);
    energies.push_back(250	); response.push_back( 0.01*35.845);
    energies.push_back(300	); response.push_back( 0.01*27.068);
    energies.push_back(350	); response.push_back( 0.01*23.910);
    energies.push_back(400	); response.push_back( 0.01*21.365);
    energies.push_back(450	); response.push_back( 0.01*19.188);
    energies.push_back(500	); response.push_back( 0.01*17.253);
    energies.push_back(550	); response.push_back( 0.01*16.496);
    energies.push_back(600	); response.push_back( 0.01*15.838);
    energies.push_back(650	); response.push_back( 0.01*15.250);
    energies.push_back(700	); response.push_back( 0.01*14.713);
    energies.push_back(750	); response.push_back( 0.01*13.892);
    energies.push_back(800	); response.push_back( 0.01*13.108);
    energies.push_back(850	); response.push_back( 0.01*12.357);
    energies.push_back(900	); response.push_back( 0.01*11.633);
    energies.push_back(950	); response.push_back( 0.01*11.136);
    energies.push_back(1000	); response.push_back( 0.01*10.657);
    energies.push_back(1100	); response.push_back( 0.01*9.7457);
    energies.push_back(1200	); response.push_back( 0.01*8.8890);
    energies.push_back(1300	); response.push_back( 0.01*8.4316);
    energies.push_back(1400	); response.push_back( 0.01*8.0034);
    energies.push_back(1500	); response.push_back( 0.01*7.6014);
    energies.push_back(1600	); response.push_back( 0.01*7.2255);
    energies.push_back(1700	); response.push_back( 0.01*6.9829);
    energies.push_back(1800	); response.push_back( 0.01*6.7548);
    energies.push_back(1900	); response.push_back( 0.01*6.5397);
    energies.push_back(2000	); response.push_back( 0.01*6.3365);
    energies.push_back(2200	); response.push_back( 0.01*6.0041);
    energies.push_back(2400	); response.push_back( 0.01*5.7056);
    energies.push_back(2600	); response.push_back( 0.01*5.4360);
    energies.push_back(2800	); response.push_back( 0.01*5.1911);
    energies.push_back(3000	); response.push_back( 0.01*4.9996);
    energies.push_back(3500	); response.push_back( 0.01*4.5891);
    energies.push_back(4000	); response.push_back( 0.01*4.2523);
    energies.push_back(4500	); response.push_back( 0.01*3.9668);
    energies.push_back(5000	); response.push_back( 0.01*3.7538);
    energies.push_back(6000	); response.push_back( 0.01*3.3865);
    energies.push_back(7000	); response.push_back( 0.01*3.0566);
    energies.push_back(8000	); response.push_back( 0.01*2.7352);
    energies.push_back(10000	); response.push_back( 0.01*2.7384);
    energies.push_back(12000	); response.push_back( 0.01*2.7081);
    energies.push_back(15000	); response.push_back( 0.01*2.6241);
    energies.push_back(20000	); response.push_back( 0.01*2.4366);
    energies.push_back(30000	); response.push_back( 0.01*2.5484);
    energies.push_back(40000	); response.push_back( 0.01*2.7018);
    energies.push_back(60000	); response.push_back( 0.01*3.0920);
    energies.push_back(100000); response.push_back( 0.01*4.0711 );
  }//if( energies.e,mpty() )


  vector<double>::const_iterator pos = lower_bound( energies.begin(), energies.end(), energy );

  if( pos == energies.end() )
    return 0.0;

  const size_t index = pos - energies.begin();
  return response[index];
}//double det_response( double energy )



