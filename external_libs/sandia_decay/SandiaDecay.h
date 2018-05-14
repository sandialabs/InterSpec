#ifndef SandiaDecay_h
#define SandiaDecay_h
/* SandiaDecay: a library that provides nuclear decay info and calculations.
 
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

#include <vector>
#include <string>

//Forward Declaration
namespace rapidxml { template<class Ch> class xml_node; }


/** SandiaDecay is a library that provides nuclear decay info and calculations.
 An example use is:
 \code{.cpp}
 using namespace SandiaDecay;
 SandiaDecayDataBase database( "sandia.decay.xml" );
 
 const Nuclide *nuclide = database.nuclide( "U238" );
 //Strings equiv to "U238" are "U-238", "92-U-238", "uranium_238", "238U", etc
 //Or you can retrieve nuclides using atomic number, mass, and isometric level
 assert( nuclide == database.nuclide( 92, 238, 0 ) );

 cout << nuclide->symbol << "has Half Life " << u238->halfLife
      << " seconds, has activity " << nuclide->activityPerGram()/curie
      << " currie/gram, with" << nuclide->decaysToChildren.size() << " decays, "
      << " and " << nuclide->decaysFromParents.size() << " parents.";
 
 //Example of how to create a mixture of nuclides, decay them, and then
 //  retrieve what gammas there will be at a given time; you can also get
 //  the alphas and betas as well.  All daughter products included in correct
 //  proportions
 NuclideMixture mixture;
 const Nuclide *np235 = database.nuclide( "Np-235" );
 const Nuclide *co60 = database.nuclide( "cobalt-60" );
 mixture.addNuclideByActivity( co60, 0.01*SandiaDecay::curie );
 mixture.addAgedNuclideByActivity( np235, 0.001*SandiaDecay::curie, 20*SandiaDecay::year );
 
 const double time = 5.7*SandiaDecay::year;
 
 //Get gammas emmitted by the nucleuses during decay, with results sorted
 //  by energy; the 'true' argument says include gammas from positron
 //  annihilations as well.
 //  The Co60 (t_{1/2}=5.7y) emmisions will be that from 5.7-year old Co60 (i.e.,
 //  0.005 ci), but the emmisions from Np237 (t_{1/2}=1.08y) will be equivalent
 //  from a sample aged 25.7y where the Np235 isotope currenly has an activity
 //  of 26.15 uCi (1mCi aged 5.27 half lives).
 const vector<EnergyRatePair> gammas = mixture.gammas( time, NuclideMixture::kByEnergy, true );
 
 //Get the xrays produced during decays
 const vector<EnergyRatePair> xrays = mixture.xrays( time, NuclideMixture::kByEnergy );
 
 //Get all xrays, gammas, and anhilation gammas, this time sorted by abundance
 const vector<EnergyRatePair> photons = mixture.photons( time, NuclideMixture::kByAbundance );
 \endcode
 */


//Known Issues:
// -Convertions between number of atoms and activity (happens durring
//  decay calculation, and other places) might cause numerical errors for
//  large activities or elements with large half lifes; superficially this
//  appears to not be an issue, but should be explicitly checked at some point.
// -The XML input file, sandia.decay.xml is known to have issues
// -Coincident gamma rays (like in Co60) are not indicated as such.
// -If there are multiple daughters in a decay chain with the same half-life
//  decay calculations will fail.
// -SPontaneous fission products are not included in this library.
// -A generall code clean up is needed

//On 20140426 I converted a number member variables of ints to short ints, and
//  doubles to floats, however I did not systematically go through and make the
//  useage of floats and doubles consistent, or even check that the correct
//  results are still gotten 

//A few design decitions:
// -Code should compile using any C++03 compliant compiler
// -A unit system based on seconds, keV, and becquerel has been used, with
//  converstions to other units provided in this namespace SandiaDecay
// -All necessary functions are included in one header and coresponding cpp
//  file; this is an attempt to make it easier to include in projects, and
//  allow the programmer to see all available functions in one file.
// -To avoid unneccassary copying of nuclide and transition object, const
//  pointers to these objects will be used by SandiaDecayDataBase; additionally
//  the user should never delete these objects, as they are always owened by
//  SandiaDecayDataBase (and in fact can only be created by SandiaDecayDataBase)
// -Lots of struct's are defined to aid in reading the code (as oposed to using
//  something like std::pair<>).  It is intended the user of this library should
//  access the member variables of structs directly. For the structs, I try to
//  follow the convention exampled by: for the struct NuclideAbundancePair,
//  there is a member variable name nuclide, and a member variable named
//  abundance.  Member variables should not be directly accesed for class object
//  but for structs they should be
// -Multithreaded initialization (parsing of xml) was tried, but added no
//  performance gain, so was removed.


namespace SandiaDecay
{
  /** An enumeration of possible decay products. */
  enum ProductType
  {
    BetaParticle,
    GammaParticle,
    AlphaParticle,
    PositronParticle,
    CaptureElectronParticle,  //neutrino
    XrayParticle
  };//enum ProductType

  
  /** An enumeration of possible decay modes. */
  enum DecayMode
  {
    AlphaDecay,
    BetaDecay,
    BetaPlusDecay,                    /**< Positron Emission */
    ProtonDecay,
    IsometricTransitionDecay,         /**< Isomeric Transition */
    BetaAndNeutronDecay,              //TODO (20151218): These look to be a subset of beta (e.g. beta br includes this, this could be a issue as BR should always sum to 1.0)
    BetaAndTwoNeutronDecay,
    ElectronCaptureDecay,
    ElectronCaptureAndProtonDecay,
    ElectronCaptureAndAlphaDecay,
    ElectronCaptureAndTwoProtonDecay,
    BetaAndAlphaDecay,
    BetaPlusAndProtonDecay,
    BetaPlusAndTwoProtonDecay,
    BetaPlusAndThreeProtonDecay,
    BetaPlusAndAlphaDecay,
    DoubleBetaDecay,
    DoubleElectronCaptureDecay,
    Carbon14Decay,
    SpontaneousFissionDecay,
    ClusterDecay,
    //"2b+"  //Double Positron Emission
    //"n" Neutron Emission
    UndefinedDecay
  };//enum DecayMode

  /** An enumeration of possible forbiddeness of decays.
      Applies to beta, positron, and capture electron decays.
   */
  enum ForbiddennessType
  {
    NoForbiddenness,
    FirstForbidden,
    FirstUniqueForbidden,
    SecondForbidden,
    SecondUniqueForbidden,
    ThirdForbidden,
    ThirdUniqueForbidden,
    FourthForbidden
  };//enum ForbiddennessType

  
  //Some forward declarations
  struct Element;
  struct Nuclide;
  struct Transition;
  struct RadParticle;
  struct EnergyRatePair;
  struct EnergyIntensityPair;
  struct NuclideAbundancePair;
  struct NuclideNumAtomsPair;
  struct NuclideActivityPair;
  struct NuclideTimeEvolution;

  
  /** Unit system used by this library, examples of use are:
  \code{.cpp}
     const double time = 10.0*day;
     cout << "time=" << time/second << " seconds\n";
     cout << "U238 HalfLife=" << db.nuclide("U238")->halfLife/year << "y\n";
  \endcode
   */
  static const double second      = 1.0;
  static const double hour        = 3600 * second;
  static const double day         = 24.0 * hour;
  static const double year        = 365.2425 * day;
  static const double month       = year / 12.0;

  static const double becquerel   = 1.0;
  static const double curie       = 3.7E10 * becquerel;
  static const double Bq          = becquerel;
  static const double MBq         = 1.0E6 * becquerel;
  static const double Ci          = curie;

  static const double keV         = 1.0;
  static const double eV          = keV / 1000.0;
  static const double MeV         = 1000.0 * keV;

  static const double mm          = 1.0;
  static const double cm          = 10.0*mm;
  static const double meter       = 1000.0*mm;
  static const double m           = meter;
  static const double cm2         = cm*cm;
  static const double cm3         = cm*cm*cm;

//A potential problem: to properly define mass, we whould do the following:
//  static const double e_SI   = 1.602176487e-19;// positron charge in coulomb
//  static const double joule  = eV/e_SI;// joule = 6.24150 e+12 * MeV
//  static const double kilogram = joule*second*second/m2, kg = joule*second*second/m2;
//  static const double g          = 1.0e-3*kilogram, gram = 1.0e-3*kilogram;
//However, it seems more natural to just define the gram as the base unit:
//  static const double g          = 1.0;
//So... just skipping defining mass right now and instead will explicitly add
//  mass units of variables or member functions in their respective names.

  
  /** Class that parses the input "sandia.decay.xml" file and "owns" the decay
   information.  Once initialized, all member functions (besides reset()) are
   const, so you should be able to use the same object from multiple threads.
   */
  class SandiaDecayDataBase
  {
  public:
    /** Constructs, and initializes the object, by parsing the specified
        "sandia.decay.xml".
        Will throw exception on error.
     */
    SandiaDecayDataBase( const std::string &sandia_decay_xml_filename );
    
    /** Default constructor: you must call initialize(string) before using the
        database.  Useful for global static variables where the name of the
        "sandia.decay.xml" file can be decided at runtime.
     */
    SandiaDecayDataBase();

    /** Destructor. */
    ~SandiaDecayDataBase();

    /** Returns true if either the SandiaDecayDataBase(string) contructor
        was used, or #initialize(string) has been called.
        I.e. If "sandia.decay.xml" has been parsed yet.
     */
    bool initialized() const;
    
    /** Parses the specified "sandia.decay.xml" file.
        Will throw exception on error.
     */
    void initialize( const std::string &sandia_decay_xml_filename );
    
    /** Discards the parsed in nuclear decay information from memory.  You must
        call #initialize() again before further using.
     */
    void reset();
    
    /** Check if the XML file contained decay x-ray information (e.g., the
        x-rays that are given off during nuclear decays).
     */
    bool xmlContainedDecayXRayInfo() const;
    
    /** Check if the XML file contained elemental x-ray information (e.g.,
       xrays that are caused by flouresence)
     */
    bool xmlContainedElementalXRayInfo() const;
    
    /** Provides an _estimate_ of the memory taken up by this object, not what
        is actually allocated for it
     */
    size_t memsize() const;
    
    /** Function to retrieve a nuclide specified by a string such as: "Co60",
        "Uranium 238", "Ho166m", "Th-232", "201 Tl", etc.
        If nuclide is not found, will return a null pointer.
     */
    //TODO: Implementation could use optimization
    const Nuclide *nuclide( const std::string &label ) const;
    
    /** Function to retrieve a nuclide specified atomic, mass, and isomeric
        numbers. If nuclide is not found, will return a null pointer.
     */
    //TODO: Implementation could use optimization
    const Nuclide *nuclide( int z, int massNumber, int iso = 0 ) const;

    /** Retireves nuclide information from a string that you know is in a
        conforming format:
        The nuclide symbol should be an uppercase letter followed by lower cases
        if it is longer than 1 letter (eg IUPAC names).
        The mass number should not be padded in any way, and come directly after
        the nuclide symbol
        The meta-stable state should be indicated by a lower case 'm' followed
        by a 2 or 3 if a second or third state (do not indicate a 1)
        Some examples: "Co60" "U258" "U238m" "Zn73m2"
    
        Will return null if Nuclide is not found.
        When in doubt use the #nuclide(string) function; this function is
        provided as a slightly more efficient option when you control the format
        of the input.
    */
    const Nuclide *nuclideUnchecked( const std::string &symbol ) const;


    /** Get the nuclides corresponding to a given element */
    //  ToDO - implement test.
    std::vector<const Nuclide *> nuclides( const Element *element ) const;
    
    /** Get the nuclides corresponding to a given element */
    //  ToDO - implement test.
    std::vector<const Nuclide *> nuclides( const std::string &element ) const;
 
    /** Get Element object for a given atomic number.
        Returns NULL on failure.
     */
    const Element *element( const int atomicNumber ) const;
    
    /** Get the Element object from either the elemental symbol (ex., U, Pu, Ba)
        or the name (uranium, plutonium, barium).
     
        Ignores non alphabetical characters, but not inteligently. E.g., Co60,
        co, cobalt will all work, but not Co60m
     
        Returns a null pointer on failure.
     */
    //ToDo: Currently very inefficient with no tests setup.
    const Element *element( const std::string &label ) const;
    
    /** Trys to turn common inputs, such as 92-U-238, uranium238, 238U, 238m,
        meta238, cobalt60, In 112-m3, etc into the normalized forms of, for
        instance, U238, U238m, Co60, In112m3.  This function throws an exception
        if it cant identify a valid nuclide.
     */
    //ToDo: This function could also use some serious optimizations and tests
    //      implemented.
    std::string toNormalSymbolForm( const std::string &symbol ) const;
    
    /** Returns the underlying store of Nuclide objects.
        Currently sorted alphabetically by the nuclides symbols, but not
        garunteed to stay that way in the future.
     */
    const std::vector<const Nuclide *> &nuclides() const;
    
    /** Returns the underlying store of Element objects.
        Currently ordered by atomic number, but could change.
     */
    const std::vector<const Element *> &elements() const;
    
    /** Returns the underlying store of Transition objects. */
    const std::vector<Transition> &transitions() const;

    
    static
    std::vector<NuclideActivityPair> decay( const Nuclide *parent,
                                            const double orignal_activity,
                                            const double time_in_seconds );
    static
    std::vector<NuclideActivityPair> decay(
                               const std::vector<NuclideNumAtomsPair> &parents,
                               const double time );
    static
    std::vector<NuclideActivityPair> decay(
                               const std::vector<NuclideActivityPair> &parents,
                               const double time );


    static
    std::vector<NuclideTimeEvolution> getTimeEvolution( const Nuclide *parent,
                                              const double orignal_activity );
    static
    std::vector<NuclideTimeEvolution> getTimeEvolution(
                              const std::vector<NuclideNumAtomsPair> &parents );
    static
    std::vector<NuclideTimeEvolution> getTimeEvolution(
                              const std::vector<NuclideActivityPair> &parents );

  private:
    /** Performs the actual work of parsing sandia.decay.xml */
    static void parseSandiaDecayXml( const std::string &filename,
                                    std::vector<const Nuclide *> &sorted_nucs,
                                    std::vector<Nuclide> &stored_nuclides,
                                    std::vector<const Element *> &elements,
                                    std::vector<Element> &stored_elements,
                                    std::vector<Transition> &transitionStore );
    
    /** Helper function that sets m_xmlFileContainedDecayXrays and
        m_xmlFileContainedElementalXrays.
     */
    void checkIfContainedXrays();
    
    /** Sorted alphabetically by their symbol */
    std::vector<const Nuclide *> m_nuclides;

    /** sorted by atomic number */
    std::vector<const Element *> m_elements;
  
    /** Storage of Element, Nuclide, and Transition objects.  The idea is to avoid
       individual allocations to keep things closer in memorry, and hence
       maybye a little bit faster (or even be able to better take advantage
       of the memory getting mmap'ed for faster swapping).
       There doesnt seem to be a big effect in decay_database_benchmark(...),
       but keeping anyway since there should be some advantage in principle.
    Note that moving to not having m_nuclides, m_elements doesnt actually add
       any time to the initialization of the database (~1ms or less), and
       it probably removes ~60 or 70 kb of memorry (this wasnt tested), but it
       has the effect of no longer ensuring each isotope has a single unique
       address, because a copy contructor becomes necassarry to make available
       for std::sort; perhaps there is some way around this (E.g. using c++11
       or specializing std::swap), but I'll leave it off to the future.
     */
    std::vector<Element> m_elementStore;
    std::vector<Nuclide> m_nuclideStore;
    std::vector<Transition> m_transitionStore;
    
    //
    bool m_xmlFileContainedDecayXrays;
    bool m_xmlFileContainedElementalXrays;
  };//class SandiaDecayDataBase



  /** Class to facilitate (potentially) mixing nuclides, and the aging them and
      getting the resulting activities of the parent and prodigy nuclides, as
      well as the decay particles produced.
   
      A single instance of this class is not safe to to use from multiple
      threads (e.g., dont share an instance across threads - make a new
      NuclideMixture instance for each thread, which is safe to do).
   */
  class NuclideMixture
  {
    //Some benchmarks for decaying isotopes via the NuclideMixture class:
    //  To decay all 3331 radioactive nuclides, and then evaluate their activity
    //  for 100 differnt times each, it only takes about 0.8 seconds total.
    //  If gammas produced are requested in addition to the activity, then the
    //  time it takes jumps up to 22 second - meaning evaluating gammas produced
    //  could use some optimizations (see code for break down of what is taking
    //  the most time).

  public:
    NuclideMixture();
    virtual ~NuclideMixture();

    /** Clear all the nuclides you have added to this mixture. */
    void clear();

    /** Add a nuclide by specifying how many nuclide atomis are initially in the
        mixture.
     */
    void addNuclideByAbundance( const Nuclide *nuclide, double num_init_atom );
    
    /** Add a nuclide by specifying the initial activity of the nuclide.
        Activity should be in SandiaDecay units.  Ex.
     \code{.cpp}
     const SandiaDecay::Nuclide *u238 = db.nuclide("U238");
     mixture.addNuclideByActivity( u238, 0.001*SandiaDecay::curie );
     \endcode
     */
    void addNuclideByActivity( const Nuclide *nuclide, double activity );
    
    /** Convenience function to add nuclide by specifying number of initial
        atoms.
     */
    void addNuclide( const NuclideNumAtomsPair &nucide );
    
    /** Convenience function to add nuclide by specifying initial activity. */
    void addNuclide( const NuclideActivityPair &nucide );

    /** Add a nuclide to the mixture that is already pre-aged.  The activity
        cooresponds to the parent nuclides activity at the mixtures t=0 age.
        For example, if you add 1uCi of U232 (t_{1/2}=68.8y) with an inital age
        of 20 years, and then ask for the gammas at time 68.8y, you will get the
        gammas a 88.8y old sample with a current U232 activity of 0.5uCi; if you
        asked for gammas at a time of 0y, you would get the gammas of a 20 year
        old sample that has an activity of 1uCi.  Ex.
     \code{.cpp}
     const SandiaDecay::Nuclide *u238 = db.nuclide("U238");
     mixture.addAgedNuclideByActivity( u238, 0.001*SandiaDecay::curie, 20.0*SandiaDecay::year );
     \endcode
     */
    //ToDo: implement test.
    void addAgedNuclideByActivity( const Nuclide *nuclide,
                                   double activity,
                                   double age_in_seconds );

    /** Add a nuclide to the mixture that has laready obtained secular
        equilibrium.  The activity specified is of the parent nuclide a the
        mixtures t=0.
     
       If the nuclide can not obtain secular equilibrium, then an error will
       be printed to stderr, the nuclide will not be added to the mixture,
       and function will return false.
       Returns true if the nuclide can obtain secular equilibrium and was
       added to mixture.
     */
    bool addNuclideInSecularEquilibrium( const Nuclide *parent,
                                         double parent_activity );

    
    /** Adds the daughter nuclides of 'parent' of whose half lives are
        monotonically decreasing, some examples are:
    
        Parent    What is Added to Mixture
        Th232     Th232, Ra228, Ac228
        U234      U234, Th230, Ra226, Rn222, Po218, At218, Rn218, Po214
        U235      U235, Th231
        U238      U238, Th234m, Pa234m
     
        The parent nuclide is always added in (unless its stable).
     */
    void addNuclideInPromptEquilibrium( const Nuclide *parent,
                                       double parent_activity );
    
    /** Returns the number of nuclides at mixture t=0. */
    int numInitialNuclides() const;
    
    /** Returns the number of nuclides given in the solution, and cooresponds
        to the index for nulides in the vector returned by
        decayedToNuclidesEvolutions().
     */
    int numSolutionNuclides() const;

    /** Get the activity of a specific nuclide in the mixture (either a parent
        or prodigeny nuclide) at a specific time after mixture t=0.
        Returned value is in units of SandiaDecay (e.x., divide result by
        #SandiaDecay::curie to determine number of Curies).
     */
    double activity( double time, const Nuclide *nuclide )    const;
    double activity( double time, const std::string &symbol ) const;
    double activity( double time, int z, int atomic_mass, int iso )  const;
    double activity( double time, int internal_index_number ) const;
    std::vector<NuclideActivityPair> activity( double time_in_seconds ) const;

    /** Get the number of atoms of a specific nuclide remaining at a given time
        after mixture.
     */
    double numAtoms( double time, const Nuclide *nuclide ) const;
    double numAtoms( double time, const std::string &symbol ) const;
    double numAtoms( double time, int z, int atomic_mass, int iso )  const;
    double numAtoms( double time, int internal_index_number ) const;
    std::vector<NuclideNumAtomsPair> numAtoms( double time_in_seconds )  const;

    /** Get the total summed activity for all nuclide in the mixture, including
        the prodigeny, at a given time.
        There are probably not to many circumstances when this function is
        useful, so make sure its what you actualy want.
     */
    double totalActivity( double time ) const;
    
    /** Get the total mass in grams of the mixture. */
    double totalMassInGrams( double time = 0.0 ) const;  //in grams

    /** How to order returned decay particles. */
    enum HowToOrder
    {
      OrderByAbundance,
      OrderByEnergy
    };//enum HowToOrder

    /** Get the rate (how many per second) of each energy of gammas at a given
        time after mixture t=0.
        If annihilation gamas is specified as true, two 511 keV gammas will
        be included for each positron emission.
        Xrays are not included in the result.
     */
    std::vector<EnergyRatePair> gammas( const double time,
                                             const HowToOrder sortType,
                                             const bool includeAnihilationGammas
                                            ) const;
    
    /** Get the rate (how many per second) of each energy of alphas at a given
        time after mixture t=0.
     */
    std::vector<EnergyRatePair> alphas( double time,
                                             HowToOrder sortType = OrderByAbundance
                                            ) const;
    
    /** Get the rate (how many per second) of each energy of betas at a given
        time after mixture t=0.
     */
    std::vector<EnergyRatePair> betas( double time,
                                            HowToOrder sortType = OrderByAbundance
                                           ) const;
    
    /** Get the rate (how many per second) of each energy of b+ at a given
        time after mixture t=0.
     */
    std::vector<EnergyRatePair> betaPlusses( double time,
                                              HowToOrder sortType = OrderByAbundance
                                                 ) const;
    
    /** Get the rate (how many per second) of each energy of x-rays at a given
        time after mixture t=0.
     
        Results will be empty if XML file didnt contain decay xrays.
        \sa SandiaDecayDataBase::xmlContainedDecayXRayInfo()
     */
    std::vector<EnergyRatePair> xrays( double time,
                                              HowToOrder sortType = OrderByAbundance
                                                 ) const;
    
    /** Returns rate of gammas, anihilation, and xrays for a given mixture age.
        X-rays only included if the XML file contained decay xrays.
     */
    std::vector<EnergyRatePair> photons( const double time,
                                              const HowToOrder sortType = OrderByAbundance
                                             ) const;
    
    /** The generic function to get rate of a specfici type of decay particle
        at a given mixture age.
     */
    std::vector<EnergyRatePair> decayParticle( double time,
                                               ProductType type,
                                               HowToOrder sortType = OrderByAbundance
                                                   ) const;
    
    /** Return a human readable summary of mixture for a given time. */
    std::string info( const double time ) const;
    

    //internal_index_number runs from 0 to numSolutionNuclides()
    int internalIndexNumber( const Nuclide *nuclide ) const;
    int internalIndexNumber( const std::string &symbol ) const;
    int internalIndexNumber( int z, int atomic_mass, int iso ) const;

    const Nuclide *initialNuclide( int index ) const;
    
    /** Internal index number spans from 0 to #numSolutionNuclides()
     */
    const Nuclide *solutionNuclide( int internal_index_number ) const;

    const std::vector<NuclideTimeEvolution> &decayedToNuclidesEvolutions() const;

  protected:
    void performTimeEvolution() const;

    std::vector<NuclideNumAtomsPair> m_originalNuclides;
    
    /** TODO: get rid of the mutable-ness of this variable, and instead do one
              of the following:
              1) Introduce a set of const member functions that require you
                 already called #performTimeEvolution
              2) Explicitly require the user to call #performTimeEvolution
                 before accessing any decay information
              3) Perform the decay everytime a nuclide is added, and keep the
                 current member functions marked as const, const
              (I originally marked this mutable to optimize a specific situation
               but it was probably not a good idea)
     */
    mutable std::vector<NuclideTimeEvolution> m_decayedToNuclides;
  };//class NuclideMixture


  /** Struct to store information about a nuclide. */
  //ToDo: check if the doubles and ints could be reduced in size with out
  //  effecting accuracy or performance
  struct Nuclide
  {
    /** The normalized ascii string symbol for this nuclide.
        Ex. "U238", "Pu237m", "Co60", "Au192m2", etc.
     */
    std::string symbol;
    
    /** The atomic number for this nuclide (number of protons in nucleus). */
    short int atomicNumber;
    /** The atomic mass number (number of protons+neutrons) */
    short int massNumber;
    /** Nuclear excitation state (isomer number) */
    short int isomerNumber;
    /** Atomic mass in AMU (e.g., in units where C12==12.0 AMU) */
    float atomicMass;
    /** Nuclide half-life in units of SandiaDecay::seconds. */
    //leaving halfLife as double for now (float~7 sig fig, double~17 sig fig)
    double halfLife;
    
    /** The nuclear transitions this nuclide decays through. */
    std::vector<const Transition *> decaysToChildren;
    /** The nuclear transitions that this nuclide can be the result of. */
    std::vector<const Transition *> decaysFromParents;

    /** The decay constant that is defined as
        0.5 = exp( -decay_const*halfLife ), or put another way ln(0.5)/halfLife.
     */
    double decayConstant() const;

    /** Number of transitions (decay channels) for this nuclide. */
    size_t numDecays() const;
    /** The decay constant for a specific transition (decay channel). */
    double partialDecayConstant( const size_t transition_num ) const;
    /** Fraction of decays that proceed by a certain transition */
    float branchingRatio( const size_t transition_num ) const;
    /** The child nuclide that is the result of a decay transition. */
    const Nuclide *child( const size_t transition_num ) const;

    /** Convert from the number of atomis of this nuclide, to its activity.
        Retuned valus is in SandiaDecay units, ex.:
        \code{.cpp}
        const SandiaDecay::Nuclide *u238 = db.nuclide("U238");
        cout << u238->numAtomsToActivity(6.022E23)/SandiaDecay::curie << "ci\n";
        \endcode
     */
    double numAtomsToActivity( const double num_atoms ) const;
    /** Convert from the activity of this nuclide, to the number of atoms.
     \code{.cpp}
     const SandiaDecay::Nuclide *u238 = db.nuclide("U238");
     cout << u238->activityToNumAtoms(0.001*SandiaDecay::curie)/ << "atoms\n";
     \endcode
     */
    double activityToNumAtoms( const double activity ) const;

    /** The fraction of decays of this nuclide that proceeds through the
        specified descendant (which may be multiple generations down the chain).
     */
    //TODO - Barely tested
    float branchRatioToDecendant( const Nuclide *descendant ) const;
    
    /** The fraction of decays of the specified ancestor (which can be mutliple
        generations above) that proceeds through this nuclide.
     */
    //TODO - Barely tested
    float branchRatioFromForebear( const Nuclide *ancestor ) const;

    /** Returns all prodigeny (descendant) isotopes (daughter, granddaughter,
        etc.) and not just immediate daughter nuclides
        Results will also include this nuclide.
     */
    std::vector<const Nuclide *> descendants() const;

    /** Returns all isotopes where this nuclide will be in their decay chain
        (E.g., parent, grandparent, great-grandparent, etc).
        Results will also include this nuclide.
     */
    std::vector<const Nuclide *> forebearers() const;

    /** Maximum half-life of all nuclide descendants; dictates time-scale for
        secular equilibrium. If this value exceeds or equals nuclide half-life,
        secular equilibrium cannot be achieved, however this function will still
        return this value.
    */
    double secularEquilibriumHalfLife() const;
    
    /** Returns true if all decendants have a shorter half life than this
        nuclide.
     */
    bool canObtainSecularEquilibrium() const;

    /** Prompt equilibrium half life is maximum half-life of nuclide descendants
        in decay series with monotonically decreasing half-lives; dictates scale
        for prompt equilibrium.  If prompt equilibrium half-life exceeds or
        equals nuclide half-life, prompt equilibrium cannot be achieved so zero
        is returned.
        For instance, an analysis applicatin may use this to define photopeaks
        for U, Pu, Th, etc, by aging a Nuclide by
        log(1000,2)*promtpEquilibriumHalfLife().
     
        If a promtpEquilibriumHalfLife() doesnt exist, then the secular
        equilibrium is used, and if this doesnt exist, only for the isotope.
     */
    double promptEquilibriumHalfLife() const;
    
    /** Returns true if the prompt equilibrium half life is defined. */
    bool canObtainPromptEquilibrium() const;

    /** The number of atoms of this nuclide it would take to equal 1 gram. */
    double atomsPerGram() const;
    
    /** The activity (in SandiaDecay units) of one gram of this nuclide. */
    double activityPerGram() const;

    /** Determine if all child isotopes are stable (e.g. aging doesnt effect
        gamma energies emmitted)
     */
    bool decaysToStableChildren() const;

    /** Provides an _estimate_ of the memory taken up by this object, not what
        is actually allocated for it.  For development purposes only.
     */
    size_t memsize() const;

    //operator< compares how we intuatively want: if lhs is a daughter of rhs
    //  then returns true, else compares based on mass, atomic, & isomer numbers
    bool operator<( const Nuclide &rhs ) const;
    bool operator==( const Nuclide &rhs ) const;
    bool operator!=( const Nuclide &rhs ) const;
    static bool lessThan( const Nuclide *lhs, const Nuclide *rhs );
    static bool greaterThan( const Nuclide *lhs, const Nuclide *rhs );
    
    ~Nuclide();
    
  protected:
    /** Constructor/destructor protected since Nuclides created are dependant
        on other Nuclides (via the Transitions pointed to by multiple Nuclides),
        and thus must be created and deleted together to avoid memory
        inconsistencies.
     */
    Nuclide( const ::rapidxml::xml_node<char> *node );

    /** Used to set information from a XML node. */
    void set( const ::rapidxml::xml_node<char> *node );


#if __cplusplus > 199711L
    //operator= declared but not implemented in order to disallow copying, thus
    //  preventing memorry leaks and pointer inconsistencies
    //For c++ pre11 we cant do this and still store m_nuclideStore in
    //  SandiaDecayDatabase...
    Nuclide &operator=( const Nuclide &rhs );
#endif

    friend class SandiaDecayDataBase;
  };//class Nuclide


  /** Information representing the nuclear transition (decay channel) for
      a nuclide.
   */
  struct Transition
  {
    /** The parent nuclide of the decay.  Will never be null when initialized
        as part of a SandiaDecayDataBase.
     */
    const Nuclide *parent;
    
    /** The resultant nuclide after the decay.
        May be null, for example for spontaneous fission decays.
     */
    const Nuclide *child;
    
    /** Decay mode represented by this transition object. */
    DecayMode mode;
    
    /** Fraction of times the parent nuclide decays through this decay channel.
     */
    float branchRatio;
    
    /** The particles (gamma, beta, alphas, etc) emmitted when this transition
        takes place.
     */
    std::vector<RadParticle> products;

    /** Provides an _estimate_ of the memory taken up by this object, not what
        is actually allocated for it.  For development purposes only.
     */
    size_t memsize() const;
    
    Transition(){}
    Transition( const ::rapidxml::xml_node<char> *node,
                const std::vector<const Nuclide *> &nuclides );
    void set( const ::rapidxml::xml_node<char> *node,
               const std::vector<const Nuclide *> &nuclides );

    friend class SandiaDecayDataBase;
  };//struct transition

  
  
  /** Represents information of a particle (of a given energy) given off during
      a nuclear transition.
  */
  struct RadParticle
  {
    /** Particle type. */
    ProductType type;
    
    /** Energy of the particle.  Applies to all ProductTypes */
    float energy;
    
    /** Intensity of this particle, for this decay channel.  E.g., between 0
        and 1.  Applies to all ProductTypes
     */
    float intensity;
    
    /** Hindrance.  Applies only to alpha decays. */
    float hindrance;
    
    /** is log-10 of fermi integral (F)*parent_half-life (T)
        Applies to beta, positron, and electronCapture decays
    */
    float logFT;
    
    /** Forbiddenss of the decay. Applies to beta, positron, and electron
        capture decays
     */
    ForbiddennessType forbiddenness;
    

    /** Revision status in the source XML data for this particle.  Can be
        ignored by users of this library.
     */
    enum RevisionType
    {
      NoRevision,
      EditedRevision,
      InsertedRevision,
      DeletedRevision
    };//enum RevisionType
    
    RevisionType revision;

    /** Other radiation particles expected to be detected when this particle
        is detected.  Useful for gamma spectroscopy where you have a chance
        of detecting two gammas as one detection event (thus you detect the
        summed energy).
        Currently this information is woefully incomplete, not well tested,
        only available for a handfull of nuclides, and only includes gammas;
        the mechanism of indicating coincident particles will be changed to
        something better in the future.
     */
    //20120804 wcjohns: I think the below should be replaced with something
    //   like a link to the other coincident gammas,
    //   e.g. struct CoincidentLink{ size_t index; float intensity; }
    //Where the index is to Transition::products of the parent of the RadParticle
    std::vector<RadParticle> coincidentRadiation; // mostly for gammas

    /** Provides an _estimate_ of the memory taken up by this object, not what
        is actually allocated for it.  For development purposes only.
     */
    size_t memsize() const;
    
  protected:
    RadParticle( const ::rapidxml::xml_node<char> *node );
    RadParticle( ProductType, float energy, float intensity );
    void set( const ::rapidxml::xml_node<char> *node );

    friend struct Transition;
  };//struct RadParticle


  /** Represents information for an element */
  struct Element
  {
    std::string symbol;
    std::string name;
    short int atomicNumber;

    /** The isotopes which make up the natural abundance of an element
        The abundance is the fractional mass of each isotope (eg add up to 1.0),
        it may be the case that abundances are all exactly 0.0 if the element
        isont naturally occuring
     */
    std::vector<NuclideAbundancePair> isotopes;

    /** Xrays that are caused by flouresence (e.g. not xrays produced as a
        result of a decay) exciting the element; intensities are relative to the
        most intense (1.0) xray, and are only approximations to be used.
        Energies are only approximate as well.  Intended as a rough refence to
        be used when xrays are wanted without a decay.
        May not be present if this data wasnt in the XML data file.
        \sa SandiaDecayDataBase::xmlContainedElementalXRayInfo
     */
    std::vector<EnergyIntensityPair> xrays;


    /** Computed from mass average of isotopes.
        returns in units of g/mole (e.g. for Iron, returns "55.8451")
        For non-naturally occuring elements, the atomic mass of the longest
        lived isotopes is returned
     */
    double atomicMass() const;

    Element(){}
    Element( ::rapidxml::xml_node<char> *node,
            const std::vector<const Nuclide *> &nuclides );
    void set( ::rapidxml::xml_node<char> *node,
              const std::vector<const Nuclide *> &nuclides );

    /** Provides an _estimate_ of the memory taken up by this object, not what
        is actually allocated for it.  For development purposes only.
     */
    size_t memsize() const;
    
    static bool lessThan( const Element *lhs, const Element *rhs );
    static bool lessThanForAtomicNum( const Element *lhs, const int atomicNum );
  };//class Element


  /** Structure to return the rate of a specific-energied decay particle.
      E.g., Give the rate for a certain energy gamma.
   */
  struct EnergyRatePair
  {
    /** Energy, in units of SandiaDecay. */
    double energy;
    
    /** The number of particles given off per second. */
    double numPerSecond;
    
    EnergyRatePair( double _numPerSecond, double _energy );
    static bool moreThanByNumPerSecond( const EnergyRatePair &lhs,
                                       const EnergyRatePair &rhs );
    static bool lessThanByEnergy( const EnergyRatePair &lhs,
                                 const EnergyRatePair &rhs );
  };//struct AbundanceEnergyPair

  
  /** Structure to return the relative (to the number of decays) intensities of
      a specific-energied decay particles.
      E.g., Specify what fraction of decay event will have a gamma of a certain
      energy.
   */
  struct EnergyIntensityPair
  {
    /** Energy, in units of SandiaDecay. */
    double energy;
    
    /** The fraction of decays with a particle of this energy.
        Typically 0 to 1.0, but for things like gammas caused by positron
        annihilation, this number can be larger than 1.
     */
    double intensity;
    
    EnergyIntensityPair( double e, double i ) : energy(e), intensity(i) {}
  };//struct EnergyIntensityPair

  
  /**  Structure to hold the information about the abundance (by mass) of a
       specific nuclide in a mixture, material, or element.
   */
  struct NuclideAbundancePair
  {
    /** Nuclide this information is for. */
    const Nuclide *nuclide;
    
    /** The fraction of this nuclide, by mass, in the mixture/material/element
     */
    double abundance;
  };//struct NuclideAbundancePair


  /**  Structure to hold the information about the number of atoms of a
       specific nuclide in a mixture, material, or element.
   */
  struct NuclideNumAtomsPair
  {
    /** Nuclide this information is for. */
    const Nuclide *nuclide;
    
    /** The number of atoms present. */
    double numAtoms;
    
    NuclideNumAtomsPair( const Nuclide *_nuclide, double _numAtoms );
  };//struct NuclideNumAtomsPair


  /**  Structure to hold the information about the activity of a specific
       nuclide in a mixture, material, or element.
   */
  struct NuclideActivityPair
  {
    /** Nuclide this information is for. */
    const Nuclide *nuclide;
    
    /** Activity of the nuclide, in units of SandiaDecay. */
    double activity;
    
    NuclideActivityPair( const Nuclide *_nuclide, double _activity );
  };//struct NuclideActivityPair

  
  /** Information to help calculate decays.
      Probably not useful other than in the impementation of this library,
      unless you are really optimizing cpu-efficiency of calculations.
   */
  struct TimeEvolutionTerm
  {
    double termCoeff;
    double exponentialCoeff;
    
    double eval( double time_in_seconds ) const;
    TimeEvolutionTerm( double mag, const Nuclide *nuc );
  };//struct TimeEvolutionTerm

  
  /** This class describes the time dependant magnitude of a single nuclide.
      Probably not useful other than in the impementation of this library,
      unless you are really optimizing cpu-efficiency of calculations.
   */
  struct NuclideTimeEvolution
  {
    const Nuclide *nuclide;
    std::vector<TimeEvolutionTerm> evolutionTerms;
    NuclideTimeEvolution( const Nuclide *_nuclide );

    double numAtoms( double time_in_seconds ) const;
    double activity( double time_in_seconds ) const;
    void addEvolutionTerm( double mag,  const Nuclide *nuc );
  };//struct TimeEvolution

  
  //some input/output operators
  const char *to_str( ProductType s );
  const char *to_str( const DecayMode &s );
  
  //Some methods to output a human readable summarry of more complex objects
  std::string human_str_summary( const Nuclide &nuclide );
  std::string human_str_summary( const Transition &transition );
  std::string human_str_summary( const RadParticle &particle );
}//namespace SandiaDecay


#endif //SandiaDecay_h
