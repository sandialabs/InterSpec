#include <set>
#include <cmath>
#include <cfloat>
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <sstream>
#include <iostream>

#include <sqlite3.h>

#include "SandiaDecay.h"

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"


using namespace std;

/** This program parses sandia.decay.xml, and removes any xray information
    present.  It then inserts the decay xray informaiton from tori.db3 into the
    xml document under each nuclides transition element.  It also parses
    sandia.xray.xml, and inserts its information under the <element> node of
    each element.
    Assumes sandia.decay.xml, tori.db3, and sandia.xray.xml are all in CWD, and
    saves sandia.decay.xray.xml to CWD.
 */


int test_sandia_decay_xrays();
int add_xray_info_to_sandia_decay_xml();
vector< pair<float,float> > test_gadras_xray_calc_for_nuc( const SandiaDecay::Nuclide *nuclide );

const std::string tori_input_filename = "old_tori.db3";
const std::string sandia_decay_only_input_filename = "sandia.decay.xml";
const std::string sandia_xray_only_input_filename = "sandia.xray.xml";
const std::string sandia_decay_xray_filename = "old_sandia.decay.xray.xml";
//const std::string sandia_decay_xray_filename = "sandia.decay.xray.xml";  //new
/* Status: 20180507:
 */


namespace
{
  //adapted from http://stackoverflow.com/questions/3152241/case-insensitive-stdstring-find
  template<typename charT>
  struct char_iequal
  {
    char_iequal( const std::locale &loc ) : m_loc(loc) {}
    bool operator()(charT ch1, charT ch2) {
      return std::toupper(ch1, m_loc) == std::toupper(ch2, m_loc);
    }
  private:
    const std::locale &m_loc;
  };
  
  bool icontains( const char *line, const size_t length,
                 const char *label, const size_t labellen )
  {
    const char *start = line;
    const char *end = start + length;
    const char *it = std::search( start, end, label, label+labellen,
                                 char_iequal<char>(std::locale()) );
    const bool answer = (it != end);
    
    return answer;
  }//icontains(...)
  
  bool icontains( const std::string &line, const char *label )
  {
    const size_t labellen = strlen(label);
    return icontains( line.c_str(), line.size(), label, labellen );
  }//icontains(...)
  
  void ireplace_all( std::string &input, const char *pattern, const char *replacement )
  {
    //This function does not handle UTF8!
    if( input.empty() )
      return;
    
    const size_t paternlen = strlen(pattern);
    if( !paternlen )
      return;
    
    //If we are replacing "XX" with "X" in the sequence "XXX" we want the
    //  result to be "X", however we have to protect against the case where
    //  the replacement string contains the search string
    const bool replace_contains_pattern = icontains( replacement, pattern );
    
    const size_t replacment_len = strlen(replacement);
    bool found = true;
    const char *start = input.c_str(), *end = input.c_str() + input.size();
    
    while( found )
    {
      const char * const it = std::search( start, end, pattern, pattern+paternlen,
                                          char_iequal<char>(std::locale()) );
      found = (it != end);
      if( found )
      {
        size_t delstart = it - input.c_str();
        input.erase( delstart, paternlen );
        input.insert( delstart, replacement );
        start = input.c_str() + delstart + (replace_contains_pattern ? replacment_len : size_t(0));
        end = input.c_str() + input.size();
      }//if( found )
    }//while( found )
  }//void ireplace_all(...)
  
  SandiaDecay::DecayMode decaymode_fromstr( const std::string &decayMode )
  {
    //Coppied from SandiaDecay.cpp
    if( decayMode == "b-" )        return SandiaDecay::BetaDecay;
    else if( decayMode == "it" )   return SandiaDecay::IsometricTransitionDecay;
    else if( decayMode == "b-n" )  return SandiaDecay::BetaAndNeutronDecay;
    else if( decayMode == "b+" )   return SandiaDecay::BetaPlusDecay;
    else if( decayMode == "ec" )   return SandiaDecay::ElectronCaptureDecay;
    else if( decayMode == "b-a" )  return SandiaDecay::BetaAndAlphaDecay;
    else if( decayMode == "b-2n" ) return SandiaDecay::BetaAndTwoNeutronDecay;
    else if( decayMode == "ecp" )  return SandiaDecay::ElectronCaptureAndProtonDecay;
    else if( decayMode == "eca" )  return SandiaDecay::ElectronCaptureAndAlphaDecay;
    else if( decayMode == "b+p" )  return SandiaDecay::BetaPlusAndProtonDecay;
    else if( decayMode == "ec2p" ) return SandiaDecay::ElectronCaptureAndTwoProtonDecay;
    else if( decayMode == "b+2p" ) return SandiaDecay::BetaPlusAndTwoProtonDecay;
    else if( decayMode == "b+3p" ) return SandiaDecay::BetaPlusAndThreeProtonDecay;
    else if( decayMode == "b+a" )  return SandiaDecay::BetaPlusAndAlphaDecay;
    else if( decayMode == "2b-" )  return SandiaDecay::DoubleBetaDecay;
    else if( decayMode == "2ec" )  return SandiaDecay::DoubleElectronCaptureDecay;
    else if( decayMode == "a" )    return SandiaDecay::AlphaDecay;
    else if( decayMode == "p" )    return SandiaDecay::ProtonDecay;
    else if( decayMode == "14c" )  return SandiaDecay::Carbon14Decay;
    else if( decayMode == "sf" )   return SandiaDecay::SpontaneousFissionDecay;
    else if( decayMode == "Undefined" ) return SandiaDecay::UndefinedDecay;
    return SandiaDecay::UndefinedDecay;
  }//DecayMode decaymode_fromstr( const char *str )
}

int main( int argc, char **argv )
{
  //add_xray_info_to_sandia_decay_xml();

  test_sandia_decay_xrays();

  return EXIT_SUCCESS;
}//int main( int argc, char **argv )


struct XrayInfo
{
  float energy, intensity;
  string assignment;
  float l1intensity, l2intensity, l3intensity;
  bool operator<( const XrayInfo &rhs ) const { return energy < rhs.energy;}
};

static int gadras_xray_info_sqlite3_callback( void *data, int argc, char **argv, char **azColName )
{
  if( argc < 6 || !argv[0] || !argv[1] )
  {
    cout << "argc=" << argc << endl;
    return 1;
  }

  XrayInfo info;
  info.l1intensity = info.l2intensity = info.l3intensity = 0.0f;

  if( !(stringstream(argv[0]) >> info.energy) || !(stringstream(argv[1]) >> info.intensity) )
  {
    cerr << "Failed to convert energy or intensity" << endl;
    return 1;
  }

  info.intensity /= 100.0f;
  if( argv[2] )
    info.assignment = argv[2];

  ireplace_all( info.assignment, "<sub><font face=\"symbol\">", "" );
  ireplace_all( info.assignment, "</font>", "" );
  ireplace_all( info.assignment, "</sub>", "" );
  ireplace_all( info.assignment, "<sub>", "" );
  ireplace_all( info.assignment, "<i>", "" );
  ireplace_all( info.assignment, "</i>", "" );

  if( argv[3] && strlen(argv[3]) > 0 )
  {
    if( (stringstream(argv[3]) >> info.l1intensity) )
    {
      info.l1intensity /= 100.f;
    }else
    {
      assert(0);
    }
  }

  if( argv[4] && strlen(argv[4]) > 0 )
  {
    if( (stringstream(argv[4]) >> info.l2intensity) )
    {
      info.l2intensity /= 100.f;
    }else
    {
        assert(0);
    }
  }

  if( argv[5] && strlen(argv[5]) > 0 )
  {
    if( (stringstream(argv[5]) >> info.l3intensity) )
    {
      info.l3intensity /= 100.f;
    }else
    {
      assert(0);
    }
  }

  if( info.intensity > 0.0f )
  {
    vector<XrayInfo> *result = (vector<XrayInfo> *)data;
    result->push_back( info );
  }

  return 0;
}//gadras_xray_info_sqlite3_callback(...)


static int tori_xray_sqlite3_callback( void *data, int argc, char **argv, char **azColName )
{
  if( argc < 2 || !argv[0] || !argv[1] )
    return 1;

  float energy, intensity;
  if( !(stringstream(argv[0]) >> energy) || !(stringstream(argv[1]) >> intensity) )
    return 1;

  std::vector< pair<float,float> > *result = (vector< pair<float,float> > *)data;
  result->push_back( make_pair(energy,intensity/100.0f) );

//  cerr << "\t{" << energy << ", " << intensity/100.0f << "}" << endl;

  return 0;
}

bool less_than_first( const std::pair<float,float> &lhs, const float rhs )
{
  return lhs.first < rhs;
}

int test_gadras_xray_calc()
{
  SandiaDecay::SandiaDecayDataBase database;
  database.initialize( sandia_decay_only_input_filename );
  
  for( const SandiaDecay::Nuclide *nuclide : database.nuclides() )
    test_gadras_xray_calc_for_nuc( nuclide );
  
  return EXIT_SUCCESS;
}//int test_gadras_xray_calc()

int add_xray_info_to_sandia_decay_xml()
{
  SandiaDecay::SandiaDecayDataBase database;
  database.initialize( sandia_decay_only_input_filename );
  const std::vector<SandiaDecay::Transition> &transitions = database.transitions();
  
  rapidxml::xml_document<char> doc;
  {
    rapidxml::file<char> xmlfile( sandia_decay_only_input_filename.c_str() );
    doc.parse<rapidxml::parse_full>( xmlfile.data() );
  }
  
  rapidxml::xml_node<char> *doc_node = doc.first_node();
  
  if( doc_node->first_attribute("version") )
    doc_node = doc_node->next_sibling();

  {//begin codeblock to add sandia.xray.xml
    rapidxml::file<char> xrayfile( sandia_xray_only_input_filename.c_str() );
    rapidxml::xml_document<char> xraydoc;
    xraydoc.parse<rapidxml::parse_full | rapidxml::parse_trim_whitespace>( xrayfile.data() );
    rapidxml::xml_node<char> *xray_base_node = xraydoc.first_node("document");
    assert( xray_base_node );
    xray_base_node = xray_base_node->first_node("xrays");
    assert( xray_base_node );
    
    //Elements that sandia.xray.xml does not have data for.
    const string known_no_data[] = { "V", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
      "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut",
      "Uuq", "Uup", "Uuh", "Uus", "Uuo"
    };
    
    
    //Now lets add elemental xrays to the XML file
    for( const SandiaDecay::Element *el : database.elements() )
    {
      rapidxml::xml_node<char> *xmlel = doc_node->first_node( "element" );
      
      bool found = false;
      do
      {
        rapidxml::xml_attribute<char> *symbol = xmlel->first_attribute( "symbol" );
        found = (symbol && (symbol->value() == el->symbol));
        if( !found )
          xmlel = xmlel->next_sibling( "element" );
      }while( xmlel && !found );
      
      assert( xmlel && found );
      
      while( rapidxml::xml_node<char> *xrayel = xmlel->first_node("xray") )
        xmlel->remove_node( xrayel );
     
      if( std::find(begin(known_no_data),end(known_no_data),el->symbol) != end(known_no_data) )
        continue;
      
      vector<double> yields;
      vector<double> energies;
      rapidxml::xml_node<char> *xrayxmlel = NULL;
      
      for( rapidxml::xml_node<char> *node = xray_base_node->first_node( "xray" );
          node;
          node = node->next_sibling("xray") )
      {
        rapidxml::xml_node<char> *n = node->first_node( "nuclide" );
        assert( n );
        if( !rapidxml::internal::compare(n->value(), n->value_size(), el->symbol.c_str(), el->symbol.size(), false ) )
          continue;
        
        xrayxmlel = n;
        for( rapidxml::xml_node<char> *en = node->first_node("energies")->first_node("energy");
            en;
            en = en->next_sibling("energy") )
        {
          double eergy;
          string val( en->value(), en->value()+en->value_size() );
          if( !(stringstream(val) >> eergy) )
            throw runtime_error( "Failed to convers yeild to double" );
          energies.push_back( eergy );
        }
        
        for( rapidxml::xml_node<char> *en = node->first_node("yields")->first_node("yield");
            en;
            en = en->next_sibling("yield") )
        {
          double yield;
          string val( en->value(), en->value()+en->value_size() );
          if( !(stringstream(val) >> yield) )
            throw runtime_error( "Failed to convers yeild to double" );
          yields.push_back( 0.01*yield );
        }
      }
      
      if( !xrayxmlel )
      {
        cerr << "Did not find '" << el->symbol << "'!" <<endl;
        continue;
      }
      
      assert( xrayxmlel );
      assert( energies.size() == yields.size() );
      
      for( size_t i = 0; i < energies.size(); ++i )
      {
        //I dont know why sandia.xray.xml has xrays with energy 1.0 that clearly
        //  arent real xrays, so lets skip em.
        if( energies[i] < 1.0001 )
          continue;
        
        rapidxml::xml_node<char> *node = doc.allocate_node( rapidxml::node_element, "xray" );
        xmlel->append_node( node );
        
        char buffer[512];
        snprintf( buffer, sizeof(buffer), "%G", energies[i] );
        const char *value = doc.allocate_string( buffer );
        node->append_attribute( doc.allocate_attribute( "energy", value ) );
        
        snprintf( buffer, sizeof(buffer), "%G", yields[i] );
        value = doc.allocate_string( buffer );
        node->append_attribute( doc.allocate_attribute( "relintensity", value ) );
      }//for( size_t i = 0; i < energies.size(); ++i )
      
    }//foreach( const SandiaDecay::Element *el, database->elements() )

  }//end codeblock to add sandia.xray.xml
  
  // SQLite3 Tori database
  sqlite3 *toriDatabase = NULL;

  if( sqlite3_open( tori_input_filename.c_str(), &toriDatabase) != 0 )
  {
    cerr << "Failed to open: " << tori_input_filename << endl;
    return EXIT_FAILURE;
  }

  
  for( const SandiaDecay::Transition &trans : transitions )
  {
    if( !trans.child )
      continue;

    const int parentZ = trans.parent->atomicNumber;
    const int parentA = trans.parent->massNumber;
    const int parentM = trans.parent->isomerNumber;
    const int daughterZ = trans.child->atomicNumber;

    bool found = false;
    rapidxml::xml_node<char> *xmltrans = doc_node->first_node("transition");

    do
    {
      try
      {
        rapidxml::xml_attribute<char> *parent = xmltrans->first_attribute( "parent" );
        rapidxml::xml_attribute<char> *child = xmltrans->first_attribute( "child" );
        rapidxml::xml_attribute<char> *mode = xmltrans->first_attribute( "mode" );
        rapidxml::xml_attribute<char> *branchRatio = xmltrans->first_attribute( "branchRatio" );

        if( !parent || !child || !mode || !branchRatio )
          throw runtime_error("");

        if( parent->value() != trans.parent->symbol )
          throw runtime_error("");

        if( child->value() != trans.child->symbol )
          throw runtime_error("");

        const SandiaDecay::DecayMode modetype = decaymode_fromstr( mode->value() );
        if( modetype != trans.mode )
          throw runtime_error("");

        found = true;
      }catch(...)
      {
        xmltrans = xmltrans->next_sibling("transition");
      }
    }while( xmltrans && !found );

    if( !found )
      xmltrans = 0;

    if( !xmltrans )
      cerr << "Failed to find xmltrans" << endl;

    if( xmltrans )
      while( rapidxml::xml_node<char> *xrayel = xmltrans->first_node("xray") )
        xmltrans->remove_node( xrayel );
    
    std::vector<XrayInfo> thisresult;


    char query_buffer[512];
    snprintf( query_buffer, sizeof(query_buffer),
             "SELECT XRays.Energy, XIntensities.Int, XRays.Assignment, XRays.L1Int, XRays.L2Int, XRays.L3Int "
             "FROM ((XRays INNER JOIN XIntensities ON XRays.XCode = XIntensities.XCode) "
             "INNER JOIN parents ON XIntensities.iZA = parents.iZA) "
             "WHERE (parents.Z = %i) AND (parents.A = %d) AND ((parents.iZA - (10000 * parents.Z + parents.A)) / 300 = %d) AND ((XRays.XCode / 100) = %d) AND (XRays.Energy >= 10)",
             parentZ, parentA, parentM, daughterZ );

    char *errorMessage = NULL;
    void *callbackArg = &thisresult;
    int returnCode = sqlite3_exec( toriDatabase, query_buffer, &gadras_xray_info_sqlite3_callback, callbackArg, &errorMessage );

    if( errorMessage || returnCode )
      cerr << "SQL error: " << errorMessage << endl;

    std::sort( thisresult.begin(), thisresult.end() );

    std::set<std::string> xraystrs;

    if( xmltrans )
    {
      for( size_t j = 0; j < thisresult.size(); ++j )
      {
        rapidxml::xml_node<char> *node = doc.allocate_node( rapidxml::node_element, "xray" );
        xmltrans->append_node( node );

        char buffer[512];
        snprintf( buffer, sizeof(buffer), "%G", thisresult[j].energy );
        const char *value = doc.allocate_string( buffer );
        node->append_attribute( doc.allocate_attribute( "energy", value ) );


        snprintf( buffer, sizeof(buffer), "%G", thisresult[j].intensity );
        value = doc.allocate_string( buffer );
        node->append_attribute( doc.allocate_attribute( "intensity", value ) );

        value = doc.allocate_string( thisresult[j].assignment.c_str() );
        node->append_attribute( doc.allocate_attribute( "assignment", value ) );

        if( thisresult[j].l1intensity > 0.0f )
        {
          snprintf( buffer, sizeof(buffer), "%G", thisresult[j].l1intensity );
          value = doc.allocate_string( buffer );
          node->append_attribute( doc.allocate_attribute( "l1intensity", value ) );
        }

        if( thisresult[j].l2intensity > 0.0f )
        {
          snprintf( buffer, sizeof(buffer), "%G", thisresult[j].l2intensity );
          value = doc.allocate_string( buffer );
          node->append_attribute( doc.allocate_attribute( "l2intensity", value ) );
        }

        if( thisresult[j].l3intensity > 0.0f )
        {
          snprintf( buffer, sizeof(buffer), "%G", thisresult[j].l3intensity );
          value = doc.allocate_string( buffer );
          node->append_attribute( doc.allocate_attribute( "l3intensity", value ) );
        }

        stringstream sxrastrstrm;
        rapidxml::internal::print_node( std::ostream_iterator<char>(sxrastrstrm), node, 0, 0);

        if( xraystrs.count( sxrastrstrm.str() ) )
          cout << "Duplicate xrayel: " << sxrastrstrm.str() << endl;
        xraystrs.insert( sxrastrstrm.str() );
      }//for( size_t j = 0; j < thisresult.size(); ++j )
    }//if( xmltrans )

//    cout << trans.parent->symbol << " --> " << trans.child->symbol << ": ";
//    for( size_t j = 0; j < thisresult.size(); ++j )
//      cout << "{" << thisresult[j].energy << ", " << thisresult[j].intensity << "}, ";
//    cout << endl;
  }//foreach( const SandiaDecay::Transition &trans, transitions )


  sqlite3_close( toriDatabase );
  toriDatabase = NULL;

  
  const char *val = 0;
  rapidxml::xml_node<char> *node = 0;
  val = "LBNL ToRI database";
  node = doc.allocate_node( rapidxml::node_element, "reference", val );
  node->append_attribute( doc.allocate_attribute( "type", "decay xray" ) );
  doc_node->prepend_node( node );
  
  val = "sandia.xray.xml";
  node = doc.allocate_node( rapidxml::node_element, "reference", val );
  node->append_attribute( doc.allocate_attribute( "type", "fluorescence xray" ) );
  doc_node->prepend_node( node );
  
  //val = "Nuclide and decay data parsed from ENSDF using DHS TRB funded code";
  //node = doc.allocate_node( rapidxml::node_comment, "", val );
  //doc_node->prepend_node( node );

  ofstream outputfile( sandia_decay_xray_filename.c_str() );
  rapidxml::print<char>( outputfile, doc );

  cout << "Wrote " << sandia_decay_xray_filename << endl;

  cout << "BTW, as of 20180503 - this program was totally untested." << endl;
  
  return EXIT_SUCCESS;
}//int add_xray_info_to_sandia_decay_xml()



vector< pair<float,float> > test_gadras_xray_calc_for_nuc( const SandiaDecay::Nuclide *nuclide )
{
//  const SandiaDecay::SandiaDecayDataBase *database = DecayDataBaseServer::database();
//  const SandiaDecay::Nuclide *nuclide = database->nuclide( "Cs137" );
  vector< pair<float,float> > result;

  const char *toriDbFileName = tori_input_filename.c_str();
  // SQLite3 Tori database
  sqlite3 *toriDatabase = NULL;

  if( sqlite3_open( toriDbFileName, &toriDatabase) != 0 )
  {
    cerr << "Failed to open: " << toriDbFileName << endl;
    return result;//EXIT_FAILURE;
  }


  const double age = 2.5*nuclide->halfLife;//3600*24*60;
  const double origParentActivity = 1.0E6 * SandiaDecay::becquerel;

  const vector<SandiaDecay::NuclideActivityPair> daughters
      = SandiaDecay::SandiaDecayDataBase::decay( nuclide, origParentActivity, age );

  double parentAgedActivity = 0.0;
  for( size_t i = 0; i < daughters.size(); ++i )
  {
    if( daughters[i].nuclide == nuclide )
      parentAgedActivity = daughters[i].activity;
  }

//  cerr << "parentAgedActivity=" << parentAgedActivity << endl;

  vector<const SandiaDecay::Transition *> transitions;
  vector<float> transitionIntensities;

  for( size_t i = 0; i < daughters.size(); ++i )
  {
//    cerr << daughters[i].nuclide->symbol << ", has activity " << daughters[i].activity << ", and parentAgedActivity=" << parentAgedActivity << endl;
    const double activityScaler = daughters[i].activity / parentAgedActivity;
//    const double activityScaler = daughters[i].activity / origParentActivity;
//    cerr << "activityScaler=" << activityScaler << endl;

    const vector<const SandiaDecay::Transition *> daughterTransitions = daughters[i].nuclide->decaysToChildren;

    for( size_t j = 0; j < daughterTransitions.size(); ++j )
    {
      const SandiaDecay::Transition *transition = daughterTransitions[j];
      if( transition->child ) // will be null for spontaneous fission
      {
        //if (transition->child->isomerNumber == 0)
        {
          transitions.push_back( transition );
          transitionIntensities.push_back( activityScaler * transition->branchRatio );  //If set activityScaler to 1, then get 5.65685
        }
      }//if( transition->child )
    }//for( size_t j = 0; j < daughterTransitions.size(); ++j )
  }//for( size_t i = 0; i < daughters.size(); ++i )


  vector<int> parentAtomicNumbers;
  vector<int> parentMassNumbers;
  vector<int> parentIsomerNumbers;
  vector<int> daughterAtomicNumbers;
  vector<float> xrayTransitionIntensities;

  // the ToRI x-ray database will return x-rays for GS->isomeric state transitions,
  // even if the x-rays originate from the isomeric de-excitation
  // so, we don't want to repeat different isomeric states of the same daughter

  for( int i = 0; i < transitions.size(); ++i )
  {
    const int parentZ = transitions[i]->parent->atomicNumber;
    const int parentA = transitions[i]->parent->massNumber;
    const int parentM = transitions[i]->parent->isomerNumber;
    const int daughterZ = transitions[i]->child->atomicNumber;

    //XXX wcjohns: I'm not sure this duplicate detection is doing what its supposed to.
    bool duplicate = false;
    const int numExistingTransitions = static_cast<int>(daughterAtomicNumbers.size());

    for (int j = 0; j < numExistingTransitions; ++j)
    {

      if ( parentZ == parentAtomicNumbers[j] &&  parentA == parentMassNumbers[j]
          &&  parentM == parentIsomerNumbers[j] && daughterZ == daughterAtomicNumbers[j])
      {
        duplicate = true;
        break;
      }
    }
    if( !duplicate )
    {
      parentAtomicNumbers.push_back(parentZ);
      parentMassNumbers.push_back(parentA);
      parentIsomerNumbers.push_back(parentM);
      daughterAtomicNumbers.push_back(daughterZ);
      xrayTransitionIntensities.push_back( transitionIntensities[i] );
    }else
    {
//      cerr << "Duplicate: " << nuclide->symbol << ": " << transitions[i]->parent->symbol << "->" << transitions[i]->child->symbol << "\n";
      /*
       Duplicate: Al34: Al34->Si34
       Duplicate: Al35: Al35->Si35
       Duplicate: As85: As85->Se85
       Duplicate: Br87: Br87->Kr87
       Duplicate: Br88: Br88->Kr88
       Duplicate: Br89: Sr89->Y89m
       Duplicate: Br89: Br89->Kr89
       Duplicate: Br90: Sr89->Y89m
       Duplicate: Br90: Br90->Kr90
       Duplicate: Br92: Br92->Kr92
       Duplicate: Br93: Rb93->Sr93
       Duplicate: Br93: Kr93->Rb93
       Duplicate: Br93: Br93->Kr93
       Duplicate: Br94: Rb93->Sr93
       Duplicate: Br94: Kr93->Rb93
       Duplicate: Cd127: In127->Sn127m
       Duplicate: Cd132: In132->Sn132m
       Duplicate: Cu76: Cu76->Zn76
       Duplicate: F25: F25->Ne25
       Duplicate: Ga79: Ga79->Ge79m
       Duplicate: Ga80: Ga80->Ge80
       Duplicate: Ga81: Ga81->Ge81m
       Duplicate: Ga82: Ga82->Ge82
       Duplicate: Ge85: As85->Se85
       Duplicate: Hg179: Hg179->Pt178
       Duplicate: Hg181: Hg181->Pt180
       Duplicate: Hg183: Hg183->Pt182
       Duplicate: I137: I137->Xe137
       Duplicate: I138: I138->Xe138
       Duplicate: I139: I139->Xe139
       Duplicate: In127: In127->Sn127m
       Duplicate: In127m: In127m->Sn127m
       Duplicate: In129: In129->Sn129m
       Duplicate: In129m: In129m->Sn129
       Duplicate: In131: In131->Sn131m
       Duplicate: In131m: In131m->Sn131
       Duplicate: In131m2: In131m2->Sn131m
       Duplicate: In132: In132->Sn132m
       Duplicate: In133: Sn133->Sb133
       Duplicate: In133: In133->Sn133
       Duplicate: In134: Sn133->Sb133
       Duplicate: K49: K49->Ca49
       Duplicate: K50: K50->Ca50
       Duplicate: K51: K51->Ca51
       Duplicate: K52: K52->Ca52
       Duplicate: Kr89: Sr89->Y89m
       Duplicate: Kr93: Rb93->Sr93
       Duplicate: Kr93: Kr93->Rb93
       Duplicate: Kr94: Rb93->Sr93
       Duplicate: Kr94: Rb94->Sr94
       Duplicate: Kr94: Kr94->Rb94
       Duplicate: Kr95: Zr95->Nb95m
       Duplicate: Kr95: Rb95->Sr95
       Duplicate: Kr99: Y98->Zr98
       Duplicate: Kr99: Sr98->Y98
       Duplicate: Kr99: Rb98->Sr97
       Duplicate: Kr99: Rb98->Sr98
       Duplicate: La149: La149->Ce149
       Duplicate: Mg33: Mg33->Al33
       Duplicate: Mg34: Al34->Si34
       Duplicate: Na27: Na27->Mg27
       Duplicate: Na29: Na29->Mg29
       Duplicate: Na30: Na30->Mg30
       Duplicate: Na31: Na31->Mg30
       Duplicate: Na31: Na31->Mg31
       Duplicate: Na32: Na32->Mg31
       Duplicate: Na32: Na32->Mg32
       Duplicate: Na33: Mg33->Al33
       Duplicate: Na33: Na33->Mg32
       Duplicate: Na33: Na33->Mg33
       Duplicate: Na34: Al34->Si34
       Duplicate: Na34: Na34->Mg34
       Duplicate: Nd133: Nd133->Pr133m
       Duplicate: Ne29: Na29->Mg29
       Duplicate: Ne30: Na30->Mg30
       Duplicate: O24: O24->F24
       Duplicate: P38: P38->S38
       Duplicate: P39: P39->S39
       Duplicate: P40: P40->S40
       Duplicate: P41: P41->S41
       Duplicate: P42: P42->S42
       Duplicate: Pb185: Hg181->Pt180
       Duplicate: Pb187: Hg183->Pt182
       Duplicate: Pb187m: Hg183->Pt182
       Duplicate: Pm133: Nd133->Pr133m
       Duplicate: Po189: Hg181->Pt180
       Duplicate: Po191: Hg183->Pt182
       Duplicate: Ra203: Ra203->Rn199m
       Duplicate: Rb100: Y98->Zr98
       Duplicate: Rb100: Y99->Zr99
       Duplicate: Rb100: Sr99->Y99
       Duplicate: Rb100: Y100->Zr100
       Duplicate: Rb100: Sr100->Y100
       Duplicate: Rb100: Rb100->Sr100
       Duplicate: Rb101: Y99->Zr99
       Duplicate: Rb101: Y100->Zr100
       Duplicate: Rb101: Sr100->Y100
       Duplicate: Rb101: Y101->Zr101
       Duplicate: Rb101: Sr101->Y101
       Duplicate: Rb101: Rb101->Sr101
       Duplicate: Rb89: Sr89->Y89m
       Duplicate: Rb93: Rb93->Sr93
       Duplicate: Rb94: Rb94->Sr94
       Duplicate: Rb95: Zr95->Nb95m
       Duplicate: Rb95: Rb95->Sr95
       Duplicate: Rb96: Zr95->Nb95m
       Duplicate: Rb96: Rb96->Sr96
       Duplicate: Rb97: Rb97->Sr97
       Duplicate: Rb98: Y98->Zr98
       Duplicate: Rb98: Sr98->Y98
       Duplicate: Rb98: Rb98->Sr97
       Duplicate: Rb98: Rb98->Sr98
       Duplicate: Rb98m: Y98->Zr98
       Duplicate: Rb98m: Sr98->Y98
       Duplicate: Rb99: Y98->Zr98
       Duplicate: Rb99: Sr98->Y98
       Duplicate: Rb99: Y99->Zr99
       Duplicate: Rb99: Sr99->Y99
       Duplicate: Rb99: Rb99->Sr99
       Duplicate: S43: S43->Cl43
       Duplicate: S44: S44->Cl44
       Duplicate: Sb135: Sb135->Te135
       Duplicate: Sb136: Te136->I136m
       Duplicate: Sb136: Sb136->Te136
       Duplicate: Se87: Br87->Kr87
       Duplicate: Se87: Se87->Br87
       Duplicate: Se88: Br88->Kr88
       Duplicate: Si36: Si36->P36
       Duplicate: Sn133: Sn133->Sb133
       Duplicate: Sn135: Sb135->Te135
       Duplicate: Sn135: Sn135->Sb135
       Duplicate: Sn136: Sb135->Te135
       Duplicate: Sr100: Y99->Zr99
       Duplicate: Sr100: Y100->Zr100
       Duplicate: Sr100: Sr100->Y100
       Duplicate: Sr101: Y100->Zr100
       Duplicate: Sr101: Y101->Zr101
       Duplicate: Sr101: Sr101->Y101
       Duplicate: Sr89: Sr89->Y89m
       Duplicate: Sr95: Zr95->Nb95m
       Duplicate: Sr98: Y98->Zr98
       Duplicate: Sr98: Sr98->Y98
       Duplicate: Sr99: Y98->Zr98
       Duplicate: Sr99: Y99->Zr99
       Duplicate: Sr99: Sr99->Y99
       Duplicate: Te136: Te136->I136m
       Duplicate: Te137: I137->Xe137
       Duplicate: Te137: Te137->I137
       Duplicate: Tm147: Tm147->Er147
       Duplicate: Y100: Y100->Zr100
       Duplicate: Y101: Y101->Zr101
       Duplicate: Y95: Zr95->Nb95m
       Duplicate: Y98: Y98->Zr98
       Duplicate: Y99: Y99->Zr99
       Duplicate: Y99m: Y99->Zr99
       Duplicate: Zn79: Ga79->Ge79m
       Duplicate: Zn79: Zn79->Ga79
       Duplicate: Zn80: Ga79->Ge79m
       Duplicate: Zn80: Ga80->Ge80
       Duplicate: Zn80: Zn80->Ga80
       Duplicate: Zn81: Ga80->Ge80
       Duplicate: Zr95: Zr95->Nb95m
       */
    }
      //cerr << "Duplicate: " << nuclide->symbol << "\n";
  }//for( int i = 0; i < transitions.size(); ++i )

  for( size_t i = 0; i < parentAtomicNumbers.size(); ++i )
  {
    vector< pair<float,float> > thisresult;

    const int parentAN = parentAtomicNumbers[i];
    const int parentIsomer = parentIsomerNumbers[i];
    const int parentMassNumber = parentMassNumbers[i];
    const int daughterAN = daughterAtomicNumbers[i];

//    string p = transitions[i]->parent->symbol, c;
//    if( transitions[i]->child )
//      c = transitions[i]->child->symbol;
//    cerr << p << "--->" << c << " gives:" << endl;

    char query_buffer[512];
    snprintf( query_buffer, sizeof(query_buffer),
              "SELECT XRays.Energy, XIntensities.Int "
              "FROM ((XRays INNER JOIN XIntensities ON XRays.XCode = XIntensities.XCode) "
              "INNER JOIN parents ON XIntensities.iZA = parents.iZA) "
              "WHERE (parents.Z = %i) AND (parents.A = %d) AND ((parents.iZA - (10000 * parents.Z + parents.A)) / 300 = %d) AND ((XRays.XCode / 100) = %d) AND (XRays.Energy >= 10)",
             parentAN, parentMassNumber, parentIsomer, daughterAN );

    char *errorMessage = NULL;
    void *callbackArg = &thisresult;
    int returnCode = sqlite3_exec( toriDatabase, query_buffer, &tori_xray_sqlite3_callback, callbackArg, &errorMessage );

    if( errorMessage )
      cerr << "SQL error: " << errorMessage << endl;

    const float br = xrayTransitionIntensities[i];

    for( size_t k = 0; k < thisresult.size(); ++k )
    {
//      if( thisresult[k].second < 1.0E-10f )
//        continue;
      if( thisresult[k].second < FLT_MIN )
        continue;

      vector< pair<float,float> >::iterator pos = std::lower_bound( result.begin(), result.end(), thisresult[k].first, &less_than_first );

      if( (pos == result.end()) || (pos->first != thisresult[k].first) )
      {
        result.insert( pos, make_pair(thisresult[k].first, br*thisresult[k].second) );
      }else
      {
        pos->second += br*thisresult[k].second;
      }
    }
  }//for( size_t i = 0; i < parentAtomicNumbers.size(); ++i )


  sqlite3_close( toriDatabase );
  toriDatabase = NULL;


//  for( size_t i = 0; i < result.size(); ++i )
//  {
//    cerr << "{ " << result[i].first << ", " << result[i].second << " }\n";
//  }


  return result;
}//int test_gadras_xray_calc_for_nuc()


int test_sandia_decay_xrays()
{
  SandiaDecay::SandiaDecayDataBase database;
  database.initialize( sandia_decay_xray_filename );

  for( const SandiaDecay::Nuclide *nuclide : database.nuclides() )
  {
    try
    {
      vector< pair<float,float> > sqlanswer = test_gadras_xray_calc_for_nuc( nuclide );

      std::sort( sqlanswer.begin(), sqlanswer.end() );

      const double origActivity = 1.0E6 * SandiaDecay::becquerel;
      const double age = 2.5*nuclide->halfLife;//3600*24*60;

      if( nuclide->decaysToChildren.empty() || nuclide->descendants().empty()
          || (nuclide->decaysToChildren.size()==1 && !nuclide->decaysToChildren[0]->child) )
        continue;

      SandiaDecay::NuclideMixture mix;
      mix.addNuclideByActivity( nuclide, origActivity );
      std::vector<SandiaDecay::EnergyRatePair> dexayanswer = mix.xrays( age, SandiaDecay::NuclideMixture::OrderByEnergy );

      const double decayedActivity = mix.activity( age, nuclide );

      for( size_t i = 0; i < dexayanswer.size(); ++i )
        dexayanswer[i].numPerSecond /= decayedActivity;

      //Combine same energy elements together (actually, this shouldnt be necassary).
      bool hasdup = true;
      while( hasdup )
      {
        hasdup = false;
        for( int i = 0; i < int(dexayanswer.size()-1); ++i )
        {
          if( dexayanswer[i].energy == dexayanswer[i+1].energy )
          {
            cerr << "Adding in " << dexayanswer[i].numPerSecond << " for " << dexayanswer[i].energy << " keV" << endl;
            dexayanswer[i].numPerSecond += dexayanswer[i+1].numPerSecond;
            dexayanswer.erase( dexayanswer.begin() + i+1 );
            hasdup = true;
            i = -1;
          }//if( dexayanswer[i].energy == dexayanswer[i+1].energy )
        }//for( int i = 0; i < int(dexayanswer.size()-1); ++i )
      }//while( hasdup )

/*
      const double minbr = -1.0E-6;
      const double minenergy = 2.0;

      {
        vector< pair<float,float> > filteredsqlanswer;
        std::vector<SandiaDecay::EnergyRatePair> dexayanswerfiltered;
        for( size_t i = 0; i < sqlanswer.size(); ++i )
        {
          if( sqlanswer[i].first > minenergy && sqlanswer[i].second > minbr )
            filteredsqlanswer.push_back( sqlanswer[i] );
        }

        for( size_t i = 0; i < dexayanswer.size(); ++i )
        {
          if( dexayanswer[i].energy > minenergy && dexayanswer[i].numPerSecond > minbr )
            dexayanswerfiltered.push_back( dexayanswer[i] );
        }

        sqlanswer.swap( filteredsqlanswer );
        dexayanswer.swap( dexayanswerfiltered );
      }
*/

      char buffer[256];

      if( dexayanswer.size() != sqlanswer.size() )
      {
        stringstream msg;
        msg << "Different number of xrays, xml vs sql gives: " << dexayanswer.size() << " vs " << sqlanswer.size() << endl;


        vector<float> xmlenergies, sqlenergies;
        for( size_t i = 0; i < dexayanswer.size(); ++i )
          xmlenergies.push_back( dexayanswer[i].energy );
        for( size_t a = 0; a < sqlanswer.size(); ++a )
          sqlenergies.push_back( sqlanswer[a].first );




        for( size_t a = 0; a < dexayanswer.size(); ++a )
        {
          if( std::find(sqlenergies.begin(),sqlenergies.end(),float(dexayanswer[a].energy)) == sqlenergies.end() )
            msg << "\tdec has: {" << dexayanswer[a].energy << " keV, " << dexayanswer[a].numPerSecond << "}\n";
        }

        for( size_t a = 0; a < sqlanswer.size(); ++a )
        {
          if( std::find(xmlenergies.begin(),xmlenergies.end(),float(sqlanswer[a].first)) == xmlenergies.end() )
            msg << "\tsql has: {" << sqlanswer[a].first << " keV, " << sqlanswer[a].second << "}\n";
        }
        throw runtime_error( msg.str() );
      }

//      for( size_t i = 0; i < dexayanswer.size(); ++i )
//        cerr << "{ " << dexayanswer[i].energy << ", " << dexayanswer[i].numPerSecond << "}\n";

      stringstream intensity_error;

      for( size_t i = 0; i < dexayanswer.size(); ++i )
      {
        const double sql_energy = sqlanswer[i].first;
        const double sql_intensity = sqlanswer[i].second;

        const double decay_energy = dexayanswer[i].energy;
        const double decay_intensity = dexayanswer[i].numPerSecond;

        if( std::fabs(decay_energy - sql_energy) > 1.0E-6 )
        {
          snprintf( buffer, sizeof(buffer), "Different energy (%f vs %f) for xray %i", decay_energy, sql_energy, (int)i );
          throw runtime_error( buffer );
        }

        const double diff = std::fabs(decay_intensity - sql_intensity);
        const double maxval = std::max( std::fabs(decay_intensity), std::fabs(sql_intensity) );

//        cerr << "For " << sql_energy << "keV xml/sql=" << (sql_intensity/decay_intensity) << endl;

        if( (diff > 1.0E-6 || diff > 0.0001*maxval) && maxval > 1.0E-15 )
        {
          snprintf( buffer, sizeof(buffer), "\tDifferent rate xml vs sql (%g vs %g) for xray at %f kev\n", decay_intensity, sql_intensity, sql_energy );
          intensity_error << buffer;
//          throw runtime_error( buffer );
        }
      }//for( size_t i = 0; i < dexayanswer.size(); ++i )

      if( intensity_error.str().size() )
        throw runtime_error( "\n" + intensity_error.str() );

    }catch( std::exception &e )
    {
      cerr << "Failed xray compare for nuclide " << nuclide->symbol << ": " << e.what() << endl;
    }
  }

  cout << "Done in test_sandia_decay_xrays()" << endl;

  return EXIT_SUCCESS;
}//int test_sandia_decay_xrays()
