#include <string>
#include <vector>
#include <fstream>

#include "SandiaDecay.h"

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

using namespace std;

/** Program to remove the edit history of sandia.decay.xml. */

int main( int argc, char **argv )
{
  if( argc != 3 )
  {
    fprintf( stderr, "Usage: %s <input sandia.decay.xml> <output file>\n", argv[0] );
    return EXIT_FAILURE;
  }//if( argc != 3 )
  
  try
  {
    rapidxml::file<char> xmlfile( argv[1] );
    
    rapidxml::xml_document<char> doc;
    doc.parse<rapidxml::parse_full>( xmlfile.data() );
    rapidxml::xml_node<char> *doc_node = doc.first_node();
  
    if( doc_node->first_attribute("version") )
      doc_node = doc_node->next_sibling();

    
    int ntrans = 0, ntrans_removed = 0;
    for( rapidxml::xml_node<char> *xmltrans = doc_node->first_node( "transition" );
        xmltrans;
        xmltrans = xmltrans->next_sibling("transition") )
    {
      vector<rapidxml::xml_node<char> *> to_delete;
      for( rapidxml::xml_node<char> *child = xmltrans->first_node();
          child;
          child = child->next_sibling() )
      {
        rapidxml::xml_attribute<char> *revision = child->first_attribute( "revision" );
        
        if( rapidxml::internal::compare(child->name(), child->name_size(), "insertion", 9, false ) )
        {
          to_delete.push_back( child );
        }else if( revision && rapidxml::internal::compare(revision->value(), revision->value_size(), "deleted", 7, false ) )
        {
          to_delete.push_back( child );
        }else
        {
          vector<rapidxml::xml_node<char> *> gchild_to_delete;
          for( rapidxml::xml_node<char> *gchild = child->first_node();
              gchild;
              gchild = gchild->next_sibling() )
          {
            if( rapidxml::internal::compare(gchild->name(), gchild->name_size(), "edit", 4, false )
               || rapidxml::internal::compare(gchild->name(), gchild->name_size(), "insertion", 9, false )
               || rapidxml::internal::compare(gchild->name(), gchild->name_size(), "deletion", 8, false ) )
              gchild_to_delete.push_back( gchild );
            else if( !rapidxml::internal::compare(gchild->name(), gchild->name_size(), "coincidentgamma", 15, false ) )
              throw runtime_error( "Unknown gchild node: "
                                   + std::string(gchild->name(), gchild->name() + gchild->name_size()) );
          }
          for( size_t i = 0; i < gchild_to_delete.size(); ++i )
            child->remove_node( gchild_to_delete[i] );
        
          if( revision )
            child->remove_attribute( revision );
        }
      }//for( loop over children )
    
      ntrans_removed += to_delete.size();
      for( size_t i = 0; i < to_delete.size(); ++i )
        xmltrans->remove_node( to_delete[i] );
      
      ++ntrans;
    }//loop over transitions.
    
    printf( "There were %i transitions; removed %i of them\n", ntrans, ntrans_removed );
    
    ofstream output( argv[2], ios::out | ios::binary );
    if( !output )
      throw runtime_error( "Could not open output file: " + string(argv[2]) );
    
    int flags = 0;
    rapidxml::print<char>( output, doc, flags );
    printf( "Saved output file '%s'\n", argv[2] );
  }catch( std::exception &e )
  {
    fprintf( stderr, "Error: %s\n", e.what() );
  }//try / catch
  
  
  
  return EXIT_SUCCESS;
}//int main( int argc, char **argv )
