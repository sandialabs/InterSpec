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

#include <string>
#include <vector>
#include <algorithm>

#include <boost/any.hpp>

#include <Wt/WString>
#include <Wt/WSuggestionPopup>
#include <Wt/WAbstractItemModel>

#include "InterSpec/PeakDef.h"
#include "SpecUtils/StringAlgo.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/IsotopeNameFilterModel.h"

using namespace std;
using namespace Wt;

/*
 * See also: http://www.webtoolkit.eu/wt/blog/2010/03/02/javascript_that_is_c__
 */
#define INLINE_JAVASCRIPT(...) #__VA_ARGS__


namespace
{
  bool less_than_by_name( const ReactionGamma::Reaction *lhs,
                          const ReactionGamma::Reaction *rhs)
  {
    return ((lhs && rhs) ? (lhs->name() < rhs->name()) : false);
  }
  
  template<class T> struct index_compare_assend
  {
    index_compare_assend(const T arr) : arr(arr) {} //pass the actual values you want sorted into here
    bool operator()(const size_t a, const size_t b) const
    {
      return arr[a] < arr[b];
    }
    const T arr;
  };//struct index_compare
  
}//namespace

IsotopeNameFilterModel::IsotopeNameFilterModel( WObject *parent )
  : WAbstractItemModel( parent ),
    m_minHalfLife( 0.0 ),
    m_includeXray( true ),
    m_includeEscape( true ),
    m_includeNuclides( true ),
    m_includeReactions( true )
{
}


IsotopeNameFilterModel::~IsotopeNameFilterModel()
{
}


void IsotopeNameFilterModel::excludeNuclides( const bool exclude )
{
  m_includeNuclides = !exclude;
}//void excludeNuclides(...);


void IsotopeNameFilterModel::excludeXrays( const bool exclude )
{
  m_includeXray = !exclude;
}//void excludeXrays(...);


void IsotopeNameFilterModel::excludeEscapes( const bool exclude )
{
  m_includeEscape = !exclude;
}//void excludeEscapes(...);


void IsotopeNameFilterModel::excludeReactions( const bool exclude )
{
  m_includeReactions = !exclude;
}//void excludeReactions(...)


Wt::WModelIndex IsotopeNameFilterModel::index( int row, int column,
                                                     const Wt::WModelIndex &parent ) const
{
  if( parent.isValid() || (column!=0) || (row<0) )
    return WModelIndex();

  const int nrows = static_cast<int>( m_candidatesNuclides.size()
                                      + m_candidatesElements.size()
                                      + m_candidatesReactions.size()
                                      + m_customSuggests.size() );

  if( row >= nrows )
    return WModelIndex();

  return createIndex( row, column, (void *)NULL );
}//index(...)


Wt::WModelIndex IsotopeNameFilterModel::parent( const Wt::WModelIndex & ) const
{
  return WModelIndex();
}//WModelIndex parent(...)


int IsotopeNameFilterModel::rowCount( const Wt::WModelIndex &parent ) const
{
  if( parent.isValid() )
    return 0;
  return static_cast<int>( m_candidatesNuclides.size()
                           + m_candidatesElements.size()
                           + m_candidatesReactions.size()
                           + m_customSuggests.size() );
}


int IsotopeNameFilterModel::columnCount( const Wt::WModelIndex &parent ) const
{
  if( parent.isValid() )
    return 0;
  return 1;
}


boost::any IsotopeNameFilterModel::data( const Wt::WModelIndex &index, int role ) const
{
  const int row = index.row();
  const int column = index.column();

  const int nnuc = static_cast<int>( m_candidatesNuclides.size() );
  const int nel = static_cast<int>( m_candidatesElements.size() );
  const int nrctn = static_cast<int>( m_candidatesReactions.size() );
  const int ncustom = static_cast<int>( m_customSuggests.size() );
  const int nrows = nnuc + nel + nrctn + ncustom;

  if( row<0 || row>=nrows || column!=0 )
    return boost::any();

//we could customize tool tip or display roles here to give more information...
//  if( role == Wt::ToolTipRole )
//    return WString( "the_tool_tip" );

  if( row < nel )
  {
    const SandiaDecay::Element *el = m_candidatesElements[row];
    return boost::any( WString( el->symbol ) );
  }else if( row < (nel+nnuc) )
  {
    const int nucnum = row - nel;
    const SandiaDecay::Nuclide *nuc = m_candidatesNuclides[nucnum];
    return boost::any( WString( m_typePrefix + nuc->symbol ) );
  }else if( row < (nel+nnuc+nrctn) )
  {
    const int rctnum = row - nel - nnuc;
    const ReactionGamma::Reaction *rctn = m_candidatesReactions[rctnum];
    return boost::any( WString( m_typePrefix + rctn->name() ) );
  }else
  {
    return boost::any( m_customSuggests[row-nel-nnuc-nrctn] );
  }
  
}//boost::any IsotopeNameFilterModel::data( const Wt::WModelIndex &index, int role = Wt::DisplayRole ) const


int IsotopeNameFilterModel::determineAndRemoveIsoLevel( std::string &label )
{
  // taken from IsotopeNameFilterModel 20121014
  //determineAndRemoveIsoLevel(...) looks for patterns such as 'Co60m', 'Co60meta',
  //  'Co 60 m', etc. to determine if the user is inputting a metts stable state.
  //  The function returns the iso level (right no just 0, 1, or 2), and
  //  removes the portion of the text indicating the meta level, from the input.

  int user_input_is_meta = 0;
  size_t len = label.size();
  size_t mpos = label.find("meta");

  if( mpos != string::npos )
  {
    if( label.find("meta2") == string::npos )
    {
      label.erase( mpos, 4 );
      user_input_is_meta = 1;
    }else
    {
      label.erase( mpos, 5 );
      user_input_is_meta = 2;
    }
  }else if( label.find('m') != string::npos )
  {
    mpos = 0;
    do
    {
      if( len==0 )
        break;

      mpos = label.find('m', mpos);

      if( mpos == string::npos )
        break;

      if( mpos == (len-1) )
      {
        size_t prev = ((mpos > 0) ? mpos-1 : mpos);
        while( prev > 0 && !isalpha(label[prev]) && !isdigit(label[prev]) )
          --prev;

        if( prev > 0 )
        {
          user_input_is_meta = (isdigit( label[prev] ) ? 1 : 0);
          if( user_input_is_meta )
          {
            label.erase( mpos, 1 );
//            len = label.length();
            break;
          }//if( user_input_is_meta )
        }//if( prev > 0 )
      }else if( mpos > 0 )
      {
        if( (!isdigit(label[mpos+1]) && !isalpha( label[mpos+1]))
            || (label[mpos+1]=='2' && ( (mpos+1)==(len-1) || (!isdigit(label[mpos+2]) && !isalpha( label[mpos+2])) ) ) )
        {
          const bool is2m = (label[mpos+1]=='2'
                             && ( (mpos+1)==(len-1)
                                  || (!isdigit(label[mpos+2]) && !isalpha( label[mpos+2])) ));
          size_t prev = mpos - 1;
          while( prev > 0 && !isalpha(label[prev]) && !isdigit(label[prev]) )
            --prev;
          if( prev > 0 )
          {
            user_input_is_meta = (isdigit( label[prev] ) ? 1 : 0);
            if( user_input_is_meta )
            {
              if( !is2m )
                label.erase( mpos, 1 );
              else
              {
                user_input_is_meta = 2;
                label.erase( mpos, 2 );
              }
//              len = label.length();
              break;
            }//if( user_input_is_meta )
          }//if( prev > 0 )
        }//if( !isdigit(label[mpos+1]) && !isalpha( label[mpos+1]) )
      }//if( mpos == (len-1) ) / else

      mpos += 1;
    }while( mpos != string::npos && mpos < len );
  }//if( has "meta" ) / else has 'm'

  return user_input_is_meta;
}//int determineAndRemoveIsoLevel( std::string &label )


void IsotopeNameFilterModel::getAlphaAndNumericSubStrs( std::string label,
                                       std::vector<std::string> &alphastrs,
                                       std::vector<std::string> &numericstrs )
{
  SpecUtils::erase_any_character( label, " -_,\t<>/?[]{}\\|!@#$%^&*();:\"'~`+=" );
  SpecUtils::to_lower_ascii( label );

  const string::size_type len = label.length();

  for( size_t i = 0; i < len; )
  {
    string strStr, numStr;
    while( (i<len) && isalpha( label[i] ) )
      strStr += label[i++];
    while( (i<len) && isdigit( label[i] ) )
      numStr += label[i++];

    if( strStr.size() )
      alphastrs.push_back( strStr );
    if( numStr.size() )
      numericstrs.push_back( numStr );
    if( numStr.empty() && strStr.empty() )
      ++i;
  }//for( size_t i = 0; i < len; )

}//void IsotopeNameFilterModel::getAlphaAndNumericSubStrs(...)

void IsotopeNameFilterModel::addCustomSuggestPossibility( const std::string &str )
{
  m_customPotentials.push_back( str );
}


void IsotopeNameFilterModel::filter( const Wt::WString &text )
{
  const int ninitialrow = rowCount();
  if( ninitialrow > 0 )
  {
    beginRemoveRows( WModelIndex(), 0, ninitialrow - 1 );
    m_candidatesElements.clear();
    m_candidatesNuclides.clear();
    m_candidatesReactions.clear();
    m_customSuggests.clear();
    m_typePrefix = "";
    endRemoveRows();
  }//if( ninitialrow > 0 )

  string testTxt = text.toUTF8();

  PeakDef::SourceGammaType srctype;
  PeakDef::gammaTypeFromUserInput( testTxt, srctype );
  
  switch( srctype )
  {
    case PeakDef::NormalGamma:
    case PeakDef::AnnihilationGamma:
      m_typePrefix = "";
    break;
      
    case PeakDef::SingleEscapeGamma:
      m_typePrefix = (m_includeEscape ? "S.E. " : "");
    break;
      
    case PeakDef::DoubleEscapeGamma:
      m_typePrefix = (m_includeEscape ? "D.E. " : "");
    break;
      
    case PeakDef::XrayGamma:
      m_typePrefix = (m_includeXray ? "xray " : "");
    break;
  }//switch( srctype )
  
  
//  bool specificyXray = false;
//  string::size_type pos = testTxt.find_first_of( ' ' );
//  if( pos != string::npos )
//  {
//    specificyXray = (testTxt.find( "xray",pos) != string::npos)
//                    || (testTxt.find( "x-ray",pos) != string::npos);
//    SpecUtils::replace_all( testTxt, "xray", "" );
//    SpecUtils::replace_all( testTxt, "x-ray", "" );
//  }
  
  //Make sure the user isnt typing in a reaction, otherwise we wont get matching
  //  isotoes
  string::size_type open_paren = testTxt.find( "(" );
  const bool reactionsOnly = (open_paren != string::npos);
  if( reactionsOnly )
    testTxt = testTxt.substr( 0, open_paren );

  const int metalevel = determineAndRemoveIsoLevel( testTxt );

  vector<string> alphastrs, numericstrs;
  getAlphaAndNumericSubStrs( testTxt, alphastrs, numericstrs );

  if( alphastrs.empty() && numericstrs.empty() )
    return;

  vector<const SandiaDecay::Nuclide *> suggestions;
  vector< const SandiaDecay::Element * > suggest_elements;
  set<const SandiaDecay::Element *> candidate_elements = possibleElements( alphastrs );
  
  suggestNuclides( alphastrs, numericstrs, metalevel,
                   candidate_elements, suggestions, suggest_elements );
  
  //Get the potential reaction such as Fe(n,g),...
  //XXX - I'm not too happy with this section of code, its sloppy, kinda a mess,
  //      and is probably pretty ineffiecnt
  //Should also be seperated out into its own funciton
  vector<const ReactionGamma::Reaction *> suggest_reactions;
  if( m_includeReactions && metalevel==0 )
    suggest_reactions = suggestReactions( text, suggestions, candidate_elements );

  //if user has entered something like "Fe(...", they only want the reactions
  if( reactionsOnly )
  {
    suggestions.clear();
    suggest_elements.clear();
  }//if( reactionsOnly )
  
  if( !m_includeNuclides )
    suggestions.clear();

  if( !m_includeXray )
    suggest_elements.clear();
  
  testTxt = text.toUTF8();
  vector<WString> customSuggests;
  for( const std::string &str : m_customPotentials )
  {
    if( SpecUtils::icontains( str, testTxt ) )
      customSuggests.push_back( str );
  }//for( const std::string &str : m_customPotentials )
  
  
  //Sort suggestions by Levenshtein distance, so that way the word closest to
  //  what the user has typed in, will be at the top of the list.  If we dont
  //  do this, then if the user types "Ra226" and hits enter, then "Rn226"
  //  will actually be selected.
  using SpecUtils::levenshtein_distance;
  testTxt = text.toUTF8();
  
  //Sort suggested nuclides
  vector<unsigned int> distance( suggestions.size() );
  vector<size_t> sort_indices( suggestions.size() );
  for( size_t i = 0; i < suggestions.size(); ++i )
  {
    sort_indices[i] = i;
    const string &symbol = suggestions[i]->symbol;
    distance[i] = levenshtein_distance( testTxt, symbol );
  }//for( size_t i = 0; i < suggestions.size(); ++i )
  std::stable_sort( sort_indices.begin(), sort_indices.end(),
                    index_compare_assend<vector<unsigned int>&>(distance) );

  //Sort suggested elements
  distance.resize( suggest_elements.size() );
  vector<size_t> element_sort_indices( suggest_elements.size() );
  for( size_t i = 0; i < suggest_elements.size(); ++i )
  {
    element_sort_indices[i] = i;
    const string &symbol = suggest_elements[i]->symbol;
    const string &name = suggest_elements[i]->name;
    distance[i] = std::min( levenshtein_distance( testTxt, symbol ),
                            levenshtein_distance( testTxt, name ) );
  }//for( size_t i = 0; i < suggest_elements.size(); ++i )
  std::stable_sort( element_sort_indices.begin(), element_sort_indices.end(),
                    index_compare_assend<vector<unsigned int>&>(distance) );
  
  //Sort suggested
  distance.resize( customSuggests.size() );
  vector<size_t> custom_sort_indices( customSuggests.size() );
  for( size_t i = 0; i < customSuggests.size(); ++i )
  {
    custom_sort_indices[i] = i;
    const string &symbol = customSuggests[i].toUTF8();
    distance[i] = levenshtein_distance( testTxt, symbol );
  }//for( size_t i = 0; i < customSuggests(); ++i )
  
  std::stable_sort( custom_sort_indices.begin(), custom_sort_indices.end(),
                    index_compare_assend<vector<unsigned int>&>(distance) );
  
  
  const int nrow = static_cast<int>( suggestions.size()
                                     + suggest_elements.size()
                                     + suggest_reactions.size()
                                     + customSuggests.size() );
  
  beginInsertRows( WModelIndex(), 0, nrow -1 );
  
  m_candidatesElements.resize( suggest_elements.size() );
  for( size_t i = 0; i < suggest_elements.size(); ++i )
    m_candidatesElements[i] = suggest_elements[ element_sort_indices[i] ];
  
  m_candidatesNuclides.resize( suggestions.size() );
  for( size_t i = 0; i < suggestions.size(); ++i )
    m_candidatesNuclides[i] = suggestions[ sort_indices[i] ];
  
  m_candidatesReactions.swap( suggest_reactions );
  
  m_customSuggests.resize( customSuggests.size() );
  for( size_t i = 0; i < customSuggests.size(); ++i )
    m_customSuggests[i] = customSuggests[ custom_sort_indices[i] ];
  
  endInsertRows();
}//void IsotopeNameFilterModel::filter( const Wt::WString &text )



std::set<const SandiaDecay::Element *> IsotopeNameFilterModel::possibleElements(
                                        const std::vector<string> &alphastrs )
{
  std::set<const SandiaDecay::Element *> candidate_elements;
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  //suggest based off of alphastrs
  const std::vector<const SandiaDecay::Element *> &elements = db->elements();
  
  for( const SandiaDecay::Element *el : elements )
  {
    const string name   = SpecUtils::to_lower_ascii_copy( el->name );
    const string symbol = SpecUtils::to_lower_ascii_copy( el->symbol );
    
    for( string str : alphastrs )
    {
      if( SpecUtils::starts_with( symbol, str.c_str() )
         || SpecUtils::starts_with( name, str.c_str() ) )
        candidate_elements.insert( el );
    }//for( const string &str : alphastrs )
  }//for( const SandiaDecay::Element *el : elements )
  
  return candidate_elements;
}//possibleElements

void IsotopeNameFilterModel::suggestNuclides(
                                          const vector<string> &alphastrs,
                                          const vector<string> &numericstrs,
                                          const int metalevel,
                                          const std::set<const SandiaDecay::Element *> &candidate_elements,
                                          vector<const SandiaDecay::Nuclide *> &suggestions,
                                          vector< const SandiaDecay::Element * > &suggest_elements )
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  for( const SandiaDecay::Element *el : candidate_elements )
  {
    bool is_exact_element = false;
    if( numericstrs.empty() )
    {
      const string name   = SpecUtils::to_lower_ascii_copy( el->name );
      const string symbol = SpecUtils::to_lower_ascii_copy( el->symbol );
      for( string str : alphastrs )
        is_exact_element |= (symbol==str || name==str);
    }//if( numericstrs.empty() )
    
    if( numericstrs.empty() && !is_exact_element )
    {
      suggest_elements.push_back( el );
    }else
    {
      //Add in just the element name, and then all the isotopes below
      if( numericstrs.empty() && is_exact_element )
        suggest_elements.push_back( el );
      
      const vector<const SandiaDecay::Nuclide *> nuclides = db->nuclides( el );
      
      for( const SandiaDecay::Nuclide *nuc : nuclides )
      {
        if( IsInf(nuc->halfLife) || nuc->decaysToChildren.empty() )
          continue;
        
        bool numeric_compat = false;
        for( const string &str : numericstrs )
          numeric_compat |= SpecUtils::contains( std::to_string(nuc->massNumber), str.c_str() );
        
        if( metalevel > 0 && metalevel!=nuc->isomerNumber )
          numeric_compat = false;
        
        if( numeric_compat || numericstrs.empty() )
          suggestions.push_back( nuc );
      }//for( const SandiaDecay::Nuclide *nuc : nuclides )
    }//if( there are no numbers, and start of an element name ) / else
  }//for( const SandiaDecay::Element *el : candidate_elements )
  
  
  if( alphastrs.empty() )
  {
    const std::vector<const SandiaDecay::Nuclide *> &nucs = db->nuclides();
    std::vector<const SandiaDecay::Nuclide *>::const_iterator pos;
    for( pos = nucs.begin(); pos != nucs.end(); ++pos )
    {
      const SandiaDecay::Nuclide *nuc = (*pos);
      if( IsInf(nuc->halfLife) || nuc->decaysToChildren.empty() )
        continue;
      
      for( string str : numericstrs )
      {
        try
        {
          if( nuc->massNumber==std::stoi(str) )
            suggestions.push_back( nuc );
        }catch(...){ cerr << "Shouldnt ever be here" << endl; }
      }//for( string str : numericstrs )
    }//for( pos = nucs.begin(); pos != nucs.end(); ++pos )
  }//if( the user has only typed in numbers )

}//suggestNuclides(...)



std::vector<const ReactionGamma::Reaction *>
        IsotopeNameFilterModel::suggestReactions( const Wt::WString &text,
              const std::vector<const SandiaDecay::Nuclide *> &suggestions,
              const std::set<const SandiaDecay::Element *> &candidate_elements )
{
  //Get the potential reaction such as Fe(n,g),...
  //XXX - I'm not too happy with this section of code, its sloppy, kinda a mess,
  //      and is probably pretty ineffiecnt
  //Should also be seperated out into its own funciton
  const ReactionGamma *reactionDb = NULL;
  vector<const ReactionGamma::Reaction *> suggest_reactions;
  try
  {
    reactionDb = ReactionGammaServer::database();
  }catch(...)
  {
    cerr << "Failed to open gamma reactions XML file" << endl;
    return vector<const ReactionGamma::Reaction *>();
  }//try / catch
  
  std::string testTxt = text.toUTF8();
  

  if( SpecUtils::iequals_ascii(testTxt, "a")
      || SpecUtils::istarts_with(testTxt, "an") )
  {
    const vector<const ReactionGamma::Reaction *> &rctns
                                 = reactionDb->reactions(AnnihilationReaction);
    
    if( rctns.size() > 0 )  //should always be the case
      suggest_reactions.push_back( rctns[0] );
  }//if( possibly annihilation )
  
  
  string::size_type open_paren = testTxt.find( "(" );
  const bool reactionsOnly = (open_paren != string::npos);
  if( reactionsOnly )
    testTxt = testTxt.substr( 0, open_paren );
  
  if( testTxt.empty() || !reactionDb )
    return suggest_reactions;
  
  for( const SandiaDecay::Element *el : candidate_elements )
    reactionDb->reactions( el, suggest_reactions );
  for( const SandiaDecay::Nuclide *nuc : suggestions )
    reactionDb->reactions( nuc, suggest_reactions );
    
  //make all rections uniqe and sorted (in not the most effiecient manor...)
  vector<const ReactionGamma::Reaction *>::iterator new_end;
  std::sort( suggest_reactions.begin(), suggest_reactions.end(),
            &less_than_by_name );
  new_end = std::unique( suggest_reactions.begin(), suggest_reactions.end() );
  suggest_reactions.erase( new_end, suggest_reactions.end() );
    
  //Now make sure if the user has specified particle types, only show those
  string rctnTxt = text.toUTF8();
  SpecUtils::to_lower_ascii( rctnTxt );
  SpecUtils::ireplace_all( rctnTxt, " ", "" );
    
  string part1_name, part2_name;
  bool got_part1 = false, got_part2 = false;
  ReactionGamma::ReactionParticle part1, part2;
    
  open_paren = rctnTxt.find( "(" );
  if( open_paren != string::npos && open_paren<(rctnTxt.size()-1) )
  {
    part1_name = rctnTxt.substr( open_paren+1, 1 );
    SpecUtils::trim( part1_name );
    try
    {
      part1 = ReactionGamma::to_particle( part1_name );
      got_part1 = true;
    }catch(...){ /*part1_name.clear();*/ }
      
    string::size_type comma_pos = rctnTxt.find( ",", open_paren );
    if( comma_pos != string::npos && comma_pos<(rctnTxt.size()-1) )
    {
      part2_name = rctnTxt.substr( comma_pos+1, 1 );
      SpecUtils::trim( part2_name );
      try
      {
        part2 = ReactionGamma::to_particle( part2_name );
        got_part2 = true;
      }catch(...) { /*part2_name.clear();*/ }
    }//if( user entered second particle )
  }//if( user enetered a "(" )
  
  
  if( part1_name.empty() && part2_name.empty() )
    return suggest_reactions;
  
  if( (part1_name.size() && !got_part1) || (part2_name.size() && !got_part2) )
  {
    cerr << "Couldnt parse reaction particle '" << part1_name << "' or '"
         << part2_name << "', not returning anything" << endl;
    return vector<const ReactionGamma::Reaction *>();
  }
  
  vector<const ReactionGamma::Reaction *> filtered_reactions;
  for( const ReactionGamma::Reaction *rctn : suggest_reactions )
  {
    if( got_part1
        && ReactionGamma::ingoing_particle(rctn->type) != part1 )
      continue;
    if( got_part2
        && ReactionGamma::outgoing_particle(rctn->type) != part2 )
      continue;
    filtered_reactions.push_back( rctn );
  }//for( const ReactionGamma::Reaction *rctn : suggest_reactions )
  
  return filtered_reactions;
}//void suggestReactions( string testTxt )




//  Wt::WFlags<Wt::ItemFlag> flags( const Wt::WModelIndex & index ) const;


void IsotopeNameFilterModel::nuclideNameMatcherJs( std::string &js )
{
  js = INLINE_JAVASCRIPT
  (
    function( edit )
    {
      try
      {
      if( !edit )
        return function(){return { match : false, suggestion : "" };};
      if( !edit.value )
        return function(){return { match : false, suggestion : "" };};
      
      var value = edit.value;
      
      return function( suggestion )
      {
        //value: what the user has typed in.
        //suggestion: what the potential suggestion in the
        try
        {
          if( !suggestion )
            return value;
          
          //XXX - TODO - right now the below ignores the case of meta-stable elements
          //             and could possibly cause results like <b><b>M</b>g13<b>m</b></b>
          //             ex. Co60m, U235m, etc
          //XXX - TODO - right now the gamma lines in parenthesis, are also
          //             matched in the regex, this should be avoided
          //replace the special characters the user may have typed in
          //value = value.replace(/[.*+?|()\\[\\]{}\\\\$^]/g, "\\$&");
          value = value.replace(/[.*+?|\\[\\]{}\\\\$^]/g, "");

          //Highlight the matching numeric parts of the answer
          var numericregex = new RegExp('(\\d+)');
          var numericmatch = numericregex.exec(value);
          if( numericmatch && numericmatch.length > 1 )
           for( var i = 1; i < numericmatch.length; ++i )
           {
             var regex = new RegExp( '(' + numericmatch[i] + ')', 'gi' );
             suggestion = suggestion.replace( regex, "<b>$1</b>" );
           }

          //Highlight the matching alpha parts of the answer
          var alpharegex = new RegExp('([a-zA-Z]+)');
          var alphamatch = alpharegex.exec(value);
          if( alphamatch && alphamatch.length > 1 )
            for( var i = 1; i < alphamatch.length; ++i )
              if( (alphamatch[i] != 'm' && alphamatch[i] != 'M')
                  || !numericmatch || numericmatch.length<2 )
              {
                var regex = new RegExp( '(' + alphamatch[i] + '(?!>))', 'gi' );  //lookahead is to keep form matching <b>
                suggestion = suggestion.replace( regex, "<b>$1</b>" );
              }

           //For right now we will filter everything server side, so all
           //  suggestions will be a match.
          return {match : true, suggestion : suggestion};
        }catch(e)
        {
          console.log( 'caught me' );
          return {match : false, suggestion : ""};
        }
      }
      }catch(e)
      {
        console.log( 'caught me here' );
      }
    }
  );
}//void nuclideNameMatcherJs( std::string &js )


void IsotopeNameFilterModel::replacerJs( std::string &js )
{
  js = INLINE_JAVASCRIPT
  (
      function (edit, suggestionText, suggestionValue)
      {
        if(!edit) return;
        try
        {
          edit.value = suggestionValue;
          if (edit.selectionStart)
            edit.selectionStart = edit.selectionEnd = suggestionValue.length;
        }catch(e)
        {
          console.log( 'caught me here' );
        }
      }//function (edit, suggestionText, suggestionValue)
   );//js = INLINE_JAVASCRIPT( ... )
}//void replacerJs( std::string &js )


void IsotopeNameFilterModel::setQuickTypeFixHackjs( Wt::WSuggestionPopup *popup )
{
#if( WT_SERIES > 3 || (WT_SERIES==3 && WT_MAJOR>3) || (WT_SERIES==3 && WT_MAJOR==3 && WT_MINOR>4) )
#warning "You can probably remove IsotopeNameFilterModel::setQuickTypeFixHackjs(...), but you should check this"
#endif
  
  string js = INLINE_JAVASCRIPT(
    var addTryCatch = function( elid ){
      var dofix = function(elid){
        var el = Wt.WT.getElement(elid);
        var self = el ? jQuery.data(el, 'obj') : null;
        if( !self ){ //Apparently not immediately available even though m_nuclideSuggest
                     //  should be in the DOM by the time this JS gets executed.
          setTimeout( function(){dofix(elid);}, 500 );
          return;
        }
                                    
        var oldfcn = self.refilter;
        self.refilter = function(value){ try{ oldfcn(value); }catch(e){ console.log('My refilter caught: ' + e ); } };
      };
      dofix(elid);
    };
  );
  
  popup->doJavaScript( js + " addTryCatch('" + popup->id() + "');" );
}//void setQuickTypeFixHackjs( Wt::WSuggestionPopup *popup )
