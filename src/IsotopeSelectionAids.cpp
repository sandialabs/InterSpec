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

#include <boost/any.hpp>
#include <boost/regex.hpp>

#include <Wt/WText>
#include <Wt/WLineEdit>
#include <Wt/WHBoxLayout>
#include <Wt/WModelIndex>
#include <Wt/WApplication>
#include <Wt/WEnvironment>
#include <Wt/WSuggestionPopup>
#include <Wt/WContainerWidget>
#include <Wt/WAbstractItemModel>
#include <Wt/WAbstractItemDelegate>

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakModel.h"
#include "SpecUtils/StringAlgo.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/IsotopeSelectionAids.h"
#include "InterSpec/IsotopeNameFilterModel.h"


using namespace std;
using namespace Wt;

/*
 * See also: http://www.webtoolkit.eu/wt/blog/2010/03/02/javascript_that_is_c__
 */
#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

namespace
{
struct NSameCharSortHelper
{
  int nsame;
  size_t index;

  NSameCharSortHelper( int _nsame, size_t _index )
    : nsame(_nsame), index( _index ) {}
  NSameCharSortHelper( const NSameCharSortHelper&rhs )
    : nsame(rhs.nsame), index( rhs.index ) {}
  bool operator<( const NSameCharSortHelper &rhs ) const
  { return nsame > rhs.nsame; } //Sort in acsending order!
};//struct NSameCharSortHelper
  
  //a struct so that we can sort an array using an extanal array of indexes,
  //  with out modifying the data array.
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


namespace IsotopeSelectionAids
{

bool distanceLessThan( const NearGammaInfo &lhs, const NearGammaInfo &rhs )
{
  return (lhs.distance < rhs.distance);
}//bool distanceLessThan( const NearGammaInfo &lhs, const NearGammaInfo &rhs )
  
vector< NearGammaInfo >
      equilibriumGammasByNearestEnergy( const SandiaDecay::Nuclide *nuclide,
                                        const double energy,
                                        const double minRelativeBr,
                                        const bool includeEscapePeaks,
                                        const bool includeXrays )
{
  vector< NearGammaInfo > energies;

  if( !nuclide )
    return energies;

  //energies.first == distance of photopeak from fit energy
  //energies.second.first == energy of photopeak
  //energies.second.second == relative intensity of photopeak

  SandiaDecay::NuclideMixture mixture;
  mixture.addNuclide( SandiaDecay::NuclideActivityPair(nuclide,1.0) );
  
  const double decaytime = PeakDef::defaultDecayTime( nuclide );
  vector<SandiaDecay::EnergyRatePair> gammas
       = mixture.gammas( decaytime, SandiaDecay::NuclideMixture::OrderByAbundance, true );

  const double max_intensity = gammas.size() ? gammas[0].numPerSecond : 1.0;

  for( const SandiaDecay::EnergyRatePair &energy_pair : gammas )
  {
    if( energy_pair.numPerSecond > minRelativeBr*max_intensity )
    {
      NearGammaInfo info;
      info.distance = static_cast<float>( fabs( energy_pair.energy - energy ) );
      info.gamma_energy = static_cast<float>( energy_pair.energy );
      info.relative_intensity = static_cast<float>( energy_pair.numPerSecond / max_intensity );
      info.gamma_type = PeakDef::NormalGamma;
      
      energies.push_back( info );
      
      if( includeEscapePeaks && (energy_pair.energy > 1022.0) )
      {
        info.distance = static_cast<float>( fabs( energy_pair.energy - energy - 510.99891 ) );
        info.gamma_type = PeakDef::SingleEscapeGamma;
        energies.push_back( info );
        
        info.distance = static_cast<float>( fabs( energy_pair.energy - energy - 2.0*510.99891 ) );
        info.gamma_type = PeakDef::DoubleEscapeGamma;
        energies.push_back( info );
      }//if( includeEscapePeaks && (energy_pair.energy > 1022.0) )
    }//if( energy_pair.numPerSecond >  0.0001*max_intensity )
  }//for( const SandiaDecay::AbundanceEnergyPair &energy_pair : gammas )

  
  if( includeXrays )
  {
    const vector<SandiaDecay::EnergyRatePair> xrays = mixture.xrays( decaytime );

    for( const SandiaDecay::EnergyRatePair &energy_pair : xrays )
    {
      if( energy_pair.numPerSecond > minRelativeBr*max_intensity )
      {
        NearGammaInfo info;
        info.distance = static_cast<float>( fabs( energy_pair.energy - energy ) );
        info.gamma_energy = static_cast<float>( energy_pair.energy );
        info.relative_intensity = static_cast<float>( energy_pair.numPerSecond / max_intensity );
        info.gamma_type = PeakDef::XrayGamma;
        
        energies.push_back( info );
      }//if( energy_pair.numPerSecond >  0.0001*max_intensity )
    }//for( const SandiaDecay::AbundanceEnergyPair &energy_pair : gammas )
  }//if( includeXrays )
  
  
  std::sort( energies.begin(), energies.end(), &distanceLessThan );

  return energies;
}//equilibriumGammasByNearestEnergy(...)

}//namespace IsotopeSelectionAids



void PhotopeakDelegate::EditWidget::nuclideNameMatcherJs( std::string &js )
{
  js = INLINE_JAVASCRIPT
  (
    function( edit )
    {
      var value = edit.value;

      return function( suggestion )
      {
        if( !suggestion )
          return value;

        //XXX - TODO - right now the below ignores the case of meta-stable elements
        //             and could possibly cause results like <b><b>M</b>g13<b>m</b></b>
        //             ex. Co60m, U235m, etc
        //XXX - TODO - right now the gamma lines in parenthesis, are also
        //             matched in the regex, this should be avoided

       //replace the special characters the user may have typed in
       //Note that . and () are not in these special characters
       value = value.replace(/[*+?|\\[\\]{}\\\\$^]/g, "\\$&");
        
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
              var regex = new RegExp( '(' + alphamatch[i] + ')', 'gi' );
              suggestion = suggestion.replace( regex, "<b>$1</b>" );
            }

         //For right now we will filter everything server side, so all
         //  suggestions will be a match.
        return { match : true, suggestion : suggestion };
      }
    }
  );
}//void nuclideNameMatcherJs( std::string &js )


void PhotopeakDelegate::EditWidget::gammaEnergyMatcherJs( std::string &js )
{
  /*
   * INLINE_JAVASCRIPT is a macro which allows entry of JavaScript
   * directly in a C++ file - awesome!
   */
  js = INLINE_JAVASCRIPT
  (
    function( edit )
    {
      var value = edit.value;

      return function( suggestion )
      {
        if( !suggestion )
          return value;

        if( suggestion.indexOf( value ) == 0 )
        {
          var highval = suggestion.replace( value, '<b>' + value + '</b>');
          return { match: true, suggestion: highval };
        }//if( pos == 0 )

        return { match: false, suggestion: suggestion };
      }//return function( suggestion )
    }//function( edit )
  );//js = INLINE_JAVASCRIPT( ... )
}//static void gammaEnergyMatcherJs( std::string &js );


void PhotopeakDelegate::EditWidget::replacerJs( std::string &js )
{
  js = INLINE_JAVASCRIPT
  (
      function (edit, suggestionText, suggestionValue)
      {
        edit.value = suggestionValue;

        if (edit.selectionStart)
          edit.selectionStart = edit.selectionEnd = suggestionValue.length;
      }//function (edit, suggestionText, suggestionValue)
   );//js = INLINE_JAVASCRIPT( ... )
}//void replacerJs( std::string &js )


PhotopeakDelegate::EditWidget::EditWidget( const Wt::WModelIndex& index,
              const Wt::WFlags<Wt::ViewItemRenderFlag> flags,
                                          const bool closeOnBlur,
              const PhotopeakDelegate::DelegateType delegateType,
              PhotopeakDelegate *parent )
  : WContainerWidget(),
    m_edit( NULL ),
    m_suggestions( NULL )
{
  char value[64], text[64];
  string matcherJS, replaceJS;

  switch( delegateType )
  {
    case NuclideDelegate:
      nuclideNameMatcherJs( matcherJS );
    break;

    case GammaEnergyDelegate:
      gammaEnergyMatcherJs( matcherJS );
    break;
  }//switch( delegateType )

  replacerJs( replaceJS );

  m_suggestions = new WSuggestionPopup( matcherJS, replaceJS );
#if( WT_VERSION < 0x3070000 ) //I'm not sure what version of Wt "wtNoReparent" went away.
  m_suggestions->setJavaScriptMember("wtNoReparent", "true");
#endif
  
  m_suggestions->setMaximumSize( WLength::Auto, WLength(15, WLength::FontEm) );

  m_edit = new WLineEdit();
  m_edit->setTextSize( 7 );

  m_edit->enterPressed().connect
    (boost::bind(&PhotopeakDelegate::doCloseEditor, parent, this, true, false ));
  m_edit->escapePressed().connect
    (boost::bind(&PhotopeakDelegate::doCloseEditor, parent, this, false, false ));
  m_edit->escapePressed().preventPropagation();
  if( closeOnBlur )
    m_edit->blurred().connect
      (boost::bind(&PhotopeakDelegate::doCloseEditor, parent, this, true, true ));

  if( flags & RenderFocused )
    m_edit->setFocus();

  m_suggestions->forEdit( m_edit,
                   WSuggestionPopup::Editing | WSuggestionPopup::DropDownIcon );

  // We use a layout so that the line edit fills the entire cell.
  // Somehow, this does not work with konqueror, but it does respond
  // properly to width, height being set to 100% !
  WApplication *app = WApplication::instance();
  if( app->environment().agent() != WEnvironment::Konqueror )
  {
    setLayout( new WHBoxLayout() );
    layout()->setContentsMargins(1, 1, 1, 1);
    layout()->addWidget( m_edit );
  }else
  {
    m_edit->resize( WLength(100, WLength::Percentage),
                    WLength(100, WLength::Percentage) );
    addWidget( m_edit );
  }//if( a non-Konqueror browser ) / else

  const PeakModel *peakModel = dynamic_cast<const PeakModel *>( index.model() );
  if( !peakModel )
    return;

  const PeakModel::PeakShrdPtr &peak = peakModel->peak( index );

  if( !peak )
    return;

  switch( delegateType )
  {
    case NuclideDelegate:
    {
      PeakIsotopeNameFilterModel *filterModel = new PeakIsotopeNameFilterModel( peak, this );
//      IsotopeNameFilterModel *filterModel = new IsotopeNameFilterModel( this );
      
      filterModel->filter( "" );
      m_suggestions->setFilterLength( -1 );
      m_suggestions->setModel( filterModel );
      m_suggestions->filterModel().connect( filterModel, &PeakIsotopeNameFilterModel::filter );
//      m_suggestions->filterModel().connect( filterModel, &IsotopeNameFilterModel::filter );
      
      break;
    }//case NuclideDelegate:

    case GammaEnergyDelegate:
    {
      const int row = index.row();
      double fitEnergy = 0.0;
      const PeakModel::PeakShrdPtr &peak = peakModel->peak( peakModel->index(row,0) );
      if( !peak )
        return;
      const SandiaDecay::Nuclide *nuclide = peak->parentNuclide();
      
      /*
       boost::any iso_any = peakModel->data( peakModel->index( row, PeakModel::kIsotope ) );
      boost::any fitEnergy_any = peakModel->data( peakModel->index( row, PeakModel::kMean ) );
      string isotope;
      try { isotope = boost::any_cast<WString>( iso_any ).narrow(); }
      catch(...){ return; }
      try{ fitEnergy = boost::any_cast<double>( fitEnergy_any ); }catch(...){}

      const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
      if( !db )
        return;

      const SandiaDecay::Nuclide *nuclide = db->nuclide( isotope );
       */

      if( nuclide )
      {
        vector< IsotopeSelectionAids::NearGammaInfo > energies
                = IsotopeSelectionAids::equilibriumGammasByNearestEnergy( nuclide, fitEnergy, 1.0E-10, true, true );

        for( const IsotopeSelectionAids::NearGammaInfo &info : energies )
        {
          const char *prefix = "", *postfix = "";
          const double energy = floor(100.0*info.gamma_energy + 0.5)/100.0;
          const double intensity = info.relative_intensity;
          
          switch( info.gamma_type )
          {
            case PeakDef::NormalGamma: case PeakDef::AnnihilationGamma: break;
            case PeakDef::SingleEscapeGamma:          prefix = "S.E. "; break;
            case PeakDef::DoubleEscapeGamma:          prefix = "D.E. "; break;
            case PeakDef::XrayGamma:                  postfix = " xray"; break;
          }//switch( info.gamma_type )
              
          snprintf( value, sizeof(value), "%s%.2f keV%s", prefix, energy, postfix );
          snprintf( text, sizeof(text), "%s%.2f keV%s I=%.1e", prefix, energy, postfix, intensity );
          m_suggestions->addSuggestion( text, value );
        }//for( const IsotopeSelectionAids::NearGammaInfo &info : energies )
      }else if( peak->xrayElement() )
      {
        for( size_t i = 0; i < peak->xrayElement()->xrays.size(); ++i )
        {
          const double energy = floor(100.0*peak->xrayElement()->xrays[i].energy + 0.5)/100.0;
          const double intensity = peak->xrayElement()->xrays[i].intensity;
          snprintf( value, sizeof(value), "%.2f keV", energy );
          snprintf( text, sizeof(text), "%.2f keV I=%.1e", energy, intensity );
          m_suggestions->addSuggestion( text, value );
        }//for( size_t i = 0; i < peak->xrayElement()->xrays.size(); ++i )
      }else if( peak->reaction() )
      {
        for( size_t i = 0; i < peak->reaction()->gammas.size(); ++i )
        {
          const double energy = floor(100.0*peak->reaction()->gammas[i].energy + 0.5)/100.0;
          const double intensity = peak->reaction()->gammas[i].abundance;
          snprintf( value, sizeof(value), "%.2f keV", energy );
          snprintf( text, sizeof(text), "%.2f keV I=%.1e", energy, intensity );
          m_suggestions->addSuggestion( text, value );
        }//for( size_t i = 0; i < peak->xrayElement()->xrays.size(); ++i )        
      }//if( nuclide ) / else if( xray ) / else if( reaction )
    }//case GammaEnergyDelegate:
  }//switch( delegateType )

}//PhotopeakDelegate::EditWidget constructor

PhotopeakDelegate::EditWidget::~EditWidget()
{
  delete m_suggestions;
}

WLineEdit *PhotopeakDelegate::EditWidget::edit()
{
  return m_edit;
}



PhotopeakDelegate::PhotopeakDelegate(  PhotopeakDelegate::DelegateType delegateType,
                                       bool closeOnBlur,
                                       Wt::WObject *parent )
  : WAbstractItemDelegate( parent ),
    m_closeOnBlur( closeOnBlur ),
    m_delegateType( delegateType ),
    m_suggestionPopup( NULL )
{
}//PhotopeakDelegate constructor


PhotopeakDelegate::~PhotopeakDelegate()
{
}//~PhotopeakDelegate()



void PhotopeakDelegate::doCloseEditor( WWidget *editor, bool save, bool isBlurr ) const
{
  PhotopeakDelegate::EditWidget *edit = dynamic_cast<PhotopeakDelegate::EditWidget *>( editor );
  WLineEdit *lineEdit = edit ? edit->edit() : nullptr;

  if( !lineEdit || (isBlurr && lineEdit->text().empty()) )
    save = false;

  closeEditor().emit( editor, save );
}//void doCloseEditor(WWidget *editor, bool save) const



WWidget *PhotopeakDelegate::update( WWidget *widget,
                                    const WModelIndex &index,
                                    WFlags< ViewItemRenderFlag > flags )
{
  bool isNew = true;
  const bool editing = (widget && (widget->find("t") == 0));

  if( flags & RenderEditing )
  {
    if( !editing )
    {
      widget = new EditWidget( index, flags, m_closeOnBlur, m_delegateType, this );
      WInteractWidget *iw = dynamic_cast<WInteractWidget *>( widget );
      if( iw ) // Disable drag & drop and selection behaviour
      {
        iw->mouseWentDown().preventPropagation();
        iw->clicked().preventPropagation();
      }//if( iw )
    }//if( !editing )
  }else
  {
    if( editing )
      widget = 0;
  }//if( flags & RenderEditing ) / lese

  if( !(flags & RenderEditing) )
  {
    WText *text = dynamic_cast<WText *>( widget );

    if( !text )
    {
      isNew = true;
      text = new WText();
      text->setObjectName( "t" );
      if( !index.isValid() || (index.isValid() && !(index.flags() & ItemIsXHTMLText)) )
        text->setTextFormat(PlainText);
      text->setWordWrap(true);
      widget = text;
    }else if( !index.isValid() )
      text->setText( "" );

    if( !index.isValid() )
      return widget;

    text->setText( asString( index.data() ) );
  }//if( !(flags & RenderEditing) )

  WString tooltip = asString( index.data(ToolTipRole) );
  if( !tooltip.empty() || !isNew )
    widget->setToolTip( tooltip );

  WT_USTRING sc = asString( index.data(StyleClassRole) );

  if( flags & RenderSelected )
    sc += WT_USTRING::fromUTF8( " Wt-selected" );

  widget->setStyleClass( sc );

  return widget;
}//WWidget *update(...)


boost::any PhotopeakDelegate::editState( WWidget *editor ) const
{
  WContainerWidget *w = dynamic_cast<WContainerWidget *>(editor);
  if( !w )
  {
    cerr << "PhotopeakDelegate::editState(...)\n\tLogic error - fix me!" << endl;
    return boost::any();
  }//if( !w )

  WLineEdit *lineEdit = dynamic_cast<WLineEdit *>(w->widget(0));

  if( !lineEdit )
  {
    cerr << "PhotopeakDelegate::editState(...)\n\tLogic error - fix me!" << endl;
    return boost::any();
  }//if( !lineEdit )

  return boost::any( lineEdit->text() );
}//boost::any editState( WWidget *editor ) const


void PhotopeakDelegate::setEditState( WWidget *editor,
                                    const boost::any& value ) const
{
  WContainerWidget *w = dynamic_cast<WContainerWidget *>(editor);
  if( !w )
  {
    cerr << "PhotopeakDelegate::setEditState(...)\n\tLogic error - fix me!" << endl;
    return;
  }

  WLineEdit *lineEdit = dynamic_cast<WLineEdit *>(w->widget(0));

  if( !lineEdit )
  {
    cerr << "PhotopeakDelegate::setEditState(...)\n\tLogic error - fix me!" << endl;
    return;
  }

  try
  {
//    lineEdit->setText( boost::any_cast<WT_USTRING>(value) );
    lineEdit->setText( boost::any_cast<WString>(value) );
  }catch(...)
  {
    cerr << "PhotopeakDelegate::setEditState(...)\n\tPossible Logic error - fix me!" << endl;
    lineEdit->setText( "" );
  }//try / catch
}//void setEditState( WWidget *editor, const boost::any& value ) const


void PhotopeakDelegate::setModelData( const boost::any& editState,
                                    WAbstractItemModel *model,
                                    const WModelIndex& index) const
{
  if( model )
    model->setData( index, editState, EditRole );
}//void setModelData(...)




PeakIsotopeNameFilterModel::PeakIsotopeNameFilterModel(
                                  const std::shared_ptr<const PeakDef> &peak,
                                  Wt::WObject *parent  )
  : WAbstractItemModel( parent ),
    m_minHalfLife( 60.0*15.0 ),
    m_filter( "" ),
    m_peak( peak )
{
}


PeakIsotopeNameFilterModel::~PeakIsotopeNameFilterModel()
{
}


WModelIndex PeakIsotopeNameFilterModel::index( int row, int column,
                                               const WModelIndex &parent ) const
{
  if( parent.isValid() || (column!=0) || (row<0) )
    return WModelIndex();

  const int nrows = static_cast<int>( m_displayData.size() );

  if( row >= nrows )
    return WModelIndex();

  return createIndex( row, column, (void *)NULL );
}//WModelIndex index(...)


WModelIndex PeakIsotopeNameFilterModel::parent( const WModelIndex &/*index*/ ) const
{
  return WModelIndex();
}//WModelIndex parent(...)


int PeakIsotopeNameFilterModel::rowCount( const WModelIndex &parent ) const
{
  if( parent.isValid() )
    return 0;
  return static_cast<int>( m_displayData.size() );
}//int rowCount(...)


int PeakIsotopeNameFilterModel::columnCount( const WModelIndex &parent ) const
{
  if( parent.isValid() )
    return 0;
  return 1;
}//int columnCount(...)


WString PeakIsotopeNameFilterModel::displayText( const SandiaDecay::Nuclide *nuclide,
                    double peakMean, PeakDef::SourceGammaType type, int role )
{
  string typePrefix = "", typePostFix = "";
  switch( type )
  {
    case PeakDef::NormalGamma:
    case PeakDef::AnnihilationGamma:
      break;
      
    case PeakDef::SingleEscapeGamma:
      typePrefix = "S.E. ";
      peakMean += 510.99891;
    break;
    case PeakDef::DoubleEscapeGamma:
      typePrefix = "D.E. ";
      peakMean += 2.0*510.99891;
    break;
      
    case PeakDef::XrayGamma:
      typePostFix = " xray";
      break;
  }//switch( type )
  
  
  if( (role==DisplayRole) || (role==ToolTipRole) )
  {
    string descript;
    descript = nuclide->symbol;
    
    const vector< IsotopeSelectionAids::NearGammaInfo > energies
             = IsotopeSelectionAids::equilibriumGammasByNearestEnergy( nuclide,
                                                     peakMean, 1.0E-10, true, true );
    const size_t nlines = energies.size();
    
    if( nlines )
      descript += " (";
    
    const size_t max_lines = (role!=ToolTipRole) ? 3 : 10;
    for( size_t line = 0; line < max_lines && line < nlines; ++line )
    {
      const IsotopeSelectionAids::NearGammaInfo &info = energies[line];
      if( line )
        descript += ", ";
      
      switch( info.gamma_type )
      {
        case PeakDef::NormalGamma:
        case PeakDef::AnnihilationGamma:
        case PeakDef::XrayGamma:
          break;
        case PeakDef::SingleEscapeGamma:       descript += "S.E. "; break;
        case PeakDef::DoubleEscapeGamma:       descript += "D.E. "; break;
      }//switch( energies[line].gamma_type )
      
      const double rounded = (floor(10.0*info.gamma_energy+0.5)/10.0);
      char buffer[64];
      snprintf( buffer, sizeof(buffer), "%.1f keV", rounded );
      descript += buffer;
      
      if( info.gamma_type == PeakDef::XrayGamma )
        descript += " xray";
    }//for( size_t line = 0; line < max_lines && line < nlines; ++line )
    
    if( nlines > max_lines )
      descript += "...";
    if( nlines )
      descript += ")";
    return typePrefix + descript;
  }//if( role == DisplayRole )
  
  if( role == UserRole )
    return typePrefix + nuclide->symbol + typePostFix;
  
  return "";
}//WString displayText( const SandiaDecay::Nuclide *txt, int role )


WString PeakIsotopeNameFilterModel::displayText( const SandiaDecay::Element *element,
                                                 int role, double minHalfLife )
{
  stringstream answer;
  answer << element->symbol;
  
  //if role==UserRole we only want the element symbol, and not all the
  //  possible isotopes so this way if the user selects this option, only
  //  the symbol will be placed into the form
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  if( db && (role==DisplayRole) )
  {
    int niso = 0;
    std::vector<const SandiaDecay::Nuclide *> nuclides = db->nuclides( element );
    for( size_t i = 0; i < nuclides.size(); ++i )
    {
      if( (nuclides[i]->halfLife > minHalfLife)
         && !IsInf( nuclides[i]->halfLife)
         && (niso<5 || role==Wt::ToolTipRole) )
      {
        ++niso;
        answer << " " << nuclides[i]->symbol.substr( element->symbol.length() );
      }else if( niso >= 5 && (role!=Wt::ToolTipRole) )
      {
        answer << "...";
        break;
      }
    }//for( const Nuclide *nuclide : nuclides )
  }//if( db )
  
  if(role == UserRole || role == DisplayRole || role == ToolTipRole )
    return answer.str();
  
  return "";
}//WString displayText( const SandiaDecay::Element *el, int role )


boost::any PeakIsotopeNameFilterModel::data( const Wt::WModelIndex &index, int role ) const
{
  const int row = index.row();
  const int column = index.column();
  const int nrows = static_cast<int>( m_displayData.size() );

  if( row < 0 || column!=0 || row >= nrows  )
    return boost::any();

  if( role == DisplayRole )
    return m_displayData[row];
  else if( role == UserRole )
    return m_userData[row];
  
  return boost::any();
}//boost::any data(..)



int PeakIsotopeNameFilterModel::determineAndRemoveIsoLevel( std::string &label )
{
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


void PeakIsotopeNameFilterModel::getAlphaAndNumericSubStrs( std::string text,
                                       std::vector<std::string> &alphastrs,
                                       std::vector<std::string> &numericstrs )
{
  SpecUtils::erase_any_character( text, "-_,\t<>/?[]{}\\|!@#$%^&*;:\"'~`+=" );
  vector<string> fields;
  SpecUtils::split( fields, text, " ()" );

  for( string &label : fields )
  {
    SpecUtils::to_lower_ascii( label );
    SpecUtils::trim( label );

    if( label.empty() )
      continue;
    
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
  }
}//void PeakIsotopeNameFilterModel::getAlphaAndNumericSubStrs(...)



void PeakIsotopeNameFilterModel::filter( const Wt::WString &text )
{
  //TODO: filter(...) is very poorly coded (it's a mess!) and probably quite
  //      inefficient - this should be cleaned up at some point.
  //      On my MacbookPro 2.3 GHz Intel Core i7 execution time was always
  //      less than 0.001 seconds wall time (debug compilation) on 20120607.

//  boost::timer::cpu_timer timer;

//  static int reminder = 0;
//  if( !reminder++ )
//    cerr << "\nPeakIsotopeNameFilterModel::filter(...): Reminder, you still need to"
//         << " implement using PeakDef::CandidateNuclide::nuclide!\n"
//         << endl;

  //m_minHalfLife
  beginRemoveRows( WModelIndex(), 0, rowCount() );
  m_filter = "";
  m_displayData.clear();
  m_userData.clear();
  endRemoveRows();

  if( !m_peak )
    return;

  string testTxt = text.toUTF8();
  const vector<PeakDef::CandidateNuclide> &candidates = m_peak->candidateNuclides();
  vector<PeakDef::CandidateNuclide> fitleredCandidates;

  PeakDef::SourceGammaType srcType;
  PeakDef::gammaTypeFromUserInput( testTxt, srcType );
  
  
  //If we find a number folowed by a 'm' or "meta" which ends the string or
  //  folowed by something other than a letter, than we have a meta.
  //  Also take into account the possiblity of a 'm2' state.
  //  - remove the m or 'meta', so it wont interfere with later on, mark it as
  //    a meta, and move on.
  //  XXX - TODO - this is a pretty shitty implementation and is messy and
  //               probably bug ridden, should be cleaned up and verified correct.
  //  TODO - Should also check to see if the user is typing in the name of
  //         the element instead of symbol (e.x. 'gold' instead of 'au')
  //         --done for the candidate nuclides
  //  TODO - Even if there are matches in m_suggestedNuclides, should also
  //         match on the symbols or atomic masses of all elements
  //         --done I thinks
  //  TODO - Possibly add in the energy of the nearest photopeak to text of the
  //         isotope name, for user convience
  //  TODO - In addition to sorting results by the number of matching characters,
  //         should first sort my photopeak-closeness as determined by SandiaDecay
  //  TODO - Support x-rays of a specific isotope.


  const int user_input_is_meta = determineAndRemoveIsoLevel( testTxt );

  vector<string> alphastrs, numericstrs;
  getAlphaAndNumericSubStrs( testTxt, alphastrs, numericstrs );


  if( text.empty() )
  {
    fitleredCandidates = candidates;
  }else
  {
    //XXX - this is at best a hack for a proof of concept implementation
    vector<int> fitleredCandidatesNumCharMatch;

    for( const PeakDef::CandidateNuclide &candidate : candidates )
    {
      int nsame = 0;
      bool numeric_candidate = false, alpha_candidate = false;
//      std::string symbol = candidate.transition->parent->symbol;
      std::string symbol = candidate.nuclide->symbol;

      std::string elementName;
      {//begin code block to find element name
        const string::size_type start = symbol.find_first_of( "0123456789" );
        if( start != string::npos )
        {
          std::string elsymbol = symbol.substr( 0, start );
          const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
          if( db )
          {
            const SandiaDecay::Element *elPtr = db->element( elsymbol );
            if( elPtr )
              elementName = elPtr->name;
          }//if( db )
        }//if( start != string::npos )
      }//end code block to find element name

      SpecUtils::to_lower_ascii( symbol );
      for( const string &userstr : alphastrs )
      {
        if( SpecUtils::starts_with( symbol, userstr.c_str() )
            || (elementName.size() && SpecUtils::starts_with( elementName, userstr.c_str() )) )
        {
          for( size_t i = 0; i < userstr.size(); ++i )
            nsame += (int)std::count( symbol.begin(), symbol.end(), userstr[i] );
          alpha_candidate = true;
//          break;
        }//if( starts_with( symbol, userstr ) )
      }//for( const string &userstr : alphastrs )

      string numberstr;
      const string::size_type start = symbol.find_first_of( "0123456789" );
      if( start != string::npos )
      {
        string::size_type end = symbol.find_first_not_of( "0123456789", start );
        if( end == string::npos )
          end = symbol.length();

        numberstr = string( symbol.begin()+start, symbol.begin()+end );
      }//if( start != string::npos )

      if( numberstr.length() )
      {
        for( const string &userstr : numericstrs )
        {
          if( SpecUtils::starts_with( numberstr, userstr.c_str() ) )
          {
            nsame += static_cast<int>( userstr.length() );
            numeric_candidate = true;
            break;
          }//if( starts_with( numberstr, userstr ) )
        }//for( const string &userstr : numericstrs )
      }//if( numberstr.length() )

//      const bool candidate_is_meta = (candidate.transition->parent->isomerNumber > 0);
      const bool candidate_is_meta = (candidate.nuclide->isomerNumber > 0);

      if( candidate_is_meta && user_input_is_meta )
        ++nsame;
      if( candidate.nuclide->isomerNumber==2 && user_input_is_meta==2 )
        ++nsame;


      //if alphastrs and numberstr are both populated, then the nuclide
      //  must be a candidate for both alpha_candidate and numeric_candidate.
      //  If only one of numberstr or alphastrs is populated, than it must
      //  be a candidate for one of them
      if( ((alpha_candidate || numeric_candidate)
          && ( (alphastrs.empty() || numericstrs.empty())
               || (alpha_candidate && numeric_candidate) ))
          || (candidate_is_meta && user_input_is_meta && alphastrs.empty() && numericstrs.empty()) )
      {
        fitleredCandidates.push_back( candidate );
        fitleredCandidatesNumCharMatch.push_back( nsame );
      }//if( good_candidate )
    }//for( const PeakDef::CandidateNuclide &candidate : candidates )

    //Now order the candidates by how many of the characters match
    vector<NSameCharSortHelper> sorterv;
    for( size_t i = 0; i < fitleredCandidates.size(); ++i )
      sorterv.push_back( NSameCharSortHelper( fitleredCandidatesNumCharMatch[i], i ) );

    std::stable_sort( sorterv.begin(), sorterv.end() );

    vector<PeakDef::CandidateNuclide> fitleredCandidatesSorted( sorterv.size() );
    for( size_t i = 0; i < sorterv.size(); ++i )
      fitleredCandidatesSorted[i] = fitleredCandidates[sorterv[i].index];
    fitleredCandidatesSorted.swap( fitleredCandidates );
  }//if( text.empty() )

  vector<const SandiaDecay::Element *> elements;
  vector<const SandiaDecay::Nuclide *> nuclides;

  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  //XXX - TODO - the below is under develoment...
  //       maybe I should get rid of m_elements and place everything into
  //       into m_isotopes
  if( db )
  {
    if( alphastrs.size() && numericstrs.empty() )
    {
      
      const vector<const SandiaDecay::Element *> &dbelements = db->elements();
      for( const string &userstr : alphastrs )
      {
        for( const SandiaDecay::Element *el : dbelements )
          if( (SpecUtils::istarts_with( el->symbol, userstr.c_str() )
               || SpecUtils::istarts_with( el->name, userstr.c_str() )) )
            elements.push_back( el );
      }
    }else if( alphastrs.empty() && numericstrs.size() )
    {
      for( const string &numstr : numericstrs )
      {
        try
        {
          const int massNumber = std::stoi( numstr );
          const vector<const SandiaDecay::Nuclide *> &dbnuclides = db->nuclides();

          for( const SandiaDecay::Nuclide *nuclide : dbnuclides )
          {
            if( nuclide->halfLife > m_minHalfLife
                && nuclide->massNumber == massNumber
                && nuclide->isomerNumber >= user_input_is_meta )
            {
              bool isSuggested = false;
              for( const PeakDef::CandidateNuclide &cand : fitleredCandidates )
//                isSuggested |= (nuclide == cand.transition->parent);
                isSuggested |= (nuclide == cand.nuclide);
              if( !isSuggested )
                nuclides.push_back( nuclide );
            }//if( we have a match )
          }//for( const SandiaDecay::Nuclide *nuclide : nuclides )
        }catch(...){}
      }//for( const string &numstr : numericstrs )
    }else if( alphastrs.size() && numericstrs.size() )
    {
      //XXX - there are _way_ too many nested loops in this section of code!
      const vector<const SandiaDecay::Element *> &dbelements = db->elements();
      for( const string &alphastr : alphastrs )  //typically alphastrs.size()==1
      {
        for( const SandiaDecay::Element *el : dbelements )
        {
          if( (SpecUtils::istarts_with( el->symbol, alphastr.c_str() )
               || SpecUtils::istarts_with( el->name, alphastr.c_str() )) )
          {
            vector<const SandiaDecay::Nuclide *> dbnuclides = db->nuclides( el );

            for( const SandiaDecay::Nuclide *nuclide : dbnuclides )
            {
              if( nuclide->isomerNumber < user_input_is_meta )
                continue;

              for( const string &numstr : numericstrs ) //typically numericstrs.size()==1
              {
                try
                {
                  const string nucstr = std::to_string( nuclide->massNumber );
                  if( SpecUtils::starts_with( nucstr, numstr.c_str() )
                      && (std::find( nuclides.begin(), nuclides.end(), nuclide ) == nuclides.end()) )
                  {
                    bool isSuggested = false;
                    for( const PeakDef::CandidateNuclide &cand : fitleredCandidates )
//                      isSuggested |= (nuclide == cand.transition->parent);
                      isSuggested |= (nuclide == cand.nuclide);

                    if( !isSuggested )
                      nuclides.push_back( nuclide );
                    break;
                  }
                }catch(...){}
              }//for( const string &numstr : numericstrs )
            }//for( const SandiaDecay::Nuclide *nuclide : nuclides )
          }//if( the symbol matches )
        }//for( const SandiaDecay::Element *el : dbelements )
      }//for( const string &userstr : alphastrs )
    }//if / else if / else if ( figure out what the user has typed in )
  }//if( db )

  
  
  vector<const ReactionGamma::Reaction *> suggest_reactions;
  set<const SandiaDecay::Element *> candidate_elements = IsotopeNameFilterModel::possibleElements( alphastrs );
  if( user_input_is_meta==0 )
    suggest_reactions = IsotopeNameFilterModel::suggestReactions( text, nuclides, candidate_elements );

  
//  cout << "filter took: " << timer.format() << endl;

  
  //Sort suggestions by Levenshtein distance, so that way the word closest to
  //  what the user has typed in, will be at the top of the list.  If we dont
  //  do this, then if the user types "Ra226" and hits enter, then "Rn226"
  //  will actually be selected.
  using SpecUtils::levenshtein_distance;
  testTxt = text.narrow();
  
  //Sort suggested nuclides
  vector<unsigned int> distance( nuclides.size() );
  vector<size_t> sort_indices( nuclides.size() );
  for( size_t i = 0; i < nuclides.size(); ++i )
  {
    sort_indices[i] = i;
    const string &symbol = nuclides[i]->symbol;
    distance[i] = levenshtein_distance( testTxt, symbol );
  }//for( size_t i = 0; i < suggestions.size(); ++i )
  std::stable_sort( sort_indices.begin(), sort_indices.end(),
                    index_compare_assend<vector<unsigned int>&>(distance) );
  
  //Sort suggested elements
  distance.resize( elements.size() );
  vector<size_t> element_sort_indices( elements.size() );
  for( size_t i = 0; i < elements.size(); ++i )
  {
    element_sort_indices[i] = i;
    const string &symbol = elements[i]->symbol;
    const string &name = elements[i]->name;
    distance[i] = std::min( levenshtein_distance( testTxt, symbol ),
                           levenshtein_distance( testTxt, name ) );
  }//for( size_t i = 0; i < suggest_elements.size(); ++i )
  std::stable_sort( element_sort_indices.begin(), element_sort_indices.end(),
                    index_compare_assend<vector<unsigned int>&>(distance) );
  
//  int nrow = static_cast<int>( fitleredCandidates.size()
//                                     + nuclides.size() + 2*elements.size() );

//  beginInsertRows( WModelIndex(), 0, nrow );
  m_filter = text;
  const double peakMean = m_peak ? m_peak->mean() : -0.0;
  vector<WString> displayData, userData;
  
  for( const PeakDef::CandidateNuclide &cand : fitleredCandidates )
  {
    displayData.push_back( displayText(cand.nuclide, peakMean, srcType, Wt::DisplayRole) );
    userData.push_back( displayText(cand.nuclide, peakMean, srcType, Wt::UserRole) );
  }
  
  for( size_t i = 0; i < nuclides.size(); ++i )
  {
    const SandiaDecay::Nuclide *nuc = nuclides[ sort_indices[i] ];
    displayData.push_back( displayText(nuc, peakMean, srcType, Wt::DisplayRole) );
    userData.push_back( displayText(nuc, peakMean, srcType, Wt::UserRole) );
  }//for( size_t i = 0; i < nuclides.size(); ++i )
  
  
  //If the user has explicitly type "xray", dont bother to show other suggestions
  if( elements.size() && srcType==PeakDef::XrayGamma )
  {
    displayData.clear();
    userData.clear();
    suggest_reactions.clear();
  }//if( the user probably only wanst reactions suggested )
  
  for( size_t i = 0; i < elements.size(); ++i )
  {
    const SandiaDecay::Element *el = elements[ element_sort_indices[i] ];
    displayData.push_back( displayText(el, Wt::DisplayRole, m_minHalfLife) );
    userData.push_back( displayText(el, Wt::UserRole, m_minHalfLife) );
    
    if( el->xrays.size() )
    {
      displayData.push_back( el->symbol + " xray" );
      userData.push_back( el->symbol + " xray" );
    }//if( el->xrays.size() )
  }
  
  if( suggest_reactions.size() && text.toUTF8().find( '(')!=string::npos )
  {
    displayData.clear();
    userData.clear();
  }//if( the user probably only wanst reactions suggested )
  
  string typePrefix = "";
  switch( srcType )
  {
    case PeakDef::NormalGamma:
    case PeakDef::AnnihilationGamma:
    case PeakDef::XrayGamma:
      break;
    case PeakDef::SingleEscapeGamma: typePrefix = "S.E. ";    break;
    case PeakDef::DoubleEscapeGamma: typePrefix = "D.E. ";    break;
  }//switch( srcType )
  
  for( size_t i = 0; i < suggest_reactions.size(); ++i )
  {
    displayData.push_back( typePrefix + suggest_reactions[i]->name() );
    userData.push_back( typePrefix + suggest_reactions[i]->name() );
  }
  
  const int nrow = static_cast<int>( displayData.size() );
  beginInsertRows( WModelIndex(), 0, nrow -1);
  m_userData.swap( userData );
  m_displayData.swap( displayData );
  endInsertRows();
}//void filter(...)



void PeakIsotopeNameFilterModel::setPeak( const std::shared_ptr<const PeakDef> &peak )
{
  m_peak = peak;
  filter( m_filter );
}//setPeak( const std::shared_ptr<const PeakDef> &peak )








