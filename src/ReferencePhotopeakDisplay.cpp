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

#include <vector>
#include <algorithm>

#include <boost/tuple/tuple.hpp>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WServer>
#include <Wt/WCheckBox>
#include <Wt/WLineEdit>
#include <Wt/WIOService>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WStringStream>
#include <Wt/WDoubleSpinBox>
#include <Wt/WContainerWidget>
#include <Wt/WSuggestionPopup>
#include <Wt/WRegExpValidator>

#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/ColorSelect.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/MoreNuclideInfo.h"
#include "InterSpec/ShieldingSelect.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/FeatureMarkerWidget.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeSelectionAids.h"
#include "InterSpec/IsotopeNameFilterModel.h"
#include "InterSpec/MoreNuclideInfoDisplay.h"
#include "InterSpec/PhysicalUnitsLocalized.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"

using namespace std;
using namespace Wt;

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__


#if( ANDROID )
// Defined in target/android/android.cpp
extern void android_download_workaround( Wt::WResource *resource, std::string description );
#endif

const int ReferencePhotopeakDisplay::sm_xmlSerializationVersion = 0;

const int DecayParticleModel::RowData::XRayDecayMode = 1000;
const int DecayParticleModel::RowData::ReactionToGammaMode = 1001;
const int DecayParticleModel::RowData::NormGammaDecayMode = 1002;
const int DecayParticleModel::RowData::CascadeSumMode = 1003;

namespace
{
  //See: https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
  //  for a pallet of distint colors
  //Or http://artshacker.com/wp-content/uploads/2014/12/Kellys-22-colour-chart.jpg
  const static vector<Wt::WColor> ns_def_line_colors{
    {"#0000FF"}, {"#006600"}, {"#0099FF"}, {"#9933FF"},
    {"#FF66FF"}, {"#CC3333"}, {"#FF6633"}, {"#FFFF99"},
    {"#CCFFCC"}, {"#0000CC"}, {"#666666"}, {"#003333"}
  };
  
  
  //struct on_scope_exit
  //{
  //  on_scope_exit( std::function<void( void )> f ) : m_function( f ) {}
  //  ~on_scope_exit( void ) { m_function(); }
  //private:
  //  std::function<void( void )> m_function;
  //};//on_scope_exit
  
  struct UpdateGuard
  {
    bool &m_guard;
    UpdateGuard( bool &guard ) : m_guard( guard ) { m_guard = true; }
    ~UpdateGuard(){ m_guard = false; }
  };

  class RefGammaCsvResource : public Wt::WResource
  {
  protected:
    ReferencePhotopeakDisplay *m_display;
    Wt::WApplication *m_app;
    
  public:
    RefGammaCsvResource( ReferencePhotopeakDisplay *parent )
    : WResource( parent ),
    m_display( parent ),
    m_app( WApplication::instance() )
    {
      assert( m_app );
      assert( m_display );
    }
    
    virtual ~RefGammaCsvResource()
    {
      beingDeleted();
    }
    
    
  private:
    virtual void handleRequest( const Wt::Http::Request &, Wt::Http::Response &response )
    {
      WApplication::UpdateLock lock( m_app );
      
      if( !lock )
      {
        log("error") << "Failed to WApplication::UpdateLock in RefGammaCsvResource.";
        
        response.out() << "Error grabbing application lock to form RefGammaCsvResource resource; please report to InterSpec@sandia.gov.";
        response.setStatus(500);
        assert( 0 );
        
        return;
      }//if( !lock )
      
    
      if( !m_display )
        return;
      
      const DecayParticleModel *particleModel = m_display->particleModel();
      if( !particleModel )
        return;
      
      char buffer[128] = { '\0' };
      
      const string eol_char = "\r\n"; //for windows - could potentially cosutomize this for the users operating system
      
      
      vector<DecayParticleModel::RowData> row_data = particleModel->rowData();
      
      stable_sort( begin(row_data), end(row_data),
        [](const DecayParticleModel::RowData &lhs,const DecayParticleModel::RowData &rhs) -> bool {
          return lhs.energy < rhs.energy;
      } );
      
      
      // Right now we'll just download the current nuclide/reaction/x-ray, and not any
      //  of the "persisted" lines
      const ReferenceLineInfo &refinfo = m_display->currentlyShowingNuclide();
      
      string filename = refinfo.m_input.m_input_txt;
      
      switch( refinfo.m_source_type )
      {
        case ReferenceLineInfo::SourceType::Nuclide:
          filename += "_lines";
          break;
          
        case ReferenceLineInfo::SourceType::FluorescenceXray:
          filename += "_xrays";
          break;
          
        case ReferenceLineInfo::SourceType::Reaction:
          filename += "_reaction";
          break;
          
        case ReferenceLineInfo::SourceType::Background:
          filename =  "background_lines";
          break;
          
        case ReferenceLineInfo::SourceType::CustomEnergy:
          filename += "_custom_energy";
          break;
          
        case ReferenceLineInfo::SourceType::NuclideMixture:
          filename += "_mixture";
          break;
          
        case ReferenceLineInfo::SourceType::OneOffSrcLines:
          filename += "_one_off_src";
          break;
          
        case ReferenceLineInfo::SourceType::None:
          filename = "empty";
          break;
      }//switch( refinfo.m_source_type )
      
      SpecUtils::ireplace_all( filename, "(", "_" );
      SpecUtils::ireplace_all( filename, ")", "_" );
      SpecUtils::ireplace_all( filename, " ", "_" );
      SpecUtils::ireplace_all( filename, ",", "-" );
      if( !filename.empty() && ((filename.back() == '_') || (filename.back() == '-')) )
         filename = filename.substr( 0, filename.size()-1 );
      
      filename += ".csv";
      suggestFileName( filename, WResource::Attachment );
      response.setMimeType( "text/csv" );
      
      std::ostream &out = response.out();
      
      if( refinfo.m_ref_lines.empty()
         || (refinfo.m_source_type == ReferenceLineInfo::SourceType::None))
      {
        assert( refinfo.m_ref_lines.empty() );
        
        out << "No displayed reference lines to output" << eol_char;
        return;
      }
      
      switch( refinfo.m_source_type )
      {
        case ReferenceLineInfo::SourceType::Nuclide:
        {
          const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
          const SandiaDecay::Nuclide * const nuc = db->nuclide(refinfo.m_input.m_input_txt);
          assert( nuc );
          const string name = nuc ? nuc->symbol : string("null");
          const string hl = nuc ? PhysicalUnits::printToBestTimeUnits(nuc->halfLife,6) : "null";
          out << "Nuclide," << (nuc ? nuc->symbol : string()) << eol_char;
          out << "HalfLife," << hl << eol_char;
          out << "AgeDecayedTo," << refinfo.m_input.m_age;
          if( refinfo.m_input.m_promptLinesOnly )
            out << ",PromptEquilibriumNuclidesOnly";
          out << eol_char;
          break;
        }
          
        case ReferenceLineInfo::SourceType::FluorescenceXray:
          out << "Element," << refinfo.m_input.m_input_txt << eol_char;
          out << "Florescent x-rays" << eol_char;
          break;
          
        case ReferenceLineInfo::SourceType::Reaction:
          out << "Reaction," << refinfo.m_input.m_input_txt << eol_char;
          break;
          
        case ReferenceLineInfo::SourceType::Background:
        case ReferenceLineInfo::SourceType::CustomEnergy:
        case ReferenceLineInfo::SourceType::OneOffSrcLines:
          out << "Source,Gammas" << eol_char;
          break;
          
        case ReferenceLineInfo::SourceType::NuclideMixture:
        {
          out << "NuclideMixture," << refinfo.m_input.m_input_txt << eol_char;
          out << "AgeDecayedTo," << refinfo.m_input.m_age << eol_char;
          break;
        }//
          
        case ReferenceLineInfo::SourceType::None:
          assert( 0 );
          break;
      }//switch( refinfo.m_source_type )
      
      
      std::function<double( float )> att_fcn = refinfo.m_input.m_shielding_att;
      
      try
      {
        if( !refinfo.m_input.m_shielding_name.empty() && !refinfo.m_input.m_shielding_thickness.empty() )
        {
          assert( refinfo.m_input.m_shielding_an.empty() );
          assert( refinfo.m_input.m_shielding_ad.empty() );
          
          const Material *material = nullptr;
          const MaterialDB *matDB = m_display->materialDB();
          if( matDB )
            material = matDB->material( refinfo.m_input.m_shielding_name );
          
          if( !material )
            throw runtime_error( "Invalid shielding '" + refinfo.m_input.m_shielding_name + "'" );
          
          out << "Shielding Material," << material->name << eol_char;
          
          const static double cm3PerG = PhysicalUnits::cm3 / PhysicalUnits::g;
          snprintf( buffer, sizeof(buffer), "%.6g", (material->density * cm3PerG) );
          out << "Shielding Density (g/cm3)," << buffer << eol_char;
          
          out << "Shielding Thickness," << refinfo.m_input.m_shielding_thickness << eol_char;
          
          snprintf( buffer, sizeof(buffer), "%1.6g", material->massWeightedAtomicNumber() );
          out << "Shielding Mass Weighted Atomic Number," << buffer << eol_char;
          
          out << "Shielding Chemical Formula," << material->chemicalFormula() << eol_char;
        }else if( !refinfo.m_input.m_shielding_an.empty() && !refinfo.m_input.m_shielding_ad.empty() )
        {
          assert( refinfo.m_input.m_shielding_name.empty() );
          assert( refinfo.m_input.m_shielding_thickness.empty() );
          out << "Shielding Atomic Number," << refinfo.m_input.m_shielding_an << eol_char;
          out << "Shielding Areal Density (g/cm2)," << refinfo.m_input.m_shielding_ad << eol_char;
        }else
        {
          assert( refinfo.m_input.m_shielding_an.empty() );
          assert( refinfo.m_input.m_shielding_ad.empty() );
          assert( refinfo.m_input.m_shielding_name.empty() );
          
          // We may have a shielding thickness, but no material
          //assert( refinfo.m_input.m_shielding_thickness.empty() );
          
          out << "Shielding,None" << eol_char;
        }//
      }catch( std::exception &e )
      {
        out << "Shielding,None" << eol_char;
      }//try / catch
      
      out << "Detector Response Function (DRF),";
      auto det = refinfo.m_input.m_det_intrinsic_eff;
      if( det && !refinfo.m_input.m_detector_name.empty() )
      {
        string name = refinfo.m_input.m_detector_name;
        SpecUtils::ireplace_all( name, ",", "-" );
        out << name << eol_char;
      }else
      {
        out << "None" << eol_char;
      }
      
      out << eol_char;
      if( (refinfo.m_input.m_showBetas || refinfo.m_input.m_showAlphas)
         && (refinfo.m_input.m_det_intrinsic_eff || refinfo.m_input.m_shielding_att) )
      {
        out << eol_char << "Note,Alphas and Betas do not include effects of shielding or DRF" << eol_char;
      }
      
      const char *rel_amp_note = "Note,The Rel. Amp. column does not include effects of shielding or DRF";
      
      switch( refinfo.m_source_type )
      {
        case ReferenceLineInfo::SourceType::Nuclide:
          out << eol_char
          << "Note,The g/Bq/second column is rate of gammas emitted per becquerel of "
          << refinfo.m_input.m_input_txt << " and does not include effects of shielding or DRF"
          << eol_char
          << eol_char;
          
          out << "Energy (keV),g/Bq/second";
          break;
          
        case ReferenceLineInfo::SourceType::FluorescenceXray:
        case ReferenceLineInfo::SourceType::Reaction:
        case ReferenceLineInfo::SourceType::Background:
        case ReferenceLineInfo::SourceType::CustomEnergy:
        case ReferenceLineInfo::SourceType::NuclideMixture:
        case ReferenceLineInfo::SourceType::OneOffSrcLines:
          out << eol_char << rel_amp_note << eol_char << eol_char;
          out << "Energy (keV),Rel. Yield";
          break;
          
        case ReferenceLineInfo::SourceType::None:
          assert( 0 );
          break;
      }//switch( refinfo.m_source_type )
      
      out << ",Parent,Mode,Particle";
      
      if( att_fcn )
      {
        out << ",Shielding Transmission";
        if( !det )
          out << ",Yield*ShieldTrans";
      }
      
      if( det )
      {
        out << ",DRF Intrinsic Efficiency";
        if( !att_fcn )
          out << ",Yield*DRF";
      }
      
      if( att_fcn && det )
        out << ",Yield*ShieldTrans*DRF";
      
      out << ",Normalized Intensity" << eol_char;
      
      for( const ReferenceLineInfo::RefLine &line : refinfo.m_ref_lines )
      {
        snprintf( buffer, sizeof(buffer), "%.3f,%1.6g", line.m_energy, line.m_decay_intensity );
        out << buffer;
        
        if( line.m_parent_nuclide
           && (refinfo.m_source_type == ReferenceLineInfo::SourceType::NuclideMixture) )
        {
          out << "," << line.m_parent_nuclide->symbol;
        }else
        {
          if( line.m_transition && line.m_transition->parent )
            out << "," << line.m_transition->parent->symbol;
          else
            out << ",";
        }
        
        out << ",";
        switch( line.m_source_type )
        {
          case ReferenceLineInfo::RefLine::RefGammaType::Normal:
            break;
          case ReferenceLineInfo::RefLine::RefGammaType::Annihilation:
            out << "Annihilation";
            break;
          case ReferenceLineInfo::RefLine::RefGammaType::SingleEscape:
            out << "Single Escape";
            break;
          case ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape:
            out << "Double Escape";
            break;
          case ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak:
            out << "cascade sum";
            break;
          case ReferenceLineInfo::RefLine::RefGammaType::SumGammaPeak:
            out << "sum";
            break;
        }//switch( line.m_source_type )
        
        if( line.m_transition )
        {
          out << " ";
          if( line.m_transition->parent )
            out << line.m_transition->parent->symbol;
          
          if( line.m_transition->child )
            out << " to " << line.m_transition->child->symbol;
          out << " via " << SandiaDecay::to_str(line.m_transition->mode);
        }//if( line.m_transition )
        
        out << ",";
        switch( line.m_particle_type )
        {
          case ReferenceLineInfo::RefLine::Particle::Alpha:
            out << "alpha";
            break;
            
          case ReferenceLineInfo::RefLine::Particle::Beta:
            out << "beta-";
            break;
            
          case ReferenceLineInfo::RefLine::Particle::Gamma:
            out << "gamma";
            break;
            
          case ReferenceLineInfo::RefLine::Particle::Xray:
            out << "xray";
            break;
        }//switch( dataRow.particle )
        
        double shield_eff = 1.0, drf_eff = 1.0;
        if( att_fcn )
        {
          shield_eff = att_fcn( line.m_energy );
          snprintf( buffer, sizeof(buffer), "%1.7g", shield_eff );
          out << "," << buffer;
        }//if( att_coef_fcn )
        
        
        if( det )
        {
          drf_eff = det( line.m_energy );
          snprintf( buffer, sizeof(buffer), "%1.7g", drf_eff );
          out << "," << buffer;
        }//if( det )
        
        if( att_fcn || det )
        {
          snprintf( buffer, sizeof(buffer), "%1.7g", line.m_decay_intensity*shield_eff*drf_eff );
          out << "," << buffer;
        }
        
        snprintf( buffer, sizeof(buffer), "%1.7g", line.m_normalized_intensity );
        out << "," << buffer;
        
        out << eol_char;
      }//for( const ReferenceLineInfo::RefLine &line : refinfo.m_ref_lines )
      
      out << eol_char;
    }//handleRequest(...)
    
  };//class RefGammaCsvResource
  
}//namespace

bool DecayParticleModel::less_than( const DecayParticleModel::RowData &lhs,
                                    const DecayParticleModel::RowData &rhs,
                                    const DecayParticleModel::Column c,
                                    const SortOrder order )
{

  bool less = false;

  switch( c )
  {
    case kEnergy:         less = (lhs.energy < rhs.energy);           break;
    case kBranchingRatio: less = (lhs.branchRatio < rhs.branchRatio); break;
    case kResponsibleNuc:
      if( !lhs.responsibleNuc || !rhs.responsibleNuc )
        less = false;
      else
        less = SandiaDecay::Nuclide::lessThan( lhs.responsibleNuc, rhs.responsibleNuc );
    break;
    case kDecayMode:      less = (lhs.decayMode < rhs.decayMode);     break;
    case kParticleType:   less = (lhs.particle < rhs.particle);       break;
    case kNumColumn:      return false;
  }//switch( r )

  if( order == Wt::AscendingOrder )
    return less;
  return !less;
}//less_than(...)


DecayParticleModel::DecayParticleModel( Wt::WObject *parent )
  : WAbstractItemModel( parent ),
    m_sortColumn( DecayParticleModel::kEnergy ),
    m_sortOrder( Wt::AscendingOrder )
{
}

DecayParticleModel::~DecayParticleModel()
{
}


WFlags<ItemFlag> DecayParticleModel::flags( const WModelIndex &p ) const
{
  if( p.isValid() && (p.column()==kDecayMode || p.column()==kParticleType) )
    return ItemIsXHTMLText;
  return WFlags<ItemFlag>();
}


boost::any DecayParticleModel::data( const WModelIndex &index, int role ) const
{
  using namespace SandiaDecay;

  if( role != DisplayRole )
    return boost::any();

  const int row = index.row();
  const int nrow = static_cast<int>( m_data.size() );
  if( row < 0 || row >= nrow )
    return boost::any();

  const int column = index.column();
  if( column < 0 || column >= kNumColumn )
    return boost::any();

  const RowData &dataRow = m_data[row];
  
  char buffer[64];
  
  switch( column )
  {
    case kEnergy:
      snprintf( buffer, sizeof(buffer), "%.4f", dataRow.energy );
      return WString( buffer );
      
    case kBranchingRatio:
      snprintf( buffer, sizeof(buffer), "%.4g", dataRow.branchRatio );
      return WString( buffer );

    case kResponsibleNuc:
      if( dataRow.responsibleNuc )
        return WString( dataRow.responsibleNuc->symbol );
      return boost::any();

    case kDecayMode:
    {
      switch( dataRow.decayMode )
      {
#ifndef WT_NO_STD_WSTRING
        case SandiaDecay::AlphaDecay:               return WString( L"\x03B1" );
        case SandiaDecay::BetaDecay:                return WString( L"\x03B2<sup>-</sup>" );
        case SandiaDecay::BetaPlusDecay:            return WString( L"\x03B2<sup>+</sup>" );
        case SandiaDecay::DoubleBetaDecay:          return WString( L"double \x03B2" );
#else
        case SandiaDecay::AlphaDecay:               return WString( "&alpha;" );
        case SandiaDecay::BetaDecay:                return WString( "&beta;<sup>-</sup>" );
        case SandiaDecay::BetaPlusDecay:            return WString( "&beta;<sup>+</sup>" );
        case SandiaDecay::DoubleBetaDecay:          return WString( "double &beta;" );
#endif
        case SandiaDecay::IsometricTransitionDecay: return WString::tr("rpd-tbl-iso");
        case SandiaDecay::ElectronCaptureDecay:     return WString::tr( "rpd-tbl-el-cap" );
        case SandiaDecay::ProtonDecay:              return WString::tr( "rpd-tbl-proton" );
        case SandiaDecay::SpontaneousFissionDecay:  return WString::tr( "rpd-tbl-spont-fis" );
        case SandiaDecay::Carbon14Decay:            return WString( "C14" );
        case RowData::XRayDecayMode:                return WString::tr( "rpd-tbl-xray" );
        case RowData::ReactionToGammaMode:          return WString::tr( "rpd-tbl-reaction" );
        case RowData::NormGammaDecayMode:           return WString::tr( "rpd-tbl-norm" );
        case RowData::CascadeSumMode:               return WString::tr( "rpd-tbl-cascade-sum" );
      }//switch( dataRow.decayMode )

      return boost::any();
    }//case kDecayMode:

    case kParticleType:
    {
      switch( dataRow.particle )
      {
#ifndef WT_NO_STD_WSTRING
        case BetaParticle:   return WString( L"\x03B2<sup>-</sup>" );
        case GammaParticle:  return WString( L"\x0263" );
        case AlphaParticle:          return WString( L"\x03B1" );
#else
        case BetaParticle:   return WString( "&beta<sup>-</sup>" );
        case GammaParticle:  return WString( "&gamma;" );
        case AlphaParticle:  return WString( "&alpha;" );
#endif
        case PositronParticle:        return WString( "e<sup>+</sup>" );
        case CaptureElectronParticle: return WString( "ec" );
        case XrayParticle:            return WString::tr( "rpd-tbl-xray" );
      }//switch( dataRow.particle )
      return boost::any();
    }//case kParticleType:
  }//switch( column )

  return boost::any();
}//boost::any data( const Wt::WModelIndex &index, int role  )


boost::any DecayParticleModel::headerData( int column,
                                           Orientation orientation,
                                           int role ) const
{
  if( role == LevelRole )
    return 0;

  if( (orientation != Horizontal)
      || ((role != DisplayRole) && (role != ToolTipRole)) )
    return WAbstractItemModel::headerData( column, orientation, role );

  //If we are here, we want the column title
  if( role == DisplayRole )
  {
    switch( column )
    {
      case kEnergy:         return WString::tr( "rpd-tbl-hdr-energy" );
      case kBranchingRatio: return WString::tr( "rpd-tbl-hdr-br" ); //"Phot/Decay" &gamma;/Decay  \u03B3/Decay  &#947;/Decay
      case kResponsibleNuc: return WString::tr( "rpd-tbl-hdr-parent" );
      case kDecayMode:      return WString::tr( "rpd-tbl-hdr-mode" );
      case kParticleType:   return WString::tr( "rpd-tbl-hdr-particle" );
      case kNumColumn:      return boost::any();
    }//switch( column )
  }else if( role == ToolTipRole )
  {
    switch( column )
    {
      case kEnergy:
        return WString::tr( "rpd-tbl-hdr-tt-energy" );
      case kBranchingRatio:
        return WString::tr( "rpd-tbl-hdr-tt-br" );
      case kResponsibleNuc:
        return WString::tr( "rpd-tbl-hdr-tt-trans" );
      case kDecayMode:
        return WString::tr( "rpd-tbl-hdr-tt-trans-mode" );
      case kParticleType:
        return WString::tr( "rpd-tbl-hdr-tt-type" );
      case kNumColumn:
      return boost::any();
    }//switch( column )
  }//if( role == DisplayRole ) / else

  return boost::any();
}//headerData(...)


int DecayParticleModel::columnCount( const Wt::WModelIndex &parent ) const
{
  if( parent.isValid() )
    return 0;

  return kNumColumn;
}


int DecayParticleModel::rowCount( const Wt::WModelIndex &parent ) const
{
  if( parent.isValid() )
    return 0;
  return static_cast<int>( m_data.size() );
}


WModelIndex DecayParticleModel::parent( const WModelIndex & ) const
{
  return WModelIndex();
}


WModelIndex DecayParticleModel::index( int row, int column,
                                       const WModelIndex &parent ) const
{
  if( parent.isValid() )
    return WModelIndex();

  if( row >= rowCount() || column >= columnCount() )
    return WModelIndex();

  return createIndex( row, column, (void *)NULL );
}//WModelIndex index( int row, int column, const WModelIndex & ) const


void DecayParticleModel::sort( int column, Wt::SortOrder order )
{
  m_sortOrder = order;
  m_sortColumn = Column(column);

  if( m_data.empty() )
    return;
  
  vector<RowData> data = m_data;

  boost::function<bool(const RowData&,const RowData&)> comparer;
  comparer = boost::bind(&DecayParticleModel::less_than, boost::placeholders::_1,
                         boost::placeholders::_2, m_sortColumn, m_sortOrder);

  stable_sort( data.begin(), data.end(), comparer );
  m_data.swap( data );
  dataChanged().emit( index(0,0), index(rowCount()-1,columnCount()-1) );
//  clear();
//  if( data.size() )
//  {
//    beginInsertRows( WModelIndex(), 0, static_cast<int>(data.size())-1 );
//    m_data = data;
//    endInsertRows();
//  }//if( data.size() )
}//sort(...)

void DecayParticleModel::clear()
{
  const int lastRow = rowCount();
  if( lastRow <= 0 )
    return;
  
  beginRemoveRows( WModelIndex(), 0, lastRow -1);
  m_data.clear();
  endRemoveRows();
}//clear()



void DecayParticleModel::setRowData( const std::vector<RowData> &newData )
{
  if( m_data.size() )
    clear();
  beginInsertRows( WModelIndex(), 0, static_cast<int>(newData.size()) -1);
  m_data = newData;
  endInsertRows();

  sort( m_sortColumn, m_sortOrder );
  //reset();
  layoutAboutToBeChanged().emit();
  layoutChanged().emit();
}//setRowData(...)


const std::vector<DecayParticleModel::RowData> &DecayParticleModel::rowData() const
{
  return m_data;
}


ReferencePhotopeakDisplay::ReferencePhotopeakDisplay(
                                            D3SpectrumDisplayDiv *chart,
                                            MaterialDB *materialDB,
                                            WSuggestionPopup *materialSuggest,
                                            InterSpec *specViewer,
                                            WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_chart( chart ),
    m_spectrumViewer( specViewer ),
    m_currently_updating( false ),
    m_undo_redo_sentry(),
    m_nuclideEdit( NULL ),
    m_nuclideSuggest( NULL ),
    m_ageEdit( NULL ),
    m_lowerBrCuttoff( NULL ),
    m_promptLinesOnly( NULL ),
    m_halflife( NULL ),
    m_moreInfoBtn( NULL ),
    m_persistLines( NULL ),
    m_clearLines( NULL ),
    //m_fitPeaks( NULL ),
    m_showGammas( NULL ),
    m_options_icon( NULL ),
    m_options( NULL ),
    m_optionsContent( NULL ),
    m_showXrays( NULL ),
    m_showAlphas( NULL ),
    m_showBetas( NULL ),
    m_showCascadeSums( NULL ),
    m_showEscapes( NULL ),
    m_cascadeWarn( NULL ),
    m_showRiidNucs( NULL ),
    m_showPrevNucs( NULL ),
    m_showAssocNucs( NULL ),
    m_showFeatureMarkers( nullptr ),
    m_otherNucsColumn( nullptr ),
    m_otherNucs( nullptr ),
    m_prevNucs{},
    m_external_ids{},
    m_featureMarkerColumn( nullptr ),
    m_detectorDisplay( NULL ),
    m_materialDB( materialDB ),
    m_materialSuggest( materialSuggest ),
    m_shieldingSelect( NULL ),
    m_particleView( NULL ),
    m_particleModel( NULL ),
    m_currentlyShowingNuclide(),
    m_colorSelect( nullptr ),
    m_csvDownload( nullptr ),
    m_userHasPickedColor( false ),
    m_peaksGetAssignedRefLineColor( false ),
    m_lineColors{ ns_def_line_colors },
    m_specificSourcelineColors{},
    m_displayingNuclide( this ),
    m_nuclidesCleared( this ),
    m_nucInfoWindow( nullptr ),
    m_featureMarkers( nullptr )
{
  auto app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  assert( app );
  if( app )
  {
    app->useMessageResourceBundle( "ReferencePhotopeakDisplay" );
    app->useStyleSheet("InterSpec_resources/ReferencePhotopeakDisplay.css");
  }//if( app )

  m_currentlyShowingNuclide.reset();
  
  if( !chart )
    throw runtime_error( "ReferencePhotopeakDisplay: a valid chart must be passed in" );

  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", specViewer );
  
  //The inputDiv/Layout is the left side of the widget that holds all the
  //  nuclide input,age, color picker, DRF, etc
  WContainerWidget *inputDiv = new WContainerWidget();
  WGridLayout *inputLayout = new WGridLayout();
  inputLayout->setContentsMargins( 0, 2, 0, 0 );
  inputDiv->setLayout( inputLayout );
    
    
  addStyleClass( "ReferencePhotopeakDisplay" );
  
  const bool isPhone = m_spectrumViewer->isPhone();
  if( isPhone )
    addStyleClass( "RefDispMobile" );
  
  const WLength labelWidth(3.5,WLength::FontEm), fieldWidth(4,WLength::FontEm);
  
  WLabel *nucInputLabel = new WLabel( WString("{1}:").arg( WString::tr("Nuclide") ) );
  nucInputLabel->setMinimumSize( labelWidth, WLength::Auto );
  m_nuclideEdit = new WLineEdit( "" );
  m_nuclideEdit->setMargin( 1 );
  m_nuclideEdit->setMargin( 2, Wt::Side::Top );
//  m_nuclideEdit->setMinimumSize( WLength(10,WLength::FontEx), WLength::Auto );
  m_nuclideEdit->setMinimumSize( fieldWidth, WLength::Auto );
  
  m_nuclideEdit->setAutoComplete( false );
  m_nuclideEdit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_nuclideEdit->setAttributeValue( "autocorrect", "off" );
  m_nuclideEdit->setAttributeValue( "spellcheck", "off" );
#endif
  nucInputLabel->setBuddy( m_nuclideEdit );
  
  
//  m_nuclideEdit->changed().connect( boost::bind( &ReferencePhotopeakDisplay::handleIsotopeChange, this, false ) );
//  m_nuclideEdit->blurred().connect( boost::bind( &ReferencePhotopeakDisplay::handleIsotopeChange, this, false ) );
//  m_nuclideEdit->enterPressed().connect( boost::bind( &ReferencePhotopeakDisplay::handleIsotopeChange, this, false ) );
  m_nuclideEdit->changed().connect( boost::bind( &ReferencePhotopeakDisplay::handleIsotopeChange, this, false ) );
  
//  m_nuclideEdit->selected().connect( boost::bind( &ReferencePhotopeakDisplay::handleIsotopeChange, this, false ) );
  m_persistLines = new WPushButton( WString::tr("rpd-add-another-btn") );
  HelpSystem::attachToolTipOn( m_persistLines, WString::tr("rpd-tt-add-another-btn"), showToolTips );
  m_persistLines->clicked().connect( this, &ReferencePhotopeakDisplay::persistCurentLines );
  m_persistLines->disable();
  
  inputLayout->addWidget( nucInputLabel, 0, 0, AlignMiddle );
  inputLayout->addWidget( m_nuclideEdit, 0, 1 );
  inputLayout->addWidget( m_persistLines, 0, 2 );
  
  
  // If we are typing in this box, we want to let app-hotkeys propogate up, but not arrow keys and
  //  stuff
  const string jsAppKeyDownFcn = wApp->javaScriptClass() + ".appKeyDown";
  const string keyDownJs = "function(s1,e1){"
  "if(e1 && e1.ctrlKey && e1.key && " + jsAppKeyDownFcn + ")"
    + jsAppKeyDownFcn + "(e1);"
  "}";
  m_nuclideEdit->keyWentDown().connect( keyDownJs );
  
  HelpSystem::attachToolTipOn( m_nuclideEdit, WString::tr("rpd-tt-nuc-edit"), showToolTips );
  
  string replacerJs, matcherJs;
  IsotopeNameFilterModel::replacerJs( replacerJs );
  IsotopeNameFilterModel::nuclideNameMatcherJs( matcherJs );
  IsotopeNameFilterModel *isoSuggestModel = new IsotopeNameFilterModel( this );
  isoSuggestModel->addCustomSuggestPossibility( "background" );
  
  for( const string &name : ReferenceLineInfo::additional_ref_line_sources() )
    isoSuggestModel->addCustomSuggestPossibility( name );
  
  m_nuclideSuggest = new WSuggestionPopup( matcherJs, replacerJs, this );
#if( WT_VERSION < 0x3070000 ) //I'm not sure what version of Wt "wtNoReparent" went away.
  m_nuclideSuggest->setJavaScriptMember("wtNoReparent", "true");
#endif
  m_nuclideSuggest->setMaximumSize( WLength::Auto, WLength(15, WLength::FontEm) );
  m_nuclideSuggest->setWidth( WLength(70, Wt::WLength::Unit::Pixel) );

  IsotopeNameFilterModel::setQuickTypeFixHackjs( m_nuclideSuggest );
  
  isoSuggestModel->filter( "" );
  m_nuclideSuggest->setFilterLength( -1 );
  m_nuclideSuggest->setModel( isoSuggestModel );
  m_nuclideSuggest->filterModel().connect( isoSuggestModel, &IsotopeNameFilterModel::filter );
  m_nuclideSuggest->forEdit( m_nuclideEdit, WSuggestionPopup::Editing );  // | WSuggestionPopup::DropDownIcon


  WLabel *ageInputLabel = new WLabel( WString("{1}:").arg( WString::tr("Age") ) );
  m_ageEdit = new WLineEdit( "" );
  WRegExpValidator *validator = new WRegExpValidator( PhysicalUnitsLocalized::timeDurationHalfLiveOptionalRegex(), this );
  validator->setFlags(Wt::MatchCaseInsensitive);
  m_ageEdit->setValidator(validator);
  
  m_ageEdit->setAutoComplete( false );
  m_ageEdit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_ageEdit->setAttributeValue( "autocorrect", "off" );
  m_ageEdit->setAttributeValue( "spellcheck", "off" );
#endif
  ageInputLabel->setBuddy( m_ageEdit );

  m_ageEdit->changed().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  m_ageEdit->blurred().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  m_ageEdit->enterPressed().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  
  
  //Well make the "Clear All" button a little bit wider for phones so that
  //  "Gammas" will be on same line as its check box, a little bit hacky
  if( specViewer->isPhone() )
    m_clearLines = new WPushButton( WString::tr("Remove") );
    else
      m_clearLines = new WPushButton( WString::tr("Clear") );
      m_clearLines->disable();
      
      if( specViewer->isMobile() )
      {
        //Get rid of the software keyboard on mobile devices.  It would be nice to
        //  simulate a screen tap on android devices as well, to get rid of system
        //  navigation UI
        m_nuclideEdit->enterPressed().connect( boost::bind( &WPushButton::setFocus, m_clearLines, true ) );
      }
  
  m_clearLines->clicked().connect( this, &ReferencePhotopeakDisplay::clearAllLines );
  
  inputLayout->addWidget( ageInputLabel, 1, 0, AlignMiddle );
  inputLayout->addWidget( m_ageEdit, 1, 1 );
  inputLayout->addWidget( m_clearLines, 1, 2 );
  
  HelpSystem::attachToolTipOn( m_ageEdit, WString::tr("rpd-tt-age"), showToolTips );
  HelpSystem::attachToolTipOn( m_clearLines, WString::tr("rpd-tt-clear"), showToolTips );
  
  
  //If we use a single layout for all the input elements, it seems when we enter
  //  a shielding, then the "Add Another" button will be shortened and layout
  //  messed up.  This must somehow be the layout for the shielding widget
  //  interacting with the layout for this element, however a quick attempt to
  //  fix didnt yield results. Sooo, the solution to this looks to be to have
  //  enough layouts, such that none of them need any spanning of columns.
  WContainerWidget *lowerInput = new WContainerWidget();
  WGridLayout *lowerInputLayout = new WGridLayout();
  lowerInputLayout->setContentsMargins(0, 0, 0, 0);
  lowerInput->setLayout( lowerInputLayout );
  
  WContainerWidget *hlRow = new WContainerWidget();
  hlRow->addStyleClass("HlOptRow");
  lowerInputLayout->addWidget( hlRow, 0, 0 );
  
  m_halflife = new WText( hlRow );
  m_halflife->addStyleClass("Hl");

  m_moreInfoBtn = new WPushButton( WString::tr("rpd-more-info"), hlRow );
  m_moreInfoBtn->addStyleClass( "LinkBtn MoreInfoBtn" );
  m_moreInfoBtn->clicked().connect( this, &ReferencePhotopeakDisplay::showMoreInfoWindow );
  m_moreInfoBtn->hide();

  // Couls add prompt and 'bare' only lines options
  //label = new WLabel( "Lowest I:" );
  //m_lowerBrCuttoff = new WDoubleSpinBox();
  //HelpSystem::attachToolTipOn( m_lowerBrCuttoff, tooltip, showToolTips );
  //m_lowerBrCuttoff->setValue( 0.0 );
  //m_lowerBrCuttoff->setSingleStep( 0.01 );
  //m_lowerBrCuttoff->setRange( 0.0, 1.0 );
  //m_lowerBrCuttoff->valueChanged().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  //m_layout->addWidget( minAmpLabel, 3, 0, AlignMiddle );
  //m_layout->addWidget( m_lowerBrCuttoff, 3, 1 );
  
  SpectraFileModel *specFileModel = specViewer->fileManager()->model();
  m_detectorDisplay = new DetectorDisplay( specViewer, specFileModel );
  
  specViewer->detectorChanged().connect( boost::bind( &ReferencePhotopeakDisplay::handleDrfChange, this, boost::placeholders::_1 ) );
  specViewer->detectorModified().connect( boost::bind( &ReferencePhotopeakDisplay::handleDrfChange, this, boost::placeholders::_1 ) );
  
  // If foreground spectrum _file_ changes, then the RIID analysis results of the suggested
  //  nuclides may need updating.  However, for simplicity, we'll update suggested nuclides
  //  whenever the foreground gets updated
  specViewer->displayedSpectrumChanged().connect(this, &ReferencePhotopeakDisplay::handleSpectrumChange);

  lowerInputLayout->addWidget( m_detectorDisplay, 1, 0 );

  m_shieldingSelect = new ShieldingSelect( m_materialDB, m_materialSuggest );
  m_shieldingSelect->materialEdit()->setEmptyText( WString("<{1}>").arg( WString::tr("rpd-shield-mat") ) );
  m_shieldingSelect->materialChanged().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  m_shieldingSelect->materialModified().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  lowerInputLayout->addWidget( m_shieldingSelect, 2, 0 );

  //m_fitPeaks = new WPushButton( "Fit Peaks" );
  //tooltip = "Fits dominant peaks for primary nuclide";
  //HelpSystem::attachToolTipOn( m_fitPeaks, tooltip, showToolTips );
  //m_fitPeaks->clicked().connect( this, &ReferencePhotopeakDisplay::fitPeaks );
  //inputLayout->addWidget( m_fitPeaks, 2, 2 );
  //m_fitPeaks->disable();
  
  if( m_lineColors.empty() )
    m_lineColors = ns_def_line_colors;
    
  m_colorSelect = new ColorSelect(ColorSelect::PrefferNative);
  m_colorSelect->setColor( m_lineColors[0] );
  
  if( ColorSelect::willUseNativeColorPicker() )
  {
    //m_colorSelect->setFloatSide( Wt::Right );
    m_colorSelect->addStyleClass("RefLinColorPick");
    hlRow->addWidget( m_colorSelect );
  }else
  {
    WContainerWidget *w = new WContainerWidget();
    w->addWidget( m_colorSelect );
    hlRow->addWidget( w );
    //w->setFloatSide( Wt::Right );
    w->addStyleClass("RefLinColorPick");
  }
  
  m_colorSelect->cssColorChanged().connect( boost::bind(
                                               &ReferencePhotopeakDisplay::userColorSelectCallback,
                                               this, boost::placeholders::_1 ) );
  m_currentlyShowingNuclide.m_input.m_color = m_lineColors[0];
  
  m_options_icon = new WPushButton();
  m_options_icon->setStyleClass("RoundMenuIcon InvertInDark RefLinesOptMenu");
  m_options_icon->clicked().preventPropagation();
  //m_options_icon->setFloatSide(Wt::Right);
  m_options_icon->clicked().connect(this, &ReferencePhotopeakDisplay::toggleShowOptions);

  hlRow->addWidget(m_options_icon);

  m_options = new WContainerWidget();
  m_options->addStyleClass("RefLinesOptions ToolTabSection ToolTabTitledColumn");
  m_options->hide();

  WContainerWidget* closerow = new WContainerWidget(m_options);
  closerow->addStyleClass( "ToolTabColumnTitle" );
  WText *txt = new WText( WString::tr("rpd-options"), closerow );
  WContainerWidget* closeIcon = new WContainerWidget(closerow);
  closeIcon->addStyleClass("closeicon-wtdefault");
  closeIcon->clicked().connect(this, &ReferencePhotopeakDisplay::toggleShowOptions);

  m_optionsContent = new WContainerWidget( m_options );
  m_optionsContent->addStyleClass( "ToolTabTitledColumnContent" );

  m_promptLinesOnly = new WCheckBox( WString::tr("rpd-opt-prompt"), m_optionsContent );  //ɣ
  HelpSystem::attachToolTipOn(m_promptLinesOnly, WString::tr("rpd-opt-tt-prompt"), showToolTips);
  m_promptLinesOnly->checked().connect(this, &ReferencePhotopeakDisplay::updateDisplayChange);
  m_promptLinesOnly->unChecked().connect(this, &ReferencePhotopeakDisplay::updateDisplayChange);
  m_promptLinesOnly->hide();

  m_showGammas = new WCheckBox( WString::tr("rpd-opt-gamma"), m_optionsContent );
  m_showXrays = new WCheckBox( WString::tr("rpd-opt-xray"), m_optionsContent );
  m_showAlphas = new WCheckBox( WString::tr("rpd-opt-alphas"), m_optionsContent );
  m_showBetas = new WCheckBox( WString::tr("rpd-opt-betas"), m_optionsContent );
  m_showCascadeSums = new WCheckBox( WString::tr("rpd-opt-cascade"), m_optionsContent );
  m_showCascadeSums->hide();
  m_showEscapes = new WCheckBox( WString::tr("rpd-opt-escapes"), m_optionsContent );
      
  m_showPrevNucs = new WCheckBox( WString::tr("rpd-prev-nucs"), m_optionsContent );
  m_showRiidNucs = new WCheckBox( WString::tr("rpd-det-nucs"), m_optionsContent );
  m_showAssocNucs = new WCheckBox( WString::tr("rpd-assoc-nucs"), m_optionsContent );
  m_showFeatureMarkers = new WCheckBox( WString::tr("rpd-feature-markers"), m_optionsContent );
      
  m_showGammas->setWordWrap( false );
  m_showXrays->setWordWrap( false );
  m_showAlphas->setWordWrap( false );
  m_showBetas->setWordWrap( false );
  m_showCascadeSums->setWordWrap( false );
  m_showEscapes->setWordWrap( false );
  m_showPrevNucs->setWordWrap( false );
  m_showRiidNucs->setWordWrap( false );
  m_showAssocNucs->setWordWrap( false );
  m_showFeatureMarkers->setWordWrap( false );

  m_showPrevNucs->checked().connect( this, &ReferencePhotopeakDisplay::updateOtherNucsDisplay );
  m_showPrevNucs->unChecked().connect( this, &ReferencePhotopeakDisplay::updateOtherNucsDisplay );
  m_showRiidNucs->checked().connect( this, &ReferencePhotopeakDisplay::updateOtherNucsDisplay );
  m_showRiidNucs->unChecked().connect( this, &ReferencePhotopeakDisplay::updateOtherNucsDisplay );
  m_showAssocNucs->checked().connect( this, &ReferencePhotopeakDisplay::updateOtherNucsDisplay );
  m_showAssocNucs->unChecked().connect( this, &ReferencePhotopeakDisplay::updateOtherNucsDisplay );
  m_showFeatureMarkers->checked().connect( this, &ReferencePhotopeakDisplay::featureMarkerCbToggled );
  m_showFeatureMarkers->unChecked().connect( this, &ReferencePhotopeakDisplay::featureMarkerCbToggled );
      

  //const bool showToolTips = InterSpecUser::preferenceValue<bool>("ShowTooltips", this);
  //HelpSystem::attachToolTipOn(m_showPrevNucs, "Show ", showToolTips);
  InterSpecUser::associateWidget( specViewer->m_user, "RefLineShowPrev", m_showPrevNucs, specViewer );
  InterSpecUser::associateWidget( specViewer->m_user, "RefLineShowRiid", m_showRiidNucs, specViewer );
  InterSpecUser::associateWidget( specViewer->m_user, "RefLineShowAssoc", m_showAssocNucs, specViewer );


  //HelpSystem::attachToolTipOn(m_options, "If checked, selection will be shown.",
  //                            showToolTips );

  m_showGammas->setChecked();
  m_showXrays->setChecked();
  
  m_showGammas->checked().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  m_showGammas->unChecked().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  
  m_showXrays->checked().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  m_showXrays->unChecked().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  
  m_showAlphas->checked().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  m_showAlphas->unChecked().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  
  m_showBetas->checked().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  m_showBetas->unChecked().connect( this, &ReferencePhotopeakDisplay::updateDisplayChange );
  
  m_showCascadeSums->checked().connect(this, &ReferencePhotopeakDisplay::updateDisplayChange);
  m_showCascadeSums->unChecked().connect(this, &ReferencePhotopeakDisplay::updateDisplayChange);
  
  m_showEscapes->checked().connect(this, &ReferencePhotopeakDisplay::updateDisplayChange);
  m_showEscapes->unChecked().connect(this, &ReferencePhotopeakDisplay::updateDisplayChange);
      
  m_otherNucsColumn = new WContainerWidget();

  m_otherNucsColumn->addStyleClass("OtherNucs ToolTabSection ToolTabTitledColumn");

  WText *otherNucTitle = new WText( WString::tr("rpd-suggestions"), m_otherNucsColumn);
  otherNucTitle->addStyleClass("ToolTabColumnTitle");

  m_otherNucs = new WContainerWidget(m_otherNucsColumn);
  m_otherNucs->addStyleClass( "OtherNucsContent ToolTabTitledColumnContent" );

  m_featureMarkerColumn = new WContainerWidget();
  m_featureMarkerColumn->addStyleClass("FeatureLines ToolTabSection ToolTabTitledColumn");
  m_featureMarkerColumn->hide();
  
  WContainerWidget *featureMarkerTitleRow = new WContainerWidget( m_featureMarkerColumn );
  featureMarkerTitleRow->addStyleClass( "ToolTabColumnTitle" );
  WText *featureMarkerTitle = new WText( WString::tr("rpd-feature-markers"), featureMarkerTitleRow );
  WContainerWidget *featureMarkerCloseIcon = new WContainerWidget(featureMarkerTitleRow);
  featureMarkerCloseIcon->addStyleClass("closeicon-wtdefault");
  //A little convoluted, but we will have the InterSpec class tell us to close feature marker widget,
  //  so it can save its state, and update the app menu-item, and handle undo/redo.
  featureMarkerCloseIcon->clicked().connect( boost::bind( &InterSpec::displayFeatureMarkerWindow, m_spectrumViewer, false) );
      
  m_particleView = new RowStretchTreeView();
  
  m_particleView->setRootIsDecorated(	false ); //makes the tree look like a table! :)
  
  m_particleView->addStyleClass( "ParticleViewTable ToolTabSection" );
  
  m_particleModel = new DecayParticleModel( this );
  m_particleView->setModel( m_particleModel );
  m_particleView->setAlternatingRowColors( true );
  m_particleView->setSortingEnabled( true );
  m_particleView->setColumnWidth( DecayParticleModel::kEnergy,
                                  WLength(6.5,WLength::FontEm) );
  m_particleView->setColumnWidth( DecayParticleModel::kBranchingRatio,
                                  WLength(6,WLength::FontEm) );
  m_particleView->setColumnWidth( DecayParticleModel::kResponsibleNuc,
                                 WLength(5,WLength::FontEm) );
  m_particleView->setColumnWidth( DecayParticleModel::kDecayMode,
                                 WLength(4,WLength::FontEm) );
  m_particleView->setColumnWidth( DecayParticleModel::kParticleType,
                                 WLength(5,WLength::FontEm) );
  
  auto bottomRow = new WContainerWidget();
  auto helpBtn = new WContainerWidget(bottomRow);
  helpBtn->addStyleClass("Wt-icon ContentHelpBtn RefGammaHelp");
  helpBtn->clicked().connect(boost::bind(&HelpSystem::createHelpWindow, "reference-gamma-lines-dialog"));

  
  RefGammaCsvResource *csv = new RefGammaCsvResource( this );
  csv->setTakesUpdateLock( true );
  
#if( BUILD_AS_OSX_APP || IOS )
  WAnchor *csvButton = new WAnchor( WLink(csv), bottomRow );
  csvButton->setTarget( AnchorTarget::TargetNewWindow );
  csvButton->setStyleClass( "LinkBtn DownloadLink RefGammaCsv" );
#else
  WPushButton *csvButton = new WPushButton( bottomRow );
  csvButton->setIcon( "InterSpec_resources/images/download_small.svg" );
  csvButton->setLink( WLink(csv) );
  csvButton->setLinkTarget( Wt::TargetNewWindow );
  csvButton->setStyleClass( "LinkBtn DownloadBtn RefGammaCsv" );
  
#if( ANDROID )
  // Using hacked saving to temporary file in Android, instead of via network download of file.
  csvButton->clicked().connect( std::bind([csv](){ android_download_workaround(csv, "photopeak_ref_info.csv"); }) );
#endif //ANDROID
  
#endif // BUILD_AS_OSX_APP / else
  
  csvButton->clicked().connect( std::bind([](){
    passMessage( WString::tr("rpd-csv-export-msg"), WarningWidget::WarningMsgInfo );
    //TODO: check about calling WarningWidget::displayPopupMessageUnsafe( msg, level, 5000 ); directly with a longer time for the message to hang around
  }));
  
  csvButton->setText( WString::tr("CSV") );
  csvButton->disable();
  m_csvDownload = csvButton;


  WGridLayout *overallLayout = new WGridLayout();
  overallLayout->setContentsMargins( 0, 0, 0, 0 );
  setLayout( overallLayout );

  overallLayout->addWidget( inputDiv,              0, 0 );
  overallLayout->addWidget( lowerInput,            1, 0 );
  overallLayout->addWidget( m_options,             0, 1, 3, 1 );
  overallLayout->addWidget( m_otherNucsColumn,     0, 2, 3, 1 );
  overallLayout->addWidget( m_featureMarkerColumn, 0, 3, 3, 1 );
  overallLayout->addWidget( m_particleView,        0, 4, 3, 1 );
  overallLayout->addWidget( bottomRow,             2, 0 );

  overallLayout->setRowStretch( 2, 1 );
  overallLayout->setColumnStretch( 4, 1 );
}//ReferencePhotopeakDisplay constructor


ReferencePhotopeakDisplay::~ReferencePhotopeakDisplay()
{
  //I think the DOM root should take care of deleting m_nuclideSuggest
  if( m_nuclideSuggest )
    delete m_nuclideSuggest;
}//~ReferencePhotopeakDisplay()


const ReferenceLineInfo &ReferencePhotopeakDisplay::currentlyShowingNuclide() const
{
  return m_currentlyShowingNuclide;
}


const std::vector<ReferenceLineInfo> &
                                 ReferencePhotopeakDisplay::persistedNuclides() const
{
  return m_persisted;
}

std::vector<ReferenceLineInfo> ReferencePhotopeakDisplay::showingNuclides() const
{
  std::vector<ReferenceLineInfo> answer;
  
  if( m_currentlyShowingNuclide.m_validity == ReferenceLineInfo::InputValidity::Valid )
    answer.push_back( m_currentlyShowingNuclide );
  answer.insert( answer.end(), m_persisted.begin(), m_persisted.end() );
  
  return answer;
}//std::vector<ReferenceLineInfo> showingNuclides() const;


const DetectorDisplay *ReferencePhotopeakDisplay::detectorDisplay() const
{
  return m_detectorDisplay;
};


const ShieldingSelect *ReferencePhotopeakDisplay::shieldingSelect() const
{
  return m_shieldingSelect;
};


const DecayParticleModel *ReferencePhotopeakDisplay::particleModel() const
{
  return m_particleModel;
}


void ReferencePhotopeakDisplay::setFocusToIsotopeEdit()
{
  m_nuclideEdit->setFocus();
  
#if( WT_VERSION >= 0x3030400 )
  InterSpecApp *app = dynamic_cast<InterSpecApp *>( wApp );
  const bool isMobile = (app && app->isMobile());
  const int nchar = static_cast<int>(m_nuclideEdit->text().narrow().length());
  if( nchar && !isMobile )
    m_nuclideEdit->setSelection( 0, nchar );
#endif
}//void setFocusToIsotopeEdit()


void ReferencePhotopeakDisplay::handleIsotopeChange( const bool useCurrentAge )
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const string isotopeLabel = m_nuclideEdit->text().toUTF8();
  const SandiaDecay::Nuclide *nuc = db->nuclide( isotopeLabel );
  const string agestr = m_ageEdit->text().toUTF8();
  const SandiaDecay::Element *el = NULL;

  if( nuc )
  {
    bool prompt = (nuc->canObtainPromptEquilibrium()
                   && m_promptLinesOnly->isChecked());
    if( (prompt || nuc->decaysToStableChildren()) && m_ageEdit->isEnabled() )
      m_ageEdit->disable();
    else if( !m_ageEdit->isEnabled() )
      m_ageEdit->enable();
  }//if( nuc )
  
  try
  {
    if( nuc && !nuc->isStable() && !nuc->decaysToChildren.empty() && !useCurrentAge )
    {
      if( nuc->decaysToStableChildren() )
      {
        m_ageEdit->setText( "0y" );
      }else if( nuc->canObtainPromptEquilibrium() && m_promptLinesOnly->isChecked() )
      {
        WString hlstr = PhysicalUnitsLocalized::printToBestTimeUnits(
                                          5.0*nuc->promptEquilibriumHalfLife(),
                                          2, SandiaDecay::second );
        m_ageEdit->setText( hlstr );
      }else if( m_currentlyShowingNuclide.m_input.m_input_txt != nuc->symbol )
      {
        string defagestr;
        PeakDef::defaultDecayTime( nuc, &defagestr );
        m_ageEdit->setText( defagestr );
      }else
      {
        const double hl = (nuc ? nuc->halfLife : -1.0);
        double age = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( agestr, hl );
        if( age > 100.0*nuc->halfLife || age < 0.0 )
          throw std::runtime_error( "" );
      }//if( nuc->decaysToStableChildren() ) / else
    }else if( nuc && useCurrentAge )
    {
      if( agestr.empty() && nuc->decaysToStableChildren() )
      {
        // We dont need an age - we will set to zero - so dont throw exception from empty string
      }else
      {
        double age = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( agestr, nuc->halfLife );
        if( age > 100.0*nuc->halfLife || age < 0.0 )
          throw std::runtime_error( "" );
      }
    }
  }catch(...)
  {
    if( nuc )
    {
      string defagestr;
      PeakDef::defaultDecayTime( nuc, &defagestr );
      passMessage( WString::tr("rpd-changed-age").arg(nuc->symbol).arg(agestr).arg(defagestr),
                   WarningWidget::WarningMsgLow );
      m_ageEdit->setText( defagestr );
    }else
    {
      m_ageEdit->setText( "0y" );
    }//if( nuc ) / else
  }//try / catch
  
  updateDisplayChange();
  
  if( nuc || el )
    m_displayingNuclide.emit();
  else if( m_persisted.empty() )
    m_nuclidesCleared.emit();
}//void handleIsotopeChange();


void ReferencePhotopeakDisplay::handleDrfChange( std::shared_ptr<DetectorPeakResponse> det )
{
  //ToDo: When the detector is changed, should maybe update all the persisted gamma lines
  //      as well - needs more thought
  
  // Block the undoing this, as the
  UpdateGuard guard( m_currently_updating );
  
  // The detector display may not have been updated yet (since it is connected to the
  //  `InterSpec::detectorChanged()` signal earlier, the call to `DetectorDisplay::setDetector(det)`
  //  actually happens later).
  m_detectorDisplay->setDetector( det );
  
  const RefLineInput user_input = userInput();
  updateDisplayFromInput( user_input );
}//void handleDrfChange()


std::map<std::string,std::vector<Wt::WColor>> ReferencePhotopeakDisplay::currentlyUsedPeakColors()
{
  std::map<string,vector<WColor>> answer;
  
  shared_ptr<const deque<PeakModel::PeakShrdPtr>> peaks;
  PeakModel *peakModel = m_spectrumViewer->peakModel();
  if( peakModel )
    peaks = peakModel->peaks();
  
  if( peaks )
  {
    for( const auto &p : *peaks )
    {
      string src;
      const WColor &color = p->lineColor();
      if( color.isDefault() )
        continue;
      if( p->parentNuclide() )
        src = p->parentNuclide()->symbol;
      else if( p->xrayElement() )
        src = p->xrayElement()->symbol;
      else if( p->reaction() )
        src = p->reaction()->name();
      else
        continue;
      vector<WColor> &colors = answer[src];
      if( std::find(begin(colors), end(colors), color) == end(colors) )
        colors.push_back( color );
    }//for( const auto &p : *peaks )
  }//if( peaks )

  for( const auto &p : m_persisted )
  {
    if( p.m_input.m_color.isDefault() )
      continue;
    
    vector<WColor> &colors = answer[p.m_input.m_input_txt];
    if( std::find(begin(colors), end(colors), p.m_input.m_color) == end(colors) )
      colors.push_back( p.m_input.m_color );
  }//for( const auto &p : m_persisted )
  
  return answer;
}//currentlyUsedPeakColors()


void ReferencePhotopeakDisplay::toggleShowOptions()
{
  if( m_options->isHidden() )
  {
    m_options->show();
    //m_options->animateShow(WAnimation(WAnimation::AnimationEffect::Pop, WAnimation::TimingFunction::Linear, 250) );
    m_options_icon->addStyleClass("active");
    
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
    if( m_particleView->isHidden() ) //Narrow phone display
    {
      m_otherNucsColumn->setHidden( true );
      m_featureMarkerColumn->setHidden( true );
    }
#endif
  }else
  {
    m_options->hide();
    //m_options->animateHide(WAnimation(WAnimation::AnimationEffect::Pop, WAnimation::TimingFunction::Linear, 250));
    m_options_icon->removeStyleClass("active");
    
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
    if( m_otherNucsColumn->isHidden() ) //Narrow phone display
      m_otherNucsColumn->setHidden( m_featureMarkerColumn->isVisible() );
#endif
  }
}//void toggleShowOptions()


void ReferencePhotopeakDisplay::handleSpectrumChange(SpecUtils::SpectrumType type)
{
  // For simplicity we'll update the suggested nuclides whenever foreground gets
  //  changed, even if its only because of sample number change (which will have
  //  same RIID results, as they are file-specific, not sample specific).
  //  We will also clear external RID results (e.g., if user has setup to
  //  automatically use a web-service on spectrum load).
  if( type == SpecUtils::SpectrumType::Foreground )
  {
    m_external_ids.clear();
    updateOtherNucsDisplay();
  }
}//handleSpectrumChange();


void ReferencePhotopeakDisplay::updateAssociatedNuclides()
{
  if( m_currentlyShowingNuclide.m_validity != ReferenceLineInfo::InputValidity::Valid )
    return;
  
  const string &currentInput = m_currentlyShowingNuclide.m_input.m_input_txt;
  
  const MoreNuclideInfo::InfoStatus status = MoreNuclideInfo::more_nuc_info_db_status();
  if( status == MoreNuclideInfo::InfoStatus::FailedToInit )
    return;

  if( status != MoreNuclideInfo::InfoStatus::Inited )
    cerr << "\n\n\nWarning, updateAssociatedNuclides called before MoreNuclideInfo initualized.\n\n";

  const auto infoDb = MoreNuclideInfo::MoreNucInfoDb::instance();
  if( !infoDb )
  {
    assert( MoreNuclideInfo::more_nuc_info_db_status() == MoreNuclideInfo::InfoStatus::FailedToInit );
    return;
  }//if( !infoDb )

  const MoreNuclideInfo::NucInfo *info = nullptr;
  if( m_currentlyShowingNuclide.m_nuclide )
    info = infoDb->info( m_currentlyShowingNuclide.m_nuclide );
  else
    info = infoDb->info( currentInput );

  if( !info || info->m_associated.empty() )
    return;

  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();

  WText *header = new WText( WString::tr("rpd-assoc-nucs") );
  header->addStyleClass( "OtherNucTypeHeader" );

  m_otherNucs->insertWidget( 0, header );

  for( size_t index = 0; index < info->m_associated.size(); ++index )
  {
    const string &nucstr = info->m_associated[index];

    WPushButton *btn = new WPushButton( nucstr );
    btn->addStyleClass( "LinkBtn" );

    m_otherNucs->insertWidget( static_cast<int>(1 + index), btn );

    // If text is a valid Nuclide, Element, or Reaction, we'll make this
    //  button clickable to display; otherwise we'll show the text, but
    //  disable the button.
    const SandiaDecay::Nuclide *const nuc = db->nuclide( nucstr );
    const SandiaDecay::Element *const el = nuc ? nullptr : db->element( nucstr );

    std::vector<ReactionGamma::ReactionPhotopeak> reactions;
    if( !nuc && !el )
    {
      try
      {
        const ReactionGamma *rctnDb = ReactionGammaServer::database();
        if( rctnDb )
          rctnDb->gammas( nucstr, reactions );
      }catch( std::exception &e)
      {
      }
    }//if( !nuc && !el )

    if( nuc || el || !reactions.empty() )
    {
      RefLineInput input = userInput();
      input.m_input_txt = nucstr;
      
      btn->clicked().connect( boost::bind( &ReferencePhotopeakDisplay::updateDisplayFromInput, this, input ) );
    }else
    {
      btn->disable();
    }
  }//for( const auto &riid : riid_nucs )
}//void updateAssociatedNuclides()


void ReferencePhotopeakDisplay::programmaticallyCloseMoreInfoWindow()
{
  if( m_nucInfoWindow )
  {
    UndoRedoManager::BlockUndoRedoInserts undo_blocker;
    m_nucInfoWindow->done(Wt::WDialog::DialogCode::Accepted);
  }
  assert( !m_nucInfoWindow );
  m_nucInfoWindow = nullptr;
}//void programmaticallyCloseMoreInfoWindow()


void ReferencePhotopeakDisplay::handleMoreInfoWindowClose( MoreNuclideInfoWindow *window )
{
  if( window == m_nucInfoWindow )
  {
    m_nucInfoWindow = nullptr;
    
    UndoRedoManager *undo_manager = UndoRedoManager::instance();
    if( undo_manager && undo_manager->canAddUndoRedoNow() )
    {
      // I *think* calling `window->currentNuclide()` would be valid, but lets not risk it
      auto undo = [](){
        InterSpec *interspec = InterSpec::instance();
        ReferencePhotopeakDisplay *disp = interspec ? interspec->referenceLinesWidget() : nullptr;
        if( disp )
          disp->showMoreInfoWindow();
      };//
      
      auto redo = [](){
        InterSpec *interspec = InterSpec::instance();
        ReferencePhotopeakDisplay *disp = interspec ? interspec->referenceLinesWidget() : nullptr;
        if( disp )
          disp->programmaticallyCloseMoreInfoWindow();
      };
      
      undo_manager->addUndoRedoStep( undo, redo, "Close nuclide more info window." );
    }//if( undo_manager && undo_manager->canAddUndoRedoNow() )
  }else
  {
    cerr << "ReferencePhotopeakDisplay::handleMoreInfoWindowClose: Received pointer (" << window
         << "), not matching m_nuclidesCleared (" << m_nucInfoWindow << ")" << endl;
  }
}//void handleMoreInfoWindowClose( MoreNuclideInfoWindow *window );


MoreNuclideInfoWindow *ReferencePhotopeakDisplay::moreInfoWindow()
{
  return m_nucInfoWindow;
}

#if( InterSpec_PHONE_ROTATE_FOR_TABS )
void ReferencePhotopeakDisplay::setNarrowPhoneLayout( const bool narrow )
{
  if( m_particleView->isHidden() == narrow )
    return;
  
  m_particleView->setHidden( narrow );
  
  const char *add_key = narrow ? "rpd-add-another-btn-narrow" : "rpd-add-another-btn";
  m_persistLines->setText( WString::tr(add_key) );
  
  WGridLayout *lay = dynamic_cast<WGridLayout *>( layout() );
  assert( lay );
  if( lay )
  {
    if( narrow )
    {
      lay->setColumnStretch( 4, 0 );
      lay->setColumnStretch( 0, 1 );
    }else
    {
      lay->setColumnStretch( 4, 1 );
      lay->setColumnStretch( 0, 0 );
    }
  }//if( lay )
}//void setNarrowPhoneLayout( const bool narrow )
#endif //InterSpec_PHONE_ROTATE_FOR_TABS


FeatureMarkerWidget *ReferencePhotopeakDisplay::featureMarkerTool()
{
  return m_featureMarkers;
}


FeatureMarkerWidget *ReferencePhotopeakDisplay::showFeatureMarkerTool()
{
  if( m_featureMarkers )
    return m_featureMarkers;
  
  m_showFeatureMarkers->setChecked( true );
  m_featureMarkerColumn->setHidden( false );
  m_featureMarkers = new FeatureMarkerWidget( m_spectrumViewer, m_featureMarkerColumn );
  
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
  if( m_particleView->isHidden() ) //Narrow phone display
  {
    if( !m_options->isHidden() )
      toggleShowOptions();
    m_otherNucsColumn->setHidden( true );
  }//
#endif
  
  return m_featureMarkers;
}//FeatureMarkerWidget *showFeatureMarkerTool()


void ReferencePhotopeakDisplay::removeFeatureMarkerTool()
{
  if( !m_featureMarkers )
    return;
  
  delete m_featureMarkers;
  m_featureMarkers = nullptr;
  m_featureMarkerColumn->hide();
  m_showFeatureMarkers->setChecked( false );
  
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
  if( m_particleView->isHidden() ) //Narrow phone display
    m_otherNucsColumn->setHidden( !m_options->isHidden() );
#endif
}//void removeFeatureMarkerTool()


void ReferencePhotopeakDisplay::featureMarkerCbToggled()
{
  // This is a little convoluted, but we will call back to the InterSpec
  //  class, which will then call to the appropriate `ReferencePhotopeakDisplay`
  //  functions to create/remove the display.
  //  This is to allow the InterSpec class to restore the widget to the correct
  //  state, and also update its menu items, and handle undo/redo.
  m_spectrumViewer->displayFeatureMarkerWindow( m_showFeatureMarkers->isChecked() );
}//void featureMarkerCbToggled()


void ReferencePhotopeakDisplay::emphasizeFeatureMarker()
{
  if( m_featureMarkers )
    m_featureMarkers->doJavaScript( "$('#" + m_featureMarkers->id() + "')"
                              ".fadeIn(100).fadeOut(100).fadeIn(100).fadeOut(100).fadeIn(100);" );
}//void emphasizeFeatureMarker()


void ReferencePhotopeakDisplay::showMoreInfoWindow()
{
  const SandiaDecay::Nuclide * const nuc = m_currentlyShowingNuclide.m_nuclide;
  assert( nuc );
  if( !nuc )
    return;

  const SandiaDecay::Nuclide *prev_orig_nuc = nullptr, *prev_current_nuc = nullptr;
  if( m_nucInfoWindow )
  {
    UndoRedoManager::BlockUndoRedoInserts undo_blocker;
    
    prev_orig_nuc = m_nucInfoWindow->originalNuclide();
    prev_orig_nuc = m_nucInfoWindow->currentNuclide();
    m_nucInfoWindow->done(Wt::WDialog::DialogCode::Accepted);
    assert( m_nucInfoWindow == nullptr );
    m_nucInfoWindow = nullptr;
  }//if( m_nucInfoWindow )
  
  m_nucInfoWindow = new MoreNuclideInfoWindow( nuc );
  m_nucInfoWindow->finished().connect( 
                               boost::bind( &ReferencePhotopeakDisplay::handleMoreInfoWindowClose,
                                 this, m_nucInfoWindow )
  );
  
  // All of this undo/redo stuff is a little over the top since we will only ever show one more-info
  //  window at a time, but oh well.
  UndoRedoManager *undo_manager = UndoRedoManager::instance();
  if( undo_manager && undo_manager->canAddUndoRedoNow() )
  {
    auto undo = [this, prev_orig_nuc, prev_current_nuc](){
      InterSpec *interspec = InterSpec::instance();
      ReferencePhotopeakDisplay *disp = interspec ? interspec->referenceLinesWidget() : nullptr;
      assert( disp );
      MoreNuclideInfoWindow *window = disp ? disp->moreInfoWindow() : nullptr;
      if( !window )
        return;
      
      window->done(Wt::WDialog::DialogCode::Accepted);
      
      if( prev_orig_nuc )
      {
        assert( !m_nucInfoWindow );
        m_nucInfoWindow = new MoreNuclideInfoWindow( prev_orig_nuc );
        m_nucInfoWindow->finished().connect( 
                  boost::bind( &ReferencePhotopeakDisplay::handleMoreInfoWindowClose,
                              this, m_nucInfoWindow )
        );
        
        if( prev_current_nuc && (prev_orig_nuc != prev_current_nuc) )
        {
          // TODO: use prev_current_nuc to kinda track history or whatever
        }
      }//if( prev_orig_nuc )
    };//undo
    
    auto redo = [this, nuc](){
      InterSpec *interspec = InterSpec::instance();
      ReferencePhotopeakDisplay *disp = interspec ? interspec->referenceLinesWidget() : nullptr;
      assert( disp );
      if( !disp )
        return;
      
      if( disp->m_nucInfoWindow )
        disp->m_nucInfoWindow->done(Wt::WDialog::DialogCode::Accepted);
      
      assert( !disp->moreInfoWindow() );
      disp->m_nucInfoWindow = new MoreNuclideInfoWindow( nuc );
      disp->m_nucInfoWindow->finished().connect(
                              boost::bind( &ReferencePhotopeakDisplay::handleMoreInfoWindowClose,
                              this, m_nucInfoWindow )
      );
    };//redo
    
    undo_manager->addUndoRedoStep( undo, redo, "Show " + nuc->symbol + " more info window." );
  }//if( undo_manager && undo_manager->canAddUndoRedoNow() )
}//void showMoreInfoWindow()


void ReferencePhotopeakDisplay::updateOtherNucsDisplay()
{
  m_otherNucs->clear();

  const bool showRiid = m_showRiidNucs->isChecked();
  const bool showPrev = m_showPrevNucs->isChecked();
  const bool showAssoc = m_showAssocNucs->isChecked();

  if( !showRiid && !showPrev && !showAssoc )
  {
    m_otherNucsColumn->hide();
    return;
  }

#if( InterSpec_PHONE_ROTATE_FOR_TABS )
  if( m_particleView->isHidden() ) //Narrow phone display
    m_otherNucsColumn->setHidden( m_options->isVisible() || m_featureMarkerColumn->isVisible() );
#else
  m_otherNucsColumn->show();
#endif
  
  vector<RefLineInput> prev_nucs;
  const string &currentInput = m_currentlyShowingNuclide.m_input.m_input_txt;

  if( showPrev )
  {
    for( const auto &prev : m_prevNucs )
    {
      if( prev.m_input_txt != currentInput )
        prev_nucs.push_back(prev);
    }
  }//if( showPrev )


  vector<pair<string, string>> riid_nucs;  //<description, nuclide>
  shared_ptr<const SpecUtils::DetectorAnalysis> riid_ana;
  shared_ptr<const SpecMeas> m = m_spectrumViewer->measurment(SpecUtils::SpectrumType::Foreground);
  if( m )
    riid_ana = m->detectors_analysis();

  if( showRiid && riid_ana )
  {
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();

    for( const SpecUtils::DetectorAnalysisResult &res : riid_ana->results_ )
    {
      if( res.nuclide_.empty() || res.isEmpty() )
        continue;

      bool already_have = false;
      for( const auto &v : riid_nucs )
        already_have |= (v.first == res.nuclide_);
      if( already_have )
        continue;

      string nuc_name = res.nuclide_;

      const SandiaDecay::Nuclide *nuc = db->nuclide(nuc_name);
      if( !nuc )
      {
        vector<string> fields;
        SpecUtils::split(fields, nuc_name, " \t,");
        for( const auto &v : fields )
        {
          nuc = db->nuclide(v);
          if( nuc )
          {
            nuc_name = v;
            break;
          }
        }//for( const auto &v : fields )
      }//if( !nuc )

      if( nuc )
        riid_nucs.push_back( {nuc->symbol, res.nuclide_} );
      else 
        riid_nucs.push_back( {res.nuclide_, ""} );
    }//for( loop over RIID results )
  }//if( riid_ana )


  // TODO: Need to implement option to show or not show suggestion catagories, and if show none, then hide the column
  if( (m_currentlyShowingNuclide.m_validity == ReferenceLineInfo::InputValidity::Valid)
     && !currentInput.empty()
     && showAssoc )
  {
    const MoreNuclideInfo::InfoStatus status = MoreNuclideInfo::more_nuc_info_db_status();

    switch( status )
    {
      case MoreNuclideInfo::InfoStatus::Inited:
        updateAssociatedNuclides();
        break;

      case MoreNuclideInfo::InfoStatus::FailedToInit:
        break;

      case MoreNuclideInfo::InfoStatus::NotInited:
      {
        const string sessionid = wApp->sessionId();
        boost::function<void()> update_gui = wApp->bind( boost::bind(&ReferencePhotopeakDisplay::updateAssociatedNuclides, this) );
       
        auto worker = [update_gui, sessionid](){
          const auto infoDb = MoreNuclideInfo::MoreNucInfoDb::instance();
          if( !infoDb )
          {
            assert( MoreNuclideInfo::more_nuc_info_db_status() == MoreNuclideInfo::InfoStatus::FailedToInit );
            return;
          }

          auto inner_worker = [update_gui]{
            update_gui();
            wApp->triggerUpdate();
          };

          Wt::WServer::instance()->post( sessionid, inner_worker );
        };//worker

        Wt::WServer::instance()->ioService().boost::asio::io_service::post( worker );
        break;
      }//case MoreNuclideInfo::InfoStatus::NotInited:
    }//switch( status )
  }//if( !currentInput.empty() )


  // TODO: Need to deal the "more info" button

  auto displayDetectorOrExternal = [=]( const WString &title, vector<pair<string,string>> nucs ){
    if( nucs.empty() )
      return;
    
    // Protect against a detector having a ton of results
    const size_t max_riid_res = 10;
    if( nucs.size() > max_riid_res )
      nucs.resize( max_riid_res );

    WText *header = new WText( title, m_otherNucs);
    header->addStyleClass("OtherNucTypeHeader");
    
    for( const auto &riid : nucs )
    {
      WPushButton *btn = new WPushButton(riid.first, m_otherNucs);
      btn->addStyleClass( "LinkBtn" );
      
      if( !riid.second.empty() )
      {
        RefLineInput input = userInput();
        input.m_input_txt = riid.first;
        
        btn->clicked().connect(boost::bind(&ReferencePhotopeakDisplay::updateDisplayFromInput, this, input));
      }else
      {
        btn->disable();
      }
    }//for( const auto &riid : nucs )
  };//displayDetectorOrExternal lamda

  if( !riid_nucs.empty() )
    displayDetectorOrExternal( WString::tr("rpd-det-id"), riid_nucs );
  
  if( !m_external_ids.empty() )
  {
    WString name = m_external_algo_name.empty() ? WString::tr("rpd-ext-rid") : WString::fromUTF8(m_external_algo_name);
    displayDetectorOrExternal( name, m_external_ids );
  }//if( !m_external_ids.empty() )
  
  
  if( !prev_nucs.empty() )
  {
    WText *header = new WText( WString::tr("rpd-prev"), m_otherNucs);
    header->addStyleClass("OtherNucTypeHeader");
    for( const auto &prev : prev_nucs )
    {
      WPushButton *btn = new WPushButton(prev.m_input_txt, m_otherNucs);
      btn->addStyleClass( "LinkBtn" );

      btn->clicked().connect(boost::bind(&ReferencePhotopeakDisplay::updateDisplayFromInput, this, prev));
    }
  }//if( !prev_nucs.empty() )
}//void updateOtherNucsDisplay()


RefLineInput ReferencePhotopeakDisplay::userInput() const
{
  RefLineInput input;

  input.m_input_txt = m_nuclideEdit->text().toUTF8();
  input.m_age = m_ageEdit->text().toUTF8();
  input.m_color = m_colorSelect->color();

  if( m_lowerBrCuttoff && (m_lowerBrCuttoff->validate() == WValidator::Valid) )
    input.m_lower_br_cutt_off = m_lowerBrCuttoff->value();

  input.m_promptLinesOnly = (m_promptLinesOnly && m_promptLinesOnly->isChecked());
  input.m_showGammas = (!m_showGammas || m_showGammas->isChecked());
  input.m_showXrays = (!m_showXrays || m_showXrays->isChecked());
  input.m_showAlphas = (m_showAlphas && m_showAlphas->isChecked());
  input.m_showBetas = (m_showBetas && m_showBetas->isChecked());
  input.m_showCascades = (m_showCascadeSums && m_showCascadeSums->isChecked());
  input.m_showEscapes = (m_showEscapes && m_showEscapes->isChecked());
  
  if( m_detectorDisplay->detector() )
  {
    input.m_det_intrinsic_eff = m_detectorDisplay->detector()->intrinsicEfficiencyFcn();
    if( input.m_det_intrinsic_eff )
      input.m_detector_name = m_detectorDisplay->detector()->name();
  }//if( m_detectorDisplay->detector() )


  try
  {
    if( m_shieldingSelect->isGenericMaterial() )
    {
      const float atomic_number = static_cast<float>(m_shieldingSelect->atomicNumber());
      const float areal_density = static_cast<float>(m_shieldingSelect->arealDensity());
      
      if( areal_density > 0.0 )
      {
        input.m_shielding_att = [atomic_number, areal_density]( float energy ) -> double {
          const double att_coef = GammaInteractionCalc::transmition_coefficient_generic( atomic_number, areal_density, energy );
          return exp( -1.0 * att_coef );
        };
        
        NativeFloatSpinBox *anEdit = m_shieldingSelect->atomicNumberEdit();
        NativeFloatSpinBox *adEdit = m_shieldingSelect->arealDensityEdit();
        assert( adEdit && adEdit );
        
        if( adEdit && adEdit )
        {
          const string an = SpecUtils::trim_copy( anEdit->text().toUTF8() );
          const string ad = SpecUtils::trim_copy( adEdit->text().toUTF8() );
          if( !an.empty() && !ad.empty() )
          {
            input.m_shielding_an = an;
            input.m_shielding_ad = ad;
          }
        }//if( adEdit && adEdit )
      }//if( areal_density > 0.0 )
    }else
    {
      std::shared_ptr<const Material> material = m_shieldingSelect->material();
      if( material )
        input.m_shielding_name = material->name;
      
      float thick = 0.0f;
      input.m_shielding_thickness = m_shieldingSelect->thicknessEdit()->text().toUTF8();
      if( !input.m_shielding_thickness.empty() )
        thick = static_cast<float>( m_shieldingSelect->thickness() );
      
      if( material && (thick > 0.0) )
      {
        input.m_shielding_att = [material, thick]( float energy ) -> double {
          const double att_coef = GammaInteractionCalc::transmition_coefficient_material( material.get(), energy, thick );
          return exp( -1.0 * att_coef );
        };
      }//if( material && (thick > 0.0) )
    }//if( isGenericMaterial ) / else
  }catch( std::exception &e )
  {
    cerr << "Exception getting shielding: " << e.what() << endl;
  }//try / catch to get shielding

  SpecUtils::trim( input.m_age );
  SpecUtils::trim( input.m_input_txt );

  return input;
}//RefLineInput userInput();


std::vector<DecayParticleModel::RowData> ReferencePhotopeakDisplay::createTableRows( const ReferenceLineInfo &refLine )
{
  vector<DecayParticleModel::RowData> inforows;

  for( const ReferenceLineInfo::RefLine &r : refLine.m_ref_lines )
  {
    DecayParticleModel::RowData row;

    row.energy = r.m_energy;
    row.branchRatio = r.m_decay_intensity;
    row.responsibleNuc = nullptr;

    row.decayMode = SandiaDecay::DecayMode::UndefinedDecay; //JIC
    row.particle = SandiaDecay::ProductType::GammaParticle; //JIC

    if( r.m_transition )
    {
      row.responsibleNuc = r.m_transition->parent;
      row.decayMode = r.m_transition->mode;
    }

    switch( r.m_particle_type )
    {
      case ReferenceLineInfo::RefLine::Particle::Alpha:
        row.particle = SandiaDecay::ProductType::AlphaParticle;
        break;

      case ReferenceLineInfo::RefLine::Particle::Beta:
        row.particle = SandiaDecay::ProductType::BetaParticle;
        break;
      
      case ReferenceLineInfo::RefLine::Particle::Gamma:
        row.particle = SandiaDecay::ProductType::GammaParticle;
        break;

      case ReferenceLineInfo::RefLine::Particle::Xray:
        row.responsibleNuc = r.m_parent_nuclide; // TODO: this is only for comparison - not correct!
        row.particle = SandiaDecay::ProductType::XrayParticle;
        row.decayMode = DecayParticleModel::RowData::XRayDecayMode;
        break;
    }//switch( r.m_particle_type )


    switch( r.m_source_type )
    {
      case ReferenceLineInfo::RefLine::RefGammaType::Normal:
      case ReferenceLineInfo::RefLine::RefGammaType::Annihilation:
        if( r.m_reaction )
          row.decayMode = DecayParticleModel::RowData::ReactionToGammaMode;
        break;

      case ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak:
      case ReferenceLineInfo::RefLine::RefGammaType::SumGammaPeak:
        row.decayMode = DecayParticleModel::RowData::CascadeSumMode;
        break;

      case ReferenceLineInfo::RefLine::RefGammaType::SingleEscape:
      case ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape:
        // We dont want these making it into the table
        continue;
    }//switch( r.m_source_type )

    inforows.push_back( row );
  }//for( const ReferenceLineInfo::RefLine &r : refLine.m_ref_lines )


  return inforows;
}//vector<DecayParticleModel::RowData> createTableRows( const ReferenceLineInfo &refLine );


Wt::WColor ReferencePhotopeakDisplay::colorForNewSource( const std::string &src )
{
  WColor color = m_colorSelect->color();
  
  if( m_peaksGetAssignedRefLineColor )
  {
#ifndef _MSC_VER
#warning "Need to test string comparison for sources always work.  E.g., need case-insensitive, etc."
#endif
    const map<string,vector<WColor>> usedColors = currentlyUsedPeakColors();
    
    auto hasBeenUsed = [&usedColors](const WColor &color)->bool{
      for( const auto &s : usedColors )
        for( const auto &c : s.second )
          if( c == color )
            return true;
      return false;
    };
    
    const auto usedIter = usedColors.find(src);
    const auto previter = m_previouslyPickedSourceColors.find(src);
    const auto specificiter = m_specificSourcelineColors.find(src);
    
    if( m_userHasPickedColor && !hasBeenUsed(color)
       && (specificiter == end(m_specificSourcelineColors)) )
    {
      //If the user picked a color, but didnt fit any peaks, or persist, the
      //  lines, and the color theme doesnt call out this current source
      //  explicitly, then use the previous color.
      //  This maybe seems a little more intuitive from the users perspective.
      
      //We also need to propagate m_userHasPickedColor==true forward for when
      //  user keeps entering different sources.
      m_userHasPickedColor = true;
    }else if( usedIter != end(usedColors) && !usedIter->second.empty() )
    {
      color = usedIter->second[0];
      //cout << "Source " << isotxt << " has color " << color.cssText() << " in the spectrum already" << endl;
    }else if( previter != end(m_previouslyPickedSourceColors)
             && !hasBeenUsed(previter->second) )
    {
      color = previter->second;
      //cout << "Source " << isotxt << " has previously had color " << color.cssText() << " picked" << endl;
    }else if( specificiter != end(m_specificSourcelineColors) )
    {
      color = specificiter->second;
      //cout << "Source " << isotxt << " has a specific color, " << color.cssText() << ", in the ColorTheme" << endl;
    }else
    {
      auto colorcopy = m_lineColors;
      for( const auto &p : usedColors )
      {
        for( const auto &c : p.second )
        {
          auto pos = std::find( begin(colorcopy), end(colorcopy), c );
          if( pos != end(colorcopy) )
            colorcopy.erase(pos);
        }
      }//for( loop over colors used for peaks already )
      
      for( const auto &p : m_persisted )
      {
        auto pos = std::find( begin(colorcopy), end(colorcopy), p.m_input.m_color );
        if( pos != end(colorcopy) )
          colorcopy.erase(pos);
      }//for( const auto &p : m_persisted )
      
      if( colorcopy.empty() )
        color = m_lineColors[m_persisted.size() % m_lineColors.size()];
      else
        color = colorcopy[0];
      //cout << "Source " << isotxt << " will select color, " << color.cssText()
      //     << ", from default list (len=" << colorcopy.size() << " of "
      //     << m_lineColors.size() << ")" << endl;
    }//if( src has been seen ) / else (user picked color previously)
  }else
  {
    color = m_lineColors[m_persisted.size() % m_lineColors.size()];
  }
  
  return color;
}//Wt::WColor colorForNewSource( const std::string &src );


void ReferencePhotopeakDisplay::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  // If we are here, we are done doing any updates, so lets make sure
  //  #m_currently_updating is set to false (even though it already should be)
  assert( !m_currently_updating );
  m_currently_updating = false;
  
  WContainerWidget::render( flags );
}


std::shared_ptr<void> ReferencePhotopeakDisplay::getDisableUndoRedoSentry()
{
  shared_ptr<void> answer = m_undo_redo_sentry.lock();
  if( answer )
    return answer;
  
  int *dummy = new int(0);
  auto deleter = []( void *obj ){
    int *sentry = (int *)obj;
    if( sentry )
      delete sentry;
  };
  
  answer = shared_ptr<void>( dummy, deleter );
  m_undo_redo_sentry = answer;
  return answer;
}//std::shared_ptr<void> getDisableUndoRedoSentry()


void ReferencePhotopeakDisplay::addUndoRedoPoint( const ReferenceLineInfo &starting_showing,
                      const std::vector<ReferenceLineInfo> &starting_persisted,
                      const std::deque<RefLineInput> &starting_prev_nucs,
                      const bool starting_user_color,
                      const RefLineInput &current_user_input )
{
  const bool ending_user_color = m_userHasPickedColor;
  
  UndoRedoManager *undo_manager = UndoRedoManager::instance();
  if( undo_manager
     && !(starting_showing == m_currentlyShowingNuclide)
     && !m_undo_redo_sentry.lock() )
  {
    auto undo = [this, starting_showing, starting_persisted, starting_prev_nucs, starting_user_color](){
      InterSpec *interspec = InterSpec::instance();
      ReferencePhotopeakDisplay *display = interspec ? interspec->referenceLinesWidget() : nullptr;
      if( !display )
        return;
      
      display->m_currentlyShowingNuclide = starting_showing;
      display->m_persisted = starting_persisted;
      display->m_prevNucs = starting_prev_nucs;
      
      display->updateDisplayFromInput( starting_showing.m_input );
      display->m_userHasPickedColor = starting_user_color;
    };//undo
    
    auto redo = [this, starting_showing, starting_persisted, starting_prev_nucs, current_user_input, ending_user_color](){
      InterSpec *interspec = InterSpec::instance();
      ReferencePhotopeakDisplay *display = interspec ? interspec->referenceLinesWidget() : nullptr;
      if( !display )
        return;
      
      display->m_currentlyShowingNuclide = starting_showing;
      display->m_persisted = starting_persisted;
      display->m_prevNucs = starting_prev_nucs;
      
      display->updateDisplayFromInput( current_user_input );
      display->m_userHasPickedColor = ending_user_color;
    };//undo
    
    
    undo_manager->addUndoRedoStep( undo, redo, "Update ref-lines." );
  }//if( undo_manager )
}//void addUndoRedoPoint(...)




void ReferencePhotopeakDisplay::updateDisplayChange()
{
  if( m_currently_updating )
    return;
  
  const ReferenceLineInfo starting_showing = m_currentlyShowingNuclide;
  const vector<ReferenceLineInfo> starting_persisted = m_persisted;
  const deque<RefLineInput> starting_prev_nucs = m_prevNucs;
  const bool starting_user_color = m_userHasPickedColor;
  
  const RefLineInput user_input = userInput();
  updateDisplayFromInput( user_input );
  
  addUndoRedoPoint( starting_showing, starting_persisted, starting_prev_nucs,
                   starting_user_color, user_input );
}//void updateDisplayChange()


void ReferencePhotopeakDisplay::updateDisplayFromInput( RefLineInput user_input )
{
  UpdateGuard guard( m_currently_updating );
  
  const ReferenceLineInfo starting_showing = m_currentlyShowingNuclide;
  const vector<ReferenceLineInfo> starting_persisted = m_persisted;
  const deque<RefLineInput> starting_prev_nucs = m_prevNucs;
  const bool starting_user_color = m_userHasPickedColor;
  
  shared_ptr<ReferenceLineInfo> ref_lines = ReferenceLineInfo::generateRefLineInfo( user_input );
  
  if( ref_lines )
  {
    assert( (ref_lines->m_source_type == ReferenceLineInfo::SourceType::Nuclide) == (ref_lines->m_nuclide != nullptr) );
    assert( (ref_lines->m_source_type == ReferenceLineInfo::SourceType::FluorescenceXray) == (ref_lines->m_element != nullptr) );
    assert( (ref_lines->m_source_type == ReferenceLineInfo::SourceType::Reaction) == (!ref_lines->m_reactions.empty()) );
    //ReferenceLineInfo::SourceType::Background
    //ReferenceLineInfo::SourceType::CustomEnergy
    //ReferenceLineInfo::SourceType::NuclideMixture:
    //ReferenceLineInfo::SourceType::OneOffSrcLines:
    //ReferenceLineInfo::SourceType::None
  }//if( ref_lines - do some quick sanity checks )
  
  
  const SandiaDecay::Nuclide * const nuclide = ref_lines ? ref_lines->m_nuclide : nullptr;
  const SandiaDecay::Element * const element = ref_lines ? ref_lines->m_element : nullptr;
  const std::string srcstr = ref_lines ? ref_lines->m_input.m_input_txt : string();
  const ReferenceLineInfo::InputValidity validity = ref_lines ? ref_lines->m_validity
  : ReferenceLineInfo::InputValidity::InvalidSource;
  const ReferenceLineInfo::SourceType src_type = ref_lines ? ref_lines->m_source_type
  : ReferenceLineInfo::SourceType::None;
  
  // Pass any warnings received from parsing the input
  if( ref_lines && !ref_lines->m_input_warnings.empty() )
  {
    for( const string &warning : ref_lines->m_input_warnings )
      passMessage( warning, WarningWidget::WarningMsgHigh );
  }
  
  const bool show_lines = ( (validity == ReferenceLineInfo::InputValidity::Valid)
                           && (src_type != ReferenceLineInfo::SourceType::None) );
  
  m_moreInfoBtn->setHidden( !nuclide );
  
  if( (validity == ReferenceLineInfo::InputValidity::Valid)
     && (m_nuclideEdit->text().toUTF8() != ref_lines->m_input.m_input_txt ) )
  {
    m_nuclideEdit->setText( WString::fromUTF8(ref_lines->m_input.m_input_txt) );
  }else if( validity == ReferenceLineInfo::InputValidity::Blank )
  {
    m_nuclideEdit->setText( "" );
  }
  
  const bool hasPromptEquilib = (nuclide && nuclide->canObtainPromptEquilibrium());
  m_promptLinesOnly->setHidden( !hasPromptEquilib );
  m_promptLinesOnly->setChecked( ref_lines && ref_lines->m_input.m_promptLinesOnly );
  if( !hasPromptEquilib )
    m_promptLinesOnly->setUnChecked();
  
  // For SourceType::NuclideMixture, we could look through all the nuclides and check if any of them
  //  are candidates for aging - but for now we'll just assume we can age them.
  const bool enable_aging = ((nuclide && !ref_lines->m_input.m_promptLinesOnly
                              && !nuclide->decaysToStableChildren())
                             || (ref_lines->m_source_type == ReferenceLineInfo::SourceType::NuclideMixture) ) ;
  
  string agestr = (!enable_aging || !ref_lines) ? string() : ref_lines->m_input.m_age;
  
  if( !agestr.empty() )
  {
    try
    {
      const Wt::WLocale &locale = Wt::WLocale::currentLocale();
      if( !locale.name().empty() && !SpecUtils::istarts_with(locale.name(), "en" ) )
      {
        const double hl = nuclide ? nuclide->halfLife : -1.0;
        const double duration = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( agestr, hl );
        agestr = PhysicalUnitsLocalized::printToBestTimeUnits( duration );
      }
    }catch( std::exception &e )
    {
      cerr << "Error converting age ('" << agestr << "') to localized time-span: " << e.what() << endl;
    }//
  }//if( !agestr.empty() )
  
  m_ageEdit->setText( WString::fromUTF8(agestr) );
  m_ageEdit->setEnabled( enable_aging );
  
  
  
  const string hl_str = !nuclide ? string() : PhysicalUnitsLocalized::printToBestTimeUnits( nuclide->halfLife, 2 );
  
  const WString hlstr = !nuclide 
                        ? WString()
                        : WString("{1}={2}").arg( WString::tr("T1/2") ).arg( hl_str );
  m_halflife->setText( hlstr );
  
  m_persistLines->setEnabled( show_lines );
  m_clearLines->setDisabled( m_persisted.empty() && !show_lines );
  
  if( m_csvDownload )
    m_csvDownload->setDisabled( !show_lines );
  
  
  const bool isPhone = ( m_spectrumViewer && m_spectrumViewer->isPhone() );
  const WString clearLineTxt = isPhone ? WString::tr(m_persisted.empty() ? "Remove" : "RemoveAll")
                                       : WString::tr( m_persisted.empty() ? "Clear" : "ClearAll" );
  if( clearLineTxt != m_clearLines->text() )
    m_clearLines->setText( clearLineTxt );
  
  
  bool showGammaCB = true, showXrayCb = true, showAplhaCb = true, showBetaCb = true, showEscapeCb = true;
  switch( src_type )
  {
    case ReferenceLineInfo::SourceType::Nuclide:
      // Show everything
      showGammaCB = showXrayCb = showAplhaCb = showBetaCb = showEscapeCb = true;
      break;
      
    case ReferenceLineInfo::SourceType::FluorescenceXray:
      // Only show x-ray option - but really we shouldnt show any of them
      showGammaCB = showAplhaCb = showBetaCb = showEscapeCb = false;
      break;
      
    case ReferenceLineInfo::SourceType::Reaction:
      // For reactions only show Show Gamma option
      showXrayCb = showAplhaCb = showBetaCb = false;
      break;
      
    case ReferenceLineInfo::SourceType::Background:
    case ReferenceLineInfo::SourceType::NuclideMixture:
      showAplhaCb = showBetaCb = showEscapeCb = false;
      break;
      
    case ReferenceLineInfo::SourceType::CustomEnergy:
    case ReferenceLineInfo::SourceType::OneOffSrcLines:
      // We treat custom energy and one-off sources as a gammas, so only show it - but really we shouldnt show any of them
      showXrayCb = showAplhaCb = showBetaCb = false;
      break;
      
    case ReferenceLineInfo::SourceType::None:
      // Show all options, otherwise whole area will be blank
      break;
  }//switch( src_type )
  
  m_showXrays->setHidden( !showXrayCb );
  m_showGammas->setHidden( !showGammaCB );
  m_showAlphas->setHidden( !showAplhaCb );
  m_showBetas->setHidden( !showBetaCb );
  m_showCascadeSums->setHidden( !ref_lines || !ref_lines->m_has_coincidences );
  m_showEscapes->setHidden( !showEscapeCb );
  
  if( ref_lines && (ref_lines->m_validity == ReferenceLineInfo::InputValidity::Valid) )
  {
    m_showXrays->setChecked( ref_lines->m_input.m_showXrays );
    m_showGammas->setChecked( ref_lines->m_input.m_showGammas );
    m_showAlphas->setChecked( ref_lines->m_input.m_showAlphas );
    m_showBetas->setChecked( ref_lines->m_input.m_showBetas );
    m_showCascadeSums->setChecked( ref_lines->m_input.m_showCascades );
    m_showEscapes->setChecked( ref_lines->m_input.m_showEscapes );
  }
  
  
  // Add current nuclide to lis of previous nuclides
  if( m_currentlyShowingNuclide.m_validity == ReferenceLineInfo::InputValidity::Valid )
  {
    RefLineInput prev = m_currentlyShowingNuclide.m_input;
    
    // Remove any other previous nuclides that have same `prev.m_nuclide` as what
    //  we are about to push on
    m_prevNucs.erase(std::remove_if(begin(m_prevNucs), end(m_prevNucs),
                                    [&prev](const RefLineInput &val) -> bool {
      return (val.m_input_txt == prev.m_input_txt);
    }), end(m_prevNucs));
    
    // Push new value onto front of history
    m_prevNucs.push_front( std::move(prev) );
    
    // Check history length, and truncate if needed
    if( m_prevNucs.size() > m_max_prev_nucs )
      m_prevNucs.resize( m_max_prev_nucs );
  }//if( !m_currentlyShowingNuclide.labelTxt.empty() )
  
  
  //Default to using the old color, unless we find a reason to change it.
  WColor color = ref_lines->m_input.m_color;
  
  //if( ref_lines->m_source_type == ReferenceLineInfo::SourceType::NuclideMixture )
  //{
   // TODO: maybe make each nuclide a different color... but this opens a can of worms
  //}//if( ref_lines->m_source_type == ReferenceLineInfo::SourceType::NuclideMixture )
  
  // Note: this is strictly the source string, and not age, shielding, or other options.
  const bool isSameSrc = (!srcstr.empty() && (srcstr == m_currentlyShowingNuclide.m_input.m_input_txt));
  if( !isSameSrc )
  {
    assert( ref_lines );
    color = colorForNewSource( srcstr );
  }//if( isSameSrc ) / else
  
  if( !isSameSrc )
    m_userHasPickedColor = false;
  ref_lines->m_input.m_color = color;
  
  if( color != m_colorSelect->color() )
    m_colorSelect->setColor( color );
  
  
  {// Begin handling shielding
    if( !user_input.m_shielding_name.empty() || !user_input.m_shielding_thickness.empty() )
    {
      const string new_mat = SpecUtils::trim_copy( user_input.m_shielding_name );
      const string new_thick = SpecUtils::trim_copy( user_input.m_shielding_thickness );
      
      WLineEdit *materialEdit = m_shieldingSelect->materialEdit();
      WLineEdit *thicknessEdit = m_shieldingSelect->thicknessEdit();
      string curr_mat = materialEdit ? materialEdit->text().toUTF8() : string();
      string curr_thick = thicknessEdit ? thicknessEdit->text().toUTF8() : string();
      SpecUtils::trim( curr_mat );
      SpecUtils::trim( curr_thick );
      
      if( !SpecUtils::iequals_ascii(curr_mat, new_mat)
          || !SpecUtils::iequals_ascii(curr_thick, new_thick) )
      {
        m_shieldingSelect->setMaterialNameAndThickness( new_mat, new_thick );
      }
    }else if( !user_input.m_shielding_an.empty() || !user_input.m_shielding_ad.empty() )
    {
#ifndef NDEBUG
      double dummy;
      assert( (std::stringstream(user_input.m_shielding_an) >> dummy) );
      assert( (std::stringstream(user_input.m_shielding_ad) >> dummy) );
#endif
      const string new_an = SpecUtils::trim_copy(user_input.m_shielding_an);
      const string new_ad = SpecUtils::trim_copy(user_input.m_shielding_ad);
      NativeFloatSpinBox *anEdit = m_shieldingSelect->atomicNumberEdit();
      NativeFloatSpinBox *adEdit = m_shieldingSelect->arealDensityEdit();
      string curr_an = anEdit ? anEdit->text().toUTF8() : string();
      string curr_ad = adEdit ? adEdit->text().toUTF8() : string();
      SpecUtils::trim( curr_an );
      SpecUtils::trim( curr_ad );
      
      if( !SpecUtils::iequals_ascii(curr_an, new_an)
         || !SpecUtils::iequals_ascii(curr_ad, new_ad) )
      {
        try
        {
          m_shieldingSelect->setAtomicNumberAndArealDensity( new_an, new_ad );
        }catch( std::exception &e )
        {
          cerr << "Error setting shielding: " << e.what() << endl;
          assert( 0 );
          m_shieldingSelect->setToNoShielding();
        }
      }//if( AN or AD needs changing )
    }else
    {
      m_shieldingSelect->setToNoShielding();
    }
  }// End handling shielding
  
  // Now check if m_options_icon should have a red background (non-default options) or not
  bool nonDefaultOpts = false;
  
  if( show_lines )
  {
    nonDefaultOpts |= ref_lines->m_input.m_showAlphas;
    nonDefaultOpts |= ref_lines->m_input.m_showBetas;
    
    switch( src_type )
    {
      case ReferenceLineInfo::SourceType::Nuclide:
        nonDefaultOpts |= ref_lines->m_input.m_promptLinesOnly;
        nonDefaultOpts |= !ref_lines->m_input.m_showXrays;
        nonDefaultOpts |= !ref_lines->m_input.m_showGammas;
        nonDefaultOpts |= ref_lines->m_input.m_showCascades;
        break;
        
      case ReferenceLineInfo::SourceType::FluorescenceXray:
        nonDefaultOpts |= !ref_lines->m_input.m_showXrays;
        break;
        
      case ReferenceLineInfo::SourceType::Reaction:
        nonDefaultOpts |= !ref_lines->m_input.m_showGammas;
        break;
        
      case ReferenceLineInfo::SourceType::Background:
      case ReferenceLineInfo::SourceType::NuclideMixture:
        nonDefaultOpts |= !ref_lines->m_input.m_showXrays;
        nonDefaultOpts |= !ref_lines->m_input.m_showGammas;
        break;
        
      case ReferenceLineInfo::SourceType::CustomEnergy:
      case ReferenceLineInfo::SourceType::OneOffSrcLines:
        nonDefaultOpts |= !ref_lines->m_input.m_showGammas;
        break;
        
      case ReferenceLineInfo::SourceType::None:
        nonDefaultOpts |= ref_lines->m_input.m_promptLinesOnly;
        nonDefaultOpts |= !ref_lines->m_input.m_showXrays;
        nonDefaultOpts |= !ref_lines->m_input.m_showGammas;
        nonDefaultOpts |= ref_lines->m_input.m_showCascades;
        break;
    }//switch( src_type )
    
    const bool showCascades = (ref_lines->m_has_coincidences && ref_lines->m_input.m_showCascades);
    if( showCascades && !m_cascadeWarn )
    {
      m_cascadeWarn = new WText( WString::tr("rpd-warn-cascade-xrays") );
      m_cascadeWarn->addStyleClass("CascadeGammaWarn");
      m_optionsContent->insertWidget( m_optionsContent->indexOf(m_showCascadeSums) + 1, m_cascadeWarn);
    }//if( show coincidences )
    
    if( m_cascadeWarn )
      m_cascadeWarn->setHidden( !showCascades );
  }//if( we are actually showing any lines )
  
  const bool hasNonDefaultStyle = m_options_icon->hasStyleClass("non-default");
  // Add or remove the .non-default style class, if necessary
  if( nonDefaultOpts != hasNonDefaultStyle )
    m_options_icon->toggleStyleClass("non-default", nonDefaultOpts);
  
  if( !show_lines )
  {
    m_currentlyShowingNuclide.reset();
    m_chart->setReferncePhotoPeakLines( {} );
    m_particleModel->clear();
  }else
  {
    assert( ref_lines );
    m_currentlyShowingNuclide = *ref_lines;
    m_chart->setReferncePhotoPeakLines( m_currentlyShowingNuclide );
    const vector<DecayParticleModel::RowData> table_rows = createTableRows( *ref_lines );
    m_particleModel->setRowData( table_rows );
  }

  updateOtherNucsDisplay();
  
  addUndoRedoPoint( starting_showing, starting_persisted, starting_prev_nucs,
                   starting_user_color, user_input );
}//void updateDisplayFromInput()



void ReferencePhotopeakDisplay::persistCurentLines()
{
  if( m_currentlyShowingNuclide.m_validity != ReferenceLineInfo::InputValidity::Valid )
    return;

  m_chart->persistCurrentReferncePhotoPeakLines();
  
  const string current_input = m_currentlyShowingNuclide.m_input.m_input_txt;
  
  ReferenceLineInfo *prev_ref = nullptr;
  for( size_t i = 0; !prev_ref && (i < m_persisted.size()); ++i )
  {
    if( m_persisted[i].m_input.m_input_txt == current_input )
      prev_ref = &(m_persisted[i]);
  }
  
  if( prev_ref )
    *prev_ref = m_currentlyShowingNuclide;
  else if( m_currentlyShowingNuclide.m_validity == ReferenceLineInfo::InputValidity::Valid )
    m_persisted.push_back( m_currentlyShowingNuclide );

  m_currentlyShowingNuclide.reset();
  m_userHasPickedColor = false;
  
  m_nuclideEdit->setText( "" );
  m_persistLines->disable();
  m_clearLines->enable();
  if( m_spectrumViewer && m_spectrumViewer->isPhone() )
    m_clearLines->setText( WString::tr("RemoveAll") );
  else
    m_clearLines->setText( WString::tr("ClearAll") );
  
  updateDisplayChange();
}//void persistCurentLines()


void ReferencePhotopeakDisplay::setColors( const std::vector<Wt::WColor> &referenceLineColor )
{
  m_lineColors.clear();
  m_lineColors.reserve( referenceLineColor.size() );
  for( const auto &i : referenceLineColor )
    if( !i.isDefault() )
      m_lineColors.push_back( i );
  
  if( m_lineColors.empty() )
    m_lineColors = ns_def_line_colors;
  
  //should update currently displayed line colors here, but there are bigger fish to fry ATM
}//void setColors( const std::vector<Wt::WColor> &referenceLineColor )


void ReferencePhotopeakDisplay::setColorsForSpecificSources( const std::map<std::string,Wt::WColor> &referenceLineColorForSources )
{
  m_specificSourcelineColors.clear();
  
  for( const auto &i : referenceLineColorForSources )
  {
    if( !i.second.isDefault() )
    m_specificSourcelineColors[i.first] = i.second;
  }
  
  //should update currently displayed line colors here, but there are bigger fish to fry ATM
}//void setColorsForSpecificSources( const std::map<std::string,Wt::WColor> &referenceLineColorForSources )


void ReferencePhotopeakDisplay::setPeaksGetAssignedRefLineColor( const bool theydo )
{
  m_peaksGetAssignedRefLineColor = theydo;
}


void ReferencePhotopeakDisplay::setExternalRidResults( const string &algo_name,
                                                      const vector<pair<string,string>> &nucs )
{
  const bool is_diff = (nucs != m_external_ids);
  m_external_ids = nucs;
  m_external_algo_name = algo_name;
  if( is_diff )
    updateOtherNucsDisplay();
}//void setExternalRidResults( const std::vector<std::string> &isotopes );


Wt::Signal<> &ReferencePhotopeakDisplay::displayingNuclide()
{
  return m_displayingNuclide;
}


Wt::Signal<> &ReferencePhotopeakDisplay::nuclidesCleared()
{
  return m_nuclidesCleared;
}

const MaterialDB *ReferencePhotopeakDisplay::materialDB() const
{
  return m_materialDB;
}


void ReferencePhotopeakDisplay::serialize( std::string &xml_data  ) const
{
  rapidxml::xml_document<char> doc;
  serialize( &doc );
  xml_data.clear();
  rapidxml::print(std::back_inserter(xml_data), doc, 0);
}//void serialize( std::string &xml_data  ) const


void ReferencePhotopeakDisplay::serialize(
                                  rapidxml::xml_node<char> *parent_node ) const
{
  rapidxml::xml_document<char> *doc = parent_node->document();
  
  const char *name, *value;
  rapidxml::xml_node<char> *base_node, *node, *element;
  rapidxml::xml_attribute<char> *attr;
  
  name = "ReferencePhotopeakDisplay";
  base_node = doc->allocate_node( rapidxml::node_element, name );
  parent_node->append_node( base_node );
  
  //If you change the available options or formatting or whatever, increment the
  //  version field of the XML!
  value = doc->allocate_string( std::to_string(sm_xmlSerializationVersion).c_str() );
  attr = doc->allocate_attribute( "version", value );
  base_node->append_attribute( attr );

  if( m_currentlyShowingNuclide.m_input.m_input_txt.empty() )
  {
    node = doc->allocate_node( rapidxml::node_element, "CurrentGuiState" );
    base_node->append_node( node );

    name = "Nuclide";
    value = doc->allocate_string( m_nuclideEdit->text().toUTF8().c_str() );
    element = doc->allocate_node( rapidxml::node_element, name, value );
    node->append_node( element );
    
    name = "Age";
    string age_str = m_ageEdit->text().toUTF8();
    value = doc->allocate_string( age_str.c_str() );
    element = doc->allocate_node( rapidxml::node_element, name, value );
    node->append_node( element );

    name = "LowestBranchRatio";
    value = doc->allocate_string( (m_lowerBrCuttoff ? m_lowerBrCuttoff->text().toUTF8().c_str() : "0.0") );
    element = doc->allocate_node( rapidxml::node_element, name, value );
    node->append_node( element );
    
    name = "PromptLinesOnly";
    value = (m_promptLinesOnly->isChecked() ? "1" : "0");
    element = doc->allocate_node( rapidxml::node_element, name, value );
    node->append_node( element );
    
    name = "ShowGammas";
    value = ((!m_showGammas || m_showGammas->isChecked()) ? "1" : "0");
    element = doc->allocate_node( rapidxml::node_element, name, value );
    node->append_node( element );

    name = "ShowXrays";
    value = ((!m_showXrays || m_showXrays->isChecked()) ? "1" : "0");
    element = doc->allocate_node( rapidxml::node_element, name, value );
    node->append_node( element );
    
    name = "ShowAlphas";
    value = (m_showAlphas->isChecked() ? "1" : "0");
    element = doc->allocate_node( rapidxml::node_element, name, value );
    node->append_node( element );
    
    name = "ShowBetas";
    value = (m_showBetas->isChecked() ? "1" : "0");
    element = doc->allocate_node( rapidxml::node_element, name, value );
    node->append_node( element );

    name = "ShowCascades";
    value = (m_showCascadeSums->isChecked() ? "1" : "0");
    element = doc->allocate_node(rapidxml::node_element, name, value);
    node->append_node(element);
    
    name = "ShowEscapes";
    value = (m_showEscapes->isChecked() ? "1" : "0");
    element = doc->allocate_node(rapidxml::node_element, name, value);
    node->append_node(element);
  }else
  {
    node = doc->allocate_node( rapidxml::node_element, "CurrentLines" );
    base_node->append_node( node );
    m_currentlyShowingNuclide.m_input.serialize( node );
  }//if( !m_currentlyShowingNuclide.empty() )
  
  if( !m_persisted.empty() )
  {
    node = doc->allocate_node( rapidxml::node_element, "PersistedLines" );
    base_node->append_node( node );
    for( const ReferenceLineInfo &n : m_persisted )
      n.m_input.serialize( node );
  }//if( !m_persisted.empty() )
  
  WColor color;
  if( m_currentlyShowingNuclide.m_input.m_input_txt.empty() )
    color = m_colorSelect->color();
  else
    color = m_currentlyShowingNuclide.m_input.m_color;
  
  if( !color.isDefault() )
  {
    name = "CurrentLineColor";
    value = doc->allocate_string( color.cssText(false).c_str() );
    element = doc->allocate_node( rapidxml::node_element, name, value );
    attr = doc->allocate_attribute( "UserSelected", (m_userHasPickedColor ? "1" : "0") );
    element->append_attribute( attr );
    base_node->append_node( element );
  }//if( !color.isDefault() )

  //I guess we dont actually need to serialize m_lineColors,
  //  m_peaksGetAssignedRefLineColor, or m_specificSourcelineColors.size as
  //  these should be captured in the ColorTheme.
  /*
  rapidxml::xml_node<char> *colorsToUse = doc->allocate_node( rapidxml::node_element, "ColorsToUse" );
  attr = doc->allocate_attribute( "PeaksGetAssignedRefLineColor", (m_peaksGetAssignedRefLineColor ? "1" : "0") );
  colorsToUse->append_attribute( attr );
  base_node->append_node( colorsToUse );
  
  for( const auto &c : m_lineColors )
  {
    value = doc->allocate_string( c.cssText(false).c_str() );
    element = doc->allocate_node( rapidxml::node_element, "Color", value );
    colorsToUse->append_node( element );
  }
  
  if( m_specificSourcelineColors.size() )
  {
    rapidxml::xml_node<char> *specificSrcColor = doc->allocate_node( rapidxml::node_element, "SpecificSrcColors" );
    base_node->append_node( specificSrcColor );
    
    for( const auto &p : m_specificSourcelineColors )
    {
      auto src = doc->allocate_string( p.first.c_str() );
      auto color = doc->allocate_string( p.second.cssText(false).c_str() );
      element = doc->allocate_node( rapidxml::node_element, "SrcColor" );
      attr = doc->allocate_attribute( "src", src );
      element->append_attribute( attr );
      attr = doc->allocate_attribute( "color", color );
      element->append_attribute( attr );
      specificSrcColor->append_node( element );
    }
  }//if( m_specificSourcelineColors.size() )
  
  if( m_previouslyPickedSourceColors.size() )
  {
    rapidxml::xml_node<char> *prevSrcColor = doc->allocate_node( rapidxml::node_element, "PreviousSrcColors" );
    base_node->append_node( prevSrcColor );
    
    for( const auto &p : m_previouslyPickedSourceColors )
    {
      auto src = doc->allocate_string( p.first.c_str() );
      auto color = doc->allocate_string( p.second.cssText(false).c_str() );
      element = doc->allocate_node( rapidxml::node_element, "SrcColor" );
      attr = doc->allocate_attribute( "src", src );
      element->append_attribute( attr );
      attr = doc->allocate_attribute( "color", color );
      element->append_attribute( attr );
      prevSrcColor->append_node( element );
    }
  }//if( m_previouslyPickedSourceColors.size() )
  */
  
  m_shieldingSelect->serialize( base_node );
}//void serialize( rapidxml::xml_document<char> &doc )


void ReferencePhotopeakDisplay::deSerialize( std::string &xml_data  )
{
  clearAllLines();
  
  try
  {
    rapidxml::xml_document<char> doc;
    const int flags = rapidxml::parse_normalize_whitespace
                      | rapidxml::parse_trim_whitespace;

    //cout << "xml_data=" << xml_data << endl;
    
    if( xml_data.size() )
      doc.parse<flags>( &(xml_data[0]) );
    
    rapidxml::xml_attribute<char> *attr;
    rapidxml::xml_node<char> *base_node, *node, *showing_node,
                             *gui_node, *persisted_node;
    
    base_node = doc.first_node( "ReferencePhotopeakDisplay", 25 );
    if( !base_node )
      base_node = doc.first_node( "PhotopeakLineDisplay", 20 );
    
    if( !base_node )
      throw runtime_error( "Couldnt get base node, ReferencePhotopeakDisplay or PhotopeakLineDisplay" );
    
    int version = 0;
    attr = base_node->first_attribute( "version", 7 );
    if( !attr || !attr->value()
        || !(stringstream(attr->value()) >> version)
        || (version != sm_xmlSerializationVersion) )
      throw runtime_error( "Missing or invalid ReferencePhotopeakDisplay version" );
    
    gui_node = base_node->first_node( "CurrentGuiState", 15 );
    showing_node = base_node->first_node( "CurrentLines", 12 );
    
    if( (gui_node && showing_node) || !(gui_node || showing_node) )
      throw runtime_error( "Inconsistent saving of photopeak lines state" );
    
    if( showing_node && showing_node->first_node("DisplayedSource",15) )
      gui_node = showing_node->first_node("DisplayedSource",15);
    
    if( gui_node )
    {
      const SandiaDecay::Nuclide *nuc = nullptr;
      node = gui_node->first_node( "Nuclide", 7 );
      if( node && node->value() )
      {
        const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
        if( db && node->value_size() )
          nuc = db->nuclide( node->value() );
        m_nuclideEdit->setText( node->value() );
      }
      
      node = gui_node->first_node( "Age", 3 );
      if( node && node->value() && nuc )
        m_ageEdit->setText( node->value() );
      else
        m_ageEdit->setText( "" );
      
      node = gui_node->first_node( "LowestBranchRatio", 17 );
      if( node && node->value() && m_lowerBrCuttoff )
        m_lowerBrCuttoff->setText( node->value() );
      
      node = gui_node->first_node( "PromptLinesOnly", 15 );
      if( node && node->value() && strlen(node->value()))
        m_promptLinesOnly->setChecked( (node->value()[0] == '1') );
      
      node = gui_node->first_node( "ShowGammas", 10 );
      if( node && node->value() && strlen(node->value()) && m_showGammas )
        m_showGammas->setChecked( (node->value()[0] == '1') );
      
      node = gui_node->first_node( "ShowXrays", 9 );
      if( node && node->value() && strlen(node->value()) && m_showXrays )
        m_showXrays->setChecked( (node->value()[0] == '1') );
      
      node = gui_node->first_node( "ShowAlphas", 10 );
      if( node && node->value() && strlen(node->value()))
        m_showAlphas->setChecked( (node->value()[0] == '1') );
      
      node = gui_node->first_node( "ShowBetas", 9 );
      if( node && node->value() && strlen(node->value()))
        m_showBetas->setChecked( (node->value()[0] == '1') );

      node = gui_node->first_node("ShowCascades", 12);
      if (node && node->value() && strlen(node->value()))
        m_showCascadeSums->setChecked((node->value()[0] == '1'));
      
      node = gui_node->first_node("ShowEscapes", 11);
      if (node && node->value() && strlen(node->value()))
        m_showEscapes->setChecked((node->value()[0] == '1'));
    }//if( gui_node )
    
    if( showing_node )
    {
      m_currentlyShowingNuclide.reset();
      node = showing_node->first_node( "DisplayedSource", 15 );
      
      if( !node )
        node = showing_node->first_node( "RefLineInput", 12 );
      
      if( node )
      {
        RefLineInput input;
        input.deSerialize( node );
        
        // Now we need to set DRF, and shielding att function
        const shared_ptr<const DetectorPeakResponse> drf = m_detectorDisplay
                                                             ? m_detectorDisplay->detector()
                                                             : nullptr;
        if( drf && drf->isValid() )
        {
          input.m_detector_name = drf->name();
          input.m_det_intrinsic_eff = drf->intrinsicEfficiencyFcn();
        }else
        {
          input.m_detector_name = "";
          input.m_det_intrinsic_eff = nullptr;
        }
        
        if( (!input.m_shielding_name.empty() && !input.m_shielding_thickness.empty())
           || (!input.m_shielding_an.empty() && !input.m_shielding_ad.empty()) )
        {
          input.setShieldingAttFcn( m_materialDB );
        }
        
        //shared_ptr<ReferenceLineInfo> ref_line = ReferenceLineInfo::generateRefLineInfo( input );
        //if( !ref_line || (ref_line->m_validity != ReferenceLineInfo::InputValidity::Valid) )
        //  throw runtime_error( "Couldn't generate reference lines from source '" + input.m_input_txt + "'" );
        //
        //m_currentlyShowingNuclide = *ref_line;
        
        updateDisplayFromInput( input );
      }
    }//if( showing_node )
    
    m_persisted.clear();
    persisted_node = base_node->first_node( "PersistedLines", 14 );
    
    if( persisted_node )
    {
      node = persisted_node->first_node( "DisplayedSource", 15 );
      if( !node )
        node = persisted_node->first_node( "RefLineInput", 12 );
      
      for( ; node; node = node->next_sibling( node->name(), node->name_size() ) )
      {
        RefLineInput input;
        input.deSerialize( node );
        
        shared_ptr<ReferenceLineInfo> ref_line = ReferenceLineInfo::generateRefLineInfo( input );
        if( !ref_line || (ref_line->m_validity != ReferenceLineInfo::InputValidity::Valid) )
          throw runtime_error( "Couldn't generate persisted reference lines from source '" + input.m_input_txt + "'" );
        m_persisted.push_back( *ref_line );
      }//for( loop over DisplayedSource nodes )
    }//if( persisted_node )
    
    node = base_node->first_node( "Shielding", 9 );
    if( node )
    {
      m_shieldingSelect->materialChanged().setBlocked( true );
      m_shieldingSelect->materialModified().setBlocked( true );
      const bool is_fixed_geom = false; //Shouldnt have an effect either way
      m_shieldingSelect->deSerialize( node, is_fixed_geom );
      m_shieldingSelect->materialChanged().setBlocked( false );
      m_shieldingSelect->materialModified().setBlocked( false );
    }
//    else
//      m_shieldingSelect-reset();
    
    auto currentColor = base_node->first_node( "CurrentLineColor", 16 );
    if( currentColor && currentColor->value_size() )
    {
      const string value( currentColor->value(), currentColor->value() + currentColor->value_size() );
      try{ m_colorSelect->setColor( WColor(value) ); }catch(...){}
      
      auto pickedAttrib = currentColor->first_attribute( "UserSelected", 12 );
      if( pickedAttrib && pickedAttrib->value_size() )
        m_userHasPickedColor = (pickedAttrib->value()[0] == '1');
    }else if( !m_currentlyShowingNuclide.m_input.m_color.isDefault() )
    {
      m_colorSelect->setColor( m_currentlyShowingNuclide.m_input.m_color );
    }
    
    //The quantities m_lineColors, m_peaksGetAssignedRefLineColor,
    //  m_previouslyPickedSourceColors, and m_specificSourcelineColors should
    //  get set by the color theme.
    
    refreshLinesDisplayedToGui();
    
//    updateDisplayChange();
    if( (m_currentlyShowingNuclide.m_validity != ReferenceLineInfo::InputValidity::Valid)
       && m_persisted.empty() )
    {
      m_nuclidesCleared.emit();
    }else
    {
      m_displayingNuclide.emit();
    }
  }catch( std::exception &e )
  {
    cerr << "ReferencePhotopeakDisplay::deSerialize() caught: " << e.what() << endl;
    stringstream msg;
    msg << "Error opening displayed photopeaks from database for display: " << e.what();
    passMessage( msg.str(), WarningWidget::WarningMsgHigh );
  }//try / catch
}//void deSerialize( std::string &xml_data  )


//std::string ReferencePhotopeakDisplay::jsonReferenceLinesArray()
//{
//  string answer = "[";
//  const bool include_primary = ( (m_currentlyShowingNuclide.m_validity == ReferenceLineInfo::InputValidity::Valid)
//                                && !m_currentlyShowingNuclide.m_ref_lines.empty() );
//  if( include_primary )
//    m_currentlyShowingNuclide.toJson(answer);
//
//  for( size_t i = 0; i < m_persisted.size(); ++i )
//  {
//    const ReferenceLineInfo &ref = m_persisted[i];
//
//    if( include_primary
//        || ref.m_input.m_input_txt != m_currentlyShowingNuclide.m_input.m_input_txt )
//    {
//      answer += (answer.size() > 1) ? "," : "";
//      ref.toJson(answer);
//    }
//  }
//  answer += "]";
//  return answer;
//}//jsonReferenceLines()


std::map<std::string,std::string> ReferencePhotopeakDisplay::jsonReferenceLinesMap()
{
  std::map<std::string,std::string> answer;
  
  if( m_currentlyShowingNuclide.m_input.m_input_txt != "" )
    m_currentlyShowingNuclide.toJson( answer[m_currentlyShowingNuclide.m_input.m_input_txt] );
    
  for( size_t i = 0; i < m_persisted.size(); ++i )
  {
    const ReferenceLineInfo &ref = m_persisted[i];
    if( ref.m_input.m_input_txt != m_currentlyShowingNuclide.m_input.m_input_txt )
      ref.toJson( answer[ref.m_input.m_input_txt] );
  }
  
  return answer;
}//std::map<std::string,std::string> jsonReferenceLinesMap();


void ReferencePhotopeakDisplay::refreshLinesDisplayedToGui()
{
  //Note that the millisecond delay _may_ be vestigulal (untested), from when
  //  this->doJavaScript() was being called instead of wApp->doJavaScript().
  
  m_chart->clearAllReferncePhotoPeakLines();
  
  for( size_t i = 0; i < m_persisted.size(); ++i )
  {
    m_chart->setReferncePhotoPeakLines( m_persisted[i] );
    m_chart->persistCurrentReferncePhotoPeakLines();
  }
  
  m_chart->setReferncePhotoPeakLines( m_currentlyShowingNuclide );
}//void refreshLinesDisplayedToGui()


void ReferencePhotopeakDisplay::setIsotope( const SandiaDecay::Nuclide *nuc,
                                       double age )
{
  if( nuc )
  {
    m_nuclideEdit->setText( nuc->symbol );
    if( age >= 0.0 )
      m_ageEdit->setText( PhysicalUnitsLocalized::printToBestTimeUnits(age) );
  }else
  {
    if( m_nuclideEdit->valueText().empty() )
      return;
    
    m_nuclideEdit->setText( "" );
  }//if( nuc ) / else
  
  handleIsotopeChange( true );
}//void setIsotope(...)


void ReferencePhotopeakDisplay::setNuclideAndAge( const string &name,
                                             const string &age )
{
  m_nuclideEdit->setText( name );
  m_ageEdit->setText( age );
  handleIsotopeChange( !age.empty() );
}//void setNuclideAndAge( const std::string &name, const std::string &age )


void ReferencePhotopeakDisplay::setShieldingMaterialAndThickness( const string &name,
                                                     const string &thickness )
{
  m_shieldingSelect->setMaterialNameAndThickness( name, thickness );
  updateDisplayChange();
}//void setShieldingMaterialAndThickness(...)


const ShieldingSelect *ReferencePhotopeakDisplay::shieldingSelect()
{
  return m_shieldingSelect;
}//const ShieldingSelect *shieldingSelect()


void ReferencePhotopeakDisplay::setElement( const SandiaDecay::Element *el )
{
  if( el )
  {
    m_nuclideEdit->setText( el->symbol );
  }else
  {
    m_nuclideEdit->setText( "" );
  }//if( nuc ) / else
  
  handleIsotopeChange( false );
}//void setIsotope(...)


void ReferencePhotopeakDisplay::setReaction( const ReactionGamma::Reaction *rctn )
{
  if( rctn )
  {
    m_nuclideEdit->setText( rctn->name() );
  }else
  {
    m_nuclideEdit->setText( "" );
  }//if( nuc ) / else
  
  handleIsotopeChange( false );
}//void setReaction(...)


void ReferencePhotopeakDisplay::fitPeaks()
{
  passMessage( "ReferencePhotopeakDisplay::fitPeaks() not implemented yet,"
               " but it will be really cool once it is", 2 );
  
  //-Get nuclide being displayed, including persisted, then do a whole lot of work...
  //-Get the existing peaks, and (temporarily fix their means).
  //-If the existing detector response (m_detectorDisplay->detector()) does not
  //   have resolution information, do a pre-fit of the highest amplitude
  //   candiates (above 100 keV, highest amp candidate in each 50 keV span
  //   maybye), as use this as a starting width response.  If this fales, than
  //   just use a generic NaI, HPGe, or LaBr resonse.
  //-Use this detector response to do the 'skyline' type filtering of potential
  //   peaks.  E.g. loop over each gamma line, and use the same guessing
  //   algorithm that asigning a nuclide line to a peak does, to see if that
  //   gamma line would be selected, if so, fit for that peak, if not, dont
  //   (btw, Ba133 window identifiaction method fails for NaI Ba133, so check
  //    that out and improve it)
  //-Now try to figure out which peaks should be fit in the same ROIs, and
  //   co-fit them, while indiviually fitting the others.  The existing peaks
  //   will have to be be tossed into the fit as well, but with fixed means.
  //   This will have to be re-evaluted after the initial fit.
  //-To actually assign the isotopes to the peaks, should consider the case
  //   of a single peak for multipe isotopes being shown: at first just identify
  //   candate nuclides, then disregard shielding and check the amplitude
  //   relative to the nearest candidate peak of each isotope (at an assumed
  //   age of course), and assign it to the peak that is closest.  Could also
  //   consider trying to make sure to assign at least some peaks to each
  //   nuclide.
  //-
  
  
//  m_spectrumViewer
  
}//void fitPeaks()

void ReferencePhotopeakDisplay::clearAllLines()
{
  //m_fitPeaks->disable();
  m_persisted.clear();
  m_currentlyShowingNuclide.reset();
  m_particleModel->clear();
  
  m_currentlyShowingNuclide.m_input.m_color = m_lineColors[0];
  m_colorSelect->setColor( m_lineColors[0] );

  m_nuclideEdit->setText( "" );
  m_persistLines->disable();
  m_clearLines->disable();
  
  // Reset the shielding as well
  if( m_shieldingSelect->isGenericMaterial() )
    m_shieldingSelect->setAtomicNumberAndArealDensity( m_shieldingSelect->atomicNumber(), 0.0 );
  else
    m_shieldingSelect->setSphericalThickness( 0.0 );
  
  if( m_spectrumViewer && m_spectrumViewer->isPhone() )
    m_clearLines->setText( WString::tr("Clear") );
  else
    m_clearLines->setText( WString::tr("Remove") );

  m_moreInfoBtn->hide();
  updateOtherNucsDisplay();

  m_chart->clearAllReferncePhotoPeakLines();
  
  m_nuclidesCleared.emit();
}//void clearAllLines()




void ReferencePhotopeakDisplay::userColorSelectCallback( const WColor &color )
{
  if( color.isDefault() )
  {
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, "Was passed an invalid color." );
#endif
    m_userHasPickedColor = false;
    m_colorSelect->setColor( m_currentlyShowingNuclide.m_input.m_color );
  }else
  {
    m_userHasPickedColor = true;
    m_currentlyShowingNuclide.m_input.m_color = color;
    m_previouslyPickedSourceColors[m_currentlyShowingNuclide.m_input.m_input_txt] = color;
    refreshLinesDisplayedToGui();
  }
}//void userColorSelectCallback( const std::string &color )



