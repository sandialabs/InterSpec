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
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeSelectionAids.h"
#include "InterSpec/IsotopeNameFilterModel.h"
#include "InterSpec/MoreNuclideInfoDisplay.h"
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
          out << "Source,Gammas" << eol_char;
          break;
          
        case ReferenceLineInfo::SourceType::None:
          assert( 0 );
          break;
      }//switch( refinfo.m_source_type )
      
      
      boost::function<double(float)> att_fcn = refinfo.m_input.m_shielding_att;
      
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
          assert( refinfo.m_input.m_shielding_thickness.empty() );
          
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
          out << eol_char << rel_amp_note << eol_char << eol_char;
          out << "Energy (keV),Rel. Yield" << eol_char;
          break;
          
        case ReferenceLineInfo::SourceType::Reaction:
          out << eol_char << rel_amp_note << eol_char << eol_char;
          out << "Energy (keV),Rel. Yield" << eol_char;
          break;
          
        case ReferenceLineInfo::SourceType::Background:
          out << eol_char << rel_amp_note << eol_char << eol_char;
          out << "Energy (keV),Rel. Yield";
          break;
          
        case ReferenceLineInfo::SourceType::CustomEnergy:
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
        
        if( line.m_transition && line.m_transition->parent )
          out << "," << line.m_transition->parent->symbol;
        else
          out << ",";
        
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
          snprintf( buffer, sizeof(buffer), "%1.7g,", line.m_decay_intensity*shield_eff*drf_eff );
          out << "," << buffer;
        }
        
        snprintf( buffer, sizeof(buffer), "%1.7g,", line.m_normalized_intensity );
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
        case SandiaDecay::IsometricTransitionDecay: return WString( "Iso" );
        case SandiaDecay::ElectronCaptureDecay:     return WString( "e.c." );
        case SandiaDecay::ProtonDecay:              return WString( "proton" );
        case SandiaDecay::SpontaneousFissionDecay:  return WString( "s.f." );
        case SandiaDecay::Carbon14Decay:            return WString( "C14" );
        case RowData::XRayDecayMode:                return WString( "xray" );
        case RowData::ReactionToGammaMode:          return WString( "Reaction" );
        case RowData::NormGammaDecayMode:           return WString( "NORM" );
        case RowData::CascadeSumMode:               return WString( "Cascade Sum" );
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
       case XrayParticle:            return WString( "xray" );
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
      case kEnergy:         return WString( "Energy (keV)" );
      case kBranchingRatio: return WString( "B.R." ); //"Phot/Decay" &gamma;/Decay  \u03B3/Decay  &#947;/Decay
      case kResponsibleNuc: return WString( "Parent" );
      case kDecayMode:      return WString( "Mode" );
      case kParticleType:   return WString( "Particle" );
      case kNumColumn:      return boost::any();
    }//switch( column )
  }else if( role == ToolTipRole )
  {
    switch( column )
    {
      case kEnergy:
        return WString( "Energy of the particle produced" );
      case kBranchingRatio:
        return WString( "Intensity of the particle, relative to the highest"
                        " intensity of that particle. For gammas, this is after"
                        " the optional shielding and detector effects are"
                        " applied." );
      case kResponsibleNuc:
        return WString( "Actual nuclide which decayed to give this particle" );
      case kDecayMode:
        return WString( "Decay mode of the parent nuclide, which produced this"
                        " particle" );
      case kParticleType:
        return WString( "The type of particle produced" );
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
    m_cascadeWarn( NULL ),
    m_showRiidNucs( NULL ),
    m_showPrevNucs( NULL ),
    m_showAssocNucs( NULL ),
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
    m_lineColors{ ns_def_line_colors }
{
  wApp->useStyleSheet("InterSpec_resources/ReferencePhotopeakDisplay.css");
  

  const char *tooltip = nullptr;
  
  m_currentlyShowingNuclide.reset();
  
  if( !chart )
    throw runtime_error( "ReferencePhotopeakDisplay: a valid chart"
                         " must be passed in" );

  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", specViewer );
  
  //The inputDiv/Layout is the left side of the widget that holds all the
  //  nuclide input,age, color picker, DRF, etc
  WContainerWidget *inputDiv = new WContainerWidget();
  WGridLayout *inputLayout = new WGridLayout();
  inputLayout->setContentsMargins( 0, 2, 0, 0 );
  inputDiv->setLayout( inputLayout );
    
    
  const bool isPhone = m_spectrumViewer->isPhone();
    
  if( isPhone )
  {
    addStyleClass( "ReferencePhotopeakDisplayMobile" );
    //m_layout->setHorizontalSpacing( 5 );
  }else
  {
    addStyleClass( "ReferencePhotopeakDisplay" );
  }//if( specViewer->isPhone() ) / else
  
  
  const WLength labelWidth(3.5,WLength::FontEm), fieldWidth(4,WLength::FontEm);
  const WLength optionWidth(5.25,WLength::FontEm), buttonWidth(5.25,WLength::FontEm);
  
  WLabel *nucInputLabel = new WLabel( "Nuclide:" );
  nucInputLabel->setMinimumSize( labelWidth, WLength::Auto );
  m_nuclideEdit = new WLineEdit( "" );
  m_nuclideEdit->setMargin( 1 );
  m_nuclideEdit->setMargin( 2, Wt::Side::Top );
//  m_nuclideEdit->setMinimumSize( WLength(10,WLength::FontEx), WLength::Auto );
  m_nuclideEdit->setMinimumSize( fieldWidth, WLength::Auto );
  m_nuclideEdit->setAutoComplete( false );
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
  m_persistLines = new WPushButton( "Add Another" );
  tooltip = "Keep the currently displayed lines and add a new nuclide/source to display.";
  HelpSystem::attachToolTipOn( m_persistLines, tooltip, showToolTips );
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
  
  tooltip = "ex. <b>U235</b>, <b>235 Uranium</b>, <b>U</b> (x-rays only)"
            ", <b>Uranium</b> (x-rays), <b>U-235m</b> (meta stable state)"
            ", <b>Cs137</b>, <b>background</b>, <b>H(n,g)</b>, etc.";
  HelpSystem::attachToolTipOn( m_nuclideEdit, tooltip, showToolTips );
  
  string replacerJs, matcherJs;
  IsotopeNameFilterModel::replacerJs( replacerJs );
  IsotopeNameFilterModel::nuclideNameMatcherJs( matcherJs );
  IsotopeNameFilterModel *isoSuggestModel = new IsotopeNameFilterModel( this );
  isoSuggestModel->addCustomSuggestPossibility( "background" );
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


  WLabel *ageInputLabel = new WLabel( "Age:" );
  m_ageEdit = new WLineEdit( "" );
  WRegExpValidator *validator = new WRegExpValidator( PhysicalUnits::sm_timeDurationHalfLiveOptionalRegex, this );
  validator->setFlags(Wt::MatchCaseInsensitive);
  m_ageEdit->setValidator(validator);
  m_ageEdit->setAutoComplete( false );
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
    m_clearLines = new WPushButton( "Remove" );
    else
      m_clearLines = new WPushButton( "Clear" );
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
  
  
  tooltip = "<div>Age can be specified using a combination of time units, "
  "similar to '<b>5.3y 8d 22m</b>' or in half lives like "
  "'<b>2.5 HL</b>'.</div>"
  "<div>"
  "Acceptible time units: <b>year</b>, <b>yr</b>, <b>y</b>, <b>day</b>, <b>d</b>, <b>hrs</b>, <b>hour</b>, <b>h</b>, <b>minute</b>, "
  "<b>min</b>, <b>m</b>, <b>second</b>, <b>s</b>, <b>ms</b>, <b>microseconds</b>, <b>us</b>, <b>nanoseconds</b>, <b>ns</b>, or "
  "you can specify time period by <b>hh:mm:ss</b>. Half life units can be "
  "specified using <b>hl</b>, <b>halflife</b>, <b>halflives</b>, <b>half-life</b>, <b>half-lives</b>, "
  "<b>half lives</b>, or <b>half life</b>."
  "</div>"
  "<div>"
  "Half life units or time periods can not be mixed with "
  "other units. When multiple time periods are "
  "specified, they are summed, e.x. '1y6months 3m' is interpreted as "
  "18 months and 3 minutes"
  "</div>";
  
  HelpSystem::attachToolTipOn( m_ageEdit, tooltip, showToolTips );
  
  
  tooltip = "Clears all persisted lines, as well as the current non-persisted"
  " lines.";
  HelpSystem::attachToolTipOn( m_clearLines, tooltip, showToolTips );
  
  
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

  m_moreInfoBtn = new WPushButton( "more info", hlRow );
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
  
  //ToDo: When the detector is changed, should actually call a specialized function
  //      so all of the persisted gamma lines will be changed...
  //      Also, it doesnt look like changing the DRF changes current showing primary lines...
  
  specViewer->detectorChanged().connect( boost::bind( &ReferencePhotopeakDisplay::handleIsotopeChange, this, true ) );
  specViewer->detectorModified().connect( boost::bind( &ReferencePhotopeakDisplay::handleIsotopeChange, this, true ) );
  
  // If foreground spectrum _file_ changes, then the RIID analysis results of the siggested  
  //  nuclides may need updating.  However, for simplicity, we'll update suggested nuclides 
  //  whenever the foreground gets updated
  specViewer->displayedSpectrumChanged().connect(this, &ReferencePhotopeakDisplay::handleSpectrumChange);

  lowerInputLayout->addWidget( m_detectorDisplay, 1, 0 );

  m_shieldingSelect = new ShieldingSelect( m_materialDB, m_materialSuggest );
  m_shieldingSelect->materialEdit()->setEmptyText( "<shielding material>" );
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
  WText *txt = new WText( "Options", closerow );
  WContainerWidget* closeIcon = new WContainerWidget(closerow);
  closeIcon->addStyleClass("closeicon-wtdefault");
  closeIcon->clicked().connect(this, &ReferencePhotopeakDisplay::toggleShowOptions);

  m_optionsContent = new WContainerWidget( m_options );
  m_optionsContent->addStyleClass( "ToolTabTitledColumnContent" );

  m_promptLinesOnly = new WCheckBox("Prompt Only", m_optionsContent );  //ɣ
  
  tooltip = "Gammas from only the original nuclide, and the descendants until one"
    " of them has a longer half-life than the original nuclide; the"
    " decay chain is in equilibrium till that point.";
  HelpSystem::attachToolTipOn(m_promptLinesOnly, tooltip, showToolTips);
  m_promptLinesOnly->checked().connect(this, &ReferencePhotopeakDisplay::updateDisplayChange);
  m_promptLinesOnly->unChecked().connect(this, &ReferencePhotopeakDisplay::updateDisplayChange);
  m_promptLinesOnly->hide();

  m_showGammas = new WCheckBox( "Show Gammas", m_optionsContent );
  m_showXrays = new WCheckBox( "Show X-rays", m_optionsContent );
  m_showAlphas = new WCheckBox( "Show Alphas", m_optionsContent );
  m_showBetas = new WCheckBox( "Show Betas", m_optionsContent );
  m_showCascadeSums = new WCheckBox("Cascade Sums", m_optionsContent );
  m_showCascadeSums->hide();
  
  m_showPrevNucs = new WCheckBox("Prev Nucs", m_optionsContent );
  m_showRiidNucs = new WCheckBox("Det RID Nucs", m_optionsContent );
  m_showAssocNucs = new WCheckBox("Assoc. Nucs", m_optionsContent );

  m_showGammas->setWordWrap( false );
  m_showXrays->setWordWrap( false );
  m_showAlphas->setWordWrap( false );
  m_showBetas->setWordWrap( false );
  m_showCascadeSums->setWordWrap( false );
  m_showPrevNucs->setWordWrap( false );
  m_showRiidNucs->setWordWrap( false );
  m_showAssocNucs->setWordWrap( false );

  m_showPrevNucs->checked().connect( this, &ReferencePhotopeakDisplay::updateOtherNucsDisplay );
  m_showPrevNucs->unChecked().connect( this, &ReferencePhotopeakDisplay::updateOtherNucsDisplay );
  m_showRiidNucs->checked().connect( this, &ReferencePhotopeakDisplay::updateOtherNucsDisplay );
  m_showRiidNucs->unChecked().connect( this, &ReferencePhotopeakDisplay::updateOtherNucsDisplay );
  m_showAssocNucs->checked().connect( this, &ReferencePhotopeakDisplay::updateOtherNucsDisplay );
  m_showAssocNucs->unChecked().connect( this, &ReferencePhotopeakDisplay::updateOtherNucsDisplay );


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
  
  m_otherNucsColumn = new WContainerWidget();

  m_otherNucsColumn->addStyleClass("OtherNucs ToolTabSection ToolTabTitledColumn");

  WText *otherNucTitle = new WText("Suggestions", m_otherNucsColumn);
  otherNucTitle->addStyleClass("ToolTabColumnTitle");

  m_otherNucs = new WContainerWidget(m_otherNucsColumn);
  m_otherNucs->addStyleClass( "OtherNucsContent ToolTabTitledColumnContent" );



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
    passMessage( "The 'Nuclide Decay Info' tool provides additional capabilities for"
                 " exporting nuclide, decay products, and decay information to CSV files.",
                 WarningWidget::WarningMsgInfo );
    //TODO: check about calling WarningWidget::displayPopupMessageUnsafe( msg, level, 5000 ); directly with a longer time for the message to hang around
  }));
  
  csvButton->setText( "CSV" );
  csvButton->disable();
  m_csvDownload = csvButton;


  WGridLayout *overallLayout = new WGridLayout();
  overallLayout->setContentsMargins( 0, 0, 0, 0 );
  setLayout( overallLayout );

  overallLayout->addWidget( inputDiv,          0, 0 );
  overallLayout->addWidget( lowerInput,        1, 0 );
  overallLayout->addWidget( m_options,         0, 1, 3, 1 );
  overallLayout->addWidget( m_otherNucsColumn, 0, 2, 3, 1 );
  overallLayout->addWidget( m_particleView,    0, 3, 3, 1 );
  overallLayout->addWidget( bottomRow,         2, 0 );

  overallLayout->setRowStretch( 2, 1 );
  overallLayout->setColumnStretch( 3, 1 );
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
        WString hlstr = PhysicalUnits::printToBestTimeUnits(
                                          5.0*nuc->promptEquilibriumHalfLife(),
                                          2, SandiaDecay::second );
        m_ageEdit->setText( hlstr );
      }else if( m_currentlyShowingNuclide.m_input.m_input_txt != nuc->symbol )
      {
        string agestr;
        PeakDef::defaultDecayTime( nuc, &agestr );
        m_ageEdit->setText( agestr );
      }else
      {
        const double hl = (nuc ? nuc->halfLife : -1.0);
        double age = PhysicalUnits::stringToTimeDurationPossibleHalfLife( agestr, hl );
        if( age > 100.0*nuc->halfLife || age < 0.0 )
          throw std::runtime_error( "" );
      }//if( nuc->decaysToStableChildren() ) / else
    }else if( nuc && useCurrentAge )
    {
      double age = PhysicalUnits::stringToTimeDurationPossibleHalfLife( agestr, nuc->halfLife );
      if( age > 100.0*nuc->halfLife || age < 0.0 )
        throw std::runtime_error( "" );
    }
  }catch(...)
  {
    if( nuc )
    {
      string defagestr;
      PeakDef::defaultDecayTime( nuc, &defagestr );
      passMessage( "Changed age to a more reasonable value for " + nuc->symbol
                   + " from '" + agestr + "' to '" + defagestr + "'",
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
  if (m_options->isHidden())
  {
    m_options->show();
    //m_options->animateShow(WAnimation(WAnimation::AnimationEffect::Pop, WAnimation::TimingFunction::Linear, 250) );
    m_options_icon->addStyleClass("active");
  }else
  {
    m_options->hide();
    //m_options->animateHide(WAnimation(WAnimation::AnimationEffect::Pop, WAnimation::TimingFunction::Linear, 250));
    m_options_icon->removeStyleClass("active");
  }
}//void toggleShowOptions()


void ReferencePhotopeakDisplay::handleSpectrumChange(SpecUtils::SpectrumType type)
{
  // For simplicity we'll update the suggested nuclides whenever foreground gets
  //  changed, even if its only because of sample number change (which will have
  //  same RIID results, as they are file-spoecific, not sample specific)
  if( type == SpecUtils::SpectrumType::Foreground )
    updateOtherNucsDisplay();
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

  WText *header = new WText( "Assoc. Nucs" );
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
      OtherNuc nuc;
      nuc.m_input = userInput();
      nuc.m_input.m_input_txt = nucstr;
      nuc.m_nuclide = nucstr;
      
      btn->clicked().connect( boost::bind( &ReferencePhotopeakDisplay::setFromOtherNuc, this, nuc ) );
    }else
    {
      btn->disable();
    }
  }//for( const auto &riid : riid_nucs )
}//void updateAssociatedNuclides()


void ReferencePhotopeakDisplay::showMoreInfoWindow()
{
  const SandiaDecay::Nuclide * const nuc = m_currentlyShowingNuclide.m_nuclide;
  assert( nuc );
  if( !nuc )
    return;

  new MoreNuclideInfoWindow( nuc );
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

  m_otherNucsColumn->show();

  vector<OtherNuc> prev_nucs;
  const string &currentInput = m_currentlyShowingNuclide.m_input.m_input_txt;

  if( showPrev )
  {
    for( const auto &prev : m_prevNucs )
    {
      if( prev.m_input.m_input_txt != currentInput )
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
        riid_nucs.push_back( {res.nuclide_, nuc->symbol} );
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

  if( !riid_nucs.empty() )
  {
    // Protect against a detector having a ton of results
    const size_t max_riid_res = 10;
    if( riid_nucs.size() > max_riid_res )
      riid_nucs.resize( max_riid_res );

    WText *header = new WText("Detector ID", m_otherNucs);
    header->addStyleClass("OtherNucTypeHeader");
    
    for( const auto &riid : riid_nucs )
    {
      WPushButton *btn = new WPushButton(riid.first, m_otherNucs);
      btn->addStyleClass( "LinkBtn" );
      
      if( !riid.second.empty() )
      {
        OtherNuc nuc;
        nuc.m_nuclide = riid.second;
        nuc.m_input = userInput();
        btn->clicked().connect(boost::bind(&ReferencePhotopeakDisplay::setFromOtherNuc, this, nuc));
      } else
      {
        btn->disable();
      }
    }//for( const auto &riid : riid_nucs )
  }//if( !prev_nucs.empty() )


  if( !prev_nucs.empty() )
  {
    WText *header = new WText("Previous", m_otherNucs);
    header->addStyleClass("OtherNucTypeHeader");
    for( const auto &prev : prev_nucs )
    {
      WPushButton *btn = new WPushButton(prev.m_nuclide, m_otherNucs);
      btn->addStyleClass( "LinkBtn" );

      btn->clicked().connect(boost::bind(&ReferencePhotopeakDisplay::setFromOtherNuc, this, prev));
    }
  }//if( !prev_nucs.empty() )
}//void updateOtherNucsDisplay()


void ReferencePhotopeakDisplay::setFromOtherNuc(const ReferencePhotopeakDisplay::OtherNuc &nuc)
{
  updateDisplayFromInput( nuc.m_input );
}//void setFromOtherNuc(const OtherNuc &nuc)



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
      const float thick = static_cast<float>( m_shieldingSelect->thickness() );
      
      if( material && (thick > 0.0) )
      {
        input.m_shielding_name = material->name;
        input.m_shielding_thickness = m_shieldingSelect->thicknessEdit()->text().toUTF8();
        
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


std::shared_ptr<ReferenceLineInfo> ReferencePhotopeakDisplay::generateRefLineInfo( RefLineInput input )
{
  // The gamma or xray energy below which we wont show lines for.
  //  x-rays for nuclides were limited at above 10 keV, so we'll just impose this
  //  as a lower limit to show to be consistent.
  const float lower_photon_energy = 10.0f;

  auto answer_ptr = make_shared<ReferenceLineInfo>();
  ReferenceLineInfo &answer = *answer_ptr;
  

  // We want to set the final _modified_ version of input to the answer, before
  // returning, so we'll just use a helper for this.
  //on_scope_exit on_exit( [&answer_ptr, &input](){
  //  answer_ptr->m_input = input;
  //  answer_ptr->lineColor = input.m_color;
  //  } );
  answer_ptr->m_input = input;

  if( input.m_input_txt.empty() )
  {
    answer.m_validity = ReferenceLineInfo::InputValidity::Blank; //Should already be this value, but being explicit

    return answer_ptr;
  }

  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  double age = 0.0;
  const SandiaDecay::Nuclide * const nuc = db->nuclide( input.m_input_txt );

  if( !nuc )
  {
    input.m_age = "";
  }else
  {
    input.m_input_txt = nuc->symbol;
    input.m_promptLinesOnly = (input.m_promptLinesOnly && nuc->canObtainPromptEquilibrium());
   
    answer.m_nuclide = nuc;
    
    if( input.m_promptLinesOnly )
    {
      age = 5.0*nuc->promptEquilibriumHalfLife();
      input.m_age = PhysicalUnits::printToBestTimeUnits(age, 2);
    }else if( input.m_age == "" )
    {
      age = PeakDef::defaultDecayTime( nuc, &input.m_age );
    }else if( nuc->isStable() )
    {
      age = 0;
      input.m_age = "";
      answer.m_input_warnings.push_back( nuc->symbol + " is a stable isotope." );
      answer.m_validity = ReferenceLineInfo::InputValidity::InvalidSource;
      answer_ptr->m_input = input;
      
      return answer_ptr;
    }else
    {
      try
      {
        age = PhysicalUnits::stringToTimeDurationPossibleHalfLife( input.m_age, nuc->halfLife );
      } catch( std::exception & )
      {
        answer.m_input_warnings.push_back( "Invalid nuclide age input." );
        answer.m_validity = ReferenceLineInfo::InputValidity::InvalidAge;

        answer_ptr->m_input = input;

        return answer_ptr;
      }//try /catch to get the age

      if( age > 100.0 * nuc->halfLife || age < 0.0 )
      {
        const string old_age_str = input.m_age;
        age = PeakDef::defaultDecayTime( nuc, &input.m_age );
        answer.m_input_warnings.push_back( "Changed age to a more reasonable value for " 
          + nuc->symbol + " from '" + old_age_str + "' to " + input.m_age );
      }
    }//if( prompt only ) / else
  }//if( nuc )

  
  const bool check_element = (!nuc && (input.m_input_txt.find_first_of( "0123456789" ) == string::npos));
  const SandiaDecay::Element * const el = check_element ? db->element( input.m_input_txt ) : nullptr;
  
  if( el )
  {
    input.m_input_txt = el->symbol;
    input.m_age = ""; //JIC
    input.m_showXrays = true;
    
    answer.m_element = el;
  }//if( !nuc )


  string reaction_txt; //CSV list of reactions - I think for our context, only ever a single reaction
  vector<ReactionGamma::ReactionPhotopeak> rctn_gammas;
  if( !nuc && !el )
  {
    const size_t open_paren = input.m_input_txt.find( "(" );
    const size_t close_paren = (open_paren == string::npos) ? string::npos 
                                            : input.m_input_txt.find( ")", open_paren );

    if( close_paren != string::npos )
    {
      try
      {
        const ReactionGamma *rctnDb = ReactionGammaServer::database();
        if( rctnDb )
        {
          reaction_txt = rctnDb->gammas( input.m_input_txt, rctn_gammas );
          SpecUtils::ireplace_all( reaction_txt, "'", "" );

          // We will fill in answer.m_reactions later on.
          
          // Note: we are-not setting the input to the reaction_txt, as we dont need to bother
          //       the user that the underlying data is for the isotopics of the element, and also
          //       we dont want to change the user input, because we want to keep contributions
          //       normalized to natural abundance.
          //input.m_input_txt = reaction_txt;
          
          input.m_age = "";
          input.m_showGammas = true;
        }
      }catch( std::exception &e )
      {
        // Not a reaction
      }//try / catch
    }//if( (open_paren != string::npos) && (close_paren != string::npos) )
  }//if( !nuc && !el )


  vector<OtherRefLine> otherRefLinesToShow;
  const bool is_background = (nuc || el || !rctn_gammas.empty()) ? false
                                : SpecUtils::icontains( input.m_input_txt, "background" );
  if( is_background )
  {
    input.m_input_txt = "background";

    for( const OtherRefLine &bl : BackgroundLines )
    {
      const bool isXray = (std::get<3>(bl) == OtherRefLineType::BackgroundXRay);
      if( (isXray && input.m_showXrays) || (!isXray && input.m_showGammas) )
        otherRefLinesToShow.push_back( bl );
    }//for( const BackgroundLine &bl : BackgroundLines )

    input.m_age = "";
  }//if( is_background )

  bool is_custom_energy = false;
  if( !nuc && !el && rctn_gammas.empty() && !is_background )
  {
    try
    {
      const float energy = static_cast<float>(PhysicalUnits::stringToEnergy( input.m_input_txt ));

      //BackgroundLine<Energy, RelBranchRatio, "Symbol", OtherRefLineType, "Description">
      OtherRefLine line{energy, 1.0f, "", OtherRefLineType::OtherBackground, input.m_input_txt};
      otherRefLinesToShow.push_back( line );
      is_custom_energy = true;
      input.m_age = "";
    }catch( std::exception & )
    {
    }
  }//if( !nuc && !el && rctnGammas.empty() )
  

  if( !nuc && !el && rctn_gammas.empty() && !is_background && !is_custom_energy )
  {
    answer.m_validity = ReferenceLineInfo::InputValidity::InvalidSource;
    answer.m_input_warnings.push_back( input.m_input_txt + " is not a valid isotope, element, reaction, or energy." );
    
    answer_ptr->m_input = input;
    
    return answer_ptr;
  }//if( we couldnt match input text to a source )

  answer.m_validity = ReferenceLineInfo::InputValidity::Valid;

  input.m_showGammas = input.m_showGammas;
  input.m_showXrays = input.m_showXrays;
  input.m_showAlphas = input.m_showAlphas;
  input.m_showBetas  = input.m_showBetas;
  // We will also update showing cascades later based on answer.m_has_coincidences
  input.m_showCascades = input.m_showCascades;
  
  answer_ptr->m_input = input;

  
  if( nuc )
    answer.m_source_type = ReferenceLineInfo::SourceType::Nuclide;
  else if( el )
    answer.m_source_type = ReferenceLineInfo::SourceType::FluorescenceXray;
  else if( !rctn_gammas.empty() )
    answer.m_source_type = ReferenceLineInfo::SourceType::Reaction;
  else if( is_background )
    answer.m_source_type = ReferenceLineInfo::SourceType::Background;
  else if( is_custom_energy )
    answer.m_source_type = ReferenceLineInfo::SourceType::CustomEnergy;
  else
    answer.m_source_type = ReferenceLineInfo::SourceType::None;

  bool use_particle[SandiaDecay::ProductType::XrayParticle + 1] = {false};

  // We'll loop + switch over SandiaDecay::ProductType so we'll at least get
  //  a compiler warning if SandiaDecay::ProductType changes.
  for( auto type = SandiaDecay::ProductType(0); 
    type <= SandiaDecay::ProductType::XrayParticle; 
    type = SandiaDecay::ProductType(type+1) )
  {
    switch( type )
    {
      case SandiaDecay::ProductType::BetaParticle:
        use_particle[type] = input.m_showBetas;
        break;

      case SandiaDecay::ProductType::GammaParticle:
        use_particle[type] = (input.m_showGammas || input.m_showCascades);
        break;

      case SandiaDecay::ProductType::AlphaParticle:
        use_particle[type] = input.m_showAlphas;
        break;

      case SandiaDecay::ProductType::PositronParticle:
        use_particle[type] = input.m_showGammas;
        break;

      case SandiaDecay::ProductType::CaptureElectronParticle:
        break;

      case SandiaDecay::ProductType::XrayParticle:
        use_particle[type] = (el || input.m_showXrays);
        break;
    }//switch(type)
  }//for( loop over SandiaDecay::ProductType )


  vector<ReferenceLineInfo::RefLine> lines;

  //transition, first gamma BR, first gamma energy, second gamma energy, coincidence fraction, second gamma BR (just for debug)
  vector<tuple<const SandiaDecay::Transition *, double, float, float, float, double>> gamma_coincidences;

  if( nuc )
  {
    SandiaDecay::NuclideMixture mixture;

    if( input.m_promptLinesOnly )
    {
      age = 0.0;
      mixture.addNuclideInPromptEquilibrium( nuc, 1.0E-3 * SandiaDecay::curie );
    } else
    {
      mixture.addNuclideByActivity( nuc, 1.0E-3 * SandiaDecay::curie );
    }//if( we want promt only ) / else


    const vector<SandiaDecay::NuclideActivityPair> activities = mixture.activity( age );

    // x-rays are slightly problematic - we can *almost* treat them like gammas, but for 
    //  some decays they essentually get duplicated - so instead we'll be a little ineffient
    //  and track them seperately.
    vector<SandiaDecay::EnergyRatePair> xrays = use_particle[SandiaDecay::ProductType::XrayParticle]
      ? mixture.xrays( age )
      : vector<SandiaDecay::EnergyRatePair>{};

    const double parent_activity = nuc ? mixture.activity( age, nuc ) : 0.0;

    // We will accumulate positrons as a single line, and just assign the transition
    //  as the first one we run into
    //  TODO: should make sure the first transition we run into is the heaviest one.
    ReferenceLineInfo::RefLine positron_line;
    positron_line.m_energy = 510.9989 * PhysicalUnits::keV;
    positron_line.m_decay_intensity = 0.0;
    positron_line.m_parent_nuclide = nuc;
    positron_line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
    positron_line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Annihilation;


    for( const SandiaDecay::NuclideActivityPair &nap : activities )
    {
      const SandiaDecay::Nuclide *nuclide = nap.nuclide;
      const double activity = nap.activity;

      for( const SandiaDecay::Transition *transition : nuclide->decaysToChildren )
      {
        for( const SandiaDecay::RadParticle &particle : transition->products )
        {
          assert( particle.type <= SandiaDecay::ProductType::XrayParticle );

          if( !use_particle[particle.type] )
            continue;

          if( ((particle.type == SandiaDecay::GammaParticle)
            || (particle.type == SandiaDecay::XrayParticle))
            && (particle.energy < lower_photon_energy) )
          {
            continue;
          }

          if( particle.type == SandiaDecay::PositronParticle )
          {
            if( !positron_line.m_transition )
              positron_line.m_transition = transition;

            const double br = activity * particle.intensity
              * transition->branchRatio / parent_activity;

            positron_line.m_decay_intensity += 2.0 * br;

            continue;
          }//if( particle.type == SandiaDecay::PositronParticle )


          ReferenceLineInfo::RefLine line;
          line.m_parent_nuclide = nuc;
          line.m_energy = particle.energy;
          line.m_transition = transition;

          if( particle.type == SandiaDecay::XrayParticle )
          {
            size_t index = 0;
            for( ; index < xrays.size(); ++index )
            {
              if( fabs( xrays[index].energy - particle.energy ) < 1.0E-6 )
                break;
            }

            if( index < xrays.size() )
            {
              line.m_decay_intensity = xrays[index].numPerSecond / parent_activity;
              line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Xray;
              line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;

              lines.push_back( line );

              // Erase this x-ray so we dont double-count it
              xrays.erase( begin( xrays ) + index );
            } else
            {
              // We've already accounted for this energy.
            }

            continue;
          }//if( particle.type == SandiaDecay::XrayParticle )


          if( !answer.m_has_coincidences && (particle.type == SandiaDecay::GammaParticle) )
            answer.m_has_coincidences = !particle.coincidences.empty();

          const double br = activity * particle.intensity
            * transition->branchRatio / parent_activity;

          if( input.m_showCascades && (particle.type == SandiaDecay::GammaParticle) )
          {
            for( size_t coinc_index = 0; coinc_index < particle.coincidences.size(); ++coinc_index )
            {
              const unsigned short int part_ind = particle.coincidences[coinc_index].first;
              const float fraction = particle.coincidences[coinc_index].second;
              assert( part_ind < transition->products.size() );
              if( part_ind < transition->products.size() )
              {
                const SandiaDecay::RadParticle &coinc_part = transition->products[part_ind];

                // The BR of second gamma is just for debugging
                const double second_br = activity * coinc_part.intensity
                  * transition->branchRatio / parent_activity;

                if( coinc_part.type == SandiaDecay::ProductType::GammaParticle )
                  gamma_coincidences.emplace_back( transition, br, particle.energy, coinc_part.energy, fraction, second_br );
              }//if (part_ind < transition->products.size())
            }//for( loop over coincidences )
          }//if( show cascade gammas )


          // If type is GammaParticle, we could be here if user selected to show cascade, but
          //  not actual gammas, so we need to check for this, and if so not add this gamma 
          //  in to be shown
          if( (particle.type == SandiaDecay::GammaParticle) && !input.m_showGammas )
            continue;

          line.m_decay_intensity = br;
          line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;

          switch( particle.type )
          {
            case SandiaDecay::ProductType::BetaParticle:
              line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Beta;
              break;

            case SandiaDecay::ProductType::GammaParticle:
              line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
              break;

            case SandiaDecay::ProductType::AlphaParticle:
              line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Alpha;
              break;

            case SandiaDecay::ProductType::PositronParticle:
            case SandiaDecay::ProductType::CaptureElectronParticle:
            case SandiaDecay::ProductType::XrayParticle:
              assert( 0 );
              continue;
              break;
          }//switch( particle.type )

          lines.push_back( line );
        }//for( const SandiaDecay::RadParticle &particle : transition->products )
      }//for( const SandiaDecay::Transition *transition : nuclide->decaysToChildren )
    }//for( const SandiaDecay::NuclideActivityPair &nap : activities )

    if( positron_line.m_decay_intensity > 0.0 )
      lines.push_back( positron_line );
  }//if( nuc )

  // Update showing cascades based on if there are actually any present
  input.m_showCascades = (answer.m_has_coincidences && input.m_showCascades);
  answer.m_input.m_showCascades = (answer.m_has_coincidences && input.m_showCascades);
  
  
  
  if( el )
  {
    for( const SandiaDecay::EnergyIntensityPair &eip : el->xrays )
    {
      if( eip.energy < lower_photon_energy )
        continue;

      ReferenceLineInfo::RefLine line;
      line.m_element = el;
      line.m_energy = eip.energy;
      line.m_decay_intensity = eip.intensity;
      line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Xray;
      line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;

      lines.push_back( line );
    }//for( const SandiaDecay::EnergyIntensityPair &eip : element->xrays )
  }//if( m_showXrays->isChecked() )


  if( !rctn_gammas.empty() )
  {
    for( const ReactionGamma::ReactionPhotopeak &eip : rctn_gammas )
    {
      if( eip.reaction )
        answer.m_reactions.insert( eip.reaction );
      
      if( eip.energy < lower_photon_energy )
        continue;

      ReferenceLineInfo::RefLine line;
      line.m_reaction = eip.reaction;
      line.m_energy = eip.energy;
      line.m_decay_intensity = eip.abundance;
      line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
      line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;

      lines.push_back( line );
    }//for( const SandiaDecay::EnergyIntensityPair &eip : element->xrays )
  }//if( !rctn_gammas.empty() )


  if( !otherRefLinesToShow.empty() )
  {
    assert( is_background || is_custom_energy );
    assert( !nuc && !el && rctn_gammas.empty() );

    const SandiaDecay::Nuclide *u238 = is_background ? db->nuclide( "U238" ) : nullptr;
    const SandiaDecay::Nuclide *u235 = is_background ? db->nuclide( "U235" ) : nullptr;
    const SandiaDecay::Nuclide *th232 = is_background ? db->nuclide( "Th232" ) : nullptr;
    const SandiaDecay::Nuclide *ra226 = is_background ? db->nuclide( "Ra226" ) : nullptr;
    const SandiaDecay::Nuclide *k40 = is_background ? db->nuclide( "K40" ) : nullptr;

    for( const OtherRefLine &bl : otherRefLinesToShow )
    {
      if( std::get<0>( bl ) < lower_photon_energy )
        continue;

      ReferenceLineInfo::RefLine line;
      line.m_energy = std::get<0>( bl );
      line.m_decay_intensity = std::get<1>( bl );
      line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
      line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;

      if( !std::get<2>(bl).empty() )
        line.m_decaystr = std::get<2>(bl) + ", ";

      switch( std::get<3>( bl ) )
      {
        case OtherRefLineType::U238Series:
          line.m_parent_nuclide = u238;
          line.m_decaystr += "U238 series";
          break;

        case OtherRefLineType::U235Series:
          line.m_parent_nuclide = u235;
          line.m_decaystr += "U235 series";
          break;

        case OtherRefLineType::Th232Series:
          line.m_parent_nuclide = th232;
          line.m_decaystr += "Th232 series";
          break;

        case OtherRefLineType::Ra226Series:
          line.m_parent_nuclide = ra226;
          line.m_decaystr += "U238 (Ra226) series";
          break;

        case OtherRefLineType::K40Background:
          line.m_parent_nuclide = k40;
          line.m_decaystr += "Primordial";
          break;

        case OtherRefLineType::BackgroundXRay:
        {
          //std::get<2>( bl ) will be like "Pb xray"
          line.m_decaystr = std::get<2>( bl );
          if( !std::get<4>( bl ).empty() )
            line.m_decaystr += (line.m_decaystr.empty() ? "" : ", ") + std::get<4>( bl );
          line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Xray;
          vector<string> parts;
          SpecUtils::split( parts, std::get<2>( bl ), " " );
          if( !parts.empty() )
            line.m_element = db->element( parts[0] );
          break;
        }//case OtherRefLineType::BackgroundXRay:

        case OtherRefLineType::OtherBackground:
        case OtherRefLineType::BackgroundReaction:
        {
          line.m_decaystr = std::get<2>( bl );
          
          line.m_parent_nuclide = db->nuclide( std::get<2>( bl ) );
          // TODO: try to get reaction if didnt get nuclide - also, nuclide list may be CSV, could split that
          //const ReactionGamma *rctnDb = ReactionGammaServer::database();
          //if( rctnDb )
          //{
          //  rctnDb->gammas( input.m_input_txt, rctn_gammas );
            //...
          //}
          
          if( !std::get<4>( bl ).empty() )
            line.m_decaystr += (line.m_decaystr.empty() ? "" : ", ") + std::get<4>( bl );
          break;
        }
      }//switch( get<3>(*bl) )


      lines.push_back( line );
    }//for( otherRefLinesToShow )
  }//if( !otherRefLinesToShow.empty() )


  // Now calc detector response and shielding
  //  Up to now, we shouldnt have any escape or sum gammas in answer.m_ref_lines
  double max_alpha_br = 0.0, max_beta_br = 0.0, max_photon_br = 0.0;
  for( ReferenceLineInfo::RefLine &line : lines )
  {
    assert( (line.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::Normal)
        || (line.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::Annihilation) );

    switch( line.m_particle_type )
    {
      case ReferenceLineInfo::RefLine::Particle::Alpha:
        max_alpha_br = std::max( max_alpha_br, line.m_decay_intensity );
        break;

      case ReferenceLineInfo::RefLine::Particle::Beta:
        max_beta_br = std::max( max_beta_br, line.m_decay_intensity );
        break;

      case ReferenceLineInfo::RefLine::Particle::Gamma:
      case ReferenceLineInfo::RefLine::Particle::Xray:
      {
        double energy = line.m_energy;
        switch( line.m_source_type )
        {
          case ReferenceLineInfo::RefLine::RefGammaType::Normal:
          case ReferenceLineInfo::RefLine::RefGammaType::Annihilation:
            break;
            
          case ReferenceLineInfo::RefLine::RefGammaType::SingleEscape:
            // TODO: Need to put in S.E. DRF factor here
            energy += 510.998950;
            break;
            
          case ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape:
            // TODO: Need to put in D.E. DRF factor here
            energy += 2.0*510.998950;
            break;
            
          case ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak:
          case ReferenceLineInfo::RefLine::RefGammaType::SumGammaPeak:
            assert( 0 );
            break;
        }//switch( line.m_source_type )
        
        if( input.m_det_intrinsic_eff )
          line.m_drf_factor = input.m_det_intrinsic_eff( energy );
        
        if( input.m_shielding_att )
          line.m_shield_atten = input.m_shielding_att( energy );
        
        max_photon_br = std::max( max_photon_br,
                                 line.m_decay_intensity * line.m_drf_factor * line.m_shield_atten );
        break;
      }
    }//switch( line.m_particle_type )
  }//for( ReferenceLineInfo::RefLine &line : lines )
 

  const double alpha_sf = ((max_alpha_br > 0.0) && !IsNan(max_alpha_br)) ? (1.0 / max_alpha_br) : 1.0;
  const double beta_sf = ((max_beta_br > 0.0) && !IsNan(max_beta_br)) ? (1.0 / max_beta_br) : 1.0;
  const double photon_sf = ((max_photon_br > 0.0) && !IsNan(max_photon_br)) ? (1.0 / max_photon_br) : 1.0;
    
  for( ReferenceLineInfo::RefLine &line : lines )
  {
    switch( line.m_particle_type )
    {
      case ReferenceLineInfo::RefLine::Particle::Alpha: 
        line.m_particle_sf_applied = alpha_sf;
        line.m_normalized_intensity = line.m_decay_intensity * alpha_sf;
        break;

      case ReferenceLineInfo::RefLine::Particle::Beta:
        line.m_particle_sf_applied = beta_sf;
        line.m_normalized_intensity = line.m_decay_intensity * beta_sf;
        break;

      case ReferenceLineInfo::RefLine::Particle::Gamma:
      case ReferenceLineInfo::RefLine::Particle::Xray:
        line.m_particle_sf_applied = photon_sf;
        line.m_normalized_intensity = photon_sf * line.m_decay_intensity * line.m_drf_factor * line.m_shield_atten;
        break;
    }//switch( line.m_particle_type )

    // We wont filter out lines smaller than wanted here
    if( (line.m_transition || is_background)
      && (line.m_normalized_intensity <= input.m_lower_br_cutt_off
        || IsInf( line.m_normalized_intensity )
        || IsNan( line.m_normalized_intensity )) )
    {
      continue;
    }


    // Now lets fill out line.m_decaystr, line.m_particlestr, and line.m_elementstr
    if( line.m_decaystr.empty() )
    {
      if( line.m_transition )
      {
        if( line.m_transition->parent )
          line.m_decaystr = line.m_transition->parent->symbol;
        if( line.m_transition->child )
          line.m_decaystr += " to " + line.m_transition->child->symbol;

        // TODO: for alphas and betas its pretty rudundant to have this next line (I guess its redundant no matter what actually)
        line.m_decaystr += string( " via " ) + SandiaDecay::to_str( line.m_transition->mode );
      }//if( line.m_transition )

      if( line.m_reaction )
        line.m_decaystr = line.m_reaction->name();
    }//if( line.m_decaystr.empty() )


    if( line.m_particlestr.empty() )
    {
      switch( line.m_particle_type )
      {
        case ReferenceLineInfo::RefLine::Particle::Alpha:
          line.m_particlestr = "alpha";
          break;

        case ReferenceLineInfo::RefLine::Particle::Beta:
          line.m_particlestr = "beta";
          break;

        case ReferenceLineInfo::RefLine::Particle::Gamma:
          line.m_particlestr = "gamma";
          break;

        case ReferenceLineInfo::RefLine::Particle::Xray:
          line.m_particlestr = "xray";
          break;
      }//switch( line.m_particle_type )
    }//if( line.m_particlestr.empty() )

    if( line.m_elementstr.empty() )
    {
      const SandiaDecay::Element *element = el;
      if( !element && nuc && (line.m_particle_type == ReferenceLineInfo::RefLine::Particle::Xray) )
        element = db->element( nuc->atomicNumber );
     
      if( element )
        line.m_elementstr = element->name;
    }//if( line.m_elementstr.empty() )

    answer.m_ref_lines.push_back( line );
  }//for( ReferenceLineInfo::RefLine &line : answer.m_ref_lines )


  // If we add in escape peaks - we could put them in here

  // Add in coincident gammas
  if( !gamma_coincidences.empty() )
  {
    double max_coincidence_br = 0.0;
    vector<ReferenceLineInfo::RefLine> coinc_ref_lines;
    for( const auto &casc : gamma_coincidences )
    {
      const SandiaDecay::Transition *const &trans = get<0>( casc );
      const double &first_br = get<1>( casc );
      const float &first_energy = get<2>( casc );
      const float &second_energy = get<3>( casc );
      const float &coinc_frac = get<4>( casc );
      const double &second_br = get<5>( casc );

      const float energy = first_energy + second_energy;

      ReferenceLineInfo::RefLine line;
      line.m_energy = energy;
      line.m_decay_intensity = first_br * coinc_frac;
      line.m_parent_nuclide = nuc;
      line.m_transition = trans;
      line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
      line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak;

      if( input.m_det_intrinsic_eff )
        line.m_drf_factor = input.m_det_intrinsic_eff( first_energy ) * input.m_det_intrinsic_eff( second_energy );

      if( input.m_shielding_att )
        line.m_shield_atten = input.m_shielding_att( first_energy ) * input.m_shielding_att( second_energy );

      const double amp = line.m_decay_intensity * line.m_drf_factor * line.m_shield_atten;
      assert( !IsNan( amp ) && !IsInf( amp ) );
      if( IsNan( amp ) || IsInf( amp ) )
      {
        cerr << "Unexpected NaN or Inf coincidence amp." << endl;
        continue;
      }

      line.m_decaystr = "Cascade sum";
      if( trans && trans->parent )
        line.m_decaystr += " " + trans->parent->symbol;
      if( trans && trans->child )
        line.m_decaystr += " to " + trans->child->symbol;

      char buffer[128];
      snprintf( buffer, sizeof( buffer ),
        " (%.1f + %.1f keV, coinc=%.3g)",
        first_energy, second_energy, coinc_frac );

      line.m_decaystr += buffer;
      line.m_particlestr = "cascade-sum";

      coinc_ref_lines.push_back( std::move( line ) );

      max_coincidence_br = std::max( max_coincidence_br, amp );
    }//for( loop over cascades )

    assert( coinc_ref_lines.empty()
      || ((max_coincidence_br > 0.0) && !IsNan( max_coincidence_br )) );

    // Scale the coincidence line amplitudes to be between 0
    for( ReferenceLineInfo::RefLine &line : coinc_ref_lines )
    {
      const double sf = 1.0 / max_coincidence_br;
      line.m_particle_sf_applied = sf;
      const double amp = line.m_decay_intensity * line.m_drf_factor * line.m_shield_atten * sf;
      line.m_normalized_intensity = amp;
    }//for( ReferenceLineInfo::RefLine &line : coinc_ref_lines )

     // There can be tons of cascade sums (4834 for U238), we'll limit the number 
    //   we draw to an arbitrary 350, because this is even more than I expect to 
    //   be relevant (although I didnt actually check this).
    //  TODO: limit based on importance, and not a flat limit, e.g., use something like
    //        yield(i)*sqrt(energy(i))/sum(yield*sqrt(energy))
    const size_t max_cascade_sums = 350;
    if( coinc_ref_lines.size() > max_cascade_sums )
    {
      std::sort( begin( coinc_ref_lines ), end( coinc_ref_lines ),
        []( const ReferenceLineInfo::RefLine &lhs, const ReferenceLineInfo::RefLine &rhs ) -> bool {
          if( lhs.m_normalized_intensity == rhs.m_normalized_intensity )
            return lhs.m_energy > rhs.m_energy;
          return lhs.m_normalized_intensity > rhs.m_normalized_intensity;
        } );

      cout << "Resizing cascade sums from " << coinc_ref_lines.size() << " to " << max_cascade_sums << endl;
      coinc_ref_lines.resize( max_cascade_sums );
    }//if( coinc_ref_lines.size() > 350 )

    answer.m_ref_lines.reserve( answer.m_ref_lines.size() + coinc_ref_lines.size() );

    for( const ReferenceLineInfo::RefLine &line : coinc_ref_lines )
    {
      const double &amp = line.m_normalized_intensity;
      if( !IsNan( amp ) && !IsInf( amp ) 
          && (amp >= std::numeric_limits<float>::min()) // numeric_limits<float>::min()==1.17549e-38
          && (amp > input.m_lower_br_cutt_off)
      )
       answer.m_ref_lines.push_back( line );
    }//for( ReferenceLineInfo::RefLine &line : coinc_ref_lines )
  }//if( !gamma_coincidences.empty() )


  //Clientside javascript currently doesnt know about this garuntee that gamma
  //  lines will be sorted by energy.
  answer.sortByEnergy();

  return answer_ptr;
}//std::shared_ptr<ReferenceLineInfo> generateRefLineInfo()


std::shared_ptr<ReferenceLineInfo> ReferencePhotopeakDisplay::refLineForUserInput()
{
  RefLineInput user_input = userInput();
  shared_ptr<ReferenceLineInfo> ref_lines = refLineForUserInput();
  if( ref_lines )
  {
    //check if user label or age changed in ref_lines->m_input, and if so, alter GUI state
    //also, check on color - it may need changing, both in ref_lines, and the gui.
  }

  //Call this function from the old function and check values.
  
  //bool showGammaCB = true, showXrayCb = true, showAplhaCb = true, showBetaCb = true;
  //if( nuc ) // Show everything
  //{
  //} else if( el ) // For x-ray only show Show x-ray option  
  //  showGammaCB = showAplhaCb = showBetaCb = false;
  //else if( !rctnGammas.empty() ) // For reactions only show Show Gamma option
  //  showXrayCb = showAplhaCb = showBetaCb = false;
  //else if( isBackground )
  //  showAplhaCb = showBetaCb = false;
  //else // Show all options, otherwise whole area will be blank 
  //{
  //}

  //m_showXrays->setHidden( !showXrayCb );
  //m_showGammas->setHidden( !showGammaCB );
  //m_showAlphas->setHidden( !showAplhaCb );
  //m_showBetas->setHidden( !showBetaCb );


  //if( show )
  //{
  //  if( m_persistLines->isDisabled() )
  //    m_persistLines->enable();
  //  if( m_clearLines->isDisabled() )
  //    m_clearLines->enable();
  //  if( m_csvDownload && !m_csvDownload->isEnabled() )
  //    m_csvDownload->enable();

  //  if( m_spectrumViewer && m_spectrumViewer->isPhone() )
  //    m_clearLines->setText( m_persisted.empty() ? "Remove" : "Remove All" );
  //  else
  //    m_clearLines->setText( m_persisted.empty() ? "Clear" : "Clear All" );


  //updateOtherNucsDisplay();

  return ref_lines;
}//refLineForUserInput()


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
}// vector<DecayParticleModel::RowData> createTableRows( const ReferenceLineInfo &refLine );


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


void ReferencePhotopeakDisplay::updateDisplayChange()
{
  if( m_currently_updating )
    return;
  
  RefLineInput user_input = userInput();
  updateDisplayFromInput( user_input );
}


void ReferencePhotopeakDisplay::updateDisplayFromInput( RefLineInput user_input )
{
  UpdateGuard guard( m_currently_updating );
  
  shared_ptr<ReferenceLineInfo> ref_lines = generateRefLineInfo( user_input );
  
  if( ref_lines )
  {
    assert( (ref_lines->m_source_type == ReferenceLineInfo::SourceType::Nuclide) == (ref_lines->m_nuclide != nullptr) );
    assert( (ref_lines->m_source_type == ReferenceLineInfo::SourceType::FluorescenceXray) == (ref_lines->m_element != nullptr) );
    assert( (ref_lines->m_source_type == ReferenceLineInfo::SourceType::Reaction) == (!ref_lines->m_reactions.empty()) );
    //ReferenceLineInfo::SourceType::Background
    //ReferenceLineInfo::SourceType::CustomEnergy
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
  
  const bool hasPromptEquilib = (nuclide && nuclide->canObtainPromptEquilibrium());
  m_promptLinesOnly->setHidden( !hasPromptEquilib );
  m_promptLinesOnly->setChecked( ref_lines && ref_lines->m_input.m_promptLinesOnly );
  if( !hasPromptEquilib )
    m_promptLinesOnly->setUnChecked();
  
  const bool enable_aging = (nuclide && !ref_lines->m_input.m_promptLinesOnly && !nuclide->decaysToStableChildren());
  const string agestr = (!enable_aging || !ref_lines) ? string() : ref_lines->m_input.m_age;
  
  m_ageEdit->setText( WString::fromUTF8(agestr) );
  m_ageEdit->setEnabled( enable_aging );
  
  const string hlstr = !nuclide ? string()
  : ("&lambda;<sub>&frac12;</sub>="
     +  PhysicalUnits::printToBestTimeUnits( nuclide->halfLife, 2 ));
  m_halflife->setText( WString::fromUTF8(hlstr) );
  
  m_persistLines->setEnabled( show_lines );
  m_clearLines->setDisabled( m_persisted.empty() && !show_lines );
  
  if( m_csvDownload )
    m_csvDownload->setDisabled( !show_lines );
  
  
  const bool isPhone = ( m_spectrumViewer && m_spectrumViewer->isPhone() );
  const WString clearLineTxt = isPhone ? (m_persisted.empty() ? "Remove" : "Remove All")
  : ( m_persisted.empty() ? "Clear" : "Clear All" );
  if( clearLineTxt != m_clearLines->text() )
    m_clearLines->setText( clearLineTxt );
  
  
  bool showGammaCB = true, showXrayCb = true, showAplhaCb = true, showBetaCb = true;
  switch( src_type )
  {
    case ReferenceLineInfo::SourceType::Nuclide:
      // Show everything
      showGammaCB = showXrayCb = showAplhaCb = showBetaCb = true;
      break;
      
    case ReferenceLineInfo::SourceType::FluorescenceXray:
      // Only show x-ray option - but really we shouldnt show any of them
      showGammaCB = showAplhaCb = showBetaCb = false;
      break;
      
    case ReferenceLineInfo::SourceType::Reaction:
      // For reactions only show Show Gamma option
      showXrayCb = showAplhaCb = showBetaCb = false;
      break;
      
    case ReferenceLineInfo::SourceType::Background:
      showAplhaCb = showBetaCb = false;
      break;
      
    case ReferenceLineInfo::SourceType::CustomEnergy:
      // We treat custom energy as a gamma, so only show it - but really we shouldnt show any of them
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
  
  if( ref_lines && (ref_lines->m_validity == ReferenceLineInfo::InputValidity::Valid) )
  {
    m_showXrays->setChecked( ref_lines->m_input.m_showXrays );
    m_showGammas->setChecked( ref_lines->m_input.m_showGammas );
    m_showAlphas->setChecked( ref_lines->m_input.m_showAlphas );
    m_showBetas->setChecked( ref_lines->m_input.m_showBetas );
    m_showCascadeSums->setChecked( ref_lines->m_input.m_showCascades );
  }
  
  
  // Add current nuclide to lis of previous nuclides
  if( m_currentlyShowingNuclide.m_validity == ReferenceLineInfo::InputValidity::Valid )
  {
    OtherNuc prev;
    prev.m_input = m_currentlyShowingNuclide.m_input;
    prev.m_nuclide = m_currentlyShowingNuclide.m_input.m_input_txt;
    
    // Remove any other previous nuclides that have same `prev.m_nuclide` as what
    //  we are about to push on
    m_prevNucs.erase(std::remove_if(begin(m_prevNucs), end(m_prevNucs),
                                    [&prev](const OtherNuc &val) -> bool {
      return (val.m_input.m_input_txt == prev.m_input.m_input_txt);
    }), end(m_prevNucs));
    
    // Push new value onto front of history
    m_prevNucs.push_front( std::move(prev) );
    
    // Check history length, and truncate if needed
    if( m_prevNucs.size() > m_max_prev_nucs )
      m_prevNucs.resize( m_max_prev_nucs );
  }//if( !m_currentlyShowingNuclide.labelTxt.empty() )
  
  
  //Default to using the old color, unless we find a reason to change it.
  WColor color = ref_lines->m_input.m_color;
  
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
    if( !user_input.m_shielding_name.empty() && !user_input.m_shielding_thickness.empty() )
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
    }else if( !user_input.m_shielding_an.empty() && !user_input.m_shielding_ad.empty() )
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
        nonDefaultOpts |= !ref_lines->m_input.m_showXrays;
        nonDefaultOpts |= !ref_lines->m_input.m_showGammas;
        break;
        
      case ReferenceLineInfo::SourceType::CustomEnergy:
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
      m_cascadeWarn = new WText("x-rays are not included in cascades");
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
}//void updateDisplayChange()



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
    m_clearLines->setText( "Remove All" );
  else
    m_clearLines->setText( "Clear All" );
  
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
    value = doc->allocate_string( m_ageEdit->text().toUTF8().c_str() );
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
        
        shared_ptr<ReferenceLineInfo> ref_line = generateRefLineInfo( input );
        if( !ref_line || (ref_line->m_validity != ReferenceLineInfo::InputValidity::Valid) )
          throw runtime_error( "Couldn't generate reference lines from source '" + input.m_input_txt + "'" );
        
        m_currentlyShowingNuclide = *ref_line;
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
        
        shared_ptr<ReferenceLineInfo> ref_line = generateRefLineInfo( input );
        if( !ref_line || (ref_line->m_validity != ReferenceLineInfo::InputValidity::Valid) )
          throw runtime_error( "Couldn't generate persisted reference lines from source '" + input.m_input_txt + "'" );
        m_persisted.push_back( *ref_line );
      }//for( loop over DisplayedSource nodes )
    }//if( persisted_node )
    
    node = base_node->first_node( "Shielding", 9 );
    if( node )
      m_shieldingSelect->deSerialize( node );
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
      m_ageEdit->setText( PhysicalUnits::printToBestTimeUnits(age) );
  }else
  {
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
    m_clearLines->setText( "Clear" );
  else
    m_clearLines->setText( "Remove" );

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



