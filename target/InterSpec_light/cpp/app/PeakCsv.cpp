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

#include <cstdio>
#include <string>
#include <ostream>
#include <algorithm>

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"

#include "PeakCsv.h"

using namespace std;


namespace PeakCsv
{

void write( std::ostream &out, const std::deque<std::shared_ptr<const PeakDef>> &peaks,
            const std::shared_ptr<const SpecUtils::Measurement> &data )
{
  const string eol_char = "\r\n";

  out << "Centroid,  Net_Area,   Net_Area,      Peak, FWHM,   FWHM,Reduced, ROI_Total,ROI, File"
         ",         ,     ,     , Nuclide, Photopeak_Energy, ROI_Lower_Energy, ROI_Upper_Energy, Color, User_Label, Continuum_Type, "
         "Skew_Type, Continuum_Coefficients, Skew_Coefficients,         , Peak_Type" << eol_char
      << "     keV,    Counts,Uncertainty,       CPS,  keV,Percent,Chi_Sqr,    Counts,ID#, Name"
         ", LiveTime, Date, Time,        ,              keV,              keV,              keV, (css),           ,               , "
         "         ,                       ,                  , RealTime,          " << eol_char;

  const auto pad_left = []( string s, const size_t len ){ while( s.size() < len ) s = " " + s; return s; };
  const auto csv_escape = []( string &s ){
    if( s.empty() )
      return;
    SpecUtils::ireplace_all( s, "\"", "\"\"" );
    if( (s.find_first_of( ",\n\r\"\t" ) != string::npos) || (s.front() == ' ') || (s.back() == ' ') )
      s = "\"" + s + "\"";
  };

  for( size_t peakn = 0; peakn < peaks.size(); ++peakn )
  {
    const PeakDef &peak = *peaks[peakn];
    const double xlow = peak.lowerX(), xhigh = peak.upperX();

    float live_time = 1.0f, real_time = 0.0f;
    double region_area = 0.0;
    SpecUtils::time_point_t meastime{};
    if( data )
    {
      live_time = data->live_time();
      real_time = data->real_time();
      meastime = data->start_time();
      region_area = SpecUtils::gamma_integral( data, static_cast<float>(xlow), static_cast<float>(xhigh) );
    }

    string nuclide;
    char energy[32] = { '\0' };
    if( peak.hasSourceGammaAssigned() )
    {
      const PeakDef::Source &src = peak.source();
      snprintf( energy, sizeof(energy), "%.2f", peak.gammaParticleEnergy() );
      switch( src.kind )
      {
        case PeakDef::Source::Kind::None:     break;
        case PeakDef::Source::Kind::Nuclide:  nuclide = src.name; break;
        case PeakDef::Source::Kind::Reaction: nuclide = src.name; SpecUtils::ireplace_all( nuclide, ",", " " ); break;
        case PeakDef::Source::Kind::Xray:     nuclide = src.name + "-xray"; break;
      }

      if( (src.kind == PeakDef::Source::Kind::Nuclide) || (src.kind == PeakDef::Source::Kind::Reaction) )
      {
        switch( src.gamma_type )
        {
          case PeakDef::NormalGamma:
          case PeakDef::AnnihilationGamma: break;
          case PeakDef::XrayGamma:         nuclide += " (x-ray)"; break;
          case PeakDef::SingleEscapeGamma: nuclide += " (S.E.)";  break;
          case PeakDef::DoubleEscapeGamma: nuclide += " (D.E.)";  break;
        }
      }
    }//if( peak.hasSourceGammaAssigned() )

    char buffer[64];
    snprintf( buffer, sizeof(buffer), "%.2f", peak.mean() );
    const string meanstr = pad_left( buffer, 8 );

    snprintf( buffer, sizeof(buffer), "%.1f", peak.peakArea() );
    const string areastr = pad_left( buffer, 10 );

    snprintf( buffer, sizeof(buffer), "  %.1f", peak.peakAreaUncert() );
    string areauncertstr = buffer;
    while( areauncertstr.size() < 11 )
      areauncertstr += " ";

    string cpsstr;
    if( (live_time > 0.0f) && data )
    {
      snprintf( buffer, sizeof(buffer), "%1.4e", (peak.peakArea()/live_time) );
      cpsstr = buffer;
      const size_t epos = cpsstr.find( "e" );
      if( epos != string::npos )
      {
        if( (cpsstr[epos+1] != '-') && (cpsstr[epos+1] != '+') )
          cpsstr.insert( cpsstr.begin() + epos + 1, '+' );
        while( (cpsstr.size() - epos) < 5 )
          cpsstr.insert( cpsstr.begin() + epos + 2, '0' );
      }
    }//if( live time )

    const double width = peak.gausPeak() ? (2.35482*peak.sigma()) : 0.5*peak.roiWidth();
    snprintf( buffer, sizeof(buffer), "%.2f", width );
    const string widthstr = pad_left( buffer, 5 );

    snprintf( buffer, sizeof(buffer), "%.2f%%", (100.0*width/peak.mean()) );
    const string widthpercentstr = pad_left( buffer, 7 );

    snprintf( buffer, sizeof(buffer), "  %.2f", peak.chi2dof() );
    const string chi2str = pad_left( buffer, 7 );

    snprintf( buffer, sizeof(buffer), "  %.1f", region_area );
    const string roiareastr = pad_left( buffer, 10 );

    snprintf( buffer, sizeof(buffer), "  %i", static_cast<int>(peakn + 1) );
    const string numstr = pad_left( buffer, 3 );

    const string specfilename = "  ";  //InterSpec leaves this blank

    string live_time_str, real_time_str;
    if( live_time > 0.0f )
    {
      snprintf( buffer, sizeof(buffer), "%.3f", live_time );
      live_time_str = buffer;
    }
    if( real_time > 0.0f )
    {
      snprintf( buffer, sizeof(buffer), "%.3f", real_time );
      real_time_str = buffer;
    }

    string datestr, timestr;
    if( !SpecUtils::is_special( meastime ) )
    {
      const string tstr = SpecUtils::to_common_string( meastime, true );
      const size_t pos = tstr.find( ' ' );
      if( pos != string::npos )
      {
        datestr = tstr.substr( 0, pos );
        timestr = tstr.substr( pos + 1 );
      }
    }

    string color_str = peak.lineColor().cssText();
    string user_label = peak.userLabel();
    csv_escape( color_str );
    csv_escape( user_label );

    const shared_ptr<const PeakContinuum> continuum = peak.continuum();
    const PeakContinuum::OffsetType cont_type = continuum->type();
    // BiLinearStepCDF coefficients follow serialization version 3; InterSpec tags the type with this
    static_assert( PeakContinuum::sm_xmlSerializationVersion == 3, "Check the BiLinearStepCDF CSV tag" );
    const string continuum_type = string( PeakContinuum::offset_type_str( cont_type ) )
                                  + ((cont_type == PeakContinuum::BiLinearStepCDF) ? "(v3)" : "");

    string cont_coefs;
    if( (cont_type != PeakContinuum::NoOffset) && (cont_type != PeakContinuum::External) )
    {
      cont_coefs = SpecUtils::printCompact( continuum->referenceEnergy(), 6 );
      const vector<double> &pars = continuum->parameters();
      const size_t npar = std::min( pars.size(), PeakContinuum::num_parameters( cont_type ) );
      for( size_t i = 0; i < npar; ++i )
        cont_coefs += " " + SpecUtils::printCompact( pars[i], 7 );
    }

    string skew_coefs;
    const size_t num_skew = PeakDef::num_skew_parameters( peak.skewType() );
    for( size_t i = 0; i < num_skew; ++i )
    {
      const auto par = PeakDef::CoefficientType( PeakDef::SkewPar0 + i );
      skew_coefs += (skew_coefs.empty() ? "" : " ") + SpecUtils::printCompact( peak.coefficient( par ), 7 );
    }

    out << meanstr << ',' << areastr << ',' << areauncertstr << ',' << cpsstr << ',' << widthstr
        << ',' << widthpercentstr << ',' << chi2str << ',' << roiareastr << ',' << numstr
        << ',' << specfilename << ',' << live_time_str << ',' << datestr << ',' << timestr
        << ',' << nuclide << ',' << energy << ',' << xlow << ',' << xhigh << ',' << color_str
        << ',' << user_label << ',' << continuum_type << ',' << PeakDef::to_string( peak.skewType() )
        << ',' << cont_coefs << ',' << skew_coefs << ',' << real_time_str
        << ',' << PeakDef::to_str( peak.type() ) << eol_char;
  }//for( loop over peaks )
}//void write(...)

}//namespace PeakCsv
