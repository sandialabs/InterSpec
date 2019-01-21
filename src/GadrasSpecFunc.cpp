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

#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <exception>
#include <algorithm>

#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/GadrasSpecFunc.h"
#include "SpecUtils/SpectrumDataStructs.h"
#include "InterSpec/MassAttenuationTool.h"


using namespace std;

namespace
{
  /**
 Compute the fraction of transmitted gamma rays given the gamma-ray energy,
 atomic number and areal density of the attenuator.
 Account for hydrogen, which exhibits twice the scatter cross-section per mass.
 
 \param E energy of gamma ray in keV
 \param AN atomic number of shield
 \param AD total areal density (g/cm^2) of shield
 \param AHIN fractional areal density of hydrogen in shield
 
 \returns Fraction of gammas that will be transmitted through shielding
 */
float TransmissionH( const float E, const float AN, const float AD, const float AHIN )
{
  if( AN <= 0.5f || AD <= 0.0 )
    return 1.0f; //transmission is approx unity for low shielding amounts
  
  if( E < 1.0f )
    return 0.0f; //transmission is approx zero for energies < 1 keV
  
  //excess hydrogen AD above 1 g/cm^2 ?
  const float FR = max( 0.0f, min( 1.0f, AHIN*AD-1.0f) );
  float AH = FR * AHIN;
  const float gcm2 = static_cast<float>( PhysicalUnits::g / PhysicalUnits::cm2 );
  
  //only apply at high hydrogen content and shielding
  if( AH > 0.0f )
  {
    // FUDGE: AH less important at high E
    if( E > 400.0f )
      AH *= exp( -((E-400.0f)/1000.0f) );
    
    const float ANother = (AN*AD - AH*AD)/((1.0f-AH)*AD);
    const float mu = MassAttenuation::massAttenuationCoeficientFracAN( ANother, E );
    const float mu_hydrogen = MassAttenuation::massAttenuationCoeficient( 1, E );
    
    return exp(-(1.0f-AH)*AD*gcm2*mu) * exp(-AH*AD*gcm2*mu_hydrogen);
  }//if( AH > 0.0f )
  
  
  const float mu = MassAttenuation::massAttenuationCoeficientFracAN( AN, E );
  return exp( -max( 0.0f, AD*gcm2*mu ) );
}//float TransmissionH( E,AN,AD,AHIN )

}//namespace


GadrasScatterTable::GadrasScatterTable( const std::string &datafile )
{
  ifstream input( datafile.c_str(), ios::in | ios::binary );
  
  if( !input.read( (char *)m_data, sizeof(m_data) ) )
    throw std::runtime_error( "Error reading " + datafile );
  
/*
 //Above read is equivalent to:
  for( int l = 0; l < sm_num_areal_density; ++l )
    for( int k = 0; k < sm_num_atomic_number; ++k )
      for( int j = 0; j < sm_num_input_energy; ++j )
        for( int i = 0; i < sm_num_out_energy; ++i )
          if( !input.read( (char *)&(m_data[l][k][j][i]), 4 ) )
            throw std::runtime_error( "Error reading " + datafile );
*/
  
  //Now seek to the end and see how much is left - it should be zero
  const std::streampos ending_pos = input.tellg();
  input.seekg( 0, ios::end );
  const std::streampos end_loc = input.tellg();
  
  //Should maybe throw and exception or something if the following fails...
  if( ending_pos != end_loc )
    cout << (end_loc - ending_pos) << " bytes left in " << datafile
         << ", there may be something wrong!!" << endl;
}//GadrasScatterTable( constructor )


GadrasScatterTable::~GadrasScatterTable()
{
  //nothing to do here
}


float GadrasScatterTable::getContinuum( std::vector<float> &answer,
                                       const float sourceEnergy,
                                       const float sourceIntensity,
                                       const float atomicNumber,
                                       const float arealDensity,
                                       const float fracHydrogen,
                                       const std::vector<float> &output_binning ) const
{
  const float energyMesh[sm_num_input_energy] = { 60.0f, 87.0f, 89.0f, 115.0f, 116.0f, 121.0f, 122.0f, 200.0f, 300.0f, 400.0f, 600.0f, 1000.0f, 1600.0f, 2600.0f, 4000.0f, 9000.0f };
  const float anMesh[sm_num_atomic_number] = { 5.28f, 13.0f, 26.0f, 82.0f, 92.0f, 94.0f }; // J (PE,Al,Fe,Pb,U,Pu)
  const float adMesh[sm_num_areal_density] = { 1.0f, 2.0f, 4.0f, 8.0f, 16.0f, 32.0f, 64.0f, 128.0f, 256.0f };

  answer.resize( output_binning.size() );
  for( size_t i = 0; i < output_binning.size(); ++i )
    answer[i] = 0.0f;
  
  
//Greg T. C# translation used the following to find i
//  int i = 6;
//  for( i = 6; i && (energyMesh[i] > sourceEnergy); --i ){}

  int i;
  if( sourceEnergy < 88.04f ) // Pb K-edge
		i = 0;
  else if( sourceEnergy < 115.5f ) // U K-edge
		i = 2;
  else if( sourceEnergy < 121.8f ) // Pu K-edge
		i = 4;
  else
  {
    i = sm_num_input_energy - 2;
    
    while( energyMesh[i] > sourceEnergy )
    {
      i -= 1;
      if( i == 6 )
        break;
    }
  }
  
  const int ip = i + 1;
  const float fi = std::max( 0.0f, (sourceEnergy-energyMesh[i])/(energyMesh[ip]-energyMesh[i]) );
  
  
//Greg T. C# translation used the following to find j
//  int j = 3;
//  for( j = 3; j && (anMesh[j] > atomicNumber); --j ){}

  int j = 4;
  while( anMesh[j] > atomicNumber )
  {
    j -= 1;
    if( j == 0 )
      break;
  }
  
  int jp = j + 1;
  
  
  if( j > 1 )//             ! get materials from emitted x-rays
  {
    //For meaning of iU, iPb, and iPu see MaterialsPresent in SpecFunc.f.
    //  I think they are supposed to represent if these materials are present
    //  in the source (or shielding) and if so
    const int iU = 0;
    const int iPb = 0;
    const int iPu = 0;
    
    if( j==4 && iU==0 ) j = 3;
    if( j==3 && iPb==0 ) j = 2;
    if( jp==3 && (atomicNumber > 82.0f || iU==1 || iPu==1) ) jp = 4;
    if( jp==4 && (atomicNumber > 92 || iPu==1) ) jp = 5;
  }//if( j > 3 )
  
  
  const float fj = (atomicNumber-anMesh[j])/(anMesh[jp]-anMesh[j]);
  
  int k = 7;
  while( (k > -1) && (adMesh[k] > arealDensity) )
    --k;

  int kp = k + 1;
  
  float fk;
  if( k > -1 )
    fk = (arealDensity - adMesh[k]) / (adMesh[kp]-adMesh[k]);
  else
    fk = arealDensity / adMesh[kp];

  const float ri = 1.0f - fi;
  const float rj = 1.0f - fj;
  const float rk = 1.0f - fk;
  
  vector<float> f( 32 );
  // get the flux distribution by matrix interpolation.
  if( k == -1 )
  {
    //No shielding
    for (int l = 0; l < 32; ++l)
    {
      f[l] = fk*(ri*(rj*m_data[kp][j][i][l] + fj*m_data[kp][jp][i][l]) +
                 fi*(rj*m_data[kp][j][ip][l] + fj*m_data[kp][jp][ip][l]));
      f[l] = max( 0.0f, f[l] );
    }
  }else
  {
    for( int l = 0; l < 32; ++l )
    {
      f[l] = ri*(rj*(rk*m_data[k][j][i][l] + fk*m_data[kp][j][i][l]) +
                 fj*(rk*m_data[k][jp][i][l] + fk*m_data[kp][jp][i][l]))
            + fi*(rj*(rk*m_data[k][j][ip][l] + fk*m_data[kp][j][ip][l]) +
               fj*(rk*m_data[k][jp][ip][l] + fk*m_data[kp][jp][ip][l]));
      f[l] = std::max( 0.0f, f[l] );
    }
  }

  
  //Next section of code following has not been checked.
  //C	Add continuum from PE if hydrogen fraction exceeds interpolation.
  float frac_pe = 0.0;
  const float hydrogen_weight_frac_pe = 0.144f;
  
  if( j == 0 )
    frac_pe = (fracHydrogen / hydrogen_weight_frac_pe) - rj;
  else
    frac_pe = fracHydrogen / hydrogen_weight_frac_pe;

  if( frac_pe > 0.0f )
  {
    cerr << "\nWarning: Correction for fraction of Hydrogen has not been checked\n";
    
    frac_pe = min( 1.0f, frac_pe );
    j = 0;
		
    if( k == -1 )
    {
      for( int l = 1; l < 32; ++l )
      {
        float fl = fk*(ri*m_data[0][j][i][l] + fi*m_data[0][j][ip][l]);
        fl = max( 0.0f, fl );
        f[l] = (1.0f - frac_pe)*f[l]+frac_pe*fl;
      }//for( int l = 1; l < 32; ++l )
    }else
    {
      for( int l = 1; l < 32; ++l )
      {
        float fl = ri*(rk*m_data[k][j][i][l]+fk*m_data[kp][j][i][l])
        	         + fi*(rk*m_data[k][j][ip][l]+fk*m_data[kp][j][ip][l]);
        fl = max(0.0f,fl);
        f[l] = (1.0f - frac_pe)*f[l] + frac_pe*fl;
      }
    }//if( k == -1 ) / else
  }//if( frac_pe > 0.0f )

  
  // rebin the flux to the proper energy group structure.
  const float sc = sourceEnergy / 31.0f;
//  double edge = sourceEnergy - sourceEnergy*(sourceEnergy/255.5)/(1 + sourceEnergy/255.5);
  
  vector<float> orig_binning(32);
  for (int l = 0; l < 32; ++l)
    orig_binning[l] = sc*l;
  
  //Section of code following next comment has not been checked since it might
  //  not be applicable.
  //	Scale scatter flux if less than 60 keV (minimum for transport).
  if( sourceEnergy < 60.0f )
  {
    cerr << "\nWarning: Continuum calculation for source energies bellow 60 keV not checked yet!\n";
    
    const float SC60 = 60.0f / 31.0f;
    for( int i = 0; i < 31; ++i )
    {
      float e_bar = SC60 * (i + 0.5f);
      const float trans60 = TransmissionH( e_bar, atomicNumber, arealDensity, fracHydrogen );
      if( trans60 < 1.0e-15)
      {
//        cout << "trans60=" << trans60 << ", for i=" << i << ", e_bar=" << e_bar << endl;
        f[i] = 0.0f;
      }else
      {
        e_bar = (sc*i + sc*(i+1)) / 2.0;
        const float transE = TransmissionH( e_bar, atomicNumber, arealDensity, fracHydrogen );
        f[i] *= transE / trans60;
//        cout << "e_bar=" << e_bar << ", transE=" << transE << ", trans60=" << trans60 << " f[" << i << "]=" << f[i] << endl;
      }
    }//for( int i = 0; i < 31; ++i )
  }//if( sourceEnergy < 60.0f )
  
  //Greg T. C# translation had the following line...
  //  f[30] = 2.0f*f[29] - f[28]; // what?

  // Project flat continuum to zero scatter angle (FUDGE).
  const float r28 = f[28]/(orig_binning[29] - orig_binning[28]);  //? equiv to: f[28]/SC
  f[29] = max( f[29], r28 * (orig_binning[30] - orig_binning[29]) );  //? equiv to: max( f[29], f[28] )
  f[30] = max( f[30], r28 * (orig_binning[31] - orig_binning[30]) );  //? equiv to: max( f[30], f[28] )
  
  
  //Section of code following next comment not implemented
  //C	Limit total scattering.
  const float trans = TransmissionH( sourceEnergy, atomicNumber, arealDensity, fracHydrogen );

  float sumf = 0.0f;
  for( size_t i = 0; i < f.size(); ++i)
  {
    f[i] *= sourceIntensity*trans;
    sumf += f[i];
  }
  
  /*
  const float ctot = MassAttenuation::massAttenuationCoeficientFracAN( atomicNumber, sourceEnergy );
  ScatToTotal=CKN/ctot;
  IF (AH1.GT.0) THEN
		CTOT=AttCoef(Energy,1.)
		ScatToTotal=(1.-AH1)*ScatToTotal+AH1*CKN/CTOT
  END IF
  Smax=3.7E10*(1.-Trans)*ScatToTotal
  IF (SumF.GT.Smax) THEN
		SC=Smax/SumF
		F(1:31)=SC*F(1:31)
  END IF
  */
  
  //Section of code following next comment not implemented
  //C	Attenuate the continuum for HPGe with thin dead layer.
  //
  
  rebin_by_lower_edge( orig_binning, f, output_binning, answer );
  
  for( size_t i = 0; i < answer.size(); ++i)
    answer[i] /= 3.7E10f;
  
  return trans * sourceIntensity;
}//void getContinuum(...)
