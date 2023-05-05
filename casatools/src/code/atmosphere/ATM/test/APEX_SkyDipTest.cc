/*******************************************************************************
 * ALMA - Atacama Large Millimeter Array
 * (c) Instituto de Estructura de la Materia, 2011
 * (in the framework of the ALMA collaboration).
 * All rights reserved.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 *******************************************************************************/

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;



#include <atmosphere/ATM/ATMPercent.h>
#include <atmosphere/ATM/ATMPressure.h>
#include <atmosphere/ATM/ATMNumberDensity.h>
#include <atmosphere/ATM/ATMMassDensity.h>
#include <atmosphere/ATM/ATMTemperature.h>
#include <atmosphere/ATM/ATMLength.h>
#include <atmosphere/ATM/ATMInverseLength.h>
#include <atmosphere/ATM/ATMOpacity.h>  
#include <atmosphere/ATM/ATMAngle.h>
#include <atmosphere/ATM/ATMHumidity.h>
#include <atmosphere/ATM/ATMFrequency.h>
#include <atmosphere/ATM/ATMWaterVaporRadiometer.h>
#include <atmosphere/ATM/ATMWVRMeasurement.h>
#include <atmosphere/ATM/ATMProfile.h>
#include <atmosphere/ATM/ATMSpectralGrid.h>
#include <atmosphere/ATM/ATMRefractiveIndex.h>
#include <atmosphere/ATM/ATMRefractiveIndexProfile.h>
#include <atmosphere/ATM/ATMSkyStatus.h>
#define STRLEN  200        // Max length of a row in a tpoint file
using namespace atm;

int main()
{  
  cout << " ApexTest:" << endl;
  cout << " ApexTest: THIS PROPOSED TEST OF THE ATM INTERFACE SOFTWARE IS BASED ON APEX SKY DIP DATA PROVIDED BY HEIKO" << endl;
  cout << " " << endl;
  
  cout << " ApexTest: STEP 1: CREATES REFERENCE ATMOSPHERIC PROFILE CORRESTONDING TO THE FOLLOWING BASIC PARAMETERS:" << endl;
  
  unsigned int atmType = 1; // TROPICAL
  Temperature      T( 270.114,"K");     // Ground temperature
  Pressure         P( 553.8,"mb");     // Ground Pressure
  Humidity         H(   9.06,"%" );     // Ground Relative Humidity (indication)
  Length         Alt(  5105,"m" );     // Altitude of the site 
  Length         WVL(   2.0,"km");     // Water vapor scale height
  double         TLR=  -6.5      ;     // Tropospheric lapse rate (must be in K/km)
  Length      topAtm(  48.0,"km");     // Upper atm. boundary for calculations
  Pressure     Pstep(  5.0,"mb");     // Primary pressure step (5.0 mb)
  double   PstepFact=         1.1;     // Pressure step ratio between two consecutive layers
  
  AtmProfile apex_AtmProfile( Alt, P, T, TLR, H, WVL, Pstep, PstepFact,  topAtm, atmType );
  
  cout << " ApexTest: First guess precipitable water vapor content: " << apex_AtmProfile.getGroundWH2O().get("mm") << "mm" << endl;
  cout << " ApexTest:  " << endl;

  cout << " ApexTest: STEP 2: CREATES SpectralGrid and RefractiveIndexProfile objects" << endl;

  unsigned int numchan=25;  unsigned int refchan=13;
  Frequency astrofreq(691.51570990969995,"GHz"); Frequency chansep(  0.04,"GHz"); Frequency intfreq(  6.0,"GHz");
  SidebandSide sidebandside=LSB; SidebandType sidebandtype=DSB;
  SpectralGrid apex_SpectralGrid(numchan, refchan, astrofreq, chansep, intfreq, sidebandside, sidebandtype); 
 
  // Creates RefractiveIndexProfile for the current atm profile (apex_AtmProfile) and spectral grid 
  // (apex_SpectralGrid, so far with only the first WVR channel). Later on we will add new spectral 
  // windows directly at the level of this RefractiveIndexProfile object.
  RefractiveIndexProfile apex_RefractiveIndexProfile(apex_SpectralGrid, apex_AtmProfile);


    
return 0;

}
