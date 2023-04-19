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
#include <stdlib.h>


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

// using namespace telcal::atm;
using namespace atm;
using namespace std;

int main(int argc, char *argv[]) {
  
  if (argc<3) {
    cout<<"ATMDiscontinuityTest.cpp fminGhz fmaxGhz"<<endl;
    exit(0);
  }

  // Read arguments
  double freqmin = atof(argv[1]);
  double freqmax = atof(argv[2]);
  
  if (freqmin<50. || freqmax>1000. || freqmax<50. || freqmax>1000.) {
    cout<<"Error : freq value should be in the interval [50Ghz,1000Ghz]"<<endl;
    exit(0);
  }
   
  // Define spectralGrid
  unsigned int numchan=1000; 
  unsigned int refchan=500;  

  double freq = (freqmin+freqmax)/2.; 
  Frequency reffreq(freq,"GHz");

  double dfreq = (freqmax-freqmin)/double(numchan); 

  Frequency chansep(  dfreq,"GHz");
  //cout<<freqmin<<" "<<freqmax<<" "<<freq<<" chansep="<<dfreq<<endl;

  SpectralGrid band(numchan, refchan, reffreq, chansep);


  // Atmospheric profile
  unsigned int atmType = 1; // TROPICAL
  Temperature      T( 273.0,"K" );     // Ground temperature
  Pressure         P( 550.0,"mb");     // Ground Pressure
  Humidity         H(  20.0,"%" );     // Ground Relative Humidity (indication)
  Length         Alt(  5000,"m" );     // Altitude of the site 
  Length         WVL(   1.0,"km");     // Water vapor scale height
  double         TLR=  -6.5      ;     // Tropospheric lapse rate (must be in K/km)
  Length      topAtm(  48.0,"km");     // Upper atm. boundary for calculations
  Pressure     Pstep(   5.0,"mb");     // Primary pressure step (10.0 mb)
  double   PstepFact=         1.1;     // Pressure step ratio between two consecutive layers

  AtmProfile myProfile( Alt, P, T, TLR, H, WVL, Pstep, PstepFact,  topAtm, atmType );

  RefractiveIndexProfile abs_band(band, myProfile); 
  SkyStatus mySky_1band(abs_band);
  mySky_1band.setAirMass(1.0);

  FILE * fp = fopen("tsky.dat","w");
  for(unsigned int i=0; i<mySky_1band.getNumChan(0); i++){
    fprintf(fp, "%f %f\n",
	    mySky_1band.getChanFreq(i).get("GHz"),
	    mySky_1band.getTebbSky(i).get("K"));
  }
  fclose(fp);

  return 0;
}

