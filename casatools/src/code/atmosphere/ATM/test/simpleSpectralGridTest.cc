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
using namespace atm;



int main()

{

  SpectralGrid *pSpectralGrid;

  int numChan(11);

  int refChan(11/2+1);

  Frequency refFreq(100.,"Hz");

  Frequency chanSep(1.0,"Hz");

  pSpectralGrid = new SpectralGrid(numChan,refChan,refFreq,chanSep);

 

  cout << "Setting SPW:" << endl;

  cout << "- numChan = " << numChan << endl;

  cout << "- refFreq = " << refFreq.get("Hz") << " [Hz]" << endl;

  cout << "- refChan = " << refChan << endl;

  cout << "- chanSep = " << chanSep.get("Hz") << " [Hz]" << endl;

 

  cout << "\nSpectral Grid defined:" << endl;

  cout << "- frequency of refChan = "

       << pSpectralGrid->getChanFreq(refChan).get("Hz") << " [Hz]" << endl;

  for(int i = 0; i < numChan; i++) {

    cout << "Frequency of channel " << i << ": " << pSpectralGrid->getChanFreq(i).get("Hz") << " [Hz]" << endl;

  }
  
  delete pSpectralGrid;

  
}
