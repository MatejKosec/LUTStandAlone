//============================================================================
// Name        : TableReader.cpp
// Author      : Matej Kosec
// Version     :
// Copyright   : 
// Description :
//============================================================================

#include <string>
#include <iostream>
#include <stdlib.h>
#include "stdio.h"
#include <cmath>
#include <memory>
#include <cassert>
#include "LUT.hpp"


using namespace std;

int main() {

	char* tablefile= "CO2.rgp";
	CLookUpTable LUT2 = CLookUpTable(tablefile);
	//LUT2.reset_Restart();
	//LUT2.SetTDState_rhoe(90, 0.545e+06);
	//LUT2.reset_Restart();
	cout<<"rhoe \n";
	//LUT2.SetTDState_PT(5520052, 270);
	//LUT2.reset_Restart();
	cout<<"PT"<<endl;
	//LUT2.SetTDState_Prho(5521052, 80);
	cout<<"Prho \n";
	LUT2.SetTDState_rhoT(90, 80);
	LUT2.reset_Restart();
	cout<<"rhoT"<<endl;
	//LUT2.SetTDState_Ps(5400000, 2200);
	//LUT2.reset_Restart();
	cout<<"PS"<<endl;
	//LUT2.SetTDState_hs(531782, 2092);
	LUT2.reset_Restart();
	cout<<"HS"<<endl;

	return 0;
	}

