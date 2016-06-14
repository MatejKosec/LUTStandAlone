//============================================================================
// Name        : TableReader.cpp
// Author      : Matej Kosec
// Version     :
// Copyright   : 
// Description :
//============================================================================

#include <string>
#include <stdlib.h>
#include <cmath>
#include <memory>
#include <cassert>
#include "LUT.hpp"


using namespace std;

int main() {

	char* tablefile= "CO2.rgp";
	CLookUpTable LUT2;
	LUT2 = CLookUpTable(tablefile);
	//LUT2.SetTDState_rhoe(90, 0.545e+06);
	//LUT2.SetTDState_Prho(5521052, 80);
	//LUT2.SetTDState_PT(5420052, 260);
	//LUT2.SetTDState_rhoT(90, 80);
	//LUT2.SetTDState_Ps(5400000, 2200);
	LUT2.SetTDState_hs(771000, 2500);

	return 0;
	}

