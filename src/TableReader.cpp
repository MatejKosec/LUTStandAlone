//============================================================================
// Name        : TableReader.cpp
// Author      : Matej Kosec
// Version     :
// Copyright   : 
// Description :
//============================================================================

#include <string>
#include <cmath>
#include <cassert>
#include "LUT.hpp"


using namespace std;

int main() {

	char* tablefile= "/home/matej/SU2/DEV/cfxlut/Code/CO2.rgp";
	CLookUpTable LUT2;
	LUT2 = CLookUpTable(tablefile);
	LUT2.SetTDState_rhoe(90, 0.535e+06);
	LUT2.SetTDState_PT(100000, 260);
	return 0;
	}

