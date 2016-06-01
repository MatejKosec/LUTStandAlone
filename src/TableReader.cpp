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
	LUT2.SetTDState_rhoe(66,-2.89275e+08);
	return 0;
	}



