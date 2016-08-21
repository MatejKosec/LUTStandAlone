#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cmath>
#include <cassert>
#include <algorithm>
#include "stdlib.h"
#include "stdio.h"
#include "LUT.hpp"
#include <iomanip>
using namespace std;

CLookUpTable::CLookUpTable(string Filename) {
	LUT_Debug_Mode = false;
	rank = 12201;

#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	Pressure_Reference_Value = 1;
	Temperature_Reference_Value = 1;
	Density_Reference_Value = 1;
	Velocity_Reference_Value = 1;
	Energy_Reference_Value = 1;

	//Detect LuT filetype
	if ((Filename).find(".tec") != string::npos) {
		if (rank == 12201) {
			cout << ".tec type LUT found" << endl;
		}
		LookUpTable_Load_TEC(Filename);
	} else {
		if (rank == 12201) {
			cout << "No recognized LUT format found, exiting!" << endl;
		}
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			Interpolation_Coeff[i][j] = -1.0;
		}
	}
	if (rank == 12201) {
		// Give the user some information on the size of the table
		cout << "Number of stations  in zone 0: " << nTable_Zone_Stations[0]
				<< endl;
		cout << "Number of triangles in zone 0: " << nTable_Zone_Triangles[0]
				<< endl;
		cout << "Number of stations  in zone 1: " << nTable_Zone_Stations[1]
				<< endl;
		cout << "Number of triangles in zone 1: " << nTable_Zone_Triangles[1]
				<< endl;

		cout << "Print LUT errors? (LUT_Debug_Mode):  " << LUT_Debug_Mode << endl;
		// Building an KD_tree for the HS thermopair
		cout << "Building trapezoidal map for rhoe..." << endl;
		cout << "Building trapezoidal map for Prho..." << endl;
		cout << "Building trapezoidal map for hs..." << endl;
		cout << "Building trapezoidal map for Ps..." << endl;
		cout << "Building trapezoidal map for rhoT..." << endl;
		//cout << "Building trapezoidal map for PT..."<<endl;

	}
}

CLookUpTable::~CLookUpTable(void) {
	// Delete the lut
	for (int i = 0; i < 2; i++) {
		delete[] ThermoTables_StaticEnergy[i];
		delete[] ThermoTables_Entropy[i];
		delete[] ThermoTables_Enthalpy[i];
		delete[] ThermoTables_Density[i];
		delete[] ThermoTables_Pressure[i];
		delete[] ThermoTables_SoundSpeed2[i];
		delete[] ThermoTables_Temperature[i];
		delete[] ThermoTables_dPdrho_e[i];
		delete[] ThermoTables_dPde_rho[i];
		delete[] ThermoTables_dTdrho_e[i];
		delete[] ThermoTables_dTde_rho[i];
		delete[] ThermoTables_Cp[i];
		delete[] ThermoTables_Mu[i];
		delete[] ThermoTables_dmudrho_T[i];
		delete[] ThermoTables_dmudT_rho[i];
		delete[] ThermoTables_Kt[i];
		delete[] ThermoTables_dktdrho_T[i];
		delete[] ThermoTables_dktdT_rho[i];
		//Delete the triangulations
		for (int j = 0; j < nTable_Zone_Triangles[i]; j++) {
			delete[] Table_Zone_Triangles[i][j];
		}
		delete[] Table_Zone_Triangles[i];
	}
}

void CLookUpTable::Search_NonEquispaced_Rho_Index(su2double rho) {
	{
		su2double grad, x00, y00;
		//  Determine the I index with binary search (rho is not assumed equispaced)
		while (UpperI - LowerI > 1) {
			middleI = (UpperI + LowerI) / 2;
			x00 = ThermoTables_Density[middleI][LowerJ];
			grad = ThermoTables_Density[middleI + 1][LowerJ] - x00;
			if (x00 * grad > rho * grad) {
				UpperI = middleI;
			} else if (x00 < rho) {
				LowerI = middleI;
			} else if (x00 == rho) {
				LowerI = middleI;
				UpperI = LowerI + 1;
				break;
			}
		}
	}
}
void CLookUpTable::Search_NonEquispaced_P_Index(su2double P) {
	su2double grad, x00, y00, y01, x01, RunVal;
	//Determine the J index using a binary search, and not assuming P is equispaced
	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2;
		x00 = ThermoTables_Pressure[LowerI][middleJ];
		grad = ThermoTables_Pressure[LowerI][middleJ + 1] - x00;
		if (x00 * grad > P * grad) {
			UpperJ = middleJ;
		} else if (x00 < P) {
			LowerJ = middleJ;
		} else if (x00 == P) {
			LowerJ = middleJ;
			UpperJ = LowerJ + 1;
			break;
		}
	}
}

void CLookUpTable::Search_i_for_X_given_j(su2double x, su2double y,
su2double **ThermoTables_X, su2double **ThermoTables_Y) {
	su2double grad, x00, y00, y01, x01, RunVal;
	while (UpperI - LowerI > 1) {
		middleI = (UpperI + LowerI) / 2;
		//Use interpolated T as the running variable for the search (RunVal)
		y00 = ThermoTables_Y[middleI][LowerJ];
		y01 = ThermoTables_Y[middleI][UpperJ];
		x00 = ThermoTables_X[middleI][LowerJ];
		x01 = ThermoTables_X[middleI][UpperJ];
		grad = ThermoTables_X[UpperI][LowerJ] - x00;
		RunVal = x00 + (x01 - x00) / (y01 - y00) * (y - y00);
		if (RunVal * grad > x * grad) {
			UpperI = middleI;
		} else if (RunVal * grad < x * grad) {
			LowerI = middleI;
		} else if (RunVal == x) {
			LowerI = middleI;
			UpperI = LowerI + 1;
			break;
		}
	}
}

void CLookUpTable::Search_j_for_Y_given_i(su2double x, su2double y,
su2double **ThermoTables_X, su2double **ThermoTables_Y) {
	su2double RunVal;
	su2double grad, x00, y00, y10, x10;
	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2; /*!< \brief Splitting index for the search */

		// The variable names is composed of a (i,j) pair
		y00 = ThermoTables_Y[LowerI][middleJ];
		y10 = ThermoTables_Y[UpperI][middleJ];
		x00 = ThermoTables_X[LowerI][middleJ];
		x10 = ThermoTables_X[UpperI][middleJ];
		//The search variable in j should be interpolated in i as well
		RunVal = y00 + (y10 - y00) / (x10 - x00) * (x - x00);
		grad = ThermoTables_Y[LowerI][middleJ + 1] - y00;
		if (RunVal * grad > y * grad) {
			UpperJ = middleJ;
		} else if (RunVal * grad < y * grad) {
			LowerJ = middleJ;
		} else if (RunVal == y) {
			LowerJ = middleJ;
			UpperJ = LowerJ + 1;
			break;
		}
	}
}

void CLookUpTable::Search_Linear_Skewed_Table(su2double x, su2double P,
su2double **ThermoTables_X) {
	su2double RunVal, rho;
	su2double grad, x00, y00, y10, x10;
	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2; /*!< \brief Splitting index for the search */

		// The variable names is composed of a (i,j) pair
		y00 = ThermoTables_Pressure[0][middleJ];
		y10 = ThermoTables_Pressure[nTable_Zone_Stations[1] - 1][middleJ];
		x00 = ThermoTables_Density[0][middleJ];
		x10 = ThermoTables_Density[nTable_Zone_Stations[1] - 1][middleJ];
		//Using the input pressure and a given middleJ value,
		//the corresponding density and iIndex may be determined
		rho = (P - y00) * (x10 - x00) / (y10 - y00) + x00;
		//If density is greater than the limits of the table,
		//then search the upper part, if it is lower than the limits,
		//search the lower part of the table
		if (rho > Density_Table_Limits[1]) {
			if (y10 > y00) {
				LowerJ = middleJ;
			} else if (y10 < y00) {
				UpperJ = middleJ;
			}
		} else if (rho < Density_Table_Limits[0]) {
			if (y10 > y00) {
				UpperJ = middleJ;
			} else if (y10 < y00) {
				LowerJ = middleJ;
			}
		} else {
			//If the density is within limits, calculate the iIndex that
			//it corresponds to
			UpperI = nTable_Zone_Stations[1] - 1;
			LowerI = 0;
			Search_NonEquispaced_Rho_Index(rho);
			//Now an i,j index pair has been found, and can be used to
			//guide the recursion further
			y00 = ThermoTables_X[LowerI][middleJ];
			y10 = ThermoTables_X[UpperI][middleJ];
			x00 = ThermoTables_Density[LowerI][middleJ];
			x10 = ThermoTables_Density[UpperI][middleJ];

			RunVal = y00 + (y10 - y00) / (x10 - x00) * (rho - x00);
			grad = ThermoTables_X[LowerI][middleJ + 1] - y00;
			if (RunVal * grad > x * grad) {
				UpperJ = middleJ;
			} else if (RunVal * grad < x * grad) {
				LowerJ = middleJ;
			} else if (RunVal == x) {
				LowerJ = middleJ;
				UpperJ = LowerJ + 1;
				break;
			}
		}
	}
}

void CLookUpTable::SetTDState_rhoe(su2double rho, su2double e) {
	// Check if inputs are in total range (necessary but not sufficient condition)
	if (rank == 12201 and LUT_Debug_Mode) {
		if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
			cerr << "RHOE Input Density out of bounds\n";
		}
		if ((e > StaticEnergy_Table_Limits[1])
				or (e < StaticEnergy_Table_Limits[0])) {
			cerr << "RHOE Input StaticEnergy out of bounds\n";
		}
	}

	// Starting values for the search
	UpperJ = nTable_Zone_Stations[0] - 1;
	LowerJ = 0;
	UpperI = nTable_Zone_Stations[1] - 1;
	LowerI = 0;

	//Searches involving density are independent of the skewness of the table
	//Fix the i index first
	Search_NonEquispaced_Rho_Index(rho);
	//Having found the i index, search in j
	Search_j_for_Y_given_i(rho, e, ThermoTables_Density,
			ThermoTables_StaticEnergy);

	//Now use the quadrilateral which contains the point to interpolate
	//Determine the interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(rho, e, ThermoTables_Density,
			ThermoTables_StaticEnergy, "RHOE");

	//Interpolate the fluid properties
	StaticEnergy = e;
	Density = rho;
	Entropy = Interpolate_2D_Bilinear(ThermoTables_Entropy);
	Pressure = Interpolate_2D_Bilinear(ThermoTables_Pressure);
	Enthalpy = Interpolate_2D_Bilinear(ThermoTables_Enthalpy);
	SoundSpeed2 = Interpolate_2D_Bilinear(ThermoTables_SoundSpeed2);
	Temperature = Interpolate_2D_Bilinear(ThermoTables_Temperature);
	dPdrho_e = Interpolate_2D_Bilinear(ThermoTables_dPdrho_e);
	dPde_rho = Interpolate_2D_Bilinear(ThermoTables_dPde_rho);
	dTdrho_e = Interpolate_2D_Bilinear(ThermoTables_dTdrho_e);
	dTde_rho = Interpolate_2D_Bilinear(ThermoTables_dTde_rho);
	Cp = Interpolate_2D_Bilinear(ThermoTables_Cp);
	Mu = Interpolate_2D_Bilinear(ThermoTables_Mu);
	Kt = Interpolate_2D_Bilinear(ThermoTables_Kt);

	//Check that the interpolated density and pressure are within LUT limits
	Check_Interpolated_PRHO_Limits("RHOE");
}

void CLookUpTable::SetTDState_PT(su2double P, su2double T) {
	// Check if inputs are in total range (necessary but not sufficient condition)
	if (rank == 12201) {
		if ((P > Pressure_Table_Limits[1]) or (P < Pressure_Table_Limits[0])) {
			cerr << "PT Input Pressure out of bounds\n";
		}
		if ((T > Temperature_Table_Limits[1])
				or (T < Temperature_Table_Limits[0])) {
			cerr << "PT Input Temperature out of bounds\n";
		}
	}

	UpperJ = nTable_Zone_Stations[0] - 1;
	LowerJ = 0;
	UpperI = nTable_Zone_Stations[1] - 1;
	LowerI = 0;

	//Determine the j index first
	Search_NonEquispaced_P_Index(P);
	//Having fixes the j index, find the i index
	Search_i_for_X_given_j(T, P, ThermoTables_Temperature, ThermoTables_Pressure);

	Search_Linear_Skewed_Table(T, P, ThermoTables_Temperature);
	//Finish of the search to check that the correct quad has indeed been selected
	//Although brief inspection shows correct element to be always selected before zig zag.

	//Determine interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(T, P, ThermoTables_Temperature,
			ThermoTables_Pressure, "PT");

	//Interpolate the fluid properties
	Pressure = P;
	Temperature = T;
	StaticEnergy = Interpolate_2D_Bilinear(ThermoTables_StaticEnergy);
	Enthalpy = Interpolate_2D_Bilinear(ThermoTables_Enthalpy);
	Entropy = Interpolate_2D_Bilinear(ThermoTables_Entropy);
	Density = Interpolate_2D_Bilinear(ThermoTables_Density);
	SoundSpeed2 = Interpolate_2D_Bilinear(ThermoTables_SoundSpeed2);
	dPdrho_e = Interpolate_2D_Bilinear(ThermoTables_dPdrho_e);
	dPde_rho = Interpolate_2D_Bilinear(ThermoTables_dPde_rho);
	dTdrho_e = Interpolate_2D_Bilinear(ThermoTables_dTdrho_e);
	dTde_rho = Interpolate_2D_Bilinear(ThermoTables_dTde_rho);
	Cp = Interpolate_2D_Bilinear(ThermoTables_Cp);
	Mu = Interpolate_2D_Bilinear(ThermoTables_Mu);
	Kt = Interpolate_2D_Bilinear(ThermoTables_Kt);

	//Check that the interpolated density and pressure are within LUT limits
	Check_Interpolated_PRHO_Limits("PT");
}

void CLookUpTable::SetTDState_Prho(su2double P, su2double rho) {
	// Check if inputs are in total range (necessary and sufficient condition)
	if (rank == 12201 and LUT_Debug_Mode) {
		if ((P > Pressure_Table_Limits[1]) or (P < Pressure_Table_Limits[0])) {
			cerr << "PRHO Input Pressure out of bounds\n";
		}
		if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
			cerr << "PRHO Input Density out of bounds\n";
		}
	}

	UpperJ = nTable_Zone_Stations[0] - 1;
	LowerJ = 0;
	UpperI = nTable_Zone_Stations[1] - 1;
	LowerI = 0;

	Search_NonEquispaced_Rho_Index(rho);
	Search_NonEquispaced_P_Index(P);

	Search_NonEquispaced_Rho_Index(rho);
	Search_j_for_Y_given_i(rho, P, ThermoTables_Density, ThermoTables_Pressure);

	//Determine interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(rho, P, ThermoTables_Density,
			ThermoTables_Pressure, "PRHO");
	//Interpolate the fluid properties
	Pressure = P;
	Density = rho;
	StaticEnergy = Interpolate_2D_Bilinear(ThermoTables_StaticEnergy);
	Enthalpy = Interpolate_2D_Bilinear(ThermoTables_Enthalpy);
	Entropy = Interpolate_2D_Bilinear(ThermoTables_Entropy);
	SoundSpeed2 = Interpolate_2D_Bilinear(ThermoTables_SoundSpeed2);
	Temperature = Interpolate_2D_Bilinear(ThermoTables_Temperature);
	dPdrho_e = Interpolate_2D_Bilinear(ThermoTables_dPdrho_e);
	dPde_rho = Interpolate_2D_Bilinear(ThermoTables_dPde_rho);
	dTdrho_e = Interpolate_2D_Bilinear(ThermoTables_dTdrho_e);
	dTde_rho = Interpolate_2D_Bilinear(ThermoTables_dTde_rho);
	Cp = Interpolate_2D_Bilinear(ThermoTables_Cp);
	Mu = Interpolate_2D_Bilinear(ThermoTables_Mu);
	Kt = Interpolate_2D_Bilinear(ThermoTables_Kt);

	//Check that the interpolated density and pressure are within LUT limits
	Check_Interpolated_PRHO_Limits("PRHO");
}

void CLookUpTable::SetEnergy_Prho(su2double P, su2double rho) {
	// Check if inputs are in total range (necessary and sufficient condition)
	if (rank == 12201 and LUT_Debug_Mode) {
		if ((P > Pressure_Table_Limits[1]) or (P < Pressure_Table_Limits[0])) {
			cerr << "PRHO Input Pressure out of bounds\n";
		}
		if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
			cerr << "PRHO Input Density out of bounds\n";
		}
	}

	UpperJ = nTable_Zone_Stations[0] - 1;
	LowerJ = 0;
	UpperI = nTable_Zone_Stations[1] - 1;
	LowerI = 0;

	Search_NonEquispaced_Rho_Index(rho);
	Search_NonEquispaced_P_Index(P);

	Search_NonEquispaced_Rho_Index(rho);
	Search_j_for_Y_given_i(rho, P, ThermoTables_Density, ThermoTables_Pressure);

	//Determine interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(rho, P, ThermoTables_Density,
			ThermoTables_Pressure, "PRHO");
	StaticEnergy = Interpolate_2D_Bilinear(ThermoTables_StaticEnergy);
	Pressure = P;
	Density = rho;

	Check_Interpolated_PRHO_Limits("PRHO (energy)");
}

void CLookUpTable::SetTDState_hs(su2double h, su2double s) {
	// Check if inputs are in total range (necessary and sufficient condition)
	if (rank == 12201 and LUT_Debug_Mode) {
		if ((h > Enthalpy_Table_Limits[1]) or (h < Enthalpy_Table_Limits[0])) {
			cerr << "HS Input Enthalpy out of bounds\n";
		}
		if ((s > Entropy_Table_Limits[1]) or (s < Entropy_Table_Limits[0])) {
			cerr << "HS Input Entropy out of bounds\n";
		}
	}

	//Preset the distance variables to something large, so they can be subsituted
	//by any point in the table.

	//Determine interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(h, s, ThermoTables_Enthalpy,
			ThermoTables_Entropy, "HS");

	//Interpolate the fluid properties
	Enthalpy = h;
	Entropy = s;
	StaticEnergy = Interpolate_2D_Bilinear(ThermoTables_StaticEnergy);
	Pressure = Interpolate_2D_Bilinear(ThermoTables_Pressure);
	Density = Interpolate_2D_Bilinear(ThermoTables_Density);
	SoundSpeed2 = Interpolate_2D_Bilinear(ThermoTables_SoundSpeed2);
	Temperature = Interpolate_2D_Bilinear(ThermoTables_Temperature);
	dPdrho_e = Interpolate_2D_Bilinear(ThermoTables_dPdrho_e);
	dPde_rho = Interpolate_2D_Bilinear(ThermoTables_dPde_rho);
	dTdrho_e = Interpolate_2D_Bilinear(ThermoTables_dTdrho_e);
	dTde_rho = Interpolate_2D_Bilinear(ThermoTables_dTde_rho);
	Cp = Interpolate_2D_Bilinear(ThermoTables_Cp);
	Mu = Interpolate_2D_Bilinear(ThermoTables_Mu);
	Kt = Interpolate_2D_Bilinear(ThermoTables_Kt);

	Check_Interpolated_PRHO_Limits("HS");

}

void CLookUpTable::SetTDState_Ps(su2double P, su2double s) {
	// Check if inputs are in total range (necessary and sufficient condition)
	if (rank == 12201 and LUT_Debug_Mode) {
		if ((P > Pressure_Table_Limits[1]) or (P < Pressure_Table_Limits[0])) {
			cerr << "PS Input Pressure out of bounds\n";
		}
		if ((s > Entropy_Table_Limits[1]) or (s < Entropy_Table_Limits[0])) {
			cerr << "PS Input Entropy  out of bounds\n";
		}
	}

	UpperJ = nTable_Zone_Stations[0] - 1;
	LowerJ = 0;
	UpperI = nTable_Zone_Stations[1] - 1;
	LowerI = 0;

	//Determine interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(s, P, ThermoTables_Entropy,
			ThermoTables_Pressure, "PS");

	//Interpolate the fluid properties
	Entropy = s;
	Pressure = P;
	StaticEnergy = Interpolate_2D_Bilinear(ThermoTables_StaticEnergy);
	Enthalpy = Interpolate_2D_Bilinear(ThermoTables_Enthalpy);
	Density = Interpolate_2D_Bilinear(ThermoTables_Density);
	SoundSpeed2 = Interpolate_2D_Bilinear(ThermoTables_SoundSpeed2);
	Temperature = Interpolate_2D_Bilinear(ThermoTables_Temperature);
	dPdrho_e = Interpolate_2D_Bilinear(ThermoTables_dPdrho_e);
	dPde_rho = Interpolate_2D_Bilinear(ThermoTables_dPde_rho);
	dTdrho_e = Interpolate_2D_Bilinear(ThermoTables_dTdrho_e);
	dTde_rho = Interpolate_2D_Bilinear(ThermoTables_dTde_rho);
	Cp = Interpolate_2D_Bilinear(ThermoTables_Cp);
	Mu = Interpolate_2D_Bilinear(ThermoTables_Mu);
	Kt = Interpolate_2D_Bilinear(ThermoTables_Kt);

	//Check that the interpolated density and pressure are within LUT limits
	Check_Interpolated_PRHO_Limits("PS");
}

void CLookUpTable::SetTDState_rhoT(su2double rho, su2double T) {
	// Check if inputs are in total range (necessary and sufficient condition)
	if (rank == 12201 and LUT_Debug_Mode) {
		if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
			cerr << "RHOT Input Density out of bounds\n";
		}
		if ((T > Temperature_Table_Limits[1])
				or (T < Temperature_Table_Limits[0])) {
			cerr << "RHOT Input Temperature out of bounds\n";
		}
	}
	// Linear interpolation requires 4 neighbors to be selected from the LUT

	UpperJ = nTable_Zone_Stations[0] - 1;
	LowerJ = 0;
	UpperI = nTable_Zone_Stations[1] - 1;
	LowerI = 0;

	Search_NonEquispaced_Rho_Index(rho);
	Search_j_for_Y_given_i(rho, T, ThermoTables_Density,
			ThermoTables_Temperature);

	//Determine the interpolation coefficients
	Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(rho, T, ThermoTables_Density,
			ThermoTables_Temperature, "RHOT");

	//Interpolate the fluid properties
	Temperature = T;
	Density = rho;
	StaticEnergy = Interpolate_2D_Bilinear(ThermoTables_StaticEnergy);
	Enthalpy = Interpolate_2D_Bilinear(ThermoTables_Enthalpy);
	Entropy = Interpolate_2D_Bilinear(ThermoTables_Entropy);
	Pressure = Interpolate_2D_Bilinear(ThermoTables_Pressure);
	SoundSpeed2 = Interpolate_2D_Bilinear(ThermoTables_SoundSpeed2);
	dPdrho_e = Interpolate_2D_Bilinear(ThermoTables_dPdrho_e);
	dPde_rho = Interpolate_2D_Bilinear(ThermoTables_dPde_rho);
	dTdrho_e = Interpolate_2D_Bilinear(ThermoTables_dTdrho_e);
	dTde_rho = Interpolate_2D_Bilinear(ThermoTables_dTde_rho);
	Cp = Interpolate_2D_Bilinear(ThermoTables_Cp);
	Mu = Interpolate_2D_Bilinear(ThermoTables_Mu);
	Kt = Interpolate_2D_Bilinear(ThermoTables_Kt);

	//Check that the interpolated density and pressure are within LUT limits
	Check_Interpolated_PRHO_Limits("RHOT");
}

void CLookUpTable::Check_Interpolated_PRHO_Limits(string interpolation_case) {
	//Check that the interpolated density and pressure are within LUT limits
	if (rank == 12201 and LUT_Debug_Mode) {
		if ((Density > Density_Table_Limits[1])
				or (Density < Density_Table_Limits[0])) {
			cerr << interpolation_case << " Interpolated Density out of bounds\n";
		}
		if ((Pressure > Pressure_Table_Limits[1])
				or (Pressure < Pressure_Table_Limits[0])) {
			cerr << interpolation_case << " Interpolated Pressure out of bounds\n";
		}
	}
}

inline void CLookUpTable::Gaussian_Inverse(int nDim) {
	//A temporary matrix to hold the inverse, dynamically allocated
	su2double **temp = new su2double*[nDim];
	for (int i = 0; i < nDim; i++) {
		temp[i] = new su2double[2 * nDim];
	}

	//Copy the desired matrix into the temporary matrix
	for (int i = 0; i < nDim; i++) {
		for (int j = 0; j < nDim; j++) {
			temp[i][j] = Interpolation_Matrix[i][j];
			temp[i][nDim + j] = 0;
		}
		temp[i][nDim + i] = 1;
	}

	su2double max_val;
	int max_idx;
	//Pivot each column such that the largest number possible divides the oter rows
	//The goal is to avoid zeros or small numbers in division.
	for (int k = 0; k < nDim - 1; k++) {
		max_idx = k;
		max_val = abs(temp[k][k]);
		//Find the largest value (pivot) in the column
		for (int j = k; j < nDim; j++) {
			if (abs(temp[j][k]) > max_val) {
				max_idx = j;
				max_val = abs(temp[j][k]);
			}
		}
		//Move the row with the highest value up
		for (int j = 0; j < (nDim * 2); j++) {
			su2double d = temp[k][j];
			temp[k][j] = temp[max_idx][j];
			temp[max_idx][j] = d;
		}
		//Subtract the moved row from all other rows
		for (int i = k + 1; i < nDim; i++) {
			su2double c = temp[i][k] / temp[k][k];
			for (int j = 0; j < (nDim * 2); j++) {
				temp[i][j] = temp[i][j] - temp[k][j] * c;
			}
		}
	}

	//Back-substitution
	for (int k = nDim - 1; k > 0; k--) {
		if (temp[k][k] != 0) {
			for (int i = k - 1; i > -1; i--) {
				su2double c = temp[i][k] / temp[k][k];
				for (int j = 0; j < (nDim * 2); j++) {
					temp[i][j] = temp[i][j] - temp[k][j] * c;
				}
			}
		}
	}
	//Normalize the inverse
	for (int i = 0; i < nDim; i++) {
		su2double c = temp[i][i];
		for (int j = 0; j < nDim; j++) {
			temp[i][j + nDim] = temp[i][j + nDim] / c;
		}
	}
	//Copy the inverse back to the main program flow
	for (int i = 0; i < nDim; i++) {
		for (int j = 0; j < nDim; j++) {
			Interpolation_Coeff[i][j] = temp[i][j + nDim];
		}
	}
	//Delete dynamic template
	for (int i = 0; i < nDim; i++) {
		delete[] temp[i];
	}
	delete[] temp;
	return;
}

void CLookUpTable::Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(su2double x,
su2double y, su2double **ThermoTables_X, su2double **ThermoTables_Y,
		std::string grid_var) {
	//The x,y coordinates of the quadrilateral
	su2double x00, y00, x10, x01, x11, y10, y01, y11;

	x00 = ThermoTables_X[LowerI][LowerJ];
	y00 = ThermoTables_Y[LowerI][LowerJ];
	x01 = ThermoTables_X[LowerI][UpperJ];
	y01 = ThermoTables_Y[LowerI][UpperJ];
	x10 = ThermoTables_X[UpperI][LowerJ];
	y10 = ThermoTables_Y[UpperI][LowerJ];
	x11 = ThermoTables_X[UpperI][UpperJ];
	y11 = ThermoTables_Y[UpperI][UpperJ];

	//Check if x, y is indeed in the quad
	//The (true and not false) type of logic is needed as the both monotonically
	//increasing and monotonically decreasing functions need to pass the same test
	bool BOTTOM, TOP, LEFT, RIGHT, OUT_OF_BOUNDS;
	su2double dy, dx, dx10, dy10, dx01, dy01, dx11, dy11;
	dx = x - x00;
	dy = y - y00;
	dx10 = x10 - x00;
	dy10 = y10 - y00;
	dx01 = x01 - x00;
	dy01 = y01 - y00;
	dx11 = x11 - x00;
	dy11 = y11 - y00;
	if (rank == 12201 and LUT_Debug_Mode) {
		BOTTOM = (dy * dx10) < (dx * dy10);
		TOP = ((dy - dy01) * (dx11 - dx01)) > ((dy11 - dy01) * (dx - dx01));
		RIGHT = ((dx - dx10) * (dy11 - dy10)) > ((dx11 - dx10) * (dy - dy10));
		LEFT = (dx * dy01) < (dx01 * dy);
		OUT_OF_BOUNDS = false;
		//Check BOTTOM quad boundary
		if (BOTTOM and !TOP) {
			OUT_OF_BOUNDS = true;
			//Check if the point is also outside the LUT
			if (LowerJ == 0) {
				cerr << grid_var << ' ' << LowerI << ", " << LowerJ
						<< " interpolation point lies below the LUT\n";
			} else {
				cerr << grid_var << ' ' << LowerI << ", " << LowerJ
						<< " interpolation point lies below bottom boundary of selected quad\n";
			}
		}
		//Check RIGHT quad boundary
		if (RIGHT and !LEFT) {
			OUT_OF_BOUNDS = true;
			//Check if the point is also outside the LUT
			if (LowerI == (nTable_Zone_Stations[1] - 2)) {
				cerr << grid_var << ' ' << LowerI << ", " << LowerJ
						<< " interpolation point lies right of the LUT\n";
			} else {
				cerr << grid_var << ' ' << LowerI << ", " << LowerJ
						<< " interpolation point lies to the right of the boundary of selected quad\n";
			}
		}
		//Check TOP quad boundary
		if (TOP and !BOTTOM) {
			OUT_OF_BOUNDS = true;
			//Check if the point is also outside the LUT
			if (LowerJ == (nTable_Zone_Stations[0] - 2)) {
				cerr << grid_var << ' ' << LowerI << ", " << LowerJ
						<< " interpolation point lies above the LUT\n";
			} else {
				cerr << grid_var << ' ' << LowerI << ", " << LowerJ
						<< +" interpolation point lies above the boundary of selected quad\n";
			}
		}
		//Check LEFT quad boundary
		if (LEFT and !RIGHT) {
			OUT_OF_BOUNDS = true;
			//Check if the point is also outside the LUT
			if (LowerI == 0) {
				cerr << grid_var << ' ' << LowerI << ", " << LowerJ
						<< " interpolation point lies left of the LUT\n";
			} else {
				cerr << grid_var << ' ' << LowerI << ", " << LowerJ
						<< +" interpolation point lies to the left of the boundary of selected quad\n";
			}
		}
	}

	//Setup the LHM matrix for the interpolation (Vandermonde)
	Interpolation_Matrix[0][0] = 1;
	Interpolation_Matrix[0][1] = 0;
	Interpolation_Matrix[0][2] = 0;
	Interpolation_Matrix[0][3] = 0;

	Interpolation_Matrix[1][0] = 1;
	Interpolation_Matrix[1][1] = x10 - x00;
	Interpolation_Matrix[1][2] = y10 - y00;
	Interpolation_Matrix[1][3] = (x10 - x00) * (y10 - y00);

	Interpolation_Matrix[2][0] = 1;
	Interpolation_Matrix[2][1] = x01 - x00;
	Interpolation_Matrix[2][2] = y01 - y00;
	Interpolation_Matrix[2][3] = (x01 - x00) * (y01 - y00);

	Interpolation_Matrix[3][0] = 1;
	Interpolation_Matrix[3][1] = x11 - x00;
	Interpolation_Matrix[3][2] = y11 - y00;
	Interpolation_Matrix[3][3] = (x11 - x00) * (y11 - y00);

	//Invert the Interpolation matrix using Gaussian elimination with pivoting
	Gaussian_Inverse(4);
	su2double d;

	//Transpose the inverse
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 4; j++) {
			d = Interpolation_Coeff[i][j];
			Interpolation_Coeff[i][j] = Interpolation_Coeff[j][i];
			Interpolation_Coeff[j][i] = d;
		}
	}
	//The transpose allows the same coefficients to be used
	// for all Thermo variables (need only 4 coefficients)
	for (int i = 0; i < 4; i++) {
		d = 0;
		d = d + Interpolation_Coeff[i][0] * 1;
		d = d + Interpolation_Coeff[i][1] * (x - x00);
		d = d + Interpolation_Coeff[i][2] * (y - y00);
		d = d + Interpolation_Coeff[i][3] * (x - x00) * (y - y00);
		Interpolation_Coeff[i][0] = d;
	}
	return;
}

su2double CLookUpTable::Interpolate_2D_Bilinear(su2double * *ThermoTables_Z) {
	//The function values at the 4 corners of the quad
	su2double func_value_at_i0j0, func_value_at_i1j0, func_value_at_i0j1,
			func_value_at_i1j1;

	func_value_at_i0j0 = ThermoTables_Z[LowerI][LowerJ];
	func_value_at_i1j0 = ThermoTables_Z[UpperI][LowerJ];
	func_value_at_i0j1 = ThermoTables_Z[LowerI][UpperJ];
	func_value_at_i1j1 = ThermoTables_Z[UpperI][UpperJ];
	//The Interpolation_Coeff values depend on location alone
	//and are the same regardless of function values
	su2double result = 0;
	result = result + Interpolation_Coeff[0][0] * func_value_at_i0j0;
	result = result + Interpolation_Coeff[1][0] * func_value_at_i1j0;
	result = result + Interpolation_Coeff[2][0] * func_value_at_i0j1;
	result = result + Interpolation_Coeff[3][0] * func_value_at_i1j1;

	return result;
}

void CLookUpTable::RecordState(char* file) {
	//Record the state of the fluid model to a file for
	//verificaiton purposes
	fstream fs;
	fs.open(file, fstream::app);
	fs.precision(17);
	assert(fs.is_open());
	fs << Temperature << ", ";
	fs << Density << ", ";
	fs << Enthalpy << ", ";
	fs << StaticEnergy << ", ";
	fs << Entropy << ", ";
	fs << Pressure << ", ";
	fs << SoundSpeed2 << ", ";
	fs << dPdrho_e << ", ";
	fs << dPde_rho << ", ";
	fs << dTdrho_e << ", ";
	fs << dTde_rho << ", ";
	fs << Cp << ", ";
	fs << Mu << ", ";
	//fs << dmudrho_T << ", ";
	//fs << dmudT_rho << ", ";
	fs << Kt << " ";
	//fs << dktdrho_T << ", ";
	//fs << dktdT_rho << ", ";
	fs << "\n";
	fs.close();
}

void CLookUpTable::LookUpTable_Print_To_File(char* filename) {
	//Print the entire table to a file such that the mesh can be plotted
	//externally (for verification purposes)
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nTable_Zone_Stations[i]; j++) {
			Temperature = ThermoTables_Temperature[i][j];
			Density = ThermoTables_Density[i][j];
			Enthalpy = ThermoTables_Enthalpy[i][j];
			StaticEnergy = ThermoTables_StaticEnergy[i][j];
			Entropy = ThermoTables_Entropy[i][j];
			Pressure = ThermoTables_Pressure[i][j];
			SoundSpeed2 = ThermoTables_SoundSpeed2[i][j];
			dPdrho_e = ThermoTables_dPdrho_e[i][j];
			dPde_rho = ThermoTables_dPde_rho[i][j];
			dTdrho_e = ThermoTables_dTdrho_e[i][j];
			dTde_rho = ThermoTables_dTde_rho[i][j];
			Cp = ThermoTables_Cp[i][j];
			Kt = ThermoTables_Kt[i][j];
			Mu = ThermoTables_Mu[i][j];
			RecordState(filename);
		}
	}

}

void CLookUpTable::LookUpTable_Load_TEC(std::string filename) {
	string line;
	string value;
	int found;
	int zone_scanned;

	ifstream table(filename.c_str());
	if (!table.is_open()) {
		if (rank == 12201) {
			cout << "The LUT file appears to be missing!! " << filename << endl;
		}
		exit(EXIT_FAILURE);
	}
	zone_scanned = 0;
	//Go through all lines in the table file.
	getline(table, line);	//Skip the header
	while (getline(table, line)) {
		found = line.find("ZONE");
		if (found != -1) {
			cout << line << endl;
			istringstream in(line);
			//Note down the dimensions of the table
			int nPoints_in_Zone, nTriangles_in_Zone;
			string c1, c2, c3, c4;
			in >> c1 >> c2 >> nPoints_in_Zone >> c3 >> c4 >> nTriangles_in_Zone;
			cout << nPoints_in_Zone << "  " << nTriangles_in_Zone << endl;
			//Create the actual LUT of CThermoLists which is used in the FluidModel
			nTable_Zone_Stations[zone_scanned] = nPoints_in_Zone;
			nTable_Zone_Triangles[zone_scanned] = nTriangles_in_Zone;
			//Allocate the memory for the table
			LookUpTable_Malloc(zone_scanned);

			//Load the values of the themordynamic properties at each table station
			for (int j = 0; j < nTable_Zone_Stations[zone_scanned]; j++) {
				getline(table, line);
				istringstream in(line);
				in >> ThermoTables_Density[zone_scanned][j];
				in >> ThermoTables_Pressure[zone_scanned][j];
				in >> ThermoTables_SoundSpeed2[zone_scanned][j];
				in >> ThermoTables_Cp[zone_scanned][j];
				in >> ThermoTables_Entropy[zone_scanned][j];
				in >> ThermoTables_Mu[zone_scanned][j];
				in >> ThermoTables_Kt[zone_scanned][j];
				in >> ThermoTables_dPdrho_e[zone_scanned][j];
				in >> ThermoTables_dPde_rho[zone_scanned][j];
				in >> ThermoTables_dTdrho_e[zone_scanned][j];
				in >> ThermoTables_dTde_rho[zone_scanned][j];
				in >> ThermoTables_Temperature[zone_scanned][j];
				in >> ThermoTables_StaticEnergy[zone_scanned][j];
				in >> ThermoTables_Enthalpy[zone_scanned][j];
			}
			//Skip empty line
			getline(table, line);
			//Load the triangles i.e. how the data point in each zone are connected
			for (int j = 0; j < nTable_Zone_Triangles[zone_scanned]; j++) {
				getline(table, line);
				istringstream in(line);
				in >> Table_Zone_Triangles[zone_scanned][j][0]
						>> Table_Zone_Triangles[zone_scanned][j][1]
						>> Table_Zone_Triangles[zone_scanned][j][2];
			}
			zone_scanned++;
		}
	}

	table.close();
	//NonDimensionalise and find limits
	NonDimensionalise_Table_Values();
	Find_Table_Limits();
}

void CLookUpTable::LookUpTable_Malloc(int Index_of_Zone) {
	ThermoTables_StaticEnergy[Index_of_Zone] =
			new su2double[nTable_Zone_Stations[Index_of_Zone]];
	ThermoTables_Entropy[Index_of_Zone] =
			new su2double[nTable_Zone_Stations[Index_of_Zone]];
	ThermoTables_Enthalpy[Index_of_Zone] =
			new su2double[nTable_Zone_Stations[Index_of_Zone]];
	ThermoTables_Density[Index_of_Zone] =
			new su2double[nTable_Zone_Stations[Index_of_Zone]];
	ThermoTables_Pressure[Index_of_Zone] =
			new su2double[nTable_Zone_Stations[Index_of_Zone]];
	ThermoTables_SoundSpeed2[Index_of_Zone] =
			new su2double[nTable_Zone_Stations[Index_of_Zone]];
	ThermoTables_Temperature[Index_of_Zone] =
			new su2double[nTable_Zone_Stations[Index_of_Zone]];
	ThermoTables_dPdrho_e[Index_of_Zone] =
			new su2double[nTable_Zone_Stations[Index_of_Zone]];
	ThermoTables_dPde_rho[Index_of_Zone] =
			new su2double[nTable_Zone_Stations[Index_of_Zone]];
	ThermoTables_dTdrho_e[Index_of_Zone] =
			new su2double[nTable_Zone_Stations[Index_of_Zone]];
	ThermoTables_dTde_rho[Index_of_Zone] =
			new su2double[nTable_Zone_Stations[Index_of_Zone]];
	ThermoTables_Cp[Index_of_Zone] =
			new su2double[nTable_Zone_Stations[Index_of_Zone]];
	ThermoTables_Mu[Index_of_Zone] =
			new su2double[nTable_Zone_Stations[Index_of_Zone]];
	ThermoTables_dmudrho_T[Index_of_Zone] =
			new su2double[nTable_Zone_Stations[Index_of_Zone]];
	ThermoTables_dmudT_rho[Index_of_Zone] =
			new su2double[nTable_Zone_Stations[Index_of_Zone]];
	ThermoTables_Kt[Index_of_Zone] =
			new su2double[nTable_Zone_Stations[Index_of_Zone]];
	ThermoTables_dktdrho_T[Index_of_Zone] =
			new su2double[nTable_Zone_Stations[Index_of_Zone]];
	ThermoTables_dktdT_rho[Index_of_Zone] =
			new su2double[nTable_Zone_Stations[Index_of_Zone]];
	Table_Zone_Triangles[Index_of_Zone] =
			new int*[nTable_Zone_Triangles[Index_of_Zone]];
	for (int j = 0; j < nTable_Zone_Triangles[Index_of_Zone]; j++) {
		Table_Zone_Triangles[Index_of_Zone][j] = new int[3];
	}
}

void CLookUpTable::Find_Table_Limits() {
	Density_Table_Limits[0] = HUGE_VAL;
	Density_Table_Limits[1] = -HUGE_VAL;
	Pressure_Table_Limits[0] = HUGE_VAL;	//lower limit
	Pressure_Table_Limits[1] = -HUGE_VAL;	//upper limit
	Enthalpy_Table_Limits[0] = HUGE_VAL;	//lower limit
	Enthalpy_Table_Limits[1] = -HUGE_VAL;	//upper limit
	SoundSpeed2_Table_Limits[0] = HUGE_VAL;	//lower limit
	SoundSpeed2_Table_Limits[1] = -HUGE_VAL;	//upper limit
	Cp_Table_Limits[0] = HUGE_VAL;	//lower limit
	Cp_Table_Limits[1] = -HUGE_VAL;	//upper limit
	Entropy_Table_Limits[0] = HUGE_VAL;	//lower limit
	Entropy_Table_Limits[1] = -HUGE_VAL;	//upper limit
	dPdrho_e_Table_Limits[0] = HUGE_VAL;	//lower limit
	dPdrho_e_Table_Limits[1] = -HUGE_VAL;	//upper limit
	dPde_rho_Table_Limits[0] = HUGE_VAL;	//lower limit
	dPde_rho_Table_Limits[1] = -HUGE_VAL;	//upper limit
	dTdrho_e_Table_Limits[0] = HUGE_VAL;	//lower limit
	dTdrho_e_Table_Limits[1] = -HUGE_VAL;	//upper limit
	dTde_rho_Table_Limits[0] = HUGE_VAL;	//lower limit
	dTde_rho_Table_Limits[1] = -HUGE_VAL;	//upper limit
	Temperature_Table_Limits[0] = HUGE_VAL;	//lower limit
	Temperature_Table_Limits[1] = -HUGE_VAL;	//upper limit
	StaticEnergy_Table_Limits[0] = HUGE_VAL;	//lower limit
	StaticEnergy_Table_Limits[1] = -HUGE_VAL;	//upper limit
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nTable_Zone_Stations[i]; j++) {
			//The table limits are stored for later checks of values becoming inconsistent
			if (ThermoTables_Density[i][j] > Density_Table_Limits[1]) {
				Density_Table_Limits[1] = ThermoTables_Density[i][j];
			}
			if (ThermoTables_Density[i][j] < Density_Table_Limits[0]) {
				Density_Table_Limits[0] = ThermoTables_Density[i][j];
			}

			if (ThermoTables_Pressure[i][j] > Pressure_Table_Limits[1]) {
				Pressure_Table_Limits[1] = ThermoTables_Pressure[i][j];
			}
			if (ThermoTables_Pressure[i][j] < Pressure_Table_Limits[0]) {
				Pressure_Table_Limits[0] = ThermoTables_Pressure[i][j];
			}
			if (ThermoTables_Enthalpy[i][j] > Enthalpy_Table_Limits[1]) {
				Enthalpy_Table_Limits[1] = ThermoTables_Enthalpy[i][j];
			}
			if (ThermoTables_Enthalpy[i][j] < Enthalpy_Table_Limits[0]) {
				Enthalpy_Table_Limits[0] = ThermoTables_Enthalpy[i][j];
			}
			if (ThermoTables_SoundSpeed2[i][j] > SoundSpeed2_Table_Limits[1]) {
				SoundSpeed2_Table_Limits[1] = ThermoTables_SoundSpeed2[i][j];
			}
			if (ThermoTables_SoundSpeed2[i][j] < SoundSpeed2_Table_Limits[0]) {
				SoundSpeed2_Table_Limits[0] = ThermoTables_SoundSpeed2[i][j];
			}
			if (ThermoTables_Cp[i][j] > Cp_Table_Limits[1]) {
				Cp_Table_Limits[1] = ThermoTables_Cp[i][j];
			}
			if (ThermoTables_Cp[i][j] < Cp_Table_Limits[0]) {
				Cp_Table_Limits[0] = ThermoTables_Cp[i][j];
			}
			if (ThermoTables_Entropy[i][j] > Entropy_Table_Limits[1]) {
				Entropy_Table_Limits[1] = ThermoTables_Entropy[i][j];
			}
			if (ThermoTables_Entropy[i][j] < Entropy_Table_Limits[0]) {
				Entropy_Table_Limits[0] = ThermoTables_Entropy[i][j];
			}
			if (ThermoTables_dPdrho_e[i][j] > dPdrho_e_Table_Limits[1]) {
				dPdrho_e_Table_Limits[1] = ThermoTables_dPdrho_e[i][j];
			}
			if (ThermoTables_dPdrho_e[i][j] < dPdrho_e_Table_Limits[0]) {
				dPdrho_e_Table_Limits[0] = ThermoTables_dPdrho_e[i][j];
			}
			if (ThermoTables_dPde_rho[i][j] > dPde_rho_Table_Limits[1]) {
				dPde_rho_Table_Limits[1] = ThermoTables_dPde_rho[i][j];
			}
			if (ThermoTables_dPde_rho[i][j] < dPde_rho_Table_Limits[0]) {
				dPde_rho_Table_Limits[0] = ThermoTables_dPde_rho[i][j];
			}
			if (ThermoTables_dTdrho_e[i][j] > dTdrho_e_Table_Limits[1]) {
				dTdrho_e_Table_Limits[1] = ThermoTables_dTdrho_e[i][j];
			}
			if (ThermoTables_dTdrho_e[i][j] < dTdrho_e_Table_Limits[0]) {
				dTdrho_e_Table_Limits[0] = ThermoTables_dTdrho_e[i][j];
			}
			if (ThermoTables_dTde_rho[i][j] > dTde_rho_Table_Limits[1]) {
				dTde_rho_Table_Limits[1] = ThermoTables_dTde_rho[i][j];
			}
			if (ThermoTables_dTde_rho[i][j] < dTde_rho_Table_Limits[0]) {
				dTde_rho_Table_Limits[0] = ThermoTables_dTde_rho[i][j];
			}
			if (ThermoTables_Temperature[i][j] > Temperature_Table_Limits[1]) {
				Temperature_Table_Limits[1] = ThermoTables_Temperature[i][j];
			}
			if (ThermoTables_Temperature[i][j] < Temperature_Table_Limits[0]) {
				Temperature_Table_Limits[0] = ThermoTables_Temperature[i][j];
			}
			if (ThermoTables_StaticEnergy[i][j] > StaticEnergy_Table_Limits[1]) {
				StaticEnergy_Table_Limits[1] = ThermoTables_StaticEnergy[i][j];
			}
			if (ThermoTables_StaticEnergy[i][j] < StaticEnergy_Table_Limits[0]) {
				StaticEnergy_Table_Limits[0] = ThermoTables_StaticEnergy[i][j];
			}
		}
	}
}
void CLookUpTable::NonDimensionalise_Table_Values() {
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < nTable_Zone_Stations[i]; j++) {
			ThermoTables_Density[i][j] /= Density_Reference_Value;
			ThermoTables_Pressure[i][j] /= Pressure_Reference_Value;
			ThermoTables_SoundSpeed2[i][j] = pow(ThermoTables_SoundSpeed2[i][j], 2);
			ThermoTables_SoundSpeed2[i][j] /= pow(Velocity_Reference_Value, 2);
			ThermoTables_Cp[i][j] *= (Temperature_Reference_Value
					/ Energy_Reference_Value);
			ThermoTables_Entropy[i][j] *= (Temperature_Reference_Value
					/ Energy_Reference_Value);
			ThermoTables_dPdrho_e[i][j] *= (Density_Reference_Value
					/ Pressure_Reference_Value);
			ThermoTables_dPde_rho[i][j] *= (Energy_Reference_Value
					/ Pressure_Reference_Value);
			ThermoTables_dTdrho_e[i][j] *= (Density_Reference_Value
					/ Temperature_Reference_Value);
			ThermoTables_dTde_rho[i][j] *= (Energy_Reference_Value
					/ Temperature_Reference_Value);
			ThermoTables_Temperature[i][j] /= Temperature_Reference_Value;
			ThermoTables_StaticEnergy[i][j] /= Energy_Reference_Value;
			ThermoTables_Enthalpy[i][j] /= Energy_Reference_Value;
		}
	}
}

