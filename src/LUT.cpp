#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <string>
#include <cstring>
#include <cmath>
#include <cassert>
#include <algorithm>
#include "stdlib.h"
#include "stdio.h"
#include "LUT.hpp"
#include <iomanip>
#include <vector>
using namespace std;

CTrapezoidalMap::CTrapezoidalMap() {
}
CTrapezoidalMap::CTrapezoidalMap(vector< su2double > const &x_samples,
		vector< su2double > const &y_samples,
		vector<vector<int> > const &unique_edges) {
	rank = 12201;

#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	clock_t build_start = clock();

	Unique_X_Bands = x_samples; //copy the x_values in
	//Sort the x_bands and make them unique
	sort(Unique_X_Bands.begin(), Unique_X_Bands.end());
	vector<su2double>::iterator iter;
	iter = unique(Unique_X_Bands.begin(), Unique_X_Bands.end());
	Unique_X_Bands.resize(distance(Unique_X_Bands.begin(), iter));
	X_Limits_of_Edges.resize(unique_edges.size(), vector<su2double>(2, 0));
	Y_Limits_of_Edges.resize(unique_edges.size(), vector<su2double>(2, 0));

	//Store the x and y values of each edge into a vector for a slight speed up as it
	//prevents some uncoalesced accesses
	for (unsigned int j = 0; j < unique_edges.size(); j++) {
		X_Limits_of_Edges[j][0] = x_samples[unique_edges[j][0]];
		X_Limits_of_Edges[j][1] = x_samples[unique_edges[j][1]];
		Y_Limits_of_Edges[j][0] = y_samples[unique_edges[j][0]];
		Y_Limits_of_Edges[j][1] = y_samples[unique_edges[j][1]];
	}

	//How many bands to search?
	int b_max = Unique_X_Bands.size() - 1;
	//Start with band 0, obviously
	int b = 0;
	//How many edges to check for intersection with the band?
	int e_max = unique_edges.size();
	//Start with edge indexes as 0.
	int i = 0;
	//Count the how many edges intersect a band
	int k = 0;
	//The high and low x value of each band
	su2double x_low = 0;
	su2double x_hi = 0;

	//Store the y_values of edges as required for searching
	Y_Values_of_Edge_Within_Band_And_Index.resize(Unique_X_Bands.size() - 1);

	//Check which edges intersect the band
	while (b < (b_max)) {
		x_low = Unique_X_Bands[b];
		x_hi = Unique_X_Bands[b + 1];
		i = 0;
		k = 0;
		//This while loop determined which edges appear in a paritcular band
		//The index of the edge being tested is 'i'
		while (i < e_max) {
			//Check if edge intersects the band (vertical edges are automatically discared)
			if (((X_Limits_of_Edges[i][0] <= x_low)
					and (X_Limits_of_Edges[i][1] >= x_hi))
					or ((X_Limits_of_Edges[i][1] <= x_low)
							and (X_Limits_of_Edges[i][0] >= x_hi))) {
				Y_Values_of_Edge_Within_Band_And_Index[b].push_back(make_pair(0.0, 0));
				//Save the edge index so it can latter be recalled (when searching)
				Y_Values_of_Edge_Within_Band_And_Index[b][k].second = i;
				//Determine what y value the edge takes in the middle of the band
				Y_Values_of_Edge_Within_Band_And_Index[b][k].first =
						Y_Limits_of_Edges[i][0]
								+ (Y_Limits_of_Edges[i][1] - Y_Limits_of_Edges[i][0])
										/ (X_Limits_of_Edges[i][1] - X_Limits_of_Edges[i][0])
										* ((x_low + x_hi) / 2.0 - X_Limits_of_Edges[i][0]);
				//k counts the number of edges which have been found to intersect with the band
				k++;
			}
			//increment i, which  moves the algorithm along to the next edge
			i++;
		}
		//Sort the edges in the band depending on the y values they were found to have
		//It is worth noting that these y values are unique (i.e. edges cannot intersect in a band)
		sort(Y_Values_of_Edge_Within_Band_And_Index[b].begin(),
				Y_Values_of_Edge_Within_Band_And_Index[b].end());
		//Move on to the next band of x values
		b++;
	}
	//Initialize the search to table limits
	UpperI = Unique_X_Bands.size() - 1;
	LowerI = 0;

	su2double duration = ((su2double) clock() - (su2double) build_start)
			/ ((su2double) CLOCKS_PER_SEC);
	if (rank == 12201)
		cout << duration << " seconds\n";
}

void CTrapezoidalMap::Find_Containing_Simplex(su2double x, su2double y) {
	//Find the x band in which the current x value sits
	Search_Bands_For(x);
	//Within that band find the two containing edges
	Search_Band_For_Edge(x, y);
}

void CTrapezoidalMap::Search_Bands_For(su2double x) {
	su2double x_middle, x_lower, x_upper;

	do {
		middleI = (UpperI + LowerI) / 2;
		x_middle = Unique_X_Bands[middleI];
		x_lower = Unique_X_Bands[LowerI];
		x_upper = Unique_X_Bands[UpperI];
		//Step used for restarting the search on the low end
		if (x < x_lower) {
			UpperI = LowerI;
			LowerI = LowerI / 2;
			//Step used for restarting the search on the upper end
		} else if (x > x_upper) {
			LowerI = UpperI;
			UpperI = (UpperI + (Unique_X_Bands.size() - 1)) / 2;
			//After the restart is cleared, do the normal binary search
		} else if (x < x_middle) {
			UpperI = middleI;
		} else if (x > x_middle) {
			LowerI = middleI;
		} else if (x_middle == x) {
			LowerI = middleI;
			UpperI = LowerI + 1;
			break;
		}

	} while (UpperI - LowerI > 1);
}

void CTrapezoidalMap::Search_Band_For_Edge(su2double x, su2double y) {

	su2double RunVal, y00, y10, x00, x10;
	int RunEdge;
	UpperJ = Y_Values_of_Edge_Within_Band_And_Index[LowerI].size() - 1;
	LowerJ = 0;

	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2;
		//Select the edge associated with the current x band (LowerI)
		//Search for the RunEdge in the middleJ direction (second value is index of edge)
		RunEdge = Y_Values_of_Edge_Within_Band_And_Index[LowerI][middleJ].second;
		y00 = Y_Limits_of_Edges[RunEdge][0];
		y10 = Y_Limits_of_Edges[RunEdge][1];
		x00 = X_Limits_of_Edges[RunEdge][0];
		x10 = X_Limits_of_Edges[RunEdge][1];
		//The search variable in j should be interpolated in i as well
		RunVal = y00 + (y10 - y00) / (x10 - x00) * (x - x00);
		if (RunVal > y) {
			UpperJ = middleJ;
		} else if (RunVal < y) {
			LowerJ = middleJ;
		} else if (RunVal == y) {
			LowerJ = middleJ;
			UpperJ = LowerJ + 1;
			break;
		}

	}
	if (UpperJ < (Y_Values_of_Edge_Within_Band_And_Index[LowerI].size() - 1)) {
		UpperEdge =
				Y_Values_of_Edge_Within_Band_And_Index[LowerI][UpperJ + 1].second;
		MiddleEdge = Y_Values_of_Edge_Within_Band_And_Index[LowerI][UpperJ].second;
		LowerEdge = Y_Values_of_Edge_Within_Band_And_Index[LowerI][LowerJ].second;
	} else {
		UpperEdge = Y_Values_of_Edge_Within_Band_And_Index[LowerI][UpperJ].second;
		MiddleEdge = Y_Values_of_Edge_Within_Band_And_Index[LowerI][LowerJ].second;
		LowerEdge =
				Y_Values_of_Edge_Within_Band_And_Index[LowerI][LowerJ - 1].second;
	}
}

CLookUpTable::CLookUpTable(string Filename) {
	LUT_Debug_Mode = false;
	rank = 12201;
	CurrentPoints.resize(4, 0);
	CurrentZone = 1;		//The vapor region

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
		cout << "Detecting all unique edges and setting edge to face connectivity..." << endl;
	}
	Get_Unique_Edges();
	if (rank == 12201) {
		cout << "Number of edges in zone 0: " << Table_Zone_Edges[0].size() << endl;
		cout << "Number of edges in zone 1: " << Table_Zone_Edges[1].size() << endl;

	}
//	Get_Edge_To_Face_Connectivty();

	if (rank == 12201) {
		// Building an KD_tree for the HS thermopair
		cout << "Building trapezoidal map for rhoe..." << endl;
	}

	rhoe_map[0] = CTrapezoidalMap(ThermoTables_Density[0],
			ThermoTables_StaticEnergy[0], Table_Zone_Edges[0]);
	rhoe_map[1] = CTrapezoidalMap(ThermoTables_Density[1],
			ThermoTables_StaticEnergy[1], Table_Zone_Edges[1]);

	if (rank == 12201) {
		cout << "Building trapezoidal map for Prho..." << endl;
	}
	Prho_map[0] = CTrapezoidalMap(ThermoTables_Pressure[0],
			ThermoTables_Density[0], Table_Zone_Edges[0]);
	Prho_map[1] = CTrapezoidalMap(ThermoTables_Pressure[1],
			ThermoTables_Density[1], Table_Zone_Edges[1]);

	if (rank == 12201) {
		cout << "Building trapezoidal map for hs..." << endl;
	}
	hs_map[0] = CTrapezoidalMap(ThermoTables_Enthalpy[0], ThermoTables_Entropy[0],
			Table_Zone_Edges[0]);
	hs_map[1] = CTrapezoidalMap(ThermoTables_Enthalpy[1], ThermoTables_Entropy[1],
			Table_Zone_Edges[1]);

	if (rank == 12201) {
		cout << "Building trapezoidal map for Ps..." << endl;
	}
	Ps_map[0] = CTrapezoidalMap(ThermoTables_Pressure[0], ThermoTables_Entropy[0],
			Table_Zone_Edges[0]);
	Ps_map[1] = CTrapezoidalMap(ThermoTables_Pressure[1], ThermoTables_Entropy[1],
			Table_Zone_Edges[1]);

	if (rank == 12201) {
		cout << "Building trapezoidal map for rhoT..." << endl;
	}
	rhoT_map[0] = CTrapezoidalMap(ThermoTables_Density[0],
			ThermoTables_Temperature[0], Table_Zone_Edges[0]);
	;
	rhoT_map[1] = CTrapezoidalMap(ThermoTables_Density[1],
			ThermoTables_Temperature[1], Table_Zone_Edges[1]);
	;

	if (rank == 12201) {
		cout << "Building trapezoidal map for PT (in vapor region only)..." << endl;
	}
	PT_map[1] = CTrapezoidalMap(ThermoTables_Pressure[1],
			ThermoTables_Temperature[1], Table_Zone_Edges[1]);
	;
	PT_map[0] = PT_map[1];

}

CLookUpTable::~CLookUpTable(void) {
// Using vectors so no need to deallocate

}

void CLookUpTable::Get_Unique_Edges() {
//Import all potential edges into a vector
	for (int j = 0; j < 2; j++) {
		Table_Zone_Edges[j].resize(3 * nTable_Zone_Triangles[j], vector<int>(3, 0));
		//Fill with edges (based on triangulation
		for (int i = 0; i < nTable_Zone_Triangles[j]; i++) {
			int smaller_point, larger_point;
			smaller_point = Table_Zone_Triangles[j][i][0];
			larger_point = Table_Zone_Triangles[j][i][1];
			Table_Zone_Edges[j][3 * i + 0][0] = min(smaller_point, larger_point);
			Table_Zone_Edges[j][3 * i + 0][1] = max(smaller_point, larger_point);
			Table_Zone_Edges[j][3 * i + 0][2] = i;
			smaller_point = Table_Zone_Triangles[j][i][1];
			larger_point = Table_Zone_Triangles[j][i][2];
			Table_Zone_Edges[j][3 * i + 1][0] = min(smaller_point, larger_point);
			Table_Zone_Edges[j][3 * i + 1][1] = max(smaller_point, larger_point);
			Table_Zone_Edges[j][3 * i + 1][2] = i;
			smaller_point = Table_Zone_Triangles[j][i][2];
			larger_point = Table_Zone_Triangles[j][i][0];
			Table_Zone_Edges[j][3 * i + 2][0] = min(smaller_point, larger_point);
			Table_Zone_Edges[j][3 * i + 2][1] = max(smaller_point, larger_point);
			Table_Zone_Edges[j][3 * i + 2][2] = i;
		}
		//Sort the edges to enable selecting unique entries only
		sort(Table_Zone_Edges[j].begin(), Table_Zone_Edges[j].end());
		//Make the list of edges unique and set connectivities at the same time
		int k_final = 0;
		int k_temp = 0;
		while (k_temp  < Table_Zone_Edges[j].size() - 1) {
			Table_Edge_To_Face_Connectivity[j].push_back(vector<int>(1,0));
			Table_Edge_To_Face_Connectivity[j][k_final][0] = Table_Zone_Edges[j][k_temp][2];
			if ((Table_Zone_Edges[j][k_temp][0] = Table_Zone_Edges[j][k_temp + 1][0])
					and (Table_Zone_Edges[j][k_temp][1] = Table_Zone_Edges[j][k_temp + 1][1])) {
						Table_Edge_To_Face_Connectivity[j][k_final].push_back(Table_Zone_Edges[j][k_temp+1][2]);
						k_temp++;
			}
			k_temp++;
			k_final++;
		}
		for (int i = 0; i < nTable_Zone_Triangles[j]; i++) {
			Table_Zone_Edges[j][3 * i + 0].erase(Table_Zone_Edges[j][3 * i + 0].begin()+2);
			Table_Zone_Edges[j][3 * i + 1].erase(Table_Zone_Edges[j][3 * i + 1].begin()+2);
			Table_Zone_Edges[j][3 * i + 2].erase(Table_Zone_Edges[j][3 * i + 2].begin()+2);
		}


		vector<vector<int> >::iterator iter;
		iter = unique(Table_Zone_Edges[j].begin(), Table_Zone_Edges[j].end());
		Table_Zone_Edges[j].resize(distance(Table_Zone_Edges[j].begin(), iter));
	}

//Filter out all the edges which have been imported twice
}


void CLookUpTable::Get_Current_Points_From_TrapezoidalMap(
		CTrapezoidalMap *t_map, su2double x, su2double y) {
	CurrentPoints.resize(6, 0);
	t_map[CurrentZone].Find_Containing_Simplex(x, y);

	CurrentPoints[0] =
			Table_Zone_Edges[CurrentZone][t_map[CurrentZone].getLowerEdge()][0];
	CurrentPoints[1] =
			Table_Zone_Edges[CurrentZone][t_map[CurrentZone].getLowerEdge()][1];
	CurrentPoints[2] =
			Table_Zone_Edges[CurrentZone][t_map[CurrentZone].getMiddleEdge()][0];
	CurrentPoints[3] =
			Table_Zone_Edges[CurrentZone][t_map[CurrentZone].getMiddleEdge()][1];
	CurrentPoints[4] =
			Table_Zone_Edges[CurrentZone][t_map[CurrentZone].getUpperEdge()][0];
	CurrentPoints[5] =
			Table_Zone_Edges[CurrentZone][t_map[CurrentZone].getUpperEdge()][1];
}

void CLookUpTable::SetTDState_rhoe(su2double rho, su2double e) {
// Check if inputs are in total range (necessary but not sufficient condition)

	Get_Current_Points_From_TrapezoidalMap(rhoe_map, rho, e);

	//Now use the quadrilateral which contains the point to interpolate
	//Determine the interpolation coefficients
	Interpolate_2D_Bilinear(rho, e, ThermoTables_Density,
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

}

void CLookUpTable::SetTDState_PT(su2double P, su2double T) {
// Check if inputs are in total range (necessary but not sufficient condition)

	// Check if inputs are in total range (necessary but not sufficient condition)
	Get_Current_Points_From_TrapezoidalMap(PT_map, P, T);
	//Determine interpolation coefficients
	Interpolate_2D_Bilinear(P, T, ThermoTables_Pressure, ThermoTables_Temperature,
			"PT");
	//Interpolate the fluid properties
	Pressure = P;
	Density = Interpolate_2D_Bilinear(ThermoTables_Density);
	StaticEnergy = Interpolate_2D_Bilinear(ThermoTables_StaticEnergy);
	Enthalpy = Interpolate_2D_Bilinear(ThermoTables_Enthalpy);
	Entropy = Interpolate_2D_Bilinear(ThermoTables_Entropy);
	SoundSpeed2 = Interpolate_2D_Bilinear(ThermoTables_SoundSpeed2);
	Temperature = T;
	dPdrho_e = Interpolate_2D_Bilinear(ThermoTables_dPdrho_e);
	dPde_rho = Interpolate_2D_Bilinear(ThermoTables_dPde_rho);
	dTdrho_e = Interpolate_2D_Bilinear(ThermoTables_dTdrho_e);
	dTde_rho = Interpolate_2D_Bilinear(ThermoTables_dTde_rho);
	Cp = Interpolate_2D_Bilinear(ThermoTables_Cp);
	Mu = Interpolate_2D_Bilinear(ThermoTables_Mu);
	Kt = Interpolate_2D_Bilinear(ThermoTables_Kt);

}

void CLookUpTable::SetTDState_Prho(su2double P, su2double rho) {

	Get_Current_Points_From_TrapezoidalMap(Prho_map, P, rho);

//Determine interpolation coefficients
	Interpolate_2D_Bilinear(rho, P, ThermoTables_Density, ThermoTables_Pressure,
			"PRHO");
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

}

void CLookUpTable::SetEnergy_Prho(su2double P, su2double rho) {

	Get_Current_Points_From_TrapezoidalMap(Prho_map, P, rho);

//Determine interpolation coefficients
	Interpolate_2D_Bilinear(rho, P, ThermoTables_Density, ThermoTables_Pressure,
			"PRHO");
	StaticEnergy = Interpolate_2D_Bilinear(ThermoTables_StaticEnergy);
	Pressure = P;
	Density = rho;

}

void CLookUpTable::SetTDState_hs(su2double h, su2double s) {

	Get_Current_Points_From_TrapezoidalMap(hs_map, h, s);

//Determine interpolation coefficients
	Interpolate_2D_Bilinear(h, s, ThermoTables_Enthalpy, ThermoTables_Entropy,
			"HS");

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

}

void CLookUpTable::SetTDState_Ps(su2double P, su2double s) {

	Get_Current_Points_From_TrapezoidalMap(Ps_map, P, s);

//Determine interpolation coefficients
	Interpolate_2D_Bilinear(P, s, ThermoTables_Pressure, ThermoTables_Entropy,
			"PS");

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

}

void CLookUpTable::SetTDState_rhoT(su2double rho, su2double T) {

	Get_Current_Points_From_TrapezoidalMap(rhoT_map, rho, T);
//Determine the interpolation coefficients
	Interpolate_2D_Bilinear(rho, T, ThermoTables_Density,
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

}

inline void CLookUpTable::Gaussian_Inverse(int nDim) {
	//A temporary matrix to hold the inverse
	vector<vector<su2double> > temp;
	temp.resize(nDim, vector<su2double>(2 * nDim, 0));

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
	return;
}

void CLookUpTable::Interpolate_2D_Bilinear(su2double x, su2double y,
		vector<su2double> *ThermoTables_X, vector<su2double> *ThermoTables_Y,
		std::string grid_var) {
	//The x,y coordinates of the quadrilateral
	su2double x0, y0, x1, x2, y1, y2, x3, y3;
	sort(CurrentPoints.begin(), CurrentPoints.end());
	unique(CurrentPoints.begin(), CurrentPoints.end());

	x0 = ThermoTables_X[CurrentZone][CurrentPoints[0]];
	y0 = ThermoTables_Y[CurrentZone][CurrentPoints[0]];
	x1 = ThermoTables_X[CurrentZone][CurrentPoints[1]];
	y1 = ThermoTables_Y[CurrentZone][CurrentPoints[1]];
	x2 = ThermoTables_X[CurrentZone][CurrentPoints[2]];
	y2 = ThermoTables_Y[CurrentZone][CurrentPoints[2]];
	x3 = ThermoTables_X[CurrentZone][CurrentPoints[3]];
	y3 = ThermoTables_Y[CurrentZone][CurrentPoints[3]];

	//cout<<CurrentPoints[0]<<" "<<CurrentPoints[1]<<" "<<CurrentPoints[2]<<" "<<CurrentPoints[3]<<endl;
	//cout<<x0<<" "<<x1<<" "<<x2<<" "<<x3<<endl;
	//cout<<y0<<" "<<y1<<" "<<y2<<" "<<y3<<endl;

	//Setup the LHM matrix for the interpolation (Vandermonde)

	Interpolation_Matrix[0][0] = 1;
	Interpolation_Matrix[0][1] = 0;
	Interpolation_Matrix[0][2] = 0;
	Interpolation_Matrix[0][3] = 0;

	Interpolation_Matrix[1][0] = 1;
	Interpolation_Matrix[1][1] = (x1 - x0);
	Interpolation_Matrix[1][2] = (y1 - y0);
	Interpolation_Matrix[1][3] = (y1 - y0) * (x1 - x0);

	Interpolation_Matrix[2][0] = 1;
	Interpolation_Matrix[2][1] = (x2 - x0);
	Interpolation_Matrix[2][2] = (y2 - y0);
	Interpolation_Matrix[2][3] = (y2 - y0) * (x2 - x0);

	Interpolation_Matrix[3][0] = 1;
	Interpolation_Matrix[3][1] = (x3 - x0);
	Interpolation_Matrix[3][2] = (y3 - y0);
	Interpolation_Matrix[3][3] = (x3 - x0) * (y3 - y0);
//	  Interpolation_Matrix[0][0] = 1;
//		Interpolation_Matrix[0][1] = x0;
//		Interpolation_Matrix[0][2] = y0;
//		Interpolation_Matrix[0][3] = x0*y0;
//
//		Interpolation_Matrix[1][0] = 1;
//		Interpolation_Matrix[1][1] = (x1);
//		Interpolation_Matrix[1][2] = (y1);
//		Interpolation_Matrix[1][3] = (y1) * (x1);
//
//		Interpolation_Matrix[2][0] = 1;
//		Interpolation_Matrix[2][1] = (x2);
//		Interpolation_Matrix[2][2] = (y2);
//		Interpolation_Matrix[2][3] = (y2)* (x2);
//
//		Interpolation_Matrix[3][0] = 1;
//		Interpolation_Matrix[3][1] = (x3);
//		Interpolation_Matrix[3][2] = (y3);
//		Interpolation_Matrix[3][3] = (x3) * (y3);

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
		d = d + Interpolation_Coeff[i][1] * (x - x0);
		d = d + Interpolation_Coeff[i][2] * (y - y0);
		d = d + Interpolation_Coeff[i][3] * (x - x0) * (y - y0);
		Interpolation_Coeff[i][0] = d;
	}
}

su2double CLookUpTable::Interpolate_2D_Bilinear(
		vector<su2double> *ThermoTables_Z) {
	//The function values at the 4 corners of the quad
	su2double func_value_0, func_value_1, func_value_2, func_value_3;

	func_value_0 = ThermoTables_Z[CurrentZone][CurrentPoints[0]];
	func_value_1 = ThermoTables_Z[CurrentZone][CurrentPoints[1]];
	func_value_2 = ThermoTables_Z[CurrentZone][CurrentPoints[2]];
	func_value_3 = ThermoTables_Z[CurrentZone][CurrentPoints[3]];

	//The Interpolation_Coeff values depend on location alone
	//and are the same regardless of function values
	su2double result = 0;
	result = result + Interpolation_Coeff[0][0] * func_value_0;
	result = result + Interpolation_Coeff[1][0] * func_value_1;
	result = result + Interpolation_Coeff[2][0] * func_value_2;
	result = result + Interpolation_Coeff[3][0] * func_value_3;

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
	//for (int i = 0; i < 2; i++) {
	int i = CurrentZone;
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
	//}

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
				//Triangles in .tec file are indexed from 1
				//In cpp it is more convenient to start with 0.
				Table_Zone_Triangles[zone_scanned][j][0]--;
				Table_Zone_Triangles[zone_scanned][j][1]--;
				Table_Zone_Triangles[zone_scanned][j][2]--;
			}

			zone_scanned++;
		}
	}

	table.close();
//NonDimensionalise
	NonDimensionalise_Table_Values();
}

void CLookUpTable::LookUpTable_Malloc(int Index_of_Zone) {
	ThermoTables_StaticEnergy[Index_of_Zone] = vector< su2double >(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_Entropy[Index_of_Zone] = vector< su2double >(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_Enthalpy[Index_of_Zone] = vector< su2double >(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_Density[Index_of_Zone] = vector< su2double >(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_Pressure[Index_of_Zone] = vector< su2double >(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_SoundSpeed2[Index_of_Zone] = vector< su2double >(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_Temperature[Index_of_Zone] = vector< su2double >(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_dPdrho_e[Index_of_Zone] = vector< su2double >(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_dPde_rho[Index_of_Zone] = vector< su2double >(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_dTdrho_e[Index_of_Zone] = vector< su2double >(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_dTde_rho[Index_of_Zone] = vector< su2double >(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_Cp[Index_of_Zone] = vector< su2double >(
			nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_Mu[Index_of_Zone] = vector< su2double >(
			nTable_Zone_Stations[Index_of_Zone], 0);
	//ThermoTables_dmudrho_T[Index_of_Zone] = vector< su2double >(
	//	nTable_Zone_Stations[Index_of_Zone], 0);
	//ThermoTables_dmudT_rho[Index_of_Zone] = vector< su2double >(
	//	nTable_Zone_Stations[Index_of_Zone], 0);
	ThermoTables_Kt[Index_of_Zone] = vector< su2double >(
			nTable_Zone_Stations[Index_of_Zone], 0);
	//ThermoTables_dktdrho_T[Index_of_Zone] = vector< su2double >(
	//		nTable_Zone_Stations[Index_of_Zone], 0);
	//ThermoTables_dktdT_rho[Index_of_Zone] = vector< su2double >(
	//		nTable_Zone_Stations[Index_of_Zone], 0);
	Table_Zone_Triangles[Index_of_Zone] = vector<vector<int> >(
			nTable_Zone_Triangles[Index_of_Zone]);
	for (int j = 0; j < nTable_Zone_Triangles[Index_of_Zone]; j++) {
		Table_Zone_Triangles[Index_of_Zone][j] = vector<int>(3, 0);
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

