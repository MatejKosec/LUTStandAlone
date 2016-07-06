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

CThermoList::CThermoList() {

	StaticEnergy = 0.0;
	Entropy = 0.0;
	Enthalpy = 0.0;
	Density = 0.0;
	Pressure = 0.0;
	SoundSpeed2 = 0.0;
	Temperature = 0.0;
	dPdrho_e = 0.0;
	dPde_rho = 0.0;
	dTdrho_e = 0.0;
	dTde_rho = 0.0;
	Cp = 0.0;
	Mu = 0.0;
	dmudrho_T = 0.0;
	dmudT_rho = 0.0;
	Kt = 0.0;
	dktdrho_T = 0.0;
	dktdT_rho = 0.0;

}

CThermoList::~CThermoList() {

}

//void CThermoList::CTLprint()
//{
//	cout<<"StaticEnergy:"<<StaticEnergy<<endl;
//	cout<<"Enthalpy    :"<<Enthalpy<<endl;
//	cout<<"Entropy     :"<<Entropy<<endl;
//	cout<<"Density     :"<<Density<<endl;
//	cout<<"Pressure    :"<<Pressure<<endl;
//	cout<<"SoundSpeed2 :"<<SoundSpeed2<<endl;
//	cout<<"Temperature :"<<Temperature<<endl;
//	cout<<"dPdrho_e    :"<<dPdrho_e<<endl;
//	cout<<"dPde_rho    :"<<dPde_rho<<endl;
//	cout<<"dTdrho_e    :"<<dTdrho_e<<endl;
//	cout<<"dTde_rho    :"<<dTde_rho<<endl;
//	cout<<"Cp          :"<<Cp<<endl;
//	cout<<"Mu          :"<<Mu<<endl;
//	cout<<"dmudrho_T   :"<<dmudrho_T<<endl;
//	cout<<"dmudT_rho   :"<<dmudT_rho<<endl;
//	cout<<"Kt          :"<<Kt<<endl;
//	cout<<"dktdrho_T   :"<<dktdrho_T<<endl;
//	cout<<"dktdT_rho   :"<<dktdT_rho<<endl;
//}

CLookUpTable::CLookUpTable() {

	ThermoTables = NULL;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Interpolation_Coeff[i][j] = 0.0;
		}
	}
	HS_tree = NULL;
	iIndex = -1;
	jIndex = -1;
	Table_Pressure_Stations = 0;
	Table_Density_Stations = 0;
}

CLookUpTable::CLookUpTable(char* Filename) {
	ThermoTables = NULL;

	TableLoadCFX(Filename);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Interpolation_Coeff[i][j] = -1.0;
		}
	}
	iIndex = -1; //negative number means it hasn't been preset yet
	jIndex = -1; //same


	cout << "Table_Pressure_Stations  : " << Table_Pressure_Stations << endl;
	cout << "Table_Density_Stations: " << Table_Density_Stations << endl;
	cout << "Building HS_tree" << endl;
	su2double* xtemp = new su2double[Table_Density_Stations
																	 * Table_Pressure_Stations];
	su2double* ytemp = new su2double[Table_Density_Stations
																	 * Table_Pressure_Stations];
	int* itemp = new int[Table_Density_Stations * Table_Pressure_Stations];
	for (int i = 0; i < Table_Density_Stations; i++) {
		for (int j = 0; j < Table_Pressure_Stations; j++) {
			xtemp[Table_Pressure_Stations * i + j] = ThermoTables[i][j].Enthalpy;
			ytemp[Table_Pressure_Stations * i + j] = ThermoTables[i][j].Entropy;
			itemp[Table_Pressure_Stations * i + j] = Table_Pressure_Stations * i + j;
		}
	}
	HS_tree = KD_Tree(xtemp, ytemp, itemp,
			Table_Pressure_Stations * Table_Density_Stations, 0);
	cout << "HS_tree built" << endl;
}

CLookUpTable::~CLookUpTable(void) {
	delete (ThermoTables);
	free_KD_tree(HS_tree);
}

void CLookUpTable::free_KD_tree(KD_node* root) {
	if (root->Branch_Dimension > 1) {
		free_KD_tree(root->upper);
		free_KD_tree(root->lower);
	}
	su2double* test = root->x_values;
	delete root;
}

struct KD_node* CLookUpTable::KD_Tree(su2double* x_values, su2double* y_values,
		int* Flattened_Point_Index, int dim, int depth) {

	struct KD_node *kdn = new KD_node;	//(KD_node*) malloc(sizeof(KD_node));;
	kdn->x_values = x_values;
	kdn->y_values = y_values;
	kdn->Flattened_Point_Index = Flattened_Point_Index;
	kdn->Branch_Splitting_Direction = depth;
	kdn->Branch_Dimension = dim;
	if (dim > 1) {
		su2double temp;
		int itemp = 0;
		int swaps = 0; //for bubblesort
		bool sorted = false;

		if (depth % 2 == 0) {
			//bubble sort by xvalues
			while (not sorted) {
				swaps = 0;

				for (int i = 0; i < dim - 1; i++) {
					if (x_values[i] > x_values[i + 1]) {
						temp = x_values[i];
						x_values[i] = x_values[i + 1];
						x_values[i + 1] = temp;

						temp = y_values[i];
						y_values[i] = y_values[i + 1];
						y_values[i + 1] = temp;

						itemp = Flattened_Point_Index[i];
						Flattened_Point_Index[i] = Flattened_Point_Index[i + 1];
						Flattened_Point_Index[i + 1] = itemp;
						//keep a record of the number of swaps performed
						swaps++;
					}
				}
				if (swaps == 0)
					sorted = true;
			}
		} else if (depth % 2 == 1) {
			//bubble sort by yvalues
			while (not sorted) {
				swaps = 0;

				for (int i = 0; i < dim - 1; i++) {
					if (y_values[i] > y_values[i + 1]) {
						temp = y_values[i];
						y_values[i] = y_values[i + 1];
						y_values[i + 1] = temp;

						temp = x_values[i];
						x_values[i] = x_values[i + 1];
						x_values[i + 1] = temp;

						itemp = Flattened_Point_Index[i];
						Flattened_Point_Index[i] = Flattened_Point_Index[i + 1];
						Flattened_Point_Index[i + 1] = itemp;
						//keep a record of the number of swaps performed
						swaps++;
					}
				}
				if (swaps == 0)
					sorted = true;
			}
		}
		//Create the new upper and lower arrays
		su2double* upperx = new su2double[dim / 2];
		su2double* uppery = new su2double[dim / 2];
		int* upperi = new int[dim / 2];
		su2double* lowerx = new su2double[dim - dim / 2];
		su2double* lowery = new su2double[dim - dim / 2];
		int* loweri = new int[dim - dim / 2];
		for (int i = dim / 2; i < dim; i++) {
			upperx[i - dim / 2] = x_values[i];
			uppery[i - dim / 2] = y_values[i];
			upperi[i - dim / 2] = Flattened_Point_Index[i];
		}
		for (int i = 0; i < dim / 2; i++) {
			lowerx[i] = x_values[i];
			lowery[i] = y_values[i];
			loweri[i] = Flattened_Point_Index[i];
		}

		kdn->upper = KD_Tree(upperx, uppery, upperi, dim / 2, depth + 1);
		kdn->lower = KD_Tree(lowerx, lowery, loweri, dim - dim / 2, depth + 1);
	}

	return kdn;
}

su2double CLookUpTable::Dist_KD_Tree(su2double x, su2double y,
		KD_node *branch) {
	su2double dist;
	dist = pow((branch->x_values[branch->Branch_Dimension / 2] - x) / x, 2)\

			+ pow((branch->y_values[branch->Branch_Dimension / 2] - y) / y, 2);
	return dist;
}

void CLookUpTable::NN_N_KD_Tree(int N, su2double thermo1, su2double thermo2,
		KD_node *root, su2double *best_dist) {
	//Recursive search for the nearest neighbors of the KD tree, with unwinding
	su2double dist;
	dist = Dist_KD_Tree(thermo1, thermo2, root);
	//cout<<"Depth "<<root->Branch_Splitting_Direction;
	int i = 0;
	while (i < N) {
		if (dist == best_dist[i])
			i = N+1;
		if (dist < best_dist[i]) {
			//cout<<" i:"<<i;
			for (int j = N-1; j > i; j--) {
				//cout<<" j: "<<j;
				best_dist[j] = best_dist[j - 1];
				Nearest_Neighbour_iIndex[j] = Nearest_Neighbour_iIndex[j - 1];
				Nearest_Neighbour_jIndex[j] = Nearest_Neighbour_jIndex[j - 1];
			}
			best_dist[i] = dist;
			Nearest_Neighbour_iIndex[i] =
					root->Flattened_Point_Index[root->Branch_Dimension / 2]
																			/ Table_Pressure_Stations;
			Nearest_Neighbour_jIndex[i] =
					root->Flattened_Point_Index[root->Branch_Dimension / 2]
																			% Table_Pressure_Stations;
			i = N+1;
		}
		i++;

	}

	if ((root->Branch_Dimension > 1)) {
		if (root->Branch_Splitting_Direction % 2 == 0) {
			if (root->x_values[root->Branch_Dimension / 2] <= thermo1) {
				NN_N_KD_Tree(N,thermo1, thermo2, root->upper, best_dist);
				if (dist < best_dist[N-1]) {

					NN_N_KD_Tree(N, thermo1, thermo2, root->lower, best_dist);
				}
			} else if (root->x_values[root->Branch_Dimension / 2] > thermo1) {
				NN_N_KD_Tree(N,thermo1, thermo2, root->lower, best_dist);
				if (dist < best_dist[N-1]) {
					NN_N_KD_Tree(N, thermo1, thermo2, root->upper, best_dist);
				}
			}
		} else if (root->Branch_Splitting_Direction % 2 == 1) {
			if (root->y_values[root->Branch_Dimension / 2] <= thermo2) {
				NN_N_KD_Tree(N,thermo1, thermo2, root->upper, best_dist);
				if (dist < best_dist[N-1]) {
					NN_N_KD_Tree(N,thermo1, thermo2, root->lower, best_dist);
				}
			} else if (root->y_values[root->Branch_Dimension / 2] > thermo2) {
				NN_N_KD_Tree(N,thermo1, thermo2, root->lower, best_dist);
				if (dist < best_dist[N-1]) {
					NN_N_KD_Tree(N,thermo1, thermo2, root->upper, best_dist);
				}
			}
		}
	}
}


void CLookUpTable::SetTDState_rhoe(su2double rho, su2double e) {

	//Check if inputs are in total range (necessary but not sufficient condition)
	if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
		cerr << "RHOE Input Density out of bounds\n";
	}
	if ((e > StaticEnergy_Table_Limits[1])
			or (e < StaticEnergy_Table_Limits[0])) {
		cerr << "RHOE Input StaticEnergy out of bounds\n";
	}
	Nearest_Neighbour_iIndex = new int[4];
	Nearest_Neighbour_jIndex = new int[4];
	su2double RunVal;
	unsigned int LowerI, UpperI, LowerJ, UpperJ, middleI, middleJ;

	UpperJ = Table_Pressure_Stations - 1;
	LowerJ = 0;
	su2double grad, x00, y00, y10, x10;

	//Determine the I index: rho is not equispaced
	UpperI = Table_Density_Stations - 1;
	LowerI = 0;

	while (UpperI - LowerI > 1) {
		middleI = (UpperI + LowerI) / 2;
		//Check current value
		x00 = ThermoTables[middleI][LowerJ].Density;
		grad = ThermoTables[middleI + 1][LowerJ].Density - x00;
		if (x00 * grad > rho * grad) {
			UpperI = middleI;
		} else if (x00 < rho) {
			LowerI = middleI;
		}
		if (x00 == rho) {
			UpperI = LowerI + 1;
			break;
		}
	}

	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2;
		//Check current value
		y00 = ThermoTables[LowerI][middleJ].StaticEnergy;
		y10 = ThermoTables[UpperI][middleJ].StaticEnergy;
		x00 = ThermoTables[LowerI][middleJ].Density;
		x10 = ThermoTables[UpperI][middleJ].Density;
		RunVal = y00 + (y10 - y00) / (x10 - x00) * (rho - x00);
		grad = ThermoTables[LowerI][middleJ + 1].StaticEnergy - y00;
		if (RunVal * grad > e * grad) {
			UpperJ = middleJ;
		} else if (RunVal * grad < e * grad) {
			LowerJ = middleJ;
		}
		if (RunVal == e) {
			UpperJ = LowerJ + 1;
			break;
		}
	}

	iIndex = LowerI;
	jIndex = LowerJ;

	//Now use the closest fit box to interpolate
	su2double x, y;
	x = (rho - ThermoTables[iIndex][jIndex].Density) / Density_Table_Limits[1];
	y = (e - ThermoTables[iIndex][jIndex].StaticEnergy)
							/ StaticEnergy_Table_Limits[1];
	//Set the nearest neigbours to the adjacent i and j vertexes
	Nearest_Neighbour_iIndex[0] = iIndex;
	Nearest_Neighbour_iIndex[1] = iIndex + 1;
	Nearest_Neighbour_iIndex[2] = iIndex;
	Nearest_Neighbour_iIndex[3] = iIndex + 1;
	Nearest_Neighbour_jIndex[0] = jIndex;
	Nearest_Neighbour_jIndex[1] = jIndex;
	Nearest_Neighbour_jIndex[2] = jIndex + 1;
	Nearest_Neighbour_jIndex[3] = jIndex + 1;
	//Determine interpolation coefficients
	Interp2D_ArbitrarySkewCoeff(x, y, "RHOE");

	StaticEnergy = e;
	Density = rho;
	Entropy = Interp2D_lin(x, y, "Entropy");
	Pressure = Interp2D_lin(x, y, "Pressure");
	Enthalpy = Interp2D_lin(x, y, "Enthalpy");
	SoundSpeed2 = Interp2D_lin(x, y, "SoundSpeed2");
	Temperature = Interp2D_lin(x, y, "Temperature");
	dPdrho_e = Interp2D_lin(x, y, "dPdrho_e");
	dPde_rho = Interp2D_lin(x, y, "dPde_rho");
	dTdrho_e = Interp2D_lin(x, y, "dTdrho_e");
	dTde_rho = Interp2D_lin(x, y, "dTde_rho");
	Cp = Interp2D_lin(x, y, "Cp");
	Mu = Interp2D_lin(x, y, "Mu");
	dmudrho_T = Interp2D_lin(x, y, "dmudrho_T");
	dmudT_rho = Interp2D_lin(x, y, "dmudT_rho");
	Kt = Interp2D_lin(x, y, "Kt");
	dktdrho_T = Interp2D_lin(x, y, "dktdrho_T");
	dktdT_rho = Interp2D_lin(x, y, "dktdT_rho");

	if ((Density > Density_Table_Limits[1])
			or (Density < Density_Table_Limits[0])) {
		cerr << "RHOE Interpolated Density out of bounds\n";
	}
	if ((Pressure > Pressure_Table_Limits[1])
			or (Pressure < Pressure_Table_Limits[0])) {
		cerr << "RHOE Interpolated Pressure out of bounds\n";
	}
	delete [] Nearest_Neighbour_iIndex;
	delete [] Nearest_Neighbour_jIndex;
}

void CLookUpTable::SetTDState_PT(su2double P, su2double T) {
	//Check if inputs are in total range (necessary but not sufficient condition)
	if ((P > Pressure_Table_Limits[1]) or (P < Pressure_Table_Limits[0])) {
		cerr << "PT Input Pressure out of bounds\n";
	}
	if ((T > Temperature_Table_Limits[1]) or (T < Temperature_Table_Limits[0])) {
		cerr << "PT Input Temperature out of bounds\n";
	}

	Nearest_Neighbour_iIndex = new int[4];
	Nearest_Neighbour_jIndex = new int[4];
	unsigned int LowerI, UpperI, LowerJ, UpperJ, middleI, middleJ;

	UpperJ = Table_Pressure_Stations - 1;
	LowerJ = 0;
	su2double grad, x00, y00, y01, x01, RunVal;

	//Determine the I index: rho is not equispaced
	UpperI = Table_Density_Stations - 1;
	LowerI = 0;

	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2;
		//Check current value
		x00 = ThermoTables[LowerI][middleJ].Pressure;
		grad = ThermoTables[LowerI][middleJ + 1].Pressure - x00;
		if (x00 * grad > P * grad) {
			UpperJ = middleJ;
		} else if (x00 < P) {
			LowerJ = middleJ;
		}
		if (x00 == P) {
			UpperJ = LowerJ + 1;
			break;
		}
	}

	//Determine the I index (for T)
	while (UpperI - LowerI > 1) {
		middleI = (UpperI + LowerI) / 2;
		//Check current value
		y00 = ThermoTables[middleI][LowerJ].Pressure;
		y01 = ThermoTables[middleI][UpperJ].Pressure;
		x00 = ThermoTables[middleI][LowerJ].Temperature;
		x01 = ThermoTables[middleI][UpperJ].Temperature;
		grad = ThermoTables[UpperI][LowerJ].Temperature - x00;
		RunVal = x00 + (x01 - x00) / (y01 - y00) * (P - y00);
		if (RunVal * grad > T * grad) {
			UpperI = middleI;
		} else if (RunVal * grad < T * grad) {
			LowerI = middleI;
		}
		if (RunVal == T) {
			UpperI = LowerI + 1;
			break;
		}
	}

	iIndex = LowerI;
	jIndex = LowerJ;

	su2double x, y;
	x = (T - ThermoTables[iIndex][jIndex].Temperature)
							/ Temperature_Table_Limits[1];
	y = (P - ThermoTables[iIndex][jIndex].Pressure) / Pressure_Table_Limits[1];
	//Set the nearest neigbours to the adjacent i and j vertexes
	Nearest_Neighbour_iIndex[0] = iIndex;
	Nearest_Neighbour_iIndex[1] = iIndex + 1;
	Nearest_Neighbour_iIndex[2] = iIndex;
	Nearest_Neighbour_iIndex[3] = iIndex + 1;
	Nearest_Neighbour_jIndex[0] = jIndex;
	Nearest_Neighbour_jIndex[1] = jIndex;
	Nearest_Neighbour_jIndex[2] = jIndex + 1;
	Nearest_Neighbour_jIndex[3] = jIndex + 1;
	//Determine interpolation coefficients
	Interp2D_ArbitrarySkewCoeff(x, y, "PT");

	Temperature = T;
	Pressure = P;
	StaticEnergy = Interp2D_lin(x, y, "StaticEnergy");
	Enthalpy = Interp2D_lin(x, y, "Enthalpy");
	Entropy = Interp2D_lin(x, y, "Entropy");
	Density = Interp2D_lin(x, y, "Density");
	SoundSpeed2 = Interp2D_lin(x, y, "SoundSpeed2");
	dPdrho_e = Interp2D_lin(x, y, "dPdrho_e");
	dPde_rho = Interp2D_lin(x, y, "dPde_rho");
	dTdrho_e = Interp2D_lin(x, y, "dTdrho_e");
	dTde_rho = Interp2D_lin(x, y, "dTde_rho");
	Cp = Interp2D_lin(x, y, "Cp");
	Mu = Interp2D_lin(x, y, "Mu");
	dmudrho_T = Interp2D_lin(x, y, "dmudrho_T");
	dmudT_rho = Interp2D_lin(x, y, "dmudT_rho");
	Kt = Interp2D_lin(x, y, "Kt");
	dktdrho_T = Interp2D_lin(x, y, "dktdrho_T");
	dktdT_rho = Interp2D_lin(x, y, "dktdT_rho");

	if ((Density > Density_Table_Limits[1])
			or (Density < Density_Table_Limits[0])) {
		cerr << "PT Interpolated Density out of bounds\n";
	}
	if ((Pressure > Pressure_Table_Limits[1])
			or (Pressure < Pressure_Table_Limits[0])) {
		cerr << "PT Interpolated Pressure out of bounds\n";
	}
	delete [] Nearest_Neighbour_iIndex;
	delete [] Nearest_Neighbour_jIndex;
}

void CLookUpTable::SetTDState_Prho(su2double P, su2double rho) {
	//Check if inputs are in total range (necessary but not sufficient condition)
	if ((P > Pressure_Table_Limits[1]) or (P < Pressure_Table_Limits[0])) {
		cerr << "PRHO Input Pressure out of bounds\n";
	}
	if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
		cerr << "PRHO Input Density out of bounds\n";
	}

	Nearest_Neighbour_iIndex = new int[4];
	Nearest_Neighbour_jIndex = new int[4];
	unsigned int LowerI, UpperI, LowerJ, UpperJ, middleI, middleJ;

	UpperJ = Table_Pressure_Stations - 1;
	LowerJ = 0;
	su2double grad, x00, y00;

	//Determine the I index: rho is not equispaced
	UpperI = Table_Density_Stations - 1;
	LowerI = 0;

	//Determine the I index: rho is not equispaced
	UpperI = Table_Density_Stations - 1;
	LowerI = 0;

	while (UpperI - LowerI > 1) {
		middleI = (UpperI + LowerI) / 2;
		//Check current value
		x00 = ThermoTables[middleI][LowerJ].Density;
		grad = ThermoTables[middleI + 1][LowerJ].Density - x00;
		if (x00 * grad > rho * grad) {
			UpperI = middleI;
		} else if (x00 < rho) {
			LowerI = middleI;
		}
		if (x00 == rho) {
			UpperI = LowerI + 1;
			break;
		}
	}

	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2;
		//Check current value
		y00 = ThermoTables[LowerI][middleJ].Pressure;
		grad = ThermoTables[LowerI][middleJ + 1].Pressure - y00;
		if (y00 * grad > P * grad) {
			UpperJ = middleJ;
		} else if (y00 < P) {
			LowerJ = middleJ;
		}
		if (y00 == P) {
			UpperJ = LowerJ + 1;
			break;
		}
	}

	iIndex = LowerI;
	jIndex = LowerJ;

	su2double x, y;
	x = (rho - ThermoTables[iIndex][jIndex].Density) / Density_Table_Limits[1];
	y = (P - ThermoTables[iIndex][jIndex].Pressure) / Pressure_Table_Limits[1];
	//Set the nearest neigbours to the adjacent i and j vertexes
	Nearest_Neighbour_iIndex[0] = iIndex;
	Nearest_Neighbour_iIndex[1] = iIndex + 1;
	Nearest_Neighbour_iIndex[2] = iIndex;
	Nearest_Neighbour_iIndex[3] = iIndex + 1;
	Nearest_Neighbour_jIndex[0] = jIndex;
	Nearest_Neighbour_jIndex[1] = jIndex;
	Nearest_Neighbour_jIndex[2] = jIndex + 1;
	Nearest_Neighbour_jIndex[3] = jIndex + 1;
	//Determine interpolation coefficients
	Interp2D_ArbitrarySkewCoeff(x, y, "PRHO");

	Pressure = P;
	Density = rho;
	StaticEnergy = Interp2D_lin(x, y, "StaticEnergy");
	Enthalpy = Interp2D_lin(x, y, "Enthalpy");
	Entropy = Interp2D_lin(x, y, "Entropy");
	SoundSpeed2 = Interp2D_lin(x, y, "SoundSpeed2");
	Temperature = Interp2D_lin(x, y, "Temperature");
	dPdrho_e = Interp2D_lin(x, y, "dPdrho_e");
	dPde_rho = Interp2D_lin(x, y, "dPde_rho");
	dTdrho_e = Interp2D_lin(x, y, "dTdrho_e");
	dTde_rho = Interp2D_lin(x, y, "dTde_rho");
	Cp = Interp2D_lin(x, y, "Cp");
	Mu = Interp2D_lin(x, y, "Mu");
	dmudrho_T = Interp2D_lin(x, y, "dmudrho_T");
	dmudT_rho = Interp2D_lin(x, y, "dmudT_rho");
	Kt = Interp2D_lin(x, y, "Kt");
	dktdrho_T = Interp2D_lin(x, y, "dktdrho_T");
	dktdT_rho = Interp2D_lin(x, y, "dktdT_rho");

	if ((Density > Density_Table_Limits[1])
			or (Density < Density_Table_Limits[0])) {
		cerr << "PRHO Interpolated Density out of bounds\n";
	}
	if ((Pressure > Pressure_Table_Limits[1])
			or (Pressure < Pressure_Table_Limits[0])) {
		cerr << "PRHO Interpolated Pressure out of bounds\n";
	}
	delete [] Nearest_Neighbour_iIndex;
	delete [] Nearest_Neighbour_jIndex;

}

void CLookUpTable::SetEnergy_Prho(su2double P, su2double rho) {

}

void CLookUpTable::SetTDState_hs(su2double h, su2double s) {
	//Check if inputs are in total range (necessary but not sufficient condition)
	if ((h > Enthalpy_Table_Limits[1]) or (h < Enthalpy_Table_Limits[0])) {
		cerr << "HS Input Enthalpy out of bounds\n";
	}
	if ((s > Entropy_Table_Limits[1]) or (s < Entropy_Table_Limits[0])) {
		cerr << "HS Input Entropy out of bounds\n";
	}

	iIndex = HS_tree->Flattened_Point_Index[HS_tree->Branch_Dimension / 2]
																					/ Table_Pressure_Stations;
	jIndex = HS_tree->Flattened_Point_Index[HS_tree->Branch_Dimension / 2]
																					% Table_Pressure_Stations;
	int N = 10;
	Nearest_Neighbour_iIndex = new int[N];
  Nearest_Neighbour_jIndex = new int[N];
  su2double *best_dist = new su2double[N];
	for (int i = 0; i < N; i++) {
		Nearest_Neighbour_iIndex[i] = -1;
		Nearest_Neighbour_jIndex[i] = -1;
	}

	for (int i = 0; i < N; i++) {
		best_dist[i] = 1E10;
	}

	//Search the HS_tree
	NN_N_KD_Tree(N,h, s, HS_tree, best_dist);

	Entropy = s;
	Enthalpy = h;
	StaticEnergy = Interp2D_Inv_Dist(N,"StaticEnergy", best_dist);
	Density = Interp2D_Inv_Dist(N,"Density", best_dist);
	Pressure = Interp2D_Inv_Dist(N,"Pressure", best_dist);
	SoundSpeed2 = Interp2D_Inv_Dist(N,"SoundSpeed2", best_dist);
	Temperature = Interp2D_Inv_Dist(N,"Temperature", best_dist);
	dPdrho_e = Interp2D_Inv_Dist(N,"dPdrho_e", best_dist);
	dPde_rho = Interp2D_Inv_Dist(N,"dPde_rho", best_dist);
	dTdrho_e = Interp2D_Inv_Dist(N,"dTdrho_e", best_dist);
	dTde_rho = Interp2D_Inv_Dist(N,"dTde_rho", best_dist);
	Cp = Interp2D_Inv_Dist(N,"Cp", best_dist);
	Mu = Interp2D_Inv_Dist(N,"Mu", best_dist);
	dmudrho_T = Interp2D_Inv_Dist(N,"dmudrho_T", best_dist);
	dmudT_rho = Interp2D_Inv_Dist(N,"dmudT_rho", best_dist);
	Kt = Interp2D_Inv_Dist(N,"Kt", best_dist);
	dktdrho_T = Interp2D_Inv_Dist(N,"dktdrho_T", best_dist);
	dktdT_rho = Interp2D_Inv_Dist(N,"dktdT_rho", best_dist);


	if ((Density > Density_Table_Limits[1])
			or (Density < Density_Table_Limits[0])) {
		cerr << "HS Interpolated Density out of bounds\n";
	}
	if ((Pressure > Pressure_Table_Limits[1])
			or (Pressure < Pressure_Table_Limits[0])) {
		cerr << "HS Interpolated Pressure out of bounds\n";
	}
	delete [] best_dist;
	delete [] Nearest_Neighbour_iIndex;
	delete [] Nearest_Neighbour_jIndex;

}

void CLookUpTable::SetTDState_Ps(su2double P, su2double s) {
	//Check if inputs are in total range (necessary but not sufficient condition)
	if ((P > Pressure_Table_Limits[1]) or (P < Pressure_Table_Limits[0])) {
		cerr << "PS Input Pressure out of bounds\n";
	}
	if ((s > Entropy_Table_Limits[1]) or (s < Entropy_Table_Limits[0])) {
		cerr << "PS Input Entropy  out of bounds\n";
	}

	unsigned int LowerI, UpperI, LowerJ, UpperJ, middleI, middleJ;
	Nearest_Neighbour_iIndex = new int[4];
	Nearest_Neighbour_jIndex = new int[4];
	UpperJ = Table_Pressure_Stations - 1;
	LowerJ = 0;
	su2double grad, x00, y00, y01, x01, RunVal;

	//Determine the I index: rho is not equispaced
	UpperI = Table_Density_Stations - 1;
	LowerI = 0;

	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2;
		//Check current value
		y00 = ThermoTables[LowerI][middleJ].Pressure;
		grad = ThermoTables[LowerI][middleJ + 1].Pressure - y00;
		if (y00 * grad > P * grad) {
			UpperJ = middleJ;
		} else if (y00 < P) {
			LowerJ = middleJ;
		}
		if (y00 == P) {
			UpperJ = LowerJ + 1;
			break;
		}
	}

	//Determine the I index (for s)
	while (UpperI - LowerI > 1) {
		middleI = (UpperI + LowerI) / 2;
		//Check current value
		y00 = ThermoTables[middleI][LowerJ].Pressure;
		y01 = ThermoTables[middleI][UpperJ].Pressure;
		x00 = ThermoTables[middleI][LowerJ].Entropy;
		x01 = ThermoTables[middleI][UpperJ].Entropy;
		grad = ThermoTables[UpperI][LowerJ].Entropy - x00;
		RunVal = x00 + (x01 - x00) / (y01 - y00) * (P - y00);
		grad = ThermoTables[middleI + 1][LowerJ].Entropy - x00;
		if (RunVal * grad > s * grad) {
			UpperI = middleI;
		} else if (RunVal * grad < s * grad) {
			LowerI = middleI;
		}
		if (RunVal == s) {
			UpperI = LowerI + 1;
			break;
		}
	}

	iIndex = LowerI;
	jIndex = LowerJ;

	su2double x, y;
	y = (P - ThermoTables[iIndex][jIndex].Pressure) / Pressure_Table_Limits[1];
	x = (s - ThermoTables[iIndex][jIndex].Entropy) / Entropy_Table_Limits[1];
	//Set the nearest neigbours to the adjacent i and j vertexes
	Nearest_Neighbour_iIndex[0] = iIndex;
	Nearest_Neighbour_iIndex[1] = iIndex + 1;
	Nearest_Neighbour_iIndex[2] = iIndex;
	Nearest_Neighbour_iIndex[3] = iIndex + 1;
	Nearest_Neighbour_jIndex[0] = jIndex;
	Nearest_Neighbour_jIndex[1] = jIndex;
	Nearest_Neighbour_jIndex[2] = jIndex + 1;
	Nearest_Neighbour_jIndex[3] = jIndex + 1;
	//Determine interpolation coefficients
	Interp2D_ArbitrarySkewCoeff(x, y, "PS");

	Entropy = s;
	Pressure = P;
	StaticEnergy = Interp2D_lin(x, y, "StaticEnergy");
	Enthalpy = Interp2D_lin(x, y, "Enthalpy");
	Density = Interp2D_lin(x, y, "Density");
	SoundSpeed2 = Interp2D_lin(x, y, "SoundSpeed2");
	Temperature = Interp2D_lin(x, y, "Temperature");
	dPdrho_e = Interp2D_lin(x, y, "dPdrho_e");
	dPde_rho = Interp2D_lin(x, y, "dPde_rho");
	dTdrho_e = Interp2D_lin(x, y, "dTdrho_e");
	dTde_rho = Interp2D_lin(x, y, "dTde_rho");
	Cp = Interp2D_lin(x, y, "Cp");
	Mu = Interp2D_lin(x, y, "Mu");
	dmudrho_T = Interp2D_lin(x, y, "dmudrho_T");
	dmudT_rho = Interp2D_lin(x, y, "dmudT_rho");
	Kt = Interp2D_lin(x, y, "Kt");
	dktdrho_T = Interp2D_lin(x, y, "dktdrho_T");
	dktdT_rho = Interp2D_lin(x, y, "dktdT_rho");

	//cout<<"Interpolated fit:"<<endl;
	//	CTLprint ();
	if ((Density > Density_Table_Limits[1])
			or (Density < Density_Table_Limits[0])) {
		cerr << "PS Interpolated Density out of bounds\n";
	}
	if ((Pressure > Pressure_Table_Limits[1])
			or (Pressure < Pressure_Table_Limits[0])) {
		cerr << "PS Interpolated Pressure out of bounds\n";
	}
	delete [] Nearest_Neighbour_iIndex;
	delete [] Nearest_Neighbour_jIndex;
}

void CLookUpTable::SetTDState_rhoT(su2double rho, su2double T) {
	//Check if inputs are in total range (necessary but not sufficient condition)
	if ((rho > Density_Table_Limits[1]) or (rho < Density_Table_Limits[0])) {
		cerr << "RHOT Input Density out of bounds\n";
	}
	if ((T > Temperature_Table_Limits[1]) or (T < Temperature_Table_Limits[0])) {
		cerr << "RHOT Input Temperature out of bounds\n";
	}

	unsigned int LowerI, UpperI, LowerJ, UpperJ, middleI, middleJ;
	Nearest_Neighbour_iIndex = new int[4];
	Nearest_Neighbour_jIndex = new int[4];
	UpperJ = Table_Pressure_Stations - 1;
	LowerJ = 0;
	su2double grad, x00, y00, y10, x10, RunVal;

	//Determine the I index: rho is not equispaced
	UpperI = Table_Density_Stations - 1;
	LowerI = 0;

	while (UpperI - LowerI > 1) {
		middleI = (UpperI + LowerI) / 2;
		//Check current value
		x00 = ThermoTables[middleI][LowerJ].Density;
		grad = ThermoTables[middleI + 1][LowerJ].Density - x00;
		if (x00 * grad > rho * grad) {
			UpperI = middleI;
		} else if (x00 < rho) {
			LowerI = middleI;
		}
		if (x00 == rho) {
			UpperI = LowerI + 1;
			break;
		}
	}

	//Determine the J index (for T)
	while (UpperJ - LowerJ > 1) {
		middleJ = (UpperJ + LowerJ) / 2;
		//Check current value
		y00 = ThermoTables[LowerI][middleJ].Temperature;
		y10 = ThermoTables[UpperI][middleJ].Temperature;
		x00 = ThermoTables[LowerI][middleJ].Density;
		x10 = ThermoTables[UpperI][middleJ].Density;
		RunVal = y00 + (y10 - y00) / (x10 - x00) * (rho - x00);
		grad = ThermoTables[LowerI][UpperJ].Temperature - y00;
		if (RunVal * grad > T * grad) {
			UpperJ = middleJ;
		} else if (RunVal * grad < T * grad) {
			LowerJ = middleJ;
		}
		if (RunVal == T) {
			UpperJ = LowerJ + 1;
			break;
		}
	}

	iIndex = LowerI;
	jIndex = LowerJ;

	su2double x, y;
	x = (rho - ThermoTables[iIndex][jIndex].Density) / Density_Table_Limits[1];
	y = (T - ThermoTables[iIndex][jIndex].Temperature)
							/ Temperature_Table_Limits[1];
	//Set the nearest neigbours to the adjacent i and j vertexes
	Nearest_Neighbour_iIndex[0] = iIndex;
	Nearest_Neighbour_iIndex[1] = iIndex + 1;
	Nearest_Neighbour_iIndex[2] = iIndex;
	Nearest_Neighbour_iIndex[3] = iIndex + 1;
	Nearest_Neighbour_jIndex[0] = jIndex;
	Nearest_Neighbour_jIndex[1] = jIndex;
	Nearest_Neighbour_jIndex[2] = jIndex + 1;
	Nearest_Neighbour_jIndex[3] = jIndex + 1;
	//Determine interpolation coefficients
	Interp2D_ArbitrarySkewCoeff(x, y, "RHOT");

	Temperature = T;
	Density = rho;
	StaticEnergy = Interp2D_lin(x, y, "StaticEnergy");
	Enthalpy = Interp2D_lin(x, y, "Enthalpy");
	Entropy = Interp2D_lin(x, y, "Entropy");
	Pressure = Interp2D_lin(x, y, "Pressure");
	SoundSpeed2 = Interp2D_lin(x, y, "SoundSpeed2");
	dPdrho_e = Interp2D_lin(x, y, "dPdrho_e");
	dPde_rho = Interp2D_lin(x, y, "dPde_rho");
	dTdrho_e = Interp2D_lin(x, y, "dTdrho_e");
	dTde_rho = Interp2D_lin(x, y, "dTde_rho");
	Cp = Interp2D_lin(x, y, "Cp");
	Mu = Interp2D_lin(x, y, "Mu");
	dmudrho_T = Interp2D_lin(x, y, "dmudrho_T");
	dmudT_rho = Interp2D_lin(x, y, "dmudT_rho");
	Kt = Interp2D_lin(x, y, "Kt");
	dktdrho_T = Interp2D_lin(x, y, "dktdrho_T");
	dktdT_rho = Interp2D_lin(x, y, "dktdT_rho");

	//cout<<"Interpolated fit:"<<endl;
	if ((Density > Density_Table_Limits[1])
			or (Density < Density_Table_Limits[0])) {
		cerr << "RHOT Interpolated Density out of bounds\n";
	}
	if ((Pressure > Pressure_Table_Limits[1])
			or (Pressure < Pressure_Table_Limits[0])) {
		cerr << "RHOT Interpolated Pressure out of bounds\n";
	}
	delete [] Nearest_Neighbour_iIndex;
	delete [] Nearest_Neighbour_jIndex;
}

void CLookUpTable::Interp2D_ArbitrarySkewCoeff(su2double x, su2double y,
		std::string grid_var) {
	//Distances in along x and y axis are taken relative to i,j point (x00, y00).
	//This reduces the interpolation to a 3by3 system rather than 4by4
	//x and y are not strictrly necessary for the calculation of the coefficients.
	//However, they do allow for checking whether the point of interest is contained in
	//the quad under consideration.
	su2double x_in, y_in, x00, y00, dx10, dx01, dx11, dy10, dy01, dy11;
	//Interpolation LHM
	su2double Vandermonde[3][3];
	//Helper variable for Gaussian elimination
	su2double c;
	//Load in the coordinates of the qudrilateral (values relative to i,j)
	if (grid_var == "RHOE") {
		x00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Density
				/ Density_Table_Limits[1];
		y00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].StaticEnergy
				/ StaticEnergy_Table_Limits[1];
		dx01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Density
				/ Density_Table_Limits[1] - x00;
		dy01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].StaticEnergy
				/ StaticEnergy_Table_Limits[1] - y00;
		dx10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Density
				/ Density_Table_Limits[1] - x00;
		dy10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].StaticEnergy
				/ StaticEnergy_Table_Limits[1] - y00;
		dx11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Density
				/ Density_Table_Limits[1] - x00;
		dy11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].StaticEnergy
				/ StaticEnergy_Table_Limits[1] - y00;
	} else if (grid_var == "PT") {
		y00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Pressure
				/ Pressure_Table_Limits[1];
		x00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Temperature
				/ Temperature_Table_Limits[1];
		dy01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Pressure
				/ Pressure_Table_Limits[1] - y00;
		dx01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Temperature
				/ Temperature_Table_Limits[1] - x00;
		dy10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Pressure
				/ Pressure_Table_Limits[1] - y00;
		dx10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Temperature
				/ Temperature_Table_Limits[1] - x00;
		dy11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Pressure
				/ Pressure_Table_Limits[1] - y00;
		dx11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Temperature
				/ Temperature_Table_Limits[1] - x00;
	} else if (grid_var == "PRHO") {
		y00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Pressure
				/ Pressure_Table_Limits[1];
		x00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Density
				/ Density_Table_Limits[1];
		dy01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Pressure
				/ Pressure_Table_Limits[1] - y00;
		dx01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Density
				/ Density_Table_Limits[1] - x00;
		dy10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Pressure
				/ Pressure_Table_Limits[1] - y00;
		dx10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Density
				/ Density_Table_Limits[1] - x00;
		dy11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Pressure
				/ Pressure_Table_Limits[1] - y00;
		dx11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Density
				/ Density_Table_Limits[1] - x00;
	} else if (grid_var == "RHOT") {
		x00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Density
				/ Density_Table_Limits[1];
		y00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Temperature
				/ Temperature_Table_Limits[1];
		dx01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Density
				/ Density_Table_Limits[1] - x00;
		dy01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Temperature
				/ Temperature_Table_Limits[1] - y00;
		dx10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Density
				/ Density_Table_Limits[1] - x00;
		dy10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Temperature
				/ Temperature_Table_Limits[1] - y00;
		dx11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Density
				/ Density_Table_Limits[1] - x00;
		dy11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Temperature
				/ Temperature_Table_Limits[1] - y00;
	} else if (grid_var == "PS") {
		y00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Pressure
				/ Pressure_Table_Limits[1];
		x00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Entropy
				/ Entropy_Table_Limits[1];
		dy01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Pressure
				/ Pressure_Table_Limits[1] - y00;
		dx01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Entropy
				/ Entropy_Table_Limits[1] - x00;
		dy10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Pressure
				/ Pressure_Table_Limits[1] - y00;
		dx10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Entropy
				/ Entropy_Table_Limits[1] - x00;
		dy11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Pressure
				/ Pressure_Table_Limits[1] - y00;
		dx11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Entropy
				/ Entropy_Table_Limits[1] - x00;
	} else if (grid_var == "HS") {
		x00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Enthalpy
				/ Enthalpy_Table_Limits[1];
		y00 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Entropy
				/ Entropy_Table_Limits[1];
		dx01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Enthalpy
				/ Enthalpy_Table_Limits[1] - x00;
		dy01 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Entropy
				/ Entropy_Table_Limits[1] - y00;
		dx10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Enthalpy
				/ Enthalpy_Table_Limits[1] - x00;
		dy10 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Entropy
				/ Entropy_Table_Limits[1] - y00;
		dx11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Enthalpy
				/ Enthalpy_Table_Limits[1] - x00;
		dy11 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Entropy
				/ Entropy_Table_Limits[1] - y00;
	}
	//Check if x, y is indeed in the quad
	//Some extra logic is needed as the both monotonically increasing and monotonically decreasing functions
	//have to be anticipated
	bool BOTTOM, TOP, LEFT, RIGHT;
	BOTTOM = (y * dx10) < (x * dy10);
	TOP = ((y - dy01) * (dx11 - dx01)) > ((dy11 - dy01) * (x - dx01));
	RIGHT = ((x - dx10) * (dy11 - dy10)) > ((dx11 - dx10) * (y - dy10));
	LEFT = (x * dy01) < (dx01 * y);
	//Check BOTTOM quad boundary
	if (BOTTOM and !TOP) {
		//added table limit detection
		if (jIndex == 0) {
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< " interpolation point lies below the LUT\n";
		} else {
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< " interpolation point lies below bottom boundary of selected quad\n";
		}
	}
	//Check RIGHT quad boundary
	if (RIGHT and !LEFT) {
		//added table limit detection
		if (iIndex == (Table_Density_Stations - 1)) {
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< " interpolation point lies right of the LUT\n";
		} else {
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< " interpolation point lies to the right of the boundary of selected quad\n";
		}
	}
	//Check TOP quad boundary
	if (TOP and !BOTTOM) {
		//added table limit detection
		if (jIndex == Table_Pressure_Stations - 1) {
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< " interpolation point lies above the LUT\n";
		} else {
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< +" interpolation point lies above the boundary of selected quad\n";
		}
	}
	//Check LEFT quad boundary

	if (LEFT and !RIGHT) {
		//added table limit detection
		if (iIndex == 0) {
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< " interpolation point lies left of the LUT\n";
		} else {
			cerr << grid_var << ' ' << Nearest_Neighbour_iIndex[0] << ", "
					<< Nearest_Neighbour_jIndex[0]
																			<< +" interpolation point lies to the left of the boundary of selected quad\n";
		}
	}

	//Setup the LHM matrix for the interpolation
	Vandermonde[0][0] = dx10;
	Vandermonde[0][1] = dy10;
	Vandermonde[0][2] = dx10 * dy10;
	Vandermonde[1][0] = dx01;
	Vandermonde[1][1] = dy01;
	Vandermonde[1][2] = dx01 * dy01;
	Vandermonde[2][0] = dx11;
	Vandermonde[2][1] = dy11;
	Vandermonde[2][2] = dx11 * dy11;

	if ((dx01 == 0) and (dy10 == 0)) {
		Interpolation_Coeff[0][0] = 1.0 / dx10;
		Interpolation_Coeff[0][1] = 0;
		Interpolation_Coeff[0][2] = 0;
		Interpolation_Coeff[1][0] = 0;
		Interpolation_Coeff[1][1] = 1.0 / dy01;
		Interpolation_Coeff[1][2] = 0;
		Interpolation_Coeff[2][0] = -1 / (dx10 * dy11);
		Interpolation_Coeff[2][1] = -1 / (dx11 * dy01);
		Interpolation_Coeff[2][2] = 1 / (dx11 * dy11);
	}

	else if (dx01 == 0) {
		Interpolation_Coeff[0][0] = dx11 * dy01 * dy11
				/ (pow(dx10, 2) * dx11 * dy10 + dx10 * dx11 * dy01 * dy11
						+ dx10 * dx11 * dy10 * (-dx10 - dy01));
		Interpolation_Coeff[0][1] = (dx10 * dy10 * dy11 - dx11 * dy10 * dy11)
								/ (pow(dx10, 2) * dx11 * dy10 + dx10 * dx11 * dy01 * dy11
										+ dx10 * dx11 * dy10 * (-dx10 - dy01));
		Interpolation_Coeff[0][2] = -dx10 * dy01 * dy10
				/ (pow(dx10, 2) * dx11 * dy10 + dx10 * dx11 * dy01 * dy11
						+ dx10 * dx11 * dy10 * (-dx10 - dy01));
		Interpolation_Coeff[1][0] = 0;
		Interpolation_Coeff[1][1] = (-dx10 * dx11 * dy10 + dx10 * dx11 * dy11)
								/ (pow(dx10, 2) * dx11 * dy10 + dx10 * dx11 * dy01 * dy11
										+ dx10 * dx11 * dy10 * (-dx10 - dy01));
		Interpolation_Coeff[1][2] = 0;
		Interpolation_Coeff[2][0] = -dx11 * dy01
				/ (pow(dx10, 2) * dx11 * dy10 + dx10 * dx11 * dy01 * dy11
						+ dx10 * dx11 * dy10 * (-dx10 - dy01));
		Interpolation_Coeff[2][1] = (-dx10 * dy11 + dx11 * dy10)
								/ (pow(dx10, 2) * dx11 * dy10 + dx10 * dx11 * dy01 * dy11
										+ dx10 * dx11 * dy10 * (-dx10 - dy01));
		Interpolation_Coeff[2][2] = dx10 * dy01
				/ (pow(dx10, 2) * dx11 * dy10 + dx10 * dx11 * dy01 * dy11
						+ dx10 * dx11 * dy10 * (-dx10 - dy01));
	} else if (dy10 == 0) {
		Interpolation_Coeff[0][0] = (-dx01 * dy01 * dy11 + dx11 * dy01 * dy11)
								/ (dx01 * pow(dy01, 2) * dy11 + dx01 * dy01 * dy11 * (-dx10 - dy01)
										+ dx10 * dx11 * dy01 * dy11);
		Interpolation_Coeff[0][1] = 0;
		Interpolation_Coeff[0][2] = 0;
		Interpolation_Coeff[1][0] = (dx01 * dx11 * dy01 - dx01 * dx11 * dy11)
								/ (dx01 * pow(dy01, 2) * dy11 + dx01 * dy01 * dy11 * (-dx10 - dy01)
										+ dx10 * dx11 * dy01 * dy11);
		Interpolation_Coeff[1][1] = dx10 * dx11 * dy11
				/ (dx01 * pow(dy01, 2) * dy11 + dx01 * dy01 * dy11 * (-dx10 - dy01)
						+ dx10 * dx11 * dy01 * dy11);
		Interpolation_Coeff[1][2] = -dx01 * dx10 * dy01
				/ (dx01 * pow(dy01, 2) * dy11 + dx01 * dy01 * dy11 * (-dx10 - dy01)
						+ dx10 * dx11 * dy01 * dy11);
		Interpolation_Coeff[2][0] = (dx01 * dy11 - dx11 * dy01)
								/ (dx01 * pow(dy01, 2) * dy11 + dx01 * dy01 * dy11 * (-dx10 - dy01)
										+ dx10 * dx11 * dy01 * dy11);
		Interpolation_Coeff[2][1] = -dx10 * dy11
				/ (dx01 * pow(dy01, 2) * dy11 + dx01 * dy01 * dy11 * (-dx10 - dy01)
						+ dx10 * dx11 * dy01 * dy11);
		Interpolation_Coeff[2][2] = dx10 * dy01
				/ (dx01 * pow(dy01, 2) * dy11 + dx01 * dy01 * dy11 * (-dx10 - dy01)
						+ dx10 * dx11 * dy01 * dy11);

	} else {
		Interpolation_Coeff[0][0] = (-dx01 * dy01 * dy11 + dx11 * dy01 * dy11)
								/ (dx11 * dy11 * (-dx01 * dy10 + dx10 * dy01)
										+ dx11 * (dx01 * dy01 * dy10 + pow(dx10, 2) * dy10)
										+ dy11 * (dx01 * dx10 * dy10 + dx01 * pow(dy01, 2))
										- (-dx10 - dy01) * (-dx01 * dy01 * dy11 - dx10 * dx11 * dy10));
		Interpolation_Coeff[0][1] = (dx10 * dy10 * dy11 - dx11 * dy10 * dy11)
								/ (dx11 * dy11 * (-dx01 * dy10 + dx10 * dy01)
										+ dx11 * (dx01 * dy01 * dy10 + pow(dx10, 2) * dy10)
										+ dy11 * (dx01 * dx10 * dy10 + dx01 * pow(dy01, 2))
										- (-dx10 - dy01) * (-dx01 * dy01 * dy11 - dx10 * dx11 * dy10));
		Interpolation_Coeff[0][2] = (dx01 * dy01 * dy10 - dx10 * dy01 * dy10)
								/ (dx11 * dy11 * (-dx01 * dy10 + dx10 * dy01)
										+ dx11 * (dx01 * dy01 * dy10 + pow(dx10, 2) * dy10)
										+ dy11 * (dx01 * dx10 * dy10 + dx01 * pow(dy01, 2))
										- (-dx10 - dy01) * (-dx01 * dy01 * dy11 - dx10 * dx11 * dy10));
		Interpolation_Coeff[1][0] = (dx01 * dx11 * dy01 - dx01 * dx11 * dy11)
								/ (dx11 * dy11 * (-dx01 * dy10 + dx10 * dy01)
										+ dx11 * (dx01 * dy01 * dy10 + pow(dx10, 2) * dy10)
										+ dy11 * (dx01 * dx10 * dy10 + dx01 * pow(dy01, 2))
										- (-dx10 - dy01) * (-dx01 * dy01 * dy11 - dx10 * dx11 * dy10));
		Interpolation_Coeff[1][1] = (-dx10 * dx11 * dy10 + dx10 * dx11 * dy11)
								/ (dx11 * dy11 * (-dx01 * dy10 + dx10 * dy01)
										+ dx11 * (dx01 * dy01 * dy10 + pow(dx10, 2) * dy10)
										+ dy11 * (dx01 * dx10 * dy10 + dx01 * pow(dy01, 2))
										- (-dx10 - dy01) * (-dx01 * dy01 * dy11 - dx10 * dx11 * dy10));
		Interpolation_Coeff[1][2] = (-dx01 * dx10 * dy01 + dx01 * dx10 * dy10)
								/ (dx11 * dy11 * (-dx01 * dy10 + dx10 * dy01)
										+ dx11 * (dx01 * dy01 * dy10 + pow(dx10, 2) * dy10)
										+ dy11 * (dx01 * dx10 * dy10 + dx01 * pow(dy01, 2))
										- (-dx10 - dy01) * (-dx01 * dy01 * dy11 - dx10 * dx11 * dy10));
		Interpolation_Coeff[2][0] = (dx01 * dy11 - dx11 * dy01)
								/ (dx11 * dy11 * (-dx01 * dy10 + dx10 * dy01)
										+ dx11 * (dx01 * dy01 * dy10 + pow(dx10, 2) * dy10)
										+ dy11 * (dx01 * dx10 * dy10 + dx01 * pow(dy01, 2))
										- (-dx10 - dy01) * (-dx01 * dy01 * dy11 - dx10 * dx11 * dy10));
		Interpolation_Coeff[2][1] = (-dx10 * dy11 + dx11 * dy10)
								/ (dx11 * dy11 * (-dx01 * dy10 + dx10 * dy01)
										+ dx11 * (dx01 * dy01 * dy10 + pow(dx10, 2) * dy10)
										+ dy11 * (dx01 * dx10 * dy10 + dx01 * pow(dy01, 2))
										- (-dx10 - dy01) * (-dx01 * dy01 * dy11 - dx10 * dx11 * dy10));
		Interpolation_Coeff[2][2] = (-dx01 * dy10 + dx10 * dy01)
								/ (dx11 * dy11 * (-dx01 * dy10 + dx10 * dy01)
										+ dx11 * (dx01 * dy01 * dy10 + pow(dx10, 2) * dy10)
										+ dy11 * (dx01 * dx10 * dy10 + dx01 * pow(dy01, 2))
										- (-dx10 - dy01) * (-dx01 * dy01 * dy11 - dx10 * dx11 * dy10));
	}

	//cout<<"Interpolation LHM matrix \n"<<"[";

	//	for (int j=0; j<3; j++)
	//	{
	//		cout<<setw(15)<<"["<<Vandermonde[j][0]<<" ,  "<<Vandermonde[j][1]<<"  , "<<Vandermonde[j][2]<<"]"<<endl;
	//	}
	//	cout<<"]\n";

	/*//Store the inverse of the LHM matrix as coeff
	 Interpolation_Coeff[0][0] = 1;
	 Interpolation_Coeff[0][1] = 0;
	 Interpolation_Coeff[0][2] = 0;
	 Interpolation_Coeff[1][0] = 0;
	 Interpolation_Coeff[1][1] = 1;
	 Interpolation_Coeff[1][2] = 0;//solved interpolation bug
	 Interpolation_Coeff[2][0] = 0;
	 Interpolation_Coeff[2][1] = 0;
	 Interpolation_Coeff[2][2] = 1;

	 //Compute inverse of LHM using Gaussian elimination
	 //Reduced Echelon form of the LHM
	 if (Vandermonde[0][0] != 0)
	 {
	 c = Vandermonde[1][0]/Vandermonde[0][0];
	 Interpolation_Coeff[1][0] = Interpolation_Coeff[1][0] -Interpolation_Coeff[0][0]*c;
	 for (int i=0; i<3; i++)
	 {
	 Vandermonde[1][i] = Vandermonde[1][i] -Vandermonde[0][i]*c;
	 }
	 c =Vandermonde[2][0]/Vandermonde[0][0];
	 Interpolation_Coeff[2][0] = Interpolation_Coeff[2][0] -Interpolation_Coeff[0][0]*c;

	 for (int i=0; i<3; i++)
	 {
	 Vandermonde[2][i] = Vandermonde[2][i] -Vandermonde[0][i]*c;
	 }
	 }

	 if (Vandermonde[1][1] != 0)
	 {
	 c = Vandermonde[2][1]/Vandermonde[1][1];
	 for (int i=0; i<2; i++)
	 {
	 Interpolation_Coeff[2][i] = Interpolation_Coeff[2][i] -Interpolation_Coeff[1][i]*c;
	 }
	 for (int i=0; i<3; i++)
	 {
	 Vandermonde[2][i] = Vandermonde[2][i] -Vandermonde[1][i]*c;
	 }
	 }
	 //Reduced reduced Echelon form of LHM
	 if (Vandermonde[2][2] != 0)
	 {

	 for (int i=0; i<3; i++)
	 {
	 Interpolation_Coeff[1][i] = Interpolation_Coeff[1][i] -Interpolation_Coeff[2][i]*Vandermonde[1][2]/Vandermonde[2][2];
	 }
	 Vandermonde[1][2] = 0;
	 for (int i=0; i<3; i++)
	 {
	 Interpolation_Coeff[0][i] = Interpolation_Coeff[0][i] -Interpolation_Coeff[2][i]*Vandermonde[0][2]/Vandermonde[2][2];
	 }
	 Vandermonde[0][2] = 0;
	 }
	 if (Vandermonde[1][1] != 0)
	 {
	 for (int i=0; i<3; i++)
	 {
	 Interpolation_Coeff[0][i] = Interpolation_Coeff[0][i] -Interpolation_Coeff[1][i]*Vandermonde[0][1]/Vandermonde[1][1];
	 }
	 }
	 Vandermonde[0][1] = 0;

	 //Normalize the RR Echelon form
	 if (Vandermonde[0][0] != 0)
	 {
	 for (int i=0; i<3; i++)
	 {
	 Interpolation_Coeff[0][i] = Interpolation_Coeff[0][i]/Vandermonde[0][0];
	 }
	 if (Vandermonde[1][1] != 0)
	 {
	 for (int i=0; i<3; i++)
	 {
	 Interpolation_Coeff[1][i] = Interpolation_Coeff[1][i]/Vandermonde[1][1];
	 }
	 }
	 if (Vandermonde[2][2] != 0)
	 {
	 for (int i=0; i<3; i++)
	 {
	 Interpolation_Coeff[2][i] = Interpolation_Coeff[2][i]/Vandermonde[2][2];
	 }

	 }
	 }*/
	return;
}
su2double CLookUpTable::Interp2D_Inv_Dist(int N, std::string interpolant_var,
		su2double* dist) {
	su2double interp_result = 0;
	//The function values to interpolate from
	su2double *Interpolation_RHS= new	su2double[N];
	//For each case the values are filled differently
	if (interpolant_var == "StaticEnergy") {
		for (int i=0;i<N;i++)
		{
			Interpolation_RHS[i] = ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].StaticEnergy;
		}

	} else if (interpolant_var == "Entropy") {
		for (int i=0;i<N;i++)
		{
			Interpolation_RHS[i] = ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].Entropy;
		}
	} else if (interpolant_var == "Density") {
		for (int i=0;i<N;i++)
		{
			Interpolation_RHS[i] = ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].Density;
		}
	} else if (interpolant_var == "Pressure") {
		for (int i=0;i<N;i++)
		{
			Interpolation_RHS[i] = ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].Pressure;
		}
	} else if (interpolant_var == "SoundSpeed2") {
		for (int i=0;i<N;i++)
		{
			Interpolation_RHS[i] = ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].SoundSpeed2;
		}
	} else if (interpolant_var == "Temperature") {
		for (int i=0;i<N;i++)
		{
			Interpolation_RHS[i] = ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].Temperature;
		}
	} else if (interpolant_var == "dPdrho_e") {
		for (int i=0;i<N;i++)
		{
			Interpolation_RHS[i] = ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].dPdrho_e;
		}


	} else if (interpolant_var == "dPde_rho") {
		for (int i=0;i<N;i++)
		{
			Interpolation_RHS[i] = ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].dPde_rho;
		}
	} else if (interpolant_var == "dTdrho_e") {
		for (int i=0;i<N;i++)
		{
			Interpolation_RHS[i] = ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].dTdrho_e;
		}
	} else if (interpolant_var == "dTde_rho") {
		for (int i=0;i<N;i++)
		{
			Interpolation_RHS[i] = ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].dTde_rho;
		}
	} else if (interpolant_var == "Cp") {
		for (int i=0;i<N;i++)
		{
			Interpolation_RHS[i] = ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].Cp;
		}
	} else if (interpolant_var == "Mu") {
		for (int i=0;i<N;i++)
		{
			Interpolation_RHS[i] = ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].Mu;
		}

	} else if (interpolant_var == "dmudrho_T") {
		for (int i=0;i<N;i++)
		{
			Interpolation_RHS[i] = ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].dmudrho_T;
		}

	} else if (interpolant_var == "dmudT_rho") {
		for (int i=0;i<N;i++)
		{
			Interpolation_RHS[i] = ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].dmudT_rho;
		}

	} else if (interpolant_var == "Kt") {
		for (int i=0;i<N;i++)
		{
			Interpolation_RHS[i] = ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].Kt;
		}


	} else if (interpolant_var == "dktdrho_T") {
		for (int i=0;i<N;i++)
		{
			Interpolation_RHS[i] = ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].dktdrho_T;
		}

	} else if (interpolant_var == "dktdT_rho") {
		for (int i=0;i<N;i++)
		{
			Interpolation_RHS[i] = ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].dktdT_rho;
		}
	} else if (interpolant_var == "Enthalpy") {
		for (int i=0;i<N;i++)
		{
			Interpolation_RHS[i] = ThermoTables[Nearest_Neighbour_iIndex[i]][Nearest_Neighbour_jIndex[i]].Enthalpy;
		}

	}

	su2double dist_sum = 0;
	for (int i = 0; i < N; i++) {
		interp_result += (1 / dist[i]) * Interpolation_RHS[i];
		dist_sum += 1 / dist[i];
	}

	interp_result = interp_result / dist_sum;
	delete [] Interpolation_RHS;
	return interp_result;
}

su2double CLookUpTable::Interp2D_lin(su2double x, su2double y,
		string interpolant_var) {
	//F is the RHS part of the interpolation equation
	su2double Interpolation_RHS[3];
	//The solution vector for the interpolation equation
	su2double Interpolation_Sol_Vector[3];
	//The function values
	su2double func_value_at_i0j0, func_value_at_i1j0, func_value_at_i0j1,
	func_value_at_i1j1;
	//For each case the values are filled differently
	if (interpolant_var == "StaticEnergy") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].StaticEnergy;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].StaticEnergy;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].StaticEnergy;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].StaticEnergy;
	} else if (interpolant_var == "Entropy") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Entropy;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Entropy;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Entropy;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Entropy;
	} else if (interpolant_var == "Density") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Density;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Density;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Density;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Density;
	} else if (interpolant_var == "Pressure") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Pressure;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Pressure;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Pressure;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Pressure;
	} else if (interpolant_var == "SoundSpeed2") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].SoundSpeed2;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].SoundSpeed2;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].SoundSpeed2;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].SoundSpeed2;
	} else if (interpolant_var == "Temperature") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Temperature;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Temperature;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Temperature;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Temperature;
	} else if (interpolant_var == "dPdrho_e") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dPdrho_e;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dPdrho_e;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dPdrho_e;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dPdrho_e;
	} else if (interpolant_var == "dPde_rho") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dPde_rho;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dPde_rho;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dPde_rho;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dPde_rho;
	} else if (interpolant_var == "dTdrho_e") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dTdrho_e;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dTdrho_e;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dTdrho_e;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dTdrho_e;
	} else if (interpolant_var == "dTde_rho") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dTde_rho;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dTde_rho;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dTde_rho;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dTde_rho;
	} else if (interpolant_var == "Cp") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Cp;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Cp;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Cp;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Cp;
	} else if (interpolant_var == "Mu") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Mu;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Mu;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Mu;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Mu;
	} else if (interpolant_var == "dmudrho_T") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dmudrho_T;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dmudrho_T;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dmudrho_T;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dmudrho_T;
	} else if (interpolant_var == "dmudT_rho") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dmudT_rho;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dmudT_rho;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dmudT_rho;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dmudT_rho;
	} else if (interpolant_var == "Kt") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Kt;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Kt;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Kt;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Kt;
	} else if (interpolant_var == "dktdrho_T") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dktdrho_T;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dktdrho_T;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dktdrho_T;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dktdrho_T;
	} else if (interpolant_var == "dktdT_rho") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].dktdT_rho;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].dktdT_rho;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].dktdT_rho;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].dktdT_rho;
	} else if (interpolant_var == "Enthalpy") {
		func_value_at_i0j0 =
				ThermoTables[Nearest_Neighbour_iIndex[0]][Nearest_Neighbour_jIndex[0]].Enthalpy;
		func_value_at_i1j0 =
				ThermoTables[Nearest_Neighbour_iIndex[1]][Nearest_Neighbour_jIndex[1]].Enthalpy;
		func_value_at_i0j1 =
				ThermoTables[Nearest_Neighbour_iIndex[2]][Nearest_Neighbour_jIndex[2]].Enthalpy;
		func_value_at_i1j1 =
				ThermoTables[Nearest_Neighbour_iIndex[3]][Nearest_Neighbour_jIndex[3]].Enthalpy;
	}

	//Using offset relative to i,j point yields a 3by3 system rather than 4by4
	Interpolation_RHS[0] = func_value_at_i1j0 - func_value_at_i0j0;
	Interpolation_Sol_Vector[0] = 0;
	Interpolation_RHS[1] = func_value_at_i0j1 - func_value_at_i0j0;
	Interpolation_Sol_Vector[1] = 0;
	Interpolation_RHS[2] = func_value_at_i1j1 - func_value_at_i0j0;
	Interpolation_Sol_Vector[2] = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Interpolation_Sol_Vector[i] = Interpolation_Sol_Vector[i]
																														 + Interpolation_RHS[i] * Interpolation_Coeff[i][j];
		}
	}

	return func_value_at_i0j0 + Interpolation_Sol_Vector[0] * x
			+ Interpolation_Sol_Vector[1] * y + Interpolation_Sol_Vector[2] * x * y;
}

void CLookUpTable::LUTprint(void) {
	//	for (int i=0; i<Table_Density_Stations; i++)
	//	{
	//		for (int j=0; j<Table_Pressure_Stations; j++)
	//		{
	////			ThermoTables[i][j].CTLprint();
	//		}
	//	}
}

void CLookUpTable::RecordState(char* file) {

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
	fs << dmudrho_T << ", ";
	fs << dmudT_rho << ", ";
	fs << Kt << ", ";
	fs << dktdrho_T << ", ";
	fs << dktdT_rho << ", ";
	fs << "\n";
	fs.close();
}

void CLookUpTable::TableDump(char* filename) {
	for (int i = 0; i < Table_Density_Stations; i++) {
		for (int j = 0; j < Table_Pressure_Stations; j++) {
			iIndex = i;
			jIndex = j;
			Temperature = ThermoTables[iIndex][jIndex].Temperature;
			Density = ThermoTables[iIndex][jIndex].Density;
			Enthalpy = ThermoTables[iIndex][jIndex].Enthalpy;
			StaticEnergy = ThermoTables[iIndex][jIndex].StaticEnergy;
			Entropy = ThermoTables[iIndex][jIndex].Entropy;
			Pressure = ThermoTables[iIndex][jIndex].Pressure;
			SoundSpeed2 = ThermoTables[iIndex][jIndex].SoundSpeed2;
			dPdrho_e = ThermoTables[iIndex][jIndex].dPdrho_e;
			dPde_rho = ThermoTables[iIndex][jIndex].dPde_rho;
			dTdrho_e = ThermoTables[iIndex][jIndex].dTdrho_e;
			dTde_rho = ThermoTables[iIndex][jIndex].dTde_rho;
			Cp = ThermoTables[iIndex][jIndex].Cp;
			Mu = ThermoTables[iIndex][jIndex].Mu;
			dmudrho_T = ThermoTables[iIndex][jIndex].dmudrho_T;
			dmudT_rho = ThermoTables[iIndex][jIndex].dmudT_rho;
			Kt = ThermoTables[iIndex][jIndex].Kt;
			dktdrho_T = ThermoTables[iIndex][jIndex].dktdrho_T;
			dktdT_rho = ThermoTables[iIndex][jIndex].dktdT_rho;
			RecordState(filename);
		}
	}

}

void CLookUpTable::TableLoadCFX(string filename) {
	int N_PARAM = 0;
	int set_x = 0;
	int set_y = 0;
	int var_steps = 0;
	int var_scanned = 0;

	string line;
	string value;

	ifstream table(filename.c_str());
	assert(table.is_open());
	//cout<<"Looking for number of parameters"<<endl;
	while (getline(table, line)) {
		unsigned int found;
		found = line.find("$$PARAM");
		if (found < 10) {
			getline(table, line);
			istringstream in(line);
			in >> N_PARAM;
			N_PARAM++;
			//	cout<<"Number of parameters "<<N_PARAM<<endl;
		}
		for (int var = var_scanned; var < N_PARAM + 1; var++) {
			string svar =
					static_cast<ostringstream*>(&(ostringstream() << var))->str();
			found = line.find("$TABLE_" + svar);
			if (found < 10) {
				var_scanned = var;
				//cout<<found<<' '<<line<<endl;
				getline(table, line);
				istringstream in(line);
				int x, y;
				in >> x >> y;
				if (var == 1) {
					ThermoTables = new CThermoList*[x];
					for (int i = 0; i < x; i++) {
						ThermoTables[i] = new CThermoList[y];
					}
					set_x = x;
					set_y = y;
					//		cout<<"Tables have been allocated"<<var<<endl;

					//Fill in the densities
					Density_Table_Limits[0] = 10E15;
					Density_Table_Limits[1] = 0;
					var_steps = 10;
					su2double* vD = new su2double[set_x];

					for (int k = 0; k < ceil(float(set_x) / 10.0); k++) {
						getline(table, line);
						//		//cout<<line<<endl;
						istringstream inD(line);
						if ((set_x - k * 10) < 10)
							var_steps = (set_x - k * 10);
						for (int i = 0; i < var_steps; i++) {
							inD >> vD[10 * k + i];
							if (vD[10 * k + i] > Density_Table_Limits[1]) {
								Density_Table_Limits[1] = vD[10 * k + i];
							}
							if (vD[10 * k + i] < Density_Table_Limits[0]) {
								Density_Table_Limits[0] = vD[10 * k + i];
							}
						}
					}
					for (int i = 0; i < set_x; i++) {
						for (int j = 0; j < set_y; j++) {
							ThermoTables[i][j].Density = vD[i];
						}
					}
					delete vD;

					//Fill in the pressures
					su2double* vP = new su2double[set_y];
					var_steps = 10; //solved pressure reading bug
					Pressure_Table_Limits[0] = 10E15; //lower limit
					Pressure_Table_Limits[1] = 0; //upper limit
					//Each line contains at most 10 pressure values
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++) {

						getline(table, line);
						//cout<<line<<endl;
						istringstream inP(line);
						//Check if line contains less than 10 values
						if ((set_y - k * 10) < 10)
							var_steps = (set_y - k * 10);
						for (int j = 0; j < var_steps; j++) {
							inP >> vP[10 * k + j];
							if (vP[10 * k + j] > Pressure_Table_Limits[1]) {
								Pressure_Table_Limits[1] = vP[10 * k + j];
							}
							if (vP[10 * k + j] < Pressure_Table_Limits[0]) {
								Pressure_Table_Limits[0] = vP[10 * k + j];
							}
						}
					}
					for (int i = 0; i < set_x; i++) {
						for (int j = 0; j < set_y; j++) {
							ThermoTables[i][j].Pressure = vP[j];
						}
					}
					delete vP;
					//	//cout<<"Tables have been filled with D and P values "<<var<<endl;

				}
				// Check that additional tables all adhere to the same x,y dimensions, otherwise throw an error
				else if (x != set_x && y != set_y) {
					cerr
					<< "The encountered dimensions of the CFX table are not the same throughout. They should be; for this to work.\n";

				}
				//Go through each one of the variables of interest
				if (var == 16) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure

					StaticEnergy_Table_Limits[0] = 10E20; //lower limit
					StaticEnergy_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								//cout<<line<<endl;
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i)); //bug fixed: detect end of table
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
								}
							}
							ThermoTables[i][j].StaticEnergy = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > StaticEnergy_Table_Limits[1]) {
								StaticEnergy_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < StaticEnergy_Table_Limits[0]) {
								StaticEnergy_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}

					}
					// //cout<<"Tables have been filled with Static Energy values "<<var<<endl;
				}
				if (var == 1) {
					//Fixed a bug: lines already skipped for var==1
					//for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					//for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure
					//cout<<"set_x "<<set_x<<endl;
					//cout<<"set_y "<<set_y<<endl;
					Enthalpy_Table_Limits[0] = 10E20;					//lower limit
					Enthalpy_Table_Limits[1] = 0;					//upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								//cout<<line<<endl;
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i)); //bug fixed: detect end of table
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
								}
							}
							ThermoTables[i][j].Enthalpy = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > Enthalpy_Table_Limits[1]) {
								Enthalpy_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < Enthalpy_Table_Limits[0]) {
								Enthalpy_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}

					}
					//cout<<"Tables have been filled with speed of specific enthalpy values "<<var<<endl;
				}
				if (var == 2) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure

					SoundSpeed2_Table_Limits[0] = 10E20; //lower limit
					SoundSpeed2_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								//cout<<line<<endl;
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i)); //bug fixed: detect end of table
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
									inp[z] = pow(inp[z], 2); //bug fixed: table features speed of sound; should be squared
								}
							}
							ThermoTables[i][j].SoundSpeed2 = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > SoundSpeed2_Table_Limits[1]) {
								SoundSpeed2_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < SoundSpeed2_Table_Limits[0]) {
								SoundSpeed2_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}

					}
					//cout<<"Tables have been filled with speed of sound values "<<var<<endl;
				}
				if (var == 5) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure

					Cp_Table_Limits[0] = 10E20; //lower limit
					Cp_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								//cout<<line<<endl;
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i)); //bug fixed: detect end of table
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
								}
							}
							ThermoTables[i][j].Cp = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > Cp_Table_Limits[1]) {
								Cp_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < Cp_Table_Limits[0]) {
								Cp_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}

					}
					//cout<<"Tables have been filled with isobaric heat capacity values "<<var<<endl;
				}
				if (var == 7) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure

					Entropy_Table_Limits[0] = 10E20; //lower limit
					Entropy_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								//cout<<line<<endl;
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i)); //bug fixed: detect end of table
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
								}
							}
							ThermoTables[i][j].Entropy = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > Entropy_Table_Limits[1]) {
								Entropy_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < Entropy_Table_Limits[0]) {
								Entropy_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}

					}
					//cout<<"Tables have been filled with entropy values "<<var<<endl;

				}
				if (var == 8) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure

					Mu_Table_Limits[0] = 10E20; //lower limit
					Mu_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								//cout<<line<<endl;
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i)); //bug fixed: detect end of table
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
								}
							}
							ThermoTables[i][j].Mu = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > Mu_Table_Limits[1]) {
								Mu_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < Mu_Table_Limits[0]) {
								Mu_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}

					}
					//cout<<"Tables have been filled with viscosity values "<<var<<endl;

				}
				if (var == 9) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure

					Kt_Table_Limits[0] = 10E20; //lower limit
					Kt_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								//cout<<line<<endl;
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i)); //bug fixed: detect end of table
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
								}
							}
							ThermoTables[i][j].Kt = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > Kt_Table_Limits[1]) {
								Kt_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < Kt_Table_Limits[0]) {
								Kt_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}

					}
					//cout<<"Tables have been filled with thermal conductivity values "<<var<<endl;

				}
				if (var == 10) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure

					dPdrho_e_Table_Limits[0] = 10E20; //lower limit
					dPdrho_e_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								//cout<<line<<endl;
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i)); //bug fixed: detect end of table
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
								}
							}
							ThermoTables[i][j].dPdrho_e = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > dPdrho_e_Table_Limits[1]) {
								dPdrho_e_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < dPdrho_e_Table_Limits[0]) {
								dPdrho_e_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}

					}
					//cout<<"Tables have been filled with specific dPdrho_e values "<<var<<endl;

				}
				if (var == 11) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure

					dPde_rho_Table_Limits[0] = 10E20; //lower limit
					dPde_rho_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								//cout<<line<<endl;
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i)); //bug fixed: detect end of table
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
								}
							}
							ThermoTables[i][j].dPde_rho = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > dPde_rho_Table_Limits[1]) {
								dPde_rho_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < dPde_rho_Table_Limits[0]) {
								dPde_rho_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}

					}
					//cout<<"Tables have been filled with specific dPde_rho values "<<var<<endl;
				}
				if (var == 12) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure

					dTdrho_e_Table_Limits[0] = 10E20; //lower limit
					dTdrho_e_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								//cout<<line<<endl;
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i)); //bug fixed: detect end of table
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
								}
							}
							ThermoTables[i][j].dTdrho_e = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > dTdrho_e_Table_Limits[1]) {
								dTdrho_e_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < dTdrho_e_Table_Limits[0]) {
								dTdrho_e_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}

					}
					// //cout<<"Tables have been filled with specific dTdrho_e values "<<var<<endl;
				}
				if (var == 13) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure

					dTde_rho_Table_Limits[0] = 10E20; //lower limit
					dTde_rho_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								//cout<<line<<endl;
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i)); //bug fixed: detect end of table
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
								}
							}
							ThermoTables[i][j].dTde_rho = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > dTde_rho_Table_Limits[1]) {
								dTde_rho_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < dTde_rho_Table_Limits[0]) {
								dTde_rho_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}

					}
					// //cout<<"Tables have been filled with specific dTde_rho values "<<var<<endl;
				}
				if (var == 15) {
					for (int k = 0; k < ceil(float(set_x) / 10.0); k++)
						getline(table, line); //skip density
					for (int k = 0; k < ceil(float(set_y) / 10.0); k++)
						getline(table, line); //skip pressure

					Temperature_Table_Limits[0] = 10E20; //lower limit
					Temperature_Table_Limits[1] = 0; //upper limit

					su2double inp[10];

					for (int j = 0; j < set_y; j++) {
						for (int i = 0; i < set_x; i++) {
							if ((j * set_x + i) % 10 == 0) {
								getline(table, line);
								//cout<<line<<endl;
								istringstream in(line);
								var_steps = 10;
								if (((set_x * set_y) - (j * set_x + i)) < 10)
									var_steps = ((set_x * set_y) - (j * set_x + i)); //bug fixed: detect end of table
								for (int z = 0; z < var_steps; z++) {
									in >> inp[z];
								}
							}
							ThermoTables[i][j].Temperature = inp[(j * set_x + i) % 10];
							if (inp[(j * set_x + i) % 10] > Temperature_Table_Limits[1]) {
								Temperature_Table_Limits[1] = inp[(j * set_x + i) % 10];
							}
							if (inp[(j * set_x + i) % 10] < Temperature_Table_Limits[0]) {
								Temperature_Table_Limits[0] = inp[(j * set_x + i) % 10];
							}
						}

					}
					// //cout<<"Tables have been filled with Temperature values "<<var<<endl;
				}

			}
		}
	}
	Table_Density_Stations = set_x;
	Table_Pressure_Stations = set_y;
	table.close();
}

