#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cassert>
#include "stdlib.h"
#include "stdio.h"
#include "LUT.hpp"

using namespace std;
CThermoList::CThermoList(){

	StaticEnergy = 0.0; //yes
	Entropy      = 0.0; //yes
	Density      = 0.0; //yes
	Pressure     = 0.0; //yes
	SoundSpeed2  = 0.0; //yes
	Temperature  = 0.0; //yes
	dPdrho_e     = 0.0; //yes
	dPde_rho     = 0.0; //yes
	dTdrho_e     = 0.0; //yes
	dTde_rho     = 0.0; //yes
	Cp           = 0.0; //yes
	Mu 			 = 0.0; //yes
	dmudrho_T    = 0.0; //no
	dmudT_rho    = 0.0; //no
	Kt           = 0.0; //yes
	dktdrho_T    = 0.0; //no
	dktdT_rho    = 0.0; //no

}
void CThermoList::CTLprint()
{
	cout<<"StaticEnergy:"<<StaticEnergy<<endl;
	cout<<"Entropy     :"<<Entropy<<endl;
	cout<<"Density     :"<<Density<<endl;
	cout<<"Pressure    :"<<Pressure<<endl;
	cout<<"SoundSpeed2 :"<<SoundSpeed2<<endl;
	cout<<"Temperature :"<<Temperature<<endl;
	cout<<"dPdrho_e    :"<<dPdrho_e<<endl;
	cout<<"dPde_rho    :"<<dPde_rho<<endl;
	cout<<"dTdrho_e    :"<<dTdrho_e<<endl;
	cout<<"dTde_rho    :"<<dTde_rho<<endl;
	cout<<"Cp          :"<<Cp<<endl;
	cout<<"Mu          :"<<Mu<<endl;
	cout<<"dmudrho_T   :"<<dmudrho_T<<endl;
	cout<<"dmudT_rho   :"<<dmudT_rho<<endl;
	cout<<"Kt          :"<<Kt<<endl;
	cout<<"dktdrho_T   :"<<dktdrho_T<<endl;
	cout<<"dktdT_rho   :"<<dktdT_rho<<endl;
}

CThermoList::~CThermoList(){

}

CLookUpTable::CLookUpTable() {
	ThermoTables = NULL;
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			coeff[i][j] = 0;
		}
	}

	iIndex = -1;
	jIndex = -1;
	p_dim  = 0;
	rho_dim  = 0;

}


CLookUpTable::CLookUpTable(char *Filename) {
	ThermoTables = NULL;
	try
	{
		TableLoadCFX(Filename);
		for (int i=0; i<3; i++)
		{
			for (int j=0; j<3; j++)
			{
				coeff[i][j] = 0;
			}
		}
		iIndex = -1;//negative number means it hasn't been preset yet
		jIndex = -1;//same
		cout<<"p_dim  : "<<p_dim<<endl;
		cout<<"rho_dim: "<<rho_dim<<endl;
	}
	catch (exception& e)
	{
		cerr<< e.what() << '\n';
	}
}



CLookUpTable::~CLookUpTable(void) {

}
void CLookUpTable::SearchKD_Tree (su2double thermo1, su2double thermo2,  char* thermoPair){

	su2double BestDist;
	su2double RunVal;
	unsigned int LowerI = 0;
	unsigned int UpperI = ceil(rho_dim/2);
	unsigned int LowerJ = 0;
	unsigned int UpperJ = ceil(p_dim/2);
	unsigned int leafsize = rho_dim*p_dim;
	unsigned int depth = 0;
	cout<<"Here"<<endl;
	while (leafsize>2)
	{
		leafsize = (UpperI-LowerI)*(UpperJ-LowerJ);
		depth++;
		cout<<RunVal<<endl;

		if (thermoPair == "RHOE")
		{
			cout<<"Here"<<endl;
			if ((depth%2)==0 and ((UpperI-LowerI)>1))
			{
				RunVal = ThermoTables[UpperI][UpperJ].Density;

			}
			else if ((depth%2)!=0 and ((UpperJ-LowerJ)>1))
			{
				RunVal = ThermoTables[UpperI][UpperJ].StaticEnergy;

			}
		}
		else if (thermoPair=="PT") cout<<endl;


		else if (thermoPair=="PRHO")cout<<endl;
		cout<<"I "<<LowerI<<" "<<UpperI<<endl;
		cout<<"J "<<LowerJ<<" "<<UpperJ<<endl;
		cout<<"Here---"<<RunVal<<endl;
		if ((depth%2)==0 and ((UpperI-LowerI)>1))
		{
			if (RunVal>thermo1)
			{
				UpperI = LowerI + ceil((UpperI-LowerI)/2);
			}
			else if (RunVal<thermo1)
			{
				LowerI = UpperI;
				UpperI = rho_dim;
			}

		}
		else if ((depth%2)!=0 and ((UpperJ-LowerJ)>1))
		{
			if (RunVal>thermo2)
			{
				UpperJ = LowerJ + ceil((UpperJ-LowerJ)/2);
			}
			else if (RunVal<thermo2)
			{
				LowerJ = UpperJ;
				UpperJ = p_dim;

			}

		}

		cout<<"I "<<LowerI<<" "<<UpperI<<endl;
		cout<<"J "<<LowerJ<<" "<<UpperJ<<endl;

	}


	//If lower than leafsize, transition to brute force
	int xmax, ymax;
	su2double D_rho, g_rho;
	su2double D_e, g_e;
	su2double e = thermo2;
	su2double rho = thermo1;
	D_rho = ThermoTables[LowerI][LowerJ].Density-rho;
	D_e = ThermoTables[LowerI][LowerJ].StaticEnergy-e;

	cout<<"Here4"<<endl;
	for (int j=0; j<(UpperJ-LowerJ); j++)
	{
		for (int i=0; i<(UpperI-LowerI); i++)
		{
			g_rho = (ThermoTables[i][j].Density-rho);
			g_e = (ThermoTables[i][j].StaticEnergy-e);
			if ((pow(g_rho/rho,2)+pow(g_e/e,2))<(pow(D_rho/rho,2)+pow(D_e/e,2)))
			{
				xmax = i;
				D_rho = g_rho;
				ymax = j;
				D_e   = g_e;
			}
		}
	}

	//Detect if it is in on a boundary:
	if(ymax==(0)) ymax++;
	if(xmax==(0)) xmax++;

	//Detect if point is not bottom corner of simplex
	if(D_e<0) ymax++;
	if(D_rho<0) xmax++;
	cout<<'KD search found'<<endl;
	cout<<xmax<<endl;
	cout<<ymax<<endl;

}


void CLookUpTable::SearchZigZag (su2double thermo1, su2double thermo2,  unsigned long thermoPair ){

	switch(thermoPair)
	{
	case 'RHOE':


		break;
	case 'PT':


		break;
	case 'PRHO':


		break;
	}

}

void CLookUpTable::SetTDState_rhoe (su2double rho, su2double e ) {
	try
	{
		//Check if inputs are in range (not a complete 2D check)
		if ((rho>Density_limits[1]) or (rho<Density_limits[0]))
		{
			throw runtime_error(string("RHOE Input Density out of bounds, using closest fit"));
		}
		if ((e>StaticEnergy_limits[1]) or (e<StaticEnergy_limits[0]))
		{	cout<<StaticEnergy_limits[0]<<" "<<StaticEnergy_limits[1]<<" "<<e<<endl;
		throw runtime_error(string("RHOE Input StaticEnergy out of bounds, using closest fit"));
		}
		cout<<endl<<"Rho desired : "<<rho<<std::endl;
		cout<<"E desired   : "<<e<<std::endl;

		su2double RunVal;
		unsigned  int CFL    = 2;
		unsigned int LowerI;
		unsigned int UpperI;
		unsigned int LowerJ;
		unsigned int UpperJ;
		//Restart search from previously used index if it exists
		if (jIndex<0) LowerJ = 0; UpperJ = ceil(p_dim/2);
		if (jIndex>=0) LowerJ = jIndex; UpperJ = jIndex+floor(p_dim/4.0);


		//Detemine the I index: assume rho is equispaced
		LowerI = floor((rho-Density_limits[0])/(Density_limits[1]-Density_limits[0])*(rho_dim-1));
		UpperI = LowerI + 1;

		//Determine the J index (for StaticEnergy), (Density invariant with j)
		while(UpperJ-LowerJ>1)
		{
			//Load the value at the upper bound
			RunVal = ThermoTables[UpperI][UpperJ].StaticEnergy;
			if (RunVal>e)
			{
				UpperJ = LowerJ + ceil((UpperJ-LowerJ)/CFL);
			}
			else if (RunVal<e)
			{
				int dif;
				dif   = UpperJ + ceil((UpperJ-LowerJ)/CFL);
				LowerJ = UpperJ;
				UpperJ = dif;
			}
			cout<<LowerJ<<"  "<<UpperJ<<endl;
		}

		iIndex = LowerI;
		jIndex = LowerJ;


		cout<<"Closest fit box :"<<endl;
		cout<<"Point i j :"<<endl;
		ThermoTables[iIndex][jIndex].CTLprint();
		cout<<"Point i+1 j :"<<endl;
		ThermoTables[iIndex+1][jIndex].CTLprint();
		cout<<"Point i j+1 :"<<endl;
		ThermoTables[iIndex][jIndex+1].CTLprint();
		cout<<"Point i+1 j+1 :"<<endl;
		ThermoTables[iIndex+1][jIndex+1].CTLprint();
		//Now use the closest fit box to interpolate



		su2double StaticEnergy = e;
		su2double Density      = rho;
		su2double Entropy      = Interp2D_lin(ThermoTables[iIndex-1][jIndex-1].Entropy,
				ThermoTables[iIndex-1][jIndex].Entropy,ThermoTables[iIndex][jIndex-1].Entropy,
				ThermoTables[iIndex][jIndex].Entropy );
		su2double Pressure     = Interp2D_lin(ThermoTables[iIndex-1][jIndex-1].Pressure,
				ThermoTables[iIndex-1][jIndex].Pressure,ThermoTables[iIndex][jIndex-1].Pressure,
				ThermoTables[iIndex][jIndex].Pressure );
		su2double SoundSpeed2  = Interp2D_lin(ThermoTables[iIndex-1][jIndex-1].SoundSpeed2,
				ThermoTables[iIndex-1][jIndex].SoundSpeed2,ThermoTables[iIndex][jIndex-1].SoundSpeed2,
				ThermoTables[iIndex][jIndex].SoundSpeed2 );
		su2double Temperature  = Interp2D_lin(ThermoTables[iIndex-1][jIndex-1].Temperature,
				ThermoTables[iIndex-1][jIndex].Temperature,ThermoTables[iIndex][jIndex-1].Temperature,
				ThermoTables[iIndex][jIndex].Temperature );
		su2double dPdrho_e     = Interp2D_lin(ThermoTables[iIndex-1][jIndex-1].dPdrho_e,
				ThermoTables[iIndex-1][jIndex].dPdrho_e,ThermoTables[iIndex][jIndex-1].dPdrho_e,
				ThermoTables[iIndex][jIndex].dPdrho_e );
		su2double dPde_rho     = Interp2D_lin(ThermoTables[iIndex-1][jIndex-1].dPde_rho,
				ThermoTables[iIndex-1][jIndex].dPde_rho,ThermoTables[iIndex][jIndex-1].dPde_rho,
				ThermoTables[iIndex][jIndex].dPde_rho );
		su2double dTdrho_e     =  Interp2D_lin(ThermoTables[iIndex-1][jIndex-1].dTdrho_e,
				ThermoTables[iIndex-1][jIndex].dTdrho_e,ThermoTables[iIndex][jIndex-1].dTdrho_e,
				ThermoTables[iIndex][jIndex].dTdrho_e );
		su2double dTde_rho     =  Interp2D_lin(ThermoTables[iIndex-1][jIndex-1].dTde_rho,
				ThermoTables[iIndex-1][jIndex].dTde_rho,ThermoTables[iIndex][jIndex-1].dTde_rho,
				ThermoTables[iIndex][jIndex].dTde_rho );
		su2double Cp           =  Interp2D_lin(ThermoTables[iIndex-1][jIndex-1].Cp,
				ThermoTables[iIndex-1][jIndex].Cp,ThermoTables[iIndex][jIndex-1].Cp,
				ThermoTables[iIndex][jIndex].Cp );
		su2double Mu 		   =  Interp2D_lin(ThermoTables[iIndex-1][jIndex-1].Mu,
				ThermoTables[iIndex-1][jIndex].Mu,ThermoTables[iIndex][jIndex-1].Mu,
				ThermoTables[iIndex][jIndex].Mu );
		su2double dmudrho_T    =  Interp2D_lin(ThermoTables[iIndex-1][jIndex-1].dmudrho_T,
				ThermoTables[iIndex-1][jIndex].dmudrho_T,ThermoTables[iIndex][jIndex-1].dmudrho_T,
				ThermoTables[iIndex][jIndex].dmudrho_T );
		su2double dmudT_rho    =  Interp2D_lin(ThermoTables[iIndex-1][jIndex-1].dmudT_rho,
				ThermoTables[iIndex-1][jIndex].dmudT_rho,ThermoTables[iIndex][jIndex-1].dmudT_rho,
				ThermoTables[iIndex][jIndex].dmudT_rho );
		su2double Kt           =  Interp2D_lin(ThermoTables[iIndex-1][jIndex-1].Kt,
				ThermoTables[iIndex-1][jIndex].Kt,ThermoTables[iIndex][jIndex-1].Kt,
				ThermoTables[iIndex][jIndex].Kt );
		su2double dktdrho_T    =  Interp2D_lin(ThermoTables[iIndex-1][jIndex-1].dktdrho_T,
				ThermoTables[iIndex-1][jIndex].dktdrho_T,ThermoTables[iIndex][jIndex-1].dktdrho_T,
				ThermoTables[iIndex][jIndex].dktdrho_T );
		su2double dktdT_rho    =  Interp2D_lin(ThermoTables[iIndex-1][jIndex-1].dktdT_rho,
				ThermoTables[iIndex-1][jIndex].dktdT_rho,ThermoTables[iIndex][jIndex-1].dktdT_rho,
				ThermoTables[iIndex][jIndex].dktdT_rho );

		if ((Density>Density_limits[1]) or (Density<Density_limits[0]))
		{

			throw runtime_error(string("RHOE Interpolated Density out of bounds, using closest fit"));

		}
		if ((Pressure>Pressure_limits[1]) or (Pressure<Pressure_limits[0]))
		{
			throw runtime_error(string("RHOE Interpolated Pressure out of bounds, using closest fit"));
		}
		cout<<"Interpolated fit:"<<endl;
		cout<<"StaticEnergy:"<<StaticEnergy<<endl;
		cout<<"Entropy     :"<<Entropy<<endl;
		cout<<"Density     :"<<Density<<endl;
		cout<<"Pressure    :"<<Pressure<<endl;
		cout<<"SoundSpeed2 :"<<SoundSpeed2<<endl;
		cout<<"Temperature :"<<Temperature<<endl;
		cout<<"dPdrho_e    :"<<dPdrho_e<<endl;
		cout<<"dPde_rho    :"<<dPde_rho<<endl;
		cout<<"dTdrho_e    :"<<dTdrho_e<<endl;
		cout<<"dTde_rho    :"<<dTde_rho<<endl;
		cout<<"Cp          :"<<Cp<<endl;
		cout<<"Mu          :"<<Mu<<endl;
		cout<<"dmudrho_T   :"<<dmudrho_T<<endl;
		cout<<"dmudT_rho   :"<<dmudT_rho<<endl;
		cout<<"Kt          :"<<Kt<<endl;
		cout<<"dktdrho_T   :"<<dktdrho_T<<endl;
		cout<<"dktdT_rho   :"<<dktdT_rho<<endl;
	}
	catch (exception& e)
	{
		cerr<<'\n'<< e.what() << '\n';
	}

}


void CLookUpTable::SetTDState_PT (su2double P, su2double T ) {

}

void CLookUpTable::SetTDState_Prho (su2double P, su2double rho ) {


}

void CLookUpTable::SetEnergy_Prho (su2double P, su2double rho ) {


}

void CLookUpTable::SetTDState_hs (su2double h, su2double s ) {

}

void CLookUpTable::SetTDState_Ps (su2double P, su2double s ){



}

void CLookUpTable::SetTDState_rhoT (su2double rho, su2double T ) {



}
void CLookUpTable::Interp2D_SingleSkewCoeff(std::string s)
{

	return;
}

void CLookUpTable::Interp2D_ArbitrarySkewCoeff(std::string s)
{
	su2double dx10,	dx01, dx11, dy10, dy01, dy11, f00, f10,	f11, f01;
	su2double A[3][3];
	su2double I[3][3];
	su2double c;
	//Load int the coordinates of the qudrilateral (values relative to i,j)

	//Setup the LHM matrix for the interpolation
	A[0][0] = dx10;
	A[0][1] = dy10;
	A[0][2] = dx10*dy10;
	A[1][0] = dx01;
	A[1][1] = dy01;
	A[1][2] = dx01*dy01;
	A[2][0] = dx11;
	A[2][1] = dy11;
	A[2][2] = dx11*dy11;
	I[0][0] = 1;
	I[1][1] = 1;
	I[2][2] = 1;

	//Compute inverse of LHM using Gaussian elimination
	//Reduced Echelon form of the LHM
	if (A[0,0] != 0)
	{
		c = A[1][0]/A[0][0];
		I[1][0] = I[1][0] -I[0][0]*c;
		for (int i=0; i<3; i++)
		{
			A[1][i] = A[1][i] -A[0][i]*c;
		}
		c =A[2][0]/A[0][0];
		I[2][0] = I[2][0] -I[0][0]*c;

		for (int i=0; i<3; i++)
		{
			A[2][i] = A[2][i] -A[0][i]*c;
		}
	}

	if (A[1][1] != 0)
	{
		c = A[2][1]/A[1][1];
		for (int i=0; i<2; i++)
		{
			I[2][i] = I[2][i] -I[1][i]*c;
		}
		for (int i=0; i<3; i++)
			A[2][i] = A[2][i] -A[1][i]*c;
	}
	//Reduced reduced Echelon form of LHM
	if (A[2][2] != 0)
	{

		for (int i=0; i<3; i++)
		{
			I[1][i] = I[1][i] -I[2][i]*A[1][2]/A[2][2];
		}
		A[1][2] = 0;
		for (int i=0; i<3; i++)
		{
			I[0][i] = I[0][i] -I[2][i]*A[0][2]/A[2][2];
		}
		A[0][2] = 0;
	}
	if (A[1][1] != 0)
	{
		for (int i=0; i<3; i++)
			I[0][i] = I[0][i] -I[1][i]*A[0][1]/A[1][1];
	}
	A[0][1] = 0;


	if (A[0][0] != 0)
	{
		for (int i=0; i<3; i++)
		{
			I[0][i] = I[0][i]/A[0][0];
		}
		if (A[1][1] != 0)
		{
			for (int i=0; i<3; i++)
			{
				I[1][i] = I[1][i]/A[1][1];
			}
		}
		if (A[2][2] != 0)
		{
			for (int i=0; i<3; i++)
			{
				I[2][i] = I[2][i]/A[2][2];
			}

		}
	}
	return;
}

su2double CLookUpTable::Interp2D_lin(su2double aa, su2double ab, su2double ba, su2double bb){

}

void CLookUpTable::LUTprint(void)
{
	for (int i=0; i<rho_dim; i++)
	{
		for (int j=0; j<p_dim; j++)
		{
			ThermoTables[i][j].CTLprint();
		}
	}
}


void CLookUpTable::TableLoadCFX(char *filename){
	int N_PARAM = 0;
	int set_x = 0;
	int set_y = 0;
	int var_steps = 0;
	int var_scanned=0;

	string line;
	string value;

	ifstream table (filename);
	assert(table.is_open());
	cout<<"Looking for number of parameters"<<endl;
	while ( getline(table,line) )
	{
		unsigned int found;
		found = line.find("$$PARAM");
		if (found<10)
		{
			getline(table,line);
			istringstream in(line);
			in>>N_PARAM;
			N_PARAM++;
			cout<<"Number of parameters "<<N_PARAM<<endl;
		}
		for (int var=var_scanned; var<N_PARAM+1; var++)
		{
			string svar = static_cast<ostringstream*>( &(ostringstream() << var) )->str();
			found = line.find("$TABLE_"+svar);
			if (found<10)
			{
				var_scanned = var;
				cout<<found<<' '<<line<<endl;
				getline(table,line);
				istringstream in(line);
				int x,y;
				in>>x>>y;
				if (var==1)
				{
					ThermoTables = new CThermoList*[x];
					for (int i=0; i<x; i++)
					{
						ThermoTables[i] = new CThermoList[y];
					}
					set_x = x;
					set_y = y;
					cout<<"Tables have been allocated"<<var<<endl;

					//Fill in the densities
					Density_limits[0] = 10E15;
					Density_limits[1] = 0;

					su2double* vD = new su2double[set_x];

					getline(table,line);
					cout<<line<<endl;
					for (int k =0; k<ceil(float(set_x)/10.0);k++)
					{
						istringstream inD(line);
						var_steps = 10;
						if ((set_x-k*10)<10) var_steps = (set_x-k*10);
						for (int i =0; i<var_steps; i++)
						{
							inD>>vD[10*k+i];
							if (vD[10*k+i]>Density_limits[1])
							{
								Density_limits[1]=vD[10*k+i];
							}
							if (vD[10*k+i]<Density_limits[0])
							{
								Density_limits[0]=vD[10*k+i];
							}
						}
					}
					for(int i =0; i<set_x; i++)
					{
						for (int j =0; j<set_y; j++)
						{
							ThermoTables[i][j].Density = vD[i];
						}
					}
					delete vD;


					//Fill in the pressures
					su2double* vP = new su2double[set_y];
					Pressure_limits[0] = 10E15; //lower limit
					Pressure_limits[1] = 0; //upper limit

					for (int k =0; k<ceil(float(set_y)/10.0);k++)
					{

						getline(table,line);
						cout<<line<<endl;
						istringstream inP(line);
						if ((set_y-k*10)<10) var_steps = (set_y-k*10);
						for (int j =0; j<var_steps; j++)
						{
							inP>>vP[10*k+j];
							if (vP[10*k+j]>Pressure_limits[1])
							{
								Pressure_limits[1]=vP[10*k+j];
							}
							if (vP[10*k+j]<Pressure_limits[0])
							{
								Pressure_limits[0]=vP[10*k+j];
							}
						}
					}
					for(int i =0; i<set_x; i++)
					{
						for (int j =0; j<set_y; j++)
						{
							ThermoTables[i][j].Pressure = vP[j];
						}
					}
					delete vP;
					cout<<"Tables have been filled with D and P values "<<var<<endl;

				}
				// Check that additional tables all adhere to the same x,y dimensions, otherwise throw an error
				else if (x != set_x && y!=set_y)
				{
					throw "The encountered dimensions of the CFX table are not the same throughout. They should be; for this to work.";

				}
				//Go through each one of the variables of interest
				if(var==16)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure

					StaticEnergy_limits[0] = 10E20;//lower limit
					StaticEnergy_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].StaticEnergy = inp[i];
							if (inp[i]>StaticEnergy_limits[1])
							{
								StaticEnergy_limits[1]= inp[i];
							}
							if (inp[i]<StaticEnergy_limits[0])
							{
								StaticEnergy_limits[0]=inp[i];
							}
						}

					}
					cout<<"Tables have been filled with Static Energy values "<<var<<endl;
				}
				if(var==2)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure

					SoundSpeed2_limits[0] = 10E20;//lower limit
					SoundSpeed2_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].SoundSpeed2 = inp[i];
							if (inp[i]>SoundSpeed2_limits[1])
							{
								SoundSpeed2_limits[1]= inp[i];
							}
							if (inp[i]<SoundSpeed2_limits[0])
							{
								SoundSpeed2_limits[0]=inp[i];
							}
						}

					}
					cout<<"Tables have been filled with speed of sound values "<<var<<endl;
				}
				if(var==5)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure


					Cp_limits[0] = 10E20;//lower limit
					Cp_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].Cp = inp[i];
							if (inp[i]>Cp_limits[1])
							{
								Cp_limits[1]= inp[i];
							}
							if (inp[i]<Cp_limits[0])
							{
								Cp_limits[0]=inp[i];
							}
						}

					}
					cout<<"Tables have been filled with isobaric heat capacity values "<<var<<endl;
				}
				if(var==7)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure

					Entropy_limits[0] = 10E20;//lower limit
					Entropy_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].Entropy = inp[i];
							if (inp[i]>Entropy_limits[1])
							{
								Entropy_limits[1]= inp[i];
							}
							if (inp[i]<Entropy_limits[0])
							{
								Entropy_limits[0]=inp[i];
							}
						}

					}
					cout<<"Tables have been filled with entropy values "<<var<<endl;

				}
				if(var==8)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure

					Mu_limits[0] = 10E20;//lower limit
					Mu_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].Mu = inp[i];
							if (inp[i]>Mu_limits[1])
							{
								Mu_limits[1]= inp[i];
							}
							if (inp[i]<Mu_limits[0])
							{
								Mu_limits[0]=inp[i];
							}
						}

					}
					cout<<"Tables have been filled with viscosity values "<<var<<endl;

				}
				if(var==9)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure

					Kt_limits[0] = 10E20;//lower limit
					Kt_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].Kt = inp[i];
							if (inp[i]>Kt_limits[1])
							{
								Kt_limits[1]= inp[i];
							}
							if (inp[i]<Kt_limits[0])
							{
								Kt_limits[0]=inp[i];
							}
						}

					}
					cout<<"Tables have been filled with thermal conductivity values "<<var<<endl;

				}
				if(var==10)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure


					dPdrho_e_limits[0] = 10E20;//lower limit
					dPdrho_e_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].dPdrho_e = inp[i];
							if (inp[i]>dPdrho_e_limits[1])
							{
								dPdrho_e_limits[1]= inp[i];
							}
							if (inp[i]<dPdrho_e_limits[0])
							{
								dPdrho_e_limits[0]=inp[i];
							}
						}

					}
					cout<<"Tables have been filled with specific dPdrho_e values "<<var<<endl;

				}
				if(var==11)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure


					dPde_rho_limits[0] = 10E20;//lower limit
					dPde_rho_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].dPde_rho = inp[i];
							if (inp[i]>dPde_rho_limits[1])
							{
								dPde_rho_limits[1]= inp[i];
							}
							if (inp[i]<dPde_rho_limits[0])
							{
								dPde_rho_limits[0]=inp[i];
							}
						}

					}
					cout<<"Tables have been filled with specific dPde_rho values "<<var<<endl;
				}
				if(var==12)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure

					dTdrho_e_limits[0] = 10E20;//lower limit
					dTdrho_e_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].dTdrho_e = inp[i];
							if (inp[i]>dTdrho_e_limits[1])
							{
								dTdrho_e_limits[1]= inp[i];
							}
							if (inp[i]<dTdrho_e_limits[0])
							{
								dTdrho_e_limits[0]=inp[i];
							}
						}

					}
					cout<<"Tables have been filled with specific dTdrho_e values "<<var<<endl;
				}
				if(var==13)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure


					dTde_rho_limits[0] = 10E20;//lower limit
					dTde_rho_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].dTde_rho = inp[i];
							if (inp[i]>dTde_rho_limits[1])
							{
								dTde_rho_limits[1]= inp[i];
							}
							if (inp[i]<dTde_rho_limits[0])
							{
								dTde_rho_limits[0]=inp[i];
							}
						}

					}
					cout<<"Tables have been filled with specific dTde_rho values "<<var<<endl;
				}
				if(var==15)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure


					Temperature_limits[0] = 10E20;//lower limit
					Temperature_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].Temperature = inp[i];
							if (inp[i]>Temperature_limits[1])
							{
								Temperature_limits[1]= inp[i];
							}
							if (inp[i]<Temperature_limits[0])
							{
								Temperature_limits[0]=inp[i];
							}
						}

					}
					cout<<"Tables have been filled with Temperature values "<<var<<endl;
				}

			}
		}
	}
	rho_dim = set_x;
	p_dim = set_y;
	table.close();
}










