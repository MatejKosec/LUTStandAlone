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

	StaticEnergy = 0.0; //no
	Entropy      = 0.0; //yes
	Density      = 0.0; //yes
	Pressure     = 0.0; //yes
	SoundSpeed2  = 0.0; //yes
	Temperature  = 0.0; //no
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
	coeff = NULL;
	iIndex = 0;
	jIndex = 0;
	p_dim  = 0;
	rho_dim  = 0;

}


CLookUpTable::CLookUpTable(char *Filename) {
	ThermoTables = NULL;
	TableLoadCFX(Filename);
	coeff = new su2double[4];
	iIndex = 0;
	jIndex = 0;
	cout<<"p_dim  : "<<p_dim<<endl;
	cout<<"rho_dim: "<<rho_dim<<endl;
}



CLookUpTable::~CLookUpTable(void) {

}

void CLookUpTable::SearchThermoPair (su2double thermo1, su2double thermo2,  unsigned short thermoPair ){

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
	cout<<endl<<"Rho desired : "<<rho<<std::endl;
	cout<<"E desired   : "<<e<<std::endl;
	cout<<"Closest fit :"<<endl;

	SearchThermoRHOE(rho, e);
	ThermoTables[iIndex][jIndex].CTLprint();




	cout<<endl;

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

void CLookUpTable::SearchThermoRHOE (su2double rho, su2double e){
	int xmax= 0;
	int ymax= 0;
	su2double D_rho, g_rho;
	su2double D_e, g_e;
	su2double f1, f2, y1, y2;
	D_rho = ThermoTables[0][0].Density-rho;
	D_e = ThermoTables[0][0].StaticEnergy-e;

	for (int j=0; j<p_dim; j++)
	{
		for (int i=0; i<rho_dim; i++)
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

	f1 = (rho-ThermoTables[xmax][ymax].Density)/(ThermoTables[xmax][ymax].Density - ThermoTables[xmax][ymax-1].Density);
	y1 = ThermoTables[xmax-1][ymax-1].StaticEnergy + f1*(ThermoTables[xmax][ymax-1].StaticEnergy-ThermoTables[xmax-1][ymax-1].StaticEnergy);
	y2 = ThermoTables[xmax-1][ymax].StaticEnergy   + f1*(ThermoTables[xmax][ymax].StaticEnergy  -ThermoTables[xmax-1][ymax].StaticEnergy);
	f2 = (e-y1)/(y2-y1);

	coeff[0]  = f1;
	coeff[1]  = f2;
	coeff[2]  = y1;
	coeff[3]  = y2;
	cout<<"f1 "<<f1<<" f2 "<<f2<<" y1 "<<y1<<" y2 "<<y2<<endl;
	cout<<"jIndex "<<ymax<<" iIndex "<<xmax<<endl;
	jIndex = ymax;
	iIndex = xmax;
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

su2double CLookUpTable::Interp2D_lin(su2double aa, su2double ab, su2double ba, su2double bb){
	su2double z, z1,z2;
	z1 = aa + coeff[0]*(ba - aa);
	z2 = ab + coeff[0]*(bb - ab);
	z  = z1 + coeff[1]*(z2 - z1);
	return z;
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
				cout<<found<<' '<<' '<<line<<endl;
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
						}
					}
					for(int i =0; i<set_x; i++)
					{
						for (int j =0; j<set_y; j++)
						{
							ThermoTables[i][j].Density = vD[j];
						}
					}
					delete vD;


					//Fill in the pressures
					su2double* vP = new su2double[set_y];

					for (int k =0; k<ceil(float(set_y)/10.0);k++)
					{

						getline(table,line);
						cout<<line<<endl;
						istringstream inP(line);
						if ((set_y-k*10)<10) var_steps = (set_y-k*10);
						for (int j =0; j<var_steps; j++)
						{
							inP>>vP[10*k+j];
						}
					}
					for(int i =0; i<set_x; i++)
					{
						for (int j =0; j<set_y; j++)
						{
							ThermoTables[i][j].Pressure = vP[i];
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
				if(var==1)
				{

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							su2double enthalpy;
							su2double pressure;
							su2double density;

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
							pressure = ThermoTables[i][j].Pressure;
							density = 1.0/ThermoTables[i][j].Density;
							ThermoTables[i][j].StaticEnergy = inp[i]-pressure/density;
						}

					}
					cout<<"Tables have been filled with Static Energy values "<<var<<endl;
				}
				if(var==2)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure

					su2double inp[10];

					for(int j =0; j<set_y; j++)
					{
						for (int i =0; i<set_x; i++)
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
							ThermoTables[i][j].SoundSpeed2=inp[i];

						}
					}
					cout<<"Tables have been filled with speed of sound values "<<var<<endl;
				}
				if(var==5)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure


					su2double inp[10];

					for(int j =0; j<set_y; j++)
					{
						for (int i =0; i<set_x; i++)
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
							ThermoTables[i][j].Cp=inp[i];
						}

					}
					cout<<"Tables have been filled with isobaric heat capacity values "<<var<<endl;
				}
				if(var==7)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure

					su2double inp[10];

					for(int j =0; j<set_y; j++)
					{
						for (int i =0; i<set_x; i++)
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

							ThermoTables[i][j].Entropy=inp[i];

						}
					}
					cout<<"Tables have been filled with entropy values "<<var<<endl;

				}
				if(var==8)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure

					su2double inp[10];

					for(int j =0; j<set_y; j++)
					{
						for (int i =0; i<set_x; i++)
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

						ThermoTables[i][j].Mu=inp[i];

					}
				}
				cout<<"Tables have been filled with viscosity values "<<var<<endl;

			}
			if(var==9)
			{
				for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
				for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure

				su2double inp[10];

				for(int j =0; j<set_y; j++)
				{
					for (int i =0; i<set_x; i++)
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

						ThermoTables[i][j].Kt=inp[i];


					}
				}
				cout<<"Tables have been filled with thermal conductivity values "<<var<<endl;

			}
			if(var==10)
			{
				for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
				for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure


				su2double inp[10];

				for(int j =0; j<set_y; j++)
				{
					for (int i =0; i<set_x; i++)
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

						ThermoTables[i][j].dPdrho_e=inp[i];


					}
				}
				cout<<"Tables have been filled with specific dPdrho_e values "<<var<<endl;

			}
			if(var==11)
			{
				for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
				for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure


				su2double inp[10];

				for(int j =0; j<set_y; j++)
				{
					for (int i =0; i<set_x; i++)
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

						ThermoTables[i][j].dPde_rho=inp[i];


					}

				}
				cout<<"Tables have been filled with specific dPde_rho values "<<var<<endl;
			}
			if(var==12)
			{
				for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
				for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure

				su2double inp[10];

				for(int j =0; j<set_y; j++)
				{
					for (int i =0; i<set_x; i++)
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

						ThermoTables[i][j].dTdrho_e=inp[i];


					}
				}
				cout<<"Tables have been filled with specific dTdrho_e values "<<var<<endl;
			}
			if(var==13)
			{
				for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
				for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure


				su2double inp[10];

				for(int j =0; j<set_y; j++)
				{
					for (int i =0; i<set_x; i++)
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

						ThermoTables[i][j].dTde_rho=inp[i];

					}
				}
				cout<<"Tables have been filled with specific dTde_rho values "<<var<<endl;
			}

		}
	}
}
rho_dim = set_x;
p_dim = set_y;
table.close();
}










