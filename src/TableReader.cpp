//============================================================================
// Name        : TableReader.cpp
// Author      : Matej Kosec
// Version     :
// Copyright   : 
// Description :
//============================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cassert>

#include "TableReader.hpp"
#include "LUT.cpp"
#include "LUT.hpp"

using namespace std;

int main() {
	double* p=NULL;
	int N_PARAM = 0;
	double* r=NULL;
	int set_x = 0;
	int set_y = 0;
	int var_scanned=0;


	CThermoList**ThermoTables;
    string line;
	string value;

	ifstream table ("/home/matej/SU2/DEV/cfxlut/Code/CO2.rgp");
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
					getline(table,line);
					istringstream inD(line);
					double* vD = new double[set_x];
					for (int i =0; i<set_x; i++)
					{
						inD>>vD[i];
					}
					for(int j =0; j<set_y; j++)
									{
						for (int i =0; i<set_x; i++)
							{
							ThermoTables[i][j].Density = vD[j];
							}
						}

				    delete vD;

				    //Fill in the pressures
				    getline(table,line);
				    istringstream inP(line);
					double* vP = new double[set_y];
					for (int j =0; j<set_y; j++)
					{
						inP>>vP[j];
					}
					for(int i =0; i<set_y; i++)
									{
						for (int j =0; j<set_y; j++)
							{
							ThermoTables[i][j].Pressure = vP[i];
							}
						}
					delete vP;
					cout<<"Tables have been filled with D and P"<<var<<endl;

				}
				// Check that additional tables all adhere to the same x,y dimensions, otherwise throw an error
				else if (x != set_x && y!=set_y)
				{
					throw "The encountered dimensions of the CFX table are not the same throughout. They should be for this to work.";
				}
				//Go through each one of the variables of interest
				if(var==2)
				{
					for (int i =0; i<set_x; i++)
					{
						getline(table,line);
						istringstream in(line);
						for(int j =0; j<set_y; j++)
						{
						in>>ThermoTables[i][j].SoundSpeed2;
						}

					}
				cout<<"Tables have been filled with speed of sound values"<<var<<endl;
				}
				if(var==5)
				{
					for (int i =0; i<set_x; i++)
					{
						getline(table,line);
						istringstream in(line);
						for(int j =0; j<set_y; j++)
						{
						in>>ThermoTables[i][j].Cp;
						}
					}
				cout<<"Tables have been filled with isobaric heat capacity values"<<var<<endl;
				}
				if(var==7)
				{
					for (int i =0; i<set_x; i++)
					{
						getline(table,line);
						istringstream in(line);
						for(int j =0; j<set_y; j++)
						{
						in>>ThermoTables[i][j].Entropy;
						}
					}
				cout<<"Tables have been filled with entropy values "<<var<<endl;
				}
				if(var==8)
				{
					for (int i =0; i<set_x; i++)
					{
						getline(table,line);
						istringstream in(line);
						for(int j =0; j<set_y; j++)
						{
						in>>ThermoTables[i][j].Mu;
						}
					}
				cout<<"Tables have been filled with viscosity values "<<var<<endl;
				}
				if(var==9)
				{
					for (int i =0; i<set_x; i++)
					{
						getline(table,line);
						istringstream in(line);
						for(int j =0; j<set_y; j++)
						{
						in>>ThermoTables[i][j].Kt;
						}

					}
				cout<<"Tables have been filled with thermal conductivity values "<<var<<endl;
				}

				if(var==10)
				{
					for (int i =0; i<set_x; i++)
					{
						getline(table,line);
						istringstream in(line);
						for(int j =0; j<set_y; j++)
						{
						in>>ThermoTables[i][j].dTdrho_e;
						}
					}
					cout<<"Tables have been filled with specific dPdrho_e values "<<var<<endl;
				}
				if(var==11)
				{
					for (int i =0; i<set_x; i++)
					{
						getline(table,line);
						istringstream in(line);
						for(int j =0; j<set_y; j++)
						{
						in>>ThermoTables[i][j].dPde_rho;
						}
					}
					cout<<"Tables have been filled with specific dPde_rho values "<<var<<endl;
				}
				if(var==12)
				{
					for (int i =0; i<set_x; i++)
					{
						getline(table,line);
						istringstream in(line);
						for(int j =0; j<set_y; j++)
						{
						in>>ThermoTables[i][j].dTdrho_e;
						}
					}
					cout<<"Tables have been filled with specific dTdrho_e values "<<var<<endl;
								}
				if(var==13)
				{
					for (int i =0; i<set_x; i++)
					{
						getline(table,line);
						istringstream in(line);
						for(int j =0; j<set_y; j++)
						{
						in>>ThermoTables[i][j].dTde_rho;
						}
					}
					cout<<"Tables have been filled with specific dTde_rho values "<<var<<endl;
				}

			}
			}
	    }
	    table.close();
	delete p;
	delete r;
	delete ThermoTables;
	return 0;
	}



