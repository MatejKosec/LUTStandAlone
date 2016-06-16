//============================================================================
// Name        : TableReader.cpp
// Author      : Matej Kosec
// Version     :
// Copyright   : 
// Description :
//============================================================================

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include "stdio.h"
#include <cmath>
#include <memory>
#include <cassert>
#include "LUT.hpp"


using namespace std;

int main() {

	char* tablefile= (char*) "CO2.rgp";
	char* outfile;
	su2double rho[500];
	su2double e[500];
	su2double h[500];
	su2double s[500];
	su2double P[500];
	su2double T[500];
	string line;
	fstream fs;
	int i;

	CLookUpTable LUT2 = CLookUpTable(tablefile);

	cout<<"--------------------------------------------------------------\n";
	cout<<"RHOE \n";
	fs.open("rhoe_in.dat", fstream::in);
	i=0;
	while(getline(fs,line))
	{
		istringstream in(line);
		in>>rho[i];
		in>>e[i];
		i++;
	}
	fs.close();
	//wipe the outfile before writing to it
	outfile = (char*) "rhoe_out.dat";
	fs.open(outfile, fstream::trunc);
	for (int j=0; j<i;j++)
	{
		cout<<j<<endl;
		LUT2.SetTDState_rhoe(rho[j], e[j]);
		LUT2.RecordState(outfile);
	}

	cout<<"--------------------------------------------------------------\n";
	cout<<"PT"<<endl;

	fs.open("PT_in.dat", fstream::in);
	i=0;
	while(getline(fs,line))
	{
		istringstream in(line);
		in>>P[i];
		in>>T[i];
		i++;
	}
	fs.close();
	//wipe the outfile before writing to it
	outfile = (char*) "PT_out.dat";
	fs.open(outfile, fstream::trunc);
	for (int j=0; j<i;j++)
	{
		cout<<j<<endl;
		LUT2.SetTDState_PT(P[j], T[j]);
		LUT2.RecordState(outfile);
	}

	//LUT2.SetTDState_PT(5520052, 270);
	//LUT2.reset_Restart();

	cout<<"--------------------------------------------------------------\n";
	cout<<"Prho \n";
	fs.open("Prho_in.dat", fstream::in);
	i=0;
	while(getline(fs,line))
	{
		istringstream in(line);
		in>>P[i];
		in>>rho[i];
		i++;
	}
	fs.close();
	//wipe the outfile before writing to it
	outfile = (char*) "Phoe_out.dat";
	fs.open(outfile, fstream::trunc);
	for (int j=0; j<i;j++)
	{
		cout<<j<<endl;
		LUT2.SetTDState_Prho(P[j], rho[j]);
		LUT2.RecordState(outfile);
	}
	//LUT2.SetTDState_Prho(5521052, 80);

	cout<<"--------------------------------------------------------------\n";

	cout<<"rhoT"<<endl;
	fs.open("rhoT_in.dat", fstream::in);
	i=0;
	while(getline(fs,line))
	{
		istringstream in(line);
		in>>rho[i];
		in>>T[i];
		i++;
	}
	fs.close();
	//wipe the outfile before writing to it
	outfile = (char*) "rhoT_out.dat";
	fs.open(outfile, fstream::trunc);
	for (int j=0; j<i;j++)
	{
		cout<<j<<endl;
		LUT2.SetTDState_rhoT(rho[j], T[j]);
		LUT2.RecordState(outfile);
	}
	//LUT2.SetTDState_rhoT(90, 80);
	//LUT2.reset_Restart();

	cout<<"--------------------------------------------------------------\n";

	cout<<"Ps"<<endl;
	fs.open("Ps_in.dat", fstream::in);
	i=0;
	while(getline(fs,line))
	{
		istringstream in(line);
		in>>P[i];
		in>>s[i];
		i++;
	}
	fs.close();
	//wipe the outfile before writing to it
	outfile = (char*) "Ps_out.dat";
	fs.open(outfile, fstream::trunc);
	for (int j=0; j<i;j++)
	{
		cout<<j<<endl;
		LUT2.SetTDState_Ps(P[j], s[j]);
		LUT2.RecordState(outfile);
	}
	//LUT2.SetTDState_Ps(5400000, 2200);
	//LUT2.reset_Restart();

	cout<<"--------------------------------------------------------------\n";

	cout<<"hs"<<endl;
	fs.open("hs_in.dat", fstream::in);
	i=0;
	while(getline(fs,line))
	{
		istringstream in(line);
		in>>h[i];
		in>>s[i];
		i++;
	}
	fs.close();
	//wipe the outfile before writing to it
	outfile = (char*) "hs_out.dat";
	fs.open(outfile, fstream::trunc);
	for (int j=0; j<i;j++)
	{
		cout<<j<<endl;
		LUT2.SetTDState_hs(h[j], s[j]);
		LUT2.RecordState(outfile);
	}
	//LUT2.SetTDState_hs(531782, 2092);
	//LUT2.reset_Restart();

	cout<<"--------------------------------------------------------------\n";

	return 0;
}

