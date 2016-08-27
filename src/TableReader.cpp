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
#include "time.h"


using namespace std;

int main() {

	//char* tablefile= (char*) "TableFile.dat";
	//char* tablefile= (char*) "AIR/air.rgp";
	char* gridfile = (char*) "mesh.dat";
	char* timefile = (char*) "time.dat";
	char* tablefile= (char*) "new_mesh/dT/Complete/lutmesh.tec";
	char* outfile;
	double duration;
	double *rho= new double[400000];
	double *e  = new double[400000];
	double *h  = new double[400000];
	double *s  = new double[400000];
	double *P  = new double[400000];
	double *T  = new double[400000];
	string line;
	fstream ft;
	fstream fs;
	int i;

	ft.open(timefile, fstream::out |fstream::trunc);
	ft<<"";
	ft.close();
	ft.open(timefile, fstream::app);
	ft.precision(17);
	assert(ft.is_open());
//
 	CLookUpTable LUT2 = CLookUpTable(tablefile);
//
	fs.open(gridfile, fstream::out |fstream::trunc);
	fs<<"";
	fs.close();
	LUT2.LookUpTable_Print_To_File(gridfile);
////
////	//LUT2.SetTDState_rhoe(85, 570000);
////
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
//
//wipe the outfile before writing to it
	outfile = (char*) "rhoe_out.dat";
	fs.open(outfile, fstream::out |fstream::trunc);
	fs.close();
	for (int j=0; j<i;j++)
	{
			LUT2.SetTDState_rhoe(rho[j], e[j]);
			LUT2.RecordState(outfile);
	}
	clock_t rhoe_start = clock();
	for (int j=0; j<i;j++)
	{
		LUT2.SetTDState_rhoe(rho[j], e[j]);
	}
	duration = ((su2double)clock()-(su2double)rhoe_start)/((su2double)CLOCKS_PER_SEC);
	ft<<"rhoe time: "<<duration<<endl;
////
////
////
//cout<<"--------------------------------------------------------------\n";
//	cout<<"PT"<<endl;
//
//	fs.open("PT_in.dat", fstream::in);
//	i=0;
//	while(getline(fs,line))
//	{
//		istringstream in(line);
//		in>>P[i];
//		in>>T[i];
//		i++;
//	}
//	fs.close();
//	//wipe the outfile before writing to it
//	outfile = (char*) "PT_out.dat";
//	fs.open(outfile, fstream::out |fstream::trunc);
//	fs.close();
//	for (int j=0; j<i;j++)
//		{
//			cout<<j<<endl;
//			LUT2.SetTDState_PT(P[j], T[j]);
//			LUT2.RecordState(outfile);
//		}
//	clock_t PT_start = clock();
//	for (int j=0; j<i;j++)
//	{
//		LUT2.SetTDState_PT(P[j], T[j]);
//	}
//	duration = ((su2double)clock()-(su2double)PT_start)/((su2double)CLOCKS_PER_SEC);
//	ft<<"PT time: "<<duration<<endl;
////	//LUT2.SetTDState_PT(5520052, 270);
////	//LUT2.reset_Restart();
////
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
	outfile = (char*) "Prho_out.dat";
	fs.open(outfile, fstream::out |fstream::trunc);
	fs.close();
	for (int j=0; j<i;j++)
		{
			LUT2.SetTDState_Prho(P[j], rho[j]);
			LUT2.RecordState(outfile);
		}
//	clock_t prho_start = clock();
//	for (int j=0; j<i;j++)
//	{
//		LUT2.SetTDState_Prho(P[j], rho[j]);
//	}
//	duration = ((su2double)clock()-(su2double)prho_start)/((su2double)CLOCKS_PER_SEC);
//	ft<<"Prho time: "<<duration<<endl;
//
//	//LUT2.SetTDState_Prho(5521052, 80);
//
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
	fs.open(outfile, fstream::out |fstream::trunc);
	fs.close();
	for (int j=0; j<i;j++)
		{
			LUT2.SetTDState_rhoT(rho[j], T[j]);
			LUT2.RecordState(outfile);
		}
//	clock_t rhoT_start = clock();
//	for (int j=0; j<i;j++)
//	{
//		LUT2.SetTDState_rhoT(rho[j], T[j]);
//	}
//	duration = ((su2double)clock()-(su2double)rhoT_start)/((su2double)CLOCKS_PER_SEC);
//	ft<<"rhoT time: "<<duration<<endl;
//	//LUT2.SetTDState_rhoT(90, 80);
//	//LUT2.reset_Restart();
//
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
	fs.open(outfile, fstream::out |fstream::trunc);
	fs.close();
	for (int j=0; j<i;j++)
		{
			LUT2.SetTDState_Ps(P[j], s[j]);
			LUT2.RecordState(outfile);
		}
//	clock_t Ps_start = clock();
//	for (int j=0; j<i;j++)
//	{
//		LUT2.SetTDState_Ps(P[j], s[j]);
//	}
//	//LUT2.SetTDState_Ps(5400000, 2200);
//	duration = ((su2double)clock()-(su2double)Ps_start)/((su2double)CLOCKS_PER_SEC);
//	ft<<"Ps time: "<<duration<<endl;
////
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
	fs.open(outfile, fstream::out |fstream::trunc);
	fs.close();
	for (int j=0; j<i;j++)
		{
			LUT2.SetTDState_hs(h[j], s[j]);
			LUT2.RecordState(outfile);
		}
//	clock_t hs_start = clock();
//	for (int j=0; j<i;j++)
//	{
//		LUT2.SetTDState_hs(h[j], s[j]);
//	}
//	duration = ((su2double)clock()-(su2double)hs_start)/((su2double)CLOCKS_PER_SEC);
//	ft<<"hs time: "<<duration<<endl;
//	//LUT2.SetTDState_hs(531782, 2092);
//	//LUT2.reset_Restart();
//
	cout<<"------------------------------END--------------------------------\n";
	ft.close();
	delete[] rho;
	delete[] e;
	delete[] h;
	delete[] s;
	delete[] P;
	delete[] T;
	return 0;
	}

