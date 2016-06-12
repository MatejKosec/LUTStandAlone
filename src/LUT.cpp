#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cassert>
#include "stdlib.h"
#include "stdio.h"
#include "LUT.hpp"
#include <iomanip>

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


CLookUpTable::CLookUpTable(char* Filename) {
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
void CLookUpTable::SearchKD_Tree (su2double thermo1, su2double thermo2,  string thermoPair){

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
	cout<<"KD search found"<<endl;
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
		//Check if inputs are in total range (necessary but not sufficient condition)
		if ((rho>Density_limits[1]) or (rho<Density_limits[0]))
		{
			throw runtime_error("RHOE Input Density out of bounds");
		}
		if ((e>StaticEnergy_limits[1]) or (e<StaticEnergy_limits[0]))
		{
			throw runtime_error("RHOE Input StaticEnergy out of bounds");
		}
		cout<<endl<<"rho desired : "<<rho<<endl;
		cout<<"e desired   : "<<e<<endl;

		su2double RunVal;
		unsigned int CFL    = 2;
		unsigned int LowerI;
		unsigned int UpperI;
		unsigned int LowerJ;
		unsigned int UpperJ;
		//Restart search from previously used index if it exists, else go to middle
		if (jIndex<0) LowerJ = 0;
		if (jIndex>=0) LowerJ = jIndex;
		if (jIndex<ceil(p_dim/2))
		{
			UpperJ = ceil(p_dim/2);
		}
		else UpperJ=p_dim;

		//Determine the I index: rho is equispaced (no restart)
		LowerI = floor((rho-Density_limits[0])/(Density_limits[1]-Density_limits[0])*(rho_dim-1));
		UpperI = LowerI + 1;

		//Determine the J index (for e), (e invariant with j)
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
		cout<<"i "<<LowerI<<"  "<<UpperI<<endl;
		cout<<"j "<<LowerJ<<"  "<<UpperJ<<endl;


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


		su2double x, y;
		x = rho - ThermoTables[iIndex][jIndex].Density;
		y = e - ThermoTables[iIndex][jIndex].StaticEnergy;

		//Determine interpolation coefficients
		Interp2D_ArbitrarySkewCoeff(x,y,"RHOE");
		cout<<"Interpolation matrix inverse \n";
		for (int j=0; j<3; j++)
		{
			cout<<setw(15)<<coeff[j][0]<<"   "<<coeff[j][1]<<"   "<<coeff[j][2]<<endl;
		}
		CThermoList interpolated;
		interpolated.StaticEnergy      = e;
		interpolated.Density           = rho ;
		interpolated.Entropy           = Interp2D_lin(x, y, "Entropy" );
		interpolated.Pressure          = Interp2D_lin(x, y, "Pressure" );
		interpolated.SoundSpeed2       = Interp2D_lin(x, y, "SoundSpeed2" );
		interpolated.Temperature       = Interp2D_lin(x, y, "Temperature" );
		interpolated.dPdrho_e          = Interp2D_lin(x, y, "dPdrho_e" );
		interpolated.dPde_rho          = Interp2D_lin(x, y, "dPde_rho" );
		interpolated.dTdrho_e          = Interp2D_lin(x, y, "dTdrho_e" );
		interpolated.dTde_rho          = Interp2D_lin(x, y, "dTde_rho" );
		interpolated.Cp                = Interp2D_lin(x, y, "Cp" );
		interpolated.Mu                = Interp2D_lin(x, y, "Mu" );
		interpolated.dmudrho_T         = Interp2D_lin(x, y, "dmudrho_T" );
		interpolated.dmudT_rho         = Interp2D_lin(x, y, "dmudT_rho" );
		interpolated.Kt                = Interp2D_lin(x, y, "Kt" );
		interpolated.dktdrho_T         = Interp2D_lin(x, y, "dktdrho_T" );
		interpolated.dktdT_rho         = Interp2D_lin(x, y, "dktdT_rho" );
		//Intermediate variables only needed for StandAlone version
		su2double Density = interpolated.Density;
		su2double Pressure = interpolated.Pressure;
		cout<<"Interpolated fit:"<<endl;
		interpolated.CTLprint ();
		if ((Density>Density_limits[1]) or (Density<Density_limits[0]))
		{
			throw runtime_error("RHOE Interpolated Density out of bounds");
		}
		if ((Pressure>Pressure_limits[1]) or (Pressure<Pressure_limits[0]))
		{
			throw runtime_error("RHOE Interpolated Pressure out of bounds");
		}

	}
	catch (exception& e)
	{
		cerr<<"\n"<< e.what() << "\n";
	}
}

void CLookUpTable::SetTDState_PT (su2double P, su2double T ) {
	try
	{
		//Check if inputs are in total range (necessary but not sufficient condition)
		if ((P>Pressure_limits[1]) or (P<Pressure_limits[0]))
		{
			throw runtime_error("PT Input Pressure out of bounds");
		}
		if ((T>Temperature_limits[1]) or (T<Temperature_limits[0]))
		{
			throw runtime_error("PT Input Temperature out of bounds");
		}
		cout<<endl<<"P desired : "<<P<<endl;
		cout<<"T desired   : "<<T<<endl;

		su2double RunVal;
		unsigned int CFL = 2;
		unsigned int LowerI;
		unsigned int UpperI;
		unsigned int LowerJ;
		unsigned int UpperJ;
		//Restart search from previously used index if it exists, else go to middle
		if (iIndex<0) LowerI = 0;
		if (iIndex>=0) LowerI = iIndex;
		if  (iIndex<ceil(rho_dim/2))
		{
			UpperI = ceil(rho_dim/2); //probably can be made more efficient (smaller square)
		}
		else UpperI = rho_dim;


		//Determine the J index: P is equispaced (no restart)
		LowerJ = floor((P-Pressure_limits[0])/(Pressure_limits[1]-Pressure_limits[0])*(p_dim-1));
		UpperJ = LowerJ + 1;

		//Determine the I index (for T)
		//Temperature is not necessarily monotonic (so use a variable to take care of overall sign)
		int grad;
		while(UpperI-LowerI>1)
		{
			grad = ThermoTables[UpperI][UpperJ].Temperature - ThermoTables[LowerI][UpperJ].Temperature;
			//Load the value at the upper bound
			RunVal = ThermoTables[UpperI][UpperJ].Temperature;
			if (grad*RunVal>grad*T)
			{
				UpperI = LowerI + ceil((UpperI-LowerI)/CFL);
			}
			else if (grad*RunVal<grad*T)
			{
				int dif;
				dif   = UpperI + ceil((UpperI-LowerI)/CFL);
				LowerI = UpperI;
				UpperI = dif;
			}
			cout<<LowerI<<"  "<<UpperI<<endl;
		}

		iIndex = LowerI;
		jIndex = LowerJ;
		cout<<"i "<<LowerI<<"  "<<UpperI<<endl;
		cout<<"j "<<LowerJ<<"  "<<UpperJ<<endl;

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

		su2double x, y;
		x = T - ThermoTables[iIndex][jIndex].Temperature;
		y = P - ThermoTables[iIndex][jIndex].Pressure;
		//Determine interpolation coefficients
		Interp2D_ArbitrarySkewCoeff(x,y,"PT");
		cout<<"Interpolation matrix inverse \n";
		for (int j=0; j<3; j++)
		{
			cout<<setw(15)<<coeff[j][0]<<"   "<<coeff[j][1]<<"   "<<coeff[j][2]<<endl;
		}
		CThermoList interpolated;
		interpolated.Temperature       = T;
		interpolated.Pressure          = P ;
		interpolated.StaticEnergy      = Interp2D_lin(x, y, "StaticEnergy" );
		interpolated.Entropy           = Interp2D_lin(x, y, "Entropy" );
		interpolated.Density           = Interp2D_lin(x, y, "Density" );
		interpolated.SoundSpeed2       = Interp2D_lin(x, y, "SoundSpeed2" );
		interpolated.dPdrho_e          = Interp2D_lin(x, y, "dPdrho_e" );
		interpolated.dPde_rho          = Interp2D_lin(x, y, "dPde_rho" );
		interpolated.dTdrho_e          = Interp2D_lin(x, y, "dTdrho_e" );
		interpolated.dTde_rho          = Interp2D_lin(x, y, "dTde_rho" );
		interpolated.Cp                = Interp2D_lin(x, y, "Cp" );
		interpolated.Mu                = Interp2D_lin(x, y, "Mu" );
		interpolated.dmudrho_T         = Interp2D_lin(x, y, "dmudrho_T" );
		interpolated.dmudT_rho         = Interp2D_lin(x, y, "dmudT_rho" );
		interpolated.Kt                = Interp2D_lin(x, y, "Kt" );
		interpolated.dktdrho_T         = Interp2D_lin(x, y, "dktdrho_T" );
		interpolated.dktdT_rho         = Interp2D_lin(x, y, "dktdT_rho" );
		//Intermediate variables only needed for StandAlone version
		su2double Density = interpolated.Density;
		su2double Pressure = interpolated.Pressure;
		cout<<"Interpolated fit:"<<endl;
		interpolated.CTLprint ();
		if ((Density>Density_limits[1]) or (Density<Density_limits[0]))
		{
			throw runtime_error("PT Interpolated Density out of bounds");
		}
		if ((Pressure>Pressure_limits[1]) or (Pressure<Pressure_limits[0]))
		{
			throw runtime_error("PT Interpolated Pressure out of bounds");
		}
	}
	catch (exception& e)
	{
		cerr<<"\n"<< e.what() << "\n";
	}
}


void CLookUpTable::SetTDState_Prho (su2double P, su2double rho ) {
	try
	{
		//Check if inputs are in total range (necessary but not sufficient condition)
		if ((P>Pressure_limits[1]) or (P<Pressure_limits[0]))
		{
			throw runtime_error("PRHO Input Pressure out of bounds");
		}
		if ((rho>Density_limits[1]) or (rho<Density_limits[0]))
		{
			throw runtime_error("PRHO Input Density out of bounds");
		}
		cout<<endl<<"rho desired : "<<rho<<endl;
		cout<<"P desired   : "<<P<<endl;

		unsigned int LowerI;
		unsigned int UpperI;
		unsigned int LowerJ;
		unsigned int UpperJ;

		//Determine the I index: RHO is equispaced
		LowerI = floor((rho-Density_limits[0])/(Density_limits[1]-Density_limits[0])*(rho_dim-1));
		UpperI = LowerI + 1;

		//Determine the J index: P is equispaced (no restart)
		LowerJ = floor((P-Pressure_limits[0])/(Pressure_limits[1]-Pressure_limits[0])*(p_dim-1));
		UpperJ = LowerJ + 1;
		cout<<"i "<<LowerI<<"  "<<UpperI<<endl;
		cout<<"j "<<LowerJ<<"  "<<UpperJ<<endl;

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


		su2double x, y;
		x = rho - ThermoTables[iIndex][jIndex].Density;
		y = P - ThermoTables[iIndex][jIndex].Pressure;
		//Determine interpolation coefficients
		Interp2D_ArbitrarySkewCoeff(x,y,"PRHO");
		cout<<"Interpolation matrix inverse \n";
		for (int j=0; j<3; j++)
		{
			cout<<setw(15)<<coeff[j][0]<<"   "<<coeff[j][1]<<"   "<<coeff[j][2]<<endl;
		}
		CThermoList interpolated;
		interpolated.Pressure           = P;
		interpolated.Density          = rho ;
		interpolated.StaticEnergy      = Interp2D_lin(x, y, "StaticEnergy" );
		interpolated.Entropy           = Interp2D_lin(x, y, "Entropy" );
		interpolated.SoundSpeed2       = Interp2D_lin(x, y, "SoundSpeed2" );
		interpolated.Temperature       = Interp2D_lin(x, y, "Temperature" );
		interpolated.dPdrho_e          = Interp2D_lin(x, y, "dPdrho_e" );
		interpolated.dPde_rho          = Interp2D_lin(x, y, "dPde_rho" );
		interpolated.dTdrho_e          = Interp2D_lin(x, y, "dTdrho_e" );
		interpolated.dTde_rho          = Interp2D_lin(x, y, "dTde_rho" );
		interpolated.Cp                = Interp2D_lin(x, y, "Cp" );
		interpolated.Mu                = Interp2D_lin(x, y, "Mu" );
		interpolated.dmudrho_T         = Interp2D_lin(x, y, "dmudrho_T" );
		interpolated.dmudT_rho         = Interp2D_lin(x, y, "dmudT_rho" );
		interpolated.Kt                = Interp2D_lin(x, y, "Kt" );
		interpolated.dktdrho_T         = Interp2D_lin(x, y, "dktdrho_T" );
		interpolated.dktdT_rho         = Interp2D_lin(x, y, "dktdT_rho" );
		//Intermediate variables only needed for StandAlone version
		su2double Density = interpolated.Density;
		su2double Pressure = interpolated.Pressure;
		cout<<"Interpolated fit:"<<endl;
		interpolated.CTLprint ();
		if ((Density>Density_limits[1]) or (Density<Density_limits[0]))
		{
			throw runtime_error("PRHO Interpolated Density out of bounds");
		}
		if ((Pressure>Pressure_limits[1]) or (Pressure<Pressure_limits[0]))
		{
			throw runtime_error("PRHO Interpolated Pressure out of bounds");
		}
	}
	catch (exception& e)
	{
		cerr<<"\n"<< e.what() << "\n";
	}

}

void CLookUpTable::SetEnergy_Prho (su2double P, su2double rho ) {


}

void CLookUpTable::SetTDState_hs (su2double h, su2double s ) {

}

void CLookUpTable::SetTDState_Ps (su2double P, su2double s )
{
	try
	{
		//Check if inputs are in total range (necessary but not sufficient condition)
		if ((P>Pressure_limits[1]) or (P<Pressure_limits[0]))
		{
			throw runtime_error("PS Input Pressure out of bounds");
		}
		if ((s>Entropy_limits[1]) or (s<Entropy_limits[0]))
		{
			throw runtime_error("PS Input Entropy  out of bounds");
		}
		cout<<endl<<"P desired : "<<P<<endl;
		cout<<"s desired   : "<<s<<endl;

		su2double RunVal;
		unsigned int CFL = 2;
		unsigned int LowerI;
		unsigned int UpperI;
		unsigned int LowerJ;
		unsigned int UpperJ;
		//Restart search from previously used index if it exists, else go to middle
		if (iIndex<0) LowerI = 0;
		if (iIndex>=0) LowerI = iIndex;
		if  (iIndex<ceil(rho_dim/2))
		{
			UpperI = ceil(rho_dim/2); //probably can be made more efficient (smaller square)
		}
		else UpperI = rho_dim;

		//Determine the I index: RHO is equispaced (no restart)
		LowerJ = floor((P-Pressure_limits[0])/(Pressure_limits[1]-Pressure_limits[0])*(p_dim-1));
		UpperJ = LowerJ + 1;

		//Determine the J index (for s), (s invariant with j)
		int grad;
		while(UpperI-LowerI>1)
		{
			//Load the value at the upper bound
			grad = ThermoTables[UpperI][UpperJ].Entropy - ThermoTables[LowerI][UpperJ].Entropy;
			RunVal = ThermoTables[UpperI][UpperJ].Entropy ;
			if (grad*RunVal>grad*s)
			{
				UpperI = LowerI + ceil((UpperI-LowerI)/CFL);
			}
			else if (grad*RunVal<grad*s)
			{
				int dif;
				dif   = UpperI + ceil((UpperI-LowerI)/CFL);
				LowerI = UpperI;
				UpperI = dif;
			}
			cout<<LowerI<<"  "<<UpperI<<endl;
		}

		iIndex = LowerI;
		jIndex = LowerJ;
		cout<<"i "<<LowerI<<"  "<<UpperI<<endl;
		cout<<"j "<<LowerJ<<"  "<<UpperJ<<endl;


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


		su2double x, y;
		y = P - ThermoTables[iIndex][jIndex].Pressure;
		x = s - ThermoTables[iIndex][jIndex].Entropy ;
		//Determine interpolation coefficients
		Interp2D_ArbitrarySkewCoeff(x,y,"PS");
		cout<<"Interpolation matrix inverse \n";
		for (int j=0; j<3; j++)
		{
			cout<<setw(15)<<coeff[j][0]<<"   "<<coeff[j][1]<<"   "<<coeff[j][2]<<endl;
		}
		CThermoList interpolated;
		interpolated.Entropy           = s;
		interpolated.Pressure          = P ;
		interpolated.StaticEnergy      = Interp2D_lin(x, y, "StaticEnergy" );
		interpolated.Density           = Interp2D_lin(x, y, "Density" );
		interpolated.SoundSpeed2       = Interp2D_lin(x, y, "SoundSpeed2" );
		interpolated.Temperature       = Interp2D_lin(x, y, "Temperature" );
		interpolated.dPdrho_e          = Interp2D_lin(x, y, "dPdrho_e" );
		interpolated.dPde_rho          = Interp2D_lin(x, y, "dPde_rho" );
		interpolated.dTdrho_e          = Interp2D_lin(x, y, "dTdrho_e" );
		interpolated.dTde_rho          = Interp2D_lin(x, y, "dTde_rho" );
		interpolated.Cp                = Interp2D_lin(x, y, "Cp" );
		interpolated.Mu                = Interp2D_lin(x, y, "Mu" );
		interpolated.dmudrho_T         = Interp2D_lin(x, y, "dmudrho_T" );
		interpolated.dmudT_rho         = Interp2D_lin(x, y, "dmudT_rho" );
		interpolated.Kt                = Interp2D_lin(x, y, "Kt" );
		interpolated.dktdrho_T         = Interp2D_lin(x, y, "dktdrho_T" );
		interpolated.dktdT_rho         = Interp2D_lin(x, y, "dktdT_rho" );
		//Intermediate variables only needed for StandAlone version
		su2double Density = interpolated.Density;
		su2double Pressure = interpolated.Pressure;
		cout<<"Interpolated fit:"<<endl;
		interpolated.CTLprint ();
		if ((Density>Density_limits[1]) or (Density<Density_limits[0]))
		{
			throw runtime_error("PS Interpolated Density out of bounds");
		}
		if ((Pressure>Pressure_limits[1]) or (Pressure<Pressure_limits[0]))
		{
			throw runtime_error("PS Interpolated Pressure out of bounds");
		}
	}
	catch (exception& e)
	{
		cerr<<"\n"<< e.what() << "\n";
	}
}

void CLookUpTable::SetTDState_rhoT (su2double rho, su2double T ) {
	try
	{
		//Check if inputs are in total range (necessary but not sufficient condition)
		if ((rho>Density_limits[1]) or (rho<Density_limits[0]))
		{
			throw runtime_error("RHOT Input Density out of bounds");
		}
		if ((T>Temperature_limits[1]) or (T<Temperature_limits[0]))
		{
			throw runtime_error("RHOT Input Temperature out of bounds");
		}
		cout<<endl<<"rho desired : "<<rho<<endl;
		cout<<"T desired   : "<<T<<endl;

		su2double RunVal;
		unsigned int CFL= 2;
		unsigned int LowerI;
		unsigned int UpperI;
		unsigned int LowerJ;
		unsigned int UpperJ;
		//Restart search from previously used index if it exists, else go to middle
		if (jIndex<0) LowerJ = 0;
		if (jIndex>=0) LowerJ = jIndex;
		if (jIndex<ceil(p_dim/2))
		{
			UpperJ = ceil(p_dim/2);
		}
		else UpperJ=p_dim;

		//Determine the I index: RHO is equispaced (no restart)
		LowerI = floor((rho-Density_limits[0])/(Density_limits[1]-Density_limits[0])*(rho_dim-1));
		UpperI = LowerI + 1;

		//Determine the I index (for T)
		//Temperature is not necessarily monotonic (so use a variable to take care of overall sign)
		while(UpperI-LowerI>1)
		{
			int grad = ThermoTables[UpperI][UpperJ].Temperature - ThermoTables[LowerI][UpperJ].Temperature;
			//Load the value at the upper bound
			RunVal = ThermoTables[UpperI][UpperJ].Temperature;
			if (grad*RunVal>grad*T)
			{
				UpperI = LowerI + ceil((UpperI-LowerI)/CFL);
			}
			else if (grad*RunVal<grad*T)
			{
				int dif;
				dif   = UpperI + ceil((UpperI-LowerI)/CFL);
				LowerI = UpperI;
				UpperI = dif;
			}
			cout<<LowerI<<"  "<<UpperI<<endl;
		}

		iIndex = LowerI;
		jIndex = LowerJ;
		cout<<"i "<<LowerI<<"  "<<UpperI<<endl;
		cout<<"j "<<LowerJ<<"  "<<UpperJ<<endl;


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


		su2double x, y;
		x = rho - ThermoTables[iIndex][jIndex].Density;
		y = T - ThermoTables[iIndex][jIndex].Temperature;
		//Determine interpolation coefficients
		Interp2D_ArbitrarySkewCoeff(x,y,"RHOT");
		cout<<"Interpolation matrix inverse \n";
		for (int j=0; j<3; j++)
		{
			cout<<setw(15)<<coeff[j][0]<<"   "<<coeff[j][1]<<"   "<<coeff[j][2]<<endl;
		}
		CThermoList interpolated;
		interpolated.Temperature       = T;
		interpolated.Density           = rho ;
		interpolated.StaticEnergy      = Interp2D_lin(x, y, "StaticEnergy" );
		interpolated.Entropy           = Interp2D_lin(x, y, "Entropy" );
		interpolated.Pressure          = Interp2D_lin(x, y, "Pressure" );
		interpolated.SoundSpeed2       = Interp2D_lin(x, y, "SoundSpeed2" );
		interpolated.dPdrho_e          = Interp2D_lin(x, y, "dPdrho_e" );
		interpolated.dPde_rho          = Interp2D_lin(x, y, "dPde_rho" );
		interpolated.dTdrho_e          = Interp2D_lin(x, y, "dTdrho_e" );
		interpolated.dTde_rho          = Interp2D_lin(x, y, "dTde_rho" );
		interpolated.Cp                = Interp2D_lin(x, y, "Cp" );
		interpolated.Mu                = Interp2D_lin(x, y, "Mu" );
		interpolated.dmudrho_T         = Interp2D_lin(x, y, "dmudrho_T" );
		interpolated.dmudT_rho         = Interp2D_lin(x, y, "dmudT_rho" );
		interpolated.Kt                = Interp2D_lin(x, y, "Kt" );
		interpolated.dktdrho_T         = Interp2D_lin(x, y, "dktdrho_T" );
		interpolated.dktdT_rho         = Interp2D_lin(x, y, "dktdT_rho" );
		//Intermediate variables only needed for StandAlone version
		su2double Density = interpolated.Density;
		su2double Pressure = interpolated.Pressure;
		cout<<"Interpolated fit:"<<endl;
		interpolated.CTLprint ();
		if ((Density>Density_limits[1]) or (Density<Density_limits[0]))
		{
			throw runtime_error("RHOT Interpolated Density out of bounds");
		}
		if ((Pressure>Pressure_limits[1]) or (Pressure<Pressure_limits[0]))
		{
			throw runtime_error("RHOT Interpolated Pressure out of bounds");
		}
	}
	catch (exception& e)
	{
		cerr<<"\n"<< e.what() << "\n";
	}


}


void CLookUpTable::Interp2D_ArbitrarySkewCoeff(su2double x, su2double y, std::string grid_var)
{
	//Distances in along x and y axis are taken relative to i,j point (x00, y00).
	//This reduces the interpolation to a 3by3 system rather than 4by4
	//x and y are not strictrly necessary for the calculation of the coefficients.
	//However, they do allow for checking whether the point of interest is contained in
	//the quad under consideration.
	su2double x00, y00, dx10, dx01, dx11, dy10, dy01, dy11;
	//Interpolation LHM
	su2double A[3][3];
	//Helper variable for Gaussian elimination
	su2double c;
	//Load in the coordinates of the qudrilateral (values relative to i,j)
	if(grid_var=="RHOE")
	{
		x00  = ThermoTables[iIndex  ][jIndex  ].Density     ;
		y00  = ThermoTables[iIndex  ][jIndex  ].StaticEnergy;
		dx01 = ThermoTables[iIndex  ][jIndex+1].Density      -x00;
		dy01 = ThermoTables[iIndex  ][jIndex+1].StaticEnergy -y00;
		dx10 = ThermoTables[iIndex+1][jIndex  ].Density      -x00;
		dy10 = ThermoTables[iIndex+1][jIndex  ].StaticEnergy -y00;
		dx11 = ThermoTables[iIndex+1][jIndex+1].Density      -x00;
		dy11 = ThermoTables[iIndex+1][jIndex+1].StaticEnergy -y00;
	}
	else if(grid_var=="PT")
	{
		y00  = ThermoTables[iIndex  ][jIndex  ].Pressure   ;
		x00  = ThermoTables[iIndex  ][jIndex  ].Temperature;
		dy01 = ThermoTables[iIndex  ][jIndex+1].Pressure    -y00;
		dx01 = ThermoTables[iIndex  ][jIndex+1].Temperature -x00;
		dy10 = ThermoTables[iIndex+1][jIndex  ].Pressure    -y00;
		dx10 = ThermoTables[iIndex+1][jIndex  ].Temperature -x00;
		dy11 = ThermoTables[iIndex+1][jIndex+1].Pressure    -y00;
		dx11 = ThermoTables[iIndex+1][jIndex+1].Temperature -x00;
	}
	else if(grid_var=="PRHO")
	{
		y00  = ThermoTables[iIndex  ][jIndex  ].Pressure;
		x00  = ThermoTables[iIndex  ][jIndex  ].Density ;
		dy01 = ThermoTables[iIndex  ][jIndex+1].Pressure -y00;
		dx01 = ThermoTables[iIndex  ][jIndex+1].Density  -x00;
		dy10 = ThermoTables[iIndex+1][jIndex  ].Pressure -y00;
		dx10 = ThermoTables[iIndex+1][jIndex  ].Density  -x00;
		dy11 = ThermoTables[iIndex+1][jIndex+1].Pressure -y00;
		dx11 = ThermoTables[iIndex+1][jIndex+1].Density  -x00;
	}
	else if(grid_var=="RHOT")
	{
		x00  = ThermoTables[iIndex  ][jIndex  ].Density    ;
		y00  = ThermoTables[iIndex  ][jIndex  ].Temperature;
		dx01 = ThermoTables[iIndex  ][jIndex+1].Density     -x00;
		dy01 = ThermoTables[iIndex  ][jIndex+1].Temperature -y00;
		dx10 = ThermoTables[iIndex+1][jIndex  ].Density     -x00;
		dy10 = ThermoTables[iIndex+1][jIndex  ].Temperature -y00;
		dx11 = ThermoTables[iIndex+1][jIndex+1].Density     -x00;
		dy11 = ThermoTables[iIndex+1][jIndex+1].Temperature -y00;
	}
	else if(grid_var=="PS")
	{
		y00  = ThermoTables[iIndex  ][jIndex  ].Pressure;
		x00  = ThermoTables[iIndex  ][jIndex  ].Entropy ;
		dy01 = ThermoTables[iIndex  ][jIndex+1].Pressure -y00;
		dx01 = ThermoTables[iIndex  ][jIndex+1].Entropy  -x00;
		dy10 = ThermoTables[iIndex+1][jIndex  ].Pressure -y00;
		dx10 = ThermoTables[iIndex+1][jIndex  ].Entropy  -x00;
		dy11 = ThermoTables[iIndex+1][jIndex+1].Pressure -y00;
		dx11 = ThermoTables[iIndex+1][jIndex+1].Entropy  -x00;
	}
	else if(grid_var=="HS")
	{
		x00  = ThermoTables[iIndex  ][jIndex  ].StaticEnergy;
		y00  = ThermoTables[iIndex  ][jIndex  ].Entropy     ;
		dx01 = ThermoTables[iIndex  ][jIndex+1].StaticEnergy -x00;
		dy01 = ThermoTables[iIndex  ][jIndex+1].Entropy      -y00;
		dx10 = ThermoTables[iIndex+1][jIndex  ].StaticEnergy -x00;
		dy10 = ThermoTables[iIndex+1][jIndex  ].Entropy      -y00;
		dx11 = ThermoTables[iIndex+1][jIndex+1].StaticEnergy -x00;
		dy11 = ThermoTables[iIndex+1][jIndex+1].Entropy      -y00;
	}
	//Check if x, y is indeed in the quad
	//Some extra logic is needed as the both monotonically increasing and monotonically decreasing functions
	//have to be anticipated
	bool BOTTOM, TOP, LEFT, RIGHT;
	BOTTOM = (y*dx10)<(x*dy10);
	TOP = ((y-dy01)*(dx11-dx01))>((dy11-dy01)*(x-dx01));
	RIGHT = ((x-dx10)*(dy11-dy10))>((dx11-dx10)*(y-dy10));
	LEFT = (x*dy01)<(dx01*y);
	//Check BOTTOM quad boundary
	if(BOTTOM and !TOP)
	{
		throw runtime_error(grid_var+" interpolation point lies below bottom boundary of selected quad");
	}
	//Check RIGHT quad boundary
	if(RIGHT and !LEFT)
	{
		throw runtime_error(grid_var+" interpolation point lies to the right of the boundary of selected quad");
	}
	//Check TOP quad boundary
	if(TOP and !BOTTOM)
	{
		throw runtime_error(grid_var+" interpolation point lies above the boundary of selected quad");
	}
	//Check LEFT quad boundary

	if(LEFT and !RIGHT)
	{
		throw runtime_error(grid_var+" interpolation point lies to the left of the boundary of selected quad");
	}
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
	cout<<"Interpolation LHM matrix \n"<<"[";

	for (int j=0; j<3; j++)
	{
		cout<<setw(15)<<"["<<A[j][0]<<" ,  "<<A[j][1]<<"  , "<<A[j][2]<<"]"<<endl;
	}
	cout<<"]\n";

	//Store the inverse of the LHM matrix as coeff
	coeff[0][0] = 1;
	coeff[0][1] = 0;
	coeff[0][2] = 0;
	coeff[1][0] = 0;
	coeff[1][1] = 1;
	coeff[1][2] = 0;//solved interpolation bug
	coeff[2][0] = 0;
	coeff[2][1] = 0;
	coeff[2][2] = 1;

	//Compute inverse of LHM using Gaussian elimination
	//Reduced Echelon form of the LHM
	if (A[0][0] != 0)
	{
		c = A[1][0]/A[0][0];
		coeff[1][0] = coeff[1][0] -coeff[0][0]*c;
		for (int i=0; i<3; i++)
		{
			A[1][i] = A[1][i] -A[0][i]*c;
		}
		c =A[2][0]/A[0][0];
		coeff[2][0] = coeff[2][0] -coeff[0][0]*c;

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
			coeff[2][i] = coeff[2][i] -coeff[1][i]*c;
		}
		for (int i=0; i<3; i++)
			A[2][i] = A[2][i] -A[1][i]*c;
	}
	//Reduced reduced Echelon form of LHM
	if (A[2][2] != 0)
	{

		for (int i=0; i<3; i++)
		{
			coeff[1][i] = coeff[1][i] -coeff[2][i]*A[1][2]/A[2][2];
		}
		A[1][2] = 0;
		for (int i=0; i<3; i++)
		{
			coeff[0][i] = coeff[0][i] -coeff[2][i]*A[0][2]/A[2][2];
		}
		A[0][2] = 0;
	}
	if (A[1][1] != 0)
	{
		for (int i=0; i<3; i++)
			coeff[0][i] = coeff[0][i] -coeff[1][i]*A[0][1]/A[1][1];
	}
	A[0][1] = 0;

	//Normalize the RR Echelon form
	if (A[0][0] != 0)
	{
		for (int i=0; i<3; i++)
		{
			coeff[0][i] = coeff[0][i]/A[0][0];
		}
		if (A[1][1] != 0)
		{
			for (int i=0; i<3; i++)
			{
				coeff[1][i] = coeff[1][i]/A[1][1];
			}
		}
		if (A[2][2] != 0)
		{
			for (int i=0; i<3; i++)
			{
				coeff[2][i] = coeff[2][i]/A[2][2];
			}

		}
	}
	return;
}

su2double CLookUpTable::Interp2D_lin(su2double x, su2double y, string interpolant_var)
{
	//F is the RHS part of the interpolation equation
	su2double F[3];
	//The solution vector for the interpolation equation
	su2double C[3];
	//The function values
	su2double f00, f10, f01, f11;
	//For each case the values are filled differently
	if(interpolant_var=="StaticEnergy")
	{
		f00 = ThermoTables[iIndex  ][jIndex  ].StaticEnergy;
		f10 = ThermoTables[iIndex+1][jIndex  ].StaticEnergy;
		f01 = ThermoTables[iIndex  ][jIndex+1].StaticEnergy;
		f11 = ThermoTables[iIndex+1][jIndex+1].StaticEnergy;
	}
	else if(interpolant_var=="Entropy")
	{
		f00 = ThermoTables[iIndex  ][jIndex  ].Entropy;
		f10 = ThermoTables[iIndex+1][jIndex  ].Entropy;
		f01 = ThermoTables[iIndex  ][jIndex+1].Entropy;
		f11 = ThermoTables[iIndex+1][jIndex+1].Entropy;
	}
	else if(interpolant_var=="Density")
	{
		f00 = ThermoTables[iIndex  ][jIndex  ].Density;
		f10 = ThermoTables[iIndex+1][jIndex  ].Density;
		f01 = ThermoTables[iIndex  ][jIndex+1].Density;
		f11 = ThermoTables[iIndex+1][jIndex+1].Density;
	}
	else if(interpolant_var=="Pressure")
	{
		f00 = ThermoTables[iIndex  ][jIndex  ].Pressure;
		f10 = ThermoTables[iIndex+1][jIndex  ].Pressure;
		f01 = ThermoTables[iIndex  ][jIndex+1].Pressure;
		f11 = ThermoTables[iIndex+1][jIndex+1].Pressure;
	}
	else if(interpolant_var=="SoundSpeed2")
	{
		f00 = ThermoTables[iIndex  ][jIndex  ].SoundSpeed2;
		f10 = ThermoTables[iIndex+1][jIndex  ].SoundSpeed2;
		f01 = ThermoTables[iIndex  ][jIndex+1].SoundSpeed2;
		f11 = ThermoTables[iIndex+1][jIndex+1].SoundSpeed2;
	}
	else if(interpolant_var=="Temperature")
	{
		f00 = ThermoTables[iIndex  ][jIndex  ].Temperature;
		f10 = ThermoTables[iIndex+1][jIndex  ].Temperature;
		f01 = ThermoTables[iIndex  ][jIndex+1].Temperature;
		f11 = ThermoTables[iIndex+1][jIndex+1].Temperature;
	}
	else if(interpolant_var=="dPdrho_e")
	{
		f00 = ThermoTables[iIndex  ][jIndex  ].dPdrho_e;
		f10 = ThermoTables[iIndex+1][jIndex  ].dPdrho_e;
		f01 = ThermoTables[iIndex  ][jIndex+1].dPdrho_e;
		f11 = ThermoTables[iIndex+1][jIndex+1].dPdrho_e;
	}
	else if(interpolant_var=="dPde_rho")
	{
		f00 = ThermoTables[iIndex  ][jIndex  ].dPde_rho;
		f10 = ThermoTables[iIndex+1][jIndex  ].dPde_rho;
		f01 = ThermoTables[iIndex  ][jIndex+1].dPde_rho;
		f11 = ThermoTables[iIndex+1][jIndex+1].dPde_rho;
	}
	else if(interpolant_var=="dTdrho_e")
	{
		f00 = ThermoTables[iIndex  ][jIndex  ].dTdrho_e;
		f10 = ThermoTables[iIndex+1][jIndex  ].dTdrho_e;
		f01 = ThermoTables[iIndex  ][jIndex+1].dTdrho_e;
		f11 = ThermoTables[iIndex+1][jIndex+1].dTdrho_e;
	}
	else if(interpolant_var=="dTde_rho")
	{
		f00 = ThermoTables[iIndex  ][jIndex  ].dTde_rho;
		f10 = ThermoTables[iIndex+1][jIndex  ].dTde_rho;
		f01 = ThermoTables[iIndex  ][jIndex+1].dTde_rho;
		f11 = ThermoTables[iIndex+1][jIndex+1].dTde_rho;
	}
	else if(interpolant_var=="Cp")
	{
		f00 = ThermoTables[iIndex  ][jIndex  ].Cp;
		f10 = ThermoTables[iIndex+1][jIndex  ].Cp;
		f01 = ThermoTables[iIndex  ][jIndex+1].Cp;
		f11 = ThermoTables[iIndex+1][jIndex+1].Cp;
	}
	else if(interpolant_var=="Mu")
	{
		f00 = ThermoTables[iIndex  ][jIndex  ].Mu;
		f10 = ThermoTables[iIndex+1][jIndex  ].Mu;
		f01 = ThermoTables[iIndex  ][jIndex+1].Mu;
		f11 = ThermoTables[iIndex+1][jIndex+1].Mu;
	}
	else if(interpolant_var=="dmudrho_T")
	{
		f00 = ThermoTables[iIndex  ][jIndex  ].dmudrho_T;
		f10 = ThermoTables[iIndex+1][jIndex  ].dmudrho_T;
		f01 = ThermoTables[iIndex  ][jIndex+1].dmudrho_T;
		f11 = ThermoTables[iIndex+1][jIndex+1].dmudrho_T;
	}
	else if(interpolant_var=="dmudT_rho")
	{
		f00 = ThermoTables[iIndex  ][jIndex  ].dmudT_rho;
		f10 = ThermoTables[iIndex+1][jIndex  ].dmudT_rho;
		f01 = ThermoTables[iIndex  ][jIndex+1].dmudT_rho;
		f11 = ThermoTables[iIndex+1][jIndex+1].dmudT_rho;
	}
	else if(interpolant_var=="Kt")
	{
		f00 = ThermoTables[iIndex  ][jIndex  ].Kt;
		f10 = ThermoTables[iIndex+1][jIndex  ].Kt;
		f01 = ThermoTables[iIndex  ][jIndex+1].Kt;
		f11 = ThermoTables[iIndex+1][jIndex+1].Kt;
	}
	else if(interpolant_var=="dktdrho_T")
	{
		f00 = ThermoTables[iIndex  ][jIndex  ].dktdrho_T;
		f10 = ThermoTables[iIndex+1][jIndex  ].dktdrho_T;
		f01 = ThermoTables[iIndex  ][jIndex+1].dktdrho_T;
		f11 = ThermoTables[iIndex+1][jIndex+1].dktdrho_T;
	}
	else if(interpolant_var=="dktdT_rho")
	{
		f00 = ThermoTables[iIndex  ][jIndex  ].dktdT_rho;
		f10 = ThermoTables[iIndex+1][jIndex  ].dktdT_rho;
		f01 = ThermoTables[iIndex  ][jIndex+1].dktdT_rho;
		f11 = ThermoTables[iIndex+1][jIndex+1].dktdT_rho;
	}

	//Using offset relative to i,j point yields a 3by3 system rather than 4by4
	F[0] = f10 - f00;
	C[0] = 0;
	F[1] = f01 - f00;
	C[1] = 0;
	F[2] = f11 - f00;
	C[2] = 0;
	for (int i = 0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			C[i] = C[i] + F[i]*coeff[i][j];
		}
	}

	return f00 + C[0]*x + C[1]*y + C[2]*x*y;
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


void CLookUpTable::TableLoadCFX(char* filename){
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
					var_steps = 10; //solved pressure reading bug
					Pressure_limits[0] = 10E15; //lower limit
					Pressure_limits[1] = 0; //upper limit
					//Each line contains at most 10 pressure values
					for (int k =0; k<ceil(float(set_y)/10.0);k++)
					{

						getline(table,line);
						cout<<line<<endl;
						istringstream inP(line);
						//Check if line contains less than 10 values
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










