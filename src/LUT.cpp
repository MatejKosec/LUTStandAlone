/*!
 * fluid_model_lut.cpp
 * \brief Source of the look-up table model.
 * \author S. Vitale, A. Rubino
 * \version 4.1.2 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#define su2double double
#include "LUT.hpp"
#include "stdlib.h"
#include "stdio.h"
#include "TableReader.hpp"

CLookUpTable::CLookUpTable() {

	ThermoTables = NULL;
	coeff = NULL;
	iIndex = 0;
	jIndex = 0;
}

CLookUpTable::CLookUpTable(CThermoList** TL ) {

	this->ThermoTables = TL;
	coeff = NULL;
	iIndex = 0;
	jIndex = 0;
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

	SearchThermoPair (rho, e, 'RHOE');
	CThermoList read_in;
	read_in = this->ThermoTables[iIndex][jIndex];

	su2double StaticEnergy = e;
	su2double	Density    = rho;
	su2double Entropy      = 0.0;
	su2double Pressure     = 0.0;
	su2double SoundSpeed2  = 0.0;
	su2double Temperature  = 0.0;
	su2double dPdrho_e     = 0.0;
	su2double dPde_rho     = 0.0;
	su2double dTdrho_e     = 0.0;
	su2double dTde_rho     = 0.0;
	su2double Cp           = 0.0;
	su2double Mu 		   = 0.0;
	su2double dmudrho_T    = 0.0;
	su2double dmudT_rho    = 0.0;
	su2double Kt           = 0.0;
	su2double dktdrho_T    = 0.0;
	su2double dktdT_rho    = 0.0;


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

void CLookUpTable::Interp2D_lin(su2double aa, su2double ab, su2double ba, su2double bb, su2double* coeffs) {



}










