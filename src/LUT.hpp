#define su2double double

/*!
 * \class CLookUpTable
 * \brief Child class for defining ideal gas model.
 * \author: A. Rubino, S.Vitale.
 * \version 4.1.2 "Cardinal"
 */
class CLookUpTable {

protected:
	CThermoList **ThermoTables;
	su2double *coeff; /*!< \brief Fluid derivative DktDT_rho. */
	unsigned long iIndex, jIndex;


public:

	/*!
	 * \brief default Constructor of the class.
	 */
	CLookUpTable(void);

	CLookUpTable(CThermoList **TL);


	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CLookUpTable(void);


	/*!
	 * \brief Search Thermo Pair
	 * \param[in] thermo1 - first thermodynamic variable.
	 * \param[in] thermo2 - second thermodynamic variable
	 * \param[in] input thermodynamic pair.
	 */
	void SearchThermoPair (su2double thermo1, su2double thermo2,  unsigned short thermoPair );


	/*!
	 * \brief Set the Dimensionless State using Density and Internal Energy
	 * \param[in] rho - first thermodynamic variable.
	 * \param[in] e - second thermodynamic variable.
	 */
	void SetTDState_rhoe (su2double rho, su2double e );

	/*!
	 * \brief Set the Dimensionless State using Pressure and Temperature
	 * \param[in] P - first thermodynamic variable.
	 * \param[in] T - second thermodynamic variable.
	 */
	void SetTDState_PT (su2double P, su2double T );

	/*!
	 * \brief Set the Dimensionless State using Pressure and Density
	 * \param[in] P - first thermodynamic variable.
	 * \param[in] rho - second thermodynamic variable.
	 */
	void SetTDState_Prho (su2double P, su2double rho );

	/*!
	 * \brief Set the Dimensionless Internal Energy using Pressure and Density
	 * \param[in] P - first thermodynamic variable.
	 * \param[in] rho - second thermodynamic variable.
	 */
	void SetEnergy_Prho (su2double P, su2double rho );

	/*!
	 * \brief Set the Dimensionless state using Enthalpy and Entropy
	 * \param[in] h - first thermodynamic variable (h).
	 * \param[in] s - second thermodynamic variable (s).
	 *
	 */
	void SetTDState_hs (su2double h, su2double s );


	/*!
	 * \brief Set the Dimensionless state using Density and Temperature
	 * \param[in] rho - first thermodynamic variable (rho).
	 * \param[in] T - second thermodynamic variable (T).
	 *
	 */
	void SetTDState_rhoT (su2double rho, su2double T );

	/*!
	 * \brief Set the Dimensionless State using Pressure and Entropy
	 * \param[in] P - first thermodynamic variable (P).
	 * \param[in] s - second thermodynamic variable (s).
	 */

	void SetTDState_Ps (su2double P, su2double s );

	void Interp2D_lin (su2double aa, su2double ab, su2double ba, su2double bb);
	void Interp2D_lin(su2double aa, su2double ab, su2double ba, su2double bb, su2double* coeffs);

};

class CThermoList {
public:
	su2double   	 StaticEnergy,			/*!< \brief Internal Energy. */
	Entropy,  				/*!< \brief Entropy. */
	Density,  				/*!< \brief Density. */
	Pressure, 				/*!< \brief Pressure. */
	SoundSpeed2, 		/*!< \brief SpeedSound. */
	Temperature,			/*!< \brief Temperature. */
	dPdrho_e, 				/*!< \brief Fluid derivative DpDd_e. */
	dPde_rho, 				/*!< \brief Fluid derivative DpDe_d. */
	dTdrho_e, 				/*!< \brief Fluid derivative DTDd_e. */
	dTde_rho, 				/*!< \brief Fluid derivative DTDe_d. */
	Cp,              /*!< \brief Specific Heat Capacity at constant pressure. */
	Mu,					    /*!< \brief Laminar Viscosity. */
	dmudrho_T, 			/*!< \brief Fluid derivative DmuDrho_T */
	dmudT_rho,				/*!< \brief Fluid derivative DmuDT_rho. */
	Kt,					    /*!< \brief Thermal Conductivity. */
	dktdrho_T, 			/*!< \brief Fluid derivative DktDrho_T.  */
	dktdT_rho;				/*!< \brief Fluid derivative DktDT_rho. */


	/*!
	 * \brief Constructor of the class.
	 */
	CThermoList(void);

	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CThermoList(void);


};


CThermoList::CThermoList(){

	StaticEnergy = 0.0;
	Entropy      = 0.0;
	Density      = 0.0;
	Pressure     = 0.0;
	SoundSpeed2  = 0.0;
	Temperature  = 0.0;
	dPdrho_e     = 0.0;
	dPde_rho     = 0.0;
	dTdrho_e     = 0.0;
	dTde_rho     = 0.0;
	Cp           = 0.0;
	Mu 			 = 0.0;
	dmudrho_T    = 0.0;
	dmudT_rho    = 0.0;
	Kt           = 0.0;
	dktdrho_T    = 0.0;
	dktdT_rho    = 0.0;

}

CThermoList::~CThermoList(){

}


