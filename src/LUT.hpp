#define su2double double
#include <string>

struct KD_node
{
	int depth, dim, *i_values;
	su2double * x_values, * y_values;
	KD_node* upper;
	KD_node* lower;
};

class CThermoList {
public:
	su2double StaticEnergy,			/*!< \brief Internal Energy. */
	Entropy,  				/*!< \brief Entropy. */
	Enthalpy,
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
	void CTLprint(void);

};

/*!
 * \class CLookUpTable
 * \brief Child class for defining ideal gas model.
 * \author: A. Rubino, S.Vitale.
 * \version 4.1.2 "Cardinal"
 */
class CLookUpTable {

protected:
	CThermoList **ThermoTables;
	su2double coeff[3][3]; //Interpolation coefficients
	/*!< \brief Fluid derivative DktDT_rho. */
	long iIndex, jIndex;
	int p_dim, rho_dim; /*!< \brief The pressure and density dimensions of the table */
	KD_node *HS_tree; //KD tree for HS thermoPair
	CThermoList interpolated;

	su2double StaticEnergy_limits[2];
	su2double Entropy_limits[2];
	su2double Enthalpy_limits[2];
	su2double Density_limits[2];
	su2double Pressure_limits[2];
	su2double SoundSpeed2_limits[2];
	su2double Temperature_limits[2];
	su2double dPdrho_e_limits[2];
	su2double dPde_rho_limits[2];
	su2double dTdrho_e_limits[2];
	su2double dTde_rho_limits[2];
	su2double Cp_limits[2];
	su2double Mu_limits[2];
	su2double dmudrho_T_limits[2];
	su2double dmudT_rho_limits[2];
	su2double Kt_limits[2];
	su2double dktdrho_T_limits[2];
	su2double dktdT_rho_limits[2];
	//Nearest neighbour's i and j indexes
	int NN_i[4];
	int NN_j[4];

public:

	/*!
	 * \brief default Constructor of the class.
	 */
	CLookUpTable(void);

	CLookUpTable(char* Filename);



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
	void reset_Restart();
	struct KD_node* KD_Tree(su2double* x_values, su2double* y_values, int* i_values, int dim, int depth);
	su2double Dist_KD_Tree (su2double x, su2double y, KD_node *branch);
	void free_KD_tree(KD_node* root);
	void NN_KD_Tree (su2double thermo1, su2double thermo2, KD_node *root, su2double best_dist);
	void NN4_KD_Tree (su2double thermo1, su2double thermo2, KD_node *root, su2double* best_dist);
	void SearchZigZag (su2double thermo1, su2double thermo2,  unsigned long thermoPair );
	void SearchThermoPair (su2double thermo1, su2double thermo2,  unsigned long thermoPair );


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
	void Interp2D_SingleSkewCoeff(std::string grid_var);
	su2double Interp2D_Inv_Dist(std::string interpolant_var, su2double* dist);
	void Interp2D_ArbitrarySkewCoeff(su2double x, su2double y, std::string grid_var);
	su2double Interp2D_lin(su2double x, su2double y, std::string interpolant_var);
	void TableLoadCFX(char* filename);
	void LUTprint(void);
	void TableDump(char* filename);
	void RecordState(char* file);


};




