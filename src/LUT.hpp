#define su2double double
#include <vector>
#include <string>

using namespace std;

/*!
 * \class CTrapezoidalMap
 * \brief An algorithm for finding the polygon
 * containing the query vector. Adapted from:
 * Computational Geometry: Algorithms and Applications,
 * 3rd edition, 2008, by M de Berg, et al.
 * \author: M.Kosec
 * \version 4.1.2 "Cardinal"
 */

class CTrapezoidalMap {
protected:
	int rank, UpperI, LowerI, middleI, LowerJ, UpperJ, middleJ;
	int UpperEdge, MiddleEdge, LowerEdge;
	//The unique values of x which exist in the data
	vector<su2double> Unique_X_Bands;
	vector<vector<int> > Unique_Edges;
	vector<vector<su2double> > X_Limits_of_Edges, Y_Limits_of_Edges;
	//The value that each edge which intersects the band takes within that
	//same band. Used to sort the edges
	vector<vector<pair<su2double, int> > > Y_Values_of_Edge_Within_Band_And_Index;
public:
	CTrapezoidalMap();
	CTrapezoidalMap(vector<su2double> const &x_samples,
			vector<su2double> const &y_samples,
			vector<vector<int> > const &unique_edges);
	void Find_Containing_Simplex(su2double x, su2double y);
	void Search_Bands_For(su2double x);
	void Search_Band_For_Edge(su2double x, su2double y);

	int getLowerI() const {
		return LowerI;
	}

	int getLowerJ() const {
		return LowerJ;
	}

	int getUpperI() const {
		return UpperI;
	}

	int getUpperJ() const {
		return UpperJ;
	}

	int getLowerEdge() const {
		return LowerEdge;
	}

	int getUpperEdge() const {
		return UpperEdge;
	}
	int getMiddleEdge() const {
		return MiddleEdge;
	}
};

/*!
 * \class CLookUpTable
 * \brief Class for defining a lookuptable fluid model
 * \author: A. Rubino, S.Vitale., M. Kosec
 * \version 4.1.2 "Cardinal"
 */
class CLookUpTable {

protected:
	int rank;
	int CurrentZone;
	vector<int> CurrentPoints;
	bool LUT_Debug_Mode;/*!< \brief If true, master node prints errors of points outside LUT*/
	su2double Pressure_Reference_Value;su2double Density_Reference_Value;su2double Temperature_Reference_Value;su2double Velocity_Reference_Value;su2double Energy_Reference_Value;

	//Put the trapezoidal maps into variables
	CTrapezoidalMap rhoe_map[2], Prho_map[2], hs_map[2], Ps_map[2], rhoT_map[2],
			PT_map[2];

	su2double StaticEnergy, /*!< \brief Internal Energy. */
	Entropy, /*!< \brief Entropy. */
	Enthalpy, /*!< \brief Enthalpy required as separate variable for use in HS tree. */
	Density, /*!< \brief Density. */
	Pressure, /*!< \brief Pressure. */
	SoundSpeed2, /*!< \brief The speed of sound squared. */
	Temperature, /*!< \brief Temperature. */
	dPdrho_e, /*!< \brief Fluid derivative DpDd_e. */
	dPde_rho, /*!< \brief Fluid derivative DpDe_d. */
	dTdrho_e, /*!< \brief Fluid derivative DTDd_e. */
	dTde_rho, /*!< \brief Fluid derivative DTDe_d. */
	Cp, /*!< \brief Specific Heat Capacity at constant pressure. */
	Mu, /*!< \brief Laminar Viscosity. */
	dmudrho_T, /*!< \brief Fluid derivative DmuDrho_T */
	dmudT_rho, /*!< \brief Fluid derivative DmuDT_rho. */
	Kt, /*!< \brief Thermal Conductivity. */
	dktdrho_T, /*!< \brief Fluid derivative DktDrho_T.  */
	dktdT_rho; /*!< \brief Fluid derivative DktDT_rho. */

	vector< su2double > ThermoTables_StaticEnergy[2], /*!< \brief Internal Energy look up table values. */
	ThermoTables_Entropy[2], /*!< \brief Entropy look up table values. */
	ThermoTables_Enthalpy[2], /*!< \brief Enthalpy required as separate variable for use in HS tree look up table values. */
	ThermoTables_Density[2], /*!< \brief Density look up table values. */
	ThermoTables_Pressure[2], /*!< \brief Pressure look up table values. */
	ThermoTables_SoundSpeed2[2], /*!< \brief The speed of sound squared look up table values. */
	ThermoTables_Temperature[2], /*!< \brief Temperature look up table values. */
	ThermoTables_dPdrho_e[2], /*!< \brief Fluid derivative DpDd_e look up table values. */
	ThermoTables_dPde_rho[2], /*!< \brief Fluid derivative DpDe_d look up table values. */
	ThermoTables_dTdrho_e[2], /*!< \brief Fluid derivative DTDd_e look up table values. */
	ThermoTables_dTde_rho[2], /*!< \brief Fluid derivative DTDe_d look up table values. */
	ThermoTables_Cp[2], /*!< \brief Specific Heat Capacity at constant pressure look up table values. */
	ThermoTables_Mu[2], /*!< \brief Laminar Viscosity look up table values. */
	ThermoTables_dmudrho_T[2], /*!< \brief Fluid derivative DmuDrho_T look up table values. */
	ThermoTables_dmudT_rho[2], /*!< \brief Fluid derivative DmuDT_rho look up table values. */
	ThermoTables_Kt[2], /*!< \brief Thermal Conductivity look up table values. */
	ThermoTables_dktdrho_T[2], /*!< \brief Fluid derivative DktDrho_T look up table values. */
	ThermoTables_dktdT_rho[2]; /*!< \brief Fluid derivative DktDT_rho look up table values. */

	su2double Interpolation_Matrix[4][4]; /*!< \brief The (Vandermonde) matrix for the interpolation (bilinear) */

	su2double Interpolation_Coeff[4][4]; /*!< \brief Used to hold inverse of Interpolation_Matrix, and solution vector */
	int LowerI, UpperI, middleI, LowerJ, UpperJ, middleJ;/*!< \brief The i,j indexes (rho, P) of the position of the table search. Can be used as a restart for next search.*/
	int nTable_Zone_Stations[2]; /*!< \brief Number of nodes in the '2' zones of the LuT*/
	int nTable_Zone_Triangles[2]; /*!< \brief Number of triangles in the '2' zones of the LuT (must be triangles for now)*/
	vector<vector<int> > Table_Zone_Triangles[2]; /*!< \brief The triangles in each zone are stored as three intgers (the tree defining data-points)*/
	vector<vector<int> > Table_Zone_Edges[2]; /*!< \brief The edges in the '2' zones of the LuT*/
	vector<vector<int> > Table_Edge_To_Face_Connectivity[2];/*!< \brief Number of edges in the '2' zones of the LuT*/

public:

	CLookUpTable(void);
	/*!
	 * \brief Constructor the LUT by reading it in from a file.
	 * \param[in] Filename - The name of the (.rgp) file from which to load the table
	 */
	CLookUpTable(string Filename);

	/*!
	 * \brief Destructor of the class, primarily handling the dealloc of the KD_trees and LUT itself.
	 */
	virtual ~CLookUpTable(void);

	/*!
	 * \brief Set the Dimensional Thermodynamic State using Density and Internal Energy as inputs. Uses binary search in both directions separately.
	 * \param[in] rho - input Density (must be within LUT limits)
	 * \param[in] e   - input StaticEnergy (must be within LUT limits)
	 */
	void Get_Unique_Edges();



	void Get_Current_Points_From_TrapezoidalMap(CTrapezoidalMap *t_map,
			su2double x, su2double y);

	void SetTDState_rhoe(su2double rho, su2double e);

	/*!
	 * \brief Set the Dimensional Thermodynamic State using Pressure and Temperature as inputs. Uses binary search in both directions separately.
	 * \param[in] P - input Pressure (must be within LUT limits)
	 * \param[in] T - input Temperature (must be within LUT limits)
	 */

	void SetTDState_PT(su2double P, su2double T);

	/*!
	 * \brief Set the Dimensional Thermodynamic State using Pressure and Density as inputs. Uses binary search in both directions separately.
	 * \param[in] P   - input Pressure (must be within LUT limits)
	 * \param[in] rho - input Density  (must be within LUT limits)
	 */
	void SetTDState_Prho(su2double P, su2double rho);

	/*!
	 * \brief Set the Dimensional Internal Energy using Pressure and Density as inputs (nearly identical to SetTDState_Prho). Uses binary search in both directions separately.
	 * \param[in] P   - input Pressure (must be within LUT limits)
	 * \param[in] rho - input Density (must be within LUT limits)
	 */
	void SetEnergy_Prho(su2double P, su2double rho);

	/*!
	 * \brief Set the Dimensionless state using Enthalpy and Entropy as inputs. Uses KD_tree nearest neighbour searching, followed by zigzag searching for quad containing point
	 * \param[in] h - input Enthalpy (must be within LUT limits)
	 * \param[in] s - input Entropy  (must be within LUT limits)
	 *
	 */
	void SetTDState_hs(su2double h, su2double s);

	/*!
	 * \brief Set the Dimensional Thermodynamic State using Density and Temperature as inputs. Uses binary search in both directions separately.
	 * \param[in] rho - input Density (must be within LUT limits)
	 * \param[in] T   - input Temperature (must be within LUT limits)
	 *
	 */
	void SetTDState_rhoT(su2double rho, su2double T);

	/*!
	 * \brief Set the Dimensional Thermodynamic State using Pressure and Entropy as inputs. Uses binary search in both directions separately.
	 * \param[in] P - input Pressure (must be within LUT limits)
	 * \param[in] s - input Entropy (must be within LUT limits)
	 */

	void SetTDState_Ps(su2double P, su2double s);

	/*!
	 * \brief Calculate the inverse of a square matrix (e.g. the Vandermonde matrix) with pivoting Gaussian elimination
	 * \param[in] nDim - the dimension of the square block to invert
	 */

	void Gaussian_Inverse(int nDim);

	void Interpolate_2D_Bilinear(su2double x, su2double y,
			vector< su2double > *ThermoTables_X, vector< su2double > *ThermoTables_Y,
			std::string grid_var);

	/*!
	 * \brief Use the interpolation coefficients to interpolate a given thermodynamic variable property. (Must calculate the interpolation coefficients first)
	 * \param[in] interpolant_var - the name of the variable to be interpolated e.g Density
	 */

	su2double Interpolate_2D_Bilinear(vector< su2double > *ThermoTables_Z);

	/*!
	 * \brief Load the LUT table from a CFX file format. X axis must be Density, and Y axis pressure. Equal spacing not required.
	 * \param[in] filename - the name of the CFX file containing the table
	 */

	void LookUpTable_Malloc(int Index_of_Zone);
	void LookUpTable_Load_TEC(std::string filename);
	void NonDimensionalise_Table_Values();

	/*!
	 * \brief Print the table to a text file (for external inspection)
	 * This was used during the verification to print the table in a simpler format
	 * \param[in] filename - the name of the file where the LUT should be stored
	 */

	void LookUpTable_Print_To_File(char* filename);

	/*!
	 * \brief Records the current thermodynamic state of the fuid model to a file
	 * This was used during verification to record the results for a random set of input samples
	 * \param[in] fil - the name of the file to which to append the current state to
	 */

	void RecordState(char* file);

};
