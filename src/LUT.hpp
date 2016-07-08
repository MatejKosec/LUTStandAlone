#define su2double double
#include <string>
/*!
	 * \brief The struct used for the KD_tree implementation, each KD_node represents a tree branch.
	 * The KD_tree is 2D, as all SetTD state functions are based on thermo-pairs.
	 * All other variable values at a given location in the tree can be obtained
	 * by using the stored indexes.
	 *
*/
struct KD_node
{
	int Branch_Splitting_Direction, /*!< \brief The depth of the branch within the tree, even numbers are split in x, odd in y */
	Branch_Dimension, /*!< \brief The number of points contained on the branch */
	*Flattened_Point_Index; /*!< \brief The flattened 2D index of the original LUT. Used to access variables other */
	su2double * x_values, /*!< \brief The x_values of the data. Values are sorted if splitting direction is even.*/
	* y_values; /*!< \brief  The (sorted) y_values of the data. Values are sorted if splitting direction is odd. */
	KD_node* upper; /*!< \brief The tree-branch on the next level of the tree containing upper 50 percentile. Based on x_values for branches of even depth and y_values for odd. */
	KD_node* lower; /*!< \brief The tree-branch on the next level of the tree containing lower 50 percentile. Based on x_values for branches of even depth and y_values for odd. */
};

/*!
 * \class CThermoList
 * \brief The class which holds the thermodynamic variables associated with a particular condition of the gas.
 * \author: M. Kosec, A.Rubino, S.Vitale
 * \version 4.1.2 "Cardinal"
 */
class CThermoList {
public:
	su2double StaticEnergy,	/*!< \brief Internal Energy. */
	Entropy,  				/*!< \brief Entropy. */
	Enthalpy,					/*!< \brief Enthalpy required as separate variable for use in HS tree. */
	Density,  				/*!< \brief Density. */
	Pressure, 				/*!< \brief Pressure. */
	SoundSpeed2, 		  /*!< \brief The speed of sound squared. */
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
	/*!
		 * \brief Print the thermodynamic variables contained in the instance (e.g. to read-out verification data)
  */
	void CThermoList_Print(void);

};

/*!
 * \class CLookUpTable
 * \brief Child class of CFluidModel which uses thermodynamic variables read in from a table in .rgp format
 * \author: M.Kosec, A. Rubino, S.Vitale
 * \version 4.1.2 "Cardinal"
 */
class CLookUpTable {

protected:
	CThermoList **ThermoTables;/*!< \brief The 2D array used to hold the values of thermodynamic properties from the LUT*/

	su2double Interpolation_Matrix[4][4]; /*!< \brief The (Vandermonde) matrix for the interpolation (bilinear) */
	su2double Interpolation_Coeff[4][4]; /*!< \brief Used to hold inverse of Interpolation_Matrix, and solution vector */
	short iIndex, jIndex; /*!< \brief The i,j indexes (rho, P) of the position of the table search. Can be used as a restart for next search.*/
	int Table_Pressure_Stations;/*!< \brief The pressure dimensions of the table */
	int Table_Density_Stations; /*!< \brief The density dimensions of the table */
	KD_node *HS_tree; /*!< \brief The pointer to the root of the KD tree for the HS thermo-pair.*/
	su2double StaticEnergy,			/*!< \brief Internal Energy. */
	Entropy,  				/*!< \brief Entropy. */
	Enthalpy,					/*!< \brief Enthalpy required as separate variable for use in HS tree. */
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

	su2double StaticEnergy_Table_Limits[2]; /*!< \brief The [min,max] values of the StaticEnergy values in the LUT */
	su2double Entropy_Table_Limits[2]; /*!< \brief The [min,max] values of the Entropy values in the LUT */
	su2double Enthalpy_Table_Limits[2]; /*!< \brief The [min,max] values of the Enthalpy values in the LUT */
	su2double Density_Table_Limits[2];/*!< \brief The [min,max] values of the Density values in the LUT */
	su2double Pressure_Table_Limits[2];/*!< \brief The [min,max] values of the Pressure values in the LUT */
	su2double SoundSpeed2_Table_Limits[2]; /*!< \brief The [min,max] values of the SoundSpeed squared values in the LUT */
	su2double Temperature_Table_Limits[2];/*!< \brief The [min,max] values of the Temperature values in the LUT */
	su2double dPdrho_e_Table_Limits[2];/*!< \brief The [min,max] values of the dPdrho_e  values in the LUT */
	su2double dPde_rho_Table_Limits[2];/*!< \brief The [min,max] values of the dPde_rho  values in the LUT */
	su2double dTdrho_e_Table_Limits[2];/*!< \brief The [min,max] values of the dPde_rho  values in the LUT */
	su2double dTde_rho_Table_Limits[2];/*!< \brief The [min,max] values of the dPde_rho  values in the LUT */
	su2double Cp_Table_Limits[2];/*!< \brief The [min,max] values of the dPde_rho  values in the LUT */
	su2double Mu_Table_Limits[2];/*!< \brief The [min,max] values of the dPde_rho  values in the LUT */
	su2double dmudrho_T_Table_Limits[2];/*!< \brief (UNUSED) The [min,max] values of the dmudrho_T  values in the LUT */
	su2double dmudT_rho_Table_Limits[2];/*!< \brief (UNUSED) The [min,max] values of the dmudT_rho  values in the LUT */
	su2double Kt_Table_Limits[2];/*!< \brief The [min,max] values of the Kt values in the LUT */
	su2double dktdrho_T_Table_Limits[2];/*!< \brief (UNUSED) The [min,max] values of the dktdrho_T values in the LUT */
	su2double dktdT_rho_Table_Limits[2];/*!< \brief (UNUSED) The [min,max] values of the dktdT_rho values in the LUT */
	//Nearest neighbour's i and j indexes
	int *Nearest_Neighbour_iIndex;/*!< \brief An array which holds the i (rho) indexes of the points used in the interpolation (usually Neighbours)*/
	int *Nearest_Neighbour_jIndex;/*!< \brief An array which holds the j (P) indexes of the points used in the interpolation (usually Neighbours)*/

public:

	/*!
	 * \brief default Constructor of the class.
	 */
	CLookUpTable(void);

	/*!
	 * \brief Constructor the LUT by reading it in from a file.
	 * \param[in] Filename - The name of the (.rgp) file from which to load the table
	 */
	CLookUpTable(char* Filename);

	/*!
	 * \brief Destructor of the class, primarily handling the dealloc of the KD_trees and LUT itself.
	 */
	virtual ~CLookUpTable(void);

	/*!
	 * \brief Recursively build a 2D KD_tree (used for HS pair, or unstructured approaches)
	 * \param[in] depth    - the depth of the tree is used to determine along which axis to split the values
	 * \param[in] x_values - dynamic array containing all values of the first thermodynamic variable at a given tree depth
	 * \param[in] y_values - dynamic array containing all values of the second thermodynamic variable at a given tree depth
	 * \param[in] i_values - the flattened LUT index to which each x,y combination corresponds (used to select interpolation quad)
	 * \param[in] dim      - the number of points held in this branch of the kd tree
	 */
	struct KD_node* KD_Tree(su2double* x_values, su2double* y_values, int* i_values, int dim, int depth);

	/*!
		 * \brief The squared Euclidian distance between the x,y search values and the median of the current tree branch
		 * \param[in] x - the x value (i.e. location) of the point being searched for
		 * \param[in] y - the y value (i.e. location) of the point being searched for
		 * \param[in] branch   - the branch of the KD_tree from which to take the median value
		 */
	su2double Dist2_KD_Tree (su2double x, su2double y, KD_node *branch);

	/*!
	* \brief Recursively descend through the tree and free up the x_values,y_values, and i_values dynamic arrays
	* \param[in] root - the branch of the tree into which to descend (should usually start with the highest level e.g. this->HS_tree)
	*/
	void free_KD_tree(KD_node* root);
	/*!
		* \brief Recursively search through the KD tree for the N nearest Neighbours to the search thermo-pair values. This includes reversing the search once the a branch with dimension 1 is found.
		* \param[in] N - the number of nearest neighbours desired
		* \param[in] thermo1 - the x_value of the thermo-pair for which to search
		* \param[in] thermo2 - the y_value of the thermo-pair for which to search
		* \param[in] root    - pointer of the KD_tree branch from which to start the search (recursively, typically root e.g this->HS_tree)
		* \param[in] best    - a shared dynamic array of the N best distances encountered in the search.
	*/
	void N_Nearest_Neighbours_KD_Tree(int N, su2double thermo1, su2double thermo2, KD_node *root, su2double* best_dist);

	/*!
	 * \brief Set the Dimensional Thermodynamic State using Density and Internal Energy as inputs. Uses binary search in both directions separately.
	 * \param[in] rho - input Density (must be within LUT limits)
	 * \param[in] e   - input StaticEnergy (must be within LUT limits)
	 */
	void SetTDState_rhoe (su2double rho, su2double e );

	/*!
	 * \brief Set the Dimensional Thermodynamic State using Pressure and Temperature as inputs. Uses binary search in both directions separately.
	 * \param[in] P - input Pressure (must be within LUT limits)
	 * \param[in] T - input Temperature (must be within LUT limits)
	 */

	void SetTDState_PT (su2double P, su2double T );

	/*!
	 * \brief Set the Dimensional Thermodynamic State using Pressure and Density as inputs. Uses binary search in both directions separately.
	 * \param[in] P   - input Pressure (must be within LUT limits)
	 * \param[in] rho - input Density  (must be within LUT limits)
	 */
	void SetTDState_Prho (su2double P, su2double rho );

	/*!
	 * \brief Set the Dimensional Internal Energy using Pressure and Density as inputs (nearly identical to SetTDState_Prho). Uses binary search in both directions separately.
	 * \param[in] P   - input Pressure (must be within LUT limits)
	 * \param[in] rho - input Density (must be within LUT limits)
	 */
	void SetEnergy_Prho (su2double P, su2double rho );

	/*!
	 * \brief Set the Dimensionless state using Enthalpy and Entropy as inputs. Uses KD_tree nearest neighbour searching, followed by zigzag searching for quad containing point
	 * \param[in] h - input Enthalpy (must be within LUT limits)
	 * \param[in] s - input Entropy  (must be within LUT limits)
	 *
	 */
	void SetTDState_hs (su2double h, su2double s );

	/*!
	 * \brief Set the Dimensional Thermodynamic State using Density and Temperature as inputs. Uses binary search in both directions separately.
	 * \param[in] rho - input Density (must be within LUT limits)
	 * \param[in] T   - input Temperature (must be within LUT limits)
	 *
	 */
	void SetTDState_rhoT (su2double rho, su2double T );

	/*!
	 * \brief Set the Dimensional Thermodynamic State using Pressure and Entropy as inputs. Uses binary search in both directions separately.
	 * \param[in] P - input Pressure (must be within LUT limits)
	 * \param[in] s - input Entropy (must be within LUT limits)
	 */

	void SetTDState_Ps (su2double P, su2double s );

	/*!
		 * \brief Calculate the inverse of a square matrix (e.g. the Vandermonde matrix) with pivoting Gaussian elimination
		 * \param[in] nDim - the dimension of the square block to invert
	*/

	void Gaussian_Inverse(int nDim);

	/*!
			 * \brief Calculate the bilinear interpolation coefficients for a quad with arbitrary skew.
			 *  The entails building the Vandermonde matrix, inverting it, transposing it, and dot product with the search values of x, and y.
			 *  This formulation with the transpose means that the coefficients depend only on the x,y cooridinate of the search and
			 *  not on the thermodynamic variable being interpolated. Thus, the same coefficients can be used across
			 *  the interpolation of all desired thermodynamic properties.
			 * \param[in] x - the x value used to set the thermodynamic state. (e.g. rho in rhoe)
			 * \param[in] x - the y value used to set the thermodynamic state. (e.g. rho in e)
			 * \param[in] grid_var - the pair of thermodynamic variables which define the grid i.e. the interpolation quad. (e.g. RHOE for rhoe)
		*/

	void Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(su2double x, su2double y, std::string grid_var);

	/*!
				 * \brief Use the interpolation coefficients to interpolate a given thermodynamic variable property. (Must calculate the interpolation coefficients first)
				 * \param[in] interpolant_var - the name of the variable to be interpolated e.g Density
			*/

	su2double Interpolate_2D_Bilinear(std::string interpolant_var);

	/*!
	* \brief Use the inverse distances in a Shephard's interpolation.
	* This was initially used for the HS tree but did not produce satisfactory accuracy
	* \param[in] N - the number of points involved in the interpolation, typically neighbors found using KD_tree
	* \param[in] interpolant_var - the distance between the sample points and the desired points. May be raised to desired power.
  */

	su2double Interp2D_Inv_Dist(int N, std::string interpolant_var, su2double* dist);

	/*!
		* \brief Load the LUT table from a CFX file format. X axis must be Density, and Y axis pressure. Equal spacing not required.
		* \param[in] filename - the name of the CFX file containing the table
	  */

	void LookUpTable_Load_CFX(std::string filename);

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




