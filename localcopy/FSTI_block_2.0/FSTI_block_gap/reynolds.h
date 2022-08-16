// ------------------------------------------------------------------------- //
// ------------------- class to solve the Reynolds equation ---------------- //
// ------------------------------------------------------------------------- //

# ifndef __reynolds__
# define __reynolds__

# include "./block_gap_main.h"
# include "../linear_solvers/PBiCG.h"
# include "./influ_set.h"
# include "../interpolation/interpl_grid.h"
# include "../interpolation/interpolation.h"
# include "../FSTI_block_gap/body_type.h"
# include <vector>

class Reynolds
{

	cblock_gap_main&	gap;	// reference to the main gap structure
	const gap_mesh& msh;

	// -------------- general info -------------- //
	
	double current_relax;		// current relaxation value
	double hmin;						// minimum film thickness
	double omega;						// pump speed
	double dt;							// simulation time step
	bool EHD_CB;						// flag for EHD on block
	bool EHD_VP;						// flag for EHD on valve plate
	public:
	bool hydrodynamic;			// flag for hydrostatic / hydrodynamic solution
	private:
	double IM_cb_avg;				// average value of the IM for each fluid cell of the CB (1Pa pressure)
	double IM_vp_avg;				// average value of the IM for each fluid cell of the VP (1Pa pressure)

	// ------------------------ Influence matrices ---------------------------- //
	
	influ_set* IMset;

	// --------------------- interpolation structures ------------------------ //

	interpl_grid& CBf;
	interpl_grid& VPf;
	interpolation CB_interpl;				// interpolation object for CB
	interpolation VP_interpl;				// interpolation object for VP
	
	// ------------- references to gap fields --------------- //

	scalar_field& p;										// pressure
	scalar_field& p_EHD;								// pressure field for the EHD loop
	scalar_field& p_sol;								// pressure field used for the solution (negative values are allowed)
	scalar_field& mu;										// viscosity
	scalar_field& rho;									// density
	const scalar_field& dhdt;						// squeeze
	scalar_field& dhdt_hs;							// squeeze due to the change in hydrostatic deformation
	scalar_field& dhdt_hd;							// squeeze due to the change in hydrostatic deformation
	scalar_field& contact;							// contact field
	
	std::vector<double>& pDC;						// pressure in each DC 
	const double& pHP;									// high pressure
	const double& pLP;									// low pressure

	// cb
	const scalar_field& hcb;						// block side film thickness
	const scalar_field& hcb_rigid;						// block side film thickness
	scalar_field& dhcb_EHD;							// block surface EHD deformation

	// vp
	const scalar_field& hvp;									// valve plate side film thickness
	scalar_field& dhvp_EHD;							// valve plate surface EHD deformation
	
	
	// ------------------ local fields ------------------- //

	scalar_field	cdef0_cb;					// elementary deformation from each CB element
	scalar_field	cdef_cb;					// elementary deformation from each CB element
	scalar_field	cdef_vp;					// elementary deformation from each VP element

	scalar_field dhcb_hs;						// cb hydrostatic deformation, reference angle, current timestep		
	scalar_field dhcb0_hs;					// cb hydrostatic deformation, reference angle, current timestep		
	scalar_field dhcb0_hs_prev;			// cb hydrostatic deformation, reference angle, previous timestep
	scalar_field dhdt_cb_hs;				// squeeze introduced by the change in hs def on the cb
	scalar_field dhvp_hs;						// valveplate hydrostatic deformation, current timestep
	scalar_field dhvp_hs_prev;			// valveplate hydrostatic deformation, previous timestep
	scalar_field dhdt_vp_hs;				// squeeze introduced by the change in hs def on the vp
	
	scalar_field dhcb0_hd;
	scalar_field dhcb0_hd_prev;
	scalar_field dhcb_hd_prev;
	scalar_field dhcb_hd;
	scalar_field dhvp_hd;
	scalar_field dhvp_hd_prev;
	scalar_field dhdt_cb_hd;
	scalar_field dhdt_vp_hd;

	// ------------------- coefficients ----------------- //

	double Ds;		// diffusive coeff on south face
	double Dn;		// diffusive coeff on north face
	double Dw;		// diffusive coeff on west face
	double De;		// diffusive coeff on east face
	double S;			// source

	// -------------------- linear system ----------------- //

	int size;																		// system size
	gmm::col_matrix<gmm::wsvector<double>> A;		// fv matrix
	std::vector<double> b;											// known term
	PBiCG solver;																// linear solver

	// ------------------- private member function ------------------- // 

	void get_D(int id);															// get diffusive coeffs for cell id
	void get_S(int id);															// get source for cell id
	void discretize(double alpha, bool EHD_sqz);		// discretise the reynolds equation
	void update_solution();													// update the pressure field
	void update_interpl(const scalar_field& _p);		// update the interpolated pressure field
	void get_cb_gap_def();													// get the deformation of the gap, CB side
	void get_vp_gap_def();													// get the deformation of the gap, VP side
	void get_dh_hs_cb();														// get the hydrostatic deformation for the CB
	void get_dh_hs_vp();														// get the hydrostatic deformation for the CB

public:
	
	Reynolds(cblock_gap_main&);													// constructor
	~Reynolds();																				// destructor
	void init_def_map(body_type B);											// initialize the deformation maps
	void update_dimension();														// update the linear system dimension
	int solve_rigid(bool EHD_sqz = false);							// solve reynolds, no EHD
	int solve_EHD();																		// solve reynolds, with EHD option
	bool EHD() const;																		//  return if the EHD option is on
	void delete_matrices();															// free the memory with the influence matrices
	void load_matrices();																// load in memory the influence matrices
	const interpl_grid& get_CBs() const;								// get the CB structure grid 
	const interpl_grid& get_VPs() const;								// get the VP structure grid 
	double get_IM_cb_avg() const { return IM_cb_avg; }	// return the average IM value for the CB fluid cell
	double get_IM_vp_avg() const { return IM_vp_avg; }	// return the average IM value for the VP fluid cell

	void temp();
};	

# endif