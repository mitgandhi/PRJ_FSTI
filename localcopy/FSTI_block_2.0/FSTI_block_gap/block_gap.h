// ------------------------------------------------------------------------- //
// the cblock_gap class contains:                                            //
//		                                                                       //
//		* gap_main, which contains all the gap data, structures and            //
//			functions                                                            //
//	  * reynolds, used to solve the Reynolds equation EHD and non EHD        //
//			versions                                                             // 
//		* energy, used to solve for the energy equation in the gap             //
//		* the multiroot solver, used to caluclate the block force balance      //
// ------------------------------------------------------------------------- //

# ifndef __cblockgap__
# define __cblockgap__

# include "../FSTI_input/input.h"
# include "./block_gap_main.h"
# include "./reynolds.h"
# include "./energy.h"
# include "../FSTI_block_thermal/te_solver.h"
# include <gsl/gsl_vector.h>
# include <gsl/gsl_multiroots.h>
# include <vector>
# include <iostream>
# include <fstream>


class cblock_gap
{
	
	const input& in;

	double tol;					// root finding tolerance

	// multiroot solver stuff
	const gsl_multiroot_fdfsolver_type *T;
	gsl_multiroot_fdfsolver *s;
	

public:

	// -------------------------------- members ------------------------------ //
	
	double phi;							// shaft angle
	int revolution;						// current revolution
	
	cblock_gap_main gap_main;			// main block gap object
	
	Reynolds reynolds;					// Reynolds object, to solve the Reynolds equation
	
	Energy energy;						// Energy object, to solve the Energy equation

	te_solver* te_CB;					// Thermo-Elastic solver for cylinder block
	
	te_solver* te_VP;					// Thermo-Elastic solver for valve plate

	// --------------------------- member functions -------------------------- //
	
	cblock_gap(const input& _in);														// constructor
	~cblock_gap() {}																				// destructor
	int get_block_balance(bool reset = false);							// force balance loop (multidim root finding)
	void get_contact_velocities();													// get velocity corrections due to contact
	//void get_contact_forces();													// get force corrections due to contact
	void check_gap();																				// check the thermal analysis
	void check_EHD();																				// check EHD analysis
	
	void check_thermal();																		// check the thermal analysis
	void get_external_loads();															// calculate the external loads on the block over one shaft revolution
	void get_balance_factors();															// get the balance factors hydrostatic/external 
	void get_piston_kinematics();														// calculate piston position respect the main reference system over one shaft revolution
	void get_oil_properties();															// calculate oil properties as a function of p and T
	void write_cb_influgen(const char* file);
	void write_vp_influgen(const char* file);
	void write_vp_pressure_field();
	void writeInflugenInputs();

	void writeVTK_block();																		// check the thermal analysis
	
};

# endif