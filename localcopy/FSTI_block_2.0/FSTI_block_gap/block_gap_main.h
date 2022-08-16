# ifndef __cblockgapmain__
# define __cblockgapmain__

# include "../FSTI_input/input.h"
# include "./gap_mesh.h"
# include "./scalar_field.h"
# include "./vector_field.h"
# include "../interpolation/interpl_grid.h"
# include "./p_file.h"
# include "./macro_geometry.h"
# include "./oil.h"
# include "./sforce.h"
# include "./dc_surf_mesh.h"

class cblock_gap_main
{

public:

	const input& in;			// reference to input structure
	p_file pfile;					// pressure file object

	// ---------------------------- public members --------------------------- //
	
	double t;										// simulation time
	double phi;									// shaft angle
	int revolution;							// current revolution
	double dt;									// simulation time step
	double step_angle;					// simulation step angle
	bool EHD_CB;								// EHD activated on cylinder block
	bool EHD_VP;								// EHD activated on valve plate
	bool TH_CB;									// Thermal activated on the cylinder block side
	bool TH_VP;									// Thermal activated on the valve plate side
	bool use_fsi_fb;						// use FSI term in the RE when solving for the force balance
	bool use_sqz_hs;						// use hydrostatic squeeze term in the RE
	bool use_sqz_hd;						// use hydrodynamic squeeze term in the RE
	
	interpl_grid CBf;						// fluid interpolation grid for cylinder block
	interpl_grid VPf;						// fluid interpolation grid for valve plate

	oil* lubricant;							// pointer to oil base class

	macro_geometry* macro_CB;		// macrogeometry for cylider block
	macro_geometry* macro_VP;		// macrogeometry for valve plate

	// mesh
	struct
	{
		gap_mesh reynolds;						// mesh used to solve the reynolds equation (2D)
		gap_mesh energy;							// mesh used to solve the reynolds equation (3D)
	} mesh;

	// fields
	struct
	{
		scalar_field p;					// pressure (for the force balance)
		scalar_field p_sol;			// pressure field with negative values
		scalar_field p_prev;
		scalar_field p_EHD;			// pressure after the solution of the EHD loop
		scalar_field rhoe;				// oil density east
		scalar_field rhow;				// oil density west
		scalar_field rhon;				// oil density north
		scalar_field rhos;				// oil density south 
		scalar_field mu;				// oil viscosity
		scalar_field mu2D;			// oil viscosity 2D field
		
		scalar_field pw;
		scalar_field pe;
		scalar_field ps;
		scalar_field pn;

		scalar_field vw;				// mass flows, west
		scalar_field ve;				// mass flows, east
		scalar_field vs;				// mass flows, south
		scalar_field vn;				// mass flows, north
		scalar_field g_pn;				// mass flows, north
		scalar_field g_pw;				// mass flows, west
		scalar_field g_pe;				// mass flows, east
		scalar_field g_ps;				// mass flows, south
		scalar_field mdotw;				// mass flows, west
		scalar_field mdote;				// mass flows, east
		scalar_field mdots;				// mass flows, south
		scalar_field mdotn;				// mass flows, north
		scalar_field fw;				// mass flows, west
		scalar_field fe;				// mass flows, east
		scalar_field fs;				// mass flows, south
		scalar_field fn;				// mass flows, north
		scalar_field rho;				// oil density
		scalar_field rho2D;			// oil density 2D
		scalar_field T;					// oil temperature
		scalar_field phid;			// dissipation field, instantaneous [W]
		scalar_field T_prev;
		scalar_field phid_avg;	// dissipation field, average [W]
		scalar_field Tcb;				// cylinder block gap surface temperature
		scalar_field Tvp;				// valve plate gap surface temperature
		scalar_field qcb;				// heat flow to cylinder block
		scalar_field qcb_avg;		// heat flow to cylinder block (average)
		scalar_field qcb_prog;	// heat flow to cylinder block progressive over the simulation
		scalar_field qvp;				// heat flow to valve plate
		scalar_field qvp_avg;		// heat flow to valve plate (average)
		scalar_field qvp_prog;	// heat flow to valve plate progressive over the simulation
		vector_field Vc;				// oil velocity, couette
		vector_field Vp;				// oil velocity, poiseuille
		vector_field V;					// oil velocity, total
	} fields;

	// film thickness
	struct
	{
		double h1, h2, h3;						// block position
		double dz, rx,ry;							// shift and rotation
		scalar_field h;								// film thickness
		scalar_field h_previous;								// film thickness - Rene
		scalar_field hcb;							// film thickness (cylinder block side)
		scalar_field hvp;							// film thickness (valveplate side)
		scalar_field hcb_rigid;				// block gap surface position (rigid)
		scalar_field hvp_rigid;				// valveplate gap surface position (rigid)
		scalar_field dhcb_thermal;		// thermal deformation cylinder block
		scalar_field dhvp_thermal;		// thermal deformation valve plate
		scalar_field dhcb_EHD;				// EHD deformation cylinder block
		scalar_field dhvp_EHD;				// EHD deformation valve plate
		scalar_field contact;					// possible theoretical contact between block and valve plate
		double hmin;
	} film;

	// squeeze motion
	struct
	{
		double dhdt1, dhdt2, dhdt3;					// block squeeze
		double d_dhdt1, d_dhdt2, d_dhdt3;		// block squeeze (contact corretion)
		scalar_field dhdt;									// squeeze motion field
		scalar_field dhdt_hs;								// squeeze due to change in hydrostatic deformation over time
		scalar_field dhdt_hd;								// squeeze due to change in hydrodynamic deformation over time
		scalar_field dhdt_total;								// squeeze due to change in hydrodynamic deformation over time
	} squeeze;

	// forces
	struct
	{
		dc_surf_mesh DC;																// object to store the dc geometry
		std::vector<double> p_cb;												// list of pressure in the elements defined in mesh.cb_gap_surf
		std::vector<std::vector<double>> pistonforces;	// external forces on the two reference points of piston, A and B four components for each piston: [FAx, FAy, FBx, FBy]
		double xFBz;																		// x position of the resulting axial force on block
		double yFBz;																		// y position of the resulting axial force on block
		double MKBx;																		// x moment on block due to piston loads
		double MKBy;																		// y moment on block due to piston loads
		double FfBz;																		// Axial fluid force
		double MfBx;																		// Fluid x moment
		double MfBy;																		// Fluid y moment
		double Ffluid[3];																// Equivalent fluid forces system
		double FBz;																			// Total external force on block from pistons, axial direction
		double FKBx;																		// Total external force on block from pistons, x direction
		double FKBy;																		// Total external force on block from pistons, y direction
		double MBx;																			// Extarnal x moment
		double MBy;																			// Extarnal y moment
		double MBx_spline;															// Moment  on the block in x direction due to the spline joint
		double Fext[3];																	// Equivalent external forces system
		double dF[3];																		// Force imbalance
		double Fc;																			// contact axial force
		double Mcx;																			// contact x moment
		double Mcy;																			// contact y moment
		double MTBz;																		// Block torque loss
		sforce FTK;																			// piston friction force
		sforce FTG;																			// slipper friction force
	} forces;

	// pistons
	struct
	{
		double z0;								// distance from block's reference system to the sealing land
		double phi_odp;						// position of the outer dead point
		double phi_idp;						// position of the inner dead point
		double delta_psi;					// offset angle between 0 and the odp (if !+ 0 only with cross angle)
		double sk_odp;						// piston position @ outer dead point
		double sk_idp;						// piston position @ inner dead point
		double Hk;								// piston stroke for current beta
		double Hk_betamax;				// piston stroke for max beta
		double lZ0;								// distance from the piston end to the end of the displ chamber (is a function of beta!)
		double pos_A;
		double pos_B;
		double radial_clearance;
		double T;
		std::vector<double> phi;	// vector of angle for each pistons
		std::vector<double> pDC;	// displacement chamber pressure for each cylinder
		double pLP;								// high pressure port
		double pHP;								// low pressure port
		std::vector<double> sK;		// vector of positions for each piston
		std::vector<double> vK;		// vector of velocity for each piston
		std::vector<double> aK;		// acceleration of each piston
		std::vector<double> zK;		// position of the ball jont center with respect to the block's reference system
		std::vector<double> lvar;							//Variable gap length [m] (of the piston-cylinder interface)
	} pistons;
	
	// ------------------------ public member functions ---------------------- //

	cblock_gap_main(const input& in);																	// constructor
	~cblock_gap_main();																								// destructor
	void initialize();																								// initialize the gap structure
	void initialize_pistons();																				// initialize the pistons structure
	void rotate(double phi_rad);																			// rotate mesh and update all fields
	void set_film_thickness();																				// set the film thickness, use true to exclude deformation effect
	double hmin() const;																							// get min film thickness
	double hmax() const;																							// get max film thickness
	double havg() const;																							// get average film thickness
	double calcFTK(double phi_rad);
	void set_squeeze();																								// set the squeeze motion field dhdt
	void set_pressure();																							// set pressure boundaries
	void set_temperature();																						// set the temperature boundaries
	void update_h_previous();
	void update_pistons(double phi_rad);															// update pistons structure with new angular position
	void update_mesh();																								// update the mesh with the new film thickness
	void update_density();																						// update oil density
	void update_viscosity();																					// update oil viscosity
	void update_Tcb();																								// update block surface temperature from fluid grid interpolation
	void update_Tvp();																								// update valve plate surface temperature from fluid grid interpolation
	void update_dhcb_thermal();																				// update the block thermal deflection from fluid grid interpolation
	void update_dhvp_thermal();																				// update the block thermal deflection from fluid grid interpolation
	void calc_fluid_forces();																					// calculate pressure force and moments on block
	void calc_piston_forces();																				// calculate forces from piston to cylider bore
	void calc_dhdt_total();	
	void calc_external_forces();																			// calculate external loads on block
	scalar_field calc_cb_heatflux();																	// calculate the heat flux towards the block
	scalar_field calc_vp_heatflux();																	// calculate the heat flux towards the valve plate
	double get_qcb_avg();																							// calculate the average heat flux over the cylinder block
	double get_qvp_avg();																							// calculate the average heat flux over the valve plate
	void apply_heat_flux(interpl_grid* str, const char* body);				// apply q_prog on the interpolation grid coming from thermal solver
	double calc_friction_moment();																		// calculate the torque loss due to friction
	double calc_leakage();																						// calculate the leakage
	void calc_dF();																										// get the force embalance
	void calc_mass_flow();																						// calculate velocity at the faces and at the cell centroid
	void define_interpl_grids();																			// define the fluid grids used for interpolation
	void interpl_pressure();																					// interpolate pressure to 
	void write_pcb_vtk(const char* name);															// write block gap surface and interpolated pressure
	void write_light_vtk(const char* name, double zscale = 1.0);			// write light vtk output
	void write_vtk(const char* name, double zscale = 1.0);						// write full vtk output

};

# endif