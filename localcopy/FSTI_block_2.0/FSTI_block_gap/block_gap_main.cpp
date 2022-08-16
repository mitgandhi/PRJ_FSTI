# include "./block_gap_main.h"
# include "../FSTI_Block_dll/log.h"
# include "./block_limits.h"
# include "./macro_geometry.h"
# include "../FSTI_block_thermal/matrix.h"
# include "../interpolation/interpolation.h"
# include "../FSTI_block_thermal/matrix.h"
# include <ANN/ANN.h>
# include <fstream>

using namespace std;


# define VTK_3D_ELM_NDS 8
# define VTK_3D_ELM_TYPE 12
# define VTK_2D_ELM_NDS 4
# define VTK_2D_ELM_TYPE 9
# define DISSIPATION

const double pi = 3.141592653589793238;

extern gaplog Log;

// ------------------------------------------------------------------------- //
cblock_gap_main::cblock_gap_main(const input& _in) : in(_in) 
{
	// initialize the cblock_gap_main object
	initialize();
}
// ------------------------------------------------------------------------- //
cblock_gap_main::~cblock_gap_main()
{
	if(macro_CB != 0)
		delete macro_CB;
	if(macro_VP != 0)
		delete macro_VP;

	delete lubricant;

}
// ------------------------------------------------------------------------- //
void cblock_gap_main::initialize_pistons()
{
	// ------------------------ copy from caspar_input ----------------------- //

	int z = in.data.operating_conditions.npistons;
	double rB = 0.50*in.data.geometry.dB;
	double rK = 0.50*in.data.geometry.dK;
	double rZ = 0.50*in.data.geometry.dZ;
	double beta = in.data.operating_conditions.beta;					
	double betamax = in.data.operating_conditions.betamax;		
	double gamma = in.data.geometry.gamma;		
	double lZ0 = in.data.geometry.lZ0;	
	double lF = in.data.geometry.lF;	
	double lK = in.data.geometry.lK;	
	double lKG = in.data.geometry.lKG;
	double le = in.data.geometry.le;
	double lch = in.data.geometry.lch;												
	double lengthB = in.data.geometry.lengthB;	
	double lengthcanalB = in.data.geometry.lengthcanalB;
	double R = 0.5*in.data.geometry.d_spherical;
	double r_gap_out = 0.5*in.data.geometry.d_gap_out;			
	double J = in.data.geometry.offset_J;						//Offset J 
	double K = in.data.geometry.offset_K;						//Offset K 
	double T;													//Translation of piston due to J and K offsets.
	// experimental delete when done
	//double gamma_z0 = in.data.options_block.general.gamma_z0;
	// ----------------------------------------------------------------------- //

	// calculate derived pistons variables (fixed throughout the simulation)
	
	// angular position of the outer dead point
	pistons.phi_odp = (beta == 0) ? 0 : -atan(tan(gamma)/sin(beta));
	pistons.phi_idp = pi + pistons.phi_odp;
	// offset in sk(0) due to gamma angle
	pistons.delta_psi = rB*tan(beta)*(1 - cos( pistons.phi_odp)) + rB*tan(gamma)*sin( pistons.phi_odp)/cos(beta);

	// geometric stroke 
	pistons.sk_odp = -rB*tan(beta)*(1-cos( pistons.phi_odp)) - rB*tan(gamma)*sin(pistons.phi_odp)/cos(beta) + pistons.delta_psi;
	pistons.sk_idp = -rB*tan(beta)*(1-cos( pistons.phi_idp)) - rB*tan(gamma)*sin(pistons.phi_idp)/cos(beta) + pistons.delta_psi;
	pistons.Hk = fabs( pistons.sk_odp -  pistons.sk_idp);

	// stroke at max beta
	pistons.Hk_betamax = 
		fabs
		(
			-rB*tan(betamax)*(1-cos(pistons.phi_odp)) - rB*tan(gamma)*sin(pistons.phi_odp)/cos(betamax)
			+ rB*tan(betamax)*(1-cos(pistons.phi_idp)) + rB*tan(gamma)*sin(pistons.phi_idp)/cos(betamax)
		);

	// account for translation of piston due to J and K offsets 
	T= tan(beta)*(J-K*tan(0.5*beta)) - tan(betamax)*(J-K*tan(0.5*betamax));
	
	// correct lZ0 according with the current beta angle
	pistons.lZ0 = lZ0 - 0.5*(pistons.Hk_betamax - pistons.Hk);

	// distance from front flat piston surface to bushing end when piston is at ODC
	double delta_A = lengthB - lengthcanalB - pistons.lZ0 - lF - le - T;
	// distance from piston head center to cylinder block beginning edge 
	// when piston is at ODC
	double delta_B = lK - delta_A - lF - le;
	
	// distance between cylinder block origin and cylinder block beginning edge
	double delta_0 = delta_B - 0.5*fabs(pistons.sk_odp - pistons.sk_idp);

	// distance from the gap surface of the block to the block origin, which
	// is supposed to be the center of instantaneous rotation
	pistons.z0 = delta_0 + lengthB;

	//radial clearance 
	pistons.radial_clearance = rZ-rK;
	pistons.T = T;

	// correct z0 with the inclination of the swash plate due to the cross angle
	if((in.data.geometry.gamma > 0))// && (gamma_z0))
	{
		Log << "\nCalculating new z0, accounting for the ""cross-angle""... ";
		// distance between the ideal origin and the ball joint center of the
		// piston at the ODP
		double d_odp = rB*sin(fabs(pistons.phi_odp))*tan(in.data.geometry.gamma);
		
		// maximum axial distance respect the ideal position due to the gamma angle
		double d_max = rB*tan(in.data.geometry.gamma);

		// add the offset
		pistons.z0 += d_max - d_odp;
	}
	
	
	// initialize pistons
	pistons.phi.resize(z, 0);
	pistons.pDC.resize(z, 0);
	pistons.sK.resize(z, 0);
	pistons.vK.resize(z, 0);
	pistons.aK.resize(z, 0);
	pistons.zK.resize(z, 0);
	pistons.lvar.resize(z, 0);
	
	// initialize pistons structure at 0 angle
	update_pistons(0);
}
// ------------------------------------------------------------------------- //
void cblock_gap_main::initialize()
{

	Log << "// - Initialization of the Cylinder block / Valve plate interface analysis - //\n\n";
	
	// ------------------------ copy from caspar_input ----------------------- //

	int oiltype = in.data.oil.general.oiltype;
	int z = in.data.operating_conditions.npistons;
	double hB1 = in.data.options_block.position.hB1;
	double hB2 = in.data.options_block.position.hB2;
	double hB3 = in.data.options_block.position.hB3;
	double stpang = in.data.options_block.general.step_angle;
	double speed = in.data.operating_conditions.speed;
	const char* pfile_name = in.data.operating_conditions.pModuleFile.c_str();
	double Tcase = in.data.operating_conditions.T_Leak;
	EHD_CB = in.data.options_block.general.EHD_CB;
	EHD_VP = in.data.options_block.general.EHD_VP;
	TH_CB = in.data.options_block.general.Thermal_CB;
	TH_VP = in.data.options_block.general.Thermal_VP;
	use_fsi_fb = in.data.options_block.numeric.use_fsi_fb;
	use_sqz_hd = in.data.options_block.numeric.use_sqz_hd;
	use_sqz_hs = in.data.options_block.numeric.use_sqz_hs;
	string DC_mesh = in.data.geometry.DC_mesh;
	int macro_CB_opt = in.data.options_block.general.macro_CB;
	int macro_VP_opt = in.data.options_block.general.macro_VP;
	int Fext_method = in.data.options_block.numeric.Fext_method;

	// ----------------------------------------------------------------------- //

	// --------- set the lubricant ----------- //

	Log << "\nInitializing the oil model ... ";
	if(oiltype == 0)
		lubricant = new constant_oil(in);		// use constant properties
	else if(oiltype == 1)
		lubricant = new user_defined_oil(in);	// user defined
	else if(oiltype == 2)
		lubricant = new HLP32_oil(in);			// HLP32
	else if(oiltype == 3)
		lubricant = new user_defined2_oil(in);	// user defined2

	//else if(oiltype == 3)
	//	lubricant = new skydrol_oil(in);		// skydrol
	else if(oiltype == 4)
		lubricant = new red_oil(in);			// red oil triumph
	else if(oiltype == 5)
		lubricant = new SAE10W(in);				// SAE10W
	else if(oiltype == 6)
		lubricant = new milh5606_oil(in);		// milh5606 oil triumph
	else if(oiltype == 7)
		lubricant = new ExxonDTE10Excel32(in);	// ExxonDTE10Excel32 
	else if(oiltype == 8)
		lubricant = new ISO46_oil(in);			// ISO46
	else if(oiltype == 9)
		lubricant = new Parker_skydrol_oil(in);	// ParkerProject
	else if(oiltype == 50)
		lubricant = new water_oil(in);			// WATER
	else if(oiltype == 10)
		lubricant = new MILPRF87257_oil(in);		//MILPRF87257 aerocontrolex

	else
	{
		Log << "Oil type " << oiltype << " not supported" << gaplog::endl;
		exit(1);
	}
	Log << "done!\n";

	// initialize variables
	revolution = 0;
	phi = 0;
	t = 0;
	dt = stpang/speed; // dphi/omega

	// read the pressure file
	pfile.initialize(pfile_name, speed);

	// build the 2D mesh and the gap surfaces
	Log << "\nBuilding the 2D gap mesh and the gap surfaces... \n\n";
	mesh.reynolds.build(in, 2);
	mesh.reynolds.define_gap_surface(CB);
	mesh.reynolds.define_gap_surface(VP);
	Log <<"\ndone!" <<gaplog::endl;
	
	// build the 3D mesh
	Log << "\nBuilding the 3D gap mesh ... \n\n";
	mesh.energy.build(in, 3);
	Log <<"\ndone!" <<gaplog::endl;

	// initialize fields
	fields.p.initialize(&mesh.reynolds, 1e5);
	fields.p_prev.initialize(&mesh.reynolds, 1e5);
	fields.p_sol.initialize(&mesh.reynolds, 1e5);
	fields.p_EHD.initialize(&mesh.reynolds, 1e5);
	fields.mu.initialize(&mesh.energy, lubricant->get_mu(1e5, 20));
	fields.mu2D.initialize(&mesh.reynolds, lubricant->get_mu(1e5, 20));
	fields.fw.initialize(&mesh.energy, 0);
	fields.fe.initialize(&mesh.energy, 0);
	fields.fs.initialize(&mesh.energy, 0);
	fields.fn.initialize(&mesh.energy, 0);
	fields.pw.initialize(&mesh.energy, 0);
	fields.pe.initialize(&mesh.energy, 0);
	fields.ps.initialize(&mesh.energy, 0);
	fields.pn.initialize(&mesh.energy, 0);
	fields.vw.initialize(&mesh.energy, 0);
	fields.ve.initialize(&mesh.energy, 0);
	fields.vs.initialize(&mesh.energy, 0);
	fields.vn.initialize(&mesh.energy, 0);
	fields.g_pw.initialize(&mesh.energy, 0);
	fields.g_pe.initialize(&mesh.energy, 0);
	fields.g_ps.initialize(&mesh.energy, 0);
	fields.g_pn.initialize(&mesh.energy, 0);
	fields.mdotw.initialize(&mesh.energy, 0);
	fields.mdote.initialize(&mesh.energy, 0);
	fields.mdots.initialize(&mesh.energy, 0);
	fields.mdotn.initialize(&mesh.energy, 0);
	fields.rhow.initialize(&mesh.energy, lubricant->get_rho(1e5, 20));
	fields.rhoe.initialize(&mesh.energy, lubricant->get_rho(1e5, 20));
	fields.rhos.initialize(&mesh.energy, lubricant->get_rho(1e5, 20));
	fields.rhon.initialize(&mesh.energy, lubricant->get_rho(1e5, 20));
	fields.rho.initialize(&mesh.energy, lubricant->get_rho(1e5, 20));
	fields.rho2D.initialize(&mesh.reynolds, lubricant->get_rho(1e5, 20));
	fields.T.initialize(&mesh.energy, Tcase);
	fields.T_prev.initialize(&mesh.energy, Tcase);
	fields.phid.initialize(&mesh.energy, 0);
	fields.phid_avg.initialize(&mesh.energy, 0);
	fields.Tcb.initialize(&mesh.reynolds, Tcase);
	fields.Tvp.initialize(&mesh.reynolds, Tcase);
	fields.qcb.initialize(&mesh.reynolds, 0);
	fields.qcb_avg.initialize(&mesh.reynolds, 0);
	fields.qcb_prog.initialize(&mesh.reynolds, 0);
	fields.qvp.initialize(&mesh.reynolds, 0);
	fields.qvp_avg.initialize(&mesh.reynolds, 0);
	fields.qvp_prog.initialize(&mesh.reynolds, 0);
	fields.Vc.initialize(&mesh.energy, cy_vector(0,0,0));
	fields.Vp.initialize(&mesh.energy, cy_vector(0,0,0));
	fields.V.initialize(&mesh.energy, cy_vector(0,0,0));	
	
	// initialize film
	film.h1 = hB1, film.h2 = hB2, film.h3 = hB3;
	film.h.initialize(&mesh.reynolds, 10e-6);
	film.h_previous.initialize(&mesh.reynolds, 10e-6);
	film.hcb.initialize(&mesh.reynolds, 10e-6);
	film.hvp.initialize(&mesh.reynolds, 0);
	film.hcb_rigid.initialize(&mesh.reynolds, 10e-6);
	film.hvp_rigid.initialize(&mesh.reynolds, 0);
	film.dhcb_thermal.initialize(&mesh.reynolds, 0);
	film.dhvp_thermal.initialize(&mesh.reynolds, 0);
	film.dhcb_EHD.initialize(&mesh.reynolds, 0);
	film.dhvp_EHD.initialize(&mesh.reynolds, 0);
	film.contact.initialize(&mesh.reynolds, 0);

	// initialize squeeze
	squeeze.dhdt1 = 0, squeeze.dhdt2 = 0, squeeze.dhdt3 = 0;
	squeeze.d_dhdt1 = 0, squeeze.d_dhdt2 = 0, squeeze.d_dhdt3 = 0;
	squeeze.dhdt.initialize(&mesh.reynolds, 0);
	squeeze.dhdt_total.initialize(&mesh.reynolds, 0);
	squeeze.dhdt_hs.initialize(&mesh.reynolds, 0);
	squeeze.dhdt_hd.initialize(&mesh.reynolds, 0);

	// initialize forces
	forces.FTG.initialize("./output/slipper/ftg.txt");
	forces.FTK.initialize("./output/piston/ftk.txt");
	forces.DC.initialize(DC_mesh.c_str());
	if(!forces.DC.is_available() && in.data.geometry.ADC == 0)
	{
		Log << "\nDC_mesh file not found and ADC seems to be 0.\n"
				<< "Please check the geometry.txt file.\n\n";
		exit(1);
	}

	forces.p_cb.resize(mesh.reynolds.cb_gap_surf.elms_0.size(), 0);
	forces.pistonforces.resize(z, std::vector<double>(4));
	forces.FfBz = 0, forces.MfBx = 0, forces.MfBy = 0;
	forces.Ffluid[0] = 0, forces.Ffluid[1] = 0, forces.Ffluid[2] = 0;
	forces.FBz = 0, forces.MBx = 0, forces.MBy = 0;
	forces.Fext[0] = 0, forces.Fext[1] = 0, forces.Fext[2] = 0;
	forces.dF[0] = 0, forces.dF[1] = 0, forces.dF[2] = 0;
	forces.Fc = 0, forces.Mcx = 0, forces.Mcy = 0;
	forces.MTBz = 0;
	forces.FKBx = 0, forces.FKBy = 0;
	forces.MKBx = 0, forces.MKBy = 0;
	forces.MBx_spline = 0;

	// initialize the piston structure
	Log << "\nInitializing pistons structure ... ";
	initialize_pistons();
	Log << "done!\n";

	// initialize the fluid interpolation grids
	Log << "\nInitializing the fluid/structure interpolation grids ... ";
	CBf.initialize(mesh.reynolds.mn, mesh.reynolds.mn + mesh.reynolds.n, 4);
	VPf.initialize(mesh.reynolds.mn, mesh.reynolds.mn + mesh.reynolds.n, 4);
	define_interpl_grids();
	Log << "done!\n";

	// macrogeometry initialization

	// block
	macro_CB = 0;
	
	if(macro_CB_opt == 1)		// axisymmetric
	{
		macro_CB = new macro_axisymmetric(film.hcb_rigid, "./input/block/macrogeometry_block.txt");
	}
	else if(macro_CB_opt == 2)	// waved
	{
		macro_CB = new macro_waved(film.hcb_rigid, "./input/block/macrogeometry_block.txt");
	}
	else if(macro_CB_opt == 3)	// spherical
	{
		macro_CB = new macro_spherical(film.hcb_rigid, "./input/block/macrogeometry_block.txt");
	}
	else if(macro_CB_opt == 4)	// from grid
	{
		macro_CB = new macro_grid(film.hcb_rigid, "./input/block/macrogeometry_block.txt");
	}
	else if(macro_CB_opt == 5)	// notches defined through stl file
	{
		macro_CB = new macro_stl(film.hcb_rigid, "./input/block/macrogeometry_block.txt");
	}
	else if(macro_CB_opt == 6)	// defined using a fluid grid file
	{
		macro_CB = new macro_fluid_grid(film.hcb_rigid, "./input/block/macrogeometry_block.txt");
	}

	// valve plate
	macro_VP = 0;
	
	if(macro_VP_opt == 1)	// axisymmetric
	{
		macro_VP = new macro_axisymmetric(film.hvp_rigid, "./input/block/macrogeometry_valveplate.txt");
	}
	else if(macro_VP_opt == 2)	// waved
	{
		macro_VP = new macro_waved(film.hvp_rigid, "./input/block/macrogeometry_valveplate.txt");
	}
	else if(macro_VP_opt == 3)	// spherical
	{
		macro_VP = new macro_spherical(film.hvp_rigid, "./input/block/macrogeometry_valveplate.txt");
	}
	else if(macro_VP_opt == 4)	// from grid
	{
		macro_VP = new macro_grid(film.hvp_rigid, "./input/block/macrogeometry_valveplate.txt");
	}
	else if(macro_VP_opt == 5)	// notches defined through stl file
	{
		macro_VP = new macro_stl(film.hvp_rigid, "./input/block/macrogeometry_valveplate.txt");
	}
	else if(macro_VP_opt == 6)	// defined using a fluid grid file
	{
		macro_VP = new macro_fluid_grid(film.hvp_rigid, "./input/block/macrogeometry_valveplate.txt");
	}
	
	
}
// ------------------------------------------------------------------------- //
void cblock_gap_main::update_pistons(double phi_rad)
{
	
	// ------------------------ copy from caspar_input ----------------------- //

	int z = in.data.operating_conditions.npistons;						// number of pistons
	double omega = in.data.operating_conditions.speed;				// speed
	double beta = in.data.operating_conditions.beta;					// swash plate angle [rad]
	double betamax = in.data.operating_conditions.betamax;		// maximum swash plate angle [rad]
	double lengthB = in.data.geometry.lengthB;								// lenght of the block (no neck!)
	double lengthcanalB = in.data.geometry.lengthcanalB;			// lenght of the kidney canal
	double lF = in.data.geometry.lF;													// bushing length
	double rB = 0.50*in.data.geometry.dB;											// pitch radius of pistons in the cylinder bore
	double lK = in.data.geometry.lK;													// piston length (from ball center to bottom, no chamfer!)
	double lch = in.data.geometry.lch;												// lenght of piston chamfer
	double lKG = in.data.geometry.lKG;												// length of piston gap surface
	double le = in.data.geometry.le;													// bushing beginning position (from top of the block)
	double gamma = in.data.geometry.gamma;										// cross angle
	//double pos_A = 0.0;			//Position of the front flat piston edge from the bottom of the cylinder (displacement chamber) [m] 
	double pos_B = 0.0;			//Position of the back flat piston edge from beginning of bushing [m] 
	double lvar = 0.0;			//Variable gap length of the piston-cylinder interface[m] 
	// ----------------------------------------------------------------------- //

		
	// ------------------- update new pistons information ------------------- //

	// first piston -> reference angle
	pistons.phi[0] = phi_rad;	
	// update pressure in DCs first member, pressure @ reference angle
	pistons.pDC[0] = pfile.getp_rad(phi_rad)[0]; 
	pistons.pLP = pfile.getp_rad(phi_rad)[1];
	pistons.pHP = pfile.getp_rad(phi_rad)[2];

	// update all the others
	for(int i=1; i < z; i++) 
	{
		pistons.phi[i] = phi_rad + ((2.0*pi)/z)*static_cast<double>(i);
		if(pistons.phi[i] >= 2.0*pi)
			pistons.phi[i] -= 2.0*pi;
		pistons.pDC[i] = pfile.getp_rad(pistons.phi[i])[0];
	}

	// distance from front flat piston surface to bushing end when piston is at ODC
	double delta_A = lengthB - lengthcanalB - pistons.lZ0 - lF - le - pistons.T;
	// distance from piston head center to cylinder block beginning edge 
	// when piston is at ODC
	double delta_B = lK - delta_A - lF - le;

	// loop over all pistons
	for(int i = 0; i < z; i++) 
	{
		
		// --------------------- update pistons stroke ------------------------- //
		
		// angle
		double phi = pistons.phi[i];

		// piston stroke (negative)
		double sK = -rB*tan(beta)*(1-cos(phi)) - rB*tan(gamma)*sin(phi)/cos(beta) + pistons.delta_psi;
		pistons.sK[i] = sK - 0.5*(pistons.Hk_betamax - pistons.Hk);

		// update piston position
		pistons.zK[i] = 0.5*pistons.Hk_betamax + pistons.sK[i];

		// update pistons velocity
		pistons.vK[i] = -omega*rB*(sin(phi)*tan(beta) + tan(gamma)*cos(phi)/cos(beta));

		// update pistons acceleration
		pistons.aK[i] = -omega*omega*rB*(tan(beta)*cos(phi) + tan(gamma)*sin(phi)/cos(beta));



		//Position of the front flat piston edge from the bottom of the cylinder
		pistons.pos_A = pistons.zK[i]+pistons.lZ0 + pistons.T;
		//Distance from piston head center to cylinder block beginning edge when piston is at ODC [m]
		pos_B = pistons.zK[i] + delta_B + le - (lK - lKG);
	
		if (pos_B < 0.0)
		{
			pos_B = 0.0;
		}
		pistons.pos_B = pos_B;
	
		//Variable gap length [m]
		lvar = lengthB - lengthcanalB - pistons.pos_A - pos_B - le ;      
		if ( lvar >= lF - pos_B )
		{
			lvar = lF  - pos_B ;
		}
		pistons.lvar[i] = lvar;

	}

}
// ------------------------------------------------------------------------- //
void cblock_gap_main::rotate(double phi_rad)
{
	// check the angle
	if(phi_rad < 0 || phi_rad > 2.0*pi)
	{
		Log << "\ncblock_gap_main::rotate: angle must be within [0,2pi[" 
				<< gaplog::endl;
		exit(1);
	}

	// update the angle
	phi = phi_rad;
	// update the pistons structure
	update_pistons(phi_rad);
	// rotate the mesh
	mesh.reynolds.rotate(phi_rad);
	mesh.energy.rotate(phi_rad);

}
// ------------------------------------------------------------------------- //
void cblock_gap_main::set_pressure()
{

	// ------------------------ copy from caspar_input ----------------------- //
	
	int z = in.data.operating_conditions.npistons;
	double pcase = in.data.operating_conditions.pCase;

	// ----------------------------------------------------------------------- //
	
	// loop to all the elements
	for(int i = 0; i < mesh.reynolds.mn; i++) 
	{
		
		bool defined = false;

		// displacement chamber
		for (int j=1; j<=z; j++) 
		{
			int ty = mesh.reynolds.elements[i].ty;
			if(ty == j || ty - LP == j || ty - HP == j) 
			{
				// set to DC pressure
				fields.p.in[i] = pfile.getp_rad(pistons.phi[j-1])[0]; 
				defined = true;
			}
		}
		// low pressure
		if(mesh.reynolds.elements[i].ty == LP) 
		{
			fields.p.in[i] = pfile.getp_rad(pistons.phi[0])[1]; 
			defined = true;
		}
		// high pressure
		else if (mesh.reynolds.elements[i].ty == HP)
		{
			fields.p.in[i] = pfile.getp_rad(pistons.phi[0])[2]; 
			defined = true;
		}
		// case
		else if (mesh.reynolds.elements[i].ty >= CASE)
		{
			fields.p.in[i] = pcase; 
			defined = true;
		}
		else if(mesh.reynolds.elements[i].ty == FLUID) 
		{
			// set to reference value
			//fields.p.in[i] = 1e5;
			//fields.p.in[i] = 1e5;
			defined = true;
		}
		
		if(defined = false)
		{
			Log << "\ncblock_gap_main::setPressure() Error! "
					 << "Undefined gap_elm type! index = " << i 
					 << " type index = " << mesh.reynolds.elements[i].ty << gaplog::endl;
			exit(1);
		}
	}

	for(int j=0; j<mesh.reynolds.n; j++) 
	{
		// set to fields.p case;
		fields.p.bound.inner[j] = pcase;
		fields.p.bound.outer[j] = pcase;
	}
}
// ------------------------------------------------------------------------- //
void cblock_gap_main::set_temperature()
{
	
	// ------------------------ copy from caspar_input ----------------------- //
	
	int z = in.data.operating_conditions.npistons;
	double TLP = in.data.operating_conditions.T_LP;
	double THP = in.data.operating_conditions.T_HP;
	double Tcase = in.data.operating_conditions.T_Leak;

	// ----------------------------------------------------------------------- //
	
	// rotate the Tcb field according with the blockposition
	scalar_field Tcb = fields.Tcb.cshift(mesh.reynolds.rotation_steps);

	double pLP = pfile.getp_rad(pistons.phi[0])[1];
	double pHP = pfile.getp_rad(pistons.phi[0])[2];

	for(int k = 0; k < mesh.energy.q; k++) 
	{
		for(int i = 0; i < mesh.energy.mn; i++) 
		{
			
			bool defined = false;

			for(int j=1; j<=z; j++) 
			{
				if (mesh.reynolds.elements[i].ty == j) 
				{
					double pDCi = fabs(pfile.getp_rad(pistons.phi[j - 1])[0]);
					// set to DC temperature
					if(fabs(pDCi - pHP) < 0.1*pHP)
						fields.T.in[k*(mesh.energy.mn) + i] = THP; 
					else if(fabs(pi - pLP) < 0.1*pLP)
						fields.T.in[k*(mesh.energy.mn) + i] = TLP; 
					else
					{
						double alpha = (pi/pHP);
						fields.T.in[k*(mesh.energy.mn) + i] = alpha*THP + (1 - alpha)*TLP; 
					}

					defined = true;
				}
			}

			if(mesh.reynolds.elements[i].ty >= LP && mesh.reynolds.elements[i].ty < HP) 
			{
				fields.T.in[k*(mesh.energy.mn) + i] = TLP; // set to LP T
				defined = true;
			}

			else if(mesh.reynolds.elements[i].ty >= HP && mesh.reynolds.elements[i].ty < CASE)
			{
				fields.T.in[k*(mesh.energy.mn) + i] = THP; // set to HP T
				defined = true;
			}
			else if(mesh.reynolds.elements[i].ty == CASE)
			{
				fields.T.in[k*(mesh.energy.mn) + i] = Tcase; // set to T case
				defined = true;
			}
			else if(mesh.reynolds.elements[i].ty == FLUID) 
			{
				fields.T.in[k*(mesh.energy.mn) + i] = 
					0.5*(Tcb.in[i] + fields.Tvp.in[i]); // set to reference value
				defined = true;
			}
			if(defined = false)
			{
				Log << "\ncblock_gap_main::setTemperature() Warning!"
						 << "Undefined gap_elm type! index = " << i 
					 << " type index = " << mesh.reynolds.elements[i].ty << gaplog::endl;
			}
		}
	}

	// top and bottom boundaries
	// this field accounts for the block rotation!
	for(int i = 0; i < mesh.energy.mn; i++) 
	{
		if(mesh.reynolds.cb_elms[i] == FLUID)
			fields.T.bound.top[i] = Tcb.in[i];
		else
			fields.T.bound.top[i] = 0.5*(THP + TLP);
		
		if(mesh.energy.vp_elms[i] == FLUID)
			fields.T.bound.bottom[i] = fields.Tvp.in[i];
		else
			fields.T.bound.bottom[i] = 0.5*(THP + TLP);
	}

	// inner and outer boundaries
	for(int i=0; i<mesh.energy.n*mesh.energy.q; i++)
	{
		fields.T.bound.inner[i] = Tcase;
		fields.T.bound.outer[i] = Tcase;
	}
}
// ------------------------------------------------------------------------- //
scalar_field cblock_gap_main::calc_cb_heatflux()
{

	// ------------------------ copy from caspar_input ----------------------- //

	double Tcase = in.data.operating_conditions.T_Leak;
	double lambda = in.data.oil.general.oillambda;
	double q_min_limit = in.data.options_block.numeric.q_min_limit;
	double q_max_limit = in.data.options_block.numeric.q_max_limit;

	// ----------------------------------------------------------------------- //

	scalar_field q_cb(&mesh.reynolds);
	q_cb = 0.0;
	
	int n = mesh.energy.n;
	int m = mesh.energy.m;
	int q = mesh.energy.q;
	int mn = mesh.energy.mn;

	// ------------------------------ dissipation ---------------------------- //
	
	# ifdef DISSIPATION

	for(int id2D=0; id2D<mn; id2D++)
	{
		// exclude the valve plate ports
		if(mesh.reynolds.elements[id2D].ty != FLUID)
		{
			q_cb.in[id2D] = 0;
		}
		// fluid domain
		else
		{
			// full film lubrication
			if(film.h.in[id2D] > 0.1e-6)
			{
				
				double q_conv = 0;
				double q_gen = 0;
				
				for(int k=0; k<mesh.energy.q; k++)
				{
					int id3D = mn*k + id2D;
					gap_elm e = mesh.energy.elements[id3D];
					double cp = lubricant->get_cp(fields.T.in[id3D]);
					
					// heat generated
					q_gen += fields.phid.in[id3D];

					// heat transfert by convection
					
					// WEST
					if(fields.fw.in[id3D] > 0) // flow in (positive fw)
						q_conv -= fabs(fields.fw.in[id3D])*cp*fields.T.in[e.w];
					else	// flow out (negative fw)
						q_conv += fabs(fields.fw.in[id3D])*cp*fields.T.in[id3D];
					
					// EAST
					if(fields.fe.in[id3D] < 0) // flow in (negative fe)
						q_conv -= fabs(fields.fe.in[id3D])*cp*fields.T.in[e.e];
					else	// flow out (positive fe)
						q_conv += fabs(fields.fe.in[id3D])*cp*fields.T.in[id3D];

					// SOUTH
					if(fields.fs.in[id3D] > 0 && e.s > -1) // flow in (positive fs)
						q_conv -= fabs(fields.fs.in[id3D])*cp*fields.T.in[e.s];
					else if(fields.fs.in[id3D] > 0 && e.s == -1)
						q_conv -= fabs(fields.fs.in[id3D])*cp*Tcase;
					else	// flow out (negative fs)
						q_conv += fabs(fields.fs.in[id3D])*cp*fields.T.in[id3D];

					// NORTH
					if(fields.fn.in[id3D] < 0 && e.n > -1) // flow in (negative fn)
						q_conv -= fabs(fields.fn.in[id3D])*cp*fields.T.in[e.n];
					else if(fields.fn.in[id3D] < 0 && e.n == -1)
						q_conv -= fabs(fields.fn.in[id3D])*cp*Tcase;
					else	// flow out (positive fn)
						q_conv += fabs(fields.fn.in[id3D])*cp*fields.T.in[id3D];
				}

				q_cb.in[id2D] = (q_gen - q_conv)/(2.0*mesh.reynolds.elements[id2D].A);
			}
			// who the hell knows what happens here in terms of dissipation ...
			else
			{
				double omega = in.data.operating_conditions.speed;
				double U = omega*mesh.reynolds.elements[id2D].r;
				double f = 0.15;
				q_cb.in[id2D] = f*fields.p.in[id2D]*U;
			}
		}
	}	
	
	# endif

	// ------------------------------- gradient ------------------------------ //

	# ifndef DISSIPATION

	for(int i=0; i<m; i++)
	{
		for(int j=0; j<n; j++)
		{
			int id = n*i + j;
			// exclude the valve plate ports
			if(mesh.reynolds.elements[id].ty != FLUID)
			{
				q_cb.in[id] = 0;
			}
			// fluid domain
			else
			{
				q_cb.in[id] = -1.0*lambda*fields.T.get_zgrad_top(i,j);
			}
		}
	}

	# endif

	// ----------------------------------------------------------------------- //

	// set the boundary
	for(int j = 0; j < n; j++)
	{
		q_cb.bound.inner[j] = q_cb.in[j];
		q_cb.bound.outer[j] = q_cb.in[n*(m-1) + j];
	}

	// shift to the zero position
	q_cb = q_cb.cshift(-mesh.reynolds.rotation_steps);


	return q_cb;

}
// ------------------------------------------------------------------------- //
scalar_field cblock_gap_main::calc_vp_heatflux()
{
	// ------------------------ copy from caspar_input ----------------------- //

	double Tcase = in.data.operating_conditions.T_Leak;
	double lambda = in.data.oil.general.oillambda;
	double q_min_limit = in.data.options_block.numeric.q_min_limit;
	double q_max_limit = in.data.options_block.numeric.q_max_limit;

	// ----------------------------------------------------------------------- //

	scalar_field q_vp(&mesh.reynolds);
	q_vp = 0.0;
	
	int n = mesh.energy.n;
	int m = mesh.energy.m;
	int q = mesh.energy.q;
	int mn = mesh.energy.mn;

	// ------------------------------ dissipation ---------------------------- //

	# ifdef DISSIPATION

	for(int id2D=0; id2D<mn; id2D++)
	{
		// exclude the valve plate ports
		if(mesh.reynolds.elements[id2D].ty != FLUID)
		{
			q_vp.in[id2D] = 0;
		}
		// fluid domain
		else
		{
			// full film lubrication
			if(film.h.in[id2D] > 0.1e-6)
			{
				double q_conv = 0;
				double q_gen = 0;
				
				for(int k=0; k<mesh.energy.q; k++)
				{
					int id3D = mn*k + id2D;
					gap_elm e = mesh.energy.elements[id3D];
					double cp = lubricant->get_cp(fields.T.in[id3D]);
					
					// heat generated
					q_gen += fields.phid.in[id3D];

					// heat transfert by convection
					
					// WEST
					if(fields.fw.in[id3D] > 0) // flow in (positive fw)
						q_conv -= fabs(fields.fw.in[id3D])*cp*fields.T.in[e.w];
					else	// flow out (negative fw)
						q_conv += fabs(fields.fw.in[id3D])*cp*fields.T.in[id3D];
					
					// EAST
					if(fields.fe.in[id3D] < 0) // flow in (negative fe)
						q_conv -= fabs(fields.fe.in[id3D])*cp*fields.T.in[e.e];
					else	// flow out (positive fe)
						q_conv += fabs(fields.fe.in[id3D])*cp*fields.T.in[id3D];

					// SOUTH
					if(fields.fs.in[id3D] > 0 && e.s > -1) // flow in (positive fs)
						q_conv -= fabs(fields.fs.in[id3D])*cp*fields.T.in[e.s];
					else if(fields.fs.in[id3D] > 0 && e.s == -1)
						q_conv -= fabs(fields.fs.in[id3D])*cp*Tcase;
					else	// flow out (negative fs)
						q_conv += fabs(fields.fs.in[id3D])*cp*fields.T.in[id3D];

					// NORTH
					if(fields.fn.in[id3D] < 0 && e.n > -1) // flow in (negative fn)
						q_conv -= fabs(fields.fn.in[id3D])*cp*fields.T.in[e.n];
					else if(fields.fn.in[id3D] < 0 && e.n == -1)
						q_conv -= fabs(fields.fn.in[id3D])*cp*Tcase;
					else	// flow out (positive fn)
						q_conv += fabs(fields.fn.in[id3D])*cp*fields.T.in[id3D];
				}

				q_vp.in[id2D] = (q_gen - q_conv)/(2.0*mesh.reynolds.elements[id2D].A);	
			}
			else
			{
				double omega = in.data.operating_conditions.speed;
				double U = omega*mesh.reynolds.elements[id2D].r;
				double f = 0.15;
				q_vp.in[id2D] = f*fields.p.in[id2D]*U;
			}
		}
	}			

	# endif

	// ------------------------------- gradient ------------------------------ //

	# ifndef DISSIPATION

	for(int i=0; i<m; i++)
	{
		for(int j=0; j<n; j++)
		{
			int id = n*i + j;
			// exclude the valve plate ports
			if(mesh.reynolds.elements[id].ty != FLUID)
			{
				q_vp.in[id] = 0;
			}
			// fluid domain
			else
			{
				q_vp.in[id] = -1.0*lambda*fields.T.get_zgrad_top(i,j);
			}
		}
	}

	# endif

	// ----------------------------------------------------------------------- //

	// set the boundary
	for(int j = 0; j < n; j++)
	{
		q_vp.bound.inner[j] = q_vp.in[j];
		q_vp.bound.outer[j] = q_vp.in[n*(m-1) + j];
	}

	return q_vp;
}
// ------------------------------------------------------------------------- //
double cblock_gap_main::get_qcb_avg()
{

	double q_avg = 0;
	double A = 0;
	for(int id2D=0; id2D<mesh.reynolds.mn; id2D++)
	{
		// exclude the block's openings
		if(mesh.reynolds.cb_elms_0[id2D] == FLUID)
		{
			q_avg += fields.qcb_prog.in[id2D]*mesh.reynolds.elements[id2D].A;
			A += mesh.reynolds.elements[id2D].A;
		}
	}

	return q_avg/A;
}
// ------------------------------------------------------------------------- //
double cblock_gap_main::get_qvp_avg()
{

	double q_avg = 0;
	double A = 0;
	for(int id2D=0; id2D<mesh.reynolds.mn; id2D++)
	{
		// exclude the valve plate ports
		if(mesh.reynolds.vp_elms[id2D] == FLUID)
		{
			q_avg += fields.qvp_prog.in[id2D]*mesh.reynolds.elements[id2D].A;
			A += mesh.reynolds.elements[id2D].A;
		}
	}

	return q_avg/A;
}
// ------------------------------------------------------------------------- //
void cblock_gap_main::apply_heat_flux(interpl_grid* structure, const char* body)
{
	scalar_field* q_prog = 0;
	interpl_grid* fluid_grid = 0;
	
	if(string(body).compare("CB") == 0)
	{
		q_prog = &fields.qcb_prog;
		fluid_grid = &CBf;
	}
	else
	{
		q_prog = &fields.qvp_prog;
		fluid_grid = &VPf;
	}

	// fill the fluid grid interpl_grid cells_data with q_prog
	for(unsigned int i=0; i<fluid_grid->ne(); i++)
	{
		if(fluid_grid->active[i])
			fluid_grid->cells_data[i] = q_prog->in[i];
		else
			fluid_grid->cells_data[i] = 0;
	}

	// interpolate q from fluid to structure grid
	interpolation fluid_to_structure;
	fluid_to_structure.cellsTocells(fluid_grid, structure);	

}
	// ------------------------------------------------------------------------- //
void cblock_gap_main::calc_fluid_forces()
{
	
	// ------------------------ copy from caspar_input ----------------------- //

	double R = 0.5*in.data.geometry.d_spherical;	

	// ----------------------------------------------------------------------- //

	// ----- interpolate pressure from mesh.reynolds to mesh.cb_gap_surf ----- //
	
	// ann stuff
	ANNpointArray dataPts;		// data points
	ANNpoint queryPt;					// query point
	ANNidxArray nnIdx;				// near neighbor indices
	ANNdistArray dists;				// near neighbor distances
	ANNkd_tree* kdTree;				// search structure

	// define dataPts using cb_gap_surf.elms center
	dataPts = annAllocPts(mesh.reynolds.mn, 2);
	for(int i=0; i<mesh.reynolds.mn; i++)
	{
		dataPts[i][0] = mesh.reynolds.elements[i].r*sin(mesh.reynolds.elements[i].theta);
		dataPts[i][1] = mesh.reynolds.elements[i].r*cos(mesh.reynolds.elements[i].theta);
	}
	// build search structure
	kdTree = new ANNkd_tree(dataPts, mesh.reynolds.mn, 2);

	// allocate query point
	queryPt = annAllocPt(2);
	nnIdx = new ANNidx[1];
	dists = new ANNdist[1];

	for(unsigned int i=0; i<mesh.reynolds.cb_gap_surf.elms.size(); i++)
	{
		point c = mesh.reynolds.cb_gap_surf.elms[i].center();
		queryPt[0] = c.x(), queryPt[1] = c.y();
		kdTree->annkSearch(queryPt, 1, nnIdx, dists);	
		forces.p_cb[i] = fields.p.in[nnIdx[0]];
	}

	// free memory
	delete [] nnIdx;
	delete [] dists;
	annDeallocPts(dataPts);
	annDeallocPt(queryPt);
	delete kdTree;
	annClose();

	// ----------------------------------------------------------------------- //

	// distance between cylinder block reference system and intersection point
	// between spherical surface and block axis
	double z0 = pistons.z0;
	
	// additional moments in x and y direction to the the spherical gap
	double MfBx_add = 0;
	double MfBy_add = 0;

	// clean
	forces.FfBz = 0;
	forces.MfBx = 0;
	forces.MfBy = 0;

	double A = 0;	// this is going to be the total cylinder block gap surface
	
	// loop through all the polygons defining the cb gap surface
	for(unsigned int i=0; i<mesh.reynolds.cb_gap_surf.elms.size(); i++)
	{
		double da = mesh.reynolds.cb_gap_surf.ai[i];
		double dF = forces.p_cb[i]*da;
		A += da;

		point c = mesh.reynolds.cb_gap_surf.elms[i].center();

		// ------------------------- spherical gap --------------------------- //
		if(R > 0)
		{
			// element center
			double xc = c.x();
			double yc = c.y();
			double rc = sqrt(xc*xc + yc*yc);
			
			// angle between element position on the spherical surface and z axis
			double alpha = asin(rc/R);

			// angle in the gap plane (projection of the spherical surface in the x-y plane)
			double phi = atan2(xc, yc);
			phi = (phi >= 0) ? phi : 2.0*pi + phi;

			// pressure force components
			double dFz = dF;
			double dFxy = dFz*tan(alpha);
			double dFx = dFxy*sin(phi);
			double dFy = dFxy*cos(phi);
				
			// axial force
			forces.FfBz += dFz;
								
			// additional length in axial direction to be added to z0 which 
			// accounts for position of the cell on the spherical surface
			double delta_z = R - sqrt(pow(R, 2.0) - pow(rc, 2.0));

			// additional moments due to the spherical gap
			MfBx_add = dFy*(z0 + delta_z);
			MfBy_add = dFx*(z0 + delta_z);

			// add to the original moments
			forces.MfBx += dFz*yc + MfBx_add;
			forces.MfBy += -1.0*(dFz*xc + MfBy_add); // a positive x yields to a negative momentum in y direction

		}
		// ---------------------------- flat gap ----------------------------- //
		else
		{
			forces.FfBz += dF;
			forces.MfBx += dF*c.y();
			forces.MfBy += -1.0*dF*c.x(); // a positive x yields to a negative momentum in y direction
		}
	}

	//cout << "Total area: " << A << endl;
	//cout << "FfBz " << forces.FfBz << "[N]" << endl;
	//cout << "MfBx " << forces.MfBx << "[N]" << endl;
	//cout << "MfBy " << forces.MfBy << "[N]" << endl;

}
// ------------------------------------------------------------------------- //
void cblock_gap_main::calc_dhdt_total()
{
	
	if (film.h_previous.in[0] > 0)
	{
		for(int i = 0; i <  mesh.reynolds.mn; i++) 
		{
			squeeze.dhdt_total.in[i] = (film.h_previous.in[i] - film.h.in[i])/dt;
			fields.p_prev.in[i]  = fields.p.in[i];
		}
	}
	else
	{
		for(int i = 0; i <  mesh.reynolds.mn; i++) 
		{
			squeeze.dhdt_total.in[i] = squeeze.dhdt.in[i] + squeeze.dhdt_hd.in[i] + squeeze.dhdt_hs.in[i] ;
			fields.p_prev.in[i]  = fields.p.in[i];
		}
	}
}
// ------------------------------------------------------------------------- //
void cblock_gap_main::calc_piston_forces()
{
	// ------------------------ copy from caspar_input ----------------------- //

	int z = in.data.operating_conditions.npistons;
	double omega = in.data.operating_conditions.speed;
	double beta = in.data.operating_conditions.beta;
	double rB = 0.50*in.data.geometry.dB;
	double rK = 0.50*in.data.geometry.dK;
	double rDK = 0.50*in.data.geometry.dDK;
	double lSK = in.data.geometry.lSK;
	double mK = in.data.geometry.mK;
	double pCase = in.data.operating_conditions.pCase;
	double gamma = in.data.geometry.gamma;

	// ----------------------------------------------------------------------- //

	// derived input

	// rotation matrix around x axis
	matrix Rx(3,3);
	Rx[0][0] = 1, Rx[0][1] = 0, Rx[0][2] = 0;
	Rx[1][0] = 0, Rx[1][1] = cos(beta), Rx[1][2] = -sin(beta);
	Rx[2][0] = 0, Rx[2][1] = sin(beta), Rx[2][2] = cos(beta);

	// rotation matrix around y axis
	matrix Ry(3,3);
	Ry[0][0] = cos(gamma), Ry[0][1] = 0, Ry[0][2] = sin(gamma);
	Ry[1][0] = 0, Ry[1][1] = 1, Ry[1][2] = 0;
	Ry[2][0] = -sin(gamma), Ry[2][1] = 0, Ry[2][2] = cos(gamma);

	// x y and z versors
	matrix nx(3,1), ny(3,1), nz(3,1);
	nx[0][0] = 1, nx[1][0] = 0, nx[2][0] = 0;
	ny[0][0] = 0, ny[1][0] = 1, ny[2][0] = 0;
	nz[0][0] = 0, nz[1][0] = 0, nz[2][0] = 1;

	// normal to the swash plate
	matrix n = Rx*(Ry*nz);

	// ----------------------------------------------------------------------- //

	// piston front face area
	double AreaK = pi*(pow(rK,2.0) - pow(rDK,2.0));

	// update friction forces files
	forces.FTG.update(2.0*pi/omega);
	forces.FTK.update(2.0*pi/omega);
	
	// update moments on block
	forces.MKBx = forces.MKBy = 0;
	forces.FKBx = forces.FKBy = 0;
	
	// total force on block, piston reference system
	//double FkAx = 0, FkAy = 0;	// point A
	//double FkBx = 0, FkBy = 0;  // point B
	//double FkA = 0, FkB = 0;
	// total force on block, block reference system
	//double FkA_xB = 0, FkA_yB = 0;	// point A
	//double FkB_xB = 0, FkB_yB = 0;  // point B
	
	// point of application of the resultant radial force on the block
	//double zFkA = 0, zFkB = 0;

	// loop to all the z pistons
	for(int i=0; i < z; i++) 
	{

		double phii = pistons.phi[i];
		double zKi = pistons.zK[i];

		// pressure force
		double FDK = AreaK*(pistons.pDC[i] - pCase);
		
		// inertia force z directions
		double FaKz = -mK*pistons.aK[i];
		
		// centifiugal force
		double FwK = mK*(omega*omega)*rB;
		
		// slipper friction force tangential direction - add fucntion 
		double FTGx = -forces.FTG.getf(phii); 
		
		// slipper friction force y direction
		double FTGy = 0.0;	
		
		// forces on the piston acting in z (axial) direction
		double FTK = forces.FTK.getf(phii);
		
		if (FTK == 0)
		{
			double area;
			
			area = 2*pi*rK*pistons.lvar[i];

			double mu = lubricant->get_mu(in.data.operating_conditions.HP,in.data.operating_conditions.T_Leak);
			FTK = mu * pistons.vK[i] * area / pistons.radial_clearance;
		}

		double FAKz = FDK + FaKz + FTK;//forces.FTK.getf(phii);
		
		// force piston perpendicular to the swashplate
		// if the force points away from the swash plate, there is no reaction of the swashplate
		// (the retainer will stop the slipper from moving away)
		double FSK = (FAKz > 0) ? FAKz/fabs(n[2][0]) : 0;
		
		// transverse piston force x direction (this should be != 0 only when gamma != 0)
		double FSKx = FSK*fabs(n[0][0]);
		// positive gamma makes the normal to point in the x direction of the 
		// main reference system (which has x in the opposite direction of the block's reference system!)
		FSKx = (gamma > 0) ? FSKx : -1.0*FSKx;	

		// transverse piston force y direction
		double FSKy = FSK*fabs(n[1][0]); 

		// update MBx
		forces.MKBx -= ( FSKy*zKi + FwK*cos(phii)*(zKi - lSK) );
		forces.MKBy += ( FSKx*zKi + FwK*sin(phii)*(zKi - lSK) );

		// update the total piston forces in x and y directions
		forces.FKBx += FSKx + FwK*sin(phii);
		forces.FKBy += FSKy + FwK*cos(phii);
		
	}


}
// ------------------------------------------------------------------------- //
void cblock_gap_main::calc_external_forces()
{
	
	// ------------------------ copy from caspar_input ----------------------- //

	//int Fext_method = in.data.options_block.numeric.Fext_method;
	int z = in.data.operating_conditions.npistons;
	double rB = 0.50*in.data.geometry.dB;
	double beta = in.data.operating_conditions.beta;
	double lengthcanalB = in.data.geometry.lengthcanalB;
	double Fblock = in.data.geometry.Fblock;
	double ADC = in.data.geometry.ADC;
	std::string DC_mesh = in.data.geometry.DC_mesh;
	double delta_z0 = in.data.geometry.delta_z0;
	double lSK = in.data.geometry.lSK;
	double lK = in.data.geometry.lK;
	
	// the following are used just by the simplified method
	double rK = 0.50*in.data.geometry.dK;
	double rZ = 0.50*in.data.geometry.dZ;
	double mK = in.data.geometry.mK;
	double omega = in.data.operating_conditions.speed;
	const double HK = 2.0*rB*tan(beta);	// maximum piston stroke [m]
	const double AK = pi*pow(rK,2);			// piston section
	const double Ac = pi*pow(rZ,2);			// cylinder bore area
	
	// ----------------------------------------------------------------------- //

	// total loads
	double FBan = -1.0*Fblock;		// total axial force (initialize to the spring value)
	double MBx = 0;								// total x-moment
	double MBy = 0;								// total y-moment

	// update the friction force
	forces.FTK.update(2.0*pi/omega);
	forces.FTG.update(2.0*pi/omega);

	// --------------- forces trasmitted from piston to block ---------------- //

	// calc piston forces with advanced calculation method
	calc_piston_forces();

	// MKBx and MKBy are calculated in calc_piston_forces
	MBx += forces.MKBx;	
	MBy += forces.MKBy;

	// add the moments associated with the spline reaction if delta z0 is
	// not on the origin of the block's reference system
	// the equation is Mx_sp = -delta_z0(-FkBy):
	// -FkBy is the reaction of the shaft ON the block
	forces.MBx_spline = delta_z0*forces.FKBy;
	MBx += forces.MBx_spline;

	// --------------- forces on block due to pressure in DCs ---------------- //

	double F0z;								// resultant pressure force on the displacement chamber
	double M0x, M0y, M0z;			// point of application of F0 respect the block reference system
				
	if(forces.DC.is_available())	// use mesh file
	{
		M0x = forces.DC.M0x, M0y = forces.DC.M0y;
		F0z = forces.DC.F0z;
	}
	else	// use simplified approach
	{
		M0x = -1.0*rB*ADC, M0y = 0, M0z = 0;
		F0z = -1.0*ADC;	
	}

	// get axial force and moments

	forces.xFBz = forces.yFBz = 0;
	
	// summing for all pistons
	for(int i = 0; i < z; i++) 
	{
		double phi = pistons.phi[i];					// angular position for piston i
		double pDC = pfile.getp_rad(phi)[0];	// displacement chamber pressure
	
		// add to the total axial force the contribution due to p in DC
		double FBDi = pDC*F0z;	
		double FTBi = -1.0*forces.FTK.getf(phi); // this is FTB
		FBan += (FBDi + FTBi);
		forces.xFBz += (FBDi + FTBi)*(rB*sin(phi));
		forces.yFBz += (FBDi + FTBi)*(rB*cos(phi));
	
		// add to the total moments
		MBx += pDC*(M0x*cos(phi) + M0y*sin(phi));
		MBy += pDC*(-1.0*M0x*sin(phi) + M0y*cos(phi));

	}
	
	// ------------- transfert results to forceblockgap structure ------------ //
	
	forces.FBz = FBan;
	forces.MBx = MBx;
	forces.MBy = MBy;
	forces.xFBz /= FBan;
	forces.yFBz /= FBan;
	
	//cout << "FBz " << forces.FBz << "[N]" << endl;
	//cout << "MBx " << forces.MBx << "[N]" << endl;
	//cout << "MBy " << forces.MBy << "[N]" << endl;

}
// ------------------------------------------------------------------------- //
void cblock_gap_main::calc_dF()
{
	// ------------------------ copy from caspar_input ----------------------- //
	
	const double RBa = 0.5*in.data.geometry.dBa;

	// ----------------------------------------------------------------------- //

	// transformation matrix 
	double val[9] = 
	{
		1,			1,										1,
		-RBa,		0.5*RBa,							0.5*RBa,
		0,			(sqrt(3.0)/2.0)*RBa,	-(sqrt(3.0)/2.0)*RBa
	};

	matrix A(3,3,val);
		
	// --------------------------- External forces  -------------------------- //

	// clean
	forces.Fext[0] = 0.0;
	forces.Fext[1] = 0.0;
	forces.Fext[2] = 0.0;
	
	// external forces vector
	double be[3] = { forces.FBz, forces.MBx, forces.MBy };
	
	matrix b = matrix(3, 1, be);

	matrix transformed = matrix(A.inv()*b);
	
	forces.Fext[0] = transformed[0][0];
	forces.Fext[1] = transformed[1][0];
	forces.Fext[2] = transformed[2][0];

	
	// --------------------------- Fluid forces  ----------------------------- //
	
	// clean
	forces.Ffluid[0] = 0.0;
	forces.Ffluid[1] = 0.0;
	forces.Ffluid[2] = 0.0;

	// fluid forces vector
	double bf[3] = { forces.FfBz, forces.MfBx, forces.MfBy };

	b = matrix(3, 1, bf);
	transformed = matrix(A.inv()*b);

	forces.Ffluid[0] = transformed[0][0];
	forces.Ffluid[1] = transformed[1][0];
	forces.Ffluid[2] = transformed[2][0];

	// update dF
	forces.dF[0] = forces.Fext[0] + forces.Ffluid[0];
	forces.dF[1] = forces.Fext[1] + forces.Ffluid[1];
	forces.dF[2] = forces.Fext[2] + forces.Ffluid[2];
}
// ------------------------------------------------------------------------- //
double cblock_gap_main::calc_friction_moment()
{

	// ------------------------ copy from caspar_input ----------------------- //
	
	double omega = in.data.operating_conditions.speed;

	// ----------------------------------------------------------------------- //
	
	// get the friction moment using the dissipation
	double phitot = 0;
	for(int i=0; i<mesh.energy.mnq; i++)
	{
		if(mesh.energy.elements[i].ty == FLUID)
			phitot += fields.phid.in[i];
	}
	
	return (phitot/omega);
}
// ------------------------------------------------------------------------- //
double cblock_gap_main::calc_leakage()
{

	// ------------------------ copy from caspar_input ----------------------- //
	
	double r_groove_in = 0.50*in.data.geometry.d_groove_in;
	double r_groove_out = 0.50*in.data.geometry.d_groove_out;
	double r_gap_in = 0.50*in.data.geometry.d_gap_in;
	double r_gap_out = 0.50*in.data.geometry.d_gap_out;

	// ----------------------------------------------------------------------- //

	double leak = 0;
	
	// ------------------------- inner boundary ------------------------------ //

	// integrate the velocity over the inner boundary
	for(int k = 0, id = 0; k < mesh.energy.q; k++)
	{
		for(int j = 0; j < mesh.energy.n; j++, id++)
		{
			int id3D = mesh.energy.mn*k + j;
			// if velocity is pointing inward is a leakage through the inner bound
			if(fields.fs.in[id3D] < 0)
				leak += fabs(fields.fs.in[id3D]/fields.rho.in[id3D]);
		}
	}

	// ------------------------- outer boundary ------------------------------ //
	
	// ---------------------- pump with outer bearing ------------------------ //
	
	if(r_groove_in > 0)	// there is the outer bearing
	{
		// get the row index (radial direction) of elements closest to the groove
		int row = 0;
		for(int i=0; i<mesh.energy.m; i++)
		{
			if(mesh.energy.elements[mesh.energy.n*i].r < r_groove_in)
				row = i;
			else
				break;
		}
		
		// make sure we have a fluid element
		//while(mesh.energy.elements[row*mesh.energy.ne].ty != FLUID)
		//	row--;

		// integrate the velocity over the groove
		for(int k = 0; k < mesh.energy.q; k++)
		{
			for(int j = 0; j < mesh.energy.n; j++)
			{
				int id3D = mesh.energy.mn*k + row*mesh.energy.n + j;
				
				// if velocity is pointing outward is a leakage through the groove
				if(fields.fn.in[id3D] > 0)
					leak += fabs(fields.fn.in[id3D]/fields.rho.in[id3D]);
			}
		}
	}
	// --------------------- pump with no outer bearing ---------------------- //
	else
	{
		// integrate the velocity over the outer boundary
		for(int k = 0, id = 0; k < mesh.energy.q; k++)
		{
			for(int j = 0; j < mesh.energy.n; j++, id++)
			{
				int id3D = mesh.energy.mn*k + (mesh.energy.m-1)*mesh.energy.n + j;
				// if velocity is pointing outward is a leakage through the outer bound
				if(fields.fn.in[id3D] > 0)
					leak += fabs(fields.fn.in[id3D]/fields.rho.in[id3D]);
			}
		}
	}

	return leak;
}
// ------------------------------------------------------------------------- //
void cblock_gap_main::calc_mass_flow()
{
	// for each element calculate:
	//
	//		- the mass flow at the faces
	//		- the velocity at the cell centroid
	//

	double omega = in.data.operating_conditions.speed;

	int m = mesh.energy.m;
	int n = mesh.energy.n;
	int q = mesh.energy.q;
	
	// loop to all the element
	for(int k = 0, id3D = 0; k < q; k++)
	{
		for(int i = 0, id2D = 0; i < m; i++)
		{
			for(int j = 0; j < n; j++, id2D++, id3D ++)
			{

				gap_elm e3D = mesh.energy.elements[id3D];
				gap_elm e2D = mesh.reynolds.elements[id2D];

				double dr = e3D.dr;
				double r = e3D.r;
				double dtheta = e3D.dtheta;
				double h = film.hcb.in[id2D] - film.hvp.in[id2D];
				double hcb = film.hcb.in[id2D];
				double hvp = film.hvp.in[id2D];

				double dz = h/q;
				double z = film.hvp.in[id2D] + 0.5*dz + dz*k;
				
				// gradient at the faces
				double gpw, gpe, gps, gpn;
				// west
				if(mesh.reynolds.elements[e2D.w].ty == FLUID)
					gpw = (fields.p.in[id2D] - fields.p.in[e2D.w])/(r*dtheta);
				else
					gpw = (fields.p.in[id2D] - fields.p.in[e2D.w])/(0.5*r*dtheta);
				// east
				if(mesh.reynolds.elements[e2D.e].ty == FLUID)
					gpe = (fields.p.in[e2D.e] - fields.p.in[id2D])/(r*dtheta);
				else
					gpe = (fields.p.in[e2D.e] - fields.p.in[id2D])/(0.5*r*dtheta);

				// south
				if(e2D.s > -1)
				{
					if(mesh.reynolds.elements[e2D.s].ty == FLUID)
						gps = (fields.p.in[id2D] - fields.p.in[e2D.s])/dr;
					else
						gps = (fields.p.in[id2D] - fields.p.in[e2D.s])/(0.5*dr);
				}
				else
					gps = (fields.p.in[id2D] - fields.p.bound.inner[j])/(0.5*dr);
				// north
				if(e2D.n > -1)
					if(mesh.reynolds.elements[e2D.n].ty == FLUID)
						gpn = (fields.p.in[e2D.n] - fields.p.in[id2D])/dr;
					else
						gpn = (fields.p.in[e2D.n] - fields.p.in[id2D])/(0.5*dr);
				else
					gpn = (fields.p.bound.outer[j] - fields.p.in[id2D])/(0.5*dr);
					

				//extract pressure gradient at the faces
				fields.g_pw.in[id3D] = (fields.p.in[id2D] - fields.p.in[e2D.w]);
				fields.g_pe.in[id3D] = (fields.p.in[id2D] - fields.p.in[e2D.e]);
				fields.g_ps.in[id3D] = (fields.p.in[id2D] - fields.p.in[e2D.s]);
				fields.g_pn.in[id3D] = (fields.p.in[id2D] - fields.p.in[e2D.n]);

				//extract pressure at the faces
				fields.pw.in[id2D] = (fields.p.in[id2D] + fields.p.in[e2D.w])/2;
				fields.pe.in[id2D] = (fields.p.in[id2D] + fields.p.in[e2D.e])/2;
				fields.pn.in[id2D] = (fields.p.in[id2D] + fields.p.in[e2D.n])/2;
				fields.ps.in[id2D] = (fields.p.in[id2D] + fields.p.in[e2D.s])/2;

				// viscosity at the faces
				double mu_w, mu_e, mu_s, mu_n;
				double rho_w, rho_e, rho_s, rho_n;
				
				mu_w = 0.5*(fields.mu.in[id3D] + fields.mu.in[e3D.w]);
				mu_e = 0.5*(fields.mu.in[id3D] + fields.mu.in[e3D.e]);
				rho_w = 0.5*(fields.rho.in[id3D] + fields.rho.in[e3D.w]);
				rho_e = 0.5*(fields.rho.in[id3D] + fields.rho.in[e3D.e]);

				if(e3D.s > -1)
				{
					mu_s = 0.5*(fields.mu.in[id3D] + fields.mu.in[e3D.s]);
					rho_s = 0.5*(fields.rho.in[id3D] + fields.rho.in[e3D.s]);
				}
				else
				{
					mu_s = fields.mu.in[id3D];
					rho_s = fields.rho.in[id3D];
				}
				if(e3D.n > -1)
				{
					mu_n = 0.5*(fields.mu.in[id3D] + fields.mu.in[e3D.n]);
					rho_n = 0.5*(fields.rho.in[id3D] + fields.rho.in[e3D.n]);
				}
				else
				{
					mu_n = fields.mu.in[id3D];
					rho_n = fields.rho.in[id3D];
				}

				// ------------------ film thickness at the faces ------------------ //

				double hcb_e, hcb_w, hcb_s, hcb_n;
				double hvp_e, hvp_w, hvp_s, hvp_n;
				double he, hw, hs, hn;
				double ze, zw, zs, zn;

				// west
				if(e2D.w == FLUID)
				{
					hcb_w = 0.5*(film.hcb.in[id2D] + film.hcb.in[e2D.w]);
					hvp_w = 0.5*(film.hvp.in[id2D] + film.hvp.in[e2D.w]);
				}
				else
				{
					hcb_w = film.hcb.in[id2D];
					hvp_w = film.hvp.in[id2D];
				}
				// east
				if(e2D.e == FLUID)
				{
					hcb_e = 0.5*(film.hcb.in[id2D] + film.hcb.in[e2D.e]);
					hvp_e = 0.5*(film.hvp.in[id2D] + film.hvp.in[e2D.e]);
				}
				else
				{
					hcb_e = film.hcb.in[id2D];
					hvp_e = film.hvp.in[id2D];
				}
				// south
				if(e2D.s == FLUID)
				{
					hcb_s = 0.5*(film.hcb.in[id2D] + film.hcb.in[e2D.s]);
					hvp_s = 0.5*(film.hvp.in[id2D] + film.hvp.in[e2D.s]);
				}
				else
				{
					hcb_s = film.hcb.in[id2D];
					hvp_s = film.hvp.in[id2D];
				}
				// north
				if(e2D.n == FLUID)
				{
					hcb_n = 0.5*(film.hcb.in[id2D] + film.hcb.in[e2D.n]);
					hvp_n = 0.5*(film.hvp.in[id2D] + film.hvp.in[e2D.n]);
				}
				else
				{
					hcb_n = film.hcb.in[id2D];
					hvp_n = film.hvp.in[id2D];
				}

				hw = hcb_w - hvp_w;
				he = hcb_e - hvp_e;
				hs = hcb_s - hvp_s;
				hn = hcb_n - hvp_n;
			
				// ----------------- cell centroid axial position ------------------ //

				zw = hvp_w + (k + 0.5)*(hw/q);
				ze = hvp_e + (k + 0.5)*(he/q);
				zs = hvp_s + (k + 0.5)*(hs/q);
				zn = hvp_n + (k + 0.5)*(hn/q);
								
				// ---------------------- velocity at the faces -------------------- //
				
				double Vsr, Vnr, Vwt, Vet;
				double Vsr_p, Vnr_p, Vwt_p, Vet_p, Vwt_c, Vet_c;

				// poiseuille
				Vsr_p = (1.0/(2.0*mu_s))*gps*(pow(zs, 2) - zs*(hcb_s + hvp_s) + hcb_s*hvp_s);
				Vnr_p = (1.0/(2.0*mu_n))*gpn*(pow(zn, 2) - zn*(hcb_n + hvp_n) + hcb_n*hvp_n);
				Vwt_p = (1.0/(2.0*mu_w))*gpw*(pow(zw, 2) - zw*(hcb_w + hvp_w) + hcb_w*hvp_w);
				Vet_p = (1.0/(2.0*mu_e))*gpe*(pow(ze, 2) - ze*(hcb_e + hvp_e) + hcb_e*hvp_e);
							
				// couette
				Vwt_c = omega*r*(zw - hvp_w)/hw;
				Vet_c = omega*r*(ze - hvp_e)/he;

				// total
				Vsr = Vsr_p;											// just poiseuille
				Vnr = Vnr_p;											// just poiseuille
				Vwt = Vwt_p + Vwt_c;		// poiseuille + couette - mesh
				Vet = Vet_p + Vet_c;		// poiseuille + couette - mesh

				// --------------------- mass flow at the faces -------------------- //

				double fs = rho_s*Vsr*(r - 0.5*dr)*dtheta*(hs/q);	
				double fn = rho_n*Vnr*(r + 0.5*dr)*dtheta*(hn/q);	
				double fw = rho_w*Vwt*dr*(hw/q);										
				double fe = rho_e*Vet*dr*(he/q);											
				
							// assign to the fields 
				
				fields.mdotw.in[id3D] = fw;
				fields.mdote.in[id3D] = fe;
				fields.mdots.in[id3D] = fs;
				fields.mdotn.in[id3D] = fn;

				// ------------- ensure the continuity to be satisfied ------------- //

				double delta_f = (fn - fs + fe - fw);
				double f_tot = fabs(fs) + fabs(fn) + fabs(fw) + fabs(fe);

				fs += (fabs(fs)/f_tot)*delta_f;
				fn -= (fabs(fn)/f_tot)*delta_f;
				fw += (fabs(fw)/f_tot)*delta_f;
				fe -= (fabs(fe)/f_tot)*delta_f;

				// assign to the fields 
				
				fields.fw.in[id3D] = fw;
				fields.fe.in[id3D] = fe;
				fields.fs.in[id3D] = fs;
				fields.fn.in[id3D] = fn;

				// velocity fields at faces 
				fields.ve.in[id3D] = Vet;
				fields.vw.in[id3D] = Vwt;
				fields.vn.in[id3D] = Vnr;
				fields.vs.in[id3D] = Vsr;


				// density fields at the faces
				fields.rhoe.in[id3D] = rho_e;
				fields.rhow.in[id3D] = rho_w;
				fields.rhon.in[id3D] = rho_n;
				fields.rhos.in[id3D] = rho_s;



				// poiseuille
				fields.Vp.in[id3D][0] = 0.5*(Vsr_p + Vnr_p);	// r
				fields.Vp.in[id3D][1] = 0.5*(Vwt_p + Vet_p);	// theta
				fields.Vp.in[id3D][2] = 0;						// z
				
				//couette
				fields.Vc.in[id3D][0] = 0.0;					// r
				fields.Vc.in[id3D][1] = 0.5*(Vwt_c + Vet_c);	// theta
				fields.Vc.in[id3D][2] = 0;						// z

				//total
				fields.V.in[id3D][0] = 0.5*(Vsr + Vnr);			// r
				fields.V.in[id3D][1] = 0.5*(Vwt + Vet);			// theta
				fields.V.in[id3D][2] = 0;						// z

				// set inner/outer boundaries
				if(e3D.n == -1)
				{
					fields.V.bound.outer[j] = fields.V.in[id3D];
					fields.Vc.bound.outer[j] = fields.Vc.in[id3D];
					fields.Vp.bound.outer[j] = fields.Vp.in[id3D];
				}
				if(e3D.s == -1)
				{
					fields.V.bound.inner[j] = fields.V.in[id3D];
					fields.Vc.bound.inner[j] = fields.Vc.in[id3D];
					fields.Vp.bound.inner[j] = fields.Vp.in[id3D];
				}
				// set top/bottom boundaries
				if(e3D.b == -1)
				{
					fields.V.bound.bottom[id2D] = 0;
					fields.Vc.bound.bottom[id2D] = 0;
					fields.Vp.bound.bottom[id2D] = 0;
				}
				if(e3D.t == -1)
				{
					fields.V.bound.top[id2D] = omega*r*(z - hvp)/h;
					fields.Vc.bound.top[id2D] = omega*r*(z - hvp)/h;
					fields.Vp.bound.top[id2D] = 0;
				}
		
			}
		}
	}

}
// ------------------------------------------------------------------------- //
void cblock_gap_main::update_density()
{

	// mesh dimensions
	int m = mesh.energy.m;
	int n = mesh.energy.n;
	int q = mesh.energy.q;
	int mn = m*n;

	// -------------- set the 3D rho ----------------- //
	
	// set the internal field
	for(int k=0, id3d = 0; k<mesh.energy.q; k++)
	{
		for(int ij=0; ij<mesh.reynolds.mn; ij++, id3d++)
		{		
			fields.rho.in[id3d] = lubricant->get_rho(fields.p.in[ij], fields.T.in[id3d]);
		}
	}
	// set top and bottom boundaries
	for(int ij=0; ij<mesh.reynolds.mn; ij++)
	{
		fields.rho.bound.top[ij] = lubricant->get_rho(fields.p.in[ij], fields.T.bound.top[ij]);
		fields.rho.bound.bottom[ij] = lubricant->get_rho(fields.p.in[ij], fields.T.bound.bottom[ij]);
	}
	// set inner and outer boundaries
	for(int k=0, id = 0; k<mesh.energy.q; k++)
	{
		for(int j=0; j<mesh.reynolds.n; j++, id++)
		{
			fields.rho.bound.inner[k*n + j] = fields.rho.in[mn*k + j];
			fields.rho.bound.outer[k*n + j] = fields.rho.in[mn*k + (m-1)*n + j];	
		}
	}

	// ----------------- get the 2D rho from the 3D rho ------------------ //
			
	// internal domain
	for(int ij=0; ij<mesh.reynolds.mn; ij++)
	{
		fields.rho2D.in[ij] = 0;
		for(int k=0; k<mesh.energy.q; k++)
			fields.rho2D.in[ij] += fields.rho.in[mesh.reynolds.mn*k + ij];
		fields.rho2D.in[ij] = fields.rho2D.in[ij]/mesh.energy.q;
	}
	// boundaries
	for(int j=0; j<mesh.reynolds.n; j++)
	{
		fields.rho2D.bound.inner[j] = 0;
		fields.rho2D.bound.outer[j] = 0;
		for(int k=0; k<mesh.energy.q; k++)
		{
			fields.rho2D.bound.inner[j] += fields.rho.bound.inner[mesh.reynolds.n*k + j];
			fields.rho2D.bound.outer[j] += fields.rho.bound.outer[mesh.reynolds.n*k + j];
		}
		fields.rho2D.bound.inner[j] /= mesh.energy.q;
		fields.rho2D.bound.outer[j] /= mesh.energy.q;
	}

}
// ------------------------------------------------------------------------- //
void cblock_gap_main::update_h_previous()
{

	film.h_previous = film.h;
	fields.T_prev = fields.T;

}

// ------------------------------------------------------------------------- //
void cblock_gap_main::update_viscosity()
{

	// mesh dimensions
	int m = mesh.energy.m;
	int n = mesh.energy.n;
	int q = mesh.energy.q;
	int mn = m*n;

	// -------------- set the 3D mu ----------------- //
	
	// set the internal field
	for(int k=0, id3d = 0; k<q; k++)
	{
		for(int ij=0; ij<mn; ij++, id3d++)
		{		
			fields.mu.in[id3d] = lubricant->get_mu(fields.p.in[ij], fields.T.in[id3d]);
		}
	}
	// set top and bottom boundaries
	for(int ij=0; ij<mn; ij++)
	{
		fields.mu.bound.top[ij] = lubricant->get_mu(fields.p.in[ij], fields.T.bound.top[ij]);
		fields.mu.bound.bottom[ij] = lubricant->get_mu(fields.p.in[ij], fields.T.bound.bottom[ij]);
	}
	// set inner and outer boundaries
	for(int k=0; k<q; k++)
	{
		for(int j=0; j<n; j++)
		{
			fields.mu.bound.inner[k*n + j] = fields.mu.in[mn*k + j];
			fields.mu.bound.outer[k*n + j] = fields.mu.in[mn*k + (m-1)*n + j];	
		}
	}

	// ----------------- get the 2D mu from the 3D mu ------------------ //
			
	// internal domain
	for(int ij=0; ij<mesh.reynolds.mn; ij++)
	{
		fields.mu2D.in[ij] = 0;
		for(int k=0; k<mesh.energy.q; k++)
			fields.mu2D.in[ij] += fields.mu.in[mesh.reynolds.mn*k + ij];
		fields.mu2D.in[ij] = fields.mu2D.in[ij]/mesh.energy.q;
	}
	// boundaries
	for(int j=0; j<mesh.reynolds.n; j++)
	{
		fields.mu2D.bound.inner[j] = 0;
		fields.mu2D.bound.outer[j] = 0;
		for(int k=0; k<mesh.energy.q; k++)
		{
			fields.mu2D.bound.inner[j] += fields.mu.bound.inner[mesh.reynolds.n*k + j];
			fields.mu2D.bound.outer[j] += fields.mu.bound.outer[mesh.reynolds.n*k + j];
		}
		fields.mu2D.bound.inner[j] /= mesh.energy.q;
		fields.mu2D.bound.outer[j] /= mesh.energy.q;
	}
}
// ------------------------------------------------------------------------- //
void cblock_gap_main::set_film_thickness()
{
	
	// ------------------------ copy from caspar_input ----------------------- //

	double RBa = 0.5*in.data.geometry.dBa;

	// ----------------------------------------------------------------------- //

	double h1 = film.h1;
	double h2 = film.h2;
	double h3 = film.h3;

	// internal field
	for(int i = 0; i < mesh.reynolds.mn; i++) 
	{
		gap_elm e = mesh.reynolds.elements[i];
		double r = e.r;
		double theta = e.theta;
		film.hcb_rigid.in[i] = 
			r*sin(theta)*(sqrt(1.0/3.0)/RBa)*(h2 - h3) + 
			r*cos(theta)*(1.0/(3.0*RBa))*(2*h1 - h2 - h3) + 
			(1.0/3.0)*(h1 + h2 + h3);
	}

	// inner boundary
	for(int j = 0; j < mesh.reynolds.n; j++) 
	{
		gap_elm e = mesh.reynolds.elements[j];
		double r = e.r;
		double dr = e.dr;
		double theta = e.theta;

		film.hcb_rigid.bound.inner[j] = 
			(r - dr)*sin(theta)*(sqrt(1.0/3.0)/RBa)*(h2 - h3) + 
			r*cos(theta)*(1.0/(3.0*RBa))*(2*h1 - h2 - h3) + 
			(1.0/3.0)*(h1 + h2 + h3);
	}

	// outer boundary
	for(int j = 0; j < mesh.reynolds.n; j++) 
	{
		gap_elm e = mesh.reynolds.elements[(mesh.reynolds.m - 2)*mesh.reynolds.n + j];
		double r = e.r;
		double dr = e.dr;
		double theta = e.theta;

		film.hcb_rigid.bound.outer[j] = 
			(r + dr)*sin(theta)*(sqrt(1.0/3.0)/RBa)*(h2 - h3) + 
			r*cos(theta)*(1.0/(3.0*RBa))*(2*h1 - h2 - h3) + 
			(1.0/3.0)*(h1 + h2 + h3);
	}

	// set the valve plate rigid position to zero
	film.hvp_rigid = 0;

	// apply macro geometry if required
	if(macro_CB != 0)
		macro_CB->apply(mesh.reynolds.rotation_steps);
	
	// apply macro geometry if required
	if(macro_VP != 0)
		macro_VP->apply();

	// block translation
	film.dz = ((h1 + h2 + h3)/3.0);
	
	// block rotations
	double dz_y = RBa*(1.0/(3.0*RBa))*(2*h1 - h2 - h3) + (1.0/3.0)*(h1 + h2 + h3);
	double dz_x = RBa*(sqrt(1.0/3.0)/RBa)*(h2 - h3) + (1.0/3.0)*(h1 + h2 + h3);
	film.rx = atan2(dz_y - film.dz, RBa);
	film.ry = atan2(film.dz - dz_x, RBa);

	film.hcb = film.hcb_rigid;
	film.hvp = film.hvp_rigid;


	// ------------------------ include the deformation ----------------------- //

	if(EHD_CB || EHD_VP || TH_CB || TH_VP)
	{
		// minimum film thickness
		double hmin = in.data.options_block.numeric.hmin;
		hmin = hmin > 0.1e-6 ? hmin : 0.1e-6;	// no less than 0.1 um

		// cylinder block
		scalar_field dhcb_thermal = film.dhcb_thermal.cshift(mesh.reynolds.rotation_steps);
		film.hcb = film.hcb_rigid + dhcb_thermal + film.dhcb_EHD;

		// valve plate
		film.hvp = film.hvp_rigid + film.dhvp_thermal + film.dhvp_EHD;
		
		// ----------------------- get the contact field ----------------------- //
		
		film.contact = 0.0;
		vector<int> cnt_cells(0);
		
		for(int id=0; id<film.h.mn; id++)
		{
			if(film.h.mesh->elements[id].ty == FLUID)
			{
				// start to consider contat when h < 2*hmin
				if(film.hcb.in[id] - film.hvp.in[id] < 2*hmin)
				{
					film.contact.in[id] = 2*hmin - (film.hcb.in[id] - film.hvp.in[id]);
					cnt_cells.push_back(id);
				}
				// changed in 07/14/2016 saturate the film @ 0.1 um --> saturate @ defined hmin in the inputs 
				if(film.hcb.in[id] - film.hvp.in[id] < hmin)
				{
					film.hcb.in[id] = film.hvp.in[id] + hmin;
					film.h.in[id] = film.hcb.in[id] - film.hvp.in[id];
				}
			}
		}

	
		int n = mesh.reynolds.n;
		int m = mesh.reynolds.m;

		// set the boundaries
		for(int j=0; j<n; j++)
		{
			film.hcb.bound.inner[j] = film.hcb.in[j];
			film.hcb.bound.outer[j] = film.hcb.in[(m-1)*n + j];
			film.hvp.bound.inner[j] = film.hvp.in[j];
			film.hvp.bound.outer[j] = film.hvp.in[(m-1)*n + j];
		}
	}

	// ------------------------------------------------------------------------ //

	film.h = film.hcb - film.hvp;

	// update film.hmin
	film.hmin = 1e10;
	for(int i=0; i<mesh.reynolds.mn; i++)
	{
		if(mesh.reynolds.elements[i].ty == FLUID)
		{
			if(film.h.in[i] < film.hmin)
				film.hmin = film.h.in[i];
		}
	}

	// update the mesh shape and elements dimension
	update_mesh();

}
// ------------------------------------------------------------------------- //
double cblock_gap_main::hmin() const
{
	double hmin = 1e10;
	
	for(int id = 0; id<mesh.reynolds.mn; id++)
	{
		if(mesh.reynolds.elements[id].ty == FLUID)
		{
			if(film.h.in[id] < hmin)
				hmin = film.h.in[id];
		}
	}

	return hmin;
}

double cblock_gap_main::calcFTK(double phi_rad)
{
	// ------------------------ copy from caspar_input ----------------------- //
	
	double rK = 0.50*in.data.geometry.dK;
	double phi = phi_rad;
	double FTK;
	// ----------------------------------------------------------------------- //
	
	FTK = forces.FTK.getf(phi);
	
	if (FTK == 0)
	{
			double area;
		
			area = 2 * pi * rK * pistons.lvar[0];

			double mu = lubricant->get_mu(in.data.operating_conditions.HP,in.data.operating_conditions.T_Leak);
			FTK = mu * pistons.vK[0] * area / pistons.radial_clearance;
			
			/*
			Log << "\nViscosity: ";
			Log << mu << "\n";
			Log << "\nTemperature: ";
			Log << in.data.operating_conditions.T_Leak << "\n";
			Log << "\nPressure: ";
			Log << in.data.operating_conditions.HP << "\n";
			Log << "\nRadial clearance: ";
			Log << pistons.radial_clearance << "\n";
			Log << "\nArea: ";
			Log << area << "\n";
			Log << "\nFTK: ";
			Log << FTK << "\n";
			*/
	}

	
	return FTK;
}
// ------------------------------------------------------------------------- //
double cblock_gap_main::hmax() const
{
	double hmax = -1e10;
	
	for(int id = 0; id<mesh.reynolds.mn; id++)
	{
		if(mesh.reynolds.elements[id].ty == FLUID)
		{
			if(film.h.in[id] > hmax)
				hmax = film.h.in[id];
		}
	}

	return hmax;
}
// ------------------------------------------------------------------------- //
double cblock_gap_main::havg() const
{

	double havg = 0;
	int n = 0;
	
	for(int id = 0; id<mesh.reynolds.mn; id++)
	{
		if(mesh.reynolds.elements[id].ty == FLUID)
		{
			havg += film.h.in[id];
			n++;
		}
	}

	return havg/n;
}
// ------------------------------------------------------------------------- //
void cblock_gap_main::set_squeeze()
{

	// ------------------------ copy from caspar_input ----------------------- //

	double RBa = 0.5*in.data.geometry.dBa;

	// ----------------------------------------------------------------------- //

	double dhdt1 = squeeze.dhdt1;
	double dhdt2 = squeeze.dhdt2;
	double dhdt3 = squeeze.dhdt3;
	
	for(int i = 0; i <  mesh.reynolds.mn; i++) 
	{
		gap_elm e = mesh.reynolds.elements[i];
		double r = e.r;
		double theta = e.theta;
		squeeze.dhdt.in[i] = 
			r*sin(theta)*(sqrt(1.0/3.0)/RBa)*(dhdt2 - dhdt3) + 
			r*cos(theta)*(1.0/(3.0*RBa))*(2*dhdt1 - dhdt2 - dhdt3) + 
			(1.0/3.0)*(dhdt1 + dhdt2 + dhdt3);
	}

	// inner and outer boundaries
	for(int j = 0; j < mesh.reynolds.n; j++) 
	{
		gap_elm e = mesh.reynolds.elements[j];
		double r = e.r;
		double theta = e.theta;
		double dr = e.dr;

		squeeze.dhdt.bound.inner[j] = 
			(r - dr)*sin(theta)*(sqrt(1.0/3.0)/RBa)*(dhdt2 - dhdt3) + 
			r*cos(theta)*(1.0/(3.0*RBa))*(2*dhdt1 - dhdt2 - dhdt3) + 
			(1.0/3.0)*(dhdt1 + dhdt2 + dhdt3);

		// go to outer
		e = mesh.reynolds.elements[(mesh.reynolds.m - 2)*mesh.reynolds.n + j];
		squeeze.dhdt.bound.outer[j] = 
			(r + dr)*sin(theta)*(sqrt(1.0/3.0)/RBa)*(dhdt2 - dhdt3) + 
			r*cos(theta)*(1.0/(3.0*RBa))*(2*dhdt1 - dhdt2 - dhdt3) + 
			(1.0/3.0)*(dhdt1 + dhdt2 + dhdt3);

	}

}
// ------------------------------------------------------------------------- //
void cblock_gap_main::update_mesh()
{
	mesh.reynolds.update_thickness(film.hcb, film.hvp);
	mesh.energy.update_thickness(film.hcb, film.hvp);
}
// ------------------------------------------------------------------------- //
void cblock_gap_main::define_interpl_grids()
{
	
	for(int i=0; i<CBf.nn(); i++)
	{
		CBf.nodes[i] = mesh.reynolds.nodes[i];
		VPf.nodes[i] = mesh.reynolds.nodes[i];
	}
	// define intepolationgrid for cylinder block
	for(int i=0; i<CBf.ne(); i++)
	{
		for(int j = 0; j<4; j++)
		{
			CBf.elements[i][j] = mesh.reynolds.elements[i].nds[j];
			
			// define the active elements for block
			if(mesh.reynolds.cb_elms_0[i] == FLUID)
				CBf.active[i] = true;
			else
				CBf.active[i] = false;
		}
	}
	// define intepolationgrid for valve plate
	for(int i=0; i<VPf.ne(); i++)
	{
		for(int j = 0; j<4; j++)
		{
			VPf.elements[i][j] = mesh.reynolds.elements[i].nds[j];
			// define the active elements for valveplate
			if(mesh.reynolds.vp_elms[i] == FLUID)
				VPf.active[i] = true;
			else
				VPf.active[i] = false;
		}
	}
}
// ------------------------------------------------------------------------- //
void cblock_gap_main::interpl_pressure()
{

	// define block interpolation object with current pressure field pressure
	if(EHD_CB)
	{
		// pressure scalar field @ the 0 position
		scalar_field p = fields.p.cshift(-mesh.reynolds.rotation_steps);

		for(int i=0; i<CBf.ne(); i++)
		{
			if(CBf.active[i])
				CBf.cells_data[i] = p.in[i];
			else
				CBf.cells_data[i] = 0;
		}

	}
	// define valve plate interpolation object with current pressure field pressure
	if(EHD_VP)
	{
		for(int i=0; i<VPf.ne(); i++)
		{
			if(VPf.active[i])
				VPf.cells_data[i] = fields.p.in[i];
			else
				VPf.cells_data[i] = 0;
		}
	}
}
// ------------------------------------------------------------------------- //
void cblock_gap_main::update_Tcb()
{

	double relax = in.data.options_block.numeric.relax_CB;

	for(int i=0; i<CBf.ne(); i++)
	{
		if(CBf.active[i])
			fields.Tcb.in[i] = fields.Tcb.in[i] + relax*(CBf.cells_data[i] - fields.Tcb.in[i]);
		else
			fields.Tcb.in[i] = 0;
	}
	
	// set the boundaries
	int m = mesh.reynolds.m, n = mesh.reynolds.n;
	for(int j=0; j<n; j++)
	{
		fields.Tcb.bound.inner[j] = fields.Tcb.in[j];
		fields.Tcb.bound.outer[j] = fields.Tcb.in[(m-1)*n + j];	
	}
}
// ------------------------------------------------------------------------- //
void cblock_gap_main::update_Tvp()
{

	double relax = in.data.options_block.numeric.relax_VP;

	for(int i=0; i<VPf.ne(); i++)
	{
		if(VPf.active[i])
			fields.Tvp.in[i] = fields.Tvp.in[i] + relax*(VPf.cells_data[i] - fields.Tvp.in[i]);
		else
			fields.Tvp.in[i] = 0;
	}
	
	// set the boundaries
	int m = mesh.reynolds.m, n = mesh.reynolds.n;
	for(int j=0; j<n; j++)
	{
		fields.Tvp.bound.inner[j] = fields.Tvp.in[j];
		fields.Tvp.bound.outer[j] = fields.Tvp.in[(m-1)*n + j];	
	}
}
// ------------------------------------------------------------------------- //
void cblock_gap_main::update_dhcb_thermal()
{
	double relax = in.data.options_block.numeric.relax_CB;

	// copy to a temporary scalar field
	scalar_field dh(&mesh.reynolds);
	for(int i=0; i<CBf.ne(); i++)
		dh.in[i] = (CBf.active[i]) ? CBf.cells_data[i] : 0;
	
	// get the minimum value
	int m = mesh.reynolds.m, n = mesh.reynolds.n, mn = m*n;
	double reference = 1e10;
	for(int i=0; i<mn; i++)
	{
		if(mesh.reynolds.cb_elms_0[i] == FLUID && dh.in[i] < reference)
			reference = dh.in[i];
	}
	// get the relative deformation
	for(int i=0; i<mn; i++)
	{
		if(mesh.reynolds.cb_elms_0[i] == FLUID)
			dh.in[i] = fabs(reference - dh.in[i]);
	}
	
	// copy to dhcb_thermal
	for(int i=0; i<mn; i++)
	{
		if(mesh.reynolds.cb_elms_0[i] == FLUID)
			film.dhcb_thermal.in[i] = film.dhcb_thermal.in[i] + relax*(dh.in[i] - film.dhcb_thermal.in[i]);
	}
	
	// set the boundaries
	for(int j=0; j<n; j++)
	{
		film.dhcb_thermal.bound.inner[j] = film.dhcb_thermal.in[j];
		film.dhcb_thermal.bound.outer[j] = film.dhcb_thermal.in[(m-1)*n + j];	
	}

}
// ------------------------------------------------------------------------- //
void cblock_gap_main::update_dhvp_thermal()
{
	double relax = in.data.options_block.numeric.relax_VP;

	// copy to a temporary scalar field
	scalar_field dh(&mesh.reynolds);
	for(int i=0; i<VPf.ne(); i++)
		dh.in[i] = (VPf.active[i]) ? VPf.cells_data[i] : 0;
	
	// get the minimum value
	int m = mesh.reynolds.m, n = mesh.reynolds.n, mn = m*n;
	double reference = -1e10;
	for(int i=0; i<mn; i++)
	{
		if(mesh.reynolds.vp_elms[i] == FLUID && dh.in[i] > reference)
			reference = dh.in[i];
	}
	// get the relative deformation
	for(int i=0; i<mn; i++)
	{
		if(mesh.reynolds.vp_elms[i] == FLUID)
			dh.in[i] = -fabs(reference - dh.in[i]);
	}
	
	// copy to dhvp_thermal
	for(int i=0; i<mn; i++)
	{
		if(mesh.reynolds.vp_elms[i] == FLUID)
			film.dhvp_thermal.in[i] = film.dhvp_thermal.in[i] + relax*(dh.in[i] - film.dhvp_thermal.in[i]);
	}
	
	// set the boundaries
	for(int j=0; j<n; j++)
	{
		film.dhvp_thermal.bound.inner[j] = film.dhvp_thermal.in[j];
		film.dhvp_thermal.bound.outer[j] = film.dhvp_thermal.in[(m-1)*n + j];	
	}

}
// ------------------------------------------------------------------------- //
static void force_big_endian(unsigned char *bytes)
{
	
	// this function forces the big endian binary format
  static int doneTest = 0;
  static int shouldSwap = 0;
  if (!doneTest)
  {
      int tmp1 = 1;
      unsigned char *tmp2 = (unsigned char *) &tmp1;
      if (*tmp2 != 0)
          shouldSwap = 1;
      doneTest = 1;
  }

  if (shouldSwap)
  {
      unsigned char tmp = bytes[0];
      bytes[0] = bytes[3];
      bytes[3] = tmp;
      tmp = bytes[1];
      bytes[1] = bytes[2];
      bytes[2] = tmp;
  }
}
// ------------------------------------------------------------------------- //
void cblock_gap_main::write_pcb_vtk(const char* name)
{
	vector<int> nodes_idx(mesh.reynolds.cb_gap_surf.elms.size()); // this is the staring index vector
	nodes_idx[0] = 0;
	int ntot = mesh.reynolds.cb_gap_surf.elms[0].vertices.size();

	for(unsigned int i=1; i<mesh.reynolds.cb_gap_surf.elms.size(); i++)
	{
		nodes_idx[i] = nodes_idx[i-1] + mesh.reynolds.cb_gap_surf.elms[i-1].vertices.size();
		ntot += mesh.reynolds.cb_gap_surf.elms[i].vertices.size();
	}

	ofstream vtk(name);
	if (!vtk.is_open()) 
	{
		cout << "Error opening " << name << ".vtk" << endl;
		exit(1);
	}

	// write VTK header
	vtk << "# vtk DataFile Version 2.0" << endl 
			<< "vtk output" << endl
			<< "ASCII" << endl 
			<< "DATASET UNSTRUCTURED_GRID" << endl 
			<< "POINTS " << ntot << " double" << endl;

	int cell_size = 0;

	// spherical valve plate //
	//double R = 0.5*in.data.geometry.d_spherical;
	//double xmax = 0.5*in.data.geometry.d_gap_out;
	//double z0 = sqrt(R*R - xmax*xmax);
	// ~sperical valve plate //

	for(unsigned int i=0; i<mesh.reynolds.cb_gap_surf.elms.size(); i++)
	{
		cell_size += (1 + 1 + mesh.reynolds.cb_gap_surf.elms[i].vertices.size());
		for(unsigned j=0; j<mesh.reynolds.cb_gap_surf.elms[i].vertices.size(); j++)
		{
			double x = mesh.reynolds.cb_gap_surf.elms[i].vertices[j].x();
			double y = mesh.reynolds.cb_gap_surf.elms[i].vertices[j].y();
			double z = 0;
			vtk << x << "\t" << y << "\t" << z << endl;
			
			// spherical valveplate //
			//double z = sqrt(R*R - (x*x + y*y)) - z0;
			//vtk << 1e3*x << "\t" << 1e3*y << "\t" << 1e3*z - 3.81 << endl;
			// ~spherical valveplate //
		}
	}

	vtk << "CELLS " << mesh.reynolds.cb_gap_surf.elms.size() << "\t" << cell_size << endl;

	for(unsigned int i = 0; i<mesh.reynolds.cb_gap_surf.elms.size(); i++) 
	{
		vtk << mesh.reynolds.cb_gap_surf.elms[i].vertices.size() + 1 << "\t";
		for(unsigned int j = 0; j<mesh.reynolds.cb_gap_surf.elms[i].vertices.size(); j++)
			vtk << nodes_idx[i] + j << "\t"; // index start from 0
		vtk << nodes_idx[i] << endl; //repeat the first
	}
	
	vtk << endl;

	vtk << "CELL_TYPES " << mesh.reynolds.cb_gap_surf.elms.size() << endl;
	for(unsigned int i=0; i<mesh.reynolds.cb_gap_surf.elms.size(); i++) 
	{
		vtk << 7 << endl;
	}

	vtk << "\nCELL_DATA " << mesh.reynolds.cb_gap_surf.elms.size() << endl;
	vtk << "\nSCALARS p[bar] float" << endl;
	vtk << "\nLOOKUP_TABLE default" << endl;

	for(unsigned int i=0; i<mesh.reynolds.cb_gap_surf.elms.size(); i++)
		vtk << forces.p_cb[i]/1e5 << endl;
		
	vtk.close();
}
// ------------------------------------------------------------------------- //
void cblock_gap_main::write_light_vtk(const char* name, double zscale)
{

	// ------------------------ copy from caspar_input ----------------------- //
	
	double R = 0.50*in.data.geometry.d_spherical;

	// ----------------------------------------------------------------------- //

	int m = mesh.reynolds.m;
	int n = mesh.reynolds.n;
	int q = mesh.energy.q;
	int mn_e = m*n;

	int mn = mesh.reynolds.m + 1;
	int nn = mesh.reynolds.n;
	int qn = mesh.energy.q + 1;
	int mn_n = mn*nn;

	int nfluid = mesh.reynolds.mn_f;
	
	int elmType = VTK_3D_ELM_TYPE;
	int elmVtx = VTK_3D_ELM_NDS;

	std::ofstream vtk;

	vtk.open(name, std::ios::binary | std::ios::out);

	// write VTK header
	vtk << 
		"# vtk DataFile Version 3.0" << std::endl <<
		"blockgap output" << std::endl <<
		"BINARY" << std::endl <<
		"DATASET UNSTRUCTURED_GRID" << std::endl << 
		"POINTS " << 2*mn_n << "\tfloat" << std::endl;

	float x,y,z;

	// bottom layer
	for(int i=0; i<mn_n; i++) 
	{
		x = static_cast<float>(mesh.reynolds.nodes[i].x());
		y = static_cast<float>(mesh.reynolds.nodes[i].y());
		z = zscale*static_cast<float>(mesh.energy.nodes[i].z());
		force_big_endian(reinterpret_cast<unsigned char*>(&x));
		vtk.write(reinterpret_cast<char*>(&x),sizeof(float));
		force_big_endian(reinterpret_cast<unsigned char*>(&y));
		vtk.write(reinterpret_cast<char*>(&y),sizeof(float));
		force_big_endian(reinterpret_cast<unsigned char*>(&z));
		vtk.write(reinterpret_cast<char*>(&z),sizeof(float));
	}
	// top layer
	for(int i=0; i<mn_n; i++) 
	{
		x = static_cast<float>(mesh.reynolds.nodes[i].x());
		y = static_cast<float>(mesh.reynolds.nodes[i].y());
		z = zscale*static_cast<float>(mesh.energy.nodes[(qn-1)*mn_n + i].z());
		force_big_endian(reinterpret_cast<unsigned char*>(&x));
		vtk.write(reinterpret_cast<char*>(&x),sizeof(float));
		force_big_endian(reinterpret_cast<unsigned char*>(&y));
		vtk.write(reinterpret_cast<char*>(&y),sizeof(float));
		force_big_endian(reinterpret_cast<unsigned char*>(&z));
		vtk.write(reinterpret_cast<char*>(&z),sizeof(float));
	}

	// -------------------------- cells connections -------------------------- //

	vtk << "\nCELLS " << nfluid << "\t" 
		  << (elmVtx + 1)*nfluid << std::endl;
	
	for(int i = 0; i < mn_e; i++) 
	{
		if(mesh.reynolds.elements[i].ty == FLUID)
		{
			int tmp = elmVtx;
			force_big_endian(reinterpret_cast<unsigned char*>(&tmp));
			vtk.write(reinterpret_cast<char*>(&tmp), sizeof(int));
			// bottom nodes
			for(int j=0; j<4; j++) 
			{
				int ndsIdx = mesh.reynolds.elements[i].nds[j];
				force_big_endian(reinterpret_cast<unsigned char*>(&ndsIdx));
				vtk.write(reinterpret_cast<char*>(&ndsIdx),sizeof(int));	
			}
			// top nodes
			for(int j=0; j<4; j++) 
			{
				int ndsIdx = mn_n + mesh.reynolds.elements[i].nds[j];
				force_big_endian(reinterpret_cast<unsigned char*>(&ndsIdx));
				vtk.write(reinterpret_cast<char*>(&ndsIdx),sizeof(int));	
			}
		}
	}

	// ------------------------------ cells type ----------------------------- //

	vtk << "\nCELL_TYPES " << nfluid << std::endl;
				
	force_big_endian(reinterpret_cast<unsigned char*>(&elmType));

	for(int i = 0; i < mn_e; i++) 
	{
		if(mesh.reynolds.elements[i].ty == FLUID)
		{
			int tmp = elmType;
			vtk.write(reinterpret_cast<char*>(&tmp),sizeof(int));
		}
	}	


	// --------------------------------- fields ------------------------------ //

	vtk << "\nCELL_DATA " << nfluid << std::endl;

	// --------------------------------- p ----------------------------------- //

	vtk << "SCALARS p[bar] float 1 " << std::endl 
			<< "LOOKUP_TABLE default" << std::endl;
		
	for(int id = 0; id<mesh.reynolds.mn; id++)
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(fields.p.in[id]/1e5);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	// --------------------------------- mu ----------------------------------- //

	vtk << "\nSCALARS mu[Pas] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
	
	for(int id = 0; id<mesh.reynolds.mn; id++)
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(fields.mu2D.in[id]);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}
	
	// --------------------------------- h ----------------------------------- //

	vtk << "\nSCALARS h[microns] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
	
	for(int id = 0; id<mesh.reynolds.mn; id++)
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(film.h.in[id]*1e6);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	// ------------------------------ contact -------------------------------- //

	vtk << "\nSCALARS contact[um] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
	
	for(int id = 0; id<mesh.reynolds.mn; id++)
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(1e6*film.contact.in[id]);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	// ------------------------------ dhdt -------------------------------- //

	vtk << "\nSCALARS dhdt[m/s] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
	
	for(int id = 0; id<mesh.reynolds.mn; id++)
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = squeeze.dhdt.in[id];
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	// ------------------------------ dhdt -------------------------------- //

	vtk << "\nSCALARS dhdt_total[m/s] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
	
	for(int id = 0; id<mesh.reynolds.mn; id++)
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = squeeze.dhdt_total.in[id];
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	// ------------------------------ dhdt_ehd -------------------------------- //

	vtk << "\nSCALARS dhdt_ehd[m/s] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
	
	for(int id = 0; id<mesh.reynolds.mn; id++)
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			double dhdt_ehd = squeeze.dhdt_hs.in[id] + squeeze.dhdt_hd.in[id];
			float val = dhdt_ehd;
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	

}
// ------------------------------------------------------------------------- //
void cblock_gap_main::write_vtk(const char* name, double zscale)
{

	// spherical gap radius (0 if the gap is flat)
	double R = 0.50*in.data.geometry.d_spherical;

	int ntot = static_cast<int>(mesh.energy.nodes.size());
	int etot = mesh.energy.mnq;
	int nfluid = mesh.energy.mnq_f;

	int m = mesh.energy.m;
	int n = mesh.energy.n;
	int mn = m*n;
	int q = mesh.energy.q;

	int elmType = VTK_3D_ELM_TYPE;
	int elmVtx = VTK_3D_ELM_NDS;
	
	std::ofstream vtk
	(
		string(name + string("_gap.vtk")).c_str(), 
		std::ios::binary | std::ios::out
	);
	
	// write VTK header
	vtk << 
		"# vtk DataFile Version 3.0" << std::endl <<
		"blockgap output" << std::endl <<
		"BINARY" << std::endl <<
		"DATASET UNSTRUCTURED_GRID" << std::endl << 
		"POINTS " << ntot << "\tfloat" << std::endl;
		
	// ----------------------------- mesh3D points --------------------------- //

	float x,y,z;

	// spherical valve plate //
	//double xmax = 0.5*in.data.geometry.d_gap_out;
	//double z0 = sqrt(R*R - xmax*xmax);
	// ~sperical valve plate //

	for(int i=0; i<ntot; i++) 
	{
		x = static_cast<float>(mesh.energy.nodes[i].x());
		y = static_cast<float>(mesh.energy.nodes[i].y());
		z = zscale*static_cast<float>(mesh.energy.nodes[i].z());
		//z = sqrt(R*R - (x*x + y*y)) - z0;
		force_big_endian(reinterpret_cast<unsigned char*>(&x));
		vtk.write(reinterpret_cast<char*>(&x),sizeof(float));
		force_big_endian(reinterpret_cast<unsigned char*>(&y));
		vtk.write(reinterpret_cast<char*>(&y),sizeof(float));
		force_big_endian(reinterpret_cast<unsigned char*>(&z));
		vtk.write(reinterpret_cast<char*>(&z),sizeof(float));
	}

	// -------------------------- cells connections -------------------------- //

	vtk << "\nCELLS " << nfluid << "\t" 
		  << (elmVtx + 1)*nfluid << std::endl;
	
	for(int i = 0; i < etot; i++) 
	{
		if(mesh.energy.elements[i].ty == FLUID)
		{
			int tmp = elmVtx;
			force_big_endian(reinterpret_cast<unsigned char*>(&tmp));
			vtk.write(reinterpret_cast<char*>(&tmp),sizeof(int));
			for(int j=0; j<elmVtx; j++) 
			{
				int ndsIdx = mesh.energy.elements[i].nds[j];
				force_big_endian(reinterpret_cast<unsigned char*>(&ndsIdx));
				vtk.write(reinterpret_cast<char*>(&ndsIdx),sizeof(int));	
			}
		}
	}

	// ------------------------------ cells type ----------------------------- //

	vtk << "\nCELL_TYPES " << nfluid << std::endl;
				
	force_big_endian(reinterpret_cast<unsigned char*>(&elmType));

	for(int i = 0; i < etot; i++) 
	{
		if(mesh.energy.elements[i].ty == FLUID)
		{
			int tmp = elmType;
			vtk.write(reinterpret_cast<char*>(&tmp),sizeof(int));
		}
	}	
	
	// --------------------------------- fields ------------------------------ //

	vtk << "\nCELL_DATA " << nfluid << std::endl;

	// --------------------------------- p ----------------------------------- //

	vtk << "SCALARS p[bar] float 1 " << std::endl 
			<< "LOOKUP_TABLE default" << std::endl;
		
	for(int k = 0, id = 0; k < q; k++) 
	{
		for(int id = 0; id<mesh.reynolds.mn; id++)
		{
			if(mesh.energy.elements[id].ty == FLUID)
			{
				float val = static_cast<float>(fields.p.in[id]/1e5);
				force_big_endian(reinterpret_cast<unsigned char*>(&val));
				vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
			}
		}
	}

	// ------------------------------- rho2D --------------------------------- //

	vtk << "SCALARS rho2D[kg/m3] float 1 " << std::endl 
			<< "LOOKUP_TABLE default" << std::endl;
		
	for(int k = 0; k < q; k++) 
	{
		for(int id = 0; id<mesh.reynolds.mn; id++)
		{
			if(mesh.energy.elements[id].ty == FLUID)
			{
				float val = static_cast<float>(fields.rho2D.in[id]);
				force_big_endian(reinterpret_cast<unsigned char*>(&val));
				vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
			}
		}
	}

	// --------------------------------- h ----------------------------------- //

	vtk << "\nSCALARS h[microns] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int k = 0, id = 0; k < q; k++) 
	{
		for(int id = 0; id<mesh.reynolds.mn; id++)
		{
			if(mesh.energy.elements[id].ty == FLUID)
			{
				float val = static_cast<float>(film.h.in[id]*1e6);
				force_big_endian(reinterpret_cast<unsigned char*>(&val));
				vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
			}
		}
	}

	// -------------------------- compenetration ----------------------------- //

	vtk << "\nSCALARS contact[microns] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int k = 0, id = 0; k < q; k++) 
	{
		for(int id = 0; id<mesh.reynolds.mn; id++)
		{
			if(mesh.energy.elements[id].ty == FLUID)
			{
				float val = static_cast<float>(film.contact.in[id]*1e6);
				force_big_endian(reinterpret_cast<unsigned char*>(&val));
				vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
			}
		}
	}

	
	// ------------------------------- dhdt ---------------------------------- //

	vtk << "\nSCALARS dhdt[mm/s] float 1 " << std::endl 
			<< "LOOKUP_TABLE default" << std::endl;
		
	for(int k = 0, id = 0; k < q; k++) 
	{
		for(int id = 0; id<mesh.reynolds.mn; id++)
		{
			if(mesh.energy.elements[id].ty == FLUID)
			{
				float val = static_cast<float>(squeeze.dhdt.in[id]*1e3);
				force_big_endian(reinterpret_cast<unsigned char*>(&val));
				vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
			}
		}
	}

	// -------------------- separate top and bottom half --------------------- //

	vtk << "\nSCALARS top_bottom float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int id = 0; id < mesh.energy.mnq; id++) 
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			int k = (id/mesh.energy.mn);
			float val = 0;
			if(k <= 0.5*mesh.energy.q)
				val = static_cast<float>(0);
			else
				val = static_cast<float>(1);

			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	// ------------------------------- fw ------------------------------------ //

	vtk << "\nSCALARS fw[kg/s] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int id = 0; id < mesh.energy.mnq; id++) 
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(fields.fw.in[id]);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	// ------------------------------- fe ------------------------------------ //

	vtk << "\nSCALARS fe[kg/s] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int id = 0; id < mesh.energy.mnq; id++) 
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(fields.fe.in[id]);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	// ------------------------------- fs ------------------------------------ //

	vtk << "\nSCALARS fs[kg/s] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int id = 0; id < mesh.energy.mnq; id++) 
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(fields.fs.in[id]);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	// ------------------------------- fn ------------------------------------ //

	vtk << "\nSCALARS fn[kg/s] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int id = 0; id < mesh.energy.mnq; id++) 
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(fields.fn.in[id]);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}


	// ------------------------------- mu ------------------------------------ //

	vtk << "\nSCALARS mu[Pas] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int id = 0; id < mesh.energy.mnq; id++) 
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(fields.mu.in[id]);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	// ------------------------------- rho ------------------------------------ //

	vtk << "\nSCALARS rho[kg/m^3] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int id = 0; id < mesh.energy.mnq; id++) 
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(fields.rho.in[id]);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}
	
	// --------------------------------- T ----------------------------------- //

	vtk << "\nSCALARS T[C] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int id = 0; id < mesh.energy.mnq; id++) 
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(fields.T.in[id]);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	// -------------------------------- phid --------------------------------- //

	vtk << "\nSCALARS phid[W] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int id = 0; id < mesh.energy.mnq; id++) 
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(fields.phid.in[id]);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	// -------------------------------- Phid --------------------------------- //

	vtk << "\nSCALARS phid_avg[W] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int id = 0; id < mesh.energy.mnq; id++) 
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(fields.phid_avg.in[id]);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

		
	// -------------------------------- mdotn --------------------------------- //

	vtk << "\nSCALARS mdotn[bar] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int id = 0; id < mesh.energy.mnq; id++) 
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(fields.mdotn.in[id]);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	// -------------------------------- mdots --------------------------------- //

	vtk << "\nSCALARS mdots[bar] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int id = 0; id < mesh.energy.mnq; id++) 
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(fields.mdots.in[id]);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	// -------------------------------- mdotw --------------------------------- //

	vtk << "\nSCALARS mdotwbar] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int id = 0; id < mesh.energy.mnq; id++) 
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(fields.mdotw.in[id]);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	// -------------------------------- mdote --------------------------------- //

	vtk << "\nSCALARS mdote[bar] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int id = 0; id < mesh.energy.mnq; id++) 
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(fields.mdote.in[id]);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}


	// -------------

	// -------------------------------- dp_n --------------------------------- //

	vtk << "\nSCALARS dp_n[bar] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int id = 0; id < mesh.energy.mnq; id++) 
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(fields.g_pn.in[id]);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	// -------------------------------- dp_s --------------------------------- //

	vtk << "\nSCALARS dp_s[bar] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int id = 0; id < mesh.energy.mnq; id++) 
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(fields.g_ps.in[id]);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	// -------------------------------- dp_w --------------------------------- //

	vtk << "\nSCALARS dp_w[bar] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int id = 0; id < mesh.energy.mnq; id++) 
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(fields.g_pw.in[id]);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	// -------------------------------- dp_e --------------------------------- //

	vtk << "\nSCALARS dp_e[bar] float 1 " << std::endl 
		<< "LOOKUP_TABLE default" << std::endl;
		
	for(int id = 0; id < mesh.energy.mnq; id++) 
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			float val = static_cast<float>(fields.g_pe.in[id]);
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}

	// --------------------------------- V ----------------------------------- //
		
	vtk << "\nVECTORS V[m/s] float " << std::endl;
			
	for(int id = 0; id < mesh.energy.mnq; id++) 
	{
		if(mesh.energy.elements[id].ty == FLUID)
		{
			gap_elm el = mesh.energy.elements[id];
			// first component
			float val = 
				static_cast<float>(fields.V.in[id].r()*sin(el.theta) 
				+ fields.V.in[id].theta()*cos(el.theta));
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
			// second component
			val = 
				static_cast<float>(fields.V.in[id].r()*cos(el.theta) 
				- fields.V.in[id].theta()*sin(el.theta));
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
			// third component
			val = static_cast<float>(fields.V.in[id].z());
			force_big_endian(reinterpret_cast<unsigned char*>(&val));
			vtk.write(reinterpret_cast<char*>(&val),sizeof(float));
		}
	}
	
	vtk.close();

	// -------------------- write cylinder block surface -------------------- //

	vtk.clear();

	vtk.open(string(name + string("_block.vtk")).c_str());

	int nelms_block = 0;
	for(int i=0; i<mesh.reynolds.mn; i++)
	{
		if(mesh.reynolds.cb_elms_0[i] == FLUID)
			nelms_block++;
	}

	int nn = mesh.reynolds.nodes.size()/2;
	
	// write VTK header
	vtk << 
		"# vtk DataFile Version 3.0" << std::endl <<
		"blockgap output" << std::endl <<
		"ASCII" << std::endl <<
		"DATASET UNSTRUCTURED_GRID" << std::endl << 
		"POINTS " << nn << "\tfloat" << std::endl;

	for(unsigned int i=0; i<nn; i++)
	{
		vtk << mesh.reynolds.nodes[i].x() << "\t"
				<< mesh.reynolds.nodes[i].y() << "\t"
				<< mesh.reynolds.nodes[i].z() << endl;
	}
	
	vtk << "\nCELLS " << nelms_block << "\t" << (VTK_2D_ELM_NDS + 1)*nelms_block << endl;

	for(int i=0; i<mesh.reynolds.mn; i++)
	{
		if(mesh.reynolds.cb_elms_0[i] == FLUID)
		{
			vtk << VTK_2D_ELM_NDS << "\t";
			for(int j=0; j<VTK_2D_ELM_NDS; j++) 
				vtk <<  mesh.reynolds.elements[i].nds[j] << "\t";
			
			vtk << endl;
		}
	}

	vtk << "\nCELL_TYPES " << nelms_block << endl;

	for(int i=0; i<mesh.reynolds.mn; i++)
	{
		if(mesh.reynolds.cb_elms_0[i] == FLUID)
				vtk << VTK_2D_ELM_TYPE << endl;
	}

	vtk << "\nCELL_DATA " << nelms_block << endl;

	vtk << "SCALARS p[bar] float 1 " << std::endl 
			<< "LOOKUP_TABLE default" << std::endl;
		
	for(int i = 0; i < mesh.reynolds.mn; i++) 
	{
		if(mesh.reynolds.cb_elms_0[i] == FLUID)
			vtk << fields.p.in[i]/1e5 << endl;
	}

	vtk << "SCALARS T[C] float 1 " << std::endl 
			<< "LOOKUP_TABLE default" << std::endl;
		
	for(int i = 0; i < mesh.reynolds.mn; i++) 
	{
		if(mesh.reynolds.cb_elms_0[i] == FLUID)
			vtk << fields.Tcb.in[i] << endl;
	}

	vtk << "SCALARS q[W/m^2] float 1 " << std::endl 
			<< "LOOKUP_TABLE default" << std::endl;
		
	for(int i = 0; i < mesh.reynolds.mn; i++) 
	{
		if(mesh.reynolds.cb_elms_0[i] == FLUID)
			vtk << fields.qcb.in[i] << endl;
	}

	vtk << "SCALARS q_avg[W/m^2] float 1 " << std::endl 
			<< "LOOKUP_TABLE default" << std::endl;
		
	for(int i = 0; i < mesh.reynolds.mn; i++) 
	{
		if(mesh.reynolds.cb_elms_0[i] == FLUID)
			vtk << fields.qcb_avg.in[i] << endl;
	}

	vtk << "SCALARS hcb[microns] float 1 " << endl 
			<< "LOOKUP_TABLE default" << std::endl;

	for(int i = 0; i < mesh.reynolds.mn; i++) 
	{
		if(mesh.reynolds.cb_elms_0[i] == FLUID)
			vtk << film.hcb.in[i]*1e6 << endl;
	}

	vtk << "SCALARS dh_ther[microns] float 1 " << endl 
			<< "LOOKUP_TABLE default" << std::endl;
	
	for(int i = 0; i < mesh.reynolds.mn; i++) 
	{
		if(mesh.reynolds.cb_elms_0[i] == FLUID)
			vtk << film.dhcb_thermal.in[i]*1e6 << endl;
	}

	vtk << "SCALARS dh_EHD[microns] float 1 " << endl 
			<< "LOOKUP_TABLE default" << std::endl;

	for(int i = 0; i < mesh.reynolds.mn; i++) 
	{
		//if(mesh.reynolds.cb_elms_0[i] == FLUID) //Marcos code
		if(mesh.reynolds.cb_elms[i] == FLUID) //removed stationary cylinder block mesh
			vtk << film.dhcb_EHD.in[i]*1e6 << endl;
	}

	vtk.close();

	// ------------------ write the valve plate surface --------------------- //

	vtk.clear();

	vtk.open(string(name + string("_valveplate.vtk")).c_str());

	int nelms_valveplate = 0;
	for(int i=0; i<mesh.reynolds.mn; i++)
	{
		if(mesh.reynolds.vp_elms[i] == FLUID)
			nelms_valveplate++;
	}

	
	// write VTK header
	vtk << 
		"# vtk DataFile Version 3.0" << std::endl <<
		"blockgap output" << std::endl <<
		"ASCII" << std::endl <<
		"DATASET UNSTRUCTURED_GRID" << std::endl << 
		"POINTS " << nn << "\tfloat" << std::endl;

	for(unsigned int i=0; i<nn; i++)
	{
		vtk << mesh.reynolds.nodes[i].x() << "\t"
				<< mesh.reynolds.nodes[i].y() << "\t"
				<< mesh.reynolds.nodes[i].z() << endl;
	}
	
	vtk << "\nCELLS " << nelms_valveplate << "\t" 
			<< (VTK_2D_ELM_NDS + 1)*nelms_valveplate << endl;

	for(int i=0; i<mesh.reynolds.mn; i++)
	{
		if(mesh.reynolds.vp_elms[i] == FLUID)
		{
			vtk << VTK_2D_ELM_NDS << "\t";
			for(int j=0; j<VTK_2D_ELM_NDS; j++) 
				vtk <<  mesh.reynolds.elements[i].nds[j] << "\t";
			
			vtk << endl;
		}
	}

	vtk << "\nCELL_TYPES " << nelms_valveplate << endl;

	for(int i=0; i<mesh.reynolds.mn; i++)
	{
		if(mesh.reynolds.vp_elms[i] == FLUID)
				vtk << VTK_2D_ELM_TYPE << endl;
	}

	vtk << "\nCELL_DATA " << nelms_valveplate << endl;

	vtk << "SCALARS p[bar] float 1 " << std::endl 
			<< "LOOKUP_TABLE default" << std::endl;
		
	for(int i = 0; i < mesh.reynolds.mn; i++) 
	{
		if(mesh.reynolds.vp_elms[i] == FLUID)
			vtk << fields.p.in[i]/1e5 << endl;
	}

	vtk << "SCALARS T[C] float 1 " << std::endl 
			<< "LOOKUP_TABLE default" << std::endl;

	for(int i = 0; i < mesh.reynolds.mn; i++) 
	{
		if(mesh.reynolds.vp_elms[i] == FLUID)
			vtk << fields.Tvp.in[i] << endl;
	}

	vtk << "SCALARS q[W/m^2] float 1 " << std::endl 
			<< "LOOKUP_TABLE default" << std::endl;
		
	for(int i = 0; i < mesh.reynolds.mn; i++) 
	{
		if(mesh.reynolds.vp_elms[i] == FLUID)
			vtk << fields.qvp.in[i] << endl;
	}

	vtk << "SCALARS q_avg[W/m^2] float 1 " << std::endl 
			<< "LOOKUP_TABLE default" << std::endl;
		
	for(int i = 0; i < mesh.reynolds.mn; i++) 
	{
		if(mesh.reynolds.vp_elms[i] == FLUID)
			vtk << fields.qvp_avg.in[i] << endl;
	}

	vtk << "SCALARS hvp[microns] float 1 " << endl 
			<< "LOOKUP_TABLE default" << std::endl;

	for(int i = 0; i < mesh.reynolds.mn; i++) 
	{
		if(mesh.reynolds.vp_elms[i] == FLUID)
			vtk << film.hvp.in[i]*1e6 << endl;
	}

	vtk << "SCALARS dh_ther[microns] float 1 " << endl 
			<< "LOOKUP_TABLE default" << std::endl;
	
	for(int i = 0; i < mesh.reynolds.mn; i++) 
	{
		if(mesh.reynolds.vp_elms[i] == FLUID)
			vtk << film.dhvp_thermal.in[i]*1e6 << endl;
	}

	vtk << "SCALARS dh_EHD[microns] float 1 " << endl 
			<< "LOOKUP_TABLE default" << std::endl;

	for(int i = 0; i < mesh.reynolds.mn; i++) 
	{
		if(mesh.reynolds.vp_elms[i] == FLUID)
			vtk << film.dhvp_EHD.in[i]*1e6 << endl;
	}


	vtk.close();

}
// ------------------------------------------------------------------------- //