# include "../FSTI_Block_dll/log.h"
# include "./block_gap.h"
# include <fstream>
# include "../FSTI_block_thermal/te_solver.h"

extern gaplog Log;

using namespace std;
//
# define pi 3.14159265359

// ------------------------------------------------------------------------- //
cblock_gap::cblock_gap(const input& _in) 
	: in(_in),  gap_main(_in), reynolds(gap_main), energy(gap_main)
{
	tol = in.data.options_block.numeric.epsilonB;
}
// ------------------------------------------------------------------------- //
// calculate the force imbalance (through the gap_main object)
static int calcdF(const gsl_vector* x, void* vp_gap, gsl_vector* dF) 
{
	
	// convert the pointer from void to cblock_gap
	// (void pointer was required by gsl)
	cblock_gap* pgap = static_cast<cblock_gap*>(vp_gap);

	// copy the squeeze velocities from the solver
	pgap->gap_main.squeeze.dhdt1 = gsl_vector_get(x,0);
	pgap->gap_main.squeeze.dhdt2 = gsl_vector_get(x,1);
	pgap->gap_main.squeeze.dhdt3 = gsl_vector_get(x,2);

	
	// calculate the force difference
	pgap->gap_main.set_squeeze();															// set the new squeeze field
	pgap->reynolds.solve_rigid(pgap->gap_main.use_fsi_fb);		// solver reynolds 
	pgap->gap_main.calc_fluid_forces();												// get fluid forces
	pgap->gap_main.calc_dF();																	// get the new force imbalance

	// copy to the solver dF vector
	gsl_vector_set (dF, 0, pgap->gap_main.forces.dF[0]);
	gsl_vector_set (dF, 1, pgap->gap_main.forces.dF[1]);
	gsl_vector_set (dF, 2, pgap->gap_main.forces.dF[2]);
	
	     		
 	return GSL_SUCCESS;
}
// ------------------------------------------------------------------------- //
// calculate the Jacobian
static int calcJ(const gsl_vector* x, void* vp_gap, gsl_matrix* J) 
{

	// convert the pointer from void to cblock_gap
	// (void pointer was required by gsl)
	cblock_gap* pgap = static_cast<cblock_gap*>(vp_gap);

	// small perturbation in the velocity (to calculate the jacobian)
	//double delta = pgap->gap_main.in.data.options_block.numeric.delta_v;
	
	gsl_vector* v = gsl_vector_alloc (3);
	gsl_vector* dF = gsl_vector_alloc (3); // [dF0, dF1, dF2]
	gsl_vector* F0 = gsl_vector_alloc (3); // [F0, F1, F2]

	//       | dF1/dv1 dF1/dv2  dF1/dv3 |
	//       |                          |
	//  J =  | dF2/dv1 dF2/dv2  dF2/dv3 |
	//       |                          |
	//       | dF3/dv1 dF3/dv2  dF3/dv3 |

	Log << "\nCalculating the Jacobian ... ";
	
	// fill the matrix column by column
	for(int j=0; j<3; j++) 
	{
		// define the perturbation
		double delta = fabs(gsl_vector_get(x,j));
		delta = delta > 0 ? 1e-2*delta : 1e-7;

		// add the perturbation
		for(int k = 0; k < 3; k++) 
		{
			if(k == j)
				gsl_vector_set(v,k, gsl_vector_get(x,k) + delta);
			else
				gsl_vector_set(v,k, gsl_vector_get(x,k));
		}
		
		// get the force imbalance with the perturbation
		// if for ex j == 0 here we get dF1(v0+dv,v1,v2) dF2(v0+dv,v1,v2) dF3(v0+dv,v1,v2)
		calcdF(v, pgap, dF);
		
		// remove the perturbation
		for(int k = 0; k < 3; k++) 
		{
			gsl_vector_set(v,k, gsl_vector_get(x,k));
		}
		
		// get the force imbalance without the perturbation
		// here we get dF1(v0,v1,v2) dF2(v0,v1,v2) dF3(v0,v1,v2)
		calcdF(v,pgap,F0);


		// populate the column j of the J matrix
		for(int i=0; i<3; i++) 
		{
			double Jij = (gsl_vector_get(dF,i) - gsl_vector_get(F0,i))/delta;
			gsl_matrix_set( J, i, j, Jij);
		}

	}


	Log << "done!" << gaplog::endl;
	
	return GSL_SUCCESS;
}
// ------------------------------------------------------------------------- //
// call both the dF and J functions
static int calcFJ(const gsl_vector* x, void* vp_gap, gsl_vector* f, gsl_matrix * J)
{
	
	// convert the pointer from void to cblock_gap
	// (void pointer was required by gsl)
	cblock_gap* pgap = static_cast<cblock_gap*>(vp_gap);
	
	calcdF (x, pgap, f);
	calcJ (x, pgap, J);
     
  return GSL_SUCCESS;
}
// ------------------------------------------------------------------------- //
void cblock_gap::get_contact_velocities()
{
	// maximum contact in the gap 
	double cnt_max = gap_main.film.contact.max();
	
	// no contact
	if(cnt_max == 0)
	{
		gap_main.forces.Fc = 0;
		gap_main.forces.Mcx = 0;
		gap_main.forces.Mcy = 0;
		
		gap_main.squeeze.d_dhdt1 = 0;
		gap_main.squeeze.d_dhdt2 = 0;
		gap_main.squeeze.d_dhdt3 = 0;

		return;
	}

	// ------------------------ copy from caspar_input ----------------------- //

	double hmin = in.data.options_block.numeric.hmin;
	double M_block = in.data.geometry.mB;
	double I_block = in.data.geometry.IMB;
	double RBa = 0.50*in.data.geometry.dBa;

	// ----------------------------------------------------------------------- //

	double dt = gap_main.dt;

	double Fc_new = 0;		// contact force
	double Mcx_new = 0;		// contact moment, x direction
	double Mcy_new = 0;		// contact moment, y direction
	
	// average value of fluid cell deformation when 1 Pa p is applied
	double IM_cb_avg = fabs(reynolds.get_IM_cb_avg());
	double IM_vp_avg = fabs(reynolds.get_IM_vp_avg());
	
	// if the info is not available just go with a rough estimation
	IM_cb_avg = (IM_cb_avg > 0) ? IM_cb_avg : 10e-3/207e9; 
	IM_vp_avg = (IM_vp_avg > 0) ? IM_vp_avg : 10e-3/207e9;
	

	int max_cnt_id = 0;
	double max_cnt = -1e-10;
	
	// loop over all elements
	int mn = gap_main.mesh.reynolds.mn;
	for(int e = 0; e < mn; e++)
	{
		gap_elm elm = gap_main.mesh.reynolds.elements[e];
	
		// start to consider contact counting any value lower than hmin
		if(gap_main.mesh.reynolds.elements[e].ty == FLUID)
		{
			double he = gap_main.film.hcb.in[e] - gap_main.film.hvp.in[e];
			if(he < 2*hmin)
			{
				double cnt_e = 2*hmin - he;
				if(cnt_e > max_cnt)
				{
					max_cnt = cnt_e;
					max_cnt_id = e;
				}
			
				double dFc = cnt_e/(IM_cb_avg + IM_vp_avg);
				double dMcx = dFc*elm.r*cos(elm.theta);
				double dMcy = -1.0*dFc*elm.r*sin(elm.theta);
			
				Fc_new += dFc;
				Mcx_new += dMcx; 
				Mcy_new += dMcy;
			}
		}
	
	}

	// assign to new forces to the variables
	gap_main.forces.Fc = Fc_new;
	gap_main.forces.Mcx = Mcx_new;
	gap_main.forces.Mcy = Mcy_new;

	Log	<< "\nContact: Fc = " << gap_main.forces.Fc << " [N]" 
				<< " Mcx = " << gap_main.forces.Mcx << " [Nm]" 
				<< " Mcy = " << gap_main.forces.Mcy << " [Nm]\n" 
				<< gaplog::endl;

	// ---------------- step 2: calculate the block velocity ----------------- //

	// -------------- get the center of instantaneous rotation --------------- //
	
	// distance from the gap surface of the block to the block origin, which
	// is supposed to be the center of instantaneous rotation
	double z0 = gap_main.pistons.z0;

	// ----------------------------------------------------------------------- //

	// F dt = m v

	// axial velocity of the block
	point VC(0, 0, (gap_main.forces.Fc*(dt))/M_block);
	// rotational velocity of the block
	double omega_x = (gap_main.forces.Mcx*(dt))/(I_block);
	double omega_y = (gap_main.forces.Mcy*(dt))/(I_block);
	point omega(omega_x, omega_y, 0);
	
	// center of istantaneous rotation
	point C(0, 0, z0);
	
	// control points
	point P1(0, RBa , 0); // P1
	point P2(cos(pi/6.0)*RBa, -sin(pi/6.0)*RBa, 0); // P2
	point P3(-cos(pi/6.0)*RBa, -sin(pi/6.0)*RBa , 0); // P3

	// velocity at the control points: VC = VP + omega ^ (P - C)
	
	point V1 = 0.5*VC + (omega^(P1 - C));
	point V2 = 0.5*VC + (omega^(P2 - C));
	point V3 = 0.5*VC + (omega^(P3 - C));

	// assign to the contact velocity variables

	double dhdt1 = gap_main.squeeze.dhdt1;
	double dhdt2 = gap_main.squeeze.dhdt2;
	double dhdt3 = gap_main.squeeze.dhdt3;
	double& d_dhdt1 = gap_main.squeeze.d_dhdt1;
	double& d_dhdt2 = gap_main.squeeze.d_dhdt2;
	double& d_dhdt3 = gap_main.squeeze.d_dhdt3;
	
	d_dhdt1 = V1.z();
	d_dhdt2 = V2.z();
	d_dhdt3 = V3.z();

	// check whether the magnitude of the correction is crazy

	double dhdt_sum = (fabs(dhdt1) + fabs(dhdt2) + fabs(dhdt3));
	double d_dhdt_sum = (fabs(d_dhdt1) + fabs(d_dhdt2) + fabs(d_dhdt3));
	
	// saturation the correction if the magnitude of the corection is bigger than
	// the magnitude of the current squeeze motion

	if(d_dhdt_sum > 2*dhdt_sum)
	{
		double R21 = fabs(d_dhdt2)/(fabs(d_dhdt1));
		double R31 = fabs(d_dhdt3)/(fabs(d_dhdt1));
		double d_dhdt1_s = 2*dhdt_sum/(1.0 + R21 + R31);
		double d_dhdt2_s = R21*d_dhdt1_s;
		double d_dhdt3_s = R31*d_dhdt1_s;
		
		d_dhdt1 = d_dhdt1_s*(d_dhdt1/fabs(d_dhdt1));
		d_dhdt2 = d_dhdt2_s*(d_dhdt2/fabs(d_dhdt2));
		d_dhdt3 = d_dhdt3_s*(d_dhdt3/fabs(d_dhdt3));

	}
}
/*
void cblock_gap::get_contact_forces()
{
	// maximum contact in the gap 
	double cnt_max = gap_main.film.contact.max();
	
	// no contact
	if(cnt_max == 0)
	{
		gap_main.forces.Fc = 0;
		gap_main.forces.Mcx = 0;
		gap_main.forces.Mcy = 0;
		
		gap_main.squeeze.d_dhdt1 = 0;
		gap_main.squeeze.d_dhdt2 = 0;
		gap_main.squeeze.d_dhdt3 = 0;

		return;
	}

	// ------------------------ copy from caspar_input ----------------------- //

	double hmin = in.data.options_block.numeric.hmin;
	double M_block = in.data.geometry.mB;
	double I_block = in.data.geometry.IMB;
	double RBa = 0.50*in.data.geometry.dBa;

	// ----------------------------------------------------------------------- //

	double dt = gap_main.dt;

	double Fc_new = 0;		// contact force
	double Mcx_new = 0;		// contact moment, x direction
	double Mcy_new = 0;		// contact moment, y direction
	
	// average value of fluid cell deformation when 1 Pa p is applied
	double IM_cb_avg = fabs(reynolds.get_IM_cb_avg());
	double IM_vp_avg = fabs(reynolds.get_IM_vp_avg());
	
	// if the info is not available just go with a rough estimation
	IM_cb_avg = (IM_cb_avg > 0) ? IM_cb_avg : 10e-3/207e9; 
	IM_vp_avg = (IM_vp_avg > 0) ? IM_vp_avg : 10e-3/207e9;
	

	int max_cnt_id = 0;
	double max_cnt = -1e-10;
	
	// loop over all elements
	int mn = gap_main.mesh.reynolds.mn;
	for(int e = 0; e < mn; e++)
	{
		gap_elm elm = gap_main.mesh.reynolds.elements[e];
	
		// start to consider contact counting any value lower than hmin
		if(gap_main.mesh.reynolds.elements[e].ty == FLUID)
		{
			double he = gap_main.film.hcb.in[e] - gap_main.film.hvp.in[e];
			if(he < 2*hmin)
			{
				double cnt_e = 2*hmin - he;
				if(cnt_e > max_cnt)
				{
					max_cnt = cnt_e;
					max_cnt_id = e;
				}
			
				double dFc = cnt_e/(IM_cb_avg + IM_vp_avg);
				double dMcx = dFc*elm.r*cos(elm.theta);
				double dMcy = -1.0*dFc*elm.r*sin(elm.theta);
			
				Fc_new += dFc;
				Mcx_new += dMcx; 
				Mcy_new += dMcy;
			}
		}
	
	}

	// assign to new forces to the variables
	gap_main.forces.Fc = Fc_new;
	gap_main.forces.Mcx = Mcx_new;
	gap_main.forces.Mcy = Mcy_new;

	Log	<< "\nContact: Fc = " << gap_main.forces.Fc << " [N]" 
				<< " Mcx = " << gap_main.forces.Mcx << " [Nm]" 
				<< " Mcy = " << gap_main.forces.Mcy << " [Nm]\n" 
				<< gaplog::endl;

	// ---------------- step 2: calculate the block velocity ----------------- //

	// -------------- get the center of instantaneous rotation --------------- //
	
	// distance from the gap surface of the block to the block origin, which
	// is supposed to be the center of instantaneous rotation
	double z0 = gap_main.pistons.z0;

	// ----------------------------------------------------------------------- //

	// F dt = m v

	// axial velocity of the block
	point VC(0, 0, (gap_main.forces.Fc*(dt))/M_block);
	// rotational velocity of the block
	double omega_x = (gap_main.forces.Mcx*(dt))/(I_block);
	double omega_y = (gap_main.forces.Mcy*(dt))/(I_block);
	point omega(omega_x, omega_y, 0);
	
	// center of istantaneous rotation
	point C(0, 0, z0);
	
	// control points
	point P1(0, RBa , 0); // P1
	point P2(cos(pi/6.0)*RBa, -sin(pi/6.0)*RBa, 0); // P2
	point P3(-cos(pi/6.0)*RBa, -sin(pi/6.0)*RBa , 0); // P3

	// velocity at the control points: VC = VP + omega ^ (P - C)
	
	point V1 = 0.5*VC + (omega^(P1 - C));
	point V2 = 0.5*VC + (omega^(P2 - C));
	point V3 = 0.5*VC + (omega^(P3 - C));

	// assign to the contact velocity variables

	double dhdt1 = gap_main.squeeze.dhdt1;
	double dhdt2 = gap_main.squeeze.dhdt2;
	double dhdt3 = gap_main.squeeze.dhdt3;
	double& d_dhdt1 = gap_main.squeeze.d_dhdt1;
	double& d_dhdt2 = gap_main.squeeze.d_dhdt2;
	double& d_dhdt3 = gap_main.squeeze.d_dhdt3;
	
	d_dhdt1 = V1.z();
	d_dhdt2 = V2.z();
	d_dhdt3 = V3.z();

	// check whether the magnitude of the correction is crazy

	double dhdt_sum = (fabs(dhdt1) + fabs(dhdt2) + fabs(dhdt3));
	double d_dhdt_sum = (fabs(d_dhdt1) + fabs(d_dhdt2) + fabs(d_dhdt3));
	
	// saturation the correction if the magnitude of the corection is bigger than
	// the magnitude of the current squeeze motion

	if(d_dhdt_sum > 2*dhdt_sum)
	{
		double R21 = fabs(d_dhdt2)/(fabs(d_dhdt1));
		double R31 = fabs(d_dhdt3)/(fabs(d_dhdt1));
		double d_dhdt1_s = 2*dhdt_sum/(1.0 + R21 + R31);
		double d_dhdt2_s = R21*d_dhdt1_s;
		double d_dhdt3_s = R31*d_dhdt1_s;
		
		d_dhdt1 = d_dhdt1_s*(d_dhdt1/fabs(d_dhdt1));
		d_dhdt2 = d_dhdt2_s*(d_dhdt2/fabs(d_dhdt2));
		d_dhdt3 = d_dhdt3_s*(d_dhdt3/fabs(d_dhdt3));

	}
}

*/

// ------------------------------------------------------------------------- //
int cblock_gap::get_block_balance(bool reset)
{
	
	int status;					// solver status
	size_t iter = 0;			// iteration
	const size_t n = 3;			// system dimension

	// solver function
	gsl_multiroot_function_fdf f = 
	{
		&calcdF,
	  &calcJ,
	  &calcFJ,
	  n, 
		this
	};

	// allocate solution vector
	gsl_vector *x = gsl_vector_alloc (n);

	// initialize the solution vector with the current squeeze velocities
	if(reset)
	{
		gsl_vector_set (x, 0, 0);
		gsl_vector_set (x, 1, 0);
		gsl_vector_set (x, 2, 0);
	}
	else
	{
		gsl_vector_set (x, 0, gap_main.squeeze.dhdt1);
		gsl_vector_set (x, 1, gap_main.squeeze.dhdt2);
		gsl_vector_set (x, 2, gap_main.squeeze.dhdt3);
	}

	// define and allocate the solver
	T = gsl_multiroot_fdfsolver_hybridsj;
	s = gsl_multiroot_fdfsolver_alloc (T, n);
	gsl_multiroot_fdfsolver_set (s, &f, x);


	Log << "\nIteration: " << iter << "\n" << gaplog::endl;
	
	std::streamsize default_val = cout.precision();
	cout.precision(3);
	cout << showpos;

	Log	<< "  dhdt[1] = " << std::scientific << gsl_vector_get (s -> x, 0) 
			<< "  dhdt[2] = " << std::scientific << gsl_vector_get (s -> x, 1) 
			<< "  dhdt[3] = " << std::scientific << gsl_vector_get (s -> x, 2) 
			<< gaplog::endl;

	Log	<< "  dF[1]   = " << std::scientific << gsl_vector_get (s -> f, 0) 
			<< "  dF[2]   = " << std::scientific << gsl_vector_get (s -> f, 1) 
			<< "  dF[3]   = " << std::scientific << gsl_vector_get (s -> f, 2) 
			<< gaplog::endl;

	// --------------------- main force balance loop ------------------------- //
	
	do 
	{
		
		iter++;

		status = gsl_multiroot_fdfsolver_iterate (s);

		Log << "\nIteration: " << iter << "\n" << gaplog::endl;
	
		Log	<< "  dhdt[1] = " << std::scientific << gsl_vector_get (s -> x, 0) 
				<< "  dhdt[2] = " << std::scientific << gsl_vector_get (s -> x, 1) 
				<< "  dhdt[3] = " << std::scientific << gsl_vector_get (s -> x, 2) 
				<< gaplog::endl;

		Log	<< "  dF[1]   = " << std::scientific << gsl_vector_get (s -> f, 0) 
				<< "  dF[2]   = " << std::scientific << gsl_vector_get (s -> f, 1) 
				<< "  dF[3]   = " << std::scientific << gsl_vector_get (s -> f, 2) 
				<< gaplog::endl;

		if(status)
			break;

		status = gsl_multiroot_test_residual (s -> f, tol);
	
	} while (status == GSL_CONTINUE && iter < 100);

	// ------------------------------------------------------------------------ //

	// store the result in the structure
	gap_main.squeeze.dhdt1 = gsl_vector_get(s -> x, 0);
	gap_main.squeeze.dhdt2 = gsl_vector_get(s -> x, 1);
	gap_main.squeeze.dhdt3 = gsl_vector_get(s -> x, 2);


	std::cout.precision(default_val);
	Log << noshowpos;

	Log << "\nStatus = " << gsl_strerror (status) << gaplog::endl;

	gsl_multiroot_fdfsolver_free (s);
	gsl_vector_free (x);


	return status;
}
// ------------------------------------------------------------------------- //
void cblock_gap::check_gap()
{
	Log << "\n* ------------------------ Check of the gap definition ---------------------- *\n" 
			<< gaplog::endl;

	if(in.data.lubrication_module.solve_block < 1)
	{
		Log << "\n\nPlease enable \"solve_block\" in lubricating_module\n\n";
		exit(1);
	}

	double cb_surf = 0;
	double vp_surf = 0;

	for(unsigned int i=0; i<gap_main.mesh.reynolds.cb_gap_surf.ai.size(); i++)
		cb_surf += gap_main.mesh.reynolds.cb_gap_surf.ai[i];
	for(unsigned int i=0; i<gap_main.mesh.reynolds.vp_gap_surf.ai.size(); i++)
		vp_surf += gap_main.mesh.reynolds.vp_gap_surf.ai[i];

	Log << "Cylinder block gap surface area: " << fixed << 1e6*cb_surf
			<< " [mm2]" << gaplog::endl;
	Log << "Valve plate gap surface area:    " << fixed << 1e6*vp_surf
			 << " [mm2]" << gaplog::endl;
	
	Log << "\nUpdating mesh ... ";
	gap_main.set_film_thickness();
	gap_main.update_mesh();
	gap_main.mesh.energy.rotate(0);
	gap_main.mesh.reynolds.rotate(0);
	Log << "done!" << gaplog::endl;	

	Log << "\nSolving Reynolds and Energy equation ... ";
	
	
	Log << "\nSetting boundary conditions ... ";
	gap_main.set_pressure();
	gap_main.set_temperature();
	Log << "done!" << gaplog::endl;	

	Log << "\nSolving Reynolds equation ... ";
	reynolds.solve_rigid();
	Log << "done!" << gaplog::endl;	

	Log << "\nCalculating mass flow ... ";
	gap_main.calc_mass_flow();
	Log << "done!" << gaplog::endl;	

	Log << "\nSolving Energy equation ... ";
	energy.solve();
	Log << "done!" << gaplog::endl;	

	Log << "\nUpdating viscosity ... ";
	gap_main.update_viscosity();
	Log << "done!" << gaplog::endl;	

	Log << "\nCalculating fluid forces ... ";
	gap_main.calc_fluid_forces();
	Log << "done!" << gaplog::endl;	


	// ------------------------- write output  ----------------------------- //

	system("IF NOT EXIST output\\block\\debug\\gap (mkdir output\\block\\debug\\gap)");

	Log << "\nWriting ... \n\n" 
			<< "   ./output/block/debug/gap/cb_surf.vtk\n" 
			<< "   ./output/block/debug/gap/test.vtk\n" 
			<< "   ./output/block/debug/gap/test_block.vtk\n" 
			<< "   ./output/block/debug/gap/test_valveplate.vtk\n" 
			<< gaplog::endl;
	
	gap_main.write_pcb_vtk("./output/block/debug/gap/cb_surf.vtk");
	gap_main.write_vtk("./output/block/debug/gap/test.vtk", 1000);
	
	Log << "done!\n" << gaplog::endl;

	gap_main.update_pistons(0);
	gap_main.calc_piston_forces();

}
// ------------------------------------------------------------------------- //
void cblock_gap::check_thermal()
{
	
	Log << "\n* -------------------------- Check Thermal calculation ------------------------ *\n" 
			<< gaplog::endl;

	if(in.data.lubrication_module.solve_block < 1)
	{
		Log << "\n\nPlease enable \"solve_block\" in lubricating_module\n\n";
		exit(1);
	}

	if(in.data.options_block.general.Thermal_CB)
		Log << "\n*  Cylinder block Thermal option: ON" << gaplog::endl;
	else
		Log << "\n*  Cylinder block Thermal option: OFF" << gaplog::endl;

	if(in.data.options_block.general.Thermal_VP)
		Log << "\n*  Valve plate Thermal option: ON" << gaplog::endl;
	else
		Log << "\n*  Valve plate Thermal option: OFF" << gaplog::endl;

	if(in.data.options_block.general.Thermal_CB || in.data.options_block.general.Thermal_VP)
	{
		system("IF NOT EXIST output\\block\\debug\\thermal (mkdir output\\block\\debug\\thermal)");

		Log << "\nSolving Reynolds and Energy equation ... ";
		gap_main.set_film_thickness();
		gap_main.update_mesh();
		gap_main.set_pressure();
		gap_main.set_temperature();
		gap_main.update_viscosity();
		reynolds.solve_rigid();
		gap_main.calc_mass_flow();
		gap_main.update_h_previous();
		energy.solve();
		Log << "done!" << gaplog::endl;	
		Log << "\nCalculating heat fluxes ... ";
		gap_main.fields.qcb_avg = gap_main.fields.qcb_prog = gap_main.calc_cb_heatflux();
		gap_main.fields.qvp_avg = gap_main.fields.qvp_prog = gap_main.calc_vp_heatflux();
		Log << "done!" << gaplog::endl;
	}
	else
	{
		Log << "\nNo thermal option activated in the options_block.txt file\n" 
				<< gaplog::endl;
		return;
	}

	// test cylinder block thermal
	if(in.data.options_block.general.Thermal_CB)
	{
		Log << "\n\nSolving for cylinder block thermo-elastic analysis" 
				<< gaplog::endl;

		// thermoelastic solver
		te_CB = new te_solver(in, "CB");

		// pointer to the interpolation grid
		interpl_grid* gap_str = &(te_CB->gap);

		// apply the heat flux from gap
		gap_main.apply_heat_flux(gap_str, "CB");	
		// solve for heat transfer
		te_CB->solve_htr();
		// interpolation of temperature from structure to fluid
		interpolation str_to_fl;
		str_to_fl.cellsTocells(gap_str, &gap_main.CBf);
		// update the surface temperature
		gap_main.update_Tcb();	

		// solve for thermal deflection
		te_CB->solve_thdfl();
		// interpolation of deformation from structure to fluid
		str_to_fl.pointsTocells(gap_str, &gap_main.CBf);
		// update the gap deformation
		gap_main.update_dhcb_thermal();

		Log << "\nWriting ./output/block/debug/thermal/cylinderblock.vtk ..." ;
		te_CB->write_vtk("./output/block/debug/thermal");
		/*
		te_CB->write_fset_vtk("./output/block/debug/thermal", "gap_block");
		te_CB->write_fset_vtk("./output/block/debug/thermal", "DCs");
		//te_CB->write_fset_vtk("./output/block/debug/thermal", "dc");
		te_CB->write_fset_vtk("./output/block/debug/thermal", "case_outer");
		te_CB->write_fset_vtk("./output/block/debug/thermal", "case_inner");
		*/
		for(unsigned int i = 0; i < 9; i++) 
		{
			ostringstream oss;
			oss << "gap_piston_" << i+1 ;
			te_CB->write_fset_vtk("./output/block/debug/thermal", oss.str().c_str());
		}
		Log << "done!" << gaplog::endl;

		delete te_CB;

	}

	// test valve plate thermal
	if(in.data.options_block.general.Thermal_VP)
	{
		Log << "\n\nSolving for valve plate thermo-elastic analysis" 
				<< gaplog::endl;

		// thermoelastic solver
		te_VP = new te_solver(in, "VP");

		// pointer to the interpolation grid
		interpl_grid* gap_str = &(te_VP->gap);

		// apply the heat flux from gap
		gap_main.apply_heat_flux(gap_str, "VP");	
		// solve for heat transfer
		te_VP->solve_htr();
		// interpolation of temperature from structure to fluid
		interpolation str_to_fl;
		str_to_fl.cellsTocells(gap_str, &gap_main.VPf);
		// update the surface temperature
		gap_main.update_Tvp();	

		// solve for thermal deflection
		te_VP->solve_thdfl();
		// interpolation of deformation from structure to fluid
		str_to_fl.pointsTocells(gap_str, &gap_main.VPf);
		// update the gap deformation
		gap_main.update_dhvp_thermal();

		Log << "\nWriting ./output/block/debug/thermal/valveplate.vtk ..." ;
		te_VP->write_vtk("./output/block/debug/thermal");
		/*
		te_VP->write_fset_vtk("./output/block/debug/thermal", "gap");
		te_VP->write_fset_vtk("./output/block/debug/thermal", "top_const");
		te_VP->write_fset_vtk("./output/block/debug/thermal", "bottom_const");
		te_VP->write_fset_vtk("./output/block/debug/thermal", "bottom_const2");
		te_VP->write_fset_vtk("./output/block/debug/thermal", "HP");
		te_VP->write_fset_vtk("./output/block/debug/thermal", "LP");
		te_VP->write_fset_vtk("./output/block/debug/thermal", "case_inner");
		te_VP->write_fset_vtk("./output/block/debug/thermal", "case_outer");
		te_VP->write_fset_vtk("./output/block/debug/thermal", "air");
		*/
		Log << "done!" << gaplog::endl;

		delete te_VP;
	}

	Log << "\nWriting ./output/block/debug/thermal/fluid.vtk ..." ;
	gap_main.write_vtk("./output/block/debug/thermal/fluid.vtk", 1000);
	Log << "done!" << gaplog::endl;

}

// ------------------------------------------------------------------------- //
void cblock_gap::check_EHD()
{
	Log << "\n* -------------------------- Check EHD calculation ------------------------ *\n" 
			<< gaplog::endl;

	if(in.data.lubrication_module.solve_block < 1)
	{
		Log << "\n\nPlease enable \"solve_block\" in lubricating_module\n\n";
		exit(1);
	}

	if(in.data.options_block.general.EHD_CB)
		Log << "\n*  Cylinder block EHD option: ON" << gaplog::endl;
	else
		Log << "\n*  Cylinder block EHD option: OFF" << gaplog::endl;

	if(in.data.options_block.general.EHD_VP)
		Log << "\n*  Valve plate EHD option: ON" << gaplog::endl;
	else
		Log << "\n*  Valve plate EHD option: OFF" << gaplog::endl;

	// load influece matrices
	if(in.data.options_block.general.EHD_CB || in.data.options_block.general.EHD_VP)
	{
		Log << "\nReading influence matrices ... " << gaplog::endl;
		reynolds.load_matrices();
		Log << "... matrices loaded succesfully!" << gaplog::endl;
	}
	else
	{
		Log << "\nNo EHD option activated in the options_block.txt file\n" 
				<< gaplog::endl;
		return;
	}

	Log << "\n* Rigid thickness set with: "
			<< "(h1 = " << gap_main.film.h1*1e6  
			 << "; h2 = " << gap_main.film.h2*1e6  
			 << "; h3 = " << gap_main.film.h3*1e6  
			 << ") [um]\n"  
			 << gaplog::endl;

	Log << "* Low pressure: " << gap_main.pfile.getp_deg(0)[1]/1e5 << " [bar]\n" 
			 << gaplog::endl;

	Log << "* High pressure: " << gap_main.pfile.getp_deg(0)[2]/1e5 << " [bar]\n" 
			<< gaplog::endl;

	// solve EHD
	gap_main.set_film_thickness();
	gap_main.update_mesh();
	gap_main.set_pressure();
	gap_main.set_temperature();
	gap_main.update_viscosity();
	reynolds.solve_EHD();
	gap_main.calc_mass_flow();

	// ------------------------- write output  ----------------------------- //

	system("IF NOT EXIST output\\block\\debug\\EHD (mkdir output\\block\\debug\\EHD)");

	Log << "\nWriting ... \n" 
			<< gaplog::endl
			<< "   ./output/block/debug/EHD/test.vtk" 
			<< gaplog::endl
			<< "   ./output/block/debug/EHD/test_block.vtk" 
			<< gaplog::endl
			<< "   ./output/block/debug/EHD/test_valveplate.vtk\n" 
			<< gaplog::endl;
	
	gap_main.update_mesh();
	gap_main.write_vtk("./output/block/debug/EHD/test", 1000);
	
	Log << "done!\n" << gaplog::endl;

	Log << "\nWriting interpolation grids ... \n" << gaplog::endl;

	// block
	if(in.data.options_block.general.EHD_CB)
	{
		Log << "   ./output/block/debug/EHD/CBs.vtk" << gaplog::endl; 
		Log << "   ./output/block/debug/EHD/CBf.vtk" << gaplog::endl; 
		gap_main.CBf.writeVTK("./output/block/debug/EHD/CBf.vtk");
		reynolds.get_CBs().writeVTK("./output/block/debug/EHD/CBs.vtk");
	}

	// valve plate
	if(in.data.options_block.general.EHD_VP)
	{
		Log << "   ./output/debug/EHD/VPs.vtk" << gaplog::endl; 
		Log << "   ./output/debug/EHD/VPf.vtk" << gaplog::endl; 
		gap_main.VPf.writeVTK("./output/block/debug/EHD/VPf.vtk");
		reynolds.get_VPs().writeVTK("./output/block/debug/EHD/VPs.vtk");
	}

	Log << "\ndone!\n" << gaplog::endl;
}
// ------------------------------------------------------------------------- //
void cblock_gap::write_cb_influgen(const char* file)
{
	
	// ------------------------ copy from caspar_input ----------------------- //
	
	double Fspring = in.data.geometry.Fblock;
	int z = in.data.operating_conditions.npistons;

	// ----------------------------------------------------------------------- //
		
	ofstream out(file);

	if(!out.is_open())
	{
		Log << "\nCould not open " << file << " for writing" << gaplog::endl;
		exit(1);
	}

	Log << "\nWriting " << file << "\n" << gaplog::endl;

	Log << "  *\tWriting block position ... ";
	out << endl << "// Block position" << endl;
	out << "phi\t" << gap_main.phi << endl;
	out << "dz\t" << gap_main.film.dz << endl;
	out << "rx\t" << gap_main.film.rx << endl;
	out << "ry\t" << gap_main.film.ry << endl;
	Log << "done!" << gaplog::endl;

	Log << "  *\tWriting deformation scaling factor ... ";
	out << "// Deformation scaling factors" << endl;
	out << "xscale\t" << 1.0 << endl;
	out << "yscale\t" << 1.0 << endl;
	out << "zscale\t" << 1000.0 << endl;
	Log << "done!" << gaplog::endl;

	// define the interpolation grid
	if(gap_main.EHD_CB)
	{

		Log << "  *\tWriting spring force ... ";
		out << endl << "// Spring force" << endl;
		out << "Fspring\t" << Fspring << endl;
		Log << "done!" << gaplog::endl;
		
		Log << "  *\tWriting displacement chambers pressures ... ";
		out << endl << "// Displacement chamber pressures" << endl;
		out << "pDC\t" << z << endl;
		for(int i=0; i<z; i++)
			out << gap_main.pistons.pDC[i] << endl;
		Log << "done!" << gaplog::endl;

		Log << "  *\tWriting spring force ... ";
		out << endl << "// Piston loads" << endl;
		out << "MKBx\t" << gap_main.forces.MKBx << endl;
		out << "MKBy\t" << gap_main.forces.MKBy << endl;
		Log << "done!" << gaplog::endl;	

		
		interpl_grid gap;
		
		// work with the internal structure grid
		Log << "  *\tGetting CB interpolation grid ...";
		const interpl_grid* tmp;
		tmp = &reynolds.get_CBs();
		tmp->copy(gap);
		Log << "done!" << gaplog::endl;
		
		// define the interpolation object and interpolate from fluid to structure
		Log << "  *\tInterpolating current pressure field to structure grid ... ";
		interpolation I;
		gap_main.interpl_pressure(); // --> update CBf
		I.cellsTocells(&gap_main.CBf, &gap);
		Log << "done!" << gaplog::endl;

		Log << "  *\tWriting gap grid and pressure field ...";
		out << endl << "// Gap grid and pressure" << endl;
		out << "NODES " << gap.nodes.size() << endl;
		for(unsigned int i = 0; i < gap.nodes.size(); i++) 
		{
			out << gap.nodes[i].x() << "\t" 
					<< gap.nodes[i].y() << "\t" 
					<< gap.nodes[i].z() << endl; 
		}
		out << "ELEMENTS " << gap.elements.size() << endl;
		for(int i = 0; i < gap.ne(); i++) 
		{
			for(int j = 0; j < gap.ELM_NDS; j++)
				out << gap.elements[i][j] << "\t"; // index start from 0
			out << endl;
		}
		out << "ACTIVES " << gap.elements.size() << endl;
		for(int i = 0; i < gap.ne(); i++) 
		{
			for(int j = 0; j < gap.ELM_NDS; j++)
				out << gap.active[i] << endl;
		}
		out << "CELLS_DATA" << endl;
		for(int i = 0; i < gap.ne(); i++) 
		{
			out << gap.cells_data[i] << endl;
		}
		out << "POINTS_DATA" << endl;
		for(int i = 0; i < gap.nn(); i++) 
		{
			out << gap.points_data[i] << endl;
		}
		
		Log << "done!" << gaplog::endl;

	}

	Log << "\ndone writing " << file << gaplog::endl;

	out.close();
	
}
// ------------------------------------------------------------------------- //
void cblock_gap::write_vp_influgen(const char* file)
{
	// ------------------------ copy from caspar_input ----------------------- //
	
	double LP = in.data.operating_conditions.LP;
	double HP = in.data.operating_conditions.HP;

	// ----------------------------------------------------------------------- //


	ofstream out(file);

	if(!out.is_open())
	{
		Log << "\nCould not open " << file << " for writing" << gaplog::endl;
		exit(1);
	}

	Log << "\nWriting " << file << "\n" << gaplog::endl;

	Log << "  *\tWriting deformation scaling factor ... ";
	out << "// Deformation scaling factors" << endl;
	out << "xscale\t" << 1.0 << endl;
	out << "yscale\t" << 1.0 << endl;
	out << "zscale\t" << 1000.0 << endl;
	Log << "done!" << gaplog::endl;

	// define the interpolation grid

	if(gap_main.EHD_VP)
	{

		interpl_grid gap;

		Log << "  *\tWriting ports pressures ... ";
		out << endl << "// Ports pressures" << endl;
		out << "LP\t" << LP << endl;
		out << "HP\t" << HP << endl;
		Log << "done!" << gaplog::endl;

		// work with the internal structure grid
		Log << "  *\tGetting VB interpolation grid ...";
		const interpl_grid* tmp;
		tmp = &reynolds.get_VPs();
		tmp->copy(gap);
		Log << "done!" << gaplog::endl;		

		// define the interpolation object and interpolate from fluid to structure
		Log << "  *\tInterpolating current pressure field to structure grid ... ";
		interpolation I;
		gap_main.interpl_pressure(); // --> update VPf
		I.cellsTocells(&gap_main.VPf, &gap);
		Log << "done!" << gaplog::endl;

		Log << "  *\tWriting gap grid and pressure field ...";
		out << endl << "// Gap grid and pressure" << endl;
		out << "NODES " << gap.nodes.size() << endl;
		for(unsigned int i = 0; i < gap.nodes.size(); i++) 
		{
			out << gap.nodes[i].x() << "\t" 
					<< gap.nodes[i].y() << "\t" 
					<< gap.nodes[i].z() << endl; 
		}
		out << "ELEMENTS " << gap.elements.size() << endl;
		for(int i = 0; i < gap.ne(); i++) 
		{
			for(int j = 0; j < gap.ELM_NDS; j++)
				out << gap.elements[i][j] << "\t"; // index start from 0
			out << endl;
		}
		out << "ACTIVES " << gap.elements.size() << endl;
		for(int i = 0; i < gap.ne(); i++) 
		{
			for(int j = 0; j < gap.ELM_NDS; j++)
				out << gap.active[i] << endl;
		}
		out << "CELLS_DATA" << endl;
		for(int i = 0; i < gap.ne(); i++) 
		{
			out << gap.cells_data[i] << endl;
		}
		out << "POINTS_DATA" << endl;
		for(int i = 0; i < gap.nn(); i++) 
		{
			out << gap.points_data[i] << endl;
		}
		
		Log << "done!" << gaplog::endl;

	}

	Log << "\ndone writing " << file << gaplog::endl;

	out.close();

}
// ------------------------------------------------------------------------- //
void cblock_gap::writeInflugenInputs()
{
	
	if(in.data.options_block.general.EHD_CB*in.data.options_block.general.EHD_VP == 0)
	{
		Log << "\nEHD_CB and EHD_VP options have to be both activated\n" << gaplog::endl;
		return;
	}
	
	reynolds.load_matrices();
	gap_main.define_interpl_grids();
	gap_main.set_pressure();
	gap_main.set_film_thickness();
	gap_main.calc_external_forces();
	get_block_balance();
	
	reynolds.solve_rigid();
	gap_main.interpl_pressure();
	gap_main.VPf.write("VPf.txt");
	gap_main.VPf.writeVTK("VPf.vtk");
	gap_main.CBf.write("CBf.txt");
	gap_main.CBf.writeVTK("CBf.vtk");
	reynolds.get_VPs().write("VPS.txt");
	reynolds.get_VPs().writeVTK("VPS.vtk");
	reynolds.get_CBs().write("CBS.txt");
	reynolds.get_CBs().writeVTK("CBS.vtk");
	write_cb_influgen("CBInflugen.txt");
	write_vp_influgen("VPInflugen.txt");

	Log << "Written:\n\n";
	Log << " * VPf.txt\n";
	Log << " * VPf.vtk\n";
	Log << " * CBf.vtk\n";
	Log << " * CBf.txt\n";
	Log << " * VPS.txt\n";
	Log << " * VPS.vtk\n";
	Log << " * CBS.txt\n";
	Log << " * CBS.vtk\n";
	Log << " * VPInflugen.txt\n";
	Log << " * CBInflugen.vtk\n";

	
}
// ------------------------------------------------------------------------- //
void cblock_gap::get_piston_kinematics()
{
		system("IF NOT EXIST output\\block\\debug (mkdir output\\block\\debug)");
	
	ofstream out("./output/block/debug/piston_position.txt");
	
	// phi angle
	// sk piston position respect ODP
	// sk_O piston position respect block origin
	
	out << "phi\tsK(ODP)\tsK(Block Origin)" << endl;

	Log << "\n* --------------------- Calculation of piston kinematics -------------------- *\n" 
			<< gaplog::endl;

	Log << " * Total piston stroke = " << 1e3*gap_main.pistons.Hk << " [mm]" << gaplog::endl;
		
	if(gap_main.in.data.geometry.d_spherical > 0)
	{
		Log << " * At ODP: derived distance from the the piston end to the point\n"
				<< "   defined by the intersection of the sherical surface with the\n"
				<< "   shaft axis = " 
				<< 1e3*(in.data.geometry.lZ0 + in.data.geometry.lengthcanalB)
				<< " [mm]" << gaplog::endl;
		
		Log << " * Distance from the cylinder block's reference system to the point defined\n"
				<< "   by the intersection of the spherical surface with the shaft axis = " 
				<< 1e3*gap_main.pistons.z0  << " [mm]" << gaplog::endl;
	}
	else
	{
		Log << " * At ODP: derived distance from the the piston end to the \n"
				<< "   cylinder block's sealing land = "
				<< 1e3*(in.data.geometry.lZ0 + in.data.geometry.lengthcanalB)
				<< " [mm]" << gaplog::endl;
		Log << " * Distance from the cylinder block's reference system to cylinder block's\n"
				<< "   sealing land = " << 1e3*gap_main.pistons.z0  << " [mm]" << gaplog::endl;
	}

	
	Log << "\nCalculating piston position over one shaft revolution ... " << gaplog::endl;

	
	for(unsigned int i=0; i<360; i++)
	{
		double phi = i*pi/180.0;
		gap_main.rotate(phi);
		out << scientific << phi << "\t" 
				<< scientific << gap_main.pistons.sK[0] << "\t" 
				<< scientific << gap_main.pistons.Hk + gap_main.pistons.sK[0]  
				<< endl;
	}

	out.close();

	Log << "\nWritten file ./output/block/debug/piston_position.txt" << gaplog::endl;

	
	
}
// ------------------------------------------------------------------------- //
void cblock_gap::get_external_loads()
{
	system("IF NOT EXIST output\\block\\debug\\loads (mkdir output\\block\\debug\\loads)");
	
	ofstream out_total("./output/block/debug/loads/total_loads.txt");
	ofstream out_pistons("./output/block/debug/loads/loads_from_pistons.txt");

	out_total << "phi\tFBz\txFBz\tyFBz\tMBx\tMBy" << endl;
	out_pistons << "phi\tFKBx\tFKBy\tMKBx\tMKBy" << endl;

	Log << "\n* ----------- Calculation of the external loads on cylinder block ----------- *\n" 
			<< gaplog::endl;

	if(gap_main.in.data.geometry.d_spherical > 0)
	{
		Log << "\nAxial distance between the cylinder block's reference system and the \n"
				<< "intersection point between the shaft axis and the spherical gap surface\n"
				<< "is " << 1e3*gap_main.pistons.z0 << " [mm]" << gaplog::endl;
	}
	else
	{
		Log << "\nAxial distance between the cylinder block's reference system and the \n"
				<< "sealing land surface is " << 1e3*gap_main.pistons.z0 << " [mm]" << gaplog::endl;
	}

	Log << "\nDistance between the point of application of the spline's reaction\n"
			<< "force and the cylinder block's reference system is " 
			<< 1e3*gap_main.in.data.geometry.delta_z0 << " [mm]" << gaplog::endl;

	if(gap_main.forces.DC.is_available())
	{
		Log << "\n * The loads associated with the DC will be calculated with the provided mesh file\n";
		Log << "   DC surface analysis report:\n" << gaplog::endl;
	
		Log << "     * Total surface area = " << 1e6*gap_main.forces.DC.Atot 
				<< "\ - [mm2]" << gaplog::endl;
		Log << "     * F0 = " << "(" 
				<< scientific << gap_main.forces.DC.F0x << ", " 
				<< scientific << gap_main.forces.DC.F0y << ", " 
				<< scientific << gap_main.forces.DC.F0z << ") - [N]" 
				<< gaplog::endl;

		Log << "     * Resultant force position: (" 
				<< gap_main.forces.DC.xR << ", " 
				<< gap_main.forces.DC.yR << ", "
				<< gap_main.forces.DC.zR << ")" 
				<< gaplog::endl; 

		Log << "     * M0 = " << "(" 
				<< scientific << gap_main.forces.DC.M0x << ", " 
				<< scientific << gap_main.forces.DC.M0y << ", " 
				<< scientific << gap_main.forces.DC.M0z << ") - [Nm]\n" 
				<< gaplog::endl;
	}
	

	Log << "\nCalculating external forces ... " << gaplog::endl;

	double FBz_avg = 0;
	double MBx_avg = 0;
	double MBy_avg = 0;

	for(unsigned int i=0; i<360; i++)
	{
		double phi = i*pi/180.0;
		gap_main.rotate(phi);
		gap_main.calc_external_forces();
		
		
		out_total << scientific << phi << "\t"
							<< scientific << gap_main.forces.FBz << "\t"
							<< scientific << gap_main.forces.xFBz << "\t"
							<< scientific << gap_main.forces.yFBz << "\t"
							<< scientific << gap_main.forces.MBx << "\t"
							<< scientific << gap_main.forces.MBy << endl;
		
		FBz_avg += gap_main.forces.FBz;
		MBx_avg += gap_main.forces.MBx;
		MBy_avg += gap_main.forces.MBy;

		out_pistons << scientific << phi << "\t"
								<< scientific << gap_main.forces.FKBx << "\t"
								<< scientific << gap_main.forces.FKBy << "\t"
								<< scientific << gap_main.forces.MKBx << "\t"
								<< scientific << gap_main.forces.MKBy << endl;

	}

	Log << "\nAverage total axial force  [N]  = " << FBz_avg/360.0 << " [N]" << gaplog::endl;
	Log << "Average total x moment     [Nm] = " << MBx_avg/360.0 << " [Nm]" << gaplog::endl;
	Log << "Average total y moment     [Nm] = " << MBy_avg/360.0 << " [Nm]" << gaplog::endl;

	Log << "\n\n * Written file ./output/block/debug/loads/total_loads.txt\n";
	Log << "\n * Written file ./output/block/debug/loads/loads/loads_from_pistons.txt\n\n";


	out_total.close();
	out_pistons.close();
}
// ------------------------------------------------------------------------- //
void cblock_gap::get_balance_factors()
{

	system("IF NOT EXIST output\\block\\debug\\bfactor (mkdir output\\block\\debug\\bfactor)");
	ofstream out("./output/block/debug/bfactor/forces.txt");
	ofstream out_balanceFactor("./output/block/debug/bfactor/balanceFactor.txt");

	out << "phi\tFfBz\tFBz\tMfBx\tMBx\MfBy\tMBy" << endl;

	Log << "\nCalculating external and hydrostatic loads ... " << gaplog::endl;

	double FBz_avg = 0;
	double MBx_avg = 0;
	double MBy_avg = 0;

	double FfBz_avg = 0;
	double MfBx_avg = 0;
	double MfBy_avg = 0;

	gap_main.film.h1 = 5e-6;
	gap_main.film.h2 = 5e-6;
	gap_main.film.h3 = 5e-6;
	gap_main.set_film_thickness();
	gap_main.set_pressure();
	reynolds.hydrodynamic = false;
	reynolds.solve_rigid();
	gap_main.update_viscosity();

	for(unsigned int i=0; i<360; i++)
	{
		double phi = i*pi/180.0;

		Log << "\rCalculating external and hydrostatic for angle " << phi*180.0/pi;

		gap_main.rotate(phi);
		gap_main.calc_external_forces();
						
		FBz_avg += gap_main.forces.FBz;
		MBx_avg += gap_main.forces.MBx;
		MBy_avg += gap_main.forces.MBy;

		// calculate fluid forces
		gap_main.update_mesh();
		gap_main.set_pressure();
		reynolds.hydrodynamic = false;
		reynolds.solve_rigid();
		gap_main.calc_fluid_forces();
		FfBz_avg += gap_main.forces.FfBz;
		MfBx_avg += gap_main.forces.MfBx;
		MfBy_avg += gap_main.forces.MfBy;

		out << phi << "\t" << gap_main.forces.FfBz << "\t" << gap_main.forces.FBz << "\t"
				<< gap_main.forces.MfBx << "\t" << gap_main.forces.MBx << "\t"
				<< gap_main.forces.MfBy << "\t" << gap_main.forces.MBy << endl;

		//ostringstream oss;
		//oss << "gap." << i << ".vtk";
		//gap_main.write_pcb_vtk(oss.str().c_str());
		

	}

	cout.precision(1);

	Log << "\n\nAverage total axial external force  [N]  = " << fixed << FBz_avg/360.0 << gaplog::endl;
	Log << "Average total x external moment      [Nm] = " << fixed << MBx_avg/360.0 << gaplog::endl;
	Log << "Average total y external moment      [Nm] = " << fixed << MBy_avg/360.0 << gaplog::endl;

	Log << "\nAverage total axial hydrostatic force  [N]  = " << fixed << FfBz_avg/360.0 << gaplog::endl;
	Log << "Average total x hydrostatic moment      [Nm] = " << fixed << MfBx_avg/360.0 << gaplog::endl;
	Log << "Average total y hydrostatic moment      [Nm] = " << fixed << MfBy_avg/360.0 << gaplog::endl;

	Log << "\nAverage Fz balance factor  [%] = " << fixed << 100*fabs(FfBz_avg/FBz_avg) << gaplog::endl;
	Log << "Average Mx balance factor    [%] = " << fixed << 100*fabs(MfBx_avg/MBx_avg) << gaplog::endl;
	Log << "Average My balance factor    [%] = " << fixed << 100*fabs(MfBy_avg/MBy_avg) << gaplog::endl << gaplog::endl;


	
	out_balanceFactor << fixed << 100*fabs(FfBz_avg/FBz_avg) << "\t"
		<< fixed << 100*fabs(MfBx_avg/MBx_avg) << "\t" << fixed << 100*fabs(MfBy_avg/MBy_avg) << endl;
	

	Log << "\n\n * Written file ./output/block/debug/bfactor/balanceFactor.txt\n";



}
// ------------------------------------------------------------------------- //
void cblock_gap::get_oil_properties()
{
	system("IF NOT EXIST output\\block\\debug\\oil (mkdir output\\block\\debug\\oil)");


	Log << "\n* -------------------- Calculation of the oil properties -------------------- *\n" 
			<< gaplog::endl;
	
	if(in.data.oil.general.oiltype == 0)
		Log << "\n\nOil type: constant properties\n\n";
	else if(in.data.oil.general.oiltype == 1)
		Log << "\n\nOil type: user defined\n\n";
	else if(in.data.oil.general.oiltype == 2)
		Log << "\n\nOil type: HLP32\n\n";
	else if(in.data.oil.general.oiltype == 3)
		Log << "\n\nOil type: user defined2\n\n";
	//else if(in.data.oil.general.oiltype == 3)
	//	Log << "\n\nOil type: Skydrol\n\n";
	else if(in.data.oil.general.oiltype == 4)
		Log << "\n\nOil type: Red Oil\n\n";
	else if(in.data.oil.general.oiltype == 5)
		Log << "\n\nOil type: SAEW10\n\n";
	else if(in.data.oil.general.oiltype == 6)
		Log << "\n\nOil type: mil 5606\n\n";
	else if(in.data.oil.general.oiltype == 7)
		Log << "\n\nOil type: Exxon DTE Excel 32\n\n";
	else if(in.data.oil.general.oiltype == 8)
		Log << "\n\nOil type: ISO 46\n\n";
	else if(in.data.oil.general.oiltype == 9)
		Log << "\n\nOil type: Skydrol - LD4\n\n";
	else if(in.data.oil.general.oiltype == 10)
		Log << "\n\nOil type: MILPRF87257\n\n"; //MILPRF87257 aerocontrolex
	
	else
	{
		Log << "\n\nOil type not supported\n\n";
		exit(1);
	}




	ofstream out_mu("./output/block/debug/oil/mu_p_T.txt");
	ofstream out_rho("./output/block/debug/oil/rho_p_T.txt");
	ofstream out_K("./output/block/debug/oil/K_p_T.txt");

	out_mu << "% LINES are mu(p=const, T), COLUMNS are mu(p, T=const)" << endl;
	out_mu << "% p range [" << gap_main.lubricant->get_p_min_limit() << ", " 
				 <<   gap_main.lubricant->get_p_max_limit() << "]" << endl;
	out_mu << "% T range [" << gap_main.lubricant->get_T_min_limit() << ", " 
				 <<   gap_main.lubricant->get_T_max_limit() << "]" << endl;

	out_rho << "% LINES are rho(p=const, T), COLUMNS are rho(p, T=const)" << endl;
	out_rho << "% p range [" << gap_main.lubricant->get_p_min_limit() << ", " 
				  <<   gap_main.lubricant->get_p_max_limit() << "]" << endl;
	out_rho << "% T range [" << gap_main.lubricant->get_T_min_limit() << ", " 
				  <<   gap_main.lubricant->get_T_max_limit() << "]" << endl;

	out_K << "% LINES are K(p=const, T), COLUMNS are K(p, T=const)" << endl;
	out_K << "% p range [" << gap_main.lubricant->get_p_min_limit() << ", " 
				<<   gap_main.lubricant->get_p_max_limit() << "]" << endl;
	out_K << "% T range [" << gap_main.lubricant->get_T_min_limit() << ", " 
				<<   gap_main.lubricant->get_T_max_limit() << "]" << endl;
	

	if(in.data.oil.general.oiltype < 1)
	{
		Log << "\nFor oil type " << in.data.oil.general.oiltype << " nothing to do.\n"
				<< gaplog::endl;
	}

	Log << "\nMin pressure " << gap_main.lubricant->get_p_min_limit()/1e5 << " [bar], ";
	Log << "Max pressure " << gap_main.lubricant->get_p_max_limit()/1e5 << " [bar]"	<< gaplog::endl;
	Log << "\nMin temperature " << gap_main.lubricant->get_T_min_limit() << " [C], ";
	Log << "Max temperature " << gap_main.lubricant->get_T_max_limit() << " [C]" 	<< gaplog::endl;

	int p_start = static_cast<int>(floor(gap_main.lubricant->get_p_min_limit()/1e5));
	int p_end = static_cast<int>(floor(gap_main.lubricant->get_p_max_limit()/1e5));
	int T_start = static_cast<int>(floor(gap_main.lubricant->get_T_min_limit()));
	int T_end = static_cast<int>(floor(gap_main.lubricant->get_T_max_limit()));

	// first column is the p value
	out_mu << "\t";
	out_rho << "\t";
	out_K << "\t";
	// write T values
	for(int j=T_start; j<T_end; j++)
	{
		double T = static_cast<double>(j);
		out_mu << T << "\t";
		out_rho << T << "\t";
		out_K << T << "\t";
	}
	out_mu << endl;
	out_rho << endl;
	out_K << endl;

	for(int i=p_start; i<p_end; i+=5)
	{
		double p = 1e5*static_cast<double>(i);
		out_mu << p << "\t";
		out_rho << p << "\t";
		out_K << p << "\t";
		
		for(int j=T_start; j<T_end; j++)
		{
			double T = static_cast<double>(j);
			out_mu << gap_main.lubricant->get_mu(p, T) << "\t";
			out_rho << gap_main.lubricant->get_rho(p, T) << "\t";
			out_K << gap_main.lubricant->get_K(p, T) << "\t";
		}
		out_mu << endl;
		out_rho << endl;
		out_K << endl;
	}

	out_mu.close();
	out_rho.close();
	out_K.close();

	Log << "\n\n * Written file ./output/block/debug/oil/mu_p_T.txt\n";
	Log << " * Written file ./output/block/debug/oil/rho_p_T.txt\n";
	Log << " * Written file ./output/block/debug/oil/K_p_T.txt\n";
	cout << "mu = " <<  gap_main.lubricant->get_mu(100e5, 50) << endl;

}
// ------------------------------------------------------------------------- //
void cblock_gap::write_vp_pressure_field()
{
	gap_main.set_film_thickness();
	gap_main.update_mesh();
	gap_main.rotate(0);
	gap_main.set_pressure();
	gap_main.set_temperature();
	reynolds.solve_rigid();
	gap_main.interpl_pressure();

	system("IF NOT EXIST output\\block\\radioss (mkdir output\\block\\radioss)");

	gap_main.VPf.writeVTK("./output/block/radioss/VPf.vtk");
	gap_main.VPf.write("./output/block/radioss/VPf.txt");
	write_vp_influgen("./output/block/radioss/VPInflugen.txt");

}
// ------------------------------------------------------------------------- //
void cblock_gap::writeVTK_block()
{
	
	Log << "\n* -------------------------- writeVTKblock ------------------------ *\n" 
			<< gaplog::endl;

	if(in.data.lubrication_module.solve_block < 1)
	{
		Log << "\n\nPlease enable \"solve_block\" in lubricating_module\n\n";
		exit(1);
	}

	if(in.data.options_block.general.Thermal_CB)
		Log << "\n*  Cylinder block Thermal option: ON" << gaplog::endl;
	else
		Log << "\n*  Cylinder block Thermal option: OFF" << gaplog::endl;

	if(in.data.options_block.general.Thermal_CB || in.data.options_block.general.Thermal_VP)
	{
		system("IF NOT EXIST output\\block\\debug\\thermal (mkdir output\\block\\debug\\thermal)");

		Log << "\nSolving Reynolds and Energy equation ... ";
		gap_main.update_mesh();
		Log << "done!" << gaplog::endl;	
		
	}
	else
	{
		Log << "\nNo thermal option activated in the options_block.txt file\n" 
				<< gaplog::endl;
		return;
	}

	// test cylinder block thermal
	if(in.data.options_block.general.Thermal_CB)
	{
		Log << "\n\nSolving for cylinder block thermo-elastic analysis" 
				<< gaplog::endl;

		// thermoelastic solver
		te_CB = new te_solver(in, "CB");

		// pointer to the interpolation grid
		interpl_grid* gap_str = &(te_CB->gap);

		Log << "\nWriting ./output/block/debug/thermal/cylinderblock.vtk ..." ;
		te_CB->write_vtk("./output/block/debug/thermal");
		te_CB->write_fset_vtk("./output/block/debug/thermal", "gap_block");
		te_CB->write_fset_vtk("./output/block/debug/thermal", "DC");
		te_CB->write_fset_vtk("./output/block/debug/thermal", "dc_mesh");
		te_CB->write_fset_vtk("./output/block/debug/thermal", "case_outer");
		te_CB->write_fset_vtk("./output/block/debug/thermal", "case_inner");
		te_CB->write_fset_vtk("./output/block/debug/thermal", "gap_case_outer");
		te_CB->write_fset_vtk("./output/block/debug/thermal", "gap_case_inner");
		for(unsigned int i = 0; i < 9; i++) 
		{
			ostringstream oss;
			oss << "gap_piston_" << i+1 ;
			te_CB->write_fset_vtk("./output/block/debug/thermal", oss.str().c_str());
		}
		Log << "done!" << gaplog::endl;

		delete te_CB;

	}
	Log << "done!" << gaplog::endl;

}
// ------------------------------------------------------------------------- //