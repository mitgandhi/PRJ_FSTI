# include "../FSTI_Block_dll/log.h"
# include "./block_ode.h"
# include <ppl.h>
# include <sstream>

# define pi 3.14159265358979323846

extern gaplog Log;

using namespace std;


//# define ANIMATION

//extern variables used in coupled caspar
extern Concurrency::event * coupled_caspar_local_event;
extern Concurrency::event * coupled_caspar_global_event;
extern bool coupled_caspar_simulation;

// ------------------------------------------------------------------------- //
cblock_ode::cblock_ode(const input& In) : in(In), block_gap(In)
{

	// ------------------------ copy from caspar_input ----------------------- //

	step_angle = in.data.options_block.general.step_angle;
	omega = in.data.operating_conditions.speed;
	t_end = in.data.lubrication_module.n_lubrication_revolutions*(2*pi)/omega;
	hmin = in.data.options_block.numeric.hmin;
	EHD_CB = in.data.options_block.general.EHD_CB;
	EHD_VP = in.data.options_block.general.EHD_VP;
	EHD = (EHD_CB || EHD_VP);
	Thermal_CB = in.data.options_block.general.Thermal_CB;
	Thermal_VP = in.data.options_block.general.Thermal_VP;
	Thermal = (Thermal_CB || Thermal_VP);
	EHD_debug = in.data.options_block.general.EHD_loop_debug;
	dense_vtk = in.data.options_block.general.dense_vtk_out;
	dense_light_vtk = in.data.options_block.general.dense_vtk_out_light;
	contact_opt = in.data.options_block.numeric.contact;

	// ----------------------------------------------------------------------- //

	phi = 0.0;
	t = 0.0;
	revcounter = 0;
	T = 2.0*pi/omega;
	dphi = step_angle;
	dt = dphi/omega;

	output.leak = 0;
	output.Mfr = 0;
	output.leak_avg = 0;
	output.Qin_leak;				
	output.Qout_leak;				
	output.phid;						
	output.Qin_leak_avg;		
	output.Qout_leak_avg;		
	output.phid_avg;				
	output.Pl;							
	output.Plf;							
	output.Pll;	
	output.FTK = 0;

}
// ------------------------------------------------------------------------- //
cblock_ode::~cblock_ode()
{

}
// ------------------------------------------------------------------------- //
void cblock_ode::initialize()
{
	
	system("IF NOT EXIST output (mkdir output)");
	system("IF NOT EXIST output\\block (mkdir output\\block)");
	system("IF NOT EXIST output\\block\\rev.0 (mkdir output\\block\\rev.0)");
	system("IF NOT EXIST output\\block\\rev.0 (mkdir output\\block\\influgen)");
	if(dense_vtk || dense_light_vtk)
		system("IF NOT EXIST output\\block\\rev.0\\dense_vtk (mkdir output\\block\\rev.0\\dense_vtk)");
	if(EHD_debug)
		system("IF NOT EXIST output\\block\\rev.0\\EHD (mkdir output\\block\\rev.0\\EHD)");

	out.open("./output/block/block.txt");
	
	if(!out.is_open())
	{
		Log << "\ncblock_ode::initialize() Error! Could not open ./output/block/block.txt"
				 << gaplog::endl;
		exit(1);
	}

	if(EHD) // load the influence matrices
	{

		// load the influence matrices
		block_gap.reynolds.load_matrices();

		if(in.data.options_block.general.EHD_CB)
			block_gap.reynolds.init_def_map(CB);
		if(in.data.options_block.general.EHD_VP)
			block_gap.reynolds.init_def_map(VP);

		block_gap.gap_main.set_pressure();
		block_gap.reynolds.solve_rigid();
		block_gap.gap_main.interpl_pressure();
	}

	// -------------------------------------

	# ifdef ANIMATION
	
	system("IF NOT EXIST output\\block\\amination (mkdir output\\block\\animation)");

	# endif


	// ------------------------------------
	
	// print the column headers in the output file

	out << "       t     " << "\t"
			 << "      rev    " << "\t"
			 << "      phi    " << "\t"
			 << "      h1     " << "\t"
			 << "      h2     " << "\t"
			 << "      h3     " << "\t"
			 << "     dhdt1   " << "\t"
			 << "     dhdt2   " << "\t"
			 << "     dhdt3   " << "\t"
			 << "      vc1    " << "\t"
			 << "      vc2    " << "\t"
			 << "      vc3    " << "\t"
			 << "     hmin    " << "\t"
			 << "     hmax    " << "\t"
			 << "     havg    " << "\t"
			 << "      MBx    " << "\t"
			 << "      MBy    " << "\t"
			 << "      FBz    " << "\t"
			 << "     MfBx    " << "\t"
			 << "     MfBy    " << "\t" 
			 << "     FfBz    " << "\t" 
			 << "     Leak    " << "\t" 
			 << "     Ploss   " << "\t" 
			 << "    PlossL   " << "\t" 
			 << "    PlossF   " << "\t"
			 << "      FTK    " << "\t" 
			 << "      FTG    " << endl;

	out.close();

	
	Log << "\n// ------------------------- Analysis options summary --------------------- //\n";
	
	// oil type
	if(in.data.oil.general.oiltype == 0)
	{
		Log << "\n * Constant oil properties oil model will be used.\n";
	}
	else if(in.data.oil.general.oiltype == 1)
	{
		Log << "\n * User defined oil properties oil model will be used.\n";
	}
	else if(in.data.oil.general.oiltype == 2)
	{
		Log << "\n * HLP32 oil model will be used.\n";
	}

	// macrogeometry initialization

	// block
	if(in.data.options_block.general.macro_CB == 0)
	{
		Log << "\n * No macro geometry will be applied to the block surface.\n";
	}
	else if(in.data.options_block.general.macro_CB == 1)
	{
		Log << "\n * Axisymmetric micro geometry profile will be applied to the block surface.\n";
	}
	else if(in.data.options_block.general.macro_CB == 2)
	{
		Log << "\n * Waved Surface micro geometry profile will be applied to the block surface.\n";
	}
	else if(in.data.options_block.general.macro_CB == 3)
	{
		Log << "\n * Spherical micro geometry profile will be applied to the block surface.\n";
	}
	else if(in.data.options_block.general.macro_CB == 4)
	{
		Log << "\n * A geometry profile defined through a grid will be applied to the block surface.\n";
	}
	else
	{
		Log << "\n * WARNING: Micro geometry option " <<in.data.options_block.general.macro_CB
			<< " on the cylinder block surface is invalid; it will be disregarded.\n";
	}

	// valve plate
	if(in.data.options_block.general.macro_VP == 0)
	{
		Log << "\n * No macro geometry will be applied to the valve plate surface.\n";
	}
	else if(in.data.options_block.general.macro_VP == 1)
	{
		Log << "\n * Axisymmetric micro geometry profile will be applied to the valve plate surface.\n";
	}
	else if(in.data.options_block.general.macro_VP == 2)
	{
		Log << "\n * Waved Surface micro geometry profile will be applied to the valve plate surface.\n";
	}
	else if(in.data.options_block.general.macro_CB == 3)
	{
		Log << "\n * Spherical micro geometry profile will be applied to the valve plate surface.\n";
	}
	else if(in.data.options_block.general.macro_CB == 4)
	{
		Log << "\n * A geometry profile defined through a grid will be applied to the valve plate surface.\n";
	}
	else
	{
		Log << "\n * WARNING: Micro geometry option " << in.data.options_block.general.macro_VP
			<< " on the valve plate surface is invalid; it will be disregarded\.n";
	}

	// check the external forces calculation method
	if(in.data.options_block.numeric.Fext_method == 0)
	{
		Log << "\n * The Simplified calculation method will be used for the loads transferred \n"
			  << "   to the cylinder block body from the piston/slipper assembly.\n";
	}
	else
	{
		Log << "\n * The Advanced calculation method will be used for the loads transferred \n"
			  << "   to the cylinder block body from the piston/slipper assembly.\n";
	}

	if(block_gap.gap_main.forces.DC.is_available())
	{
		Log << "\n * The loads associated with the DC will be calculated with the provided mesh file\n";
		Log << "   DC surface analysis report:\n" << gaplog::endl;
	
		Log << "     * Total surface area = " << 1e6*block_gap.gap_main.forces.DC.Atot 
				<< "\ - [mm2]" << gaplog::endl;
		Log << "     * F0 = " << "(" 
				<< scientific << block_gap.gap_main.forces.DC.F0x << ", " 
				<< scientific << block_gap.gap_main.forces.DC.F0y << ", " 
				<< scientific << block_gap.gap_main.forces.DC.F0z << ") - [N]" 
				<< gaplog::endl;

		Log << "     * Resultant force position: (" 
				<< block_gap.gap_main.forces.DC.xR << ", " 
				<< block_gap.gap_main.forces.DC.yR << ", "
				<< block_gap.gap_main.forces.DC.zR << ")" 
				<< gaplog::endl; 

		Log << "     * M0 = " << "(" 
				<< scientific << block_gap.gap_main.forces.DC.M0x << ", " 
				<< scientific << block_gap.gap_main.forces.DC.M0y << ", " 
				<< scientific << block_gap.gap_main.forces.DC.M0z << ") - [Nm]\n" 
				<< gaplog::endl;
	}
	else
	{
		Log << "\n * DC loads will be calculated with a simplified approach because \n"
			  << "   the DC mesh was speficied or could not be read\n";

		double AK = 0.6*(0.25*pi*in.data.geometry.dK*in.data.geometry.dK);
		if( fabs(in.data.geometry.ADC - AK)/AK > 0.9)
		{
			Log << "   WARNING: It seems that the provided ADC value is unreasonable."	<< gaplog::endl;
			Log << "  Are you sure that the value " << in.data.geometry.ADC*1e6 
					<< " [mm^2] is correct?" << gaplog::endl << gaplog::endl;
		}
	}

	if(in.data.geometry.d_spherical > 0)
	{
		Log << "\n * The reference system used for the loads calculation on the cylinder block body\n"
				<< "   will be at a distance of " << 1e3*block_gap.gap_main.pistons.z0 
				<< " [mm] from the point defined by the\n"
				<< "   intersection of the spherical surface of the cylinder block with the shaft axis." 
				<< "\n";
		Log	<< "\n * The reaction force of the spline joint on the cylider block body will be \n"
				<< "   applied at a distance " << 1e3*block_gap.gap_main.in.data.geometry.delta_z0 << " [mm]"
				<< " from the reference system."
				<< "\n";
	}
	else
	{
		Log << "\n * The reference system used for the loads calculation on the cylinder block body\n"
				<< "   will be at a distance of " << 1e3*block_gap.gap_main.pistons.z0 
				<< " [mm] from the cylinder block's sealing land surface\n"
				<< "\n";
	}

}
// ------------------------------------------------------------------------- //
void cblock_ode::print_state()
{
	Log << "\n// ------------------------------------------------------------"
			<< "------------- //\n" << gaplog::endl;

	//time(&overall_time);
	//timeinfo = localtime(&overall_time);

	// get elapsed time
	int hours, minutes, seconds;
	time_t elapsed = (clock() - overall_time)/CLOCKS_PER_SEC;
	hours = (elapsed / 3600);
	minutes = (elapsed / 60) - (hours * 60);
	seconds = elapsed%60;

	Log << "Elapsed time: " 
			<< hours << " h, " 
			<< minutes << " m, " 
			<< seconds << " s.\n\n";

	Log << "Current status:\n\n";
	
	Log << "Revolution " << revcounter << ", ";

	cout.precision(1);

	Log << "Shaft angle: " << fixed << phi*(180.0/pi) 
			<< ", pDC: " << block_gap.gap_main.pistons.pDC[0]/1e5 << " [bar]"
			<< "\n\n";

	//string now(asctime(timeinfo));
	//now.erase(std::remove(now.begin(), now.end(), '\n'), now.end());

	cout.precision(3);

	cout << showpos;
	
	Log << "  dhdt[1] = " << scientific << block_gap.gap_main.squeeze.dhdt1 
			<< "  dhdt[2] = " << scientific << block_gap.gap_main.squeeze.dhdt2 
			<< "  dhdt[3] = " << scientific << block_gap.gap_main.squeeze.dhdt3 
			<< gaplog::endl;
	
	Log << "  h[1]    = " << scientific << block_gap.gap_main.film.h1 
			<< "  h[2]    = " << scientific << block_gap.gap_main.film.h2 
			<< "  h[3]    = " << scientific << block_gap.gap_main.film.h3 
			<< gaplog::endl;

	cout << noshowpos;
	
	Log << "  hmin    =  " << scientific << block_gap.gap_main.hmin() 
			<< "  hmax    =  " << scientific << block_gap.gap_main.hmax() 
			<< "  havg    =  " << scientific << block_gap.gap_main.havg() 
			<< gaplog::endl;

}
// ------------------------------------------------------------------------- //
void cblock_ode::run()
{

	clock_start = clock();
	overall_time = clock();

	double next_step = 0.0;

	ostringstream oss;

	std::streamsize default_val = cout.precision();

	// ----------------------- Run the thermal analysis ----------------------- //

	if(in.data.options_block.general.StartWithTH)
	{

		if(in.data.options_block.general.Thermal_CB || in.data.options_block.general.Thermal_VP)
		{

			Log << "\nEstimating heat fluxes over one shaft revolution ... ";

			block_gap.gap_main.fields.qcb_avg = 0;
			block_gap.gap_main.fields.qvp_avg = 0;
			block_gap.gap_main.film.h1 = 10e-6;
			block_gap.gap_main.film.h2 = 10e-6;
			block_gap.gap_main.film.h3 = 10e-6;

			for(int i=0; i<360; i+=18)
			{
			
				double phi = i*pi/180.0;
				block_gap.gap_main.rotate(phi);
				block_gap.reynolds.update_dimension();
				block_gap.energy.update_dimension();
				block_gap.gap_main.set_film_thickness();	
				block_gap.gap_main.update_mesh();
		
				block_gap.gap_main.set_pressure();
				block_gap.gap_main.set_temperature();
				block_gap.gap_main.update_viscosity();
				block_gap.reynolds.solve_rigid();
				block_gap.gap_main.calc_mass_flow();
				block_gap.gap_main.update_h_previous();
				block_gap.energy.solve();
			
				block_gap.gap_main.fields.qcb_avg = 
					block_gap.gap_main.fields.qcb_avg  + block_gap.gap_main.calc_cb_heatflux();
			
				block_gap.gap_main.fields.qvp_avg = 
					block_gap.gap_main.fields.qvp_avg + block_gap.gap_main.calc_vp_heatflux();
			}

			block_gap.gap_main.fields.qcb_avg = 0.5*block_gap.gap_main.fields.qcb_avg/20.0;
			block_gap.gap_main.fields.qvp_avg = 0.5*block_gap.gap_main.fields.qvp_avg/20.0;
			block_gap.gap_main.fields.qcb_prog = block_gap.gap_main.fields.qcb_avg;
			block_gap.gap_main.fields.qvp_prog = block_gap.gap_main.fields.qvp_avg;

			Log << "done!" << gaplog::endl;

			if(in.data.options_block.general.Thermal_CB)
			{
			
				Log << "\n\nSolving for the cylinder block thermo-elastic analysis" 
						<< gaplog::endl;

				
				// run the thermal analysis

				block_gap.te_CB = new te_solver(in, "CB");
				
				// pointer to the interpolation grid
				interpl_grid* gap_str = &(block_gap.te_CB->gap);

				// apply the heat flux from gap
				block_gap.gap_main.apply_heat_flux(gap_str, "CB");	
				// solve for heat transfer
				block_gap.te_CB->solve_htr();
				// interpolation of temperature from structure to fluid
				interpolation str_to_fl;
				str_to_fl.cellsTocells(gap_str, &block_gap.gap_main.CBf);
				// update the surface temperature
				block_gap.gap_main.update_Tcb();	

				// solve for thermal deflection
				block_gap.te_CB->solve_thdfl();
				// interpolation of deformation from structure to fluid
				str_to_fl.pointsTocells(gap_str, &block_gap.gap_main.CBf);
				// update the gap deformation
				block_gap.gap_main.update_dhcb_thermal();

				// write output
				system("IF NOT EXIST output\\block\\starting (mkdir output\\block\\starting)");
				block_gap.te_CB->write_vtk("./output/block/starting");
				block_gap.te_CB->write_fset_vtk("./output/block/starting", "gap_block");

				delete block_gap.te_CB;
			}

			if(in.data.options_block.general.Thermal_VP)
			{
			
				Log << "\n\nSolving for valve plate thermo-elastic analysis" 
						<< gaplog::endl;

				block_gap.te_VP = new te_solver(in, "VP");

				// pointer to the interpolation grid
				interpl_grid* gap_str = &(block_gap.te_VP->gap);

				// apply the heat flux from gap
				block_gap.gap_main.apply_heat_flux(gap_str, "VP");	
				// solve for heat transfer
				block_gap.te_VP->solve_htr();
				// interpolation of temperature from structure to fluid
				interpolation str_to_fl;
				str_to_fl.cellsTocells(gap_str, &block_gap.gap_main.VPf);
				// update the surface temperature
				block_gap.gap_main.update_Tvp();	

				// solve only if required
				if(in.data.thermal.valveplate.calc_deflect_on[0].compare("NONE") != 0)
				{

					// solve for thermal deflection
					block_gap.te_VP->solve_thdfl();
					// interpolation of deformation from structure to fluid
					str_to_fl.pointsTocells(gap_str, &block_gap.gap_main.VPf);
					// update the gap deformation
					block_gap.gap_main.update_dhvp_thermal();

				}

				system("IF NOT EXIST output\\block\\starting (mkdir output\\block\\starting)");
				block_gap.te_VP->write_vtk("./output/block/starting");
				block_gap.te_VP->write_fset_vtk("./output/block/starting", "gap");

				delete block_gap.te_VP;
			}

			// initialize the fields for the gap calculation
			block_gap.gap_main.fields.qcb_prog = 0;
			block_gap.gap_main.fields.qcb_avg = 0;
			block_gap.gap_main.fields.qvp_prog = 0;
			block_gap.gap_main.fields.qvp_avg = 0;

			block_gap.gap_main.film.h1 = in.data.options_block.position.hB1;
			block_gap.gap_main.film.h2 = in.data.options_block.position.hB2;
			block_gap.gap_main.film.h3 = in.data.options_block.position.hB3;
			block_gap.gap_main.set_film_thickness();

			Log << "\nWriting ./output/block/debug/thermal/fluid.vtk ..." ;
			block_gap.gap_main.write_vtk("./output/block/starting/fluid.vtk", 1000);
			Log << "done!" << gaplog::endl;
		
		}

	}

	// ------------------------------------------------------------------------ //
	// ------------------------------------------------------------------------ //
	// ------------------------------------------------------------------------ //

		
	Log << "\n\n\n* ---- Starting the gap flow calculation ---- *\n\n\n";

	while(t < t_end)
	{
		
		print_state();

		// rotate the block

		Log << "\nRotating block ... ";

		block_gap.gap_main.t = t;
		block_gap.gap_main.rotate(phi);
		block_gap.reynolds.update_dimension();
		block_gap.energy.update_dimension();
		
		Log << "done (steps: " << block_gap.gap_main.mesh.reynolds.rotation_steps 
				<< ")" << gaplog::endl;

		Log << "Seting relative position between cylinder block and valve plate ... ";
		
		block_gap.gap_main.set_film_thickness();	
		
		Log << "done!" << gaplog::endl;

		//double _cnt_max = block_gap.gap_main.film.contact.max();
		//if(_cnt_max > 0)
		//{
		//	cout.precision(1);
		//	Log << "WARNING: Possible contact detected, theoretical penetration: " << fixed 
		//			<< 1e6*(_cnt_max) << "[um]\n"; 
		//}

		Log << "Updating p and T boundary conditions ... ";
		
		block_gap.gap_main.set_pressure();
		block_gap.gap_main.set_temperature();

		Log << "done!" << gaplog::endl;


		//if(EHD_CB || EHD_VP)
		//{
		//	Log << "\nCalculating external forces before FSI loop ... ";
		//	block_gap.gap_main.calc_external_forces();
		//	Log << "done!" << gaplog::endl;

		//	block_gap.gap_main.calc_fluid_forces();
		//	block_gap.gap_main.calc_dF();

		//	Log << "\nForce imbalance: \n" << gaplog::endl
		//			 << "  "
		//			 << "dF[0] = " << block_gap.gap_main.forces.dF[0]
		//			 << "  dF[1] = " << block_gap.gap_main.forces.dF[1]
		//			 << "  dF[2] = " << block_gap.gap_main.forces.dF[2] 
		//			 << gaplog::endl;

		//	cout.precision(default_val);
		//
		//	Log << "\nForce balance loop ... " << gaplog::endl;
		//
		//	int fbalance = block_gap.get_block_balance();

		//	if(fbalance != GSL_SUCCESS)
		//	{
		//		Log << "\nWARNING: something wrong in the force balance loop!\n" 
		//				 << gaplog::endl;

		//		Log << "\nWARNING: something wrong in the force balance loop!\n" 
		//			 << gaplog::endl;
		//	
		//	}
		//}

		// ------------------ solve the Reynolds equation ------------------------ //

		// set the squeeze
		block_gap.gap_main.set_squeeze();
		
		if(EHD)
		{
			Log << gaplog::endl << gaplog::endl;
			block_gap.reynolds.solve_EHD();

			// update the field which stores the EHD solution
			block_gap.gap_main.fields.p_EHD = block_gap.gap_main.fields.p;
		}
		else
		{
			Log << "\nSolving the Reynolds equation ... ";
			block_gap.reynolds.solve_rigid();
			Log << "done!" << gaplog::endl;
		}

		block_gap.gap_main.update_mesh();

		// ------------------------ solve the energy equation -------------------- //
		
		Log << "\nSolving the Energy equation ... ";

		block_gap.gap_main.calc_mass_flow();
		block_gap.energy.solve();

		Log << "done!" << gaplog::endl;

		Log << "\nUpdating oil properties ... ";

		block_gap.gap_main.update_density();
		block_gap.gap_main.update_viscosity();

		// store previous film thickness - Rene
		block_gap.gap_main.update_h_previous();
		
		Log << "done!" << gaplog::endl;

		Log << "\nCalculating heat fluxes ... ";

		// heat flux towards cylinder block
		block_gap.gap_main.fields.qcb = block_gap.gap_main.calc_cb_heatflux();
		block_gap.gap_main.fields.qcb_avg = 
			block_gap.gap_main.fields.qcb_avg + block_gap.gap_main.fields.qcb;
		// heat flux towards valve plate
		block_gap.gap_main.fields.qvp = block_gap.gap_main.calc_vp_heatflux();
		block_gap.gap_main.fields.qvp_avg = 
			block_gap.gap_main.fields.qvp_avg + block_gap.gap_main.fields.qvp;

		// update average dissipation field
		block_gap.gap_main.fields.phid_avg = 
			block_gap.gap_main.fields.phid_avg + block_gap.gap_main.fields.phid;

		Log << "done!" << gaplog::endl;

		
		// -------------------------- write VTK output ------------------------- //
		
		block_gap.gap_main.update_mesh();
		
		if(dense_light_vtk)
		{
			oss.str("");
			oss.precision(2);
			oss << "./output/block/rev." << revcounter << "/dense_vtk/phi." << fixed 
					<< (180/pi)*phi << "_gap.vtk";
			Log << "\nWriting " << oss.str();
			block_gap.gap_main.write_light_vtk(oss.str().c_str(), 1000);
			Log << "done!\n" << gaplog::endl;
			
			//oss.str("");
			//oss << "./output/block/rev." << revcounter << "/dense_vtk/p." << fixed 
			//	<< (180/pi)*phi << "_cb.vtk";
			//block_gap.gap_main.write_pcb_vtk(oss.str().c_str());

		}
		if(dense_vtk)
		{
			oss.str("");
			oss.precision(2);
			oss << "./output/block/rev." << revcounter << "/dense_vtk/phi." << fixed 
					<< (180/pi)*phi << "_gap.vtk";
			Log << "\nWriting " << oss.str();
			block_gap.gap_main.write_vtk(oss.str().c_str(), 1000);
			Log << "done!\n" << gaplog::endl;
		}
		
		// ----------------------- write stuff for making animations ----------------------- //

		# ifdef ANIMATION

		ostringstream oss;
		oss << "IF NOT EXIST output\\block\\amination\\"
				<< "rev." << revcounter << " (mkdir output\\block\\animation\\"
				<< "rev." << revcounter << ")";
		system(oss.str().c_str());
		oss.str("");
		oss.precision(0);
		oss << "./output/block/animation/rev." << revcounter << "/block." 
				<< fixed << (180/pi)*phi << ".txt";
		block_gap.write_cb_influgen(oss.str().c_str());
		oss.str("");
		oss << "./output/block/animation/rev." << revcounter << "/valveplate." 
				<< fixed << (180/pi)*phi << ".txt";
		block_gap.write_vp_influgen(oss.str().c_str());
		oss.str("");
		
		# endif



		
		// write influgen inputs for postprocessing
		//{
			//oss.str("");
			//oss.precision(2);
			//oss << "IF NOT EXIST output\\block\\rev." << revcounter 
			//		<< "\\influgen\\" << (180/pi)*phi << " (mkdir output\\block\\rev." << revcounter 
			//		<< "\\influgen\\" << (180/pi)*phi << ")";
			//system(oss.str().c_str());
			//oss.str("");
			//oss.precision(2);
			//oss << "./output/block/rev." << revcounter << "/influgen/" << (180/pi)*phi << "/cb_influgen.txt";
			//block_gap.write_cb_influgen(oss.str().c_str());
			//oss.str("");
			//oss.precision(2);
			//oss << "./output/block/rev." << revcounter << "/influgen/" << (180/pi)*phi << "/vp_influgen.txt";
			//block_gap.write_vp_influgen(oss.str().c_str());
		//}


		// ----------------------- force balance loop -------------------------- //

		//bool balance = false;
		//int fiters = 0;
		//double norm = 0;

		Log << "\nCalculating external forces ... ";
		block_gap.gap_main.calc_external_forces();
		Log << "done!" << gaplog::endl;

		block_gap.gap_main.calc_fluid_forces();
		block_gap.gap_main.calc_dF();

		Log << "\nForce imbalance: \n" << gaplog::endl
				 << "  "
				 << "dF[0] = " << block_gap.gap_main.forces.dF[0]
				 << "  dF[1] = " << block_gap.gap_main.forces.dF[1]
				 << "  dF[2] = " << block_gap.gap_main.forces.dF[2] 
				 << gaplog::endl;

		cout.precision(default_val);
		
		Log << "\nForce balance loop ... " << gaplog::endl;
		
		int fbalance = block_gap.get_block_balance();

		if(fbalance != GSL_SUCCESS)
		{
			Log << "\nWARNING: something wrong in the force balance loop!\n" 
					 << gaplog::endl;
			
			if(block_gap.gap_main.use_fsi_fb)
			{
				Log << "Trying with no fsi term in the force balance ... \n\n";
				
				block_gap.gap_main.use_fsi_fb = false;
				fbalance = block_gap.get_block_balance();
				block_gap.gap_main.use_fsi_fb = true;
			}
			
		}

		// ------------------- calculate performance parameters ------------------ //

		double pHP = block_gap.gap_main.pfile.getp_rad(phi)[2];
		double pLP = block_gap.gap_main.pfile.getp_rad(phi)[1];
		double dp = pHP - pLP;
		
		output.leak = block_gap.gap_main.calc_leakage();
		output.phid = block_gap.gap_main.fields.phid.sum();
		
		//output.Q_in_leak = block_gap.gap_main.calcQin();
		//output.Q_out_leak = block_gap.gap_main.calcQout();
		
		output.leak_avg += output.leak;
		output.phid_avg += output.phid;
		
		//output.Qin_leak_avg += output.Q_in_leak;
		//output.Qout_leak_avg += output.Q_out_leak;
		
		// total power loss
		output.Pl = output.phid;
		// leakage power loss
		output.Pll = dp*output.leak;
		// friction loss
		output.Plf = output.Pl - output.Pll;
		output.FTK = block_gap.gap_main.calcFTK(phi);

		write_output();

		// -------------- integration (simple euler integration) --------------- //

		// do nothing, just integrate (negative film may occur)
		if(contact_opt == 0)	
		{
			Log << "\nIntegrating film thickness:\n";
			Log << "h[1](t) = " << block_gap.gap_main.film.h1 << "\t";
			Log << "h[2](t) = " << block_gap.gap_main.film.h2 << "\t";
			Log << "h[3](t) = " << block_gap.gap_main.film.h3 << "\n";
			block_gap.gap_main.film.h1 += dt*block_gap.gap_main.squeeze.dhdt1;
			block_gap.gap_main.film.h2 += dt*block_gap.gap_main.squeeze.dhdt2;
			block_gap.gap_main.film.h3 += dt*block_gap.gap_main.squeeze.dhdt3;
			Log << "h[1](t+dt) = " << block_gap.gap_main.film.h1 << "\t";
			Log << "h[2](t+dt) = " << block_gap.gap_main.film.h2 << "\t";
			Log << "h[3](t+dt) = " << block_gap.gap_main.film.h3 << "\n\n";
		}
		// integrate just positive velocities
		else
		{
			// no contact
			if(block_gap.gap_main.film.contact.max() == 0)
			{
				block_gap.gap_main.film.h1 += dt*block_gap.gap_main.squeeze.dhdt1;
				block_gap.gap_main.film.h2 += dt*block_gap.gap_main.squeeze.dhdt2;
				block_gap.gap_main.film.h3 += dt*block_gap.gap_main.squeeze.dhdt3;
				
				block_gap.gap_main.squeeze.d_dhdt1 = 0;
				block_gap.gap_main.squeeze.d_dhdt2 = 0;
				block_gap.gap_main.squeeze.d_dhdt3 = 0;
			}
			else // contact
			{
				// get the contact velocities
				block_gap.get_contact_velocities();
			
				Log << "\ndh/dt correction: "
					<< "P1: " << scientific << block_gap.gap_main.squeeze.d_dhdt1
					<< " P2: " << scientific << block_gap.gap_main.squeeze.d_dhdt2
					<< " P3: " << scientific << block_gap.gap_main.squeeze.d_dhdt3
					<< "\n" << gaplog::endl;

				block_gap.gap_main.film.h1 += dt*block_gap.gap_main.squeeze.d_dhdt1;
				block_gap.gap_main.film.h2 += dt*block_gap.gap_main.squeeze.d_dhdt2;
				block_gap.gap_main.film.h3 += dt*block_gap.gap_main.squeeze.d_dhdt3;
			}
		}		
		
				
		// ----------------------- move forward in time ------------------------ //
		
		t += dt;
		phi += dphi;

		// check for and of revolution
		end_revolution();

	}
}
// ------------------------------------------------------------------------- //
void cblock_ode::end_revolution()
{
	
	if(phi <= (2.0*pi - step_angle*pi/180.0))
	{
		return;
	}
		
	// ------------------------ copy from caspar_input ----------------------- //

	double alpha_CB = in.data.options_block.numeric.relax_CB;
	double alpha_VP = in.data.options_block.numeric.relax_VP;
	alpha_CB = (alpha_CB <= 1.0) ? alpha_CB : 1.0;
	alpha_VP = (alpha_VP <= 1.0) ? alpha_VP : 1.0;
	
	// ----------------------------------------------------------------------- //

	// Handle coupled caspar blocking events
	if(coupled_caspar_simulation)
	{
		//this is a new revolution
		coupled_caspar_local_event->set();			//signal that the slipper has finished a revolution
		coupled_caspar_global_event->wait();		//wait for the all interfaces to finish
		coupled_caspar_global_event->reset();		//reset the global event flag
	}

	// get elapsed time
	int hours, minutes, seconds;
	clock_end = clock();
	time_t elapsed = (clock_end - clock_start)/CLOCKS_PER_SEC;
	hours = (elapsed / 3600);
	minutes = (elapsed / 60) - (hours * 60);
	seconds = elapsed%60;

	Log << "\n// ------------------------------------------------------------------------- //\n" 
			<< gaplog::endl;

	Log << "Revolution " << revcounter << " completed in "
			<< hours << " h, " << minutes << " min and " << seconds << " s\n\n"; 
					
	// ----------------------------------------------------------------------- //

	phi = 0.0;
	block_gap.gap_main.rotate(phi);

	// ---------------------------- average fields --------------------------- //

	output.leak_avg = output.leak_avg*(dt/T);
	
	Log << "Average leakage [l/min]: " << 60000*output.leak_avg << gaplog::endl;

	output.leak_avg = 0;
	
	// --------------- delete the influence matrices to save memory ---------- //
	if(EHD)
	{
		Log << "\nDeleting influence matrices to save memory ... ";
		block_gap.reynolds.delete_matrices();
		Log << "done!\n" << gaplog::endl;
	}

	// --------------------- thermal analysis stuff ------------------------ //

	if(Thermal)
	{
		// ------------------- get the average heat flux ------------------- //

		Log << "Averaging the heat fluxes ... ";

		// calculate the average heat fluxes
		block_gap.gap_main.fields.qcb_avg = (dt/T)*block_gap.gap_main.fields.qcb_avg;
		block_gap.gap_main.fields.qvp_avg = (dt/T)*block_gap.gap_main.fields.qvp_avg;

		// limit the heat fluxes
		double q_min_limit = in.data.options_block.numeric.q_min_limit;
		double q_max_limit = in.data.options_block.numeric.q_max_limit;
		block_gap.gap_main.fields.qcb_avg.limit(q_min_limit, q_max_limit);
		block_gap.gap_main.fields.qvp_avg.limit(q_min_limit, q_max_limit);

		Log << "done!" << gaplog::endl;

		
		// calculate the average dissipation field
		output.phid_avg = output.phid_avg*(dt/T);
		Log << "\nAverage dissipation: " << output.phid_avg << " [W]" << gaplog::endl;

		// update the heat fluxes: new = old + alpha(new - old)
		
		// cylinder block
		block_gap.gap_main.fields.qcb_prog =	
			block_gap.gap_main.fields.qcb_prog
			+ alpha_CB*(block_gap.gap_main.fields.qcb_avg - block_gap.gap_main.fields.qcb_prog);
		
		// valve plate
		block_gap.gap_main.fields.qvp_prog =	
			block_gap.gap_main.fields.qvp_prog
			+ alpha_VP*(block_gap.gap_main.fields.qvp_avg - block_gap.gap_main.fields.qvp_prog);

		// write the average heat flux on cylinder block to the exchange file 
		ofstream flux_out("./output/block/block_flux.txt");
		flux_out << block_gap.gap_main.get_qcb_avg() << endl;
		flux_out.close();
		

		// thermal analysis on cylinder block
		if(Thermal_CB)
		{

			block_gap.te_CB = new te_solver(in, "CB");
				
			// pointer to the interpolation grid
			interpl_grid* gap_str = &(block_gap.te_CB->gap);

			// apply the heat flux from gap
			block_gap.gap_main.apply_heat_flux(gap_str, "CB");	
			// solve for heat transfer
			block_gap.te_CB->solve_htr();
			// interpolation of temperature from structure to fluid
			interpolation str_to_fl;
			str_to_fl.cellsTocells(gap_str, &block_gap.gap_main.CBf);
			// update the surface temperature
			block_gap.gap_main.update_Tcb();	

			if(in.data.thermal.block.calc_deflect_on[0].compare("NONE") != 0)
			{
				// solve for thermal deflection
				block_gap.te_CB->solve_thdfl();
				// interpolation of deformation from structure to fluid
				str_to_fl.pointsTocells(gap_str, &block_gap.gap_main.CBf);
				// update the gap deformation
				block_gap.gap_main.update_dhcb_thermal();
			}
			ostringstream oss;
			oss.str("");
			oss << "./output/block/rev." << revcounter;
			block_gap.te_CB->write_vtk(oss.str().c_str());
			block_gap.te_CB->write_fset_vtk(oss.str().c_str(), "gap_block");

			// delete the solver structure 
			delete block_gap.te_CB;
			block_gap.te_CB = 0;

		}
		// thermal analysis on valve plate
		if(Thermal_VP)
		{

			block_gap.te_VP = new te_solver(in, "VP");

			// pointer to the interpolation grid
			interpl_grid* gap_str = &(block_gap.te_VP->gap);

			// apply the heat flux from gap
			block_gap.gap_main.apply_heat_flux(gap_str, "VP");	
			// solve for heat transfer
			block_gap.te_VP->solve_htr();
			// interpolation of temperature from structure to fluid
			interpolation str_to_fl;
			str_to_fl.cellsTocells(gap_str, &block_gap.gap_main.VPf);
			// update the surface temperature
			block_gap.gap_main.update_Tvp();	

			// solve only if required
			if(in.data.thermal.valveplate.calc_deflect_on[0].compare("NONE") != 0)
			{
				// solve for thermal deflection
				block_gap.te_VP->solve_thdfl();
				// interpolation of deformation from structure to fluid
				str_to_fl.pointsTocells(gap_str, &block_gap.gap_main.VPf);
				// update the gap deformation
				block_gap.gap_main.update_dhvp_thermal();
			}


			ostringstream oss;
			oss.str("");
			oss << "./output/block/rev." << revcounter;
			block_gap.te_VP->write_vtk(oss.str().c_str());
			block_gap.te_VP->write_fset_vtk(oss.str().c_str(), "gap");


			// delete the solver structure 
			delete block_gap.te_VP;
			block_gap.te_VP = 0;

		}
	}
		
	// reload the matrices
	if(EHD)
	{
		Log << "\nReloading influece matrices ... " << gaplog::endl;
		block_gap.reynolds.load_matrices();
		Log << "matrices reloaded!\n" << gaplog::endl;
	}

	// solve the reynolds equation
	block_gap.gap_main.set_film_thickness();
	block_gap.gap_main.update_mesh();
	block_gap.gap_main.set_pressure();

	if(EHD)
		block_gap.reynolds.solve_EHD();
	else
		block_gap.reynolds.solve_rigid();

	// calculate the mass flow and solve for the enrgy equation
	block_gap.gap_main.calc_mass_flow();
	block_gap.gap_main.set_temperature();
	block_gap.energy.solve();
	block_gap.gap_main.update_viscosity();

	// write gap output
	ostringstream oss;
	oss.str("");
	oss << "./output/block/rev." << revcounter << "/fluid.vtk";
	block_gap.gap_main.write_vtk(oss.str().c_str(), 1000);

	// reset to zero the averaged fields
	output.phid_avg = 0;
	block_gap.gap_main.fields.phid_avg = 0;
	block_gap.gap_main.fields.qcb_avg = 0;
	block_gap.gap_main.fields.qvp_avg = 0;


	// ----------------------------------------------------------------------- //
		
	revcounter++;
	block_gap.gap_main.revolution = revcounter;

	if(revcounter < in.data.lubrication_module.n_lubrication_revolutions)
	{
		// create the folders structure for the next revolution
		oss.str("");
		oss << "IF NOT EXIST output\\block\\rev." << revcounter 
				<< " (mkdir output\\block\\rev." << revcounter << ")";
		system(oss.str().c_str());
		
		oss.str("");
		oss << "IF NOT EXIST output\\block\\rev." << revcounter 
				<< "\\influgen (mkdir output\\block\\rev." << revcounter << "\\influgen)";

		if(dense_vtk || dense_light_vtk)
		{
			oss.str("");
			oss << "IF NOT EXIST output\\block\\rev." << revcounter 
					<< "\\dense_vtk (mkdir output\\block\\rev." << revcounter << "\\dense_vtk)";
			system(oss.str().c_str());
		}
	
		if(EHD_debug)
		{
			oss.str("");
			oss << "IF NOT EXIST output\\block\\rev." << revcounter 
					<< "\\EHD (mkdir output\\block\\rev." << revcounter << "\\EHD)";
			system(oss.str().c_str());
		}
	}
		

	clock_start = clock();	// reset the clock

}
// ------------------------------------------------------------------------- //
void cblock_ode::write_output()
{

	out.open("./output/block/block.txt",  fstream::out | fstream::app);

	out << scientific << t << "\t"
			<< scientific << t/T << "\t"
			<< scientific << (180/pi)*phi << "\t"
			// block position
			<< scientific << block_gap.gap_main.film.h1 << "\t"
			<< scientific << block_gap.gap_main.film.h2 << "\t"
			<< scientific << block_gap.gap_main.film.h3 << "\t"
			// squeeze
			<< scientific << block_gap.gap_main.squeeze.dhdt1 << "\t"
			<< scientific << block_gap.gap_main.squeeze.dhdt2 << "\t"
			<< scientific << block_gap.gap_main.squeeze.dhdt3 << "\t"
			// contact velocities
			<< scientific << block_gap.gap_main.squeeze.d_dhdt1 << "\t"
			<< scientific << block_gap.gap_main.squeeze.d_dhdt2 << "\t"
			<< scientific << block_gap.gap_main.squeeze.d_dhdt3 << "\t"
			// film thickness
			<< scientific << block_gap.gap_main.hmin() << "\t"
			<< scientific << block_gap.gap_main.hmax() << "\t"
			<< scientific << block_gap.gap_main.havg() << "\t"
			// external forces
			<< scientific << block_gap.gap_main.forces.MBx << "\t"
			<< scientific << block_gap.gap_main.forces.MBy << "\t"
			<< scientific << block_gap.gap_main.forces.FBz << "\t"
			// fluid forces
			<< scientific << block_gap.gap_main.forces.MfBx << "\t"
			<< scientific << block_gap.gap_main.forces.MfBy << "\t"
			<< scientific << block_gap.gap_main.forces.FfBz << "\t"
			// leakage
			<< scientific << output.leak << "\t"
			// Total losses
			<< scientific << output.Pl << "\t"
			// Leakage loss
			<< scientific << output.Pll << "\t"
			// Friction loss
			<< scientific << output.Plf << "\t"
			// friction forces from the other interfaces
			<< scientific << output.FTK << "\t"
			//<< scientific << block_gap.gap_main.forces.FTK.getf(phi) << "\t"
			<< scientific << block_gap.gap_main.forces.FTG.getf(phi)
			<< endl;
	
	out.close();
}
// ------------------------------------------------------------------------- //
void cblock_ode::testing()
{

	initialize();

	
	block_gap.gap_main.set_film_thickness();
	block_gap.gap_main.set_pressure();
	block_gap.gap_main.set_temperature();
	block_gap.reynolds.solve_EHD();
	block_gap.gap_main.set_film_thickness();
	

	/*

	block_gap.gap_main.set_film_thickness();
	block_gap.gap_main.set_pressure();
	block_gap.gap_main.set_temperature();
	block_gap.reynolds.solve_EHD();
	//block_gap.gap_main.write_light_vtk("phi.CB.vtk", 1e3);

	block_gap.gap_main.calc_mass_flow();
	block_gap.energy.solve();
	block_gap.gap_main.fields.qcb_prog = block_gap.gap_main.calc_cb_heatflux();
	block_gap.gap_main.fields.qvp_prog = block_gap.gap_main.calc_vp_heatflux();

	
	// thermoelastic solver
	block_gap.te_CB = new te_solver(in, "CB");

	// pointer to the interpolation grid
	interpl_grid* gap_str = &(block_gap.te_CB->gap);

	// apply the heat flux from gap
	block_gap.gap_main.apply_heat_flux(gap_str, "CB");	
	// solve for heat transfer
	block_gap.te_CB->solve_htr();
	// interpolation of temperature from structure to fluid
	interpolation str_to_fl_cb;
	str_to_fl_cb.cellsTocells(gap_str, &block_gap.gap_main.CBf);
	// update the surface temperature
	block_gap.gap_main.update_Tcb();	

	// solve for thermal deflection
	block_gap.te_CB->solve_thdfl();
	// interpolation of deformation from structure to fluid
	str_to_fl_cb.pointsTocells(gap_str, &block_gap.gap_main.CBf);
	// update the gap deformation
	block_gap.gap_main.update_dhcb_thermal();

	block_gap.te_CB->write_vtk("./");
	
	
	// thermoelastic solver
	block_gap.te_VP = new te_solver(in, "VP");

	// pointer to the interpolation grid
	gap_str = &(block_gap.te_VP->gap);

	// apply the heat flux from gap
	block_gap.gap_main.apply_heat_flux(gap_str, "VP");	
	// solve for heat transfer
	block_gap.te_VP->solve_htr();
	// interpolation of temperature from structure to fluid
	interpolation str_to_fl_vp;
	str_to_fl_vp.cellsTocells(gap_str, &block_gap.gap_main.VPf);
	// update the surface temperature
	block_gap.gap_main.update_Tvp();	

	// solve for thermal deflection
	block_gap.te_VP->solve_thdfl();
	// interpolation of deformation from structure to fluid
	str_to_fl_vp.pointsTocells(gap_str, &block_gap.gap_main.VPf);
	// update the gap deformation
	block_gap.gap_main.update_dhvp_thermal();

	block_gap.te_VP->write_vtk("./");
	
	block_gap.gap_main.set_film_thickness();
	block_gap.gap_main.write_vtk("gap_th.vtk", 1e3);
	
	

	//block_gap.write_cb_influgen("cb.txt");
	//block_gap.write_vp_influgen("vp.txt");

	*/

}
// ------------------------------------------------------------------------- //