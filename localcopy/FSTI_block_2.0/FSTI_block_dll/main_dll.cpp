# include "../FSTI_block_dll/FSTI_block_dll.h"
# include "../FSTI_input/input.h"
# include "../FSTI_block_gap/block_ode.h"
# include "../FSTI_block_gap/block_gap.h"
# include "./log.h"
# include "./fsti_netlic_client.h"

# include <iostream>

using namespace std;

gaplog Log("./block_log.txt");

//The PPL concurency events used by coupled caspar
Concurrency::event * coupled_caspar_local_event;
Concurrency::event * coupled_caspar_global_event;
bool coupled_caspar_simulation = false;

// ------------------------------------------------------------------------- //

void testing()
{
	input in;
	
	cblock_ode block_ode(in);

	block_ode.testing();

}
// ------------------------------------------------------------------------- //
int run_block_standalone(int argc, char* argv[])
{

	if(argc < 2)
	{
		Log << "Usage CasparBlock -<option>" << gaplog::endl;
		Log << "\n\n";
		Log << "   -gap (run gap calculation)" << gaplog::endl;
		Log << "   -testing (development testing)" << gaplog::endl;
		Log << "   -pModule (run the pressure module)" << gaplog::endl;
		Log << "   -getExtLoads (calculate the external loads)" << gaplog::endl;
		Log << "   -getKinematics (calculate the piston position over one shaft rev)" << gaplog::endl;
		Log << "   -getOilProperties (write a table of the oil properties)" << gaplog::endl; 
		Log << "   -checkGap (check gap geometry and boundary definition)" 
				<< gaplog::endl;
		Log << "   -checkEHD (check influece matrices and EHD calculation)" << gaplog::endl;
		Log << "   -checkThermal (check thermal inputs and thermal solution)" << gaplog::endl;
		Log << "   -writeInflugenInputs" << gaplog::endl;
		Log << "\n\n";

		return 0;
	}
	else if(string(argv[1]).compare("-getExtLoads") == 0)
	{
		
		// check for the authorization //
		void * lic = fsti_netlic::checkout_license(fsti_netlic::NETWORK, "block");
		// --------------------------- //

		input in;
		cblock_gap block_gap(in);
		block_gap.get_external_loads();
		
		// --------------------------- //
		fsti_netlic::release_license(lic);
		return 0;
	}
	else if(string(argv[1]).compare("-getBalanceFactor") == 0)
	{
		
		// check for the authorization //
		void * lic = fsti_netlic::checkout_license(fsti_netlic::NETWORK, "block");
		// --------------------------- //

		input in;
		cblock_gap block_gap(in);
		block_gap.get_balance_factors();

		// --------------------------- //
		fsti_netlic::release_license(lic);
		return 0;
	}
	else if(string(argv[1]).compare("-getKinematics") == 0)
	{
		
		// check for the authorization //
		void * lic = fsti_netlic::checkout_license(fsti_netlic::NETWORK, "block");
		// --------------------------- //

		input in;
		cblock_gap block_gap(in);
		block_gap.get_piston_kinematics();


		// --------------------------- //
		fsti_netlic::release_license(lic);
		return 0;
	}
	else if(string(argv[1]).compare("-getOilProperties") == 0)
	{
		
		// check for the authorization //
		void * lic = fsti_netlic::checkout_license(fsti_netlic::NETWORK, "block");
		// --------------------------- //

		input in;
		cblock_gap block_gap(in);
		block_gap.get_oil_properties();
		
		// --------------------------- //
		fsti_netlic::release_license(lic);
		return 0;
	}
	else if(string(argv[1]).compare("-checkGap") == 0)
	{
		
		// check for the authorization //
		void * lic = fsti_netlic::checkout_license(fsti_netlic::NETWORK, "block");
		// --------------------------- //

		input in;
		cblock_gap block_gap(in);
		block_gap.check_gap();
		
		// --------------------------- //
		fsti_netlic::release_license(lic);
		return 0;
	}
	else if(string(argv[1]).compare("-checkEHD") == 0)
	{

		// check for the authorization //
		void * lic = fsti_netlic::checkout_license(fsti_netlic::NETWORK, "block");
		// --------------------------- //

		input in;
		cblock_gap block_gap(in);
		block_gap.check_EHD();
		
		// --------------------------- //
		fsti_netlic::release_license(lic);
		return 0;
	}
	else if(string(argv[1]).compare("-checkThermal") == 0)
	{
		// check for the authorization //
		void * lic = fsti_netlic::checkout_license(fsti_netlic::NETWORK, "block");
		// --------------------------- //

		input in;
		cblock_gap block_gap(in);
		block_gap.check_thermal();
		
		// --------------------------- //
		fsti_netlic::release_license(lic);
		return 0;
	}
	else if(string(argv[1]).compare("-writeVTKblock") == 0)
	{
		// check for the authorization //
		void * lic = fsti_netlic::checkout_license(fsti_netlic::NETWORK, "block");
		// --------------------------- //

		input in;
		cblock_gap block_gap(in);
		block_gap.writeVTK_block();
		
		// --------------------------- //
		fsti_netlic::release_license(lic);
		return 0;
	}
	else if(string(argv[1]).compare("-gap") == 0)
	{
		
		// check for the authorization //
		void * lic = fsti_netlic::checkout_license(fsti_netlic::NETWORK, "block");
		// --------------------------- //

		input in;

		if(in.data.lubrication_module.solve_block == 0)
		{
			cout << "Please enable block calculation in lubrication_module.txt" << endl;
			system("pause");
			exit(1);
		}

		cblock_ode block_ode(in);
		block_ode.initialize();
		block_ode.run();

		// --------------------------- //
		fsti_netlic::release_license(lic);
		return 0;
	}
	
	else if(string(argv[1]).compare("-writeInflugenInputs") == 0)
	{
		input in;
		cblock_gap block_gap(in);
		block_gap.writeInflugenInputs();
		return 0;
	}
	else if(string(argv[1]).compare("-writeVP") == 0)
	{
		input in;
		cblock_gap block_gap(in);
		block_gap.write_vp_pressure_field();
		return 0;
	}
	else 
	{
		if(string(argv[1]).compare("-testing") == 0)
		{
			testing();
			return 0;
		}
	}
}
// ------------------------------------------------------------------------- //
void run_block_gap(Concurrency::event & local_event, Concurrency::event & global_event, const input& in)
{
	// authorization //

	// check for the authorization //
	//void * lic = fsti_netlic::checkout_license(fsti_netlic::NODE_THEN_NETWORK, "block");
	
	// ------------- //

	coupled_caspar_simulation = true;

	coupled_caspar_local_event = &local_event;
	coupled_caspar_global_event = &global_event;

	cblock_ode block_ode(in);
	block_ode.initialize();
	block_ode.run();

}
// ------------------------------------------------------------------------- //