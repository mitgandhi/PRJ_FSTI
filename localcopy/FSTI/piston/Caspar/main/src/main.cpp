#include "CGapODE.h"
#include "CGapInput.h"
#include "CGapUtils.h"
#include "CPressure.h"
#include "main.h"
#include "logger.h"
#include "../../caspar_input/input.h"
//#include <iostream>
//#include <iomanip>
#include <cmath>
#include "fsti_netlic_client.h"
#pragma once

//The PPL concurency events used by coupled caspar
Concurrency::event * coupled_caspar_local_event;
Concurrency::event * coupled_caspar_global_event;
bool coupled_caspar_simulation;
bool output_piston_window;
double caspar_autoconverge;
int min_revs;
bool force_balance_iterative;
double mixed_friction;
int hd_revs;
int hd_levels;

//classes objects declaration
CGapInput myGapInput;
sGapResult myGapResult;
CMesh myMesh;
CGapUtils myGapUtils;
CThermal myThermal;
CFEMThermal myFEMThermal;
input myinput;
//function prototypes
void StartPressureModule(void);
void StartGapModule(void);



int piston_standalone(int argc,char *argv[])
{
	//cout << argv[0] << "\n" << argv[1] << "\n";
	output_piston_window = 1;
	caspar_autoconverge = 0;
	min_revs = 0;
	force_balance_iterative = 0;
	mixed_friction = 0.0;
	hd_revs = 0;
	int i = 1;
	if(argc > 1){

		
		for(int a = 1;a<argc ;a++){

			if(!strcmp(argv[a],"-q"))
				output_piston_window = 0;
			else if(!strcmp(argv[a],"-p"))
				i = 0;
			else if(!strcmp(argv[a],"-a")){
				caspar_autoconverge = atof(argv[a+1]);
				a++;
			}
			else if(!strcmp(argv[a],"-m")){
				min_revs = atoi(argv[a+1]);
				a++;
			}
			else if(!strcmp(argv[a],"-i"))
				force_balance_iterative = 1;
			else if(!strcmp(argv[a],"-f")){
				a++;
				mixed_friction = atof(argv[a]);
			}
			else if(!strcmp(argv[a],"-d")){
				a++;
				hd_revs = atof(argv[a]);
				a++;
				hd_levels = atof(argv[a]);
				//cout << hd_revs << "\t" << hd_levels << endl;
			}
			else{
				cout << "Usage: fsti_piston.dll [Options]" << "\n";
				cout << "     -a #     Stops simulation if consecutive revolutions converge to within #%." << "\n";
				cout << "     -d # #   Begin HD grid refinement at revolution # with # levels of refinement.\n";
				cout << "     -f #     Viscosity in contact areas calculated with additional # fraction of contact stress.\n";
				cout << "     -i       Uses force balance iterative method." << "\n";
				cout << "     -m #     Runs at least # revolutions before stopping (used with autoconverge option)." << "\n";
				cout << "     -p       DEPRECATED - Runs piston pressure module." << "\n";
				cout << "     -q       No console output." << "\n";
				
				exit(1);
			}
		}
	}

	
	double pi = 4.0 * atan(1.0);

	
	if(i)
	{
		cout << "Starting the Gap Module!" << "\n";
		cout << "\n";
		coupled_caspar_simulation = 0;
		//create necessary folders for outputs
		char CurrentPath[_MAX_PATH];
		_getcwd(CurrentPath, _MAX_PATH);
		string timeroot = CurrentPath;
		string path;
		//main output folders
		path = timeroot + "\\output\\";					CreateDirectory(path.c_str(),NULL);
		path = timeroot + "\\output\\piston\\";			CreateDirectory(path.c_str(),NULL);
		path = timeroot + "\\output\\piston\\matlab\\";	CreateDirectory(path.c_str(),NULL);
		path = timeroot + "\\output\\piston\\vtk\\";	CreateDirectory(path.c_str(),NULL);
		//create log file
		Log.open("./piston_log.txt");
		StartGapModule();
	}
	else
	{
		cout << "You selected to run Pressure Module." << "\n";
		cout << "Warning: This pressure module is for internal research only." << "\n";
		cout << "It may not be approved for use by the Maha lab.  For best results" << "\n";
		cout << "use the pressure module included with the FSTI GUI." << "\n";
		cout << "\n";
		system("pause");
		StartPressureModule();
	}

	//system("PAUSE");
	Log << "Piston/Cylinder Completed Successfully!" << "\n";


	return 0;

};
void StartPressureModule(void)
{

	//read input files (dwm)
	//cout << "Reading input files... " << "\t";
	//myGapInput.readGeneral();
	//cout << "done!" << "\n";
	//cout << "\n";
	//cout << "\n";


	//output main simulation paramters to console
	cout << "**************************************************" << "\n";
	cout << "*           Main Simulation Parameters           *" << "\n";
	cout << "**************************************************" << "\n";
	cout << "\n";
	cout << "Revs. = " << myinput.data.lubrication_module.n_lubrication_revolutions << " [-]" << "\n";  //Use ONLY for internal research, this runs the Piston/Cylinder simplified pressure module.  For coupled simulation, there is a more detailed pressure module that should be run automatically!!!
	cout << "pHP = " << myinput.data.operating_conditions.HP*1e-5 << " [bar]" << "\n";
	cout << "pLP = " << myinput.data.operating_conditions.LP*1e-5 << " [bar]" << "\n";
	cout << "pCase = " << myinput.data.operating_conditions.pCase*1e-5 << " [bar]" << "\n";
	cout << "n = " << myinput.data.operating_conditions.speed * 30/PI << " [rpm]" << "\n";
	cout << "Beta = " << myinput.data.operating_conditions.beta * (180/PI) << " [deg]" << "\n";
	cout << "Gamma = " << myinput.data.geometry.gamma * (180/PI) << " [deg]" << "\n"; //dwm implement gamma angle.
	cout << "\n";
	cout << "\n";

	//system("PAUSE");
	

	/*/create necessary folders for outpurs
	char CurrentPath[_MAX_PATH];
	_getcwd(CurrentPath, _MAX_PATH);
	string timeroot = CurrentPath;
	string path;
	//main output folders
	path = timeroot + "\\outputs_piston\\";		CreateDirectory(path.c_str(),NULL); dwm */

	//launch calculations
	CPressure myPressure;
	myPressure.Pressure();

};
void StartGapModule(void)
{
	
	//read input files
	//Log << "Reading input files... " << "\n";

	//myGapInput.readGeneral();
	//myGapInput.readGeometry();
	//myGapInput.readBoundary();
	//myGapInput.readOptions();
	//myGapInput.readGrids();
	//input();
	//Log << "done!" << "\n";
	//Log << "\n";
	//Log << "\n";

	Log << "Checking Registration... " ;

	void * license = fsti_netlic::checkout_license(fsti_netlic::NETWORK, "piston");
	//for network licenses only use fsti_netlic::NETWORK

	Log << "Success!" << "\n";

	//Output Internal Use Header
	/*
	Log << "" << "\n";
	Log << "" << "\n";
	Log << "                             _" << "\n";
	Log << "    _._ _..._ .-',     _.._(`))" << "\n";
	Log << "   '-. `     '  /-._.-'    ',/" << "\n";
	Log << "      )         \\            '." << "\n";
	Log << "     / _    _    |  F S T I    \\ " << "\n";
	Log << "    |  a    a    /  Internal    |" << "\n";
	Log << "    \\   .-.         Version     ; " << "\n";
	Log << "     '-('' ).-'       ,'       ;" << "\n";
	Log << "        '-;           |      .'" << "\n";
	Log << "           \\           \\    /" << "\n";
	Log << "           | 7  .__  _.-\\   \\ " << "\n";
	Log << "           | |  |  ``/  /`  /" << "\n";
	Log << "          /,_|  |   /,_/   /" << "\n";
	Log << "             /,_/      '`-'" << "\n";
	Log << "Version: 3 - " << __DATE__ << " " << __TIME__ << "\n";
	Log << "-Autoconvergence Option" << "\n";
	Log << "-Adaptive Multithreading" << "\n";
	Log << "-Compatible with Influgen 2" << "\n";
	Log << "-Linear Force Balance Default" << "\n";
	Log << " -Iterative Force Balance Option" << "\n";
	Log << "-Mixed Friction Option" << "\n";
	Log << "-Reynolds boundary for pressure calculations." << "\n";
	Log << "" << "\n";
	Log << "" << "\n";
	*/
	//Output Common Header
	Log << "\n\n\n\n" << "\n";
	Log << "*******************************************************************" << "\n";
	Log << "FSTI Gap Design" << "\n";
	Log << "Piston/Cylinder Lubrication Module" << "\n";
	Log << "Version: See FSTI Gap Design Info" << "\n";
	Log << "\n";
	Log << "Purdue University" << "\n";
	Log << "Maha Fluid Power Research Center" << "\n";
	Log << "https://engineering.purdue.edu/Maha" << "\n";
	Log << "\n";
	Log << "PROPRIETARY AND CONFIDENTIAL." << "\n";
	Log << "INFORMATION CONTAINED HEREIN SHALL NOT BE DECOMPILED, DISASSEMBLED," << "\n";
	Log << "DUPLICATED, OR DISCLOSED IN WHOLE OR IN PART FOR ANY PURPOSE." << "\n";
	Log << "UNAUTHORIZED USE IS STRICTLY PROHIBITED." << "\n";
	Log << "*******************************************************************" << "\n";
	Log << "\n\n\n\n" << "\n";
	Log << "\n";

	//output main simulation paramters to console
	double tempdouble;
	Log << "**************************************************" << "\n";
	Log << "*           Main Simulation Parameters           *" << "\n";
	Log << "**************************************************" << "\n";
	Log << "\n";
	Log << "Revs. = " << myinput.data.lubrication_module.n_lubrication_revolutions << " [-]" << "\n";
	tempdouble = myinput.data.operating_conditions.HP*1e-5;
	Log << "pHP = " << tempdouble << " [bar]" << "\n";
	tempdouble = myinput.data.operating_conditions.LP*1e-5;
	Log << "pLP = " << tempdouble << " [bar]" << "\n";
	tempdouble = myinput.data.operating_conditions.pCase*1e-5;
	Log << "pCase = " << tempdouble << " [bar]" << "\n";
	tempdouble = myinput.data.operating_conditions.speed * 30/PI;
	Log << "n = " << tempdouble << " [rpm]" << "\n";
	tempdouble = myinput.data.operating_conditions.beta * (180/PI);
	Log << "Beta = " << tempdouble << " [deg]" << "\n";
	//lizhi check his input here:
	//tempdouble = myinput.data.options_piston.numeric.stgv * 1e3;
	//Log << "stgv = " << tempdouble << " [mm]" << "\n";
	//tempdouble = myinput.data.options_piston.numeric.wgv * 1e3;
	//Log << "wgv = " << tempdouble << " [mm]" << "\n";
	//tempdouble = myinput.data.options_piston.numeric.ngv;
	//Log << "ngv = " << tempdouble << " [-]" << "\n";
	//tempdouble = myinput.data.options_piston.numeric.spgv * 1e3;
	//Log << "spgv = " << tempdouble << " [mm]" << "\n";
	//tempdouble = myinput.data.options_piston.numeric.stploc * 1e3;
	//Log << "stploc = " << tempdouble << " [mm]" << "\n";
	//tempdouble = myinput.data.options_piston.numeric.stpdep * 1e6;
	//Log << "stpdep = " << tempdouble << " [micron]" << "\n";
	Log << "\n";
	Log << "\n";
	
	//Can't do force balance iterative without IM's
	if(myinput.data.options_piston.general.PressureDeformation == 0)
		force_balance_iterative = 0;

	//read optional piston macro geometry files
	if(myinput.data.options_piston.general.McrK)
	{
		Log << "Reading piston macrogeometry file... " << "\t";
		myGapInput.readMacroGeometryPiston();
		Log << "done!" << "\n";
		Log << " " << "\n";
	};


	//read optional cylinder macro geometry files
	if(myinput.data.options_piston.general.McrB)
	{
		//Log << "Cylinder Macrogeometry File is: " << myinput.data.options_piston.general.McrB_file << "\n";
		Log << "Reading cylinder macrogeometry file... " << "\t";
		myGapInput.readMacroGeometryCylinder();
		Log << "done!" << "\n";
		Log << " " << "\n";
	};


	//read optional dc pressure file (not required dwm)
	if(myinput.data.options_piston.general.ReadpFile)
	{
		Log << "Reading DC pressure file... " << "\t";
		myGapInput.readpFile();
		Log << "done!" << "\n";
		Log << " " << "\n";
	};

	//Inertial Relief for thermal Calculations
	if(myinput.data.options_piston.general.HeatTransfer || myinput.data.options_piston.general.ThermalDeformation)
	{
		myMesh.IR_c = myinput.data.thermal.block.IR;
		myMesh.IR_p = myinput.data.thermal.piston.IR;
		myMesh.d_tol_c = 1e-4;
		myMesh.d_tol_p = 1e-4;
	};

	//read main file for bodies mesh definition (not required dwm)
	//if(myGapInput.PistonOptions.PressureDeformation || myGapInput.PistonOptions.ThermalDeformation || myGapInput.PistonOptions.HeatTransfer)
	//{
	//	Log << "Reading mesh main.dat file... " << "\t";
	//	myMesh.ReadMeshMain();
	//	Log << "done!" << "\n";
	//	Log << " " << "\n";
	//};


	//read influence matrix sets
	/*if(myinput.data.options_piston.general.PressureDeformation)
	{
		//read surface coordinates
		Log << "Reading pressure mesh coordinates file... " << "\t";
		myGapInput.readBodySurfacexyzPressure();
		Log << "done!" << "\n";
		Log << " " << "\n";
		//read matrices
		Log << "Reading influence matrices... " << "\t";
		myGapInput.readInfluenceMatricesPistonCylinder(myinput.data.options_piston.general.IM_piston_path,myinput.data.options_piston.general.IM_bushing_path);
		Log << "done!" << "\n";
		Log << " " << "\n";
	};*/


	//Choose logical multigrid levels if user did not specify
	if((myinput.data.options_piston.fluid_grid.MG.nL == 0)&&(myinput.data.options_piston.general.ReynoldsMultiGrid)){
		while(myinput.data.options_piston.fluid_grid.MG.MG_M.size() > 1)
			myinput.data.options_piston.fluid_grid.MG.MG_M.pop_back();
		while(myinput.data.options_piston.fluid_grid.MG.MG_N.size() > 1)
			myinput.data.options_piston.fluid_grid.MG.MG_N.pop_back();
		double M = myinput.data.options_piston.fluid_grid.MG.MG_M[0];
		double N = myinput.data.options_piston.fluid_grid.MG.MG_N[0];
		while((N > 2) || (M > 2)){
			M /= 2;
			if(M<2)
				myinput.data.options_piston.fluid_grid.MG.MG_M.push_back(2);
			else
				myinput.data.options_piston.fluid_grid.MG.MG_M.push_back(ceil(M));
			N /= 2;
			if(N<2)
				myinput.data.options_piston.fluid_grid.MG.MG_N.push_back(2);
			else
				myinput.data.options_piston.fluid_grid.MG.MG_N.push_back(ceil(N));
		}
		myinput.data.options_piston.fluid_grid.MG.nL = myinput.data.options_piston.fluid_grid.MG.MG_M.size();
	}


	//make sure a boundary for gap_block is defined.
	bool found = false;
	for(int i = 0;i<myinput.data.thermal.block.dirichlet_bc.size();i++)
		if(myinput.data.thermal.block.dirichlet_bc[i].setname.compare("gap_block") == 0)
			found = true;
	for(int i = 0;i<myinput.data.thermal.block.mixed_bc.size();i++)
		if(myinput.data.thermal.block.mixed_bc[i].setname.compare("gap_block") == 0)
			found = true;
	for(int i = 0;i<myinput.data.thermal.block.neumann_bc.size();i++)
		if(myinput.data.thermal.block.neumann_bc[i].setname.compare("gap_block") == 0)
			found = true;
	if(!found){
		myinput.data.thermal.block.dirichlet_bc.resize(myinput.data.thermal.block.dirichlet_bc.size()+1);
		myinput.data.thermal.block.dirichlet_bc[myinput.data.thermal.block.dirichlet_bc.size()-1].setname = "gap_block";
		myinput.data.thermal.block.dirichlet_bc[myinput.data.thermal.block.dirichlet_bc.size()-1].Tp = myinput.data.operating_conditions.T_Leak;
	}

	//vtk output folder for surface temperature
	/*if(myGapInput.PistonOptions.HeatTransfer || myGapInput.PistonOptions.ThermalDeformation)
	{
		path = timeroot + "\\outputs_piston\\vtk\\";		CreateDirectory(path.c_str(),NULL);
	}*/
	//fem folder for thermal fem
	if(myinput.data.options_piston.general.ThermalDeformation || myinput.data.options_piston.general.HeatTransfer)
	{
		//create necessary folders for outputs
		char CurrentPath[_MAX_PATH];
		_getcwd(CurrentPath, _MAX_PATH);
		string timeroot = CurrentPath;
		string path;
		path = timeroot + "\\temp\\";							CreateDirectory(path.c_str(),NULL);
		path = timeroot + "\\temp\\piston\\";					CreateDirectory(path.c_str(),NULL);
		path = timeroot + "\\temp\\piston\\input\\";			CreateDirectory(path.c_str(),NULL);
		path = timeroot + "\\temp\\piston\\output\\";			CreateDirectory(path.c_str(),NULL);
		path = timeroot + "\\temp\\piston\\solver\\";			CreateDirectory(path.c_str(),NULL);
		path = timeroot + "\\temp\\piston\\piston\\";			CreateDirectory(path.c_str(),NULL);
		path = timeroot + "\\temp\\piston\\bushing\\";			CreateDirectory(path.c_str(),NULL);
	}

	//Initialize random shit
	srand(0);

	//launch calculations
	CGapODE myGapODE;
	myGapODE.ODEmain();

	Log << "Releasing Licence" << "\n";

	fsti_netlic::release_license(license);

};
int piston_GUI(Concurrency::event & local_event, Concurrency::event & global_event)//Called by the GUI to run the program
{
	//Assign a pointer to the events to global extern variables
	coupled_caspar_local_event = &local_event;
	coupled_caspar_global_event = &global_event;

	coupled_caspar_simulation = 1;
	//create necessary folders for outputs
	char CurrentPath[_MAX_PATH];
	_getcwd(CurrentPath, _MAX_PATH);
	string timeroot = CurrentPath;
	string path;
	//main output folders
	path = timeroot + "\\output\\";					CreateDirectory(path.c_str(),NULL);
	path = timeroot + "\\output\\piston\\";			CreateDirectory(path.c_str(),NULL);
	path = timeroot + "\\output\\piston\\matlab\\";	CreateDirectory(path.c_str(),NULL);
	path = timeroot + "\\output\\piston\\vtk\\";	CreateDirectory(path.c_str(),NULL);
	//create header for piston.txt file
	std::ofstream fout;
	fout.open("./output/piston/piston.txt",ios::app);
	fout << "%time\trevolution\tshaft_angle\te_1\te_2\te_3\te_4\te_1_velocity\te_2_velocity\te_3_velocity\te_4_velocity\t" 
		<< "variable_gap_length\tlout\tpos_A\tpos_B\tzRK\tStroke\tlA\tlB\tz_A\tz_B\tPiston_Velocity\tPiston_Leakage\t" 
		<< "Poiseuille_Leakage\tCouette_Leakage\tMTK\tMTKtan\tAxial_Friction\tFTKy_p\tFTKy_c\tCircumferential_Friction\tFTKx_p\t" 
		<< "FTKx_c\tFTG\tpDC\tFTby\tFTby_p\tFTby_c\tFTbx\tFTbx_p\tFTbx_c\tHigh_Pressure\tLow_Pressure\tF_1\tF_2\tF_3\tF_4\t" 
		<< "F_Fluid_1\tF_Fluid_2\tF_Fluid_3\tF_Fluid_4\tF_Contact_1\tF_Contact_2\tF_Contact_3\tF_Contact_4\tdF_1\tdF_2\t" 
		<< "dF_3\tdF_4\tTotal_Leakage\tTotal_Torque_Loss\tMFT\tTotal_Power_Loss\tTotal_Volumetric_Power_Loss\t" 
		<< "Power_Loss\tVolumetric_Power_Loss" << "\n";
	fout.close();
	fout.clear();
	//create log file
	Log.open("./piston_log.txt");
	StartGapModule();
	return 0;
}