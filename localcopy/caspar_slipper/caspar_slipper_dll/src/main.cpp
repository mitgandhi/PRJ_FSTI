#include "../../caspar_input/caspar_input.h"
#ifdef _WIN32
#include "caspar_slipper_dll.h"
#endif
#include "main.h"
#include "CGapLog.h"
#include "../../fstiAuthorization/fsti_netlic_client.h"

#include "CGapODE.h"	//Gap Module

//Main Functions
void Gap(const bool Resume, const input & I);
void Help();
void CreateDirs();

CGapLog GapLog;			//The master log file

#ifdef _WIN32
//The PPL concurency events used by coupled caspar
Concurrency::event * coupled_caspar_local_event;
Concurrency::event * coupled_caspar_global_event;
bool coupled_caspar_simulation = false;
#endif

int slipper_standalone(int argc, char *argv[], const input & I)
{
	//Check if FSTI is authorized
	void * license = fsti_netlic::checkout_license(fsti_netlic::NODE_THEN_NETWORK, "slipper");

	if(argc >= 2)
	{
		if(strcmp(argv[1],"-gap") == 0)
		{
			bool Resume = false;
			if (argc == 3)
			{
				if(strcmp(argv[2],"-resume") == 0)
				{
					Resume = true;
				}
			}

			Gap(Resume, I);
		} else {
			Help();	
		}
	} else {
		//Go interactive mode
		
		//CASPAR intro banner
		GapLog << "********************************************************************************" << endl;
		GapLog << "*" << endl;
		GapLog << "*" << "\t" << "FSTI Gap Design" << endl;
		GapLog << "*" << "\t"  << "Slipper/Swashplate Lubrication Module" << endl;
		string datetime = __DATE__;
		datetime += " ";
		datetime += __TIME__;
		GapLog << "*" << "\t" << "Version: " << datetime << endl;
		GapLog << "*" << "\t" << "Version: " << "Run \"About FSTI Gap Design\"" << endl;
		GapLog << "*" << endl;
		GapLog << "*" << "\t" << "Purdue University" << endl;
		GapLog << "*" << "\t" << "Maha Fluid Power Research Center" << endl;
		GapLog << "*" << "\t" << "https://engineering.purdue.edu/Maha" << endl;
		GapLog << "*" << endl;
		GapLog << "*" << "\t" << "PROPRIETARY AND CONFIDENTIAL." << endl;
		GapLog << "*" << "\t" << "INFORMATION CONTAINED HEREIN SHALL NOT BE DECOMPILED, DISASSEMBLED," << endl;
		GapLog << "*" << "\t" << "DUPLICATED, OR DISCLOSED IN WHOLE OR IN PART FOR ANY PURPOSE." << endl;
		GapLog << "*" << "\t" << "UNAUTHORIZED USE IS STRICTLY PROHIBITED." << endl;
		GapLog << "*" << endl;
		GapLog << "********************************************************************************" << endl;
		GapLog << endl;

		//module selection loop
		bool quit = false;
		while(!quit)
		{
			GapLog << endl;
			GapLog << "Please select an option:" << endl;
			GapLog << endl;
			GapLog << "\t" << "2   Run Slipper Lubrication Module" << endl;
			
			GapLog << "\t" << "4   Quit FSTI" << endl;
			GapLog << endl;
			GapLog << "Enter Option Number: ";
			string input;
			int option;
			cin >> input;
			option = atoi(input.c_str());

			system("cls");
			
			GapLog << endl << "********************************************************************************" << endl;
			switch(option)
			{
			case 2:
				GapLog << "Running Gap Module..." << endl;
				Gap(false, I);
				break;
			case 4:
				GapLog << "Quitting FSTI..." << endl;
				quit = true;
				break;
			default:
				GapLog << "ERROR: Invalid Selection!" << endl;
				break;
			}
			GapLog << endl << "********************************************************************************" << endl;
			

		}
		GapLog << endl;

	}
	
	fsti_netlic::release_license(license);
	return 0;
}

void CreateDirs()
{
	//Hmm.. I don't want to use this for now
	//return;

#if OS == WINDOWS
	GapLog << "Creating Directories for Windows...";

	//system(string("IF EXIST .\\Logs (rd /S/Q .\\Logs )").c_str());	
	//system(string("IF EXIST .\\Outputs (rd /S/Q .\\Outputs )").c_str());	
	//system(string("IF EXIST .\\vtk (rd /S/Q .\\vtk)").c_str());	

	system(string("IF NOT EXIST .\\output (mkdir .\\output)").c_str());	
	system(string("IF NOT EXIST .\\output\\slipper (mkdir .\\output\\slipper)").c_str());	
	system(string("IF NOT EXIST .\\output\\slipper\\vtk (mkdir .\\output\\slipper\\vtk)").c_str());	
	system(string("IF NOT EXIST .\\output\\slipper\\matlab (mkdir .\\output\\slipper\\matlab)").c_str());	

		
	GapLog << " Finished!" << endl;
#endif

#if OS == LINUX
	GapLog << "Creating Directories for Linux..." << endl;

	system(string("rm -fr ./Logs").c_str());	
	system(string("rm -fr ./Outputs").c_str());	
	
	system(string("mkdir ./Logs").c_str());
	system(string("mkdir ./Outputs").c_str());

	GapLog << " Finished!" << endl;
#endif
}

void Help(void)
{
	GapLog << "Usage: [-pressure | -gap] [-resume]" << endl;
}

void Gap(const bool Resume, const input & I)
{
	if (!Resume)
	{
		CreateDirs();
	}
	
	caspar_input gapinput(I);

	CGapODE GapModule(&gapinput);
	GapModule.ODEmain(Resume);
}

#ifdef _WIN32
void run_slipper_gap(Concurrency::event & local_event, Concurrency::event & global_event, const input & I)
{
	//disable the cout becuase we are running coupled
	GapLog.enable_cout = false;

	//Check if FSTI is authorized
	void * license = fsti_netlic::checkout_license(fsti_netlic::NODE_THEN_NETWORK, "slipper");

	//Assign a pointer to the events to global extern variables
	coupled_caspar_local_event = &local_event;
	coupled_caspar_global_event = &global_event;

	//Set the coupled caspar bool to true
	coupled_caspar_simulation = true;

	Gap(false, I);

	//release the license
	fsti_netlic::release_license(license);
}
#endif

string slipper_version()
{
	string datetime = __DATE__;
	datetime += " ";
	datetime += __TIME__;

	string ver = "FSTI VERSION ";

	string version = ver + " " + datetime;
	return version;
}

