# include "./influ_set.h"
# include "../FSTI_Block_dll/log.h"
# include <string>
# include <fstream>
# include <sstream>
# include <iostream>
# include <ctime>

extern gaplog Log;

using namespace std;

// ------------------------------------------------------------------------- //
influ_set::influ_set(const input& In, const gap_mesh* _mesh)
: in(In), msh(_mesh)
{

	// ------------------------ copy from caspar_input ----------------------- //

	EHD_CB = in.data.options_block.general.EHD_CB;
	EHD_VP = in.data.options_block.general.EHD_VP;
	IM_CB_path = in.data.options_block.general.IM_CB;
	IM_VP_path = in.data.options_block.general.IM_VP;

	// ----------------------------------------------------------------------- //

	CB.gap = 0;
	CB.DC = 0;
	CB.SPRING = 0;
	VP.gap = 0;
	VP.LP = 0;
	VP.HP = 0;

	if(EHD_CB)
	{
		// read the body structure surface geometry and build the kdtree
		string gapfile = (IM_CB_path + "/gap.txt");
		CBs.read(gapfile.c_str());
	}

	if(EHD_VP)
	{
		// read the body structure surface geometry and build the kdtree
		string gapfile = (IM_VP_path + "/gap.txt");
		VPs.read(gapfile.c_str());
	}

}
// ------------------------------------------------------------------------- //
influ_set::~influ_set()
{
	delete_matrices();
}
// ------------------------------------------------------------------------- //
void influ_set::read_matrices()
{

	// ------------------------ copy from caspar_input ----------------------- //

	int z = in.data.operating_conditions.npistons;

	// ----------------------------------------------------------------------- //

	ifstream infile;
	ostringstream oss;

	// -------------------------- cylinder block ----------------------------- //
	if(EHD_CB)
	{

		// allocate storage for the matrices
		
		CB.gap = new double*[CBs.ne()];
		for(int i=0; i<CBs.ne(); i++)
			CB.gap[i] = new double[CBs.nn()];
	
		CB.DC = new double*[z];
		for(int i=0; i<z; i++)
			CB.DC[i] = new double[CBs.nn()];
		
		CB.SPRING = new double[CBs.nn()];

		Log << "\nReading cylinder block matrices ...\n" << gaplog::endl;


		// ------------------------ read gap matrices -------------------------- //

		Log << "   * reading gap matrices ... ";

		oss.str("");
		oss << IM_CB_path << "/IM.gap.bin";
		infile.open(oss.str().c_str(), ios::in | ios::binary);
		if(!infile.is_open())
		{
			Log << "\ninflu_set::read_matrices(): ERROR file IM.gap.bin not found!" 
					 << gaplog::endl;
			exit(1);
		}

		for(int i=0; i<CBs.ne(); i++)
		{	
			infile.read
			(
				reinterpret_cast<char*>(CB.gap[i]),
				CBs.nn()*sizeof(double)
			);
		}

		infile.close();

		Log << "done!" << gaplog::endl;

		// ------------------------ read DC matrices --------------------------- //
		
		for(int i=0; i<z; i++)
		{
			Log << "   * reading matrix DC " << i+1 << " ... ";

			oss.str("");
			oss << IM_CB_path << "/IM.DC" << i+1 << ".txt";
			infile.clear();
			infile.open(oss.str().c_str());
			if(!infile.is_open())
			{
				Log << "\ninflu_set::read_matrices(): ERROR file IM.DC" << i+1 << ".txt not found!" 
						<< gaplog::endl;
				exit(1);
			}
			for(int j=0; j<CBs.nn(); j++)
				infile >> CB.DC[i][j];
			

			infile.close();
			
			Log << "done!" << gaplog::endl;
		}
	

		// ----------------------- read SPRING matrix -------------------------- //

		Log << "   * reading matrix SPRING ... ";

		oss.str("");
		infile.clear();
		infile.open(string(IM_CB_path + "/IM.SPRING.txt").c_str());
		if(!infile.is_open())
		{
			Log << " (WARNING file IM.SPRING.txt not found!) ";
			for(int i=0; i<CBs.nn(); i++)
			{
				CB.SPRING[i] = 0.0;
			}
			
		}
		for(int i=0; i<CBs.nn(); i++)
		{
			infile >> CB.SPRING[i];
		}
		
		Log << "done!" << gaplog::endl;
		
		infile.close();

		Log << "\ndone!" << gaplog::endl;
	}

	// --------------------------- valve plate ------------------------------- //
	
	if(EHD_VP)
	{

		// allocate storage for the matrices
		
		VP.gap = new double*[VPs.ne()];
		for(int i=0; i<VPs.ne(); i++)
			VP.gap[i] = new double[VPs.nn()];

		VP.LP = new double[VPs.nn()];
		VP.HP = new double[VPs.nn()];
		
		// ------------------------ read gap matrices -------------------------- //

		Log << "\nReading valve plate matrices ...\n" << gaplog::endl;

		Log << "   * reading gap matrices ... ";

		oss.str("");	
		oss << IM_VP_path << "/IM.gap.bin";
		infile.clear();
		infile.open(oss.str().c_str(), ios::in | ios::binary);
		if(!infile.is_open())
		{
			Log << "\ninflu_set::read_matrices(): ERROR file IM.gap.bin not found!" 
					 << gaplog::endl;
			exit(1);
		}

		for(int i=0; i<VPs.ne(); i++)
		{
			infile.read
			(
				reinterpret_cast<char*>
				(
					VP.gap[i]),
					VPs.nn()*sizeof(double
				)
			);
		}

		infile.close();

		Log << "done!" << gaplog::endl;

		// ------------------------- read LP matrix ---------------------------- //

		Log << "   * reading matrix LP ... ";

		oss.str("");
		infile.clear();
		infile.open(string(IM_VP_path + "/IM.LP.txt").c_str());

		if(!infile.is_open())
		{
			Log << "\ninflu_set::read_matrices(): ERROR file IM.LP.txt not found!" 
					 << gaplog::endl;
			exit(1);
		}
		else
		{
			for(int i=0; i<VPs.nn(); i++)
				infile >> VP.LP[i];
		}

		infile.close();

		Log << "done!" << gaplog::endl;

		// ------------------------- read HP matrix ---------------------------- //

		Log << "   * reading matrix HP ... ";

		oss.str("");
		infile.clear();
		infile.open(string(IM_VP_path + "/IM.HP.txt").c_str());
		
		if(!infile.is_open())
		{
			Log << "\ninflu_set::read_matrices(): ERROR file IM.HP.txt not found!" 
					 << gaplog::endl;
			exit(1);
		}
		else
		{
			for(int i=0; i<VPs.nn(); i++)
				infile >> VP.HP[i];
		}
		
		Log << "done!" << gaplog::endl;
		
		infile.close();

		Log << "\ndone!" << gaplog::endl;
	}

}
// ------------------------------------------------------------------------- //
void influ_set::delete_matrices()
{
	// ------------------------ copy from caspar_input ----------------------- //

	int z = in.data.operating_conditions.npistons;

	// ----------------------------------------------------------------------- //

	if(EHD_CB)
	{
		// delete gap matrices
		if(CB.gap != 0)
		{
			for(int i=0; i<CBs.ne(); i++)
				delete [] CB.gap[i];
			delete [] CB.gap;
			// set the pointer to zero
			CB.gap = 0;
		}
		// delete spring matrix
		if(CB.SPRING != 0)
		{
			delete [] CB.SPRING;
			// set the pointer to zero
			CB.SPRING = 0;
		}
		// delete displacement chambers matrices
		if(CB.DC != 0)
		{
			for(int i=0; i<z; i++)
				delete [] CB.DC[i];
			delete [] CB.DC;
			// set the pointer to zero
			CB.DC = 0;
		}
	}
	if(EHD_VP)
	{
		// delete gap matrices
		if(VP.gap != 0)
		{
			for(int i=0; i<VPs.ne(); i++)
				delete [] VP.gap[i];
			delete [] VP.gap;
			// set the pointer to zero
			VP.gap = 0;
		}
		if(VP.LP != 0)
		{
			if(VP.LP != 0)
				delete [] VP.LP;
			// set the pointer to zero
			VP.LP = 0;
		}
		if(VP.HP != 0)
		{
			if(VP.HP != 0)
				delete [] VP.HP;
			// set the pointer to zero
			VP.HP = 0;
		}
	}
}
