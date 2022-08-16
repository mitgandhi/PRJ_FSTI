# ifndef __influSet__
# define __influSet__

# include <vector>
# include <string>
# include "../FSTI_input/input.h"
# include "./gap_mesh.h"
# include "./scalar_field.h"
# include "../interpolation/interpl_grid.h"

// this class read and stores and delete the influence matrices

class influ_set
{

	const input& in;
	const gap_mesh* msh;

	bool EHD_CB;
	bool EHD_VP;
	std::string IM_CB_path;
	std::string IM_VP_path;	
	
public:

	interpl_grid CBs;			// structure interpolation grid for CB
	interpl_grid VPs;			// structure interpolation grid for VP

	// block influ set
	struct
	{
		double** gap;			// gap 
		double** DC;			// displacement chambers
		double* SPRING;		// spring
	} CB;
	// valve plate influ set
	struct
	{
		double** gap;		// gap
		double* LP;			// low pressure port
		double* LPB;		// beneath the valve plate, low pressure side
		double* HP;			// high pressure port
		double HPB;			// beneath the valve plate, high pressure side
	} VP;

	// ----------------------------- public functions ------------------------ //
	
	influ_set(const input& input,	const gap_mesh* mesh);	// constructor
	~influ_set();																					// destructor
	void read_matrices();																	// read influence matrices
	void delete_matrices();																// delete influence matrices
	
};

# endif