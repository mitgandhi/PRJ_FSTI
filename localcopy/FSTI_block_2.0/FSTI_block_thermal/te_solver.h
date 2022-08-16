# ifndef __te_solver__
# define __te_solver__

# include <vector>
# include "../FSTI_input/input.h"
# include "./mesh.h"
# include "./htr_analysis.h"
# include "./thdfl_analysis.h"

# include "../interpolation/interpl_grid.h"


class te_solver
{
	// body: CB or VP
	const char* body;

	// reference to input structure
	const input& in;
	
	// mesh structure
	tetra_mesh msh;
	
	// heat transfer analysis object
	htr_analysis htr;
	// structural analysis object
	thdfl_analysis thdfl;

	// temperature @ each element
	std::vector<double> Te;
	// temperature @ each node
	std::vector<double> Tn;
	// displacement @ each node
	std::vector<double> displ;

	// define the interpolation grid for the gap surface
	void define_gap_interpl();

public:

	// interpolation grid for the gap surface
	interpl_grid gap;	

	// constructor
	te_solver(const input& in, const char* body);
	// solve for the heat transfer analysis
	void solve_htr();
	// solve for the thermal deflection
	void solve_thdfl();
	// write vtk output
	void write_vtk(const char* path);
	// write face set in vtk format
	void write_fset_vtk(const char* path, const char* set);


};

# endif