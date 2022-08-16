# ifndef __thdfl_analysis__
# define __thdfl_analysis__

# include <array>
# include <vector>
# include <map>
# include "../linear_solvers/PCG.h"
# include "../linear_solvers/PBiCG.h"
# include "./mesh.h"
# include "./th_tetra_3dof.h"
# include "./th_load.h"
# include "../FSTI_input/input.h"

// thermal deflection analysis over a mesh subset
class thdfl_analysis
{
	
	//friend class th_load;
	//friend class expansion;

	// body: CB or VP
	const char* body;
	
	const input& inp;

	// reference to the mesh
	const tetra_mesh& msh;	
		
	// list of finite elements
	std::vector<th_tetra_3dof> elements;	

	gmm::col_matrix<gmm::wsvector<double>> K;
	double K_avg;
	std::vector<double> b;
	PCG solver;
	std::vector<int> gdof;
	int ndof;
	int ndof_eff;

	std::vector<double> node_mass;
	std::vector<double> ir_loads;
	
public:

	// temperature @ each element
	const std::vector<double>* Te;


	// mesh subset where the deflection is calculated
	tetra_sub_mesh smsh;
	// thermal load
	std::vector<double> th_load;	
	// solution vector (T @ each node)
	std::vector<double> sol;

	// constructor
	thdfl_analysis(const input& _in, const tetra_mesh&	_msh);
	// initialize the analysis
	void initialize(const char* body);
	//
	void interpl_T();
	// calculate the thermal loads
	void calc_thloads();
	// inertial loads
	int calc_inertial_load();
	// apply intertia relief
	int apply_ir();
	// define the thermal stiffness matrix
	void get_K();
	// solve the linear system
	void solve();

};

# endif