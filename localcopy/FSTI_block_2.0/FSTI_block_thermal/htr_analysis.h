# ifndef __htr_analysis__
# define __htr_analysis__

# include "./mesh.h"
# include "../linear_solvers/PCG.h"
# include "./th_tetra_1dof.h"
# include "./th_load.h"
# include "./th_boundary.h"
# include "../FSTI_input/input.h"

// heat transfer analysis
class htr_analysis
{
	
	// reference to input structure
	const input& in;

	// body CB or VP
	const char* body;

	friend class th_load;
	friend class heat_flux;
	friend class convection;
	
	std::vector<th_boundary> bcs; 

	std::vector<th_tetra_1dof> elements;	
	gmm::col_matrix<gmm::wsvector<double>> K;
	std::vector<double> b;
	PCG solver;
	std::vector<int> gdof;
	int ndof;
	int ndof_eff;
		
	// returnt the bc associated with the face fid
	const th_boundary& get_face_bnd(int fid);

public:

	// heat flux from gap
	std::vector<double> qgap;

	// reference to the mesh
	const tetra_mesh& msh;	

	// list of thermal loads
	std::vector<th_load*> loads;
	// solution vector (T @ each node)
	std::vector<double> sol;
	
	// constructor
	htr_analysis(const input& in, const tetra_mesh& msh);
	// initialize the analysis
	void initialize(const char* body);
	// define the thermal stiffness matrix
	void get_K();
	// solve the linear system
	void solve();
	// clear the data structure
	void clear();
};

# endif