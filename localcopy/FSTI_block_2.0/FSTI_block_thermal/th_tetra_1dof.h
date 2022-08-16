# ifndef __th_tetra_1dof__
# define __th_tetra_1dof__

# include "./mesh.h"
# include "./matrix.h"
# include "../FSTI_input/inputdata.h"  // material definition

class th_tetra_1dof
{
	const tetra_mesh& msh;

	// get nodes coordinates (nidx is the global node index)
	double x(int nidx) const { return msh.nodes[msh.elements[id].nodes[nidx]].x(); }
	double y(int nidx) const { return msh.nodes[msh.elements[id].nodes[nidx]].y(); }
	double z(int nidx) const { return msh.nodes[msh.elements[id].nodes[nidx]].z(); }
	
public:

	int id;												// element id
	static const int ndof = 1;		// dof for each node
	std::vector<int> gnum;				// global numbering
	const material* m;						// reference to the associated material

	th_tetra_1dof(int _id, const tetra_mesh& _msh);
	~th_tetra_1dof() {}

	// calculate the K matrix
	matrix calc_K();
	// calculate the integral of [N]'[N] over the element face indicated by side
	matrix calc_intNTN(int side) const;
	// calculate the integral of [N] over the element face indicated by side
	matrix calc_intN(int side) const;

};


# endif