# ifndef __th_tetra_3dof__
# define __th_tetra_3dof__

# include "./mesh.h"
# include "./matrix.h"
# include "../FSTI_input/inputdata.h"  // material definition

class th_tetra_3dof
{
	// structural analysis can be limited to a mesh subset
	const tetra_sub_mesh& smsh;

	// get nodes coordinates (nidx is the global node index)
	double x(int nidx) const { return smsh.msh.nodes[gm.nodes[nidx]].x(); }
	double y(int nidx) const { return smsh.msh.nodes[gm.nodes[nidx]].y(); }
	double z(int nidx) const { return smsh.msh.nodes[gm.nodes[nidx]].z(); }
	
public:

	int id;												// element id (local to the submesh)
	const tetra& gm;							// element geometry
	static const int ndof = 3;		// dof for each node
	std::vector<int> gnum;				// global numbering
	const material* m;						// reference to the associated material
	matrix D;											// constitution matrix
	matrix B;											// strain displacement matrix

	th_tetra_3dof(int _id, const tetra_sub_mesh& _smsh);
	~th_tetra_3dof() {}

	matrix calc_B();		// calculate the strain-displacement matrix
	matrix calc_D();		// calculate the consitutive matrix
	matrix calc_K();		// calculate the stiffness matrix

};



# endif