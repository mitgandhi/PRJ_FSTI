// ------------------------------------------------------------------------- //
// --- this class handles the DC surface definition, which is used to ------ //
// ------------ calculate the external forces on the block ----------------- //
// ------------------------------------------------------------------------- //

# ifndef __dcsurfmesh__
# define __dcsurfmesh__

# include "../interpolation/point.h"
# include <vector>

class dc_surf_mesh
{

	bool available;
	std::vector<point> nodes;							// surface mesh nodes
	std::vector<std::vector<int>> faces;	// surface mesh elements (faces)
	std::vector<double> ai;								// face area
	std::vector<point> ni;								// face normals

	void analyze();												// get resultant force and point of application
	int read(const char* mesh_name);			// read in abaqus format

public:

	double Atot;													// total surface area
	double M0x, M0y, M0z;									// coordinates of the resultant force
	double F0x, F0y, F0z;									// components of the resultant force
	double xR, yR, zR;										// coordinates of the resultant force				

	dc_surf_mesh() {}
	~dc_surf_mesh() {}

	bool is_available() const { return available; }
	void initialize(const char* mesh_name);

};

# endif