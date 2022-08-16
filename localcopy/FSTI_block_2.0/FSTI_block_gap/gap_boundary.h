// ------------------------------------------------------------------------- //
// ----- class to handle cylinder block and valve plate boundary shape ----- //
// ------------------------------------------------------------------------- //

# ifndef __gapboundary__
# define __gapboundary__

# include "./stl.h"
# include "../interpolation/polygon.h"
# include <vector>

class gap_boundary
{

	// -------------------------- private members ---------------------------- //

	// ------- geometric data ------- //

	// reference to the stl solid
	const stl* solid;
	// top and bottom gap plane in the solid
	double bottom, top;
	// surface nodes
	std::vector<point> nodes;
	// triangles (defined through the three nodes index)
	std::vector<std::vector<int>> triangles;
	// triangle neighbors
	std::vector<std::vector<int>> neigh;
	
	// -------------------------- private member functions ------------------- //
	
	// return true if is a top node
	bool is_top_node(const point& node) const;
	// select triangles with two nodes in the top plane
	void select_triangles();
	// define the neighbors for each triangle
	void define_neighbors();
	// define the boundaries of the gap surface and populates the boundaries vector
	void define_boundaries();

public:

	// --------------------------- public members  --------------------------- //

	// boundary definition
	std::vector<polygon> boundaries;

	// -------------------------- public member functions -------------------- //

	// constructor
	gap_boundary() {}
	// destructor
	~gap_boundary() {}
	// build the 
	void build(const stl* solid);
	// determine if a point is inside or outside the gap 
	bool is_inside(const point& P) const;

	
};

 # endif