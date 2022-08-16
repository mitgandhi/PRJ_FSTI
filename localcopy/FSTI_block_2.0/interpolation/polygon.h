// ------------------------------------------------------------------------- //
// -- simple polygon class to wrap the gpc library (which is a C library) -- //
//																																					 //
//				this class is used to interpolate from the gap grid to a					 //
//													trianglulated surface														 //
// ------------------------------------------------------------------------- //

# ifndef __polygon__
# define __polygon__

# include <vector>
# include "./point.h"

extern "C" 
{
	# include "gpc.h"
}

struct polygon
{
	// number of vertices
	int nv;	
	// list of vertices
	std::vector<point> vertices;

	// constructor
	polygon(int nv = 0);
	// destructor
	~polygon();
	// operator =
	polygon& operator=(const polygon& p);

	// return the intersection with another polygon
	// GPC_DIFF			Difference
  // GPC_INT			Intersection
  // GPC_XOR,     Exclusive or
  // GPC_UNION		Union
	std::vector<polygon> clip(const polygon& p, gpc_op = GPC_INT);
	
	// return polygon area
	double area();
	// determine if a point lies inside or outside the polygon
	bool is_inside(const point& P);
	// get the polygon center
	point center();
	// print the polygon to the standard output
	void print();
	// write the polygon into a vtk file
	void write_vtk(const char* name);
	// write a polygon surface defined by a vector of polys into a vtk file
	static void write_vtk(const char* filename, const std::vector<polygon>& surf);
};

# endif