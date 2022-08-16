// ------------------------------------------------------------------------- //
// -------------- class to handle CAD geometry in stl format --------------- //
// ------------------------------------------------------------------------- //

#	ifndef __stl__
#	define __stl__

#	include <iostream>
#	include <string>
#	include <vector>
#	include "../interpolation/point.h"


// class to handle a STL geometry
struct stl 
{

	// solid name
	std::string name;		
	// triangles nodes
	std::vector<point> nodes;	
	// triangles list
	std::vector<std::vector<int>> triangles;	

	// determine is a point is inside or outside a triangle (2D)
	int is_inside_2D(int triangleidx, const point& P);
		
public:

	// construct by STL geometry
	stl();
	~stl();
	
	// read from file
	void read(const std::string& STLfile, double scale = 1.0);
	// write vtk file for the stl
	void writeVTK();

};


#endif
