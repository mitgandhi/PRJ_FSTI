# ifndef __interpl_grid__
# define __interpl_grid__

# include <vector>
# include "../interpolation/point.h"

struct interpl_grid
{
	
	// ---------------------------- public members --------------------------- //

	int ELM_NDS;																// nodes per element
	std::vector<point> nodes;										// list of nodes
	std::vector<std::vector<int>> elements;			// list of elements
	std::vector<bool> active;										// list of active elements

	std::vector<double> points_data;						// field value @ grid nodes
	std::vector<double> cells_data;							// field value @ element center
	
	// ------------------------ public member functions ---------------------- //
	
	interpl_grid();																		// default constructor
	interpl_grid(int ne, int nn, int elmnds);					// construct with given size
	void initialize(int ne, int nn, int elmnds);			// initialize an empty grid
	int ne() const { return elements.size(); }				// return the number of elements
	int nn() const { return nodes.size(); }						// return the number of nodes
	int isInside(int elmidx, const point& P) const;		// point is inside/outside

	void cellsToPoints();															// map cells data into points data
	void read(const char* name);											// read from grid file
	void write(const char* name) const;								// write grid file
	void writeVTK(const char* name) const;						// write vtk grid
	void copy(interpl_grid& newgrid) const;						// copy from another grid object

};

# endif