# ifndef __interpolation__
# define __interpolation__

# include "./interpl_grid.h"
# include <ANN/ANN.h>

class interpolation
{

	// ---------------------------- private members -------------------------- //

	// kd tree of interpl_grid "from" element centers
	struct
	{
		static const int dim = 2;		// dimension
		double eps;									// error bound
		ANNpointArray dataPts;			// data points
		ANNpoint queryPt;						// query point
		ANNidxArray nnIdx;					// near neighbor indices
		ANNdistArray dists;					// near neighbor distances
		ANNkd_tree* kdTree;					// search structure
	} ANN;

	// data to interpolate from cells to cells
	struct
	{
		bool ready;																			// flag to understand is elms_idx and elms_weights are already defined
		std::vector<std::vector<int>> elms_idx;					// index of elements of interpl_grid "from" intersecting interpl_grid "to"
		std::vector<std::vector<double>> elms_weights;	// values of intersection areas of interpl_grid "from" for each element of interpl_grid "to"
		std::vector<double> elms_weight;								
	} cells2cells;

	// data to interpolate from nodes to cells
	struct
	{
		bool ready;										// flag to understand is elms_idx and nds_weights are already defined
		std::vector<int> elm_idx;			// index of "from" elements containing each element of interpl_grid "to"
		std::vector<point> centers;		// interpl_grid "to" centers
	} points2cells;

	// ------------------------ private member functions --------------------- //

	void define_cells2cells(const interpl_grid* from, const interpl_grid* to);	// initialize the cells2cells structure
	void define_points2cells(const interpl_grid* from, const interpl_grid* to);	// initialize the points2cells structure

public:

	// ------------------------ public member functions ---------------------- //

	interpolation();																									// constructor
	~interpolation();																									// destructor
	void clear();																											// clear the structures and kdtree
	void buildTree(const interpl_grid* from);													// build kdtree from interpl_grid from
	void cellsTocells(const interpl_grid* from, interpl_grid* to);		// interpolate cells data from one interpl_grid to onother interpl_grid
	void pointsTocells(const interpl_grid* from, interpl_grid* to);		// interpolate points data from one interpl_grid to cell data of onother interpl_grid

};

# endif