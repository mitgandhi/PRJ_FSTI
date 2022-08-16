// ------------------------------------------------------------------------- //
// -------------------- class to handle gap fluid mesh --------------------- //
// ------------------------------------------------------------------------- //

# ifndef __gapmesh__
# define __gapmesh__

# include "../interpolation/point.h"
# include "../interpolation/polygon.h"
# include "./gap_boundary.h"
# include "./scalar_field.h"
# include "../FSTI_input/input.h"
# include "./body_type.h"
# include <vector>

// forward class declarations
class scalar_field;

// labelling used to distinguish fluid elements to boundaries
enum
{
	FLUID = 0,
	LP = 100,
	HP = 200,
	CASE = 300
};

// fluid element
struct gap_elm
{
	double x;				// x coord of the centroid
	double y;				// x coord of the centroid
	double r;				// radius of the centeroid
	double theta;		// centroid angular position
	double z;				// centroid axial position
	double dr;			// radial dimension
	double dtheta;	// circumferential dimension
	double dz;			// axial dimension
	double A;				// element area
	double V;				// element volume
	int nds[8];			// nodes defining the element 2D elm just first 4

	int ty;					// element type FLUID, HP, LP, ...
	int fluid_id;		// element id, considering just fluid domain

	int w;					// west neighbor
	int e;					// east neighbor
	int s;					// south neighbor
	int n;					// north neighbor
	int b;					// bottom neighbor
	int t;					// top neighbor

};

struct gap_mesh
{
	
	// -------------------- geometry ------------------------- //

	int dim;								// mesh dimension 2D or 3D
	
	double r_in;						// gap inner radius
	double r_openin;				// kidney opening inner radius
	double r_openout;				// kidney opening outer radius
	double r_groovein;			// groove inner radius
	double r_grooveout;			// groove outer radius
	double r_out;						// gap outer radius
	bool outerbearing;			// outer bearing (hydrodynamic bearing) present or not
	
	int rotation_steps;			// rotation steps

	stl cb_stl;							// stl geometry for the block gap
	stl vp_stl;							// stl geometry for the valve plate gap
	gap_boundary cb_bound;	// cylinder block geometry
	gap_boundary vp_bound;	// valve plate geometry
		
	// -------------------- elements ------------------------- //

	// radial sectors in the mesh
	static const int sectors = 5;
		
	int m;						// elements in radial direction
	int mi[sectors];	// elements in radial direction for each sector
										// m[0]: from r_in to r_opein
										// m[1]: from r_openin to r_openout
										// m[2]: from r_openout to r_groovein
										// m[3]: from r_groovein to r_grooveout
										// m[4]: from r_grooveout to r_out
	int n;						// elements in circumferential direction
	int q;						// elements in axial direction
	
	int mn;						// number of elements in one layer (nominal size, no openings)
	int mn_f;					// actual number of fluid elements in one layer
	int mnq;					// total number of elements (nominal size, no openings)
	int mnq_f;				// actual total number of fluid elements
	
	double dr[sectors];			// elements radial dimension
	double dtheta;					// elements circumferential dimension
	std::vector<double> dz;	// dz depends on the film thickness, size is mn

	std::vector<gap_elm> elements;	// element list
	std::vector<point> nodes;				// node list
	int cb_fluid;										// number of fluid elements on the block surface
	std::vector<int> cb_elms;				// element definition for cylinder block
	std::vector<int> cb_elms_0;			// element definition for cylinder block (fixed to zero position)
	int vp_fluid;										// number of fluid elements on the valve plate surface
	std::vector<int> vp_elms;				// element definition for valve plate

	// cylinder block gap surface, as list of polygons
	struct
	{
		std::vector<polygon> elms;
		std::vector<polygon> elms_0;	// fixed to zero position
		std::vector<double> ai;
	} cb_gap_surf;
	
	// valve plate gap surface, as list of polygons
	struct
	{
		std::vector<polygon> elms;
		std::vector<double> ai;
	} vp_gap_surf;
	

	// -------------------- member functions -------------------- //

	gap_mesh() {}																						// default constructor
	~gap_mesh() {}																					// default destructor
	void build(const input& in, int dim);										// build the mesh
	void define_elements_type();														// define the element types
	void define_gap_surface(body_type which);								// define the two gap surfaces
	int	rotate(double angle_rad);														// rotate the mesh by angle in radiants
	void update_thickness																		// update the mesh with the current film thickness
	(
		const scalar_field& top,															// block surface
		const scalar_field& bottom														// valve plate surface
	);
	void write_vtk(const char* file, double scale = 1.0);		// write the mesh in vtk file format
	void write_gap_surface_vtk(body_type which);						// write the gap surfaces

};

# endif