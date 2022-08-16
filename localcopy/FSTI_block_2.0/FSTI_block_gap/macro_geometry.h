# ifndef __macrogeometry__
# define __macrogeometry__

# include "./scalar_field.h"
# include "../interpolation/interpl_grid.h"
# include "../FSTI_input/input.h"
# include "./block_gap_main.h"

class cblock_gap_main;	// forward declaration

class macro_geometry
{

protected:

	scalar_field& h;	
	scalar_field*	dh; 

public:
	
	macro_geometry(scalar_field& h);
	virtual ~macro_geometry();

	void apply(int rotation_steps = 0);
	const scalar_field* getdh() const { return dh; }
	virtual void read(const char* filename);

};

// waved surface
class macro_waved : public macro_geometry
{
	
	double f;
	double A;
	double offset;

	void read(const char* filename);

public:

	macro_waved(scalar_field& h, const char* filename);

};

// axisymmetric radial profile
class macro_axisymmetric : public macro_geometry
{

	std::vector<double> ri;
	std::vector<double> dhi;

	void read(const char* filename);

public:

	macro_axisymmetric(scalar_field& h, const char* filename);

};

// two different spherical diameter for block and valve plate
class macro_spherical : public macro_geometry
{
	
	// original spherical diameter
	double D;
	// modified spherical diameter
	double Dm;

	void read(const char* filename);

public:

	macro_spherical(scalar_field& h, const char* filename);

};

// used to apply a fixed deformation through the simulation
// the deformation is provided through a grid (see interpl_grid)
class macro_grid : public macro_geometry
{
	
	string grid_file;
	string data_at;
	string body;

	// grid that specify the deformation
	interpl_grid from;

	void read(const char* filename);

public:

	macro_grid(scalar_field& h, const char* filename);
};

// use an stl geometry to define notches of a given depth
class macro_stl : public macro_geometry
{
	
	string stl_file;
	string body_type;
	double notch_depth;

	
	stl body;								// stl geometry for the block gap
	gap_boundary boundary;	// cylinder block geometry

	void read(const char* filename);

public:

	macro_stl(scalar_field& h, const char* filename);

};

// the deformation is given using fluid grid indicies
class macro_fluid_grid : public macro_geometry
{
	
	vector<int> i;
	vector<int> j;
	vector<double> macro_height;

	void read(const char* filename);

public:

	macro_fluid_grid(scalar_field& h, const char* filename);

};

# endif