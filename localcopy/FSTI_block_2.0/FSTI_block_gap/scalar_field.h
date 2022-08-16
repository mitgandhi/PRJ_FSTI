# ifndef __scalarField__
# define __scalarField__

# include "./gap_mesh.h"
# include <vector>

// forward declarations
struct gap_mesh;

class scalar_field
{

public:

	// pointer to the mesh
	const gap_mesh* mesh;

	// dimensions
	int m;		// radial
	int n;		// circumferential
	int q;		// axial
	int mn;		// layer size
	int nq;		// inner boundary size
	int mnq;	// total size

	// internal field
	std::vector<double>	in;
	
	// boundary field
	struct
	{	
		std::vector<double> inner;
		std::vector<double> outer;
		std::vector<double> bottom;
		std::vector<double> top;
	} bound;

	// constructors
	scalar_field();
	scalar_field(const gap_mesh* msh, double val = 0);
	scalar_field(const scalar_field& s);
	void initialize(const gap_mesh* msh, double val = 0);
	~scalar_field();

	// rotate the field (step > 0 -> clockwise)
	scalar_field cshift(int steps) const;

	
	// calculate the gradient of the scalar field in the film thickness direction
	double get_zgrad_ijk(int i, int j, int k) const;		// get z gradient for cell i,j,k
	double get_zgrad_bottom(int i, int j) const;	// get z gradient at the bottom wall
	double get_zgrad_top(int i, int j) const;			// get z gradient at the top wall

	double norm2() const;								// get the the squared norm
	double min() const;									// get the minimum value
	double max() const;									// get the maximum value
	double avg() const;									// get the average value
	double sum() const;									// get the internal field summation


	// operator overloading

	// operator =
	scalar_field& operator=(const scalar_field& rhs);
	scalar_field& operator=(double rhs);

	// operator +
	friend scalar_field operator+(const scalar_field& l, const scalar_field& r);
	friend scalar_field operator+(double l, const scalar_field& r);
	friend scalar_field operator+(const scalar_field& l, double r);
	// operator -
	friend scalar_field operator-(const scalar_field& l, const scalar_field& r);
	friend scalar_field operator-(double l, const scalar_field& r);
	friend scalar_field operator-(const scalar_field& l, double r);
	// operator *
	friend scalar_field operator*(const scalar_field& l, const scalar_field& r);
	friend scalar_field operator*(double l, const scalar_field& r);
	friend scalar_field operator*(const scalar_field& l, double r);
	// operator /
	friend scalar_field operator/(const scalar_field& l, const scalar_field& r);
	friend scalar_field operator/(const scalar_field& l, double r);

	// limit numeric values to [min max
	void limit(double min, double max);

	// save the scalar field
	void save(const char* filename);
	// load the scalar_field
	void load(const char* filename);


};

# endif