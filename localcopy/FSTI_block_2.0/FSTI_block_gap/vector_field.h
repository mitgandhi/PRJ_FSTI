# ifndef __vectorField__
# define __vectorField__

# include <vector>
# include "./gap_mesh.h"
# include "./scalar_field.h"
# include "./gap_vectors.h"


// forward declarations
//struct gap_mesh;
//class scalar_field;

class vector_field
{

public:

	// pointer to the mesh
	const gap_mesh* mesh;

	// dimensions
	int m;		// radial
	int n;		// circumferential
	int q;		// axial
	int mn;		// layer size
	int nq;		// inner outer boundary size
	int mnq;	// total size

	// internal field
	std::vector<cy_vector> in;
	
	// boundary field
	struct
	{	
		std::vector<cy_vector> inner;
		std::vector<cy_vector> outer;
		std::vector<cy_vector> bottom;
		std::vector<cy_vector> top;
	} bound;

	vector_field();
	vector_field(const gap_mesh* msh, cy_vector val = cy_vector(0,0,0));
	void initialize(const gap_mesh* msh, cy_vector val = cy_vector(0,0,0));
	vector_field(const vector_field& s);
	~vector_field();

	scalar_field r() const;
	scalar_field theta() const;
	scalar_field z() const;

	// operator overloading

	// operator =
	vector_field& operator=(const vector_field& rhs);
	vector_field& operator=(const cy_vector& rhs);

	// operator +
	friend vector_field operator+(const vector_field& l, const vector_field& r);
	friend vector_field operator+(const scalar_field&, const vector_field& r);
	friend vector_field operator+(const vector_field& l, const scalar_field& r);
	// operator -
	friend vector_field operator-(const vector_field& l, const vector_field& r);
	friend vector_field operator-(const scalar_field&, const vector_field& r);
	friend vector_field operator-(const vector_field& l, const scalar_field& r);
	// operator *
	friend vector_field operator*(const vector_field& l, const vector_field& r);
	friend vector_field operator*(const scalar_field& l, const vector_field& r);
	friend vector_field operator*(const vector_field& l, const scalar_field& r);
	// operator /
	friend vector_field operator/(const vector_field& l, const vector_field& r);
	friend vector_field operator/(const vector_field& l, const scalar_field& r);
};

# endif