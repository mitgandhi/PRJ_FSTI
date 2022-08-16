// ------------------------------------------------------------------------- //
// ----------------- Matrix class written by Marco Zecchi ------------------ //
// ------------------------------------------------------------------------- //
//																																					 //
// This matrix class is very simple and not optimized in speed!							 //
// It can be used to work with small, dense matrices.												 //
//																																					 //
// ------------------------------------------------------------------------- //


# ifndef __matrix__
# define __matrix__

#include <iostream>

class matrix
{

protected:

	int _m, _n;
	double** val;

public:

	// default constructor
	matrix() { val = 0; _m = _n = 0; }
	// square matrix
	matrix(int sz, double diag_val = 0.0);
	// non square matrix
	matrix(int m, int n);
	// construc by providing size and values
	matrix(int m, int n, double d[]);
	// construct by components (default to column vector)
	matrix(double x, double y, double z);
	
	// copy constructor
	matrix(const matrix& M);
	// operator=
	matrix& operator=(const matrix& M);
	matrix& operator=(double v);
	
	// destructor
	~matrix();

	// resize the matrix
	void resize(int m, int n);

	// return number of rows
	int m() const { return _m; }
	// return number of columns
	int n() const { return _n; }
	// return element by position
	double operator()(int idx) const { return val[(idx/_n)][idx%_n]; }
	// return reference to element by position
	double& operator()(int idx) { return val[(idx/_n)][idx%_n]; }
	// operator [] to access with the [][] sintax
	double* operator[](int row) { return val[row]; }
	// operator [] to access with the [][] sintax, const version
	const double* operator[](int row) const { return val[row]; }

	// matrix-matrix multiplication
	friend matrix operator*(const matrix& left, const matrix& right);
	// matrix-scalar multiplication
	friend matrix operator*(const matrix& left, double v);
	// scalar-matrix multiplication
	friend matrix operator*(double v, const matrix& right);
	// matrix-scalar division
	friend matrix operator/(const matrix& left, double v);

	// matrix-matrix sum
	friend matrix operator+(const matrix& left, const matrix& right);
	// matrix-scalar sum
	friend matrix operator+(const matrix& left, double v);
	// scalar-matrix sum
	friend matrix operator+(double v, const matrix& right);
	// matrix-matrix subtraction
	friend matrix operator-(const matrix& left, const matrix& right);
	// matrix-scalar subtraction
	friend matrix operator-(const matrix& left, double v);
	// scalar-matrix subctration
	friend matrix operator-(double v, const matrix& right);

	friend matrix cross (const matrix& left, const matrix& right);
	
	// get a matrix minor 
	matrix getMinor(int row, int col);
	// get the matrix determinant
	double det();
	// get the matrix transpose
	matrix T();
	// get the matrix inverse
	matrix inv();

	// print
	friend std::ostream& operator<<(std::ostream& os, const matrix& t); 
	
};

# endif