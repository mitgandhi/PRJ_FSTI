# include "matrix.h"
# include "../FSTI_Block_dll/log.h"
# include <iostream>
# include <cmath>

using namespace std;

extern class gaplog Log;

// ------------------------------------------------------------------------- //
matrix::matrix(int m, int n) : _m(m), _n(n)
{
	val = new double*[_m];
	for(int i=0; i<_m; i++)
		val[i] = new double[_n];

	for(int i=0; i<_m; i++)
		for(int j=0; j<_n; j++)
			val[i][j] = 0;

}
// ------------------------------------------------------------------------- //
matrix::matrix(int m, int n, double d[]) : _m(m), _n(n)
{
	val = new double*[_m];
	for(int i=0; i<_m; i++)
		val[i] = new double[_n];

	for(int i=0; i<_m; i++)
		for(int j=0; j<_n; j++)
			val[i][j] = d[_n*i + j];
}
// ------------------------------------------------------------------------- //
matrix::matrix(double x, double y, double z) : _m(3), _n(1)
{
	val = new double*[_m];
	for(int i=0; i<_m; i++)
		val[i] = new double[_n];

	val[0][0] = x, val[1][0] = y, val[2][0] = z;
}
// ------------------------------------------------------------------------- //
matrix::matrix(int sz, double diag_val) : _m(sz), _n(sz)
{
	val = new double*[_m];
	for(int i=0; i<_m; i++)
		val[i] = new double[_n];

	for(int i=0; i<_m; i++)
	{
		for(int j=0; j<_n; j++)
		{
			if(i == j)
				val[i][j] = diag_val;
			else
				val[i][j] = 0.0;
		}
	}
}
// ------------------------------------------------------------------------- //
matrix::matrix(const matrix& M)
{

	_m = M.m();
	_n = M.n();

	val = new double*[_m];
	for(int i=0; i<_m; i++)
		val[i] = new double[_n];

	for(int i=0; i<_m; i++)
		for(int j=0; j<_n; j++)
			val[i][j] = M[i][j];

}
// ------------------------------------------------------------------------- //
matrix& matrix::operator=(const matrix& M)
{
	if(val != 0)
	{
		for(int i=0; i<_m; i++)
		{
			if(val[i] != 0)
				delete [] val[i];
		}
		delete [] val;

		val = 0;
	}

	_m = M.m();
	_n = M.n();

	val = new double*[_m];
	for(int i=0; i<_m; i++)
		val[i] = new double[_n];

	for(int i=0; i<_m; i++)
		for(int j=0; j<_n; j++)
			val[i][j] = M[i][j];

	return *this;
}
// ------------------------------------------------------------------------- //
matrix& matrix::operator=(double v)
{
	for(int i=0; i<_m; i++)
		for(int j=0; j<_n; j++)
			val[i][j] = v;

	return *this;
}
// ------------------------------------------------------------------------- //
matrix::~matrix()
{
	for(int i=0; i<_m; i++)
		delete [] val[i];

	delete [] val;
}
// ------------------------------------------------------------------------- //
void matrix::resize(int m, int n)
{
	// update dimension
	_m = m;
	_n = n;


	// deallocate storage
	if(val != 0)
	{
		for(int i=0; i<_m; i++)
			delete [] val[i];
		delete [] val;
		val = 0;
	}

	// reallocate the new storage
	val = new double*[_m];
	for(int i=0; i<_m; i++)
		val[i] = new double[_n];

	// set all the values to zero
	for(int i=0; i<_m; i++)
		for(int j=0; j<_n; j++)
			val[i][j] = 0;

}
// -------------------------  operator * -------------------------- //
matrix operator*(const matrix& left, const matrix& right)
{
	if(left.n() != right.m())
	{
		Log << "operator* -> can't multiply a " 
			 << left.m() << " x " << left.n()
			 << " matrix with a "
			 << right.m() << " x " << right.n()
			 << " matrix!" 
			 << gaplog::endl;
		
		return matrix();
	}

	else
	{
		matrix tmp(left.m(),right.n());

		for(int i=0; i<tmp.m(); i++)
		{
			for(int j=0; j<tmp.n(); j++)
			{	
				for(int h=0; h<left.n(); h++)
					tmp[i][j] += left[i][h]*right[h][j];
			}
		}		
		return (tmp);
	}

}
// ------------------------------------------------------------------------- //
matrix operator*(const matrix& left, double v)
{
	matrix tmp(left);
	for(int i=0; i<tmp.m(); i++)
		for(int j=0; j<tmp.n(); j++)
			tmp[i][j] *= v;
		
	return tmp;

}
// ------------------------------------------------------------------------- //
matrix operator*(double v, const matrix& right)
{
	matrix tmp(right);
	for(int i=0; i<tmp.m(); i++)
		for(int j=0; j<tmp.n(); j++)
			tmp[i][j] *= v;
		
	return tmp;
}
// -------------------------  operator / -------------------------- //
matrix operator/(const matrix& left, double v)
{
	if(v != 0)
	{
		matrix tmp(left);
		for(int i=0; i<tmp.m(); i++)
			for(int j=0; j<tmp.n(); j++)
				tmp[i][j] /= v;
			
		return tmp;
	}
	else
	{
		Log << "matrix::operator/ -> attempt to divide by zero!" << gaplog::endl;
		return left;
	}

}
// -------------------------  operator + -------------------------- //
matrix operator+(const matrix& left, const matrix& right)
{
	if(left.n() != right.n() && left.m() != right.m())
	{
		Log << "operator+ -> can't sum a " 
			 << left.m() << " x " << left.n()
			 << " matrix with a "
			 << right.m() << " x " << right.n()
			 << " matrix!" 
			 << gaplog::endl;
		
		return matrix();
	}

	else
	{
		matrix tmp(left.m(),right.n());

		for(int i=0; i<tmp.m(); i++)
			for(int j=0; j<tmp.n(); j++)
				tmp[i][j] = left[i][j] + right[i][j];
					
		return tmp;
	}

}
// ------------------------------------------------------------------------- //
matrix operator+(const matrix& left, double v)
{
	matrix tmp(left);
	for(int i=0; i<tmp.m(); i++)
		for(int j=0; j<tmp.n(); j++)
			tmp[i][j] += v;
		
	return tmp;
}
// ------------------------------------------------------------------------- //
matrix operator+(double v, const matrix& right)
{
	matrix tmp(right);
	for(int i=0; i<tmp.m(); i++)
		for(int j=0; j<tmp.n(); j++)
			tmp[i][j] += v;
		
	return tmp;
}
// -------------------------  operator - -------------------------- //
matrix operator-(const matrix& left, const matrix& right)
{
	if(left.n() != right.n() && left.m() != right.m())
	{
		Log << "operator- -> can't subtract a " 
			 << left.m() << " x " << left.n()
			 << " matrix with a "
			 << right.m() << " x " << right.n()
			 << " matrix!" 
			 << gaplog::endl;
		
		return matrix();
	}

	else
	{
		matrix tmp(left.m(),right.n());

		for(int i=0; i<tmp.m(); i++)
			for(int j=0; j<tmp.n(); j++)
				tmp[i][j] = left[i][j] - right[i][j];
					
		return tmp;
	}

}
// ------------------------------------------------------------------------- //
matrix operator-(const matrix& left, double v)
{
	matrix tmp(left);
	for(int i=0; i<tmp.m(); i++)
		for(int j=0; j<tmp.n(); j++)
			tmp[i][j] -= v;
		
	return tmp;
}
// ------------------------------------------------------------------------- //
matrix operator-(double v, const matrix& right)
{
	matrix tmp(right);
	for(int i=0; i<tmp.m(); i++)
		for(int j=0; j<tmp.n(); j++)
			tmp[i][j] -= v;
		
	return tmp;
}
// ------------------------------------------------------------------------- //
matrix cross (const matrix& left, const matrix& right)
{
	matrix crossprod;
	if(left.m() == 3 && left.n() == 1)	// check if it is a column vector
	{
		if(right.m() == 3 && right.n() == 1)
		{
			crossprod = matrix(3,1);
			crossprod[0][0] = left[1][0]*right[2][0] - left[2][0]*right[1][0];
			crossprod[1][0] = left[2][0]*right[0][0] - left[0][0]*right[2][0];
			crossprod[2][0] = left[0][0]*right[1][0] - left[1][0]*right[0][0];
		}
		else
			Log << "cross product allowed only beetwen 1x3 or 3x1 matrices!" << gaplog::endl;
	}
	else if(left.m() == 1 && left.n() == 3)	// check if it is a row vector
	{
		if(right.m() == 1 && right.n() == 3)
		{
			crossprod = matrix(1,3);
			crossprod[0][0] = left[0][1]*right[0][2] - left[0][2]*right[0][1];
			crossprod[0][1] = left[0][2]*right[0][0] - left[0][0]*right[0][2];
			crossprod[0][2] = left[0][0]*right[0][1] - left[0][1]*right[0][0];
		}
		else
			Log << "cross product allowed only beetwen 1x3 or 3x1 matrices!" << gaplog::endl;
	}
	else
		Log << "cross product allowed only beetwen 1x3 or 3x1 matrices!" << gaplog::endl;


	return crossprod;
}
// ------------------------------------------------------------------------- //
// get the minor for position (row,col)
matrix matrix::getMinor(int row, int col)
{
	matrix tmp(_n-1,_m-1);

	for(int i=0; i<row; i++)
	{
		for(int j=0; j<col; j++)
			tmp[i][j] = val[i][j];
		for(int j=col+1; j<_n; j++)
		
			tmp[i][j-1] = val[i][j];
	}
	for(int i=row+1; i<_m; i++)
	{
		for(int j=0; j<col; j++)
			tmp[i-1][j] = val[i][j];
		for(int j=col+1; j<_n; j++)
		
			tmp[i-1][j-1] = val[i][j];
	}


	return tmp;
   
}
// ------------------------------------------------------------------------- //
// calculate the determinant (algebric method)
double matrix::det()
{
	double d = 0;

	if(_m == _n)
	{

		if(_m == 1)
			return val[0][0];
		else
		{
			for(int i=0; i<_m; i++)
			{
				matrix M(getMinor(i,0));
				d += pow(-1.0,i)*val[i][0]*M.det();		
			}
		}
	}
	else
		Log << "matrix::det() -> can't compute the determinant of a non-square matrix!" << gaplog::endl;

	return d;
}
// ------------------------------------------------------------------------- //
// get the transpose of the matrix
matrix matrix::T()
{
	matrix tmp(_n,_m);
	for(int i=0; i<_m; i++)
		for(int j=0; j<_n; j++)
			tmp[j][i] = val[i][j];

	return tmp;
}
// ------------------------------------------------------------------------- //
matrix matrix::inv()
{
	if(_m != _n)
	{
		Log << "matrix::inv -> can't invert a non square matrix!" << gaplog::endl;
		return *this;
	}
	else
	{

		double d = this -> det();

		if(d != 0)
		{

			matrix inv(_n,_m);
			matrix tmp = this -> T();
			
			for(int i=0; i<_m; i++)
			{
				for(int j=0; j<_n; j++)
				{
					matrix minor = tmp.getMinor(i,j);
					inv[i][j] = (1/d)*pow(-1.0,i+j)*minor.det();
				}
			}

			return inv;
		}
		else
		{
			Log << "matrix::inv -> the matrix is singular to the machine precision!" << gaplog::endl;
			return *this;
		}
	}
}
// ------------------------------------------------------------------------- //
std::ostream& operator<<(std::ostream& os, const matrix& t) 
{

	os  << t.m() << " x " << t.n() << gaplog::endl;
	os	<< "[" <<gaplog::endl;

	for(int i=0; i<t.m(); i++)
	{
		os << "\t";
		for(int j=0; j<t.n(); j++)
			os << t[i][j] << "\t";
		os <<gaplog::endl;
	}

	os	<< "]" <<gaplog::endl;

	return os;
}
// ------------------------------------------------------------------------- //



