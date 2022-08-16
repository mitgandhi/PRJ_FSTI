//This implements a "dense" style matrix and limited matrix operations using std::vector
namespace CasparSlipperFEM
{
	template <class T>
	class dmatrix
	{
		vector<vector<T> > mat;

	public:
		int m;
		int n;

		//Constructor
		dmatrix(const int i, const int j)
		{
			m = i;
			n = j;
			
			mat.resize(m);
			for(int r=0; r<m; r++)
			{
				mat[r].resize(n);
				for(int c=0; c<n; c++)
				{
					mat[r][c] = 0;
				}
			}
		};

		dmatrix()
		{
		};

		void resize(const int i, const int j)
		{
			//clear the current matrix
			mat.clear();

			m = i;
			n = j;
			
			mat.resize(m);
			for(int r=0; r<m; r++)
			{
				mat[r].resize(n);
				for(int c=0; c<n; c++)
				{
					mat[r][c] = 0;
				}
			}
		};

		//Direct matrix access (i,j)
		T & operator() (const int i, const int j) 
		{
			return mat[i][j];
		}

		const T & operator() (const int i, const int j) const
		{
			return mat[i][j];
		}

		//Direct matrix access [i]
		T & operator[] (const int i) 
		{
			return mat[i/m][i%m];
		}

		//Scalar multiplication
		dmatrix<T> & operator*= (const double s)
		{
			for(int i=0; i<m; i++)
			{
				for(int j=0; j<n; j++)
				{
					mat[i][j] *= s;
				}
			}
			return *this;
		}

		//Matrix multiplication
		dmatrix<T> operator* (const dmatrix<T> & M) const
		{
			if(n != M.m)
			{
				//Invalid matrix sizes for multiplication
				return dmatrix<T> (0,0);
			}
			dmatrix<T> r(m,M.n);

			for(int i=0; i<m; i++)
			{
				for(int j=0; j<M.n; j++)
				{
					for(int k=0; k<n; k++)
					{
						r(i,j) += mat[i][k]*M(k,j);
					}
				}
			} 

			return r;
		}

		//Matrix addition
		dmatrix<T> & operator+= (const dmatrix<T> M)
		{
			if(m != M.m || n != M.n)
			{
				//Invalid matrix sizes for addition so don't do anything
				return *this;
			}

			for(int i=0; i<m; i++)
			{
				for(int j=0; j<n; j++)
				{
					mat[i][j] += M(i,j);
				}
			}
			return *this;
		}

		//Transpose
		dmatrix<T> t()
		{
			dmatrix<T> trans(n,m);

			for(int i=0; i<m; i++)
			{
				for(int j=0; j<n; j++)
				{
					trans(j,i) = mat[i][j];
				}
			}

			return trans;
		}

		//Matrix Determinate
		double det()
		{
			if(m != n)
			{
				return 0;
			}

			if(m == 1)
			{
				return mat[0][0];
			}

			//We will use Laplace's formula for the determinate
			double d = 0;
			for(int j=0; j<n; j++)
			{
				//Create the 'Minor' matrix
				dmatrix M(m-1, n-1);
				for(int r=1; r<m; r++)
				{
					for(int c=0; c<n; c++)
					{
						if(c < j)
						{
							M(r-1,c) = mat[r][c];
						} else if (c > j)
						{
							M(r-1,c-1) = mat[r][c];
						}
					}
				}
				
				d += pow(-1.0,j) * mat[0][j] * M.det();
			}

			return d;
		}

	};
};
