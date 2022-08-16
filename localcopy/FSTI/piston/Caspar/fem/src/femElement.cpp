# include "femElement.h"
# include <string>
# include <iostream>
# include <fstream>
using namespace std;

femElement::femElement() :
	E(0),
	nu(0),
	x(0),
	y(0),
	z(0),
	conn(0),
	globalNumbering(0),
	stiffMat(0)
{

}

femElement::~femElement()
{
	// the virtual destructor will call the proper destructor for each object
}

void femElement::calcStiff()
{

}



// ------------------------------------------------------------------------------------------------------ //
// --------------------------------- isoparametric 8 nodes solid ---------------------------------------- //
// ------------------------------------------------------------------------------------------------------ //
brick::brick(int intOp, double E, double nu, double* x, double* y, double* z, int* conn) {
	
	// initialize the scalar values
	this -> intOp = intOp;
	this -> E = E;
	this -> nu = nu;
	jacobDet = 0.0;
	jacobDetW = 0.0;
	v = 0.0;
	r1 = r2 = r3 = 0.0;

	// initialize the vectors
	this -> conn = new int[8];
	memcpy (this -> conn, conn, 8*sizeof(int));		// initialization of constitution matrix
	this -> x = new double[8];
	memcpy (this -> x, x, 8*sizeof(double));		// initialization of x coordinate vector
	this -> y = new double[8];
	memcpy (this -> y, y, 8*sizeof(double));		// initialization of y coordinate vector
	this -> z = new double[8];
	memcpy (this -> z, z, 8*sizeof(double));		// initialization of z coordinate vector
	
	globalNumbering = new int[24];
	memset (globalNumbering,0,24*sizeof(int));		// initialization of the global numbering vector for the element

	for(int i = 0, j = 0; i < 24; i += 3,j++) {
		globalNumbering[i]= conn[j]*3;			
		globalNumbering[i+1] = conn[j]*3 + 1;	
		globalNumbering[i+2] = conn[j]*3 + 2;	
	}

	// allocate storage for the stiffness matrix
	stiffMat = new double *[24];
	for(int i = 0; i < 24; i++) {
		stiffMat[i] = new double[24];
	}
			
}

brick::~brick()
{
	for(int i = 0; i < 24; i++)
		delete [] stiffMat[i];
	delete [] stiffMat;

	delete [] x;
	delete [] y;
	delete [] z;
	delete [] conn;
	delete [] globalNumbering;
}

void brick::calcJacob() {

	memset (jacob,0,9*sizeof(double));		// initialization of jacobian matrix
	memset (jacobInv,0,9*sizeof(double));	// initialization of inverse jacobian matrix
	
	// filling the Jacobian matrix
	// the Jacobian matrix is defined by eq (12.17) at pag 389; each derivative in the matrix, can be
	// expressed in term of shape functions using equations (12.15) at pag 388 of G.Hartley book.
	for(int k = 0; k < 8; k++) {
		for(int j = 0; j < 3; j++) {
			jacob[0][j] += dndr[k][j]*x[k];
			jacob[1][j] += dndr[k][j]*y[k];
			jacob[2][j] += dndr[k][j]*z[k];
		}
	}
	
	// calculating the Jacobian determinant
	jacobDet = jacob[0][0]*jacob[1][1]*jacob[2][2] 
		+ jacob[0][1]*jacob[1][2]*jacob[2][0]
		+ jacob[0][2]*jacob[1][0]*jacob[2][1]
		- jacob[2][0]*jacob[1][1]*jacob[0][2]
		- jacob[2][1]*jacob[1][2]*jacob[0][0]
		- jacob[2][2]*jacob[1][0]*jacob[0][1];

	// calculation of jacobian inverse matrix
	jacobInv[0][0] = (jacob[1][1]*jacob[2][2] - jacob[1][2]*jacob[2][1])/jacobDet;
	jacobInv[0][1] = (jacob[1][2]*jacob[2][0] - jacob[1][0]*jacob[2][2])/jacobDet;
	jacobInv[0][2] = (jacob[1][0]*jacob[2][1] - jacob[1][1]*jacob[2][0])/jacobDet;
	jacobInv[1][0] = (jacob[0][2]*jacob[2][1] - jacob[0][1]*jacob[2][2])/jacobDet;
	jacobInv[1][1] = (jacob[0][0]*jacob[2][2] - jacob[0][2]*jacob[2][0])/jacobDet;
	jacobInv[1][2] = (jacob[0][1]*jacob[2][0] - jacob[0][0]*jacob[2][1])/jacobDet;
	jacobInv[2][0] = (jacob[0][1]*jacob[1][2] - jacob[0][2]*jacob[1][1])/jacobDet;
	jacobInv[2][1] = (jacob[0][2]*jacob[1][0] - jacob[0][0]*jacob[1][2])/jacobDet;
	jacobInv[2][2] = (jacob[0][0]*jacob[1][1] - jacob[0][1]*jacob[1][0])/jacobDet;


	// jacobInv[2][1] = (jacob[0][2]*jacob[0][1] - jacob[0][0]*jacob[1][2])/jacobDet; //wrong
}

void brick::calcD() 
{
	// filling D_mat
	memset (D,0,36*sizeof(double));				// initialization of constitution matrix
	double coef  = E/((1.0 + nu)*(1.0 - 2.0*nu));
	D[0][0] = D[1][1] = D[2][2] = coef*(1.0 - nu);
	D[3][3] = D[4][4] = D[5][5] = coef*(1.0 - 2.0*nu)/2.0;
	D[0][1] = D[1][0] = D[0][2] = D[2][0] = D[1][2] = D[2][1] = coef*nu;

}

void brick::calcB() {

	// initialization of the local derivative of interpolation functions matrix
	memset (dndr,0,24*sizeof(double));	
	
	// see pag 391 of G. Hartley book
	dndr[0][0] = -0.1250 * (1.0 - r2) * (1.0 - r3);
	dndr[0][1] = -0.1250 * (1.0 - r1) * (1.0 - r3);
	dndr[0][2] = -0.1250 * (1.0 - r1) * (1.0 - r2);
	
	dndr[1][0] = 0.1250 * (1.0 - r2) * (1.0 - r3);
	dndr[1][1] = -0.1250 * (1.0 + r1) * (1.0 - r3);
	dndr[1][2] = -0.1250 * (1.0 + r1) * (1.0 - r2);
	
	dndr[2][0] = 0.1250 * (1.0 + r2) * (1.0 - r3);
	dndr[2][1] = 0.1250 * (1.0 + r1) * (1.0 - r3);
	dndr[2][2] = -0.1250 * (1.0 + r1) * (1.0 + r2);
	
	dndr[3][0] = -0.1250 * (1.0 + r2) * (1.0 - r3);
	dndr[3][1] = 0.1250 * (1.0 - r1) * (1.0 - r3);
	dndr[3][2] = -0.1250 * (1.0 - r1) * (1.0 + r2);
	
	dndr[4][0] = -0.1250 * (1.0 - r2) * (1.0 + r3);
	dndr[4][1] = -0.1250 * (1.0 - r1) * (1.0 + r3);
	dndr[4][2] = 0.1250 * (1.0 - r1) * (1.0 - r2);
	
	dndr[5][0] = 0.1250 * (1.0 - r2) * (1.0 + r3);
	dndr[5][1] = -0.1250 * (1.0 + r1) * (1.0 + r3);
	dndr[5][2] = 0.1250 * (1.0 + r1) * (1.0 - r2);
	
	dndr[6][0] = 0.1250 * (1.0 + r2) * (1.0 + r3);
	dndr[6][1] = 0.1250 * (1.0 + r1) * (1.0 + r3);
	dndr[6][2] = 0.1250 * (1.0 + r1) * (1.0 + r2);
	
	dndr[7][0] = -0.1250 * (1.0 + r2) * (1.0 + r3);
	dndr[7][1] = 0.1250 * (1.0 - r1) * (1.0 + r3);
	dndr[7][2] = 0.1250 * (1.0 - r1) * (1.0 + r2);


	// jacobian calculation
	calcJacob();

	// calculating ndnx as dndx = (J^-1)(dndr)
	memset (dndx,0,24*sizeof(double));	// initialization of the global derivative of interpolation functions matrix
	for(int i = 0; i < 8; i++) {
		for(int j = 0; j < 3; j++) {
			dndx[i][j] = jacobInv[j][0]*dndr[i][0] + 
						+ jacobInv[j][1]*dndr[i][1]
						+ jacobInv[j][2]*dndr[i][2];
		}
	}

		
	// Construction of B matrix (B matrix rapresents the displacements in terms of global coordinate system)
	
	memset (B,0,144*sizeof(double)); // initialization of displacement-strain matrix

	for(int i = 0; i < 8; i++) {
		int i1 = 3*i, i2 = 3*i + 1, i3 = 3*i + 2; // index for the proper position in Bmatrix (see pag 392)
		B[0][i1] = dndx[i][0];
		B[1][i2] = dndx[i][1];
		B[2][i3] = dndx[i][2];
		B[3][i1] = B[1][i2];
		B[3][i2] = B[0][i1];
		B[4][i2] = B[2][i3];
		B[4][i3] = B[1][i2];
		B[5][i1] = B[2][i3];
		B[5][i3] = B[0][i1];
	}


}

void brick::calcBDB() {

	// the product is performed calculating [D]*[B], and then [B]T*([D]*[B])
	
	// definition and intialization of the DB matrix
	double DB[6][24];
	memset(DB,0,144*sizeof(double));

	// calculation of the product 
	for(int i = 0; i < 6; i++) {
		for(int k = 0; k < 24; k++) {
			for(int j = 0; j < 6; j++) {
				DB[i][k] += D[i][j]*B[j][k];
			}
		}
	}

	// calculation of the [B]T*([D]*[B]) product
	for(int i = 0; i < 24; i++) {
		for(int k = 0; k < 24; k++) {
			double tmp = 0.0;
			for(int j = 0; j < 6; j++) {
				tmp += B[j][i]*DB[j][k];
			}
			stiffMat[i][k] += tmp*jacobDetW;
		}
	}

}

void brick::calcStiff() {
	
	// Sampling point of Gaussian numerical integral
	double arg_Gauss_point[] = {		
		0.0, 0.0, 0.0, 0.0,																// one point
		-0.577350269189626, 0.577350269189626, 0.0, 0.0,								// two points
		-0.774596669241483, 0.0, 0.774596669241483, 0.0,								// three points
		-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053	// four points
	};
	
	// Weight of Gaussian numerical integral
	double arg_Gauss_weight[] = {
		2.0, 0.0, 0.0, 0.0,																// one point
		1.0, 1.0, 0.0, 0.0,																// two points
		0.555555555555556, 0.888888888888889, 0.555555555555556, 0.0,					// three points
		0.347854845137454,0.652145154862546, 0.652145154862546, 0.347854845137454		// four points
	};
	
	// fill the constitution matrix
	calcD();

		
	// Gaussian numerical integration for 3 dimensions
	v = 0.0;
	jacobDetW = 0.0;

	// initialize the stiffness matrix for the element
	for(int i=0; i<24; i++)
		memset(stiffMat[i],0,24*sizeof(double));
	
	for (int i = 0; i < intOp; i++) {
		for (int j = 0; j < intOp; j++) {
			for (int k = 0; k < intOp; k++) {
				// definition of integration points 
				r1 = arg_Gauss_point[i + 4*(intOp - 1)];
				r2 = arg_Gauss_point[j + 4*(intOp - 1)];
				r3 = arg_Gauss_point[k + 4*(intOp - 1)];

				// calculation of displacement-strain matrix
				calcB();
				// calculation of weightning factors products
				jacobDetW = jacobDet*arg_Gauss_weight[i + 4*(intOp - 1)]
						*arg_Gauss_weight[j + 4*(intOp - 1)]
						*arg_Gauss_weight[k + 4*(intOp - 1)];
				// calculation of the element volume
				v += jacobDetW;
				// calculation of the product [B]T[D][B]
				calcBDB();
			}
		}
	}
	
	
	/*
	// code to print out the element stiffness matrix

	std::ofstream out("./fem/solver/elementStiffness.txt") ;
	for(int i=0; i<24; i++) {
		for(int j=0; j<24; j++) {
			out << std::scientific << stiffMat[i][j] << "\t";
		}
		out << std::"\n";
	}
	out.close();
	*/
	
}


// -------------------------------------------------------------------------------------------------- //
// --------------------------------- four-nodes tetrahedron  ---------------------------------------- //
// -------------------------------------------------------------------------------------------------- //
tetra::tetra(double E, double nu, double *x, double *y, double *z, int *conn) {
	
	// initialize the scalar values
	this -> E = E;
	this -> nu = nu;
	
	// initialize the vectors
	this -> conn = new int[4];
	memcpy (this -> conn, conn, 4*sizeof(int));		// initialization of constitution matrix
	this -> x = new double[4];
	memcpy (this -> x, x,4*sizeof(double));			// initialization of x coordinate vector
	this -> y = new double[4];
	memcpy (this -> y, y,4*sizeof(double));			// initialization of y coordinate vector
	this -> z = new double[4];
	memcpy (this -> z, z,4*sizeof(double));			// initialization of z coordinate vector

	double a0 = (y[1] - y[3])*(z[2] - z[3]) - (y[2] - y[3])*(z[1] - z[3]);
	double b0 = (x[2] - x[3])*(z[1] - z[3]) - (x[1] - x[3])*(z[2] - z[3]);
	double c0 = (x[1] - x[3])*(y[2] - y[3]) - (x[2] - x[3])*(y[1] - y[3]);
	
	jacobDet =	a0*(x[0] - x[3]) + b0*(y[0] - y[3]) + c0*(z[0] - z[3]);

	// the element is inside-out, invert the jacobian sign
	v = jacobDet/6.0;
	if(v < 0)
	{
		v *= -1.0;
	}

	globalNumbering = new int[12];
	memset (globalNumbering,0,12*sizeof(int));	// initialization of the global numbering vector for the element
	for(int i=0,j=0; i<12; i+=3,j++) {
		globalNumbering[i] = conn[j]*3;			
		globalNumbering[i+1] = conn[j]*3 + 1;	
		globalNumbering[i+2] = conn[j]*3 + 2;	
	}

	stiffMat = new double *[12];
	for(int i=0; i<12; i++) {
		stiffMat[i] = new double[12];
	}
}

tetra::~tetra()
{
	for(int i=0; i<12; i++)
		delete [] stiffMat[i];
	delete [] stiffMat;

	delete [] x;
	delete [] y;
	delete [] z;
	delete [] conn;
	delete [] globalNumbering;
}

void tetra::calcB() {

	memset(B,0,72*sizeof(double));

	// matrix entries
	double a0 = (y[1] - y[3])*(z[2] - z[3]) - (y[2] - y[3])*(z[1] - z[3]);
	double a1 = (y[2] - y[3])*(z[0] - z[3]) - (y[0] - y[3])*(z[2] - z[3]);
	double a2 = (y[0] - y[3])*(z[1] - z[3]) - (y[1] - y[3])*(z[0] - z[3]);
	double a3 = -(a0 + a1 + a2);
	double b0 = (x[2] - x[3])*(z[1] - z[3]) - (x[1] - x[3])*(z[2] - z[3]);
	double b1 = (x[0] - x[3])*(z[2] - z[3]) - (x[2] - x[3])*(z[0] - z[3]);
	double b2 = (x[1] - x[3])*(z[0] - z[3]) - (x[0] - x[3])*(z[1] - z[3]);
	double b3 = -(b0 + b1 + b2);
	double c0 = (x[1] - x[3])*(y[2] - y[3]) - (x[2] - x[3])*(y[1] - y[3]);
	double c1 = (x[2] - x[3])*(y[0] - y[3]) - (x[0] - x[3])*(y[2] - y[3]);
	double c2 = (x[0] - x[3])*(y[1] - y[3]) - (x[1] - x[3])*(y[0] - y[3]);
	double c3 = -(c0 + c1 + c2);

	double invDet = (1/jacobDet);

	B[0][0] = invDet*a0;
	B[0][3] = invDet*a1;
	B[0][6] = invDet*a2;
	B[0][9] = invDet*a3;
	//
	B[1][1] = invDet*b0;
	B[1][4] = invDet*b1;
	B[1][7] = invDet*b2;
	B[1][10] = invDet*b3;
	//
	B[2][2] = invDet*c0;
	B[2][5] = invDet*c1;
	B[2][8] = invDet*c2;
	B[2][11] = invDet*c3;
	//
	B[3][0] = invDet*b0;
	B[3][1] = invDet*a0;
	B[3][3] = invDet*b1;
	B[3][4] = invDet*a1;
	B[3][6] = invDet*b2;
	B[3][7] = invDet*a2;
	B[3][9] = invDet*b3;
	B[3][10] = invDet*a3;
	//
	B[4][1] = invDet*c0;
	B[4][2] = invDet*b0;
	B[4][4] = invDet*c1;
	B[4][5] = invDet*b1;
	B[4][7] = invDet*c2;
	B[4][8] = invDet*b2;
	B[4][10] = invDet*c3;
	B[4][11] = invDet*b3;
	//
	B[5][0] = invDet*c0;
	B[5][2] = invDet*a0;
	B[5][3] = invDet*c1;
	B[5][5] = invDet*a1;
	B[5][6] = invDet*c2;
	B[5][8] = invDet*a2;
	B[5][9] = invDet*c3;
	B[5][11] = invDet*a3;

}

void tetra::calcD() {

	memset(D,0,36*sizeof(double));

	double coef  = E/((1.0 + nu)*(1.0 - 2.0*nu));
	D[0][0] = D[1][1] = D[2][2] = coef*(1.0 - nu);
	D[3][3] = D[4][4] = D[5][5] = coef*(1.0 - 2.0*nu)/2.0;
	D[0][1] = D[1][0] = D[0][2] = D[2][0] = D[1][2] = D[2][1] = coef*nu;

}
void tetra::calcStiff() {
	
	// make sure is initialized to zero
	for(int i = 0; i < 12; i++) {
		memset(stiffMat[i],0,12*sizeof(double));
	}

	// filling the constitutive matrix
	calcD();
	// filling the B matrix 
	calcB();
	
	// product [B]^T*[D]*[B] is performed in two steps: [D]*[B] and then B^T*(D*B)
	// calculation of the [D][B] product
	double DB[6][12];
	memset(DB,0,72*sizeof(double));
	for(int i = 0; i < 6; i++) {
		for(int k = 0; k < 12; k++) {
			for(int j = 0; j < 6; j++) {
				DB[i][k] += D[i][j]*B[j][k];
			}
		}
	}

	// calculation of the [B]T*([D]*[B]) product
	for(int i = 0; i < 12; i++) {
		for(int k = 0; k < 12; k++) {
			double tmp = 0.0;
			for(int j = 0; j < 6; j++) {
				tmp += B[j][i]*DB[j][k];
			}
			stiffMat[i][k] += tmp;
		}
	}

	for(int i = 0; i < 12; i++) {
		for(int j = 0; j < 12; j++) {
			stiffMat[i][j] *= v;
		}
	}


	/*
	// code to print out the element stiffness matrix
	std::ofstream out("./fem/solver/elementStiffness.txt") ;
	for(int i=0; i<12; i++) {
		for(int j=0; j<12; j++) {
			out << std::scientific << stiffMat[i][j] << "\t";
		}
		out << std::"\n";
	}
	out.close();
	*/
}