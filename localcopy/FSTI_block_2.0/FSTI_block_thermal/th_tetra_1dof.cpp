# include "./th_tetra_1dof.h"
# include "../FSTI_Block_dll/log.h"

extern class gaplog Log;
using namespace std;

// ------------------------------------------------------------------------- //
th_tetra_1dof::th_tetra_1dof(int _id, const tetra_mesh& _msh) : id(_id), msh(_msh)
{

	m = 0;

	// define the global node numbering
	gnum.resize(0);
	for(int n=0; n<4; n++)
	{
		int nid = msh.elements[id].nodes[n];	// node id
		for(int j=0; j<ndof; j++)
			gnum.push_back(ndof*nid + j);
	}

}
// ------------------------------------------------------------------------- //
matrix th_tetra_1dof::calc_K()
{

	if(m == 0)
	{
		Log << "\nth_tetra::calc_K() material was not initialized!" << gaplog::endl;
		exit(1);
	}

	matrix D(3,3);
	D[0][0] = m->lambda;
	D[1][1] = m->lambda;
	D[2][2] = m->lambda;

	matrix B(3,4);

	// matrix entries
	double a0 = (y(1) - y(3))*(z(2) - z(3)) - (y(2) - y(3))*(z(1) - z(3));
	double a1 = (y(2) - y(3))*(z(0) - z(3)) - (y(0) - y(3))*(z(2) - z(3));
	double a2 = (y(0) - y(3))*(z(1) - z(3)) - (y(1) - y(3))*(z(0) - z(3));
	double a3 = -(a0 + a1 + a2);
	double b0 = (x(2) - x(3))*(z(1) - z(3)) - (x(1) - x(3))*(z(2) - z(3));
	double b1 = (x(0) - x(3))*(z(2) - z(3)) - (x(2) - x(3))*(z(0) - z(3));
	double b2 = (x(1) - x(3))*(z(0) - z(3)) - (x(0) - x(3))*(z(1) - z(3));
	double b3 = -(b0 + b1 + b2);
	double c0 = (x(1) - x(3))*(y(2) - y(3)) - (x(2) - x(3))*(y(1) - y(3));
	double c1 = (x(2) - x(3))*(y(0) - y(3)) - (x(0) - x(3))*(y(2) - y(3));
	double c2 = (x(0) - x(3))*(y(1) - y(3)) - (x(1) - x(3))*(y(0) - y(3));
	double c3 = -(c0 + c1 + c2);

	double invJdet = 1.0/( a0*(x(0) - x(3)) + b0*(y(0) - y(3)) + c0*(z(0) - z(3)) );

	// (d[N]'/dx)(d[N]/dx)
	B[0][0] = invJdet*a0;
	B[0][1] = invJdet*a1;
	B[0][2] = invJdet*a2;
	B[0][3] = invJdet*a3;
	// (d[N]'/dy)(d[N]/dy)
	B[1][0] = invJdet*b0;
	B[1][1] = invJdet*b1;
	B[1][2] = invJdet*b2;
	B[1][3] = invJdet*b3;
	// (d[N]'/dz)(d[N]/dz)
	B[2][0] = invJdet*c0;
	B[2][1] = invJdet*c1;
	B[2][2] = invJdet*c2;
	B[2][3] = invJdet*c3;

	//
	matrix K = msh.elements[id].volume*((B.T()*D)*B);

	return K;
}
// ------------------------------------------------------------------------- //
matrix th_tetra_1dof::calc_intNTN(int side) const
{
	matrix NTN(4,4);

	// fill the matrix
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<4; j++)
		{
			if(i == side || j == side)
				NTN[i][j] = 0;
			else if(i==j)
				NTN[i][j] = 1.0;
			else if(i!=j)
				NTN[i][j] = 0.5;
		}
	}

	// define the face vertices
	vector<int> vtx(0);
	for(int i=0; i<4; i++)
	{
		if(i != side)
			vtx.push_back(i);
	}
	
	point v0 = msh.nodes[msh.elements[id].nodes[vtx[0]]];
	point v1 = msh.nodes[msh.elements[id].nodes[vtx[1]]];
	point v2 = msh.nodes[msh.elements[id].nodes[vtx[2]]];

	// calculate the area
	double A = 0.5*((v1 - v0) ^ (v2 - v0)).mag();
	
	// multiply by A/6 constnat (deriving from integration)
	NTN = NTN*(A/6.0);

	return NTN;
}
// ------------------------------------------------------------------------- //
matrix th_tetra_1dof::calc_intN(int side) const
{

	// resize to 4 row 1 column matrix
	matrix N(4,1);

	// fill the matrix
	for(int i=0; i<4; i++)
	{
		if(i == side)
			N[i][0] = 0;
		else
			N[i][0] = 1.0;
	}

	// define the face vertices
	vector<int> vtx(0);
	for(int i=0; i<4; i++)
	{
		if(i != side)
			vtx.push_back(i);
	}

	point v0 = msh.nodes[msh.elements[id].nodes[vtx[0]]];
	point v1 = msh.nodes[msh.elements[id].nodes[vtx[1]]];
	point v2 = msh.nodes[msh.elements[id].nodes[vtx[2]]];

	// calculate the area
	double A = 0.5*((v1 - v0) ^ (v2 - v0)).mag();
	
	// multiply by A/3 constnat (deriving from integration)
	N = N*(A/3.0);

	return N;
}
// ------------------------------------------------------------------------- //