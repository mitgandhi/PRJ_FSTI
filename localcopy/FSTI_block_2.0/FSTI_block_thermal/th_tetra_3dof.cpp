# include "./th_tetra_3dof.h"

// ------------------------------------------------------------------------- //
th_tetra_3dof::th_tetra_3dof(int _id, const tetra_sub_mesh& _smsh) 
	: id(_id), smsh(_smsh), gm(smsh.get_elm(_id))
{

	m = 0;

	// define the global node numbering
	gnum.resize(0);
	for(int n=0; n<4; n++)
	{
		// local node id node id
		int nid = smsh.g2l_nn.at(gm.nodes[n]);	
		for(int j=0; j<ndof; j++)
			gnum.push_back(ndof*nid + j);
	}

}
// ------------------------------------------------------------------------- //
matrix th_tetra_3dof::calc_D()
{
	matrix D(6,6);

	double coef  = m->E/((1.0 + m->nu)*(1.0 - 2.0*m->nu));
	D[0][0] = D[1][1] = D[2][2] = coef*(1.0 - m->nu);
	D[3][3] = D[4][4] = D[5][5] = coef*(1.0 - 2.0*m->nu)/2.0;
	D[0][1] = D[1][0] = D[0][2] = D[2][0] = D[1][2] = D[2][1] = coef*m->nu;

	return D;
}
// ------------------------------------------------------------------------- //
matrix th_tetra_3dof::calc_B()
{
	
	matrix B(6,12);

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

	// define the B matrix
	
	B[0][0] = invJdet*a0;
	B[0][3] = invJdet*a1;
	B[0][6] = invJdet*a2;
	B[0][9] = invJdet*a3;
	//
	B[1][1] = invJdet*b0;
	B[1][4] = invJdet*b1;
	B[1][7] = invJdet*b2;
	B[1][10] = invJdet*b3;
	//
	B[2][2] = invJdet*c0;
	B[2][5] = invJdet*c1;
	B[2][8] = invJdet*c2;
	B[2][11] = invJdet*c3;
	//
	B[3][0] = invJdet*b0;
	B[3][1] = invJdet*a0;
	B[3][3] = invJdet*b1;
	B[3][4] = invJdet*a1;
	B[3][6] = invJdet*b2;
	B[3][7] = invJdet*a2;
	B[3][9] = invJdet*b3;
	B[3][10] = invJdet*a3;
	//
	B[4][1] = invJdet*c0;
	B[4][2] = invJdet*b0;
	B[4][4] = invJdet*c1;
	B[4][5] = invJdet*b1;
	B[4][7] = invJdet*c2;
	B[4][8] = invJdet*b2;
	B[4][10] = invJdet*c3;
	B[4][11] = invJdet*b3;
	//
	B[5][0] = invJdet*c0;
	B[5][2] = invJdet*a0;
	B[5][3] = invJdet*c1;
	B[5][5] = invJdet*a1;
	B[5][6] = invJdet*c2;
	B[5][8] = invJdet*a2;
	B[5][9] = invJdet*c3;
	B[5][11] = invJdet*a3;

	return B;
}
// ------------------------------------------------------------------------- //
matrix th_tetra_3dof::calc_K()
{
	
	matrix D = calc_D();

	matrix B = calc_B();
	
	// get the stiffness
	matrix K = smsh.get_elm(id).volume*((B.T()*D)*B);

	return K;

}
// ------------------------------------------------------------------------- //