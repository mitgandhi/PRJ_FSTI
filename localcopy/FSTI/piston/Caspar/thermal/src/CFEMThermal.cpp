#include "CFEMThermal.h"
#include "../../main/src/logger.h"
#include <iostream>
#include <iomanip>
#include "../../caspar_input/input.h"
#pragma once


extern struct sGapResult myGapResult;
extern class CGapInput myGapInput;
extern class CMesh myMesh;
extern class CThermal myThermal;
extern class input myinput;


CFEMThermal::CFEMThermal()
{
	
};
CFEMThermal::~CFEMThermal()
{
	
};
void CFEMThermal::FEMThermalSolve(string body,Array<double,1> T_body,Array<double,1> &def_surf,string mesh_name)
{
	Log << "\n";
	Log << "Solving " + body + " thermal deformation..." << "\n";
	Log << "\n";

	//Determine constraints for inertia relief
	if( (myMesh.IR_p==1 && myMesh.body==1) || (myMesh.IR_c==1 && myMesh.body==0) )
	{
		FEMThermalInertiaConstraints( );
	};

	//Write general file
	FEMThermalGeneral();

	//Write thermal constraints file
	FEMThermalConstraints();
	
	//Write thermal loads file
	FEMThermalLoads(T_body);

	double tol;
	int maxit;
	if(body == "cylinder"){
		tol = myinput.data.thermal.block.solver.tol;
		maxit = myinput.data.thermal.block.solver.maxit;
	}
	else if(body == "piston"){
		tol = myinput.data.thermal.piston.solver.tol;
		maxit = myinput.data.thermal.piston.solver.maxit;
	}
	else{
		Log << "Body " << body << " not supported for thermal deformation!" << endl;
		exit(1);
	}

	//Call Marco's fem
	fem::stiffness(".","ascii");
	fem::loads(".","ascii");
	fem::solve(".",maxit,tol);
	fem::writeDisplacement(".","ascii");


	//read displacement
	ifstream fin("./temp/piston/output/displacement.fem");
	if (!fin.is_open()) {
		Log << "Error opening results file displacement.fem!\n" << endl;
		exit(1);
	}
	udisp.resize(nNodes,3);
	udisp = 0.0;
	for(int i = 0; i<nNodes; i++)
	{
		fin >> udisp(i,0) >> udisp(i,1) >> udisp(i,2);
	}
	fin.close();


	//assign surface deformation under-relaxing final value
	Array<double,1> def_surf_old(def_surf.extent(0));
	def_surf_old = def_surf;
	//assign
	def_surf = FEMThermalSurfaceDeformation();
	//under-relax
	double AlphaTh = myinput.data.options_piston.numeric.AlphaTh;
	def_surf = def_surf_old + AlphaTh * (def_surf - def_surf_old);


	//write vtk
	FEMThermalOutput(1.0,body);


	//clear
	udisp.free(); myMesh.cxyz_th.free();

	Log << "\n";
	Log << "Done " + body + " thermal deformation!" << "\n";
	Log << "\n";

};
void CFEMThermal::FEMThermalGeneral(void)
{
	//Element type
	Teth=myMesh.Teth;
	Hexa=myMesh.Hexa;
	//Element type
	int elmType=0;
	if(Teth)
	{
		elmType=1;
	};
	//Number of gauss quadrature points considered in the element
	int GaussPts = 2;
	//Degrees of freedom for each node in the problem
	int dof = 3;
	//Number of nodes of the mesh
	nNodes = myMesh.nNodes;
	//Number of elements of the mesh
	nCells = myMesh.nCells;
	//number of constraints
	int nConstraints = (int) myMesh.cxyz_th(0).size() + (int) myMesh.cxyz_th(1).size() + (int) myMesh.cxyz_th(2).size();	
	//inertia relief
	if( (myMesh.IR_p==1 && myMesh.body==1) || (myMesh.IR_c==1 && myMesh.body==0) )
	{
		nConstraints = (int) cxyz_IR(0).size() + (int) cxyz_IR(1).size() + (int) cxyz_IR(2).size() ;	
	}
	//Number of loaded nodes (for the logic of the program a load is considered for all the nodes that in case of no load is set to zero)
	int nLoads = 3*nNodes;

	fout.open("./temp/piston/input/general.fem");
	fout << "// ------------ fem general input file ------------ //" << "\n"
			<< "\n" << "// Element type: 0 -> brick, 1 -> tetra" << "\n"
			<< "elmType" << "\t" << elmType << "\n"
			<< "\n" << "// Number of Gauss integration points" << "\n"
			<< "intOp" << "\t" << GaussPts << "\n"
			<< "\n" << "// Problem dimension (2D or 3D)" << "\n"
			<< "nDof" << "\t" << dof << "\n"
			<< "\n" << "// Number of elements" << "\n"
			<< "nElements" << "\t" << nCells << "\n"
			<< "\n" << "// Number of nodes" << "\n"
			<< "nNodes" << "\t" << nNodes << "\n"
			<< "\n" << "// Number of constraints" << "\n"
			<< "nc" << "\t" << nConstraints << "\n"
			<< "\n" << "// Number of surface tractions" << "\n"
			<< "nst" << "\t" << nLoads << "\n";
	fout.close();
	fout.clear();

};
void CFEMThermal::FEMThermalLoads(Array<double,1> T_body)
{

	//locals
	phi.resize(T_body.extent(0));
	phi = T_body;

	//loads
	Array<double,1> loads;
	loads.resize(3*nNodes);
	Array<double,2> C(6,6);	
	Array<double,1> epsT(6);

	//Reference temperature
	double phi_ref = 25.0;

	//Tethraedron elements
	if(Teth)
	{
		int n;
		double a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4,J,C1,
			x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
		Array<double,2> B(6,12);
		Array<double,2> B1(12,6);
		Array<double,2> C(6,6);	
		Array<double,2> load_T(12);	
		//Thermal loads for each element node
		n = 0;
		loads = 0.0;
		while(n < nCells)
		{
			//coordinates
			x1 = myMesh.xyz( myMesh.conn(n,0),0);	y1 = myMesh.xyz( myMesh.conn(n,0),1);	z1 = myMesh.xyz( myMesh.conn(n,0),2);
			x2 = myMesh.xyz( myMesh.conn(n,1),0);	y2 = myMesh.xyz( myMesh.conn(n,1),1);	z2 = myMesh.xyz( myMesh.conn(n,1),2);
			x3 = myMesh.xyz( myMesh.conn(n,2),0);	y3 = myMesh.xyz( myMesh.conn(n,2),1);	z3 = myMesh.xyz( myMesh.conn(n,2),2);
			x4 = myMesh.xyz( myMesh.conn(n,3),0);	y4 = myMesh.xyz( myMesh.conn(n,3),1);	z4 = myMesh.xyz( myMesh.conn(n,3),2);
			//a
			a1 = ( y2 - y4 )*( z3 - z4 ) - ( y3 - y4 )*( z2 - z4 ) ;
			a2 = ( y3 - y4 )*( z1 - z4 ) - ( y1 - y4 )*( z3 - z4 ) ;
			a3 = ( y1 - y4 )*( z2 - z4 ) - ( y2 - y4 )*( z1 - z4 ) ;	
			a4 = -1.0*(a1+a2+a3);
			//b
			b1 = ( x3 - x4 )*( z2 - z4 ) - ( x2 - x4 )*( z3 - z4 ) ;
			b2 = ( x1 - x4 )*( z3 - z4 ) - ( x3 - x4 )*( z1 - z4 ) ;
			b3 = ( x2 - x4 )*( z1 - z4 ) - ( x1 - x4 )*( z2 - z4 ) ;
			b4 = -1.0*(b1+b2+b3);
			//c
			c1 = ( x2 - x4 )*( y3 - y4 ) - ( x3 - x4 )*( y2 - y4 ) ;
			c2 = ( x3 - x4 )*( y1 - y4 ) - ( x1 - x4 )*( y3 - y4 ) ;
			c3 = ( x1 - x4 )*( y2 - y4 ) - ( x2 - x4 )*( y1 - y4 ) ;
			c4 = -1.0*(c1+c2+c3);
			//Jacobian determinant
			J = ( x1 - x4 )*a1 + ( y1 - y4 )*b1 + ( z1 - z4 )*c1;
			//B
			B=0.0;
			B(0,0)=a1;		B(0,3)=a2;		B(0,6)=a3;		B(0,9)=a4;
			B(1,1)=b1;		B(1,4)=b2;		B(1,7)=b3;		B(1,10)=b4;
			B(2,2)=c1;		B(2,5)=c2;		B(2,8)=c3;		B(2,11)=c4;
			B(3,0)=b1;		B(3,1)=a1;		B(3,3)=b2;		B(3,4)=a2;		B(3,6)=b3;		B(3,7)=a3;		B(3,9)=b4;		B(3,10)=a4;
			B(4,1)=c1;		B(4,2)=b1;		B(4,4)=c2;		B(4,5)=b2;		B(4,7)=c3;		B(4,8)=b3;		B(4,10)=c4;		B(4,11)=b4;
			B(5,0)=c1;		B(5,2)=a1;		B(5,3)=c2;		B(5,5)=a2;		B(5,6)=c3;		B(5,8)=a3;		B(5,9)=c4;		B(5,11)=a4;
			B*=1.0/J;
			//C
			C=0.0;
			C1 = myMesh.E(n)/((1+myMesh.v(n))*(1-2*myMesh.v(n)));
			C(0,0) = C1*(1-myMesh.v(n));			C(0,1) = C1*myMesh.v(n);		C(0,2) = C(0,1);
			C(1,0) = C(0,1);						C(1,1) = C(0,0);				C(1,2) = C(0,2);
			C(2,0) = C(0,2);						C(2,1) = C(0,2);				C(2,2) = C(0,0);
			C(3,3) = C1*(1-2.0*myMesh.v(n))/2.0;	C(4,4) = C(3,3);				C(5,5) = C(3,3);
			//Thermal strain
			epsT=0.0;
			epsT(0) = myMesh.alpha(n) * (phi(n)-phi_ref);		epsT(1) = myMesh.alpha(n) * (phi(n)-phi_ref);		epsT(2) = myMesh.alpha(n) * (phi(n)-phi_ref);
			//epsT(0) = myMesh.alpha(n) * (phi(myMesh.conn(n,0))-phi_ref);	epsT(1) = myMesh.alpha(n) * (phi(myMesh.conn(n,1))-phi_ref);	epsT(2) = myMesh.alpha(n) * (phi(myMesh.conn(n,2))-phi_ref);
			//B' * C
			B1=0.0;
			for(int i=0;i<12;i++)
			{
				for(int j=0;j<6;j++)
				{
					for(int k=0;k<6;k++)
					{
						B1(i,j) += B(k,i) * C(k,j);
					}
				}
			}
			//Thermal load
			load_T = 0.0;
			for(int i=0;i<12;i++)
			{
				for(int j=0;j<6;j++)
				{
					load_T(i) += B1(i,j) * epsT(j) * (fabs(J)/6.0);
				}
			}
			//Total thermal load
			for(int i = 0; i < 4 ; i++)
			{
				loads(3*(myMesh.conn(n,i)+1)-2-1) += load_T(3*(i+1)-2-1);
				loads(3*(myMesh.conn(n,i)+1)-1-1) += load_T(3*(i+1)-1-1);
				loads(3*(myMesh.conn(n,i)+1)-1) += load_T(3*(i+1)-1);
			}
			//Next element
			n++;
		};
	};

	//Hexahedron elements
	if(Hexa)
	{
		int n;
		double a,b,c,d,e,f,g,h,hh,C1;
		Array<double,1> xi(8);				Array<double,1> a31(8);	
		Array<double,1> yi(8);				Array<double,1> a32(8);
		Array<double,1> zi(8);				Array<double,1> a33(8);
		Array<double,2> str(8,3);			Array<double,1> a11(8);	
		Array<double,2> dNds(8,8);			Array<double,1> a12(8);		
		Array<double,2> dNdt(8,8);			Array<double,1> a13(8);		
		Array<double,2> dNdr(8,8);			Array<double,1> a21(8);
		Array<double,3> B(6,24,8);			Array<double,1> a22(8);
		Array<double,3> B1(6,24,8);			Array<double,1> a23(8);
		Array<double,2> load_T1(24,8);		Array<double,1> InvJac(8);	
		Array<double,1> load_T(24);			
			
		//Gauss quadrature matrix for numerical integration
		str(0,0) = +0.5774;		str(0,1) = -0.5774;		str(0,2) = -0.5774;
		str(1,0) = +0.5774;		str(1,1) = +0.5774;		str(1,2) = -0.5774;
		str(2,0) = -0.5774;		str(2,1) = +0.5774;		str(2,2) = -0.5774;
		str(3,0) = -0.5774;		str(3,1) = -0.5774;		str(3,2) = -0.5774;
		str(4,0) = +0.5774;		str(4,1) = -0.5774;		str(4,2) = +0.5774;
		str(5,0) = +0.5774;		str(5,1) = +0.5774;		str(5,2) = +0.5774;
		str(6,0) = -0.5774;		str(6,1) = +0.5774;		str(6,2) = +0.5774;
		str(7,0) = -0.5774;		str(7,1) = -0.5774;		str(7,2) = +0.5774;

		//Shape functions and strain-displacement functions
		for(int j=0;j<8;j++)
		{
			dNds(j,0) = (1-str(j,1))*(1-str(j,2)) / 8.0;		dNdt(j,0) = (1+str(j,0))*(str(j,2)-1) / 8.0;		dNdr(j,0) = (1+str(j,0))*(str(j,1)-1) / 8.0;
			dNds(j,1) = (1+str(j,1))*(1-str(j,2)) / 8.0;		dNdt(j,1) = (1+str(j,0))*(1-str(j,2)) / 8.0;		dNdr(j,1) = (-1-str(j,0))*(1+str(j,1)) / 8.0;
			dNds(j,2) = (1+str(j,1))*(str(j,2)-1) / 8.0;		dNdt(j,2) = (1-str(j,0))*(1-str(j,2)) / 8.0;		dNdr(j,2) = (str(j,0)-1)*(1+str(j,1)) / 8.0;
			dNds(j,3) = (str(j,1)-1)*(1-str(j,2)) / 8.0;		dNdt(j,3) = (str(j,0)-1)*(1-str(j,2)) / 8.0;		dNdr(j,3) = (str(j,0)-1)*(1-str(j,1)) / 8.0;
			dNds(j,4) = (1-str(j,1))*(1+str(j,2)) / 8.0;		dNdt(j,4) = (-1-str(j,0))*(1+str(j,2)) / 8.0;		dNdr(j,4) = (1+str(j,0))*(1-str(j,1)) / 8.0;
			dNds(j,5) = (1+str(j,1))*(1+str(j,2)) / 8.0;		dNdt(j,5) = (1+str(j,0))*(1+str(j,2)) / 8.0;		dNdr(j,5) = (1+str(j,0))*(1+str(j,1)) / 8.0;
			dNds(j,6) = (-1-str(j,1))*(1+str(j,2)) / 8.0;		dNdt(j,6) = (1-str(j,0))*(1+str(j,2)) / 8.0;		dNdr(j,6) = (1-str(j,0))*(1+str(j,1)) / 8.0;
			dNds(j,7) = (str(j,1)-1)*(1+str(j,2)) / 8.0;		dNdt(j,7) = (str(j,0)-1)*(1+str(j,2)) / 8.0;		dNdr(j,7) = (1-str(j,0))*(1-str(j,1)) / 8.0;
		}
						
		//Thermal loads for each element node
		n = 0;
		loads = 0.0;
		while(n < nCells)
		{
			//Constitutiva matrix
			C=0.0;
			C1 = myMesh.E(n)/((1+myMesh.v(n))*(1-2*myMesh.v(n)));
			C(0,0) = C1*(1-myMesh.v(n));			C(0,1) = C1*myMesh.v(n);		C(0,2) = C(0,1);
			C(1,0) = C(0,1);						C(1,1) = C(0,0);				C(1,2) = C(0,2);
			C(2,0) = C(0,2);						C(2,1) = C(0,2);				C(2,2) = C(0,0);
			C(3,3) = C1*(1-2*myMesh.v(n))/2.0;		C(4,4) = C(3,3);				C(5,5) = C(3,3);
			
			//Thermal strain
			epsT=0.0;
			epsT(0) = myMesh.alpha(n);		epsT(1) = myMesh.alpha(n);		epsT(2) = myMesh.alpha(n);

			//Sizing strain-displacmenet matrix for new element
			InvJac=0.0; a11=0.0; a12=0.0; a13=0.0; a21=0.0; a22=0.0; a23=0.0; a31=0.0; a32=0.0; a33=0.0;
			B=0.0; B1=0.0; xi=0.0; yi=0.0; zi=0.0; load_T1=0.0; load_T=0.0;
			for(int i=0;i<8;i++)
			{
				//Element absolute spatial coordinates
				xi(i) = myMesh.xyz( myMesh.conn(n,i),0);
				yi(i) = myMesh.xyz( myMesh.conn(n,i),1);
				zi(i) = myMesh.xyz( myMesh.conn(n,i),2);
			}

			//B matrix 
			for(int i=0;i<8;i++)
			{
				//Inverse Jacobian calculation
				a=0; b=0; c=0; d=0;	e=0; f=0; g=0; h=0; hh=0;										
				for(int j=0;j<8;j++)						
				{	
					a = a + dNds(i,j)*xi(j);	d = d + dNdt(i,j)*xi(j);	g = g + dNdr(i,j)*xi(j);
					b = b + dNds(i,j)*yi(j);	e = e + dNdt(i,j)*yi(j);	h = h + dNdr(i,j)*yi(j);
					c = c + dNds(i,j)*zi(j);	f = f + dNdt(i,j)*zi(j);	hh = hh + dNdr(i,j)*zi(j);
				}
				a11(i) = e*hh-f*h;	a12(i) = f*g-d*hh;	a13(i) = d*h-e*g;
				a21(i) = c*h-b*hh;	a22(i) = a*hh-c*g;	a23(i) = b*g-a*h;
				a31(i) = b*f-c*e;	a32(i) = c*d-a*f;	a33(i) = a*e-b*d;
			
				//Inverse Jacobian
				InvJac(i) =  a11(i)*a + a12(i)*b + a13(i)*c ;

				//B matrix definition for each node of the element
				for(int j=0;j<8;j++)
				{
					B(0,3*(j+1)-2-1,i) = (a11(i)*dNds(i,j) + a21(i)*dNdt(i,j) + a31(i)*dNdr(i,j)) / InvJac(i);
					B(1,3*(j+1)-1-1,i) = (a12(i)*dNds(i,j) + a22(i)*dNdt(i,j) + a32(i)*dNdr(i,j)) / InvJac(i);
					B(2,3*(j+1)-1,i) = (a13(i)*dNds(i,j) + a23(i)*dNdt(i,j) + a33(i)*dNdr(i,j)) / InvJac(i);
					B(3,3*(j+1)-2-1,i) = B(1,3*(j+1)-1-1,i);
					B(3,3*(j+1)-1-1,i) = B(0,3*(j+1)-2-1,i);
					B(4,3*(j+1)-1-1,i) = B(2,3*(j+1)-1,i);
					B(4,3*(j+1)-1,i) = B(1,3*(j+1)-1-1,i);
					B(5,3*(j+1)-2-1,i) = B(2,3*(j+1)-1,i);
					B(5,3*(j+1)-1,i) = B(0,3*(j+1)-2-1,i);
				}
			}

			//Thermal load vector definition for the element
			for(int j=0;j<8;j++)
			{
				for(int i=0;i<24;i++)
				{
					for(int k=0;k<6;k++)
					{
						for(int kk=0;kk<6;kk++)
						{
							B1(k,i,j) += B(kk,i,j) * C(kk,k);
						}
						load_T1(i,j) += B1(k,i,j)*epsT(k);
					}
					load_T1(i,j) = load_T1(i,j) * InvJac(j) * ( phi(j) - phi_ref );
				}
			}

			//Summation of thermal load vectors using gauss quadrature and weight factor 1
			load_T = sum(load_T1,tensor::j);

			//Total thermal load
			for(int i = 0; i < 8 ; i++)
			{
				//Thermal stress
				loads(3*(myMesh.conn(n,i)+1)-2-1) += load_T(3*(i+1)-2-1);
				loads(3*(myMesh.conn(n,i)+1)-1-1) += load_T(3*(i+1)-1-1);
				loads(3*(myMesh.conn(n,i)+1)-1) += load_T(3*(i+1)-1);
			}

			//Next element
			n++;
		};
	};


	//inertia relief
	if( (myMesh.IR_p==1 && myMesh.body==1) || (myMesh.IR_c==1 && myMesh.body==0) )
	{
		Log << "\n";
		Log << "---------------------------------" << "\n";
		Log << "Inertia Relief Loads ... " << "\n";
		FEMThermalInertiaRelief(loads);
		Log << "Done!" << "\n";
		Log << "---------------------------------" << "\n";
		Log << "\n";
	}


	//Output piston thermal load
	fout.open("./temp/piston/input/loads.fem");
	for(int i = 0; i < 3*nNodes ; i++)
	{
		fout << i << "\t" << scientific << loads(i) << "\n";
	}
	fout.close();
	fout.clear();

};
void CFEMThermal::FEMThermalConstraints(void)
{
	vector<int> constraints;
	int cx,cy,cz;

	//inertia relief
	if( (myMesh.IR_p==1 && myMesh.body==1) || (myMesh.IR_c==1 && myMesh.body==0) )
	{
		cx = (int) cxyz_IR(0).size();
		cy = (int) cxyz_IR(1).size();
		cz = (int) cxyz_IR(2).size();
	}
	else
	{
		//number of constarints in x y z directions
		cx = (int) myMesh.cxyz_th(0).size();
		cy = (int) myMesh.cxyz_th(1).size();
		cz = (int) myMesh.cxyz_th(2).size();
	}

	
	//assign constraints
	Array<int,1> nn_cx;
	Array<int,1> nn_cy;
	Array<int,1> nn_cz;
	Array<int,1> nn_c;
	nn_cx.resize(cx);
	nn_cy.resize(cy);
	nn_cz.resize(cz);
	nn_c.resize(cx+cy+cz);
	for(int i = 0 ; i < cx ; i++)
	{
		//inertia relief
		if( (myMesh.IR_p==1 && myMesh.body==1) || (myMesh.IR_c==1 && myMesh.body==0) )
		{
			nn_cx(i) = 3*cxyz_IR(0)[i];
		}
		else
		{
			nn_cx(i) = 3*myMesh.cxyz_th(0)[i];
		}
	}
	for(int i = 0 ; i < cy ; i++)
	{
		//inertia relief
		if( (myMesh.IR_p==1 && myMesh.body==1) || (myMesh.IR_c==1 && myMesh.body==0) )
		{
			nn_cy(i) = 3*cxyz_IR(1)[i] + 1;
		}
		else
		{
			nn_cy(i) = 3*myMesh.cxyz_th(1)[i] + 1; 
		}
	}
	for(int i = 0 ; i < cz ; i++)
	{
		//inertia relief
		if( (myMesh.IR_p==1 && myMesh.body==1) || (myMesh.IR_c==1 && myMesh.body==0) )
		{
			nn_cz(i) = 3*cxyz_IR(2)[i] + 2 ;
		}
		else
		{
			nn_cz(i) = 3*myMesh.cxyz_th(2)[i] + 2 ;
		}
	}
	nn_c(Range(0,cx-1)) = nn_cx;
	nn_c(Range(cx,cx+cy-1)) = nn_cy;
	nn_c(Range(cx+cy,cx+cy+cz-1)) = nn_cz;
	for(int i = 0; i < cx+cy+cz ; i++)
	{
		constraints.push_back(nn_c(i));
	}
	sort(constraints.begin(),constraints.end());

	//Output piston constraint matrix
	fout.open("./temp/piston/input/constraints.fem");
	for(int i = 0; i < cx+cy+cz ; i++)
	{
		fout << constraints[i] << "\n";
	}
	fout.close();
	fout.clear();

};
Array<double,1> CFEMThermal::FEMThermalSurfaceDeformation(void)
{
	//variables
	double a,b,phi,x,y;
	int nNodes_gap;
	vector<int> Nodes_gap;
	Array<double,1> def_surf;


	Nodes_gap = myThermal.nodeid_qb[0];
	nNodes_gap = (int) Nodes_gap.size();
	def_surf.resize(nNodes_gap);	def_surf=0.0;
	//calculate gap shift compensation in y-direction.
	double defshift = 0.0;
	bool shiftcoord = 1;
	for(int a=0;a<nNodes_gap;a++){
		if(myMesh.xyz(Nodes_gap[a],1)<0.0)
			shiftcoord = 0;
		defshift += udisp(Nodes_gap[a],1);
	}
	//xshift /= nNodes_gap;
	defshift /= nNodes_gap;
	centershift = defshift;
	//surface deformation from nodal displacement
	for(int i=0;i<nNodes_gap;i++)
	{
		//displacment x and y 
		a = udisp(Nodes_gap[i],0);
		b = udisp(Nodes_gap[i],1) - defshift;
		//node coordinates x and y
		x = myMesh.xyz(Nodes_gap[i],0)+1.0e-20;
		if(shiftcoord)
			y = myMesh.xyz(Nodes_gap[i],1) - (myinput.data.geometry.dB/2);
		else
			y = myMesh.xyz(Nodes_gap[i],1) +1.0e-20;
		//angle respect to axis
		phi = atan2(y,x);
		if(phi<0.0) { phi+=2*PI; };
		//deformation of film thickness [m]
		def_surf(i) = a*cos(phi) + b*sin(phi);
	};

	//Output for vtk files entire body
	def_vtk.resize(nNodes);
	def_vtk = 0.0;
	def_vtk_x.resize(nNodes);
	def_vtk_x = 0.0;
	def_vtk_y.resize(nNodes);
	def_vtk_y = 0.0;
	def_vtk_z.resize(nNodes);
	def_vtk_z = 0.0;
	for(int j=0;j<nNodes;j++)
	{
		//displacment x and y 
		double a = udisp(j,0);
		double b = udisp(j,1);
		//node coordinates x and y
		double x = myMesh.xyz(j,0) + 1.0e-20;
		double y = myMesh.xyz(j,1) + 1.0e-20;
		//angle respect to axis
		double phi = atan2(y,x);
		if(phi<0.0) { phi += 2*PI; };
		//deformation of film thickness [m]
		double def = a*cos(phi) + b*sin(phi);
		//output [microns]
		def_vtk(j) = 1.0e6 * def;
		def_vtk_x(j) = 1.0e6 * udisp(j,0);
		def_vtk_y(j) = 1.0e6 * udisp(j,1);
		def_vtk_z(j) = 1.0e6 * udisp(j,2);
	}

	return def_surf;

};
void CFEMThermal::FEMThermalOutput(double scale,string body)
{

	//file name based on revolution
	string rev_n;
	char buffer[500];
	//revolution number
	_itoa_s(myGapResult.revcounter,buffer,10);
	rev_n = buffer;


	//constrained nodes inertia relief
	Array<double,1> consx;
	Array<double,1> consy;
	Array<double,1> consz;
	if( (myMesh.IR_p==1 && myMesh.body==1) || (myMesh.IR_c==1 && myMesh.body==0) )
	{
		consx.resize(nNodes);	consx = 0.0;
		consy.resize(nNodes);	consy = 0.0;
		consz.resize(nNodes);	consz = 0.0;
		for(int i=0;i<(int)cxyz_IR(0).size();i++)
		{
			consx(cxyz_IR(0)[i]) = 1.0;
		}
		for(int i=0;i<(int)cxyz_IR(1).size();i++)
		{
			consy(cxyz_IR(1)[i]) = 1.0;
		}
		for(int i=0;i<(int)cxyz_IR(2).size();i++)
		{
			consz(cxyz_IR(2)[i]) = 1.0;
		}
	}


	if(Hexa)
	{
		//------------Mesh VTK------------//
		string output = "./output/piston/vtk/" + body + "_th_rev_" + rev_n + ".vtk";
		fout.open(output.c_str());
		fout << "# vtk DataFile Version 2.0" << "\n";
		fout << "vtk output" << "\n";
		fout << "ASCII" << "\n";
		fout << "DATASET UNSTRUCTURED_GRID" << "\n";
		fout << "POINTS " << nNodes << " double" << "\n";
		fout.precision(8);
		for(int i=0;i<nNodes;i++)
		{
			fout << scientific <<  myMesh.xyz(i,0) + udisp(i,0)*scale << "\t" <<  myMesh.xyz(i,1) + udisp(i,1)*scale  << "\t" <<  myMesh.xyz(i,2) + udisp(i,2)*scale << "\n";
		}
		fout << " " << "\n";
		fout << "CELLS " << "\t" << nCells << "\t" << 9*nCells << "\n";
		for(int i = 0; i < nCells; i++)
		{
			fout << 8 << "\t" << myMesh.conn(i,0) << "\t" << myMesh.conn(i,1)  << "\t" << myMesh.conn(i,2) << "\t" << myMesh.conn(i,3) << "\t" << myMesh.conn(i,4) << "\t" << myMesh.conn(i,5) << "\t" << myMesh.conn(i,6) << "\t" << myMesh.conn(i,7) << "\n";
		}
		fout << " " << "\n";
		fout << "CELL_TYPES " << nCells << "\n";
		for(int i = 0; i < nCells; i++)
		{
			fout << 12 << "\n";
		}
		fout << " " << "\n";
		fout << "CELL_DATA " << nCells << "\n";
		fout << "SCALARS T_[C] double" << "\n";
		fout << "LOOKUP_TABLE default" << "\n";
		for(int i = 0; i < nCells; i++)
		{
			fout << scientific << phi(i) << "\n";
		}
		fout << " " << "\n";
		fout << "POINT_DATA " << nNodes << "\n";
		fout << "SCALARS Total_Deformation_[um] double 1" << "\n";
		fout << "LOOKUP_TABLE default" << "\n";
		for(int i = 0; i < nNodes; i++)
		{
			fout << def_vtk(i) << "\n";
		}
		fout << " " << "\n";
		fout << "SCALARS x_Deformation_[um] double 1" << "\n";
		fout << "LOOKUP_TABLE default" << "\n";
		for(int i = 0; i < nNodes; i++)
		{
			fout << def_vtk_x(i) << "\n";
		}
		fout << " " << "\n";
		fout << "SCALARS y_Deformation_[um] double 1" << "\n";
		fout << "LOOKUP_TABLE default" << "\n";
		for(int i = 0; i < nNodes; i++)
		{
			fout << def_vtk_y(i) << "\n";
		}
		fout << " " << "\n";
		fout << "SCALARS z_deformation_[um] double 1" << "\n";
		fout << "LOOKUP_TABLE default" << "\n";
		for(int i = 0; i < nNodes; i++)
		{
			fout << def_vtk_z(i) << "\n";
		}
		if( (myMesh.IR_p==1 && myMesh.body==1) || (myMesh.IR_c==1 && myMesh.body==0) )
		{
			fout << " " << "\n";
			fout << "SCALARS x_Constraints double 1" << "\n";
			fout << "LOOKUP_TABLE default" << "\n";
			for(int i = 0; i < nNodes; i++)
			{
				fout << consx(i) << "\n";
			}
			fout << " " << "\n";
			fout << "SCALARS y_Constraints double 1" << "\n";
			fout << "LOOKUP_TABLE default" << "\n";
			for(int i = 0; i < nNodes; i++)
			{
				fout << consy(i) << "\n";
			}
			fout << " " << "\n";
			fout << "SCALARS z_Constraints double 1" << "\n";
			fout << "LOOKUP_TABLE default" << "\n";
			for(int i = 0; i < nNodes; i++)
			{
				fout << consz(i) << "\n";
			}
		}
		fout.close();
		fout.clear();
	}
	else
	{
		//-------------------Mesh VTK-------------------//
		string output = "./output/piston/vtk/" + body + "_th_rev_" + rev_n + ".vtk";
		fout.open(output.c_str());
		fout << "# vtk DataFile Version 2.0" << "\n";
		fout << "vtk output" << "\n";
		fout << "ASCII" << "\n";
		fout << "DATASET UNSTRUCTURED_GRID" << "\n";
		fout << "POINTS " << nNodes << " double" << "\n";
		fout.precision(8);
		for(int i=0;i<nNodes;i++)
		{
			fout << scientific <<  myMesh.xyz(i,0) + udisp(i,0)*scale  << "\t" <<  myMesh.xyz(i,1) + udisp(i,1)*scale  << "\t" <<  myMesh.xyz(i,2) + udisp(i,2)*scale  << "\n";
		}
		fout << " " << "\n";
		fout << "CELLS " << "\t" << nCells << "\t" << 5*nCells << "\n";
		for(int i = 0; i < nCells; i++)
		{
			fout << 4 << "\t" << myMesh.conn(i,0) << "\t" << myMesh.conn(i,1)  << "\t" << myMesh.conn(i,2) << "\t" << myMesh.conn(i,3) << "\n";
		}
		fout << " " << "\n";
		fout << "CELL_TYPES " << nCells << "\n";
		for(int i = 0; i < nCells; i++)
		{
			fout << 10 << "\n";
		}
		fout << " " << "\n";
		fout << "CELL_DATA " << nCells << "\n";
		fout << "SCALARS T_[C] double" << "\n";
		fout << "LOOKUP_TABLE default" << "\n";
		for(int i = 0; i < nCells; i++)
		{
			fout << scientific << phi(i) << "\n";
		}
		fout << " " << "\n";
		fout << "POINT_DATA " << nNodes << "\n";
		fout << "SCALARS Total_Deformation_[um] double 1" << "\n";
		fout << "LOOKUP_TABLE default" << "\n";
		for(int i = 0; i < nNodes; i++)
		{
			fout << def_vtk(i) << "\n";
		}
		fout << " " << "\n";
		fout << "SCALARS x_Deformation_[um] double 1" << "\n";
		fout << "LOOKUP_TABLE default" << "\n";
		for(int i = 0; i < nNodes; i++)
		{
			fout << def_vtk_x(i) << "\n";
		}
		fout << " " << "\n";
		fout << "SCALARS y_Deformation_[um] double 1" << "\n";
		fout << "LOOKUP_TABLE default" << "\n";
		for(int i = 0; i < nNodes; i++)
		{
			fout << def_vtk_y(i) << "\n";
		}
		fout << " " << "\n";
		fout << "SCALARS z_Deformation_[um] double 1" << "\n";
		fout << "LOOKUP_TABLE default" << "\n";
		for(int i = 0; i < nNodes; i++)
		{
			fout << def_vtk_z(i) << "\n";
		}
		if( (myMesh.IR_p==1 && myMesh.body==1) || (myMesh.IR_c==1 && myMesh.body==0) )
		{
			fout << " " << "\n";
			fout << "SCALARS x_Constraints double 1" << "\n";
			fout << "LOOKUP_TABLE default" << "\n";
			for(int i = 0; i < nNodes; i++)
			{
				fout << consx(i) << "\n";
			}
			fout << " " << "\n";
			fout << "SCALARS y_Constraints double 1" << "\n";
			fout << "LOOKUP_TABLE default" << "\n";
			for(int i = 0; i < nNodes; i++)
			{
				fout << consy(i) << "\n";
			}
			fout << " " << "\n";
			fout << "SCALARS z_Constraints double 1" << "\n";
			fout << "LOOKUP_TABLE default" << "\n";
			for(int i = 0; i < nNodes; i++)
			{
				fout << consz(i) << "\n";
			}
		}
	}

	fout.close();
	fout.clear();

	//free output containers
	def_vtk.free(); def_vtk_x.free(); def_vtk_y.free(); def_vtk_z.free();
	consx.free(); consy.free(); consz.free(); 

};
//Inertia relief
void CFEMThermal::FEMThermalInertiaConstraints(void)
{

	//Center of gravity coordinates
	xcg = myThermal.xcg;
	ycg = myThermal.ycg;
	zcg = myThermal.zcg;


	Log << "------------------------------------------"  << "\n";
	Log << "Inertia Relief Constraints ... " << "\n";
	Log << "Searching nodes closer to CG ..." << "\n";
	double tempdouble = xcg*1.0e3;
	Log << "xCG: " << tempdouble << " [mm]" << "\n";
	tempdouble = ycg*1.0e3;
	Log << "yCG: " << tempdouble << " [mm]" << "\n";
	tempdouble = zcg*1.0e3;
	Log << "zCG: " << tempdouble << " [mm]" << "\n";

	
	int dim;
	double tol;
	ANNpointArray dataPts;	//Data points
	ANNpoint queryPt;		//Query point
	ANNidxArray nnIdx;		//Near neighbor indices
	ANNdistArray dists;		//Near neighbor distances
	ANNkd_tree* kdTree;		//Search structure
	
	//Numer of nodes to consider 
	nNodes = myMesh.nNodes;
	int ncg = (int) nNodes/10;

	//KD-tree
	dim = 3;
	tol = 1.0e-9;
	//Allocate data points
	dataPts = annAllocPts(nNodes,dim);
	//Allocate query pt
	queryPt = annAllocPt(dim); 

	//Face cooridnates of all the faces
	for(int i=0;i<nNodes;i++)
	{
		dataPts[i][0] = myMesh.xyz(i,0);
		dataPts[i][1] = myMesh.xyz(i,1);
		dataPts[i][2] = myMesh.xyz(i,2);
	}

	//Creating kd tree
	Log << "Creating KD tree ..." << " ";
	kdTree = new ANNkd_tree(dataPts,nNodes,dim);
	nnIdx = new ANNidx[ncg];		
	dists = new ANNdist[ncg];	
	Log << "done!" << "\n";

	//Searching the kd tree for 2 closest centers to query face center: one is same center other is duplicate
	Log << "Searching KD tree ..." << " ";
	//Face center under investigation
	queryPt[0] = xcg;
	queryPt[1] = ycg;
	queryPt[2] = zcg;

	//Search kd-tree
	kdTree -> annkSearch(queryPt,ncg,nnIdx,dists);
	Log << "finished searching!" << "\n";

	//Minimum distance
	double dist_min = 10e10;
	for(int i=0;i<ncg;i++)
	{
		if(dists[i]<dist_min)
		{
			dist_min = dists[i];
		}
	};
	dist_min = sqrt(dist_min);

	Log << "Assigning IR constraints...";
	
	//Look for nodes to constrain close to CG
	vector<int> nnIdx_min;
	cxyz_IR.resize(3);
	cxyz_IR(0).clear();
	cxyz_IR(1).clear();
	cxyz_IR(2).clear();
	double d_tol;
	if(myMesh.body)
	{
		 d_tol = myMesh.d_tol_p;
	}
	else
	{
		 d_tol = myMesh.d_tol_c;
	}

	//Nodes closer to CG are IR constraints
	for(int i=0;i<ncg;i++)
	{
		if( sqrt(dists[i]) < (dist_min+d_tol) )
		{
			nnIdx_min.push_back(nnIdx[i]);
		}
	}

	//Assign the minimum number of constraints to avoid rigid body motion
	double d = 0.0;
	double d_c = myThermal.d_c;
	for(int i=0;i<(int) nnIdx_min.size();i++)
	{
		//x-constraint
		d = fabs( myMesh.xyz(nnIdx_min[i],0) - xcg );
		if(d<d_c)
		{
			cxyz_IR(0).push_back(nnIdx_min[i]);
		};
		//y-constraint
		d = fabs( myMesh.xyz(nnIdx_min[i],1) - ycg );
		if(d<d_c)
		{
			cxyz_IR(1).push_back(nnIdx_min[i]);
		};
		//z-constraint
		d = fabs( myMesh.xyz(nnIdx_min[i],2) - zcg );
		if(d<d_c)
		{
			cxyz_IR(2).push_back(nnIdx_min[i]);
		};
	};

	//Clean and free kdtree
	annDeallocPt(queryPt);
	annDeallocPts(dataPts);
	delete [] nnIdx;
	delete [] dists;
	delete kdTree;
	annClose();

	Log << " Done!" << "\n";
	Log << "------------------------------------------" << "\n";

};
Array<double,2> CFEMThermal::FEMThermalCellInertia(int i)
{

	//Moment of inertia
	Array<double,2> I(3,3);	I=0.0;
	
	//Principal
	I(0,0) = M(i) * ( pow(myThermal.xyzc(i,1)-ycg,2.0) + pow(myThermal.xyzc(i,2)-zcg,2.0) );
	I(1,1) = M(i) * ( pow(myThermal.xyzc(i,0)-xcg,2.0) + pow(myThermal.xyzc(i,2)-zcg,2.0) );
	I(2,2) = M(i) * ( pow(myThermal.xyzc(i,0)-xcg,2.0) + pow(myThermal.xyzc(i,1)-ycg,2.0) );
	//Other
	I(0,1) = M(i) * ( myThermal.xyzc(i,0)-xcg ) * ( myThermal.xyzc(i,1)-ycg ) ;
	I(0,2) = M(i) * ( myThermal.xyzc(i,0)-xcg ) * ( myThermal.xyzc(i,2)-zcg ) ;
	I(1,2) = M(i) * ( myThermal.xyzc(i,1)-ycg ) * ( myThermal.xyzc(i,2)-zcg ) ;
	I(0,1) *= -1.0;
	I(0,2) *= -1.0;
	I(1,2) *= -1.0;
	I(1,0) = I(0,1);
	I(2,0) = I(0,2);
	I(2,1) = I(1,2);

	return I;

}
Array<double,2> CFEMThermal::FEMThermalInertiaInverse(Array<double,2> A)
{
	Array<double,2> B(3,3);
	B=0.0;

	//Determinant of A
	double det = A(0,0)*( A(1,1)*A(2,2)-A(2,1)*A(1,2) )
                 - A(0,1)*( A(1,0)*A(2,2)-A(1,2)*A(2,0) )
                 + A(0,2)*( A(1,0)*A(2,1)-A(1,1)*A(2,0) );
	
	//Inverse of A
	double invdet = 1.0/det;
	B(0,0) =  ( A(1,1)*A(2,2) - A(2,1)*A(1,2) ) * invdet;
	B(1,0) = -( A(0,1)*A(2,2) - A(0,2)*A(2,1) ) * invdet;
	B(2,0) =  ( A(0,1)*A(1,2) - A(0,2)*A(1,1) ) * invdet;
	B(0,1) = -( A(1,0)*A(2,2) - A(1,2)*A(2,0) ) * invdet;
	B(1,1) =  ( A(0,0)*A(2,2) - A(0,2)*A(2,0) ) * invdet;
	B(2,1) = -( A(0,0)*A(1,2) - A(1,0)*A(0,2) ) * invdet;
	B(0,2) =  ( A(1,0)*A(2,1) - A(2,0)*A(1,1) ) * invdet;
	B(1,2) = -( A(0,0)*A(2,1) - A(2,0)*A(0,1) ) * invdet;
	B(2,2) =  ( A(0,0)*A(1,1) - A(1,0)*A(0,1) ) * invdet;

	return B;
}
void CFEMThermal::FEMThermalInertiaRelief(Array<double,1> &loads)
{
	Range all = Range::all();

	//Nodes per element
	double ne;
	if(Teth)
	{
		ne=4;
	}
	else
	{
		ne=8;
	}

	//Mass of each cell array
	M.resize(nCells);
	M = myThermal.M;
	//Total mass of the body
	Mtot = myThermal.Mtot;

	//Calculate thermal loads and moments respect CG
	Array<double,1> T(3);			T = 0.0;
	Array<double,1> T_i(3);			T_i = 0.0;
	Array<double,1> F(3);			F = 0.0;
	Array<double,1> F_i(3);			F_i = 0.0;
	//Sum over each loaded mesh face
	for(int i=0;i<nNodes;i++)
	{
		//Load vector face-i
		F_i(0) = loads(3*i);
		F_i(1) = loads(3*i+1);
		F_i(2) = loads(3*i+2);
		//Moment arm for face-i
		double rx = myMesh.xyz(i,0) - xcg;
		double ry = myMesh.xyz(i,1) - ycg;
		double rz = myMesh.xyz(i,2) - zcg;
		//Torque for face-i
		T_i = 0.0;
		T_i(0) = -1.0 * F_i(1) * rz + F_i(2) * ry;
		T_i(1) = F_i(0) * rz - 1.0 * F_i(2) * rx;
		T_i(2) = -1.0 * F_i(0) * ry + F_i(1) * rx;
		//Total thermal load
		F += F_i;
		//Total thermal torque
		F += F_i;
	}
	Log << "Body Total Mass: " << Mtot << " [kg]" << "\n";
	Log << "Thermal Load:" << "\n";
	Log << "	Fx: " << F(0) << " [N]" << "\n";
	Log << "	Fy: " << F(1) << " [N]" << "\n";
	Log << "	Fz: " << F(2) << " [N]" << "\n";
	Log << "Thermal Torque:" << "\n";
	Log << "	Tx: " << T(0) << " [Nm]" << "\n";
	Log << "	Ty: " << T(1) << " [Nm]" << "\n";
	Log << "	Tz: " << T(2) << " [Nm]" << "\n";


	//Calculate inertia tensor inverse
	Array<double,2> I_Inv(3,3);	I_Inv=0.0;
	I_Inv = FEMThermalInertiaInverse(myThermal.I);


	//Angular acceleration at CG associated with thermal torques
	Array<double,1> A_ang(3);	A_ang=0.0;
	for(int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
		{
			A_ang(i) += I_Inv(i,j) * T(j);
		}
	}
	//Negative cause opposing load moments
	A_ang *= -1.0;


	//Linear acceleration associated with thermal loads
	Array<double,1> A_lin(3);	A_lin=0.0;
	for(int i=0;i<3;i++)
	{
		A_lin(i) = F(i) / Mtot;
	}
	//Negative cause opposing load forces
	A_lin *= -1.0;

	//Output
	Log << "Linear Inertial Acceleration:" << "\n";
	Log << "	Ax: " << A_lin(0) << " [m/s2]" << "\n";
	Log << "	Ay: " << A_lin(1) << " [m/s2]" << "\n";
	Log << "	Az: " << A_lin(2) << " [m/s2]" << "\n";
	Log << "Angular Inertial Acceleration:" << "\n";
	Log << "	ax: " << A_ang(0) << " [rad/s2]"<< "\n";
	Log << "	ay: " << A_ang(1) << " [rad/s2]"<< "\n";
	Log << "	az: " << A_ang(2) << " [rad/s2]"<< "\n";


	//Angular opposing force and torque for each cell
	Array<double,1> T_cell(3);		T_cell = 0.0;
	Array<double,2> I_cell(3,3);	I_cell = 0.0;
	Array<double,1> F_ang(3);		F_ang = 0.0;
	Array<double,1> F_lin(3);		F_lin = 0.0;
	Array<double,1> F_e(3);			F_e = 0.0;
	Array<double,1> T_tot(3);		T_tot=0.0;
	Array<double,1> F_tot(3);		F_tot=0.0;
	for(int i=0;i<nCells;i++)
	{
		//Cell linear inertia force
		F_lin = 0.0;
		F_lin(0) = M(i) * A_lin(0);
		F_lin(1) = M(i) * A_lin(1);
		F_lin(2) = M(i) * A_lin(2);

		//Cell inertia tensor
		I_cell = 0.0;
		I_cell = FEMThermalCellInertia(i);
		//Cell inertia torque
		T_cell=0.0;
		for(int j=0;j<3;j++)
		{
			for(int k=0;k<3;k++)
			{
				T_cell(j) += I_cell(j,k) * A_ang(k) ;
			}
		}
		//Cell center moment arms
		double rx = myThermal.xyzc(i,0) - xcg;
		double ry = myThermal.xyzc(i,1) - ycg;
		double rz = myThermal.xyzc(i,2) - zcg;
		//Cell inertia torque components
		double Tx = T_cell(0);
		double Ty = T_cell(1);
		double Tz = T_cell(2);
		//Cell angular inertia force from vector triple product and angular force-position vector orthogonality
		F_ang = 0.0;
		double norm_r = sqrt(rx*rx+ry*ry+rz*rz);
		F_ang(0) = 1.0/pow(norm_r,2.0) * (Ty*rz-Tz*ry) ;
		F_ang(1) = 1.0/pow(norm_r,2.0) * (Tz*rx-Tx*rz) ;
		F_ang(2) = 1.0/pow(norm_r,2.0) * (Tx*ry-Ty*rx) ;


		//Total cell inertia forces
		F_e = 0.0;
		F_e = ( F_lin + F_ang ) / ne;

		//Correct loads array to consider inertial contribution
		for(int j=0;j<ne;j++)
		{
			//x
			loads(3*myMesh.conn(i,j)) += F_e(0);
			//y
			loads(3*myMesh.conn(i,j)+1) += F_e(1);
			//z
			loads(3*myMesh.conn(i,j)+2) += F_e(2);
		}

		//Total inertia torque back calculations to check balance
		T_tot(0) += -1.0 * F_ang(1) * rz + F_ang(2) * ry;
		T_tot(1) += F_ang(0) * rz - 1.0 * F_ang(2) * rx;
		T_tot(2) += -1.0 * F_ang(0) * ry + F_ang(1) * rx;
		//Total inertia force back calculations to check balance
		F_tot(0) += F_lin(0);
		F_tot(1) += F_lin(1);
		F_tot(2) += F_lin(2);
	};


	//Correct geometry file for IR calculations using CG as origin
	ofstream fout("./temp/piston/input/geometry.fem");
	for(int i = 0; i < nNodes ; i++)
	{
		fout << scientific << myMesh.xyz(i,0) - xcg << "\t" << myMesh.xyz(i,1) - ycg << "\t" << myMesh.xyz(i,2) - zcg << "\n";
	}
	fout.close();
	fout.clear();

	//Check balance of external and inertia torques and forces
	fout.open("./temp/piston/IR_th_check.dat",ios::app);
	fout << T_tot(0)+T(0) << " " << T_tot(1)+T(1) << " " << T_tot(2)+T(2) << "\n";
	fout << F_tot(0)+F(0) << " " << F_tot(1)+F(1) << " " << F_tot(2)+F(2) << "\n";
	fout.close();
	fout.clear();


	myThermal.M.free(); M.free(); myThermal.xyzc.free();		

	//Array<double,1> T_ang(3);		T_ang=0.0;
	/*//Check balance
	T_ang(0) += -1.0 * F_ang(1) * rz + F_ang(2) * ry;
	T_ang(1) += F_ang(0) * rz - 1.0 * F_ang(2) * rx;
	T_ang(2) += -1.0 * F_ang(0) * ry + F_ang(1) * rx;*/
	/*ofstream fout("./delta.dat",ios::app);
	fout << T_ang(0)+T(0) << " " << T_ang(1)+T(1) << " " << T_ang(2)+T(2) << "\n";
	fout.close();*/

};