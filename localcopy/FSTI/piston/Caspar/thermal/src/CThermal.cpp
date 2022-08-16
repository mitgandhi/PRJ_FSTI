#include "CThermal.h"
//#include "../../main/src/CPistonGap.h"
#include "../../main/src/logger.h"
#include <iostream>
#include <iomanip>
#include "../../caspar_input/input.h"
#include <cmath>
#include <vector>
#include <string>

#pragma once

extern class CMesh myMesh;
extern struct sGapResult myGapResult;
extern class CGapInput myGapInput;
extern class CGapUtils myGapUtils;
extern class input myinput;

using namespace std;

CThermal::CThermal()
{
	
};
CThermal::~CThermal()
{
	
};
void CThermal::readMeshThermal(string body)
{
	string mesh_path, mesh_name, mesh_type;

	//Mesh type
	myMesh.Teth = 0;
	myMesh.Hexa = 0;
	//Solid body: piston=1,cylinder=0
	myMesh.body = 0;

	//Define paths according to solid body
	if(body=="piston")
	{
		mesh_path = myinput.data.thermal.piston.meshfile; //myMesh.piston_path_th;
		//mesh_name = myMesh.piston_name;
		//mesh_type = myMesh.piston_type;
		myMesh.body = 1;
	}
	else if(body=="cylinder")
	{
		mesh_path = myinput.data.thermal.block.meshfile; // myMesh.cylinder_path_th;
		//mesh_name = myMesh.cylinder_name;
		//mesh_type = myMesh.cylinder_type;
		myMesh.body = 0;
	}
	
	//Read mesh
	//if(mesh_type=="Teth")
	//{
	myMesh.Teth = 1;	myMesh.Hexa = 0;
	Teth = myMesh.Teth; Hexa = myMesh.Hexa;
	myMesh.ReadABAQUSTeth(mesh_path);

	//Shift coordinate system so origin is on plane defined by bottom of bushing gap surface (closest to valve plate)
	double shift = 0.0;
	if(body=="cylinder"){
		shift = myMesh.xyz(myMesh.qbABQS[0][0],2);
		for(int j = 0;j<myMesh.qbABQS[0].size();j++){
			double zval = myMesh.xyz(myMesh.qbABQS[0][j],2);
			if(zval<shift)
				shift = zval;
		}
		for(int i = 0;i<myMesh.xyz.extent(0);i++)
			myMesh.xyz(i,2) = myMesh.xyz(i,2) - shift;
	}
	//}
	//else
	//{
	//	myMesh.Hexa = 1;	myMesh.Teth = 0;
	//	Teth = myMesh.Teth; Hexa = myMesh.Hexa;
	//	myMesh.ReadABAQUSHexa(mesh_path,mesh_name);
	//};

	//Assign mesh parameters
	nNodes=myMesh.nNodes;
	nCells=myMesh.nCells;
	if(Teth)
	{
		nn=3;
		nf=4;
	}
	else
	{
		nn=4;
		nf=6;
	}

	//Define faces centers for mesh
	FaceCenters();

	//Determine net faces and list of faces and shared cells
	Face2Cells();

	//Check solid mesh is fully bounded
	CheckBoundaryFaces(body);

	//Define boundary faces
	BoundaryFaces();

	//Calculate cell centers
	CellCenters();

	//Determine faces normal area vector
	FaceAf();

	//Write xyz coordinates gap surface
	if(myGapResult.revcounter==0)
	{
		WriteGapSurfacexyz(myMesh.body);
	};
};
void CThermal::ThermalSolve(string body,vector<Array<double,1>> qbi,Array<double,1> &T_body,Array<double,1> &T_surface)
{

	Log << "\n";
	Log << "Solving " + body + " temperature distribution... " << "\n";
	Log << "\n";

	/*for(int a = 0; a < qbi.size(); a++){
		Log << "Boundary " << a << "\n";
		for(int b = 0; b < qbi[a].size(); b++){
			Log << qbi[a](b) << "\n";
		}
	}*/

	/*for(short i = 0;i<myinput.data.thermal.piston.mixed_bc.size();i++){
		Log << myinput.data.thermal.piston.mixed_bc[i].setname << "\n";
	}*/
	

	string mesh_name;

	//Define paths according to solid body
	if(body=="piston")
	{
		mesh_name = myMesh.piston_name;
	}
	else if(body=="cylinder")
	{
		mesh_name = myMesh.cylinder_name;
	}
	

	//Size and initialize scalar field to be solved
	phi.resize(nNodes,0);
	
	//Log << nCells << "\n";
	//Log << myMesh.conn.size() << "\n";
	for(int a = 0; a < nNodes; a++){
		double temp = T_body(a);
		//phi[a] = temp;
	//	Log << temp << "\n";
		for(int b = 0; b < 4; b++){
			phi[myMesh.conn(a,b)] = temp;
		}
	}

	//vector<short> boundaries;
	boundaries.resize(nNodes,0);
	//Set boundaries
	SetBoundaries(qbi);
	

	//Develop Stiffness Matrix
	Log << "Calculating Stiffness Matrix... ";
	gmm::row_matrix<gmm::wsvector<double > > matK(nNodes,nNodes);
	gmm::clear(matK);
	for(int cell = 0;cell<nCells;cell++){
		/*Log << "Cell: " << cell << "\n";
		Log << "\n";*/
		double J = SetBMatrix(cell);
		/*Log << "J: " << J << "\n";
		Log << "\n";*/
		vector<vector<double> > Kcd;
		Kcd = CalcKcd(matB,myMesh.Tf(cell),J);//calculate the conduction matrix
		for(int a = 0;a<4;a++){
			for(int b = 0;b<4;b++){
				//if(myMesh.conn(cell,a)>=myMesh.conn(cell,b))
					matK(myMesh.conn(cell,a),myMesh.conn(cell,b)) = matK(myMesh.conn(cell,a),myMesh.conn(cell,b)) + Kcd[a][b];
			}
		}
		//Log << "\n";
	}


	//add convective terms for mixed boundary faces.
	for(int a = 0;a<faceid_mixb.size();a++){//for each mixed boundary face
		//Log << "Applying mixed boundary face " << a << "\n";
		//double phi = myMesh.phibinfi[a];
		double h = myMesh.hbi[a];
		for(int b = 0;b<faceid_mixb[a].size();b++){
			int faceid = faceid_mixb[a][b];
			//cout << "Face: " << faceid << "\n";
			//cout << Af.size() << "\n";
			double dA = FaceArea(connf(faceid,0),connf(faceid,1),connf(faceid,2));
			//cout << "Area: " << dA << "\n";
			//cout << connallf.size() << "\n";
			for(int c = 0;c<3;c++){
				int nodeid = connf(faceid,c);
				//cout << "Node: " << nodeid << "\n";
				matK(nodeid,nodeid) = matK(nodeid,nodeid) + h*dA/6;
				if(c != 0){
					//if(connallf(faceid,0) <= nodeid)
						matK(nodeid,connf(faceid,0)) = matK(nodeid,connf(faceid,0)) + h*dA/12;
				}
				if(c != 1){
					//if(connallf(faceid,1) <= nodeid)
						matK(nodeid,connf(faceid,1)) = matK(nodeid,connf(faceid,1)) + h*dA/12;
				}
				if(c != 2){
					//if(connallf(faceid,2) <= nodeid)
						matK(nodeid,connf(faceid,2)) = matK(nodeid,connf(faceid,2)) + h*dA/12;
				}
			}
		}
	}



	Log << "Done!" << "\n";

	//Develop Load Matrix
	Log << "Setting Load Vector... ";
	vecq.clear();
	vecq.resize(nNodes,0);
	
	
	for(int a = 0;a<faceid_mixb.size();a++){//for each mixed boundary face
		double phiint = myMesh.phibinfi[a];
		double h = myMesh.hbi[a];
		for(int b = 0;b<faceid_mixb[a].size();b++){
			int faceid = faceid_mixb[a][b];
			//Log << "Face: " << faceid << "\t" ;
			double dA = FaceArea(connf(faceid,0),connf(faceid,1),connf(faceid,2));
			//Log << "Area: " << dA << "\tPhiinfinite: " << phiint << "\th: " << h << "\n";;
			for(int c = 0;c<3;c++){
				vecq[connf(faceid,c)] = vecq[connf(faceid,c)] + h*phiint*dA/3;
				boundaries[connf(faceid,c)] = 2;
			}
		}
	}

	for(int a = 0;a<faceid_qb.size();a++){//for each flux boundary face
		for(int b = 0;b<faceid_qb[a].size();b++){
			int faceid = faceid_qb[a][b];
			double q = qb(faceid);
			double dA = FaceArea(connf(faceid,0),connf(faceid,1),connf(faceid,2));
			//Log << scientific << "Face: " << faceid << "\tArea: " << dA << "\tq: " << q << "\n";
			for(int c = 0;c<3;c++){
				vecq[connf(faceid,c)] = vecq[connf(faceid,c)] + q*dA/3;
				boundaries[connf(faceid,c)] = 1;
			}
		}
	}
	//dwm debug
	/*for(int i = 0; i<nNodes;i++){
		for( int j = 1+i; j < nNodes ; j++ ){
			if(matK(i,j) != 0.0){
				Log << "K(" << i << "," << j << ") = " << matK(i,j) << "\n";
				system("PAUSE");
			}
		}
	}*/
	/*Log << "0: " << matK(0,0) << "\t" << matK(0,1) << "\n";
	for( int i = 1; i<nNodes-1; i++){
		Log << i << ": " << matK(i,i-1) << "\t" << matK(i,i) << "\t" << matK(i,i+1) << "\n";
		system("PAUSE");
	}
	Log << nNodes << ": " << matK(nNodes,nNodes-1) << "\t" << matK(nNodes,nNodes);
	system("PAUSE");*/
	//vecq[1] = 0.05;
	
	Log << "Done!" << "\n";

	//Set Constraints (phib)
	Log << "Setting Temperature Constraints... " ;

	for(int a = 0;a<faceid_phib.size();a++){//for each fixed temperature face
		for(int b = 0;b<faceid_phib[a].size();b++){
			int faceid = faceid_phib[a][b];
			double phiint = phib(faceid);
			for(int c = 0;c<3;c++){
				int node = connf(faceid,c);
				//phi[node] = phiint; // Set Temperature
				vecq[node] = phiint;
				T_body(node) = phiint;//eliminate under-relaxation for this node.
				for(int i = 0;i<nNodes;i++){
					//matK(i,node) = 0;
					if(i == node)
						matK(node,i) = 1;
					else
						matK(node,i) = 0;
				}
			}
		}
	}

	/*Log << "global K = " << "\n";
	for(int a = 0;a<nNodes;a++){
		for(int b = 0;b<nNodes;b++){
			Log << matK[a][b] << "\t";
		}
		Log << "\n";
	}
	Log << "\n";

	Log << "q vector" << "\n";
	for(int a = 0;a<nNodes;a++){
		Log << vecq[a] << "\n";
	}*/
	/*for(int a = 0;a<faceid_phib.size();a++){//for each fixed temperature face
		for(int b = 0;b<faceid_phib[a].size();b++){
			int faceid = faceid_phib[a][b];
			for(int c = 0;c<3;c++){
				int node = connf(faceid,c);
				matK(node,node) = 1; //pivot this column and row to prevent changes.
			}
		}
	}*/
	
	connallf.free();
	Log << "Done!" << "\n";
	gmm::csr_matrix<double> Kint;
	//gmm::clean(Kint, 1E-12);
	gmm::copy(matK,Kint);
	gmm::clear(matK);
	
	//Solve System
	double tol;
	if(body == "cylinder")
		tol = myinput.data.thermal.block.solver.tol;
	else
		tol = myinput.data.thermal.piston.solver.tol;

	gmm::iteration iter(tol);
	if(body == "cylinder"){
		iter.set_maxiter(myinput.data.thermal.block.solver.maxit);
	}
	else if(body == "piston"){
		iter.set_maxiter(myinput.data.thermal.piston.solver.maxit);
	}
	else{
		Log << body << " is an invalid entry!" << "\n";
		system("PAUSE");
	}
	if(coupled_caspar_simulation)
		iter.set_noisy(0);
	else
		iter.set_noisy(1);
	gmm::identity_matrix BS;
	gmm::diagonal_precond<gmm::csr_matrix<double > >P(Kint);
	gmm::cg(Kint, phi, vecq, BS, P, iter);

	if(iter.converged())
		Log << "Temperature Distribution Complete!" << "\n";
	else{
		Log << "Temperature Distribution Failed!" << endl;
		exit(1);
	}

	
	//---------------------------------------------------------------------------------------------------------
	//end dwm

	//Assign results
	T_body.resize(nCells);
	Array<double,1> T_body_old = T_body;
	for(int i = 0;i<nCells;i++)
		T_body(i) = 0.25*(phi[myMesh.conn(i,0)] + phi[myMesh.conn(i,1)] + phi[myMesh.conn(i,2)] + phi[myMesh.conn(i,3)]);
	T_body = T_body_old + myinput.data.options_piston.numeric.AlphaTh * (T_body - T_body_old);
	T_body_old.free();
	
	T_surface.resize(faceid_qb[0].size());
	Array<double,1> T_surf_old = T_surface;
	for(int i = 0;i<faceid_qb[0].size();i++){
		//double temp = phi[connf(i,0)]/3 + phi[connf(i,1)]/3 + phi[connf(i,2)]/3;
		int id = faceid_qb[0][i];
		T_surface(i) = (phi[connf(id,0)] + phi[connf(id,1)] + phi[connf(id,2)]) / 3.0;
		//Log << "T_surface(" << i << ") = " << temp << "\n";
	}
	T_surface = T_surf_old + myinput.data.options_piston.numeric.AlphaTh * (T_surface - T_surf_old);
	T_surf_old.free();


	//Output VTK and others
	ThermalOutput(body);

	

	//free arrays
	connf.free(); Afn.free(); phi.clear(); phi_surf.free(); faceid_qb.clear(); faceid_mixb.clear();


	Log << "\n";
	Log << "Done " + body + " temperature distribution!" << "\n";
	Log << "\n";


};
// add conduction matrix for piston body
void CThermal::ThermalSolve_Piston(vector<Array<double,1>> qbi, vector<Array<double,1>> wtd_h, vector<double> Tb, Array<double,1> &T_body,Array<double,1> &T_surface)
{

	string mesh_name;
	mesh_name = myMesh.piston_name;
	
	//Size and initialize scalar field to be solved
	phi_both.resize(nNodes,0);
	
	//Log << nCells << "\n";
	//Log << myMesh.conn.size() << "\n";
	for(int a = 0; a < nNodes; a++){
		double temp = T_body(a);
		//phi[a] = temp;
	//	Log << temp << "\n";
		for(int b = 0; b < 4; b++){
			phi_both[myMesh.conn(a,b)] = temp;
		}
	}

	//vector<short> boundaries;
	boundaries.resize(nNodes,0);
	//Set boundaries
	SetBoundaries(qbi);
	

	//Develop Stiffness Matrix
	Log << "Calculating Stiffness Matrix... ";
	gmm::resize(matK,nNodes,nNodes);
	//gmm::row_matrix<gmm::wsvector<double > > matK(nNodes,nNodes);
	gmm::clear(matK);
	for(int cell = 0;cell<nCells;cell++){
		/*Log << "Cell: " << cell << "\n";
		Log << "\n";*/
		double J = SetBMatrix(cell);
		/*Log << "J: " << J << "\n";
		Log << "\n";*/
		vector<vector<double> > Kcd;
		Kcd = CalcKcd(matB,myMesh.Tf(cell),J);//calculate the conduction matrix
		for(int a = 0;a<4;a++){
			for(int b = 0;b<4;b++){
				//if(myMesh.conn(cell,a)>=myMesh.conn(cell,b))
					matK(myMesh.conn(cell,a),myMesh.conn(cell,b)) = matK(myMesh.conn(cell,a),myMesh.conn(cell,b)) + Kcd[a][b];
			}
		}
		//Log << "\n";
	}


	//add convective terms for mixed boundary faces.
	for(int a = 0;a<faceid_mixb.size();a++){//for each mixed boundary face
		//Log << "Applying mixed boundary face " << a << "\n";
		//double phi = myMesh.phibinfi[a];
		double h = myMesh.hbi[a];
		for(int b = 0;b<faceid_mixb[a].size();b++){
			int faceid = faceid_mixb[a][b];
			//cout << "Face: " << faceid << "\n";
			//cout << Af.size() << "\n";
			double dA = FaceArea(connf(faceid,0),connf(faceid,1),connf(faceid,2));
			//cout << "Area: " << dA << "\n";
			//cout << connallf.size() << "\n";
			for(int c = 0;c<3;c++){
				int nodeid = connf(faceid,c);
				//cout << "Node: " << nodeid << "\n";
				matK(nodeid,nodeid) = matK(nodeid,nodeid) + h*dA/6;
				if(c != 0){
					//if(connallf(faceid,0) <= nodeid)
						matK(nodeid,connf(faceid,0)) = matK(nodeid,connf(faceid,0)) + h*dA/12;
				}
				if(c != 1){
					//if(connallf(faceid,1) <= nodeid)
						matK(nodeid,connf(faceid,1)) = matK(nodeid,connf(faceid,1)) + h*dA/12;
				}
				if(c != 2){
					//if(connallf(faceid,2) <= nodeid)
						matK(nodeid,connf(faceid,2)) = matK(nodeid,connf(faceid,2)) + h*dA/12;
				}
			}
		}
	}

	//add convective terms for mixed heat transfer on flux boundary faces.
	for(int a = 0;a<2;a++){//for each mixed boundary face
		//Log << "Applying mixed boundary face " << a << "\n";
		//double phi = myMesh.phibinfi[a];
		for(int b = 0;b<faceid_qb[0].size();b++){
			int faceid = faceid_qb[0][b];
			double h = wtd_h[a](b);
			//cout << "Face: " << faceid << "\n";
			//cout << Af.size() << "\n";
			double dA = FaceArea(connf(faceid,0),connf(faceid,1),connf(faceid,2));
			//cout << "Area: " << dA << "\n";
			//cout << connallf.size() << "\n";
			for(int c = 0;c<3;c++){
				int nodeid = connf(faceid,c);
				//cout << "Node: " << nodeid << "\n";
				matK(nodeid,nodeid) = matK(nodeid,nodeid) + h*dA/6;
				if(c != 0){
					//if(connallf(faceid,0) <= nodeid)
						matK(nodeid,connf(faceid,0)) = matK(nodeid,connf(faceid,0)) + h*dA/12;
				}
				if(c != 1){
					//if(connallf(faceid,1) <= nodeid)
						matK(nodeid,connf(faceid,1)) = matK(nodeid,connf(faceid,1)) + h*dA/12;
				}
				if(c != 2){
					//if(connallf(faceid,2) <= nodeid)
						matK(nodeid,connf(faceid,2)) = matK(nodeid,connf(faceid,2)) + h*dA/12;
				}
			}
		}
	}

	Log << "Done!" << "\n";

	//Develop Load Matrix
	Log << "Setting Load Vector... ";
	vecq.clear();
	vecq.resize(nNodes,0);
	
	
	for(int a = 0;a<faceid_mixb.size();a++){//for each mixed boundary face
		double phiint = myMesh.phibinfi[a];
		double h = myMesh.hbi[a];
		for(int b = 0;b<faceid_mixb[a].size();b++){
			int faceid = faceid_mixb[a][b];
			//Log << "Face: " << faceid << "\t" ;
			double dA = FaceArea(connf(faceid,0),connf(faceid,1),connf(faceid,2));
			//Log << "Area: " << dA << "\tPhiinfinite: " << phiint << "\th: " << h << "\n";;
			for(int c = 0;c<3;c++){
				vecq[connf(faceid,c)] = vecq[connf(faceid,c)] + h*phiint*dA/3;
				boundaries[connf(faceid,c)] = 2;
			}
		}
	}

	for(int a = 0;a<2;a++){//for each mixed boundary face
		double phiint = Tb[a];
		for(int b = 0;b<faceid_qb[0].size();b++){
			int faceid = faceid_qb[0][b];
			double h = wtd_h[a](b);
			//Log << "Face: " << faceid << "\t" ;
			double dA = FaceArea(connf(faceid,0),connf(faceid,1),connf(faceid,2));
			//Log << "Area: " << dA << "\tPhiinfinite: " << phiint << "\th: " << h << "\n";;
			for(int c = 0;c<3;c++){
				vecq[connf(faceid,c)] = vecq[connf(faceid,c)] + h*phiint*dA/3;
				boundaries[connf(faceid,c)] = 2;
			}
		}
	}

	area_piston.resize(faceid_qb.size());

	for(int a = 0;a<faceid_qb.size();a++){//for each flux boundary face
		area_piston[a].resize(faceid_qb[a].size());
		for(int b = 0;b<faceid_qb[a].size();b++){
			int faceid = faceid_qb[a][b];
			double q = qb(faceid);
			double dA = FaceArea(connf(faceid,0),connf(faceid,1),connf(faceid,2));
			area_piston[a][b] = dA;
			//Log << scientific << "Face: " << faceid << "\tArea: " << dA << "\tq: " << q << "\n";
			for(int c = 0;c<3;c++){
				vecq[connf(faceid,c)] = vecq[connf(faceid,c)] + q*dA/3;
				boundaries[connf(faceid,c)] = 1;
			}
		} 
	}
	
	Log << "Done!" << "\n";

	//Set Constraints (phib)
	Log << "Setting Temperature Constraints... " ;

	for(int a = 0;a<faceid_phib.size();a++){//for each fixed temperature face
		for(int b = 0;b<faceid_phib[a].size();b++){
			int faceid = faceid_phib[a][b];
			double phiint = phib(faceid);
			for(int c = 0;c<3;c++){
				int node = connf(faceid,c);
				//phi[node] = phiint; // Set Temperature
				vecq[node] = phiint;
				T_body(node) = phiint;//eliminate under-relaxation for this node.
				for(int i = 0;i<nNodes;i++){
					//matK(i,node) = 0;
					if(i == node)
						matK(node,i) = 1;
					else
						matK(node,i) = 0;
				}
			}
		}
	}

	connf_piston.resize(connf.extent(0),connf.extent(1));
	connf_piston = connf;
	faceid_qb_piston = faceid_qb;


	connallf.free();
	Log << "Done!" << "\n";

	connf.free(); Afn.free(); phi.clear(); phi_surf.free(); faceid_qb.clear(); faceid_mixb.clear();

};
// add conduction matrix for cylinder body
void CThermal::ThermalSolve_Cylinder(vector<Array<double,1>> qbi, vector<Array<double,1>> wtd_h, vector<double> Tb, Array<double,1> &T_body,Array<double,1> &T_surface)
{
	int n_piston_Nodes;
	n_piston_Nodes = matK.nrows();


	string mesh_name;
	mesh_name = myMesh.cylinder_name;
	
	//Size and initialize scalar field to be solved
	gmm::resize(phi_both,n_piston_Nodes+nNodes);
	//phi.resize(nNodes,0);
	
	//Log << nCells << "\n";
	//Log << myMesh.conn.size() << "\n";
	for(int a = 0; a < nNodes; a++){
		double temp = T_body(a);
		//phi[a] = temp;
	//	Log << temp << "\n";
		for(int b = 0; b < 4; b++){
			phi_both[myMesh.conn(a,b) + n_piston_Nodes] = temp;
		}
	}

	//vector<short> boundaries;
	boundaries.resize(nNodes,0);
	//Set boundaries
	SetBoundaries(qbi);
	

	//Develop Stiffness Matrix
	Log << "Calculating Stiffness Matrix... ";
	gmm:resize(matK,nNodes+n_piston_Nodes,nNodes+n_piston_Nodes);
	//gmm::row_matrix<gmm::wsvector<double > > matK(nNodes,nNodes);
	//gmm::clear(matK);
	for(int cell = 0;cell<nCells;cell++){
		/*Log << "Cell: " << cell << "\n";
		Log << "\n";*/
		double J = SetBMatrix(cell);
		/*Log << "J: " << J << "\n";
		Log << "\n";*/
		vector<vector<double> > Kcd;
		Kcd = CalcKcd(matB,myMesh.Tf(cell),J);//calculate the conduction matrix
		for(int a = 0;a<4;a++){
			for(int b = 0;b<4;b++){
				//if(myMesh.conn(cell,a)>=myMesh.conn(cell,b))
					matK(myMesh.conn(cell,a)+n_piston_Nodes,myMesh.conn(cell,b)+n_piston_Nodes) = matK(myMesh.conn(cell,a)+n_piston_Nodes,myMesh.conn(cell,b)+n_piston_Nodes) + Kcd[a][b];
			}
		}
		//Log << "\n";
	}


	//add convective terms for mixed boundary faces.
	for(int a = 0;a<faceid_mixb.size();a++){//for each mixed boundary face
		//Log << "Applying mixed boundary face " << a << "\n";
		//double phi = myMesh.phibinfi[a];
		double h = myMesh.hbi[a];
		for(int b = 0;b<faceid_mixb[a].size();b++){
			int faceid = faceid_mixb[a][b];
			//cout << "Face: " << faceid << "\n";
			//cout << Af.size() << "\n";
			double dA = FaceArea(connf(faceid,0),connf(faceid,1),connf(faceid,2));
			//cout << "Area: " << dA << "\n";
			//cout << connallf.size() << "\n";
			for(int c = 0;c<3;c++){
				int nodeid = connf(faceid,c) + n_piston_Nodes;
				//cout << "Node: " << nodeid << "\n";
				matK(nodeid,nodeid) = matK(nodeid,nodeid) + h*dA/6;
				if(c != 0){
					//if(connallf(faceid,0) <= nodeid)
						matK(nodeid,connf(faceid,0)+ n_piston_Nodes) = matK(nodeid,connf(faceid,0)+ n_piston_Nodes) + h*dA/12;
				}
				if(c != 1){
					//if(connallf(faceid,1) <= nodeid)
						matK(nodeid,connf(faceid,1)+ n_piston_Nodes) = matK(nodeid,connf(faceid,1)+ n_piston_Nodes) + h*dA/12;
				}
				if(c != 2){
					//if(connallf(faceid,2) <= nodeid)
						matK(nodeid,connf(faceid,2)+ n_piston_Nodes) = matK(nodeid,connf(faceid,2)+ n_piston_Nodes) + h*dA/12;
				}
			}
		}
	}

	//add convective terms for mixed heat transfer on flux boundary faces.
	for(int a = 0;a<2;a++){//for each mixed boundary face
		//Log << "Applying mixed boundary face " << a << "\n";
		//double phi = myMesh.phibinfi[a];
		for(int b = 0;b<faceid_qb[0].size();b++){
			int faceid = faceid_qb[0][b];
			double h = wtd_h[a](b);
			//cout << "Face: " << faceid << "\n";
			//cout << Af.size() << "\n";
			double dA = FaceArea(connf(faceid,0),connf(faceid,1),connf(faceid,2));
			//cout << "Area: " << dA << "\n";
			//cout << connallf.size() << "\n";
			for(int c = 0;c<3;c++){
				int nodeid = connf(faceid,c);
				//cout << "Node: " << nodeid << "\n";
				matK(nodeid,nodeid) = matK(nodeid,nodeid) + h*dA/6;
				if(c != 0){
					//if(connallf(faceid,0) <= nodeid)
						matK(nodeid,connf(faceid,0)) = matK(nodeid,connf(faceid,0)) + h*dA/12;
				}
				if(c != 1){
					//if(connallf(faceid,1) <= nodeid)
						matK(nodeid,connf(faceid,1)) = matK(nodeid,connf(faceid,1)) + h*dA/12;
				}
				if(c != 2){
					//if(connallf(faceid,2) <= nodeid)
						matK(nodeid,connf(faceid,2)) = matK(nodeid,connf(faceid,2)) + h*dA/12;
				}
			}
		}
	}

	Log << "Done!" << "\n";

	//Develop Load Matrix
	Log << "Setting Load Vector... ";
	//vecq.clear();
	vecq.resize(nNodes+n_piston_Nodes,0);
	
	
	for(int a = 0;a<faceid_mixb.size();a++){//for each mixed boundary face
		double phiint = myMesh.phibinfi[a];
		double h = myMesh.hbi[a];
		for(int b = 0;b<faceid_mixb[a].size();b++){
			int faceid = faceid_mixb[a][b];
			//Log << "Face: " << faceid << "\t" ;
			double dA = FaceArea(connf(faceid,0),connf(faceid,1),connf(faceid,2));
			//Log << "Area: " << dA << "\tPhiinfinite: " << phiint << "\th: " << h << "\n";;
			for(int c = 0;c<3;c++){
				vecq[connf(faceid,c)+n_piston_Nodes] = vecq[connf(faceid,c)+n_piston_Nodes] + h*phiint*dA/3;
				boundaries[connf(faceid,c)] = 2;
			}
		}
	}

	for(int a = 0;a<2;a++){//for each mixed boundary face
		double phiint = Tb[a];
		for(int b = 0;b<faceid_qb[0].size();b++){
			int faceid = faceid_qb[0][b];
			double h = wtd_h[a](b);
			//Log << "Face: " << faceid << "\t" ;
			double dA = FaceArea(connf(faceid,0),connf(faceid,1),connf(faceid,2));
			//Log << "Area: " << dA << "\tPhiinfinite: " << phiint << "\th: " << h << "\n";;
			for(int c = 0;c<3;c++){
				vecq[connf(faceid,c)] = vecq[connf(faceid,c)] + h*phiint*dA/3;
				boundaries[connf(faceid,c)] = 2;
			}
		}
	}

	area_cylinder.resize(faceid_qb.size());

	for(int a = 0;a<faceid_qb.size();a++){//for each flux boundary face
		area_cylinder[a].resize(faceid_qb[a].size());
		for(int b = 0;b<faceid_qb[a].size();b++){
			int faceid = faceid_qb[a][b];
			double q = qb(faceid);
			double dA = FaceArea(connf(faceid,0),connf(faceid,1),connf(faceid,2));
			area_cylinder[a][b] = dA;
			//Log << scientific << "Face: " << faceid << "\tArea: " << dA << "\tq: " << q << "\n";
			for(int c = 0;c<3;c++){
				vecq[connf(faceid,c)+n_piston_Nodes] = vecq[connf(faceid,c)+n_piston_Nodes] + q*dA/3;
				boundaries[connf(faceid,c)] = 1;
			}
		}
	}
	
	Log << "Done!" << "\n";

	//Set Constraints (phib)
	Log << "Setting Temperature Constraints... " ;

	for(int a = 0;a<faceid_phib.size();a++){//for each fixed temperature face
		for(int b = 0;b<faceid_phib[a].size();b++){
			int faceid = faceid_phib[a][b];
			double phiint = phib(faceid);
			for(int c = 0;c<3;c++){
				int node = connf(faceid,c);
				//phi[node] = phiint; // Set Temperature
				vecq[node+n_piston_Nodes] = phiint;
				T_body(node+n_piston_Nodes) = phiint;//eliminate under-relaxation for this node.
				for(int i = n_piston_Nodes;i<nNodes+n_piston_Nodes;i++){
					//matK(i,node) = 0;
					if(i == node+n_piston_Nodes)
						matK(node+n_piston_Nodes,i) = 1;
					else
						matK(node+n_piston_Nodes,i) = 0;
				}
			}
		}
	}

	
	connf_cylinder.resize(connf.extent(0),connf.extent(1));
	connf_cylinder = connf;
	faceid_qb_cylinder = faceid_qb;

	connallf.free();
	Log << "Done!" << "\n";

	connf.free(); Afn.free(); phi.clear(); phi_surf.free(); faceid_qb.clear(); faceid_mixb.clear();

};
// add conduction matrix between piston and cylinder body
void CThermal::ThermalSolve_interbody_conduction(vector<vector<int>> &Mid_K2B, vector<vector<double>> &Mwt_K2B, vector<vector<int>> &Mid_B2K, vector<vector<double>> &Mwt_B2K, double C)
{
	//double lambda_oil = my_oil -> get_lambda();//oilpistongap.oillambda;
	//double C = operatingpistongap.speed * lambda_oil / 60 / 9;
	// from piston to cylinder
	for(int j = 0; j < Mid_K2B.size(); j++)
	{
		int faceid = faceid_qb_cylinder[0][j];
		for(int i = 0; i < Mid_K2B[j].size(); i++)
		{
			int faceid_p = faceid_qb_piston[0][Mid_K2B[j][i]];
			double dA = area_cylinder[0][j];
			double kj = +1 * C * dA * Mwt_K2B[j][i];
			double ki = -1 * C * dA * Mwt_K2B[j][i];
			matK(connf_cylinder(faceid,0),connf_cylinder(faceid,0)) += kj;
			matK(connf_cylinder(faceid,0),connf_cylinder(faceid,1)) += kj;
			matK(connf_cylinder(faceid,0),connf_cylinder(faceid,2)) += kj;
			matK(connf_cylinder(faceid,1),connf_cylinder(faceid,0)) += kj;
			matK(connf_cylinder(faceid,1),connf_cylinder(faceid,1)) += kj;
			matK(connf_cylinder(faceid,1),connf_cylinder(faceid,2)) += kj;
			matK(connf_cylinder(faceid,2),connf_cylinder(faceid,0)) += kj;
			matK(connf_cylinder(faceid,2),connf_cylinder(faceid,1)) += kj;
			matK(connf_cylinder(faceid,2),connf_cylinder(faceid,2)) += kj;
			matK(connf_cylinder(faceid,0),connf_piston(faceid_p,0)) += ki;
			matK(connf_cylinder(faceid,0),connf_piston(faceid_p,1)) += ki;
			matK(connf_cylinder(faceid,0),connf_piston(faceid_p,2)) += ki;
			matK(connf_cylinder(faceid,1),connf_piston(faceid_p,0)) += ki;
			matK(connf_cylinder(faceid,1),connf_piston(faceid_p,1)) += ki;
			matK(connf_cylinder(faceid,1),connf_piston(faceid_p,2)) += ki;
			matK(connf_cylinder(faceid,2),connf_piston(faceid_p,0)) += ki;
			matK(connf_cylinder(faceid,2),connf_piston(faceid_p,1)) += ki;
			matK(connf_cylinder(faceid,2),connf_piston(faceid_p,2)) += ki;
		}
	}
	// from cylinder to piston
	for(int i = 0; i < Mid_B2K.size(); i++)
	{
		int faceid_p = faceid_qb_piston[0][i];
		for(int j = 0; j < Mid_B2K[i].size(); j++)
		{
			int faceid = faceid_qb_cylinder[0][Mid_B2K[i][j]];
			double dA = area_piston[0][i];
			double ki = +1 * C * dA * Mwt_B2K[i][j];
			double kj = -1 * C * dA * Mwt_B2K[i][j];
			matK(connf_piston(faceid_p,0),connf_piston(faceid_p,0)) += ki;
			matK(connf_piston(faceid_p,0),connf_piston(faceid_p,1)) += ki;
			matK(connf_piston(faceid_p,0),connf_piston(faceid_p,2)) += ki;
			matK(connf_piston(faceid_p,1),connf_piston(faceid_p,0)) += ki;
			matK(connf_piston(faceid_p,1),connf_piston(faceid_p,1)) += ki;
			matK(connf_piston(faceid_p,1),connf_piston(faceid_p,2)) += ki;
			matK(connf_piston(faceid_p,2),connf_piston(faceid_p,0)) += ki;
			matK(connf_piston(faceid_p,2),connf_piston(faceid_p,1)) += ki;
			matK(connf_piston(faceid_p,2),connf_piston(faceid_p,2)) += ki;
			matK(connf_piston(faceid_p,0),connf_cylinder(faceid,0)) += kj;
			matK(connf_piston(faceid_p,0),connf_cylinder(faceid,1)) += kj;
			matK(connf_piston(faceid_p,0),connf_cylinder(faceid,2)) += kj;
			matK(connf_piston(faceid_p,1),connf_cylinder(faceid,0)) += kj;
			matK(connf_piston(faceid_p,1),connf_cylinder(faceid,1)) += kj;
			matK(connf_piston(faceid_p,1),connf_cylinder(faceid,2)) += kj;
			matK(connf_piston(faceid_p,2),connf_cylinder(faceid,0)) += kj;
			matK(connf_piston(faceid_p,2),connf_cylinder(faceid,1)) += kj;
			matK(connf_piston(faceid_p,2),connf_cylinder(faceid,2)) += kj;

		}
	}
	connf_piston.free(); connf_cylinder.free(); area_piston.clear(); area_cylinder.clear(); faceid_qb_piston.clear(); faceid_qb_cylinder.clear();
}
// solving heat transfer linear system
void CThermal::ThermalSolve_new(void)
{

	Log << "\n";
	Log << "Solving both body temperature distribution... " << "\n";
	Log << "\n";

	gmm::csr_matrix<double> Kint;
	//gmm::clean(Kint, 1E-12);
	gmm::copy(matK,Kint);
	gmm::clear(matK);
	
	//Solve System
	double tol;
	tol = myinput.data.thermal.piston.solver.tol;

	gmm::iteration iter(tol);
	iter.set_maxiter(myinput.data.thermal.piston.solver.maxit);


	if(coupled_caspar_simulation)
		iter.set_noisy(0);
	else
		iter.set_noisy(1);

	gmm::identity_matrix BS;
	gmm::diagonal_precond<gmm::csr_matrix<double > >P(Kint);
	gmm::cg(Kint, phi_both, vecq, BS, P, iter);

	if(iter.converged())
		Log << "Temperature Distribution Complete!" << "\n";
	else{
		Log << "Temperature Distribution Failed!" << endl;
		exit(1);
	}
};
// assign heat transfer linear system result to piston/cylinder body/surface temperature array
void CThermal::ThermalSolve_post(string body,vector<Array<double,1>> qbi,Array<double,1> &T_body,Array<double,1> &T_surface)
{

	Log << "\n";
	Log << "Solving " + body + " temperature distribution... " << "\n";
	Log << "\n";

	string mesh_name;
	phi.resize(nNodes,0);
	//Define paths according to solid body
	if(body=="piston")
	{
		mesh_name = myMesh.piston_name;
		for(int i=0;i<nNodes;i++)
			phi[i] = phi_both[i];
	}
	else if(body=="cylinder")
	{
		mesh_name = myMesh.cylinder_name;
		for(int i=0;i<nNodes;i++)
			phi[i] = phi_both[phi_both.size()-nNodes+i];
	}
		

	//Assign results
	T_body.resize(nCells);
	Array<double,1> T_body_old = T_body;
	for(int i = 0;i<nCells;i++)
		T_body(i) = 0.25*(phi[myMesh.conn(i,0)] + phi[myMesh.conn(i,1)] + phi[myMesh.conn(i,2)] + phi[myMesh.conn(i,3)]);
	T_body = T_body_old + myinput.data.options_piston.numeric.AlphaTh * (T_body - T_body_old);
	T_body_old.free();
	
	T_surface.resize(faceid_qb[0].size());
	Array<double,1> T_surf_old = T_surface;
	for(int i = 0;i<faceid_qb[0].size();i++){
		//double temp = phi[connf(i,0)]/3 + phi[connf(i,1)]/3 + phi[connf(i,2)]/3;
		int id = faceid_qb[0][i];
		T_surface(i) = (phi[connf(id,0)] + phi[connf(id,1)] + phi[connf(id,2)]) / 3.0;
		//Log << "T_surface(" << i << ") = " << temp << "\n";
	}
	T_surface = T_surf_old + myinput.data.options_piston.numeric.AlphaTh * (T_surface - T_surf_old);
	T_surf_old.free();


	//Output VTK and others
	ThermalOutput(body);

	

	//free arrays
	connf.free(); Afn.free(); phi.clear(); phi_surf.free(); faceid_qb.clear(); faceid_mixb.clear();


	Log << "\n";
	Log << "Done " + body + " temperature distribution!" << "\n";
	Log << "\n";


};

void CThermal::CheckBoundaryFaces(string body)
{
	
	//Count boundary faces from solid mesh
	int nf_bd_mesh=0;
	for(int i=0;i<nFaces;i++)
	{
		//Log << Face2Cell(i,0) << "\t" << Face2Cell(i,1) << "\n";
		for(int j=0;j<2;j++)
		{
			
			if(Face2Cell(i,j)==-1)
			{
				nf_bd_mesh++;
				//face = i
				int found;
				for(int a = 0;a<myMesh.mixbABQS.size();a++){//index of the boundary set
					found = 0;
					for(int b = 0;b<3;b++){//search for the three faces
						for(int c = 0;c<myMesh.mixbABQS[a].size();c++){
							if(myMesh.mixbABQS[a][c]==connf(i,b)){
								found++;
								break;
							}
						}
					}
					if(found==3)
						break;
				}
				if(found!=3){
					
					for(int a = 0;a<myMesh.phibABQS.size();a++){//index of the boundary set
						found = 0;
						for(int b = 0;b<3;b++){//search for the three faces
							for(int c = 0;c<myMesh.phibABQS[a].size();c++){
								if(myMesh.phibABQS[a][c]==connf(i,b)){
									found++;
									break;
								}
							}
						}
						if(found==3)
							break;
					}
				}
				if(found!=3){
					
					for(int a = 0;a<myMesh.qbABQS.size();a++){//index of the boundary set
						found = 0;
						for(int b = 0;b<3;b++){//search for the three faces
							for(int c = 0;c<myMesh.qbABQS[a].size();c++){
								if(myMesh.qbABQS[a][c]==connf(i,b)){
									found++;
									break;
								}
							}
						}
						if(found==3)
							break;
					}
				}
				if(found!=3){
					Log << "Error: Face " << i << " not defined in thermal boundary conditions!" << "\n";
					Log << "Nodes: " << connf(i,0) << ", " << connf(i,1) << ", " << connf(i,2) << endl;
					exit(1);
				}
				else
					myMesh.nf_bd_user++;
			} // for each surface node you know 3 node id's from connf.  Search for these in the user defined node sets.  If they are all there, then the face is defined.  If not, send an error.
		}		// Marco recommends using an ANN kd tree.
		
	}

	//free variables
	//connallf.free(); moved dwm
	

	//Stop if number of boundary faces is not the same as user defined
	
	if(nf_bd_mesh!=(myMesh.nf_bd_user / 2))
	{
		int temp = myMesh.nf_bd_user / 2;
		Log << "\n";
		Log << "Warning: Some of " + body + " boundary faces may be missing!" << "\n";
		Log << "Solid mesh boundary faces are: " << nf_bd_mesh << "\n";
		Log << "User defined boundary faces are: " << temp << "\n";
		Log << "\n";
	}
};
void CThermal::BoundaryFaces(void)
{

	//Assign value to boundary faces 
	FaceBoundId.resize(nFaces);
	FaceBoundId = 0;

	//Dirichlet
	Facephib();

	//Mixed
	Facemixb();

	//Neumann
	Faceqb();

	//clear mesh inputs
	myMesh.qbABQS.clear();	myMesh.phibABQS.clear();	myMesh.mixbABQS.clear();

};
void CThermal::Facephib(void)
{

	//-----------phi value at node sets from boundaries Dirichlet----------//
	Array<double,2> phib0(nNodes,(int) myMesh.phibABQS.size());
	phib0 = 0.0;
	if((int) myMesh.phibABQS.size())
	{
		for(int i=0;i<(int) myMesh.phibABQS.size();i++)
		{
			for(int j=0;j<(int) myMesh.phibABQS[i].size();j++)
			{
				phib0(myMesh.phibABQS[i][j],i) = 1.0;
			}
		};
	};

	//Identify boundary faces
	int nphib;
	nphib = (int) phib0.extent(1);
	faceid_phib.resize(nphib);
	for(int j=0;j<nphib;j++)
	{
		for(int i=0;i<nFaces;i++)
		{
			if(Face2Cell(i,0)==-1 || Face2Cell(i,1)==-1)
			{
				if(Teth)
				{
					if( phib0(connf(i,0),j)!=0.0 && phib0(connf(i,1),j)!=0.0 && phib0(connf(i,2),j)!=0.0 )
					{
						//identify face boundary
						FaceBoundId(i) = 1;
						//face id dirichlet
						faceid_phib[j].push_back(i);
					}
				};
				if(Hexa)
				{
					if( phib0(connf(i,0),j)!=0.0 && phib0(connf(i,1),j)!=0.0 && phib0(connf(i,2),j)!=0.0 && phib0(connf(i,3),j)!=0.0 )
					{
						//identify face boundary
						FaceBoundId(i) = 1;
						//face id dirichlet
						faceid_phib[j].push_back(i);
					}
				};
			};
		};
	};

};
void CThermal::Faceqb(void)
{
	//-----------flux value at node sets from boundaries Neumann-----------//
	Array<double,2> qb0(nNodes,(int) myMesh.qbABQS.size());
	qb0 = 0.0;
	nodeid_qb.resize((int) myMesh.qbABQS.size());	
	for(int i=0;i<(int) myMesh.qbABQS.size();i++)
	{
		for(int j=0;j<(int) myMesh.qbABQS[i].size();j++)
		{
			qb0(myMesh.qbABQS[i][j],i) = 1.0;
			//nodes coordinates gap surface
			nodeid_qb[i].push_back( myMesh.qbABQS[i][j] );
		}
	};

	unsigned short tempint = myMesh.qbABQS.size();
	Log << "Applying Heat Flux Boundaries to " << tempint << " Surfaces." << "\n";

	//Identify boundary faces
	int nqb;
	nqb = (int) qb0.extent(1);
	vector<int> free; free.resize(0,0);
	faceid_qb.resize(nqb,free);
	for(int j=0;j<nqb;j++)
	{
		for(int i=0;i<nFaces;i++)
		{
			if( Face2Cell(i,0)==-1 || Face2Cell(i,1)==-1 )
			{
				if(Teth)
				{
					if( qb0(connf(i,0),j)!=0.0 && qb0(connf(i,1),j)!=0.0 && qb0(connf(i,2),j)!=0.0 )
					{
						//identify face boundary
						FaceBoundId(i) = 2;
						//face id neumann
						faceid_qb[j].push_back(i);
					}
				};
				if(Hexa)
				{
					if( qb0(connf(i,0),j)!=0.0 && qb0(connf(i,1),j)!=0.0 && qb0(connf(i,2),j)!=0.0 && qb0(connf(i,3),j)!=0.0 )
					{
						//identify face boundary
						FaceBoundId(i) = 2;
						//face id neumann
						faceid_qb[j].push_back(i);
					}
				};
			};
		};
	}

};
void CThermal::Facemixb(void)
{

	//-----------mixb value at node sets from boundaries Mixed-----------//
	Array<double,2> mixb0(nNodes,(int) myMesh.mixbABQS.size());
	mixb0 = 0.0;
	if((int) myMesh.mixbABQS.size())
	{
		for(int i=0;i<(int) myMesh.mixbABQS.size();i++)
		{
			for(int j=0;j<(int) myMesh.mixbABQS[i].size();j++)
			{
				mixb0(myMesh.mixbABQS[i][j],i) = 1.0;
			}
		};
	};

	//Identify boundary faces
	int nmixb;
	nmixb = (int) mixb0.extent(1);
	faceid_mixb.resize(nmixb);
	for(int j=0;j<nmixb;j++)
	{
		for(int i=0;i<nFaces;i++)
		{
			if(Face2Cell(i,0)==-1 || Face2Cell(i,1)==-1)
			{
				if(Teth)
				{
					if( mixb0(connf(i,0),j)!=0.0 && mixb0(connf(i,1),j)!=0.0 && mixb0(connf(i,2),j)!=0.0 )
					{
						//identify face boundary
						FaceBoundId(i) = 3;
						//face id mixed
						faceid_mixb[j].push_back(i);
					}
				};
				if(Hexa)
				{
					if( mixb0(connf(i,0),j)!=0.0 && mixb0(connf(i,1),j)!=0.0 && mixb0(connf(i,2),j)!=0.0 && mixb0(connf(i,3),j)!=0.0 )
					{
						//identify face boundary
						FaceBoundId(i) = 3;
						//face id mixed
						faceid_mixb[j].push_back(i);
					}
				};
			};
		};
	};

};
void CThermal::WriteGapSurfacexyz(int body)
{

	//Write Gap Surface Face Centers Coordinates
	int nFaces_gap;
	vector<int> Faces_gap;
	nFaces_gap = (int) faceid_qb[0].size();
	Faces_gap = faceid_qb[0];
	string output;
	if(body == 1)
		output = "./temp/piston/piston/xyz_faces.dat";//mesh_path + "/" + "xyz_faces.dat";
	else if(body == 0)
		output = "./temp/piston/bushing/xyz_faces.dat";
	fout.open(output.c_str());
	fout << nFaces_gap << "\n";
	for(int i=0;i<nFaces_gap;i++)
	{
		
		fout << scientific << xyzf(Faces_gap[i],0) << "\t" << 
			xyzf(Faces_gap[i],1) << "\t" << xyzf(Faces_gap[i],2) << "\n";
	}
	fout << "*" << "\n";
	fout.close();
	fout.clear();

	//white Gap Surface Face Centers Coordinate for other cylinder bores at ODC
	/*double dtheta = 2.0*PI/(double) faceid_qb.size();
	double theta = dtheta;
	if(body == 1)
	{
		for(int j=1;j<faceid_qb.size();j++)
		{
			int nFaces_gap;
			vector<int> Faces_gap;
			nFaces_gap = (int) faceid_qb[j].size();
			Faces_gap = faceid_qb[j];
			string output;
			string bore_n;
			char buffer[500];
			_itoa_s(j,buffer,10);
			bore_n = buffer;
			output = "./temp/piston/bushing/xyz_faces" + bore_n + ".dat";
			fout.open(output.c_str());
			fout << nFaces_gap << "\n";
			for(int i=0;i<nFaces_gap;i++)
			{				
				fout << scientific << xyzf(Faces_gap[i],0) * cos(theta) - xyzf(Faces_gap[i],1) * sin(theta) << "\t" << 
					xyzf(Faces_gap[i],0) * sin(theta) + xyzf(Faces_gap[i],1) * cos(theta)<< "\t" << xyzf(Faces_gap[i],2) << "\n";
			}
			fout << "*" << "\n";
			fout.close();
			fout.clear();
			theta += dtheta;
		}
	}*/


	//Write Gap Surface Nodes Coordinates
	int nNodes_gap;
	vector<int> Nodes_gap;
	//number of nodes in surface
	nNodes_gap = (int) nodeid_qb[0].size();
	//nodes ids in surface
	Nodes_gap = nodeid_qb[0];
	//write file
	if(body == 1)
		output = "./temp/piston/piston/xyz_nodes.dat";
	else if(body == 0)
		output = "./temp/piston/bushing/xyz_nodes.dat";
	//output = mesh_path + "/" + "xyz_nodes.dat";
	fout.open(output.c_str());
	fout << nNodes_gap << "\n";
	for(int i=0;i<nNodes_gap;i++)
	{
		fout << scientific << myMesh.xyz(Nodes_gap[i],0) << "\t" << 
			myMesh.xyz(Nodes_gap[i],1) << "\t" << myMesh.xyz(Nodes_gap[i],2) << "\n";
	}
	fout << "*" << "\n";
	fout.close();
	fout.clear();

}

void CThermal::SetBoundaries(vector<Array<double,1>> qbi)
{
	//Boundary faces values Dirichlet
	phib.resize(nFaces);	phib = 0.0;
	for(int i=0;i<(int) faceid_phib.size();i++)
	{
		for(int j=0;j<(int) faceid_phib[i].size();j++)
		{
			int id = faceid_phib[i][j];
			phib(id) = myMesh.phibi[i];
		}
	}


	//Boundary faces values Mixed
	phibinf.resize(nFaces); phibinf=0.0;
	hb.resize(nFaces);		hb=0.0;
	for(int i=0;i<(int) faceid_mixb.size();i++)
	{
		for(int j=0;j<(int) faceid_mixb[i].size();j++)
		{
			int id = faceid_mixb[i][j];
			phibinf(id) = myMesh.phibinfi[i];
			hb(id) = myMesh.hbi[i];
		}
	}


	//Boundary faces values Neumann
	qb.resize(nFaces);	qb=0.0;
	for(int i=0;i<(int) faceid_qb.size();i++)
	{
		//cout<<i<<endl;
		for(int j=0;j<(int) faceid_qb[i].size();j++)
		{
			int id = faceid_qb[i][j];
			qb(id) = qbi[i](j);
			//Log << "Face " << id << " loaded with " << qb(id) << " from boundary " << i << "\n";
		}
	}

}
void CThermal::CellCenters(void)
{
	Range all = Range::all();

	//Cell centers coordinates
	xyzc.resize(nCells,3);
	xyzc=0.0;

	//Cells Centroids coordinate for teth mesh
	if(Teth)
	{
		for(int i = 0; i < nCells; i++)
		{
			for(int j=0;j<4;j++)
			{
				xyzc(i,0) +=  myMesh.xyz(myMesh.conn(i,j),0) ;
				xyzc(i,1) +=  myMesh.xyz(myMesh.conn(i,j),1) ;
				xyzc(i,2) +=  myMesh.xyz(myMesh.conn(i,j),2) ;
			};
			xyzc(i,all) *= 0.25;
		};
		
	}

	//Cells Centroids coordinate for hexa mesh
	if(Hexa)
	{
		for(int i = 0; i < nCells; i++)
		{
			for(int j=0;j<8;j++)
			{
				xyzc(i,0) +=  myMesh.xyz(myMesh.conn(i,j),0);
				xyzc(i,1) +=  myMesh.xyz(myMesh.conn(i,j),1);
				xyzc(i,2) +=  myMesh.xyz(myMesh.conn(i,j),2);
			};
			xyzc(i,all) *= 0.125;
		};
	}

	//Calculate cell volumes
	if( (myMesh.IR_p==1 && myMesh.body==1) || (myMesh.IR_c==1 && myMesh.body==0) )
	{
		CellMassVolume();
	};

};
void CThermal::CellMassVolume(void)
{

	//Volume
	Vol.resize(nCells);
	Vol=0.0;
	Array<double,1> cross(3);
	Array<double,1> Vol3(nCells);
	Vol3=0.0;
	cross=0.0;
	if(Teth)
	{
		for(int i=0;i<nCells;i++)
		{
			//x component vector 1 and 2
			double x1 =  myMesh.xyz(myMesh.conn(i,1),0) -  myMesh.xyz(myMesh.conn(i,0),0);
			double x2 =  myMesh.xyz(myMesh.conn(i,2),0) -  myMesh.xyz(myMesh.conn(i,0),0);
			double x3 =  myMesh.xyz(myMesh.conn(i,3),0) -  myMesh.xyz(myMesh.conn(i,0),0);
			//y component vector 1 and 2
			double y1 =  myMesh.xyz(myMesh.conn(i,1),1) -  myMesh.xyz(myMesh.conn(i,0),1);
			double y2 =  myMesh.xyz(myMesh.conn(i,2),1) -  myMesh.xyz(myMesh.conn(i,0),1);
			double y3 =  myMesh.xyz(myMesh.conn(i,3),1) -  myMesh.xyz(myMesh.conn(i,0),1);
			//z component vector 1 and 2
			double z1 =  myMesh.xyz(myMesh.conn(i,1),2) -  myMesh.xyz(myMesh.conn(i,0),2);
			double z2 =  myMesh.xyz(myMesh.conn(i,2),2) -  myMesh.xyz(myMesh.conn(i,0),2);
			double z3 =  myMesh.xyz(myMesh.conn(i,3),2) -  myMesh.xyz(myMesh.conn(i,0),2);
			//Cross product
			cross(0) = y1*z2 - z1*y2;
			cross(1) = z1*x2 - x1*z2;
			cross(2) = x1*y2 - y1*x2;
			//Volume
			Vol(i) = 1.0/6.0 * fabs( x3 * cross(0) + y3 * cross(1) + z3 * cross(2) );
			Vol3(i) = pow(Vol(i),1.0/3.0);
		}
	}
	else
	{
		//5 teth in hexa edge combinations
		Array<int,2> perm(5,4);
		perm(0,0)=1; perm(0,1)=2; perm(0,2)=5; perm(0,3)=0;
		perm(1,0)=0; perm(1,1)=2; perm(1,2)=7; perm(1,3)=3;
		perm(2,0)=2; perm(2,1)=7; perm(2,2)=5; perm(2,3)=6;
		perm(3,0)=5; perm(3,1)=7; perm(3,2)=0; perm(3,3)=4;
		perm(4,0)=2; perm(4,1)=5; perm(4,2)=7; perm(4,3)=0;
		for(int i=0;i<nCells;i++)
		{
			//5 teth in a hexa
			for(int j=0;j<5;j++)
			{
				//x component vector 1 and 2
				double x1 =  myMesh.xyz(myMesh.conn(i,perm(j,0)),0) -  myMesh.xyz(myMesh.conn(i,0),perm(j,3));
				double x2 =  myMesh.xyz(myMesh.conn(i,perm(j,1)),0) -  myMesh.xyz(myMesh.conn(i,0),perm(j,3));
				double x3 =  myMesh.xyz(myMesh.conn(i,perm(j,2)),0) -  myMesh.xyz(myMesh.conn(i,0),perm(j,3));
				//y component vector 1 and 2
				double y1 =  myMesh.xyz(myMesh.conn(i,perm(j,0)),1) -  myMesh.xyz(myMesh.conn(i,0),perm(j,3));
				double y2 =  myMesh.xyz(myMesh.conn(i,perm(j,0)),1) -  myMesh.xyz(myMesh.conn(i,0),perm(j,3));
				double y3 =  myMesh.xyz(myMesh.conn(i,perm(j,0)),1) -  myMesh.xyz(myMesh.conn(i,0),perm(j,3));
				//z component vector 1 and 2
				double z1 =  myMesh.xyz(myMesh.conn(i,perm(j,0)),2) -  myMesh.xyz(myMesh.conn(i,0),perm(j,3));
				double z2 =  myMesh.xyz(myMesh.conn(i,perm(j,0)),2) -  myMesh.xyz(myMesh.conn(i,0),perm(j,3));
				double z3 =  myMesh.xyz(myMesh.conn(i,perm(j,0)),2) -  myMesh.xyz(myMesh.conn(i,0),perm(j,3));
				//Cross product
				cross(0) = y1*z2 - z1*y2;
				cross(1) = z1*x2 - x1*z2;
				cross(2) = x1*y2 - y1*x2;
				//Volume
				Vol(i) += 1.0/6.0 * fabs( x3 * cross(0) + y3 * cross(1) + z3 * cross(2) );
			}
			Vol3(i) = pow(Vol(i),1.0/3.0);
		}
	};


	//Cell Mass
	Range all=Range::all();

	M.resize(nCells);
	M = 0.0;
	
	//Steel density
	rho.resize(nCells);
	rho = myMesh.rho;

	//Cells mass
	M = rho * Vol;

	//Total mass
	Mtot = sum(M);

	//Log << "Total Mass: " << Mtot << "\n";

	//Center of gravity coordinates
	xcg = sum(xyzc(all,0)*M)/Mtot;
	ycg = sum(xyzc(all,1)*M)/Mtot;
	zcg = sum(xyzc(all,2)*M)/Mtot;

	//Moment of inertia repect to CG
	I.resize(3,3);	I=0.0;
	for(int i=0;i<nCells;i++)
	{
		//Principal
		I(0,0) += M(i) * ( pow(xyzc(i,1)-ycg,2.0) + pow(xyzc(i,2)-zcg,2.0) );
		I(1,1) += M(i) * ( pow(xyzc(i,0)-xcg,2.0) + pow(xyzc(i,2)-zcg,2.0) );
		I(2,2) += M(i) * ( pow(xyzc(i,0)-xcg,2.0) + pow(xyzc(i,1)-ycg,2.0) );
		//Other
		I(0,1) += M(i) * (xyzc(i,0)-xcg) * (xyzc(i,1)-ycg) ;
		I(0,2) += M(i) * (xyzc(i,0)-xcg) * (xyzc(i,2)-zcg) ;
		I(1,2) += M(i) * (xyzc(i,1)-ycg) * (xyzc(i,2)-zcg) ;
	};
	I(0,1) *= -1.0;
	I(0,2) *= -1.0;
	I(1,2) *= -1.0;
	I(1,0) = I(0,1);
	I(2,0) = I(0,2);
	I(2,1) = I(1,2);

	//Elements characteristic dimension as cubic root of average cell volume
	d_c = 0.5 * ( max(Vol3) + min(Vol3) );

	myMesh.rho.free();  Vol.free();  Vol3.free();

};
void CThermal::FaceCenters(void)
{
	Range all = Range::all();

	if(Teth)
	{
		//Total number of faces
		nFacesTot = nCells*nf;
		//Coordinates of all faces
		xyzallf.resize(nFacesTot,3);
		xyzallf=0.0;
		//Connectivity of all faces
		connallf.resize(nFacesTot,nf);
		connallf=0;
		//Nodes permutation to define faces
		Array<int,2> perm(nf,nn);
		perm(0,0)=0; perm(0,1)=1; perm(0,2)=2;
		perm(1,0)=0; perm(1,1)=3; perm(1,2)=1;
		perm(2,0)=0; perm(2,1)=2; perm(2,2)=3;
		perm(3,0)=3; perm(3,1)=2; perm(3,2)=1;
		for(int i=0; i<nCells; i++)
		{
			//Face connectivity and centers coordinates
			for(int j=0; j<nn; j++)
			{
				//Face 1
				connallf(nf*i,j) = myMesh.conn(i,perm(0,j));
				//Face 2
				connallf(nf*i+1,j) = myMesh.conn(i,perm(1,j));
				//Face 3
				connallf(nf*i+2,j) = myMesh.conn(i,perm(2,j));
				//Face 4
				connallf(nf*i+3,j) = myMesh.conn(i,perm(3,j));
				//Face centers
				for(int k=0; k<3;k++)
				{
					//Face 1
					xyzallf(nf*i,k) +=  myMesh.xyz(myMesh.conn(i,perm(0,j)),k);
					//Face 2
					xyzallf(nf*i+1,k) +=  myMesh.xyz(myMesh.conn(i,perm(1,j)),k); 
					//Face 3
					xyzallf(nf*i+2,k) +=  myMesh.xyz(myMesh.conn(i,perm(2,j)),k); 
					//Face 4
					xyzallf(nf*i+3,k) +=  myMesh.xyz(myMesh.conn(i,perm(3,j)),k);
				}
			};
			xyzallf(nf*i,all) *= 1.0/3.0; xyzallf(nf*i+1,all) *= 1.0/3.0;
			xyzallf(nf*i+2,all) *= 1.0/3.0; xyzallf(nf*i+3,all) *= 1.0/3.0;
		};
	};
	//int temp = connallf.size();
	//Log << "Size of connallf: " << temp << "\n";

	//Define faces centers for hexa mesh
	if(Hexa)
	{
		//Total number of faces
		nFacesTot = nCells*nf;
		//Coordinates of all faces
		xyzallf.resize(nFacesTot,3);
		xyzallf=0.0;
		//Connectivity of all faces
		connallf.resize(nFacesTot,nf);
		connallf=0;
		//Nodes permutation to define faces
		Array<int,2> perm(nf,nn);
		perm(0,0)=0; perm(0,1)=3; perm(0,2)=2; perm(0,3)=1;
		perm(1,0)=4; perm(1,1)=5; perm(1,2)=6; perm(1,3)=7;
		perm(2,0)=4; perm(2,1)=7; perm(2,2)=3; perm(2,3)=0;
		perm(3,0)=1; perm(3,1)=2; perm(3,2)=6; perm(3,3)=5;
		perm(4,0)=0; perm(4,1)=1; perm(4,2)=5; perm(4,3)=4;
		perm(5,0)=2; perm(5,1)=3; perm(5,2)=7; perm(5,3)=6;
		for(int i=0; i<nCells; i++)
		{
			//Face connectivity and centers coordinates
			for(int j=0; j<nn; j++)
			{
				//Face 1
				connallf(nf*i,j) = myMesh.conn(i,perm(0,j)); 
				//Face 2
				connallf(nf*i+1,j) = myMesh.conn(i,perm(1,j)); 
				//Face 3
				connallf(nf*i+2,j) = myMesh.conn(i,perm(2,j)); 
				//Face 4
				connallf(nf*i+3,j) = myMesh.conn(i,perm(3,j)); 
				//Face 5
				connallf(nf*i+4,j) = myMesh.conn(i,perm(4,j)); 
				//Face 6
				connallf(nf*i+5,j) = myMesh.conn(i,perm(5,j));
				//Face centers
				for(int k=0; k<3;k++)
				{
					//Face 1
					xyzallf(nf*i,k) +=  myMesh.xyz(myMesh.conn(i,perm(0,j)),k); 
					//Face 2
					xyzallf(nf*i+1,k) +=  myMesh.xyz(myMesh.conn(i,perm(1,j)),k); 
					//Face 3
					xyzallf(nf*i+2,k) +=  myMesh.xyz(myMesh.conn(i,perm(2,j)),k); 
					//Face 4
					xyzallf(nf*i+3,k) +=  myMesh.xyz(myMesh.conn(i,perm(3,j)),k); 
					//Face 5
					xyzallf(nf*i+4,k) +=  myMesh.xyz(myMesh.conn(i,perm(4,j)),k);
					//Face 6
					xyzallf(nf*i+5,k) +=  myMesh.xyz(myMesh.conn(i,perm(5,j)),k);
				}
 			}
			xyzallf(nf*i,all) *= 0.25; xyzallf(nf*i+1,all) *= 0.25; xyzallf(nf*i+2,all) *= 0.25;
			xyzallf(nf*i+3,all) *= 0.25; xyzallf(nf*i+4,all) *= 0.25; xyzallf(nf*i+5,all) *= 0.25;
		}
	};
		
};
void CThermal::Face2Cells(void)
{

	Log << "Defining Mesh Faces and Neighbours..." << "\n";

	//Defining Cell2Face array, faces 0 to nf are shared by cell 0 etc --> Cell2Face(0:nf) = 0 ...
	Array<int,1> Cell2Face;
	Cell2Face.resize(nFacesTot);
	int n=-1;
	for(int i=0;i<nFacesTot;i++)
	{
		if(i%nf==0)
		{
			n++;
		}
		Cell2Face(i)=n;
	}

	//Search for duplicate face in the list based on face center coordinates
	int dim;
	double tol;
	ANNpointArray dataPts;	//Data points
	ANNpoint queryPt;		//Query point
	ANNidxArray nnIdx;		//Near neighbor indices
	ANNdistArray dists;		//Near neighbor distances
	ANNkd_tree* kdTree;		//Search structure

	Log << "Searching duplicate faces: " << "\n";
	dim = 3;
	tol = 1.0e-9;
	//Allocate data points
	dataPts = annAllocPts(nFacesTot, dim);
	//Allocate query pt
	queryPt = annAllocPt(dim); 

	//Face cooridnates of all the faces
	for (int i=0;i<nFacesTot;i++)
	{
		dataPts[i][0] = xyzallf(i,0);
		dataPts[i][1] = xyzallf(i,1);
		dataPts[i][2] = xyzallf(i,2);
	}

	//Creating kd tree
	Log << "Creating KD tree ..." << " ";
	kdTree = new ANNkd_tree(dataPts,nFacesTot,dim);
	nnIdx = new ANNidx[2];		
	dists = new ANNdist[2];	
	Log << "done!" << "\n";

	//Searching the kd tree for 2 closest centers to query face center: one is same center other is duplicate
	Log << "Searching KD tree ..." << " ";
	vector<int> CellId1;
	vector<int> CellId2;
	vector<int> FaceId;
	for(int i=0;i<nFacesTot;i++)
	{
		//This means that the face center has not been considered yet
		if( Cell2Face(i) != -1 ) 
		{
			//Face Id
			FaceId.push_back(i);

			//Face center under investigation
			queryPt[0] = xyzallf(i,0);
			queryPt[1] = xyzallf(i,1);
			queryPt[2] = xyzallf(i,2);

			//Search kd-tree
			kdTree -> annkSearch(queryPt,2,nnIdx,dists);
			
			//Assign cellid: if face is boundary then cellid is -1 elseif face is shared then 2 cellids
			CellId1.push_back(-1);
			CellId2.push_back(-1);
			if(sqrt(dists[0]) < tol)
			{
				CellId1.pop_back();
				CellId1.push_back(Cell2Face(nnIdx[0]));
				//If face center is considered id is set to -1 to avoid repetition
				Cell2Face(nnIdx[0]) = -1;
			}
			if(sqrt(dists[1]) < tol)
			{
				CellId2.pop_back();
				CellId2.push_back(Cell2Face(nnIdx[1]));
				//If face center is considered id is set to -1 to avoid repetition
				Cell2Face(nnIdx[1]) = -1;
			}
		}
	};
	Log << "finished searching!" << "\n";


	//Clean and free kdtree
	annDeallocPt(queryPt);
	annDeallocPts(dataPts);
	delete [] nnIdx;
	delete [] dists;
	delete kdTree;
	annClose();


	Log << "Populate faces and neighbours matrix ..." << " ";
	//Face structure and face coordinates
	nFaces = (int) CellId1.size();
	Face2Cell.resize(nFaces,2);
	xyzf.resize(nFaces,3);
	connf.resize(nFaces,nn);
	nb.resize(nCells);
	for(int i=0;i<nFaces;i++)
	{
		//Face to cells structure
		Face2Cell(i,0) = CellId1[i];
		Face2Cell(i,1) = CellId2[i];
		//Face centers coordinates
		xyzf(i,0) = xyzallf(FaceId[i],0);
		xyzf(i,1) = xyzallf(FaceId[i],1);
		xyzf(i,2) = xyzallf(FaceId[i],2);
		//Face connectivity
		for(int j=0;j<nn;j++)
		{
			connf(i,j) = connallf(FaceId[i],j);
		}
		//Neighbors of cell
		int c0 = Face2Cell(i,0);
		int c1 = Face2Cell(i,1);
		//Two cells shared
		if(c0!=-1 && c1!=-1)
		{
			//Cell id as boundary
			nb(c0).push_back( c1 );
			nb(c1).push_back( c0 );
		};
		//Boundary cell
		if(c0==-1 || c1==-1)
		{
			int cb = (c0>c1) ? c0 : c1; 
			//(Faceid - 1) as boundary
			nb(cb).push_back( -i-1 );
		};
	};
	Log << "done!" << "\n";
	
	Log << "Faces & Neighbours Definition Completed!" << "\n";
	Log << "\n";

	//free variables
	xyzallf.free();	
	
}
void CThermal::FaceAf(void)
{
	Range all = Range::all();

	//Find area, normal vector and area normal vector	
	double x1,x2,y1,y2,z1,z2;
	xyzfn.resize(nFaces,3);
	Af.resize(nFaces);
	Afn.resize(nFaces,3);
	for(int i=0;i<nFaces;i++)
	{
		//x component vector 1 and 2
		x1 =  myMesh.xyz(connf(i,1),0) -  myMesh.xyz(connf(i,0),0);
		x2 =  myMesh.xyz(connf(i,2),0) -  myMesh.xyz(connf(i,0),0);
		//y component vector 1 and 2
		y1 =  myMesh.xyz(connf(i,1),1) -  myMesh.xyz(connf(i,0),1);
		y2 =  myMesh.xyz(connf(i,2),1) -  myMesh.xyz(connf(i,0),1);
		//z component vector 1 and 2
		z1 =  myMesh.xyz(connf(i,1),2) -  myMesh.xyz(connf(i,0),2);
		z2 =  myMesh.xyz(connf(i,2),2) -  myMesh.xyz(connf(i,0),2);
		//Cross product
		xyzfn(i,0) = y1*z2 - z1*y2;
		xyzfn(i,1) = z1*x2 - x1*z2;
		xyzfn(i,2) = x1*y2 - y1*x2;
		//Cross product magnitude (Area of paralellogram or 2xArea of triangle)
		Af(i) = sqrt( pow(xyzfn(i,0),2.0) + pow(xyzfn(i,1),2.0) + pow(xyzfn(i,2),2.0) );
		//Normal vector
		xyzfn(i,all) /= Af(i);
		
		//Teth mesh: face area is half the area of the parallelogram obtained with cross product
		if(Teth)
		{
			//Face area
			Af(i) *= 0.5;
			//Area vector
			Afn(i,all) = Af(i) * xyzfn(i,all);
		};
		//Hexa mesh: face area is sum of area of two triangles obtained by cross product
		if(Hexa)
		{
			//x component vector 1 and 2
			x1 =  myMesh.xyz(connf(i,1),0) -  myMesh.xyz(connf(i,3),0);
			x2 =  myMesh.xyz(connf(i,2),0) -  myMesh.xyz(connf(i,3),0);
			//y component vector 1 and 2
			y1 =  myMesh.xyz(connf(i,1),1) -  myMesh.xyz(connf(i,3),1);
			y2 =  myMesh.xyz(connf(i,2),1) -  myMesh.xyz(connf(i,3),1);
			//z component vector 1 and 2
			z1 =  myMesh.xyz(connf(i,1),2) -  myMesh.xyz(connf(i,3),2);
			z2 =  myMesh.xyz(connf(i,2),2) -  myMesh.xyz(connf(i,3),2);
			//Cross product
			Array<double,1> xyzfn2(3);
			xyzfn2(0) = y1*z2 - z1*y2;
			xyzfn2(1) = z1*x2 - x1*z2;
			xyzfn2(2) = x1*y2 - y1*x2;
			//Cross product magnitude (area of second triangle)
			double Af2 = 0.5 * (sqrt( pow(xyzfn2(0),2.0) + pow(xyzfn2(1),2.0) + pow(xyzfn2(2),2.0) ));
			//Face area
			Af(i) = 0.5*Af(i) + Af2;
			//Area vector
			Afn(i,all) = Af(i) * xyzfn(i,all);
		}
	};

};
Array<double,2> CThermal::MTMInverse(Array<double,2> A)
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

};
void CThermal::ThermalOutput(string body)
{
	string rev_n;
	char buffer[500];
	//revolution number
	_itoa_s(myGapResult.revcounter,buffer,10);
	rev_n = buffer;

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
			fout << scientific <<  myMesh.xyz(i,0) << "\t" <<  myMesh.xyz(i,1) << "\t" <<  myMesh.xyz(i,2) << "\n";
		}
		fout << "CELLS " << "\t" << nCells << "\t" << 9*nCells << "\n";
		for(int i = 0; i < nCells; i++)
		{
			fout << 8 << "\t" << myMesh.conn(i,0) << "\t" << myMesh.conn(i,1)  << "\t" << myMesh.conn(i,2) << "\t" << myMesh.conn(i,3) << "\t" << myMesh.conn(i,4) << "\t" << myMesh.conn(i,5) << "\t" << myMesh.conn(i,6) << "\t" << myMesh.conn(i,7) << "\n";
		}
		fout << "CELL_TYPES " << nCells << "\n";
		for(int i = 0; i < nCells; i++)
		{
			fout << 12 << "\n";
		}
		fout << "CELL_DATA " << nCells << "\n";
		fout << "SCALARS T_[C] double" << "\n";
		fout << "LOOKUP_TABLE default" << "\n";
		for(int i = 0; i < nCells; i++)
		{
			fout << scientific << phi[i] << "\n";
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
			fout << scientific <<  myMesh.xyz(i,0) << "\t" <<  myMesh.xyz(i,1) << "\t" <<  myMesh.xyz(i,2) << "\n";
		}
		fout << "CELLS " << "\t" << nCells << "\t" << 5*nCells << "\n";
		for(int i = 0; i < nCells; i++)
		{
			fout << 4 << "\t" << myMesh.conn(i,0) << "\t" << myMesh.conn(i,1)  << "\t" << myMesh.conn(i,2) << "\t" << myMesh.conn(i,3) << "\n";
		}
		fout << "CELL_TYPES " << nCells << "\n";
		for(int i = 0; i < nCells; i++)
		{
			fout << 10 << "\n";
		}
		fout << "POINT_DATA " << nNodes << "\n";
		fout << "SCALARS T_[C] double" << "\n";
		fout << "LOOKUP_TABLE default" << "\n";
		for(int i = 0; i < nNodes; i++)
		{
			fout << scientific << phi[i] << "\n";
		}
		fout << "\n";
		//fout << "POINT_DATA " << nNodes << "\n";
		fout << "SCALARS Boundary int 1" << "\n";
		fout << "LOOKUP_TABLE default" << "\n";
		for(int i = 0; i < nNodes; i++)
		{
			fout << scientific << boundaries[i] << "" << "\n";
		}
		fout << "SCALARS q_[W/m2] double" << "\n";
		fout << "LOOKUP_TABLE default" << "\n";
		for(int i = 0; i < nNodes; i++)
		{
			fout << scientific << vecq[i] << "" << "\n";
		}
		fout.close();
		fout.clear();
	}


	
};
double CThermal::SetBMatrix(int cell){
	//Load data about the cell in question
	//node 1
	int node = myMesh.conn(cell,0);
	double x1 = myMesh.xyz(node,0);
	double y1 = myMesh.xyz(node,1);
	double z1 = myMesh.xyz(node,2);
	//Log << "Node: " << node << "\t" << x1 << "\t" << y1 << "\t" << z1 << "\n";
	//node 2
	node = myMesh.conn(cell,1);
	double x2 = myMesh.xyz(node,0);
	double y2 = myMesh.xyz(node,1);
	double z2 = myMesh.xyz(node,2);
	//Log << "Node: " << node << "\t" << x2 << "\t" << y2 << "\t" << z2 << "\n";
	//node 3
	node = myMesh.conn(cell,2);
	double x3 = myMesh.xyz(node,0);
	double y3 = myMesh.xyz(node,1);
	double z3 = myMesh.xyz(node,2);
	//Log << "Node: " << node << "\t" << x3 << "\t" << y3 << "\t" << z3 << "\n";
	//node 4
	node = myMesh.conn(cell,3);
	double x4 = myMesh.xyz(node,0);
	double y4 = myMesh.xyz(node,1);
	double z4 = myMesh.xyz(node,2);
	//Log << "Node: " << node << "\t" << x4 << "\t" << y4 << "\t" << z4 << "\n";
	//Log << "\n";

	//Clear and resize the B matrix
	matB.clear();
	vector<double> temp;temp.resize(4);
	matB.resize(3,temp);
	temp.clear();
	
	//Fill in B matrix
	matB[0][0] = (y2-y4)*(z3-z4)-(y3-y4)*(z2-z4);
	matB[0][1] = (y3-y4)*(z1-z4)-(y1-y4)*(z3-z4);
	matB[0][2] = (y1-y4)*(z2-z4)-(y2-y4)*(z1-z4);
	matB[1][0] = (x3-x4)*(z2-z4)-(x2-x4)*(z3-z4);
	matB[1][1] = (x1-x4)*(z3-z4)-(x3-x4)*(z1-z4);
	matB[1][2] = (x2-x4)*(z1-z4)-(x1-x4)*(z2-z4);
	matB[2][0] = (x2-x4)*(y3-y4)-(x3-x4)*(y2-y4);
	matB[2][1] = (x3-x4)*(y1-y4)-(x1-x4)*(y3-y4);
	matB[2][2] = (x1-x4)*(y2-y4)-(x2-x4)*(y1-y4);

	//Calculate Jacobian Determinant
	double J = (x1-x4)*matB[0][0]+(y1-y4)*matB[1][0]+(z1-z4)*matB[2][0]; //determinant of the Jacobian is consistently negative.  Could be a problem to check dwm.

	//Log << "Element " << cell << " J = " << J << "\n";
	
	for(int a = 0;a<3;a++){
		for( int b = 0; b<3; b++)
			matB[a][b] *= (1/J);
		double temp = matB[a][0]+matB[a][1]+matB[a][2];
		temp = temp * -1.0;
		matB[a][3] = temp;
	};

	/*Log << "B =" << "\n";
	for(int a = 0;a<3;a++){
		for(int b = 0;b<4;b++){
			Log << matB[a][b] << "\t";
		}
		Log << "\n";
	}
	Log << "\n";*/

	/*Debug check
	cout << "For Cell " << cell << ", the B matrix is:" << "\n";
	for (int a = 0;a<3;a++){
		for (int b = 0;b<3;b++){
			cout << matB[a][b] << '\t';
		}
		cout << "\n";
	}*/
	J = abs(J);
	return(J);
};
vector<vector<double > > CThermal::CalcKcd(vector<vector<double> > B,double lambda,double J){
	//multiply BtB
	vector<vector<double > > Kcd;
	vector<double> temp; temp.resize(4,0);
	Kcd.resize(4,temp);
	temp.clear();

	//Log << "Kcd" << "\n";
	for(int a = 0; a < 4; a++){
		for(int b = 0; b < 4; b++){
			for(int c = 0; c < 3; c++){
				double temp = B[c][a] * B[c][b];
				temp = temp * lambda * J / ( 6 );
				Kcd[a][b] = Kcd[a][b] + temp;
			}
			//Log << Kcd[a][b] ;
		}
		//Log << "\n";
	}
	//Log << "\n";*/

	/*Log << "Local Kcd =" << "\n";
	for(int a = 0;a<4;a++){
		for(int b = 0;b<4;b++){
			Log << Kcd[a][b] << "\t";
		}
		Log << "\n";
	}
	Log << "\n";*/

	return(Kcd);
};
double CThermal::FaceArea(int node1, int node2, int node3){
	double x1 = myMesh.xyz(node1,0);
	double y1 = myMesh.xyz(node1,1);
	double z1 = myMesh.xyz(node1,2);
	double x2 = myMesh.xyz(node2,0);
	double y2 = myMesh.xyz(node2,1);
	double z2 = myMesh.xyz(node2,2);
	double x3 = myMesh.xyz(node3,0);
	double y3 = myMesh.xyz(node3,1);
	double z3 = myMesh.xyz(node3,2);

	//make two vectors
	double dx1 = x3-x1;
	double dy1 = y3-y1;
	double dz1 = z3-z1;
	double dx2 = x3-x2;
	double dy2 = y3-y2;
	double dz2 = z3-z2;

	

	double ax1 = pow(dy1*dz2-dz1*dy2,2);
	double ay1 = pow(dz1*dx2-dx1*dz2,2);
	double az1 = pow(dx1*dy2-dy1*dx2,2);

	

	double area = 0.5 * sqrt(ax1+ay1+az1);

	

	return(area);
};