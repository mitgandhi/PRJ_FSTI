#include "../../mesh/src/CMesh.h"
#include "../../main/src/sGapResult.h"
#include "../../main/src/CGapInput.h"
#include "../../main/src/CGapUtils.h"
#include "../../gmm/gmm_kernel.h"
#include "../../gmm/gmm_iter_solvers.h"
#include "../../gmm/gmm_inoutput.h"
#pragma once


class CThermal
{
	public:

	ofstream fout;

	//Default constructor-destructor
	CThermal(void);
	~CThermal(void);


	int nNodes;
	int nCells;
	bool Teth,Hexa;
	double Mtot,xcg,ycg,zcg,d_c;
	//Total face number
	int nFacesTot;
	//Net face number
	int nFaces;
	//Number of nodes per face
	int nn;
	//Number of faces per cell
	int nf;
	//id neumann faces
	vector<vector<int>> faceid_qb;
	//id dirichlet faces
	vector<vector<int>> faceid_phib;
	//id mixed faces
	vector<vector<int>> faceid_mixb;
	//id gap nodes
	vector<vector<int>> nodeid_qb;

	//Scalar fields boundary
	Array<double,1> phib;		//Boundary faces values Dirichlet 
	Array<double,1> qb;			//Boundary faces values Neumann
	Array<double,1> phibinf;	//Boundary faces convection value mixed
	Array<double,1> hb;			//Boundary faces convection value mixed

	//Cell values
	//Array<double,1> phi;	
	//Surface value
	Array<double,1> phi_surf;
	//Face boundary type id
	Array<int,1> FaceBoundId;
	//Coordinates of all faces
	Array<double,2> xyzallf;
	//Connectivity of all faces
	Array<int,2> connallf;
	//List of faces and shared cells
	Array<int,2> Face2Cell;
	//Cell centers coordinates
	Array<double,2> xyzc;
	//Net faces coordinates
	Array<double,2> xyzf;
	//Net faces connectivity
	Array<int,2> connf;
	//Faces normal vector
	Array<double,2> xyzfn;
	//Faces area
	Array<double,1> Af;
	//Faces normal area vector
	Array<double,2> Afn;
	//Cell volumes
	Array<double,1> Vol;
	Array<double,1> M;
	Array<double,1> rho;
	Array<double,2> I;
	//Neighbors primary diffusion coeffcients
	Array<vector<double>,1> anb;
	//Neighbors id
	Array<vector<int>,1> nb;
	//Source term
	Array<double,1> b;
	//Old secondary gradient
	Array<double,1> SecGrad;
	//Cell diffusion coefficient
	Array<double,1> ap;
	//Cell geometric matrix for cell gradient calculations
	vector<Array<double,2>> C;
	//Cell gradient
	vector<Array<double,1>> Gradc;
	//Face gradient
	vector<Array<double,1>> Gradf;



	//Read mesh
	void readMeshThermal(string body);
	//Solve thermal problem
	void ThermalSolve(string body,vector<Array<double,1>> qbi,Array<double,1> &T_body,Array<double,1> &T_surface);	
	void ThermalSolve_Piston(vector<Array<double,1>> qbi, vector<Array<double,1>> wtd_h, vector<double> Tb, Array<double,1> &T_body,Array<double,1> &T_surface);	
	void ThermalSolve_Cylinder(vector<Array<double,1>> qbi, vector<Array<double,1>> wtd_h, vector<double> Tb, Array<double,1> &T_body,Array<double,1> &T_surface);	
	void ThermalSolve_interbody_conduction(vector<vector<int>> &Mid_K2B, vector<vector<double>> &Mwt_K2B, vector<vector<int>> &Mid_B2K, vector<vector<double>> &Mwt_B2K, double C);	
	void ThermalSolve_new(void);	
	void ThermalSolve_post(string body,vector<Array<double,1>> qbi,Array<double,1> &T_body,Array<double,1> &T_surface);	
	//Check all part si bounded and no boundary faces are missing
	void CheckBoundaryFaces(string body);		
	//Set boundary faces from mesh analysis
	void BoundaryFaces(void);		
	//Set boundary conditons arrays
	void SetBoundaries(vector<Array<double,1>> qbi);
	//Define Dirichlet phi at boundaries 
	void Facephib(void);	
	//Define Neumann q at boundaries 
	void Faceqb(void);				
	//Define mixed condition at boundaries 
	void Facemixb(void);	
	//Write xyz cooridnates for gap surface
	void WriteGapSurfacexyz(int body);	

	//Calculate cell centers
	void CellCenters(void);	
	//Calculate cell valume and mass
	void CellMassVolume(void);	
	//Calculate face centers
	void FaceCenters(void);	
	//Calculate face 2 cells list: array containig cells shared by each face
	void Face2Cells(void);
	//Calculate normal to face and area
	void FaceAf(void);	

	//Calculate cells sec gradient			
	Array<double,2> MTMInverse(Array<double,2> A);		

	//VTK output
	void ThermalOutput(string body);	

	double SetBMatrix(int cell);
	double FaceArea(int node1, int node2, int node3);
	vector<vector<double > > CalcKcd(vector<vector<double> > B,double lambda,double J);

	vector<vector<double>> matN	;//Natural Coordinate Matrix
	vector<vector<double>> matB;	//Gradient of shape functions
	
	vector<double> vecq;			//Load Vector
	vector<double> phi_both;
	vector<double> phi;
	gmm::row_matrix<gmm::wsvector<double > > matK;

	int n_piston_Nodes;
	Array<int,2> connf_piston;
	vector<vector<double>> area_piston;
	vector<vector<int>> faceid_qb_piston;
	Array<int,2> connf_cylinder;
	vector<vector<double>> area_cylinder;
	vector<vector<int>> faceid_qb_cylinder;

	

	//debug
	vector<short> boundaries;

};