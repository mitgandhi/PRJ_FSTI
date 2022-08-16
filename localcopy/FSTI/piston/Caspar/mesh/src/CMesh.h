#include "../../main/src/CGapInput.h"
#pragma once


class CMesh
{
	public:

	ofstream fout;

	//Default constructor-destructor
	CMesh(void);
	~CMesh(void);

	//Nodes number
	int nNodes;
	//Cells number
	int nCells;
	//Hexa or Teth mesh boolean
	bool Hexa,Teth,IR_p,IR_c,body;
	double d_tol_p,d_tol_c;
	//Boundary faces defined by user
	int nf_bd_user;
	//Boundary nodes defined by user
	int nn_bd_user;
	
	//string body mesh type
	string piston_type;
	string cylinder_type;
	//string body mesh name
	string piston_name;
	string cylinder_name;
	//folder_paths
	string piston_path_th;
	string cylinder_path_th;
	//influence matrices path
	string IM_piston_path;
	string IM_cylinder_path;

	//Nodes coordinates
	Array<double,2> xyz;
	//Cell Connectivity
	Array<int,2> conn;
	//Node constraints thermal load
	Array<vector<int>,1> cxyz_th;
	//Thermal conductivity
	Array<double,1> Tf;
	//Elastic modulus 
	Array<double,1> E;
	//Poisson ratio
	Array<double,1> v;
	//Density
	Array<double,1> rho;
	//linear expansion coefficient
	Array<double,1> alpha;
	//input parameters materials
	vector<double> namei;
	vector<double> Tfi;		
	vector<double> Ei;		
	vector<double> vi;		
	vector<double> rhoi;
	vector<double> alphai;
	//input parameters boundaires
	int qbi;
	vector<double> phibi;		
	vector<double> hbi;		
	vector<double> phibinfi;		

	//nodeset neumann boundary abaqus
	vector<vector<int>> qbABQS;
	//nodeset dirichelt boundary abaqus
	vector<vector<int>> phibABQS;
	//nodeset mixed boundary abaqus
	vector<vector<int>> mixbABQS;
	
	//functions
	void ReadABAQUSTeth(string mesh_path);		//Read teth mesh
	void ReadABAQUSHexa(string mesh_path,string mesh_name);		//Read hexa mesh
	void ReadMaterials(void);						//Read material properties
	void ReadBoundaries(void);						//Read boundaries
	//void ReadMeshMain(void);									//Read main


};