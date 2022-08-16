# include <fstream>
# include <iostream>
# include <sstream>
# include <vector>
# include <string>
# include "femFunctions.h"
# include "../../main/src/logger.h"
# include <iomanip>

using namespace std;

// General file can be only in ASCII mode
// the structure is:
//
// intOp	nDof	nElements	nNodes	nc	nst; 
//
// where
//
// intOp:		is the integral operator 
// nDof:		is the number of degree of freedom 1,2 or 3
// nElements:	is the number of element in the system
// nNodes:		is the number of nodes in the system
// nc:			is the number of the total dof (nDof*nNodes) constrined
// nst:			is the number of surface traction (is a portion of nDof*nNodes)
//

void fem::readGeneral(const char* path, int &elmType, int &intOp, int &nDof, int &nElements, int &nNodes, int &nc, int &nst) {
	
	ifstream general(string(string(path) + string("/temp/piston/input/general.fem")).c_str());
	if (!general.is_open()) {
		Log << "Error opening general.fem file" << endl;
		//system("pause");
		exit(1);
	}

	string tmp;
	string line;

	while(getline(general,line)) {
		istringstream iss (line,istringstream::in);
		while(iss) {
			iss >> tmp;
			size_t comment  = tmp.find("//"); // check if there is a comment
			if (comment!=string::npos) {
				break;
			}
			if (tmp == "elmType") {
				iss >> elmType;
				break;
			}
			if (tmp == "intOp") {
				iss >> intOp;
				break;
			}
			if (tmp == "nDof") {
				iss >> nDof;
				break;
			}
			if (tmp == "nElements") {
				iss >> nElements;
				break;
			}
			if (tmp == "nNodes") {
				iss >> nNodes;
				break;
			}
			if (tmp == "nc") {
				iss >> nc;
				break;
			}
			if (tmp == "nst") {
				iss >> nst;
				break;
			}
		}
	}
	
	general.close(); // close file
}

// Geometry file
// Coordinates of all the nodes of the system
// Can be either in ASCII format or in BINARY format
//
// ASCII version
// list of coordinates in the following format
//
// x_i	y_i	z_i
//
// the list size if nNodes
//
// BINARY format
// the coordinates are stored with the format {x} {y} {z}
// ex of code
// file.write(reinterpret_cast<char*>(x),nNodes*sizeof(double));
// file.write(reinterpret_cast<char*>(y),nNodes*sizeof(double));
// file.write(reinterpret_cast<char*>(z),nNodes*sizeof(double));
//

void fem::readGeometry(const char* path, int nNodes, vector<double>& x, vector<double>& y, vector<double>& z, const char* mode) {

	x.resize(nNodes);
	y.resize(nNodes);
	z.resize(nNodes);

	if(!strcmp(mode,"ascii")) {
		// open
		ifstream geometry(string(string(path) + string("/temp/piston/input/geometry.fem")).c_str());
		if (!geometry.is_open()) {
			Log << "Error opening geometry file" << endl;
			//system("pause");
			exit(1);
		}
		// read
		for (int i=0; i < nNodes; i++)
			geometry >> x[i] >> y[i] >> z[i];
		// close
		geometry.close();
	}

	else if(!strcmp(mode,"binary")) {
		// storage allocation
		double* _x = new double[nNodes];
		double* _y = new double[nNodes];
		double* _z = new double[nNodes];
		// open
		ifstream geometry(string(string(path) + string("/temp/piston/input/geometry.fem")).c_str(), std::ios::in | std::ios::binary );
		if (!geometry.is_open()) {
			Log << "Error opening geometry file" << endl;
			//system("pause");
			exit(1);
		}
		// read
		geometry.read(reinterpret_cast<char*>(_x), nNodes*sizeof(double));
		geometry.read(reinterpret_cast<char*>(_y), nNodes*sizeof(double));
		geometry.read(reinterpret_cast<char*>(_z), nNodes*sizeof(double));
		// close
		geometry.close();
		
		for(int i=0; i<nNodes; i++) {
			x[i] = _x[i], y[i] = _y[i], z[i] = _z[i];
		}

		delete[] _x;
		delete[] _y;
		delete[] _z;
	}
	else {
		Log << "Wrong input mode: use ascii or binary mode!" << endl;
		//system("pause");
		exit(1);
	}
}

// Connettivity file
// List of nodes defining the elements (1 element -> 8 nodes)
// Can be either in ASCII format or in BINARY format
// 
// ASCII version
// List (size = nElements), each line 8 entries (nodes defining the element)
//
// n_0	n_1	n_2	n_3	n_4 n_5	n_6	n_7
//
// BINARY version
// one array of size 8*nElements with all the connection
// { ... [n_0, n_1, ... n_7], [n_0, n_1, ... n_7], ... }
//


void fem::readConnectivity(const char* path, int elmType, int nElements, vector<int>& connect, const char* mode) {

	int elementNodes = 0;

	if(elmType == BRICK)
		elementNodes = 8;
	else if(elmType == TETRA)
		elementNodes = 4;

	connect.resize(elementNodes*nElements);

	if(!strcmp(mode,"ascii")) {
		// open
		ifstream connectivity(string(string(path) + string("/temp/piston/input/connectivity.fem")).c_str());
		if (!connectivity.is_open()) {
			Log << "Error opening connectivity.fem file" << endl;
			//system("pause");
			exit(1);
		}
		// read
		for (int i=0, c = 0; i < nElements; i++)
			for(int j = 0; j < elementNodes; j++, c++)
				connectivity >> connect[c];
		// close
		connectivity.close();
	}
	else if(!strcmp(mode,"binary")) {
		// open
		ifstream connectivity(string(string(path) + string("/temp/piston/input/connectivity.fem")).c_str(), std::ios::in | std::ios::binary );
		if (!connectivity.is_open()) {
			Log << "Error opening connectivity file" << endl;
			//system("pause");
			exit(1);
		}

		// read
		int* _connect = new int[elementNodes*nElements];
		connectivity.read(reinterpret_cast<char*>(_connect), elementNodes*nElements*sizeof(int));
		// close
		connectivity.close();

		for(int i=0; i<elementNodes*nElements; i++) {
			connect[i] = _connect[i];
		}

		delete[] _connect;

	}
	else {
		Log << "Wrong input mode: use ascii or binary mode!" << endl;
		//system("pause");
		exit(1);
	}
}

// Constrains file
// List of degrees of freedom constrained
// The degrees of freedom of the entire system is represented by a list from 0 to nDof*nNodes
// Each nodes has nDof degree of freedom, so the node "i" has is represented by the entries
//		i*nDof, i*nDof + 1, ... , i*nDof + nDof - 1
// 
// Can be either in ASCII format or in BINARY format
//
// ASCII and BINARY version
// list of dof constrained
// if for example nDof = 3, the system has 3 elements (24 nodes) and the list is
// 3,4,6,7,8,13
// it means that:
//		element 1 is constrained in y and z direction
//		element 3 is constrained in x, y and z direction
//		element 5 is constrained in y direction
//

void fem::readConstraints(const char* path, int nc, vector<int>& iConstrains, const char* mode) {
	
	iConstrains.resize(nc);

	if(!strcmp(mode,"ascii")) {
		// open
		ifstream constrains(string(string(path) + string("/temp/piston/input/constraints.fem")).c_str());
		if (!constrains.is_open()) {
			Log << "Error opening constrains.fem file" << endl;
			//system("pause");
			exit(1);
		}
		// read
		for (int i=0; i < nc; i++)
			constrains >> iConstrains[i];
		// close
		constrains.close();
	}
	else if(!strcmp(mode,"binary")) {
		// open
		ifstream constrains(string(string(path) + string("/temp/piston/input/constraints.fem")).c_str(), std::ios::in | std::ios::binary );
		if (!constrains.is_open()) {
			Log << "Error opening connectivity file" << endl;
			//system("pause");
			exit(1);
		}

		// read
		int* _iConstrains = new int[nc];
		constrains.read(reinterpret_cast<char*>(_iConstrains), nc*sizeof(int));
		// close
		constrains.close();

		for(int i = 0; i < nc; i++) {
			iConstrains[i] = _iConstrains[i];
		}

		delete[] _iConstrains;

	}
	else {
		Log << "Wrong input mode: use ascii or binary mode!" << endl;
		//system("pause");
		exit(1);
	}
}

// Material file
// For each element are indicated the Elastic module and the poisson ratio
// Can be either in ASCII format or in BINARY format
//
// ASCII version
// List (size = nElements), for each entry
//
// E	nu
//
// BINARY version
// Array for E followed by array for nu
// ex of code
// file.write(reinterpret_cast<char*>(E),nElements*sizeof(double));
// file.write(reinterpret_cast<char*>(nu),nElements*sizeof(double));

void fem::readMaterials(const char* path, int nElements, vector<double>& E, vector<double>& nu, const char* mode) {
	
	E.resize(nElements);
	nu.resize(nElements);

	if(!strcmp(mode,"ascii")) {
		// open
		ifstream materials(string(string(path) + string("/temp/piston/input/materials.fem")).c_str());
		if (!materials.is_open()) {
			Log << "Error opening materials.fem file" << endl;
			//system("pause");
			exit(1);
		}
		// read
		for (int i=0; i < nElements; i++)
			materials >> E[i] >> nu[i];
		// close
		materials.close();
	}
	else if(!strcmp(mode,"binary")) {
		// open
		ifstream materials(string(string(path) + string("/temp/piston/input/materials.fem")).c_str(), std::ios::in | std::ios::binary );
		if (!materials.is_open()) {
			Log << "Error opening materials file" << endl;
			//system("pause");
			exit(1);
		}
		// read E
		double _E;
		double _nu;
		for(int i=0; i<nElements; i++) {
			materials.read(reinterpret_cast<char*>(&_E), sizeof(double));	
			E[i] = _E;
		}
		for(int i=0; i<nElements; i++) {
			materials.read(reinterpret_cast<char*>(&_nu), sizeof(double));	
			nu[i] = _nu;
		}
				
		// close
		materials.close();

	}
	else {
		Log << "Wrong input mode: use ascii or binary mode!" << endl;
		//system("pause");
		exit(1);
	}
}

// Loads file
// each node can be load in nDof directions
// Can be either in ASCII format or in BINARY format
// 
// ASCII version
// List (size = nst) each entry
//
// iForce	force
//
// where
//		iForce is the label of the degree of freedom that is loaded (see constrains file)
//		force is the magnitude (in N) if the force applied
//
// BINARY version
// array of iForces followed by array of forces

void fem::readLoads(const char* path, int nst, vector<int>& iForces, vector<double>& forces, const char* mode) {
	
	iForces.resize(nst);
	forces.resize(nst);

	if(!strcmp(mode,"ascii")) {
		// open
		ifstream loads(string(string(path) + string("/temp/piston/input/loads.fem")).c_str());
		if (!loads.is_open()) {
			Log << "Error opening loads.fem file" << endl;
			//system("pause");
			exit(1);
		}
		// read
		for (int i=0; i < nst; i++) {
			loads >> iForces[i] >> forces[i];
		}
		// close
		loads.close();
	}
	else if(!strcmp(mode,"binary")) {
		// open
		ifstream loads(string(string(path) + string("/temp/piston/input/loads.fem")).c_str(), std::ios::in | std::ios::binary );
		if (!loads.is_open()) {
			Log << "Error opening loads.fem file" << endl;
			//system("pause");
			exit(1);
		}
		// read
		int* _iForces = new int[nst];
		double* _forces = new double[nst];
		loads.read(reinterpret_cast<char*>(_iForces), nst*sizeof(int));
		loads.read(reinterpret_cast<char*>(_forces), nst*sizeof(double));
		// close
		loads.close();

		for(int i = 0; i < nst; i++) {
			iForces[i] = _iForces[i];
			forces[i] = _forces[i];
		}

		delete[] _iForces;
		delete[] _forces;
	}
	else {
		Log << "Wrong input mode: use ascii or binary mode!" << endl;
		///system("pause");
		exit(1);
	}	
}