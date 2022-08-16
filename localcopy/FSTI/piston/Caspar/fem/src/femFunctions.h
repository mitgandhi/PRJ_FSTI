# include <vector>

namespace fem {
	
	// element type
	enum {
		BRICK,
		TETRA
	};

	// functions to read from input files
	void readGeneral(const char* path, int &elmType, int &intOp, int &nDof, int &nElements, int &nNodes, int &nc, int &nst);
	void readGeometry(const char* path, int nNodes, std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, const char* mode);
	void readConnectivity(const char* path, int elmType, int nElements, std::vector<int>& connect, const char* mode);
	void readConstraints(const char* path, int nc, std::vector<int>& iConstrains, const char* mode);
	void readMaterials(const char* path, int nElements, std::vector<double>& E, std::vector<double>& nu, const char* mode);
	void readLoads(const char* path, int nst, std::vector<int>& iForces, std::vector<double>& forces, const char* mode);

	void stiffness(const char* path, const char* inputType);
	void loads(const char* path, const char* inputType);
	int solve(const char* path, const int iterMax, const double convergeTol);
	void writeDisplacement(const char* path, const char* inputType);

}
