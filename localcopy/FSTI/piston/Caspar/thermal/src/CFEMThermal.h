#include "CThermal.h"
#include "../../main/src/sGapResult.h"
#include "../../fem/src/femFunctions.h"
#pragma once


class CFEMThermal
{
	public:
	
	ofstream fout;

	//Default constructor-destructor
	CFEMThermal(void);
	~CFEMThermal(void);

	//locals
	int nNodes;
	int nCells;
	bool Teth,Hexa;
	double Mtot;
	double xcg,ycg,zcg;


	//Scalar fields boundary
	Array<double,1> phi;
	Array<double,2> udisp;
	//vtk check outputs
	Array<double,1> def_vtk;
	//vtk check outputs
	Array<double,1> def_vtk_x;
	//vtk check outputs
	Array<double,1> def_vtk_y;
	//vtk check outputs
	Array<double,1> def_vtk_z;
	//mass
	Array<double,1> M;
	//mass
	Array<vector<int>,1> cxyz_IR;

	//Thermal load
	void FEMThermalGeneral(void);
	void FEMThermalLoads(Array<double,1> T_body);
	void FEMThermalConstraints(void);
	void FEMThermalSolve(string body,Array<double,1> T_body,Array<double,1> &def_surf,string mesh_name);
	Array<double,1> FEMThermalSurfaceDeformation(void);
	void FEMThermalOutput(double scale,string body);
	//Inertia relief
	void FEMThermalInertiaRelief(Array<double,1> &loads);
	Array<double,2> FEMThermalInertiaInverse(Array<double,2> A);
	Array<double,2> FEMThermalCellInertia(int i);
	void FEMThermalInertiaConstraints(void);

	//Contact Algorithm
	double centershift;

};