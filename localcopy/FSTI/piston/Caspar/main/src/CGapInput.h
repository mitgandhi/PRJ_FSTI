#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <direct.h>
#include <time.h>
#include <windows.h>
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#include <blitz\array.h>
#include <ANN.h>
#pragma once

using namespace std;
using namespace blitz;
using namespace blitz::tensor;

const double PI = 4.0*atan(1.0);

class CGapInput
{

	public:


	//Constructor destructor
	CGapInput(void);
	~CGapInput(void);


	//----------------------------------------------------------------------------------------------------//
	//-------------------------DECLARATIONS FOR GENERAL INPUT FILE .GEN-----------------------------------//
	//----------------------------------------------------------------------------------------------------//
	/*struct sGeneral
	{
		int mode;										//Working mode
		int npistons;									//Number of pistons
		int nrevolutions;								//Number of revolutions
		double speed;									//Shaft speed
		double dB;										//Cylinder block pitch diameter
		double dK;										//Piston diameter
		double beta;									//Swashplate angle
		double gamma;									//Swashplate cross angle
		double betamax;									//Maximum swashplate angle
		double pHP;										//High pressure
		double pLP;										//Low pressure
		double pCase;									//Case pressure
		double THP;									//Temperature high pressure
		double TLP;									//Temperature low pressure
		double TCase;									//Temperature leakage
		int oiltype;									//Oil type
		double oildensity;								//Oil density at reference pressure
		double oilbetaP;								//Density pressure coefficient
		double oilbetaT;								//Density thermal expansion coefficient
		double oilK;									//Oil bulk modulus at reference point
		double oilbetaKP;								//Oil bulk modulus pressure coefficient
		double oilbetaKT;								//Oil bulk modulus thermal coefficient
		double oilviscosity;							//Oil dynamic viscosity
		double oilW;									//Oil kinematic viscosity weightening factor
		double oilTc1;									//Kinematic viscosity temperature coefficient 1
		double oilPc1;									//Kinematic viscosity pressure coefficent 2
		double oilTc2;									//Kinematic viscosity tempearature coefficent 2
		double oilPc2;									//Kinematic viscosity pressure coefficient 2
		double oillambda;								//Oil thermal conductivity
		double oilC;									//Oil heat capacity
		double alpha1;									//Alpha 1 used for oil
		double alpha2;									//Alpha 2 used for oil
		double alpha3;									//Alpha 3 used for oil
	};
	sGeneral General;*/
	

	//----------------------------------------------------------------------------------------------------//
	//-------------------------DECLARATIONS FOR GEOMETRY INPUT FILE .GEO----------------------------------//
	//----------------------------------------------------------------------------------------------------//
	struct sGeometry
	{
		//double dZ;										//Bushing diameter
		//double lK;										//Piston length
		//double lF;										//Bushing length
		//double le;										//Bushing beginning position
		//double lZ0;										//Displacement chamber lenght outer dead point
		//double dDK;										//Orifice piston head diameter
		//double doutG;									//Slipper outer diameter
		//double dinG;									//Slipper inner diameter
		//double lSK;										//Distance to the center of mass of piston/slipper assembly from piston head
		///double mK;										//Mass of piston slipper assembly
		//double lKG;										//Lenght of piston surface
		//double lch;										//Lenght piston chamfer
		//double hmin;									//Minimum gap height
		//double lengthB;									//Length cylinder block
		//double lengthcanalB;							//Length canal cylinder block
		//double speedK;									//Relative rotation speed piston
		//Piston macro geometry inputs
		//int PistonMacroGeometry;
		//double rK_red;						//Reduced piston radius for sperical gap
		//double lK_hs;						//Piston cylindrical length for halfspherical gap
		vector<double> polygap_coeff;		//Polynomial coeffcients for gap definition
		vector<double> stepwisegap_d_K;		//Stepwise gap diameter for stepwise-linear gap
		vector<double> stepwisegap_l_K;		//Stepwise gap length for stepwise-linear gap
		//Cylinder macro geometry inputs
		//int CylinderMacroGeometry;
		vector<double> stepwisegap_d_B;		//Stepwise gap diameter for stepwise-linear gap
		vector<double> stepwisegap_l_B;		//Stepwise gap length for stepwise-linear gap
		vector<double> bushingsurface;		//2D surface macrogeometry for bushing.
		vector<double> pistonsurface;		//2D surface macrogeometry for piston.
		unsigned short BushingCirc	;					//number of entries in circumferential directon in McrB file
		unsigned short BushingAx	;					//number of entries in axial direction in McrB file.
		unsigned short PistonCirc	;					//number of entries in circumferential directon in McrK file
		unsigned short PistonAx	;					//number of entries in axial direction in McrK file.
	};
	sGeometry Geometry;


	//---------------------------------------------------------------------------------------------------//
	//------------------------------DECLARATIONS FOR PRESSURE INPUT FILE .PRS----------------------------//
	//---------------------------------------------------------------------------------------------------//
	struct sPressure
	{
		double Vdead;								    //Dead volume displacement chamber
		double alphaD_LP;								//Low pressure flow coefficient
		double alphaD_HP;								//High pressure flow coefficient
		double V_LP;									//Low pressure line volume
		double V_HP;									//High pressure line volume
		double AD_LP;									//Cross section area high pressure
		double AD_HP;									//Cross section area low pressure
		double P1;										//Upstream pressure low pressure throttle
		double P2;										//Downstream pressure high pressure throttle
		unsigned short leakageoption;					//External leakage option
		double QSK;										//Constant leakage flow
	};
	sPressure Pressure;


	//---------------------------------------------------------------------------------------------------//
	//------------------------------DECLARATIONS FOR PRESSURE OUTPUT PFILE .DAT--------------------------//
	//---------------------------------------------------------------------------------------------------//
	struct spFile
	{
		vector<double> time;
		vector<double> pDC;
		vector<double> pHP;
		vector<double> pLP;
	};
	spFile pFile; 

	//added by dwm
	struct sFTGFile
	{
		vector<double> time;
		//vector<double> phi;
		vector<double> FTG;
	};
	sFTGFile FTGFile;

	struct sVPFluxFile
	{
		double qvp;
	};
	sVPFluxFile VPFluxFile;


	//---------------------------------------------------------------------------------------------------------//
	//-----------------------------DECLARATIONS FOR BOUNDARY PARAMETERS INPUT FILE .BPD------------------------//
	//---------------------------------------------------------------------------------------------------------//
	/*struct sBoundary
	{
		double Tmax;							//Max temperature gap
		double Simalphastep;					//Maximum step size
		double Simalphaplot;					//Step size for plotting
		//Parts position vector declaration
		double xA;										//Initial postion piston xA
		double yA;										//Initial position piston yA
		double xB;										//Initial position piston xB
		double yB;										//Initial position piston yB
		//Convection coefficients
		double AlphaDC;									//Convection coeffcient DC	
		double AlphaCase;								//Convection coefficent Case
		//Materials
		double EmodK;								
		double EmodB;	
		double vK;								
		double vB;	
		//Pressure-Deformation Damping Coefficients
		double AlphaP;
		double AlphaDef;									
		double AlphaMu;
		double AlphaTh;
		double Rmin_h;									
		double Rmin_p;
		double epsilonK;
		int nmax;
		int jmax;
		int kmax;
		double delta_v;
		double p_max;
		double Rmin_R;
		double Rmin_E;
		int penCells;						//Number of rows of cells on which to apply contact forces
	};
	sBoundary Boundary;*/


	//--------------------------------------------------------------------------------------------------------//
	//-----------------------------DECLARATIONS FOR OPTIONS PARAMETERS INPUT FILE .OPT------------------------//
	//--------------------------------------------------------------------------------------------------------//
	/*struct sPistonOption
	{
		int ReadpFile;
		int ReynoldsMultiGrid;
		int EnergyEquation;
		int HeatTransfer;
		int PressureDeformation;
		int PressureDeformationOMP;
		int ThermalDeformation;
		int EHDTestRig;
		int TriboTestRig;
	} PistonOptions;*/


	//------------------------------------------------------------------------------------------------------//
	//-----------------------------DECLARATIONS FOR GRIDS PARAMETERS INPUT FILE .GRS------------------------//
	//------------------------------------------------------------------------------------------------------//
	/*struct sPistonGapGrid
	{
		int N;
		int M;
		int Q;
	} PistonGapGrid;

	struct sPistonGapMultiGrid
	{
		int nL,v1,v2,MGInt,VW,Q;
		vector<int> N;
		vector<int> M;
	} PistonGapMultiGrid;*/

	/*struct sSlipperGapGrid
	{
		int M;
		int N;
		int Q;
	} SlipperGapGrid;*/

	/*struct sBlockGapGrid
	{
		int N;
		int M[5];
		int Q;
	} BlockGapGrid;*/


	//--------------FEM PRESSURE MESH-------------//
	//faces coordinates piston and cylinder mesh
	Array<double,2> xyzfK_p;
	Array<double,2> xyzfB_p;
	//nodes coordinates piston and cylinder mesh
	Array<double,2> xyznK_p;
	Array<double,2> xyznB_p;

	//--------------FEM THERMAL MESH-------------//
	//faces coordinates piston and cylinder mesh
	Array<double,2> xyzfK_th;
	Array<double,2> xyzfB_th;
	//nodes coordinates piston and cylinder mesh
	Array<double,2> xyznK_th;
	Array<double,2> xyznB_th;

	//influence matrices
	Array<double,1> IM_piston;
	Array<double,1> IM_cylinder;

	//-----------------------------------------------------------------------------------//
	//-----------------------------DECLARATIONS FOR READING FUNCTIONS--------------------//
	//-----------------------------------------------------------------------------------//
	//void readGeneral(void);
	//void readGeometry(void);
	//void readBoundary(void);
	//void readGrids(void);
	//void readOptions(void);
	void readFTG(void); //dwm
	void readVPFlux(void);
	void readpFile(void);
	void readMacroGeometryPiston(void);
	void readMacroGeometryCylinder(void);

	void readBodySurfacexyzPressure(void);
	void readBodySurfacexyzThermal(string body,string body_path_th);
	void readInfluenceMatricesPistonCylinder(string IM_piston_path,string IM_cylinder_path);

};