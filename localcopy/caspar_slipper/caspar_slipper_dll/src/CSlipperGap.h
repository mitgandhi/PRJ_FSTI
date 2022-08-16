#pragma once
#include "CGapLog.h"
#include "../../caspar_input/caspar_input.h"
#include <blitz/array.h>
#include <ANN.h>
#include <gmm/gmm.h>
#include "CGapResult.h"
#include "oil.h"

//Logging object
extern CGapLog GapLog;

struct ANN
{
	// ANN stuff
	int	dim;							// problem dimension
	double eps;						// tolerance
	ANNpointArray points;	// elements centers
	ANNkd_tree* kdtree;		// kd tree

	ANN()
	{	
		eps = 0;
		dim = 2;
		kdtree = NULL;
	}
};

class CSlipperGap
{

public:

	//Gap inputs
	caspar_input * gapinput;

	CSlipperGap(caspar_input * gapinputs, CGapResult &Result);
	~CSlipperGap(void);

	CGapResult* GapResult;

//--------------------DECLARATION OF VARIABLES-----------------------------//
	
	//General variables
	Range all;

	// Old gap position / velocity vector
	// Needed for the implicit method
	vector<double> pold;
	vector<double> vold;

	//Current position / velocity vectors
	vector<double> control_point_position;
	vector<double> control_point_velocity;

	//what was the FullLoop mode SlipperGap() was called with
	int curFullLoop;

//----------------------DECLARATION OF STRUCTURES-------------------------//
	struct sFluid
	{
		//Fluid Grid Size
		int M;			//Number of elements in gap radial direction (r) (M)
		int N;			//Number of elements in gap circumferential direction (theta) (N)
		int Q;			//Number of elements in gap height direction (z) (Q)

		//In the case of non-uniform radial grids
		int Groove1Location;
		int Groove2Location;
		double Groove1r;
		double Groove2r;
		double Groove1dr;
		double Groove2dr;
		
		//Fluid Grid Spacing
		Array<double,2> dr;
		Array<double,2> dtheta;

		//Finite volume area
		Array<double,2> dA;

		//Rigid Gap Height 2D [m]
		Array<double,2> hrigid;

		//Deformed Gap Height 2D [m]
		Array<double,2> h;

		Array<double,2> hgroove;

		//3D Gap Height Position Vector
		Array<double,3> z;
		
		//Squeeze Velocity 2D [m/s]
		Array<double,2> dht;

		//Old ehdsqz
		Array<double,2> ehdsqzOld;


		//Local coordinate system (LCS)
		//The origin is located at the slipper center and slipper gap surface
		//The orientation is always such that
			//x-axis points in slipper tangental travel direction
			//y-axis points radially outward
			//z-axis Right Hand Rule

			//Polar Coordinates of Fluid Grid [m, rad]
			Array<double,2> r;
			Array<double,2> theta;

			//Cartesian Coordinates of Fluid Grid [m]
			Array<double,2> Lx;
			Array<double,2> Ly;

		//Local rotating coordinate system (LRCS)
		//The origin is located at the slipper center and slipper gap surface
		//The orientation depends on speedK. 
			//if speedk = 1, LRCS will be aligned with the GCS
			//if speedk = 0, LRCS will be aligned with the LCS
			//z-axis Right Hand Rule

			//Cartesian Coordinates
			Array<double,2> LRx;
			Array<double,2> LRy;

		//Global coordinate system (GCS)
		//This origin is located at pump shaft center and swashplate gap surface
		//The orientation is always such that
			//y-axis aligned with LCS y-axis at phi = 0
			//z-axis perpendicular and away from swashplate gap surface
			//x-axis Right Hand Rule
		
			//Cartesian Coordinates of Fluid Grid [m]
			Array<double,2> Gx;
			Array<double,2> Gy;


		//Pressure 2D [Pa]
		Array<double,2> p;

		//Pressure, full loop 2D [Pa] (Used to carry the full FSI field over to the vtk)
		Array<double,2> pFullLoop;

		//Uncut Pressure 2D [Pa]
		Array<double,2> p_uncut;

		//3D Fluid temperature [C] 
		Array<double,3> T;

		//Heat flux [W/m^2]
		Array<double,2> Qflux;

		//Sum of the energy source
		double Qphisum;

		//Old pressure field
		Array<double,2> pold;
				
		//3D Fluid Velocity
		Array<double,3> vr;
		Array<double,3> vtheta;
		
		//3D Fluid Viscosity / Density
		Array<double,3> oilviscosity;
		Array<double,3> oildensity;

		//Contact ?: 0 No (Reynods). 1 Yes (Contact)
		Array<int,2> contact;

		//What is the contact "penetration"
		Array<double,2> h_contact;

		//What is the contact pressure
		Array<double,2> p_contact;

		//What is the height contact should be considered at
		double ContactHeight;

		//Boundary value of the cell
		//-1 Solve reynolds
		//-2 Pocket pressure
		//-3 Case pressure
		//+val [Pa]
		Array<double,2> boundary;

		//A value that varies radially attempting to aprox the IM stiffness of the slipper
		Array<double,1> newt_ehdsqzapprox;

		//A KD tree of the fluid gap to be used with the slipper IMs
		//Note this consideres the rotation of the slipper in GCS, but not rotation about the shaft axia
		ANN KDslip;

		//The alpha value
		double alphaPold;
		double alphaPoldCP;

		Array<double,2> ap;
		Array<double,2> b;

		Array<double,3> ReyVals;
		
	}; sFluid Fluid;

	class sSolid
	{
	public:
		CGapLog* GapLog;

		struct face
		{
			int nodes[3];
			double x;
			double y;
			face(const int n1, const int n2, const int n3, const double X, const double Y)
			{
				nodes[0] = n1;
				nodes[1] = n2;
				nodes[2] = n3;
				x = X;
				y = Y;
			}
		};

		struct node
		{
			double x;
			double y;
			node(const double X, const double Y)
			{
				x = X;
				y = Y;
			}
		};

		string path;
		string suffix;

		vector<face> faces;
		vector<node> nodes;

		vector<double> nodeDeformation;
		vector<double> facePressure;

		int M,N;

		int nodecnt;
		int facecnt;

		//A KD tree of the gap face centroids
		ANN KDfaces;

		double ** IM;
		double * IMpocket;
		double * IMcase;
		double * IMsocket;

		double avgIM;

		Array<double,2> ehd;
		Array<double,2> oldehd;
		Array<double,2> fakeehd;

		vector< vector<int> > quickF2S;

		void loadIM();
		void loadIM_ig1();
		void Tokenize(const string& str,vector<string>& tokens,const string& delimiters)
		{
			 // Skip delimiters at beginning.
			 string::size_type lastPos = str.find_first_not_of(delimiters, 0);
			 // Find first "non-delimiter".
			 string::size_type pos = str.find_first_of(delimiters, lastPos);

			tokens.clear();

			while (string::npos != pos || string::npos != lastPos)
			 {
				  // Found a token, add it to the vector.
				  tokens.push_back(str.substr(lastPos, pos - lastPos));
				  // Skip delimiters.  Note the "not_of"
				  lastPos = str.find_first_not_of(delimiters, pos);
				  // Find next "non-delimiter"
				  pos = str.find_first_of(delimiters, lastPos);
			 }
		}
		void calcDeform(const double pG, const double pcase, const double psocket);

		void writevtk(const string file);
		void writevtk(const string file, const vector<double> &qflux);

		sSolid(CGapLog * gaplog) : GapLog(gaplog)
		{
		}
	};
	
	//A class used for the thermo-elastic solid body
	class sTSolid
	{
	public:
		CGapLog* GapLog;

		struct face
		{
			int nodes[3];
			double x;
			double y;
			face(const int n1, const int n2, const int n3, const double X, const double Y)
			{
				nodes[0] = n1;
				nodes[1] = n2;
				nodes[2] = n3;
				x = X;
				y = Y;
			}
		};

		struct node
		{
			int Gid; //used when referring to the thermal solver

			double x;
			double y;
			node(const double X, const double Y, const int GID)
			{
				x = X;
				y = Y;
				Gid = GID;
			}
		};

		string path;
		string suffix;

		vector<face> faces;
		vector<node> nodes;

		vector<double> nodeDeformation;
		vector<double> nodeTemperature;
		vector<double> nodeTemperature_old_FullMesh;
		vector<double> faceFlux;
		vector<double> faceFlux_old;

		//Because faceFlux needs to be an average over a whole revolution, what is the count
		double faceFluxcnt;

		int M,N;

		int nodecnt;
		int facecnt;

		//A KD tree of the gap face centroids
		ANN KDfaces;

		Array<double,2> temp;
		Array<double,2> deform;
		Array<double,2> pdeform;	//deform is shifted so only positive

		//last step thermal analysis ran
		int last_run_step;

		//Options file
		string option_file;

		//Methods
		void load(const string Option_File);

		sTSolid(CGapLog * gaplog) : GapLog(gaplog)	//constructer
		{
			last_run_step = -1;	//forces thermal analysis to run the first timestep
		};
	};

	//The slipper & swashplate structs for EHD (pressure) deformation
	sSolid slipper;
	sSolid swashplate;

	//The slipper & swashplate structs for thermoelastic analysis
	sTSolid t_slipper;
	sTSolid t_swashplate;

	//A structure containing the 'full' rotating group of slippers
	struct sFullGroup
	{
		//The Gx, Gy, and p at every alpha step
		Array<double,3> Gx;
		Array<double,3> Gy;
		Array<double,3> p;

		//A KD tree of the 9 slipper fluid volume centroids at the current time step
		ANN KDslip;
	}; sFullGroup FullGroup;

	struct sOperatingSlipperGap
	{

		double T_Leak;			//Temperature of leakage flow [C]

		double speed;			//Pump rotational speed	[rpm]
		double omega;			//Pump rotational speed [rad/s]
		double betamax_deg;		//Maximum displacement angle [deg]
		double betamax_rad;		//Maximum displacement angle [rad]
		double beta_deg;		//Swashplate angle [deg]
		double beta_rad;		//Swashplate angle [rad]
		double gamma_rad;		//Swashplate cross angle [rad]
		double phi_rad;			//Slipper angular position [rad]
		double phi_deg;			//Slipper angular position [deg]
		double pHP;				//High pressure [bar]
		double pLP;				//Low pressure	[bar]
		double pDC;				//Pressure Dicplacement Chamber [Pa]
		double pcase;			//Case pressure	[Pa]
		double speedK;			//Piston relative rotational speed

		double pG;				//slippper pocket pressure [Pa]
		double pGprevious;		//previous timestep slippper pocket pressure [Pa]
		double QSG;				//flow from under slipper (calculated through the gap)
		double QSG_orifice;		//flow through slipper orifice (considered the 'proper' slipper leakage)
		double Q_S_pois;		//Poiseuille leakage flow [m^3/s]
		double Q_S_couette;		//Couette leakage flow [m^3/s]
		
		double Ploss;			//Total power loss due to friction & leakage [W]
		double PlossMech;		//Total power loss due to ONLY friction [W]

		//Old model approach - depreciated
		Array<double,2> vgr;			//Slipper radial velocity
		Array<double,2> vgtheta;	//Slipper tangental velocity

		//slipper velocities, useful for plotting
		Array<double,2> vgx;
		Array<double,2> vgy;

		//Slipper surface velocity
		Array<double,2> tvr;
		Array<double,2> tvtheta;

		//Swash/wobble plate velocity
		Array<double,2> bvr;
		Array<double,2> bvtheta;

	};
	sOperatingSlipperGap operatingslippergap;

	struct sGeometrySlipperGap
	{

		double lKG;										//Length of piston gap surface
		double dZ;										//Bore Diameter
		
		int npistons;									//number of pistons / slippers
		double doutG;									//Slipper outer diameter
		double routG;									//Slipper outer radius
		double dinG;									//Slipper inner diameter
		double rinG;									//Slipper inner radius
		double mG;										//Slipper mass
		double mK;										//Piston / Slipper assembly mass
		double lSG;										//Distance from piston head to slipper center of mass
		double lG;										//Distance from piston head to gap area
		double rB;										//Cylinder block pitch radius
		double dDG;										//Slipper orifice diameter
		double lDG;										//Slipper orifice length
		double dK;										//Piston Diameter
		double dDK;										//Orifice piston head diameter
		double lDK;										//Orifice piston head length
		double Fslipper;								//Slipper hold down force
		double vPocket;									//Volume of the slipper pocket
		double hmaxG;									//Maximum height of slipper holder
		double alphaD_KG;								//Orifice coefficient piston-slipper assembly
		int flowtype;									//Type of flow option piston/slipper 0, 1, 2, 3

		Array<double, 2> SlipperMacro;
	};
	sGeometrySlipperGap geometryslippergap;

	struct sForcesSlipperGap
	{
		//----------Fluid Forces------------//
		double FfGz;			//Fluid force slipper z direction [N]
		double MfGx;			//Fluid moment slipper x direction [Nm]
		double MfGy;			//Fluid moment slipper y direction [Nm]
		Array<double,2> F;		//Fluid slipper pressure force array , Z direction, [N]
		vector<double> F_fluid;
		/* 
		//----------Thermal Wedge Fluid Forces------------//
		double FfK_TWx;			//Fluid force piston x direction
		double FfK_TWy;			//Fluid force piston y direction
		double MfK_TWx;			//Fluid moment piston x direction
		double MfK_TWy;			//Fluid moment piston y direction
		double FfK_TW;				//Fluid force piston
		double MfK_TW;				//Fluid moment piston
		*/
		
		//----------Friction Forces------------//
		double M_FTGx;				//Friction moment to slipper
		double M_FTGy;				//Friction moment to slipper

		double F_FTGx;				//Friction force to slipper
		double F_FTGy;				//Friction force to slipper
		double FTG;				//Friction force to slipper

		//----------External Forces------------//
		double MGx_centrifugal; //Centrifugal force in Mx direction
		double FSK;				//Reaction force from piston, z direction
		double FTK;				//Piston friction force
		double psocket;		//Socket pressure	[Pa]
		vector<double> F_external;
		
		//----------Contact Forces------------//
		vector<double> F_contact;

		//----------Force Balance----------//
		vector<double> dF;

		//Coupled caspar - Vector of piston friction forces//
		vector<vector<double> > FTK_vector;
	};
	sForcesSlipperGap forcesslippergap;

	struct sTemperatureSlipperGap
	{
		double T_HP;		//Temperature high pressure port [°C]
		double T_LP;		//Temperature low pressure port [°C]
		double T_Leak;		//Leakage temperature [°C]
	};
	sTemperatureSlipperGap temperatureslippergap;

	/*
	struct sOilSlipperGap
	{
		int oiltype;
		double oildensityaverage;
		double oilviscosityaverage;
		double oilviscosity;	//Oil dynamic viscosity
		double oildensity;		//Oil density 
		double oilbetaP;		//Oil density pressure coefficient
		double oilbetaT;		//Oil desity thermal coefficient
		double alpha1;			//Alpha 1 used for oil
		double alpha2;			//Alpha 2 used for oil
		double alpha3;			//Alpha 3 used for oil
		double oilW;			//Oil kinematic viscosity weightening factor
		double oilTc1;			//Kinematic viscosity temperature coefficient 1
		double oilPc1;			//Kinematic viscosity pressure coefficent 2
		double oilTc2;			//Kinematic viscosity tempearature coefficent 2
		double oilPc2;			//Kinematic viscosity pressure coefficient 2
		double oilK;				//Oil Bulk modulus
		double oilbetaKP;			//Oil bulk modulus pressure coefficient
		double oilbetaKT;			//Oil bulk modulus thermal coefficient
		double oillambda;		//Oil thermal conductivity
		double oilC;			//Oil heat capacity
	};
	sOilSlipperGap oilslippergap;
	*/
	oil* oil_properties;


	//-----------------------DECLARATION OF FUNCTIONS------------------------//
	void SlipperGap(vector<double> &xg,vector<double> &vg,vector<double> &dF, int FullLoop = 0);
	void SlipperStartTimestep(vector<double> &Newtonpositionslipper, vector<double> &Newtonvelocityslipper, vector<double> &oldPosition, vector<double> &oldVelocity);
	void SlipperEndTimestep(void);
	void SlipperStartRevolution(void);
	void SlipperEndRevolution(void);
	void SlipperSetCoordSys(void);
	void SlipperDefinePressureBounds(void);
	int SlipperReynoldsGS(double alphaPold);
	int SlipperReynolds(double alphaPold);
	int SlipperReynoldsPC(double alphaPold);
	int SlipperReynoldsFkSqz(double alphaPold);
	int SlipperReynoldsCut(double alphaPold);
	void alphaP(const double init_residual, const double resid, double & alphaPold);
	double presid();
	double presidPC();
	int SlipperEnergy();
	int SlipperEnergyCG();
	void SlipperCalcFluidForces(void);
	void SlipperCalcExternalForces(void);
	void SlipperCalcViscosity(const double alpha);
	void SlipperCalcDensity(void);
	double CalcK(const double& T_C, const double& p_Pa);
	void SlipperInitializeTemperature(void);
	void SlipperCalcFluidV(void);
	void SlipperCalch(vector<double> &xp);
	void SlipperCalch(void);
	void SlipperCalcdht(void);
	void SlipperCalcdF(vector<double> &dF, const vector<double> &xg);
	void SlipperCalcHolder(const vector<double> &xg, vector<double> &F, const Array<double,2>& A);
	bool SlipperCalcContactPressure(const double alphaEHD, const int FullLoop, const double alphaContactP);
	void SlipperpG(const double alpha);
	void SlipperAssignResults(const vector<double> &xg);
	void GetFSK(void);
	void GetFTK(void);
	void SlipperVTK(const string file);
	void LoadSlipperIM(void);
	void LoadSwashplateIM(void);
	void SlipperInitEHD(void);
	void SlipperCalcEHD(const double alpha);
	void SwashplateInitEHD(const double alpha);
	void SwashplateCalcEHD(const double alpha);
	void SwashplateInitPressure(vector<double> & pDC);
	void SlipperPressureBounds(void);
	TinyVector<double,3> Solve3(Array<double,2> A, TinyVector<double,3> b);
	double det(Array<double,2> m,int i, int j);
	static void Tokenize(const string& str,vector<string>& tokens,const string& delimiters);
	void testDeform(void);
	void testDeform(vector<double> &xg,vector<double> &vg);
	void cell2face(Array<double,2> & p, Array<double,2> & pe, Array<double,2> & pw, Array<double,2> & pn, Array<double,2> & ps);
	void cell2face(Array<double,3> & p, Array<double,3> & pe, Array<double,3> & pw, Array<double,3> & pn, Array<double,3> & ps, Array<double,3> & pt, Array<double,3> & pb);
	Array<double,2> North(Array<double,2> & phi);
	Array<double,2> South(Array<double,2> & phi);
	Array<double,2> East(Array<double,2> & phi);
	Array<double,2> West(Array<double,2> & phi);
	void Fluid2Slipper();
	void Fluid2Swashplate();
	void Slipper2Fluid(const double alpha);
	void Swashplate2Fluid(const double alpha);
	void writevtk(string fileName);
	void testGap(void);
	void debug_makeEHDfluidpressurefile();
	void CalcSlipperThermal(const double alpha = 1.0);
	void CalcSwashplateThermal(const double alpha = 1.0);
	void Fluid2Slipper_Thermal(void);
	void Slipper2Fluid_Thermal(const double alpha = 1.0);
	void Fluid2Swashplate_Thermal(void);
	void Swashplate2Fluid_Thermal(const double alpha = 1.0);
	
};

