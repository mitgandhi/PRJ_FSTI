#include "CGapInput.h"
#include "CGapUtils.h"
#include "sGapResult.h"
#include <CFEMThermal.h>
#include "oil.h"
#pragma once

class CPistonGap
{
	public:


	//Default constructor-destructor
	CPistonGap(void);
	~CPistonGap(void);

	//dwm check load
	vector<int> processorload;
	vector<int> nprocessors;

	//Resume File
	string ResumeFile;
	streampos ResumePosition;
	
	//General variables
	ofstream fout;
	int N,M,Q,counter,nmax,revcounter,revtime;
	double dx,dy,dphi,dAz,rZ,rB,rK,AreaG,AreaK,timeold,timenew,AlphaMu,
		AlphaDef,AlphaP,AlphaTh,R_h,Rold_h,R_p,Rold_p,R_mu,Rold_mu,Rmin_p,
		Rmin_h,Rmin_R,Rmin_E,epsilonK,lvarold,EmodK,EmodB,vK,vB,Eprime;
	int EHDTestRig;
	int TriboTestRig;
	int ReadpFile;
	int PressureDeformation;
	int ReynoldsMultiGrid;
	int EnergyEquation;
	int HeatTransfer;
	int ThermalDeformation;
	int PistonMacroGeometry;
	int CylinderMacroGeometry;
	int PressureDeformationOMP;

	//Processor Load dwm
	FILETIME IdleTime;
	FILETIME KernelTime;
	FILETIME UserTime;

	//Groove LIZHI
	Array<double,1> pgvlizhi;//pressure in the grooves
	Array<int,1> cgvlizhi;//condition of the grooves: cgvlizhi = 1 or 3, 1 = normal, 3 = DC
	Array<double,1> lgvlizhi;//distance between the center of the grooves and DC
	Array<int,1> zerolizhi;//use to be mutipled by the velocity in order to cancel groove fluid velocity
	Array<int,1> slimitlizhi;
	Array<int,1> nlimitlizhi;//slimit and nlimit refer to the lubrication section, not the groove
	int jlizhi;
	int FLAGlizhi;
	//Array<double,1> pplizhi;//Pressure profile
	double strposilizhi;//distance between the first groove to the end of the piston near DC
	double wgvlizhi;//Width of the groove
	double stpgvlizhi;//distance between each groove
	int numgvlizhi;//number of the grooves
	int poblizhi;//groove on piston or bushing? 1: piston linear; 2: bushing linear; 3:piston constant; 4:bushing constant

	//Step LIZHI
	Array<double,1> stepfield;//stepfield show step depth at each positon
	int stepboundary;//interpolate step location into fluid grid

	//PistonGap fields definition
	Array<double,1> p;					//Pressure 2D vector
	Array<double,1> r_p;				//Pressure 2D vector
	Array<double,1> h;					//Height 2D vector
	Array<double,1> h_pre;					//Height 2D vector
	Array<double,1> hT;					//Height 2D vector
	Array<double,1> xm;					//Height 2D vector
	Array<double,1> ym;					//Height 2D vector
	Array<double,1> hold;				//Height 2D vector
	Array<double,1> h1;					//Height 2D vector
	Array<double,1> hK;					//Height 2D vector piston
	Array<double,1> dht;				//Squeeze motion
	Array<double,1> dht_total;				//Squeeze motion
	Array<double,1> sigma;				//Conatct stress
	Array<double,1> T;					//Gap fluid temeprature
	Array<double,1> T_old;					//Gap fluid temeprature
	Array<double,1> T_2d;					//Gap fluid temeprature
	Array<double,1> r_T;				//Pressure 2D vector
	Array<double,1> TK_surf_gap;		//Gap fluid temeprature
	Array<double,1> TB_surf_gap;		//Gap fluid temeprature
	Array<double,1> vx;					//Fluid velocity x
	Array<double,1> vy;					//Fluid velocity y
	Array<double,1> vx_edge;
	Array<double,1> vy_edge;
	Array<double,1> vy_edge_case;	
	Array<double,1> phiD;				//energy dissipation
	Array<double,1> phiD_2d;				//energy dissipation
	Array<double,1> vx_p;				//Fluid velocity x poiseuille
	Array<double,1> vy_p;				//Fluid velocity y poiseuille
	Array<double,1> vx_c;				//Fluid velocity x couette
	Array<double,1> vy_c;				//Fluid velocity y couette
	Array<double,1> dz3;				//Height increment 3D vector
	Array<double,1> dz2;				//Height incrment 2D vector
	Array<double,1> oilviscosity;		//Viscoity field
	Array<double,1> oilviscosity_old;	//Viscoity field previous step
	Array<double,1> mu_n;
	Array<double,1> mu_s;
	Array<double,1> mu_e;
	Array<double,1> mu_w;
	Array<double,1> oildensity;			//Density field
	Array<double,1> oildensity_old;
	Array<double,1> rho_n;
	Array<double,1> rho_s;
	Array<double,1> rho_e;
	Array<double,1> rho_w;
	Array<double,1> dAy;				//Element area y direction
	//Array<double,1> Tfluid;				//Fluid properties temperature field
	
	//PistonReynoldsGS fields definition
	Array<double,1> ploop;		//Convergence loop pressure field
	Array<double,1> ploop_old;		//Convergence loop pressure field
	Array<double,1> pgvold;		//Convergence loop groove pressure
	Array<double,1> pfluid;		//Fluid properties pressure field
	Array<double,1> pold;		//Previuos pressure field 2D vector
	Array<double,1> pnew;		//Actual pressure field 2D vector
	Array<double,1> pEHD;		//EHD Test rig pressure
	Array<double,1> mu;			//2D viscosity field
	Array<double,1> rho2d;		//2D density field --- Lizhi 03/17/2015
	Array<double,1> density_expansion;
	Array<double,1> muT;		//2D viscosity field from friction forces
	Array<double,1> dpx;		//2D circumferential pressure gradient
	Array<double,1> dpy;		//2D axial pressure gradient
	Array<double,1> ap;			//Finite volume coefficient ap
	Array<double,1> an;			//Finite volume coefficent an
	Array<double,1> as;			//Finite volume coeffcient as
	Array<double,1> ae;			//Finite volume coefficent ae
	Array<double,1> aw;			//Finite volume coefficent aw
	Array<double,1> b;			//Reynolds equation right side group
	Array<double,1> As_g;		//groove pressure coeffcient As --- Lizhi 03/19/2015
	Array<double,1> An_g;		//groove pressure coeffcient An
	Array<double,1> Cs_g;		//groove pressure coeffcient Cs
	Array<double,1> Cn_g;		//groove pressure coeffcient Cn
	vector<vector<double>> Bs_g;	//groove pressure coeffcient Bs
	vector<vector<double>> Bn_g;	//groove pressure coeffcient Bn
	Array<double,1> Bsps;
	Array<double,1> Bnpn;

	
	//PistonEnergyGS fields definition
	Array<double,1> Tnew;		//Actual temperature field 3D vector
	Array<double,1> apE;		//Finite volume coefficient apE
	Array<double,1> anE;		//Finte volume coefficent anE
	Array<double,1> asE;		//Finite volume coeffcient asE
	Array<double,1> aeE;		//Finite volume coefficent aeE
	Array<double,1> awE;		//Finite volume coefficent awE
	Array<double,1> atE;		//Finite volume coefficent atE
	Array<double,1> abE;		//Finite volume coefficent abE
	Array<double,1> bE;			//Dissipation term right side of equation
	Array<double,1> dvxz;		//Fluid velocity gradient dxdz
	Array<double,1> dvyz;		//Fluid velocity gradietn dydz
	Array<double,1> Dx;			//Diffusive coefficnet energy equation
	Array<double,1> Dy;			//Diffusive coefficnet energy equation
	Array<double,1> Dz;			//Diffusive coefficnet energy equation
	Array<double,1> Fx;			//Convective coefficnet energy equation
	Array<double,1> Fy;			//Convective coefficnet energy equation
	Array<double,1> Px;			//Peclet number
	Array<double,1> Py;			//Peclet number
	Array<double,1> Ax;			//A coefficient power law scheme
	Array<double,1> Ay;			//A coefficient power law scheme
	
	//PistonForces fields definition
	Array<double,1> dvxT_p;		//velocity gradient friction forces
	Array<double,1> dvyT_p;		//velocity gradient friction forces
	Array<double,1> dvxT_c;		//velocity gradient friction forces
	Array<double,1> dvyT_c;		//velocity gradient friction forces
	Array<double,1> taux;		//tau friction forces
	Array<double,1> tauy;		//tau friction forces
	Array<double,1> FfK;		//fluid force piston gap
	Array<double,1> FfKx;		//fluid force piston gap x
	Array<double,1> FfKy;		//fluid force piston gap y
	Array<double,1> MfKx;		//fluid force moment piston gap x
	Array<double,1> MfKy;		//fluid force moment piston gap y
	Array<double,1> FcK;		//contact force piston gap
	Array<double,1> FcKx;		//contact force piston gap x
	Array<double,1> FcKy;		//contact force piston gap y
	Array<double,1> McKx;		//contact force moment piston gap x
	Array<double,1> McKy;		//contact force moment piston gap y
	Array<double,1> phi;		//angular increment gap
	Array<double,1> zKj;		//moment arm each volume
	
	//PistonThermal fields definition
	Array<double,1> QgapK;		//energy dissipation heat flux piston [W/m2]
	Array<double,1> QgapB;		//energy dissipation heat flux block [W/m2]
	double PhiD_mech;			//total energy dissipation viscous friciton [W] 
	double PhiD_vol;			//total energy dissipation leakage [W] 
	
	//Piston Solids Thermal fields definition
	vector<Array<double,1>> qbi_piston;
	vector<Array<double,1>> qbi_cylinder;
	Array<double,1> EbodyK;			//energy flux to piston body [J/m2]
	Array<double,1> tbodyK_DC;		//length of the period when gap surface is in DC or case volume
	Array<double,1> tbodyK_case;	
	Array<double,1> tbodyB_DC;		
	Array<double,1> tbodyB_case;	
	Array<double,1> EbodyB;			//energy flux to cylinder body [J/m2]
	Array<double,1> EbodyK_old;		//energy flux to piston body previous revolution [J/m2]
	Array<double,1> EbodyB_old;		//energy flux to piston body previous revolution [J/m2]	

	//Solid bodies surface and body temperatures
	Array<double,1> TK_surf;
	Array<double,1> TB_surf;
	Array<double,1> TK_body;
	Array<double,1> TB_body;
	
	//Solids surface deformation fields definition
	Array<double,1> defK_th;		//piston body surface deformation thermal
	Array<double,1> defB_th;		//bushing body surface deformation thermal
	Array<double,1> defK_p;			//piston body surface deformation pressure
	Array<double,1> defB_p;			//bushing body surface deformation thermal
	Array<double,1> defK_th_gap;	//piston gap surface deformation thermal
	Array<double,1> defB_th_gap;	//bushing gap surface deformation thermal
	Array<double,1> defK_p_gap;		//piston gap surface deformation pressure
	Array<double,1> defB_p_gap;		//bushing gap surface deformation pressure
	Array<double,1> defK_gap;			//piston gap surface deformation total
	Array<double,1> defB_gap;			//bushing gap surface deformation total
	Array<double,1> defK_p_gap_old;		//piston gap surface deformation total
	Array<double,1> defB_p_gap_old;		//bushing gap surface deformation total
	Array<double,1> defK_p_gap_squeeze;
	Array<double,1> defB_p_gap_squeeze;

	//Solids gap surface macro-geometry fields definition
	Array<double,1> McrK;	//piston gap surface macro geomnetry
	Array<double,1> McrB;	//cylinder gap surface macro geomnetry


	//-------------------FLUID FVM MESH------------------//
	Array<double,2> xyzf_gap;


	//--------------------------FEM PRESSURE MESH------------------//
	//axis coordinate surface faces
	Array<double,1> zfK_p;
	Array<double,1> zfB_p;
	//node ids neighbours structure to fluid interpolation - body nodes dispalcement to fluid face centers
	Array<int,2> NodeIdK_s2f_p;
	Array<int,2> NodeIdB_s2f_p;
	Array<double,2> NodeDistK_s2f_p;
	Array<double,2> NodeDistB_s2f_p;
	//face ids neighbours fluid to structure interpolation - fluid face centers pressure to body faces
	Array<int,2> FaceIdK_f2s_p;
	Array<int,2> FaceIdB_f2s_p;
	Array<double,2> FaceDistK_f2s_p;
	Array<double,2> FaceDistB_f2s_p;


	//--------------------------FEM THERMAL MESH------------------//
	//axis coordinate surface faces
	Array<double,1> zfK_th;
	Array<double,1> zfB_th;
	//face ids neighbours structure to fluid interpolation - body face centers temperature to fluid face centers
	Array<int,2> FaceIdK_s2f_th;
	Array<int,2> FaceIdB_s2f_th;
	Array<double,2> FaceDistK_s2f_th;
	Array<double,2> FaceDistB_s2f_th;
	//node ids neighbours structure to fluid interpolation - body nodes displacement to fluid face centers
	Array<int,2> NodeIdK_s2f_th;
	Array<int,2> NodeIdB_s2f_th;
	Array<double,2> NodeDistK_s2f_th;
	Array<double,2> NodeDistB_s2f_th;
	//face ids neighbours fluid to structure interpolation - fluid face centers heat flux to body face centers
	Array<int,2> FaceIdK_f2s_th;
	Array<int,2> FaceIdB_f2s_th;
	Array<double,2> FaceDistK_f2s_th;
	Array<double,2> FaceDistB_f2s_th;
	//face ids structure to structure
	Array<int,2> FaceId_K2B_th;
	Array<int,2> FaceId_B2K_th;
	Array<double,2> FaceDist_K2B_th;
	Array<double,2> FaceDist_B2K_th;
	//
	vector<vector<int>> Mid_K2B;
	vector<vector<int>> Mid_B2K;
	vector<vector<double>> Mwt_K2B;
	vector<vector<double>> Mwt_B2K;


	//Convergence dwm
	Array<double,1> pcon;

	//HD dwm
	int levelrequest;
	vector<bool> refinememory;
	int refinecounter;
	vector<int> boundarynodes;

	//----------------------DECLARATION OF STRUCTURES-------------------------//
	struct sGridPistonGap
	{
		int N;			//Number of elements in gap circumferential direction
		int M;			//Number of elements in gap direction
		int Q;			//Number of elements in gap height
	};
	sGridPistonGap gridpistongap;

	struct sGeometryPistonGap
	{
		int npistons;			//Number of pistons
		double dK;				//Piston diameter [mm]
		double doutB;			//Block outer diamater [mm]
		double dinB;			//Block inner diamater [mm]
		double lK;				//Piston length [mm]
		double lch;				//Piston chamfer length [mm]
		double dB;				//Cylinder block pitch diameter [mm]
		double dZ;				//Bushing diameter [mm]
		double lF;				//Bushing length	[mm]
		double lSK;				//Distance piston center of mass/slipper assembly [mm] 
		double le;				//Bushing beginning position [mm]
		double lZ0;				//Length displacement chamber ODC [mm]
		double lKG;				//Length pisto surface [mm]
		double dDK;				//Diameter orifice piston head [mm]
		double lengthB;			//Length cylinder block [mm]
		double lengthcanalB;	//Length canal cylinder block [mm]
		double mK;				//Mass piston/slipper assembly [g]
		double hmin;			//Minimum gap height [microns]
		double doutG;			//External diameter slipper [mm]
		double dinG;			//Internal diameter slipper	[mm]
		//------Derived variables--------//
		double lvar;			//Variable gap length due to piston movement [m]
		double lout;			//Variable distance between piston head and beginning edge of cylinder block [m]
		double zRK;				//Variable distance between piston head and far end of the gap [m]
		double lA;				//Eventual portion of the piston exceeding gap in displacment chamber [m]
		double lB;				//Eventual portion of the cylinder exceeding gap in displacment chamber [m]
		double pos_A;			//Position of the front flat piston edge from the bottom of the cylinder [m]
		double pos_B;			//Position of the back flat piston edge from beginning of bushing [m]
		double delta_A;			//Distance from front flat piston surface to bushing end when piston is at ODC [m]
		double delta_B;			//Distance from piston head center to cylinder block beginning edge when piston is at ODC [m]
		double delta_0;			//Distance from cylinder block origin to cylinder block beginning edge [m]
		double sK;				//Piston stroke [m]
		//------Piston macro geometry variables-----//
		double rK_red;					//Reduced piston radius for sperical gap [microns]
		double lK_hs;					//Piston cylindrical length for halfspherical gap [mm]
		//-------Moment arms for cylinder block balance---//
		double z_A;						//Distance from piston DC coordinate system origin to cylinder block origin [m]
		double z_B;						//Distance from piston Case coordinate system origin to cylinder block origin [m]
		double z_A_old;						//Distance from piston DC coordinate system origin to cylinder block origin [m]
		double z_B_old;						//Distance from piston Case coordinate system origin to cylinder block origin [m]
	};
	sGeometryPistonGap geometrypistongap;

	struct sOperatingPistonGap
	{
		int mode;				//Working mode
		int nrevolutions;		//Number of revolutions
		double speed;			//Pump rotational speed	[rpm]
		double omega;			//Pump rotational speed [rad/s]
		double vK;				//Piston axial velocity [m/s]
		double betamax_deg;		//Maximum displacement angle [deg]
		double betamax_rad;		//Maximum displacement angle [rad]
		double beta_deg;		//Swashplate angle [deg]
		double beta_rad;		//Swashplate angle [rad]
		double gamma_deg;		//Swashplate cross angle [deg]
		double gamma_rad;		//Swashplate cross angle [rad]
		double phi_rad;			//Piston angular position [rad]
		double phi_deg;			//Piston angular position [deg]
		double pHP;				//High pressure [bar]
		double pLP;				//Low pressure	[bar]
		double pDC;				//Pressure Dicplacement Chamber [Pa]
		double pCase;			//Case pressure	[bar]
		double speedK;			//Piston relative rotational speed
		double QSK;				//Piston gap leakage [m^3/s]
		double QSK_p;			//Piston gap leakage poiseuille[m^3/s]
		double QSK_c;			//Piston gap leakage couette [m^3/s]
		
	};
	sOperatingPistonGap operatingpistongap;

	struct sTemperaturePistonGap
	{
		double Tmax;		//Maximum temperature in the gap [°C]
		double THP;		//Temperature high pressure port [°C]
		double TLP;		//Temperature low pressure port [°C]
		double TCase;		//Leakage temperature [°C]
	};
	sTemperaturePistonGap temperaturepistongap;
	
	struct sOilPistonGap
	{
		//int oiltype;
		/*double oildensityaverage;
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
		double oillambda;		//Oil thermal conductivity
		double oilC;			//Oil heat capacity*/
		double AlphaDC;			//Convection coeffcient displacmetn chamber
		double AlphaCase;		//Convection coefficient case
	};
	sOilPistonGap oilpistongap;

	

	struct sForcesPistonGap
	{
		//----------Fluid Forces------------//
		double FfKx;			//Fluid force piston x direction
		double FfKy;			//Fluid force piston y direction
		double MfKx;			//Fluid moment piston x direction
		double MfKy;			//Fluid moment piston y direction
		double FfK;				//Fluid force piston
		double MfK;				//Fluid moment piston
		vector<double> F_fluid;
		//----------Thermal Wedge Fluid Forces------------//
		double FfK_TWx;			//Fluid force piston x direction
		double FfK_TWy;			//Fluid force piston y direction
		double MfK_TWx;			//Fluid moment piston x direction
		double MfK_TWy;			//Fluid moment piston y direction
		double FfK_TW;			//Fluid force piston
		double MfK_TW;			//Fluid moment piston
		//----------Friction Forces------------//
		Array<double,1> FTKx;			//Friction force piston circumferential direction
		Array<double,1> FTKy;			//Friction force piston axial direction
		Array<double,1> FTbx;			//Friction force bushing circumferential direction
		Array<double,1> FTby;			//Friction force bushing circumferential direction
		Array<double,1> FTKx_p;			//Friction force piston circumferential direction poiseuille
		Array<double,1> FTKy_p;			//Friction force piston axial direction poiseuille
		Array<double,1> FTbx_p;			//Friction force bushing circumferential direction poiseuille
		Array<double,1> FTby_p;			//Friction force bushing circumferential direction poiseuille
		Array<double,1> FTKx_c;			//Friction force piston circumferential direction couette
		Array<double,1> FTKy_c;			//Friction force piston axial direction couette
		Array<double,1> FTbx_c;			//Friction force bushing circumferential direction couette
		Array<double,1> FTby_c;			//Friction force bushing circumferential direction couette
		double FTG;				//Friction force from slipper
		double MTK;				//Axial friction force moment
		double MTKtan;			//Circumferential friction force moment
		//----------External Forces------------//
		double FKx;				//External force piston x direction
		double FKy;				//External force piston y direction
		double FK;				//Total external force piston
		double MKx;				//Moment from external force piston around x axis
		double MKy;				//Moment from external force piston around y axis
		double MK;				//Total moment from external force piston
		double Fsk;				//Piston normal force with swashplate
		vector<double> F_external;
		//----------Contact Forces------------//
		double FcKx;			//Contact force piston x direction
		double FcKy;			//Contact force piston y direction
		double McKx;			//Contact moment piston x direction
		double McKy;			//Contact moment piston y direction
		double FcK;				//Contact force piston
		double McK;				//Contact moment piston
		vector<double> F_contact;
		//----------Force Balance----------//
		vector<double> dF;
	};
	sForcesPistonGap forcespistongap;

	//------------Reynolds Multigrid level variables-----------//
	typedef struct
	{
		//Mesh size
		int N,M,nL,v1,v2,MGInt,VW;
		//Geometric increments
		double dx,dy,dphi,Rx;
		//Field coordinates
		Array<double,2> xyz;
		//Centroid spatial coefficients
		Array<double,1> ap;
		Array<double,1> an;
		Array<double,1> as;
		Array<double,1> ae;
		Array<double,1> aw;
		Array<double,1> b;
		Array<double,1> phi;
		//Residual
		Array<double,1> r;
		//Error
		Array<double,1> e;
		//Scalar variable
		Array<double,1> x;
		//Film thickness
		Array<double,1> h;
		//Fluid viscosity
		Array<double,1> mu;
		//Fluid density
		Array<double,1> rho;
		//FaceId prolongation
		Array<int,2> FaceId_c2f;
		Array<double,2> K_bilin;
		//FaceId and distance restriction
		Array<int,2> FaceId_f2c;
		Array<double,2> FaceDist_f2c;
	} LEVEL;
	vector<LEVEL> levels;
	Array<double,1> xold;

	//------------------HD Pressure Calculation-------------------------------//
	typedef struct
	{
		vector<double> GridD;	//FSTI fluid density [kg/m3]
		vector<double> GridH;	//FSTI film thickness of each grid point [m]
		vector<double> GridP;	//FSTI fluid pressure of each grid point [Pa]
		vector<double> GridV;	//FSTI fluid viscosity [Pa-s]
		vector<double> GridX;	//X coordinate of each grid point. Axial from DC [m]
		vector<double> GridY;	//Y coordinate of each grid point. Circumferentially around bushing, same direction as FSTI Phi_k [m]
		vector<double> GridQ;	//Squeeze [m/s]
		vector<double> GridA;	//Expansion
		vector<int> GridE;		//Number of the grid point neighboring to the East (axially positive)
		vector<int> GridN;		//Number of the grid point neighboring to the North (circumferentially positive)
		vector<int> GridS;		//Number of the grid point neighboring to the South (circumferentially negative)
		vector<int> GridW;		//Number of the grid point neighboring to the West (axially negative)
		//vector<int> GridC;		//Section ID number
		//vector<vector<int > > GridT; //Sections under each section
	} HDLEVEL;
	vector<HDLEVEL> HDlevels;
	vector<int> reshere;			//integer for each grid point indicating which level contains pressure values (coarser points in same location contain residual values).
	vector<double> updatepress;		//intermediate storage for new pressure values during calculations


	//------------------DECLARATION OF MULTIGRID FUNCTIONS--------------------//
	void PistonSetMG();
	void PistonGaussSeidelMG(int l,int leg);
	Array<double,1> PistonRestrictionMG(Array<double,1> field_f,Array<int,2> FaceId_f2c,Array<double,2> FaceDist_f2c);
	Array<double,1> PistonProlongationMG(Array<double,1> field_c,Array<int,2> FaceId_c2f);
	Array<double,1> PistonProlongationBilinearMG(int l_c,Array<double,1> field_c,Array<int,2> FaceId_c2f);
	void PistonVcycleMG(void);
	void PistonWcycleMG(void);
	void PistonReynoldsCalcCoefficientsMG(double dt);
	void PistonReynoldsMG(void);
	void PistonReynoldsCalcResidualMG(int l);
	

	//----------------------DECLARATION OF GAP FUNCTIONS----------------------//
	void PistonGap(vector<double> &xp,vector<double> &vp,double dt);
	//Main function to call Piston gap calculations
	void PistonGapJacobian(vector<double> &xp,vector<double> &vp,vector<double> &dF,double dt);
	//Calculate the shaft angle based on the time input
	void PistonCalcPhi(double time);
	//Set DC pressure for piston simulation ideally or reading from .dat file given by pressure module simulation
	void PistonSetPressure(void);
	//Get ideal ramp of DC pressure based on pHP and pLP
	vector<double> PistonGetIdealPressure(double deltaangle);
	//Get pressure from .dat file given by pressure module simulation
	vector<double> PistonGetFilePressure(double deltaangle);
	//Calculate intermediate variables for piston simulation
	void PistonInitializeVariables(double time);
	//Calculate piston axial velocity
	void PistonCalcvK(void);
	//Calculate piston stroke	
	void PistonCalcsK(void);
	//Initialize a 2D pressure vector
	void PistonInitializePressure(void);				
	//Initialize gap 3D temperature vector
	void PistonInitializeTemperature(void);
	//Calculate variable gap length
	void PistonCalclvar(void);				
	//Calculate fluid density as a function of pressure and viscosity in the gap
	void PistonCalcDensity(void);	
	//Initialize fluid vicosity as a function of pressure and visocsity in the gap
	void PistonInitializeViscosity(double plim,double Tlim);	
	//Calculate fluid vicosity as a function of pressure and visocsity in the gap
	void PistonCalcViscosity(double plim,double Tlim);	
	//Guess a fluid viscosity for initialization and slipper approximated friction force
	void PistonGuessViscosity(double &oilviscosityguess);
	//Calculate fluid velocity in x direction
	void PistonCalcVedge(void);	
	void PistonCalcV(void);	
	//Calculate piston gap height
	void PistonCalch(vector<double> &xp,vector<double>&vp);
	void PistonCalch_FSI(vector<double> &xp,vector<double>&vp);
	void PistonCalch_pre(void);
	//Calculate piston gap height
	void PistonCalcRigidh(vector<double> &xp);
	//Calculate piston gap height
	void PistonCalcdht(vector<double> xp,vector<double> &vp);	
	//Calculate piston macro geometry gap surface
	void PistonCalcMcrK(void);
	//Calculate cyinder macro geometry gap surface
	void PistonCalcMcrB(void);
	//Solve piston gap Reynolds Equation
	void PistonReynoldsGS(void);
	void Averagepressure(void);
	//Solve piston gap Reynolds Equation
	void PistonReynoldsCalcCoefficientsGS(double dt);
	//Solve piston gap Reynolds Equation
	void PistonReynoldsCalcResidualGS(double &R);
	//Calculate residual based on fluid pressure
	int PistonCalcLoopResidual(void); 
	//Solve piston gap Energy Equation
	void PistonEnergyGS(double dt); 
	void PistonEnergyGSedge(double dt); 
	void Pistonenergydissipation(void); 
	//Solve piston gap Energy Equation
	void PistonEnergyCalcResidualGS(double &R);
	//Calculate piston leakage in axial direction
	void PistonCalcLeakage(void);
	//Calculate groove stuff outside the loop
	void CPistonGap::Pistongroove(void);
	//Calculate step location
	void CPistonGap::Pistonsteplocation(void);
	//Refine the HD Pressure Grid
	bool CPistonGap::RefineMesh(int gridref,bool forcebalance);
	//Check for duplicate HD mesh points
	int CPistonGap::CheckPoint(double nX,double nY, double tol, int gridref);
	//Updates the pressure field using the Reynolds Equation
	void CPistonGap::SolveReynoldsHD(int gridref,int maxlev,double alphapress,double dt);
	void CPistonGap::ReynoldsResidualHD(int gridref, int maxlev,double dt);
	double CPistonGap::CalcRes(vector<double> pold,vector<double> pnew);
	void CPistonGap::Recursion(int gridref, int maxlev, double dt,int section,bool updown,double relax);
	void CPistonGap::CalcPressureHD(double dt,bool forcebalance);


	//dwm check load
	//Number of parallel processes
	void CPistonGap::IMParallel(int i,int load,int n);



	//-------------DECLARATION OF FORCE CALCULATION FUNCTIONS----------------//
	void PistonCalcFluidForces(void);			//Calculate fluid forces
	void PistonCalcContactForces(void);			//Calculate piston contact forces
	void PistonCalcFrictionForces(void);		//Calculate piston friction forces
	void PistonRelaxFields(void);				//Relax viscosity and thickness fields to smooth friction
	void PistonGuessFTG(void);					//Calculate an approximated friction force from the slipper						
	void PistonReadFTG(void);
	void PistonCalcExternalForces(void);		//Calculate piston external forces	
	void PistonCalcdF(vector<double> &dF);		//Calculate piston force balance considering all the forces
	

	//-------------DECLARATION OF MESH READING AND SIZING FIELDS FUNCTIONS----------------//
	void SizeFieldsPiston(void);
	void SizeFieldsCylinder(void);


	//----------------------DECLARATION OF THERMAL FUNCTIONS-------------------------//
	void PistonCylinderCalcGapThermalFlux(void);						//Gap Thermal heat fluxes
	void PistonCylinderCalcBodyThermalFlux(void);						//Gap Thermal heat fluxes calculation to solid bodies
	void PistonCylinderCalcGapThermalFlux_new(void);						//Gap Thermal heat fluxes
	void PistonCylinderCalcBodyThermalFlux_new(void);						//Gap Thermal heat fluxes calculation to solid bodies
	void PistonCylinderConstructWeightedMatrix(void);
	void PistonCylinderBlockFlux(vector<Array<double,1>> &qbi);			//Cylinder block bore fluxes
	void PistonCylinderSolveBodyThermal(void);							//Solve bodeis thermal problem
	void PistonCylinderSolveBodyThermal_new(void);							//Solve bodeis thermal problem
	void PistonCylinderSolveBodyThermal_iter(void);							//Solve bodeis thermal problem
	void Heatfluxrecovering(int option);

	
	//----------------------DECLARATION OF GAP INTERPOLATION FUNCTIONS-------------------------//
	void PistonCylinderInterpolateSurfaceTemperatures(void);
	void PistonCylinderInterpolatePressureDeformations(void);
	void PistonCylinderSmoothField(Array<double,1> &field, int iter);
	void PistonCylinderStraightenField(Array<double,1> &field);
	void PistonCylinderInterpolateThermalDeformations(void);


	//----------------------DECLARATION OF FEM ANALYSIS FUNCTIONS-------------------------//
	void PistonFEMPressureSurfaceDeformation(void);				//Calculate pressure deformation from influence matrices
	void PistonCylinderFEMPressureSurfaceDeformation(void);		//Calculate pressure deformation from influence matrices


	//----------------------DECLARATION OF UTILITY FUNCTIONS-------------------------//
	void PistonCylinderSearchNodesNeighbours(int nb);			//Search neighbouring nodes for piston cylinder interface interpolation
	void PistonCylinderGapSurfaceCoordinates(void);				//Define surface coordinates for gap,piston and cylinder


	//-----------------------CALCULATE EHD PRESSURE---------------------------//
	void PistonCalcEHDTestRigPressureField(void);


	//-----------------------CALCULATE MESH MOTION FOR .VTK FILE SEQUENCE------------------------//
	void PistonMeshMotion(void);

	oil* my_oil;
	

	

};