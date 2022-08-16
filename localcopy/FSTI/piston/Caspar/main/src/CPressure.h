#include "CGapInput.h"
#pragma once


class CPressure
{

public:
	CPressure(void);
	~CPressure(void);
	void Pressure(void);		//Main pressure module function that calls the integration
	ofstream fout;

private:
 //-------------DECLARATION OF STRUCTURES---------------//
	struct sPump
	{	
		//From General Input File .gen
		int mode;							//Working mode
		int npistons;						//Number of pistons
		int nrevolutions;					//Number of revolutions
		double speed;						//Shaft speed
		double dB;							//External diameter cylinder block
		double dK;							//Piston diameter
		double beta_deg;					//Swash plate angle
		double betamax_deg;					//Max swash plate angle
		double gamma_deg;					//Swash plote cross angle
		double pHP;							//High pressure
		double pLP;							//Low pressure
		double pCase;						//Case pressure
		double THP;						//Temperature in high pressure
		double TLP;						//Temperature at low pressure
		double TCase;						//Leakage temperature
		//From Pressure Input File .prm
		double Vdead;						//Dead volume
		double alphaD_LP;					//Low pressure flow coefficient
		double alphaD_HP;					//HIgh pressure flow coefficent
		double V_LP;						//Volume low pressure
		double V_HP;						//Volume high pressure
		double AD_LP;						//Cross section area low pressure
		double AD_HP;						//Cross section area high pressure
		double P1;							//Upstream pressure low pressure throttle
		double P2;							//Downstream pressure high pressure
		int leakageoption;					//Leakage option
		double QSK;							//Constant leakage value
		//Derived variables
		double beta_rad;					//Swash plate angle
		double betamax_rad;					//Max swash plate angle
		double gamma_rad;					//Swash plate cross angle
		double Ap;							//Piston area [m^2]
		double hK;							//Piston stroke [mm]
		double omega;						//Angular velocity [rad/s]
		double tr;							//Time per revolution
		double Qr_theorSI;					//Theoretical flow rate [m^3/s]
		double Qr_theor;					//Theoretical flow rate [L/min]
		vector<double> arhp;				//Pump opening HP area all pistons one time [mm^2]
		vector<double> arlp;				//Pump opening LP area all pistons one time [mm^2]
		vector<double> phi_rad;				//Pump phi all pistons one time [rad]
		double errorHP;
		double errorLP;
		double errorHP2;
		double errorLP2;
		double errorHPsum;
		double errorLPsum;
		double errorHPsum2;
		double errorLPsum2;
		double pHPsum;
		double pLPsum;
		double AD_HPsum;
		double AD_LPsum;
	};

	struct sPumpVPArea
	{
		vector<double> arphi;			//Pump opening phi [rad]
		vector<double> arhp;			//Pump opening HP area [mm^2]
		vector<double> arlp;			//Pump opening LP area [mm^2]
	};
	struct sOil
	{
		//From General Input File .gen
		int oiltype;				//Type oil chosen
		double oildensity;			//Oil density
		double oilbetaP;			//Oil density pressure coefficient
		double oilbetaT;		    //Oil desity thermal coefficient
		double alpha1;				//Alpha 1 used for oil
		double alpha2;				//Alpha 2 used for oil
		double alpha3;				//Alpha 3 used for oil
		//From Pressure Input File .prm
		double oilK;				//Oil Bulk modulus
		double oilbetaKP;			//Oil bulk modulus pressure coefficient
		double oilbetaKT;			//Oil bulk modulus thermal coefficient
	};

	struct sSimVariables
	{
	vector<double> pD;				//Pressure states vector pD
	vector<double> dpD;				//Differential pressure states vector dpD
	vector<double> QrHPi;			//Flow rate from cylinder to HP for each piston
	vector<double> QrLPi;			//Flow rate from low pressure to one cylinder
	vector<double> QHPitheo;		//Theoretical flow rate HP side each piston
	vector<double> Qrleaki;			//Constant leakage for each piston
	vector<double> QsKi;			//Gap flow between piston and cylinder for each piston
	vector<double> QsGi;			//Gap flow between slipper and swashplate for each piston
	vector<double> QsBi;			//Gap flow between cylinder block and valveplate  for each piston
	vector<double> dVdT;			//Derivation of the volume of the displacement chamber
	double Vo;						//Current volume displacement chamber
	double QrHP;					//Total flow rate HP side
	double QrLP;					//Total flow rate LP side
	double Q2;						//Total flow rate from pump to downstream load
	double Q1;						//Total flow rate from upstream inlet to pump
	double QHPtheo;					//Total theoretical flow rate HP side
	double QLPtheo;					//Total theoretical flow rate LP side
	};

	struct sSim
	{	
		double tStart;					// Simulation start time [s]
		double DegStep;					// Simulation angular step [deg]
		double resume;					// Simulation Resume option 0:no 1:yes [-]
		double tEnd;					// Simulation end time [s];
		double tDelta;					// Time step [s];
		double IndexEnd;				// Simulation time step Index at end of 1 revolution [-]
		int Index;						// Index of time steps in one revolution [-]
		int rev;						// Current revolution during simulation
		int stepcounter;
		int stepcounterint;
		int revcounter;
		int dispcounter;
	};
	struct sResult
	{
		vector<double> t,phi,pD,dPD,QrHPi,QrLPi,Qri,QsK,QsG,QsB,dVdT,Vo,arhp,arlp,pHP,pLP,QrHP,QrLP,QHPtheo,QLPtheo;
	};

	
	//----------------STRUCTURE VARIABLES------------------------
	sPump Pump;
	sOil Oil;
	sPumpVPArea PumpVPArea;
	sSimVariables SV;
	sSim Sim;
	sResult Result;


	//-------------DECLARATION OF STRINGS----------------------
	string Area_file_pump;


	//-----------DECLARATION OF FUNCTIONS--------------------
	void readPressure(void);
	string Input_file;
	//Read pump area from .ar file
	void ReadPumpArea(void);						
	//Calculate the valve plate area depending on the time step
	void CalcVPArea(void);								
	//Initialize states vector elements	
	void InitializeValues(void);
	//Calculate variables in SI units and derived variables
	void CalcIntermediateVariables(void);
	//Determine pistons angular positions
	void CalcPistonPositions(double time);	
	//Calculate theoretical flow rates HP and LP side
	void CalcTheoreticalFlow(void);
	//Calculate the sign of the input a
	int CalcSgn(double a);
	//Calculate fluid density [kg/m^3] based on the input option. Inputs: T_C [C], p_Pa [Pa]
	double CalcDensity(const double& T_C,const double& p_Pa);
	//Calculate fluid bulk modulus [Pa] on the input option. Inputs:  T_C [C], p_Pa [Pa]
	double CalcK(const double& T_C,const double& p_Pa);
	//Main solver function containing the pressure built up
	//Inputs: time [s], states pressure vector y, differential pressure vector yp
	void odef(double time,vector<double> &y,vector<double> &yp);	
																																
	

    //------------DECLARATION OF OFSTREAM----------------
	ofstream output, errorout, fsolutionwrite;


	//------------SOLVER VARIABLES------------------
	int n;				// Dimension of the system
	double x;			// Initial value of dependent variable (usually time)
	double xbeg;		// Beginning x-value
	double xend;		// Final x value (xend-x may be positive or negative).
	double dx;			// time step for intermediate output
	int nrdens;
	// switch for rtol and atol; if itol = 0, rtol and atol are scalars
	// if itol = 1, rtol and atol are vectors
	int itoler;
	double rtoler_s, atoler_s;
	vector<double> rtoler,atoler;
	//double rtoler; //Relative and absolute error tolerances
	//double atoler; //They can be both scalars or vectors of length n 
	//(in the scalar case pass the addresses of	variables 
	//where you have placed the tolerance values). If set as
	//NULL on input, they are set to 1.0e-7 and itoler is set to 0.
	int iout;			// routine for dense output at every time step is called if iout = 1
	double h;			// integration step length
	// Derived variables
	double hmax;		// maximal step size
	double hinit;
	int nmax;			// maximal number of steps
	double uround;		// smallest number satisfying 1.0 + uround > 1.0
	double safe;		// safety factor in step size prediction
	double facl;		// facl, facr--parameters for step size selection
	double facr;		// facl, facr--parameters for step size selection
	double beta;		// for stabilized step-size control
	int nstiff;			// test for stiffness
	int meth;			// switch for the choice of the coefficients
	int *indir;			// array used for dense output
	unsigned *icont;	// indices for components for which dense output is required
	// number of function evaluations (not counting those in numerical Jacobian calculations)
	int nfcn;
	// number of attempted steps
	int nstep;
	// number of accepted steps
	int naccpt;
	// number of rejected steps
	int nrejct;
	// stores past value of x
	double xold;
	// stores past value of h
	double hold;
	// x at discrete points specified by dx interval
	double xd;
	// vectors used in integration steps
	vector<double> rcont1;
	vector<double> rcont2;
	vector<double> rcont3;
	vector<double> rcont4;
	vector<double> rcont5;
	vector<double> y; //Initial y values.
	vector<double> yp; 
	vector<double> f;
	vector<double> yy1;
	vector<double> k1;
	vector<double> k2;
	vector<double> k3;
	vector<double> k4;
	vector<double> k5;
	vector<double> k6;
	vector<double> ysti;

//------------SOLVER METHODS------------------
	double ODEsign(double a, double b);
	void ODEInitialize();
	void ODEIntegrate();
	double ODEhinit();
	int ODECoreIntegrator();
	double ODEContinuousOutput(unsigned i);
	void ODESolutionWrite();
	
	int ODESolutionOutput();///Function that controls the output of the results. Modify this routine according to your needs
	// get number of function evaluations
	int ODENumFunction() const { return nfcn; }
	// get number of attempted steps
	int ODENumStep() const { return nstep; }
	// get number of accepted steps
	int ODENumAccept() const { return naccpt; }
	// get number of rejected steps
	int ODENumReject() const { return nrejct; }
};


