#include "CGapLog.h"
#include "../../caspar_input/caspar_input.h"
#include <blitz/array.h>
#pragma once

using namespace std;
using namespace blitz;

//Foward decleration required for the plot method to call the slipper vtk method
class CSlipperGap;

class CGapResult
{
public:

	//what directory should the .txt (matlab) files go
	string txt_dir;

	//what directory should the .vtk (paraview) files go
	string vtk_dir;

	//Gap inputs
	caspar_input * gapinput;

	//General results
	double time;					//time [s] - REPLACED WITH GAPINPUT.COMMON
	double phi;						//shaft angle [deg] - REPLACED WITH GAPINPUT.COMMON

	//Pressures
	double pDC;						//Displacement chamber pressure [bar] - REPLACED WITH GAPINPUT.COMMON
	double pHP;						//High pressure port pressure [bar] - REPLACED WITH GAPINPUT.COMMON
	double pLP;						//Low pressure port pressure [bar] - REPLACED WITH GAPINPUT.COMMON
	
	//Positions and velocities
	vector<double> gappositions;	//parts positions vector [m]
	vector<double> gapvelocities;	//parts velocities vector [m/s]
	vector<double> QSKTot;			//Total leakage from pistons [l/min]
	vector<double> MSKTot;			//Total torque loss from pistons [Nm]
	vector<double> MFT;				//Total torque loss due to pistons spin in chamber [Nm]
	vector<double> PhiDK;			//Total energy dissipated [W]



	//Slipper results
	struct sSlipperResult
	{
		Array<double,2> P_slipper;
		Array<double,2> P_ehdsq;
		Array<double,3> T_slipper;
		Array<double,2> h_slipper;
		Array<double,2> slipperEHD;
		Array<double,2> swashplateEHD;
		Array<double,3> vr_slipper;
		Array<double,3> vtheta_slipper;
		Array<double,2> vgx;
		Array<double,2> vgy;


		//----------Fluid Forces------------//
		double FfGz;			//Fluid force slipper z direction [N]
		double MfGx;			//Fluid moment slipper x direction [Nm]
		double MfGy;			//Fluid moment slipper y direction [Nm]
		Array<double,2> F;		//Fluid slipper pressure force array , Z direction, [N]
		vector<double> F_fluid;

		//----------Friction Forces------------//
		double M_FTGx;				//Friction moment to slipper
		double M_FTGy;				//Friction moment to slipper
		double F_FTGx;				//Friction force to slipper
		double F_FTGy;				//Friction force to slipper
		double FTG;					//Friction force to slipper
		double FTK;					//Piston friction force
		
		//----------External Forces------------//
		double MGx_centrifugal; //Centrifugal force in Mx direction
		double FSK;				//Reaction force from piston, z direction
		vector<double> F_external;
		
		//----------Contact Forces------------//
		vector<double> F_contact;	//These forces are not "contact" in the gap, but rather hold down forces acting on the case end of the slipper.
		double avg_p_contact;		//Contact pressure in the gap, averaged over all fluid cells. 

		//----------Force Balance----------//
		vector<double> dF;
		vector<double> dFcomp;

		double minh;		//Minimum gap height
		double meanh;		//Average gap height
		double maxh;		//Maximum gap height

		double pG;			//Pocket pressure [Pa]
		double Ploss;		//Total power loss [W]
		double PlossMech;	//Mechanical power loss [W]
		double QSG;			//External leakage [m^3/s]
		double TorqueLoss;	//Torque loss [N*m]
		double xg0;			
		double xg1;	
		double xg2;	

		//socket friction
		double Msock_maxfric;
		double Msock_sum;
		bool Sock_fixed;
		double M_TJx;		//x-friction moment in the ball joint [N*m]
		double M_TJy;		//y-friction moment in the ball joint [N*m]
		double visc_sock;	//estimated viscosity of oil in the socket [Pa*s]
		double tilt_speed;	//rotational velocity of socket about ball joint [rps]
		double psocket;		//estimated pressure in the socket [Pa]
		double mu_coef;		//friction coefficient [-]
		double F_fric;		//Force of friction acting between socket and ball joint [N]

		//Leakage components
		double Q_S_pois;		//Poiseuille leakage flow [m^3/s]
		double Q_S_couette;		//Couette leakage flow [m^3/s]
		
	};
	sSlipperResult SlipperResult;

	//Methods
	void plot1d(void);
	void plot2d(CSlipperGap * SlipperGap);
	void make_summed(vector<double> &sums);
	CGapResult(caspar_input * gapinputs);	
};
