#include "CPressure.h"
#include <iostream>
#include <iomanip>
#include "../../caspar_input/input.h"
#pragma once


extern class CGapInput myGapInput;

extern class input myinput;

CPressure::CPressure()//Constructor
{

	readPressure();

	Pump.mode = myinput.data.operating_conditions.mode;
	Pump.npistons = myinput.data.operating_conditions.npistons;
	Pump.nrevolutions = myinput.data.lubrication_module.n_lubrication_revolutions;//ONLY FOR INTERNAL USE Coupled Caspar will have its own, separate pressure module!
	Pump.speed = myinput.data.operating_conditions.speed * PI/30;
	Pump.dB = myinput.data.geometry.dB;
	Pump.dK = myinput.data.geometry.dK;
	Pump.beta_deg = myinput.data.operating_conditions.beta * (180/PI);
	Pump.betamax_deg = myinput.data.operating_conditions.betamax * (180/PI);
	Pump.gamma_deg = myinput.data.geometry.gamma * (180/PI);
	Pump.pHP = myinput.data.operating_conditions.HP;
	Pump.pLP = myinput.data.operating_conditions.LP;
	Pump.pCase = myinput.data.operating_conditions.pCase;
	Pump.THP = myinput.data.operating_conditions.T_HP;
	Pump.TLP = myinput.data.operating_conditions.T_LP;
	Pump.TCase = myinput.data.operating_conditions.T_Leak;

	Pump.Vdead = myGapInput.Pressure.Vdead;
	Pump.alphaD_LP = myGapInput.Pressure.alphaD_LP;
	Pump.alphaD_HP = myGapInput.Pressure.alphaD_HP;
	Pump.V_LP = myGapInput.Pressure.V_LP;
	Pump.V_HP = myGapInput.Pressure.V_HP;
	Pump.AD_LP = myGapInput.Pressure.AD_LP;
	Pump.AD_HP = myGapInput.Pressure.AD_HP;
	Pump.P1 = myGapInput.Pressure.P1;
	Pump.P2 = myGapInput.Pressure.P2;
	Pump.leakageoption = myGapInput.Pressure.leakageoption;
	Pump.QSK = myGapInput.Pressure.QSK;

	//Oil.oiltype = myinput.data.oil.oiltype;
	//Oil.oildensity = myinput.data.oil.oildensity;
	//Oil.oilbetaP = myinput.data.oil.oilbetaP;
	//Oil.oilbetaT = myinput.data.oil.oilbetaT;
	//Oil.oilK = myinput.data.oil.oilK;
	//Oil.oilbetaKP = myinput.data.oil.oilbetaKP;
	//Oil.oilbetaKT = myinput.data.oil.oilbetaKT;
	//Oil.alpha1 = myinput.data.oil.alpha1;
	//Oil.alpha2 = myinput.data.oil.alpha2;
	//Oil.alpha3 = myinput.data.oil.alpha3;

	Pump.errorHP = 0.0;
	Pump.errorLP = 0.0;
	Pump.errorHP2 = 0.0;
	Pump.errorLP2 = 0.0;
	Pump.errorHPsum = 0.0;
	Pump.errorLPsum = 0.0;
	Pump.errorHPsum2 = 0.0;
	Pump.errorLPsum2 = 0.0;
	Pump.pHPsum = 0.0;
	Pump.pLPsum = 0.0;
	Pump.AD_HPsum = 0.0;
	Pump.AD_LPsum = 0.0;
	Sim.stepcounter = 0;
	Sim.stepcounterint = 0;
	Sim.revcounter = 0;
	Sim.dispcounter = 0;



};
CPressure::~CPressure()//Desctructor
{


};
void CPressure::Pressure(void)
{
	//Reading Pump Valve Plate Opening Area
	ReadPumpArea();
	
	//Adjusting the values to SI units
	CalcIntermediateVariables();
	
	//Define vectors sizes
	Pump.phi_rad.resize(Pump.npistons);
	Pump.arhp.resize(Pump.npistons);
	Pump.arlp.resize(Pump.npistons);

	SV.dVdT.resize(Pump.npistons);
	SV.QrHPi.resize(Pump.npistons);
	SV.QrLPi.resize(Pump.npistons);
	SV.QHPitheo.resize(Pump.npistons);
	SV.Qrleaki.resize(Pump.npistons);
	SV.QsKi.resize(Pump.npistons);
	SV.QsGi.resize(Pump.npistons);
	SV.QsBi.resize(Pump.npistons);
	
	// Total number of states to be calculated for solver
	n = Pump.npistons+2;
	// Vector containing change in value of all states
	SV.dpD.resize(n);	
	// Vector containing value of all states
	y.resize(n);					

	//Initialize the states of vector y
	InitializeValues();
	//Calculate initial piston positions
	CalcPistonPositions(0);
	//Initialize time to zero when writing first line of Results.time vector
	x = 0;
	
	//Solver variables initializations
	//Setting Relative and Absolute tolerances
	itoler = 1;
	rtoler_s = 1;
	atoler_s = 1;
	rtoler.resize(n);
	atoler.resize(n);
	for(int z=0;z<n;z++)
	{
		rtoler[z] = 1e-6;
		atoler[z] = 1e-6;
	}
	//Initial step size
	hinit = 0.0;
	//Use SolutionOutput routine
	iout = 2;
	safe=0.0;
	beta= 0.0;
	nstiff =-1;
	nrdens= n;
	icont = NULL;
	uround = 0.0;
	nstep=0;
	naccpt=0;
	nrejct=0;
	//Use default values (see header files) for these parameters:
	hmax= 0.0;
	nmax= 0;
	meth = 1;
	facl= 0.0;
	facr= 0.0;
	nfcn=0;
	h =0.0;
	hold=hinit;
	//====sizing vectors used ============
	yp.resize(n);
	f.resize(n);
	yy1.resize(n);
	k1.resize(n);
	k2.resize(n);
	k3.resize(n);
	k4.resize(n);
	k5.resize(n);
	k6.resize(n);
	ysti.resize(n);
	//==============================================
	// initial value for x
	xbeg = 0;
	// final value for x
	xend = Sim.tEnd;
	// interval of x for printing output
	dx = Sim.tDelta;
	xold=xbeg;
	x=xbeg;
	xd=xbeg;
	//====ODE solver ============
	ODEInitialize();
	ODEIntegrate();
	ODESolutionWrite();

};


void CPressure::readPressure(void)
{
	// open the input file
	ifstream fpressure("./inputs/Pressure.prm");
	
	if (!fpressure.is_open()) 	{
		cout << "CPressure::readPressure -> Unable to open Pressure input file!" << endl;
		exit(1);
	}

	string tmp;
	string line;

	while(fpressure) {
		getline(fpressure,line);
		istringstream iss (line,istringstream::in);
		while(iss) {
			iss >> tmp;
			size_t comment  = tmp.find("//"); // check if there is a comment
			if (comment!=string::npos) {
				while(iss) // go to the end of the line
					iss >> tmp;
			}
			if (tmp == "Vdead")
				iss >> myGapInput.Pressure.Vdead;
			if (tmp == "alphaD_LP")
				iss >> myGapInput.Pressure.alphaD_LP;
			if (tmp == "alphaD_HP")
				iss >> myGapInput.Pressure.alphaD_HP;
			if (tmp == "V_LP")
				iss >> myGapInput.Pressure.V_LP;	
			if (tmp == "V_HP")
				iss >> myGapInput.Pressure.V_HP;				
			if (tmp == "AD_LP")
				iss >> myGapInput.Pressure.AD_LP;
			if (tmp == "AD_HP")
				iss >> myGapInput.Pressure.AD_HP;
			if (tmp == "P1")
				iss >> myGapInput.Pressure.P1;
			if (tmp == "P2")
				iss >> myGapInput.Pressure.P2;	
			if (tmp == "leakageoption")
				iss >> myGapInput.Pressure.leakageoption;
			if (tmp == "QSK")
				iss >> myGapInput.Pressure.QSK;
		}
	}

	fpressure.close();


};
void CPressure::ReadPumpArea(void)
{
	char str[256];
	char cbuffer;
	int i=0;
	double xarea;
	double yarea;
	double zarea;
	Area_file_pump = "./inputs/Pump_Area.ar";

	fstream file_PumpArea;
	file_PumpArea.open(Area_file_pump.c_str(),ios::in);   
	if (file_PumpArea.is_open()) 
	{
		if (file_PumpArea.good()) 
		{
			cout << "Reading PumpArea file...";
		}
		else 
		{
			cout << "Unable to read PumpArea file!" << "\n";
		}
		while(!file_PumpArea.eof() && i < 721) 
		{
			if( file_PumpArea.get(cbuffer) && (cbuffer =='%'))
			{
				file_PumpArea.getline(str,sizeof(str));continue;
			}
			else
			{
				file_PumpArea.putback(cbuffer);
			}
			file_PumpArea >> xarea >> yarea >> zarea;
			PumpVPArea.arphi.push_back(xarea*(PI/180));
			PumpVPArea.arhp.push_back(yarea);
			PumpVPArea.arlp.push_back(zarea);
			i++;
		}
		file_PumpArea.close();
		cout << " Done!";
	}
	else
	{
		cout << "Unable to open PumpArea file!'" << "\n";
	}
	cout << " " << "\n";
	cout << " " << "\n";

}
void CPressure::InitializeValues(void)
{
	//  Note that y is a vector containing all the states to be integrated.
	//  # of states = # pistons pump + 2
	// y[1...npistons] = piston displacement pressure
	// y[npistons+1] = high pressure port pressure
	// y[npistons+2] = low pressure port pressure
	
	//-------------------------------------------------------------------------------------------------
	//STEP 1: Initialize Pump And Motor Displacement Chamber Pressures
	//-------------------------------------------------------------------------------------------------
	//Pistons counted in clockwise direction
	//From second piston to half piston -1 high pressure assigned
	for(int j=1;j<(int)(Pump.npistons/2)+1;j++)
	{
		y[j] = Pump.pHP;
	}
	//From half piston to final piston low pressure is assigned
	for(int j =(int)(Pump.npistons/2)+1; j < Pump.npistons;j++)
	{
		y[j] = Pump.pLP;
	}
	//First piston is assumed low pressure
	y[0] = Pump.pLP;

	//-------------------------------------------------------------------------------------
	//STEP2: Initialize Pump And Motor HP and LP Port Pressures
	//-------------------------------------------------------------------------------------
	//Last two states are assigned high pressure and low pressure ports respectively
	y[Pump.npistons] = Pump.pHP;
	y[Pump.npistons+1] = Pump.pLP;

};
void CPressure::CalcVPArea(void)
{
	
	int i=0;
	double slopehp=0,slopelp=0;
	for (int zi=0; zi<Pump.npistons;zi++)  
	{
		while((PumpVPArea.arphi[i] > Pump.phi_rad[zi]) || (Pump.phi_rad[zi] > PumpVPArea.arphi[i+1]))
		{
			i++;
		}
		slopehp = (PumpVPArea.arhp[i+1]-PumpVPArea.arhp[i])/(PumpVPArea.arphi[i+1]-PumpVPArea.arphi[i]); 
		slopelp = (PumpVPArea.arlp[i+1]-PumpVPArea.arlp[i])/(PumpVPArea.arphi[i+1]-PumpVPArea.arphi[i]); 
		Pump.arhp[zi] = ((Pump.phi_rad[zi]-PumpVPArea.arphi[i])*slopehp+PumpVPArea.arhp[i])*0.000001;		
		Pump.arlp[zi] = ((Pump.phi_rad[zi]-PumpVPArea.arphi[i])*slopelp+PumpVPArea.arlp[i])*0.000001;
		i=0;
	}
}


void CPressure::CalcPistonPositions(double time)
{
	for (int zi=0; zi<Pump.npistons;zi++)  
	{
		Pump.phi_rad[zi] = fmod(Pump.omega*(time)+(zi)*(2*PI/Pump.npistons)+2*PI,2*PI);
	}
	
}
void CPressure::CalcIntermediateVariables(void)
{

//-----------------------------------------------------------------------------------------
// STEP 2: Pump Variables
//-----------------------------------------------------------------------------------------	
Pump.beta_rad = Pump.beta_deg*(PI/180);							//Displacement angle [rad]
Pump.betamax_rad = Pump.betamax_deg*(PI/180);			    //Displacement angle max [rad]
Pump.gamma_rad = Pump.gamma_deg*(PI/180);					//Cross angle [rad]
if (Pump.mode == 2)																	//Motoring mode: negative displacement angle
{
	Pump.beta_rad *= -1.0;
}
Pump.dK=Pump.dK/1000;														//Piston diameter conversion [mm] to [m]
Pump.dB=Pump.dB/1000;														//Cylinder block pitch diameter conversion [mm] to [m]
Pump.pHP=Pump.pHP*100000;													//Pump high pressure conversion [bar] to [Pa]							
Pump.pLP=Pump.pLP*100000;													//Pump low pressure conversion [bar] to [Pa]	
Pump.P1=Pump.P1*100000;													//Pump upstream pressure conversion [bar] to [Pa]	
Pump.P2=Pump.P2*100000;													//Pump downstream pressure conversion [bar] to [Pa]	
Pump.pCase=Pump.pCase*100000;										//Pump case pressure conversion [bar] to [Pa]	
Pump.Ap = PI*(pow(Pump.dK,2)/4);										// Piston area [m^2]
Pump.omega = Pump.speed*(PI/30);									    // Angular velocity [rad/s]

//-----------------------------------------------------------------------------------------
// STEP 2: Simulation Variables
//-----------------------------------------------------------------------------------------	
Sim.DegStep = 0.2;																	// Constant angular step [deg]
Sim.tDelta = Sim.DegStep/(6.0*Pump.speed);					//  Simulation time step [s]
Sim.tEnd = Pump.nrevolutions*(60/Pump.speed);				//  Simulation end time [s]

}

void CPressure::CalcTheoreticalFlow(void)
{
	SV.QHPtheo = 0.0;
	SV.QLPtheo = 0.0;
	for (int zi=0; zi<Pump.npistons;zi++)  
	{
		SV.QHPitheo[zi] = Pump.omega*Pump.Ap*Pump.dB/2*tan(Pump.beta_rad)*sin(Pump.phi_rad[zi]);
		if 	(SV.QHPitheo[zi] < 0)
		{
			SV.QHPitheo[zi] = 0;
		};
		SV.QHPtheo+= SV.QHPitheo[zi];
	};
	SV.QLPtheo = -SV.QHPtheo;

};


int CPressure::CalcSgn(double a)  
{
	if (a>0) 
		return 1; 
	else if (a<0) 
		return -1;  
	else  return 0; 
}
double CPressure::CalcDensity(const double& T_C,const double& p_Pa)
{
	double RhoT=0.0, Rho=0.0;
	double T_K = T_C + 273.15;
	
	//Oiltype option 0: User defined constant density 
	if (Oil.oiltype == 0)
	{
		Rho= Oil.oildensity; //[kg/m3]
	}
	//Oiltype option 1: User defined linear formula density
	if (Oil.oiltype == 1) 
	{
		double Rhozero = Oil.oildensity; //[kg/m3]
		double betaT = Oil.oilbetaT; //[1/K]
		double betaP = Oil.oilbetaP; //[1/Pa]
		Rho = Rhozero*(1 + betaP*(p_Pa) - betaT*(T_C - 20)); // Rho(T,p)
	}
	//Oiltype option 2: Oil is HLP 32
	if (Oil.oiltype == 2) 
	{
		double rs = 1047.03; // [kg/m3]
		double als = 5.761668e-4; // [1/K]
		double a1 = 0.0732965;
		double a2 = 1965.02; // [bar]
		double a3 = -2.96813; // [bar/K]
		RhoT = rs*(1-als*T_K); 
		Rho= RhoT/(1-a1*log((a2+a3*T_K+p_Pa*1.0e-5)/(a2+a3*T_K))); // Rho(T,p)
	}
	if (Oil.oiltype == 3)
	{
		//Skydrol - model developed using AMEsim & SAE specs by Marco Zecchi
		double pref = 1.00000000000000e005;
		double Tref = 1.00000000000000e001;
		double rhoref = 1.06690035585407e003;
		double muref = 4.78435154853258e-002;
		double cpref = 1.67717744042304e003;
		double lambdaref = 1.36331614988375e-001;

		double T1 = 7.80474256715724e-004;
		double T2 = 8.33673780452984e-007;
		double P1 = -8.76030682143087e-010;
		double P2 = 3.76207238974272e-018;
		double PT = -5.44707978002713e-012;
		double PT2 = 1.17343267760062e-014;
		double P2T = 2.78987659880380e-020;
		double P2T2 = -7.01924728092232e-023;
		
		double dp = p_Pa - pref;
		double dT = T_C - Tref;

		double vs = ( 
							(1.0/rhoref)*
							( 
								1 + P1*dp + P2*pow(dp, 2.0) + T1*dT + T2*pow(dT, 2.0) +  
								PT*dp*dT + PT2*dp*pow(dT,2) + P2T*dT*pow(dp,2) +
								P2T2*pow(dp, 2)*pow(dT, 2)
							)
						);
			
		Rho = 1.0/vs;

	}
	return Rho;
}

double CPressure::CalcK(const double& T_C,const double& p_Pa)
{
	double RhoT =0.0, RhoTp=0.0, K=0.0;
	double T_K = T_C + 273.15;
	//Oiltype option 0: User defined constant bulk modulus
	if (Oil.oiltype == 0)
	{
		K = Oil.oilK; // [Pa]
	}
	//Oiltype option 1: User defined linear formula bulk modulus
	if (Oil.oiltype == 1)
	{
		double Kzero = Oil.oilK; //[Pa]
		double betaKT = Oil.oilbetaKT; //[1/K]
		double betaKP = Oil.oilbetaKP; //[1/Pa]
		K= Kzero*(1 + betaKP*(p_Pa-1e5) - betaKT*(T_C - 20)); // Rho(T,p)
	}
	//Oiltype option 2: Oil is HLP 32
	if (Oil.oiltype == 2) //Oil is HLP 32
	{
		double rs = 1047.03; // [kg/m3]
		double als = 5.761668e-4; // [1/K]
		double a1 = 0.0732965;
		double a2 = 1965.02; // [bar]
		double a3 = -2.96813; // [bar/K]
		RhoT = rs*(1-als*T_K); 
		RhoTp= RhoT/(1-a1*log((a2+a3*T_K+p_Pa*1.0e-5)/(a2+a3*T_K)));
		K = 1.0e5*(RhoT*(a2+a3*T_K+p_Pa*1.0e-5)/(a1*RhoTp)); //[Pa]
	}
	if (Oil.oiltype == 3) //Oil is Skydrol
	{
		//Skydrol - model developed using AMEsim & SAE specs by Marco Zecchi
		//Using definition of Bulk Modulus = rho / (drho/dp)
		//drho/dp is the analytical expression from the rho Skydrol model
		
		double pref = 1.00000000000000e005;
		double Tref = 1.00000000000000e001;
		double rhoref = 1.06690035585407e003;
		double muref = 4.78435154853258e-002;
		double cpref = 1.67717744042304e003;
		double lambdaref = 1.36331614988375e-001;

		double T1 = 7.80474256715724e-004;
		double T2 = 8.33673780452984e-007;
		double P1 = -8.76030682143087e-010;
		double P2 = 3.76207238974272e-018;
		double PT = -5.44707978002713e-012;
		double PT2 = 1.17343267760062e-014;
		double P2T = 2.78987659880380e-020;
		double P2T2 = -7.01924728092232e-023;
		
		double dp = p_Pa - pref;
		double dT = T_C - Tref;

		double drhodp =
		-(rhoref*(P1 + 2*P2*dp + PT*dT + PT2*dT*dT + 2*P2T*dT*dp + 2*P2T2*dT*dT*dp))
			/
		pow(P2T2*dT*dT*dp*dp + PT2*dT*dT*dp + T2*dT*dT + P2T*dT*dp*dp + PT*dT*dp + T1*dT + P2*dp*dp + P1*dp + 1,2.0);

		K = CalcDensity(T_C, p_Pa)/drhodp;
	}

	return K;
}

void CPressure::odef(double time,vector<double> &y,vector<double> &yp)
{
	//States vector with pressure elements
	SV.pD = y;

	//----------------------------------------------------------------------------------------
	//STEP 1: Evaluate Pump Pistons Angular Positions (phi) and Corresponding VP Areas
	//----------------------------------------------------------------------------------------
    CalcPistonPositions(time);
	CalcVPArea();

	//--------------------------------------------------------------------
	//STEP 2: Set Pressure Limits to Avoid Very Low Pressure Values
	//--------------------------------------------------------------------
	for(int j=0;j<Pump.npistons+2;j++)
	{
		if( y[j] < -0.8e5)
		{
			y[j] = -0.8e5;
		}
	};

	//-----------------------------------------------------------
	//STEP 3: Evaluate the variables inside each piston
	//-----------------------------------------------------------
	double Vp=0,Vo=0,tempHP=0,tempLP=0,temp=0,C=0;
	double P2 = Pump.P2;
	double P1 = Pump.P1;
	double pHP = y[Pump.npistons];
	double pLP = y[Pump.npistons+1];
	SV.QrHP = 0;
	SV.QrLP = 0;

	for(int z=0;z<Pump.npistons;z++)
	{
		double pD = SV.pD.at(z);
		//Current displacement volume
		Vo =  Pump.Vdead + Pump.dB/2.0*Pump.Ap*(tan(Pump.beta_rad) + tan(Pump.betamax_rad));
		SV.Vo = Vo;
		//Change in Volume over time 
		Vp =  -Pump.omega*Pump.dB/2.0*(tan(Pump.beta_rad)*sin(Pump.phi_rad[z])+tan(Pump.gamma_rad)*cos(Pump.phi_rad[z])/cos(Pump.beta_rad));
		SV.dVdT[z] =Vp*Pump.Ap;
		//Flow in/out of one piston
		tempHP = Pump.alphaD_HP*pow((2/CalcDensity(Pump.THP,pD)),0.5);
		SV.QrHPi[z] = tempHP*Pump.arhp[z]*pow((fabs(-pHP+pD)),0.5)*CalcSgn(pHP-pD);
		tempLP = Pump.alphaD_LP*pow((2/CalcDensity(Pump.TLP,pD)),0.5);
		SV.QrLPi[z] = tempLP*Pump.arlp[z]*pow((fabs(-pLP+pD)),0.5)*CalcSgn(pLP-pD); 
		//Leakage through piston/cylinder and block/valve plate gap
		//Add options
		SV.Qrleaki[z]= Pump.QSK;
		SV.QsKi[z]=0.0;
		SV.QsGi[z]=0.0;
		SV.QsBi[z]=0.0;
		//Hydraulic Capacitance
		double psi = -1.0 * atan(tan(Pump.gamma_rad)/sin(Pump.beta_rad));
		C = Vo - Pump.Ap*
			(
			Pump.dB/2*tan(Pump.beta_rad)*(1-cos(Pump.phi_rad[z]))
			+ Pump.dB/2*tan(Pump.gamma_rad)*sin(Pump.phi_rad[z])/cos(Pump.beta_rad)
			- Pump.dB/2*tan(Pump.beta_rad)*(1-cos(psi))
			- Pump.dB/2 * tan(Pump.gamma_rad) * sin(psi) / cos(Pump.beta_rad)
			)
			;
		//Pressure Built up Equation: change in pressure inside each chamber
		SV.dpD[z] = ((CalcK(Pump.THP,pD))/C)*(-SV.dVdT[z]+SV.QrHPi[z]+SV.QrLPi[z]-SV.Qrleaki[z]);
		//Summation of in/out flow for all pistons
		SV.QrHP+= SV.QrHPi[z];
		SV.QrLP+= SV.QrLPi[z];
	};

	//-------------------------------------------------------------
	//STEP 4: Evaluate dpHP and dpLP 
	//-------------------------------------------------------------
	// Downstream flow rate Q2 calculation from pump to load
	temp = 2*fabs(pHP - P2)/CalcDensity(Pump.THP,pHP);
	temp = sqrt(temp);
	SV.Q2 = Pump.alphaD_HP*Pump.AD_HP*temp*CalcSgn(P2 - pHP);

	// Upstream flow rate Q1 calculation from inlet to pump
	temp = 2*fabs(pLP - P1)/CalcDensity(Pump.TLP,pLP);
	temp = sqrt(temp);
	SV.Q1 = Pump.alphaD_LP*Pump.AD_LP*temp*CalcSgn(P1-pLP);

	//Pressure Built Up for HP and LP ports based on Q2 and Q1 respectively
	SV.dpD[Pump.npistons] = ((CalcK(Pump.THP,pHP))/Pump.V_HP)*(SV.Q2 - SV.QrHP);
	SV.dpD[Pump.npistons+1] = ((CalcK(Pump.TLP,pLP))/Pump.V_LP)*(SV.Q1 - SV.QrLP);

	//---------------------------------------------------------------
	//STEP 5: Store the calculated dpD into the ouput container		
	//---------------------------------------------------------------
	yp = SV.dpD;

	return;

};

void CPressure::ODEInitialize()
{
	//Checking the inputs and setting appropriate values
	// n, the dimension of the system
	if (n == UINT_MAX) {
		cout << "System too big, max. n = " << UINT_MAX - 1 << "\n";
		throw -1;
	}

	// rtoler, the relative tolerance of the integration
	if (!rtoler_s) {
		itoler = 0;
		rtoler_s = 1.0e-7;
		//rtolerNULL = true;
	}

	// atoler, the absolute tolerance of the integration
	if (!atoler_s) {
		itoler = 0;
		atoler_s = 1.0e-7;
		//atolerNULL = true;
	}

	// -------- maximal step size
	if (hmax == 0.0) hmax = xend - x;

	// -------- nmax--maximal number of steps
	if (nmax == 0) nmax = 10000000;
	if (nmax <= 0) {
		cout << " wrong input, nmax = " << nmax << "\n";
		throw -1;
	}

	// -------- uround--smallest number satisfying 1.0 + uround > 1.0
	if (uround == 0.0) uround = 1.0e-16;
	if ((uround <= 1.0e-19) || (uround >= 1.0)) {
		cout << " coefficients have 20 digits, uround = " << uround << "\n";
		throw -1;
	}

	// --------- safe--safety factor in step size prediction
	if (safe == 0.0) safe = 0.9;
	if ((safe <= 0.001) || (safe >= 1.0)) {
		cout << " curious input for safety factor, safe = " << safe << "\n";
		throw -1;
	}
	// iout, switch for calling SolutionOutput
	if ((iout < 0) || (iout > 2)) {
		cout << "Wrong input, iout = " << iout << "\n";
		throw -1;
	}

	// facl, facr, parameters for step size selection
	if (facl == 0.0) facl = 0.2;
	if (facr == 0.0) facr = 10.0;

	// beta for step control stabilization
	if (beta == 0.0)
		beta = 0.04;
	else if (beta < 0.0)
		beta = 0.0;
	else if (beta > 0.2) {
		cout << "Curious input for beta : beta = " << beta << "\n";
		throw -1;
	}

	// nstiff, parameter for stiffness detection
	if (nstiff == 0) nstiff = 1000;
	else if (nstiff < 0) nstiff = nmax + 10;

	indir = NULL;
	//	rcont1 = rcont2 = rcont3 = rcont4 = rcont5 = NULL;

	// nrdens, number of dense output components
	if (nrdens > n) {
		cout << "Curious input, nrdens = " << nrdens << "\n";
		throw -1;
	}
	else if (nrdens != 0) {
		rcont1.resize(nrdens); 
		rcont2.resize(nrdens);
		rcont3.resize(nrdens);
		rcont4.resize(nrdens);
		rcont5.resize(nrdens);
		if (nrdens < n) indir = new int[n];

		// control of length of icont
		if (nrdens == n) {
			if (icont) cout << "Warning : when nrdens = n there is no need " <<
				"allocating memory for icont\r\n";
		}
		else {
			if (iout < 2) cout << "Warning : put iout = 2 for dense output\r\n";
			for (int i = 0; i < n; i++) indir[i] = UINT_MAX;
			for (int i = 0; i < nrdens; i++) indir[icont[i]] = i;
		}
	}

}

double CPressure::ODEhinit()
{
	int iord = 5;
	double posneg = ODEsign(1.0, xend-x);
	double sk, sqr;
	double dnf = 0.0;
	double dny = 0.0;

	vector<double> atoli = atoler;
	vector<double> rtoli = rtoler;
	double atoli_s = atoler_s;
	double rtoli_s = rtoler_s;


	if (itoler == 0)
	{
		for (int i = 0; i < n; i++) 
		{
			sk = atoli_s + rtoli_s * fabs(y[i]);
			sqr = k1[i]/sk;
			dnf += sqr*sqr;
			sqr = y[i]/sk;
			dny += sqr*sqr;
		}
	}
	else
	{

		for (int i = 0; i < n; i++) 
		{
			sk = atoler[i] + rtoler[i] * fabs(y[i]);
			sqr = k1[i]/sk;
			dnf += sqr*sqr;
			sqr = y[i]/sk;
			dny += sqr*sqr;
		}
	}
	if ((dnf <= 1.0e-10) || (dny <= 1.0e-10)) h = 1.0e-6;
	else h = sqrt(dny/dnf)*0.01;

	h = min(h, hmax);
	h = ODEsign(h, posneg);

	// perform an explicit Euler step
	for (int i = 0; i < n; i++) k3[i] = y[i] + h * k1[i];

	odef(x+h, k3, k2);

	// estimate the second derivative of the solution
	double der2 = 0.0;
	if (itoler == 0)
		for (int i = 0; i < n; i++) {
			sk = atoli_s + rtoli_s * fabs(y[i]);
			sqr = (k2[i] - k1[i])/sk;
			der2 += sqr*sqr;
		}
	else
		for (int i = 0; i < n; i++) 
		{
			sk = atoler[i] + rtoler[i] * fabs(y[i]);
			sqr = (k2[i] - k1[i])/sk;
			der2 += sqr*sqr;
		}

		der2 = sqrt(der2)/h;

		// step size is computed such that
		// h**iord * max(norm(k1), norm(der2)) = 0.01
		double der12 = max(fabs(der2), sqrt(dnf));

		double h1;
		if (der12 <= 1.0e-15) h1 = max(1.0e-6, fabs(h)*1.0e-3);
		else h1 = pow(0.01/der12, 1.0/(double)iord);

		h = min(100.0*h, min(h1, hmax));

		return ODEsign(h, posneg);

} // hinit
void CPressure::ODEIntegrate()
{
	//system("cls");
	int idid = ODECoreIntegrator();

	if (idid < 0) 
	{
		cout << " Computation failed " << "\n";
		return;
	}
	return;
} // Integrate


double CPressure::ODEsign(double a, double b)
{
	return (b > 0.0) ? fabs(a) : -fabs(a);
} // sign

int CPressure::ODECoreIntegrator()
{
	double c2, c3, c4, c5, e1, e3, e4, e5, e6, e7, d1, d3, d4, d5, d6, d7;
	double a21, a31, a32, a41, a42, a43, a51, a52, a53, a54;
	double a61, a62, a63, a64, a65, a71, a73, a74, a75, a76;

	// initialisations
	switch (meth)
	{
	case 1:

		c2=0.2, c3=0.3, c4=0.8, c5=8.0/9.0;
		a21=0.2, a31=3.0/40.0, a32=9.0/40.0;
		a41=44.0/45.0, a42=-56.0/15.0; a43=32.0/9.0;
		a51=19372.0/6561.0, a52=-25360.0/2187.0;
		a53=64448.0/6561.0, a54=-212.0/729.0;
		a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0;
		a64=49.0/176.0, a65=-5103.0/18656.0;
		a71=35.0/384.0, a73=500.0/1113.0, a74=125.0/192.0;
		a75=-2187.0/6784.0, a76=11.0/84.0;
		e1=71.0/57600.0, e3=-71.0/16695.0, e4=71.0/1920.0;
		e5=-17253.0/339200.0, e6=22.0/525.0, e7=-1.0/40.0;
		d1=-12715105075.0/11282082432.0, d3=87487479700.0/32700410799.0;
		d4=-10690763975.0/1880347072.0, d5=701980252875.0/199316789632.0;
		d6=-1453857185.0/822651844.0, d7=69997945.0/29380423.0;

		break;
	}

	double posneg = ODEsign(1.0, xend-x);
	double facold = 1.0e-4;
	double expo1 = 0.2 - beta*0.75;
	double facc1 = 1.0/facl;
	double facc2 = 1.0/facr;

	// initial preparations
	double atoli_s = atoler_s;
	double rtoli_s = rtoler_s;

	bool last = false;
	double hlamb = 0.0;
	int iasti = 0;

	odef(x, y, k1);

	hmax = fabs(hmax);
	if (h == 0.0) h = ODEhinit();

	nfcn += 2;
	bool reject = false;
	int nonsti=0;

	if (iout != 0) {
		hold = h;
		int irtrn = ODESolutionOutput();
		if (irtrn < 0) {
			cout << "Exit of dopri5 at x = " << x << "\n";
			return 2;
		}
	}

	// basic integration step
	while (true) {
		if (nstep > nmax) {
			cout << "Exit of dopri5 at x = " << x << ", more than nmax = "
				<< nmax << " are needed." << "\n";
			hold = h;
			return -2;
		}
		if (0.1*fabs(h) <= fabs(x)*uround) {
			cout << "Exit of dopri5 at x = " << x << ", step size too small h = "
				<< h << "\n";
			hold = h;
			return -3;
		}
		if ((x + 1.01*h - xend)*posneg > 0.0) {
			h = xend - x;
			last = true;
		}
		nstep++;

		// the first 6 stages
		for (int i = 0; i < n; i++)
			yy1[i] = y[i] + h*a21*k1[i];

		odef(x+c2*h, yy1, k2);

		for (int i = 0; i < n; i++)
			yy1[i] = y[i] + h*(a31*k1[i] + a32*k2[i]);

		odef(x+c3*h, yy1, k3);

		for (int i = 0; i < n; i++)
			yy1[i] = y[i] + h*(a41*k1[i] + a42*k2[i] + a43*k3[i]);

		odef(x+c4*h, yy1, k4);

		for (int i = 0; i <n; i++)
			yy1[i] = y[i] + h*(a51*k1[i] + a52*k2[i] + a53*k3[i] + a54*k4[i]);

		odef(x+c5*h, yy1, k5);

		for (int i = 0; i < n; i++)
			ysti[i] = y[i] + h*(a61*k1[i] + a62*k2[i] + a63*k3[i] +
			a64*k4[i] + a65*k5[i]);

		odef(x+h, ysti, k6);

		for (int i = 0; i < n; i++)
			yy1[i] = y[i] + h*(a71*k1[i] + a73*k3[i] + a74*k4[i] +
			a75*k5[i] + a76*k6[i]);

		odef(x+h, yy1, k2);

		if (iout == 2) {
			if (nrdens == n) {
				for (int i = 0; i < n; i++) {
					rcont5[i] = h*(d1*k1[i] + d3*k3[i] + d4*k4[i] +
						d5*k5[i] + d6*k6[i] + d7*k2[i]);
				}
			}
			else {
				for (int j = 0; j < nrdens; j++) {
					unsigned i = icont[j];
					rcont5[j] = h*(d1*k1[i] + d3*k3[i] + d4*k4[i] +
						d5*k5[i] + d6*k6[i] + d7*k2[i]);
				}
			}
		}

		for (int i = 0; i < n; i++)
			k4[i] = h*(e1*k1[i] + e3*k3[i] + e4*k4[i] + e5*k5[i] +
			e6*k6[i] + e7*k2[i]);
		nfcn += 6;

		// error estimation
		double err = 0.0, sk, sqr;
		if (!itoler)
			for (int i = 0; i < n; i++) {
				sk = atoli_s + rtoli_s*max(fabs(y[i]), fabs(yy1[i]));
				sqr = k4[i]/sk;
				err += sqr*sqr;
			}
		else
			for (int i = 0; i < n; i++) 
			{
				sk = atoler[i] + rtoler[i]*max(fabs(y[i]), fabs(yy1[i]));
				sqr = k4[i]/sk;
				err += sqr*sqr;
			}

			err = sqrt(err/(double)n);

			// computation of hnew
			double fac11 = pow(err, expo1);
			// Lund-stabilization
			double fac = fac11/pow(facold,beta);
			// we require facl <= hnew/h <= facr
			fac = max(facc2, min(facc1, fac/safe));
			double hnew = h/fac;

			if (err <= 1.0) 
			{
				// step accepted
				facold = max(err, 1.0e-4);
				naccpt++;

				// stiffness detection
				if (!(naccpt % nstiff) || (iasti > 0)) 
				{
					double stnum = 0.0, stden = 0.0;
					for (int i = 0; i < n; i++) 
					{
						sqr = k2[i] - k6[i];
						stnum += sqr*sqr;
						sqr = yy1[i] - ysti[i];
						stden += sqr*sqr;
					}
					if (stden > 0.0) hlamb = h*sqrt(stnum/stden);
					if (hlamb > 3.25) 
					{
						nonsti = 0;
						iasti++;
						if (iasti == 15) {
							cout << "The problem seems to become stiff at x = " << x << "\n";
							hold = h;
							return -4;
						}
					}
					else {
						nonsti++;
						if (nonsti == 6) iasti = 0;
					}
				}

				if (iout == 2) {
					double yd0, ydiff, bspl;
					if (nrdens == n)
						for (int i = 0; i < n; i++) {
							yd0 = y[i];
							ydiff = yy1[i] - yd0;
							bspl = h*k1[i] - ydiff;
							rcont1[i] = y[i];
							rcont2[i] = ydiff;
							rcont3[i] = bspl;
							rcont4[i] = -h*k2[i] + ydiff - bspl;
						}
					else
						for (int j = 0; j < nrdens; j++) {
							unsigned i = icont[j];
							yd0 = y[i];
							ydiff = yy1[i] - yd0;
							bspl = h * k1[i] - ydiff;
							rcont1[j] = y[i];
							rcont2[j] = ydiff;
							rcont3[j] = bspl;
							rcont4[j] = -h * k2[i] + ydiff - bspl;
						}
				}

				//memcpy(k1, k2, n*sizeof(double));
				//memcpy(y, yy1, n*sizeof(double));
				copy(k2.begin(),k2.end(),k1.begin());
				copy(yy1.begin(),yy1.end(),y.begin());
				xold = x;
				x += h;

				if (iout) {
					hold = h;
					int irtrn = ODESolutionOutput();
					if (irtrn < 0) {
						cout << "Exit of dopri5 at x = " << x << "\n";
						return 2;
					}
				}

				// normal exit
				if (last) {
					hold = hnew;
					return 1;
				}

				if (fabs(hnew) > hmax) hnew = posneg*hmax;
				if (reject) hnew = posneg*min(fabs(hnew), fabs(h));
				reject = false;
			}
			else {
				// step rejected
				hnew = h/min(facc1, fac11/safe);
				reject = true;
				if (naccpt >= 1) nrejct++;
				last = false;
			}
			h = hnew;
	}

} // CoreIntegrator

double CPressure::ODEContinuousOutput(unsigned i)
{
	unsigned ii = UINT_MAX;

	if (!indir) ii = i;
	else ii = indir[i];

	if (ii == UINT_MAX) {
		cout << "No dense output available for %uth component" << i << "\n";
		return 0.0;
	}

	double theta = (xd - xold)/hold;
	double theta1 = 1.0 - theta;

	return rcont1[ii] + theta*(rcont2[ii] +
		theta1*(rcont3[ii] + theta*(rcont4[ii] + theta1*rcont5[ii])));

} // contd5


int CPressure::ODESolutionOutput()
{
	output << setiosflags(ios::showpoint);

	if (naccpt == 0) xd = xold;

	while (xd < x) {
		if ((xold <= xd) && (x >= xd)) 
		{
			//PI Controller on the inlet and outlet areas
			double pHPsp,pLPsp,KpHP,KpLP,KiHP,KiLP,pHP,pLP,steprev;
			if(Sim.stepcounter==(int)360/Sim.DegStep)
			{
				steprev = Sim.stepcounter;
				Pump.AD_HP = Pump.AD_HPsum/steprev;
				Pump.AD_LP = Pump.AD_LPsum/steprev;
				Sim.stepcounter=0;
				Sim.revcounter++;
				
			}
			//First revolutions chasing right areas until error is less then 3 bars
			if(Sim.revcounter==0)
			{
				//Set point pressures
				pHPsp = Pump.pHP*1e-5;
				pLPsp = Pump.pLP*1e-5;
				pHP = y[Pump.npistons]*1e-5;
				pLP = y[Pump.npistons+1]*1e-5;
				Pump.errorHP = pHPsp-pHP;
				Pump.errorLP = pLPsp-pLP;
				Pump.errorHPsum += Pump.errorHP;
				Pump.errorLPsum += Pump.errorLP;
				//Initialize gains
				if(Sim.stepcounter==0)
				{
					KpHP = 1e-9;
					KiHP = 1e-10;
					KiLP = 1e-9;
					KpLP = 1e-10;
				}
				//Gain Scheduling Large speeds and displacements
				if(abs(Pump.errorHP) >= 50)
				{
					KpHP = 1e-6; KiHP = 1e-13;
				}
				if(abs(Pump.errorHP) >= 5 && abs(Pump.errorHP) < 50)
				{
					KpHP = 1e-7; KiHP = 1e-14;
				}
				if(abs(Pump.errorHP) >= 1 && abs(Pump.errorHP) < 5)
				{
					KpHP = 1e-8; KiHP = 1e-15;
				}
				if(abs(Pump.errorHP) >= 0.0 && abs(Pump.errorHP) < 1)
				{
					KpHP = 1e-9; KiHP = 1e-16;
				}

				if(abs(Pump.errorLP) >= 10)
				{
					KpLP = 1e-6; KiLP = 1e-13;
				}
				if(abs(Pump.errorLP) >= 5 && abs(Pump.errorLP) < 10)
				{
					KpLP = 1e-7; KiLP = 1e-14;
				}
				if(abs(Pump.errorLP) >= 1 && abs(Pump.errorLP) < 5)
				{
					KpLP = 1e-8; KiLP = 1e-15;
				}
				if(abs(Pump.errorLP) >= 0.0 && abs(Pump.errorLP) < 1)
				{
					KpLP = 1e-9; KiLP = 1e-16;
				}


				//Gain Scheduling Small displacements
				/*if(abs(Pump.errorHP) >= 50)
				{
					KpHP = 1e-9; KiHP = 1e-13;
				}
				if(abs(Pump.errorHP) >= 5 && abs(Pump.errorHP) < 50)
				{
					KpHP = 1e-10; KiHP = 1e-14;
				}
				if(abs(Pump.errorHP) >= 1 && abs(Pump.errorHP) < 5)
				{
					KpHP = 1e-11; KiHP = 1e-15;
				}
				if(abs(Pump.errorHP) >= 0.1 && abs(Pump.errorHP) < 1)
				{
					KpHP = 1e-12; KiHP = 1e-16;
				}

				if(abs(Pump.errorLP) >= 10)
				{
					KpLP = 1e-9; KiLP = 1e-13;
				}
				if(abs(Pump.errorLP) >= 5 && abs(Pump.errorLP) < 10)
				{
					KpLP = 1e-10; KiLP = 1e-14;
				}
				if(abs(Pump.errorLP) >= 1 && abs(Pump.errorLP) < 5)
				{
					KpLP = 1e-11; KiLP = 1e-15;
				}
				if(abs(Pump.errorLP) >= 0.1 && abs(Pump.errorLP) < 1)
				{
					KpLP = 1e-12; KiLP = 1e-16;
				}*/
				
				//Control
				if(Pump.mode == 2)
				{
					Pump.AD_HP += (KpHP*Pump.errorHP + KiHP*Pump.errorHPsum);
					Pump.AD_LP -= (KpLP*Pump.errorLP + KiLP*Pump.errorLPsum);
				} else {
					Pump.AD_HP -= (KpHP*Pump.errorHP + KiHP*Pump.errorHPsum);
					Pump.AD_LP += (KpLP*Pump.errorLP + KiLP*Pump.errorLPsum);
				}
				Pump.AD_HPsum += Pump.AD_HP;
				Pump.AD_LPsum += Pump.AD_LP;
			}
			//Further revolutions for precise adjustment to match deltap
			else if(Sim.revcounter>0 && Sim.revcounter < Pump.nrevolutions-1)
			{
				//Set point pressures
				pHPsp = Pump.pHP*1e-5;
				pLPsp = Pump.pLP*1e-5;
				pHP = y[Pump.npistons]*1e-5;
				pLP = y[Pump.npistons+1]*1e-5;
				Pump.pHPsum+=pHP;
				Pump.pLPsum+=pLP;
				if(Sim.stepcounterint == (int) 360/Pump.npistons/Sim.DegStep)
				{
					Pump.pHPsum=Pump.pHPsum/Sim.stepcounterint;
					Pump.pLPsum=Pump.pLPsum/Sim.stepcounterint;
					Pump.errorHP2 = pHPsp - Pump.pHPsum;
					Pump.errorLP2 = pLPsp - Pump.pLPsum;

					Pump.errorHPsum2 += Pump.errorHP2;
					Pump.errorLPsum2 += Pump.errorLP2;
					if(abs(Pump.errorHP2)>1.0)
					{
						KpHP = 2e-9;
						KiHP = 2e-13;
						if(Pump.mode == 2)
						{
							Pump.AD_HP += (KpHP*Pump.errorHP2 + KiHP*Pump.errorHPsum2);
						} else {
							Pump.AD_HP -= (KpHP*Pump.errorHP2 + KiHP*Pump.errorHPsum2);
						}
					}
					if(abs(Pump.errorLP2)>1.0)
					{
						KpLP = 2e-9;
						KiLP = 2e-13;
						if(Pump.mode == 2)
						{
							Pump.AD_LP -= (KpLP*Pump.errorLP2 + KiLP*Pump.errorHPsum2);
						} else {
							Pump.AD_LP += (KpLP*Pump.errorLP2 + KiLP*Pump.errorHPsum2);
						}
					}
					Sim.stepcounterint = 0;
					Pump.pHPsum = 0;
					Pump.pLPsum = 0;
				}
				Sim.stepcounterint++;
			}
			//Last revolution with final calculated areas kept fixed
			else
			{
				Pump.AD_HP = Pump.AD_HP;
				Pump.AD_LP = Pump.AD_LP;
			};
			Sim.stepcounter++;

			//Outputs_Piston
			CalcTheoreticalFlow();
			Result.t.push_back(x);
			Result.phi.push_back((Pump.phi_rad[0])*180/PI);
			Result.pD.push_back(y[0]);
			Result.dPD.push_back(SV.dpD[0]);
			Result.QrHPi.push_back(SV.QrHPi[0]);
			Result.QrLPi.push_back(SV.QrLPi[0]);
			Result.Qri.push_back(SV.QrHPi[0]+SV.QrLPi[0]);
			Result.QsK.push_back(SV.QsKi[0]);
			Result.QsG.push_back(SV.QsGi[0]);
			Result.QsB.push_back(SV.QsBi[0]);
			Result.dVdT.push_back(SV.dVdT[0]);
			Result.Vo.push_back(SV.Vo);
			Result.arhp.push_back(Pump.arhp[0]);
			Result.arlp.push_back(Pump.arlp[0]);
			Result.pHP.push_back(y[Pump.npistons]);
			Result.pLP.push_back(y[Pump.npistons+1]);
			Result.QrHP.push_back(-SV.QrHP);
			Result.QrLP.push_back(-SV.QrLP);
			Result.QHPtheo.push_back(SV.QHPtheo);
			Result.QLPtheo.push_back(SV.QLPtheo);


			//Console and file writing
			cout.setf(0,ios::floatfield);
			cout.precision(5);
			if(Sim.dispcounter > 21)
			{
				Sim.dispcounter = 0;
				system("cls");
			}
			if(Sim.dispcounter == 0)
			{
				cout << "Deg" << "\t" << "HP" << "\t" << "LP" << "\t" << "dP" << "\t" << "AD_HP" << "\t" << "AD_LP";
				cout << "\t\t\t   Rev: " << Sim.revcounter+1 << " of " << Pump.nrevolutions << "\n";
				cout << "--------------------------------------------------------------------------------";
			}
			Sim.dispcounter++;
			cout << (Pump.phi_rad[0])*180/PI << "\t" << y[Pump.npistons]/100000 << "\t" << y[Pump.npistons+1]/100000 << "\t" << (y[Pump.npistons]-y[Pump.npistons+1])/100000 << "\t" << Pump.AD_HP*1e6 << "\t" << Pump.AD_LP*1e6 << "\t" <<  "\n";
			

			xd += dx;
		}
	}


	return 0;
};
void CPressure::ODESolutionWrite()
{
		//cout the final areas
		ostringstream oss (ostringstream::out);
		oss.str("");
		oss << "--------------------------------------------------" << "\n";
		oss << "Final Areas:" << "\n";
		oss << "AD_HP: " << Pump.AD_HP << "\tAD_LP: " << Pump.AD_LP << "\n";
		oss << "--------------------------------------------------";

		cout << "Writing pressure files..." << "\n";
		int tend = (int) Result.t.size();

		//write dp file
		fsolutionwrite.open("./outputs_piston/Pressure_Module_Output.dp");
		for(int i=0; i < tend;i++)
		{
			fsolutionwrite << setw(20) << Result.t[i];
			fsolutionwrite << setw(20) << Result.phi[i];
			fsolutionwrite << setw(20) << Result.pD[i]/1e5;
			fsolutionwrite << setw(20) << Result.dPD[i];
			fsolutionwrite << setw(20) << Result.QrHPi[i]*6e4;
			fsolutionwrite << setw(20) << Result.QrLPi[i]*6e4;
			fsolutionwrite << setw(20) << Result.Qri[i]*6e4;
			fsolutionwrite << setw(20) << Result.QsK[i]*6e4;
			fsolutionwrite << setw(20) << Result.QsG[i]*6e4;
			fsolutionwrite << setw(20) << Result.QsB[i]*6e4;
			fsolutionwrite << setw(20) << Result.dVdT[i];
			fsolutionwrite << setw(20) << Result.Vo[i];
			fsolutionwrite << setw(20) << Result.arhp[i];
			fsolutionwrite << setw(20) << Result.arlp[i];
			fsolutionwrite << setw(20) << Result.pHP[i]/1e5;
			fsolutionwrite << setw(20) << Result.pLP[i]/1e5;
			fsolutionwrite << setw(20) << Result.QrHP[i]*6e4;
			fsolutionwrite << setw(20) << Result.QrLP[i]*6e4;
			fsolutionwrite << setw(20) << Result.QHPtheo[i]*6e4;
			fsolutionwrite << setw(20) << Result.QLPtheo[i]*6e4 << "\n";
		};
        fsolutionwrite.close();
		
		//write pressure file for gap simulation
		fsolutionwrite.open("./inputs/pFile.dat");
		double step1;
		step1 = 360/0.2;
		int step = (int) step1;
		for(int i = tend - step; i < tend;i++)
		{
			fsolutionwrite << setw(20) << Result.t[i] - (Pump.nrevolutions-1)*60/(Pump.speed);
			fsolutionwrite << setw(20) << Result.pD[i];
			fsolutionwrite << setw(20) << Result.pHP[i];
			fsolutionwrite << setw(20) << Result.pLP[i] << "\n";
		};
        fsolutionwrite.close();
		cout << "Finished!" << "\n";

};