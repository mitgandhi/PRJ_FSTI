#include "CSlipperGap.h"
#include <omp.h>
#include <time.h>
#include <float.h>
#include <deque>
#pragma once
#define OMP

extern void matlab(const Array<double,2>& data,const string file);

CSlipperGap::CSlipperGap(caspar_input * gapinputs, CGapResult &Result) : gapinput(gapinputs), slipper(&GapLog), swashplate(&GapLog), t_slipper(&GapLog), t_swashplate(&GapLog)		//Constructor
{
	GapResult = &Result;
	
	GapLog.message("*-----------------------------------------------------------------------------*");
	GapLog.message("|                             SLIPPER GAP LOG FILE                            |");
	GapLog.message("*-----------------------------------------------------------------------------*");
	GapLog.message("\n");
	GapLog.message("Constructing object SlipperGap ... ");

	//Blitz++
	all = Range::all();

	//Grid parameters
	Fluid.M = gapinput->options_slipper.fluid_grid.M;
	Fluid.N = gapinput->options_slipper.fluid_grid.N;
	Fluid.Q = gapinput->options_slipper.fluid_grid.Q;
	Fluid.Groove1Location = gapinput->options_slipper.fluid_grid.Groove1Location;
	Fluid.Groove2Location = gapinput->options_slipper.fluid_grid.Groove2Location;
	Fluid.Groove1r = gapinput->options_slipper.fluid_grid.Groove1r;
	Fluid.Groove2r = gapinput->options_slipper.fluid_grid.Groove2r;
	Fluid.Groove1dr = gapinput->options_slipper.fluid_grid.Groove1dr;
	Fluid.Groove2dr = gapinput->options_slipper.fluid_grid.Groove2dr;

	//Geomerty parameters
	geometryslippergap.npistons = gapinput->operating_conditions.npistons;
	geometryslippergap.doutG = gapinput->geometry.doutG;							//[m]
	geometryslippergap.routG = geometryslippergap.doutG/2.0;						//radius [m]
	geometryslippergap.dinG = gapinput->geometry.dinG;								//[m]
	geometryslippergap.rinG = geometryslippergap.dinG/2.0;							//radius [m]
	geometryslippergap.mG = gapinput->geometry.mG;									//[kg]
	geometryslippergap.mK = gapinput->geometry.mK;									//[kg]
	geometryslippergap.lSG = gapinput->geometry.lSG;								//[m]
	geometryslippergap.lG = gapinput->geometry.lG;									//[m]
	geometryslippergap.rB = gapinput->geometry.dB/2.0;								//radius [m]
	geometryslippergap.dDG = gapinput->geometry.dDG;								//[m]
	geometryslippergap.lDG = gapinput->geometry.lDG;								//[m]
	geometryslippergap.vPocket = gapinput->geometry.vPocket;						//[m^3]
	geometryslippergap.hmaxG = gapinput->geometry.hmaxG;							//[m]
	geometryslippergap.Fslipper = gapinput->geometry.Fslipper;						//[N]
	geometryslippergap.alphaD_KG = gapinput->options_slipper.general.alphaD_KG;		//unitless orifice coefficient
	geometryslippergap.dDK = gapinput->geometry.dDK;								//[m]
	geometryslippergap.dK = gapinput->geometry.dK;									//[m]
	geometryslippergap.lDK = gapinput->geometry.lDK;								//[m]
	geometryslippergap.flowtype = gapinput->options_slipper.general.flowtype;
	//Added by Jeremy on 10.9.2015---------------//
	geometryslippergap.lKG = gapinput->geometry.lKG;								//[m];
	geometryslippergap.dZ = gapinput->geometry.dZ;									//[m];
	//-------------------------------------------//


	//Operating parameters
	operatingslippergap.omega = gapinput->operating_conditions.speed;
	operatingslippergap.speed = operatingslippergap.omega*(30.0/PI);
	operatingslippergap.beta_rad = gapinput->operating_conditions.beta;
	if(gapinput->operating_conditions.mode == 2)
	{
		operatingslippergap.beta_rad *= -1.0;
	}
	operatingslippergap.beta_deg = operatingslippergap.beta_rad*(180.0/PI);
	operatingslippergap.betamax_rad = gapinput->operating_conditions.betamax;
	operatingslippergap.betamax_deg = operatingslippergap.betamax_rad*(180.0/PI);
	operatingslippergap.gamma_rad = gapinput->geometry.gamma;
	operatingslippergap.speedK = gapinput->geometry.speedK;
	operatingslippergap.pHP = gapinput->operating_conditions.HP;
	operatingslippergap.pLP = gapinput->operating_conditions.LP;
	operatingslippergap.pcase = gapinput->operating_conditions.pCase;
	//Added by Jeremy on 10.9.2015---------------//
	operatingslippergap.T_Leak = gapinput->operating_conditions.T_Leak;
	//-------------------------------------------//

	//Fluid Mesh
	Fluid.dr.resize(Fluid.M,Fluid.N);	
	Fluid.dr = 0;

	//Presently a constant dtheta is used
	Fluid.dtheta.resize(Fluid.M,Fluid.N);	
	Fluid.dtheta = 0;

	//Boundary values [Pa / special]
	Fluid.boundary.resize(Fluid.M,Fluid.N);
	Fluid.boundary = -1;								//Default to solve Reynolds
	
	//Polar Coordinates of Fluid Grid [m, rad]
	Fluid.r.resize(Fluid.M,Fluid.N);
	Fluid.r = 0;

	//Set inner / outter radius
	Fluid.r(0, all) = geometryslippergap.dinG/2.0;
	Fluid.r(Fluid.M-1, all) = geometryslippergap.doutG/2.0;

	

	//Set non-uniform radial grid if applicable
	/*
	for(unsigned int i=0; i<GapInput.SlipperGapGrid.nonuniformr.size(); i++)
	{
		CGapInput::SlipperNonUniformGridR tmpr(GapInput.SlipperGapGrid.nonuniformr[i]);
		Fluid.r(tmpr.id, all) = tmpr.r/1000.0;
		Fluid.dr(tmpr.id, all) = tmpr.dr/1000.0;
	}
	*/

	if(Fluid.Groove1Location < 0 && Fluid.Groove2Location >= 0)
	{
		GapLog << CGapLog::error() << "Error: SlipperGrid Inputs. If Groove2Location is specified, Fluid.Groove1Location must be specified!" << endl;
	}

	if(Fluid.Groove1Location < 0 && Fluid.Groove2Location < 0)
	{
		//uniform dr
		Fluid.dr = 0.5*(geometryslippergap.doutG - geometryslippergap.dinG) / double(Fluid.M-1);
	} else if (Fluid.Groove2Location < 0)
	{
		//only one groove
		Fluid.dr(Range(0,Fluid.Groove1Location-1), all) = (Fluid.Groove1r-0.5*Fluid.Groove1dr - 0.5*geometryslippergap.dinG) / (double(Fluid.Groove1Location) - 1.5);
		Fluid.dr(Range(Fluid.Groove1Location+1,Fluid.M-1), all) = (0.5*geometryslippergap.doutG - (Fluid.Groove1r+0.5*Fluid.Groove1dr)) / (double(Fluid.M-Fluid.Groove1Location) - 2.5);
		Fluid.dr(Fluid.Groove1Location, all) = Fluid.Groove1dr - (Fluid.dr(Fluid.Groove1Location-1, all)+Fluid.dr(Fluid.Groove1Location+1, all));
	} else {
		//two grooves

		Fluid.dr(Range(0,Fluid.Groove1Location-1), all) = (Fluid.Groove1r-0.5*Fluid.Groove1dr - 0.5*geometryslippergap.dinG) / (double(Fluid.Groove1Location) - 1.5);
		Fluid.dr(Range(Fluid.Groove1Location+1,Fluid.Groove2Location-1), all) = ((Fluid.Groove2r-0.5*Fluid.Groove2dr) - (Fluid.Groove1r+0.5*Fluid.Groove1dr)) / double(Fluid.Groove2Location-Fluid.Groove1Location-3);
		Fluid.dr(Range(Fluid.Groove2Location+1,Fluid.M-1), all) = (0.5*geometryslippergap.doutG - (Fluid.Groove2r+0.5*Fluid.Groove2dr)) / (double(Fluid.M-Fluid.Groove2Location) - 2.5);

		Fluid.dr(Fluid.Groove1Location, all) = Fluid.Groove1dr - (Fluid.dr(Fluid.Groove1Location-1, all)+Fluid.dr(Fluid.Groove1Location+1, all));
		Fluid.dr(Fluid.Groove2Location, all) = Fluid.Groove2dr - (Fluid.dr(Fluid.Groove2Location-1, all)+Fluid.dr(Fluid.Groove2Location+1, all));
	}

	//Calculate the radius
	for(int i=1; i<Fluid.M; i++)
	{
		Fluid.r(i, all) = Fluid.r(i-1, all) + 0.5*(Fluid.dr(i-1, all)+Fluid.dr(i, all));
	}

	//Initialize theta
	Fluid.theta.resize(Fluid.M,Fluid.N);
	Fluid.theta = 0;

	Fluid.dtheta = 2*PI/Fluid.N;

	//Calculate theta
	for(int j=1; j<Fluid.N; j++)
	{
		Fluid.theta(all, j) = Fluid.theta(all, j-1) + 0.5*(Fluid.dtheta(all, j-1) + Fluid.dtheta(all, j));
	}
	//Fluid.theta = Fluid.dtheta*tensor::j;

	//Finite volume area of fluid grid
	{
		Fluid.dA.resize(Fluid.M,Fluid.N);
		Range R = Range(1,Fluid.M-1);
		Fluid.dA(R,all) = Fluid.r(R,all)*Fluid.dr(R,all)*Fluid.dtheta(R,all);
		
		//The first and last volumes are only half volumes
		Fluid.dA(0,all) = 0.5*Fluid.dtheta(0,all)*(
												pow(Fluid.r(0,all)+0.5*Fluid.dr(0,all),2.0)-pow(Fluid.r(0,all),2.0)
											);
		Fluid.dA(Fluid.M-1,all) = 0.5*Fluid.dtheta(Fluid.M-1,all)*(
												pow(Fluid.r(Fluid.M-1,all),2.0)-pow(Fluid.r(Fluid.M-1,all)-0.5*Fluid.dr(Fluid.M-1,all),2.0)
											);
	}
	
	//LCS Cartesian Coordinates of Fluid Grid [m]
	Fluid.Lx.resize(Fluid.M,Fluid.N);
	Fluid.Lx =(Fluid.r)*sin(Fluid.theta);
	Fluid.Ly.resize(Fluid.M,Fluid.N);
	Fluid.Ly =(Fluid.r)*cos(Fluid.theta);

	//Local rotating coordinate system
	Fluid.LRx.resize(Fluid.M,Fluid.N);
	Fluid.LRy.resize(Fluid.M,Fluid.N);

	//Global coordinate system
	Fluid.Gx.resize(Fluid.M,Fluid.N);
	Fluid.Gy.resize(Fluid.M,Fluid.N);

	//3D Gap Height Position Vector
	Fluid.z.resize(Fluid.M,Fluid.N,Fluid.Q);
	Fluid.z = 0;

	//Pressure 2D [Pa]
	Fluid.p.resize(Fluid.M,Fluid.N);
	Fluid.p_uncut.resize(Fluid.M,Fluid.N);
	Fluid.pold.resize(Fluid.M,Fluid.N);
	Fluid.pFullLoop.resize(Fluid.M,Fluid.N);
	Fluid.pold = 0;

	//Rigid Gap Height 2D [m]
	Fluid.hrigid.resize(Fluid.M,Fluid.N);
	Fluid.hgroove.resize(Fluid.M,Fluid.N);

	//Deformed Gap Height 2D [m]
	Fluid.h.resize(Fluid.M,Fluid.N);				

	//Squeeze Velocity 2D [m/s]
	Fluid.dht.resize(Fluid.M,Fluid.N);

	//Old ehdSqz
	Fluid.ehdsqzOld.resize(Fluid.M,Fluid.N);
	Fluid.ehdsqzOld = 0;

	//2D Fluid force on slipper [N]
	forcesslippergap.F.resize(Fluid.M,Fluid.N);
	forcesslippergap.F = 0;

	//3D Fluid Velocity
	Fluid.vr.resize(Fluid.M,Fluid.N,Fluid.Q);
	Fluid.vr = 0;
	Fluid.vtheta.resize(Fluid.M,Fluid.N,Fluid.Q);
	Fluid.vtheta = 0;

	//2D Boundary Slipper Velocity in Polar Coords
	operatingslippergap.vgr.resize(Fluid.M,Fluid.N);
	operatingslippergap.vgr = 0;
	operatingslippergap.vgtheta.resize(Fluid.M,Fluid.N);
	operatingslippergap.vgtheta = 0;

	//2D Boundary Slipper Velocity in Cartesian Coords (useful for plotting)
	operatingslippergap.vgx.resize(Fluid.M,Fluid.N);
	operatingslippergap.vgx = 0;
	operatingslippergap.vgy.resize(Fluid.M,Fluid.N);
	operatingslippergap.vgy = 0;

	//Slipper surface velocity
	operatingslippergap.tvr.resize(Fluid.M,Fluid.N);
	operatingslippergap.tvr = 0;
	operatingslippergap.tvtheta.resize(Fluid.M,Fluid.N);
	operatingslippergap.tvtheta = 0;

	//Swash/wobble plate velocity
	operatingslippergap.bvr.resize(Fluid.M,Fluid.N);
	operatingslippergap.bvr = 0;
	operatingslippergap.bvtheta.resize(Fluid.M,Fluid.N);
	operatingslippergap.bvtheta = 0;

	//Contact
	Fluid.contact.resize(Fluid.M,Fluid.N);
	Fluid.h_contact.resize(Fluid.M,Fluid.N);
	Fluid.p_contact.resize(Fluid.M,Fluid.N);
	Fluid.p_contact = 0;

	Fluid.alphaPold = 0.8;

	Fluid.ReyVals.resize(11,Fluid.M,Fluid.N);

	//Slipper EHD Deformation [m]
	slipper.ehd.resize(Fluid.M,Fluid.N);				
	slipper.oldehd.resize(Fluid.M,Fluid.N);				
	slipper.ehd = 0.0;
	slipper.oldehd = slipper.ehd;
	Fluid.newt_ehdsqzapprox.resize(Fluid.M);

	//Swashplate EHD Deformation [m]
	swashplate.ehd.resize(Fluid.M,Fluid.N);				
	swashplate.oldehd.resize(Fluid.M,Fluid.N);				
	swashplate.fakeehd.resize(Fluid.M,Fluid.N);				
	swashplate.ehd = 0.0;
	swashplate.fakeehd = 0.0;
	swashplate.oldehd = swashplate.ehd;

	//Slipper Thermal
	t_slipper.temp.resize(Fluid.M,Fluid.N);
	t_slipper.temp = 0;
	t_slipper.deform.resize(Fluid.M,Fluid.N);
	t_slipper.deform = 0;
	t_slipper.pdeform.resize(Fluid.M,Fluid.N);
	t_slipper.pdeform = 0;

	//Swashplate Thermal
	t_swashplate.temp.resize(Fluid.M,Fluid.N);
	t_swashplate.temp = 0;
	t_swashplate.deform.resize(Fluid.M,Fluid.N);
	t_swashplate.deform = 0;
	t_swashplate.pdeform.resize(Fluid.M,Fluid.N);
	t_swashplate.pdeform = 0;

	//Set EHD Squeeze Option
	if(gapinput->options_slipper.general.EHDsqueeze == 1 && gapinput->options_slipper.general.SlipperPressureDeformation != 1)
	{
		//Disable EHD squeeze if pressure deformation is not enabled
		gapinput->options_slipper.general.EHDsqueeze = 0;
	}

	//Macro Geometry
	if(gapinput->options_slipper.general.EnableSlipperMacro == 1)
	{
		geometryslippergap.SlipperMacro.resize(Fluid.M,Fluid.N);
		geometryslippergap.SlipperMacro = 0;
		
		ifstream macro(gapinput->options_slipper.general.SlipperMacroFile);
		if(!macro.is_open())
		{
			GapLog << CGapLog::error() << "Slipper Macro File: " << gapinput->options_slipper.general.SlipperMacroFile << " doesn't exist!" << endl;
		}

		int cnt = 0;
		while(!macro.eof())
		{
			int r1, r2, t1, t2;
			double h;
			macro >> r1;
			macro >> r2;
			macro >> t1;
			macro >> t2;
			macro >> h;
			//Scale from micron to meter
			h /= 1.0e+6;

			//Check for out of bound values
			if(r1 < 0 || r2 >= Fluid.M || t1 < 0 || t2 >= Fluid.N || r1 > r2 || t1 > t2)
			{
				GapLog.message("Error in the slipper macro geometry file! Cell values out of range! \n");
				GapLog.message("Compare fluid grid dimensions (options_slipper.txt -> M,N) ");
				GapLog.message("with the dimensions of the slipper macro geometry file. \n");
				exit(1);
			}

			//Set the macro to be whichever marco is higher
			//This allows for "nice" blends between steps
			geometryslippergap.SlipperMacro(Range(r1,r2),Range(t1,t2)) =  max(geometryslippergap.SlipperMacro(Range(r1,r2),Range(t1,t2)), h);
			cnt++;
		}
		macro.close();
	
		GapLog.message("Loaded " + n2s(cnt) + " slipper macro geometry values.");

	}

	//Influence Matricies
	if(gapinput->options_slipper.general.SlipperPressureDeformation == 1)
	{
		GapLog.message("Loading Slipper IM...");
		LoadSlipperIM();
	}
	if(gapinput->options_slipper.general.SwashplatePressureDeformation == 1)
	{
		GapLog.message("Loading Swashplate IM...");
		LoadSwashplateIM();

		//KD tree of the slipper fluid cell centroids for the full rotating group
		FullGroup.KDslip.points = annAllocPts(Fluid.M*Fluid.N*geometryslippergap.npistons,Fluid.KDslip.dim);

	}

	//KD Tree of the slipper cell centroids
	Fluid.KDslip.points = annAllocPts(Fluid.M*Fluid.N,Fluid.KDslip.dim);

	//Thermoelastic analysis
	if(gapinput->options_slipper.general.SlipperThermoElastic == 1)
	{
		//load the slipper thermal
		t_slipper.load("./input/thermal_slipper.txt");
	}

	if(gapinput->options_slipper.general.SwashplateThermoElastic)
	{
		//load the swashplate thermal
		t_swashplate.load("./input/thermal_swashplate.txt");
	}

	//Temperatures
	temperatureslippergap.T_HP = gapinput->operating_conditions.T_HP;
	temperatureslippergap.T_LP = gapinput->operating_conditions.T_LP;
	temperatureslippergap.T_Leak = gapinput->operating_conditions.T_Leak;

	//3D Fluid temperature [C] 
	Fluid.T.resize(Fluid.M,Fluid.N,Fluid.Q);				
	Fluid.T = temperatureslippergap.T_Leak;

	//Heat flux [W/m^2]
	Fluid.Qflux.resize(Fluid.M,Fluid.N);
	Fluid.Qflux = 0;

	//3D Fluid Viscosity
	Fluid.oilviscosity.resize(Fluid.M,Fluid.N,Fluid.Q);
	Fluid.oilviscosity = 0;

	//3D Fluid Density
	Fluid.oildensity.resize(Fluid.M,Fluid.N,Fluid.Q);
	Fluid.oildensity = 0;

	/*
	//Oil Parameters
	oilslippergap.oiltype = gapinput->oil.oiltype;
	oilslippergap.oildensity = gapinput->oil.oildensity;
	oilslippergap.oilviscosity = gapinput->oil.oilviscosity;
	oilslippergap.oilbetaP = gapinput->oil.oilbetaP;
	oilslippergap.oilbetaT = gapinput->oil.oilbetaT;
	oilslippergap.oilPc1 = gapinput->oil.oilPc1;
	oilslippergap.oilPc2 = gapinput->oil.oilPc2;
	oilslippergap.oilTc1 = gapinput->oil.oilTc1;
	oilslippergap.oilTc2 = gapinput->oil.oilTc2;
	oilslippergap.oilW = gapinput->oil.oilW;
	oilslippergap.alpha1 = gapinput->oil.alpha1;
	oilslippergap.alpha2 = gapinput->oil.alpha2;
	oilslippergap.alpha3 = gapinput->oil.alpha3;
	oilslippergap.oilC = gapinput->oil.oilC;
	oilslippergap.oillambda = gapinput->oil.oillambda;
	oilslippergap.oilK = gapinput->oil.oilK;
	oilslippergap.oilbetaKP = gapinput->oil.oilbetaKP;
	oilslippergap.oilbetaKT = gapinput->oil.oilbetaKT;
	*/
	//Oil
	if(gapinput->oil.general.oiltype == 0)
	{
		oil_properties = new constant_oil(*gapinput->parent_input);		// use constant properties
	}
	else if(gapinput->oil.general.oiltype == 1)
	{
		oil_properties = new user_defined_oil(*gapinput->parent_input);	// user defined
	}
	else if(gapinput->oil.general.oiltype == 2)
	{
		oil_properties = new HLP32_oil(*gapinput->parent_input);		// HLP32
	}
	else if(gapinput->oil.general.oiltype == 3)
	{
		oil_properties = new user_defined2_oil(*gapinput->parent_input);		// userdefined2
	}/*
	else if(gapinput->oil.general.oiltype == 3)
	{
		oil_properties = new skydrol_oil(*gapinput->parent_input);		// skydrol
	}
	else if(gapinput->oil.general.oiltype == 4)
	{
		oil_properties = new red_oil(*gapinput->parent_input);		// red oil
	}
	else if(gapinput->oil.general.oiltype == 6)
	{
		oil_properties = new milh5606_oil(*gapinput->parent_input);		// mil H5606
	}
	else if(gapinput->oil.general.oiltype == 7)
	{
		oil_properties = new ExxonDTE10Excel32(*gapinput->parent_input);		// exxon
	}
	else if(gapinput->oil.general.oiltype == 8)
	{
		oil_properties = new ISO46_oil(*gapinput->parent_input);		// iso 46
	}
	else if(gapinput->oil.general.oiltype == 9)
	{
		oil_properties = new Parker_skydrol_oil(*gapinput->parent_input);		// Parker Skydrol
	}
	else if(gapinput->oil.general.oiltype == 50)
	{
		oil_properties = new water_oil(*gapinput->parent_input);		// water
	}*/
	else
	{
		GapLog << "ERROR: Oil type " << gapinput->oil.general.oiltype << " not supported" << endl;
		exit(1);
	}

	//Results
	
	//2D Fluid Pressure [Pa]
	GapResult->SlipperResult.P_slipper.resize(Fluid.M,Fluid.N);
	GapResult->SlipperResult.P_slipper=0;

	GapResult->SlipperResult.P_ehdsq.resize(Fluid.M,Fluid.N);
	GapResult->SlipperResult.P_ehdsq=0;

	//2D Deformed Gap Height [m]
	GapResult->SlipperResult.h_slipper.resize(Fluid.M,Fluid.N);
	GapResult->SlipperResult.h_slipper=0;

	//2D EHD [m]
	GapResult->SlipperResult.slipperEHD.resize(Fluid.M,Fluid.N);
	GapResult->SlipperResult.slipperEHD=0;
	GapResult->SlipperResult.swashplateEHD.resize(Fluid.M,Fluid.N);
	GapResult->SlipperResult.swashplateEHD=0;

	//3D Fluid Velocity in Polar Coords
	GapResult->SlipperResult.vr_slipper.resize(Fluid.M,Fluid.N,Fluid.Q);
	GapResult->SlipperResult.vr_slipper=0;
	GapResult->SlipperResult.vtheta_slipper.resize(Fluid.M,Fluid.N,Fluid.Q);
	GapResult->SlipperResult.vtheta_slipper=0;

	//2D Slipper Velocity in Rec Coords
	GapResult->SlipperResult.vgx.resize(Fluid.M,Fluid.N);
	GapResult->SlipperResult.vgx=0;
	GapResult->SlipperResult.vgy.resize(Fluid.M,Fluid.N);
	GapResult->SlipperResult.vgy=0;

	//3D Fluid Temperature
	GapResult->SlipperResult.T_slipper.resize(Fluid.M,Fluid.N,Fluid.Q);
	GapResult->SlipperResult.T_slipper=0;

	//3 Point Forces [N]
	forcesslippergap.F_contact.resize(3);
	GapResult->SlipperResult.F_contact.resize(3);
	GapResult->SlipperResult.F_fluid.resize(3);
	GapResult->SlipperResult.F_external.resize(3);
	GapResult->SlipperResult.dFcomp.resize(3);

	GapResult->SlipperResult.Sock_fixed = false;

	//Hard coded min gap height
	Fluid.ContactHeight = gapinput->options_slipper.general.RoughnessRq * 1.0e-6;
	if(Fluid.ContactHeight < 0.05e-6 || Fluid.ContactHeight > 2.0e-6)
	{
		Fluid.ContactHeight = 0.1e-6;
	}

	Fluid.ap.resize(Fluid.M,Fluid.N);
	Fluid.b.resize(Fluid.M,Fluid.N);
	
	GapLog.message("SlipperGap constructed successfully.");

	//Create the Lx, Ly coords
	matlab(Fluid.Lx, "./output/slipper/matlab/Grx.txt");
	matlab(Fluid.Ly, "./output/slipper/matlab/Gry.txt");
	matlab(Fluid.dA, "./output/slipper/matlab/dA.slipper");
}
CSlipperGap::~CSlipperGap()	//Destructor
{
	if(gapinput->lubrication_module.solve_slipper == 0)
	{
		//the slipper is not calculated so don't do anything
		return;
	}

	//Clean up influence matricies
	/*
	for(int i=0; i<Fluid.M; i++)
	{
		delete [] Solid.IM[i];
	}
	delete [] Solid.IM;
	delete [] Solid.IMpocket;
	*/

	GapLog.message("Object CSlipperGap destructed successfully!");
}
void CSlipperGap::SlipperGap(vector<double> &xg,vector<double> &vg,vector<double> &dF,int FullLoop)
{
	//meaning of FullLoop
	//FullLoop == 0 : Do NOT solve the FSI problem - no iterations
	//FullLoop == 1 : Solve the full FSI problem for the 1st time this timestep
	//FullLoop == 2 : Solve the FSI problem but it has already been solved at least once this timestep
	curFullLoop = FullLoop;

	if(!gapinput->options_slipper.general.Explicit)
	{
		//Implicit Method
		//Calculate the position using the 1st order backward Euler method
		//Assume constant acceleration between timesteps
		for(int i = 0; i < 3; i++)
		{
			//xg[i] = pold[i]+gapinput->common.timestep*(vold[i]+vg[i])/2.0;
			xg[i] = pold[i]+gapinput->common.timestep*vg[i];
		}
	}

	//Just for testing purposed .. not to be used in production code
	//testDeform();
	//testDeform(xg, vg);
	//debug_makeEHDfluidpressurefile();

	//Set the current control point position and velocities
	control_point_position = xg;
	control_point_velocity = vg;

	//Calculate the 2D slipper gap height and shifting velocity using the three control points
	if(FullLoop == 1)
	{
		//Under relax the new position in a full FSI loop
		SlipperCalch(pold);
	} else {
		SlipperCalch();
	}

	//Calculate the shifting velocity
	SlipperCalcdht();

	//Logging
	int lastRprint = INT_MIN;
	if(FullLoop == 1)
	{
		GapLog.message("\t****** Initial Picard Iteration ******");
	}

	//Debug log
	ofstream dbglog;
	if(gapinput->options_slipper.general.DebugMode == 1 && FullLoop >= 0)
	{
		dbglog.open("./output/slipper/SlipperFxtPtLog.txt");
		dbglog << "%Slipper Debug Mode FullLoop Convergence..." << endl;
		dbglog << 
			"%x0: " << xg[0] << "\t" <<
			"x1: " << xg[1] << "\t" <<
			"x2: " << xg[2] << "\t" <<
			"v0: " << vg[0] << "\t" <<
			"v1: " << vg[1] << "\t" <<
			"v2: " << vg[2] << "\t" <<
			"t: " << gapinput->common.time <<
		endl;
		dbglog << "%" << endl;
		dbglog <<
			"%cnt\t" <<
			"pG\t\t" <<
			"pSocket\t\t" <<
			"min(h)\t\t" <<
			"mean(mu)\t\t" <<
			"resid\t\t" <<
			"aPold\t\t" <<
			"aEHD\t" <<
			"apG\t\t" <<
		endl;
		dbglog << "%-------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	}


	//Under-relaxation for pocket pressure
	double AlphapG = gapinput->options_slipper.numeric.AlphapG;
	
	//Under-relaxation for pressure field (This is the main FSI relaxation)
	double alphaPold = gapinput->options_slipper.numeric.AlphaReynolds;

	if(FullLoop < 1)
	{
		//This _MUST_ be forced to 1.0 if not in fullloop
		alphaPold = 1.0;
	} else if (FullLoop == 1)
	{
		alphaPold = Fluid.alphaPold;
		Fluid.p_contact = 0;

	} else {
		alphaPold = Fluid.alphaPoldCP;
		Fluid.p_contact = 0;
	}

	//A sum of each alphaPold used to construct an average
	double alphaPoldSum = 0;

	//alphaEHD is used to slowly bring the ehd from the previous FSI solution into this solution
	double alphaEHD; 
	if(sum(swashplate.oldehd) == 0 || FullLoop != 1)
	{
		//Disable if not in full loop or if there is no old swashplate deformation
		alphaEHD = 1.0;
	} else {
		//Start with a small value
		alphaEHD = 0.01;
	}

	vector<double> dpg;
	size_t last_decpg = 0;

	vector<double> dresid;
	size_t last_decresid = 0;
	
	//The start residual of the current loop	
	double resid = 0;
	
	//The end residual of the current loop	
	double Fresid = 0;

	//The old pG & pGerror
	double pGold = operatingslippergap.pG;
	double pGerror = fabs((pGold-operatingslippergap.pG)/operatingslippergap.pG);

	//Total number of fixed point loops
	int slippercnt = 0;
	
	//Used for alpha control
	int goodctr = 0;
	
	if(FullLoop > 0 && gapinput->options_slipper.general.DenseMode == 1)
	{
		matlab(Fluid.boundary,"./dense/boundary.txt");
		matlab(Fluid.p,"./dense/p-1.txt");
		matlab(Fluid.hrigid,"./dense/hrigid-1.txt");
		matlab(slipper.ehd,"./dense/slipper-1.txt");
		matlab(swashplate.ehd,"./dense/swashplate-1.txt");
		matlab(Fluid.h,"./dense/h-1.txt");
		matlab(Fluid.hrigid,"./dense/hrigid-1.txt");
		//writevtk("./dense/Initial-FluidGrid.vtk");
	}

	//Time the fixed point iteration
	clock_t start = clock();
	clock_t trap = clock();
	double Rtime = 0;
	double Etime = 0;
	double Swashtime = 0;
	double Sliptime = 0;
	double Contacttime = 0;

	int maxFullLoop = 500;
	if(gapinput->options_slipper.general.ComplexPicard == 2)
	{
		maxFullLoop = 100;
	}

	/* ------------------------------------------------------- */
	// Begin Main FSI Fixed Point Loop
	/* ------------------------------------------------------- */
	do
	{
		//The old residual
		double oldresid = resid;

		//The residual at the start of this current loop
		if(gapinput->options_slipper.general.EnableRoughness == 1)
		{
			resid = presidPC();
		} else {
			resid = presid();
		}

		//Fixed point control
		if (FullLoop > 0)		
		{
			//Logging progress to cout
			if(lastRprint+50 <= slippercnt && FullLoop == 1)
			{

				ostringstream oss (ostringstream::out);
				oss.str("");
				oss << setprecision (5) <<  scientific << resid;

				GapLog.message("\tResidual: " + oss.str() + "\tCnt: " + n2s(slippercnt));
				lastRprint = slippercnt;
			}

			double minH = min(where(Fluid.hgroove == 0,Fluid.h,1));
			
			unsigned int switch_cnt = 0;
			for(size_t i=last_decresid+1; i<dresid.size(); i++)
			{
				if((dresid[i-1] < 0 && dresid[i] > 0) || (dresid[i-1] > 0 && dresid[i] < 0))
				{
					switch_cnt++;
				}
			}

			//if(slippercnt > 0 && oldresid < resid)
			if(slippercnt > 0 && switch_cnt > 1)
			{	
				last_decresid = dresid.size();

				//We don't have good convergence
				goodctr = 0;
				
		//		if(	alphaEHD >= 1 || 
		//				(alphaPold > 0.5 && minH >= 1.0e-6) ||
		//				(alphaEHD >= 0.8 && minH < 0.8e-6) 
		//			)
				{
					if(alphaPold >= 0.14)
					{
						//Decrease our palpha
						alphaPold -= 0.05;
					} else {
						alphaPold -= 0.005;
					}

					//Limit the palpha
					if(minH > 0.5e-6 && alphaPold < 0.01)
					{
						alphaPold = 0.01;
					} else if ( alphaPold < 0.001 )
					{
						alphaPold = 0.001;
					}
				}

			}

			if(goodctr > 15 || (alphaPold < 0.1 && goodctr > 5) )
			{
				//We've had good convergence for a while

				if( alphaPold >= 0.1)
				{
					alphaPold += 0.025;
				} else {
					alphaPold += 0.0025;
				}

				//Limit the alpha
				if(alphaPold >= 0.9)
				{
					alphaPold = 0.9;
				}
			
				//Reset the goodctr
				goodctr = 0;
			}
		
			//Increase the goodctr
			goodctr++;

			//Increase the alphaEHD
			alphaEHD += 0.02;
			if(alphaEHD > 1.0)
			{
				alphaEHD = 1.0;
			}

		}
		//End fixed point control

		//Relax the rigid height according to alphaEHD
		vector<double> xgRelax(3);
		for(int i=0; i<3; i++)
		{
			xgRelax[i] = pold[i] - alphaEHD * (pold[i] - xg[i]);
		}
		SlipperCalch(xgRelax);

		//cut the gap height
		const double wearHeight = 2.0*Fluid.ContactHeight;
		Fluid.h = where(Fluid.h < wearHeight, wearHeight, Fluid.h);
		
		trap = clock();
		//Solve Reynolds
		int nRiterations;
		if(gapinput->options_slipper.general.EnableRoughness == 1)
		{
			nRiterations = SlipperReynoldsPC(alphaPold);
		} else {
			nRiterations = SlipperReynolds(alphaPold);
			//nRiterations = SlipperReynoldsGS(alphaPold);
		}	
		Rtime += clock()-trap;
		
		/*
		//Update the socket pressure assuming it is = to fluid pressure
		{
			double fsocket = sum((Fluid.p-operatingslippergap.pcase)*Fluid.dA)+
									PI*(pow(geometryslippergap.dinG/2.0,2)-pow(geometryslippergap.dDG/2.0,2))*(operatingslippergap.pG-operatingslippergap.pcase);
			double AreaSock = GapInput.Geometry.aSock * 1e-6; //convert to m^2
			forcesslippergap.psocket = fsocket / AreaSock;
		}
		*/

		//Update fluid velocities
		SlipperCalcFluidV();
		
		//Calculate fluid density
		SlipperCalcDensity();
		
		//Update pocket pressure
		//SlipperpG(AlphapG);	
		{
			{
				unsigned int switch_cnt = 0;
				for(size_t i=last_decpg+1; i<dpg.size(); i++)
				{
					if((dpg[i-1] < 0 && dpg[i] > 0) || (dpg[i-1] > 0 && dpg[i] < 0))
					{
						switch_cnt++;
					}
				}

				if(switch_cnt >= 2)
				{
					AlphapG = 0.75*AlphapG;
					if(AlphapG < 0.05)
					{
						AlphapG = 0.05;
					}

					last_decpg = dpg.size();
				} 
					
				else if (dpg.size() - last_decpg > 10 && (dpg.size() - last_decpg) % 20 == 0)
				{
					AlphapG = 1.05*AlphapG;
					if(AlphapG > gapinput->options_slipper.numeric.AlphapG)
					{
						AlphapG = gapinput->options_slipper.numeric.AlphapG;
					}
				}
					
			}
				
			double alphapg = AlphapG;

			double pgOld = operatingslippergap.pG;
			
			if(slippercnt < 10 && FullLoop > 0)
			{
				operatingslippergap.pG = operatingslippergap.pGprevious;

				 //Update pressure boundary values
				Fluid.p = where(Fluid.boundary >= 0, Fluid.boundary, Fluid.p);
				Fluid.p_uncut = where(Fluid.boundary >= 0, Fluid.boundary, Fluid.p_uncut);
				Fluid.p = where(Fluid.boundary == -2, operatingslippergap.pG, Fluid.p);
				Fluid.p_uncut = where(Fluid.boundary == -2, operatingslippergap.pG, Fluid.p_uncut);
				Fluid.p = where(Fluid.boundary == -3, operatingslippergap.pcase, Fluid.p);
				Fluid.p_uncut = where(Fluid.boundary == -3, operatingslippergap.pcase, Fluid.p_uncut);
			} else {
				SlipperpG(alphapg);
			}

			dpg.push_back(operatingslippergap.pG - pgOld);
		}	

		//Pocket pressure error
		pGerror = fabs(pGold-operatingslippergap.pG);
		
		/* ------------------------------------------------------- */
		// Full Loop Calculations
		/* ------------------------------------------------------- */
		if(FullLoop > 0)
		{
			if((slippercnt < 10 || slippercnt % 10 == 0) && FullLoop == 1) //Used to reduce computational effort
			{
				if(gapinput->options_slipper.general.CalcEnergy)
				{
					trap = clock();
					//Solve Energy
					SlipperEnergy();
					//SlipperEnergyCG();
					Etime += clock()-trap;
				}

				//Calculate viscosity
				SlipperCalcViscosity(1.0);	//No relaxation
			}

			//Calculate swashplate pressure deformation
			if(gapinput->options_slipper.general.SwashplatePressureDeformation == 1)
			{
				trap = clock();
				SwashplateCalcEHD(alphaEHD);
				Swashtime += clock()-trap;
			}

			//Calculate slipper pressure deformation
			if(gapinput->options_slipper.general.SlipperPressureDeformation == 1)
			{
				trap = clock();
				SlipperCalcEHD(alphaEHD);
				Sliptime += clock()-trap;

				if(gapinput->options_slipper.general.EHDsqueeze == 1)
				{
					//Update the ehd squeeze term
					SlipperCalcdht();
				}
				Fluid.ehdsqzOld = Fluid.dht;
				
			}

		} 
		/* ------------------------------------------------------- */
		// End Full Loop Calculations
		/* ------------------------------------------------------- */

		if(_isnan(sum(Fluid.p)))
		{
			//something has gone very wrong with the solution
			GapLog << "WARNING: Fluid pressure has reported NaN value! Resetting..." << endl;
			
			//a long list of resets, quickly hacked together, should be checked
			operatingslippergap.pG = operatingslippergap.pDC;
			operatingslippergap.pGprevious = operatingslippergap.pG;
			Fluid.p_contact = 0;
			SlipperPressureBounds();
			Fluid.pFullLoop = Fluid.p;
			Fluid.pold = Fluid.p;
			Fluid.p_uncut = Fluid.p;
			if(gapinput->options_slipper.general.SlipperPressureDeformation == 1)
			{
				slipper.ehd = 0;
				slipper.oldehd = 0;
				SlipperInitEHD();
				slipper.oldehd = slipper.ehd;
			}
			if(gapinput->options_slipper.general.SwashplatePressureDeformation == 1)
			{
				swashplate.ehd = 0;
				swashplate.oldehd = 0;
				SwashplateInitEHD(1.0);
				swashplate.oldehd = swashplate.ehd;
			}
			SlipperInitializeTemperature();
			SlipperCalcViscosity(1.0);
			Fluid.alphaPold = 0.1;
			Fluid.alphaPoldCP = Fluid.alphaPold;
			SlipperCalch();
			SlipperCalcdht();

			double lifth = 2.0e-6-min(Fluid.h);
			if(lifth > 15e-6)
			{
				lifth = 15e-6;	//limit excessive values
			}
			if(lifth >= 0)
			{
				for(int i=0; i<3; i++)
				{
					//don't think this would work correctly in explicit ode mode
					control_point_position[i] += lifth;
					xg[i] += lifth;
				}
			}
			SlipperCalch();

			if(_isnan(sum(Fluid.p)))
			{
				GapLog << "WARNING: Fluid.p still NaN" << endl;
			}

			if(_isnan(sum(Fluid.h)))
			{
				GapLog << "WARNING: Fluid.h still NaN" << endl;
			}

			if(_isnan(sum(Fluid.hrigid)))
			{
				GapLog << "WARNING: Fluid.hrigid still NaN" << endl;
			}

			if(_isnan(sum(slipper.ehd)))
			{
				GapLog << "WARNING: slipper.ehd still NaN" << endl;
			}

			if(_isnan(sum(swashplate.ehd)))
			{
				GapLog << "WARNING: swashplate.ehd still NaN" << endl;
			}

			if(_isnan(sum(t_slipper.pdeform)))
			{
				GapLog << "WARNING: t_slipper.pdeform still NaN" << endl;
			}
			
			if(_isnan(sum(geometryslippergap.SlipperMacro)))
			{
				GapLog << "WARNING: geometryslippergap.SlipperMacro still NaN" << endl;
			}

		}

		if(gapinput->options_slipper.general.DebugMode == 1)
		{
				dbglog << slippercnt << setprecision (5) <<  scientific << "\t" <<
						operatingslippergap.pG << "\t" << 
						forcesslippergap.psocket << "\t" <<
					   min(where(Fluid.hgroove == 0,Fluid.h,1)) << "\t" <<
						mean(Fluid.oilviscosity) << "\t" <<
						resid << "\t" <<
						alphaPold << "\t" <<
						alphaEHD << "\t" <<
						AlphapG << "\t" <<
					endl;
		}

		//Calculate contact pressure
		{
			if(FullLoop == 0)
			{
				//the height was cut before Reynolds
				SlipperCalch();
			}
			const double ContactHeight = 2.0*Fluid.ContactHeight;
			Fluid.contact = where(Fluid.h < ContactHeight , 1, 0);
			Fluid.h_contact = where(Fluid.h < ContactHeight , ContactHeight-Fluid.h, 0);
			//Fluid.h = where(Fluid.h < ContactHeight , ContactHeight , Fluid.h);

			trap = clock();
			SlipperCalcContactPressure(alphaEHD, FullLoop, alphaPold);
			Contacttime += clock()-trap;
		}

		if(FullLoop > 0 && gapinput->options_slipper.general.DebugMode == 1)
		{
			/* ------------------------------------------------------- */
			// Dense Mode Logging
			/* ------------------------------------------------------- */
			if(gapinput->options_slipper.general.DenseMode == 1)
			{
				/*
				matlab(Fluid.p,"./dense/p" + n2s(slippercnt) + ".txt");
				matlab(Fluid.h,"./dense/h" + n2s(slippercnt) + ".txt");
				//matlab(Fluid.hrigid,"./dense/hrigid" + n2s(slippercnt) + ".txt");
				matlab(Fluid.p_uncut,"./dense/puncut" + n2s(slippercnt) + ".txt");
				//matlab(Fluid.h_contact,"./dense/hc" + n2s(slippercnt) + ".txt");
				matlab(Fluid.p_contact,"./dense/pc" + n2s(slippercnt) + ".txt");
				matlab(slipper.ehd,"./dense/slipper" + n2s(slippercnt) + ".txt");
				matlab(swashplate.ehd,"./dense/swashplate" + n2s(slippercnt) + ".txt");
				*/
				writevtk("./dense/vtk" + n2s(slippercnt) + ".vtk");
			}
		}

		//Increment fixed point itteration loop counter
		slippercnt++;
		
		//This will be used to take a 'mean' of the alpha during the FSI loop
		alphaPoldSum += alphaPold;

		//The residual at the end of this current loop
		if(gapinput->options_slipper.general.EnableRoughness == 1)
		{
			Fresid = presidPC();
		} else {
			Fresid = presid();
		}

		dresid.push_back(Fresid-resid);

		//The old pocket pressure
		pGold = operatingslippergap.pG;

	} while ( (Fresid > gapinput->options_slipper.numeric.FSIresidTol || pGerror > 0.05e+5|| alphaEHD < 1.0 || slippercnt < 10)
					&& slippercnt < maxFullLoop && FullLoop >= 0);	//Convergence criteria
	/* ------------------------------------------------------- */
	// End Main FSI Fixed Point Loop
	/* ------------------------------------------------------- */


	/* ------------------------------------------------------- */
	// Debug Mode
	/* ------------------------------------------------------- */
	if(gapinput->options_slipper.general.DebugMode == 1 && FullLoop == 1)
	{
		dbglog.close();
		//writevtk("./vtk/" + n2s(operatingslippergap.phi_deg) + ".vtk");
	}
	/* ------------------------------------------------------- */
	// End Debug Mode
	/* ------------------------------------------------------- */

	if(FullLoop == 1)
	{
			GapLog.message(
			"\tTotal time: " + n2s(double(clock()-start)/double(CLOCKS_PER_SEC)) + " [s]" + 
			"\tTotal Cnt: " + n2s(slippercnt) +
			"\tRtime: " + n2s(Rtime/double(CLOCKS_PER_SEC)) +
			"\tEtime: " + n2s(Etime/double(CLOCKS_PER_SEC)) +
			"\tSliptime: " + n2s(Sliptime/double(CLOCKS_PER_SEC)) +
			"\tSwashtime: " + n2s(Swashtime/double(CLOCKS_PER_SEC)) +
			"\tContacttime: " + n2s(Contacttime/double(CLOCKS_PER_SEC))
			);

			Fluid.pFullLoop = Fluid.p;
	} else if (FullLoop == 2)
	{
		GapLog.message(
			n2s(slippercnt) + " ");
	} else if (FullLoop == 3)
	{
		GapLog.message(
			n2s(slippercnt) + " ");
	}

	//Calculate the fluid forces due to gap pressure field
	SlipperCalcFluidForces();

	//Calculate external forces
	SlipperCalcExternalForces();

	//Calculate force balance
	SlipperCalcdF(dF, xg);

	//Assign Results
	SlipperAssignResults(xg);

	if(FullLoop == 1)
	{
		Fluid.alphaPold = alphaPoldSum / (double) slippercnt;
		Fluid.alphaPoldCP = Fluid.alphaPold;
	} else if (FullLoop > 1)
	{
		Fluid.alphaPoldCP = alphaPoldSum / (double) slippercnt;
	}

	if(FullLoop > 0)
	{
		Fluid.pold = Fluid.p;
	}

	if(gapinput->options_slipper.general.DenseMode == 1 && FullLoop > 0)
	{
		writevtk("./dense/Final-FluidGrid.vtk");
		if(gapinput->options_slipper.general.SlipperPressureDeformation == 1)
		{
			SlipperVTK("./dense/Final-SlipperGrid.vtk");
		}

		if(FullLoop != 3)
		{
			GapLog << "fx mx my: " << GapResult->SlipperResult.dFcomp[0] << "\t"
										<< GapResult->SlipperResult.dFcomp[1] << "\t"
										<< GapResult->SlipperResult.dFcomp[2] << "\t"
										<< endl;
			GapLog << "p0 p1 p2: " << forcesslippergap.dF[0] << "\t"
										<< forcesslippergap.dF[1] << "\t"
										<< forcesslippergap.dF[2] << "\t"
										<< endl;

	
			string junk;
			cout << "Press Enter to continue to Ctrl-C to quit...";
			cin >> junk;
			cout << endl;
		}
		system("del /f/q .\\dense\\");

	}

}
void CSlipperGap::SlipperStartTimestep(vector<double> &Newtonpositionslipper, vector<double> &Newtonvelocityslipper, vector<double> &oldPosition, vector<double> &oldVelocity)
{
	//Update the old positions / velocities
	control_point_position = pold = oldPosition;
	control_point_velocity = vold = oldVelocity;

	//Set pressures
	operatingslippergap.pHP = gapinput->common.pHP;
	operatingslippergap.pLP = gapinput->common.pLP;
	operatingslippergap.pDC = gapinput->common.pDC;
	
	//Phi
	operatingslippergap.phi_deg = gapinput->common.phi_deg; //[deg] between 0 and 360
	operatingslippergap.phi_rad = operatingslippergap.phi_deg*PI/180; //[rad] between 0 and 2PI
	
	//Get / zero the piston friction force
	GetFTK();

	//Assign the piston reaction force
	GetFSK();

	//Find slipper (r,theta) velocities - pg 37-38 of Wieczorek thesis
	double re;
	Array<double,2> rm(Fluid.M,Fluid.N);
	Array<double,2> delta(Fluid.M,Fluid.N);
	Array<double,2> angle(Fluid.M,Fluid.N);
	
	re = pow(pow(geometryslippergap.rB*sin(operatingslippergap.phi_rad),2)+pow(geometryslippergap.rB/cos(operatingslippergap.beta_rad)*cos(operatingslippergap.phi_rad),2),0.5);  //ellipse radius
	rm = pow(Fluid.r*Fluid.r+re*re-2*Fluid.r*re*cos(PI-Fluid.dtheta*tensor::j),0.5);
	angle = (rm*rm+Fluid.r*Fluid.r-re*re)/(2*Fluid.r*rm);		
		angle = where(angle<1,angle,1);			//used to correct for numerical error when angle is slightly > 1
		angle = where(angle>-1,angle,-1);		//used to correct for numerical error when angle is slightly < -1
	delta = PI/2-acos(angle);
	
	//operatingslippergap.vgr = -operatingslippergap.omega*rm*cos(delta)*where(tensor::j<Fluid.N/2,1,-1);
	//operatingslippergap.vgtheta = operatingslippergap.speedK*operatingslippergap.omega*Fluid.r-operatingslippergap.omega*rm*sin(delta);

	operatingslippergap.vgr = operatingslippergap.omega*rm*cos(delta)*where(tensor::j<Fluid.N/2,1,-1);
	operatingslippergap.vgtheta = -operatingslippergap.speedK*operatingslippergap.omega*Fluid.r+operatingslippergap.omega*rm*sin(delta);
	
	//operatingslippergap.vgr = operatingslippergap.omega*re*sin(Fluid.theta);
	//operatingslippergap.vgtheta = -operatingslippergap.speedK*operatingslippergap.omega*Fluid.r+operatingslippergap.omega*re*cos(Fluid.theta);
	
	operatingslippergap.vgx = operatingslippergap.vgr * sin(Fluid.theta) + operatingslippergap.vgtheta * cos(Fluid.theta);
	operatingslippergap.vgy = operatingslippergap.vgr * cos(Fluid.theta) - operatingslippergap.vgtheta * sin(Fluid.theta);

	//New model method - top and bottom surface velocity
	operatingslippergap.tvr = operatingslippergap.vgr;
	operatingslippergap.tvtheta = operatingslippergap.vgtheta;

	operatingslippergap.bvr = 0;
	operatingslippergap.bvtheta = 0;

	//Initialize pG to pDC
	operatingslippergap.pGprevious = operatingslippergap.pG;
	operatingslippergap.pG = operatingslippergap.pDC;

	//to be safe, especially in some resume'd case
	Fluid.p_contact = 0;

	if(gapinput->common.first_step)
	{
		operatingslippergap.pGprevious = operatingslippergap.pG;
	}
	
	//Update the coordinate systems
	SlipperSetCoordSys();

	if(gapinput->options_slipper.general.SwashplatePressureDeformation == 1)
	{
		//Because the phi angle has changed we need to update which slippers are in the FullGroup

		int steps = gapinput->common.rev_steps;
		int step = int(floor((gapinput->common.phi_rev_deg+gapinput->common.phi_deg_tol) / gapinput->common.phistep_deg)) % steps;

		for(int np = 0, cnt = 0; np<geometryslippergap.npistons; np++)
		{
			int s = (step + np*steps/geometryslippergap.npistons) % steps;
			for(int i=0; i<Fluid.M; i++)
			{
				for(int j=0; j<Fluid.N; j++, cnt++)
				{
					FullGroup.KDslip.points[cnt][0] = FullGroup.Gx(i,j,s);
					FullGroup.KDslip.points[cnt][1] = FullGroup.Gy(i,j,s);
				}
			}
		}

		if(FullGroup.KDslip.kdtree == NULL)
		{
			FullGroup.KDslip.kdtree = new ANNkd_tree(FullGroup.KDslip.points, Fluid.M*Fluid.N*geometryslippergap.npistons, FullGroup.KDslip.dim);
		} else {
			delete FullGroup.KDslip.kdtree;
			FullGroup.KDslip.kdtree = new ANNkd_tree(FullGroup.KDslip.points, Fluid.M*Fluid.N*geometryslippergap.npistons, FullGroup.KDslip.dim);
		}
	}

	//Define the pressure boundaries (note that the KD tree is required for theta bounds)	
	SlipperDefinePressureBounds();

	if(sum(Fluid.pold) == 0)
	{
		//Pressure bounds need to be defined first
		//This is needed to initialize temperature (and thus visco) and deformation below
		SlipperPressureBounds();

		//Initialize pFullLoop to the interpolate pressure since we don't have anything better to use
		Fluid.pFullLoop = Fluid.p;
	} else {
		//Use the old FSI pressure field
		Fluid.p = Fluid.pold;

		//Initialize pFullLoop to the old final pressure
		//This will be used by the preNewt EHDSQZ
		Fluid.pFullLoop = Fluid.p;
	
		//Update pressure boundary values
		Fluid.p = where(Fluid.boundary >= 0, Fluid.boundary, Fluid.p);
		Fluid.p = where(Fluid.boundary == -2, operatingslippergap.pG, Fluid.p);
		Fluid.p = where(Fluid.boundary == -3, operatingslippergap.pcase, Fluid.p);

		//Set the uncut field to this field
		Fluid.p_uncut = Fluid.p;
	}

	if(gapinput->options_slipper.general.SlipperPressureDeformation == 1)
	{	
		static bool newt_ehdsqz_built = false;
		if(gapinput->options_slipper.general.EHDsqueeze && !newt_ehdsqz_built)
		{
			//generate the sqz approx for newton
			GapLog << "Building Newton EHD Squeeze IM Approximation ..." << endl;
			
			Array<double, 2> pbackup(Fluid.p.copy());
			Array<double, 2> ehdbackup(slipper.ehd.copy());
		
			int N = Fluid.N/8;
			int boxsize = 2;

			for(int i=0; i<Fluid.M; i++)
			{
				int starti = i-boxsize >= 0 ? i-boxsize : 0;
				int stopi = i+boxsize < Fluid.M ? i+boxsize : Fluid.M-1;

				Fluid.p = 0;
				Fluid.p(Range(starti, stopi), Range(N-boxsize, N+boxsize)) = 100e5;

				Fluid2Slipper();
				slipper.calcDeform(0, 0, 0);
				Slipper2Fluid(1);
				Fluid.newt_ehdsqzapprox(i) = mean(slipper.ehd(Range(starti, stopi), Range(N-boxsize, N+boxsize)))/100e5;
			}

			Fluid.p = pbackup;
			slipper.ehd = ehdbackup;
			newt_ehdsqz_built = true;
		}

		if(sum(slipper.ehd) == 0)
		{
			//Initialize EHD
			SlipperInitEHD();

			//Init oldehd
			//slipper.oldehd = slipper.ehd;
			//Setting oldehd = 0 disables the ehdsqz at the first timestep
			slipper.oldehd = 0;
		}  else {
			//Update oldehd
			slipper.oldehd = slipper.ehd;
			
			//Initialize EHD
			//SlipperInitEHD();
		}

	}

	if(gapinput->options_slipper.general.SwashplatePressureDeformation == 1)
	{
		if(sum(swashplate.ehd) == 0)
		{
			//Initialize EHD
			//SwashplateCalcEHD(1.0);
			SwashplateInitEHD(1.0);

			//Init oldehd
			//swashplate.oldehd = swashplate.ehd;
			//Setting oldehd = 0 disables the ehdsqz at the first timestep
			swashplate.oldehd = 0;
		}  else {

			//old = old
			Array<double, 2> p_swashehd_old(swashplate.ehd.copy());

			//old defomation @ old location, used in part by the modified ehdsqz
			//swashplate.oldehd = swashplate.ehd;

			//Update the EHD field using swashplate deformation at the previous timestep
			//but at the new GCS slipper position. Note that SlipperSetCoordSys() needs to be called first.
			Swashplate2Fluid(1.0);	//WHEN THIS IS DISABLED, SWASHPLATE SQZ DOES NOT MAKE SENSE

			//modified ehdsqz which considers an avg of the old & new loc
			//swashplate.oldehd += swashplate.ehd;
			//swashplate.oldehd *= 0.5;

			//the traditional ehdsqz 
			swashplate.oldehd = swashplate.ehd;

			//old = movedold
			//Array<double, 2> p_swashehd_old(swashplate.ehd.copy());
	
			//SwashplateInitEHD needs to be called because the normal calc function does NOT solve the normal IM problem.
			//The calc method only updates the 'active' slipper deformation, holding the other n-1 pressures constant
			//SwashplateInitEHD instead solves the full swashplate IM problem
			SwashplateInitEHD(0.01);
			//SwashplateInitEHD(1.0);
			
			/*
			Array<double, 2> p_swashehd_new(swashplate.ehd.copy());

			//predict a change in the new velocity based on the change of swashplate deformation
			{
				//compute the regression plane's for p_swashehd_old and p_swashehd_new.
				
				//The A matrix is the same for both
				Array<double, 2> A(3,3);
				A(0,0) = sum(Fluid.Lx*Fluid.Lx);
				A(0,1) = sum(Fluid.Lx*Fluid.Ly);
				A(0,2) = sum(Fluid.Lx);
				A(1,0) = A(0,1);
				A(1,1) = sum(Fluid.Ly*Fluid.Ly);
				A(1,2) = sum(Fluid.Ly);
				A(2,0) = A(0,2);
				A(2,1) = A(1,2);
				A(2,2) = Fluid.M*Fluid.N;

				//create the two b vectors
				TinyVector<double, 3> b_old, b_new;
				
				b_old(0) = sum(Fluid.Lx*p_swashehd_old);
				b_old(1) = sum(Fluid.Ly*p_swashehd_old);
				b_old(2) = sum(p_swashehd_old);

				b_new(0) = sum(Fluid.Lx*p_swashehd_new);
				b_new(1) = sum(Fluid.Ly*p_swashehd_new);
				b_new(2) = sum(p_swashehd_new);

				//solve for the two actual planes
				TinyVector<double, 3> plane_old, plane_new;
				plane_old = Solve3(A, b_old);
				plane_new = Solve3(A, b_new);

				//get the Lx & Ly for the 3 control points
				TinyVector<double, 2> pt1, pt2, pt3;
				pt1(0) = geometryslippergap.routG*sin(0.0);
				pt1(1) = geometryslippergap.routG*cos(0.0);

				pt2(0) = geometryslippergap.routG*sin(2.0*PI/3.0);
				pt2(1) = geometryslippergap.routG*cos(2.0*PI/3.0);

				pt3(0) = geometryslippergap.routG*sin(4.0*PI/3.0);
				pt3(1) = geometryslippergap.routG*cos(4.0*PI/3.0);

				//get the z_old for each of the points
				double z1_old = plane_old(0)*pt1(0)+plane_old(1)*pt1(1)+plane_old(2);
				double z2_old = plane_old(0)*pt2(0)+plane_old(1)*pt2(1)+plane_old(2);
				double z3_old = plane_old(0)*pt3(0)+plane_old(1)*pt3(1)+plane_old(2);

				//get the z_new for each of the points
				double z1_new = plane_new(0)*pt1(0)+plane_new(1)*pt1(1)+plane_new(2);
				double z2_new = plane_new(0)*pt2(0)+plane_new(1)*pt2(1)+plane_new(2);
				double z3_new = plane_new(0)*pt3(0)+plane_new(1)*pt3(1)+plane_new(2);
			
				//find the dh/dt for each control point
				double dhdt1 = (z1_new-z1_old)/gapinput->common.timestep;
				double dhdt2 = (z2_new-z2_old)/gapinput->common.timestep;
				double dhdt3 = (z3_new-z3_old)/gapinput->common.timestep;
				
				
				//GapLog << "z_old = " << z1_old << "\t" << z2_old << "\t" << z3_old << endl;
				//GapLog << "z_new = " << z1_new << "\t" << z2_new << "\t" << z3_new << endl;
				GapLog << "Updating the initial velocity by " << dhdt1 << "\t" << dhdt2 << "\t" << dhdt3 << endl;
				Newtonvelocityslipper[0] += dhdt1;
				Newtonvelocityslipper[1] += dhdt2;
				Newtonvelocityslipper[2] += dhdt3;
				
			}
			
			//set the swash deformation to the damped value to ease it in the new FSI loop
			Swashplate2Fluid(0.01);

			*/
		}
	}

	//Slipper thermoelastic
	if(gapinput->options_slipper.general.SlipperThermoElastic == 1)
	{
		if(gapinput->common.first_rev_step)
		{	
			//Update the thermal problem / thermal deformation
			CalcSlipperThermal(gapinput->options_slipper.numeric.AlphaTEHD);
		}

		//Update the slipper thermal deformation / temperature to the fluid grid
		Slipper2Fluid_Thermal(gapinput->options_slipper.numeric.AlphaTEHD);
	}

	//Swashplate thermoelastic
	if(gapinput->options_slipper.general.SwashplateThermoElastic)
	{
		if(gapinput->common.first_rev_step)
		{	
			//Update the thermal problem / thermal deformation
			CalcSwashplateThermal(gapinput->options_slipper.numeric.AlphaTEHD);
		}

		//Update the slipper thermal deformation / temperature to the fluid grid
		Swashplate2Fluid_Thermal(gapinput->options_slipper.numeric.AlphaTEHD);
	}

	//Requires fluid boundaries to be initialized, requires Slipper/Swashplate2Fluid_Thermal to be called first
	SlipperInitializeTemperature();

	//Requires temperatue and pressure fields to be initialized
	SlipperCalcViscosity(1.0);

	//Update the old alpha value for the complex picard relaxation loop
	Fluid.alphaPoldCP = Fluid.alphaPold;

	//This will set the fluid film thickness to the specified initial gap height even considering all the of deformations
	//This is required in the case of large pressure deformations (common when considering the swashplate)
	static bool initGapHeight = true;	//This would fail if the executable is not restarted for a new simulation (unlikely, but a bug)
	if(initGapHeight)
	{
		if(!gapinput->common.resumed)	//Don't init the gap height if a resumed simulation
		{
			//Update the gap height
			SlipperCalch();

			//update the min gap height to input value
			const double dh = oldPosition[0] - min(Fluid.h);
			for(int i=0; i<3; i++)
			{
				oldPosition[i] += dh;
			}

			Newtonpositionslipper = control_point_position = pold = oldPosition;
		}
		initGapHeight = false;
	}
}
void CSlipperGap::SlipperEndTimestep(void)
{
	if(gapinput->options_slipper.general.SlipperThermoElastic == 1)
	{
		//Add the qflux from the fluid at this timestep to the slipper
		Fluid2Slipper_Thermal();
	}
	
	if(gapinput->options_slipper.general.SwashplateThermoElastic)
	{
		//Add the qflux from the fluid at this timestep to the swashplate
		Fluid2Swashplate_Thermal();
	}
}
void CSlipperGap::SlipperStartRevolution(void)
{
	//CURRENTLY NOT EVER CALLED BY ODE
}
void CSlipperGap::SlipperEndRevolution(void)
{
	//CURRENTLY NOT EVER CALLED BY ODE
}
void CSlipperGap::SlipperSetCoordSys(void)
{
	//Set up the LRCS (Local rotating coordinate system)
	{
		//ROTATION calculations based on speedk
		
		//An offset used to set the 'initial' rotation - USE WITH CAUTION. **DO NOT COMMIT TO GIT A NON-ZERO VALUE**
		double Roffset = 0.0;

		for(int i=0; i<Fluid.M; i++)
		{
			for(int j=0; j<Fluid.N; j++)
			{
				Fluid.LRx(i,j) = Fluid.r(i,j)*sin(Fluid.theta(i,j) + operatingslippergap.speedK*operatingslippergap.phi_rad + Roffset);
				Fluid.LRy(i,j) = Fluid.r(i,j)*cos(Fluid.theta(i,j) + operatingslippergap.speedK*operatingslippergap.phi_rad + Roffset);
			}
		}
	}

	//Set up the GCS (Global coordinate system)
	{
		const double phi = operatingslippergap.phi_rad;
		const double beta = operatingslippergap.beta_rad;
		const double gamma = operatingslippergap.gamma_rad;
		const double rB = geometryslippergap.rB;

		for(int i=0; i<Fluid.M; i++)
		{
			for(int j=0; j<Fluid.N; j++)
			{
				//Move LCS -> GCS
				//Rotate the LCS about its origin to the GCS alignment
				Fluid.Gx(i,j) = Fluid.Lx(i,j)*cos(-phi) - Fluid.Ly(i,j)*sin(-phi);
				Fluid.Gy(i,j) = Fluid.Lx(i,j)*sin(-phi) + Fluid.Ly(i,j)*cos(-phi);
				
				//Translate the LCS origin to the GCS origin
				//Fluid.Gx(i,j) += rB/cos(gamma)*sin(phi) - rB*tan(beta)*sin(gamma)*cos(phi);
				//Fluid.Gy(i,j) += rB/cos(beta)*cos(phi) - rB*tan(gamma)*sin(beta)*sin(phi);
				Fluid.Gx(i,j) += rB/cos(gamma)*sin(phi);
				Fluid.Gy(i,j) += rB/cos(beta)*cos(phi);
			}
		}
	}

	if(gapinput->options_slipper.general.SlipperPressureDeformation == 1 || gapinput->options_slipper.general.SlipperThermoElastic == 1)
	{
		//Because the LRCS has changed, we need to update the KDslip tree

		for(int i=0; i<Fluid.M; i++)
		{
			for(int j=0; j<Fluid.N; j++)
			{
				Fluid.KDslip.points[i*Fluid.N+j][0] = Fluid.LRx(i,j);
				Fluid.KDslip.points[i*Fluid.N+j][1] = Fluid.LRy(i,j);
			}
		}

		if(Fluid.KDslip.kdtree == NULL)
		{
			Fluid.KDslip.kdtree = new ANNkd_tree(Fluid.KDslip.points, Fluid.M*Fluid.N, Fluid.KDslip.dim);
		} else {
			delete Fluid.KDslip.kdtree;
			Fluid.KDslip.kdtree = new ANNkd_tree(Fluid.KDslip.points, Fluid.M*Fluid.N, Fluid.KDslip.dim);
		}
	}
}
void CSlipperGap::SlipperCalcFluidForces(void)
{
	//Calculate fluid pressure force and moments
	int M = Fluid.M; int N = Fluid.N; int Q = Fluid.Q;
	Array<double,2> Tr(M,N);
	Array<double,2> Tt(M,N);
	Array<double,2> Tx(M,N);
	Array<double,2> Ty(M,N);

	forcesslippergap.F = Fluid.p*Fluid.dA; //calcualte force in ring
	
	Array<double,2> mu(M,N);
	mu = Fluid.oilviscosity(Range::all(),Range::all(),Fluid.Q-2);

	//calculate viscous friction in polar
	Tr = mu*(Fluid.vr(all,all,Q-1)-Fluid.vr(all,all,Q-2))/(Fluid.h/(Q-1));
	Tt = mu*(Fluid.vtheta(all,all,Q-1)-Fluid.vtheta(all,all,Q-2))/(Fluid.h/(Q-1));

	//transform viscous friction to cartesian
	Tx = Tr*sin(Fluid.theta)+Tt*cos(Fluid.theta);
	Ty = Tr*cos(Fluid.theta)-Tt*sin(Fluid.theta);

	//Sum for friction force
	forcesslippergap.F_FTGx = -sum(Tx*Fluid.dA);
	forcesslippergap.F_FTGy = -sum(Ty*Fluid.dA);

	//Calculate moment due to friction
	forcesslippergap.M_FTGx = forcesslippergap.F_FTGy * geometryslippergap.lG;
	forcesslippergap.M_FTGy = -forcesslippergap.F_FTGx * geometryslippergap.lG;
	
	//Friction force magnitude	
	forcesslippergap.FTG = pow(forcesslippergap.F_FTGx*forcesslippergap.F_FTGx+forcesslippergap.F_FTGy*forcesslippergap.F_FTGy,.5);
	
	//Total fluid force in Z [N]
	forcesslippergap.FfGz = sum(forcesslippergap.F)+PI*(pow(geometryslippergap.dinG/2,2)-pow(geometryslippergap.dDG/2,2))*operatingslippergap.pG;

	//Moment due to unequal fluid pressure field
	forcesslippergap.MfGx = sum(forcesslippergap.F*Fluid.Ly);	//[Nm]
	forcesslippergap.MfGy = -sum(forcesslippergap.F*Fluid.Lx);	//[Nm]

	//Calculate power loss
	//Due to friction
	if(gapinput->options_slipper.general.CalcEnergy)
	{
		operatingslippergap.Ploss = Fluid.Qphisum;
		operatingslippergap.PlossMech = sum(Tr*Fluid.dA*operatingslippergap.vgr)+sum(Tt*Fluid.dA*operatingslippergap.vgtheta);
	} else {
		operatingslippergap.Ploss = sum(Tr*Fluid.dA*operatingslippergap.vgr)+sum(Tt*Fluid.dA*operatingslippergap.vgtheta);
		operatingslippergap.PlossMech = operatingslippergap.Ploss;
		//Due to leakage
		if(operatingslippergap.QSG > 0)
		{
			operatingslippergap.Ploss += operatingslippergap.QSG*(operatingslippergap.pG - operatingslippergap.pcase);
		}
	}


}
void CSlipperGap::SlipperCalcExternalForces(void) {
	forcesslippergap.MGx_centrifugal = geometryslippergap.mG*(geometryslippergap.rB)*operatingslippergap.omega*operatingslippergap.omega*geometryslippergap.lSG*cos(operatingslippergap.beta_rad);

}
void CSlipperGap::GetFSK(void) 
{
	//guess FSK function
	//For pressure
	const double dK = geometryslippergap.dK;
	const double dDK = geometryslippergap.dDK;
	//const double AreaK = PI / 4 * ( pow(dK,2) - pow(dDK,2) );
	const double AreaK = PI / 4 * ( pow(dK,2) - pow(geometryslippergap.dDG,2) );	

	const double omega = operatingslippergap.omega;
	const double rB = geometryslippergap.rB;
	const double beta = operatingslippergap.beta_rad;
	const double phi = operatingslippergap.phi_rad;
	const double mK = geometryslippergap.mK;

	//Pressure force
	forcesslippergap.FSK = AreaK *(operatingslippergap.pDC);
	//Inertial force
	forcesslippergap.FSK +=  mK * (omega*omega) * rB * tan(beta) * cos(phi);
	//Friction force
	forcesslippergap.FSK += forcesslippergap.FTK;
	
	//Reaction
	forcesslippergap.FSK /= cos(beta)*cos(operatingslippergap.gamma_rad);

	//Socket pressure
	//first find the 'projected' socket area	
	double FHD = 0; 
	if(geometryslippergap.hmaxG < -999*1e-6)
	{
		FHD = geometryslippergap.Fslipper/(double)geometryslippergap.npistons;		
	}
	const double AreaSock = gapinput->geometry.aSock;
	forcesslippergap.psocket = (forcesslippergap.FSK-FHD) / AreaSock;

	//Case pressure pushing on the top of the slipper
	//Calculated by ensuring force balance if piston+slipper were surrounded by constant pressure
	forcesslippergap.FSK += PI/4.0*operatingslippergap.pcase * (
												(pow(geometryslippergap.doutG,2) - pow(geometryslippergap.dDG,2)) - 
												(pow(dK,2) - pow(geometryslippergap.dDG,2))/(cos(operatingslippergap.beta_rad)*cos(operatingslippergap.gamma_rad))
												);
	
	GapResult->SlipperResult.FSK = forcesslippergap.FSK;
}
void CSlipperGap::GetFTK(void) 
{
	//Coupled caspar forces - fix start deg to 0
	if(gapinput->common.first_rev_step)	//used to correct floating point error near 360 deg
	{
		//we are at the beginning of a new revolution
		//see if we can open the ftk file
		ifstream ftk("./output/piston/ftk.txt");
		if(ftk.is_open())
		{
			//the file exists, let's grab the data
			forcesslippergap.FTK_vector.clear();

			while(!ftk.eof())
			{
				string l;
				getline(ftk,l);
				vector<string> line;
				Tokenize(l, line, "\t");
				if(line.size() >= 3)
				{
					vector<double> tmp(2);
					tmp[0] = atof(line[0].c_str());
					tmp[1] = atof(line[2].c_str());
					forcesslippergap.FTK_vector.push_back(tmp);
				}
			}

		} else {
			//file doesn't exist
			forcesslippergap.FTK_vector.clear();
		}
	}

	//coupled caspar, try and set FTK
	if(forcesslippergap.FTK_vector.size() > 0)
	{
		//time for one rev
		double revtime = gapinput->common.rev_period;
		double time = gapinput->common.time;
		double searchtime = time - revtime;

		int a = (int) forcesslippergap.FTK_vector.size()-1;
		while(forcesslippergap.FTK_vector[a][0] > searchtime && a > 0)
		{
			a--;
		}

		if((a == 0 && forcesslippergap.FTK_vector[a][0] >= searchtime)  || a == forcesslippergap.FTK_vector.size()-1)
		{
			//check to see why a landed at the bound
			if(fabs(forcesslippergap.FTK_vector[a][0] - searchtime) < revtime/3600.0)
			{
				//within 0.1 degree so just use it
				forcesslippergap.FTK = forcesslippergap.FTK_vector[a][1];
			} else {
				//more than 0.1 degree away so disregard
				forcesslippergap.FTK = 0.0;
			}
		} else {
			//interpolate
			int b = a + 1;
			double m = (forcesslippergap.FTK_vector[b][1]-forcesslippergap.FTK_vector[a][1])/(forcesslippergap.FTK_vector[b][0]-forcesslippergap.FTK_vector[a][0]);
			forcesslippergap.FTK = forcesslippergap.FTK_vector[a][1] + m*(searchtime-forcesslippergap.FTK_vector[a][0]);
		}

	} else {
		//no FTK data exists
		forcesslippergap.FTK = 0.0;

		//Added by Jeremy on 10.9.2015---------------------------------------------//
		//If no FTK data exists, then slipper simulation is not coupled to a piston simulation
		//In this case, we calculate piston axial friction assuming "parallel gap" conditions. That is, the piston axis coincident with the bore axis. 
		//For derivation, see [fpworks:\Jeremy Beale\03Research\Thesis Research\Piston Couette and Poiseuille Friction\Piston Couette and Poiseuille Friction Cart Coords.docx]
		const double dZ = geometryslippergap.dZ;
		const double dK = geometryslippergap.dK;
		const double dDK = geometryslippergap.dDK;
		const double dDG = geometryslippergap.dDG;
		const double AreaK = PI / 4 * ( pow(dK,2) - pow(dDG,2) );	
		const double lKG = geometryslippergap.lKG;

		const double T_Leak = operatingslippergap.T_Leak;
		const double pDC = operatingslippergap.pDC;
		const double pcase = operatingslippergap.pcase;
		const double omega = operatingslippergap.omega;
		const double rB = geometryslippergap.rB;
		const double beta = operatingslippergap.beta_rad;
		const double phi = operatingslippergap.phi_rad;
		const double mK = geometryslippergap.mK;
		
		double U_Ky = omega*rB*tan(beta)*sin(phi);
		double gap_height = dZ-dK;
		double dpdy = (pDC-pcase)/lKG;
		double mu = oil_properties->get_mu(pDC/2.0, T_Leak);

		double F_zy = (3.1415*dK*lKG)*(-gap_height/2.0*dpdy + mu*U_Ky/gap_height);
		forcesslippergap.FTK = F_zy;
		//-------------------------------------------------------------------------//
	}
}
void CSlipperGap::SlipperCalcHolder(const vector<double> &xg, vector<double> &F, const Array<double,2>& A) 
{
	/*
	//Method used in SHS optimization
	if(geometryslippergap.hmaxG < -999*1e-6)
	{
		//spring hold down

		if(abs(xg[0] - xg[1]) < 0.5e-6 && abs(xg[1] - xg[2]) < 0.5e-6 && abs(xg[0]-xg[2]) < 0.5e-6)
		{
			//The slipper is flat so apply the holder uniformly
			for(int i=0;i<3;i++)
			{
				F[i] = geometryslippergap.Fslipper/3.0/(double)geometryslippergap.npistons;
			}
		} else {
			//Apply the spring force at the "holder" contact point

			TinyVector<double,3> b;
			TinyVector<double,3> x;
			//get the maximum gap location
			TinyVector<int,1> loc = maxIndex(Fluid.hrigid(Fluid.M-1,all));
			double HoldForce = geometryslippergap.Fslipper/(double)geometryslippergap.npistons;
			double dx = Fluid.Rx(Fluid.M-1,loc(0));
			double dy = Fluid.Ry(Fluid.M-1,loc(0));
			
			//b = Fz,Mx,My
			b(0) = HoldForce;		//Z
			b(1) = HoldForce*dy;	//Mx
			b(2) = -HoldForce*dx;	//My
			
			x = Solve3(A,b);
			for(int i=0;i<3;i++)
			{
				F[i] = x(i);
			}
		}
	} else {
		//fixed hold down	
		for(int i=0;i<3;i++)
		{
			if(xg[i] > geometryslippergap.hmaxG)
			{
				//Simple F=kx reaction force
				F[i] = (xg[i]-geometryslippergap.hmaxG)*geometryslippergap.Fslipper;
			} else {
				F[i] = 0.0;
			}	
		}
	}
	*/

	//New method
	if(geometryslippergap.hmaxG < -999*1e-6)
	{
		//spring hold down
		Array<double,2> hT (Fluid.h + swashplate.ehd - Fluid.hgroove*5.0e-6);
		double maxdist = max(hT(Fluid.M-1,all))-min(hT(Fluid.M-1,all));
		
		// maxdist < 5 um will apply force equally
		// 10 um < maxdist will apply force at max point
		// maxdist in between will use a linear interpolation
		double upper = 10e-6;
		double lower = 5e-6; //5: value used for Japan, 2 is "new" testing val //double lower = 2e-6;

		double Fu, Fp;
		double Ftot = geometryslippergap.Fslipper/(double)geometryslippergap.npistons;
		
		//H1 Wave Spring
		//Ftot += (max(Fluid.hrigid(Fluid.M-1,all)) * 0.5733) / (double)geometryslippergap.npistons;

		/*
		//TESTING :: APPLYING ALL THE FORCE UNIFORMLY
		Fu = Ftot;
		Fp = 0.0;
		*/
		
		if (maxdist < lower)
		{
			Fu = Ftot;
			Fp = 0.0;
		} else if ( maxdist > upper)
		{
			Fp = Ftot;
			Fu = 0.0;
		} else {
			//use linear interpolation
			Fp = (maxdist - lower)/(upper -  lower)*Ftot;
			Fu = Ftot - Fp;
		}
		
		
		//Apply the uniform portion
		for(int i=0;i<3;i++)
		{
			F[i] = Fu/3.0;
		}
		
		//Apply the portion of force at the "holder" contact point
		TinyVector<double,3> b;
		TinyVector<double,3> x;
		//get the maximum gap location
		TinyVector<int,1> loc = maxIndex(hT(Fluid.M-1,all));
		double HoldForce = Fp;
		double dx = Fluid.Lx(Fluid.M-1,loc(0));
		double dy = Fluid.Ly(Fluid.M-1,loc(0));
	
		//b = Fz,Mx,My
		b(0) = HoldForce;		//Z
		b(1) = HoldForce*dy;	//Mx
		b(2) = -HoldForce*dx;	//My
		
		x = Solve3(A,b);
		for(int i=0;i<3;i++)
		{
			F[i] += x(i);
		}
		
	} else {
		/*
		//fixed hold down	
		for(int i=0;i<3;i++)
		{
			if(xg[i] > geometryslippergap.hmaxG)
			{
				//Simple F=kx reaction force
				F[i] = (xg[i]-geometryslippergap.hmaxG)*geometryslippergap.Fslipper;
			} else {
				F[i] = 0.0;
			}	
		}
		*/
		//Array<double,1> hT (Fluid.h(Fluid.M-1, Range::all()) + swashplate.ehd(Fluid.M-1, Range::all()));
		Array<double,1> hT (Fluid.h(Fluid.M-1, Range::all()));
		Array<double,1> contact(hT - geometryslippergap.hmaxG);
		contact = where(contact < 0, 0, contact);

		Array<double,1> contactp ( contact*geometryslippergap.Fslipper / (double)Fluid.N );

		TinyVector<double,3> b;
		TinyVector<double,3> x;

		double Fz = sum(contactp);
		double Mx = sum(contactp*Fluid.Ly(Fluid.M-1,Range::all()));
		double My = sum(-contactp*Fluid.Lx(Fluid.M-1,Range::all()));

		b = Fz,Mx,My;
		x = Solve3(A,b);

		for(int i=0; i<3; i++)
		{
			F[i] = x(i);
		}

	}



}
void CSlipperGap::SlipperCalcdF(vector<double> &dF, const vector<double> &xg) 
{
	double Mx,My,Fz,Rg;
	vector<double> HoldDownF(3);

	Rg = geometryslippergap.routG;
	Array<double,2> A(3,3);
	TinyVector<double,3> b;
	TinyVector<double,3> x;

	//Create geometry matrix
	A = 1,1,1,Rg,-.5*Rg,-.5*Rg,0,-sqrt(3.0)/2*Rg,sqrt(3.0)/2*Rg;
	
	//Fluid Forces
	Fz = forcesslippergap.FfGz;
	Mx = forcesslippergap.MfGx; 
	My = forcesslippergap.MfGy;
	b = Fz,Mx,My;
	x = Solve3(A,b);
	GapResult->SlipperResult.F_fluid[0] = x(0);
	GapResult->SlipperResult.F_fluid[1] = x(1);
	GapResult->SlipperResult.F_fluid[2] = x(2);

	//a positive socketR value means that for pumping mode at ODC the fsk resultant force is located
	//radially outward of the slipper center. a negative socketR value means that the resultant force
	//is located radially inward at pumping ODC.

	//double socketR = 0.11/1000.0;	//s90 75cc @ 100% beta?
	//double socketR = -0.0464/1000.0;	//h1 130cc @ 47% beta
	double socketR = 0.0;
	double socketRx = -forcesslippergap.FSK*socketR*cos(operatingslippergap.phi_rad);
	double socketRy = -forcesslippergap.FSK*socketR*sin(operatingslippergap.phi_rad);


	//External Forces
	Fz = -forcesslippergap.FSK;
	Mx = forcesslippergap.MGx_centrifugal + socketRx;
	//Mx = forcesslippergap.MGx_centrifugal;
	My = forcesslippergap.M_FTGy + socketRy;
	b = Fz,Mx,My;
	x = Solve3(A,b);
	GapResult->SlipperResult.F_external[0] = x(0);
	GapResult->SlipperResult.F_external[1] = x(1);
	GapResult->SlipperResult.F_external[2] = x(2);
	
	//Contact Forces
	SlipperCalcHolder(xg,HoldDownF,A);
	forcesslippergap.F_contact = HoldDownF;
	
	//Net Forces
	Fz = forcesslippergap.FfGz-forcesslippergap.FSK;

	//Mx = forcesslippergap.MfGx+forcesslippergap.MGx_centrifugal;
	Mx = forcesslippergap.MfGx+forcesslippergap.MGx_centrifugal+socketRx;

	//My = forcesslippergap.MfGy+forcesslippergap.M_FTGy;
	My = forcesslippergap.MfGy+forcesslippergap.M_FTGy+socketRy;
	//My = forcesslippergap.MfGy+forcesslippergap.M_FTGy-2.0;
	//My = forcesslippergap.MfGy+forcesslippergap.M_FTGy + forcesslippergap.FSK*(0.0/1000.0);
	

	//Added by Jeremy -------------------------------------------------------------------------
	//This code accounts for the effect of ball-joint friction on the slipper tilting moments.
	//In short, the ball-joint friction acts to reduce (or negate) the slipper tilting velocity.
	//For details on implementation, see [fpworks:\Jeremy Beale\01Presentations\GroupMeetings\Piston_Friction.pptx].
	//For details on simulation results, see [fpworks:\Jeremy Beale\01Presentations\GroupMeetings\Slipper_Sensitivity.pptx].

	//Initially assume socket is free to move.
	GapResult->SlipperResult.Sock_fixed = false;

	int ball_joint_friction = gapinput->options_slipper.friction.ball_joint_friction;
	if (ball_joint_friction != 0)
	{
		double r_socket = gapinput->geometry.r_socket;
		double doutG = gapinput->geometry.doutG;
		double C_joint = gapinput->options_slipper.friction.C_joint;
		//double stribeck_scale = gapinput->options_slipper.friction.stribeck_scale;
		double F_fric;
		double M_fric;
		double M_ext_x;
		double M_ext_y;
		double M_ext_mag;
		double socket_fric_moment_x;
		double socket_fric_moment_y;
		double visc_sock = 0.0;  //Initialize socket oil viscosity 
		double mu_coef = 0.0;  //Initialize friction coefficient
		

		//Calculate rotational speed [rps] about the socket using control point velocities
		double v1 = control_point_velocity[0];
		double v2 = control_point_velocity[1];
		double v3 = control_point_velocity[2];

		double v_avg = (v1+v2+v3)/3.0;

		//v*_rot represents the portion of velocity that contributes to rotation
		double v1_rot = v1 - v_avg;
		double v2_rot = v2 - v_avg;
		double v3_rot = v3 - v_avg;

		//Can determine rotational velocities directly from v1_rot and v2_rot
		double omega_x = v1_rot / doutG;

		//Solve the following equation for omega_y
		//v2_rot = omega_x * doutG * cos(60) + omega_y * doutG * cos(30)
		double omega_y = (v2_rot - omega_x * doutG * cos(60.0*3.1415/180.0)) / (doutG * cos(30.0*3.1415/180.0));

		//Overall rotational speed is the RMS magnitude of components
		double omega_net = pow(pow(omega_x,2.0) + pow(omega_y,2.0), 0.5);
		
		//Convert from rad/s to rev/s
		double omega_rps = omega_net / (2.0 * 3.1415);



		if (ball_joint_friction == 1)
		{
			bool tilting;
			double mu_coef_static = gapinput->options_slipper.friction.mu_coef_static;
			double mu_coef_dyn = gapinput->options_slipper.friction.mu_coef_dyn;

			//Determine if slipper was tilting during the previous time-step
			if (control_point_velocity[0] == control_point_velocity[1] && control_point_velocity[1] == control_point_velocity[2])
			{
				tilting = false;
			}
			else
			{
				tilting = true;
			}

			//Calculate the maximum possible friction moment in the ball joint.
			//This value can later be reduced to match the magnitude of the external moment.
			if (tilting == true)  //Use dynamic friction coefficient if slipper has tilt velocity
			{
				mu_coef = mu_coef_dyn;
			}
			else  //Otherwise, use static friction coefficient
			{
				mu_coef = mu_coef_static;
			}
		} 
		else if (ball_joint_friction == 2)
		{
			double mu_coef_pressure = gapinput->options_slipper.friction.mu_coef_pressure;
			double mu_coef_speed = gapinput->options_slipper.friction.mu_coef_speed;
			
			//Convert socket pressure to bar
			double psocket_bar = pow(10.0,-5.0) * forcesslippergap.psocket;

			mu_coef = fabs(psocket_bar * mu_coef_pressure) + fabs(omega_rps * mu_coef_speed);
			/*
			cout << endl << "CSlipperGap.cpp";
			cout << endl << "mu_coef_speed";
			cout << endl << gapinput->options_slipper.friction.mu_coef_speed;
			cout << endl << "mu_coef_dyn";
			cout << endl << gapinput->options_slipper.friction.mu_coef_dyn;
			cout << endl << "psocket_bar";
			cout << endl << psocket_bar;
			cout << endl << "mu_coef_pressure";
			cout << endl << mu_coef_pressure;
			cout << endl << "omega_rps";
			cout << endl << omega_rps;
			cout << endl << "mu_coef_speed";
			cout << endl << mu_coef_speed;
			cin.get();
			*/
		} 
		else if (ball_joint_friction == 3)
		{
			//Multiply by 2/pi, so that the range for mu_coeff is between [0,1]
			mu_coef = 2.0 / 3.1415 * atan(omega_rps);
		} 
		else if (ball_joint_friction == 4)
		{
			//Pressure in the socket was calculated previously in CSlipperGap::GetFSK
			double psocket = forcesslippergap.psocket;

			//Estimate viscosity of socket oil at leakage temperature and p_socket
			const double T_Leak = operatingslippergap.T_Leak;
			visc_sock = oil_properties->get_mu(psocket, T_Leak);
			
			//Calculate duty parameter (Hersey number)
			//double duty = mu * speed / psocket;
			//The previous line looks like the correct equation in literature, but doesn't work for this application, because duty would have values around 10^-9 or lower
			double duty = pow(10.0,12.0) * visc_sock * omega_rps / psocket;
			
			//Duty is going up to 1000 right now. Experimental data doesn't really exist beyond duty = 100. So let's reduce it.
			//Square root is a little too aggressive, so I'll use duty^0.7
			duty = pow(duty,0.7);

			//Calculate friction coefficient from Stribeck curve
			//See [fpworks:\Jeremy Beale\01Presentations\GroupMeetings\Slipper_Sensor_Questions.pptx] for visual comparison of these equations with Stribeck curve.
			//Equations were manually curve-fit using Matlab script in [fpworks:\Jeremy Beale\03Research\Slipper Code\stribeck_curve]
			double asperity = 0.0005 * pow(cos(0.17 * duty) + 0.7, 10.0);
			double fluid_film = 0.0004 * pow(duty, 1.1);

			if (duty < 10.0) 
			{
				mu_coef = asperity + fluid_film;
			} 
			else  //As duty parameter increases, socket motion becomes hydrodynamic
			{
				mu_coef = fluid_film;
			}
		}

		//On first time-step, assume that no friction is present. This may prevent the socket from freezing initially over the HP stroke.
		if(gapinput->common.first_step)
		{
			mu_coef = 0.0;
		}

		F_fric = fabs(forcesslippergap.FSK * mu_coef);

		//At this step of the calculation, 'F_fric' and 'M_fric' are scalar quantities, because their orientation is unknown
		M_fric = r_socket * F_fric;


		//Calculate external moments. 
		//That is, moments due to pressure (MfG), friction (M_FTGy), and inertia (MGx_centr) acting on the slipper.
		M_ext_x = forcesslippergap.MfGx + forcesslippergap.MGx_centrifugal;
		M_ext_y = forcesslippergap.MfGy + forcesslippergap.M_FTGy;
		M_ext_mag = pow(pow(M_ext_x,2.0)+pow(M_ext_y,2.0),0.5);  //Magnitude of external friction


		//Combine external moments with friction moments
		//The friction moment is artificially reduced by the coefficient 'C_joint', which is an empirical value 
		//Reduction by 'C_joint' is sometimes necessary because the slipper tilt angle might otherwise be "frozen" during the high pressure stroke.
		if (M_ext_mag > C_joint * M_fric) //Friction moment cannot exceed the external moment
		{
			Mx = (1.0 - C_joint * M_fric / M_ext_mag) * M_ext_x;
			My = (1.0 - C_joint * M_fric / M_ext_mag) * M_ext_y;

			//The Sock_fixed parameter is not output in slipper.txt, but it is used in CNewtonIteration::NewtonCalcShiftingVelocities. 
			GapResult->SlipperResult.Sock_fixed = false;
		}
		else  //If the friction moment appears to exceed the external moment, then the moments cancel and the net moment is 0. 
		{
			
			Mx = 0.0;
			My = 0.0;
			
			//The Sock_fixed parameter is not output in slipper.txt, but it is used in CNewtonIteration::NewtonCalcShiftingVelocities. 
			GapResult->SlipperResult.Sock_fixed = true;

			/*
			// Unfortunately, simulations are exibiting film thickness instability with all friction options.
			// So let's change net moments more gradually, by using Mx,My from the previous time-step.
			double blend_factor = 0.1;
			double Mx_old = GapResult->SlipperResult.dFcomp[1];
			double My_old = GapResult->SlipperResult.dFcomp[2];

			Mx = blend_factor * Mx + (1 - blend_factor) * Mx_old;
			My = blend_factor * My + (1 - blend_factor) * My_old;
			*/
		}


		//Also calculate the friction moment's x and y components for slipper.txt
		// friction moment = external moment - net moment
		socket_fric_moment_x = M_ext_x - Mx;
		socket_fric_moment_y = M_ext_y - My;

		
		//Assign net moments to result class
		GapResult->SlipperResult.M_TJx = socket_fric_moment_x;
		GapResult->SlipperResult.M_TJy = socket_fric_moment_y;

		//Assign other quantities to results class
		GapResult->SlipperResult.visc_sock = visc_sock;
		GapResult->SlipperResult.tilt_speed = omega_rps;
		GapResult->SlipperResult.psocket = forcesslippergap.psocket;
		GapResult->SlipperResult.mu_coef = mu_coef;
		GapResult->SlipperResult.F_fric = F_fric;
	}
	//-------------------------------------------------------------------------------------------


	//Determine the primary rotation magnitude
	double rSocket = 7.8/1000;	//xxxxxx ONLY GOOD FOR THE S90-75CC xxxxxxx
	double mu_cof_static = 0.05;
	double mu_cof_dyn = 0.045;
	double socket_fric_moment = 0;

	//we're now always using a full 3dof rbd solver
	//GapResult->SlipperResult.Sock_fixed = false;
	/*
	if(forcesslippergap.FSK > 0)	//friction when the socket is in tension is probably different
	{
		//this is the rotation vector
		//this is the velocity at the control points. To convert to Rx & Ry we need the A matrix
		double Rx = (2.0*control_point_velocity[0]-control_point_velocity[1]-control_point_velocity[2])/Rg;
		double Ry = (4*control_point_velocity[0]-3.0*control_point_velocity[1]-control_point_velocity[2])/(Rg*pow(3.0, 0.5));
		double theta = 0;
		if(Rx != 0)
		{
			theta = atan(Ry/Rx);
		} else {
			if(Ry > 0)
			{
				theta = PI/2.0;
			} else if (Ry < 0)
			{
				theta = -PI/2.0;
			}
		}
		double Rmag = pow(Rx*Rx+Ry*Ry, 0.5);
		if(Rx < 0)
		{
			Rmag = -Rmag;
		}

		//moment vector
		//double thetaM = atan(Mx/My);
		//double Msum = Mx/sin(thetaM);
		double Msum = pow(Mx*Mx+My*My, 0.5);
		GapResult->SlipperResult.Msock_sum = Msum;

		//if(fabs(Msum) > fabs(max_fric_moment))
		//if(fabs(Msum) > fabs(forcesslippergap.FSK*mu_cof_static*rSocket))
		if(fabs(Rmag) > 0.025 || fabs(Msum) > fabs(forcesslippergap.FSK*mu_cof_static*rSocket))
		{
			//dynamic friction
			double max_fric_moment = forcesslippergap.FSK*mu_cof_dyn*rSocket;
			GapResult->SlipperResult.Msock_maxfric = max_fric_moment;

			//this is only used to help the solver near 0 velocity
			double Rtol = 0.000648;
			double fric_mom = 0;
			if(fabs(Rmag) < Rtol)
			{
				fric_mom = max_fric_moment*fabs(Rmag)/Rtol;
			} else {
				fric_mom = max_fric_moment;
			}

			if(Rmag > 0)
			{
				Mx -= fric_mom*cos(theta);
				My -= fric_mom*sin(theta);
			} 
			else if(Rmag <= 0)
			{
				Mx -= -fric_mom*cos(theta);
				My -= -fric_mom*sin(theta);
			}
		} else {
			if(!gapinput->options_slipper.general.HybridForceBalance)
			{
				GapResult->SlipperResult.Sock_fixed = true;
				double max_fric_moment = forcesslippergap.FSK*mu_cof_static*rSocket;
				GapResult->SlipperResult.Msock_maxfric = max_fric_moment;
			} else {			
				//the slipper is being held by friction. introduce moment penilties for rotational motion
				double max_fric_moment = forcesslippergap.FSK*mu_cof_static*rSocket;
				GapResult->SlipperResult.Msock_maxfric = max_fric_moment;

				//this is only used to help the solver near 0 velocity
				double Rtol = 0.029651;
				double fric_mom = 0;
				if(fabs(Rmag) < Rtol)
				{
					fric_mom = max_fric_moment*fabs(Rmag)/Rtol;
				} else {
					fric_mom = max_fric_moment;
				}

				if(Rmag > 0)
				{
					Mx = fric_mom*cos(theta);
					My = fric_mom*sin(theta);
				} 
				else if(Rmag <= 0)
				{
					Mx = -fric_mom*cos(theta);
					My = -fric_mom*sin(theta);
				}
			}
			
		}
	}
	*/
	
	b = Fz,Mx,My;
	x = Solve3(A,b);

	GapResult->SlipperResult.dFcomp[0] = Fz;
	GapResult->SlipperResult.dFcomp[1] = Mx;
	GapResult->SlipperResult.dFcomp[2] = My;

	dF[0] = x(0)-HoldDownF[0];
	dF[1] = x(1)-HoldDownF[1];
	dF[2] = x(2)-HoldDownF[2];

	forcesslippergap.dF = dF;

	//Assign force to results
	GapResult->SlipperResult.dF = forcesslippergap.dF;
}
TinyVector<double,3> CSlipperGap::Solve3(Array<double,2> A, TinyVector<double,3> b)
{
	//Solves a 3x3 of the from {A}[x]=[b]
	
	TinyVector<double,3> x;
	Array<double,2> Ai(3,3); //A^-1
	double dA; //det(A)

	dA = A(0,0)*det(A,0,0)+A(0,1)*det(A,0,1)+A(0,2)*det(A,0,2);
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			Ai(j,i) = det(A,i,j);
		}
	}
	Ai = (1/dA)*Ai;

	x(0) = Ai(0,0)*b(0)+Ai(0,1)*b(1)+Ai(0,2)*b(2);
	x(1) = Ai(1,0)*b(0)+Ai(1,1)*b(1)+Ai(1,2)*b(2);
	x(2) = Ai(2,0)*b(0)+Ai(2,1)*b(1)+Ai(2,2)*b(2);

	return x;
}
double CSlipperGap::det(Array<double,2> m,int i, int j)
{
	//find the 2x2 det for pos i,j in the 3x3 matrix m
	int r[2],c[2];
	int sign=1;

	switch(i)
	{
	case 0:
		r[0] = 1;r[1] = 2;break;
	case 1:
		r[0] = 0;r[1] = 2;break;
	case 2:
		r[0] = 0;r[1] = 1;
	}
	switch(j)
	{
	case 0:
		c[0] = 1;c[1] = 2;break;
	case 1:
		c[0] = 0;c[1] = 2;break;
	case 2:
		c[0] = 0;c[1] = 1;
	}
	if ((i+j)%2==1)
	{
		sign = -1;
	}

	return sign*(m(r[0],c[0])*m(r[1],c[1])-m(r[0],c[1])*m(r[1],c[0]));
}
void CSlipperGap::SlipperAssignResults(const vector<double> &xg)
{
	GapResult->SlipperResult.P_slipper = Fluid.p;
	GapResult->SlipperResult.T_slipper = Fluid.T;
	GapResult->SlipperResult.h_slipper = Fluid.h;
	GapResult->SlipperResult.slipperEHD = slipper.ehd;
	GapResult->SlipperResult.swashplateEHD = swashplate.ehd;
	GapResult->SlipperResult.vr_slipper = Fluid.vr;
	GapResult->SlipperResult.vtheta_slipper = Fluid.vtheta;
	GapResult->SlipperResult.pG = Fluid.p(0,0);
	GapResult->SlipperResult.Ploss = operatingslippergap.Ploss;
	GapResult->SlipperResult.PlossMech = operatingslippergap.PlossMech;
	GapResult->SlipperResult.vgx = operatingslippergap.vgx;
	GapResult->SlipperResult.vgy = operatingslippergap.vgy;
	GapResult->SlipperResult.F_contact = forcesslippergap.F_contact;

	//----------Fluid Forces------------//
	GapResult->SlipperResult.FfGz = forcesslippergap.FfGz;
	GapResult->SlipperResult.MfGx = forcesslippergap.MfGx;
	GapResult->SlipperResult.MfGy = forcesslippergap.MfGy;
	
	//----------Contact Forces------------//
	GapResult->SlipperResult.avg_p_contact = sum(Fluid.p_contact * Fluid.dA) / sum(Fluid.dA);

	//----------Friction Forces------------//
	GapResult->SlipperResult.M_FTGx = forcesslippergap.M_FTGx;
	GapResult->SlipperResult.M_FTGy = forcesslippergap.M_FTGy;
	GapResult->SlipperResult.F_FTGx = forcesslippergap.F_FTGx;
	GapResult->SlipperResult.F_FTGy = forcesslippergap.F_FTGy;
	GapResult->SlipperResult.FTG = forcesslippergap.FTG;
	GapResult->SlipperResult.FTK = forcesslippergap.FTK;
	GapResult->SlipperResult.TorqueLoss = -forcesslippergap.F_FTGx*geometryslippergap.rB;
	
	//----------External Forces------------//
	GapResult->SlipperResult.MGx_centrifugal = forcesslippergap.MGx_centrifugal;
	//GapResult->SlipperResult.QSG = operatingslippergap.QSG_orifice;
	GapResult->SlipperResult.QSG = operatingslippergap.QSG;
	GapResult->SlipperResult.Q_S_pois = operatingslippergap.Q_S_pois;
	GapResult->SlipperResult.Q_S_couette = operatingslippergap.Q_S_couette;
	GapResult->SlipperResult.xg0 = xg[0];
	GapResult->SlipperResult.xg1 = xg[1];
	GapResult->SlipperResult.xg2 = xg[2];

	//----------Gap Heights------------//
	GapResult->SlipperResult.maxh = 0;
	GapResult->SlipperResult.minh = 1;
	GapResult->SlipperResult.meanh = 0;
	for(int i=0; i<Fluid.M; i++)
	{
		for(int j=0; j<Fluid.N; j++)
		{
			const double h = Fluid.h(i,j);
			if(h > GapResult->SlipperResult.maxh)
			{
				GapResult->SlipperResult.maxh = h;
			}

			if(h < GapResult->SlipperResult.minh)
			{
				GapResult->SlipperResult.minh = h;
			}

			GapResult->SlipperResult.meanh += h;
		}
	}
	GapResult->SlipperResult.meanh /= (double) Fluid.M*Fluid.N;
	
}
void CSlipperGap::Tokenize(const string& str,vector<string>& tokens,const string& delimiters)
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

