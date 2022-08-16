#include "CGapODE.h"
#include "Coupled.h"
#include "CPistonGap.h"
#include "logger.h"
#include <iomanip>
#include "../../caspar_input/input.h"
#include <iostream>
#pragma once

extern Concurrency::event * coupled_caspar_local_event;
extern Concurrency::event * coupled_caspar_global_event;
extern bool coupled_caspar_simulation;

extern class CGapInput myGapInput;
extern struct sGapResult myGapResult;

//extern class CPistonGap myPistonGap;//Lizhi
extern class CGapInput myGapInput;
extern struct sGapResult myGapResult;

extern class input myinput;

CGapODE::CGapODE(void)
{
	
	//Number of states for the solver (4 piston)
	n = 4;
	//Initial time
	x = 0.0;

	//ODE vector allocation
	y.resize(n+1);
	//ODE dydt vector allocation
	yp.resize(n);

	//Initialize velocity vector
	for(int i=0;i<n;i++)
	{
		yp[i] = 0.0;
		y[i+1] = 0.0;
	}

	//Initialize the parts positon based on the .bpd file
	y[0] = x;
	y[1] = myinput.data.options_piston.position.xA;
	y[2] = myinput.data.options_piston.position.yA;
	y[3] = myinput.data.options_piston.position.xB;
	y[4] = myinput.data.options_piston.position.yB;
	
	//Initial values for time
	xbeg = 0;
	myGapResult.time = xbeg;

};
CGapODE::~CGapODE(void)
{

};
//Solver main starting function
void CGapODE::ODEmain(void)
{	
	//Setup all the solver variables
	ODESetupSolver();

	//ODE solver 
	ODEIntegrate();
};
//Solver ode core function
void CGapODE::odef(double x,vector<double> &y,vector<double> &yp)
{

	Log << " " << endl;
	Log << " " << endl;
	Log << dec << "REVOLUTION: " << myGapResult.revcounter << endl;
	Log << "TIME: " << x << " [s]" << endl;

	//Newton Iteration main call for shifting velocities calculation
	myNewtonIteration.NewtonCalcNewtonIteration(y,yp);
	
};
//Solver variables initialization
void CGapODE::ODESetupSolver(void)
{	
	//Maximum stepsize during simulation [s]: (Desired Angular Step [rad])/(Pump Speed [rad/s])
	hmax = myinput.data.options_piston.numeric.Simalphastep*PI/180.0 / (myinput.data.operating_conditions.speed);
	myNewtonIteration.dt = hmax;
	//Revolutions counter
	revcounter = 0;
	//Final value for x
	xend = myinput.data.lubrication_module.n_lubrication_revolutions / (myinput.data.operating_conditions.speed/(2*PI));
	//Interval of x for printing output
	dx = hmax;
	//Output for 2D and 3D matrixes with the frequancy of plot step step
	dx2D = myinput.data.options_piston.numeric.Simalphaplot*PI/180.0 / (myinput.data.operating_conditions.speed); 
	//Final setup
	xold = xbeg;
	x = xbeg;
	xd = xbeg;
	xd2D = xbeg;

};
//Solver Euler fixed step integration routine
void CGapODE::ODEIntegrate(void)
{

	do
	{
		
		//call function
		odef(x,y,yp);
		//increase time step
		ODEAdaptTimeStep();
		//output solution
		ODESolutionOutput();
		//calculate new position
		for(int i=0; i<4; i++)
		{
			y[i+1] += hmax*yp[i];
		};
		//assign new time
		y[0] = x;

		//Handle coupled caspar blocking events
		if(coupled_caspar_simulation)
		{
			static int revcounter_GUI = 0; //counts the number of times the simulation has been held at the end of a revolution.
			
			//Check for the end of a revolution
			
			//double deg = fmod(6.0*myGapInput.General.speed*x, 360.0);//should work if needed
			int rev = floor(myinput.data.operating_conditions.speed*x/(2*PI));//may need to add ~1/7200 to get 
			
			if(revcounter_GUI<rev)
			{
				//this is a new revolution
				coupled_caspar_local_event->set();			//signal that the slipper has finished a revolution
				coupled_caspar_global_event->wait();		//wait for the all interfaces to finish
				coupled_caspar_global_event->reset();		//reset the global event flag
				revcounter_GUI += 1;
				//Check for updated FTG file from slipper interface dwm
				myGapInput.readFTG();
				
			}

		}

	}while(x<xend);

};
void CGapODE::ODEAdaptTimeStep(void)
{

	//Update xold time step
	xold = x;

	//Phi angle [rad] between 0 and 2PI
	double omega = myinput.data.operating_conditions.speed ;
	double phi_rad = fmod(omega*x,2.0*PI) ;
	double phi_deg = phi_rad * 180.0/PI ;
	double tol = 1.0e-7;
	double gamma = myinput.data.geometry.gamma;//[rad]
	double beta = myinput.data.operating_conditions.beta;//[rad]

	//Refine angular step to 0.5*dphi [deg] when switching form pHP to pLP

	//correct for gamma angle
	phi_deg -=atan(tan(gamma)/sin(beta))*(180/PI);

	if( phi_deg>=(175.0-tol) && phi_deg<(185.0-tol) )
	{
		double dphi_min = 0.1;
		hmax = dphi_min * PI/180.0 / omega ;
		dx = hmax;
		myNewtonIteration.dt = hmax;
	}
	//Keep input time step
	else
	{
		hmax = myinput.data.options_piston.numeric.Simalphastep * PI/180.0 / omega ;
		dx = hmax;
		myNewtonIteration.dt = hmax;
	}

	//Increase time step
	x += hmax;

}
//Solver output results
int CGapODE::ODESolutionOutput(void)
{

	//Fluid grid points
	int N,M,Q;
	if(myinput.data.options_piston.general.ReynoldsMultiGrid){
		N = myinput.data.options_piston.fluid_grid.MG.MG_N[0];
		M = myinput.data.options_piston.fluid_grid.MG.MG_M[0];
		Q = myinput.data.options_piston.fluid_grid.MG.Q;
	}
	else{
		N = myinput.data.options_piston.fluid_grid.GS.N;
		M = myinput.data.options_piston.fluid_grid.GS.M;
		Q = myinput.data.options_piston.fluid_grid.GS.Q;
	}
	double tol = 1.0e-7;
	double revtime = myNewtonIteration.revtime ;


	//OUTPUT OF 1D RESULTS AT EACH SIMULATION STEP
	if ( xold==0.0 || x>(xd+tol) ) 
	{

		//Population of total loss vectors
		myGapResult.QSKTot.push_back(myNewtonIteration.myPistonGap.operatingpistongap.QSK);
		myGapResult.MSKTot.push_back(myNewtonIteration.myPistonGap.forcespistongap.MTK);
		myGapResult.MFT.push_back(myNewtonIteration.myPistonGap.forcespistongap.MTKtan);
		myGapResult.PhiD_mech.push_back(myNewtonIteration.myPistonGap.PhiD_mech);
		myGapResult.PhiD_vol.push_back(myNewtonIteration.myPistonGap.PhiD_vol);
		//Time array
		xi.push_back(xold);

		double QSKTot = 0.0;
		double MSKTot = 0.0;
		double MFT = 0.0;
		double PhiD_mech_tot = 0;
		double PhiD_vol_tot = 0;

		//Total losses after first revolution, .qst file and .Kdisst file
		double omega = myinput.data.operating_conditions.speed;
		double dt_z = 2.0*PI/myinput.data.operating_conditions.npistons/omega;
		int n_z = myinput.data.operating_conditions.npistons;
		if(x > revtime && myinput.data.options_piston.general.EHDTestRig==0)
		{
			
			//Sum (n_z-1) previous pistons contributions
			for(int i=0; i<n_z-1; i++)
			{
				double t_z = xold - (i+1)*dt_z;
				for(int j=0;j<(int) xi.size();j++)
				{
					if(xi[j]>=t_z)
					{
						QSKTot += myGapResult.QSKTot[j];
						MSKTot += myGapResult.MSKTot[j];
						MFT += myGapResult.MFT[j];
						PhiD_mech_tot += myGapResult.PhiD_mech[j];
						PhiD_vol_tot += myGapResult.PhiD_vol[j];
						break;
					}
				}
			}
			//Sum actual piston contribution
			QSKTot += myNewtonIteration.myPistonGap.operatingpistongap.QSK;
			MSKTot += myNewtonIteration.myPistonGap.forcespistongap.MTK;
			MFT += myNewtonIteration.myPistonGap.forcespistongap.MTKtan;
			PhiD_mech_tot += myNewtonIteration.myPistonGap.PhiD_mech;
			PhiD_vol_tot += myNewtonIteration.myPistonGap.PhiD_vol;
		}


		//Positions and velocities, .pos file
		fout.open("./output/piston/piston.txt",ios::app);
		fout << scientific << myGapResult.time												//1
			<< "\t" << myinput.data.operating_conditions.speed * myGapResult.time / (2*PI)	//2
			<< "\t" << myGapResult.phi														//3
			<< "\t" << myGapResult.gappositions[1]											//4
			<< "\t" << myGapResult.gappositions[2]											//5
			<< "\t" << myGapResult.gappositions[3]											//6
			<< "\t" << myGapResult.gappositions[4]											//7
			<< "\t" << myGapResult.gapvelocities[0]											//8
			<< "\t" << myGapResult.gapvelocities[1]											//9
			<< "\t" << myGapResult.gapvelocities[2]											//10
			<< "\t" << myGapResult.gapvelocities[3]											//11
			<< "\t" << myNewtonIteration.myPistonGap.geometrypistongap.lvar					//12
			<< "\t" << myNewtonIteration.myPistonGap.geometrypistongap.lout					//13
			<< "\t" << myNewtonIteration.myPistonGap.geometrypistongap.pos_A				//14
			<< "\t" << myNewtonIteration.myPistonGap.geometrypistongap.pos_B				//15
			<< "\t" << myNewtonIteration.myPistonGap.geometrypistongap.zRK					//16
			<< "\t" << myNewtonIteration.myPistonGap.geometrypistongap.sK					//17
			<< "\t" << myNewtonIteration.myPistonGap.geometrypistongap.lA					//18
			<< "\t" << myNewtonIteration.myPistonGap.geometrypistongap.lB					//19
			<< "\t" << myNewtonIteration.myPistonGap.geometrypistongap.z_A					//20
			<< "\t" << myNewtonIteration.myPistonGap.geometrypistongap.z_B					//21
			<< "\t" << myNewtonIteration.myPistonGap.operatingpistongap.vK					//22
			<< "\t" << myNewtonIteration.myPistonGap.operatingpistongap.QSK					//23
			<< "\t" << myNewtonIteration.myPistonGap.operatingpistongap.QSK_p				//24
			<< "\t" << myNewtonIteration.myPistonGap.operatingpistongap.QSK_c				//25
			<< "\t" << myNewtonIteration.myPistonGap.forcespistongap.MTK					//26
			<< "\t" << myNewtonIteration.myPistonGap.forcespistongap.MTKtan					//27
			<< "\t" << sum(myNewtonIteration.myPistonGap.forcespistongap.FTKy)				//28
			<< "\t" << sum(myNewtonIteration.myPistonGap.forcespistongap.FTKy_p)			//29
			<< "\t" << sum(myNewtonIteration.myPistonGap.forcespistongap.FTKy_c)			//30
			<< "\t" << sum(myNewtonIteration.myPistonGap.forcespistongap.FTKx	)			//31
			<< "\t" << sum(myNewtonIteration.myPistonGap.forcespistongap.FTKx_p)			//32
			<< "\t" << sum(myNewtonIteration.myPistonGap.forcespistongap.FTKx_c)			//33
			<< "\t" << myNewtonIteration.myPistonGap.forcespistongap.FTG					//34
			<< "\t" << myGapResult.pDC														//35
			<< "\t" << sum(myNewtonIteration.myPistonGap.forcespistongap.FTby	)			//36
			<< "\t" << sum(myNewtonIteration.myPistonGap.forcespistongap.FTby_p)			//37
			<< "\t" << sum(myNewtonIteration.myPistonGap.forcespistongap.FTby_c)			//38
			<< "\t" << sum(myNewtonIteration.myPistonGap.forcespistongap.FTbx	)			//39
			<< "\t" << sum(myNewtonIteration.myPistonGap.forcespistongap.FTbx_p)			//40
			<< "\t" << sum(myNewtonIteration.myPistonGap.forcespistongap.FTbx_c)			//41
			<< "\t" << myGapResult.pHP														//42
			<< "\t" << myGapResult.PLP														//43
			<< "\t" << myNewtonIteration.myPistonGap.forcespistongap.F_external[0]			//44
			<< "\t" << myNewtonIteration.myPistonGap.forcespistongap.F_external[1]			//45
			<< "\t" << myNewtonIteration.myPistonGap.forcespistongap.F_external[2]			//46
			<< "\t" << myNewtonIteration.myPistonGap.forcespistongap.F_external[3]			//47
			<< "\t" << myNewtonIteration.myPistonGap.forcespistongap.F_fluid[0]				//48
			<< "\t" << myNewtonIteration.myPistonGap.forcespistongap.F_fluid[1]				//49
			<< "\t" << myNewtonIteration.myPistonGap.forcespistongap.F_fluid[2]				//50
			<< "\t" << myNewtonIteration.myPistonGap.forcespistongap.F_fluid[3]				//51
			<< "\t" << myNewtonIteration.myPistonGap.forcespistongap.F_contact[0]			//52
			<< "\t" << myNewtonIteration.myPistonGap.forcespistongap.F_contact[1]			//53
			<< "\t" << myNewtonIteration.myPistonGap.forcespistongap.F_contact[2]			//54
			<< "\t" << myNewtonIteration.myPistonGap.forcespistongap.F_contact[3]			//55
			<< "\t" << myNewtonIteration.myPistonGap.forcespistongap.dF[0]					//56
			<< "\t" << myNewtonIteration.myPistonGap.forcespistongap.dF[1]					//57
			<< "\t" << myNewtonIteration.myPistonGap.forcespistongap.dF[2]					//58
			<< "\t" << myNewtonIteration.myPistonGap.forcespistongap.dF[3]					//59
			<< "\t" << QSKTot																//60
			<< "\t" << MSKTot																//61
			<< "\t" << MFT																	//62
			<< "\t" << PhiD_mech_tot														//63
			<< "\t" << PhiD_vol_tot															//64
			<< "\t" << myNewtonIteration.myPistonGap.PhiD_mech								//65
			<< "\t" << myNewtonIteration.myPistonGap.PhiD_vol								//66
			<< endl;
		fout.close();
		fout.clear();

		//Output FTK force for coupled simulation.
		
		if(coupled_caspar_simulation){
			fout.open("./output/piston/ftk.txt",ios::app);
			fout << scientific << myGapResult.time											//1
				<< "\t" << myGapResult.phi													//2
				<< "\t" << sum(myNewtonIteration.myPistonGap.forcespistongap.FTKy) << endl;		//3
			fout.close();
			fout.clear();
		};


		//Output bodies surface temepratures and deformations at each revlution beginning
		if( xold==0.0 || xold > (revcounter*revtime-tol) )
		{
			//Surface temperature output
			if(myinput.data.options_piston.general.HeatTransfer)
			{
				//Piston total heat flux, .QK file  [W/m2]
				fout.open("./output/piston/matlab/Piston_Body_Heat_Flux.txt",ios::app);
				fout << "%Revolution Number: " << revcounter << endl;
				for(int i=0;i<(int) myNewtonIteration.myPistonGap.qbi_piston[0].extent(0);i++)
				{
					fout << scientific << myNewtonIteration.myPistonGap.qbi_piston[0](i) << endl;
				}
				fout << "%" << endl;
				fout.close();
				fout.clear();

				//Cylinder total heat flux, .QB file  [W/m2]
				fout.open("./output/piston/matlab/Piston_Cylinder_Body_Heat_Flux.txt",ios::app);
				fout << "%Revolution Number: " << revcounter << endl;
				for(int i=0;i<(int) myNewtonIteration.myPistonGap.qbi_cylinder[0].extent(0);i++)
				{
					fout << scientific << myNewtonIteration.myPistonGap.qbi_cylinder[0](i) << endl;
				}
				fout << "%" << endl;
				fout.close();
				fout.clear();

				//Piston surface temperature, .TKs file [C]
				fout.open("./output/piston/matlab/Piston_Body_Surface_Temperature.txt",ios::app);
				fout << "%Revolution Number: " << revcounter << endl;
				for(int i=0;i<(int) myNewtonIteration.myPistonGap.TK_surf.extent(0);i++)
				{
					fout << scientific << myNewtonIteration.myPistonGap.TK_surf(i) << endl;
				}
				fout << "%" << endl;
				fout.close();
				fout.clear();

				//Bushing surface temperature, .TBs file [C]
				fout.open("./output/piston/matlab/Piston_Cylinder_Body_Surface_Temperature.txt",ios::app);
				fout << "%Revolution Number: " << revcounter << endl;
				for(int i=0;i<(int) myNewtonIteration.myPistonGap.TB_surf.extent(0);i++)
				{
					fout << scientific << myNewtonIteration.myPistonGap.TB_surf(i) << endl;
				}
				fout << "%" << endl;
				fout.close();
				fout.clear();
			}

			//Thermal deformation output
			if(myinput.data.options_piston.general.ThermalDeformation)
			{
				//Piston surface elatic deformation due to thermal stress, .defKthb file [microns]
				fout.open("./output/piston/matlab/Piston_Thermal_Deformation.txt",ios::app);
				fout << "%Revolution Number: " << revcounter << endl;
				for(int i=0;i<(int) myNewtonIteration.myPistonGap.defK_th.extent(0);i++)
				{
					fout << scientific << myNewtonIteration.myPistonGap.defK_th(i)*1.0e6 << endl;
				}
				fout << "%" << endl;
				fout.close();
				fout.clear();

				//Piston cylinder surface elatic deformation due to thermal stress, .defBthb file [microns]
				fout.open("./output/piston/matlab/Piston_Cylinder_Thermal_Deformation.txt",ios::app);
				fout << "%Revolution Number: " << revcounter << endl;
				for(int i=0;i<(int) myNewtonIteration.myPistonGap.defB_th.extent(0);i++)
				{
					fout << scientific << myNewtonIteration.myPistonGap.defB_th(i)*1.0e6 << endl;
				}
				fout << "%" << endl;
				fout.close();
				fout.clear();
			}

			

			//Revolution number increment
			revcounter++;
		};

		//Increment time
		xd += dx;

	};


	//OUTPUT OF 2D RESULTS AND VTK AT EACH PLOT STEP	
	if( xold==0.0 || x>(xd2D+tol) ) 
	{
		/*Test outputs by dan
		//Piston gap viscosity from convergence loop, .pK file [bar]
		fout.open("./outputs_piston/Piston_Gap_Viscosity.txt",ios::app);
		fout << "%PHI [deg]: " << myGapResult.phi << endl;
		for(int i = 0; i < N ; i++)
		{
			for(int j = 0; j < M ; j++)
            {
					fout << scientific << myNewtonIteration.myPistonGap.oilviscosity(j+i*M) << "\t" ;
            }
			fout << endl;
		}
		fout << "%" << endl;
		fout.close();
		fout.clear();
		
		//Piston gap density from convergence loop, .pK file [bar]
		fout.open("./outputs_piston/Piston_Gap_Density.txt",ios::app);
		fout << "%PHI [deg]: " << myGapResult.phi << endl;
		for(int i = 0; i < N ; i++)
		{
			for(int j = 0; j < M ; j++)
            {
					fout << scientific << myNewtonIteration.myPistonGap.oildensity(j+i*M) << "\t" ;
            }
			fout << endl;
		}
		fout << "%" << endl;
		fout.close();
		fout.clear();
		*/


		//Piston gap pressure from convergence loop, .pK file [bar]
		double pDC = myGapResult.pDC;
		double pCase = myinput.data.operating_conditions.pCase*1e-5;
		fout.open("./output/piston/matlab/Piston_Gap_Pressure.txt",ios::app);
		fout << "%PHI [deg]: " << myGapResult.phi << endl;
		for(int i = 0; i < N ; i++)
		{
			for(int j = 0; j < M ; j++)
            {
				if(j==0)
				{
					fout << scientific << pDC << "\t" ;
				}
				else if(j==M-1)
				{
					fout << scientific << pCase << "\t" ;
				}
				else
				{
					fout << scientific << myNewtonIteration.myPistonGap.ploop(j+i*M)*1.0e-5 << "\t" ;
				}
            }
			fout << endl;
		}
		fout << "%" << endl;
		fout.close();
		fout.clear();


		//Piston test rig gap pressure, .pKEHD file [bar]
		if(myNewtonIteration.myPistonGap.EHDTestRig)
		{
			fout.open("./output/piston/matlab/Piston_Gap_Pressure_EHD.txt",ios::app);
			fout << "%PHI [deg]: " << myGapResult.phi << endl;
			for(int i = 0; i < 180 ; i++)
			{
				for(int j = 0; j < 9 ; j++)
				{
						fout << scientific << myNewtonIteration.myPistonGap.pEHD(j+i*9)*1e-5 << "\t" ;
				}
				fout << endl;
			}
			fout << "%" << endl;
			fout.close();
			fout.clear();
		}
		
		//Piston gap height, .hK file [microns]
		fout.open("./output/piston/matlab/Piston_Gap_Height.txt",ios::app);
		fout << "%PHI [deg]: " << myGapResult.phi << endl;
		for(int i = 0; i < N ; i++)
		{
			for(int j = 0; j < M ; j++)
            {
					fout << scientific << myNewtonIteration.myPistonGap.h(j+i*M)*1.0e6 <<  "\t" ;
            }
			fout << endl;
		}
		fout << "%" << endl;
		fout.close();
		fout.clear();

		//Piston gap surface temperature, .TK file [C]
		fout.open("./output/piston/matlab/Piston_Gap_Surface_Temperature.txt",ios::app);
		fout << "%PHI [deg]: " << myGapResult.phi << endl;
		for(int i = 0; i < N ; i++)
		{
			for(int j = 0; j < M ; j++)
            {
					fout << scientific << myNewtonIteration.myPistonGap.TK_surf_gap(j+i*M) <<  "\t" ;
            }
			fout << endl;
		}
		fout << "%" << endl;
		fout.close();
		fout.clear();

		//Piston cylinder gap surface temperature, .TB file [C]
		fout.open("./output/piston/matlab/Piston_Cylinder_Gap_Surface_Temperature.txt",ios::app);
		fout << "%PHI [deg]: " << myGapResult.phi << endl;
		for(int i = 0; i < N ; i++)
		{
			for(int j = 0; j < M ; j++)
            {
					fout << scientific << myNewtonIteration.myPistonGap.TB_surf_gap(j+i*M) <<  "\t" ;
            }
			fout << endl;
		}
		fout << "%" << endl;
		fout.close();
		fout.clear();

		//Surface deformations
		if(myinput.data.options_piston.general.ThermalDeformation || myinput.data.options_piston.general.PressureDeformation)
		{
			//Piston gap height deformation, .defK file [microns]
			fout.open("./output/piston/matlab/Piston_Gap_Deformation.txt",ios::app);
			fout << "%PHI [deg]: " << myGapResult.phi << endl;
			for(int i = 0; i < N ; i++)
			{
				for(int j = 0; j < M ; j++)
				{
						fout << scientific << myNewtonIteration.myPistonGap.defK_gap(j+i*M)*1e6 <<  "\t" ;
				}
				fout << endl;
			}
			fout << "%" << endl;
			fout.close();
			fout.clear();

			//Piston cylinder gap height deformation, .defK file [microns]
			fout.open("./output/piston/matlab/Piston_Cylinder_Gap_Deformation.txt",ios::app);
			fout << "%PHI [deg]: " << myGapResult.phi << endl;
			for(int i = 0; i < N ; i++)
			{
				for(int j = 0; j < M ; j++)
				{
						fout << scientific << myNewtonIteration.myPistonGap.defB_gap(j+i*M)*1e6 <<  "\t" ;
				}
				fout << endl;
			}
			fout << "%" << endl;
			fout.close();
			fout.clear();

			//Pressure deformation
			if(myinput.data.options_piston.general.PressureDeformation)
			{
				//Piston gap height pressure deformation, .defKp file [microns]
				fout.open("./output/piston/matlab/Piston_Gap_Pressure_Deformation.txt",ios::app);
				fout << "%PHI [deg]: " << myGapResult.phi << endl;
				for(int i = 0; i < N ; i++)
				{
					for(int j = 0; j < M ; j++)
					{
							fout << scientific << myNewtonIteration.myPistonGap.defK_p_gap(j+i*M)*1e6 <<  "\t" ;
					}
					fout << endl;
				}
				fout << "%" << endl;
				fout.close();
				fout.clear();

				//Piston cylinder gap height pressure deformation, .defBp file [microns]
				fout.open("./output/piston/matlab/Piston_Cylinder_Gap_Pressure_Deformation.txt",ios::app);
				fout << "%PHI [deg]: " << myGapResult.phi << endl;
				for(int i = 0; i < N ; i++)
				{
					for(int j = 0; j < M ; j++)
					{
							fout << scientific << myNewtonIteration.myPistonGap.defB_p_gap(j+i*M)*1e6 <<  "\t" ;
					}
					fout << endl;
				}
				fout << "%" << endl;
				fout.close();
				fout.clear();

				//Piston surface elatic deformation due to pressure stress, .defKpb file [microns]
				fout.open("./output/piston/matlab/Piston_Pressure_Deformation.txt",ios::app);
				fout << "%PHI [deg]: " << myGapResult.phi << endl;
				for(int i=0;i<(int) myNewtonIteration.myPistonGap.defK_p.extent(0);i++)
				{
					fout << scientific << myNewtonIteration.myPistonGap.defK_p(i) << endl;
				}
				fout << "%" << endl;
				fout << "%" << endl;
				fout.close();
				fout.clear();

				//Piston cylinder surface elatic deformation due to pressure stress, .defBpb file [microns]
				fout.open("./output/piston/matlab/Piston_Cylinder_Pressure_Deformation.txt",ios::app);
				fout << "%PHI [deg]: " << myGapResult.phi << endl;
				for(int i=0;i<(int) myNewtonIteration.myPistonGap.defB_p.extent(0);i++)
					{
						fout << scientific << myNewtonIteration.myPistonGap.defB_p(i) << endl;
					}
					fout << "%" << endl;
				fout << "%" << endl;
				fout.close();
				fout.clear();
			}

			//Thermal deformation
			if(myinput.data.options_piston.general.ThermalDeformation)
			{
				//Piston gap height thermal deformation, .defKth file [microns]
				fout.open("./output/piston/matlab/Piston_Gap_Thermal_Deformation.txt",ios::app);
				fout << "%PHI [deg]: " << myGapResult.phi << endl;
				for(int i = 0; i < N ; i++)
				{
					for(int j = 0; j < M ; j++)
					{
							fout << scientific << myNewtonIteration.myPistonGap.defK_th_gap(j+i*M)*1e6 <<  "\t" ;
					}
					fout << endl;
				}
				fout << "%" << endl;
				fout.close();
				fout.clear();

				//Piston cylinder gap height thermal deformation, .defBth file [microns]
				fout.open("./output/piston/matlab/Piston_Cylinder_Gap_Thermal_Deformation.txt",ios::app);
				fout << "%PHI [deg]: " << myGapResult.phi << endl;
				for(int i = 0; i < N ; i++)
				{
					for(int j = 0; j < M ; j++)
					{
							fout << scientific << myNewtonIteration.myPistonGap.defB_th_gap(j+i*M)*1e6 <<  "\t" ;
					}
					fout << endl;
				}
				fout << "%" << endl;
				fout.close();
				fout.clear();
			}

		};

		//Output vtk (vectors are outside first, then inside).
		static int vtkcount = -1;
		static int vtknum = ceil(360 / myinput.data.options_piston.numeric.Simalphaplot);//Number of VTK files to retain not to exceed one revolution to minimize file size.
		short Mcells, Ncells;
		if(myinput.data.options_piston.general.ReynoldsMultiGrid){
			Mcells = myinput.data.options_piston.fluid_grid.MG.MG_M[0];
			Ncells = myinput.data.options_piston.fluid_grid.MG.MG_N[0];
		}
		else{
			Mcells = myinput.data.options_piston.fluid_grid.GS.M;
			Ncells = myinput.data.options_piston.fluid_grid.GS.N;
		}


		vtkcount += 1;
		std::ostringstream strs;
		strs << vtkcount;
		std::string str = strs.str();
		string output("./output/piston/vtk/piston_gap." + str + ".vtk");
		fout.open(output.c_str());
		fout << "# vtk DataFile Version 2.0" << endl;
		fout << "vtk output" << endl;
		fout << "ASCII" << endl;
		fout << "DATASET STRUCTURED_GRID" << endl;
		fout << "DIMENSIONS " << M << " " << N+1 << " 2" << endl;
		fout << "POINTS " << 2*Mcells*(Ncells+1) << " double" << endl;
		fout.precision(8);
		for(int c=0;c<2;c++)
		{
			for(int i=0;i<Ncells;i++)
			{
				for(int j=0;j<Mcells;j++)
				{
					fout << scientific <<  myNewtonIteration.myPistonGap.xyzf_gap(j+i*Mcells,0) << "\t" << myNewtonIteration.myPistonGap.xyzf_gap(j+i*Mcells,1) << "\t" << myNewtonIteration.myPistonGap.xyzf_gap(j+i*Mcells,2) + myNewtonIteration.myPistonGap.geometrypistongap.lB << endl;
				}
				
			}
			for(int j=0;j<Mcells;j++)
			{
				fout << scientific <<  myNewtonIteration.myPistonGap.xyzf_gap(j,0) << "\t" << myNewtonIteration.myPistonGap.xyzf_gap(j,1) << "\t" << myNewtonIteration.myPistonGap.xyzf_gap(j,2) + myNewtonIteration.myPistonGap.geometrypistongap.lB << endl;
			}
		}
		fout << " " << endl;
		fout << "POINT_DATA " << 2*Mcells*(Ncells+1)<<endl;
		fout << "VECTORS gap_height double"<<endl;
		for(int i=0;i<Ncells;i++)
		{
			for(int j=0;j<Mcells;j++)
			{
				fout << scientific << cos(i*2*PI/Ncells) * ( myNewtonIteration.myPistonGap.defB_gap(j+i*Mcells) - myNewtonIteration.myPistonGap.McrB(j+i*Mcells) ) << "\t"
				<< sin(i*2*PI/Ncells) * ( myNewtonIteration.myPistonGap.defB_gap(j+i*Mcells) - myNewtonIteration.myPistonGap.McrB(j+i*Mcells) ) << "\t0" << endl;
			}
		}
		for(int j=0;j<Mcells;j++)
		{
			fout << scientific << ( myNewtonIteration.myPistonGap.defB_gap(j) - myNewtonIteration.myPistonGap.McrB(j) ) << "\t0\t0" << endl;
		}
		for(int i=0;i<Ncells;i++)
		{
			for(int j=0;j<Mcells;j++)
			{
				fout << scientific << cos(i*2*PI/Ncells) * ( myNewtonIteration.myPistonGap.defB_gap(j+i*Mcells) - myNewtonIteration.myPistonGap.McrB(j+i*Mcells) - myNewtonIteration.myPistonGap.h(j+i*Mcells) ) << "\t"
				<< sin(i*2*PI/Ncells) * ( myNewtonIteration.myPistonGap.defB_gap(j+i*Mcells) - myNewtonIteration.myPistonGap.McrB(j+i*Mcells) - myNewtonIteration.myPistonGap.h(j+i*Mcells) ) << "\t0" << endl;
			}
		}
		for(int j=0;j<Mcells;j++)
		{
			fout << scientific << ( myNewtonIteration.myPistonGap.defB_gap(j) - myNewtonIteration.myPistonGap.McrB(j) - myNewtonIteration.myPistonGap.h(j) ) << "\t0\t0" << endl;
		}
		fout << endl;
		fout << "VECTORS pressure_deformation double" << endl;
		for(int i = 0; i < Ncells ; i++)
		{
			for(int j = 0; j < Mcells ; j++)
			{
				fout << scientific << cos(i*2*PI/Ncells) * myNewtonIteration.myPistonGap.defB_p_gap(j+i*Mcells) <<  "\t" << sin(i*2*PI/Ncells) * myNewtonIteration.myPistonGap.defB_p_gap(j+i*Mcells) << "\t0" << endl ;
			}
		}
		for(int j = 0; j < Mcells ; j++)
		{
			fout << scientific << myNewtonIteration.myPistonGap.defB_p_gap(j) <<  "\t0\t0" << endl ;
		}
		for(int i = 0; i < Ncells ; i++)
		{
			for(int j = 0; j < Mcells ; j++)
			{
				fout << scientific << cos(i*2*PI/Ncells) * myNewtonIteration.myPistonGap.defK_p_gap(j+i*Mcells) << "\t" << sin(i*2*PI/Ncells) * myNewtonIteration.myPistonGap.defK_p_gap(j+i*Mcells) <<  "\t0" << endl ;
			}
		}
		for(int j = 0; j < Mcells ; j++)
		{
			fout << scientific << myNewtonIteration.myPistonGap.defK_p_gap(j) << "\t0\t0" << endl ;
		}
		fout << endl;
		fout << "VECTORS thermal_deformation double" << endl;
		for(int i = 0; i < Ncells ; i++)
		{
			for(int j = 0; j < Mcells ; j++)
			{
				fout << scientific << cos(i*2*PI/Ncells) * myNewtonIteration.myPistonGap.defB_th_gap(j+i*Mcells) << "\t" << sin(i*2*PI/Ncells) * myNewtonIteration.myPistonGap.defB_th_gap(j+i*Mcells) << "\t0" << endl ;
			}
		}
		for(int j = 0; j < Mcells ; j++)
		{
			fout << scientific << myNewtonIteration.myPistonGap.defB_th_gap(j)<< "\t0\t0" << endl ;
		}
		for(int i = 0; i < Ncells ; i++)
		{
			for(int j = 0; j < Mcells ; j++)
			{
				fout << scientific << cos(i*2*PI/Ncells) * myNewtonIteration.myPistonGap.defK_th_gap(j+i*Mcells) << "\t" << sin(i*2*PI/Ncells) * myNewtonIteration.myPistonGap.defK_th_gap(j+i*Mcells) << "\t0" << endl ;
			}
		}
		for(int j = 0; j < Mcells ; j++)
		{
			fout << scientific << myNewtonIteration.myPistonGap.defK_th_gap(j) << "\t0\t0" << endl ;
		}
		fout << endl;
		fout << "VECTORS friction_cylinder double" << endl;
		//int temp = myNewtonIteration.myPistonGap.forcespistongap.FTbx.size();
		//Log << "size FTbx = " << temp << endl;
		for(int i = 0; i < Ncells ; i++)
		{
			for(int j = 0; j < Mcells ; j++)
			{
				double a = myNewtonIteration.myPistonGap.forcespistongap.FTbx(j+i*Mcells);
				double b = myNewtonIteration.myPistonGap.forcespistongap.FTby(j+i*Mcells);
				fout << scientific << -1.0* sin(i*2*PI/Ncells) * a << "\t" << cos(i*2*PI/Ncells) * a << "\t" << b << endl ;
			}
		}
		for(int j = 0; j < Mcells ; j++)
		{
			fout << scientific << "0\t" << myNewtonIteration.myPistonGap.forcespistongap.FTbx(j) << "\t" << myNewtonIteration.myPistonGap.forcespistongap.FTby(j) << endl ;
		}
		for(int i = 0; i < Ncells ; i++)
		{
			for(int j = 0; j < Mcells ; j++)
			{
				fout << "0\t0\t0" << endl ;
			}
		}
		for(int j = 0; j < Mcells ; j++)
		{
			fout << scientific << "0\t0\t0" << endl ;
		}
		fout << endl;
		fout << "VECTORS friction_piston double" << endl;
		for(int i = 0; i < Ncells ; i++)
		{
			for(int j = 0; j < Mcells ; j++)
			{
				fout << "0\t0\t0" << endl ;
			}
		}
		for(int j = 0; j < Mcells ; j++)
		{
			fout << scientific << "0\t0\t0" << endl ;
		}
		for(int i = 0; i < Ncells ; i++)
		{
			for(int j = 0; j < Mcells ; j++)
			{
				fout << scientific << -1.0 * sin(i*2*PI/Ncells) * myNewtonIteration.myPistonGap.forcespistongap.FTKx(j+i*Mcells) << "\t" << cos(i*2*PI/Ncells) * myNewtonIteration.myPistonGap.forcespistongap.FTKx(j+i*Mcells) << "\t" << myNewtonIteration.myPistonGap.forcespistongap.FTKy(j+i*Mcells) << endl ;
			}
		}
		for(int j = 0; j < Mcells ; j++)
		{
			fout << scientific << "0\t" << myNewtonIteration.myPistonGap.forcespistongap.FTKx(j) << "\t" << myNewtonIteration.myPistonGap.forcespistongap.FTKy(j) << endl ;
		}
		fout << endl;
		//fout << "POINT DATA" << 2*M*(N+1) << endl;
		fout << "SCALARS gap_pressure double 1" << endl;
		fout << "LOOKUP_TABLE default" << endl;
		for(int c = 0; c < 2 ; c++)
		{
			for(int i = 0; i < Ncells ; i++)
			{
				for(int j = 0; j < Mcells ; j++)
			    {
					if(j==0)
					{
						fout << scientific << pDC << "\t" ;
					}
					else if(j==Mcells-1)
					{
						fout << scientific << pCase << "\t" ;
					}
					else
					{
						fout << scientific << myNewtonIteration.myPistonGap.ploop(j+i*Mcells)*1.0e-5 << "\t" ;
					}
			    }
				fout << endl;
			}
			for(int j = 0; j < Mcells ; j++)
			{
				if(j==0)
				{
					fout << scientific << pDC << "\t" ;
				}
				else if(j==Mcells-1)
				{
					fout << scientific << pCase << "\t" ;
				}
				else
				{
					fout << scientific << myNewtonIteration.myPistonGap.ploop(j)*1.0e-5 << "\t" ;
				}
			}
			fout<<endl;
		}
		fout << "SCALARS surface_temperature double 1" << endl;
		fout << "LOOKUP_TABLE default" << endl;
		for(int i = 0; i < Ncells ; i++)
		{
			for(int j = 0; j < Mcells ; j++)
            {
					fout << scientific << myNewtonIteration.myPistonGap.TB_surf_gap(j+i*Mcells) <<  "\t" ;
            }
			fout << endl;
		}
		for(int j = 0; j < Mcells ; j++)
        {
			fout << scientific << myNewtonIteration.myPistonGap.TB_surf_gap(j) <<  "\t" ;
        }
		fout << endl;
		for(int i = 0; i < Ncells ; i++)
		{
			for(int j = 0; j < Mcells ; j++)
            {
					fout << scientific << myNewtonIteration.myPistonGap.TK_surf_gap(j+i*Mcells) <<  "\t" ;
            }
			fout << endl;
		}
		for(int j = 0; j < Mcells ; j++)
        {
			fout << scientific << myNewtonIteration.myPistonGap.TK_surf_gap(j) <<  "\t" ;
        }
		fout << endl;
		fout.close();
		fout.clear();

		if(vtkcount > (vtknum-1)) //remove prior VTK files.
		{
			std::ostringstream strs;
			strs << vtkcount - vtknum;
			std::string str = strs.str();
			string output("./output/piston/vtk/piston_gap." + str + ".vtk");
			remove(output.c_str());
		};

		xd2D += dx2D;

		/*std::ostringstream strs;
		strs << vtkcount;
		std::string str = strs.str();
		string output("./outputs_piston/vtk/piston_gap." + str + ".vtk");
		fout.open(output.c_str());*/
	};



	return 0;

};