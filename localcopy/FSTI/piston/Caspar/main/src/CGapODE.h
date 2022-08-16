#include "CNewtonIteration.h"
#include "CGapInput.h"
#include "sGapResult.h"
#include <ppl.h>
#pragma once

#include <iostream>
using namespace Concurrency;
using namespace std;


class CGapODE
{
public:
	
	//Default constructor-destructor
	CGapODE(void);
	~CGapODE(void);
	
	//Instances of other classes used in CGapODE
	CNewtonIteration myNewtonIteration;
	

	//----------------------DECLARATION OF OFSTREAM----------------------
	ofstream fout;

	//-------------------------SOLVER VARIABLES--------------------------
	unsigned short n;				// Dimension of the system
	double x;			// Initial value of dependent variable (usually time)
	double xbeg;		// Beginning x-value
	double xend;		// Final x value (xend-x may be positive or negative).
	double dx;			// time step for intermediate output
	double dx2D;		// time step for intermediate 2D and 3D output
	// Derived variables
	double hmax;		// maximal step size
	// stores past value of x
	double xold;
	// stores past value of h
	double hold;
	// x at discrete points specified by dx interval
	double xd;
	double xd2D;

	double PhiD_mech_tot;
	double PhiD_vol_tot;
	// revolution counter
	unsigned short revcounter;
	// y values - positons
	vector<double> y; 
	// yp values - velocities
	vector<double> yp; 
	// time vector
	vector<double> xi;

	//counters
	int vtkcount;
	int revcounter_GUI;
	int revcounter_Convergance;
	int resumecount;

	//=======================SOLVER METHODS===============================//
	//Start the integration process
	void ODEmain(void);
	//Setup all the ODE solver variables
	void ODESetupSolver(void);	
	//ODE definition, used directly by the Runge-Kutta solver methods (input: time [s],postions [m],velocities [m/s])
	void odef(double x,vector<double> &y,vector<double> &yp);
	//Call core integration function
	void ODEIntegrate();
	//Adap time step
	void ODEAdaptTimeStep();
	//Control the output of the results
	int ODESolutionOutput();

};

