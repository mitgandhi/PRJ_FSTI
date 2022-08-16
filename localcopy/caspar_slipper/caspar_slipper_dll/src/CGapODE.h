#include "CNewtonIteration.h"
#include "CGapLog.h"
#include "CGapResult.h"
#include "../../caspar_input/caspar_input.h"
#pragma once

using namespace std;

class CGapODE
{
public:

	//Gap inputs
	caspar_input * gapinput;
	
	//Default constructor-destructor
	CGapODE(caspar_input * gapinputs);
	~CGapODE(void);

	//Instance of the NewtonIneration
	CNewtonIteration NewtonIteration;

	//-------------------------SOLVER VARIABLES--------------------------
	int n;										// Dimension of the system
	double timeStart;							// Beginning time
	
	double time;								// Time vector
	vector<double> p;							// Gap position vector
	vector<double> v;							// Gap velocity vector

	vector<double> pold;						// Gap old position vector
	vector<double> vold;						// Gap old velocity vector


	//=======================SOLVER METHODS===============================//
	void ODEmain(const bool Resume);
	void ODESetupSolver(const bool Resume);	
	void ODECleanSolver(void);
	void ODEIntegrate();

	void ODEResumeWrite(const bool DirtyPlot);
	bool ODEResumeRead(void);

	void updateTime(const double time_s);
	
};

