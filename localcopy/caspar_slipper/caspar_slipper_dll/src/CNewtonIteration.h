#include "CSlipperGap.h"
#include "../../caspar_input/caspar_input.h"
#include "CGapResult.h"
#include "CGapLog.h"
#include "time.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#pragma once

extern CGapLog GapLog;

class CNewtonIteration
{

public:

	//Gap inputs
	caspar_input * gapinput;

	//Default constructor-destructor
	CNewtonIteration(caspar_input * gapinputs);
	~CNewtonIteration(void);
	
	//Instances of the classes used by Newton Iteration
	CGapResult GapResult;							//Result Class
	CSlipperGap* SlipperGap;						//SlipperGap Class

	//This is a structure used to pass data into HybridCallFunc
	struct hybrid_params
	{
		CNewtonIteration * newton;
		vector<double> gappositions;
		size_t n;
		int ComplexPicard;

		hybrid_params()
		{
			ComplexPicard = 0;
		}
	};
	
	//Main variables passed through Newton Iteration functions
	double epsilonG;	//Slipper force balance tolerance [N]

	//----------------------------------------------//
	//------------FUNCTIONS DECLARATION------------//
	//---------------------------------------------//
	//Main function to call Newton Iteration for the 3 gaps
	int NewtonCalcNewtonIteration(vector<double> &gappositions,vector<double> &gapvelocities,  vector<double> &pold,vector<double> &vold);
	//Newton Iteration core function for shifting velocity calculation
	int NewtonCalcShiftingVelocities(int gapindex,int n,double &epsilon,vector<double> &gappositions,vector<double> &gapvelocities, double & dFnorm);
	//Uses the GSL hybrid method instead of Newton
	int HybridCalcShiftingVelocities(int gapindex,int n,double &epsilon,vector<double> &gappositions,vector<double> &gapvelocities, double & dFnorm);
	//Used by HybridCalcShiftingVelocities to call the Slipper
	static int HybridCallFunc(const gsl_vector * x, void *params, gsl_vector * f);
	static int HybridCalldFunc(const gsl_vector * x, void *params, gsl_matrix * J);
	static int HybridCallFuncdFunc(const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J);
	//Function used to calculate the euclidean norm of a gsl vector. I'm sure the GSL library has some way to do this with BLAS, but I didn't feel like fighting with the manual.
	static double gsl_vector_norm(const gsl_vector * v);
	//Newton function to calculate euclidean norm of the dF vector
	double NewtonCalcL2norm(int n,vector<double> &dF); 
	//Newton function to calculate dF tolerance for force balance convergence
	void NewtonCalcEpsilonK(void); 
	//Jacobian matrix calculation for force balance
	int NewtonCalcJacobian(int gapindex,int n,vector<double> &gappositions,vector<double> &gapvelocities,vector< vector<double> > &Jacobian,vector<double> &dF,vector<double> &tempdF);
	//Main Gauss solver function calling the desired solution procedure
	int NewtonCalcGaussMain(int mod,int n,vector< vector<double> > &mat,vector< vector<double> > &lumat,vector<int> &perm,vector<double> &b,vector<double> &x,int *signd);
	//Function to calculate LU matrix decomposition
	int NewtonCalcGaussDecomposition(int n,vector< vector<double> > &mat,vector< vector<double> > &lumat,vector<int> &perm,int *signd);
	//Function to calculate linear system of equation solution
	int NewtonCalcGaussSolution(int n,vector< vector<double> > &lumat,vector<int> &perm,vector<double> &b,vector<double> &x);
	//Pressure Functions
	vector<double> SetPressure(const double phi_deg);
	vector<double> GetIdealPressure(double phi_deg);
	vector<double> GetFilePressure(double phi_deg);

	//A simple sign function
	int sgn(const double x);

};


