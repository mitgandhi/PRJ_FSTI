#include "CPistonGap.h"
#include "CGapInput.h"
#include "sGapResult.h"
#pragma once

class CNewtonIteration
{

public:
	
	//Default constructor-destructor
	CNewtonIteration(void);
	~CNewtonIteration(void);
	
	//Instances of the classes used by Newton Iteration
	CPistonGap myPistonGap;		//PistonGap class as variable for Newton class

	//Main variables passed through Newton Iteration functions
	double epsilonK,delta_v,dt,revtime;	//Piston force balance tolerance [N]
	int jmax,kmax;


	


	//----------------------------------------------//
	//------------FUNCTIONS DECLARATION------------//
	//---------------------------------------------//
	//Main function to call Newton Iteration for the 3 gaps
	void NewtonCalcNewtonIteration(vector<double> &gappositions,vector<double> &gapvelocities);
	//Newton Iteration core function for shifting velocity calculation
	void NewtonCalcShiftingVelocities(double epsilon,vector<double> &gappositions,vector<double> &gapvelocities);
	//Newton function to calculate euclidean norm of the dF vector
	double NewtonCalcL2norm(int n,vector<double> &dF); 
	//Newton function to calculate dF tolerance for force balance convergence
	void NewtonCalcEpsilonK(void); 
	//Jacobian matrix calculation for force balance
	int NewtonCalcJacobian(int n,vector<double> &gappositions,vector<double> &gapvelocities,vector< vector<double> > &Jacobian,vector<double> &dF,vector<double> &tempdF);
	//Main Gauss solver function calling the desired solution procedure
	int NewtonCalcGaussMain(int mod,int n,vector< vector<double> > &mat,vector< vector<double> > &lumat,vector<int> &perm,vector<double> &b,vector<double> &x,int *signd);
	//Function to calculate LU matrix decomposition
	int NewtonCalcGaussDecomposition(int n,vector< vector<double> > &mat,vector< vector<double> > &lumat,vector<int> &perm,int *signd);
	//Function to calculate linear system of equation solution
	int NewtonCalcGaussSolution(int n,vector< vector<double> > &lumat,vector<int> &perm,vector<double> &b,vector<double> &x);
	//Solve Newton loop
	void NewtonLoop(double epsilon,vector<double> &gappositions,vector<double> &gapvelocities);
};


