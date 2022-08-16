#include "CNewtonIteration.h"
#include "logger.h"
#include <iostream>
#include <iomanip>
#include "../../caspar_input/input.h"
#include <math.h>
#pragma once

extern class CGapInput myGapInput;
extern struct sGapResult myGapResult;
extern class input myinput;

CNewtonIteration::CNewtonIteration(void)
{

	epsilonK = 100.0;
	jmax = myinput.data.options_piston.numeric.jmax;
	kmax = myinput.data.options_piston.numeric.kmax;
	delta_v = myinput.data.options_piston.numeric.delta_v;

	//revolution time
	revtime = 2*PI/myinput.data.operating_conditions.speed ;

};
CNewtonIteration::~CNewtonIteration(void)
{

};
//Perform Newton iteration
void CNewtonIteration::NewtonCalcNewtonIteration(vector<double> &gappositions, vector<double> &gapvelocities)
{
	//Temporary vectors
	vector<double> positionpiston(4);
	vector<double> velocitypiston(4);

	//Assign parts positions to the respective vectors from global positions vector
	double time = gappositions[0];
	double tol = 1.0e-7;
	for(int i=0;i<4;i++)
	{
		positionpiston[i] = gappositions[i+1];
		velocitypiston[i] = gapvelocities[i];
	}

	//Piston initialize main variables for gap calculations
	myPistonGap.PistonInitializeVariables(time);
	
	//Calculate the piston shifting velocities
	NewtonCalcShiftingVelocities(epsilonK,positionpiston,velocitypiston);
	
	//Calculate force balance tolerance based on external forces
	NewtonCalcEpsilonK();
	
	//Calculate energy thermal flux to solid bodies
	if(myPistonGap.HeatTransfer)
	{
		myPistonGap.PistonCylinderCalcBodyThermalFlux();
		myPistonGap.PistonCylinderConstructWeightedMatrix();
	}

	//Solve for thermal problem at the end of each revolution
	if( time > (myGapResult.revcounter*revtime-tol) )
	{
		if(myPistonGap.HeatTransfer)
		{
			if(coupled_caspar_simulation){
				//Check for updated head flux from valveplate gap.
				myGapInput.readVPFlux();
				//modify thermal boundary conditions accordingly.
				if(myGapInput.VPFluxFile.qvp != 0.0){
					for(int i = 0;i < myinput.data.thermal.block.dirichlet_bc.size();i++){
						if(myinput.data.thermal.block.dirichlet_bc[i].setname.compare("gap_block") == 0){
							myinput.data.thermal.block.dirichlet_bc.erase(myinput.data.thermal.block.dirichlet_bc.begin()+i);//this may not work correctly dwm
						}
					}
					for(int i = 0;i < myinput.data.thermal.block.mixed_bc.size();i++){
						if(myinput.data.thermal.block.mixed_bc[i].setname.compare("gap_block") == 0){
							myinput.data.thermal.block.mixed_bc.erase(myinput.data.thermal.block.mixed_bc.begin()+i);//this may not work correctly dwm
						}
					}
					for(int i = 0;i < myinput.data.thermal.block.neumann_bc.size();i++){
						if(myinput.data.thermal.block.neumann_bc[i].setname.compare("gap_block") == 0){
							myinput.data.thermal.block.neumann_bc.erase(myinput.data.thermal.block.neumann_bc.begin()+i);//this may not work correctly dwm
						}
					}
					
					myinput.data.thermal.block.neumann_bc.resize(myinput.data.thermal.block.neumann_bc.size()+1);
					myinput.data.thermal.block.neumann_bc[myinput.data.thermal.block.neumann_bc.size()-1].setname = "gap_block";
					myinput.data.thermal.block.neumann_bc[myinput.data.thermal.block.neumann_bc.size()-1].q.push_back(myGapInput.VPFluxFile.qvp);
				}
			}

			myPistonGap.PistonCylinderSolveBodyThermal_new();

			if(coupled_caspar_simulation){
				//Compute and output bore average flux.
				double fluxout = 0.0;
				for(int i = 0;i<myPistonGap.qbi_cylinder[0].size();i++){
					fluxout += myPistonGap.qbi_cylinder[0](i);
				}
				fluxout /= myPistonGap.qbi_cylinder[0].size();

				remove("./output/piston/piston_flux.txt");

				ofstream fout;

				fout.open("./output/piston/piston_flux.txt",ios::app);
				fout << scientific << fluxout;
				fout.close();
				fout.clear();
			}
		}
		//Revolution counter
		myGapResult.revcounter++;
	};

	//Assign the velocities to the global velocities vector
	for(int i=0;i<4;i++)
	{
		gapvelocities[i] = velocitypiston[i];
	};
	
	//Results assignment
	myGapResult.time = time;								
	myGapResult.phi = myPistonGap.operatingpistongap.phi_deg;
	myGapResult.gappositions = gappositions;
	myGapResult.gapvelocities = gapvelocities;


};
//Calculate squeeze motion using Newton relaxed method
void CNewtonIteration::NewtonCalcShiftingVelocities(double epsilon,vector<double> &gappositions_in,vector<double> &gapvelocities_in)
{

	//Assign previous velocity
	//gapvelocities = gapvelocities_0;

	//Shift velocities and positions according to block rotation. dwm
	double dphi = myinput.data.geometry.speedK * myinput.data.operating_conditions.speed * dt; //rad/s * sec
	vector<double> temp(4);//temp[1] = dphi * 180 / 3.14157;
	//Log << "Shifting Velocities by " << temp[1] << " degrees" << "\n";
	temp = gapvelocities_in;
	gapvelocities_in[0] = cos(dphi) * temp[0] - sin(dphi) * temp[1];
	gapvelocities_in[1] = cos(dphi) * temp[1] + sin(dphi) * temp[0];
	gapvelocities_in[2] = cos(dphi) * temp[2] - sin(dphi) * temp[3];
	gapvelocities_in[3] = cos(dphi) * temp[3] + sin(dphi) * temp[2];
	temp = gappositions_in;
	gappositions_in[0] = cos(dphi) * temp[0] - sin(dphi) * temp[1];
	gappositions_in[1] = cos(dphi) * temp[1] + sin(dphi) * temp[0];
	gappositions_in[2] = cos(dphi) * temp[2] - sin(dphi) * temp[3];
	gappositions_in[3] = cos(dphi) * temp[3] + sin(dphi) * temp[2];

	myPistonGap.IMParallel(0,0,0);
	
	//Solve main convergence loop
	myPistonGap.PistonGap(gappositions_in,gapvelocities_in,dt);

	if(myinput.data.options_piston.general.PressureDeformation && myinput.data.options_piston.general.PressureDeformationOMP == 0)
		myPistonGap.IMParallel(2,0,0);

	//Call Newton iteration loop
	
	NewtonLoop(epsilon,gappositions_in,gapvelocities_in);

	//Assign current velocity
	//gapvelocities_0 = gapvelocities;
		
};
//Calculate problem jacobian matrix
int CNewtonIteration::NewtonCalcJacobian(int n,vector<double> &gappositions,vector<double> &gapvelocities,vector< vector<double> > &Jacobian,vector<double> &dF,vector<double> &tempdF)
{

	//Jacobian matrix calculation: each term of the matrix expressed with central difference as (dF - tempdF)/(2*h)
	//where dF and tempdF are evaluated at (gapvelocities + h) and (gapvelocities - h)
	double  xj,denom;

	for (int j = 0; j < n; j++)
	{
		//Velocity
		xj = gapvelocities[j];

		//Delta of speed
		denom = 1.0/(2.0*delta_v);

		//Forward velocity
		//if(gappositions[j] < 0.0)
			gapvelocities[j] += delta_v;      
		
		//Vector tempdF is returned
		myPistonGap.PistonGapJacobian(gappositions,gapvelocities,tempdF,dt);
		
		//Reset velocity
		gapvelocities[j] = xj;

		//Backward velocity
		//if(gappositions[j] >= 0.0)
			gapvelocities[j] -= delta_v;
		
		//Vector dF is returned
		myPistonGap.PistonGapJacobian(gappositions,gapvelocities,dF,dt);
		
		//Reset velocity
		gapvelocities[j] = xj;

		//Jacobian
		for (int i = 0; i < n; i++)
		{
			Jacobian[i][j] = (tempdF[i] - dF[i]) * denom;
		}
	
	}

	return 0;

};
//Solve linear system usign Gauss back substitution to determine delta velocity
int CNewtonIteration::NewtonCalcGaussMain(int mod,int n,vector< vector<double> > &mat,vector< vector<double> > &lumat,vector<int> &perm,vector<double> &b,vector<double> &x,int *signd)
{
	
	/*====================================================================*
	 *                                                                    *
	 *  The funktion gauss solves a linear system :  mat * x = b.         *
	 *  Here mat is the nonsingular system matrix, b the right hand side  *
	 *  of the system and x the solution vector.                          *
	 *                                                                    *
	 *  gauss uses the Gauss algorithm and computes a triangular factori- *
	 *  zation of mat and skaled column pivot search.  (Crout method with *
	 *  row swaps).                                                       *
	 *                                                                    *
	 *====================================================================*
	.BE*)
	 *                                                                    *
	 *   Application:                                                     *
	 *   ============                                                     *
	 *      Solve general linear system with a nonsingular coefficient    *
	 *      matrix.                                                       *
	 *                                                                    *
	 *====================================================================*
	 *                                                                    *
	 *   Control parameter:                                               *
	 *   ==================                                               *
	 *      mod      int mod;                                             *
	 *               calling modus for gauss:                             *
	 *       = 0     Find factorization and solve linear system           *
	 *       = 1     Find factorization only.                             *
	 *       = 2     Solve linear system only; the factorization is       *
	 *               already available in lumat. This saves work when     *
	 *               solving a linear system repeatedly for several right *
	 *               hand sides and the same system matrix such as when   *
	 *               inverting the matrix.                                *
	 *       = 3     as under 2, additionally we improve the solution     *
	 *               via iterative refinement.                            *
	 *                                                                    *
	 *   Input parameters:                                                *
	 *   ================                                                 *
	 *      n        int n;  ( n > 0 )                                    *
	 *               Dimension of mat and lumat,                          *
	 *               size of the vector b, the right hand side, the       *
	 *               solution x and the permutation vector perm.          *
	 *      mat      mat[];                                       *
	 *               Matrix of the linear system. It is stored in vector  *
	 *               form.                                                *
	 *      lumat    lumat[];      ( for mod = 2, 3 )             *
	 *               LU factors of mat                                    *
	 *               mat in eine untere und obere Dreieckmatrix ent-      *
	 *               lumat can be stored in the space of mat.             *
	 *      perm     int perm[];           ( for mod = 2, 3 )             *
	 *               Permutation vector, of the row interchangfes in lumat*
	 *      b        REAL   b[];           ( for mod = 0, 2, 3 )          *
	 *               Right hand side of the system.                       *
	 *      signd    int *signd;           ( for mod = 2, 3 )             *
	 *               sign of the permutation in perm; the determinant of  *
	 *               mat can be computed as the product of the diagonal   *
	 *               entries of lumat times signd.                        *
	 *                                                                    *
	 *   Output parameters:                                               *
	 *   ==================                                               *
	 *      lumat    REAL   *lumat[];      ( for mod = 0, 1 )             *
	 *               LU factorization of mat.                             *
	 *      perm     int perm[];           ( for mod = 0, 1 )             *
	 *               row ermutation vektor                                *
	 *      x        REAL   x[];           ( for mod = 0, 2, 3 )          *
	 *               solution vector.                                     *
	 *      signd    int *signd;           ( for mod = 0, 1 )             *
	 *               sign of perm.                                        *
	 *                                                                    *
	 *   Return value :                                                   *
	 *   ==============                                                   *
	 *      =-1      Max. number (MAXITER) of iterative refinements       *
	 *               reached (MAXITER) while mod = 3                      *
	 *      = 0      all ok                                               *
	 *      = 1      n < 1 or other invalid input                         *
	 *      = 2      lack of memory                                       *
	 *      = 3      Matrix singular                                      *
	 *      = 4      Matrix numerically singular                          *
	 *      = 5      incorrect call                                       *
	 *                                                                    *
	 *====================================================================*
	 *                                                                    *
	 *   Functions used :                                                 *
	 *   ================                                                 *
	 *                                                                    *
	 *      int gaudec  (): determines  LU decomposition                  *
	 *      int gausol  (): solves the linear system                      *
	 *                                                                    *
	 *====================================================================*/
	int  rc;
	
	if (n < 1) return (1);
	
	switch (mod)
	{
	case 0: /* Find factorization and solve system ...................*/
		rc = NewtonCalcGaussDecomposition(n,mat,lumat,perm,signd);
		if (rc == 0)
			return (NewtonCalcGaussSolution(n,lumat,perm,b,x));
		else
			return (rc);
	case 1: /* Find factorization only ...............................*/
		return (NewtonCalcGaussDecomposition(n,mat,lumat,perm,signd));
	case 2: /* Solve only ............................................*/
		return (NewtonCalcGaussSolution(n,lumat,perm,b,x));
	}
	
	return (5);                                           /* Wrong call */

};
int CNewtonIteration::NewtonCalcGaussDecomposition(int n,vector<vector<double> > &mat,vector<vector<double> > &lumat,vector<int> &perm, int *signd)
{
	/*====================================================================*
	 *                                                                    *
	 *  gaudec decomposes a nonsingular n x n matrix into a product of a  *
	 *  unit lower and an upper triangular matrix. Both triangular factors*
	 *  are stored in lumat (minus the unit diagonal, of course).         *
	 *                                                                    *
	 *====================================================================*
	 *                                                                    *
	 *   Eingabeparameter:                                                *
	 *   ================                                                 *
	 *      n        int n;  ( n > 0 )                                    *
	 *               Dimension of  mat and lumat,                         *
	 *               size of  b , x and perm                              *
	 *      mat      REAL   *mat[];                                       *
	 *               original system matrix in vector form.               *
	 *                                                                    *
	 *   Output parameters:                                               *
	 *   ==================                                               *
	 *      lumat    REAL   *lumat[];                                     *
	 *               LU factorization                                     *
	 *      perm     int perm[];                                          *
	 *               row permutation vector for lumat                     *
	 *      signd    int *signd;                                          *
	 *               sign of perm. The determinant of mat can be computed *
	 *               as the product of the diagonal entreis of lumat times*
	 *               signd.                                               *
	 *                                                                    *
	 *   Return value:                                                    *
	 *   =============                                                    *
	 *      = 0      all ok                                               *
	 *      = 1      n < 1 or invalid input                               *
	 *      = 2      lack of memory                                       *
	 *      = 3      Matrix is singular                                   *
	 *      = 4      Matrix numerically singular                          *
	 *                                                                    *
	 *====================================================================*
	 *                                                                    *
	 *   Functions in use :                                               *
	 *   ==================                                               *
	 *                                                                    *

	 *      void vmfree():    free list of vectors and matrices           *
	 *                                                                    *
	 *====================================================================*
	 *====================================================================*/
	
	int  m, j, i, j0, ci;
	double piv, tmp, zmax,epsilon, cd;
	vector<double> d(n);
	epsilon = 2.2204460492503131e-016;
	
	if (n < 1) return (1);                   /*  Invalid parameters     */
	
	/*if (mat == NULL || lumat == NULL) return (1);
	if (perm == NULL) return (1);
	
	for (i = 0; i < n; i++)
    if (mat[i] == NULL || lumat[i] == NULL) return (1);*/

	/* d = Skaling vector for pivoting allocate storage   */
	//d[n];
	for(i=0;i<n;i++){ d[i]=0; };

	/* If  lumat and  mat are distinct, copy mat to lumat*/
	if (lumat != mat)                     
	{
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				lumat[i][j] = mat[i][j];
			}
		}
	}

	/* Initializations */
	for (i = 0; i < n; i++)
	{
		perm[i] = i;                        /* Initialize perm            */
		zmax = 0;
		for (j = 0; j < n; j++)             /* find row maxima            */
		{
			tmp = fabs(lumat[i][j]);
			if (tmp > zmax) zmax = tmp;
		}
		
		if (zmax == 0)                   /* mat is singular            */
		{
			//free(d);
			return (3);
		}
		d[i] = 1.0 / zmax;
	}

	*signd = 1;                         /* initialize sign of perm      */

	for (i = 0; i < n; i++)
	{
		piv = fabs(lumat[i][i]) * d[i];
		j0 = i;                           /* Search for pivot element     */
		for (j = i + 1; j < n; j++)
		{
			tmp = fabs(lumat[j][i]) * d[j];
			if (piv < tmp)
			{
				piv = tmp;                    /* Mark pivot element and       */
				j0 = j;                       /* its location                 */
			}
		}
		
		if (piv < epsilon)               /* If piv is small, mat is      */
		{                                 /* nearly singular              */
			*signd = 0;
			//delete(d);
			return (4);
		}
		
		if (j0 != i)
		{
			*signd = - *signd;              /* update signd  */
			/* Swap pivotentries   */
			ci = perm[j0];
			perm[j0]=perm[i];
			perm[i]=ci;
			/* Swap skaling vector   */
			cd=d[j0];
			d[j0]=d[i];
			d[i]=cd;
			/* Swap j0-th and i-th row of  lumat  */
			lumat[j0].swap(lumat[i]);
		}
		
		/* Gauss elimination    */
		for (j = i + 1; j < n; j++)       
		{
			if (lumat[j][i] != 0.0)
			{
				lumat[j][i] /= lumat[i][i];
				tmp = lumat[j][i];
				for (m = i + 1; m < n; m++)
					lumat[j][m] -= tmp * lumat[i][m];
			}
		}
	} /* end i */
  //delete(d);
  return (0);

};
int CNewtonIteration::NewtonCalcGaussSolution(int n,vector< vector<double> > &lumat,vector<int> &perm,vector<double> &b,vector<double> &x)
{
	/*====================================================================*
	 *                                                                    *
	 *  gausol  finds the solution x of the linear system  lumat * x = b  *
	 *  for the product matrix lumat, that describes an LU decomposition, *
	 *  as produced by gaudec.                                            *
	 *                                                                    *
	 *====================================================================*
	 *                                                                    *
	 *   Input parameters:                                                *
	 *   ================                                                 *
	 *      n        int n;  ( n > 0 )                                    *
	 *               Dimension of lumat,                                  *
	 *      lumat    REAL   *lumat[];                                     *
	 *               LU factorization, as produced from gaudec            *
	 *      perm     int perm[];                                          *
	 *               row permutation vector for lumat                     *
	 *      b        REAL   b[];                                          *
	 *               Right hand side of the system.                       *
	 *                                                                    *
	 *   Output parameter:                                                *
	 *   ================                                                 *
	 *      x        REAL   x[];                                          *
	 *               solution vector                                      *
	 *                                                                    *
	 *   Return value :                                                   *
	 *   =============                                                    *
	 *      = 0      all ok                                               *
	 *      = 1      n < 1 or other invalid input parameter               *
	 *      = 3      improper LU decomposition ( zero diagonal entry)     *
	 *                                                                    *
	 *====================================================================*/
	int  j, k;
	double sum;
	if (n < 1) return (1);                   /* Invalid input parameter */

	//if (lumat == NULL || b == NULL || perm == NULL) return (1);

	//for (j = 0; j < n; j++)
    //if (lumat[j] == NULL) return (1);

	for (k = 0; k < n; k++)                              /* update b    */
	{
		x[k] = b[perm[k]];
		for (j = 0; j < k; j++)
			x[k] -= lumat[k][j] * x[j];
	}

	for (k = n - 1; k >= 0; k--)                    /* back substitute  */
	{
		sum = 0;
		for (j = k + 1; j < n; j++)
			sum += lumat[k][j] * x[j];
		if (lumat[k][k] == 0) return (3);
		x[k] = (x[k] - sum) / lumat[k][k];
	}
	
	return (0);

};
//Calculate actual dF as vectorial norm
double CNewtonIteration::NewtonCalcL2norm(int n,vector<double> &dF)
{

/*==========================================================================================*
 *                                                                                          *
 *				CalcL2norm computes the euclidean or L2 norm of a vector                    *
 *				dF = (dF[0],dF[1],...,dF[n-1]). It avoides underflow.						*
 *	dFnorm = dF[0]*sqrt( 1.0 + dF[1]^2/dF[0]^2+ dF[2]^2/dF[0]^2 +,...,+ dF[n-1]^2/dF[0]^2 ) *
 *																							*
 *==========================================================================================*/

  int  i, j;
  double scale, sum, tmp, xiabs, minepsilon, minepsilon2;
  minepsilon = 2.2204460492503131e-016;
  minepsilon2 = minepsilon*minepsilon;

  if (n <= 0) return (0.0);          //n <= 0 ==> Norm = 0

  for (i = 0; i < n; i++)
  {
    if (dF[i] != 0.0) break;
  }

  if (i == n) return (0.0);         //zero vector

  scale = fabs( dF[i] );
  
  if (i == n - 1) return (scale);        //only one component  != 0


  j = i + 1;
  for (sum = 1.0, i = j; i < n; i++)
  {
    xiabs = fabs( dF[i] );
    if (xiabs <= scale)                 //scale = previous max      
    {                                   //of  ABS(x[i])             
      tmp = xiabs / scale;
      if (tmp > minepsilon2)
        sum += tmp * tmp;               //sum = sum + temp*temp     
    }
    else
    {
      tmp = scale / xiabs;
      if (tmp <= minepsilon2) tmp = 0.0;
      sum *= tmp * tmp;
      sum += 1.0;                     //sum = sum * temp * temp + 1
      scale = xiabs;
    }
  }

  return (scale * sqrt(sum));

};
//Determine piston force tolerance based on displacement chamber pressure force
void CNewtonIteration::NewtonCalcEpsilonK(void)
{
	//Force balance tolerance is set percentage of external force
	epsilonK = fabs(myPistonGap.forcespistongap.F_external[0]);
	for(int i=0;i<3;i++)
	{
		if( fabs(myPistonGap.forcespistongap.F_external[i+1]) > epsilonK )
		{
			epsilonK = fabs(myPistonGap.forcespistongap.F_external[i+1]);
		}
	}
	epsilonK *= myPistonGap.epsilonK;

};
//Solve Newton Iterative Loop
void CNewtonIteration::NewtonLoop(double epsilon,vector<double> &gappositions,vector<double> &gapvelocities)
{
	
	//Parameters definition
	int rc,kk,cas,signd,kn,ks;
	double dFnorm,dFnormmean,dFnormsigma;
	vector<double> dFnormtemp;

	//Number of states for Newton itration
	int n = (int) gapvelocities.size(); 
	
	//Vectors definition for Newton Iteration, Jacobian and Gauss methods
	vector<double> dF(n);						//dF vector out of the force calculation
	vector<double> tempdF(n);					//tempdF vector for Jacobian calculation
	vector<double> tempgapvelocities(n);		//temporary gap shifting velocities for Newton
	vector< vector <double> > Jacobian(n,n);	//Jacobian matrix
	vector<int> perm(n);						//perm vector for gauss solution
	vector<double> deltagapvelocities(n);		//delta gap velocities solved with NewtonCalcGauss	


	//First call to PistonGap function with no convergence loop
	myPistonGap.PistonGapJacobian(gappositions,gapvelocities,dF,dt);
	//Calculate dF
	dFnorm = NewtonCalcL2norm(n,dF);


	//---------------MAIN FORCE BALANCE LOOP - NEWTON ITERATION----------------//
	Log << "Piston/Cylinder Newton Iteration Loop..." << "\n";
	kk = 0;
	kn = 1;
	ks = 0;
	cas = 0;
	//----------External while loop----------//
	//fout.open("./output/piston/dF.txt",ios::app);
	do
	{
		//Delta Shifting Velocity Calculation
		myPistonGap.refinememory.resize(0);
		rc = NewtonCalcJacobian(n,gappositions,gapvelocities,Jacobian,dF,tempdF);
		rc = NewtonCalcGaussMain(cas,n,Jacobian,Jacobian,perm,dF,deltagapvelocities,&signd);
		//fout << phi << "\t" << gapvelocities[0] << "\t" << gapvelocities[1] << "\t" << gapvelocities[2] << "\t" << gapvelocities[3] << "\t" << dF[0] << "\t" << dF[1] << "\t" << dF[2] << "\t" << dF[3] << "\n";
		//limit magnitude of deltagapvelocities
		//double scale = 0.001;
		//scale /= sqrt(pow(deltagapvelocities[0],2.0)+pow(deltagapvelocities[1],2.0))+sqrt(pow(deltagapvelocities[2],2.0)+pow(deltagapvelocities[3],2.0));
		//scale *= 0.01;
		//scale = min(scale,1.0);
		//for(int i = 0;i<deltagapvelocities.size();i++)
			//deltagapvelocities[i] *= scale;


		//----------Internal while loop based on shifting velocities damping----------//
		double omega = 1.0;
		int jj = 0;
		double dFnormint;
		dFnormint=dFnorm;
		do
		{
			//Determine temporary gap velocities
			for (int i = 0; i < n; i++)
			{
				tempgapvelocities[i] = gapvelocities[i] - omega * deltagapvelocities[i];
			}
			//Newton relaxed loop
			myPistonGap.PistonGapJacobian(gappositions,tempgapvelocities,dF,dt);
			//fout << phi << "\t" << tempgapvelocities[0] << "\t" << tempgapvelocities[1] << "\t" << tempgapvelocities[2] << "\t" << tempgapvelocities[3] << "\t" << dF[0] << "\t" << dF[1] << "\t" << dF[2] << "\t" << dF[3] << "\n";
		
			//dF result vector
			dFnormint = NewtonCalcL2norm(n,dF);
			//Relaxation
			omega *= 0.5;
			//Internal counter
			jj++;
		}while(dFnormint > dFnorm && jj < jmax);
		Log << "j: " << jj << "\t" << "dFnorm_j: " << dFnormint << "\t scale: 1.0\n";

		//Assign values to main variables
		dFnorm = dFnormint;
		for (int i = 0; i < n; i++)                 
		{
			gapvelocities[i] = tempgapvelocities[i];
		}
		//External counter
		kk++;
		ks++;

		//----------Extra calculations to speed up convergence-----------//
		dFnormtemp.push_back(dFnorm);
		if(ks>=8)
		{
			dFnormsigma = 0.0;
			dFnormmean = 0.0;
			//Mean
			for(int i=0;i<(int) dFnormtemp.size();i++)
			{
				dFnormmean+=dFnormtemp[i];
			}
			dFnormmean /= (dFnormtemp.size()+1.0e-20);
			//Standard deviation
			for(int i=0;i<(int) dFnormtemp.size();i++)
			{
				dFnormsigma+=pow((dFnormtemp[i]-dFnormmean),2.0);
			}
			dFnormsigma /= (dFnormtemp.size()+1.0e-20);
			dFnormsigma = sqrt(dFnormsigma);
			if(dFnormsigma < 10.0)
			{
				dFnorm = 0.0;
			}
			dFnormtemp.clear();
			ks=0;
			Log << "dFnorm Mean: " << dFnormmean << "\t" << "dFnorm Sigma: " << dFnormsigma << "\n";
		}
		//Increase progressively tolerance to reach convergence
		if(kk>kn*8)
		{
			epsilon = epsilon + 10.0;
			kn++;
		}
		//Check for value very close to tolerance
		if( fabs(dFnorm-epsilon) < 3.0 )
		{
			dFnorm = 0.0; 
		};
		//---------------------------------------------------------------------//

	}while(dFnorm > epsilon && kk < kmax);
	Log << "k: " << kk << "\t" << "dFnorm_k: " << dFnorm << "\n";
	Log << "Done!" << "\n";
	//fout.close();
	//fout.clear();
	//Perform Sanity Check (added by dwm)
	if(!((dFnorm >= 0) && (dFnorm <= 1e9))){
		Log << "Error: unusual values for dFnorm detected.  Terminating Simulation." << endl;
		exit(1);
	}

	double maxecc = max(abs(gappositions[0]),abs(gappositions[1]));
	maxecc = max(maxecc,abs(gappositions[2]));
	maxecc = max(maxecc,abs(gappositions[3]));
	//cout<<maxecc<<"\n";
	if(maxecc>myinput.data.geometry.dZ){
		Log << "Error: piston eccentricity value has exceeded bushing diameter.  Terminating Simulation." << endl;
		Log << "gappositions[1] = " <<gappositions[0]<< endl;
		Log << "gappositions[2] = " <<gappositions[1]<< endl;
		Log << "gappositions[3] = " <<gappositions[2]<< endl;
		Log << "gappositions[4] = " <<gappositions[3]<< endl;
		exit(1);
	}

};
	