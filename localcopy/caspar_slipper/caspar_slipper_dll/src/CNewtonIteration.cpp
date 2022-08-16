#include "CNewtonIteration.h"
#include <blitz/tinyvec.h>

extern void matlab(const Array<double,2>& data,const string file);

string gettime()
{
	time_t rawtime;
	struct tm * timeinfo;
	char buffer [80];
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	strftime (buffer,80,"%x %X",timeinfo);
	return buffer;
}
CNewtonIteration::CNewtonIteration(caspar_input * gapinputs) : gapinput(gapinputs), GapResult(gapinput)
{
	GapLog.message("NewtonIteration -> Constructing object NewtonIteration... ");
	
	if(gapinput->lubrication_module.solve_slipper == 1)
	{
		SlipperGap = new CSlipperGap(gapinput, GapResult);
		if(gapinput->options_slipper.general.SwashplatePressureDeformation == 1)
		{
			//We need to create a vector of all the displacement chamber pressures to init the swash deformation
			
			int steps = gapinput->common.rev_steps;
			vector<double> pDC(steps, 0);
			for(int i=0; i<steps; i++)
			{
				vector<double> p = SetPressure(gapinput->common.phistep_deg * double(i));
				pDC[i] = p[0];
			}

			SlipperGap->SwashplateInitPressure(pDC);
		}
		
	}

	GapLog.message("NewtonIteration -> Object NewtonIteration constructed successfully!");
};

CNewtonIteration::~CNewtonIteration(void)
{

	GapLog.message("NewtonIteration -> Object NewtonIteration destructed successfully!");

};
//MAIN FUNCTIONS FOR NEWTON LOOP AND SHIFTING VELOCITIES CALCULATION
int CNewtonIteration::NewtonCalcNewtonIteration(vector<double> &gappositions,vector<double> &gapvelocities,  vector<double> &pold,vector<double> &vold)
{
	int rc=0;

	vector<double> positionslipper(3);
	vector<double> velocityslipper(3);
	positionslipper = gappositions;
	velocityslipper = gapvelocities;

	
	//Calculate basic parameters
	double phi = gapinput->common.phi_deg; //current shaft angle

	GapResult.time = gapinput->common.time;						//Time
	GapResult.phi = phi;						//Shaft angle
	
	//Current pressures
	vector<double> pressure = SetPressure(gapinput->common.phi_rev_deg);
	gapinput->common.pDC = pressure[0];
	gapinput->common.pHP = pressure[1];
	gapinput->common.pLP = pressure[2];

	GapResult.pDC = pressure[0];	//Displacement chamber pressure [Pa]
	GapResult.pHP = pressure[1];	//High pressure port pressure [Pa]
	GapResult.pLP = pressure[2];	//Low pressure port pressure [Pa]

	ostringstream oss (ostringstream::out);
	
	//------------------SLIPPER SHIFTING VELOCITIES CALCULATION-------------------//
	if(gapinput->lubrication_module.solve_slipper == 1)
	{
		GapLog << "Solving for the Slipper Interface..."  << endl;
		
		int gapindex = 2;	//Gap index. Used with multiple interfaces
		int n = 3;			//Number of states for slipper positions
		vector<double> dFinit(3);

		oss.str("");
		oss << endl << gettime() << "\tPhi: " << gapinput->common.phi_deg << "\tpDC: " << gapinput->common.pDC/1.0e+5 << 
									"\tDegStep: " << gapinput->common.phistep_deg;
		GapLog.message(oss.str());

		//Call the slipper start timestep method
		SlipperGap->SlipperStartTimestep(positionslipper, velocityslipper, pold, vold);

		if(gapinput->options_slipper.general.preFSIforceBalance)
		{
			//just need some values for newton
			double dFnormInit = 100.0;
			double epsilonG = 50.0;

			GapLog << "\n\t****** Pre FSI Force Balance ******" << endl;
			if(gapinput->options_slipper.general.HybridForceBalance)
			{
				HybridCalcShiftingVelocities(gapindex,n,epsilonG,positionslipper,velocityslipper,dFnormInit);
			} else {
				NewtonCalcShiftingVelocities(gapindex,n,epsilonG,positionslipper,velocityslipper,dFnormInit);
			}
			GapLog << endl;

		}
		
		//I generally find that using 1/2 of the squeeze velocity in the FSI loop is better for convergence, especially at 
		//low film thicknesses
		vector<double> zero(3,0);
		vector<double> half(3,0);
		for (int i=0; i<3; i++)
		{
			half[i] = 0.5*velocityslipper[i];
		}

		if (gapinput->options_slipper.general.ComplexPicard == 0)
		{
			//SlipperGap->SlipperGap(positionslipper,zero,dFinit,1);
			//SlipperGap->SlipperGap(positionslipper,half,dFinit,1);
			if(gapinput->options_slipper.general.preFSIforceBalance)
			{
				SlipperGap->SlipperGap(positionslipper,velocityslipper,dFinit,1);
			} else {
				//it seems using 'half' velocity results in a "smoother" simulation @ low gap heights
				//when prenewtFSI isn't enabled
				SlipperGap->SlipperGap(positionslipper,half,dFinit,1);
			}
		} else {
			SlipperGap->SlipperGap(positionslipper,velocityslipper,dFinit,1);
		}
		
		double dFnormInit = NewtonCalcL2norm(n,dFinit);

		epsilonG = dFnormInit*gapinput->options_slipper.numeric.Newton1;

		if (epsilonG < gapinput->options_slipper.numeric.Newton3)
		{
			//epsilonG = dFnormInit*GapInput.Boundary.Slipper_Newton4;
			epsilonG = gapinput->options_slipper.numeric.Newton3;
		}
		if(epsilonG > gapinput->options_slipper.numeric.Newton2)
		{
			epsilonG = gapinput->options_slipper.numeric.Newton2;
		}


		oss.str("");
		oss << "\n\t****** Post FSI Force Balance ( dFnorm: " << dFnormInit << "\t epsG: " << epsilonG << " ) ******";
		GapLog.message(oss.str());

		//Calculate the piston shifting velocities
		if(gapinput->options_slipper.general.HybridForceBalance)
		{
			HybridCalcShiftingVelocities(gapindex,n,epsilonG,positionslipper,velocityslipper,dFnormInit);
		} else {
			NewtonCalcShiftingVelocities(gapindex,n,epsilonG,positionslipper,velocityslipper,dFnormInit);
		}
		
		// log message

		oss.str("");
		oss << endl << "\t****** ODE Integration ******" << endl;
		oss << "\th_r1 = " << positionslipper[0] << "\th_r2 = " << positionslipper[1] << "\th_r3 = " << positionslipper[2];
		GapLog.message(oss.str());
		oss.str("");
		oss << "\tdh1/dt = " << velocityslipper[0] << "\tdh2/dt = " << velocityslipper[1] << "\tdh3/dt = " << velocityslipper[2];
		GapLog.message(oss.str());
		oss.str("");
		oss << "\th_d1 = " << SlipperGap->Fluid.h(SlipperGap->Fluid.M-1,0) << "\th_d2 = " << SlipperGap->Fluid.h(SlipperGap->Fluid.M-1,SlipperGap->Fluid.N/3) << "\th_d3 = " << SlipperGap->Fluid.h(SlipperGap->Fluid.M-1,2*SlipperGap->Fluid.N/3);
		oss << endl << "\tMin h = " << min(where(SlipperGap->Fluid.hgroove==0,SlipperGap->Fluid.h,1)) << "\tMax h = " << max(SlipperGap->Fluid.h);
		GapLog.message(oss.str());

		if(gapinput->options_slipper.general.DenseMode == 2)
		{
			string junk;
			cout << "Press Enter to continue to Ctrl-C to quit...";
			cin >> junk;
			cout << endl;
			system("del /f/q .\\dense\\");
		}

		//Call the slipper end timestep method
		SlipperGap->SlipperEndTimestep();

	}
	

	gappositions = positionslipper;
	gapvelocities = velocityslipper;

	GapResult.gappositions = gappositions;		//parts positions vector
	GapResult.gapvelocities = gapvelocities;	//parts velocities vector

	return rc;
};
double CNewtonIteration::gsl_vector_norm(const gsl_vector * v)
{
	double sumsq = 0;

	for(size_t i=0; i<v->size; i++)
	{
		sumsq += pow(gsl_vector_get(v, i),2);
	}
	return pow(sumsq, 0.5);
}
int CNewtonIteration::HybridCallFunc(const gsl_vector * x, void *Params, gsl_vector * f)
{
	hybrid_params * params = (hybrid_params *) Params;
	
	vector<double> gapvelocities(params->n, 0);
	vector<double> dF(params->n, 0);

	for(size_t i=0; i<params->n; i++)
	{
		gapvelocities[i] = gsl_vector_get(x, i);
	}
	
	if(params->ComplexPicard >= 1)
	{
		params->newton->SlipperGap->SlipperGap(params->gappositions,gapvelocities,dF,2);
	} else {
		params->newton->SlipperGap->SlipperGap(params->gappositions,gapvelocities,dF,0);
	}

	for(size_t i=0; i<params->n; i++)
	{
		gsl_vector_set(f, i, dF[i]);
	}

	return GSL_SUCCESS;
}
int CNewtonIteration::HybridCalldFunc(const gsl_vector * x, void *params, gsl_matrix * J)
{
	//hopefully this isn't called a lot because it's a waste compared to dFunc

	hybrid_params * h_params = (hybrid_params *) params;
	gsl_vector *f = gsl_vector_alloc (h_params->n);
	
	HybridCallFuncdFunc(x,params,f,J);
	
	gsl_vector_free (f);

	return GSL_SUCCESS;
} 
int CNewtonIteration::HybridCallFuncdFunc(const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J)
{
	hybrid_params * h_params = (hybrid_params *) params;
	int cp = h_params->ComplexPicard;	//save the CP state
	h_params->ComplexPicard = 0;	//force the state to 0, b/c this is a Jacob calc

	HybridCallFunc (x, params, f);

	//Adaptave delta
	double mean = 0;
	for(size_t j=0; j<x->size; j++)
	{
		mean += fabs(gsl_vector_get(x, j));
	}
	mean /= double(x->size);
	
	double delta;

	if(mean < 1e-10)
	{
		delta = 1.0e-7;
	} else {
		delta = mean/20.0;
		//delta = mean/10.0;
	}

	if(delta > 1e-2)
	{
		delta = 1e-2;
	}

	
	gsl_multiroot_function func = {&HybridCallFunc, h_params->n, params};

	//Use the GSL jacob function
	gsl_multiroot_fdjacobian (&func, x, f, delta, J);

	/*
	{
		double df = 0;
		const double mindf = 50;
		const double maxdf = 2000;
		int cnt = 0;

		do
		{
			df = 0;

			//Use the GSL jacob function
			gsl_multiroot_fdjacobian (&func, x, f, delta, J);
    
			for(size_t i=0; i<x->size; i++)
			{
				df += fabs(gsl_matrix_get(J, i, i)*delta);
			}
			df /= double(x->size);

			//GapLog << "Jacob df = " << df << endl;

			if(df < mindf)
			{
				delta *= 1.25*mindf/df;
			}

			if(df > maxdf)
			{
				delta *= 0.8*maxdf/df;
			}

			cnt++;
		}
		while((df < mindf || df > maxdf) && cnt < 5);
	}
	*/

	/*
	{
		const double maxdf = 500;
		double df = maxdf*2.0;	//just to start up the while (should use a do-loop)
		while(df > maxdf)
		{
			df = 0;

			//Use the GSL jacob function
			gsl_multiroot_fdjacobian (&func, x, f, delta, J);
    
			for(size_t i=0; i<x->size; i++)
			{
				df += fabs(gsl_matrix_get(J, i, i)*delta);
			}
			df /= double(x->size);

			GapLog << "Jacob df = " << df << endl;

			if(df > maxdf)
			{
				delta *= 0.8*maxdf/df;
			}
		}
	}
	*/
	

	/*
	//Use my jacob code
	gsl_vector *x2 = gsl_vector_alloc (x->size);
	gsl_vector *f2 = gsl_vector_alloc (f->size);
	for(size_t i=0; i<x->size; i++)
	{
		for(size_t j=0; j<x->size; j++)
		{
			if(i == j)
			{
				gsl_vector_set(x2, j,  gsl_vector_get(x, j)+delta);
			} else {
				gsl_vector_set(x2, j,  gsl_vector_get(x, j));
			}
		}

		HybridCallFunc (x2, params, f2);

		for(size_t j=0; j<x->size; j++)
		{
			double df = gsl_vector_get(f2, j)-gsl_vector_get(f, j);
			gsl_matrix_set(J, j, i, df/delta);
		}

	}
	gsl_vector_free (x2);
	gsl_vector_free (f2);
	*/

	//restore the CP state
	h_params->ComplexPicard = cp;

	return GSL_SUCCESS;
}
int CNewtonIteration::HybridCalcShiftingVelocities(int gapindex,int n,double &epsilon,vector<double> &gappositions,vector<double> &gapvelocities, double & dFnorm)
{
	//const gsl_multiroot_fsolver_type *T;
	const gsl_multiroot_fdfsolver_type *T;
	//gsl_multiroot_fsolver *s;
	gsl_multiroot_fdfsolver *s;
     
	int status;
	size_t iter = 0;
     
	//const size_t n = 2;
	//struct rparams p = {1.0, 10.0};
	//void * params = this;
	   
	hybrid_params params;
	params.gappositions = gappositions;
	params.newton = this;
	params.n = n;
	if(gapinput->options_slipper.general.ComplexPicard >= 1)
	{
		params.ComplexPicard = 1;
	} else {
		params.ComplexPicard = 0;
	}

	//gsl_multiroot_function f = {&HybridCallFunc, n, &params};
	gsl_multiroot_function_fdf f = {&HybridCallFunc, &HybridCalldFunc, &HybridCallFuncdFunc, n, &params};
     
	gsl_vector *x = gsl_vector_alloc (n);
     
	for(size_t i=0; i<n; i++)
	{
		gsl_vector_set(x, i, gapvelocities[i]);
	}
    
	//T = gsl_multiroot_fsolver_hybrids;
	T = gsl_multiroot_fdfsolver_hybridsj;
	//s = gsl_multiroot_fsolver_alloc (T, n);
	s = gsl_multiroot_fdfsolver_alloc (T, n);
	//gsl_multiroot_fsolver_set(s, &f, x);
	gsl_multiroot_fdfsolver_set(s, &f, x);
     
	//print_state (iter, s);

	do
	{
		iter++;
		GapLog << "\tdFnorm = " << gsl_vector_norm(s->f);

		//status = gsl_multiroot_fsolver_iterate (s);
		status = gsl_multiroot_fdfsolver_iterate(s);
     
		//print_state (iter, s);
     
		if (status)   /* check if solver is stuck */
		{
			GapLog << "WARNING: Hybrid Force Balance Solver Stuck" << endl;
			break;
		}

		//status = gsl_multiroot_test_residual (s->f, epsilon);		
		double dFnorm = gsl_vector_norm(s->f);
		if(dFnorm > epsilon)
		{
			status = GSL_CONTINUE;
		} else {
			status = GSL_SUCCESS;
		}

		GapLog << " -> " << dFnorm << "\t( v0: " << gsl_vector_get(s->x, 0) << " v1: " <<  gsl_vector_get(s->x, 1) << " v2: " << gsl_vector_get(s->x, 2) << " )" << endl;
	}
	while (status == GSL_CONTINUE && iter < 30);
         
	gappositions = params.gappositions;
	for(size_t i=0; i<n; i++)
	{
		//since we can't limit velocities in gsl without writing the hybrid algorythm ourselves (and seeing if it makes sense)
		//we will just limit it on save
		
		double newval = gsl_vector_get(s->x, i);

		//limit crazy velocities
		const double vlimit = 5.0e-6/gapinput->common.timestep;
		if(newval > vlimit)
		{
			newval = vlimit;
		} else if (newval < -vlimit)
		{
			newval = -vlimit;
		}

		gapvelocities[i] = newval;
	}

	//gsl_multiroot_fsolver_free (s);
	gsl_multiroot_fdfsolver_free(s);
	gsl_vector_free (x);
       

	return 0;
}
int CNewtonIteration::NewtonCalcShiftingVelocities(int gapindex,int n,double &epsilon,vector<double> &gappositions,vector<double> &gapvelocities, double & dFnorm)
{
	//Parameters definition
	int rc,kk,jj,kmax,jmax,first,counter,cas,signd,kn,ks;
	double omega,dFnormint,dFnormold;
	vector<double>dFnormtemp;
	ostringstream oss (ostringstream::out);

	//Vectors definition for Newton Iteration, Jacobian and Gauss methods
	vector<double> dF(n);						//dF vector out of the force calculation
	vector<double> tempdF(n);					//tempdF vector for Jacobian calculation
	vector<double> tempgapvelocities(n);		//temporary gap shifting velocities for Newton
	vector< vector <double> > Jacobian(n);	//Jacobian matrix
	for(int i=0; i<n; i++)
	{
		Jacobian[i].resize(n);
	}
	vector<int> perm(n);						//perm vector for gauss solution
	vector<double> deltagapvelocities(n);		//delta gap velocities solved with NewtonCalcGauss

	vector<double> bestv(n);
	bestv = gapvelocities;
	double bestdFnorm = dFnorm;


	//if the socket is fixed, the three gap velocities need to be equal, preventing the slipper from
	//tipping. take the mean of the three for a first guess
	if(GapResult.SlipperResult.Sock_fixed)
	{
		double vz = gapvelocities[1]+gapvelocities[2]-gapvelocities[0];

		for(int i=0; i<n; i++)
		{
			gapvelocities[i] = vz;
		}
	}
		
	//-------------------------------------------------------------------------//
	//---------------MAIN FORCE BALANCE LOOP - NEWTON ITERATION----------------//
	//-------------------------------------------------------------------------//
	kk = 0;
	kn = 1;
	ks=0;
	kmax = 15;
		jj = 0;
		jmax = 8;

	first = 1;
	counter = first;
	//===================
	//External while loop
	//===================

	int maxjjcnt = 0;

	//the  jj!=jmax stops newton from going when it's in a local min
	// the kk==0 forces it to do at least one loop
	while((dFnorm > epsilon && kk < kmax && maxjjcnt < 15) || kk==0)	
	{
		//Jacobian Matrix Calculation
		if (counter < first)
		{
			counter = counter + 1;
			cas = 2;
		}
		else
		{
			counter = 0;
			cas = 0;
			rc = NewtonCalcJacobian(gapindex,n,gappositions,gapvelocities,Jacobian,dF,tempdF);
		}
		
		if(GapResult.SlipperResult.Sock_fixed)
		{
			double Jsum = 0;
			double dFsum = 0;
			for(int i=0; i<n; i++)
			{
				dFsum += dF[i];
				for(int j=0; j<n; j++)
				{
					Jsum += Jacobian[i][j];
				}
			}

			double deltav = dFsum/Jsum;
			for(int i=0; i<n; i++)
			{
				deltagapvelocities[i] = deltav;
			}
		} else {
			//Delta Shifting Velocity Calculation
			rc = NewtonCalcGaussMain(cas,n,Jacobian,Jacobian,perm,dF,deltagapvelocities,&signd);
		}

		//========================================================
		//Internal while loop based on shifting velocities damping
		//========================================================

		jj = 0;
		jmax = 4;
		dFnormint = dFnorm;
		dFnormold = dFnorm;
	
		oss.str("");
		if(GapResult.SlipperResult.Sock_fixed)
		{
			oss << "\t~dFnorm = " << dFnormint;
		} else {
			oss << "\tdFnorm = " << dFnormint;
		}
		
		omega = 1.0;
		{
			double maxgv =  fabs(gapvelocities[0]);
			double maxdgv =  fabs(deltagapvelocities[0]);
			for(int i=1; i<n; i++)
			{
				if(fabs(gapvelocities[i]) > maxgv)
				{
					maxgv = fabs(gapvelocities[i]);
				}
				if(fabs(deltagapvelocities[i]) > maxdgv)
				{
					maxdgv = fabs(deltagapvelocities[i]);
				}
			}
		}

		if(gapinput->options_slipper.general.ComplexPicard >= 1)
		{
			if(omega == 1.0)
			{
				GapLog.message("\t\tStep :: Picard Cnt: "); 
			} else {
				GapLog.message("\t\tStep (w=" + n2s(omega) + "):: Picard Cnt: "); 
			}
		}

		do
		{
			if(GapResult.SlipperResult.Sock_fixed)
			{
				double vz = gapvelocities[1]+gapvelocities[2]-gapvelocities[0];

				for(int i=0; i<n; i++)
				{
					gapvelocities[i] = vz;
				}
			}

			for (int i = 0; i < n; i++)
			{
				double newval = gapvelocities[i] - omega * deltagapvelocities[i];
				
				//limit crazy velocities
				const double vlimit = 5.0e-6/gapinput->common.timestep;
				if(newval > vlimit)
				{
					newval = vlimit;
				} else if (newval < -vlimit)
				{
					newval = -vlimit;
				}
				
				tempgapvelocities[i] = newval;
			}

			switch(gapindex)
			{
				case 1:
					break;
				case 2:
					if(gapinput->options_slipper.general.ComplexPicard >= 1)
					{
						SlipperGap->SlipperGap(gappositions,tempgapvelocities,dF,2);
					} else {
						SlipperGap->SlipperGap(gappositions,tempgapvelocities,dF,0);
					}
					break;
				case 3: 
					break;
			}

			if(GapResult.SlipperResult.Sock_fixed)
			{
				dFnormint = 0;
				for(int i=0; i<n; i++)
				{
					dFnormint += dF[i];		
				}
				dFnormint = fabs(dFnormint);
			} else {
				dFnormint = NewtonCalcL2norm(n,dF);
			}

			omega *= 0.5;
			jj++;
			//Internal while loop counter
		}while(dFnormint > dFnorm && jj < jmax);

		//oss << " -> " << dFnormint << " in " << jj << " iterations";
		oss << " -> " << dFnormint << " in " << jj << " iterations.\t( v0: " << tempgapvelocities[0] << " v1: " <<  tempgapvelocities[1] << " v2: " << tempgapvelocities[2] << " )";

		if(gapinput->options_slipper.general.DenseMode == 2)
		{

			matlab(SlipperGap->Fluid.p,"./dense/p" + n2s(kk) + ".txt");
			matlab(SlipperGap->Fluid.h,"./dense/h" + n2s(kk) + ".txt");
		}

		switch(gapindex)
		{
			case 1:
				break;
			case 2:
				if(gapinput->options_slipper.general.ComplexPicard >= 1)
				{
					GapLog.message(""); 
				}
				GapLog.message(oss.str());
				break;
			case 3: 
				break;
		}

		//=======================
		//Internal while loop end
		//=======================
		//Assign the current values to the main variables
		double dFimprove = dFnorm - dFnormint;
		dFnorm = dFnormint;
		oss.str("");
		oss << dFnormint;
		for (int i = 0; i < n; i++)                 
		{
			gapvelocities[i] = tempgapvelocities[i];
			oss << "\t" << gapvelocities[i];
		}
		
		//Log the best velocities
		if(dFnorm < bestdFnorm)
		{
			bestv = gapvelocities;
			bestdFnorm = dFnorm;
		}

		if(jj == jmax)
		{
			maxjjcnt++;
		}

		//External while loop counter
		kk++;
		ks++;



		//Increase progressively tolerance to reach convergence
		if(dFimprove < 1.0)
		{
			if(kk>kn*2)
			{
				if (epsilon*1.1 > epsilon+1.0)
				{
					epsilon *= 1.1;
				} else {
					epsilon += 1.0;
				}
				kn++;
			}
		}
		
	}	//outer while loop end

	if(dFnorm > epsilon)	
	{
		//newton didn't converge so just use the best velocities found
		gapvelocities = bestv;

		//update the slipper gap to this new V
		if(gapinput->options_slipper.general.ComplexPicard >= 1)
		{
			SlipperGap->SlipperGap(gappositions,gapvelocities,dF,2);
		} else {
			SlipperGap->SlipperGap(gappositions,gapvelocities,dF,0);
		}
	}
	
	if(dFnorm > epsilon && kk >= kmax)
	{
		//Newton failed to converge
		return 1;
	} else {
		return 0;
	}

};
int CNewtonIteration::NewtonCalcJacobian(int gapindex,int n,vector<double> &gappositions,vector<double> &gapvelocities,vector< vector<double> > &Jacobian,vector<double> &dF,vector<double> &tempdF)
{
	//Jacobian matrix calculation: each term of the matrix expressed with central difference as (dF - tempdF)/(2*h)
	//where dF and tempdF are evaluated at (gapvelocities + h) and (gapvelocities - h)
	double  xj;
	vector<double> delta(n);
	vector<double> denom(n);

	
	
	//Delta of speed

	//Constant
	//delta =  1.0e-6;
	
	
	//Adaptave delta
	double mean = 0;
	for(int j=0; j<n; j++)
	{
		mean += fabs(gapvelocities[j]);
	}
	mean /= double(n);
	
	/*
	if(mean < 1e-10)
	{
		delta = 1.0e-7;
	} else {
		delta = mean/30.0;
		//delta = mean/75.0;
	}
	*/

	for(int j=0; j<n; j++)
	{
		if(mean < 1e-10)
		{
			delta[j] = 1.0e-7;
		} else {
			delta[j] = mean/20.0;
			//delta = mean/10.0;
		}
		denom[j] = 1.0/(2.0*delta[j]);
	}

	

	if(gapinput->options_slipper.general.ComplexPicard == 2)
	{
		GapLog.message("\t\tJacobian (" + n2s(delta[0]) + ", " + n2s(delta[1]) + ", " + n2s(delta[2]) + ") :: Picard Cnt: ");
	}

	for (int j = 0; j < n; j++)
	{
		xj = gapvelocities[j];
		gapvelocities[j] += delta[j];      
		
		//Switch to calculate the correct gap based on the value of gapindex
		//Vector tempdF is returned
		switch(gapindex)
		{
		case 1:
		
			break;
		case 2:
			if(gapinput->options_slipper.general.ComplexPicard == 2)
			{
				SlipperGap->SlipperGap(gappositions,gapvelocities,tempdF,3);
			} else {
				SlipperGap->SlipperGap(gappositions,gapvelocities,tempdF,0);
			}
			break;
		case 3: 
		
			break;
		}
		
		gapvelocities[j] = xj;
		gapvelocities[j] -= delta[j];
		
		//Switch to calculate the correct gap based on the value of gapindex
		//Vector dF is returned
		switch(gapindex)
		{
		case 1:
			break;
		case 2:
			if(gapinput->options_slipper.general.ComplexPicard == 2)
			{
				SlipperGap->SlipperGap(gappositions,gapvelocities,dF,3);
			} else {
				SlipperGap->SlipperGap(gappositions,gapvelocities,dF,0);
			}
			break;
		case 3: 
			break;
		}
		
		gapvelocities[j] = xj;

		//Jacobian calculation
		for (int i = 0; i < n; i++)
		{
			Jacobian[i][j] = (tempdF[i] - dF[i]) * denom[j];
		}
	
	}

	if(gapinput->options_slipper.general.ComplexPicard == 2)
	{
		GapLog.message("");
	}

	return 0;

};
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
vector<double> CNewtonIteration::SetPressure(const double phi_deg)
{
	int readpressure = 1;	//User interface option
	vector<double> pressure;
	if(readpressure != 1)
	{ 
		pressure =  GetIdealPressure(phi_deg);
	}
	else 
	{
		pressure = GetFilePressure(phi_deg);
	}

	return pressure;
};
vector<double> CNewtonIteration::GetIdealPressure(double phi_deg)
{
	//the requested phi / time
	phi_deg = fmod(phi_deg, 360.0);

	vector<double> pressure;
	pressure.resize(3);
	
	const double pHP = gapinput->operating_conditions.HP;
	const double pLP = gapinput->operating_conditions.LP;

    if(phi_deg < 3.0)
	{         
		pressure[0] = pLP + (pHP - pLP)/3.0*phi_deg;      
	}
    else if(phi_deg < 179.0 && phi_deg >=3.0)
	{        
		pressure[0] = pHP;
	}
    else if(phi_deg < 182.0 && phi_deg >= 179.0)
	{
		pressure[0] = pHP - (pHP-pLP)/3.0*(phi_deg - 179.00 );
	}
    else
	{
		pressure[0] = pLP;      
	}
	pressure[1] = pHP;
	pressure[2] = pLP;

    return pressure;

};
vector<double> CNewtonIteration::GetFilePressure(double phi_deg)
{
	caspar_input::Pfile * pfile(&gapinput->pfile);	//a pointer just to shorten the local code

	//the requested phi / time
	phi_deg = fmod(phi_deg, 360.0);
	const double time = phi_deg/360.0 * gapinput->common.rev_period;
	
	double pTime = 0;	//the current pressure file time
	double pAngle = 0;	//the current pressure file phi (deg)
	int size = (int) pfile->time.size();
	vector<double> pressure;
	for(int i=0; i < size ; i++)
	{
		pTime = pfile->time[i];
		pAngle = pTime / gapinput->common.rev_period * 360.0;
		pAngle = fmod(pAngle, 360.0);
		
		if(phi_deg <= pAngle)
		{
			//adding interpolation
			if(i>0)
			{
				double slope = (pfile->pDC[i] - pfile->pDC[i-1]) / (pfile->time[i] - pfile->time[i-1]);
				pressure.push_back((time - pfile->time[i-1])*slope + pfile->pDC[i-1]);

				slope = (pfile->pHP[i] - pfile->pHP[i-1]) / (pfile->time[i] - pfile->time[i-1]);
				pressure.push_back((time - pfile->time[i-1])*slope + pfile->pHP[i-1]);
				
				slope = (pfile->pLP[i] - pfile->pLP[i-1]) / (pfile->time[i] - pfile->time[i-1]);
				pressure.push_back((time - pfile->time[i-1])*slope + pfile->pLP[i-1]);
			} else {
				pressure.push_back(pfile->pDC[i]);
				pressure.push_back(pfile->pHP[i]);
				pressure.push_back(pfile->pLP[i]);
			}
			
			return pressure;
		}
	}
	//In case simulated angle is bigger then last time in file get last pressure
	pressure.push_back(pfile->pDC[size-1]);
	pressure.push_back(pfile->pHP[size-1]);
	pressure.push_back(pfile->pLP[size-1]);

	return pressure;
		
};
int CNewtonIteration::sgn(const double x)
{
	//A simple sign function
	if(x > 0)
	{
		return 1;
	} 
	else if (x < 0)
	{
		return -1;
	} else {
		return 0;
	}
}
