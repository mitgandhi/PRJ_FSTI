#include "CPistonGap.h"
#include <omp.h>
#include "../../caspar_input/input.h"
#include "logger.h"
#pragma once


extern class CGapInput myGapInput;
extern class CGapUtils myGapUtils;
extern class input myinput;

int nproc;

//Calculate piston surface pressure deformation based on IM
void CPistonGap::PistonFEMPressureSurfaceDeformation(void)
{	
	int nfK,nnK;//,phioffset;
	double pref,lKG,lch,lA,lvar,pDC,pCase,scale;//,phi_rad;//,beta,betamax;
	double beta = operatingpistongap.beta_rad;
	double betamax = operatingpistongap.betamax_rad;

	Range all = Range::all();
	//Element pressure
	Array<double,1> p_e,p_e1,p_e2;


	//Geometry
	lKG = geometrypistongap.lKG;
	lch = geometrypistongap.lch;
	lA = geometrypistongap.lA;
	//lvar = geometrypistongap.lvar;
	//cout << "lvar is: " << lvar << "\n";
	//Pressures
	pref = 100.0e5;
	pDC = operatingpistongap.pDC;
	pCase = operatingpistongap.pCase;
	//double psocket = forcespistongap.Fsk/(myinput.data.geometry.aSock);//[Pa]
	double psocket_0 = forcespistongap.Fsk*(1-(beta/betamax))/(myinput.data.geometry.aSock);//[Pa]
	double psocket_1 = forcespistongap.Fsk*(beta/betamax)/(myinput.data.geometry.aSock);
	//Saturate socket pressure to 2x DC pressure for stability
	/*if((psocket_0 + psocket_1) > (2 * pDC)){
		scale = (2 * pDC) / (psocket_0 + psocket_1);
		psocket_0 *= scale;
		psocket_1 *= scale;
		scale *= 100;
		Log << "Scaling FSK force by " << scale << "%" << "\n";
	}*/
	//cout << "lvar is: " << lvar << "\n";
	//Pressure interpolation from fluid to structure
	nfK = myGapInput.xyzfK_p.extent(0);
	p_e.resize(nfK);
	p_e1.resize(nfK);
	p_e2.resize(nfK);
	p_e = 0.0;
	p_e1 = 0.0;
	p_e2 = 0.0;


	//Interpolation
	Array<double,1> pdef(N*M);
	pdef = p;
	//cout << "lvar is: " << lvar << "\n";
	//Add pressure field from contact forces
	pdef += sigma;
	
	//Rotate piston pressure field as seen by piston
	/*phi_rad = operatingpistongap.phi_rad;
	phioffset = (int) (floor(phi_rad/dphi));
	phioffset = N - phioffset;
	if(phioffset<(N))
	{
		Array<double,1> pdef_copy;
		pdef_copy.resize(N*M);
		pdef_copy = pdef;
		pdef(Range(phioffset*M,N*M-1)) = pdef_copy(Range(0,(N-phioffset)*M-1));
		pdef(Range(0,phioffset*M-1)) = pdef_copy(Range((N-phioffset)*M,N*M-1));
	};*/

	p_e2 = myGapUtils.InterpolateFields(FaceIdK_f2s_p,FaceDistK_f2s_p,pdef);
	//cout << "lvar is: " << lvar << "\n";

	//Boundaries
	//cout << "lvar is: " << lvar << "\n";
	double tol = 0.5e-3;
	//cout << "lvar is: " << lvar << "\n";
	p_e1 = where( zfK_p>(lch+lA+geometrypistongap.lvar+tol), pCase , p_e2 );
	p_e = where( zfK_p<(geometrypistongap.lch+geometrypistongap.lA-tol), pDC , p_e1 );
	p_e2.free(); p_e1.free();
	//cout << "lvar is: " << lvar << "\n";
	


	//Deformation due to gap pressure
	defK_p = 0.0;
	nnK = myGapInput.xyznK_p.extent(0);
	int load;
	if(PressureDeformationOMP == 0){
		//initialize arrays and variables
		int flag;
		
		ULONGLONG idla, kera, usea, idlb, kerb, useb, syst, idle;
		idla = (((ULONGLONG) IdleTime.dwHighDateTime) << 32) + IdleTime.dwLowDateTime;
		kera = (((ULONGLONG) KernelTime.dwHighDateTime) << 32) + KernelTime.dwLowDateTime;
		usea = (((ULONGLONG) UserTime.dwHighDateTime) << 32) + UserTime.dwLowDateTime;
		//Recheck system times and take difference.
		GetSystemTimes(&IdleTime, &KernelTime, &UserTime);
		idlb = (((ULONGLONG) IdleTime.dwHighDateTime) << 32) + IdleTime.dwLowDateTime;
		kerb = (((ULONGLONG) KernelTime.dwHighDateTime) << 32) + KernelTime.dwLowDateTime;
		useb = (((ULONGLONG) UserTime.dwHighDateTime) << 32) + UserTime.dwLowDateTime;
		syst = (kerb+useb)-(kera+usea);
		idle = (idlb-idla);
		if(syst == 0){
			//Log << "Cannot Determine Number of Processors to Use." << "\n";
			nproc = 1;
		}
		else{
				
			load = 100 - (100 * idle/syst);
			if(load >= 65)
				nproc = 1;
			else if(load >= 60)
				nproc = 2;
			else if(load >= 55)
				nproc = 3;
			else
				nproc = 4;
			//Log << "Using " << nproc << " Processors! (load is " << load << "%.)" << "\n";
		}
	}
	else nproc = PressureDeformationOMP;
	omp_set_num_threads(nproc);
	IMParallel(1,load,nproc);

	//Parallelized
	if(nproc != 1)
	{
		//cout << "Running Calculations on " << nproc << " Processors" << "\n";
		double** defK_p_OMP = new double*[nproc];
		for(int i=0; i<nproc; i++)
		{
			defK_p_OMP[i] = new double[nnK];
			for(int j=0; j<nnK; j++)
			{
				defK_p_OMP[i][j] = 0.0;
			}
		}
		//parallel loop openMP
		# pragma omp parallel for
		for(int i=0;i<nfK;i++)
		{
			//get thread number
			int th = omp_get_thread_num( );
			//temporary deformation array for current thread
			//double* defK_p_temp = new double[nnK];
			//calculate thread gap pressure deformation
			double scale = p_e(i) / pref;
			/*for(int j=0;j<nnK;j++)
			{
				defK_p_temp[j] = scale * myGapInput.IM_piston(i*nnK+j);
			}*/
			//assign thread gap pressure deformation to proper array
			for(int j=0;j<nnK;j++)
			{
				//defK_p_OMP[th][j] += defK_p_temp[j];
				defK_p_OMP[th][j] += scale * myGapInput.IM_piston(i*nnK+j);
			}
			//delete temporary array
			//delete [] defK_p_temp; 
		}
		//sum threads arrays of gap deformation
		for(int i=0;i<nproc;i++)
		{
			for(int j=0;j<nnK;j++)
			{
				defK_p(j) += defK_p_OMP[i][j];
			}
		}
		//delete temporary deformation threads arrays
		for(int i=0; i<nproc; i++)
		{
			delete [] defK_p_OMP[i];
		}
		delete [] defK_p_OMP;
	}
	//Not parallelized
	else
	{
		//cout << "Running Calculations on 1 Processor" << "\n";
		for(int i=0;i<nfK;i++)
		{
			//calculate gap pressure deformation
			scale = p_e(i)/pref;
			defK_p += scale * myGapInput.IM_piston(Range(i*nnK,(i+1)*nnK-1));
		}
	}

	//Deformation due to displacement chamber pressure
	scale = pDC/pref;
	defK_p += scale * myGapInput.IM_piston(Range(nnK*nfK,nnK*nfK+nnK-1));

	//Deformation due to case pressure
	pref = 1.0e5;
	scale = pCase/pref;
	defK_p += scale * myGapInput.IM_piston(Range((nfK+1)*nnK,(nfK+1)*nnK+nnK-1));
	pref = 100.0e5;

	//Deformation due to socket pressure
	//scale = psocket/pref;
	//defK_p += scale * myGapInput.IM_piston(Range((nfK+2)*nnK,(nfK+2)*nnK+nnK-1));

	//Deformation due to socket_0 pressure
	scale = psocket_0/pref;
	defK_p += scale * myGapInput.IM_piston(Range((nfK+2)*nnK,(nfK+2)*nnK+nnK-1));

	//Deformation due to socket_1 pressure
	scale = psocket_1/pref;
	defK_p += scale * myGapInput.IM_piston(Range((nfK+3)*nnK,(nfK+3)*nnK+nnK-1));

};
//Calculate cylinder pressure deformation based on IM
void CPistonGap::PistonCylinderFEMPressureSurfaceDeformation(void)
{
	
	
	int nfB,nnB;
	double pref,scale,lF,pDC,pCase,lB,lvar,shiftangle;
	Range all = Range::all();
	vector<double> pressure;
	//Element pressure
	Array<double,1> p_e;
	
	//Geometry
	lF = geometrypistongap.lF;
	lB = geometrypistongap.lB;
	lvar = geometrypistongap.lvar;
	//Pressures
	pref = 100.0e5;
	//pDC = operatingpistongap.pDC;
	pCase = operatingpistongap.pCase;


	//Pressure interpolation from fluid to structure
	nfB = myGapInput.xyzfB_p.extent(0);
	p_e.resize(nfB);
	p_e = 0.0;


	//Interpolation
	Array<double,1> pdef(N*M);
	pdef = p;
	//add contact pressure
	pdef += sigma;
	p_e = myGapUtils.InterpolateFields(FaceIdB_f2s_p,FaceDistB_f2s_p,pdef);
	//Boundaries
	double tol = 0.5e-3;
	p_e = where( zfB_p<(lB-tol), pDC , p_e );
	p_e = where( zfB_p>(lB+lvar+tol), pCase , p_e );

	
	//Deformation due to gap pressure
	defB_p = 0.0;
	nnB = myGapInput.xyznB_p.extent(0);
	/*int load;
	int nproc;
	if(PressureDeformationOMP == 0){
		//initialize arrays and variables
		int flag;
		
		ULONGLONG idla, kera, usea, idlb, kerb, useb, syst, idle;
		idla = (((ULONGLONG) IdleTime.dwHighDateTime) << 32) + IdleTime.dwLowDateTime;
		kera = (((ULONGLONG) KernelTime.dwHighDateTime) << 32) + KernelTime.dwLowDateTime;
		usea = (((ULONGLONG) UserTime.dwHighDateTime) << 32) + UserTime.dwLowDateTime;
		//Recheck system times and take difference.
		GetSystemTimes(&IdleTime, &KernelTime, &UserTime);
		idlb = (((ULONGLONG) IdleTime.dwHighDateTime) << 32) + IdleTime.dwLowDateTime;
		kerb = (((ULONGLONG) KernelTime.dwHighDateTime) << 32) + KernelTime.dwLowDateTime;
		useb = (((ULONGLONG) UserTime.dwHighDateTime) << 32) + UserTime.dwLowDateTime;
		syst = (kerb+useb)-(kera+usea);
		idle = (idlb-idla);
		if(syst == 0){
			//Log << "Cannot Determine Number of Processors to Use." << "\n";
			nproc = 1;
		}
		else{
				
			nproc = (omp_get_num_procs()*((idle-0.5*syst)/syst))+1;
			if(nproc > omp_get_num_procs()/2)
				nproc = omp_get_num_procs()/2;
			if(nproc > 4)
				nproc = 4;
			if(nproc < 1)
				nproc = 1;
			load = 100 - (100 * idle/syst);

			//Log << "Using " << nproc << " Processors! (load is " << load << "%.)" << "\n";
		}
	}
	else nproc = PressureDeformationOMP;*/

	//IMParallel(1,load,nproc);

	//Parallelized
	if(nproc != 1)
	{
		//initialize arrays and variables
		double** defB_p_OMP = new double*[nproc];
		for(int i=0; i<nproc; i++)
		{
			defB_p_OMP[i] = new double[nnB];
			for(int j=0; j<nnB; j++)
			{
				defB_p_OMP[i][j] = 0.0;
			}
		}
		//parallel loop openMP
		# pragma omp parallel for
		for(int i=0;i<nfB;i++)
		{
			//get thread number
			int th = omp_get_thread_num( );
			//temporary deformation array for current thread
			//double* defB_p_temp = new double[nnB];
			//calculate thread gap pressure deformation
			double scale = p_e(i) / pref;
			/*for(int j=0;j<nnB;j++)
			{
				defB_p_temp[j] = scale * myGapInput.IM_cylinder(i*nnB+j);
			}*/
			//assign thread gap pressure deformation to proper array
			for(int j=0;j<nnB;j++)
			{
				//defB_p_OMP[th][j] += defB_p_temp[j];
				defB_p_OMP[th][j] += scale * myGapInput.IM_cylinder(i*nnB+j);
			}
			//delete temporary array
			//delete [] defB_p_temp; 
		}
		//sum threads arrays of gap deformation
		for(int i=0;i<nproc;i++)
		{
			for(int j=0;j<nnB;j++)
			{
				defB_p(j) += defB_p_OMP[i][j];
			}
		}
		//delete temporary deformation threads arrays
		for(int i=0; i<nproc; i++)
		{
			delete [] defB_p_OMP[i];
		}
		delete [] defB_p_OMP;
	}
	//Not parallelized
	else
	{
		for(int i=0;i<nfB;i++)
		{
			//calculate gap pressure deformation
			scale = p_e(i)/pref;
			defB_p += scale * myGapInput.IM_cylinder(Range(i*nnB,(i+1)*nnB-1));
		}
	}

	//Deformation due to other gaps pressure
	for(int i = 0;i<myinput.data.operating_conditions.npistons-1;i++){
		shiftangle = 360.0 * ((double)i+1.0)/(double)myinput.data.operating_conditions.npistons;
		if(ReadpFile)
		{ 
			pressure = PistonGetFilePressure( shiftangle );
		}
		else 
		{
			pressure = PistonGetIdealPressure( shiftangle );
		}
		scale = 0.5*(pressure[0]+pCase)/pref;
		defB_p += scale * myGapInput.IM_cylinder(Range(nnB*(nfB+i),nnB*(nfB+1+i)-1));
	}

	//Deformation due to displacement chamber pressure
	for(int i = 0;i<myinput.data.operating_conditions.npistons;i++){
		shiftangle = 360.0 * ((double)i)/(double)myinput.data.operating_conditions.npistons;
		if(ReadpFile)
		{ 
			pressure = PistonGetFilePressure( shiftangle );
		}
		else 
		{
			pressure = PistonGetIdealPressure( shiftangle );
		}
		scale = pressure[0]/pref;
		defB_p += scale * myGapInput.IM_cylinder(Range(nnB*(nfB-1+myinput.data.operating_conditions.npistons+i),nnB*(nfB+myinput.data.operating_conditions.npistons+i)-1));
	}

	//Deformation due to case pressure
	pref = 1.0e5;
	scale = pCase/pref;
	defB_p += scale * myGapInput.IM_cylinder(Range(nnB*(nfB+2*myinput.data.operating_conditions.npistons-1),nnB*(nfB+2*myinput.data.operating_conditions.npistons)-1));

};