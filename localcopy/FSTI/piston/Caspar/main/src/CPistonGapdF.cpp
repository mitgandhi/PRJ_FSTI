#include "CPistonGap.h"
#include "..\..\caspar_input\input.h"
#pragma once

extern class CGapInput myGapInput;
extern struct sGapResult myGapResult;

extern class input myinput;
extern bool force_balance_iterative;

//CGapODE myGapODE; dwm

//Calculate fluid forces from fluid pressure
void CPistonGap::PistonCalcFluidForces(void)
{

	double FfKxtot, FfKytot, MfKxtot, MfKytot;
	
	//Initializations
	FfK = 0; FfKx = 0; FfKy = 0; MfKx = 0; MfKy = 0;
	FfKxtot = 0; FfKytot = 0; MfKxtot = 0; MfKytot = 0;
	
	//Fluid forces
	FfK = p * dAz ;
	
	//Fluid forces components
	FfKx = -1.0 * FfK * cos( phi ) ;
	FfKy = -1.0 * FfK * sin( phi ) ;

	//Fluid forces moments
	MfKx = -1.0 * FfKy * zKj ;  
	MfKy = FfKx * zKj ;

	//Total fluid force and moment
	FfKxtot = sum( FfKx ) ;
	FfKytot = sum( FfKy ) ;
	MfKxtot = sum( MfKx ) ;
	MfKytot = sum( MfKy ) ;

	forcespistongap.FfKx = FfKxtot ;
	forcespistongap.FfKy = FfKytot ;
	forcespistongap.MfKx = MfKxtot ;
	forcespistongap.MfKy = MfKytot ;
	forcespistongap.FfK = sqrt( (FfKxtot*FfKxtot) + (FfKytot*FfKytot) ) ;
	forcespistongap.MfK = sqrt( (MfKxtot*MfKxtot) + (MfKytot*MfKytot) ) ;

}
void CPistonGap::PistonCalcContactForces(void)
{
	double FcKxtot, FcKytot, McKxtot, McKytot, scale, sigmalimit;
	double hmin = geometrypistongap.hmin;
	int contact;

	//Calcuate contact stress from elastic compenetration
	if(force_balance_iterative){
		sigma += 0.1*Eprime*(hmin-h1)/rK;
		sigma = where( sigma < 0.0 , 0.0 , sigma ) ;
		sigma = where( sigma > 288e6 , 288e6 , sigma ) ;
	}
	else
		sigma = where( h1 < hmin, Eprime*(hmin-h1)/rK, 0.0);


	//Limit contact stress to first nl elements
	/*int nl = myinput.data.options_piston.numeric.penCells - 1; // (int) floor(dl/dy) ;
	Array<double,1> sigmadf = sigma;
	for(int i=0;i<(N*M);i++)
	{
		if( (i%M>nl) && (i%M<M-1-nl) )
		{
			sigmadf(i) = 0.0 ;
		}
	};*/

	//Saturate contact pressure to ~30% of pDC for stability
	//double pmax = 0.3 * operatingpistongap.pDC;
	//sigma = where( sigma>pmax , pmax , sigma ) ;

	//Scale contact pressure based on amount of surface in contact
	//contact = 0;
	//for(int i = 0;i<h1.size();i++)
	//	if(h1(i)<=hmin)
	//		contact++;
	//scale = contact / (h1.size() * 0.05);
	//sigmalimit = scale * max(sigma);
	//sigma = where(sigma>sigmalimit , sigmalimit , sigma);

	//Contact force field
	FcK = 0.0;
	FcK = sigma * dAz ;
	
	//Fluid forces components
	FcKx = 0.0; FcKy = 0.0;
	FcKx = -1.0 * FcK * cos( phi );
	FcKy = -1.0 * FcK * sin( phi );

	//Fluid forces moments
	McKx = 0.0; McKy = 0.0;
	McKx = -1.0 * FcKy * zKj;  
	McKy = FcKx * zKj;

	//Total fluid force and moment
	FcKxtot = 0.0; FcKytot = 0.0;
	McKxtot = 0.0; McKytot = 0.0;
	FcKxtot = sum( FcKx );
	FcKytot = sum( FcKy );
	McKxtot = sum( McKx );
	McKytot = sum( McKy );

	//Assign variables
	forcespistongap.FcKx = FcKxtot;
	forcespistongap.FcKy = FcKytot;
	forcespistongap.McKx = McKxtot;
	forcespistongap.McKy = McKytot;
	forcespistongap.FcK = sqrt( (FcKxtot*FcKxtot) + (FcKytot*FcKytot) );
	forcespistongap.McK = sqrt( (McKxtot*McKxtot) + (McKytot*McKytot) );

	//Contact forces in control points
	forcespistongap.F_contact[2] = forcespistongap.McKy / geometrypistongap.lvar;
	forcespistongap.F_contact[3] = -1.0 * forcespistongap.McKx / geometrypistongap.lvar;
	forcespistongap.F_contact[0] = -1.0 * forcespistongap.F_contact[2] + forcespistongap.FcKx;
	forcespistongap.F_contact[1] = -1.0 * forcespistongap.F_contact[3] + forcespistongap.FcKy;

};
//Calculate gap axial and circumferential friction forces
void CPistonGap::PistonCalcFrictionForces(void)
{
	Range all = Range::all();
	double omega = operatingpistongap.omega;
	double vK = operatingpistongap.vK;

	forcespistongap.FTKy.resize(N*M);
	forcespistongap.FTKy_p.resize(N*M);
	forcespistongap.FTKy_c.resize(N*M);
	forcespistongap.FTby.resize(N*M);
	forcespistongap.FTby_p.resize(N*M);
	forcespistongap.FTby_c.resize(N*M);
	forcespistongap.FTKx.resize(N*M);
	forcespistongap.FTKx_p.resize(N*M);
	forcespistongap.FTKx_c.resize(N*M);
	forcespistongap.FTbx.resize(N*M);
	forcespistongap.FTbx_p.resize(N*M);
	forcespistongap.FTbx_c.resize(N*M);

	//----------------------------------------------------------------//
	//---------------------------PISTON-------------------------------//
	//----------------------------------------------------------------//
	dvxT_p = 0.0; dvxT_c = 0.0; dvyT_p = 0.0; dvyT_c = 0.0; 
	//wall velocity gradient poiseuille - backward 3rd order
	dvxT_p = ( 11.0/6.0*vx_p(Range(N*M*(Q-1),N*M*Q-1)) - 3.0*vx_p(Range(N*M*(Q-2),N*M*(Q-1)-1)) +
		3.0/2.0*vx_p(Range(N*M*(Q-3),N*M*(Q-2)-1)) - 1.0/3.0*vx_p(Range(N*M*(Q-4),N*M*(Q-3)-1)) ) / dz2(Range(0,N*M-1)) ;
	dvyT_p = ( 11.0/6.0*vy_p(Range(N*M*(Q-1),N*M*Q-1)) - 3.0*vy_p(Range(N*M*(Q-2),N*M*(Q-1)-1)) +
		3.0/2.0*vy_p(Range(N*M*(Q-3),N*M*(Q-2)-1)) - 1.0/3.0*vy_p(Range(N*M*(Q-4),N*M*(Q-3)-1)) ) / dz2(Range(0,N*M-1)) ;
	//wall velocity gradient couette - backward 3rd order
	//Really? Not Sure about that...
	//dvxT_c = (omega*rK) / hT ;
	//dvyT_c = vK / hT ;
	//wall velocity gradient couette - backward 3rd order - Lizhi
	dvxT_c = ( 11.0/6.0*vx_c(Range(N*M*(Q-1),N*M*Q-1)) - 3.0*vx_c(Range(N*M*(Q-2),N*M*(Q-1)-1)) +
		3.0/2.0*vx_c(Range(N*M*(Q-3),N*M*(Q-2)-1)) - 1.0/3.0*vx_c(Range(N*M*(Q-4),N*M*(Q-3)-1)) ) / dz2(Range(0,N*M-1)) ;
	dvyT_c = ( 11.0/6.0*vy_c(Range(N*M*(Q-1),N*M*Q-1)) - 3.0*vy_c(Range(N*M*(Q-2),N*M*(Q-1)-1)) +
		3.0/2.0*vy_c(Range(N*M*(Q-3),N*M*(Q-2)-1)) - 1.0/3.0*vy_c(Range(N*M*(Q-4),N*M*(Q-3)-1)) ) / dz2(Range(0,N*M-1)) ;
	//Shear piston poiseuille
	taux = 0.0; tauy = 0.0;
	taux = muT * dvxT_p ;
	tauy = muT * dvyT_p ;
	//Force piston poiseuille
	forcespistongap.FTKx_p = -1.0 * taux * dAz ;
	forcespistongap.FTKy_p = -1.0 * tauy * dAz ;
	//Shear piston couette
	taux = 0.0; tauy = 0.0;
	taux = muT * dvxT_c ;
	tauy = muT * dvyT_c ;
	//Force piston couette
	forcespistongap.FTKx_c = -1.0 * taux * dAz ;
	forcespistongap.FTKy_c = -1.0 * tauy * dAz ;
	//Total force piston
	forcespistongap.FTKx = forcespistongap.FTKx_p + forcespistongap.FTKx_c;
	forcespistongap.FTKy = forcespistongap.FTKy_p + forcespistongap.FTKy_c;


	//----------------------------------------------------------------//
	//---------------------------CYLINDER-----------------------------//
	//----------------------------------------------------------------//
	dvxT_p = 0.0; dvxT_c = 0.0; dvyT_p = 0.0; dvyT_c = 0.0;
	//wall velocity gradient poiseuille - forward 3rd order
	dvxT_p = ( -11.0/6.0*vx_p(Range(0,N*M-1)) + 3.0*vx_p(Range(N*M,2*N*M-1)) - 
		3.0/2.0*vx_p(Range(2*N*M,3*N*M-1)) + 1.0/3.0*vx_p(Range(3*N*M,4*N*M-1))  ) / dz2(Range(0,N*M-1)) ;
	dvyT_p = ( -11.0/6.0*vy_p(Range(0,N*M-1)) + 3.0*vy_p(Range(N*M,2*N*M-1)) - 
		3.0/2.0*vy_p(Range(2*N*M,3*N*M-1)) + 1.0/3.0*vy_p(Range(3*N*M,4*N*M-1))  ) / dz2(Range(0,N*M-1)) ;
	//wall velocity gradient couette - forward 3rd order
	//dvxT_c = (omega*rK) / hT ;
	//dvyT_c = vK / hT ;
	//wall velocity gradient couette - forward 3rd order - Lizhi
	dvxT_c = ( -11.0/6.0*vx_c(Range(0,N*M-1)) + 3.0*vx_c(Range(N*M,2*N*M-1)) - 
		3.0/2.0*vx_c(Range(2*N*M,3*N*M-1)) + 1.0/3.0*vx_c(Range(3*N*M,4*N*M-1))  ) / dz2(Range(0,N*M-1)) ;
	dvyT_c = ( -11.0/6.0*vy_c(Range(0,N*M-1)) + 3.0*vy_c(Range(N*M,2*N*M-1)) - 
		3.0/2.0*vy_c(Range(2*N*M,3*N*M-1)) + 1.0/3.0*vy_c(Range(3*N*M,4*N*M-1))  ) / dz2(Range(0,N*M-1)) ;

	//Shear cylinder poiseuille
	taux = 0.0; tauy = 0.0;
	taux = muT * dvxT_p ;
	tauy = muT * dvyT_p ; 
	//Force cylinder poiseuille
	forcespistongap.FTbx_p = taux * dAz;
	forcespistongap.FTby_p = tauy * dAz;
	//Shear cylinder couette
	taux = 0.0; tauy = 0.0;
	taux = muT * dvxT_c ;
	tauy = muT * dvyT_c ;
	//Force cylinder couette
	forcespistongap.FTbx_c = taux * dAz;
	forcespistongap.FTby_c = tauy * dAz;
	//Total force cylinder
	forcespistongap.FTbx = forcespistongap.FTbx_p + forcespistongap.FTbx_c ;
	forcespistongap.FTby = forcespistongap.FTby_p + forcespistongap.FTby_c ;


	//---------------------------GENERAL-----------------------------//
	double beta = operatingpistongap.beta_rad;
	double phi = operatingpistongap.phi_rad;
	//Friction forces moments
	forcespistongap.MTK = sum(forcespistongap.FTKy) * tan(beta) * sin(phi) * rB;
	forcespistongap.MTKtan = sum(forcespistongap.FTKx) * rZ;
	//Mechanical power loss [W]
	PhiD_mech = 0.0;
	/*for(int i=0;i<Q;i++)
	{
		PhiD_mech += sum( muT * ( pow(dvxz(Range(i*N*M,(i+1)*N*M-1)),2.0) + pow(dvyz(Range(i*N*M,(i+1)*N*M-1)),2.0) ) * dx * dy * dz2(Range(0,N*M-1)) );
	};*/
	myGapResult.PhiD_mech_2d.resize(M*N);
	myGapResult.PhiD_mech_2d = 0.0;
	for(int i = 0;i<Q;i++){
		for(int j = 0;j<N;j++){
			for(int k = 0;k<M;k++){
				myGapResult.PhiD_mech_2d(k+M*j) += muT(k+M*j) * ( pow(dvxz(k+M*j+M*N*i),2.0) + pow(dvyz(k+M*j+M*N*i),2.0) ) * dx * dy * dz2(k+M*j+M*N*i);
			}
		}
	}
	vector<double> PhiD_sort(myGapResult.PhiD_mech_2d.size());
	for(int i = 0;i<myGapResult.PhiD_mech_2d.size();i++)
		PhiD_sort[i] = myGapResult.PhiD_mech_2d(i);
	sort(PhiD_sort.begin(),PhiD_sort.end());
	for(int i = PhiD_sort.size()-6;i<PhiD_sort.size();i++)
		PhiD_sort[i] = PhiD_sort[PhiD_sort.size()-7];
	for(int i = 0;i<PhiD_sort.size();i++)
		PhiD_mech += PhiD_sort[i];

};
//Guess slipper friction forces
void CPistonGap::PistonGuessFTG(void)
{
	double hG,oilviscosityguess, omega;
		
	hG = 5.0e-06;
	oilviscosityguess = 0.0;
	omega = operatingpistongap.omega;

	PistonGuessViscosity(oilviscosityguess);
	
	forcespistongap.FTG = oilviscosityguess * omega * rB / hG * AreaG ;

};
//Read slipper friction forces from slipper module.
void CPistonGap::PistonReadFTG(void)
{
	int searchindex = myGapInput.FTGFile.time.size() - 1;
	//cout << "Search Index: " << searchindex << "\n";
	double targettime = myGapResult.time - 2*PI / myinput.data.operating_conditions.speed;
	if(myGapInput.FTGFile.time.empty())
		PistonGuessFTG();
	else if(myGapInput.FTGFile.time[searchindex]<targettime)
		PistonGuessFTG();
	else{
		bool search = true;
		while(search && (searchindex >= 0)){
			searchindex += -1;
			if(myGapInput.FTGFile.time[searchindex]<targettime)
				search = false;
		};
		//cout << "Search Index: " << searchindex << "\n";
		if(searchindex<0){
			PistonGuessFTG();
			//cout << "Guessing FTG, " << targettime << "too small!" << "\n";
		}
		else{
			double y1 = myGapInput.FTGFile.FTG[searchindex];
			double y2 = myGapInput.FTGFile.FTG[searchindex +1];
			double t1 = myGapInput.FTGFile.time[searchindex];
			double t2 = myGapInput.FTGFile.time[searchindex + 1];
			forcespistongap.FTG = y1 + (y2 - y1) * (targettime - t1) / (t2 - t1) ;
		};
	};
};
//Calculate piston external forces
void CPistonGap::PistonCalcExternalForces(void)
{
	//Variable declaration
	double mK,lSK,beta,betamax,phi,omega,pDC,pCase,
		zRK,FDK,FaKz,FwK,FTKy,FTG,FTGx,FTGy,FAKz,FSK,FSKx,FSKy,FKx,FKy,MKx,MKy,gamma;

	//Variable initialization
	mK = geometrypistongap.mK;					//Mass piston/slipper assembly
	lSK = geometrypistongap.lSK;				//Distance piston center of mass/slipper assembly from piston head
	zRK = geometrypistongap.zRK;				//Distance between piston head and beggining of the gap (gap lenght included)
	beta = operatingpistongap.beta_rad;			//Swashplate angle
	betamax = operatingpistongap.betamax_rad;	//Max swashplate angle
	gamma = operatingpistongap.gamma_rad;		//Swashplate cross angle
	phi = operatingpistongap.phi_rad;			//Angular position
	omega = operatingpistongap.omega;			//Angular speed
	pDC = operatingpistongap.pDC;				//Pressure displacemnt chamber
	pCase = operatingpistongap.pCase;			//Case pressure
	FTG = forcespistongap.FTG;					//Slipper friction force

	//Pressure force
    FDK = AreaK * (pDC - pCase);
	
	//Center of mass inertia force z directions
    FaKz =  mK * (omega*omega) * rB * (tan(beta) * cos(phi) - tan(gamma) * sin(phi) / cos(beta));

	//Centifiugal force (set to zero to simulate EHD test rig)
    if(EHDTestRig)
	{
		FwK =  0.0;
	}
	else
	{
		FwK =  mK * (omega*omega) * rB;
	}

	//Friction force axial direction
	FTKy = sum(forcespistongap.FTKy);
    
	//Slipper friction force tangential direction  
	FTGx = - FTG;
	//Slipper friction force y direction
	FTGy = 0.0;		

	//Forces on the piston acting in z (axial) direction
	FAKz = FDK + FaKz + FTKy;

	//Force piston perpendicular to the swashplate
	//Anz = asin( pow( ( pow( sin(gamma) , 2.0 ) + pow( ( -1.0 * sin(beta) * cos(gamma) ) , 2.0 ) ) , 0.5 ) );
	//FSK = FAKz / cos(Anz);
	FSK = FAKz / (cos(beta)*cos(gamma));//check dwm
	forcespistongap.Fsk=FSK;
	//FSKxy = FAKz * tan(Anz);//pow( ( pow( FSK , 2.0 ) - pow( FAKz , 2.0 ) ) , 0.5 );

	//Transverse piston force x and y direction
	//Any = PI/2 - asin( sin(gamma) / pow( ( pow( sin(gamma) , 2.0 ) + pow( ( -1.0 * sin(beta) * cos(gamma) ) , 2.0 ) ) , 0.5 ) );
	//FSKx = FSKxy * sin(Any);
	//FSKy = FSKxy * cos(Any);
	FSKx = -1.0 * FSK * sin(gamma);
	FSKy = FSK * cos(gamma)*sin(beta);
	//FSKy = FAKz * tan(beta); //check dwm
	//FSKx = - FAKz * tan(gamma)/cos(beta); //check dwm
	
	//External force components x and y directions
	forcespistongap.FKx = -1.0 * FSKy * sin(phi) + FSKx * cos(phi) + FTGx;
	FKx = forcespistongap.FKx;
	forcespistongap.FKy = FSKy * cos(phi) + FSKx * sin(phi) + FwK + FTGy;
	FKy = forcespistongap.FKy;
	//Total external force
	forcespistongap.FK = sqrt( (FKx*FKx) + (FKy*FKy) );

	//External moments components x and y directions
	forcespistongap.MKx = - zRK * ( FSKy * cos(phi) + FSKx * sin(phi) ) - (zRK - lSK) * FwK;//dwm corrected
	MKx = forcespistongap.MKx;
	forcespistongap.MKy = zRK * FKx;
	MKy = forcespistongap.MKy;
	//Total external moment
	forcespistongap.MK = sqrt( (MKx*MKx) + (MKy*MKy) );
	
};

//Calculate force difference on control points
void CPistonGap::PistonCalcdF(vector<double> &dF)  
{
	int j;

	if(!force_balance_iterative)
		PistonCalcContactForces();

	//Fluid forces
	forcespistongap.F_fluid[2] = forcespistongap.MfKy / geometrypistongap.lvar;
	forcespistongap.F_fluid[3] = -1.0 * forcespistongap.MfKx / geometrypistongap.lvar;
	forcespistongap.F_fluid[0] = -1.0 * forcespistongap.F_fluid[2] + forcespistongap.FfKx;
	forcespistongap.F_fluid[1] = -1.0 * forcespistongap.F_fluid[3] + forcespistongap.FfKy;

	//External forces
	forcespistongap.F_external[2] = forcespistongap.MKy / geometrypistongap.lvar;
	forcespistongap.F_external[3] = -1.0 * forcespistongap.MKx / geometrypistongap.lvar;
	forcespistongap.F_external[0] = -1.0 * forcespistongap.F_external[2] + forcespistongap.FKx;
	forcespistongap.F_external[1] = -1.0 * forcespistongap.F_external[3] + forcespistongap.FKy;

	//dF calculation for force balance
	for(j=0;j<4;j++)
	{
		forcespistongap.dF[j] = - forcespistongap.F_external[j] - forcespistongap.F_fluid[j] - forcespistongap.F_contact[j];
	}
	//dF
	dF = forcespistongap.dF;

};