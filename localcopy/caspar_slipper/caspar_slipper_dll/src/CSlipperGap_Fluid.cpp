#include "CSlipperGap.h"
#include <omp.h>
#include <time.h>
#pragma once
#define OMP

extern void matlab(const Array<double,2>& data,const string file);

void CSlipperGap::SlipperCalcViscosity(const double alpha)
{
	/*
	using namespace blitz::tensor; 

	int M = Fluid.M; int N = Fluid.N; int Q = Fluid.Q;
	double oilPc1,oilPc2,oilTc1,oilTc2,oilbetaP,oilbetaT,oilW;
	
	//Old mu for under relax
	Array<double,3> oldmu(M,N,Q);
	oldmu = Fluid.oilviscosity;

	double mu_0 = 1.702e-2;
	double A0 = 2.278e-8;
	double A1 = 1.477e-10;
	double A2 = 5.127e-13;
	double B = 5.428;
	double C = 94.72;
	Fluid.oilviscosity = 0;

	oilPc1 = oilslippergap.oilPc1;
	oilPc2 = oilslippergap.oilPc2;
	oilTc1 = oilslippergap.oilTc1;
	oilTc2 = oilslippergap.oilTc2;
	oilbetaP = oilslippergap.oilbetaP;
	oilbetaT = oilslippergap.oilbetaT;
	oilW = oilslippergap.oilW;

	//Limit pressure and temperature
	Array<double,2> p(where(Fluid.p(i,j)>10000e+5,10000e+5,Fluid.p(i,j)));
	p = where(p(i,j)<0.1e+5,0.1e+5,Fluid.p(i,j));

	Array<double,3> T(where(Fluid.T(i,j,k)>100,100,Fluid.T(i,j,k)));
	T = where(T(i,j,k)<1,1,T(i,j,k));

	if (oilslippergap.oiltype == 0)
	{
		Fluid.oilviscosity = oilslippergap.oilviscosity;
	}
	if (oilslippergap.oiltype == 1)
	{
		Array<double,3> a(M,N,Q);
		Array<double,3> b(M,N,Q);
		Array<double,3> mu_kin(M,N,Q);

		a = oilPc1 * p(i,j) * pow( (T(i,j,k)+273.15),oilPc2 );        //change for viscosity second option 10.11.07 Daniel
		b = pow( 10.0, oilTc1 - oilTc2*log10(T(i,j,k)+273.15) );
		mu_kin = oilW * exp( a(i,j,k) ) * ( pow(10.0,b(i,j,k)) - 0.5 );
		Fluid.oilviscosity = 1.0e-6 * mu_kin(i,j,k) * ( oilslippergap.oildensity*(1 + oilbetaP * p(i,j) - oilbetaT * (T(i,j,k)- 20)) );
	}
	if (oilslippergap.oiltype == 2)
	{
		Array<double,3> alpha(M,N,Q);
		Array<double,3> Temp_coeff(M,N,Q);
		Array<double,3> Temp_coeff2(M,N,Q);

		alpha = A0 - A1*T + A2*( pow(T,2.0) );
		Temp_coeff =  B*( 50 - T(i,j,k) )/( C + T(i,j,k) ) + alpha(i,j,k) * p(i,j);
		
		//Limit Temp_coeff to 4
		//Temp_coeff = where(Temp_coeff(i,j,k)>4,4,Temp_coeff(i,j,k));
		
		Temp_coeff2 = exp(Temp_coeff);
		Fluid.oilviscosity = mu_0 * Temp_coeff2;
	}
	if (oilslippergap.oiltype == 3)
	{
		//Skydrol - model developed using AMEsim & SAE specs by Marco Zecchi
		double pref = 1.00000000000000e005;
		double Tref = 1.00000000000000e001;
		double rhoref = 1.06690035585407e003;
		double muref = 4.78435154853258e-002;
		double cpref = 1.67717744042304e003;
		double lambdaref = 1.36331614988375e-001;

		double P1 = 1.02059515692251e-008;
		double T1 = -1.72947118956637e-002;
		double T2 = 6.53539084574889e-005;
		double PT = -1.40663906488398e-011;
		double PT2 = -3.12090395313133e-013;
		double T3 = -9.56687972187690e-008;
		double PT3 = 1.01110563823919e-015;

		Array<double,2> dp(p - pref);
		Array<double,3> dT(T - Tref);
		
		Array<double,3> Psi ( 
			P1*dp(i,j) + T1*dT(i,j,k) + T2*pow(dT(i,j,k), 2) + PT*dp(i,j)*dT(i,j,k) + PT2*dp(i,j)*pow(T1,2) 
			+ T3*pow(dT(i,j,k), 3) + PT3*dp(i,j)*pow(T1, 3)
									);

		Fluid.oilviscosity = muref*pow(10.0, Psi);
	}

	//Underrelax
	Fluid.oilviscosity = oldmu - alpha*(oldmu-Fluid.oilviscosity);

	Fluid.oilviscosity = where(Fluid.oilviscosity>1.0,1.0,Fluid.oilviscosity);
	*/

	for(int i=0; i<Fluid.M; i++)
	{
		for(int j=0; j<Fluid.N; j++)
		{
			for(int k=0; k<Fluid.Q; k++)
			{
				Fluid.oilviscosity(i,j,k) = oil_properties->get_mu(Fluid.p(i,j), Fluid.T(i,j,k));
			}
		}
	}
}
void CSlipperGap::SlipperCalcDensity(void)
{
	/*
	using namespace blitz::tensor;

	double oilRhozero = oilslippergap.oildensity; //[kg/m3]
	double oilbetaT = oilslippergap.oilbetaT; //[1/K]
	double oilbetaP = oilslippergap.oilbetaP; //[1/Pa]
	Array<double,3> oilRhoT(Fluid.M,Fluid.N,Fluid.Q);

	double rs = 1047.03; // [kg/m3]
	double als = 5.761668e-4; // [1/K]
	double a1 = 0.0732965;
	double a2 = 1965.02; // [bar]
	double a3 = -2.96813; // [bar/K]
	Fluid.oildensity = 0;

	//Limit pressure and temperature
	Array<double,2> p(where(Fluid.p(i,j)>10000e+5,10000e+5,Fluid.p(i,j)));
	p = where(p(i,j)<0.1e+5,0.1e+5,Fluid.p(i,j));

	Array<double,3> T(where(Fluid.T(i,j,k)>100,100,Fluid.T(i,j,k)));
	T = where(T(i,j,k)<1,1,T(i,j,k));

	if (oilslippergap.oiltype == 0)
	{
		Fluid.oildensity = oilRhozero;
	}
	if (oilslippergap.oiltype == 1)
	{
		Fluid.oildensity = oilRhozero * (1 + oilbetaP * p(i,j) - oilbetaT * (T(i,j,k) - 20));
	}
	if (oilslippergap.oiltype == 2)
	{
		oilRhoT = rs*(1 - als*(Fluid.T + 273.15)); 
		Fluid.oildensity = oilRhoT / ( 1-a1*log( (a2+a3*(T(i,j,k) + 273.15) + 1e-5 * p(i,j) )/(a2+a3*(T(i,j,k) + 273.15)) ) );
	}
	if (oilslippergap.oiltype == 3)
	{
		//Skydrol - model developed using AMEsim & SAE specs by Marco Zecchi
		double pref = 1.00000000000000e005;
		double Tref = 1.00000000000000e001;
		double rhoref = 1.06690035585407e003;
		double muref = 4.78435154853258e-002;
		double cpref = 1.67717744042304e003;
		double lambdaref = 1.36331614988375e-001;

		double T1 = 7.80474256715724e-004;
		double T2 = 8.33673780452984e-007;
		double P1 = -8.76030682143087e-010;
		double P2 = 3.76207238974272e-018;
		double PT = -5.44707978002713e-012;
		double PT2 = 1.17343267760062e-014;
		double P2T = 2.78987659880380e-020;
		double P2T2 = -7.01924728092232e-023;
		
		Array<double,2> dp(p - pref);
		Array<double,3> dT(T - Tref);

		Array<double,3> vs ( 
										(1.0/rhoref)*
										( 
											1 + P1*dp(i,j) + P2*pow(dp(i,j), 2.0) + T1*dT(i,j,k) + T2*pow(dT(i,j,k), 2.0) +  
											PT*dp(i,j)*dT(i,j,k) + PT2*dp(i,j)*pow(dT(i,j,k),2) + P2T*dT(i,j,k)*pow(dp(i,j),2) +
											P2T2*pow(dp(i,j), 2)*pow(dT(i,j,k), 2)
										)
									);
			
		Fluid.oildensity = 1.0/vs;
	}
	*/
	for(int i=0; i<Fluid.M; i++)
	{
		for(int j=0; j<Fluid.N; j++)
		{
			for(int k=0; k<Fluid.Q; k++)
			{
				Fluid.oildensity(i,j,k) = oil_properties->get_rho(Fluid.p(i,j), Fluid.T(i,j,k));
			}
		}
	}
};
double CSlipperGap::CalcK(const double& T_C,const double& p_Pa)
{
	/*
	double RhoT =0.0, RhoTp=0.0, K=0.0;
	double T_K = T_C + 273.15;
	//Oiltype option 0: User defined constant bulk modulus
	if (oilslippergap.oiltype == 0)
	{
		K = oilslippergap.oilK; // [Pa]
	}
	//Oiltype option 1: User defined linear formula bulk modulus
	if (oilslippergap.oiltype == 1)
	{
		double Kzero = oilslippergap.oilK; //[Pa]
		double betaKT = oilslippergap.oilbetaKT; //[1/K]
		double betaKP = oilslippergap.oilbetaKP; //[1/Pa]
		K= Kzero*(1 + betaKP*(p_Pa-1e5) - betaKT*(T_C - 20)); // Rho(T,p)
	}
	//Oiltype option 2: Oil is HLP 32
	if (oilslippergap.oiltype == 2) //Oil is HLP 32
	{
		double rs = 1047.03; // [kg/m3]
		double als = 5.761668e-4; // [1/K]
		double a1 = 0.0732965;
		double a2 = 1965.02; // [bar]
		double a3 = -2.96813; // [bar/K]
		RhoT = rs*(1-als*T_K); 
		RhoTp= RhoT/(1-a1*log((a2+a3*T_K+p_Pa*1.0e-5)/(a2+a3*T_K)));
		K = 1.0e5*(RhoT*(a2+a3*T_K+p_Pa*1.0e-5)/(a1*RhoTp)); //[Pa]
	}
	if (oilslippergap.oiltype == 3) //Oil is Skydrol
	{
		//Skydrol - model developed using AMEsim & SAE specs by Marco Zecchi
		//Using definition of Bulk Modulus = rho / (drho/dp)
		//drho/dp is the analytical expression from the rho Skydrol model
		
		double pref = 1.00000000000000e005;
		double Tref = 1.00000000000000e001;
		double rhoref = 1.06690035585407e003;
		double muref = 4.78435154853258e-002;
		double cpref = 1.67717744042304e003;
		double lambdaref = 1.36331614988375e-001;

		double T1 = 7.80474256715724e-004;
		double T2 = 8.33673780452984e-007;
		double P1 = -8.76030682143087e-010;
		double P2 = 3.76207238974272e-018;
		double PT = -5.44707978002713e-012;
		double PT2 = 1.17343267760062e-014;
		double P2T = 2.78987659880380e-020;
		double P2T2 = -7.01924728092232e-023;
		
		double dp = p_Pa - pref;
		double dT = T_C - Tref;

		double drhodp =
		-(rhoref*(P1 + 2*P2*dp + PT*dT + PT2*dT*dT + 2*P2T*dT*dp + 2*P2T2*dT*dT*dp))
			/
		pow(P2T2*dT*dT*dp*dp + PT2*dT*dT*dp + T2*dT*dT + P2T*dT*dp*dp + PT*dT*dp + T1*dT + P2*dp*dp + P1*dp + 1,2.0);

		K = Fluid.oildensity(0,0,0)/drhodp;
	}

	return K;
	*/
	return oil_properties->get_K(p_Pa, T_C);
	
}
void CSlipperGap::SlipperCalcFluidV(void)
{
	//An update
	Fluid.z = Fluid.h(tensor::i,tensor::j)/(Fluid.Q-1)*tensor::k;

	int M = Fluid.M; int N = Fluid.N; int Q = Fluid.Q;
	Array<double,2> dpr(M,N);
	Array<double,2> dptheta(M,N);
		
	Array<double,2> pn(M,N);		pn = 0;			//Actual pressure field 2D vector
	Array<double,2> ps(M,N);		ps = 0;			//Actual pressure field 2D vector
	Array<double,2> pe(M,N);		pe = 0;			//Actual pressure field 2D vector
	Array<double,2> pw(M,N);		pw = 0;			//Actual pressure field 2D vector
	Array<double,2> mu(M,N);		mu = 0;			//Actual pressure field 2D vector
	
	//Calculate dpr
	dpr(Range(1,M-2),all) = (Fluid.p(Range(2,M-1),all)-Fluid.p(Range(0,M-3),all))/(0.5*Fluid.dr(Range(0,M-3),all)+0.5*Fluid.dr(Range(2,M-1),all)+Fluid.dr(Range(1,M-2),all));
	dpr(0,all) = (Fluid.p(1,all)-Fluid.p(0,all))/(0.5*Fluid.dr(0,all)+0.5*Fluid.dr(1,all));
	dpr(M-1,all) = (Fluid.p(M-1,all)-Fluid.p(M-2,all))/(0.5*Fluid.dr(M-1,all)+0.5*Fluid.dr(M-2,all));

	//Calculate dptheta
	pe(all,Range(0,N-2)) = Fluid.p(all,Range(1,N-1));
	pe(all,N-1) = Fluid.p(all,0);
	pw(all,Range(1,N-1)) = Fluid.p(all,Range(0,N-2));
	pw(all,0) = Fluid.p(all,N-1);
	dptheta = (pw - pe)/(2*Fluid.dtheta);

	//mu = mean(Fluid.oilviscosity(tensor::i,tensor::j,tensor::k),tensor::k);

	/*
	for(int i=0; i<Fluid.M; i++)
	{
		for(int j=0; j<Fluid.N; j++)
		{
			Fluid.vr(i,j,all) = .5/mu(i,j)*dpr(i,j)*((Fluid.z(i,j,all)*Fluid.z(i,j,all))-Fluid.h(i,j)*Fluid.z(i,j,all));
			//Fluid.vr(i,j,all) = .5/mu(i,j)*dpr(i,j)*((Fluid.z(i,j,all)*Fluid.z(i,j,all))-Fluid.h(i,j)*Fluid.z(i,j,all))+operatingslippergap.vgr(i,j)*Fluid.z(i,j,all)/Fluid.h(i,j);
			//Fluid.vtheta = .5/(mu(i,j)*Fluid.r(i,j))*dptheta(i,j)*((Fluid.z*Fluid.z)-Fluid.h(i,j)*Fluid.z)+operatingslippergap.vgtheta(i,j)*Fluid.z/Fluid.h(i,j);	
		}
	}
	*/

	Fluid.vr = .5/Fluid.oilviscosity*dpr(tensor::i,tensor::j)*((Fluid.z*Fluid.z)-Fluid.h(tensor::i,tensor::j)*Fluid.z)+operatingslippergap.vgr(tensor::i,tensor::j)*Fluid.z/Fluid.h(tensor::i,tensor::j);
	Fluid.vtheta = .5/(Fluid.oilviscosity*Fluid.r(tensor::i,tensor::j))*dptheta(tensor::i,tensor::j)*((Fluid.z*Fluid.z)-Fluid.h(tensor::i,tensor::j)*Fluid.z)+operatingslippergap.vgtheta(tensor::i,tensor::j)*Fluid.z/Fluid.h(tensor::i,tensor::j);	

	Array<double,3> dz(Fluid.M,Fluid.N,Fluid.Q);
	dz(all,all,Range(1,Fluid.Q-2)) = 0.5*(Fluid.z(all,all,Range(2,Fluid.Q-1))-Fluid.z(all,all,Range(0,Fluid.Q-3)));
	dz(all,all,0) = 0.5*(Fluid.z(all,all,1)-Fluid.z(all,all,0));
	dz(all,all,Fluid.Q-1) = 0.5*(Fluid.z(all,all,Fluid.Q-1)-Fluid.z(all,all,Fluid.Q-2));
	
	//old methods
	//operatingslippergap.QSG = sum(sum(Fluid.vr(tensor::i,tensor::j,tensor::k)*dz(tensor::i,tensor::j,tensor::k),tensor::k)*Fluid.r(Fluid.M-1,all)*Fluid.dtheta);
	//operatingslippergap.QSG = sum(sum(Fluid.vr(Fluid.M-1,all,all)*dz(Fluid.M-1,all,all),tensor::j)*Fluid.r(Fluid.M-1,all)*Fluid.dtheta);
	
	//considering flux at the 'outter' ring only
	//operatingslippergap.QSG = sum(mean(Fluid.vr(Fluid.M-1,all,all),tensor::j)*Fluid.h(Fluid.M-1,all))*geometryslippergap.routG*Fluid.dtheta;

	//considering flux at the 'inner' (#2) ring only
	//operatingslippergap.QSG = sum(mean(Fluid.vr(2,all,all),tensor::j)*Fluid.h(2,all)*Fluid.r(2,all)*Fluid.dtheta(2,all));

	if(gapinput->options_slipper.fluid_grid.SealingLand < Fluid.M-1 && gapinput->options_slipper.fluid_grid.SealingLand > 0)
	{
		const int m = gapinput->options_slipper.fluid_grid.SealingLand;
		//operatingslippergap.QSG = sum(mean(Fluid.vr(m,all,all),tensor::j)*Fluid.h(m,all)*Fluid.r(m,all)*Fluid.dtheta(m,all));
		//operatingslippergap.QSG = sum(sum(Fluid.vr(m,all,all)*dz(m,all,all), tensor::j)*Fluid.r(m,all)*Fluid.dtheta(m,all));
		Array<double, 1> tmp(sum((-1.0/(12.0*mean(Fluid.oilviscosity,tensor::k))*dpr*pow(Fluid.h,3.0)-0.5*operatingslippergap.vgr*Fluid.h)*Fluid.r*Fluid.dtheta, tensor::j));
		operatingslippergap.QSG = tmp(m);

		//Also output the individual leakage components - Poiseuille and Couette
		Array<double, 1> tmp_pois(sum((-1.0/(12.0*mean(Fluid.oilviscosity,tensor::k))*dpr*pow(Fluid.h,3.0))*Fluid.r*Fluid.dtheta, tensor::j));  //Poiseuille
		Array<double, 1> tmp_couette(sum((-0.5*operatingslippergap.vgr*Fluid.h)*Fluid.r*Fluid.dtheta, tensor::j));  //Couette
		operatingslippergap.Q_S_pois = tmp_pois(m);
		operatingslippergap.Q_S_couette = tmp_couette(m);
		

	} else {
		GapLog.message("ERROR: Slipper Grid SealingLand should be greater than 0, less than M-1, and not on a boudary.");
		exit(1);
	}
	
	

	//considering the mean flux throughout the ring
	//operatingslippergap.QSG = mean(sum(mean(Fluid.vr(tensor::i,tensor::j,tensor::k),tensor::k)*Fluid.h(tensor::i,tensor::j)*Fluid.r(tensor::i,tensor::j)*Fluid.dtheta,tensor::j));
}
void CSlipperGap::SlipperInitializeTemperature(void)
{
	double T_DC,T_Leak,phi;

	//----------Initialize gap temperature based-----------//
	//Fixed surface temperature profile for slipper and bushing
	phi = gapinput->common.phi_rev_deg;
	if(phi < 180)
	{
		T_DC = temperatureslippergap.T_HP;
	}
	else
	{
		T_DC = temperatureslippergap.T_LP;
	}
	T_Leak = temperatureslippergap.T_Leak;

	/*
	//linear method

	//Update temperature boundary values
	Fluid.T = where(Fluid.boundary(tensor::i,tensor::j) == -2, T_DC, Fluid.T);
	Fluid.T = where(Fluid.boundary(tensor::i,tensor::j) == -3, T_Leak, Fluid.T);

	
	//The first and last row MUST be a boundary or we have problems
	for(int i=0; i<Fluid.M-1;)
	{
		int j=i+1;

		//Find the next boundary value
		while(Fluid.boundary(j,0) == -1)
		{
			j++;
			//We should put an error message here if j == Fluid.M
		}
		
		//Initialize the 2D pressure using linear interpolation between the boundaries	
		for(int ii=i+1; ii<j; ii++)
		{
			//Blitz++ doesn't seem to let me write the expression I need to
			//Doing it the loop way then
			for(int n=0; n<Fluid.N; n++)
			{
				for(int q=0; q<Fluid.Q; q++)
				{
					if(Fluid.boundary(ii,n) == -1)
					{
						Fluid.T(ii,n,q) = (Fluid.T(j,n,q)-Fluid.T(i,n,q))/(Fluid.r(j,n)-Fluid.r(i,n))*(Fluid.r(ii,n)-Fluid.r(i,n))+Fluid.T(i,n,q);	
					}
				}
			}
		}

		i=j;
	}

	//Set the swashplate boundary T to just T_Leak
	Fluid.T(Range::all(),Range::all(),0) = T_Leak;
	*/

	//avg method
	Fluid.T = 0.5 * (T_DC + T_Leak);

	//Update temperature boundary values
	Fluid.T = where(Fluid.boundary(tensor::i,tensor::j) == -2, T_DC, Fluid.T);
	Fluid.T = where(Fluid.boundary(tensor::i,tensor::j) == -3, T_Leak, Fluid.T);

	//Set the swashplate boundary T to just T_Leak
	Fluid.T(Range::all(),Range::all(),0) = T_Leak;

	//update thermal analysis boundries
	if(gapinput->options_slipper.general.SlipperThermoElastic)
	{
		Fluid.T(Range::all(), Range::all(), Fluid.Q-1) = t_slipper.temp;
	}
	if(gapinput->options_slipper.general.SwashplateThermoElastic)
	{
		Fluid.T(Range::all(), Range::all(), 0) = t_swashplate.temp;
	}
}
void CSlipperGap::SlipperpG(const double alpha) 
{
	double pGold = operatingslippergap.pG;
	double QSG = operatingslippergap.QSG;

	/*
	//Does this make sense??
	double pockFlow = sum(0.5*Fluid.r(0,all)*Fluid.r(0,all)*Fluid.dtheta(0,all)*Fluid.dht(0,all));
	double expFlow = sum(0.5*Fluid.r(Fluid.M-1,all)*Fluid.r(Fluid.M-1,all)*Fluid.dtheta(Fluid.M-1,all)*Fluid.dht(Fluid.M-1,all));
	QSG += expFlow;
	*/

	double pG_conv;

	/*
	//conventional method - lumped steady-state conservation of volume
	if(geometryslippergap.flowtype == 0)
	{
		//option 0
		//turbulent flow through two orifices - piston and slipper
		double AG = PI * (geometryslippergap.dDG*geometryslippergap.dDG) / 4;
		double AK = PI * (geometryslippergap.dDK*geometryslippergap.dDK) / 4;
		double a1 = 0.5*Fluid.oildensity(0,0,0)*QSG*QSG;
		double a2 = (AG*AG+AK*AK)/(geometryslippergap.alphaD_KG*geometryslippergap.alphaD_KG*AG*AG*AK*AK);
		
		if(QSG < 0)
		{
			operatingslippergap.pG = operatingslippergap.pDC+a1*a2;
		} else {
			operatingslippergap.pG = operatingslippergap.pDC-a1*a2;
		}

		//pG_conv = operatingslippergap.pDC-a1*a2;
		
	} else if(geometryslippergap.flowtype == 1)
	{
		//option 1
		//laminar flow through two orifices - piston and slipper
		double RG = 128*Fluid.oilviscosity(0,0,0)*geometryslippergap.lDG/(PI * pow(geometryslippergap.dDG,4.0));
		double RK = 128*Fluid.oilviscosity(0,0,0)*geometryslippergap.lDK/(PI * pow(geometryslippergap.dDK,4.0));
		double a1 = (RG+RK)*QSG;
		
		operatingslippergap.pG = operatingslippergap.pDC-a1;
	}
	*/


	//internal report method

	//pressure buildup approach
	double pocket_volume = geometryslippergap.vPocket; //m^3, should be input in geometry.txt
	pocket_volume += sum(0.5*Fluid.h(0, all)*pow(Fluid.r(0, all), 2.0)*Fluid.dtheta(0, all));
	double dvdt = sum(0.5*Fluid.dht(0, all)*pow(Fluid.r(0, all), 2.0)*Fluid.dtheta(0, all));
	double K = CalcK(Fluid.T(0,0,Fluid.Q-1), operatingslippergap.pG);

	if(geometryslippergap.flowtype == 0)
	{
		double A,B,C;

		//option 0
		//turbulent flow through two orifices - piston and slipper
		double AG = PI * (geometryslippergap.dDG*geometryslippergap.dDG) / 4;
		double AK = PI * (geometryslippergap.dDK*geometryslippergap.dDK) / 4;

		B = geometryslippergap.alphaD_KG*AG*AK*pow(2.0/(Fluid.oildensity(0,0,0)*(AG*AG+AK*AK)), 0.5);

		A = operatingslippergap.pGprevious - gapinput->common.timestep*K/pocket_volume*(QSG+dvdt);
		B *= gapinput->common.timestep*K/pocket_volume;
		C = operatingslippergap.pDC;
		
		if(A > C)
		{
			double sqrt_v = B*B-4*C+4*A;
			if(sqrt_v < 0)
			{
				GapLog.message("Slipper pocket pressure Error (1). Sqrt of less than 0. QSG = " + n2s(QSG));
				operatingslippergap.pG = 0.1e+5;
			} else {
				operatingslippergap.pG = A-0.5*B*(-B+pow(sqrt_v, 0.5));
			}
		} else {
			double sqrt_v = B*B+4*C-4*A;
			if(sqrt_v < 0)
			{
				GapLog.message("Slipper pocket pressure Error (2). Sqrt of less than 0. QSG = " + n2s(QSG));
				operatingslippergap.pG = 0.1e+5;
			} else {
				operatingslippergap.pG = A+0.5*B*(-B+pow(sqrt_v, 0.5));
			}
		}

	} else if(geometryslippergap.flowtype == 1)
	{
		double A,B;

		double rG4 = pow(0.5*geometryslippergap.dDG,4.0);
		double lDG = geometryslippergap.lDG;
		double rK4 = pow(0.5*geometryslippergap.dDK,4.0);
		double lDK = geometryslippergap.lDK;
		
		B = PI/(8.0*Fluid.oilviscosity(0,0,0))*(rG4*rK4/(rK4*lDG+rG4*lDK));
		A = operatingslippergap.pGprevious - gapinput->common.timestep*K/pocket_volume*(QSG+dvdt);
		B *= gapinput->common.timestep*K/pocket_volume;

		operatingslippergap.pG = (A+B*operatingslippergap.pDC)/(1+B);

	}

	
	//end of simpler method
	
	//rate limit to n bar change
	const double dp_limit = 5.0e5;
	if(operatingslippergap.pG > pGold + dp_limit)
	{
		operatingslippergap.pG = pGold + dp_limit;
	}

	if(operatingslippergap.pG < pGold - dp_limit)
	{
		operatingslippergap.pG = pGold - dp_limit;
	}

	if(operatingslippergap.pG < 0.1e5) //limit pressure
	{
		operatingslippergap.pG = 0.1e5;
	}

	if(operatingslippergap.pG > 1000.0e5) //limit pressure
	{
		operatingslippergap.pG = 1000.0e5;
	}

	//underrelax the new pG pressure
	operatingslippergap.pG = pGold + alpha*(operatingslippergap.pG-pGold); 

	if(operatingslippergap.pG < 0.1e5) //limit pressure
	{
		operatingslippergap.pG = 0.1e5;
	}

	if(operatingslippergap.pG > 1000.0e5) //limit pressure
	{
		operatingslippergap.pG = 1000.0e5;
	}

	
	//find the leakage through the slipper orifice using this pG
	//this is what we will call 'QSG'
	if(geometryslippergap.flowtype == 0)
	{
		//option 0
		//turbulent flow through two orifices - piston and slipper
		double AG = PI * (geometryslippergap.dDG*geometryslippergap.dDG) / 4;
		double AK = PI * (geometryslippergap.dDK*geometryslippergap.dDK) / 4;

		double dp = operatingslippergap.pDC-operatingslippergap.pG;
		double sgn = 1;
		if(dp < 0)
		{
			sgn = -1;
			dp = fabs(dp);
		}

		operatingslippergap.QSG_orifice = sgn*geometryslippergap.alphaD_KG*AG*AK*pow(2.0*dp/(Fluid.oildensity(0,0,0)*(AG*AG+AK*AK)), 0.5);
	} else if(geometryslippergap.flowtype == 1)
	{
		//option 1
		//laminar flow through two throttles - piston and slipper

		double rG4 = pow(0.5*geometryslippergap.dDG,4.0);
		double lDG = geometryslippergap.lDG;
		double rK4 = pow(0.5*geometryslippergap.dDK,4.0);
		double lDK = geometryslippergap.lDK;

		double dp = operatingslippergap.pDC-operatingslippergap.pG;
		
		operatingslippergap.QSG_orifice = dp*PI/(8.0*Fluid.oilviscosity(0,0,0))*(rG4*rK4/(rK4*lDG+rG4*lDK));
	}
	

	//Update pressure boundary values
	Fluid.p = where(Fluid.boundary >= 0, Fluid.boundary, Fluid.p);
	Fluid.p_uncut = where(Fluid.boundary >= 0, Fluid.boundary, Fluid.p_uncut);
	Fluid.p = where(Fluid.boundary == -2, operatingslippergap.pG, Fluid.p);
	Fluid.p_uncut = where(Fluid.boundary == -2, operatingslippergap.pG, Fluid.p_uncut);
	Fluid.p = where(Fluid.boundary == -3, operatingslippergap.pcase, Fluid.p);
	Fluid.p_uncut = where(Fluid.boundary == -3, operatingslippergap.pcase, Fluid.p_uncut);
	
}
void CSlipperGap::SlipperDefinePressureBounds(void)
{
	Fluid.boundary = -1;								//Default to solve Reynolds
	Fluid.boundary(0,all) = -2;					//Set pocket pressure boundary
	Fluid.boundary(Fluid.M-1,all) = -3;			//Set case pressure boundary
	Fluid.hgroove = 0;

	//Set non-uniform radial grid pressure bounds if applicable
	if(Fluid.Groove1Location >= 0)
	{
		if(Fluid.Groove1Location <= gapinput->options_slipper.fluid_grid.SealingLand)
		{
			//pocket pressure
			Fluid.boundary(Range(Fluid.Groove1Location-1,Fluid.Groove1Location+1),all) = -2;
			Fluid.hgroove(Range(Fluid.Groove1Location-1,Fluid.Groove1Location+1),all) = 1;
		} else {
			//case pressure
			Fluid.boundary(Range(Fluid.Groove1Location-1,Fluid.Groove1Location+1),all) = -3;
			Fluid.hgroove(Range(Fluid.Groove1Location-1,Fluid.Groove1Location+1),all) = 1;
		}
	}

	if(Fluid.Groove2Location >= 0)
	{
		if(Fluid.Groove2Location <= gapinput->options_slipper.fluid_grid.SealingLand)
		{
			//pocket pressure
			Fluid.boundary(Range(Fluid.Groove2Location-1,Fluid.Groove2Location+1),all) = -2;
			Fluid.hgroove(Range(Fluid.Groove2Location-1,Fluid.Groove2Location+1),all) = 1;
		} else {
			//case pressure
			Fluid.boundary(Range(Fluid.Groove2Location-1,Fluid.Groove2Location+1),all) = -3;
			Fluid.hgroove(Range(Fluid.Groove2Location-1,Fluid.Groove2Location+1),all) = 1;
		}
		
	}

	//We can only begin to consider nonuniform theta bounds when we have a slipper IM set
	if(gapinput->options_slipper.general.SlipperPressureDeformation == 1)
	{
		//how many pts should be tried for interpolation boundedness
		int maxBoundTry = 10;
		ANNpoint qPt = annAllocPt(2);

		for(int m=0; m<Fluid.M; m++)
		{
			for(int n=0; n<Fluid.N; n++)
			{
				if(Fluid.boundary(m,n) != -1)
				{
					//We don't need to test for a boundary on a cell that is already a boundary
					continue;
				}
				qPt[0] = Fluid.LRx(m,n);
				qPt[1] = Fluid.LRy(m,n);


				//Interpolation points
				for(int k = 1; k<=maxBoundTry+1; k++)
				{
					if(k == maxBoundTry+1)
					{
						//This is a step pressure boundary

						if(m < gapinput->options_slipper.fluid_grid.SealingLand)
						{
							//Pocket pressure
							Fluid.boundary(m,n) = -2;
							Fluid.hgroove(m,n) = 1;
						} else {
							//Case pressure
							Fluid.boundary(m,n) = -3;
							Fluid.hgroove(m,n) = 1;
						}			
					}

					ANNidxArray	  nnIdx = new ANNidx[k];
					ANNdistArray  dists = new ANNdist[k];

					slipper.KDfaces.kdtree->annkSearch
					(
						qPt,
						k,
						nnIdx,
						dists,
						slipper.KDfaces.eps
					);	


					//Barycentric Interpolation
					//http://en.wikipedia.org/wiki/Barycentric_coordinates_(mathematics)

					const double x = qPt[0];
					const double y = qPt[1];
					const double x1 = slipper.nodes[slipper.faces[nnIdx[k-1]].nodes[0]].x;
					const double y1 = slipper.nodes[slipper.faces[nnIdx[k-1]].nodes[0]].y;
					const double x2 = slipper.nodes[slipper.faces[nnIdx[k-1]].nodes[1]].x;
					const double y2 = slipper.nodes[slipper.faces[nnIdx[k-1]].nodes[1]].y;
					const double x3 = slipper.nodes[slipper.faces[nnIdx[k-1]].nodes[2]].x;
					const double y3 = slipper.nodes[slipper.faces[nnIdx[k-1]].nodes[2]].y;

					const double detT = (y2-y3)*(x1-x3)+(x3-x2)*(y1-y3);
					const double lam1 = ((y2-y3)*(x-x3)+(x3-x2)*(y-y3))/detT;
					const double lam2 = ((y3-y1)*(x-x3)+(x1-x3)*(y-y3))/detT;
					const double lam3 = 1.0 - lam1 - lam2;

					if( (lam1 < 0 || lam2 < 0 || lam3 < 0 || detT == 0) )
					{
						//Not bounded by this cell

						delete [] nnIdx;
						delete [] dists;
					
						continue;
					}


					//not a boundary cell
					delete [] nnIdx;
					delete [] dists;
					break;
				}

			}
		}

		//clean up the ann pt
		annDeallocPt(qPt);

		//Set boundary height based on radial neighbor for edge cells
		Fluid.hgroove(0,all) = where(Fluid.hgroove(1,all) != 0,1,0);
		Fluid.hgroove(Fluid.M-1,all) = where(Fluid.hgroove(Fluid.M-2,all) != 0,1,0);
	}

	matlab(Fluid.boundary,"./Outputs/slipper.boundary");
}
void CSlipperGap::SlipperPressureBounds(void)
{
	//Update pressure boundary values
	Fluid.p = where(Fluid.boundary >= 0, Fluid.boundary, Fluid.p);
	Fluid.p = where(Fluid.boundary == -2, operatingslippergap.pG, Fluid.p);
	Fluid.p = where(Fluid.boundary == -3, operatingslippergap.pcase, Fluid.p);

	//The first and last row MUST be a boundary or we have problems
	for(int i=0; i<Fluid.M-1;)
	{
		int j=i+1;

		//Find the next boundary value
		while(Fluid.boundary(j,0) == -1)
		{
			j++;
			//We should put an error message here if j == Fluid.M
		}
		
		//Initialize the 2D pressure using linear interpolation between the boundaries	
		for(int ii=i+1; ii<j; ii++)
		{
			//The where statement preserves and special pressure boundaries set above
			Fluid.p(ii,all) = where(Fluid.boundary(ii,all) == -1,
												(Fluid.p(j,all)-Fluid.p(i,all))/
												(Fluid.r(j,all)-Fluid.r(i,all))*(Fluid.r(ii,all)-Fluid.r(i,all))+Fluid.p(i,all)
												,Fluid.p(ii,all));

		}

		i=j;
	}

	//Set the uncut field to this field
	Fluid.p_uncut = Fluid.p;
}
void CSlipperGap::cell2face(Array<double,3> & p, Array<double,3> & pe, Array<double,3> & pw, Array<double,3> & pn, Array<double,3> & ps, Array<double,3> & pt, Array<double,3> & pb)
{
	const int M = Fluid.M;
	const int N = Fluid.N;
	const int Q = Fluid.Q;

	/*
	//pn = (p*drp+pn*drn)/(drp+drn)
	pn(Range(0,M-2),all,all) = 
		(p(Range(0,M-2),all,all)*Fluid.dr(Range(0,M-2),all) +					//(p*drp +
		p(Range(1,M-1),all,all)*Fluid.dr(Range(1,M-1),all))						//pn*drn)
		/(Fluid.dr(Range(0,M-2),all)+Fluid.dr(Range(1,M-1),all));		// /(drp+drn);
	//"upwind" the boundary
	pn(M-1,all,all) = p(M-1,all,all);

	//ps = (p*drp+ps*drs)/(drp+drs)
	ps(Range(1,M-1),all,all) = 
		(p(Range(1,M-1),all,all)*Fluid.dr(Range(1,M-1),all) +					//(p*drp +
		p(Range(0,M-2),all,all)*Fluid.dr(Range(0,M-2),all))						//ps*drs)
		/(Fluid.dr(Range(1,M-1),all)+Fluid.dr(Range(0,M-2),all));		// /(drp+drs);
	//"upwind" the boundary
	ps(0,all,all) = p(0,all,all);
	
	//pe = (p*dtp+pe*dte)/(dtp+dte)
	pe(all,Range(0,N-2),all) =
		(p(all,Range(0,N-2),all)*Fluid.dtheta(all,Range(0,N-2)) +						//(p*dtp +
		p(all,Range(1,N-1),all)*Fluid.dtheta(all,Range(1,N-1)))						//pe*dte)
		/(Fluid.dtheta(all,Range(0,N-2))+Fluid.dtheta(all,Range(1,N-1)));		// /(dtp+dte);
	pe(all,N-1,all) =
		(p(all,N-1,all)*Fluid.dtheta(all,N-1) +												//(p*dtp +
		p(all,0,all)*Fluid.dtheta(all,0))														//pe*dte)
		/(Fluid.dtheta(all,N-1)+Fluid.dtheta(all,0));								// /(dtp+dte);

	//pw = (p*dtp+pw*dtw)/(dtp+dtw)
	pw(all,Range(1,N-1),all) =
		(p(all,Range(1,N-1),all)*Fluid.dtheta(all,Range(1,N-1)) +						//(p*dtp +
		p(all,Range(0,N-2),all)*Fluid.dtheta(all,Range(0,N-2)))						//pw*dtw)
		/(Fluid.dtheta(all,Range(1,N-1))+Fluid.dtheta(all,Range(0,N-2)));		// /(dtp+dtw);
	pw(all,0,all) =
		(p(all,0,all)*Fluid.dtheta(all,0) +													//(p*dtp +
		p(all,N-1,all)*Fluid.dtheta(all,N-1))												//pw*dtw)
		/(Fluid.dtheta(all,0)+Fluid.dtheta(all,N-1));								// /(dtp+dtw);
		*/
}
