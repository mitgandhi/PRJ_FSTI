# include "./oil.h"

// ------------------------------------------------------------------------- //
void limit(double &val, double min, double max)
{
	val = (val > min) ? val : min;
	val = (val < max) ? val : max;
}
// ------------------------------------------------------------------------- //
oil::oil(const input& inp) : in(inp)
{

	// read from caspar input
	lambda = in.data.oil.general.oillambda;
	cp = in.data.oil.general.oilC;

	p_min_limit = 0;
	p_max_limit = 0;
	T_min_limit = 0;
	T_max_limit = 0;

}
// ------------------------------------------------------------------------- //
constant_oil::constant_oil(const input& inp) : oil(inp)
{
	rho_0 = in.data.oil.constant_properties.oildensity;
	mu_0 = in.data.oil.constant_properties.oilviscosity;
	K_0 = in.data.oil.constant_properties.oilK;
}
// ------------------------------------------------------------------------- //
double constant_oil::get_rho(double p, double T)
{
	return rho_0;
}
// ------------------------------------------------------------------------- //
double constant_oil::get_mu(double p, double T)
{
	return mu_0;
}
// ------------------------------------------------------------------------- //
double constant_oil::get_K(double p, double T)
{
	return K_0;
}
// ------------------------------------------------------------------------- //
user_defined_oil::user_defined_oil(const input& inp) : oil(inp)
{
	
	// read from caspar input

	// copy the domain limits
	p_min_limit = in.data.oil.user_defined.pmin;
	p_max_limit = in.data.oil.user_defined.pmax;
	T_min_limit = in.data.oil.user_defined.Tmin;
	T_max_limit = in.data.oil.user_defined.Tmax;
		
	// copy the coefficients
	nuw = in.data.oil.user_defined.w;
	nupf = in.data.oil.user_defined.nupf;
	nup1 = in.data.oil.user_defined.nup1;
	nup2 = in.data.oil.user_defined.nup2;
	nuT1 = in.data.oil.user_defined.nuT1;
	nuT2 = in.data.oil.user_defined.nuT2;
	rhop1 = in.data.oil.user_defined.rhop1;
	rhop2 = in.data.oil.user_defined.rhop2;
	rhoT = in.data.oil.user_defined.rhoT;
	rhopT = in.data.oil.user_defined.rhopT;

	//

	// define the reference values
	p0 = 0;		// Pa
	T0 = 293;	// K
	rho0 = in.data.oil.user_defined.rho0;	// [kg/m3]


}
// ------------------------------------------------------------------------- //
double user_defined_oil::get_rho(double _p, double _T)
{
	
	// IMPORTANT: T must be in C, p in Pa

	double p = _p, T = _T;

	limit(p, p_min_limit, p_max_limit);
	limit(T, T_min_limit, T_max_limit);

	
	double dp = p - p0;
	double dT = (T + 273.0) - T0;
	double rho = rho0*(1 + rhoT*dT + rhop1*dp + rhop2*dp*dp + rhopT*dp*dT);
		
	return rho;
	
}
// ------------------------------------------------------------------------- //
double user_defined_oil::get_mu(double _p, double _T)
{

	// IMPORTANT: T must be in C, p in Pa
	// The coefficients of nu MUST be derived with mu expressed in cSt!

	double p = _p, T = _T;

	limit(p, p_min_limit, p_max_limit);
	limit(T, T_min_limit, T_max_limit);
		
	// nupf == 0 -> linear law, nupf == 1 -> power law
	double fp = (nupf == 0) ? nup1 + nup2*(T + 273.0) : nup1*pow(T + 273.0, nup2);
	
	double fT = nuT1 + nuT2*log10(T + 273.0);

	double nu = nuw*exp(fp*p)*(pow(10, (pow(10, fT))) - 0.5);	 // cSt
	double rho = get_rho(p, T);
	double mu = 1e-6*nu*rho; // Pa s

	return mu;
	
}
// ------------------------------------------------------------------------- //
double user_defined_oil::get_K(double _p, double _T)
{
	
	// IMPORTANT: T must be in C, p in Pa
	
	double p = _p, T = _T;

	limit(p, p_min_limit, p_max_limit);
	limit(T, T_min_limit, T_max_limit);

	double dp = p - p0;
	double dT = (T + 273.0) - T0;

	double num = 1.0 + rhoT*dT + rhop1*dp + rhop2*dp*dp + rhopT*dp*dT;
	double den = rhop1 + rhopT*dT + 2*rhop2*dp;
	double K = num/den;

	return K;
	
}

// ------------------------------------------------------------------------- //
user_defined2_oil::user_defined2_oil(const input& inp) : oil(inp)
{
	
	// read from caspar input

	// copy the domain limits
	p_min_limit = in.data.oil.user_defined.pmin;
	p_max_limit = in.data.oil.user_defined.pmax;
	T_min_limit = in.data.oil.user_defined.Tmin;
	T_max_limit = in.data.oil.user_defined.Tmax;
		
	// copy the coefficients
	/*nuw = in.data.oil.user_defined.w;
	nupf = in.data.oil.user_defined.nupf;
	nup1 = in.data.oil.user_defined.nup1;
	nup2 = in.data.oil.user_defined.nup2;
	nuT1 = in.data.oil.user_defined.nuT1;
	nuT2 = in.data.oil.user_defined.nuT2;
	nuw = in.data.oil.user_defined.w;*/
	
	C_2 = in.data.oil.user_defined.c_2;
	D_2 = in.data.oil.user_defined.d_2;
	G_0 = in.data.oil.user_defined.g_0;
	S_0 = in.data.oil.user_defined.s_0;
	
	rhop1 = in.data.oil.user_defined.rhop1;
	rhop2 = in.data.oil.user_defined.rhop2;
	rhoT = in.data.oil.user_defined.rhoT;
	rhopT = in.data.oil.user_defined.rhopT;
	rhop2T = in.data.oil.user_defined.rhop2T;

	//

	// define the reference values
	p0 = 1.0e5;		// Pa
	T0 = 25;		// C
	rho0 = in.data.oil.user_defined.rho0;	// [kg/m3]


}
// ------------------------------------------------------------------------- //
double user_defined2_oil::get_rho(double _p, double _T)
{
	
	// IMPORTANT: T must be in C, p in Pa

	double p = _p, T = _T;

	limit(p, p_min_limit, p_max_limit);
	limit(T, T_min_limit, T_max_limit);

	
	double dp = p - p0;
	double dT = T  - T0;
	// added second degree term 
	double rho = rho0*(1 + rhoT*dT + rhop1*dp + rhop2*dp*dp + rhopT*dp*dT+ rhop2T*dp*dp*dT);
		
	return rho;
	
}
// ------------------------------------------------------------------------- //
double user_defined2_oil::get_mu(double _p, double _T)
{

	// IMPORTANT: T must be in C, p in Pa
	// The coefficients of nu MUST be derived with mu expressed in cSt!

	double p = _p, T = _T;


	limit(p, p_min_limit, p_max_limit);
	limit(T, T_min_limit, T_max_limit);

	p = p/1.0e5;			// pa to bar
			
	double a = C_2*log10(1+(T/135.0))+D_2;
	double d = pow( (1.0+ (T/135.0) ) , S_0);
	double f = (G_0 * (pow( (1.0+(p/2000.0)) , a) / d) ) - 1.2 ;

	double mu = pow(10, f);	 // cSt
	
	mu = mu * 1.0e-3;		// cP to Pa s

	return mu;
	
}
// ------------------------------------------------------------------------- //
double user_defined2_oil::get_K(double _p, double _T)
{
	
	// IMPORTANT: T must be in C, p in Pa
	
	double p = _p, T = _T;

	limit(p, p_min_limit, p_max_limit);
	limit(T, T_min_limit, T_max_limit);

	double dp = p - p0;
	double dT = T - T0;

	double num = 1.0 + rhoT*dT + rhop1*dp + rhop2*dp*dp + rhopT*dp*dT + rhop2T*dp*dp*dT;
	double den = rhop1 + rhopT*dT + 2*rhop2*dp + 2*rhop2T*dp*dT;
	double K = num/den;

	return K;
	
}

// ------------------------------------------------------------------------- //
HLP32_oil::HLP32_oil(const input& inp) : oil(inp)
{

	// define the domain limits
	p_min_limit = 0.9e5;
	p_max_limit = 1000.0e5;
	T_min_limit = 0;
	T_max_limit = 130;

	// viscosity coeffs
	visco_0 = 1.702e-2;
	A0 = 2.278e-8;
	A1 = 1.477e-10;
	A2 = 5.127e-13;
	B = 5.428;
	C = 94.72;

	// density coeffs
	rs = 1047.03;			
	als = 5.761668e-4;
	a1 = 0.0732965;
	a2 = 1965.02;			
	a3 = -2.96813;		
}
// ------------------------------------------------------------------------- //
double HLP32_oil::get_cp(double _T)
{
	double T = _T;
	limit(T, T_min_limit, T_max_limit);

	double cp_25 = 1895.4; // 25C
	double cp_80 = 2118.1; // 80C

	// linear variation knowing cp @ two different temepratures
	return (cp_25 + (cp_80-cp_25)*(T - 25.0)/(80.0-25.0));

}
// ------------------------------------------------------------------------- //
double HLP32_oil::get_rho(double _p, double _T)
{
	
	double p = _p, T = _T;
	limit(p, p_min_limit, p_max_limit);
	limit(T, T_min_limit, T_max_limit);

	double TK = T + 273.15;
	double rho_T = rs*(1 - als*TK); 
	return rho_T/(1.0 - a1*log((a2 + a3*TK + p*1.0e-5)/(a2 + a3*TK)));
}
// ------------------------------------------------------------------------- //
double HLP32_oil::get_mu(double _p, double _T)
{
	
	double p = _p, T = _T;
	limit(p, p_min_limit, p_max_limit);
	limit(T, T_min_limit, T_max_limit);
	
	double alpha = A0 - A1*T + A2*(T*T);
	double Temp_coeff = B*(50.0 - T)/(C + T) + alpha*p;

	Temp_coeff = (Temp_coeff < 4.0) ? Temp_coeff : 4.0;
	return visco_0*exp(Temp_coeff);		

}
// ------------------------------------------------------------------------- //
double HLP32_oil::get_K(double _p, double _T)
{
	
	double p = _p, T = _T;
	limit(p, p_min_limit, p_max_limit);
	limit(T, T_min_limit, T_max_limit);
	double TK = T + 273.15;
	double rs = 1047.03;			// [kg/m3]
	double als = 5.761668e-4; // [1/K]
	double a1 = 0.0732965;
	double a2 = 1965.02;			// [bar]
	double a3 = -2.96813;			// [bar/K]
	double RhoT = rs*(1-als*TK); 
	double RhoTp = RhoT/(1 - a1*log((a2 + a3*TK + p*1.0e-5)/(a2 + a3*TK)));
	return 1.0e5*(RhoT*(a2 + a3*TK + p*1.0e-5)/(a1*RhoTp)); //[Pa]
}
// ------------------------------------------------------------------------- //
ExxonDTE10Excel32::ExxonDTE10Excel32(const input& inp) : oil(inp)
{

	// define the domain limits
	p_min_limit = 0.9e5;
	p_max_limit = 1400e5;
	T_min_limit = 20;
	T_max_limit = 130;
	
	// density / bulk modulus
	rho0	= 843.30;
	rho_T = -7.365E-04;
	rho_p = 6.345E-10;
	rho_p2 = -1.029E-18;
	rho_pT = 2.737E-12;
	rho_p2T = -8.683E-21;

	// viscosity
	nu_T0 = 3.983;
	nu_T1 = -1.587;
	nu_p10	= 0.8383;
	nu_p11	= -7.1136;
	nu_p12	= 1.5421;
	nu_p20	= -16.502;
	nu_p21	= 5.0636;
	nu_p22	= -2.1678;

}
// ------------------------------------------------------------------------- //
double ExxonDTE10Excel32::get_cp(double _T)
{

	double T = _T;
	limit(T, T_min_limit, T_max_limit);

	double cp_25 = 1895.4; // 25C
	double cp_80 = 2118.1; // 80C

	// linear variation knowing cp @ two different temepratures
	return (cp_25 + (cp_80-cp_25)*(T - 25.0)/(80.0-25.0));
}
// ------------------------------------------------------------------------- //
double ExxonDTE10Excel32::get_rho(double _p, double _T)
{
	double rho = 0;

	double p = _p, T = _T;
	limit(p, p_min_limit, p_max_limit);
	limit(T, T_min_limit, T_max_limit);

	// reference p and T
	double p0 = 1e5, T0 = 25;
	double dT = T - T0, dp = p - p0;

	return rho0*(1.0 + rho_T*dT + rho_p*dp + rho_p2*pow(dp,2) + rho_pT*dp*dT + rho_p2T*pow(dp,2)*dT); 

}
// ------------------------------------------------------------------------- //
double ExxonDTE10Excel32::get_mu(double _p, double _T)
{
	double mu = 0;

	// IMPORTANT: T must be in C, p in Pa
	// The coefficients of nu MUST be derived with mu expressed in cSt!

	double p = _p, T = _T;
	limit(p, p_min_limit, p_max_limit);
	limit(T, T_min_limit, T_max_limit);

	// reference p 
	double p0 = 1e5;
	double Tlog = log10(T);

	double nu_T = pow(10.0, nu_T0 + nu_T1*Tlog);
	double nu_p1 = pow(10.0, nu_p10 + nu_p11*Tlog + nu_p12*pow(Tlog,2.0));
	double nu_p2 = pow(10.0, nu_p20 + nu_p21*Tlog + nu_p22*pow(Tlog,2.0));

	// this is nu in cSt
	double nu = nu_T + nu_p1*(p - p0) + nu_p2*pow((p - p0), 2.0);

	// return the viscosity in Pas
	return 1e-6*nu*get_rho(p,T);
}
// ------------------------------------------------------------------------- //
double ExxonDTE10Excel32::get_K(double _p, double _T)
{
	double K = 0;

	double p = _p, T = _T;
	limit(p, p_min_limit, p_max_limit);
	limit(T, T_min_limit, T_max_limit);
	
	double p0 = 1e5, T0 = 25;
	double dT = T - T0, dp = p - p0;

	// K = rho/(drho/dp)_T
	double rho = get_rho(p,T);
	double drho_dp = rho0*(rho_p + 2.0*rho_p2*dp + rho_pT*dT + 2.0*rho_p2T*dp*dT);

	return (rho/drho_dp);
}
// ------------------------------------------------------------------------- //
MILPRF87257_oil::MILPRF87257_oil(const input& inp) : oil(inp)
{

	// define the domain limits
	p_min_limit = 0.9e5;
	p_max_limit = 1400e5;
	T_min_limit = -50;
	T_max_limit = 140;
	
	// density / bulk modulus
	rho0	= 835;   
	rho_T = -8.0528E-04;
	rho_p = 6.361E-10;
	rho_p2 = -2.725E-18;
	rho_pT = 1.817E-12;
	rho_p2T = -1.8480E-20;

	// viscosity coeffs
	cnu.p1 = -3.6021e-11;		
	cnu.p2 = 2.4487e-08;		
	cnu.T1 = 11.436;			
	cnu.T2 = -4.6274;	
	cnu.w = 1.294;

}
// ------------------------------------------------------------------------- //
double MILPRF87257_oil::get_cp(double _T)
{

	double T = _T;
	limit(T, T_min_limit, T_max_limit);

	double cp_25 = 1800; // 25C
	double cp_190 = 2500; // 80C

	// linear variation knowing cp @ two different temepratures
	return (cp_25 + ( cp_190 -cp_25)*(T - 25.0)/(190.0-25.0));
}
// ------------------------------------------------------------------------- //
double MILPRF87257_oil::get_rho(double _p, double _T)
{
	double rho = 0;

	double p = _p, T = _T;
	limit(p, p_min_limit, p_max_limit);
	limit(T, T_min_limit, T_max_limit);

	// reference p and T
	double p0 = 1e5, T0 = 25;
	double dT = T - T0, dp = p - p0;

	return rho0*(1.0 + rho_T*dT + rho_p*dp + rho_p2*pow(dp,2) + rho_pT*dp*dT + rho_p2T*pow(dp,2)*dT); 

}
// ------------------------------------------------------------------------- //
double MILPRF87257_oil::get_mu(double _p, double _T)
{
	//double mu = 0;

	// IMPORTANT: T must be in K, p in Pa
	// The coefficients of nu MUST be derived with mu expressed in cSt!

	double pw = _p, Tw = _T;

	limit(pw, p_min_limit, p_max_limit);
	limit(Tw, T_min_limit, T_max_limit);

	double fp = cnu.p1*(Tw + 273.15)+ cnu.p2;
	double fT = cnu.T1 + cnu.T2*log10(Tw + 273.15);

	double nu = cnu.w*exp(fp*pw)*(pow(10, (pow(10, fT))) - 0.5);	 // cSt
	double rho = get_rho(pw, Tw);
	double mu = 1e-6*nu*rho; // Pa s

	return mu;
}
// ------------------------------------------------------------------------- //
double MILPRF87257_oil::get_K(double _p, double _T)
{
	double K = 0;

	double p = _p, T = _T;
	limit(p, p_min_limit, p_max_limit);
	limit(T, T_min_limit, T_max_limit);
	
	double p0 = 1e5, T0 = 25;
	double dT = T - T0, dp = p - p0;

	// K = rho/(drho/dp)_T
	double rho = get_rho(p,T);
	double drho_dp = rho0*(rho_p + 2.0*rho_p2*dp + rho_pT*dT + 2.0*rho_p2T*dp*dT);

	return (rho/drho_dp);
}
// ------------------------------------------------------------------------- //
skydrol_oil::skydrol_oil(const input& inp) : oil(inp)
{

	// define the domain limits
	p_min_limit = 0.9e5;
	p_max_limit = 1000e5;
	T_min_limit = -50;
	T_max_limit = 140;

	// reference values
	p0 = 0;
	T0 = 293.0;					// K (20 C)
	rho0 = 1008.0;			// Kg/m3

	// viscosity coeffs
	cnu.p1 = 3.0649e-6;		
	cnu.p2 = -8.5688e-1;		
	cnu.T1 = 7.69698;			
	cnu.T2 = -3.07574;	
	cnu.w = 0.85;
	
	// density coeffs
	crho.T = -7.9365e-4;
	crho.p1 = 7.2993e-10;
	crho.p2 = -1.2951e-18;
	crho.pT = 3.9360e-12;

}
// ------------------------------------------------------------------------- //
double skydrol_oil::get_rho(double _p, double _T)
{

	// IMPORTANT: T must be in C, p in Pa

	double pw = _p, Tw = _T;
	
	limit(pw, p_min_limit, p_max_limit);
	limit(Tw, T_min_limit, T_max_limit);

	double dp = pw - p0;
	double dT = (Tw + 273.0) - T0;

	double rho = rho0*(1 + crho.T*dT + crho.p1*dp + crho.p2*dp*dp + crho.pT*dp*dT);
		
	return rho;
}
// ------------------------------------------------------------------------- //
double skydrol_oil::get_K(double _p, double _T)
{

	// IMPORTANT: T must be in C, p in Pa

	double pw = _p, Tw = _T;
	
	limit(pw, p_min_limit, p_max_limit);
	limit(Tw, T_min_limit, T_max_limit);
	
	double dp = pw - p0;
	double dT = (Tw + 273.0) - T0;

	double num = 1.0 + crho.T*dT + crho.p1*dp + crho.p2*dp*dp + crho.pT*dp*dT;
	double den = crho.p1 + crho.pT*dT + 2*crho.p2*dp;
	double K = num/den;

	return K;

}
// ------------------------------------------------------------------------- //
double skydrol_oil::get_mu(double _p, double _T)
{

	// IMPORTANT: T must be in C, p in Pa

	double pw = _p, Tw = _T;

	limit(pw, p_min_limit, p_max_limit);
	limit(Tw, T_min_limit, T_max_limit);

	double fp = cnu.p1*pow(Tw + 273.0, cnu.p2);
	double fT = cnu.T1 + cnu.T2*log10(Tw + 273.0);

	double nu = cnu.w*exp(fp*pw)*(pow(10, (pow(10, fT))) - 0.5);	 // cSt
	double rho = get_rho(pw, Tw);
	double mu = 1e-6*nu*rho; // Pa s

	return mu;

}
// ------------------------------------------------------------------------- //
red_oil::red_oil(const input& inp) : oil(inp)
{
	// define the domain limits
	p_min_limit = 0.9e5;
	p_max_limit = 1000e5;
	T_min_limit = -50;
	T_max_limit = 140;

	// reference values
	p0 = 0;
	T0 = 313.15;		// K (40 C)
	rho0 = 825;			// Kg/m3

	// viscosity coeffs
	cnu.p1 = -3.2988e-11;		
	cnu.p2 = 2.28150e-08;		
	cnu.T1 = 10.44885;			
	cnu.T2 = -4.166341;	
	cnu.w = 1.11;
	
	// density coeffs
	crho.T = -7.5e-4;
	crho.p1 = 6.6667e-10;
	crho.p2 = -6.8359e-19;
	crho.pT = 3.78e-12;
	crho.pT2 = -1.0e-22;

}
// ------------------------------------------------------------------------- //
double red_oil::get_rho(double _p, double _T)
{

	// IMPORTANT: T must be in C, p in Pa

	double pw = _p, Tw = _T;
	
	limit(pw, p_min_limit, p_max_limit);
	limit(Tw, T_min_limit, T_max_limit);

	double dp = pw - p0;
	double dT = (Tw + 273.0) - T0;

	double rho = rho0*(1 + crho.T*dT + crho.p1*dp + crho.p2*dp*dp + crho.pT*dp*dT + crho.pT2*dp*dT*dp*dT);
		
	return rho;
}
// ------------------------------------------------------------------------- //
double red_oil::get_K(double _p, double _T)
{

	// IMPORTANT: T must be in C, p in Pa

	double pw = _p, Tw = _T;
	
	limit(pw, p_min_limit, p_max_limit);
	limit(Tw, T_min_limit, T_max_limit);
	
	double dp = pw - p0;
	double dT = (Tw + 273.0) - T0;

	double num = 1.0 + crho.T*dT + crho.p1*dp + crho.p2*dp*dp + crho.pT*dp*dT + crho.pT2*dp*dT*dp*dT;
	double den = crho.p1 + crho.pT*dT + 2*crho.p2*dp + 2*crho.pT2*dp*dT*dT;
	double K = num/den;

	return K;

}
// ------------------------------------------------------------------------- //
double red_oil::get_mu(double _p, double _T)
{

	// IMPORTANT: T must be in C, p in Pa

	double pw = _p, Tw = _T;

	limit(pw, p_min_limit, p_max_limit);
	limit(Tw, T_min_limit, T_max_limit);

	double fp = cnu.p1*(Tw + 273.0)+ cnu.p2;
	double fT = cnu.T1 + cnu.T2*log10(Tw + 273.0);

	double nu = cnu.w*exp(fp*pw)*(pow(10, (pow(10, fT))) - 0.5);	 // cSt
	double rho = get_rho(pw, Tw);
	double mu = 1e-6*nu*rho; // Pa s

	return mu;

}
// ------------------------------------------------------------------------- //
SAE10W::SAE10W(const input& inp) : oil(inp)
{
	// define the domain limits
	p_min_limit = 0.9e5;
	p_max_limit = 1000e5;
	T_min_limit = 0;
	T_max_limit = 120;
}
// ------------------------------------------------------------------------- //
double SAE10W::get_rho(double _p, double _T)
{

	// IMPORTANT: T must be in C, p in Pa

	double pw = _p, Tw = _T;
	
	limit(pw, p_min_limit, p_max_limit);
	limit(Tw, T_min_limit, T_max_limit);

	double rho_T = 1000*(-0.0006*Tw + 0.8819);

	double rho = rho_T + 9.8*(pw/1e10)*rho_T/(1e-9*get_K(pw,Tw));

	return rho;
}
// ------------------------------------------------------------------------- //
double SAE10W::get_K(double _p, double _T)
{

	// IMPORTANT: T must be in C, p in Pa

	double pw = _p, Tw = _T;
	
	limit(pw, p_min_limit, p_max_limit);
	limit(Tw, T_min_limit, T_max_limit);

	double K1 = (1.51 + 7.0*(0.87 - 0.86));						// GPa
	double K2 = 1e4*exp(0.0023*(20-Tw));							// GPa
	double K3 = 5.6*(pw/1e5);														// GPa

	double K = (K1*K2 + K3)*1e5; // from GPa to Pa

	return K;
}
// ------------------------------------------------------------------------- //
double SAE10W::get_mu(double _p, double _T)
{

	// IMPORTANT: T must be in C, p in Pa

	double pw = _p, Tw = _T;

	limit(pw, p_min_limit, p_max_limit);
	limit(Tw, T_min_limit, T_max_limit);

	double fp = exp(128014.0*pow(Tw + 273.0,-3.184)*(pw/1e5));
	double fT = 9.4252 - 3.7*log10(Tw + 273.0);

	double mu = fp*(pow(10.0, pow(10.0, fT)) - 0.7)/1000.0; 

	return mu;

}
// ------------------------------------------------------------------------- //
milh5606_oil::milh5606_oil(const input& inp) : oil(inp)
{
	// define the domain limits
	p_min_limit = 0.9e5;
	p_max_limit = 1000e5;
	T_min_limit = -50;
	T_max_limit = 140;

	// reference values
	p0 = 0;
	T0 = 313.15;		// K (40 C)
	rho0 = 842.5;			// Kg/m3

	// viscosity coeffs
	cnu.p1 = -9.7961e-11;		
	cnu.p2 = 5.26414e-08;		
	cnu.T1 = 7.6245;			
	cnu.T2 = -3.03896;	
	cnu.w = 1.05;
	
	// density coeffs
	crho.T = -0.00088;
	crho.p1 = 7.7e-10;
	crho.p2 = -9.5e-19;
	crho.pT = 4.056e-12;
	crho.pT2 = -9e-23;

}
// ------------------------------------------------------------------------- //
double milh5606_oil::get_rho(double _p, double _T)
{

	// IMPORTANT: T must be in C, p in Pa

	double pw = _p, Tw = _T;
	
	limit(pw, p_min_limit, p_max_limit);
	limit(Tw, T_min_limit, T_max_limit);

	double dp = pw - p0;
	double dT = (Tw + 273.0) - T0;

	double rho = rho0*(1 + crho.T*dT + crho.p1*dp + crho.p2*dp*dp + crho.pT*dp*dT + crho.pT2*dp*dT*dp*dT);
		
	return rho;
}
// ------------------------------------------------------------------------- //
double milh5606_oil::get_K(double _p, double _T)
{

	// IMPORTANT: T must be in C, p in Pa

	double pw = _p, Tw = _T;
	
	limit(pw, p_min_limit, p_max_limit);
	limit(Tw, T_min_limit, T_max_limit);
	
	double dp = pw - p0;
	double dT = (Tw + 273.0) - T0;

	double num = 1.0 + crho.T*dT + crho.p1*dp + crho.p2*dp*dp + crho.pT*dp*dT + crho.pT2*dp*dT*dp*dT;
	double den = crho.p1 + crho.pT*dT + 2*crho.p2*dp + 2*crho.pT2*dp*dT*dT;
	double K = num/den;

	return K;

}
// ------------------------------------------------------------------------- //
double milh5606_oil::get_mu(double _p, double _T)
{

	// IMPORTANT: T must be in C, p in Pa

	double pw = _p, Tw = _T;

	limit(pw, p_min_limit, p_max_limit);
	limit(Tw, T_min_limit, T_max_limit);

	double fp = cnu.p1*(Tw + 273.0)+ cnu.p2;
	double fT = cnu.T1 + cnu.T2*log10(Tw + 273.0);

	double nu = cnu.w*exp(fp*pw)*(pow(10, (pow(10, fT))) - 0.5);	 // cSt
	double rho = get_rho(pw, Tw);
	double mu = 1e-6*nu*rho; // Pa s

	return mu;

}
// ------------------------------------------------------------------------- //
water_oil::water_oil(const input& inp) : oil(inp)
{
	// define the domain limits (in C and Pa)
	p_min_limit = 0.9e5;
	p_max_limit = 5000e5;
	T_min_limit = -50;
	T_max_limit = 140;

	//setup density coefficients

		density.p00 =        1037;
		density.p10 =      -42.23;
		density.p01 =       45.35;
		density.p11 =       3.352;
		density.p02 =      -4.862;
      

	//setup bulk modulus coefficients

       K.p00 =        3551;
       K.p10 =      -434.8;
       K.p01 =       879.1;
       K.p20 =      -145.4;
       K.p11 =       13.14;
       K.p02 =      -23.77;
       K.p30 =        48.1;
       K.p21 =       19.08;
       K.p12 =      -16.42;

	//setup viscosity coefficients

	   Mu.p00 =       306.1;
       Mu.p10 =      -172.4;
       Mu.p01 =       31.42;
       Mu.p20 =       53.81;
       Mu.p11 =      -13.05;
       Mu.p02 =      -4.493;
       Mu.p30 =      -26.78;
       Mu.p21 =       7.016;
       Mu.p12 =      0.6557;
       Mu.p03 =      0.0416;
       Mu.p40 =       55.03;
       Mu.p31 =      -5.841;
       Mu.p22 =       8.413;
       Mu.p13 =       1.029;
       Mu.p04 =      0.1655;
       Mu.p50 =      -24.63;
       Mu.p41 =       2.577;
       Mu.p32 =      -5.568;
       Mu.p23 =     -0.6933;
       Mu.p14 =     -0.2535;
       Mu.p05 =    -0.08575;

}
double water_oil::get_rho(double _p, double _T)
{
	//IMPORTANT: T must be in C, p in Pa
                
	limit(_p, p_min_limit, p_max_limit);
	limit(_T, T_min_limit, T_max_limit);

	//convert _T to K and _p to MPa
	double T = _T+273.15;
	double P = _p/1e6;
    
	density.x=(T-390)/69.24;
	density.y=(P-252.5)/144.4;

	//Finding density:
        
	double rho=density.p00 + density.p10*density.x + density.p01*density.y + density.p11*density.x*density.y + density.p02*pow(density.y,2);
	
	return rho;
}
double water_oil::get_K(double _p, double _T)
{
	// IMPORTANT: T must be in C, p in Pa

	limit(_p, p_min_limit, p_max_limit);
	limit(_T, T_min_limit, T_max_limit);

	//convert _T to K and _p to MPa
	double T = _T+273.15;
	double p = _p/1e6;

	K.x=(T-390)/69.24;
	K.y=(p-252.5)/144.4;
	
	double BulkModulus=(K.p00 + K.p10*K.x + K.p01*K.y + K.p20*pow(K.x,2) + K.p11*K.x*K.y + K.p02*pow(K.y,2) + K.p30*pow(K.x,3) + K.p21*pow(K.x,2)*K.y + K.p12*K.x*pow(K.y,2))*1e6;

	return BulkModulus;
}
double water_oil::get_mu(double _p, double _T)
{
	// IMPORTANT: T must be in C, p in Pa

	limit(_p, p_min_limit, p_max_limit);
	limit(_T, T_min_limit, T_max_limit);

	//convert _T to K and _p to MPa
	double T = _T+273.15;
	double p = _p/1e6;

     Mu.x=(T-390)/69.24;
     Mu.y=(p-252.5)/144.4;

	double DynamicViscosity=(Mu.p00 + 1.11*Mu.p10*Mu.x + 0.9*Mu.p01*Mu.y + Mu.p20*pow(Mu.x,2) + 0.95*Mu.p11*Mu.x*Mu.y + 1.15*Mu.p02*pow(Mu.y,2) + Mu.p30*pow(Mu.x,3) + Mu.p21*pow(Mu.x,2)*Mu.y + Mu.p12*Mu.x*pow(Mu.y,2) + Mu.p03*pow(Mu.y,3) + Mu.p40*pow(Mu.x,4) + Mu.p31*pow(Mu.x,3)*Mu.y + Mu.p22*pow(Mu.x,2)*pow(Mu.y,2)  + Mu.p13*Mu.x*pow(Mu.y,3) + Mu.p04*pow(Mu.y,4) + Mu.p50*pow(Mu.x,5) + Mu.p41*pow(Mu.x,4)*Mu.y + Mu.p32*pow(Mu.x,3)*pow(Mu.y,2)  + Mu.p23*pow(Mu.x,2)*pow(Mu.y,3) + Mu.p14*Mu.x*pow(Mu.y,4) + Mu.p05*pow(Mu.y,5))*1e-6;

	return DynamicViscosity;
}
//-------------------------------------------------------------------------------------------//

ISO46_oil::ISO46_oil(const input& inp) : oil(inp)
{

	// define the domain limits
	p_min_limit = 0.9e5;
	p_max_limit = 1000.0e5;
	T_min_limit = -40;
	T_max_limit = 120;

	// viscosity coeffs
	TC1 = 8.800722708112780;
	TC2 = -3.438322049472131;
	alpha = 2.475e-8;

	// density coeffs
	TR1 = 956.0191915132108;
	TR2 = -0.3538895116093;
	PR1 = 0.6/1e9;
	PR2 = 1.17/1e9;
}
// ------------------------------------------------------------------------- //
double ISO46_oil::get_rho(double _p, double _T)
{
	double p = _p, T = _T;
	limit(p, p_min_limit, p_max_limit);
	limit(T, T_min_limit, T_max_limit);
	double TK = _T + 273.15;

	double rho_T = TR1+TR2*TK; 
	return rho_T*(1+(PR1*p)/(1+PR2*p));
}
// ------------------------------------------------------------------------- //
double ISO46_oil::get_mu(double _p, double _T)
{
	double p = _p, T = _T;
	limit(p, p_min_limit, p_max_limit);
	limit(T, T_min_limit, T_max_limit);
	double TK = _T + 273.15;

	double nu_0 = pow(10, pow(10, TC1+TC2*log10(TK)));
	double nu = nu_0*exp(alpha*p);
	return 1e-6*nu*get_rho(_p,_T);
}
// ------------------------------------------------------------------------- //
double ISO46_oil::get_K(double _p, double _T)
{
	double p = _p, T = _T;
	limit(p, p_min_limit, p_max_limit);
	limit(T, T_min_limit, T_max_limit);
	double TK = _T + 273.15;

	double rho_T = TR1+TR2*TK;
	double drhodp = rho_T*PR1/pow((PR2*p+1),2.0);
	return get_rho(_p,_T)*1.0/drhodp;
}
// ------------------------------------------------------------------------- //
Parker_skydrol_oil::Parker_skydrol_oil(const input& inp) : oil(inp)
{
	// define the domain limits
	p_min_limit = 0.9e5;
	p_max_limit = 1000e5;
	T_min_limit = -50;
	T_max_limit = 140;

	// reference values
	p0 = 0;
	T0 = 293.0;					// K (20 C)
	rho0 = 1008.0;				// Kg/m3
	


	//**************************** Begin Parker Code *******************************//
	//Reference temperature and pressure
	tref = 100; //F
	pref = 0;	//psi

		//  Specific gravity coefficients for Skydrol LD-4
	crho.A = -0.000465142;
	crho.B = 1.04679;

	// Viscosity coefficients for Skydrol LD-4
	cnu.A = 8.57158850327909;
	cnu.B = 3.10174401619516;
	cnu.C = 0.60556845;

	//**************************** End Parker Code *******************************//

}
// ------------------------------------------------------------------------- //
double Parker_skydrol_oil::get_rho(double p, double T)
{

	// IMPORTANT: T must be in C, p in Pa

	double pw = p, Tw = T;

	limit(pw, p_min_limit, p_max_limit);
	limit(Tw, T_min_limit, T_max_limit);

	

	//**************************** Begin Parker Code *******************************//


	//Specified temperature and pressure
	crho.T_degF = Tw*9/5+32;
	crho.p_psi = pw*14.504*1e-5;


	// Density
	gamma_ref = 1.25 - 0.000405 * tref;											//Ratio of specific heats at reference temperature
    Beta_at_ref = get_K(pref/14.504*1e5,(tref-32)*5/9) * 1e-5 * 14.504;          //Adiabatic tangent bulk modulus at reference conditions
    Beta_it_ref = Beta_at_ref / gamma_ref;									//isothermal tangent bulk modulus at reference conditions

    gamma = 1.25 - 0.000405 * crho.T_degF;                                 //Ratio of specific heats at specified temperature
    Beta_at = get_K(pw,Tw) * 1e-5 * 14.504;									 //Adiabatic tangent bulk modulus at specified conditions
    Beta_it = Beta_at / gamma;                                   //Isothermal tangent bulk modulus at specified conditions
    Beta_as = (Beta_at_ref + Beta_at) / 2;                       //Adiabatic secant bulk modulus at specified conditions
    Beta_is = (Beta_it_ref + Beta_it) / 2;                      //Isothermal secant bulk modulus at specified conditions

    p_ref = 0;																//Reference pressure for bulk modulus
    SG = (crho.B + crho.A * crho.T_degF) * ((crho.p_psi - p_ref) / Beta_is + 1);      //Calc SG and correct for pressure

    dens_SI = SG * 1000; //kg/m^3
    dens_EN = dens_SI * 3.612729e-005; 
	
	//**************************** End Parker Code *******************************//
		
	rho = dens_SI;
	return rho;
}
// ------------------------------------------------------------------------- //
double Parker_skydrol_oil::get_K(double p, double T)
{

	// IMPORTANT: T must be in C, p in Pa

	double pw = p, Tw = T;
	
	limit(pw, p_min_limit, p_max_limit);
	limit(Tw, T_min_limit, T_max_limit);
	
	
	//**************************** Begin Parker Code *******************************//


	//Adiabatic tangent bulk modulus
	//valid for type 4 phosphate ester fluid per AS1241

	//Convert C to F
	cK.T_degF = Tw * 9/5 + 32;

	//Convert Pa to psi
	cK.p_psi = pw * 1e-5 * 14.504;

	if (cK.T_degF > 75)
	{
		BE1 = (-0.0000003593732 * cK.p_psi - 0.008160747) * pow(cK.T_degF,3);
		BE2 = (0.0002268273 * cK.p_psi + 6.395661) * pow(cK.T_degF,2);
		BE3 = (-0.06115304 * cK.p_psi - 1862.197) * cK.T_degF;
		BE4 = (24.15663 * cK.p_psi + 305663.3);
		B = BE1 + BE2 + BE3 + BE4;
	}
	else{
		slope = -0.03318917 * cK.p_psi - 1040.634;
		BATP75 = (-0.0000003593732 * cK.p_psi - 0.008160747) * pow(75.0,3) + (0.0002268273 * cK.p_psi + 6.395661) * pow(75.0,2) + (-0.06115304 * cK.p_psi - 1862.197) * 75 + (24.15663 * cK.p_psi + 305663.3);
		YINTP0 = BATP75 - slope * 75;
		B = slope * cK.T_degF + YINTP0;	//psi
	}

	K = B / 14.504 * 1e5; //Pascal

	//**************************** End Parker Code *******************************//

	return K;
}
// ------------------------------------------------------------------------- //
double Parker_skydrol_oil::get_mu(double p, double T)
{

	// IMPORTANT: T must be in C, p in Pa

	double pw = p, Tw = T;
	
	limit(pw, p_min_limit, p_max_limit);
	limit(Tw, T_min_limit, T_max_limit);

	
	
	//**************************** Begin Parker Code *******************************//

	//KINEMATIC VISCOSITY

    //Fitting function is:  LOG LOG (v + .6) = A - B*LOG(T + 459.69) (MacCoull viscosity-temperature equation)
    //where v = kinematic viscosity, cs; T = temperature in degF
    //Pressure effect on viscosity is then taken into account with the Barus Equation:
    //vp = v0 e^(2.3(alpha)^a(P)(10^4))
	
    //Convert C to F
    cnu.T_degF = Tw * 9/5 + 32;

    //Convert Pa to psi
    cnu.p_psi = pw * 1e-5 * 14.504;

    //Calculate 1 Atm Viscosity
    V0 = pow(10,pow(10,(cnu.A - cnu.B* log10(cnu.T_degF + 460)))) - 0.6;   

    //Correct for pressure
    KlauseA = 3.23523 - 11.3886 * cnu.C + 13.1735 * pow(cnu.C,2) - 4.8881 * pow(cnu.C,3);
    KlauseB = -5.33425 + 19.9521 * cnu.C - 23.9448 * pow(cnu.C,2) + 10.155 * pow(cnu.C,3);
    KlauseC = 3.35452 - 13.1273 * cnu.C + 17.1712 * pow(cnu.C,2) - 7.6551 * pow(cnu.C,3);
    alpha = KlauseA + KlauseB * log10(V0) + KlauseC * pow((log10(V0)),2);
    Kin_Vis = V0 * exp(2.3 * pow(alpha,(560 / (cnu.T_degF + 460))) * cnu.p_psi * 0.0001); //cSt
			
	//**************************** End Parker Code *******************************//
		

	rho = get_rho(pw, Tw);
	mu = (Kin_Vis/1e6)*rho;   // Pa*s
	return mu;
}

