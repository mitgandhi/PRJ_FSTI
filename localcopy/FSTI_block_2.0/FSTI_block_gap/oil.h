# ifndef __oil__
# define __oil__

# include "../FSTI_input/input.h"

// ------------------------------------------------------------------------- //
// general oil class
class oil
{

protected:

	const input& in;

	// domain of definition
	double p_min_limit;
	double p_max_limit;
	double T_min_limit;
	double T_max_limit;

	double lambda;		// thermal conductivity [W/m2]
	double cp;				// heat capacity [J/kgK]

public:

	oil(const input& in);
	
	virtual double get_lambda(double T) = 0;
	virtual double get_cp(double T) = 0;
	virtual double get_rho(double p, double T) = 0;
	virtual double get_mu(double p, double T) = 0;
	virtual double get_K(double p, double T) = 0;

	double get_p_min_limit() const { return p_min_limit; }
	double get_p_max_limit() const { return p_max_limit; }
	double get_T_min_limit() const { return T_min_limit; }
	double get_T_max_limit() const { return T_max_limit; }

};
// ------------------------------------------------------------------------- //
// oil with constant properties
class constant_oil : public oil
{

	double rho_0;			// [kg/m3]
	double mu_0;			// [Pas]
	double K_0;				// [Pa]
			
public:

	constant_oil(const input& in);

	double get_cp(double T = 0) { return cp; }
	double get_lambda(double T = 0) { return lambda; }
	double get_rho(double p = 0, double T = 0);
	double get_mu(double p = 0, double T = 0);
	double get_K(double p, double T);

};
// ------------------------------------------------------------------------- //
// oil with user defined properties
class user_defined_oil : public oil
{
	// reference values
	double p0;
	double T0;
	double rho0;

	// coefficients
	double nuw;
	int nupf;
	double nup1;
	double nup2;
	double nuT1;
	double nuT2;
	double rhop1;
	double rhop2;
	double rhoT;
	double rhopT;

public:

	user_defined_oil(const input& in);

	double get_cp(double T = 0) { return cp; }
	double get_lambda(double T = 0) { return lambda; }
	double get_rho(double p, double T);
	double get_mu(double p, double T);
	double get_K(double p, double T);

};

// ------------------------------------------------------------------------- //
// oil with user defined properties
// NEW for Daniel Hasko 09 28 2016 change made by Rene Chacon
class user_defined2_oil : public oil
{
	// reference values
	double p0;
	double T0;
	double rho0;

	// coefficients
	//double nuw;
	//int nupf;
	//double nup1;
	//double nup2;
	//double nuT1;
	//double nuT2;
	
	// Roelands eq empirical coefficients:  (added by Rene Chacon 09 28 2016)
	double C_2;
	double D_2;
	double S_0;
	double G_0;
	//
	double nuT2;
	double rhop1;
	double rhop2;
	double rhoT;
	double rhopT;
	double rhop2T;  // NEW second degree term coefficient

public:

	user_defined2_oil(const input& in);

	double get_cp(double T = 0) { return cp; }
	double get_lambda(double T = 0) { return lambda; }
	double get_rho(double p, double T);
	double get_mu(double p, double T);
	double get_K(double p, double T);

};

// HLP32 oil
class HLP32_oil : public oil
{
	
	// viscosity coeffs
	double visco_0;
	double A0;
	double A1;
	double A2;
	double B;
	double C;

	// density coeffs
	double rs;			// [kg/m3]
	double als;			// [1/K]
	double a1;		
	double a2;			// [bar]
	double a3; 			// [bar/K]

public:

	HLP32_oil(const input& in);

	double get_cp(double T = 0);
	double get_lambda(double T = 0) { return lambda; }
	double get_rho(double p, double T);
	double get_mu(double p, double T);
	double get_K(double p, double T);

};
// ------------------------------------------------------------------------- //
class ExxonDTE10Excel32 : public oil
{

	// density/bulk modulus
	double rho0;

	double rho_p, rho_p2;
	double rho_T;
	double rho_pT, rho_p2T;

	// viscosity
	double nu_T0, nu_T1;
	double nu_p10, nu_p11, nu_p12;
	double nu_p20, nu_p21, nu_p22;

public:

	ExxonDTE10Excel32(const input& in);

	double get_cp(double T = 0);
	double get_lambda(double T = 0) { return lambda; }
	double get_rho(double p, double T);
	double get_mu(double p, double T);
	double get_K(double p, double T);

};

// ------------------------------------------------------------------------- //
class MILPRF87257_oil : public oil
{

	// density/bulk modulus
	double rho0;

	double rho_p, rho_p2;
	double rho_T;
	double rho_pT, rho_p2T;

	// viscosity
	struct // coefficients for kinematic viscosity	
	{
		double w, p1, p2, T1, T2;
	} cnu;

public:

	MILPRF87257_oil(const input& in);

	double get_cp(double T = 0);
	double get_lambda(double T = 0) { return lambda; }
	double get_rho(double p, double T);
	double get_mu(double p, double T);
	double get_K(double p, double T);

};

// ------------------------------------------------------------------------- //
// skydrol oil
class skydrol_oil : public oil
{

	// reference values 
	double p0;			// [Pa]
	double T0;			// [degC]
	double rho0;		// [kg/m3]

	struct	// corfficients for density
	{
		double p1, p2, T, pT;
	} crho;
	struct // coefficients for kinematic viscosity	
	{
		double w, p1, p2, T1, T2;
	} cnu;

public:

	skydrol_oil(const input& in);

	double get_cp(double T = 0) { return cp; }
	double get_lambda(double T = 0) { return lambda; }
	double get_rho(double p, double T);
	double get_mu(double p, double T);
	double get_K(double p, double T);

};
// aerospace "red oil"
class red_oil : public oil
{

	// reference values 
	double p0;			// [Pa]
	double T0;			// [degC]
	double rho0;		// [kg/m3]

	struct	// coefficients for density
	{
		double p1, p2, T, pT, pT2;
	} crho;
	struct // coefficients for kinematic viscosity	
	{
		double w, p1, p2, T1, T2;
	} cnu;

public:

	red_oil(const input& in);

	double get_cp(double T = 0) { return cp; }
	double get_lambda(double T = 0) { return lambda; }
	double get_rho(double p, double T);
	double get_mu(double p, double T);
	double get_K(double p, double T);

};
// ------------------------------------------------------------------------- //
// SAE10W
class SAE10W : public oil
{

public:

	SAE10W(const input& in);

	double get_cp(double T = 0) { return cp; }
	double get_lambda(double T = 0) { return lambda; }
	double get_rho(double p, double T);
	double get_mu(double p, double T);
	double get_K(double p, double T);

};
//aerospace "MIL-H-5606" (Low temperatures)
class milh5606_oil : public oil
{

	// reference values 
	double p0;			// [Pa]
	double T0;			// [degC]
	double rho0;		// [kg/m3]

	struct	// coefficients for density
	{
		double p1, p2, T, pT, pT2;
	} crho;
	struct // coefficients for kinematic viscosity	
	{
		double w, p1, p2, T1, T2;
	} cnu;

public:

	milh5606_oil(const input& in);

	double get_cp(double T = 0) { return cp; }
	double get_lambda(double T = 0) { return lambda; }
	double get_rho(double p, double T);
	double get_mu(double p, double T);
	double get_K(double p, double T);

};



// ------------------------------------------------------------------------- //
class water_oil : public oil
{
		struct 
	{
		double x;
		double y;
           
		double p00;
		double p10;
		double p01;
		double p11;
		double p02;

	} density;

	struct 
	{
	   double x;
	   double y;
       double p00;
       double p10;
       double p01;
       double p20;
       double p11;
       double p02;
       double p30;
       double p21;
       double p12;
	} K;

	struct
	{
		double x;
		double y;
		double p00;
        double p10;
        double p01;
        double p20;
        double p11;
        double p02;
        double p30;
        double p21;
        double p12;
        double p03;
        double p40;
        double p31;
        double p22;
        double p13;
        double p04;
        double p50;
        double p41;
        double p32;
        double p23;
        double p14;
        double p05;
	} Mu;

public:
	water_oil(const input& in);
	double get_cp(double T = 0) { return cp; }
	double get_lambda(double T = 0) { return lambda; }
	double get_rho(double p, double T);
	double get_mu(double p, double T);
	double get_K(double p, double T);
};

class ISO46_oil : public oil
{

	// viscosity coeffs
	double TC1;
	double TC2;
	double alpha;

	// density coeffs
	double TR1;
	double TR2;
	double PR1;
	double PR2;
public:

	ISO46_oil(const input& in);
	double get_cp(double T = 0) { return cp; }
	double get_lambda(double T = 0) { return lambda; }
	double get_rho(double p, double T);
	double get_mu(double p, double T);
	double get_K(double p, double T);

};

class Parker_skydrol_oil : public oil
{
	// reference values 
	double p0;			// [Pa]
	double T0;			// [degC]
	double rho0;		// [kg/m3]
//- Begin Parker Code -//
	//Reference temperature and pressure
    double tref;	//F
    double pref;	//psi

		// Density variables and constants
	double  gamma_ref, Beta_at_ref, Beta_it_ref, gamma, 
			Beta_at, Beta_it, Beta_as, Beta_is, p_ref, SG, dens_SI, dens_EN, rho;

		// Bulk modulus variables and constants
	double p, T, K, BE1, BE2, BE3, BE4, B, slope, BATP75, YINTP0;
		
		// Kinematic viscosity variables and constants
	double  V0, KlauseA, KlauseB, KlauseC, alpha, Kin_Vis, mu;


	struct	// corfficients for density
	{
		double A, B, T_degF, p_psi;
	} crho;
	
	struct // coefficients for bulk modulus
	{
		double T_degF, p_psi;
	} cK;
	
	struct // coefficients for kinematic viscosity	
	{
		double A, B, C, T_degF, p_psi;
	} cnu;

//- End Parker Code -/


public:

	Parker_skydrol_oil(const input& in);
	double get_cp(double T = 0) { return cp; }
	double get_lambda(double T = 0) { return lambda; }
	double get_rho(double p, double T);
	double get_mu(double p, double T);
	double get_K(double p, double T);

};

# endif