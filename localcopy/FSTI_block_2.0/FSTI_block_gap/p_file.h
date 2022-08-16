# ifndef __CpFile__
# define __CpFile__

# include <gsl/gsl_errno.h>
# include <gsl/gsl_spline.h>
# include <vector>

class p_file
{
	
	static const int max_sz = 10000;
	int sz;							// file data size
	double speed;				// pump/motor speed
	
	double php_avg;			// average high pressure
	double plp_avg;			// average low pressure
	
	gsl_interp_accel* acc;	// used for interpolation
	gsl_spline* pDC_s;			// spline fr DC pressure
	gsl_spline* pLP_s;			// spline for low pressure
	gsl_spline* pHP_s;			// spline for high pressure
	
	double pi() const { return 4.0*atan(1.0); } 

public:

	p_file();																							// default constructor
	~p_file();																						// destructor
	void initialize(const char* fame, double speed);			// initialize the pressure file
	double pHP_avg() const { return php_avg; }						// get average HP
	double pLP_avg() const { return plp_avg; }						// get average LP
	std::vector<double> getp_rad(double phi_rad) const;		// get [pDC, pLP, pHP] given as input an angle in radiants
	std::vector<double> getp_deg(double phi_deg) const;		// get [pDC, pLP, pHP] given as input an angle in deg

};

# endif