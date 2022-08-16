# ifndef __sforce__
# define __sforce__

# include <vector>
# include <gsl/gsl_errno.h>
# include <gsl/gsl_spline.h>

struct sforce
{
	// file name
	std::string filename;

	gsl_interp_accel* acc;
	gsl_spline* f;	// spline 

	int last_complete_rev;

	struct
	{
		bool available;
		std::vector<double> t;
		std::vector<double> phi;
		std::vector<double> val;
	} data;

	sforce();
	~sforce();
	void initialize(const char* name);
	int read();									// read from input file
	void update(double trev);					// update the spline, input is the time for one revolution
	double getf(double phirad);					// update the reading, input is the angle in rad
	
};

# endif