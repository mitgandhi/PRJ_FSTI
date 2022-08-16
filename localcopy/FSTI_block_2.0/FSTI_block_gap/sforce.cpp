# include "./sforce.h"
# include "../FSTI_Block_dll/log.h"
# include <fstream>
# include <string>
# include <sstream>

extern gaplog Log;

using namespace std;

const double pi = 4.0*atan(1.0);

// ------------------------------------------------------------------------- //
sforce::sforce()
{
	acc = 0;
	f = 0;
	last_complete_rev = -1;
}
// ------------------------------------------------------------------------- //
sforce::~sforce()
{
	if(acc != 0)
		gsl_interp_accel_free (acc);
	if(f != 0)
		gsl_spline_free(f);
}
// ------------------------------------------------------------------------- //
void sforce::initialize(const char* name)
{
	filename = string(name);
}
// ------------------------------------------------------------------------- //
int sforce::read()
{
	
	ifstream in(filename.c_str());

	if(!in.is_open())
		return -1;

	string line;

	data.t.resize(0);
	data.phi.resize(0);
	data.val.resize(0);

	while(getline(in,line))
	{
		istringstream iss(line);
		double tmp;
		iss >> tmp;
		data.t.push_back(tmp);
		iss >> tmp;
		data.phi.push_back(tmp);
		iss >> tmp;
		data.val.push_back(tmp);
	}

	in.close();

	return 0;
}
// ------------------------------------------------------------------------- //
void sforce::update(double trev)
{

	// read the file
	int was_read = read();

	if(was_read == -1)
	{
		//Log << "\nForce definition file " << filename << " could not be read\n"
		//		<< gaplog::endl;
		
		// set to zero the spline pointer (will not be used)
		acc = 0;
		f = 0;

		return;
	}

	int last = -1;

	if(data.t.size() > 0)
	{
		// number of complete revolutions
		int nrev = static_cast<int>(floor(data.t.back()/trev));
		last = nrev - 1;
	}

	if(last >= 0 && last_complete_rev < last)
	{
		// update the last complete revolution index
		last_complete_rev = last; 

		// ------------------- update the spline definition -------------------- //

		Log << "\nUpdating force definition associated with file " 
				 << filename << " (rev. " << last_complete_rev << ") ... ";

		vector<double> t(0), phi(0), val(0);

		// start time
		double tstart = trev*last_complete_rev;
		// end time
		double tend = trev*(last_complete_rev + 1);

		// save just the last complete revolution
		for(unsigned int i=0; i<data.t.size(); i++)
		{
			if(data.t[i] >= tstart && data.t[i] <= tend)
			{
				t.push_back(data.t[i]);
				phi.push_back(2.0*pi*(data.t[i]-tstart)/trev);
				val.push_back(data.val[i]);
			}
		}

		// make sure it starts from zero and ends at 2pi
		if(phi[0] > 0)
		{
			phi.insert(phi.begin(), 0.0);
			val.insert(val.begin(), val[0]);
		}
		if(phi.back() < 2.0*pi)
		{
			phi.insert(phi.end(), 2.0*pi);
			val.insert(val.end(), val.back());
		}

		// define size and copy into the double[] containers
		int sz = phi.size();

		double* _phi = new double[sz];
		double* _val = new double[sz];

		for(int i=0; i<sz; i++)
		{
			_phi[i] = phi[i];
			_val[i] = val[i];
		}

		// clear the spline and the interpolation object
		if(f != 0)
			gsl_spline_free(f);
		if(acc != 0)
			gsl_interp_accel_free(acc);

		// allocate storage
		acc = gsl_interp_accel_alloc ();
		f = gsl_spline_alloc (gsl_interp_cspline, sz);
		// create the spline
		gsl_spline_init (f, _phi, _val, sz);

		// delete temporary objects
		delete [] _phi;
		delete [] _val;

		Log << "done!" << gaplog::endl;

	}
	
}
// ------------------------------------------------------------------------- //
double sforce::getf(double phirad)
{
	double F = 0;

	if(f == 0)
	{
		return F;	// if the force file could not be read, return 0
	}
	else
	{
		if(phirad >= 0.0 && phirad <= 2.0*pi)
		{
			if(f != 0)
			{
				F = gsl_spline_eval (f, phirad, acc);
			}
		}
		else
		{
			Log << "\nsforce::getf: Angle " << phirad << " is not within [0,2pi]\n" 
					 << gaplog::endl;
			exit(1);
		}

		return F;
	}
}
// ------------------------------------------------------------------------- //