# include <fstream>
# include <sstream>
# include <string>

# include "./p_file.h"
# include "../FSTI_Block_dll/log.h"

using namespace std;

extern gaplog Log;

// ------------------------------------------------------------------------- //
p_file::p_file()
{
}
// ------------------------------------------------------------------------- //
void p_file::initialize(const char* filename, double _speed)
{

	// update the speed
	speed = _speed;

	vector<double> _time(0);
	vector<double> _pDC(0);
	vector<double> _pHP(0);
	vector<double> _pLP(0);

	ifstream in(filename);	
	// check for correct opening
	if(!in.is_open()) 
	{
		Log << gaplog::endl << "p_file::p_file: Unable to open DC pressure input file!" 
				 << gaplog::endl;
		exit(1);
	}
	
	Log << "Reading pressure file " << filename << gaplog::endl << gaplog::endl;
	
	string line;
	double read;

	while(getline(in,line))
	{
		if(line.size() > 0)
		{
			istringstream iss (line,istringstream::in);
			vector<double> colData(0);
			while(iss) 
			{
				string tmp;
				iss >> tmp;
				if(tmp.size() > 0)
				{
					// read all data inm the line using a istringstream
					istringstream(tmp) >> read;
					colData.push_back(read);
				}
			}
		
			// store the reading
			if(colData.size() == 4)
			{
				// check that the values of time are monotonically increasing
				if(_time.size() > 0)
				{
					if(colData[0] > _time.back())
					{
						_time.push_back(colData[0]);
						_pDC.push_back(colData[1]);
						_pHP.push_back(colData[2]);
						_pLP.push_back(colData[3]);
					}
				}
				else
				{
					_time.push_back(colData[0]);
					_pDC.push_back(colData[1]);
					_pHP.push_back(colData[2]);
					_pLP.push_back(colData[3]);
				}
			}
			else
			{
				Log << "\n\nError reading " << filename << ": found " << colData.size() 
						 << " columns" << gaplog::endl << gaplog::endl;
				exit(1);
			}
		}
		
	}
	
	in.close();

	double dt = _time[_time.size()-1] - _time[_time.size()-2];

	// check that the time starts from zero, if not add a line to avoid interpolation error
	if(_time[0] > 0.0)
	{
		// insert one line
		_time.insert(_time.begin(), 0.0);
		_pDC.insert(_pDC.begin(), _pDC[0]);
		_pHP.insert(_pHP.begin(), _pHP[0]);
		_pLP.insert(_pLP.begin(), _pLP[0]);
		
	}

	// check if the end time is correct
	double T = 2.0*pi()/speed;
	int counte = 0;
	// add a line in the end to make sure the interpolation goes smooth at 2pi
	
	while(_time.back() > T)	// check the consistency of the pressure file
	{
		Log << "\nNumber of Rows before: "<<  _time.size()
					<< gaplog::endl;	
				Log << "\nNumber of Rows before: "<<  _time.size()
			<< gaplog::endl;	
		// DELETE EXTRA ROWS IN P_FILE
		_time.erase (_time.end()-1);
		_pDC.erase (_pDC.end()-1);
		_pHP.erase (_pHP.end()-1);
		_pLP.erase (_pLP.end()-1);
		counte++;
		Log << "\nNumber of Rows after: "<<  _time.size()
			<< gaplog::endl;	
		if (counte > 10)
		{
			Log << "\nERROR: the pressure file doesn't seem to match with the declared speed." 
				<< gaplog::endl;
			exit(1);
		}
	}
	if(_time.back() < T)
	{
		_time.insert(_time.end(), T);
		_pDC.insert(_pDC.end(), _pDC[_pDC.size()-1]);
		_pHP.insert(_pHP.end(), _pHP[_pHP.size()-1]);
		_pLP.insert(_pLP.end(), _pLP[_pLP.size()-1]);
	}


	sz = _time.size();
	
	double* time = new double[sz];
	double* phi = new double[sz];
	double* pHP = new double[sz];
	double* pLP = new double[sz];
	double* pDC = new double[sz];

	php_avg = 0;
	plp_avg = 0;

	
	for(int i=0; i<sz; i++)
	{
		time[i] = _time[i];
		phi[i] = speed*_time[i];
		pDC[i] = _pDC[i];
		pHP[i] = _pHP[i];
		pLP[i] = _pLP[i];
		php_avg += pHP[i];
		plp_avg += pLP[i];
	}

	

	php_avg = (php_avg/sz)/1e5;
	plp_avg = (plp_avg/sz)/1e5;

	acc = gsl_interp_accel_alloc ();
	
	// allocate storage
	pLP_s = gsl_spline_alloc (gsl_interp_cspline, sz);
	pHP_s = gsl_spline_alloc (gsl_interp_cspline, sz);
	pDC_s = gsl_spline_alloc (gsl_interp_cspline, sz);
	
	// create splines
  gsl_spline_init (pLP_s, phi, pLP, sz);
	gsl_spline_init (pHP_s, phi, pHP, sz);
	gsl_spline_init (pDC_s, phi, pDC, sz);

	delete[] time;
	delete[] pLP;
	delete[] pHP;
	delete[] pDC;
	delete[] phi;

	Log << "  * read data in the time interval [" 
			 << _time[0] << "," << _time[_time.size()-1] << "]" << gaplog::endl;

	Log << "  * mean value for HP: " << php_avg << " bar" << gaplog::endl;
	Log << "  * mean value for LP: " << plp_avg << " bar" << gaplog::endl;

	Log << gaplog::endl << "pressure file reading completed!" << gaplog::endl;
}
// ------------------------------------------------------------------------- //
p_file::~p_file()
{
	gsl_interp_accel_free (acc);
	gsl_spline_free(pLP_s);
	gsl_spline_free(pHP_s);
	gsl_spline_free(pDC_s);
}
// ------------------------------------------------------------------------- //
vector<double> p_file::getp_rad(double phi_rad) const
{
	
	double pDC = gsl_spline_eval (pDC_s, phi_rad, acc);
	double pLP = gsl_spline_eval (pLP_s, phi_rad, acc);
	double pHP = gsl_spline_eval (pHP_s, phi_rad, acc);

	vector<double> tmp(3);
	tmp[0] = pDC;
	tmp[1] = pLP;
	tmp[2] = pHP;

	return tmp;
}
// ------------------------------------------------------------------------- //
vector<double> p_file::getp_deg(double phi_deg) const
{
	double phi_rad = phi_deg*pi()/180.0;

	double pDC = gsl_spline_eval (pDC_s, phi_rad, acc);
	double pLP = gsl_spline_eval (pLP_s, phi_rad, acc);
	double pHP = gsl_spline_eval (pHP_s, phi_rad, acc);

	vector<double> tmp(3);
	tmp[0] = pDC;
	tmp[1] = pLP;
	tmp[2] = pHP;

	return tmp;
}
// ------------------------------------------------------------------------- //
