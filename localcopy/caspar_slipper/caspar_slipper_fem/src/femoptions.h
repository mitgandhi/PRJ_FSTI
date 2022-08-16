#include <string>

namespace CasparSlipperFEM
{
	struct femoptions
	{
		struct material
		{
			string name;
			double E;
			double nu;
			double alpha;
			double lambda;
			double rho;
			vector<string> sets;
		};

		struct thermal_boundary
		{
			enum types {DIRICHLET, NEUMANN, ROBBINS};
			types type;
			string set;
			double Tp;
			double Tinf;
			double h;
			double q;

			thermal_boundary()
			{
				set = ""; //just to make sure
			}
		};

		struct elastic_constraint
		{
			string set;
			double x;
			double y;
			double z;
		};

		femoptions()
		{
			//default inertia relief to 0
			inrel = 0;

			//default refTemp to 25
			reftemp = 25;

			//default scale factor to mm
			scalefactor = 1.0e-3;
		}

		//general options
		string casename;
		string meshfile;
		double scalefactor;
		double reftemp;

		//boundary conditions
		vector<material> materials;
		vector<thermal_boundary> thermal_boundries;
		vector<elastic_constraint> elastic_constraints;
		int inrel;

		//solver options
		int maxIters;
		double tolerance;

	};

};
