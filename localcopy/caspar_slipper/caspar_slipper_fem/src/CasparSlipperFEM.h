#include "fem.h"
#include "../../caspar_input/input.h"

namespace CasparSlipperFEM
{
	class thermoelastic
	{
	public:
		//inputs
		string option_file;
		string vtk_file;
		vector<double> heat_flux;
		const input_data * caspar_gap_input;

		//outputs
		vector<double> node_temp;
		vector<vector<double> > node_deform;

		//methods
		void run_analysis(const string heatfluxsetname, const double alpha = 1, const vector<double> & old_temp = vector<double> ());
		fem * open();
		void get_faceset(fem * body, const string name, vector<int> & nid, vector<vector<double> > & nodexyz, vector<vector<int> > & setfaces);
		unsigned int get_nodecnt(fem * body);
		void close(fem * body);

		thermoelastic()
		{
			caspar_gap_input = NULL;
		}
	};

}

