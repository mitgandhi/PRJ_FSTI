#include "main.h"
#include <gmm/gmm.h>
#include "tetra_element.h"
#include "femoptions.h"
#include "../../caspar_input/input.h"

namespace CasparSlipperFEM
{
	class fem
	{
	public:

		const input_data * caspar_gap_input;

		//fem options
		femoptions options;

		//protected class members
		//protected:

		//Element set definition
		typedef vector<int> elementset;
		
		//Faceset definition
		typedef vector<face> faceset;

		//Nodeset definition
		typedef vector<int> nodeset;

		//What type of analysis is this
		analysis_types analysis_type;

		//Nodes
		vector<node> nodes;
		int nodecnt;

		//Elements
		vector<element *> elements;
		int elecnt;

		//Build a map of element faces - this is time consuming to build and uses 100's MB memory .. is there another way?
		map<face, vector<int> > elementfaces;

		//Element sets
		map<string, elementset> elementsets;

		//Face sets
		map<string, faceset> facesets;

		//Node sets
		map<string, nodeset> nodesets;
		
		//A simple pair of the nid and dof used to define BC's
		typedef pair<int, int> niddof;

		//Essential boundary conditions (Dirichlet)
		map<niddof, double> ebcs;

		//Rigid body elements
		map<niddof, niddof> rbe;

		//stiffness matrix
		gmm::col_matrix<gmm::wsvector<double> > K;

		//load vector
		vector<double> b;

		//X vector (u or phi) size = uDOF
		vector<double> X;

		//Number of unconstrained DOF
		int uDOF;	

		//Number of DOF
		int DOF;

		//Utility methods
		void numberDOF();
		void setLoads();
		void clearElements();
		void checkAnalysis();
		bool checkNodeset(const string ns);
		bool checkFaceset(const string fs);
		bool checkElementset(const string es);
		void warning(const string w);
		void error(const string e);
		void apply_inertia_relief(const double avgK);

		//public class methods
	public:
		//Member methods - general
		fem();
		~fem();
		void reset();

		//Member methods - IO
		void readoptions(string file);
		void loadinp();
		void writeK(const string filename);
		void writeB(const string filename);
		void writeVTK(const string filename, const vector<double> & nodetmp = vector<double>(),
											const vector<double> & dirchlet = vector<double>(),
											const vector<double> & neumann = vector<double>(),
											const vector<double> & mixed_h = vector<double>(),
											const vector<double> & mixed_phi = vector<double>() );
		void writeFacesetVTK(const string filename, const string name);											

		//Member methods - props
		void assign_mats();
		void assign_analysis(analysis_types type);

		//primarly used internally
		void assign_E(const string elset, const double E);
		void assign_nu(const string elset, const double nu);
		void assign_rho(const string elset, const double rho);
		void assign_k(const string elset, const double k);
		void assign_alpha(const string elset, const double alpha);
		
		//Member methods - BCs
		bool inertia_relief;
		void assign_bcs();
		void get_thermal_bcs(vector<double> &dirchlet, vector<double> &neumann, vector<double> &mixed_h, vector<double> &mixed_phi, const string heatfluxsetname, const fem::faceset &fullheatfluxset, const vector<double> &heat_flux_vals);

		//primarly used internally
		void apply_ebc(const string ns, const int dof, const double val);
		void apply_faceset_nbc(const string fs, const double val);
		void apply_faceset_nbc(const string fs, vector<double> vals);
		void apply_faceset_mbc(const string fs, const double h, const double phi);
		void apply_ele_therm(const int eid, const double t);

		//Member methods - solving
		void stiffnessMatrix();
		bool solve(int SolSequence = 0);

		//Output
		int get_elecnt();
		int get_nodecnt();
		double getele_temp(const int eid);
		double getnode_temp(const int nid);
		vector<double> getnode_deform(const int nid);

	};
};
