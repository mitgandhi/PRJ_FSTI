#include "element.h"

namespace CasparSlipperFEM
{
	class tetra_element: public element
	{

	private:
		double x(const int n);
		double y(const int n);
		double z(const int n);
		
		//Determinant of J
		double detJ();

		//Elasticity
		dmatrix<double> Be();
		dmatrix<double> De();
		dmatrix<double> Ke();

		//Thermal
		dmatrix<double> Bt();
		dmatrix<double> Dt();
		dmatrix<double> Kt();

		//A first order tetra has 4 nodes
		static const int nnodes = 4;
		
		//How many degree of freedom does each node have (depends on analysis type)
		int nodeDOF;

	public:
		tetra_element();
		tetra_element(node * n0, node * n1, node * n2, node * n3);
		
		//Volume
		double V();

		//Local 2 Global DOF map
		vector<int> l2gDOF();

		//Return the niddof of a stiffness matrix entry
		pair<int, int> niddof(const int ki);

		//Stiffness matrix
		dmatrix<double> K();

		//Load vector
		vector<double> F();

		//Assign analysis type
		void assign_analysis(const analysis_types type);

		//Return a vector of element faces
		vector<face> getfaces();

		//Assign face flux
		void assign_face_nbc(const face & f, const double val);

		//Assign face mbc
		void assign_face_mbc(const face & f, const double h, const double phi);

		//Assign a element thermal load
		void assign_ele_therm(const double dt);

		//Get element center
		point get_center();

	};
};
