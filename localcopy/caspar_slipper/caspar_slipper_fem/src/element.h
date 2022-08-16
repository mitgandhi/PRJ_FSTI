#include "main.h"
#include "dmatrix.h"
#include "node.h"
#include "face.h"

#pragma once

namespace CasparSlipperFEM
{
	class element
	{
	protected:
		//These are values from the assigned from the boundary conditions
		dmatrix<double> Kbc;
		vector<double> Fbc;

	public:
		//Element type
		enum element_types { TETRA };
		const element_types element_type;

		analysis_types analysis_type;

		//Pointer to nodes
		vector<node *> nodes;

		//Constant properities
		double E;		//Elastic modulus
		double nu;		//Possion ratio
		double rho;		//Density
		double k;		//Thermal conductivity
		double alpha;	//Thermal expansion

		//Element volume
		virtual double V() = 0;

		//local DOF -> global DOF map
		virtual vector<int> l2gDOF() = 0;

		//Return the niddof of a stiffness matrix entry
		virtual pair<int, int> niddof(const int ki) = 0;

		//Stiffness matrix
		virtual dmatrix<double> K() = 0;

		//B vector
		virtual vector<double> F() = 0;

		//Return a vector of element faces
		virtual vector< face > getfaces() = 0;

		//Assign analysis type
		virtual void assign_analysis(const analysis_types type) = 0;

		//Assign face flux
		virtual void assign_face_nbc(const face & f, const double val) = 0;
		
		//Assign face mbc
		virtual void assign_face_mbc(const face & f, const double h, const double phi) = 0;

		//Assign ele temp
		virtual void assign_ele_therm(const double dt) = 0;

		//Get element center
		virtual point get_center() = 0;
		
		//Constructor
		element(const element_types e) : element_type(e)
		{
		};

	};
};