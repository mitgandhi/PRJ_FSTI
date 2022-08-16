#include "main.h"
#include "point.h"

namespace CasparSlipperFEM
{
	struct node
	{
		//id - must be the same int value as the vector index in the fem class
		int id;

		//XYZ location
		point coord;

		//global DOF numbering - size will depend on analysis type
		vector<int> DOF;

		node()
		{
			coord = point(0,0,0);
		};

		node(const point p)
		{
			coord = p;
		};

		node(const double x,const double y, const double z)
		{
			coord = point(x,y,z);
		};

		double operator[] (const int i) 
		{
			return coord[i];
		};

		double x()
		{
			return coord.x();
		};

		double y()
		{
			return coord.y();
		};

		double z()
		{
			return coord.z();
		};
	};
};
