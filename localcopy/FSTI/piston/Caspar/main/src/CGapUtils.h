#include "CGapInput.h"
#pragma once


class CGapUtils
{
	public:

	//Constructor destructor
	CGapUtils(void);
	~CGapUtils(void);

	//search nodes neighbours
	void SearchNeighbours(Array<double,2> xyz_data,Array<double,2> xyz_query,Array<int,2> &NodeId,Array<double,2> &NodeDist,int nn);
	//interpolate fields based on Shepard's interpolation
	Array<double,1> InterpolateFields(Array<int,2> NodeId,Array<double,2> NodeDist,Array<double,1> field_data);
	//calculate KD tree for the field
	void CalcKDtree(int n,Array<double,2> data,Array<double,2> query,Array<int,2> &NodeId,Array<double,2> &NodeDist);

};

