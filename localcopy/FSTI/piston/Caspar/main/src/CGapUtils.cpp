#include "CGapUtils.h"
#pragma once

CGapUtils::CGapUtils(void)
{


}
CGapUtils::~CGapUtils(void)	//Destructor
{

}; 
//Search node neighbours
void CGapUtils::SearchNeighbours(Array<double,2> xyz_data,Array<double,2> xyz_query,Array<int,2> &NodeId,Array<double,2> &NodeDist, int nn)
{

	int nPts_q,n;
	Range all = Range::all();

	//Number of points query mesh
	nPts_q = xyz_query.extent(0);
	
	//Search neighbouring nodes with KD-tree
	n = nn;
	NodeId.resize(nPts_q,n);	NodeId=0;
	NodeDist.resize(nPts_q,n);	NodeDist=0.0;
	CalcKDtree(n,xyz_data,xyz_query,NodeId,NodeDist);
	
};
//Interpolate fields based on weighted distance
Array<double,1> CGapUtils::InterpolateFields(Array<int,2> NodeId,Array<double,2> NodeDist,Array<double,1> field_data)
{
	int nPts_q,n;
	Range all = Range::all();

	//Number of points fine mesh
	nPts_q = NodeId.extent(0);
	n = NodeId.extent(1);

	//Shepard's inverse distance weighted interpolation
	Array<double,1> field_query(nPts_q);
	field_query = 0.0;
	double wi,dn;
	for(int i=0;i<nPts_q;i++)
	{
		//Simple Shepard's method
		dn = 0.0;
		wi = 0.0;
		for(int j=0;j<n;j++)
		{
			wi = 1.0 / NodeDist(i,j) ;
			field_query(i) += field_data( NodeId(i,j) ) * wi ;
			dn += wi ;
		}
		field_query(i) /= dn;
	};

	return field_query;

};
//KD tree fuction to determine nodes neighbours 
void CGapUtils::CalcKDtree(int n,Array<double,2> data,Array<double,2> query,Array<int,2> &NodeId,Array<double,2> &NodeDist)
{
	int nPts_d,nPts_q,dim;
	double tol;
	Array<int,2> nodeId;
	ANNpointArray dataPts;	//Data points
	ANNpoint queryPt;		//Query point
	ANNidxArray nnIdx;		//Near neighbor indices
	ANNdistArray dists;		//Near neighbor distances
	ANNkd_tree* kdTree;		//Search structure


	nPts_d = data.extent(0);
	nPts_q = query.extent(0);
	dim = 3;
	tol = 0.0;
	//Allocate data points
	dataPts = annAllocPts(nPts_d,dim);
	//Allocate query pt
	queryPt = annAllocPt(dim); 

	//cout << "Creating dataPts..." << "\n";
	for (int i=0;i<nPts_d;i++)
	{
		dataPts[i][0] = data(i,0);
		dataPts[i][1] = data(i,1);
		dataPts[i][2] = data(i,2);
	}

	//cout << "Creating kd tree..." << "\n";
	kdTree = new ANNkd_tree(dataPts,nPts_d,dim);
	nnIdx = new ANNidx[n];						// allocate near neigh indices
	dists = new ANNdist[n];						// allocate near neighbor dists

	//cout << "Searching kd tree..." << "\n";
	int nNodes=0;
	NodeId=0;
	NodeDist=0;
	for(int i=0;i<nPts_q;i++)
	{
		queryPt[0] = query(i,0);
		queryPt[1] = query(i,1);
		queryPt[2] = query(i,2);
		//Tree search
		kdTree->annkSearch(queryPt,n,nnIdx,dists);
		for(int j=0;j<n;j++)
		{
			//Assign to node id array
			NodeId(i,j)=nnIdx[j];
			//Assign the node distances squared
			NodeDist(i,j)=sqrt(dists[j]);
		}
		//new node
		nNodes++;

	};
	NodeDist = where(NodeDist<1.0e-6,1.0e-6,NodeDist);

	annDeallocPt(queryPt);
	annDeallocPts(dataPts);
	delete [] nnIdx;
	delete [] dists;
	delete kdTree;
	annClose();
	//done with ANN

};