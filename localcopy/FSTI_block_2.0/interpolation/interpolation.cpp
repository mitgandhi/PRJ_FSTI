# include "./interpolation.h"
# include "../interpolation/polygon.h"
# include "../FSTI_Block_dll/log.h"

# include <iostream>

using namespace std;

extern class gaplog Log;

// ------------------------------------------------------------------------- //
interpolation::interpolation()
{
	ANN.eps = 0.0;
	ANN.dataPts = 0;
	ANN.dists = 0;
	ANN.nnIdx = 0;
	ANN.kdTree = 0;
	cells2cells.ready = false;
	points2cells.ready = false;
}
// ------------------------------------------------------------------------- //
interpolation::~interpolation()
{
	clear();
	annClose();
}
// ------------------------------------------------------------------------- //
void interpolation::clear()
{
	
	if(ANN.dataPts != 0)
	{
		annDeallocPts(ANN.dataPts);
		ANN.dataPts = 0;
	}
	if(ANN.nnIdx != 0)
	{
		delete [] ANN.nnIdx;
		ANN.nnIdx = 0;
	}
	if(ANN.dists != 0)
	{
		delete [] ANN.dists;
		ANN.dists = 0;
	}
	if(ANN.kdTree != 0)
	{
		delete ANN.kdTree;
		ANN.kdTree = 0;
	}

	cells2cells.ready = false;
	cells2cells.elms_idx.resize(0);
	cells2cells.elms_weight.resize(0);
	cells2cells.elms_weights.resize(0);

	points2cells.ready = false;
	points2cells.elm_idx.resize(0);
	points2cells.centers.resize(0);
	
}
// ------------------------------------------------------------------------- //
void interpolation::buildTree(const interpl_grid* from)
{
	if(ANN.dataPts != 0)
	{
		annDeallocPts(ANN.dataPts);
		ANN.dataPts = 0;
	}
	if(ANN.nnIdx != 0)
	{
		delete [] ANN.nnIdx;
		ANN.nnIdx = 0;
	}
	if(ANN.dists != 0)
	{
		delete [] ANN.dists;
		ANN.dists = 0;
	}
	if(ANN.kdTree != 0)
	{
		delete ANN.kdTree;
		ANN.kdTree = 0;
	}

	ANN.queryPt = annAllocPt(ANN.dim);
	ANN.dataPts = annAllocPts(from->ne(), ANN.dim);
	// define elements center
	for(int i=0; i<from->ne(); i++)
	{
		point c(0,0);
		for(int j=0; j<from->ELM_NDS; j++)
			c = c + (1.0/from->ELM_NDS)*from->nodes[from->elements[i][j]];
		
		ANN.dataPts[i][0] = c.x();
		ANN.dataPts[i][1] = c.y();
	}

	// build search structure
	ANN.kdTree = new ANNkd_tree
	( 
		ANN.dataPts,	// the data points
		from->ne(),		// number of points
		ANN.dim
	);

}
// ------------------------------------------------------------------------- //
void interpolation::define_cells2cells
(
	const interpl_grid* from, 
	const interpl_grid* to
)
{
	
	buildTree(from);

	// per each body element 
	const int max_elms = 5 + 5*from->ne()/to->ne();

	cells2cells.elms_idx.resize(to->ne());
	cells2cells.elms_weights.resize(to->ne());
	cells2cells.elms_weight.resize(to->ne());

	for(int i=0; i<to->ne(); i++)
	{
		if(to->active[i])
		{

			// define polygon using a "to" interpl_grid element
			polygon p_to(to->ELM_NDS);
			for(int h=0; h<to->ELM_NDS; h++)
				p_to.vertices[h] = to->nodes[to->elements[i][h]];
		
			point c = p_to.center();
			
			double p_to_area = p_to.area();
			ANN.queryPt[0] = c.x();
			ANN.queryPt[1] = c.y();

			// flag used to stop the searching process 
			// The searching process stops when the sum of the intersection areas,
			// area_i will be approximately equal to the area of the triangle.
			// Since in the interpl_grid there are openings boudnaries (and this 
			// boundaries could not match exactly in the two grids) the searching 
			// process may be stopped is more than max_elements are involved.
			bool completed = false;

			// sum of area_i
			cells2cells.elms_weight[i] = 0.0;

			int k = 1;

			if(ANN.nnIdx != 0)
				delete [] ANN.nnIdx;
			if(ANN.dists != 0)
				delete [] ANN.dists;

			// searching loop
			do
			{

				ANN.nnIdx = new ANNidx[k];
				ANN.dists = new ANNdist[k];
				
				// searching for the k closest fluid element
				ANN.kdTree->annkSearch
				(
					ANN.queryPt, 
					k, 
					ANN.nnIdx, 
					ANN.dists, 
					ANN.eps
				);

				if(from->active[ANN.nnIdx[k-1]])
				{
					polygon p_from(from->ELM_NDS);

					for(int h=0; h<from->ELM_NDS; h++)
						p_from.vertices[h] = from->nodes[from->elements[ANN.nnIdx[k-1]][h]];

					std::vector<polygon> inters = p_to.clip(p_from);	// clip defaults to intersection
					polygon I;
					if(inters.size() > 0)
					{
						I = inters[0]; 
					}
				
					// assign the intersection area to the area_i vector, only if the value
					// is bigger than zero
					double thisarea = I.area();
					if(thisarea > 0)
					{
						cells2cells.elms_weight[i] += thisarea;
						cells2cells.elms_idx[i].push_back(ANN.nnIdx[k-1]);
						cells2cells.elms_weights[i].push_back(thisarea);
					}

					// use 1% tolerance
					if
					(
						fabs(cells2cells.elms_weight[i]  -  p_to_area) < 0.01*p_to_area 
						|| k > max_elms
					)
					{
						completed = true;
					}
				}

				delete [] ANN.nnIdx;
				delete [] ANN.dists;
				ANN.nnIdx = 0;
				ANN.dists = 0;
				
				k++;

			} while (!completed);
		}
	}

	cells2cells.ready = true;
}
// ------------------------------------------------------------------------- //
void interpolation::cellsTocells
(
	const interpl_grid* from, 
	interpl_grid *to
)
{
	
	if(cells2cells.ready == false)
	{
		define_cells2cells(from, to);
	}

	for(int i=0; i<to->ne(); i++)
	{
		if(to->active[i])
		{
			int nweights = cells2cells.elms_weights[i].size();
			double sum_weight = cells2cells.elms_weight[i];
			
			to->cells_data[i] = 0;
			for(int j=0; j<nweights; j++)
			{
				double weight = cells2cells.elms_weights[i][j];
				to->cells_data[i] += 
					from->cells_data[cells2cells.elms_idx[i][j]]*(weight/sum_weight);
			}
		}
	}

}
// ------------------------------------------------------------------------- //
void interpolation::define_points2cells
(
	const interpl_grid* from, 
	const interpl_grid* to
)
{

	if(from->ELM_NDS != 3)
	{
		Log << "\ninterpolation::define_points2cells: ERROR interpl_grid \"from\" must be " 
				 << "a triangular interpl_grid!" << gaplog::endl;
		exit(1);
	}
	
	// build the kdtree
	buildTree(from);

	points2cells.elm_idx.resize(to->ne());
	points2cells.centers.resize(to->ne(),point(0,0));

	for(int i=0; i<to->ne(); i++)
	{

		if(to->active[i])
		{			 
			for(int j=0; j<to->ELM_NDS; j++)
			{
				points2cells.centers[i] = 
					points2cells.centers[i] + (1.0/to->ELM_NDS)*to->nodes[to->elements[i][j]];
			}

			// find the triangle 
			bool found = false;
			int k = 1;
			int attempts = 0;
			ANN.queryPt[0] = points2cells.centers[i].x();
			ANN.queryPt[1] = points2cells.centers[i].y();
			
			if(ANN.nnIdx != 0)
				delete [] ANN.nnIdx;
			if(ANN.dists != 0)
				delete [] ANN.dists;

			// searching loop
			do
			{
				ANN.nnIdx = new ANNidx[k];
				ANN.dists = new ANNdist[k];
				
				// searching for the k closest fluid element
				ANN.kdTree->annkSearch
				(
					ANN.queryPt, 
					k, 
					ANN.nnIdx, 
					ANN.dists, 
					ANN.eps
				);

				if(from->active[ANN.nnIdx[k-1]])
				{
					if(from->isInside(ANN.nnIdx[k-1], points2cells.centers[i]))
					{
						points2cells.elm_idx[i] = ANN.nnIdx[k-1];
						found = true;
					}
					else
					{
						attempts++;
					}
				}

				if(attempts > 10)
				{	
					//Log << "\nWarning: fluid element with center (" 
					//		 << points2cells.centers[i].x() << ", " 
					//		 << points2cells.centers[i].y() << ")" << gaplog::endl
					//		 << "         seems to be outside any triangle." << gaplog::endl
					//		 << "         Interpolation was accomplished using the triangle" << gaplog::endl
					//		 << "         with center in ("
					//		 << ANN.dataPts[ANN.nnIdx[0]][0] 
					//		 << ", " << ANN.dataPts[ANN.nnIdx[0]][1]  << ")" << gaplog::endl;

					points2cells.elm_idx[i] = ANN.nnIdx[0];

					break;
				}

				delete [] ANN.nnIdx;
				delete [] ANN.dists;
				ANN.nnIdx = 0;
				ANN.dists = 0;
				k++;

			} while (!found);

		}
	}

	points2cells.ready = true;
}
// ------------------------------------------------------------------------- //
void interpolation::pointsTocells
(
	const interpl_grid* from, 
	interpl_grid* to
)
{
	if(points2cells.ready == false)
	{
		define_points2cells(from, to);
	}

	for(int i=0; i<to->ne(); i++)
	{
		if(to->active[i])
		{

			double x0 = from->nodes[from->elements[points2cells.elm_idx[i]][0]].x();
			double x1 = from->nodes[from->elements[points2cells.elm_idx[i]][1]].x();
			double x2 = from->nodes[from->elements[points2cells.elm_idx[i]][2]].x();
			double y0 = from->nodes[from->elements[points2cells.elm_idx[i]][0]].y();
			double y1 = from->nodes[from->elements[points2cells.elm_idx[i]][1]].y();
			double y2 = from->nodes[from->elements[points2cells.elm_idx[i]][2]].y();

			double v0 = from->points_data[from->elements[points2cells.elm_idx[i]][0]];
			double v1 = from->points_data[from->elements[points2cells.elm_idx[i]][1]];
			double v2 = from->points_data[from->elements[points2cells.elm_idx[i]][2]];

			// linear interpolation from triangle nodes to the fluid point (x,y) 
			// by solving the system:

			//	a*x0 + b*y0 + c = v0;
			//	a*x1 + b*y1 + c = v1;
			//	a*x2 + b*y2 + c = v2;

			// v0,v1,v2 displacements at the nodes of the triangle
			// (x0,y0),(x1,y1),(x2,y2) triangle coordinates
			// displacement(x,y) = ax + by + c

			// coefficients calculation (expressions found with the Matlab symbolic toolbox)

			double a = 
				(v0*(y1 - y2))/(x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1) 
				- (v1*(y0 - y2))/(x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1) 
				+ (v2*(y0 - y1))/(x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1);
		    
			double b = 
				(v1*(x0 - x2))/(x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1) 
				- (v0*(x1 - x2))/(x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1) 
				- (v2*(x0 - x1))/(x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1);

			double c = 
				(v2*(x0*y1 - x1*y0))/(x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1) 
				- (v1*(x0*y2 - x2*y0))/(x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1) 
				+ (v0*(x1*y2 - x2*y1))/(x0*y1 - x1*y0 - x0*y2 + x2*y0 + x1*y2 - x2*y1);


			to->cells_data[i] = 
				a*points2cells.centers[i].x() + b*points2cells.centers[i].y() + c;

		}
	}
}
// ------------------------------------------------------------------------- //
