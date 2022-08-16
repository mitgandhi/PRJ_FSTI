# include <algorithm>
# include <sstream>
# include "./gap_boundary.h"
# include "../FSTI_Block_dll/log.h"

# define pi 3.14159265358979323846

using namespace std;

extern class gaplog Log;

// ------------------------------------------------------------------------- //
bool gap_boundary::is_top_node(const point& node) const
{
	if(fabs(node.z() - top) < 0.1*(top - bottom))
		return true;
	else
		return false;
}
// ------------------------------------------------------------------------- //
void gap_boundary::select_triangles()
{
	const double ang_tol = 30.0;

	point dir(0,0);
	point c(0,0);

	top = bottom = 0;
	int nmax = 0;
	int nmin = 0;

	vector<int> side_triangles(0);

	// determine bottom and top
	for(unsigned int i=0; i<solid->triangles.size(); i++)
	{

		// get the triangle vertex
		point V0 = solid->nodes[solid->triangles[i][0]];
		point V1 = solid->nodes[solid->triangles[i][1]];
		point V2 = solid->nodes[solid->triangles[i][2]];

		// get the triangle normal
		point n = (V1 - V0)^(V2 - V0);
		n = n/n.mag();
		// get the triangle center
		c = (V0 + V1 + V2)/3.0;

		dir = point(0,0,1);
		
		// calculate the scalar product
		double prod = dir*n;
		
		// if the triangle normal is aligned with the c direction the triangle
		// is part of the gap surface, use to calculate the radius
		double tol = 1.0 - cos(ang_tol*pi/180.0);
		//if(fabs(prod) > 1.0 - tol && fabs(prod) < 1.0 + tol)
		if(fabs(prod) >= 1.0 - tol)
		{
			// use just the upper surface
			if(prod > 0)
			{
				top += c.z();
				nmax++;
			}
			else
			{
				bottom += c.z();
				nmin++;
			}
		}
		else
		{
			side_triangles.push_back(i);
		}
	}

	
	top /= nmax;
	bottom /= nmin;

	// select just the side triangles having two vertices on the top plane

	// list of nodes defining the surface
	vector<int> nodes_idx;
	vector<vector<int>> tmp_triangle(0);

	for(unsigned int i=0; i<side_triangles.size(); i++)
	{
		// count the number of vertices on the top surface
		int count = 0;
		for(int j=0; j<3; j++)
		{
			if(is_top_node(solid->nodes[solid->triangles[side_triangles[i]][j]]))
				count++;
		}
		if(count == 2)
		{
			vector<int> triangle(3);
			// this triangle will be part of surface
			for(int j=0; j<3; j++)
			{
				vector<int>::iterator it = 
					find(nodes_idx.begin(), nodes_idx.end(), solid->triangles[side_triangles[i]][j]);
				// new node in the list
				if(it == nodes_idx.end())
				{
					nodes_idx.push_back(solid->triangles[side_triangles[i]][j]);
					nodes.push_back(solid->nodes[solid->triangles[side_triangles[i]][j]]);
					triangle[j] = nodes_idx.size() - 1;
				}
				else
				{
					triangle[j] = it - nodes_idx.begin();
				}
			}
			triangles.push_back(triangle);
		}
	}
}
// ------------------------------------------------------------------------- //
void gap_boundary::build(const stl* s)
{
	solid = s;

	// select the triangles in the top plane
	select_triangles();

	// define the triangles neighbors
	define_neighbors();

	// define the boundaries
	define_boundaries();

}
// ------------------------------------------------------------------------- //
void gap_boundary::define_neighbors()
{
	
	neigh.resize(triangles.size());

	// loop to all the triangles
	for(unsigned int i=0; i<triangles.size(); i++)
	{
		// loop to all the triangles
		for(unsigned int j=0; j<triangles.size(); j++)
		{
			// if it is not the same triangle go
			if(j != i)
			{
				// number of shared nodes found
				int found = 0;
				// search for shared nodes
				for(int h=0; h<3;h++)
				{
					vector<int>::iterator it = 
						find(triangles[j].begin(),triangles[j].end(),triangles[i][h]);
					if(it != triangles[j].end())
					{
						if(is_top_node(nodes[triangles[i][h]]))
							found++;	
					}
				}
				if(found == 1)
				{
					neigh[i].push_back(j);
				}
			}
		}
		if(neigh[i].size() != 2)
		{
			Log << "problem defining the neighbors for triangel " << i << gaplog::endl;
			Log << "triangle " << i << " vertices:" << gaplog::endl;
			Log << nodes[triangles[i][0]] << gaplog::endl;
			Log << nodes[triangles[i][1]] << gaplog::endl;
			Log << nodes[triangles[i][2]] << gaplog::endl;
			Log << "List of neighbors: {";
			for(unsigned int h=0; h<neigh[i].size(); h++)
				Log << neigh[i][h] << ", ";
			Log << "}" << gaplog::endl;
			exit(1);
		}
	}
}
// ------------------------------------------------------------------------- //
void gap_boundary::define_boundaries()
{
	vector<bool> processed(triangles.size(),false);
	bool completed = false;
	boundaries.resize(0);
	
	do
	{
		int start = -1;

		for(unsigned int i=0; i<processed.size(); i++)
		{
			if(processed[i] == false)
			{
				start = i;
			}
		}

		if(start == -1)
		{
			completed = true;
			break;
		}
		else
		{
			// triangles already included
			vector<int> incl(0);
			incl.push_back(start);	// push back the first one
			// last triangle included
			int last_incl = start;
			
			// triangle processed
			processed[start] = true;
			
			// boundary node list
			vector<int> nodelist(0);

			// push back the first two top nodes
			for(int h=0; h<3; h++)
			{
				if(is_top_node(nodes[triangles[start][h]]))
				{
					nodelist.push_back(triangles[start][h]);
				}
			}

			while(true)
			{
				// extract two groups of nodes index (each group is composed by two nodes)
				// tn0 -> 2 top nodes for neigh[0]
				// tn1 -> 2 top nodes for neigh[1]
				vector<int> tn0, tn1;
				for(int h=0; h<3; h++)
				{
					if(is_top_node(nodes[triangles[neigh[last_incl][0]][h]]))
						tn0.push_back(triangles[neigh[last_incl][0]][h]);

					if(is_top_node(nodes[triangles[neigh[last_incl][1]][h]]))
						tn1.push_back(triangles[neigh[last_incl][1]][h]);
				}

				// get the last node index in the list defining the new boundary
				int lastnode = nodelist[nodelist.size()-1];

				// determine the next node in the list
				if(lastnode == tn0[0])
				{
					// push back the other
					last_incl = neigh[last_incl][0];
					nodelist.push_back(tn0[1]);
					processed[last_incl] = true;
					// check if the boundary is completed
					if(tn0[1] == nodelist[0])
						break;
				}
				else if(lastnode == tn0[1])
				{
					// push back the other
					last_incl = neigh[last_incl][0];
					nodelist.push_back(tn0[0]);
					processed[last_incl] = true;
					// check if the boundary is completed
					if(tn0[0] == nodelist[0])
						break;
				}
				else if(lastnode == tn1[0])
				{
					// push back the other
					last_incl = neigh[last_incl][1];
					nodelist.push_back(tn1[1]);
					processed[last_incl] = true;
					// check if the boundary is completed
					if(tn1[1] == nodelist[0])
						break;
				}
				else if(lastnode == tn1[1])
				{
					// push back the other
					last_incl = neigh[last_incl][1];
					nodelist.push_back(tn1[0]);
					processed[last_incl] = true;
					// check if the boundary is completed
					if(tn1[0] == nodelist[0])
						break;
				}
				else
				{
					// this should never happen
					Log << "something is wrong during the definition of one of the boundaries" << gaplog::endl;
					Log << lastnode << "\t" 
							 << tn0[0] << "\t" << tn0[1] << "\t" << tn1[0] << "\t" << tn1[1] << gaplog::endl;
					exit(1);
				}
			}

			//Log << gaplog::endl<< "Defined boundary " << boundaries.size() 
			//		 << "\tsize: " << nodelist.size() << gaplog::endl;

			polygon boundary;
			boundary.nv = nodelist.size();
			
			
			for(unsigned int i=0; i<nodelist.size(); i++)
				boundary.vertices.push_back(nodes[nodelist[i]]);	
			
			boundaries.push_back(boundary);

			//ostringstream oss;
			//oss << "boundary." << boundaries.size()-1 << ".vtk";
			//boundary.write_vtk(oss.str().c_str());
			
		}
			
	} while(!completed);

	Log << "  * Found " << boundaries.size() 
			 << " boundary lines defining the gap" << gaplog::endl;


}
// ------------------------------------------------------------------------- //
bool gap_boundary::is_inside(const point& P) const
{
	// crossing number method

	int total_cn = 0;
	
	for(int b=0; b<boundaries.size(); b++)
	{
		int cn = 0;    // the crossing number counter
		const vector<point>& V = boundaries[b].vertices; // vertices
		
		// loop through all edges of the polygon
		for (unsigned int i=0; i<V.size()-1; i++) 
		{  
			// edge from V[i] to V[i+1]
			if
			(
				((V[i].y() <= P.y()) && (V[i+1].y() > P.y())) ||		// an upward crossing
				((V[i].y() > P.y()) && (V[i+1].y() <= P.y())))			// a downward crossing
			{ 
				// compute the actual edge-ray intersect x-coordinate
				double vt = (P.y() - V[i].y()) / (V[i+1].y() - V[i].y());
				if (P.x() < V[i].x() + vt*(V[i+1].x() - V[i].x())) // P.x < intersect
				{
					++cn;   // a valid crossing of y=P.y right of P.x
				}
			}
		}

		total_cn += cn;
	}

	return (total_cn&1);    // 0 if even (out), and 1 if odd (in)
}
// ------------------------------------------------------------------------- //