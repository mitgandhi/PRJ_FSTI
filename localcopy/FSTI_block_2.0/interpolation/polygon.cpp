# include "polygon.h"
# include <iostream>
# include <cmath>
# include <fstream>

using namespace std;

polygon::polygon(int _nv)
{
	nv = _nv;
	vertices.resize(nv);
}
// ------------------------------------------------------------------------- //
polygon::~polygon()
{
}
// ------------------------------------------------------------------------- //
polygon& polygon::operator=(const polygon& p)
{
	nv = p.nv;
	vertices.resize(nv);
	for(int i=0; i<nv; i++)
		vertices[i] = p.vertices[i];
	return *this;
}
// ------------------------------------------------------------------------- //
void polygon::print()
{
	cout << "Polygon [Vertices: " << nv << "\tArea: " << area() << "]" << endl;
	for(int i=0; i<nv; i++)
		cout << "\t(" << vertices[i].x() << ", " << vertices[i].y() << ")" << endl;
}
// ------------------------------------------------------------------------- //
double polygon::area()
{
	double A = 0;

	if(vertices.size() > 0)
	{
		for(int i=1; i<nv; i++)
		{
			A += vertices[i-1].x()*vertices[i].y() - vertices[i].x()*vertices[i-1].y();
		}
		A += vertices[nv-1].x()*vertices[0].y() - vertices[0].x()*vertices[nv-1].y();
		
		A *= 0.5;
	}

	return fabs(A);
}
// ------------------------------------------------------------------------- //
point polygon::center()
{
	double cx = 0, cy = 0;
	const vector<point>& v = vertices;
	for(unsigned int i=0; i<vertices.size(); i++)
	{
		cx += v[i].x();
		cy += v[i].y();
	}
		
	return point(cx/vertices.size(), cy/vertices.size());
}
// ------------------------------------------------------------------------- //
bool polygon::is_inside(const point& P)
{
	int cn = 0;    // the crossing number counter
	
	vector<point> V = vertices; // vertices
	if(V[0] != V[V.size()-1])
		V.push_back(V[0]);
		
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

	return (cn&1);
}
// ------------------------------------------------------------------------- //
std::vector<polygon> polygon::clip(const polygon& p, gpc_op opt)
{
	
	// create the subject polygon data structure to interface with gpc
	int s_nc = 1;
	int* s_nv = new int[1];
	s_nv[0] = nv;
	double** s_x = new double*[1];
	s_x[0] = new double[nv];
	double** s_y = new double*[1];
	s_y[0] = new double[nv];
	for(int i=0; i<nv; i++)
	{
		s_x[0][i] = vertices[i].x();
		s_y[0][i] = vertices[i].y();
	}

	// create the clipping polygon data structure to interface with gpc
	int c_nc = 1;
	int* c_nv = new int[1];
	c_nv[0] = p.nv;
	double** c_x = new double*[1];
	c_x[0] = new double[p.nv];
	double** c_y = new double*[1];
	c_y[0] = new double[p.nv];
	for(int i=0; i<p.nv; i++)
	{
		c_x[0][i] = p.vertices[i].x();
		c_y[0][i] = p.vertices[i].y();
	}

	// gpc stuff
	gpc_polygon subject, clip, result;
	gpc_create_polygon(1,s_nv,s_x,s_y,&subject);
	gpc_create_polygon(1,c_nv,c_x,c_y,&clip);
	
	// GPC_DIFF			Difference
  // GPC_INT			Intersection
  // GPC_XOR,     Exclusive or
  // GPC_UNION		Union
	gpc_polygon_clip(opt, &subject, &clip, &result);
	
	std::vector<polygon> intersections;
	for(int i=0; i<result.num_contours; i++)
	{
		// copy the result to the resulting polygon
		polygon I;
		I.nv = result.contour[i].num_vertices;
		I.vertices.resize(I.nv);
		for(int j=0; j<I.nv; j++)
		{
			I.vertices[j]  = 
				point(result.contour[i].vertex[j].x, result.contour[i].vertex[j].y);
		}
		intersections.push_back(I);
	}
	
	// delete all the temporary stuff
	delete [] s_x[0];
	delete [] s_y[0];
	delete [] s_nv;
	delete [] s_x;
	delete [] s_y;
	delete [] c_x[0];
	delete [] c_y[0];
	delete [] c_x;
	delete [] c_y;
	delete [] c_nv;
	
	gpc_free_polygon(&subject);
	gpc_free_polygon(&clip);
  gpc_free_polygon(&result);
		
	
	return intersections;
}
// ------------------------------------------------------------------------- //
void polygon::write_vtk(const char* filename)
{
	ofstream vtk(filename);
	if (!vtk.is_open()) 
	{
		cout << "Error opening " << filename << ".vtk" << endl;
		system("pause");
		exit(1);
	}

	vtk << "# vtk DataFile Version 2.0" << endl 
			<< "vtk output" << endl 
			<< "ASCII" << endl 
			<< "DATASET UNSTRUCTURED_GRID" << endl 
			<< "POINTS " << vertices.size() << " double" 
			<< endl;

		
	// --------------------------------- mesh nodes ----------------------------------//
	
	for(unsigned int i = 0; i < vertices.size(); i++) 
	{
		vtk << scientific << vertices[i].x() << "\t" 
				<< scientific << vertices[i].y() << "\t" 
				<< scientific << vertices[i].z() << endl;
	}
	
	vtk << endl;

	// ----------------------------- elements definition ------------------------------- //

	vtk << "CELLS " << 1 << "\t" 
		<< 2 + vertices.size() << endl;

	vtk << vertices.size() + 1 << "\t";
	for(unsigned int i = 0; i < vertices.size(); i++) 
	{
			vtk << i << "\t";
	}
	vtk << 0 << "\t"; // repeat the first
	
	vtk << endl;

	// ----------------------------- elements type ------------------------------- //

	vtk << "CELL_TYPES " << 1 << endl;
			
	vtk << 4 << endl;
		
}
// ------------------------------------------------------------------------- //
void polygon::write_vtk(const char* filename, const std::vector<polygon>& surf)
{
	ofstream vtk(filename);
	if (!vtk.is_open()) 
	{
		cout << "Error opening " << filename << ".vtk" << endl;
		system("pause");
		exit(1);
	}

	int nvertices = 0;
	for(int poly = 0; poly<surf.size(); ++poly)
	{
			// --------------------------------- mesh nodes ----------------------------------//
			for(unsigned int i = 0; i < surf[poly].vertices.size(); i++) 
			{
				nvertices++;
			}
	}

	vtk << "# vtk DataFile Version 2.0" << endl 
			<< "vtk output" << endl 
			<< "ASCII" << endl 
			<< "DATASET UNSTRUCTURED_GRID" << endl 
			<< "POINTS " << nvertices << " double" 
			<< endl;

	

	for(int poly = 0; poly<surf.size(); ++poly)
	{
			// --------------------------------- mesh nodes ----------------------------------//
			for(unsigned int i = 0; i < surf[poly].vertices.size(); i++) 
			{
				vtk << scientific << surf[poly].vertices[i].x() << "\t" 
						<< scientific << surf[poly].vertices[i].y() << "\t" 
						<< scientific << surf[poly].vertices[i].z() << endl;
			}
	
			vtk << endl;
	}

	// ----------------------------- elements definition ------------------------------- //
		vtk << "CELLS " << surf.size() << "\t" 
			<< surf.size() + nvertices << endl;

	int gidctr = 0;
	for(int poly = 0; poly<surf.size(); ++poly)
	{
		vtk << surf[poly].vertices.size() << "\t";
		for(unsigned int i = 0; i < surf[poly].vertices.size(); i++) 
		{
				vtk << gidctr << "\t";
				gidctr++;
		}
		//vtk << gidctr-surf[poly].vertices.size() << "\t"; // repeat the first
	
		vtk << endl;
	}

	// ----------------------------- elements type ------------------------------- //

	vtk << "CELL_TYPES " << surf.size() << endl;
	for(int poly = 0; poly<surf.size(); ++poly)
	{
		vtk << 7 << endl;
	}
			
	
}