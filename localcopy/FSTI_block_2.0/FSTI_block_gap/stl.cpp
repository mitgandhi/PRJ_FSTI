# include <iostream>
# include <fstream>
# include <vector>
# include <string>
# include <sstream>
# include <algorithm>
# include <cstdlib>

# include "../interpolation/point.h"
# include "./stl.h"
# include "../FSTI_Block_dll/log.h"

# define pi 3.14159265358979323846
# define VTK_TRIANGLE 5

extern class gaplog Log;

using namespace std;

stl::stl() 
{

}
// ------------------------------------------------------------------------- //
void stl::read(const string& STLfile, double scale)
{
	ostringstream msg;
	// open the input file
	ifstream in;
	in.open(STLfile.c_str());
	if (!in.is_open()) 	
	{
		Log << gaplog::endl << "STL::Error opening " << STLfile << gaplog::endl;
		exit(1);
	}

	// ------------------------- Reading from STL file ------------------------//

	// stuff for the file reading
	string word; 
	string tmp;
	int nsolids = 0;
	
	int ID(0); // ID for each triangle
	
	triangles.resize(0);

	while(in.tellg() > -1) 
	{

		in >> word;

		if (word == "solid") // Beginning of solid definition
		{ 
			nsolids++;
		}

		ostringstream oss;
		oss << "solid." << nsolids;
		name = oss.str();

		if (word == "facet") // Beginning of a new triangle
		{ 
			
			// used to read coordinates		
			double c0,c1,c2;
			
			in >> tmp; 
			if (tmp == "normal") // check for the "normal" keyword
				in >> c0;
			else 
			{
				Log << gaplog::endl 
						 << "stl::stl: error in STL file, check for the correct structure!" 
						 << gaplog::endl;
				exit(1);
			}

			in >> c1;
			in >> c2;
			
			point normal(c0,c1,c2); // used for double check
			
			// -------------------------- first vertex --------------------------- //
			in >> tmp;
			if (tmp == "outer") // check for the "outer" keyword
			{ 
				in >> tmp; 
				if (tmp == "loop") // check for the "loop" keyword
				{ 
						in >> tmp;
						if (tmp == "vertex") // check for the "vertex" keyword
							in >> c0;
						else 
						{
							Log << gaplog::endl
									 << "stl::stl: error in STL file, check for the correct structure!" 
									 << gaplog::endl;
							exit(1);
						}
				}
				else 
				{
					Log << gaplog::endl
							 << "stl::stl: error in STL file, check for the correct structure!" 
							 << gaplog::endl;
					exit(1);
				}
			}
			else 
			{
				Log << gaplog::endl
						 << "stl::stl: error in STL file, check for the correct structure!" 
						 << gaplog::endl;
				exit(1);
			}
						
			in >> c1;
			in >> c2;

			c0 *= scale;
			c1 *= scale;
			c2 *= scale;

			point V(c0,c1,c2);
			int V0 = -1;

			bool add = true;

			for(unsigned int i=0; i<nodes.size(); i++)
			{
				if(V == nodes[i])
				{
					add = false;
					V0 = i;
				}
			}

			if(add)
			{
				nodes.push_back(V);
				V0 = nodes.size() - 1;
			}

			
			// -------------------------- second vertex -------------------------- //

			in >> tmp; 
			if (tmp == "vertex") // check for the "vertex" keyword
				in >> c0;
			else 
			{
				Log << gaplog::endl 
						 << "stl::stl: error in STL file, check for the correct structure!" 
						 << gaplog::endl;
				exit(1);
			}  
		 
			in >> c1;
			in >> c2;

			c0 *= scale;
			c1 *= scale;
			c2 *= scale;
			
			V = point(c0,c1,c2);

			int V1 = -1;

			add = true;

			for(unsigned int i=0; i<nodes.size(); i++)
			{
				if(V == nodes[i])
				{
					add = false;
					V1 = i;
				}
			}

			if(add)
			{
				nodes.push_back(V);
				V1 = nodes.size() - 1;
			}
				

			// --------------------------- third vertex -------------------------- //
			
			in >> tmp;
			if (tmp == "vertex") // check for the "vertex" keyword
				in >> c0;
			else 
			{
				Log << gaplog::endl
						 << "stl::stl: error in STL file, check for the correct structure!" 
						 << gaplog::endl;
				exit(1);
			}
			
			in >> c1;
			in >> c2;

			c0 *= scale;
			c1 *= scale;
			c2 *= scale;
			
			V = point(c0,c1,c2);

			int V2 = -1;

			add = true;

			for(unsigned int i=0; i<nodes.size(); i++)
			{
				if(V == nodes[i])
				{
					add = false;
					V2 = i;
				}
			}

			if(add)
			{
				nodes.push_back(V);
				V2 = nodes.size() - 1;
			}
			

			// ------------------ end of the triangle definition ----------------- //
			
			in >> tmp;

			if (tmp == "endloop") // check for the "endloop" keyword
			{ 
				in >> tmp;
				if (tmp != "endfacet") // check for the "endfacet" keyword
				{	
					Log << gaplog::endl 
							 << "stl::stl: error in STL file, check for the correct structure!" 
							 << gaplog::endl;
					exit(1);
				}
			}
			else 
			{
				Log << gaplog::endl
						 << "stl::stl: error in STL file, check for the correct structure!" 
						 << gaplog::endl;
				exit(1);
			}
			
			// define the triangle with the address of the three nodes
			vector<int> triangle(3);
			triangle[0] = V0, triangle[1] = V1, triangle[2] = V2;
									
			// push triangle in container
			triangles.push_back(triangle); 

			ID++; // increment the ID value
		}

		if (word == "endsolid") 
		{ 
			// End of solid definition
		}
		
	}

	in.close();

}
// ------------------------------------------------------------------------- //
stl::~stl() 
{
}
// ------------------------------------------------------------------------- //
void stl::writeVTK()
{
	
	ofstream vtk(string(name + string(".vtk")).c_str());
	if (!vtk.is_open()) 
	{
		Log << "Error opening " << name << ".vtk" << gaplog::endl;
		system("pause");
		exit(1);
	}

	vtk << "# vtk DataFile Version 2.0\n"
			<< "vtk output\n"
			<< "ASCII\n"
			<< "DATASET UNSTRUCTURED_GRID\n"
			<< "POINTS " << nodes.size() << " double\n\n";

		
	// --------------------------------- mesh nodes ----------------------------------//
	
	for(unsigned int i = 0; i < nodes.size(); i++) 
	{
		vtk <<	scientific << nodes[i].x() << "\t" << 
				scientific << nodes[i].y() << "\t" << 
				scientific << nodes[i].z() << endl;
	}

	vtk << endl;

	// ----------------------------- elements definition ------------------------------- //

	vtk << "CELLS " << triangles.size() << "\t" 
		<< (1 + 3)*triangles.size() << endl;

	for(unsigned int i = 0; i < triangles.size(); i++) 
	{
			vtk << 3 << "\t";
			vtk << triangles[i][0] << "\t"
				<< triangles[i][1] << "\t"
				<< triangles[i][2] << endl;
				
			vtk << endl;
	}
	
	vtk << endl;

	// ----------------------------- elements type ------------------------------- //

	vtk << "CELL_TYPES " << triangles.size() << endl;
			
	for (unsigned int i = 0; i < triangles.size(); i++) 
	{
		vtk << VTK_TRIANGLE << endl;
	}	

	vtk << "CELL_DATA" << "\t" << triangles.size() << endl;
	
}
// ------------------------------------------------------------------------- //