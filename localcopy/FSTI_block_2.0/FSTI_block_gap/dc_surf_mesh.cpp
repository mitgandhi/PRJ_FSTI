# include "./dc_surf_mesh.h"
# include "../FSTI_Block_dll/log.h"
# include <string>
# include <sstream>

using namespace std;

extern gaplog Log;

// ------------------------------------------------------------------------- //
static void removecomma(string& s) // 0:start 1:end
{
	int idx = s.rfind(',');
	if(idx != std::string::npos)
		s = s.substr(0,idx);
}
// ------------------------------------------------------------------------- //
int dc_surf_mesh::read(const char* mesh_name)
{

	string read_s;

	string inputfile;

	ifstream infile(mesh_name);
	if (!infile.is_open()) 	
	{
		return -1;
	}

	Log << "\nReading DC surface mesh " << mesh_name  << " in Abaqus format\n" 
			 << gaplog::endl;

	double scale = 1e-3;
	bool firstnode = true;
	int shift = 0;
	nodes.resize(0);
	faces.resize(0);

	while(infile) 
	{

		infile >> read_s; 

		// reading list of nodes
		if(read_s.compare("*NODE") == 0) 
		{
			Log << "   * Reading nodes ... ";
			while(true) 
			{
				// store the pointer location before reading
				int pos = infile.tellg(); 
				infile >> read_s;	// node idx
				if (read_s.find("*") != string::npos)
				{
					infile.seekg(pos); // restore the pointer location
					break;
				}
				else 
				{
					if(firstnode)
					{
						removecomma(read_s);
						// convert to a double
						istringstream iss(read_s.c_str());
						int idx;
						iss >> idx;
						if(idx != 1)
						{
							shift = idx;
							//Log << "\nNodes index start from " << idx << gaplog::endl
							//	  << "Nodes index should start from 1. Please check the mesh\n" 
							//		<< gaplog::endl;
							//exit(1);
						}
					}
				}
				int c = 0;
				double x = 0, y = 0, z = 0;
				double val;
				// read nodes	
				while(c < 3) 
				{
					infile >> read_s;	
					// get only numbers
					if(read_s.compare(",") != 0) {
						// remove comma if present
						removecomma(read_s);
						// convert to a double
						istringstream iss(read_s.c_str());
						iss >> val;
						// store the value
						if(c == 0)
							x = val*scale;
						if(c == 1)
							y = val*scale;
						if(c == 2)
							z = val*scale;
						c++;
					}
				}
				nodes.push_back(point(x,y,z));
				firstnode = false;
			}
			
			Log << "done! (Found " << nodes.size() << " nodes)" << gaplog::endl;
		}

		// reading list of elements
		else if(read_s.find("*ELEMENT") != string::npos)
		{
			// new face set
			if(read_s.find("TYPE=S3") != string::npos)
			{
				// extract the name
				string setname;
				string keyword("ELSET="); 
				size_t f = read_s.find(keyword);
				if(f != std::string::npos)
					setname = read_s.substr(f + keyword.size(),read_s.size() - 1);

				Log << "   * Reading faces set " << setname << " ... ";

				while(true) 
				{
					// store the pointer location before reading
					int pos = infile.tellg(); 
					infile >> read_s;	// node idx
					// * is the end of section
					if (read_s.find("*") != string::npos) 
					{
						infile.seekg(pos); // restore the pointer location
						break;
					}
					
					int id;
					removecomma(read_s);
					istringstream iss(read_s.c_str());
					iss >> id; // crap
					
					vector<int> facenodes(3);										
					int c = 0;
					// read nodes	
					while(c < 3) 
					{
						infile >> read_s;	
						// get only numbers
						if(read_s.compare(",") != 0) {
							// remove comma if present
							removecomma(read_s);
							// convert to a double
							istringstream iss(read_s.c_str());
							int val;
							iss >> val;
							// index 1 based so add - 1 to have a 0 based
							facenodes[c] = val - 1 - shift; 
							// increase
							c++;
						}
					}

					faces.push_back(facenodes);
				}

				Log << "done! (Found " << faces.size() << " triangular faces)." << gaplog::endl;
								
			}
			
		}
		
	}

	Log << "\ndone!" << gaplog::endl;
		
	infile.close();

	return 0;

}
// ------------------------------------------------------------------------- //
void dc_surf_mesh::analyze()
{
	
	ni.resize(faces.size());
	ai.resize(faces.size());

	Atot = 0;
	
	point r(0,0,0);
	point Fr(0,0,0);
	point Mr(0,0,0);

	xR = yR = zR = 0.0;

	for(unsigned int i=0; i<faces.size(); i++)
	{
		point c = (nodes[faces[i][0]] + nodes[faces[i][1]] + nodes[faces[i][2]])/3.0;	// face center
		point r0 = (nodes[faces[i][1]] - nodes[faces[i][0]]);
		point r1 = (nodes[faces[i][2]] - nodes[faces[i][0]]);
		ni[i] = r0 ^ r1;							// get the normal
		ai[i] = 0.5*ni[i].mag();			// get the area
		ni[i] = ni[i]/ni[i].mag();		// normalize

		point F0i(ai[i]*ni[i].x(), ai[i]*ni[i].y(), ai[i]*ni[i].z());
		point M0i = c ^ F0i;
		
		xR += c.x()*ai[i];
		yR += c.y()*ai[i];
		zR += c.z()*ai[i];

		Fr = Fr + F0i;
		Mr = Mr + M0i;

		Atot += ai[i];
	}

	if(Fr.z() > 0)	// invert the sign (flip the normals direction)
	{
		Fr = -1.0*Fr;
		Mr = -1.0*Mr;
	}

	xR = xR / Atot;
	yR = yR / Atot;
	zR = zR / Atot;

	F0x = Fr.x(), F0y = Fr.y(), F0z = Fr.z();
	M0x = Mr.x(), M0y = Mr.y(), M0z = Mr.z();

	Log << "done! \nAREA: " << Atot << "\nxR:" << xR << "\nyR:" << yR << "\nzR:" << zR<<   gaplog::endl;
	Log << "\nF0x: " << Fr.x() << "\nF0y =" <<  Fr.y() << "\nF0z = " << Fr.z() <<   gaplog::endl;
}
// ------------------------------------------------------------------------- //
void dc_surf_mesh::initialize(const char* mesh_name)
{
	int found = read(mesh_name);

	if(found < 0)
	{
		available = false;
	}
	else
	{
		available = true;
		analyze();
	}
}
// ------------------------------------------------------------------------- //

