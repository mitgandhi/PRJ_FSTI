# include "./block_gap_main.h"
# include "./macro_geometry.h"
# include "../FSTI_Block_dll/log.h"
# include "../interpolation/interpolation.h"
# include <string>
# include <sstream>
# include <fstream>
# include <iostream>


using namespace std;

extern gaplog Log;

# define pi 3.14159265358979323846

// ------------------------------------------------------------------------- //
macro_geometry::macro_geometry(scalar_field& _h)
:
	h(_h),
	dh(0)
{
	
}
// ------------------------------------------------------------------------- //
macro_geometry::~macro_geometry()
{
	if(dh!=0)
		delete dh;
}
// ------------------------------------------------------------------------- //
void macro_geometry::apply(int rotation_steps)
{
	h = h + dh->cshift(rotation_steps);
}
// ------------------------------------------------------------------------- //
void macro_geometry::read(const char* filename)
{
	
}
// ------------------------------------------------------------------------- //
macro_waved::macro_waved(scalar_field& _h, const char* filename) : macro_geometry(_h)
{
	dh = new scalar_field(h.mesh);
	
	read(filename);

	// set the internal field
	for(int i=0; i<h.m; i++) 
	{
		for(int j=0; j<h.n; j++) 
		{
			int id = i*h.n + j;
			gap_elm e = h.mesh->elements[id];
			double angle = j*2*pi/h.n;
			dh->in[id] = 1e-6*A*sin(f*(angle+offset*pi/180));
		}
	}
	for(int j=0; j<h.n; j++) 
	{
		double angle = j*2*pi/h.n;
		dh->bound.inner[j] = dh->bound.outer[j] = 1e-6*A*sin(f*(angle+offset));
	}
}
// ------------------------------------------------------------------------- //
void macro_waved::read(const char* filename)
{
	
	ifstream in(filename);

	if(!in.is_open()) 	
	{
		Log << "\nmacro_geometry::read(): Unable to open " << filename << gaplog::endl;
		exit(1);
	}
	else
	{
		Log << "\nReading macro_geometry file " << filename 
				<< gaplog::endl << gaplog::endl;
	}

	string tmp;
	string line;
	double val = 0;
	int linenumber = 1;

	while(getline(in,line)) 
	{
		if(line.size() > 0)
		{
			istringstream iss (line,istringstream::in);
			iss >> tmp;
			size_t comment = tmp.find("//"); // check if there is a comment
			
			// if is not a comment read
			if (comment == string::npos) 
			{	
				if(tmp.compare("frequency") == 0)
				{
					iss >> tmp;
					istringstream(tmp) >> f;
				}
				if(tmp.compare("amplitude") == 0)
				{
					iss >> tmp;
					istringstream(tmp) >> A;
					Log << "\nAmplitude read... " << filename
						<< gaplog::endl << gaplog::endl;
				}
				if(tmp.compare("offset") == 0)
				{
					iss >> tmp;
					istringstream(tmp) >> offset;
				}

			}

			linenumber++;
		}
	}
}
// ------------------------------------------------------------------------- //
macro_axisymmetric::macro_axisymmetric(scalar_field& _h, const char* filename) : macro_geometry(_h)
{
	
	dh = new scalar_field(h.mesh);

	read(filename);

	double r_in = h.mesh->r_in;
	double r_out = h.mesh->r_out;
	vector<double> new_ri(0);
	vector<double> new_dhi(0);

	if(ri[0] > r_out)
	{
		Log << "\nmacro_axisymmetric: WARNING first radius is bigger than outer radius!"
				 << gaplog::endl;
		return;
	}

	// copy just within the sealing land region
	for(unsigned int i=0; i<ri.size(); i++)
	{
		if(ri[i] >= r_in && ri[i] <= r_out)
		{
			new_ri.push_back(ri[i]);
			new_dhi.push_back(dhi[i]);
		}
	}

	//cout << r_in << "\t" << r_out << endl;

	// need to add an element at r_in position
	if(ri[0] < new_ri[0])
	{

		new_ri.insert(new_ri.begin(), r_in);
		int k=0;
		while(ri[k] < new_ri[0])
			k++;
		double dhi0 =  dhi[k] + (dhi[k+1] - dhi[k])*(r_in - ri[k])/(ri[k+1] - ri[k]);
		new_dhi.insert(new_dhi.begin(), dhi0);
	}
	// need to add an element at r_out position
	if(ri.back() > new_ri.back())
	{
		new_ri.insert(new_ri.end(), r_out);
		int k=0;
		while(ri[k] < new_ri.back())
			k++;
		
		double dhin = dhi[k-1] + (dhi[k] - dhi[k-1])*(r_out - ri[k-1])/(ri[k] - ri[k-1]);
		new_dhi.insert(new_dhi.end(), dhin);
	}

	ri = new_ri;
	dhi = new_dhi;
	
	// set to zero if some radial position is not defined at the beginning
	if(fabs(ri[0] - r_in) > 1e-6)
	{	
		if(ri[0] > r_in)
		{
			ri.insert(ri.begin(), ri[0]);
			dhi.insert(dhi.begin(), 0);
			ri.insert(ri.begin(), r_in);
			dhi.insert(dhi.begin(), 0);
		}
	}
	// set to zero if some radial position is not defined at the end
	if(fabs(ri.back() - r_out) > 1e-6)
	{
		if(ri.back() < r_out)
		{
			ri.insert(ri.end(), ri.back());
			dhi.insert(dhi.end(), 0);
			ri.insert(ri.end(), r_out);
			dhi.insert(dhi.end(), 0);
		}
	}
	
	// set the internal field
	for(int i=0; i<h.m; i++) 
	{
		double r = h.mesh->elements[i*h.n].r;
		int k = 0;
		while(ri[k] <= r)
			k++;
		
		for(int j=0; j<h.n; j++) 
		{
			int id = h.n*i + j;
			dh->in[id] = dhi[k-1] + (dhi[k] - dhi[k-1])*(r-ri[k-1])/(ri[k]-ri[k-1]);
			
		}
		
	}

	
}
// ------------------------------------------------------------------------- //
void macro_axisymmetric::read(const char* filename)
{
	
	ifstream in(filename);

	if(!in.is_open()) 	
	{
		Log << "\nmacro_geometry::read(): Unable to open " << filename << gaplog::endl;
		exit(1);
	}
	else
	{
		Log << "\nReading macro_geometry file " << filename 
				<< gaplog::endl << gaplog::endl;
	}

	string tmp;
	string line;
	double val = 0;
	int linenumber = 1;

	while(getline(in,line)) 
	{
		if(line.size() > 0)
		{
			istringstream iss (line,istringstream::in);
			iss >> tmp;
			size_t comment = tmp.find("//"); // check if there is a comment
			// if is not a comment read
			if (comment == string::npos) 
			{	
				double val;
				
				// read ri
				istringstream(tmp) >> val;
				ri.push_back(1e-3*val); // [mm->m]
				
				// read dhi
				iss >> tmp;
				if(tmp.size() > 0)
				{
					istringstream(tmp) >> val;
					dhi.push_back(1e-6*val);	// [um->m]
				}
				else
				{
					Log << "\nmacro_geometry::read(): Error in file " << filename 
							<< " at line " << linenumber << gaplog::endl;
					exit(1);
				}
			}

			linenumber++;
		}
	}
}
// ------------------------------------------------------------------------- //
macro_spherical::macro_spherical(scalar_field& _h, const char* filename) : macro_geometry(_h)
{

	dh = new scalar_field(h.mesh);

	read(filename);

	double R0 = 0.5*D;
	double Re = 0.5*Dm;

	cout << R0 << "\t" << Re << endl;
	
	// offset between the nominal and the effective radius
	double d = Re - R0;

	double dh_inner = 0;

	// set the internal field
	for(int i=0; i<h.m; i++) 
	{
		// position on the nominal spherical surface
		double r0 = h.mesh->elements[i*h.n].r;
		double alpha = asin(r0/R0);
		double z0 = R0*cos(alpha);

		// position on the effective spherical surface
		double ze = (d +  sqrt(pow(d,2) + (pow(Re,2) - pow(d,2))/pow(cos(alpha),2)))/(1 + pow(tan(alpha),2));
		double re = ze*tan(alpha);

		double rho0 = sqrt(pow(r0,2) + pow(z0,2));
		double rhoe = sqrt(pow(re,2) + pow(ze,2));

		// change in film thickness due to difference in diameter
		double dhr = rhoe - rho0;

		if(i==0)
			dh_inner = dhr;
			
		for(int j=0; j<h.n; j++) 
		{
			int id = h.n*i + j;
			dh->in[id] = dh_inner - dhr;
		}
	}
}
// ------------------------------------------------------------------------- //
void macro_spherical::read(const char* filename)
{
	ifstream in(filename);

	if(!in.is_open()) 	
	{
		Log << "\nmacro_geometry::read(): Unable to open " << filename << gaplog::endl;
		exit(1);
	}
	else
	{
		Log << "\nReading macro_geometry file " << filename 
				<< gaplog::endl << gaplog::endl;
	}

	string tmp;
	string line;
	double val = 0;
	int linenumber = 1;

	while(getline(in,line)) 
	{
		if(line.size() > 0)
		{
			istringstream iss (line,istringstream::in);
			iss >> tmp;
			size_t comment = tmp.find("//"); // check if there is a comment
			
			// if is not a comment read
			if (comment == string::npos) 
			{	
				if(tmp.compare("D") == 0)
				{
					iss >> tmp;
					istringstream(tmp) >> D;
					D *= 1e-3;
				}
				if(tmp.compare("Dm") == 0)
				{
					iss >> tmp;
					istringstream(tmp) >> Dm;
					Dm *= 1e-3;
				}
			}

			linenumber++;
		}
	}
}
// ------------------------------------------------------------------------- //
macro_grid::macro_grid(scalar_field& _h, const char* filename) : macro_geometry(_h)
{

	// create the dh scalar field;
	dh = new scalar_field(h.mesh);

	// read the input file
	read(filename);

	// read the grid file and data
	from.read(grid_file.c_str());

	system("IF NOT EXIST output (mkdir output)");
	system("IF NOT EXIST output\\block (mkdir output\\block)");
	from.writeVTK("./output/block/macro_grid.vtk");
	
	interpolation I;
	interpl_grid to;
	to.initialize(h.mesh->elements.size(),h.mesh->nodes.size(), 4);
	for(unsigned int i=0;i<h.mesh->nodes.size(); i++)
		to.nodes[i] = h.mesh->nodes[i];
	for(unsigned int i=0;i<h.mesh->elements.size(); i++)
	{
		if(body.compare("CB") == 0)
			to.active[i] = !h.mesh->cb_elms_0[i];
		else
			to.active[i] = !h.mesh->vp_elms[i];
		// define nodes
		for(int j=0; j<to.ELM_NDS; j++)
			to.elements[i][j] = h.mesh->elements[i].nds[j];
	}

	if(data_at.compare("CELLS") == 0)
		I.cellsTocells(&from, &to);
	else
		I.pointsTocells(&from, &to);
	
	for(unsigned int i=0; i<to.active.size(); i++)
	{
		if(to.active[i])
			dh->in[i] = to.cells_data[i];
		else
			dh->in[i] = 0;
	}
	
}
// ------------------------------------------------------------------------- //
void macro_grid::read(const char* filename)
{
	ifstream in(filename);

	if(!in.is_open()) 	
	{
		Log << "\nmacro_geometry::read(): Unable to open " << filename << gaplog::endl;
		exit(1);
	}
	else
	{
		Log << "\nReading macro_geometry file " << filename 
				<< gaplog::endl << gaplog::endl;
	}

	string tmp;
	string line;
	double val = 0;
	int linenumber = 1;

	while(getline(in,line)) 
	{
		if(line.size() > 0)
		{
			istringstream iss (line, istringstream::in);
			iss >> tmp;
			size_t comment = tmp.find("//"); // check if there is a comment
					
			// if is not a comment read
			if (comment == string::npos) 
			{	
				if(tmp.compare("grid_file") == 0)
					iss >> grid_file;
				if(tmp.compare("data_at") == 0)
					iss >> data_at;
				if(tmp.compare("body") == 0)
					iss >> body;
			}

			linenumber++;
		}
	}


	if(grid_file.size() == 0)
	{
		Log << "\nExpected keyword grid_file followed by the grid file path.\n";
		exit(1);
	}

	if(data_at.size() == 0)
	{
		Log << "\nExpected keyword data_at followed by either CELLS or NODES.\n";
		exit(1);
	}

	if(body.size() == 0)
	{
		Log << "\nExpected keyword body followed by either CB or VP.\n";
		exit(1);
	}

	in.close();
}
// ------------------------------------------------------------------------- //
void macro_stl::read(const char* filename)
{
	ifstream in(filename);

	if(!in.is_open()) 	
	{
		Log << "\nmacro_geometry::read(): Unable to open " << filename << gaplog::endl;
		exit(1);
	}
	else
	{
		Log << "\nReading macro_geometry file " << filename 
				<< gaplog::endl << gaplog::endl;
	}

	string tmp;
	string line;
	double val = 0;
	int linenumber = 1;

	while(getline(in,line)) 
	{
		if(line.size() > 0)
		{
			istringstream iss (line, istringstream::in);
			iss >> tmp;
			size_t comment = tmp.find("//"); // check if there is a comment
					
			// if is not a comment read
			if (comment == string::npos) 
			{	
				if(tmp.compare("body_type") == 0)
					iss >> body_type;
				if(tmp.compare("stl_file") == 0)
					iss >> stl_file;
				if(tmp.compare("notch_depth") == 0)
					iss >> notch_depth;
			}

			linenumber++;
		}
	}
}
// ------------------------------------------------------------------------- //
macro_stl::macro_stl(scalar_field& _h, const char* filename) : macro_geometry(_h)
{

	// create the dh scalar field;
	dh = new scalar_field(h.mesh);

	// read the input file
	read(filename);

	body.read(stl_file.c_str(), 1e-3);
	boundary.build(&body);
	
	double k = (body_type.compare("VP") == 0) ? -1e-6 : 1e-6;
	
	for(unsigned int i=0; i<h.mesh->elements.size(); i++)
	{
		point center(h.mesh->elements[i].x, h.mesh->elements[i].y);
		if(boundary.is_inside(center))
			dh->in[i] = k*notch_depth;
	}
	
}
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
macro_fluid_grid::macro_fluid_grid(scalar_field& _h, const char* filename) : macro_geometry(_h)
{
	
	dh = new scalar_field(h.mesh);
	
	read(filename);

	// set the internal field
	for(int e=0; e<macro_height.size(); e++)
	{
		if(i[e] < 0 || i[e] >= h.m || j[e] < 0 || j[e] >= h.n)
		{
			Log << "\nmacro_geometry::macro_fluid_grid(): Invalid cell i,j value of " << i[e] << ", " << j[e]  << gaplog::endl;
			exit(1);
		}
		int id = i[e]*h.n + j[e];
		dh->in[id] = macro_height[e];
	}

	//set inner bound to use the fluid circumfrence 0 values
	for(int j=0; j<h.n; j++) 
	{
		int id = 0*h.n + j;
		dh->bound.inner[j] = dh->in[id];
	}

	//set outer bound to use the fluid circumfrence m-1 values
	for(int j=0; j<h.n; j++) 
	{
		int id = (h.m-1)*h.n + j;
		dh->bound.outer[j] = dh->in[id];
	}
	
	
}
// ------------------------------------------------------------------------- //
void macro_fluid_grid::read(const char* filename)
{
	
	ifstream in(filename);

	if(!in.is_open()) 	
	{
		Log << "\nmacro_geometry::read(): Unable to open " << filename << gaplog::endl;
		exit(1);
	}
	else
	{
		Log << "\nReading macro_geometry file " << filename 
				<< gaplog::endl << gaplog::endl;
	}

	string tmp;
	string line;
	double val = 0;
	int linenumber = 1;

	while(getline(in,line)) 
	{
		if(line.size() > 0)
		{
			istringstream iss (line,istringstream::in);
			iss >> tmp;
			size_t comment = tmp.find("//"); // check if there is a comment
			// if is not a comment read
			if (comment == string::npos) 
			{	
				int val;
				
				// read i
				istringstream(tmp) >> val;
				i.push_back(val);
				
				// read j
				iss >> tmp;
				if(tmp.size() > 0)
				{
					istringstream(tmp) >> val;
					j.push_back(val);
				}
				else
				{
					Log << "\nmacro_geometry::read(): Error in file " << filename 
							<< " at line " << linenumber << gaplog::endl;
					exit(1);
				}

				// read h
				iss >> tmp;
				if(tmp.size() > 0)
				{
					double dval;
					istringstream(tmp) >> dval;
					macro_height.push_back(1e-6*dval);	// [um->m]
				}
				else
				{
					Log << "\nmacro_geometry::read(): Error in file " << filename 
							<< " at line " << linenumber << gaplog::endl;
					exit(1);
				}
			}

			linenumber++;
		}
	}
}
