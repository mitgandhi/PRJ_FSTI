# include "./interpl_grid.h"
# include "../FSTI_Block_dll/log.h"

# include <fstream>
# include <vector>
# include <string>
# include <sstream>

# define VTK_T 10

using namespace std;

extern class gaplog Log;

// ------------------------------------------------------------------------- //
interpl_grid::interpl_grid() : ELM_NDS(0) { }
// ------------------------------------------------------------------------- //
interpl_grid::interpl_grid(int _ne, int _nn, int elmnds)
{
	ELM_NDS = elmnds;
	nodes.resize(_nn,point(0,0));
	elements.resize(_ne,std::vector<int>(ELM_NDS));
	active.resize(_ne,true);
	points_data.resize(_nn, 0.0);
	cells_data.resize(_ne, 0.0);
}
// ------------------------------------------------------------------------- //
void interpl_grid::initialize(int _ne, int _nn, int elmnds)
{
	ELM_NDS = elmnds;
	nodes.resize(_nn,point(0,0));
	elements.resize(_ne,std::vector<int>(ELM_NDS));
	active.resize(_ne,true);
	points_data.resize(_nn, 0.0);
	cells_data.resize(_ne, 0.0);
}
// ------------------------------------------------------------------------- //
int interpl_grid::isInside(int elmidx, const point& P) const
{
	int cn = 0;    // the crossing number counter
	std::vector<point> V(ELM_NDS + 1);
	for(int i=0; i<ELM_NDS; i++)
	{
		V[i] = nodes[elements[elmidx][i]];
	}
	V[ELM_NDS] = V[0];
	
	// loop through all edges of the polygon
	for (int i = 0; i < ELM_NDS; i++) 
	{  
		// edge from V[i] to V[i+1]
		bool upward_crossing = ( V[i].y() <= P.y() ) && ( V[i+1].y() > P.y() );
		bool downward_crossing = ( V[i].y() > P.y() ) && ( V[i+1].y() <= P.y() );
		
		if (upward_crossing || downward_crossing)	
		{ 
			// compute the actual edge-ray intersect x-coordinate
			double vt = (P.y() - V[i].y()) / (V[i+1].y() - V[i].y());
			if (P.x() < V[i].x() + vt * (V[i+1].x() - V[i].x())) // P.x < intersect
				++cn;   // a valid crossing of y=P.y right of P.x
		}
	}
	
	return (cn&1);    // 0 if even (out), and 1 if odd (in)
}
// ------------------------------------------------------------------------- //
void interpl_grid::cellsToPoints()
{
	// determin wich elements share each node
	vector<vector<int>> shared_cells(nodes.size());
	for(unsigned int i=0; i<elements.size(); i++)
	{
		if(active[i])
		{
			// add just if it is an active element
			for(unsigned int j=0; j<elements[i].size(); j++)
				shared_cells[elements[i][j]].push_back(i);
		}
	}

	// cet the elements center
	vector<point> centers(elements.size());
	{
		for(unsigned int i=0; i<elements.size(); i++)
		{
			centers[i] = 0;
			for(unsigned int j=0; j<elements[i].size(); j++)
				centers[i] = centers[i] + nodes[elements[i][j]]/ELM_NDS;
		}
	}

	// interpolate
	for(unsigned int i=0; i<nodes.size(); i++)
	{
		// get the distances
		vector<double> wij(0);
		double sum_wij = 0;
		for(unsigned int j=0; j<shared_cells[i].size(); j++)
		{
			double d = nodes[i].distance(centers[shared_cells[i][j]]);
			double wj = 1.0/pow(d,4);
			wij.push_back(wj);
			sum_wij += wj;
		}
		// get the interpolated value
		points_data[i] = 0;
		for(unsigned int j=0; j<shared_cells[i].size(); j++)
		{
			points_data[i] += wij[j]*cells_data[shared_cells[i][j]]/sum_wij;
		}

	}
}
// ------------------------------------------------------------------------- //
void interpl_grid::read(const char* name)
{
	ifstream in(name);
	if (!in.is_open()) 
	{
		Log << "Error opening " << name << gaplog::endl;
		exit(1);
	}
	//else
	//{
	//	Log << gaplog::endl << "Reading grid " << name << gaplog::endl << gaplog::endl;
	//}

	string line;
	
	bool check_nodes = false;
	bool check_elements = false;
	bool check_actives = false;
	bool check_pointsdata = false;
	bool check_cellsdata = false;
	
	while(getline(in,line))
	{

		// read nodes
		if(line.find("NODES") != std::string::npos)	
		{
			check_nodes = true;
			istringstream iss(line);
			string keyword;
			iss >> keyword;
			int n_nodes;
			iss >> n_nodes;			
			//Log << "  *\tFound " << n_nodes << " nodes" << gaplog::endl;
			nodes.resize(n_nodes);
			for(int i=0; i<n_nodes; i++)
			{
				getline(in,line);
				istringstream _iss(line);
				for(int j=0; j<3; j++)
					_iss >> nodes[i][j];
			}
		}
		// read elements
		if(line.find("ELEMENTS") != std::string::npos && line.find("ACTIVE") == std::string::npos)	
		{
			check_elements = true;
			istringstream iss(line);
			string keyword;
			iss >> keyword;
			int n_elements;
			iss >> n_elements;
			//Log << "  *\tFound " << n_elements << " elements" << gaplog::endl;
			elements.resize(n_elements);

			ELM_NDS = 0;
			getline(in,line);
			istringstream firstentry(line);
			while(firstentry)
			{
				string val;
				firstentry >> val;
				if(val.size() > 0)
				{
					int idx;
					istringstream _iss(val);
					_iss >> idx;
					elements[0].push_back(idx);
					ELM_NDS++;
				}
			}
			for(int i=1; i<n_elements; i++)
			{
				getline(in,line);
				istringstream _iss(line);
				for(int j=0; j<ELM_NDS; j++)
				{
					int idx;
					_iss >> idx;
					elements[i].push_back(idx);
				}
			}
		}
		// read elements
		if(line.find("ACTIVE_ELEMENTS") != std::string::npos)	
		{
			check_actives = true;
			active.resize(elements.size());
			for(int i=0; i<elements.size(); i++)
			{
				getline(in,line);
				istringstream iss(line);
				int _tmp;
				iss >> _tmp;
				active[i] = _tmp;
			}
			//Log << "done!" << gaplog::endl;
			
		}
		// read cells data
		if(line.find("CELLS_DATA") != std::string::npos)	
		{
			check_cellsdata = true;
			//Log << "  *\t Reading cells data ... ";
			cells_data.resize(elements.size());
			for(int i=0; i<elements.size(); i++)
			{
				getline(in,line);
				istringstream iss(line);
				iss >> cells_data[i];
			}
			//Log << "done!" << gaplog::endl;
		}
		// read points data
		if(line.find("POINTS_DATA") != std::string::npos)	
		{
			check_pointsdata = true;
			//Log << "  *\tReading point data ... ";
			points_data.resize(nodes.size());
			for(int i=0; i<nodes.size(); i++)
			{
				getline(in,line);
				istringstream iss(line);
				iss >> points_data[i];
			}
			//Log << "done!" << gaplog::endl;
		}
	}
	
	// cells data not specified, initialize to zero
	if(!check_cellsdata)
	{
		//Log << "  *\twarning: no cell data specified in "
		//		 << name << gaplog::endl;

		cells_data.resize(elements.size(), 0.0);
	}
	// active elements not specified, initialize to 1
	if(!check_actives)
	{
		//Log << "  *\twarning: no active elements specified in "
		//		 << name << gaplog::endl;

		active.resize(elements.size(), true);
	}
	// cells data not specified, initialize to zero
	if(!check_pointsdata)
	{
		//Log << "  *\twarning: no points data specified in "
		//		 << name << gaplog::endl;
		points_data.resize(nodes.size(), 0.0);
	}

	//Log << gaplog::endl << "done reading " << name << gaplog::endl;
}
// ------------------------------------------------------------------------- //
void interpl_grid::write(const char* name) const 
{
	ofstream out(name);
	if (!out.is_open()) 
	{
		Log << "Error opening " << name << gaplog::endl;
		exit(1);
	}
	//else
	//{
	//	cou*/t << gaplog::endl << "Writing grid " << name << gaplog::endl;
	//}

	out << "NODES " << nodes.size() << std::endl;
	for(unsigned int i = 0; i < nodes.size(); i++) 
	{
		out << nodes[i].x() << "\t" 
				<< nodes[i].y() << "\t" 
				<< nodes[i].z() << std::endl; 
	}
	out << "ELEMENTS " << elements.size() << std::endl;
	for(int i = 0; i < ne(); i++) 
	{
		for(int j = 0; j < ELM_NDS; j++)
			out << elements[i][j] << "\t"; // index start from 0
		out << std::endl;
	}
	out << "ACTIVE_ELEMENTS " << elements.size() << std::endl;
	for(int i = 0; i < ne(); i++) 
	{
		out << active[i] << std::endl;
	}
	out << "CELLS_DATA" << std::endl;
	for(int i = 0; i < ne(); i++) 
	{
		out << cells_data[i] << std::endl;
	}
	out << "POINTS_DATA" << std::endl;
	for(int i = 0; i < nn(); i++) 
	{
		out << points_data[i] << std::endl;
	}
}
// ------------------------------------------------------------------------- //
void interpl_grid::writeVTK(const char* name) const
{

	ofstream vtk(name);
	if (!vtk.is_open()) 
	{
		Log << "Error opening " << name << gaplog::endl;
		exit(1);
	}

	// write VTK header
	vtk << "# vtk DataFile Version 2.0\n"
			<< "vtk output\n"
			<< "ASCII\n"
			<< "DATASET UNSTRUCTURED_GRID\n"
			<< "POINTS " << nodes.size() << " double\n\n";
		
	// ----------------------------- mesh nodes -------------------------------//
	
	for(unsigned int i = 0; i < nodes.size(); i++) 
	{
		vtk << nodes[i].x() << "\t" 
				<< nodes[i].y() << "\t" 
				<< nodes[i].z() << std::endl; 
	}

	vtk << std::endl;


	// ------------------------ elements definition -------------------------- //

	vtk << "CELLS " << ne() << "\t" 
		<< (1 + ELM_NDS)*ne() << std::endl;

	for(int i = 0; i < ne(); i++) {
		vtk << ELM_NDS << "\t";
		for(int j = 0; j < ELM_NDS; j++)
			vtk << elements[i][j] << "\t"; // index start from 0
		vtk << std::endl;
	}
	
	vtk << std::endl;

	// --------------------------- elements type ----------------------------- //

	vtk << "CELL_TYPES " << ne() << std::endl;
	
	int ELM_TYPE;
	if(ELM_NDS == 3)
		ELM_TYPE = 5; // VTK_VOL_ELM_TYPE
	else
		ELM_TYPE = 9; // VTK_QUAD

	for (int i = 0; i < ne(); i++) 
		vtk << ELM_TYPE << std::endl;

	// ----------------------------- cells data ------------------------------ //

	vtk << "\nCELL_DATA " << ne() << std::endl;
	vtk << "\nSCALARS cellsdata float" << std::endl;
	vtk << "\nLOOKUP_TABLE default" << std::endl;

	for(int i=0; i<ne(); i++)
		vtk << cells_data[i] << std::endl;

	vtk << "\nSCALARS active float" << std::endl;
	vtk << "\nLOOKUP_TABLE default" << std::endl;

	for(int i=0; i<ne(); i++)
		vtk << active[i] << std::endl;

	
	// ----------------------------- points data ----------------------------- //

	vtk << "\nPOINT_DATA " << nn() << std::endl;
	vtk << "\nSCALARS pointsdatafield float" << std::endl;
	vtk << "\nLOOKUP_TABLE default" << std::endl;

	for(int i=0; i<nn(); i++)
		vtk << points_data[i] << std::endl;

}
// ------------------------------------------------------------------------- //
void interpl_grid::copy(interpl_grid& newgrid) const
{

	newgrid.ELM_NDS = this->ELM_NDS;

	newgrid.nodes.resize(this->nodes.size());
	newgrid.elements.resize(this->elements.size());
	newgrid.active.resize(this->active.size());
	newgrid.points_data.resize(this->points_data.size());
	newgrid.cells_data.resize(this->cells_data.size());

	for(int i=0; i<newgrid.nn(); i++)
	{
		newgrid.nodes[i] = this->nodes[i];
		newgrid.points_data[i] = this->points_data[i];
	}
	for(int i=0; i<newgrid.ne(); i++)
	{
		newgrid.elements[i] = this->elements[i];
		newgrid.cells_data[i] = this->cells_data[i];
		newgrid.active[i] = this->active[i];
	}

}
// ------------------------------------------------------------------------- //