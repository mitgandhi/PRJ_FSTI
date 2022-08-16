# include "./te_solver.h"
# include "../FSTI_Block_dll/log.h"
# include <set>

# define VTK_TETRA 10
# define VTK_TRI 5

using namespace std;
extern class gaplog Log;

// -------------------------------------------------------------------------- //
te_solver::te_solver(const input& _in, const char* _body)
	: in(_in), body(_body), htr(in,msh), thdfl(in,msh)
{

	// build the mesh
	if(string(body).compare("CB") == 0)
		msh.build(in.data.thermal.block.meshfile.c_str(), 1e-3);
	else
		msh.build(in.data.thermal.valveplate.meshfile.c_str(), 1e-3);

	Tn.resize(msh.nodes.size());
	Te.resize(msh.elements.size());
	displ.resize(3*msh.nodes.size(), 0.0);

	// define the interpolation grid for the gap surface
	define_gap_interpl();

}
// -------------------------------------------------------------------------- //
void te_solver::define_gap_interpl()
{
	
	geom_set::const_iterator gapfs;
	if(string(body).compare("CB") == 0)
		gapfs = msh.face_sets.find("gap_block");
	else
		gapfs = msh.face_sets.find("gap");

	if(gapfs == msh.face_sets.end())
	{
		Log << "\nCould not find face set gap_block in the mesh ";
		
		if(string(body).compare("CB") == 0)
			Log << in.data.thermal.block.meshfile << "\n\n";
		else
			Log << in.data.thermal.valveplate.meshfile << "\n\n";

		exit(1);
	}

	geom_set::const_iterator gapns;
	if(string(body).compare("CB") == 0)
		gapns = msh.node_sets.find("gap_block");
	else
		gapns = msh.node_sets.find("gap");
	if(gapns == msh.node_sets.end())
	{
		Log << "\nCould not find node set gap_block in the mesh ";
		
		if(string(body).compare("CB") == 0)
			Log << in.data.thermal.block.meshfile << "\n\n";
		else
			Log << in.data.thermal.valveplate.meshfile << "\n\n";

		exit(1);
	}

	// initialize the gap interpolation structure
	gap.initialize(gapfs->second.size(), gapns->second.size(), 3);
	
	// define the nodes
	map<int,int> g2l;
	for(unsigned int i=0, c=0; i<gapns->second.size(); i++)
	{
		g2l[gapns->second[i]] = c++;
		gap.nodes[i] = msh.nodes[gapns->second[i]];
	}
	// define the elements
	for(unsigned int i=0; i<gapfs->second.size(); i++)
	{
		for(int j=0; j<3; j++)
			gap.elements[i][j] = g2l.at(msh.faces[gapfs->second[i]].nodes[j]);
		
		gap.active[i] = true;
	}

}
// -------------------------------------------------------------------------- //
void te_solver::solve_htr()
{
	// copy the gap heat flux from the inteprolation grid
	htr.qgap = gap.cells_data;

	// solve for the T field
	if(string(body).compare("CB") == 0)
		htr.initialize("CB");
	else
		htr.initialize("VP");
	
	// calculate the thermal stiffness matrix
	htr.get_K();

	// solve the system
	htr.solve();

	// copy the temperature @ nodes
	for(unsigned int i=0; i<msh.nodes.size(); i++)
		Tn[i] = htr.sol[i];
	
	// interpolate the temperature @ elements
	for(unsigned int i=0; i<msh.elements.size(); i++)
	{
		Te[i] = 0;
		for(int j=0; j<4; j++)
			Te[i] += 0.25*htr.sol[msh.elements[i].nodes[j]];
	}

	// update the interpolation structure
	geom_set::const_iterator it;
	if(string(body).compare("CB") == 0)
		it = msh.face_sets.find("gap_block");
	else
		it = msh.face_sets.find("gap");

	const vector<int>& fs = it->second;
	for(unsigned int i=0; i<fs.size(); i++)
	{
		double Tf = 0;
		for(int j=0; j<3; j++)
			Tf += Tn[msh.faces[fs[i]].nodes[j]]/3.0;
		// update the cells data with the calculated temperature
		gap.cells_data[i] = Tf;
	}

	// clear the htr object to save memory
	htr.clear();

}
// -------------------------------------------------------------------------- //
void te_solver::solve_thdfl()
{
	// assign pointer to the elements center
	thdfl.Te = &Te;

	// solve for the T field
	if(string(body).compare("CB") == 0)
		thdfl.initialize("CB");
	else
		thdfl.initialize("VP");

	// calculate the stiffness matrix
	thdfl.get_K();

	// now calculate and apply the thermal loads
	thdfl.calc_thloads();

	// apply ir
	thdfl.apply_ir();

	// solve the linear system
	thdfl.solve();

	// map the displacement on the original mesh nodes
	displ.resize(3*msh.nodes.size(), 0.0);
	for(int i=0; i<thdfl.smsh.nn(); i++)
	{
		displ[3*thdfl.smsh.gnid[i]] = thdfl.sol[3*i];
		displ[3*thdfl.smsh.gnid[i]+1] = thdfl.sol[3*i+1];
		displ[3*thdfl.smsh.gnid[i]+2] = thdfl.sol[3*i+2];
	}

	// update the interpolation structure
	geom_set::const_iterator it;
	if(string(body).compare("CB") == 0)
		it = msh.node_sets.find("gap_block");
	else
		it = msh.node_sets.find("gap");

	const vector<int>& ns = it->second;
	for(unsigned int i=0; i<ns.size(); i++)
		gap.points_data[i] = displ[3*ns[i]+2];

}
// -------------------------------------------------------------------------- //
void te_solver::write_vtk(const char* path)
{

	ostringstream oss;
	if(string(body).compare("CB") == 0)
		oss << path << "/cylinderblock.vtk";
	else
		oss << path << "/valveplate.vtk";
	
	ofstream vtk(oss.str());

	// write VTK header
	vtk << 
		"# vtk DataFile Version 2.0" << endl <<
		"vtk output" << endl <<
		"ASCII" << endl <<
		"DATASET UNSTRUCTURED_GRID" << endl << 
		"POINTS " << msh.nodes.size() << " float" << endl;
	
	// ----------------------------- mesh nodes -------------------------------//
	
	for(unsigned int i = 0; i < msh.nodes.size(); i++) 
	{
		vtk << msh.nodes[i].x() << "\t" 
				<< msh.nodes[i].y() << "\t" 
				<< msh.nodes[i].z() << endl; 
	}

	vtk << endl;

	// ------------------------ elements definition -------------------------- //

	vtk << "CELLS " << msh.elements.size() << "\t" 
			<< (1 + 4)*msh.elements.size() << endl;

	for(int i = 0; i < msh.elements.size(); i++) 
	{
		vtk << 4 << "\t";
		for(unsigned int j=0; j<4; j++)
			vtk << msh.elements[i].nodes[j] << "\t"; // index start from 0
		vtk << endl;
	}
	
	vtk << endl;

	// --------------------------- elements type ----------------------------- //

	vtk << "CELL_TYPES " << msh.elements.size() << endl;
			
	for(int i=0; i<msh.elements.size(); i++) 
	{
		vtk << VTK_TETRA << endl;
	}

	vtk << "CELL_DATA " << msh.elements.size() << endl;
	
	vtk << "SCALARS sets float 1" << endl;
	vtk << "LOOKUP_TABLE default" << endl;
	for(int i=0; i<msh.elements.size(); i++)
	{
		string seti = msh.elements[i].elm_set;
		geom_set::iterator it = msh.elm_sets.find(seti.c_str());
		vtk << distance(msh.elm_sets.begin(), it) << endl;
	}

	vtk << "SCALARS T float 1" << endl;
	vtk << "LOOKUP_TABLE default" << endl;
	for(int i=0; i<msh.elements.size(); i++)
			vtk << Te[i] << endl;

	vtk << "POINT_DATA " << msh.nodes.size() << endl;
	vtk << "SCALARS T float 1" << endl;
	vtk << "LOOKUP_TABLE default" << endl;

	for(int i=0; i<msh.nodes.size(); i++)
			vtk << Tn[i] << endl;

	vtk << "VECTORS displ float" << endl;

	
	for(int i=0; i<msh.nodes.size(); i++)
	{
		vtk << displ[3*i] << "\t" 
				<< displ[3*i+1] << "\t"
				<< displ[3*i+2] << endl;
	}

}
// ------------------------------------------------------------------------- //
void te_solver::write_fset_vtk(const char* path, const char* fset)
{

	geom_set::const_iterator it = msh.face_sets.find(fset);
	if(it ==  msh.face_sets.end())
	{
		Log << "\nCould not find face set " << fset << "\n\n";
		exit(1);
	}

	const vector<int>& fs = it->second;

	set<int> _nodes;
	for(unsigned int i=0; i<fs.size(); i++)
	{
		for(int j=0; j<3; j++)
			_nodes.insert(msh.faces[fs[i]].nodes[j]);
	}

	map<int,int> g2l;
	set<int>::const_iterator itt;
	int c=0;
	for(itt=_nodes.begin(); itt!=_nodes.end(); itt++)
		g2l[*itt] = c++;


	ostringstream oss;
	if(string(body).compare("CB") == 0)
		oss << path << "/cylinderblock." << fset << ".vtk";
	else
		oss << path << "/valveplate." << fset << ".vtk";

	ofstream vtk(oss.str().c_str());

	// write VTK header
	vtk << 
		"# vtk DataFile Version 2.0" << endl <<
		"vtk output" << endl <<
		"ASCII" << endl <<
		"DATASET UNSTRUCTURED_GRID" << endl << 
		"POINTS " << _nodes.size() << " float" << endl;
	
	// ----------------------------- mesh nodes -------------------------------//
	
	for(itt=_nodes.begin(); itt!=_nodes.end(); itt++)
	{
		vtk << msh.nodes[*itt].x() << "\t" 
				<< msh.nodes[*itt].y() << "\t" 
				<< msh.nodes[*itt].z() << endl; 
	}

	vtk << endl;

	// ------------------------ elements definition -------------------------- //

	vtk << "CELLS " << fs.size() << "\t" 
			<< (1 + 3)*fs.size() << endl;

	for(unsigned int i=0; i<fs.size(); i++) 
	{
		vtk << 3 << "\t";
		for(unsigned int j=0; j<3; j++)
			vtk << g2l[msh.faces[fs[i]].nodes[j]] << "\t";
		vtk << endl;
	}
	
	vtk << endl;

	// --------------------------- elements type ----------------------------- //

	vtk << "CELL_TYPES " << fs.size() << endl;
			
	for(int i=0; i<fs.size(); i++) 
	{
		vtk << VTK_TRI << endl;
	}

	vtk << "CELL_DATA " << fs.size() << endl;
	vtk << "SCALARS T float 1" << endl;
	vtk << "LOOKUP_TABLE DEFAULT" << endl;
	for(unsigned int i=0; i<fs.size(); i++)
	{
		double Tf = 0;
		for(int j=0; j<3; j++)
			Tf += Tn[msh.faces[fs[i]].nodes[j]]/3.0;
		vtk << Tf << endl;
	}

	vtk << "POINT_DATA " << _nodes.size() << endl;
	vtk << "SCALARS T float 1" << endl;
	vtk << "LOOKUP_TABLE DEFAULT" << endl;

	map<int,int>::const_iterator mit;
	for(mit=g2l.begin(); mit!=g2l.end(); mit++)
		vtk << Tn[mit->first] << endl;

	vtk << "VECTORS displ float" << endl;
	for(mit=g2l.begin(); mit!=g2l.end(); mit++)
	{
		vtk << displ[3*mit->first] << "\t"
				<< displ[3*mit->first+1] << "\t"
				<< displ[3*mit->first+2] << endl;
	}

}
// ------------------------------------------------------------------------- //