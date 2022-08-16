# include "./thdfl_analysis.h"
# include "../FSTI_Block_dll/log.h"
# include <array>
# include <set>
# include <ANN\ANN.h>

using namespace std;

extern class gaplog Log;

// ------------------------------------------------------------------------- //
thdfl_analysis::thdfl_analysis(const input& _inp, const tetra_mesh&	_msh) 
	: inp(_inp), msh(_msh), smsh(_msh), solver(&K,&b)
{
	ndof = 0;
	ndof_eff = 0;
	Te = 0;
}
// ------------------------------------------------------------------------- //
void thdfl_analysis::initialize(const char* _body)
{


	body = _body;


	// define the mesh subset

	if(string(body).compare("CB") == 0)
	{

		Log << "\nInitialization of cylinder block deflection analysis\n\n";

		// create a "subset" containing the whole mesh
		vector<string> sets;
		geom_set::const_iterator it;
		for(it=msh.elm_sets.begin(); it!=msh.elm_sets.end(); it++)
			sets.push_back(it->first);

		smsh.build(sets);

	}
	else // valveplate
	{

		Log << "\nInitialization of valve plate deflection analysis\n\n";

		vector<string> sets;
		// create a "subset" containing the whole mesh
		if(inp.data.thermal.valveplate.calc_deflect_on.front().compare("ALL") == 0)
		{
			geom_set::const_iterator it;
			for(it = msh.elm_sets.begin(); it!=msh.elm_sets.end(); it++)
				sets.push_back(it->first);
			smsh.build(sets);
		}
		else
		{
			Log << "  * Defining mesh subset at ";
			for(unsigned int i=0; i<inp.data.thermal.valveplate.calc_deflect_on.size(); i++)
				Log << inp.data.thermal.valveplate.calc_deflect_on[i] << " ";
			Log << " ... ";

			// build the mesh subset
			smsh.build(inp.data.thermal.valveplate.calc_deflect_on);

			Log << "done!" << gaplog::endl;

			Log << "  * Defined " << smsh.ne() << " elements and " 
					 << smsh.nn() << " nodes." << gaplog::endl;
		}
	}

	ndof = 3*smsh.nn();
	ndof_eff = ndof;
	gdof.resize(ndof, 0);
	sol.resize(ndof, 0.0);

	// create the list of elements
	for(unsigned int i=0; i<smsh.ne(); i++)
	{
		th_tetra_3dof t(i,smsh);
		elements.push_back(t);
	}

	// ------------------------- assign the materials ------------------------- //

	Log << "  * assigning material ... ";

	if(string(body).compare("CB") == 0)
	{
		// assign the materials
		for(unsigned int i=0; i<elements.size(); i++)
		{
			bool found = false;
			// look into set_name -> material_name list
			for(unsigned int j=0; j<inp.data.thermal.block.materials.size(); j++)
			{
				// elm set and set associated to current material are the same
				if(msh.elements[i].elm_set.compare(inp.data.thermal.block.materials[j].first) == 0)
				{
					// find the material in the materials list
					for(unsigned int k=0; k<inp.data.materials.size(); k++)
					{
						if(inp.data.materials[k].name.compare(inp.data.thermal.block.materials[j].second) == 0)
						{
							elements[i].m = &inp.data.materials[k];
							found = true;
							break;
						}
					}
				}

				if(found)
					break;
			}

			if(!found)
			{
				Log << "\nCould not find material associated with element " << i
						 << "Please fix the mesh.\n\n";
				exit(1);
			}
		}
	}
	else
	{
		// assign the materials
		for(unsigned int i=0; i<elements.size(); i++)
		{
			bool found = false;
			// look into set_name -> material_name list
			for(unsigned int j=0; j<inp.data.thermal.valveplate.materials.size(); j++)
			{
				// elm set and set associated to current material are the same
				if(msh.elements[i].elm_set.compare(inp.data.thermal.valveplate.materials[j].first) == 0)
				{
					// find the material in the materials list
					for(unsigned int k=0; k<inp.data.materials.size(); k++)
					{
						if(inp.data.materials[k].name.compare(inp.data.thermal.valveplate.materials[j].second) == 0)
						{
							elements[i].m = &inp.data.materials[k];
							found = true;
							break;
						}
					}
				}

				if(found)
					break;
			}

			if(!found)
			{
				Log << "\nCould not find material associated with element " << i
						 << "Please fix the mesh.\n\n";
				exit(1);
			}
		}
	}
	
	Log << "done!" << gaplog::endl;

	
	// define the global dof numbering 
	for(unsigned int i=0, id=0; i<ndof; i++)
	{
		if(gdof[i] > -1)
			gdof[i] = id++;
		else
			ndof_eff--;
	}

	// resize the stiffness matrix accounting for the real size of the problem
	// 3 dof are added because inertia relief will be used
	K.resize(ndof_eff+3, ndof_eff+3);
	// resize b to proper size
	b.resize(ndof_eff+3, 0.0);


}
// ------------------------------------------------------------------------- //
void thdfl_analysis::calc_thloads()
{

	Log << "  * calculating thermal loads ... ";
	
	// reference temperature
	double Tref = 20.0;

	th_load.resize(ndof_eff, 0.0);

	for(int i=0; i<smsh.ne(); i++)
	{
		// global element id and reference
		int gid = smsh.geid[i];
		const tetra& elm = smsh.get_elm(i);
		
		// element volume
		double V = elm.volume;

		// coefficient of linear expansion
		double alpha = elements[i].m->alpha;

		// deflection matrix
		matrix epsilon(6,1);
		epsilon[0][0] = alpha*(Te->at(gid) - Tref);			// epsilon x
		epsilon[1][0] = alpha*(Te->at(gid) - Tref);			// epsilon y
		epsilon[2][0] = alpha*(Te->at(gid) - Tref);			// epsilon z
		epsilon[3][0] = 0;															// gamma xy
		epsilon[4][0] = 0;															// gamma xz
		epsilon[5][0] = 0;															// gamma yz

		matrix D = elements[i].calc_D();
		matrix B = elements[i].calc_B();

		matrix Q = V*(B.T()*D)*epsilon;

		// sum the load to the nodes
		for(int j=0, id=0; j<4; j++)
		{
			for(int k=0; k<3; k++, id++)
				th_load[elements[i].gnum[id]] += Q[3*j + k][0];
		}

	}

	for(unsigned int i=0; i<th_load.size(); i++)
		b[i] += th_load[i];

	Log << "done!" << gaplog::endl;

}
// ------------------------------------------------------------------------- //
int thdfl_analysis::calc_inertial_load()
{

	// Size the nodal mass matrix
	node_mass.resize(smsh.nn(), 0.0);
	// resize IR loads
	ir_loads.resize(3*smsh.nn(), 0);

	// Find the mass, total mass, nodal mass, and COG	
	vector<double> M(smsh.ne());
	double Mtot = 0;
	double Vtot = 0;
	double xcg = 0, ycg = 0, zcg = 0;
	
	for(unsigned int c=0; c<smsh.ne(); c++)
	{
		// reference to the element
		const tetra& elm = smsh.get_elm(c);

		// Find cell mass
		double rho = elements[c].m->rho;
		M[c] = elm.volume*rho;

		// Set the nodal mass
		for(int i=0; i<4; i++)
			node_mass[ smsh.g2l_nn[ elm.nodes[i] ] ] += (M[c]/4.0);

		// COG
		xcg += M[c]*elm.center.x();
		ycg += M[c]*elm.center.y();
		zcg += M[c]*elm.center.z();
		
		// Mass sum
		Mtot += M[c];
		Vtot += elm.volume;
	}
	
	xcg /= Mtot;
	ycg /= Mtot;
	zcg /= Mtot;

	// Calculate the inertia tensor of each element and body inertia tensor
	vector<matrix> I(smsh.ne());
	matrix Itot(3, 3);
	Itot = 0;

	for(unsigned int c=0; c<smsh.ne(); c++)
	{
		// reference to the element
		const tetra& elm = smsh.get_elm(c);

		//Find the cell centroid position vector from the COG
		double rx = elm.center.x() - xcg;
		double ry = elm.center.y() - ycg;
		double rz = elm.center.z() - zcg;

		// inertial tensor for the current element
		matrix Ie(3,3);
		Ie[0][0] = M[c]*(ry*ry + rz*rz);		// Ixx
		Ie[1][1] = M[c]*(rx*rx + rz*rz);		// Iyy
		Ie[2][2] = M[c]*(rx*rx + ry*ry); 		// Izz
		Ie[0][1] = Ie[1][0] = -M[c]*rx*ry;	// Ixy = Iyx
		Ie[0][2] = Ie[2][0] = -M[c]*rx*rz;	// Ixz = Izx
		Ie[1][2] = Ie[2][1] = -M[c]*ry*rz;	// Iyz = Izy
		
		I[c] = Ie;
		Itot = Itot + Ie;
	}
	
	// Calculate net force & torque about COG
	matrix F(3, 1);
	matrix T(3, 1);

	F = 0, T = 0;

	for(unsigned int l=0; l<th_load.size(); l++)
	{
		
		int d = l % 3;		// What direction (x,y,z) does this force act in?
		matrix f(3,1);
		f[d][0] += th_load[l];			// load
		matrix r(3,1);
		r[0][0] =	smsh.get_node(l/3).x() - xcg;
		r[1][0] = smsh.get_node(l/3).y() - ycg;
		r[2][0] = smsh.get_node(l/3).z() - zcg;

		F = F + f;						// add to the applied force
		T = T + cross(r, f);	// add to the applied torque
		
	}

	Log << "  * Solid Body analysis ..." << gaplog::endl << gaplog::endl;
	Log << "     * V = " << Vtot << " m3" << gaplog::endl;
	Log << "     * M = " << Mtot << " kg" << gaplog::endl;
	Log << "     * I = [ " << Itot[0][0] << "\t" << Itot[0][1] << "\t" << Itot[0][2] << gaplog::endl
			 << "             " << Itot[1][0] << "\t" << Itot[1][1] << "\t" << Itot[1][2] << gaplog::endl
		   << "             " << Itot[2][0] << "\t" << Itot[2][1] << "\t" << Itot[2][2] 
			 << " ] kgm^2" << gaplog::endl;
	Log << "     * COG = [ " << xcg << "\t" << ycg << "\t" << zcg << " ] m" << gaplog::endl;
	
	Log << "     * Forces and moment imbalance: " << gaplog::endl;
	Log << "       * dF = [ " << scientific << F[0][0] 
			 << "\t" << scientific << F[1][0] 
			 << "\t" << scientific << F[2][0] << " ] N" << gaplog::endl;

	Log << "       * dM = [ " << scientific << T[0][0] 
			 << "\t" << scientific << T[1][0] 
			 << "\t" << scientific << T[2][0] << " ] Nm" << gaplog::endl;
	
	double dF_norm = sqrt(pow(F[0][0],2) + pow(F[1][0],2) + pow(F[2][0],2));
	double dT_norm = sqrt(pow(T[0][0],2) + pow(T[1][0],2) + pow(T[2][0],2));

	if(dF_norm < 1e-9 && dT_norm < 1e-9)	// no need to apply inertial loads
	{
		
		Log << "\n    No need to apply inertial loads." << gaplog::endl;
		
		return -1;
	}

	// Linear acceleration
	matrix A_lin = -1.0*F/Mtot;
	
	// Angular acceleration
	matrix A_ang = -1.0*Itot.inv()*T;
	
	Log << "     * Lin. acc. = [ " 
			 << A_lin[0][0] << "\t" << A_lin[1][0] << "\t" << A_lin[2][0] 
		   << " ] m/s^2" << gaplog::endl;
	
	Log << "     * Ang. acc. = [ " 
			 << A_ang[0][0] << "\t" << A_ang[1][0] << "\t" << A_ang[2][0] 
			 << " ] rad/s^2" << gaplog::endl;

	Log << "     * Applying inertial loads ... ";


	// Calculate the inertia load on each cell
	for(unsigned int c=0; c<smsh.ne(); c++)
	{
		
		// reference to the element
		const tetra& elm = smsh.get_elm(c);

		// Cell force
		matrix F_lin = M[c]*A_lin;
		matrix T_cell = I[c]*A_ang;

		// Cell centroid position vector
		matrix r(3,1);
		r[0][0] = elm.center.x() - xcg;
		r[1][0] = elm.center.y() - ycg;
		r[2][0] = elm.center.z() - zcg;
		double r_norm2 = pow(r[0][0],2) + pow(r[1][0],2) + pow(r[2][0],2);
		
		matrix F_ang = (1.0/r_norm2)*cross(T_cell, r);
		matrix F_tot = (F_lin + F_ang)/4.0;

		// Apply the inertial reaction force to the cell nodes
		for(int i=0; i<4; i++)
		{
			int lnid = smsh.g2l_nn[elm.nodes[i]];
			ir_loads[3*lnid] += F_tot[0][0];
			ir_loads[3*lnid + 1] += F_tot[1][0];
			ir_loads[3*lnid + 2] += F_tot[2][0];
		}	
	}

	// Calculate net force & torque about COG now that the inertial reaction loads have been applied
	// Just for checking, not actually needed
	
	// sum the inertial loads
	for(unsigned int i=0; i<smsh.nn(); i++)
	{
		matrix r(3,1);
		r[0][0] = smsh.get_node(i).x() - xcg;
		r[1][0] = smsh.get_node(i).y() - ycg;
		r[2][0] = smsh.get_node(i).z() - zcg;
		matrix f(3,1);
		for(int j=0; j<3; j++)
		{
			f[j][0] += ir_loads[3*i + j];
		}
		
		F = F + f;						// add to the applied force
		T = T + cross(r, f);	// add to the applied torque
	}
		
	Log << "done!" << gaplog::endl;

	Log << "     * Forces and moment imbalance: " << gaplog::endl;
	
	Log << "       * dF = [ " << scientific << F[0][0] 
			 << "\t" << scientific << F[1][0] 
			 << "\t" << scientific << F[2][0] << " ] N" << gaplog::endl;

	Log << "       * dM = [ " << scientific << T[0][0] 
			 << "\t" << scientific << T[1][0] 
			 << "\t" << scientific << T[2][0] << " ] Nm" << gaplog::endl;

	Log << gaplog::endl;


	return 0;

}
// ------------------------------------------------------------------------- //
int thdfl_analysis::apply_ir()
{
	// calculate the inertial loads
	calc_inertial_load();

	// ----------------- update the load vector ----------------- //

	for(unsigned int i=0; i<ir_loads.size(); i++)
		b[i] += ir_loads[i];

	//  ----------------- update the stiffness matrix ----------------- //

	double mass_avg = 0;
	for(int i=0; i<node_mass.size(); i++)
		mass_avg += node_mass[i];
	mass_avg /= node_mass.size();

	// Scale factor = [abs mean stiffness matrix] / [mean node mass]
	//double scale = K_avg/mass_avg;
	double scale = K_avg;

	for(int i=0; i<smsh.nn(); i++)
	{
		for(int d=0; d<3; d++)
		{
			int row = ndof_eff + d;
			int col = gdof[3*i + d];
			K(row, col) = scale;
			K(col, row) = scale;
		}
	}

	return 0;
}
// ------------------------------------------------------------------------- //
void thdfl_analysis::get_K()
{
	// set to zero the content of the stiffness matrix!
	gmm::clear(K);

	// number of non zero value in stiffness matrix
	int nnz = 0;
	K_avg = 0;

	// ----------- this is for the operation progress information ------------ //
	
	double single_step = 1.0; // percentage
	double next_step = single_step; // percentage

	std::streamsize default_val = cout.precision();

	cout.precision(1);
	Log << "  * Assembling the stiffness matrix " << std::fixed << 0.00 << " %";

	// ---------------------- loop to all the elements  ---------------------- //

	for(int e=0; e<smsh.ne(); e++)
	{
		
		// ---------------- gives information on operation progress ------------ //
		double current = (double) e/smsh.ne();
		if(current > next_step/100) 
		{
			next_step += single_step;
			Log << "\r  * Assembling the stiffness matrix " 
					 << std::fixed << 100*current << " %";
		}

		// get the element stiffness matrix
		matrix ke = elements[e].calc_K();

		// imagine to have the ideal stiffness matrix Kid without any constraint;
		// the size of Kid would be [ndof*nn x ndof*nn]
		// think about the actual stiffness matrix as a sub-matrix of Kid
		// When non zero Dirichlet th_boundary conditions are present the load vector
		// need to be corrected using the non constrained rows and the constrained 
		// columns of Kid and the known solution values. Of course, the approach works 
		// also for zero dirichlet th_boundaries 

		for(unsigned int i=0; i<elements[e].gnum.size(); i++)
		{
			if(gdof[elements[e].gnum[i]] > -1)	// row elements[e].gnum[i] not constraint
			{
				for(unsigned int j=0; j<elements[e].gnum.size(); j++)
				{
					if(gdof[elements[e].gnum[j]] > -1) // column elements[e].gnum[j] not constraint
					{
						K(gdof[elements[e].gnum[i]], gdof[elements[e].gnum[j]]) += ke[i][j];
						K_avg += fabs(ke[i][j]);
						nnz++;
					}
					// constrained columns, correct the known vector at the same position of the row id
					// (gdof[elm->gnum[i]]) withe the product between 
					// the element stiffness and the known solution in the position 
					// of the coulumn (elm->gnum[j])
					else 
					{
						b[gdof[elements[e].gnum[i]]] -= ke[i][j]*sol[elements[e].gnum[j]];
					}
				}
			}
		}
						
	}

	Log << "\r  * Assembling the stiffness matrix 100.00 % (nnz = " 
			 << nnz << ")" << gaplog::endl;



	cout.precision(default_val);
}
// ------------------------------------------------------------------------- //
void thdfl_analysis::solve()
{
	
	Log << "\nSolving the linear system ... ";
	
	solver.system.quiet = true;
	solver.system.tol = 1e-6;
	solver.system.itermax = 2000;
	solver.system.P = ILDLT;
	solver.update();

	// clear the sriffness matrix, since has been copied into
	// the solver structure 
	K.clear_mat();
	K.resize(0,0);
	
	solver.initialize_x();
	solver.preconditioner();
	int iters = solver.solve();

	Log << "done!" << gaplog::endl << "Final residual " << scientific << solver.residual() 
			 << " in " << iters << " iterations." << gaplog::endl;

	// copy from the solution vector only where is required
	// dirichlet th_boundary values are of course skipped
	for(unsigned int i=0; i<ndof; i++)
	{	
		if(gdof[i] > -1)
			sol[i] = solver.x[gdof[i]]; 
	}

	// the lambda constraint method is used
	if(ndof_eff < b.size())
	{
		
		Log << "\n\n    Lambda displ. = [ " << solver.x[solver.x.size()-3] << "\t"
									 << solver.x[solver.x.size()-2] << "\t"
									 << solver.x[solver.x.size()-1] << " ]\n" << gaplog::endl;
	}

	
}
/*
// ------------------------------------------------------------------------- //
void thdfl_analysis::write_gap_defl()
{

	ofstream out("./data_exchange/gap_defl.csv");

	geom_set::const_iterator it = msh.node_sets.find("gap");
	const vector<int>& gapnds = it->second;

	//double min = 1e10;
	//for(unsigned int i=0; i<gapnds.size(); i++)
	//	min = (thdfl[3*gapnds[i]+2] < min) ? thdfl[3*gapnds[i]+2] : min;

	for(unsigned int i=0; i<gapnds.size(); i++)
	{
		out << msh.n2o.at(gapnds[i]) << ", "
				<< 1e3*(thdfl[3*gapnds[i]+2]) << "," << gaplog::endl;	// write in mm
	}

	out.close();
}
// ------------------------------------------------------------------------- //
*/