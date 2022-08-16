# include "./htr_analysis.h"
# include "../FSTI_Block_dll/log.h"
# include <ANN/ANN.h>
# include <set>
# include <map>

using namespace std;

extern class gaplog Log;

# define VTK_TETRA 10

// ------------------------------------------------------------------------- //
htr_analysis::htr_analysis(const input& _in, const tetra_mesh& _msh) 
	: in(_in), msh(_msh), solver(&K,&b)
{
	ndof = 0;
	ndof_eff = 0;
}
// ------------------------------------------------------------------------- //
const th_boundary& htr_analysis::get_face_bnd(int fid)
{

	// get the name of the face set
	string set = msh.faces[fid].face_set;

	const geom_set::const_iterator it = msh.face_sets.find(set.c_str());
	if(it == msh.face_sets.end())
	{
		Log << "\nth_analysis::get_face_bnd() face " << fid 
					<< " is not associated with any face set in the mesh!" << gaplog::endl;
		exit(1);
	}

	for(unsigned int b=0; b<bcs.size(); b++)
	{
		if(bcs[b].setname.compare(set) == 0)
		return bcs[b];
	}
	

	Log << "\nFace " << fid << " is not associated to any th_boundary condition!"
			 << gaplog::endl;
	
	exit(1);
}
// ------------------------------------------------------------------------- //
void htr_analysis::initialize(const char* _body) 
{

	// ------------------- copy from FSTI input ---------------------- //

	int z = in.data.operating_conditions.npistons;
	double Tcase = in.data.operating_conditions.T_Leak;

	// --------------------------------------------------------------- //

	body = _body;

	ndof = msh.nodes.size();
	ndof_eff = ndof;
	gdof.resize(ndof, 0);
	sol.resize(ndof, 0.0);

	// create the list of elements
	for(unsigned int i=0; i<msh.elements.size(); i++)
	{
		th_tetra_1dof t(i,msh);
		elements.push_back(t);
	}

	// copy bc list from FSTI input for CB analysis
	if(string(body).compare("CB") == 0)
	{

		Log << "\nInitialization of cylinder block thermal analysis\n\n";

		bcs.resize(0);

		// copy dirichlet boundaries
		for(unsigned int i=0; i<in.data.thermal.block.dirichlet_bc.size(); i++)
		{
			th_boundary thbc;
			thbc.bctype = dirichlet;
			thbc.setname = in.data.thermal.block.dirichlet_bc[i].setname;
			thbc.dirichlet.Ts = in.data.thermal.block.dirichlet_bc[i].Tp;
			bcs.push_back(thbc);
		}
		// copy neumann boundaries
		for(unsigned int i=0; i<in.data.thermal.block.neumann_bc.size(); i++)
		{
			th_boundary thbc;
			thbc.bctype = neumann;
			thbc.setname = in.data.thermal.block.neumann_bc[i].setname;
			thbc.neumann.q = in.data.thermal.block.neumann_bc[i].q;
			bcs.push_back(thbc);
		}
		// copy mixed boundaries
		for(unsigned int i=0; i<in.data.thermal.block.mixed_bc.size(); i++)
		{
			th_boundary thbc;
			thbc.bctype = mixed;
			thbc.setname = in.data.thermal.block.mixed_bc[i].setname;
			thbc.mixed.h = in.data.thermal.block.mixed_bc[i].h;
			thbc.mixed.Tinf = in.data.thermal.block.mixed_bc[i].Tinf;
			bcs.push_back(thbc);
		}

		// add the gap heat flux
		th_boundary thbc;
		thbc.bctype = neumann;
		thbc.setname = "gap_block";
		thbc.neumann.q = qgap;
		bcs.push_back(thbc);

		// add the heat fluxes to gap_piston
		// the file is located in ./output/piston/piston_flux.txt
		// if the file is present, it will contain just one number:
		// the most updated average heat flux over the total piston
		// gap surface
		//
		ifstream piston_flux("./output/piston/piston_flux.txt");
		if(piston_flux.is_open())
		{
			vector<double> q_piston(1);
			piston_flux >> q_piston[0];
			for(int i=1; i<=z; i++)
			{
				ostringstream oss;
				oss << "gap_piston_" << i;
				th_boundary thbc;
				thbc.bctype = neumann;
				thbc.setname = oss.str();
				thbc.neumann.q = q_piston;
				bcs.push_back(thbc);
			}
			piston_flux.close();
		}
		else
		{
			// the file is not there, just apply the case temperature
			for(int i=1; i<=z; i++)
			{
				ostringstream oss;
				oss << "gap_piston_" << i;
				th_boundary thbc;
				thbc.bctype = mixed;
				thbc.setname = oss.str();
				thbc.mixed.h = 5000;
				thbc.mixed.Tinf = Tcase;
				bcs.push_back(thbc);
			}
		}
		
		// assign the materials
		for(unsigned int i=0; i<elements.size(); i++)
		{
			bool found = false;
			// look into set_name -> material_name list
			for(unsigned int j=0; j<in.data.thermal.block.materials.size(); j++)
			{
				// elm set and set associated to current material are the same
				if(msh.elements[i].elm_set.compare(in.data.thermal.block.materials[j].first) == 0)
				{
					// find the material in the materials list
					for(unsigned int k=0; k<in.data.materials.size(); k++)
					{
						if(in.data.materials[k].name.compare(in.data.thermal.block.materials[j].second) == 0)
						{
							elements[i].m = &in.data.materials[k];
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
	// copy bc list from FSTI input for VP analysis
	else
	{
		Log << "\nInitialization of valve plate thermal analysis\n\n";

		bcs.resize(0);

		// copy dirichlet boundaries
		for(unsigned int i=0; i<in.data.thermal.valveplate.dirichlet_bc.size(); i++)
		{
			th_boundary thbc;
			thbc.bctype = dirichlet;
			thbc.setname = in.data.thermal.valveplate.dirichlet_bc[i].setname;
			thbc.dirichlet.Ts = in.data.thermal.valveplate.dirichlet_bc[i].Tp;
			bcs.push_back(thbc);
		}
		// copy neumann boundaries
		for(unsigned int i=0; i<in.data.thermal.valveplate.neumann_bc.size(); i++)
		{
			th_boundary thbc;
			thbc.bctype = neumann;
			thbc.setname = in.data.thermal.valveplate.neumann_bc[i].setname;
			thbc.neumann.q = in.data.thermal.valveplate.neumann_bc[i].q;
			bcs.push_back(thbc);
		}
		// copy mixed boundaries
		for(unsigned int i=0; i<in.data.thermal.valveplate.mixed_bc.size(); i++)
		{
			th_boundary thbc;
			thbc.bctype = mixed;
			thbc.setname = in.data.thermal.valveplate.mixed_bc[i].setname;
			thbc.mixed.h = in.data.thermal.valveplate.mixed_bc[i].h;
			thbc.mixed.Tinf = in.data.thermal.valveplate.mixed_bc[i].Tinf;
			bcs.push_back(thbc);
		}

		// add the gap heat flux
		th_boundary thbc;
		thbc.bctype = neumann;
		thbc.setname = "gap";
		thbc.neumann.q = qgap;
		bcs.push_back(thbc);

		// assign the materials
		for(unsigned int i=0; i<elements.size(); i++)
		{
			bool found = false;
			// look into set_name -> material_name list
			for(unsigned int j=0; j<in.data.thermal.valveplate.materials.size(); j++)
			{
				// elm set and set associated to current material are the same
				if(msh.elements[i].elm_set.compare(in.data.thermal.valveplate.materials[j].first) == 0)
				{
					// find the material in the materials list
					for(unsigned int k=0; k<in.data.materials.size(); k++)
					{
						if(in.data.materials[k].name.compare(in.data.thermal.valveplate.materials[j].second) == 0)
						{
							elements[i].m = &in.data.materials[k];
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


	// -------------------- settle the dirichlet th_boundaries ------------------ //

	for(unsigned int bnd = 0; bnd<bcs.size(); bnd++)
	{
		if(bcs[bnd].bctype == dirichlet)
		{
			// find the node set
			const string& sname(bcs[bnd].setname);
			const geom_set::const_iterator nsit = msh.node_sets.find(sname.c_str());
			const vector<int>& ns = nsit->second;
			for(unsigned int i=0; i<ns.size(); i++)
			{
				gdof[ns[i]] = -1;
				sol[ns[i]] = bcs[bnd].dirichlet.Ts;
			}
		}
	}

	// define the global dof numbering 
	for(unsigned int i=0, id=0; i<ndof; i++)
	{
		if(gdof[i] > -1)
			gdof[i] = id++;
		else
			ndof_eff--;
	}

	// resize the stiffness matrix accounting for the real size of the problem
	K.resize(ndof_eff, ndof_eff);
	// resize b to proper size
	b.resize(ndof_eff, 0.0);

	// --------------------------- define the loads --------------------------- //

	Log << "  * adding thermal loads\n" << gaplog::endl;
	for(unsigned int i=0; i<bcs.size(); i++)
	{
		if(bcs[i].bctype == mixed)
		{
			Log << "      - adding mixed th_boundary on set " 
					 << bcs[i].setname << " ... ";

			th_load* thisload = 
				new convection
				(
					*this,
					bcs[i].setname.c_str(),
					bcs[i].mixed.h,
					bcs[i].mixed.Tinf
				);
			loads.push_back(thisload);
			
			Log << "done!" << gaplog::endl;
		}
		if(bcs[i].bctype == neumann)
		{

			Log << "      - adding neumann th_boundary on set " 
					 << bcs[i].setname << " ... ";

			th_load* thisload = 
				new heat_flux
				(
					*this,
					bcs[i].setname.c_str(),
					bcs[i].neumann.q
				);
			loads.push_back(thisload);
			
			Log << "done! " << gaplog::endl;
		
		}
	}
	
	Log << "\n  * applying the loads ... ";

	// apply the loads
	for(unsigned int l=0; l<loads.size(); l++)
	{
		
		// apply the load
		loads[l]->apply();

		std::map<int, double>::iterator it;
		for(it = loads[l]->lgdof_val.begin(); it!=loads[l]->lgdof_val.end(); it++)
		{
			int l_id = it->first;
			if(gdof[l_id] > -1)
				b[gdof[l_id]] += it->second;
		}

	}

	Log << "done!" << gaplog::endl;


}
// ------------------------------------------------------------------------- //
void htr_analysis::get_K()
{
	// set to zero the content of the stiffness matrix!
	gmm::clear(K);

	// number of non zero value in stiffness matrix
	int nnz = 0;

	// ----------- this is for the operation progress information ------------ //
	
	double single_step = 1.0; // percentage
	double next_step = single_step; // percentage

	std::streamsize default_val = cout.precision();

	cout.precision(1);
	Log << "  * Assembling the stiffness matrix " << std::fixed << 0.00 << " %";

	// ---------------------- loop to all the elements  ---------------------- //

	for(int e = 0; e < msh.elements.size(); e++)
	{
		
		// ---------------- gives information on operation progress ------------ //
		double current = (double) e/msh.elements.size();
		if(current > next_step/100) 
		{
			next_step += single_step;
			Log << "\r  * Assembling the stiffness matrix " 
					 << std::fixed << 100*current << " %";
		}

		// calculate the stiffness matrix
		matrix ke = elements[e].calc_K();

		// add convective stiffness matrix to thermal th_boundary elements 
		if(msh.elements[e].belm)
		{
			for(unsigned int f=0; f<msh.elements[e].ext_faces.size(); f++)
			{
				int ext_face_id = msh.elements[e].ext_faces[f];
				const th_boundary bc = get_face_bnd(ext_face_id);
				// if convection th_boundary, then add the conductive stiffness matrix
				if(bc.bctype == mixed)
				{
					int elm_fid = msh.faces[ext_face_id].elm_fid;	
					double h = bc.mixed.h;
					// calculate the int([N]'[N]) on the face
					matrix Kcv = h*(elements[e].calc_intNTN(elm_fid));
					// sum to the conductive stiffness matrix
					ke = (ke + Kcv);
				}
			}
		}

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
void htr_analysis::solve()
{
	solver.system.quiet = true;
	solver.system.tol = 1e-6;
	solver.system.itermax = 2000;
	solver.system.P = ILDLT;

	Log << "\nSolving the linear system ... ";

	solver.initialize_x();
	solver.update();

	// clear the stiffness matrix now that it has been copied
	// into the solver structure
	K.clear_mat();
	K.resize(0,0);
	
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
}
// ------------------------------------------------------------------------- //
void htr_analysis::clear()
{
	bcs.clear();
	while (!elements.empty())
		elements.pop_back();
	b.clear();
	K.clear_mat();
	solver.clear();
	gdof.clear();
}
// ------------------------------------------------------------------------- //

