# include "./th_load.h"

# include "./th_load.h"
# include "./htr_analysis.h"
# include <algorithm>
# include "../FSTI_Block_dll/log.h"


using namespace std;
extern class gaplog Log;


// ------------------------------------------------------------------------- //
void heat_flux::apply()
{

	// get the list of faces in the face set
	const geom_set::const_iterator it = fea.msh.face_sets.find(set);
	if(it == fea.msh.face_sets.end())
	{
		Log << "\nheat_flux::apply() could not find face set "
				 << set << " in the mesh!" << gaplog::endl;
		exit(1);
	}

	const vector<int>& fs = it->second;
	
	lgdof_val.clear();
			
	// adjust the heat flux vector to the proper size, if a constant heat is used
	if(q.size() == 1)
		q.resize(fs.size(), q[0]);
	
	for(unsigned int f=0; f<fs.size(); f++)
	{

		// get the id of the associated element
		int eid = fea.msh.faces[fs[f]].elm;
		// get the face id for this element
		int elm_fid = fea.msh.faces[fs[f]].elm_fid;
			
		// ndof x 1 size matrix with the load
		matrix Q = q[f]*fea.elements[eid].calc_intN(elm_fid);

		for(unsigned int i=0; i<Q.m(); i++) 
		{
			lgdof_val[fea.elements[eid].gnum[i]] += Q[i][0];
		}
	}
	
}
// ------------------------------------------------------------------------- //
void convection::apply()
{

	// get the list of faces in the face set
	const geom_set::const_iterator it = fea.msh.face_sets.find(set);
	if(it == fea.msh.face_sets.end())
	{
		Log << "\nheat_flux::apply() could not find face set "
				 << set << " in the mesh!" << gaplog::endl;
		exit(1);
	}

	const vector<int>& fs = it->second;
		
	for(unsigned int f=0; f<fs.size(); f++)
	{

		// get the id of the associated element
		int eid = fea.msh.faces[fs[f]].elm;
		// get the face id for this element
		int elm_fid = fea.msh.faces[fs[f]].elm_fid;
		
		// ndof x 1 size matrix with the load
		matrix Q = h*Tinf*fea.elements[eid].calc_intN(elm_fid);

		for(unsigned int i=0; i<Q.m(); i++)
		{
			lgdof_val[fea.elements[eid].gnum[i]] += Q[i][0];
		}
	}
}
/*
// ------------------------------------------------------------------------- //
void expansion::apply()
{
	
	const mesh& msh = fea.get_mesh();
	const vector<fem_element*> femelm = fea.get_elements();

	lgdof_val.clear();

	double Tref = 20.0;
	
	for(int e=0; e<msh.ne; e++)
	{
		matrix B = femelm[e]->B;							// strain-displ matrix
		matrix D = femelm[e]->D;							// constitution matrix
		double V = femelm[e]->V;							// elm volume
		double alpha = femelm[e]->mat.alpha;	// corff of linear expansion
		
		matrix epsilon(6,1);
		epsilon[0][0] = alpha*(T[e] - Tref);	// epsilon x
		epsilon[1][0] = alpha*(T[e] - Tref);	// epsilon y
		epsilon[2][0] = alpha*(T[e] - Tref);	// epsilon z
		epsilon[3][0] = 0;										// gamma xy
		epsilon[4][0] = 0;										// gamma xz
		epsilon[5][0] = 0;										// gamma yz

		matrix Q = V*(B.T()*D)*epsilon;

		// summ the contribution of this element on the loads vector
		for(unsigned int i=0, id = 0; i<msh.elements[e].nodes.size(); i++)
		{
			for(unsigned int j=0; j<3; j++, id++)
			{
				lgdof_val[femelm[e]->gnum[id]] += Q[3*i + j][0];
			}
		}

	}
}
// ------------------------------------------------------------------------- //
*/
