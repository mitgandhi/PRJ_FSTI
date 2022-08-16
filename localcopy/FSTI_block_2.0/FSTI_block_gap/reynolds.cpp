# include "./reynolds.h"
# include "./block_limits.h"
# include "../FSTI_Block_dll/log.h"
# include "./block_limits.h"
# include <iomanip>
# include <omp.h>
# include <ANN/ANN.h>

# define pi 3.1415926535898

extern gaplog Log;

using namespace std;

// ------------------------------------------------------------------------- //
Reynolds::Reynolds(cblock_gap_main& Gap) :
	gap(Gap),
	msh(gap.mesh.reynolds),
	p(gap.fields.p),
	p_EHD(gap.fields.p),
	p_sol(gap.fields.p_sol),
	mu(gap.fields.mu2D),
	rho(gap.fields.rho2D),
	dhdt(gap.squeeze.dhdt),
	dhdt_hs(gap.squeeze.dhdt_hs),
	dhdt_hd(gap.squeeze.dhdt_hd),
	hcb(gap.film.hcb),
	hcb_rigid(gap.film.hcb_rigid),
	hvp(gap.film.hvp),
	dhcb_EHD(gap.film.dhcb_EHD),
	dhvp_EHD(gap.film.dhvp_EHD),
	contact(gap.film.contact),
	dhcb_hs(&gap.mesh.reynolds,0),
	dhcb0_hs(&gap.mesh.reynolds,0),
	dhcb0_hs_prev(&gap.mesh.reynolds,0),
	dhvp_hs(&gap.mesh.reynolds,0),
	dhvp_hs_prev(&gap.mesh.reynolds,0),
	dhdt_cb_hs(&gap.mesh.reynolds,0),
	dhdt_vp_hs(&gap.mesh.reynolds,0),
	dhcb_hd(&gap.mesh.reynolds,0),
	dhcb0_hd(&gap.mesh.reynolds,0),
	dhcb0_hd_prev(&gap.mesh.reynolds,0),
	dhcb_hd_prev(&gap.mesh.reynolds,0),
	dhvp_hd(&gap.mesh.reynolds,0),
	dhvp_hd_prev(&gap.mesh.reynolds,0),
	dhdt_cb_hd(&gap.mesh.reynolds,0),
	dhdt_vp_hd(&gap.mesh.reynolds,0),
	cdef0_cb(&gap.mesh.reynolds,0),
	cdef_cb(&gap.mesh.reynolds,0),
	cdef_vp(&gap.mesh.reynolds,0),
	pDC(gap.pistons.pDC),
	pHP(gap.pistons.pHP),
	pLP(gap.pistons.pLP),
	solver(&A, &b),
	CBf(gap.CBf),
	VPf(gap.VPf),
	IMset(0)
{

	// ------------ copy from caspar input --------------- //
	
	hmin = gap.in.data.options_block.numeric.hmin;
	omega = gap.in.data.operating_conditions.speed;
	EHD_CB = gap.in.data.options_block.general.EHD_CB;
	EHD_VP = gap.in.data.options_block.general.EHD_VP;
			
	// --------------------------------------------------- //

	hydrodynamic = true;
	dt = gap.dt;
	size = gap.mesh.reynolds.mn_f;
	A.resize(size,size);
	b.resize(size);
	update_dimension();
	solver.update();
	current_relax = 1.0;

	// average value of the IM associated to each fluid cell
	IM_cb_avg = 0;
	IM_vp_avg = 0;

	// initialize the influence matrices structure
	if(EHD_CB || EHD_VP)
	{
		IMset = new influ_set(gap.in, gap.fields.p.mesh);
	}

}
// ------------------------------------------------------------------------- //
Reynolds::~Reynolds()
{

}
// ------------------------------------------------------------------------- //
void Reynolds::init_def_map(body_type B)
{
	

	const interpl_grid* S;
	
	if(B == CB)
	{
		// check if the matrices are loaded
		if(IMset->CB.gap == 0)
		{
			Log << "\nReynolds::init_def_map: The matrices for the CB structure are not loaded!";
			Log << "\nTerminating.\n\n";
			exit(1);
		}

		S = &IMset->CBs;
	}
	if(B == VP)
	{
		// check if the matrices are loaded
		if(IMset->VP.gap == 0)
		{
			Log << "\nReynolds::init_def_map: The matrices for the VP structure are not loaded!";
			Log << "\nTerminating.\n\n";
			exit(1);
		}

		S = &IMset->VPs;
	}

	int dim = 2;								// dimension
		
	// fluid to structure id
	vector<int> FC2SC(gap.mesh.reynolds.mn);					
	
	ANNpointArray dataPts = annAllocPts(S->elements.size(), dim);

	// fill the data points with the structure centroids
	for(unsigned int i=0; i<S->elements.size(); i++)
	{
		point C(0,0,0);
		for(int j=0; j<3; j++)
			C = C + S->nodes[S->elements[i][j]];
		dataPts[i][0] = C.x()/3.0; 
		dataPts[i][1] = C.y()/3.0;
		dataPts[i][2] = C.z()/3.0;
	}

	// build search structure
	ANNkd_tree* kdTree = new ANNkd_tree
	( 
		dataPts,	// the data points
		S->elements.size(),		// number of points
		dim
	);

	ANNpoint queryPt = annAllocPt(dim);
	ANNidxArray nnIdx = new ANNidx[1];
	ANNdistArray dists = new ANNdist[1];

	// for each fluid element, define the closest structure element
	for(unsigned int i=0; i<gap.mesh.reynolds.mn; i++)
	{

		queryPt[0] = gap.mesh.reynolds.elements[i].x;
		queryPt[1] = gap.mesh.reynolds.elements[i].y;

		kdTree->annkSearch(queryPt, 1, nnIdx, dists);
		FC2SC[i] = nnIdx[0];
	}

	// now, looking into the proper IM, define the displacement associated
	// with the node closest to the fluid element centroid

	delete kdTree;
	annDeallocPts(dataPts);

	dataPts = annAllocPts(S->nodes.size(), dim);

	// fill the data points with the structure nodes
	for(unsigned int i=0; i<S->nodes.size(); i++)
	{
		dataPts[i][0] = S->nodes[i].x(); 
		dataPts[i][1] = S->nodes[i].y();
		dataPts[i][2] = S->nodes[i].z();
	}

	// build the new kdtree
	kdTree = new ANNkd_tree
	( 
		dataPts,	// the data points
		S->nodes.size(),		// number of points
		dim
	);

	// define the pointer to the proper vector 
	vector<double>* cdef;
	const vector<int>* fluid_id;
	if(B == CB)
	{
		cdef = &(cdef0_cb.in);
		fluid_id = &gap.mesh.reynolds.cb_elms_0;
	}
	else
	{
		cdef = &(cdef_vp.in);
		fluid_id = &gap.mesh.reynolds.vp_elms;
	}

	// for each fluid element, define the closest structure element
	for(unsigned int i=0; i<gap.mesh.reynolds.mn; i++)
	{

		if(fluid_id->at(i) == FLUID)
		{

			queryPt[0] = gap.mesh.reynolds.elements[i].x;
			queryPt[1] = gap.mesh.reynolds.elements[i].y;

			// at this point I know what is the node closest to the 
			// fluid element centroid
			kdTree->annkSearch(queryPt, 1, nnIdx, dists);
			int closest_node = nnIdx[0];

			// scale the deformation using the ration between fluid and
			// structure element size

			double AF = gap.mesh.reynolds.elements[i].A;
			point PA = S->nodes[S->elements[FC2SC[i]][0]];
			point PB = S->nodes[S->elements[FC2SC[i]][1]];
			point PC = S->nodes[S->elements[FC2SC[i]][2]];
		
			double AS = 0.5*fabs( 
				(PA.x() - PC.x())*(PB.y() - PA.y()) 
				-	(PA.x() - PB.x())*(PC.y() - PA.y()) 
			);

			// the deformation is also divided by 100e5 since the IM was 
			// calculated using a pressure of 100 bar

			if(B == CB)
				cdef->at(i) = (AF/AS)*IMset->CB.gap[FC2SC[i]][closest_node]/100e5;
			
			else
				cdef->at(i) = (AF/AS)*IMset->VP.gap[FC2SC[i]][closest_node]/100e5;
			
		}
		
	}

	// update the average value
	if(B == CB)
	{
		IM_cb_avg = 0;
		for(unsigned int i=0; i < cdef->size(); i++)
			IM_cb_avg += cdef->at(i);
		IM_cb_avg /= cdef->size();
	}
	else
	{
		IM_vp_avg = 0;
		for(unsigned int i=0; i < cdef->size(); i++)
			IM_vp_avg += cdef->at(i);
		IM_vp_avg /= cdef->size();
	}

	delete kdTree;
	delete [] nnIdx;
	delete [] dists;

	annDeallocPts(dataPts);

	annClose();

}
// ------------------------------------------------------------------------- //
bool Reynolds::EHD() const
{
	return (EHD_CB || EHD_VP);
}
// ------------------------------------------------------------------------- //
const interpl_grid& Reynolds::get_CBs() const
{
	return IMset->CBs;
}
// ------------------------------------------------------------------------- //
const interpl_grid& Reynolds::get_VPs() const
{
	return IMset->VPs;
}
// ------------------------------------------------------------------------- //
void Reynolds::get_D(int id)
{

	int i = id/gap.mesh.reynolds.n;		// row id
	int j = id%gap.mesh.reynolds.n;		// column id

	// gap element
	gap_elm e = msh.elements[id];
	double dr = e.dr;
	double dtheta = e.dtheta;
	double r = e.r;
	
	// film thickness at the faces
	double hcb_w, hcb_e, hcb_s, hcb_n;
	double hvp_w, hvp_e, hvp_s, hvp_n;
	
	// viscosity at the faces
	double mu_w, mu_e, mu_s, mu_n;
	
	// density at the faces
	double rho_w, rho_e, rho_s, rho_n;

	// ---------------------------- West face -------------------------- //

	//if(msh.elements[e.w].ty == FLUID && saturation.in[e.w] == 0)
	if(msh.elements[e.w].ty == FLUID)
	{
		hcb_w = 0.50*(hcb.in[id] + hcb.in[e.w]);
		hvp_w = 0.50*(hvp.in[id] + hvp.in[e.w]);
	}
	else
	{
		hcb_w = hcb.in[id];
		hvp_w = hvp.in[id];
	}

	mu_w = 0.50*(mu.in[id] + mu.in[e.w]);
	rho_w = 0.50*(rho.in[id] + rho.in[e.w]);

	// ---------------------------- East face -------------------------- //

	//if(msh.elements[e.e].ty == FLUID && saturation.in[e.e] == 0)
	if(msh.elements[e.e].ty == FLUID)
	{
		hcb_e = 0.50*(hcb.in[id] + hcb.in[e.e]);
		hvp_e = 0.50*(hvp.in[id] + hvp.in[e.e]);
	}
	else
	{
		hcb_e = hcb.in[id];
		hvp_e = hvp.in[id];
	}

	mu_e = 0.50*(mu.in[id] + mu.in[e.e]);
	rho_e = 0.50*(rho.in[id] + rho.in[e.e]);
	
	// --------------------------- South face -------------------------- //
	
	if(e.s > -1) 
	{
		//if(msh.elements[e.s].ty == FLUID && saturation.in[e.s] == 0)
		if(msh.elements[e.s].ty == FLUID)
		{
			hcb_s = 0.50*(hcb.in[id] + hcb.in[e.s]);
			hvp_s = 0.50*(hvp.in[id] + hvp.in[e.s]);
		}
		else
		{
			hcb_s = hcb.in[id];
			hvp_s = hvp.in[id];
		}
		
		mu_s = 0.50*(mu.in[id] + mu.in[e.s]);
		rho_s = 0.50*(rho.in[id] + rho.in[e.s]);
	}
	else 
	{
		hcb_s = hcb.bound.inner[j];
		hvp_s = hvp.bound.inner[j];
		mu_s = mu.bound.inner[j];
		rho_s = rho.bound.inner[j];
	}
	
	// --------------------------- North face -------------------------- //
	
	if(e.n > -1) 
	{
		//if(msh.elements[e.n].ty == FLUID && saturation.in[e.n] == 0)
		if(msh.elements[e.n].ty == FLUID)
		{
			hcb_n = 0.50*(hcb.in[id] + hcb.in[e.n]);
			hvp_n = 0.50*(hvp.in[id] + hvp.in[e.n]);
		}
		else
		{
			hcb_n = hcb.in[id];
			hvp_n = hvp.in[id];
		}
		mu_n = 0.50*(mu.in[id] + mu.in[e.n]);
		rho_n = 0.50*(rho.in[id] + rho.in[e.n]);
	}
	else 
	{
		hcb_n = hcb.bound.outer[j];
		hvp_n = hvp.bound.outer[j];
		mu_n = mu.bound.outer[j];
		rho_n = rho.bound.outer[j];
	}

	// --------------------------------------------------------------------- //

	// double check on h_cb and h_vp values
	hvp_n = (hcb_n > hvp_n + hmin) ? hvp_n : hcb_n - hmin;
	hvp_s = (hcb_s > hvp_s + hmin) ? hvp_s : hcb_s - hmin;
	hvp_w = (hcb_w > hvp_w + hmin) ? hvp_w : hcb_w - hmin;
	hvp_e = (hcb_e > hvp_e + hmin) ? hvp_e : hcb_e - hmin;

	// get film thickness
	double hs = (hcb_s - hvp_s);
	double hn = (hcb_n - hvp_n);
	double hw = (hcb_w - hvp_w);
	double he = (hcb_e - hvp_e);

	// ---------------------------------- Ds --------------------------------- //

	if(e.s > -1) 
	{
		// fluid
		//if(msh.elements[e.s].ty == FLUID && saturation.in[e.s] == 0) 
		if(msh.elements[e.s].ty == FLUID)
			Ds = (rho_s*pow(hs, 3.0)*(r - 0.5*dr)*dtheta)/(12.0*mu_s*dr);
		// opening
		else 
			Ds = (rho_s*pow(hs, 3.0)*(r - 0.5*dr)*dtheta)/(12.0*mu_s*0.5*dr);
	}
	// inner boundary
	else 
	{
		Ds = (rho_s*pow(hs, 3.0)*(r - 0.5*dr)*dtheta)/(12.0*mu_s*0.5*dr); 
	}

	// ---------------------------------- Dn --------------------------------- //

	if(e.n > -1) 
	{
		// fluid
		//if(msh.elements[e.n].ty == FLUID && saturation.in[e.n] == 0) 
		if(msh.elements[e.n].ty == FLUID)
			Dn = (rho_n*pow(hn, 3.0)*(r + 0.5*dr)*dtheta)/(12.0*mu_n*dr);
		// opening
		else 
			Dn = (rho_n*pow(hn, 3.0)*(r + 0.5*dr)*dtheta)/(12.0*mu_n*0.5*dr);
	}
	// inner boundary
	else 
	{
		Dn = (rho_n*pow(hn, 3.0)*(r + 0.5*dr)*dtheta)/(12.0*mu_n*0.5*dr); 
	}

	// ---------------------------------- Dw --------------------------------- //

	//if(msh.elements[e.w].ty == FLUID && saturation.in[e.w] == 0) // fluid
	if(msh.elements[e.w].ty == FLUID)
		Dw = (rho_w*pow(hw, 3.0)*dr)/(12.0*mu_w*r*dtheta);
	// opening
	else 
	{
		Dw = (rho_w*pow(hw, 3.0)*dr)/(12.0*mu_w*r*0.5*dtheta);
	}

	// ---------------------------------- De --------------------------------- //

	//if(msh.elements[e.e].ty == FLUID && saturation.in[e.e] == 0) // fluid
	if(msh.elements[e.e].ty == FLUID)
		De = (rho_e*pow(he, 3.0)*dr)/(12.0*mu_e*r*dtheta);
	// opening
	else 
	{
		De = (rho_e*pow(he, 3.0)*dr)/(12.0*mu_e*r*0.5*dtheta);
	}
}
// ------------------------------------------------------------------------- //
void Reynolds::get_S(int id)
{
	
	int i = id/gap.mesh.reynolds.n;		// row id
	int j = id%gap.mesh.reynolds.n;		// column id

	// gap element
	gap_elm e = msh.elements[id];
	double dr = e.dr;
	double dtheta = e.dtheta;
	double r = e.r;

	// ----------------------------- grad_h ---------------------------- //
				
	double delta_theta_h = 0;

	// -------- east -------- //

	double hE = 0;
	double rhoE = rho.in[e.e];
	if(msh.elements[e.e].ty == FLUID)
	{
		hE = hcb.in[e.e] - hvp.in[e.e];
		delta_theta_h += dtheta;
	}
	else
	{
		hE = hcb.in[id] - hvp.in[id];
	}
	
	// ensure a positive film thickness
	// hE = (hE > hmin) ? hE : hmin;
	
	// -------- west -------- //

	double hW = 0;
	double rhoW = rho.in[e.w];
	if(msh.elements[e.w].ty == FLUID)
	{
		hW = hcb.in[e.w] - hvp.in[e.w];
		delta_theta_h += dtheta;
	}
	else
	{
		hW = hcb.in[id] - hvp.in[id];
	}
	// ensure a positive film thickness
	// hW = (hW > hmin) ? hW : hmin;

	// gradient of film thickness along the circumferential direction
	double grad_h = (delta_theta_h > 0) ? (rhoE*hE - rhoW*hW)/delta_theta_h : 0;

	// ---------------------------- grad_ht ---------------------------- //

	double delta_theta_hcb = 0;

	// -------- east -------- //

	double hcbE = 0;
	if(msh.elements[e.e].ty == FLUID)
	{
		hcbE = hcb.in[e.e] - dhcb_hd.in[e.e] - hcb_rigid.in[e.e];
		delta_theta_hcb += dtheta;
	}
	else
	{
		hcbE = hcb.in[id] - dhcb_hd.in[id] - hcb_rigid.in[id];
	}
	// ensure a positive film thickness
	//hcbE = (hcbE > hmin) ? hcbE : hmin;

	// -------- west -------- //

	double hcbW = 0;
	if(msh.elements[e.w].ty == FLUID)
	{
		hcbW = hcb.in[e.w] - dhcb_hd.in[e.w] - hcb_rigid.in[e.w];
		delta_theta_hcb += dtheta;
	}
	else
	{
		hcbW = hcb.in[id] - dhcb_hd.in[id] - hcb_rigid.in[id];
	}
	// ensure a positive film thickness
	//hcbW = (hcbW > hmin) ? hcbW : hmin;

	// gradient of h_cb along the circumferential direction
	double grad_ht = 
		(delta_theta_hcb > 0) ? (rho.in[id]*hcbE - rho.in[id]*hcbW)/delta_theta_hcb : 0;
	// ----------------------------------------------------------------- //

	// S = rho*(Vt * grad(ht) -(Vt/2) * grad(h) - dhdt)
	S = omega*grad_ht	- 0.5*omega*grad_h - rho.in[id]*dhdt.in[id];
	
	// add the EHD squeeze if specified
	if(gap.use_sqz_hs)
		S -= rho.in[id]*dhdt_hs.in[id];
	if(gap.use_sqz_hd)
		S -= rho.in[id]*dhdt_hd.in[id];
	
}
// ------------------------------------------------------------------------- //
void Reynolds::discretize(double alpha, bool EHD_sqz)
{

	// clear A and b
	gmm::clear(A);
	gmm::clear(b);

	int me = msh.m;
	int ne = msh.n;

	// rotate the elementary deformation firld for the CB if necessary
	if(EHD_sqz)
	{
		cdef_cb = cdef0_cb.cshift(gap.mesh.reynolds.rotation_steps);
	}
		
	// matrix coefficients
	double ap, aw, ae, as, an;

	// loop to all gap elements 
	for(int i = 0, f = 0, id = 0; i < me; i++) 
	{
		for(int j = 0; j < ne; j++, id++) 
		{
			
			// calculate coefficients only for the fluid elements
			//if(msh.elements[id].ty == FLUID && saturation.in[id] == 0) 
			if(msh.elements[id].ty == FLUID) 
			{

				gap_elm e = msh.elements[id];
				double dr = e.dr;
				double dtheta = e.dtheta;
				double r = e.r;
								
				b[f] = 0;
				ap = 0;
				
				// get the diffusive coeffs
				get_D(id);
				
				// ------------------------ West neighbor -------------------------- //

				// fluid
				//if(msh.elements[e.w].ty == FLUID && saturation.in[e.w] == 0) 
				if(msh.elements[e.w].ty == FLUID) 
					aw = Dw;
				// opening
				else 
				{
					aw = 0;
					double ab = Dw;
					ap += ab;
					//b[f] += ab*p.in[e.w];
					b[f] += ab*p_sol.in[e.w];
				}
				
				// ------------------------ East neighbor -------------------------- //
				
				// fluid
				//if(msh.elements[e.e].ty == FLUID && saturation.in[e.e] == 0) 
				if(msh.elements[e.e].ty == FLUID) 
					ae = De;
				// opening
				else	
				{
					ae = 0;
					double ab = De;
					ap += ab;
					//b[f] += ab*p.in[e.e];
					b[f] += ab*p_sol.in[e.e];
				}
				
				// ----------------------- South neighbor -------------------------- //

				// internal mesh
				if(e.s > -1) 
				{
					// fluid
					//if(msh.elements[e.s].ty == FLUID && saturation.in[e.s] == 0) 
					if(msh.elements[e.s].ty == FLUID) 
						as = Ds;
					// opening
					else 
					{
						as = 0;
						double ab = Ds;
						ap += ab;
						//b[f] += ab*p.in[e.s];
						b[f] += ab*p_sol.in[e.s];
					}
				}
				// inner boundary
				else 
				{
					as = 0;
					double ab = Ds; 
					ap += ab;
					//b[f] += ab*p.bound.inner[j];
					b[f] += ab*p_sol.bound.inner[j];
				}
				
				// ----------------------- North neighbor -------------------------- //

				// internal mesh
				if(e.n > -1) 
				{
					// fluid
					//if(msh.elements[e.n].ty == FLUID && saturation.in[e.n] == 0) 
					if(msh.elements[e.n].ty == FLUID) 
						an = Dn;
					// opening 
					else 
					{
						an = 0;
						double ab = Dn;
						ap += ab;
						//b[f] += ab*p.in[e.n];
						b[f] += ab*p_sol.in[e.n];
					}
				}
				// outer boundary
				else 
				{
					an = 0;
					double ab = Dn;
					ap += ab;
					//b[f] += ab*p.bound.outer[j];
					b[f] += ab*p_sol.bound.outer[j];
				}

				// --------------------------- source ------------------------------ //
				
				if(hydrodynamic) // normal solution
					get_S(id);
				else						 // hydrostatic solution
					S = 0;
	
				double Sc = 0, Sp = 0;

				if(alpha < 1.0)
				{
					//Sc = S + 1e-6*fabs(S)*p.in[id];
					Sc = S + 1e-8*fabs(S)*p_sol.in[id];
					Sp = -1e-8*fabs(S);
				}
				else
				{
					Sc = S;
					Sp = 0;
				}

				b[f] += Sc*e.A;
				ap -= Sp*e.A;

				// ---------- add the dh_sqz/dt squeeze term ---------- //
				
				if(EHD_sqz)
				{
					// cell P
					double dhsqz_P = cdef_cb.in[id] - cdef_vp.in[id];
					b[f] += rho.in[id]*p_EHD.in[id]*(dhsqz_P/dt)*e.A;
					ap += rho.in[id]*(dhsqz_P/dt)*e.A;
				}

				// ---------------------------------------------------- //
				
				// ----------------- Fill the finite volume matrix ----------------- //

				ap += (aw + ae + as  + an);	

				// add the under-relaxation term in the source
				//b[f] += ((1.0 - alpha)/alpha)*ap*p.in[id];
				b[f] += ((1.0 - alpha)/alpha)*ap*p_sol.in[id];
			
				// main diagonal
				A(f,f) = ap/alpha;	// account for under-relaxation on ap

				// neighbors
				if(aw != 0)
					A(f, msh.elements[e.w].fluid_id) = -aw;
				if(ae != 0)
					A(f, msh.elements[e.e].fluid_id) = -ae;
				if(as != 0)
					A(f, msh.elements[e.s].fluid_id) = -as;
				if(an != 0)
					A(f, msh.elements[e.n].fluid_id) = -an;

				// increment the fluid id index
				f++;
			}
		}
	}

	// update the linear system definition
	solver.update();

	// clear the matrix
	gmm::clear(A); // the data is now stored in the M matrix of the solver
}
// ------------------------------------------------------------------------- //
void Reynolds::update_dimension()
{
	// get the domain dimension, accounting for possible saturation regions
	int new_size = 0;
	int jump = 0;
	
	for(int i = 0; i<msh.mn; i++)
	{
		if(msh.elements[i].ty == FLUID)
		{
			// no saturation region
			//if(saturation.in[i] == 0) 
			new_size++;
			//else
			//	jump++;
		}
		else
		{
			jump++;
		}
		// update global fluid gapelement ID
		gap.mesh.reynolds.elements[i].fluid_id = i - jump;
	}
	
	if(new_size != size)
	{
		size = new_size;
		A.resize(new_size, new_size);
		b.resize(new_size);
		
		// update the solver with the new matrix size
		solver.update();
	}
}
// ------------------------------------------------------------------------- //
void Reynolds::update_solution()
{
	// update the p field
	for(int i = 0, f = 0; i < msh.mn; i++) 
	{
		if(msh.elements[i].ty == FLUID)
		{
			p_sol.in[i] = solver.x[f++];
			p.in[i] = p_sol.in[i];

			// saturate pressure to saturation
			if(p_sol.in[i] < p_min_limit)
			{
				p.in[i] = p_min_limit;
			}
			// limit maximim pressure value
			if(p_sol.in[i] > p_max_limit)
				p.in[i] = p_max_limit;
		}
	}
}
// ------------------------------------------------------------------------- //
void Reynolds::load_matrices()
{
	IMset->read_matrices();
}
// ------------------------------------------------------------------------- //
void Reynolds::delete_matrices()
{
	IMset->delete_matrices();
}
// ------------------------------------------------------------------------- //
void Reynolds::get_cb_gap_def()
{

	// rest the deformation field
	dhcb0_hd = 0;
	
	// don't do anything is the EHD_CB option is off
	if(!EHD_CB)
		return;

	int nn = IMset->CBs.nn();
	int ne = IMset->CBs.ne();

	// interpolate pressure from fluid to structure
	CB_interpl.cellsTocells(&CBf , &IMset->CBs);
	// pressure value at each structure element
	const vector<double>& p = IMset->CBs.cells_data;

	double* gapdef = new double[nn];
	// initialize the gap deformation array
	memset(gapdef, 0, nn*sizeof(double));	

	// get the number of processors available
	int nproc = omp_get_num_procs();
	// allocate as many deformation array as the available processors are
	double** gapdef_pr = new double*[nproc];
	// initialize each array to zero
	for(int pr=0; pr<nproc; pr++)
	{
		gapdef_pr[pr] = new double[nn];
		for(int i=0; i<nn; i++)
			gapdef_pr[pr][i] = 0.0;
	}

	// loop to all the elements using open mp
	# pragma omp parallel for
	for(int i=0; i<ne; i++)
	{
		// get the thread number
		int th = omp_get_thread_num();
		// calculate the scale factor
		double scale = p[i]/100e5;
		// get the deformation
		for(int j=0; j<nn; j++)
		{
			gapdef_pr[th][j] += scale*IMset->CB.gap[i][j];
		}
	}

	// merge all the deformation arrays into one
	for(int pr=0; pr<nproc; pr++)
	{
		for(int i=0; i<nn; i++)
			gapdef[i] += gapdef_pr[pr][i];
	}

	// interpolate from structure to fluid and copy to the scalar field

	// copy the result to CBs
	for(int i=0; i<IMset->CBs.nn(); i++)
		IMset->CBs.points_data[i] = gapdef[i];

	CB_interpl.pointsTocells(&IMset->CBs, &CBf);
	for(int i=0; i<CBf.ne(); i++)
	{
		if(CBf.active[i])
			dhcb0_hd.in[i] = CBf.cells_data[i];
		else
			dhcb0_hd.in[i] = 0;
	}

	// settle the boundaries
	for(int j=0; j<msh.n; j++)
	{
		dhcb0_hd.bound.inner[j] = dhcb0_hd.in[j];
		dhcb0_hd.bound.outer[j] = dhcb0_hd.in[(msh.m-1)*msh.n + j];
	}

	// rotate and copy to the dhcb_hd field
	dhcb_hd = dhcb0_hd.cshift(msh.rotation_steps);

	// delete temperary vectors
	for(int pr=0; pr<nproc; pr++)
		delete [] gapdef_pr[pr];
	delete [] gapdef_pr;
	delete [] gapdef;

}
// ------------------------------------------------------------------------- //
void Reynolds::get_vp_gap_def()
{

	// reset the deformation field
	dhvp_hd = 0;
	

	// don't do anything is the EHD_VP option is off
	if(!EHD_VP)
		return;

	int nn = IMset->VPs.nn();
	int ne = IMset->VPs.ne();

	// interpolate pressure from fluid to structure
	VP_interpl.cellsTocells(&VPf , &IMset->VPs);
	// pressure value at each structure element
	const vector<double>& p = IMset->VPs.cells_data;

	double* gapdef = new double[nn];
	// initialize the gap deformation array
	memset(gapdef, 0, nn*sizeof(double));	
	
	// get the number of processors available
	int nproc = omp_get_num_procs();
	// allocate as many deformation array as the available processors are
	double** gapdef_pr = new double*[nproc];
	// initialize each array to zero
	for(int pr=0; pr<nproc; pr++)
	{
		gapdef_pr[pr] = new double[nn];
		for(int i=0; i<nn; i++)
			gapdef_pr[pr][i] = 0.0;
	}

	// loop to all the elements using open mp
	# pragma omp parallel for
	for(int i=0; i<ne; i++)
	{
		// get the thread number
		int th = omp_get_thread_num();
		// calculate the scale factor
		double scale = p[i]/100e5;
		// get the deformation
		for(int j=0; j<nn; j++)
		{
			gapdef_pr[th][j] += scale*IMset->VP.gap[i][j];
		}
	}

	// merge all the deformation arrays into one
	for(int pr=0; pr<nproc; pr++)
	{
		for(int i=0; i<nn; i++)
			gapdef[i] += gapdef_pr[pr][i];
	}

	// interpolate from structure to fluid and copy to the scalar field

	// copy the result to VPs
	for(int i=0; i<IMset->VPs.nn(); i++)
		IMset->VPs.points_data[i] = gapdef[i];

	VP_interpl.pointsTocells(&IMset->VPs, &VPf);
	for(int i=0; i<VPf.ne(); i++)
	{
		if(VPf.active[i])
			dhvp_hd.in[i] = VPf.cells_data[i];
		else
			dhvp_hd.in[i] = 0;
	}

	// settle the boundaries
	for(int j=0; j<msh.n; j++)
	{
		dhvp_hd.bound.inner[j] = dhvp_hd.in[j];
		dhvp_hd.bound.outer[j] = dhvp_hd.in[(msh.m-1)*msh.n + j];
	}

	// delete temperary vectors
	for(int pr=0; pr<nproc; pr++)
		delete [] gapdef_pr[pr];
	delete [] gapdef_pr;
	delete [] gapdef;

}
// ------------------------------------------------------------------------- //
void Reynolds::get_dh_hs_cb()
{

	// copy the current into the previous
	dhcb0_hs_prev = dhcb0_hs;
	// reset the current hydrostatic deformation
	dhcb0_hs = 0;
	dhcb_hs = 0;
	dhdt_cb_hs = 0;

	if(EHD_CB)
	{
		// number of pistons
		int z = gap.in.data.operating_conditions.npistons;

		// deformation for each node
		vector<double> def(IMset->CBs.nn(), 0.0);

		// add the scaled IM for each DC
		for(int i=0; i<z; i++)
		{
			double scale = pDC[i]/100e5;	// influence matrix scale factor
			for(int j=0; j<IMset->CBs.nn(); j++)
				def[j] += scale*IMset->CB.DC[i][j];
		}

		// interpolate the deformation form the structure to the fluid grid
		IMset->CBs.points_data = def;
		CB_interpl.pointsTocells(&IMset->CBs, &CBf);
	
		// copy to the deformation scalar_field
		for(unsigned int j=0; j<dhcb0_hs.in.size(); j++)
			dhcb0_hs.in[j] = CBf.cells_data[j];

		// update the EHD squeeze
		dhdt_cb_hs = (dhcb0_hs - dhcb0_hs_prev)/dt;
		dhdt_cb_hs = dhdt_cb_hs.cshift(msh.rotation_steps);

		// return the DC deformation properly scaled
		dhcb_hs = dhcb0_hs.cshift(msh.rotation_steps);
	}

}
// ------------------------------------------------------------------------- //
void Reynolds::get_dh_hs_vp()
{

	// copy the current into the previous
	dhvp_hs_prev = dhvp_hs;
	dhvp_hs = 0;
	dhdt_vp_hs = 0;

	if(EHD_VP)
	{
		// deformation for each node
		vector<double> def(IMset->VPs.nn(), 0.0);
	
		// get the deformation of the low pressure port 
		double scale = pLP/100e5;
		for(int j=0; j<IMset->VPs.nn(); j++)
				def[j] += scale*IMset->VP.LP[j];
	
		// get the deformation of the high pressure port 
		scale = pHP/100e5;
		for(int j=0; j<IMset->VPs.nn(); j++)
				def[j] += scale*IMset->VP.HP[j];

		// interpolate the deformation form the structure to the fluid grid
		IMset->VPs.points_data = def;
		VP_interpl.pointsTocells(&IMset->VPs, &VPf);

		// copy to the deformation scalar_field
		for(unsigned int j=0; j<dhvp_hs.in.size(); j++)
			dhvp_hs.in[j] = VPf.cells_data[j];

		// update the EHD squeeze
		dhdt_vp_hs = (dhvp_hs - dhvp_hs_prev)/dt;
	}

}
// ------------------------------------------------------------------------- //
void Reynolds::update_interpl(const scalar_field& _p)
{
	// update the cylinder block fluid interpolation grid
	if(EHD_CB)
	{
		scalar_field p_ref = _p;
		
		p_ref = p_ref.cshift(-msh.rotation_steps);
				
		for(int i=0; i<CBf.ne(); i++)
		{
			double _pi = p_ref.in[i];
			
			// saturate pressure values
			_pi = (_pi > p_min_limit) ? _pi : p_min_limit;
			_pi = (_pi < p_max_def) ? _pi : p_max_def;	// use pressure limit for deformation
			
			if(CBf.active[i])
				CBf.cells_data[i] = _pi;
			else
				CBf.cells_data[i] = 0;
		}
	}
	// update the valve plate fluid interpolation grid
	if(EHD_VP)
	{
		scalar_field p_ref = _p;

		for(int i=0; i<VPf.ne(); i++)
		{
			double _pi = p_ref.in[i];

			
			// saturate pressure values
			_pi = (_pi > p_min_limit) ? _pi : p_min_limit;
			_pi = (_pi < p_max_def) ? _pi : p_max_def; // use pressure limit for deformation
	
			if(VPf.active[i])
				VPf.cells_data[i] = _pi;
			else
				VPf.cells_data[i] = 0;
		}
	}
}
// ------------------------------------------------------------------------- //
double avg(const std::vector<double>& v)
{

	int n = static_cast<int>(v.size());

	// get the average
	double avg = 0;
	for(int i=0; i<n; i++)
		avg += v[i];
	avg /= n;

	return avg;

}
// ------------------------------------------------------------------------- //
void Reynolds::temp()
{
	// get the starting deformations for CB and VP
	get_dh_hs_cb();
	gap.film.dhcb_EHD = dhcb_hs;
	gap.film.dhvp_EHD = 0;
	dhdt_hs = (dhdt_cb_hs - dhdt_vp_hs);
}
// ------------------------------------------------------------------------- //
int Reynolds::solve_rigid(bool EHD_sqz)
{
	// initialize the solution with the current pressure field
	// this is very important for the boundary conditions!
	p_sol = p;

	// reset the saturation field
	//saturation = 0.0;
	update_dimension();

	// initialize the solution vector
	solver.initialize_x(); 

	// discretize the system
	discretize(1.0, EHD_sqz);
	
	// precondition the fv matrix
	solver.preconditioner();
	
	// solve
	int iters = solver.solve();
	
	// update the pressure field
	update_solution();
	
	// update the interpolation grids
	// update_interpl();

	// cout << "Iters: " << iters << "\t" << solver.residual() << endl;


	return iters;
}
// ------------------------------------------------------------------------- //
int Reynolds::solve_EHD()
{
	// -------- initialization -------- //

	// initialize the solution field (very important for the boundary conditions!)
	p_sol = p;

	// copy the previous hydrodynamic deformation fields
	dhcb0_hd_prev = dhcb0_hd;
	dhcb_hd_prev = dhcb_hd;
	dhvp_hd_prev = dhvp_hd;

	// update the liner system dimension
	update_dimension();

	// update the interpolation structure with the current pressure field
	update_solution();
	
	// it is crucial to have update_interpl before get_cb_gap_def and get_vp_gap_def
	update_interpl(p_sol);
	get_cb_gap_def();
	get_vp_gap_def();

	// get the starting deformations for CB and VP
	get_dh_hs_cb();
	get_dh_hs_vp();

	// update the deformation squeeze
	if(gap.t > 0)
		dhdt_hs = (dhdt_cb_hs - dhdt_vp_hs);
	
	dhdt_hd = 0;	// reset the hydrodynamic squeeze
	
	// update the overall deformation, just with the hydrodynamic deformation
	// in order to avoid initial penetration
	dhcb_EHD = dhcb_hd;
	dhvp_EHD = dhvp_hd;
	
	// get the film thickness and set to hmin if neces
	gap.set_film_thickness();

	// initialize solver 
	for(int i = 0, f = 0; i < msh.mn; i++) 
	{
		if(msh.elements[i].ty == FLUID)
			solver.x[f++] = p_sol.in[i];
	}

	cout.precision(3);
	
	// convergence stuff		
	
	double relax = current_relax + 0.1;
	relax = (relax <= 0.9) ? relax : 0.9;
	
	int iters = 0;								// overall iterations
	std::vector<double> res(0);		// residual vector
	int repeat = 5;								// convergence check
	double tol = 1.0;							// [ % ]
	double min_relax = 0.1;				// minimum relaxation factor
	double res_i = 0;							// current residual
	double max_i = 0;							// current maximum pressure
	//double min_max = 1e20;				// minimum value of maximum pressure reached
	double min_res = 1e20;				// minimum residual reached
	int max_iters = 250;
	double conv_tol = 1e-3;				// overall convergence tolerance
	
	int add_hd = 0;

	double mean = 0, mean_old = 0;

	// hydrostatic deformation, progressive magnitude
	scalar_field dhcb_hs_prog(&gap.mesh.reynolds,0);
	scalar_field dhvp_hs_prog(&gap.mesh.reynolds,0);


	// ---------------------------- main loop ------------------------------ //

	Log << "FSI convergence loop ... \n\n";

	Log << "|  it   "  << "|" << "  rlx  " << "|" << " pmax[bar] " << "|" 
			<< "  residual  " << "|"  << " cnt[um] " << "| hmin[um] |\n";
	Log << "+-------+-------+-----------+------------+---------+----------+\n";

	
	while(true)
	{
			
		// ------------- discretize and solve with the relaxation -------------- //

		// discretize with relaxation, no dh_EHD squeeze
		discretize(relax, false);
		solver.preconditioner();
		solver.solve();
		update_solution();

		// update oil properties
		gap.update_viscosity();
		gap.update_density();

		// update the film thickness with the new pressure field
		update_interpl(p_sol);
		get_cb_gap_def();
		get_vp_gap_def();

		

		// update the deformation
		dhcb_EHD = dhcb_hd + dhcb_hs;
		dhvp_EHD = dhvp_hd + dhvp_hs;

		// update film thickness
		gap.set_film_thickness();

		// update the hydrodynamic squeeze 
		// make sure that it is not the very first time step
		if(gap.use_sqz_hd && gap.t > 0)
		{
			scalar_field tmp_cb = (dhcb0_hd - dhcb0_hd_prev)/dt;
			tmp_cb = tmp_cb.cshift(gap.mesh.reynolds.rotation_steps);
			tmp_cb.limit(-1e-2, 1e-2); // saturate to avoid instabilities
			dhdt_cb_hd = dhdt_cb_hd + 0.1*(tmp_cb - dhdt_cb_hd);

			scalar_field tmp_vp = (dhvp_hd - dhvp_hd_prev)/dt;
			tmp_vp.limit(-1e-2, 1e-2); // saturate to avoid instabilities
			dhdt_vp_hd = dhdt_vp_hd + 0.1*(tmp_vp - dhdt_vp_hd);

			dhdt_hd = dhdt_cb_hd - dhdt_vp_hd;
		}

		// ---------------------- residual and convergence --------------------- //
		
		// get the residual
		discretize(1.0, false);
		res_i = solver.residual();

		ostringstream oss;
		oss << "./vtk/iter." << iters << ".vtk";
		gap.write_light_vtk(oss.str().c_str(),1e3);

		// check for convergence
		if(res_i < conv_tol && iters > 30)
		{
			
			current_relax = relax;

			// update the solution
			update_solution();

			cout.precision(8);
	
			return iters;
			
		}

		// get minimum residual
		res.push_back(res_i);
		if(res_i < min_res)
			min_res = res_i;
		
		max_i = p.max(); // get the maximum pressure
		
	
		// ------------------------ tune the relaxation ------------------------ //

		if(iters%repeat == 0 && iters > 0)
		{

			/*

			// get the mean residual
			mean_old = mean;
			mean = avg(res);
			
			// resize the residual
			res.resize(0);

			// reached the best achievable with this relaxation level
			if(mean > mean_old && mean !=0 && mean_old != 0 && iters > 0)
			{
				if(relax > min_relax)
					relax = 0.75*relax;
				
				mean = mean_old = 0;
			}

			// decide if is better to stop because there is no further improvement 
			// to the solution
			if
			(
				mean != 0 && mean_old != 0	// make sure the relax has not just changed
				&& fabs(100*(mean-mean_old)/mean_old) < 1
				&& res_i < 1e-3
			)
			{
	
				current_relax = relax;
		
				Log << "\nThe iteration process is not moving to any sensible improvement after "
						 << iters << " iterations\n" << gaplog::endl
						 << "\tres: " << std::scientific << res_i
						 << " pmax: " << max_i/1e5 
						 << gaplog::endl << gaplog::endl;

				// before quitting, set the hydrostatic field
				update_solution();

				cout.precision(8);
				
				return iters;
			}

			*/

			// tune the relaxation
			if(res.front() < res.back())
			{
				// update the relaxation
				relax = 0.9*relax;
				
				if(relax < min_relax)
				{
					Log << "\nThe iteration process is not moving to any sensible improvement after "
						 << iters << " iterations\n" << gaplog::endl
						 << "\tres: " << std::scientific << res_i
						 << " pmax: " << max_i/1e5 
						 << gaplog::endl << gaplog::endl;

					update_solution();

					cout.precision(8);
				
					return iters;
				}
			}

			// resize the residual array
			res.resize(0);

			// ---------------------- print convergence information -------------------------- //

			std::streamsize default_val = cout.precision();
			
			Log << "| ";
			cout << setw(5);
			Log << iters << " | ";
			cout.precision(3), cout << setw(5);
			Log << std::fixed << relax  << " | ";
			cout.precision(1), cout << setw(9);
			Log << std::fixed << max_i/1e5 << " | ";
			cout.precision(3);
			Log << std::scientific << res_i << " | ";
			cout.precision(3), cout << setw(7);
			Log << std::fixed << 1e6*contact.max() << " | ";
			cout << setw(8);
			Log << std::fixed << 1e6*gap.film.hmin << " | ";
			cout.precision(1);

			Log << "\n";
			
			cout.precision(default_val);
		
		}
			
		// update iterations number
		iters++;		

		// ------------------- check for convergence problems ------------------ //		
		
		if(iters > max_iters)
		{

			Log << "reynolds::solve(): WARNING! Reynolds did not converged!\n"
					 << gaplog::endl;
			
			// update the solution

			update_solution();
			
			current_relax = relax;
			cout.precision(8);
		
			return -1;
		}
		
	}


	
}
// ------------------------------------------------------------------------- //