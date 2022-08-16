# include "./Energy.h"
# include <fstream>

Energy::Energy(cblock_gap_main& Gap)
:
	gap(Gap),
	T(Gap.fields.T),
	V(Gap.fields.V),
	Vc(Gap.fields.Vc),
	hcb(gap.film.hcb),
	hvp(gap.film.hvp),
	mu(Gap.fields.mu),
	phid(Gap.fields.phid),
	size(Gap.fields.T.mesh->mnq_f),
	b(size, 0.0),
	A(size, size),
	solver(&A, &b)
{
	
	// ------------------------ copy from caspar_input ----------------------- //

	TLP = Gap.in.data.operating_conditions.T_LP;
	THP = Gap.in.data.operating_conditions.T_HP;
	Tcase = Gap.in.data.operating_conditions.T_Leak;

	// ----------------------------------------------------------------------- //

}
// ------------------------------------------------------------------------- //
Energy::~Energy()
{
}
// ------------------------------------------------------------------------- //
void Energy::update_dimension()
{
	int new_size = T.mesh->mnq_f;
	
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
std::vector<double> Energy::getd(int i, int j, int k)
{

	int id2D = gap.mesh.reynolds.n*i + j;
	int id3D = k*(gap.mesh.reynolds.mn) + id2D;

	gap_elm e3D = gap.mesh.energy.elements[id3D];
	gap_elm e2D = gap.mesh.energy.elements[id2D];

	double lambda = gap.lubricant->get_lambda(T.in[id3D]);
	double cp = gap.lubricant->get_cp(T.in[id3D]);
	double gamma = lambda;    //lambda/cp;

	// element dimension
	double dr = gap.mesh.reynolds.elements[id2D].dr;
	if(e2D.s == -1 || e2D.n == -1) { dr = 0.5*dr; }
	double dtheta = gap.mesh.reynolds.elements[id2D].dtheta;
	double dz = (hcb.in[id2D] - hvp.in[id2D])/gap.mesh.energy.q;
	if(e3D.b == -1 || e3D.t == -1) { dz = 0.5*dz; }
		
	double hcb_s = (e2D.s > -1) ? 0.5*(hcb.in[id2D] + hcb.in[e2D.s]) : hcb.bound.inner[j];
	double hcb_n = (e2D.n > -1) ? 0.5*(hcb.in[id2D] + hcb.in[e2D.n]) : hcb.bound.outer[j];
	double hcb_w = 0.5*(hcb.in[id2D] + hcb.in[e2D.w]);
	double hcb_e = 0.5*(hcb.in[id2D] + hcb.in[e2D.e]);

	double hvp_s = (e2D.s > -1) ? 0.5*(hvp.in[id2D] + hvp.in[e2D.s]) : hvp.bound.inner[j];
	double hvp_n = (e2D.n > -1) ? 0.5*(hvp.in[id2D] + hvp.in[e2D.n]) : hvp.bound.outer[j];
	double hvp_w = 0.5*(hvp.in[id2D] + hvp.in[e2D.w]);
	double hvp_e = 0.5*(hvp.in[id2D] + hvp.in[e2D.e]);

	double dz_s = (hcb_s - hvp_s)/gap.mesh.energy.q;
	double dz_n = (hcb_n - hvp_n)/gap.mesh.energy.q;
	double dz_w = (hcb_w - hvp_w)/gap.mesh.energy.q;
	double dz_e = (hcb_e - hvp_e)/gap.mesh.energy.q;
	
	std::vector<double> d(6);

	// d[0] = ds
	// d[1] = dn
	// d[2] = dw
	// d[3] = de
	// d[4] = db
	// d[5] = dt

	d[0] = gamma*((e2D.r - 0.5*dr)*dtheta*dz_s)/dr;
	d[1] = gamma*((e2D.r + 0.5*dr)*dtheta*dz_n)/dr;
	d[2] = gamma*(dr*dz_w)/(e2D.r*dtheta);
	d[3] = gamma*(dr*dz_e)/(e2D.r*dtheta);
	d[4] = gamma*(e3D.A)/dz;
	d[5] = gamma*(e3D.A)/dz;


	return d;

}
// ------------------------------------------------------------------------- //
void Energy::discretize()
{

	const double hmin = gap.in.data.options_block.numeric.hmin;
	
	update_dimension();
		
	// calc the mass flow in each element
	gap.calc_mass_flow();
	gap.calc_dhdt_total();
	
	// clear A and b
	gmm::clear(A);
	gmm::clear(b);

	// mesh dimension
	int m = T.mesh->m;	// elm radial dir
	int n = T.mesh->n;	// elm circ. dir
	int q = T.mesh->q;	// elm axial dir
	int mn = m*n;				// element in one layer
	
	// velocity @ cell centroids, radial direction
	const scalar_field V_r = V.r(); // Vc.r()
	// velocity @ cell centroids, circumferential direction
	const scalar_field V_theta = V.theta();	//	Vc.theta();
	
	// diffusion coefficients
	double d_w, d_e, d_s, d_n, d_b, d_t;
	
	// convection coefficients
	double f_w, f_e, f_s, f_n;
	//pressurere gradients
	double gpw, gpe, gps, gpn;
	//pressurere gradients
	double pw, pe, ps, pn;

	// mass flow at the faces
	double mdotw,mdote,mdots,mdotn;
	// algebraic equation coefficient
	double ap, aw, ae, as, an, ab, at;

	// loop to all gap elements 
	for(int k = 0, id3D = 0, f = 0; k < q; k++) 
	{
		for(int i = 0, id2D = 0; i < m; i++) 
		{
			for(int j = 0; j < n; j++, id2D++, id3D++) 
			{
				double cp = gap.lubricant->get_cp(T.in[id3D]);
				// calculate coefficients only for the fluid elements
				if(T.mesh->elements[id3D].ty == FLUID) 
				{
					// get element i,j,k and its properties
					gap_elm e = T.mesh->elements[id3D];
					double r = e.r;
					double theta = e.theta;
					double dr = e.dr;
					double dtheta = e.dtheta;
					double dz = e.dz;
					
					// -------------------- diffusion coefficients ------------------- //

					std::vector<double> dijk = getd(i, j, k);
					d_s = dijk[0], d_n = dijk[1], d_w = dijk[2], d_e = dijk[3];
					d_b = dijk[4], d_t = dijk[5];			
															
					// ------------------- convection coefficients ------------------- //
					
					f_s = gap.fields.fs.in[id3D]; 
					f_n = gap.fields.fn.in[id3D]; 
					f_w = gap.fields.fw.in[id3D]; 
					f_e = gap.fields.fe.in[id3D]; 

					// ------------------- mass flow ------------------- //
					
					mdots = gap.fields.mdots.in[id3D]; 
					mdotn = gap.fields.mdotn.in[id3D]; 
					mdotw = gap.fields.mdotw.in[id3D]; 
					mdote = gap.fields.mdote.in[id3D]; 
					
					// ------------------- pressure gradients ------------------- //
					
					gps = gap.fields.g_ps.in[id3D]; 
					gpn = gap.fields.g_pn.in[id3D]; 
					gpw = gap.fields.g_pw.in[id3D]; 
					gpe = gap.fields.g_pe.in[id3D]; 

					// ------------------- pressure at face ------------------- //
					
					ps = gap.fields.ps.in[id2D]; 
					pn = gap.fields.pn.in[id2D]; 
					pw = gap.fields.pw.in[id2D]; 
					pe = gap.fields.pe.in[id2D]; 


					// ---------------- algebraic equation coefficient --------------- //

					b[f] = 0;
					ap = 0;
					aw = ae = as = an = ab = at = 0;

					// ----------------------------- west ---------------------------- //

					// aw
					// interior
					if(T.mesh->elements[e.w].ty == FLUID)
					{
						aw = (f_w > 0) ? f_w : 0;			// UDS
						aw += d_w;
					}
					// boundary 
					else
					{
						aw = 0;	
						if(f_w > 0) // inflow
						{
							// diffusion
							double a_bound = d_w;			
							ap += a_bound;			
							b[f] += a_bound*T.in[e.w];
							// convection
							ap += cp *f_w;	// this will make up for the -fw in ap further down
							b[f] += cp *fabs(f_w)*T.in[e.w];
						}
						// else nothing to be done, ap -= fw = ap += fabs(fw) is further down
					}
									
					// ---------------------------- east ----------------------------- //

					// ae
					// interior
					if(T.mesh->elements[e.e].ty == FLUID)
					{
						ae = (-f_e > 0) ? -f_e : 0;			// UDS
						ae += d_e;
					}
					// boundary
					else
					{
						ae = 0;
						if(f_e < 0) // inflow
						{
							// diffusion 
							double a_bound = d_e;
							ap += a_bound;
							b[f] += a_bound*T.in[e.e];
							// convection
							ap -= cp *f_e; // this will make up for the +fe in ap further down
							b[f] += cp *fabs(f_e)*T.in[e.e];
						}
						// else nothing to be done, ap += fe = ap += fabs(fe) is further down
					}				
								
					// ---------------------------- south ---------------------------- //

					// as
					// case boundary
					if(e.s == -1)	
					{
						as = 0;
						double a_bound = d_s;
						ap += a_bound;			
						b[f] += a_bound*Tcase;
					}
					// interior
					else
					{
						if(T.mesh->elements[e.s].ty == FLUID)
						{
							as = (f_s > 0) ? f_s : 0;			// UDS
							as += d_s;
						}
						// port boundary
						else
						{
							as = 0;
							if(f_s > 0) // inflow
							{
								// diffusion
								double a_bound = d_s;
								ap += a_bound;
								b[f] += a_bound*T.in[e.s];
								// convection
								ap += cp *f_s;	// this will make up for the -fs in ap further down
								b[f] += cp *fabs(f_s)*T.in[e.s];
							}
							// else nothing to be done, ap -= fs = ap += fabs(fs) is further down
						}
					}

					// --------------------------- north ----------------------------- //

					// an
					// case boundary
					if(e.n == -1) 
					{
						an = 0;
						double a_bound = d_n;
						ap += a_bound;
						b[f] += a_bound*Tcase;
					}
					// interior
					else 
					{
						if(T.mesh->elements[e.n].ty == FLUID)
						{
							an = (-f_n > 0) ? -f_n : 0;			// UDS
							an += d_n;
						}
						// port boundary
						else
						{	
							an = 0;
							if(f_n < 0) // inflow
							{
								// diffusion
								double a_bound = d_n;
								ap += a_bound;
								b[f] += a_bound*T.in[e.n];
								// convection
								ap -= cp *f_n;	// this will make up for the +fn in ap further down
								b[f] += cp *fabs(f_n)*T.in[e.n];
							}
							// else nothing to be done, ap += fn = ap += fabs(fn) is further down
						}
					}
				
					// ---------------------------- bottom --------------------------- //
					
					// ab (only diffusion)
					// interior
					if(e.b != -1)
					{	
						ab = d_b;
					}
					// boundary, valve plate side
					else
					{
						ab = 0;
						double a_bound = d_b;
						b[f] += a_bound*T.bound.bottom[id2D];
						ap += a_bound;
					}
					
					// ----------------------------- top ----------------------------- //
					
					// top
					// interior
					if(e.t != -1) 
					{
						at = d_t;
					}
					// boundary, cylinder block side
					else
					{
						at = 0;
						double a_bound = d_t;
						b[f] += a_bound*T.bound.top[id2D];
						ap += a_bound;
					}
					
					
					// --------------------------- source ---------------------------- //

					

					

					// dissipation
					double phi = 
						(mu.in[id3D])*e.V*                       //(mu.in[id3D]/cp)*e.V*
						(
							pow(V_r.get_zgrad_ijk(i,j,k), 2.0) + 
							pow(V_theta.get_zgrad_ijk(i,j,k), 2.0) + 
							(4.0/3.0)*pow(V.in[id3D].r()/r, 2.0) +
							pow(V.in[id3D].theta()/r, 2.0)
						);
					
					b[f] += phi;	// add to the source
					
					
					double cT = 9.0447e-4;
					double h_0 = 3.6159e5;
					double specific_v = 1 / gap.fields.rho.in[id3D];

					
					
					b[f] -= cT * mdote * pe;
					b[f] -= -1 * cT * mdotw * pw;
					b[f] -= cT * mdotn * pn;
					b[f] -= -1 * cT * mdots * ps;
					
					
					b[f] -= mdote * gpe * specific_v / 2;
					b[f] -= -1 * mdotw * gpw * specific_v / 2;
					b[f] -= mdotn * gpn * specific_v / 2;
					b[f] -= -1 * mdots * gps * specific_v / 2;
					
					
					double dhdt_total = gap.squeeze.dhdt_total.in[id2D]/q;

					
					double dpdt = (gap.fields.p.in[id2D]-gap.fields.p_prev.in[id2D]);

					
					// first step
					b[f] += (gap.fields.rho.in[id3D] * r* dtheta * dr * dhdt_total * gap.dt + 
								gap.fields.rho.in[id3D] * r * dtheta * dr * dz) /  gap.fields.rho.in[id3D]	* 
									dpdt / gap.dt;
					
					
					// third step
					b[f] -= gap.fields.rho.in[id3D] * r * dtheta * dr * dz * gap.fields.p.in[id2D]*cT/gap.dt;

					b[f] += (gap.fields.rho.in[id3D] *r*dtheta * dr * dz + gap.fields.rho.in[id3D] * 
								r*dtheta * dr * dhdt_total *gap.dt)*(cp*gap.fields.T_prev.in[id3D] + 
									cT * gap.fields.p_prev.in[id2D]) / gap.dt;
					
					
					
					
					double sp = -1 * ( dr *r* dtheta * dz * gap.fields.rho.in[id3D] * cp / gap.dt);
					
					phid.in[id3D] = phi;  //cp*phi;	// fill the dissipation field

					// --------------------------- cell P ---------------------------- //
					//ap += (aw + ae + an + as + ab + at) + (f_e - f_w + f_n - f_s);
					ap += (aw + ae + an + as + ab + at - sp) + cp *(f_e - f_w + f_n - f_s);

					// --------------------- fill the fv matrix ---------------------- //

					if(aw != 0)
						A(f, T.mesh->elements[e.w].fluid_id) = -aw;
					if(ae != 0)
						A(f, T.mesh->elements[e.e].fluid_id) = -ae;
					if(as != 0)
						A(f, T.mesh->elements[e.s].fluid_id) = -as;
					if(an != 0)
						A(f, T.mesh->elements[e.n].fluid_id) = -an;
					if(ab != 0)
						A(f, T.mesh->elements[e.b].fluid_id) = -ab;
					if(at != 0)
						A(f, T.mesh->elements[e.t].fluid_id) = -at;

					A(f,f) = ap;


									
					f++;	// increment the fluid id counter
				}
				else
				{
					phid.in[id3D] = 0;
				}
			}
		}
	}

	// update the linear system definition
	solver.update();
	
	// clear the matrix
	gmm::clear(A);

}
// ------------------------------------------------------------------------- //
int Energy::solve()
{

	int iters = 0;
	
	//solver.system.quiet = false;
	//solver.system.itermax = 1000;


	discretize();
	double res = solver.residual();

	if(res > solver.system.tol)
	{

		solver.initialize_x(); // initialize the solution vector
		
		solver.preconditioner();
			
		iters = solver.solve();

		// update the fields.T field
		for(int i = 0, f = 0; i < T.mesh->mnq; i++) 
		{
			if(T.mesh->elements[i].ty == FLUID)
				T.in[i] = solver.x[f++];
		}
	}

	return iters;
}
// ------------------------------------------------------------------------- //
