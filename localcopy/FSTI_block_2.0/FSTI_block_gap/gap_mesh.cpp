# include <cmath>
# include <fstream>
# include "./gap_mesh.h"
# include "./gap_vectors.h"
# include "../FSTI_Block_dll/log.h"

# define ELM_NDS 8
# define pi 3.14159265358979323846

// the stl file can be used either to define the fluid region,
// or to define the openings and ports
// 
# define USE_STL_FLUID

using namespace std;

extern class gaplog Log;

// ------------------------------------------------------------------------- //
void gap_mesh::build(const input& in, int _dim)
{

	// ------------------------ copy from caspar_input ----------------------- //

	int z = in.data.operating_conditions.npistons;
	const char* cb_stl_file = in.data.options_block.fluid_grid.stl_cb.c_str();
	const char* vp_stl_file = in.data.options_block.fluid_grid.stl_vp.c_str();
	
	r_in = 0.5*in.data.geometry.d_gap_in;
	r_openin = 0.5*in.data.geometry.d_ope_in;
	r_openout = 0.5*in.data.geometry.d_ope_out;
	r_groovein = 0.5*in.data.geometry.d_groove_in;
	r_grooveout = 0.5*in.data.geometry.d_groove_out;
	r_out = 0.5*in.data.geometry.d_gap_out;

	n = in.data.options_block.fluid_grid.N;
	m = in.data.options_block.fluid_grid.M;
	q = in.data.options_block.fluid_grid.Q;


	// ------------------------------------------------------------------------ //

	dim = _dim;	// define the dimension
	if(dim != 2 && dim != 3)
	{
		Log << "\nWrong gap mesh dimension, use 2 or 3" << gaplog::endl;
		exit(1);
	}

	// update q in the case of a 2D mesh
	if(dim == 2)
		q = 1;
	
	// update mesh size
	mn = m*n;
	mnq = mn*q;
		
	// no hydrodynamic bearing
	if(r_groovein == 0 || r_grooveout == 0)
	{
		outerbearing = false;
		mi[0] = static_cast<int>(floor(m*(r_openin - r_in)/(r_out - r_in)));
		mi[1] = static_cast<int>(floor(m*(r_openout - r_openin)/(r_out - r_in)));
		mi[2] = m - (mi[0] + mi[1]);
		mi[3] = 0;
		mi[4] = 0;
	}
	// with hydrodynamic bearing
	else
	{
		outerbearing = true;
		double groove_width = r_grooveout - r_groovein;
		mi[0] = static_cast<int>(floor(m*(r_openin - r_in)/((r_out - r_in) - groove_width)));
		mi[1] = static_cast<int>(floor(m*(r_openout - r_openin)/((r_out - r_in) - groove_width)));
		mi[2] = static_cast<int>(floor(m*(r_groovein - r_openout)/((r_out - r_in) - groove_width)));
		mi[3] = 1;
		mi[4] = m - (mi[0] + mi[1] + mi[2] + mi[3]);
	}

	// define outerbearing flag
	if(r_groovein == r_grooveout == 0)
		outerbearing = true;
	else
		outerbearing = false;
	

	// define elements dimension
	dr[0] = (r_openin - r_in)/mi[0];
	dr[1] = (r_openout - r_openin)/mi[1];
	if(outerbearing)
	{
		dr[2] = (r_groovein - r_openout)/mi[2];
		dr[3] = (r_grooveout - r_groovein);
		dr[4] = (r_out - r_grooveout)/mi[4];
	}
	else
	{
		dr[2] = (r_out - r_openout)/mi[2];
		dr[3] = dr[4] = 0;
	}
	
	dtheta = 2.0*pi/n;

	dz.resize(mn);
	for(int i=0; i<mn; i++)
		dz[i] = 10e-6/q;	// set default to 10 um

	rotation_steps = 0;

	// ----------------------- populate the nodes vector --------------------- //
	
	// loop in axial direction
	for(int k=0; k<q+1; k++)
	{
		double r = r_in;
		// loop throughout the radial sectors
		for(int im=0; im<5; im++)
		{
			// loop throughout the elements in the radial sector
			for(int i=0; i<mi[im]; i++)
			{
				// loop in circumferential direction
				for(int j=0; j<n; j++)
				{
					double theta = (j-0.5)*dtheta;
					double x = r*sin(theta);
					double y = r*cos(theta);
					double z = k*dz[0];
					nodes.push_back(point(x,y,z));
				}
				r += dr[im];
			}
		}
		// add the last circle
		for(int j=0; j<n; j++)
		{
			double theta = (j-0.5)*dtheta;
			double x = r*sin(theta);
			double y = r*cos(theta);
			double z = k*dz[0];
			nodes.push_back(point(x,y,z));
		}
	}

	// ---------------------------- define elements -------------------------- //

	// loop in axial direction
	for(int k=0, e = 0; k<q; k++)
	{
		double r = r_in;
		// loop throughout the radial sectors
		for(int im=0, id2d = 0, c = 0; im<5; im++)
		{
			// loop throughout the elements in the radial sector
			for(int i=0; i<mi[im]; i++, c++)
			{
				// loop in circumferential direction
				for(int j=0; j<n; j++, e++, id2d++)
				{
					gap_elm elm;
					elm.r = r + 0.5*dr[im];
					elm.theta = j*dtheta;
					elm.x = elm.r*sin(elm.theta);
					elm.y = elm.r*cos(elm.theta);
					//cout << id2d << "\t" << dz.size() << endl;
					elm.z = (0.5 + k)*dz[id2d];
					elm.dr = dr[im];
					elm.dtheta = dtheta;
					elm.dz = dz[id2d];
					// element area (ok, checked)
					elm.A = 0.5*(pow(elm.r + 0.5*elm.dr, 2.0) - pow(elm.r - 0.5*elm.dr, 2.0))*dtheta;
					// element volume
					if(q > 1)
						elm.V = elm.A*elm.dz;
					else
						elm.V = 0;
															
					// element nodes definition
					memset(elm.nds, 0, 8);
					if(j < n - 1) 
					{
						// bottom layer (k)
						elm.nds[0] = c*n + j + k*(m + 1)*n;
						elm.nds[1] = c*n + j + 1 + k*(m + 1)*n;
						elm.nds[2] = (c + 1)*n + j + 1 + k*(m + 1)*n;
						elm.nds[3] = (c + 1)*n + j + k*(m + 1)*n;
						// top layer (k + 1)
						elm.nds[4] = c*n + j + (k + 1)*(m + 1)*n;
						elm.nds[5] = c*n + j + 1 + (k + 1)*(m + 1)*n;
						elm.nds[6] = (c + 1)*n + j + 1 + (k + 1)*(m + 1)*n;
						elm.nds[7] = (c + 1)*n + j + (k + 1)*(m + 1)*n;
					}
					else 
					{
						// bottom layer (k)
						elm.nds[0] = c*n + j + k*(m + 1)*n;
						elm.nds[1] = c*n + k*(m + 1)*n;
						elm.nds[2] = (c + 1)*n + k*(m + 1)*n;
						elm.nds[3] = (c + 1)*n + j + k*(m + 1)*n;
						// top layer (k + 1)
						elm.nds[4] = c*n + j + (k + 1)*(m + 1)*n;
						elm.nds[5] = c*n + (k + 1)*(m + 1)*n;
						elm.nds[6] = (c + 1)*n + (k + 1)*(m + 1)*n;
						elm.nds[7] = (c + 1)*n + j + (k + 1)*(m + 1)*n;
					}

					// -------------------------- neighbors -------------------------- //

					// South
					if(c > 0) 
						elm.s = e - n;
					else
						elm.s = -1;
					// North
					if(c < m - 1)
						elm.n = e + n;
					else 
						elm.n = -1;
					// West
					if(j > 0)
						elm.w = e - 1;
					else
						elm.w = e + n - 1;
					// East
					if(j < n - 1)
						elm.e = e + 1;
					else
						elm.e = e - (n - 1);
					// Bottom
					if(k > 0)
						elm.b = e - m*n;
					else
						elm.b = -1;
					// Top
					if(k < q - 1)
						elm.t = e + m*n;
					else
						elm.t = -1;

					//cout << "id: " << e << endl;
					//cout << "w: " << elm.w << "\t" 
					//		 << "e: " << elm.e << "\t"
					//		 << "s: " << elm.s << "\t"
					//		 << "n: " << elm.n << "\t"
					//		 << "b: " << elm.b << "\t"
					//		 << "t: " << elm.t << endl;
					//system("pause");

					// push back into the list
					cb_elms.push_back(FLUID);	
					cb_elms_0.push_back(FLUID);		
					vp_elms.push_back(FLUID);
					
					// define the gap element type
					elm.ty = FLUID;

					// add the element in the list and update counter
					elements.push_back(elm);

				}
				r += dr[im];
			}
		}
	}

	// ---------------------- define element properties ---------------------- //

	// stl geometry for the block
	cb_stl.read(cb_stl_file, 1e-3);
	cb_bound.build(&cb_stl);
	// stl geometry for the valveplate
	vp_stl.read(vp_stl_file, 1e-3);
	vp_bound.build(&vp_stl);

	Log << "  * Defining elements ... ";

	mn_f = 0;
	mnq_f = 0;
	int f_id = 0;		// fluid id
	cb_fluid = 0;
	vp_fluid = 0;

	for(int k=0, e = 0; k<q; k++)
	{
		// loop throughout the radial sectors
		for(int im=0, c = 0; im<5; im++)
		{
			// loop throughout the elements in the radial sector
			for(int i=0; i<mi[im]; i++, c++)
			{
				// loop in circumferential direction
				for(int j=0; j<n; j++, e++)
				{

					// element center
					double cx = elements[e].r*sin(elements[e].theta);
					double cy = elements[e].r*cos(elements[e].theta);
					double r = elements[e].r;
					double theta = elements[e].theta;
					double dr = elements[e].dr;
					
					// cylinder block and valve plate labels
					int cbij, vpij;

					// use the stl geometry just for the first layer
					if(k == 0)
					{
						// get inside/outside information
						// cbij, vpij: 0 -> fluid, !0 -> opening/groove

						# ifdef USE_STL_FLUID
						cbij = !cb_bound.is_inside(point(cx,cy));
						vpij = !vp_bound.is_inside(point(cx,cy));
						# else
						cbij = cb_bound.is_inside(point(cx,cy));
						vpij = vp_bound.is_inside(point(cx,cy));
						# endif
																		
						// fluid label
						int label = FLUID;
						
						// opening angle
						const double aope = 2.0*pi/z; 

						// openings region
						if(r > r_openin && r < r_openout)
						{
							// opening of first piston on cylinder block
							if((theta >=0 && theta < 0.5*aope) || (theta > (2.0*pi - 0.5*aope) && theta < 2*pi)) 
							{
								label = 1;
							}
							// opening for others pistons
							for(int h=1; h<z; h++) 
							{
								// opening of j + 1_th piston on cylinder block
								double ang_min = 0.5*aope + (h - 1)*aope;
								double ang_max = 0.5*aope + h*aope;
								if((theta > ang_min) && (theta <= ang_max))
								{
									label = h + 1;
								}
							}
						}
						else // case
						{	
							label = CASE;
						}

						// define the cylinder block label
						cbij = label*(cbij);	
						
						// update the cb_fluid size
						if(cbij == FLUID)
							cb_fluid++;

						// set gapelement type for valve plate
						if(r > r_openin && r < r_openout) // openings region
						{
							if(cx >= 0)
								label = HP;
							else
								label = LP;
						}
						else // case
						{
							label = CASE;
						}

						// define the valve plate label
						vpij = label*(vpij);	
						
						// update the vp_fluid size
						if(vpij == FLUID)
							vp_fluid++;
						
					}
					else
					{
						// just copy from the first layer
						cbij = cb_elms[e - m*n];
						vpij = vp_elms[e - m*n];
					}

					// push back into the list
					cb_elms[e] = cbij;		
					cb_elms_0[e] = cbij;		
					vp_elms[e] = vpij;
					
					// define the gap element type
					elements[e].ty = cbij + vpij;

					if(elements[e].ty == FLUID)
					{
						mnq_f++;
						elements[e].fluid_id = f_id++;
					}
				}
			}
		}
	}

	mn_f = mnq_f/q;

	Log << "done!" << gaplog::endl;

}
// ------------------------------------------------------------------------- //
void gap_mesh::define_gap_surface(body_type which)
{
	
	// pointer to the gap surface
	vector<polygon>* surf;
	gap_boundary* bound;
	vector<double>* ai;

	if(which == CB)
	{
		Log << "  * Defining cylinder block gap surface ...";
		surf = &cb_gap_surf.elms;
		ai = &cb_gap_surf.ai;
		bound = &cb_bound;
	}
	else
	{
		Log << "  * Defining valve plate gap surface ...";
		surf = &vp_gap_surf.elms;
		ai = &vp_gap_surf.ai;
		bound = &vp_bound;
	}

	surf->resize(0);
	ai->resize(0);
	double atot = 0;
	
	int sz = m*n;

	for(int i=0; i<sz; i++)
	{
		// define a polygon by using the vertices of the current fluid element
		polygon this_elm;
		this_elm.nv = 4;
		this_elm.vertices.resize(this_elm.nv);
		for(int h=0; h<4; h++)
		{
			point V
			(
				nodes[elements[i].nds[h]].x(), 
				nodes[elements[i].nds[h]].y()
			);
			this_elm.vertices[h] = V;
		}

		// determine the intersections of each element with the boundaries
		vector<vector<polygon>> intersections(0);

		for(unsigned int j=0; j<bound->boundaries.size(); j++)
		{
			# ifdef USE_STL_FLUID
			vector<polygon> tmp = this_elm.clip(bound->boundaries[j], GPC_INT);
			# else
			vector<polygon> tmp = this_elm.clip(bound->boundaries[j], GPC_DIFF);
			# endif
			if(tmp.size() > 0)
				intersections.push_back(tmp);
			else
				intersections.push_back(vector<polygon>(0));
		}

		// clip out what doesn't need to be there				
		bool partial_intersection = false;
		bool full_intersection = false;
		for(unsigned int j=0; j<bound->boundaries.size(); j++)
		{
			// at leat one intersection was found
			if(intersections[j].size() > 0)
			{
				double int_area = intersections[j][0].area();
				// the element is fully included into a boundary
				if(int_area/elements[i].A > 0.999)
				{
					full_intersection = true;
				}
				// the element is partially included into a boundary
				else if(int_area/elements[i].A > 0.001 && int_area/elements[i].A <= 0.999)
				{
					partial_intersection = true;

					// find a point inside the intersection
					point d1 = 1e-3*(intersections[j][0].vertices[0] - intersections[j][0].vertices[1]);
					point d2 = 1e-3*(intersections[j][0].vertices[2] - intersections[j][0].vertices[1]);
					point c = intersections[j][0].vertices[1] + (d1 + d2);
					// adjust the position if c lies outside the intersection
					if(!intersections[j][0].is_inside(c))
						c = intersections[j][0].vertices[1] - (d1 + d2);
					
					// if the interior point of the polygon is inside the gap, then add to the list
					# ifdef USE_STL_FLUID
					if(bound->is_inside(c))
					{
						surf->push_back(intersections[j][0]);
						ai->push_back(intersections[j][0].area());
						atot += (*ai)[ai->size()-1];
					}
					# else
					if(!bound->is_inside(c))
					{
						surf->push_back(intersections[j][0]);
						ai->push_back(intersections[j][0].area());
						atot += (*ai)[ai->size()-1];
					}
					# endif
					// the intersection is not inside the gap, then try clipping out
					else
					{
						vector<polygon> diff = this_elm.clip(bound->boundaries[j], GPC_DIFF);
						if(diff.size() > 0)
						{
							for(unsigned int h=0; h<diff.size(); h++)
							{
								//d1 = 1e-2*(diff[h].vertices[0] - diff[0].vertices[1]);
								//d2 = 1e-2*(diff[h].vertices[2] - diff[0].vertices[1]);
								//c = diff[h].vertices[1] + (d1 + d2);
								//if(bound->is_inside(c))
								//{
								surf->push_back(diff[h]);
								ai->push_back(diff[h].area());
								atot += (*ai)[ai->size()-1];
								//}
							}
						}
					}
				}
			}
		}
		if(full_intersection && !partial_intersection)
		{
			# ifdef USE_STL_FLUID
			if(bound->is_inside(this_elm.center()))
			{
				surf->push_back(this_elm);
				ai->push_back(this_elm.area());
				atot += (*ai)[ai->size()-1];
			}
			# else
			if(!bound->is_inside(this_elm.center()))
			{
				surf->push_back(this_elm);
				ai->push_back(this_elm.area());
				atot += (*ai)[ai->size()-1];
			}
			# endif
		}	
		
	}


	// copy cb_gap_surf.elms to cb_gap_surf.elms_0
	if(which == CB)
	{
		cb_gap_surf.elms_0.resize(cb_gap_surf.elms.size());
		for(unsigned int i=0; i<surf->size(); i++)
			cb_gap_surf.elms_0[i] = cb_gap_surf.elms[i];
		//polygon::write_vtk("./cb.vtk", *surf);
	} else {
		//polygon::write_vtk("./vp.vtk", *surf);
	}

	Log << "done! (Total area: " << atot*1e6 << " [mm2])" << gaplog::endl;
	
	
}
// ------------------------------------------------------------------------- //
int gap_mesh::rotate(double angle_rad)
{
	
	double step_angle = 2.0*pi/n;
	int steps = static_cast<int>(floor(angle_rad/step_angle));
	double difference = angle_rad - steps*step_angle;
	
	// if the difference is greather than 0.5 deg, rotate
	if(difference*(180.0/pi) >= 0.50)
		steps += 1;

	// update the rotation steps variable
	rotation_steps = steps;

	if(steps != 0) // perform clock-wise rotation
	{
		// update the element definition
		for(int i = 0; i<mn; i++)
			cb_elms[i] = cb_elms_0[n*(i/n) + (n + i - steps)%n];
		// update the cylinder block surface
		for(unsigned int i=0; i<cb_gap_surf.elms_0.size(); i++)
		{
			for(unsigned int j=0; j<cb_gap_surf.elms_0[i].vertices.size(); j++)
			{
				double x = cb_gap_surf.elms_0[i].vertices[j].x();
				double y = cb_gap_surf.elms_0[i].vertices[j].y();
				double ang = steps*step_angle; // use the same angle
				cb_gap_surf.elms[i].vertices[j][0] = (x*cos(ang) + y*sin(ang));
				cb_gap_surf.elms[i].vertices[j][1] = (-x*sin(ang) + y*cos(ang));
			}
		}
	}
	else	// just copy from zero position
	{
		for(int i = 0; i<mn; i++)
			cb_elms[i] = cb_elms_0[i];
		for(unsigned int i=0; i<cb_gap_surf.elms_0.size(); i++)
			cb_gap_surf.elms[i] = cb_gap_surf.elms_0[i];
	}

	int jump = 0;
	mnq_f = 0;

	// loop through all the elements
	for(int i = 0, j = 0, f_id = 0; i < mnq; i++, j++) 
	{
		if(j == mn)
			j = 0;
			
		// update the gapelement type
		elements[i].ty = cb_elms[j] + vp_elms[j];
		// update the fluid gapelement number
		if(elements[i].ty == FLUID)
		{
			mnq_f++;
			elements[i].fluid_id = f_id++;
		}
		
	}

	// update the layer fluid gapelement number
	mn_f = mnq_f / q;
		
	return steps;
	
	
}
// ------------------------------------------------------------------------- //
void gap_mesh::update_thickness(const scalar_field& top, const scalar_field& bottom)
{
	// update just for 3D mesh
	if(dim == 2)
		return;

	int w, e;	
	int p1, p2, p3, p4;

	// update gapelement dz
	for(int k = 0, id = 0; k < q; k++)
	{
		for(int i = 0; i < m; i++)
		{
			for(int j = 0; j < n; j++, id++)
			{
				elements[id].dz = 
					(top.in[n*i + j] - bottom.in[n*i + j])/q;
			}
		}
	}

	// update nodes coordinates

	int mn = m+1;	// nodes in radial direction
	int nn = n;		// nodes in circumferential direction
	int qn = q+1;	// nodes in axial direction

	for(int i = 0; i < mn; i++)
	{
		for(int j = 0; j < nn; j++)
		{
			// top coordinate for the nodes in position (i,j)
			double tij;
			// bottom coordinate for nodes in position (i,j)
			double bij;
			// used for opening 1 -> node to be used; 0 -> node not to be used
			p1 = p2 = p3 = p4 = 1;
				
			// set the w and e index
			if(j == 0)
				w = n - 1, e = 0;
			else
				w = j - 1, e = j;

			// inner boundary
			if(i == 0)
			{
					
				p3 = static_cast<int>(!elements[w].ty);
				p4 = static_cast<int>(!elements[e].ty);
				p1 = p3;
				p2 = p4;
					
				if(p1 + p2 + p3 + p4 > 0)
				{
					tij = 
					(
						p1*top.bound.inner[w] + p2*top.bound.inner[e] +
						p3*top.in[w] + p4*top.in[e]
					)/(p1 + p2 + p3 + p4);
					//
					bij =
					(
						p1*bottom.bound.inner[w] + p2*bottom.bound.inner[e] +
						p3*bottom.in[w] + p4*bottom.in[e]
					)/(p1 + p2 + p3 + p4);
				}
			}
			// outer boundary
			else if(i == mn - 1)
			{
				
				p3 = static_cast<int>(!elements[(i-1)*n + w].ty);
				p4 = static_cast<int>(!elements[(i-1)*n + e].ty);
				p1 = p3;
				p2 = p4;

				if(p1 + p2 + p3 + p4 > 0)
				{
					tij = 
					(
						p1*top.bound.outer[w] + p2*top.bound.outer[e] +
						p3*top.in[(i-1)*n + w] + p4*top.in[(i-1)*n + e]				
					)/(p1 + p2 + p3 + p4);
					//
					bij = 
					(
						p1*bottom.bound.outer[w] + p2*bottom.bound.outer[e] +
						p3*bottom.in[(i-1)*n + w] + p4*bottom.in[(i-1)*n + e]
					)/(p1 + p2 + p3 + p4);
				}
			}
			// internal mesh
			else
			{
				// find if any of the four elements are inside an opening
				p1 = static_cast<int>(!elements[(i - 1)*n + w].ty);
				p2 = static_cast<int>(!elements[(i - 1)*n + e].ty);
				p3 = static_cast<int>(!elements[i*n + w].ty);
				p4 = static_cast<int>(!elements[i*n + e].ty);

				if(p1 + p2 + p3 + p4 > 0)
				{
					// top side
					tij = 
					(
						p1*top.in[(i - 1)*n + w] +
						p2*top.in[(i - 1)*n + e] +
						p3*top.in[i*n + w] +
						p4*top.in[i*n + e]
					)/(p1 + p2 + p3 + p4);
					// bottom side
					bij = 
					(
						p1*bottom.in[(i - 1)*n + w] +
						p2*bottom.in[(i - 1)*n + e] +
						p3*bottom.in[i*n + w] +
						p4*bottom.in[i*n + e]
					)/(p1 + p2 + p3 + p4);
				}
			}
				
			// update the node coordinate
			if(p1 + p2 + p3 + p4 > 0)
			{
				double dnz = (tij - bij)/q;
				for(int k = 0; k < qn; k++)
					nodes[mn*nn*k + nn*i + j][2] = bij + k*dnz;
			}
		}
	}	

	// update elements center and volume
	for(int i = 0; i < mnq; i++)
	{
		// update gapelement center
		point c(0,0,0);
		for(int j = 0; j < 8; j++)
			c = c + nodes[elements[i].nds[j]];
		c = c/8.0;
		
		// store the new center in cylindrical coordinates
		ca_vector elmc(c.x(), c.y(), c.z());
		elements[i].r = elmc.r();
		elements[i].theta = elmc.theta();
		elements[i].z = elmc.z();

		// update gap element area
		elements[i].A =
		0.5*
		(
			pow(elements[i].r + 0.5*elements[i].dr, 2.0) -
			pow(elements[i].r - 0.5*elements[i].dr, 2.0) 
		)*elements[i].dtheta;

		// update gap element volume
		elements[i].V = elements[i].A*elements[i].dz;

	}

}
// ------------------------------------------------------------------------- //
void gap_mesh::write_gap_surface_vtk(body_type which)
{

	string name;
	const vector<polygon>* surf;

	if(which == CB)
	{
	 surf = &cb_gap_surf.elms;
	 name = "filename";
	}
	else
	{
	 surf = &vp_gap_surf.elms;
	 name = "valveplate_gap_surface.vtk";
	}
	
	
	Log << "  * Writing vtk ... ";

	vector<int> nodes_idx(surf->size()); // this is the staring index vector
	nodes_idx[0] = 0;
	int ntot = (*surf)[0].vertices.size();

	for(unsigned int i=1; i<surf->size(); i++)
	{
		nodes_idx[i] = nodes_idx[i-1] + (*surf)[i-1].vertices.size();
		ntot += (*surf)[i].vertices.size();
	}

	ofstream vtk(name.c_str());
	if (!vtk.is_open()) 
	{
		cout << "Error opening " << name << ".vtk" << endl;
		exit(1);
	}

	// write VTK header
	vtk << "# vtk DataFile Version 2.0" << endl 
			<< "vtk output" << endl
			<< "ASCII" << endl 
			<< "DATASET UNSTRUCTURED_GRID" << endl 
			<< "POINTS " << ntot << " double" << endl;

	int cell_size = 0;

	for(unsigned int i=0; i<surf->size(); i++)
	{
		cell_size += (1 + 1 + (*surf)[i].vertices.size());
		for(unsigned j=0; j<(*surf)[i].vertices.size(); j++)
		{
			vtk << (*surf)[i].vertices[j].x() << "\t" 
					<< (*surf)[i].vertices[j].y() << "\t" 
					<< 0 << endl;
		}
	}

	vtk << "CELLS " << surf->size() << "\t" << cell_size << endl;

	for(unsigned int i = 0; i<surf->size(); i++) 
	{
		vtk << (*surf)[i].vertices.size() + 1 << "\t";
		for(unsigned int j = 0; j<(*surf)[i].vertices.size(); j++)
			vtk << nodes_idx[i] + j << "\t"; // index start from 0
		vtk << nodes_idx[i] << endl; //repeat the first
	}
	
	vtk << endl;

	vtk << "CELL_TYPES " << surf->size() << endl;
	for(unsigned int i=0; i<surf->size(); i++) 
	{
		vtk << 7 << endl;
	}

	//vtk << "\nCELL_DATA " << mnq << endl;
	//vtk << "\nSCALARS gap float" << endl;
	//vtk << "\nLOOKUP_TABLE default" << endl;

	//for(int i=0; i<mnq; i++)
	//	vtk << elms[i] << endl;
	//	//vtk << elements[i].ty << endl;

	Log << "done!" << gaplog::endl;

	vtk.close();
}
// ------------------------------------------------------------------------- //
void gap_mesh::write_vtk(const char* file, double scale)
{

	Log << "Writing vtk ... ";

	ofstream vtk(file);
	if (!vtk.is_open()) 
	{
		cout << "Error opening " << file << " .vtk" << endl;
		exit(1);
	}

	// write VTK header
	vtk << 
		"# vtk DataFile Version 2.0" << endl <<
		"vtk output" << endl <<
		"ASCII" << endl <<
		"DATASET UNSTRUCTURED_GRID" << endl << 
		"POINTS " << nodes.size() << " double" << endl;

	for(unsigned int i=0; i<nodes.size(); i++)
	{
		vtk << nodes[i].x() << "\t" << nodes[i].y() << "\t" << scale*nodes[i].z() << endl;
	}

	vtk << "CELLS " << mnq << "\t" << (1 + ELM_NDS)*mnq << endl;

	for(int i = 0; i <mnq; i++) 
	{
		vtk << ELM_NDS << "\t";
		for(int j = 0; j < ELM_NDS; j++)
			vtk << elements[i].nds[j] << "\t"; // index start from 0
		vtk << endl;
	}

	vtk << "CELL_TYPES " << mnq << endl;
	for(int i = 0; i <mnq; i++) 
	{
		vtk << 12 << endl;
	}

	vtk << "\nCELL_DATA " << mnq << endl;
	vtk << "\nSCALARS gap float" << endl;
	vtk << "\nLOOKUP_TABLE default" << endl;

	for(int i=0; i<mnq; i++)
		vtk << elements[i].ty << endl;

	Log << "done!" << gaplog::endl;

	vtk.close();
}
// ------------------------------------------------------------------------- //