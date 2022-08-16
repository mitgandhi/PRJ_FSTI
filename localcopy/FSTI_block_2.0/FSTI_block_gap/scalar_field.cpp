# include "./scalar_field.h"
# include "../FSTI_Block_dll/log.h"
# include <fstream>

# include <gsl/gsl_errno.h>
# include <gsl/gsl_spline.h>
# include <gsl/gsl_math.h>
# include <gsl/gsl_deriv.h>

extern class gaplog Log;

// ------------------------------------------------------------------------- //
scalar_field::scalar_field() : mesh(0)
{
	m = 0;
	n = 0;
	q = 0;
	mn = 0;
	nq = 0;
	mnq = 0;

	in.resize(0);
	bound.inner.resize(0);
	bound.outer.resize(0);
	bound.top.resize(0);
	bound.bottom.resize(0);
}
// ------------------------------------------------------------------------- //
scalar_field::scalar_field(const gap_mesh *msh, double val) : mesh(msh)
{
	m = mesh -> m;
	n = mesh -> n;
	q = mesh -> q;
	mn = m*n;
	nq = n*q;
	mnq = m*n*q;

	in.resize(mnq, val);
	bound.inner.resize(nq, val);
	bound.outer.resize(nq, val);
	bound.top.resize(mn, val);
	bound.bottom.resize(mn, val);


}
// ------------------------------------------------------------------------- //
void scalar_field::initialize(const gap_mesh* msh, double val)
{
	mesh = msh;
	m = mesh -> m;
	n = mesh -> n;
	q = mesh -> q;
	mn = m*n;
	nq = n*q;
	mnq = m*n*q;

	in.resize(mnq, val);
	bound.inner.resize(nq, val);
	bound.outer.resize(nq, val);
	bound.top.resize(mn, val);
	bound.bottom.resize(mn, val);
}
// ------------------------------------------------------------------------- //
scalar_field::scalar_field(const scalar_field& s)
{
	mesh = s.mesh;
	m = s.m;
	n = s.n;
	q = s.q;
	mn = s.mn;
	nq = s.nq;
	mnq = s.mnq;

	in.resize(mnq);
	bound.inner.resize(s.bound.inner.size());
	bound.outer.resize(s.bound.outer.size());
	bound.top.resize(mn);
	bound.bottom.resize(mn);

	for(register int i=0; i<mnq; i++)
		in[i] = s.in[i];
	for(register int i=0; i<mn; i++)
	{
		bound.bottom[i] = s.bound.bottom[i];
		bound.top[i] = s.bound.top[i];
	}
	for(register int i=0; i<nq; i++)
	{
		bound.inner[i] = s.bound.inner[i];
		bound.outer[i] = s.bound.outer[i];
	}
}
// ---------------------------- operator = --------------------------------- //
scalar_field& scalar_field::operator=(const scalar_field &rhs)
{
	mesh = rhs.mesh;
	m = rhs.m;
	n = rhs.n;
	q = rhs.q;
	mn = rhs.mn;
	nq = rhs.nq;
	mnq = rhs.mnq;

	in.resize(mnq);
	bound.inner.resize(rhs.bound.inner.size());
	bound.outer.resize(rhs.bound.outer.size());
	bound.top.resize(mn);
	bound.bottom.resize(mn);

	for(register int i=0; i<mnq; i++)
		in[i] = rhs.in[i];
	for(register int i=0; i<mn; i++)
	{
		bound.bottom[i] = rhs.bound.bottom[i];
		bound.top[i] = rhs.bound.top[i];
	}
	for(register int i=0; i<nq; i++)
	{
		bound.inner[i] = rhs.bound.inner[i];
		bound.outer[i] = rhs.bound.outer[i];
	}

	return *this;
}
// ------------------------------------------------------------------------- //
scalar_field& scalar_field::operator=(double rhs)
{
	for(register int i=0; i<mnq; i++)
		in[i] = rhs;
	for(register int i=0; i<mn; i++)
	{
		bound.bottom[i] = rhs;
		bound.top[i] = rhs;
	}
	for(register int i=0; i<nq; i++)
	{
		bound.inner[i] = rhs;
		bound.outer[i] = rhs;
	}

	return *this;

}
// ------------------------------------------------------------------------- //
scalar_field::~scalar_field() {}
// ------------------------------------------------------------------------- //
scalar_field scalar_field::cshift(int steps) const
{

	scalar_field tmp(*this);

	if(steps > 0) // clockwise rotation
	{
		// rotate the internal field 
		for(int k=0; k<q; k++)
		{
			for(int i=0; i<mn; i++)
				tmp.in[i] = in[k*mn + n*(i/n) + (i + n - steps)%n];
		}
		
		// rotate the top and bottom boundaries
		for(int i=0; i<mn; i++)
		{
			tmp.bound.bottom[i] = bound.bottom[n*(i/n) + (i + n - steps)%n];
			tmp.bound.top[i] = bound.top[n*(i/n) + (i + n - steps)%n];
		}
		
		// rotate the inner and outer boundaries
		for(int i=0; i<nq; i++)
		{
			tmp.bound.inner[i] = bound.inner[(i + n - steps)%n];
			tmp.bound.outer[i] = bound.outer[(i + n - steps)%n];
		}
	}
	else if(steps < 0) // counterclockwise rotation
	{
		// rotate the internal field 
		for(int k=0; k<q; k++)
		{
			for(int i=0; i<mn; i++)
				tmp.in[i] = in[k*mn + n*(i/n) + (i - steps)%n];
		}
		
		// rotate the top and bottom boundaries
		for(int i=0; i<mn; i++)
		{
			tmp.bound.bottom[i] = bound.bottom[n*(i/n) + (i - steps)%n];
			tmp.bound.top[i] = bound.top[n*(i/n) + (i - steps)%n];
		}
		
		// rotate the inner and outer boundaries
		for(int i=0; i<nq; i++)
		{
			tmp.bound.inner[i] = bound.inner[(i - steps)%n];
			tmp.bound.outer[i] = bound.outer[(i - steps)%n];
		}
	}

	return tmp;

}
// ------------------------------------------------------------------------- //
//// calculate the gradient of the scalar field
//vectorField scalar_field::grad()
//{
//	// The gradient in polar coordinates has this form:
//	//
//	// (grad(phi)) = grad(phi)_r r^ + grad(phi)_theta theta^ 
//	//			   = diffp(phi,r) r^ + (1/r)diffp(phi,theta) theta^
//	//
//	// where 
//	//
//	// * grad() is the gradient vector, which has 2 coordinates
//	// * diffp(phi,r) is the partial derivative of phi respect r
//	// * diffp(phi,theta) is the partial derivative of phi respect theta
//	// * r^ is the versor of r direction
//	// * theta^ is the versor of theta direction
//	//
//	// if we use the central difference scheme, at the point P we have
//	//
//	// (grad(phi))_P = (phi_N-phi_S)/(2delta_r) r^ + 
//	//  + (1/r_P)*(phi_E - phi_W)/(2*delta_theta) theta^
//	//
//
//	// on the boundary the radial component of gradient is calculated using 
//	// forward/backward difference, so
//	//
//	// (grad(phi))_P = (phi_N-phi_P)/(delta_r) r^ + 
//	//  + (1/r_P)*(phi_E-phi_W)/(2*delta_theta) theta^
//	// (grad(phi))_P = (phi_P-phi_S)/(delta_r) r^ + 
//	//  + (1/r_P)*(phi_E-phi_W)/(2*delta_theta) theta^
//
//	vectorField gradient(mesh);
//
//	for(int id = 0; id < mnq; id++)
//	{
//		
//		gradient.in[id] = 0;
//
//		gapelement e = mesh -> elements[id];
//
//		if(e.ty == FLUID)
//		{
//		
//			double r = e.center.r;
//			double dr = e.dimension.r;
//			double dz = e.dimension.z;
//			double dtheta = e.dimension.theta;
//			int j = id%n;
//
//			// circumferential direction
//			if(mesh->elements[e.e].ty == FLUID && mesh->elements[e.w].ty == FLUID)
//				gradient.in[id].theta = (in[e.e] - in[e.w])/(2.0*r*dtheta);
//			// opening on the right
//			else if(mesh->elements[e.e].ty != FLUID && mesh->elements[e.w].ty == FLUID)
//				gradient.in[id].theta = (in[id] - in[e.w])/(r*dtheta);
//			// opening on the left
//			else if (mesh->elements[e.e].ty != FLUID && mesh->elements[e.w].ty == FLUID)
//				gradient.in[id].theta = (in[e.e] - in[id])/(r*dtheta);
//
//
//			// ------------------------- radial direction ------------------------- //
//			if (e.n > -1 && e.s > -1)
//			{
//				// all fluid domain
//				if(mesh->elements[e.n].ty == FLUID && mesh->elements[e.s].ty == FLUID)
//					gradient.in[id].r = (in[e.n] - in[e.s])/(2*dr);
//				// opening above
//				else if(mesh->elements[e.n].ty != FLUID && mesh->elements[e.s].ty == FLUID)
//					gradient.in[id].r = (in[id] - in[e.s])/(dr);
//				// opening belove
//				else if(mesh->elements[e.n].ty == FLUID && mesh->elements[e.s].ty != FLUID)
//					gradient.in[id].r = (in[e.n] - in[id])/(dr);
//			}
//			else if (e.n == -1)
//			{
//				gradient.in[id].r = (in[id] - in[e.s])/(dr);
//			}
//			else if (e.s == -1)
//			{
//				gradient.in[id].r = (in[e.n] - in[id])/(dr);
//			}
//		}
//		else
//		{
//			gradient.in[id] = 0;
//		}
//		
//	}
//
//	// set the boundaries
//	for(int k=0, id = 0; k<q; k++)
//	{
//		for(int j=0; j<n; j++, id++)
//		{
//			gradient.bound.inner[id] = gradient.in[mn*k + j];
//			gradient.bound.outer[id] = gradient.in[mn*k + (m-1)*n + j];
//		}
//	}
//
//	// three dimensional field, calc grad_z
//	// fot three dimensional fields the gradient in z direction is not calculated, 
//	// because is most of the time required only at the boundary heat flux and
//	// torque losses...
//	// if grad_z is required, please use the zgrad class, which is design 
//	// calculate the gradient in z direction very accurately
//
//	return gradient;
//}
// ------------------------------- norm2 ----------------------------------- //
double scalar_field::norm2() const
{
	double norm = 0;
	for(int i=0; i<mn; i++)
		norm += in[i]*in[i];

	return sqrt(norm);
}
// ---------------------------------- min ---------------------------------- //
double scalar_field::min() const
{
	double min = in[0];
	for(int i=1; i<mnq; i++)
	{
		if(in[i] < min)
			min = in[i];
	}

	return min;
}
// ---------------------------------- max ---------------------------------- //
double scalar_field::max() const
{
	double max = in[0];
	for(int i=1; i<mnq; i++)
	{
		if(in[i] > max)
			max = in[i];
	}

	return max;
}
// ---------------------------------- avg ---------------------------------- //
double scalar_field::avg() const
{
	double avg = 0;
	for(int i=0; i<mnq; i++)
		avg += in[i];

	return avg/mnq;
}
// ---------------------------------- sum ---------------------------------- //
double scalar_field::sum() const
{
	double s = 0;
	for(int i=0; i<mnq; i++)
		s += in[i];

	return s;
}
// ---------------------------- operator + --------------------------------- //
scalar_field operator+(const scalar_field& l, const scalar_field& r)
{
	scalar_field tmp(l.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
		tmp.in[i] = l.in[i] + r.in[i];
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i] = l.bound.inner[i] + r.bound.inner[i];
		tmp.bound.outer[i] = l.bound.outer[i] + r.bound.outer[i];
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i] = l.bound.bottom[i] + r.bound.bottom[i];
		tmp.bound.top[i] = l.bound.top[i] + r.bound.top[i];
	}

	return tmp;
}
// ------------------------------------------------------------------------- //
scalar_field operator+(double l, const scalar_field& r)
{
	scalar_field tmp(r.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
		tmp.in[i] = l + r.in[i];
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i] = l + r.bound.inner[i];
		tmp.bound.outer[i] = l + r.bound.outer[i];
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i] = l + r.bound.bottom[i];
		tmp.bound.top[i] = l + r.bound.top[i];
	}

	return tmp;
}
// ------------------------------------------------------------------------- //
scalar_field operator+(const scalar_field& l, double r)
{
	scalar_field tmp(l.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
		tmp.in[i] = l.in[i] + r;
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i] = l.bound.inner[i] + r;
		tmp.bound.outer[i] = l.bound.outer[i] + r;
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i] = l.bound.bottom[i] + r;
		tmp.bound.top[i] = l.bound.top[i] + r;
	}

	return tmp;
}
// ---------------------------- operator - --------------------------------- //
scalar_field operator-(const scalar_field& l, const scalar_field& r)
{
	scalar_field tmp(l.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
		tmp.in[i] = l.in[i] - r.in[i];
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i] = l.bound.inner[i] - r.bound.inner[i];
		tmp.bound.outer[i] = l.bound.outer[i] - r.bound.outer[i];
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i] = l.bound.bottom[i] - r.bound.bottom[i];
		tmp.bound.top[i] = l.bound.top[i] - r.bound.top[i];
	}

	return tmp;
}
// ------------------------------------------------------------------------- //
scalar_field operator-(double l, const scalar_field& r)
{
	scalar_field tmp(r.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
		tmp.in[i] = l - r.in[i];
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i] = l - r.bound.inner[i];
		tmp.bound.outer[i] = l - r.bound.outer[i];
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i] = l - r.bound.bottom[i];
		tmp.bound.top[i] = l - r.bound.top[i];
	}

	return tmp;
}
// ------------------------------------------------------------------------- //
scalar_field operator-(const scalar_field& l, double r)
{
	scalar_field tmp(l.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
		tmp.in[i] = l.in[i] - r;
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i] = l.bound.inner[i] - r;
		tmp.bound.outer[i] = l.bound.outer[i] - r;
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i] = l.bound.bottom[i] - r;
		tmp.bound.top[i] = l.bound.top[i] - r;
	}

	return tmp;
}
// ---------------------------- operator * --------------------------------- //
scalar_field operator*(const scalar_field& l, const scalar_field& r)
{
	scalar_field tmp(l.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
		tmp.in[i] = l.in[i] * r.in[i];
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i] = l.bound.inner[i] * r.bound.inner[i];
		tmp.bound.outer[i] = l.bound.outer[i] * r.bound.outer[i];
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i] = l.bound.bottom[i] * r.bound.bottom[i];
		tmp.bound.top[i] = l.bound.top[i] * r.bound.top[i];
	}

	return tmp;
}
// ------------------------------------------------------------------------- //
scalar_field operator*(double l, const scalar_field& r)
{
	scalar_field tmp(r.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
		tmp.in[i] = l * r.in[i];
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i] = l * r.bound.inner[i];
		tmp.bound.outer[i] = l * r.bound.outer[i];
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i] = l * r.bound.bottom[i];
		tmp.bound.top[i] = l * r.bound.top[i];
	}

	return tmp;
}
// ------------------------------------------------------------------------- //
scalar_field operator*(const scalar_field& l, double r)
{
	scalar_field tmp(l.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
		tmp.in[i] = l.in[i] * r;
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i] = l.bound.inner[i] * r;
		tmp.bound.outer[i] = l.bound.outer[i] * r;
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i] = l.bound.bottom[i] * r;
		tmp.bound.top[i] = l.bound.top[i] * r;
	}

	return tmp;
}
// ---------------------------- operator / --------------------------------- //
scalar_field operator/(const scalar_field& l, const scalar_field& r)
{
	scalar_field tmp(l.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
		tmp.in[i] = l.in[i] / r.in[i];
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i] = l.bound.inner[i] / r.bound.inner[i];
		tmp.bound.outer[i] = l.bound.outer[i] / r.bound.outer[i];
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i] = l.bound.bottom[i] / r.bound.bottom[i];
		tmp.bound.top[i] = l.bound.top[i] / r.bound.top[i];
	}

	return tmp;
}
// ------------------------------------------------------------------------- //
scalar_field operator/(const scalar_field& l, double r)
{
	scalar_field tmp(l.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
		tmp.in[i] = l.in[i] / r;
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i] = l.bound.inner[i] / r;
		tmp.bound.outer[i] = l.bound.outer[i] / r;
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i] = l.bound.bottom[i] / r;
		tmp.bound.top[i] = l.bound.top[i] / r;
	}

	return tmp;
}
// ------------------------------------------------------------------------- //
void scalar_field::limit(double min, double max)
{
	for(register int i=0; i<mnq; i++)
	{
		if(in[i] < min)
			in[i] = min;
		if(in[i] > max)
			in[i] = max;
	}
}
// ------------------------------------------------------------------------- //
double scalar_field::get_zgrad_ijk(int i, int j, int k) const
{

	int id = mn*k + n*i + j;	// element id
	
	if(k > 0 && k < q - 1)			// use central difference
	{
		return (in[id + mn] - in[id - mn])/(2*mesh->elements[id].dz);
	}
	else if(k == 0)						// use forward difference
	{
		return (in[id + mn] - in[id])/mesh->elements[id].dz;
	}
	else if(k == q - 1)					// use backward difference
	{
		return (in[id] - in[id - mn])/mesh->elements[id].dz;
	}

}
// ------------------------------------------------------------------------- //
double scalar_field::get_zgrad_bottom(int i, int j) const
{
	int id = n*i + j;
	return (in[id] - bound.bottom[id])/(0.5*mesh->elements[id].dz);
}
// ------------------------------------------------------------------------- //
double scalar_field::get_zgrad_top(int i, int j) const
{
	int id = n*i + j;
	return (bound.top[id] - in[(q-1)*mn + id])/(0.5*mesh->elements[(q-1)*mn + id].dz);
}
// ------------------------------------------------------------------------- //
void scalar_field::save(const char *filename)
{
	std::ofstream out(filename);
	
	out << m << std::endl;
	out << n << std::endl;
	out << q << std::endl;

	out << "internal_field" << std::endl;
	for(int i=0; i<mnq; i++)
		out << in[i] << std::endl;
	out << "inner_boundary" << std::endl;
	for(int j=0; j<nq; j++)
		out << bound.inner[j] << std::endl;
	out << "outer_boundary" << std::endl;
	for(int j=0; j<nq; j++)
		out << bound.outer[j] << std::endl;
	if(q>1)
	{
		out << "bottom_boundary" << std::endl;
		for(int i=0; i<mn; i++)
			out << bound.bottom[i] << std::endl;
		out << "bottom_boundary" << std::endl;
		for(int i=0; i<mn; i++)
			out << bound.top[i] << std::endl;
	}

	out.close();
}
// ------------------------------------------------------------------------- //
void scalar_field::load(const char *filename)
{
	std::ifstream input(filename);
	
	if(!input.is_open())
	{
		Log << "scalar_field::load file name " << filename << " not found!" 
			<< gaplog::endl;
		exit(1);
	}
	
	int _m,_n,_q;
	
	input >> _m;
	input >> _n;
	input >> _q;

	if(_m != m || _n != n || _q != q)
	{
		Log << "scalar_field::load the scalar field to load has an "
						  << "inconsistent dimension!" 
							<< gaplog::endl;
		exit(1);
	}

	std::string tmp;

	// internal field
	input >> tmp;
	if(tmp.compare("internal_field") != 0)
	{
		Log << "scalar_field::load internal_field expected!"
							<< gaplog::endl;
		exit(1);
	}
	for(int i=0; i<mnq; i++)
		input >> in[i];

	// inner boundary
	input >> tmp;
	if(tmp.compare("inner_boundary") != 0)
	{
		Log << "scalar_field::load inner_boundary expected!"
							<< gaplog::endl;
		exit(1);
	}
	for(int j=0; j<nq; j++)
		input >> bound.inner[j];

	// outer boundary
	input >> tmp;
	if(tmp.compare("outer_boundary") != 0)
	{
		Log << "scalar_field::load inner_boundary expected!"
							<< gaplog::endl;
		exit(1);
	}
	for(int j=0; j<nq; j++)
		input >> bound.outer[j];

	if(q>1)
	{
		// bottom boundary
		input >> tmp;
		if(tmp.compare("bottom_boundary") != 0)
		{
			Log << "scalar_field::load inner_boundary expected!"
								<< gaplog::endl;
			exit(1);
		}
		for(int i=0; i<mn; i++)
			input >> bound.bottom[i];
		// top boundary
		input >> tmp;
		if(tmp.compare("top_boundary") != 0)
		{
			Log << "scalar_field::load inner_boundary expected!"
								<< gaplog::endl;
			exit(1);
		}
		for(int i=0; i<mn; i++)
			input >> bound.top[i];
	}

	input.close();
}
// ------------------------------------------------------------------------- //