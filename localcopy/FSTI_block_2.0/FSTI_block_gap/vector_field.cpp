# include "vector_field.h"

// ------------------------------------------------------------------------- //
vector_field::vector_field() : mesh(0)
{
	m = 0;
	n = 0;
	q = 0;
	mn = 0;
	mnq = 0;

	in.resize(0);
	bound.inner.resize(0);
	bound.outer.resize(0);
	bound.top.resize(0);
	bound.bottom.resize(0);
}
// ------------------------------------------------------------------------- //
vector_field::vector_field(const gap_mesh *msh, cy_vector val) : mesh(msh)
{
	m = mesh -> m;
	n = mesh -> n;
	q = mesh -> q;
	nq = n*q;
	mn = m*n;
	mnq = m*n*q;

	in.resize(mnq, val);
	bound.inner.resize(nq, val);
	bound.outer.resize(nq, val);
	bound.top.resize(mn, val);
	bound.bottom.resize(mn, val);
}
// ------------------------------------------------------------------------- //
void vector_field::initialize(const gap_mesh* msh, cy_vector val)
{
	
	mesh = msh;
	m = mesh -> m;
	n = mesh -> n;
	q = mesh -> q;
	nq = n*q;
	mn = m*n;
	mnq = m*n*q;

	in.resize(mnq, val);
	bound.inner.resize(nq, val);
	bound.outer.resize(nq, val);
	bound.top.resize(mn, val);
	bound.bottom.resize(mn, val);
}
// ------------------------------------------------------------------------- //
vector_field::vector_field(const vector_field& s)
{
	mesh = s.mesh;
	m = s.m;
	n = s.n;
	q = s.q;
	mn = s.mn;
	mnq = s.mnq;

	in.resize(mnq);
	bound.inner.resize(n);
	bound.outer.resize(n);
	bound.top.resize(mn);
	bound.bottom.resize(mn);

	for(register int i=0; i<mnq; i++)
		in[i] = s.in[i];
	for(register int i=0; i<mn; i++)
	{
		bound.bottom[i] = s.bound.bottom[i];
		bound.top[i] = s.bound.top[i];
	}
	for(register int i=0; i<n; i++)
	{
		bound.inner[i] = s.bound.inner[i];
		bound.outer[i] = s.bound.outer[i];
	}
}
// ------------------------------------------------------------------------- //
vector_field::~vector_field() {}
// ------------------------------------------------------------------------- //
scalar_field vector_field::r() const
{
	scalar_field tmp(mesh);

	for(register int i=0; i<mnq; i++)
		tmp.in[i] = in[i].r();
	for(register int i=0; i<mn; i++)
	{
		tmp.bound.bottom[i] = bound.bottom[i].r();
		tmp.bound.top[i] = bound.top[i].r();
	}
	for(register int i=0; i<n; i++)
	{
		tmp.bound.inner[i] = bound.inner[i].r();
		tmp.bound.outer[i] = bound.outer[i].r();
	}

	return tmp;
}
// ------------------------------------------------------------------------- //
scalar_field vector_field::theta() const
{
	scalar_field tmp(mesh);

	for(register int i=0; i<mnq; i++)
		tmp.in[i] = in[i].theta();
	for(register int i=0; i<mn; i++)
	{
		tmp.bound.bottom[i] = bound.bottom[i].theta();
		tmp.bound.top[i] = bound.top[i].theta();
	}
	for(register int i=0; i<n; i++)
	{
		tmp.bound.inner[i] = bound.inner[i].theta();
		tmp.bound.outer[i] = bound.outer[i].theta();
	}

	return tmp;
}
// ------------------------------------------------------------------------- //
scalar_field vector_field::z() const
{
	scalar_field tmp(mesh);

	for(register int i=0; i<mnq; i++)
		tmp.in[i] = in[i].z();
	for(register int i=0; i<mn; i++)
	{
		tmp.bound.bottom[i] = bound.bottom[i].z();
		tmp.bound.top[i] = bound.top[i].z();
	}
	for(register int i=0; i<n; i++)
	{
		tmp.bound.inner[i] = bound.inner[i].z();
		tmp.bound.outer[i] = bound.outer[i].z();
	}

	return tmp;
}
// ---------------------------- operator = --------------------------------- //

vector_field& vector_field::operator=(const vector_field &rhs)
{
	mesh = rhs.mesh;
	m = rhs.m;
	n = rhs.n;
	q = rhs.q;
	mn = rhs.mn;
	mnq = rhs.mnq;

	in.resize(mnq);
	bound.inner.resize(n);
	bound.outer.resize(n);
	bound.top.resize(mn);
	bound.bottom.resize(mn);

	for(register int i=0; i<mnq; i++)
		in[i] = rhs.in[i];
	for(register int i=0; i<mn; i++)
	{
		bound.bottom[i] = rhs.bound.bottom[i];
		bound.top[i] = rhs.bound.top[i];
	}
	for(register int i=0; i<n; i++)
	{
		bound.inner[i] = rhs.bound.inner[i];
		bound.outer[i] = rhs.bound.outer[i];
	}

	return *this;
}
// ------------------------------------------------------------------------- //
vector_field& vector_field::operator=(const cy_vector& rhs)
{
	for(register int i=0; i<mnq; i++)
		in[i] = rhs;
	for(register int i=0; i<mn; i++)
	{
		bound.bottom[i] = rhs;
		bound.top[i] = rhs;
	}
	for(register int i=0; i<n; i++)
	{
		bound.inner[i] = rhs;
		bound.outer[i] = rhs;
	}

	return *this;
}
 // ---------------------------- operator + --------------------------------- //
vector_field operator+(const vector_field& l, const vector_field& r)
{
	vector_field tmp(l.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
	{
		tmp.in[i][0] = l.in[i].r() + r.in[i].r();
		tmp.in[i][1] = l.in[i].theta() + r.in[i].theta();
		tmp.in[i][2] = l.in[i].z() + r.in[i].z();
	}
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i][0] = l.bound.inner[i].r() + r.bound.inner[i].r();
		tmp.bound.inner[i][1] = l.bound.inner[i].theta() + r.bound.inner[i].theta();
		tmp.bound.inner[i][2] = l.bound.inner[i].z() + r.bound.inner[i].z();
		tmp.bound.outer[i][0] = l.bound.outer[i].r() + r.bound.outer[i].r();
		tmp.bound.outer[i][1] = l.bound.outer[i].theta() + r.bound.outer[i].theta();
		tmp.bound.outer[i][2] = l.bound.outer[i].z() + r.bound.outer[i].z();
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i][0] = l.bound.bottom[i].r() + r.bound.bottom[i].r();
		tmp.bound.bottom[i][1] = l.bound.bottom[i].theta() + r.bound.bottom[i].theta();
		tmp.bound.bottom[i][2] = l.bound.bottom[i].z() + r.bound.bottom[i].z();
		tmp.bound.top[i][0] = l.bound.top[i].r() + r.bound.top[i].r();
		tmp.bound.top[i][1] = l.bound.top[i].theta() + r.bound.top[i].theta();
		tmp.bound.top[i][2] = l.bound.top[i].z() + r.bound.top[i].z();
	}

	return tmp;
}
// ------------------------------------------------------------------------- //
vector_field operator+(const scalar_field& l, const vector_field& r)
{
	vector_field tmp(r.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
	{
		tmp.in[i][0] = l.in[i] + r.in[i].r();
		tmp.in[i][1] = l.in[i] + r.in[i].theta();
		tmp.in[i][2] = l.in[i] + r.in[i].z();
	}
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i][0] = l.bound.inner[i] + r.bound.inner[i].r();
		tmp.bound.inner[i][1] = l.bound.inner[i] + r.bound.inner[i].theta();
		tmp.bound.inner[i][2] = l.bound.inner[i] + r.bound.inner[i].z();
		tmp.bound.outer[i][0] = l.bound.outer[i] + r.bound.outer[i].r();
		tmp.bound.outer[i][1] = l.bound.outer[i] + r.bound.outer[i].theta();
		tmp.bound.outer[i][2] = l.bound.outer[i] + r.bound.outer[i].z();
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i][0] = l.bound.bottom[i] + r.bound.bottom[i].r();
		tmp.bound.bottom[i][1] = l.bound.bottom[i] + r.bound.bottom[i].theta();
		tmp.bound.bottom[i][2] = l.bound.bottom[i] + r.bound.bottom[i].z();
		tmp.bound.top[i][0] = l.bound.top[i] + r.bound.top[i].r();
		tmp.bound.top[i][1] = l.bound.top[i] + r.bound.top[i].theta();
		tmp.bound.top[i][2] = l.bound.top[i] + r.bound.top[i].z();
	}

	return tmp;
}
// ------------------------------------------------------------------------- //
vector_field operator+(const vector_field& l, const scalar_field& r)
{
	vector_field tmp(l.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
	{
		tmp.in[i][0] = l.in[i].r() + r.in[i];
		tmp.in[i][1] = l.in[i].theta() + r.in[i];
		tmp.in[i][2] = l.in[i].z() + r.in[i];
	}
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i][0] = l.bound.inner[i].r() + r.bound.inner[i];
		tmp.bound.inner[i][1] = l.bound.inner[i].theta() + r.bound.inner[i];
		tmp.bound.inner[i][2] = l.bound.inner[i].z() + r.bound.inner[i];
		tmp.bound.outer[i][0] = l.bound.outer[i].r() + r.bound.outer[i];
		tmp.bound.outer[i][1] = l.bound.outer[i].theta() + r.bound.outer[i];
		tmp.bound.outer[i][2] = l.bound.outer[i].z() + r.bound.outer[i];
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i][0] = l.bound.bottom[i].r() + r.bound.bottom[i];
		tmp.bound.bottom[i][1] = l.bound.bottom[i].theta() + r.bound.bottom[i];
		tmp.bound.bottom[i][2] = l.bound.bottom[i].z() + r.bound.bottom[i];
		tmp.bound.top[i][0] = l.bound.top[i].r() + r.bound.top[i];
		tmp.bound.top[i][1] = l.bound.top[i].theta() + r.bound.top[i];
		tmp.bound.top[i][2] = l.bound.top[i].z() + r.bound.top[i];
	}

	return tmp;
}
// ---------------------------- operator - --------------------------------- //
vector_field operator-(const vector_field& l, const vector_field& r)
{
	vector_field tmp(l.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
	{
		tmp.in[i][0] = l.in[i].r() - r.in[i].r();
		tmp.in[i][1] = l.in[i].theta() - r.in[i].theta();
		tmp.in[i][2] = l.in[i].z() - r.in[i].z();
	}
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i][0] = l.bound.inner[i].r() - r.bound.inner[i].r();
		tmp.bound.inner[i][1] = l.bound.inner[i].theta() - r.bound.inner[i].theta();
		tmp.bound.inner[i][2] = l.bound.inner[i].z() - r.bound.inner[i].z();
		tmp.bound.outer[i][0] = l.bound.outer[i].r() - r.bound.outer[i].r();
		tmp.bound.outer[i][1] = l.bound.outer[i].theta() - r.bound.outer[i].theta();
		tmp.bound.outer[i][2] = l.bound.outer[i].z() - r.bound.outer[i].z();
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i][0] = l.bound.bottom[i].r() - r.bound.bottom[i].r();
		tmp.bound.bottom[i][1] = l.bound.bottom[i].theta() - r.bound.bottom[i].theta();
		tmp.bound.bottom[i][2] = l.bound.bottom[i].z() - r.bound.bottom[i].z();
		tmp.bound.top[i][0] = l.bound.top[i].r() - r.bound.top[i].r();
		tmp.bound.top[i][1] = l.bound.top[i].theta() - r.bound.top[i].theta();
		tmp.bound.top[i][2] = l.bound.top[i].z() - r.bound.top[i].z();
	}

	return tmp;
}
// ------------------------------------------------------------------------- //
vector_field operator-(const scalar_field& l, const vector_field& r)
{
	vector_field tmp(r.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
	{
		tmp.in[i][0] = l.in[i] - r.in[i].r();
		tmp.in[i][1] = l.in[i] - r.in[i].theta();
		tmp.in[i][2] = l.in[i] - r.in[i].z();
	}
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i][0] = l.bound.inner[i] - r.bound.inner[i].r();
		tmp.bound.inner[i][1] = l.bound.inner[i] - r.bound.inner[i].theta();
		tmp.bound.inner[i][2] = l.bound.inner[i] - r.bound.inner[i].z();
		tmp.bound.outer[i][0] = l.bound.outer[i] - r.bound.outer[i].r();
		tmp.bound.outer[i][1] = l.bound.outer[i] - r.bound.outer[i].theta();
		tmp.bound.outer[i][2] = l.bound.outer[i] - r.bound.outer[i].z();
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i][0] = l.bound.bottom[i] - r.bound.bottom[i].r();
		tmp.bound.bottom[i][1] = l.bound.bottom[i] - r.bound.bottom[i].theta();
		tmp.bound.bottom[i][2] = l.bound.bottom[i] - r.bound.bottom[i].z();
		tmp.bound.top[i][0] = l.bound.top[i] - r.bound.top[i].r();
		tmp.bound.top[i][1] = l.bound.top[i] - r.bound.top[i].theta();
		tmp.bound.top[i][2] = l.bound.top[i] - r.bound.top[i].z();
	}

	return tmp;
}
// ------------------------------------------------------------------------- //
vector_field operator-(const vector_field& l, const scalar_field& r)
{
	vector_field tmp(l.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
	{
		tmp.in[i][0] = l.in[i].r() - r.in[i];
		tmp.in[i][1] = l.in[i].theta() - r.in[i];
		tmp.in[i][2] = l.in[i].z() - r.in[i];
	}
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i][0] = l.bound.inner[i].r() - r.bound.inner[i];
		tmp.bound.inner[i][1] = l.bound.inner[i].theta() - r.bound.inner[i];
		tmp.bound.inner[i][2] = l.bound.inner[i].z() - r.bound.inner[i];
		tmp.bound.outer[i][0] = l.bound.outer[i].r() - r.bound.outer[i];
		tmp.bound.outer[i][1] = l.bound.outer[i].theta() - r.bound.outer[i];
		tmp.bound.outer[i][2] = l.bound.outer[i].z() - r.bound.outer[i];
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i][0] = l.bound.bottom[i].r() - r.bound.bottom[i];
		tmp.bound.bottom[i][1] = l.bound.bottom[i].theta() - r.bound.bottom[i];
		tmp.bound.bottom[i][2] = l.bound.bottom[i].z() - r.bound.bottom[i];
		tmp.bound.top[i][0] = l.bound.top[i].r() - r.bound.top[i];
		tmp.bound.top[i][1] = l.bound.top[i].theta() - r.bound.top[i];
		tmp.bound.top[i][2] = l.bound.top[i].z() - r.bound.top[i];
	}

	return tmp;
}
// ---------------------------- operator * --------------------------------- //
vector_field operator*(const vector_field& l, const vector_field& r)
{
	vector_field tmp(l.mesh);
	//// internal field
	for(register int i=0; i<tmp.mnq; i++)
	{
		tmp.in[i][0] = l.in[i].r() * r.in[i].r();
		tmp.in[i][1] = l.in[i].theta() * r.in[i].theta();
		tmp.in[i][2] = l.in[i].z() * r.in[i].z();
	}
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i][0] = l.bound.inner[i].r() * r.bound.inner[i].r();
		tmp.bound.inner[i][1] = l.bound.inner[i].theta() * r.bound.inner[i].theta();
		tmp.bound.inner[i][2] = l.bound.inner[i].z() * r.bound.inner[i].z();
		tmp.bound.outer[i][0] = l.bound.outer[i].r() * r.bound.outer[i].r();
		tmp.bound.outer[i][1] = l.bound.outer[i].theta() * r.bound.outer[i].theta();
		tmp.bound.outer[i][2] = l.bound.outer[i].z() * r.bound.outer[i].z();
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i][0] = l.bound.bottom[i].r() * r.bound.bottom[i].r();
		tmp.bound.bottom[i][1] = l.bound.bottom[i].theta() * r.bound.bottom[i].theta();
		tmp.bound.bottom[i][2] = l.bound.bottom[i].z() * r.bound.bottom[i].z();
		tmp.bound.top[i][0] = l.bound.top[i].r() * r.bound.top[i].r();
		tmp.bound.top[i][1] = l.bound.top[i].theta() * r.bound.top[i].theta();
		tmp.bound.top[i][2] = l.bound.top[i].z() * r.bound.top[i].z();
	}

	return tmp;
}
// ------------------------------------------------------------------------- //
vector_field operator*(const scalar_field& l, const vector_field& r)
{
	vector_field tmp(r.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
	{
		tmp.in[i][0] = l.in[i] * r.in[i].r();
		tmp.in[i][1] = l.in[i] * r.in[i].theta();
		tmp.in[i][2] = l.in[i] * r.in[i].z();
	}
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i][0] = l.bound.inner[i] * r.bound.inner[i].r();
		tmp.bound.inner[i][1] = l.bound.inner[i] * r.bound.inner[i].theta();
		tmp.bound.inner[i][2] = l.bound.inner[i] * r.bound.inner[i].z();
		tmp.bound.outer[i][0] = l.bound.outer[i] * r.bound.outer[i].r();
		tmp.bound.outer[i][1] = l.bound.outer[i] * r.bound.outer[i].theta();
		tmp.bound.outer[i][2] = l.bound.outer[i] * r.bound.outer[i].z();
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i][0] = l.bound.bottom[i] * r.bound.bottom[i].r();
		tmp.bound.bottom[i][1] = l.bound.bottom[i] * r.bound.bottom[i].theta();
		tmp.bound.bottom[i][2] = l.bound.bottom[i] * r.bound.bottom[i].z();
		tmp.bound.top[i][0] = l.bound.top[i] * r.bound.top[i].r();
		tmp.bound.top[i][1] = l.bound.top[i] * r.bound.top[i].theta();
		tmp.bound.top[i][2] = l.bound.top[i] * r.bound.top[i].z();
	}

	return tmp;
}
// ------------------------------------------------------------------------- //
vector_field operator*(const vector_field& l, const scalar_field& r)
{
	vector_field tmp(l.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
	{
		tmp.in[i][0] = l.in[i].r() * r.in[i];
		tmp.in[i][1] = l.in[i].theta() * r.in[i];
		tmp.in[i][2] = l.in[i].z() * r.in[i];
	}
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i][0] = l.bound.inner[i].r() * r.bound.inner[i];
		tmp.bound.inner[i][1] = l.bound.inner[i].theta() * r.bound.inner[i];
		tmp.bound.inner[i][2] = l.bound.inner[i].z() * r.bound.inner[i];
		tmp.bound.outer[i][0] = l.bound.outer[i].r() * r.bound.outer[i];
		tmp.bound.outer[i][1] = l.bound.outer[i].theta() * r.bound.outer[i];
		tmp.bound.outer[i][2] = l.bound.outer[i].z() * r.bound.outer[i];
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i][0] = l.bound.bottom[i].r() * r.bound.bottom[i];
		tmp.bound.bottom[i][1] = l.bound.bottom[i].theta() * r.bound.bottom[i];
		tmp.bound.bottom[i][2] = l.bound.bottom[i].z() * r.bound.bottom[i];
		tmp.bound.top[i][0] = l.bound.top[i].r() * r.bound.top[i];
		tmp.bound.top[i][1] = l.bound.top[i].theta() * r.bound.top[i];
		tmp.bound.top[i][2] = l.bound.top[i].z() * r.bound.top[i];
	}

	return tmp;
}
// ---------------------------- operator / --------------------------------- //
vector_field operator/(const vector_field& l, const vector_field& r)
{
	vector_field tmp(l.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
	{
		tmp.in[i][0] = l.in[i].r() / r.in[i].r();
		tmp.in[i][1] = l.in[i].theta() / r.in[i].theta();
		tmp.in[i][2] = l.in[i].z() / r.in[i].z();
	}
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i][0] = l.bound.inner[i].r() / r.bound.inner[i].r();
		tmp.bound.inner[i][1] = l.bound.inner[i].theta() / r.bound.inner[i].theta();
		tmp.bound.inner[i][2] = l.bound.inner[i].z() / r.bound.inner[i].z();
		tmp.bound.outer[i][0] = l.bound.outer[i].r() / r.bound.outer[i].r();
		tmp.bound.outer[i][1] = l.bound.outer[i].theta() / r.bound.outer[i].theta();
		tmp.bound.outer[i][2] = l.bound.outer[i].z() / r.bound.outer[i].z();
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i][0] = l.bound.bottom[i].r() / r.bound.bottom[i].r();
		tmp.bound.bottom[i][1] = l.bound.bottom[i].theta() / r.bound.bottom[i].theta();
		tmp.bound.bottom[i][2] = l.bound.bottom[i].z() / r.bound.bottom[i].z();
		tmp.bound.top[i][0] = l.bound.top[i].r() / r.bound.top[i].r();
		tmp.bound.top[i][1] = l.bound.top[i].theta() / r.bound.top[i].theta();
		tmp.bound.top[i][2] = l.bound.top[i].z() / r.bound.top[i].z();
	}

	return tmp;
}
// ------------------------------------------------------------------------- //
vector_field operator/(const vector_field& l, const scalar_field& r)
{
	vector_field tmp(l.mesh);
	// internal field
	for(register int i=0; i<tmp.mnq; i++)
	{
		tmp.in[i][0] = l.in[i].r() / r.in[i];
		tmp.in[i][1] = l.in[i].theta() / r.in[i];
		tmp.in[i][2] = l.in[i].z() / r.in[i];
	}
	// inner and outer boundaries
	for(register int i=0; i<tmp.n; i++)
	{
		tmp.bound.inner[i][0] = l.bound.inner[i].r() / r.bound.inner[i];
		tmp.bound.inner[i][1] = l.bound.inner[i].theta() / r.bound.inner[i];
		tmp.bound.inner[i][2] = l.bound.inner[i].z() / r.bound.inner[i];
		tmp.bound.outer[i][0] = l.bound.outer[i].r() / r.bound.outer[i];
		tmp.bound.outer[i][1] = l.bound.outer[i].theta() / r.bound.outer[i];
		tmp.bound.outer[i][2] = l.bound.outer[i].z() / r.bound.outer[i];
	}
	// bottom and top boundaries
	for(register int i=0; i<tmp.mn; i++)
	{
		tmp.bound.bottom[i][0] = l.bound.bottom[i].r() / r.bound.bottom[i];
		tmp.bound.bottom[i][1] = l.bound.bottom[i].theta() / r.bound.bottom[i];
		tmp.bound.bottom[i][2] = l.bound.bottom[i].z() / r.bound.bottom[i];
		tmp.bound.top[i][0] = l.bound.top[i].r() / r.bound.top[i];
		tmp.bound.top[i][1] = l.bound.top[i].theta() / r.bound.top[i];
		tmp.bound.top[i][2] = l.bound.top[i].z() / r.bound.top[i];
	}

	return tmp;
}
// ------------------------------------------------------------------------- //
