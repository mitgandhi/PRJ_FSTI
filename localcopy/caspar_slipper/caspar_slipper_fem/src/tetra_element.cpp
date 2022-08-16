#include "tetra_element.h"

namespace CasparSlipperFEM
{

	tetra_element::tetra_element() : element(TETRA)
	{
		nodes.resize(4);
	};

	tetra_element::tetra_element(node * n0, node * n1, node * n2, node * n3) : element(TETRA)
	{
		nodes.resize(4);
		nodes[0] = n0;
		nodes[1] = n1;
		nodes[2] = n2;
		nodes[3] = n3;
	};

	double tetra_element::V()
	{
		return 1.0/6.0*detJ();
	};

	vector<int> tetra_element::l2gDOF()
	{
		vector<int> m;
		for(int n=0; n<nnodes; n++)
		{
			for(int d=0; d<nodeDOF; d++)
			{
				m.push_back(nodes[n]->DOF[d]);
			}
		}

		return m;
	}

	pair<int, int> tetra_element::niddof(const int ki)
	{
		int KI = 0;

		for(int n=0; n<nnodes; n++)
		{
			for(int d=0; d<nodeDOF; d++)
			{
				if(KI == ki)
				{
					return pair<int, int> (nodes[n]->id, d);
				}
				KI++;
			}
		}

		//invalid ki value passed
		return pair<int, int> (-1, -1);
	}

	void tetra_element::assign_analysis(const analysis_types type)
	{
		analysis_type = type;

		//Update the nodeDOF
		nodeDOF = 0;

		//Update the nodal DOF depending on type
		if(analysis_type == ELASTIC)
		{
			nodeDOF = 3;
		} else if (analysis_type == THERMAL)
		{
			nodeDOF = 1;
		}

		for(int n=0; n<nnodes; n++)
		{
			nodes[n]->DOF.clear();
			nodes[n]->DOF.resize(nodeDOF);
		}

		//Also resize the BC stores
		Kbc.resize(nnodes*nodeDOF, nnodes*nodeDOF);
		Fbc.clear();
		Fbc.resize(nnodes*nodeDOF, 0);
	}

	//Stiffness matrix construction
	dmatrix<double> tetra_element::Be()
	{
		dmatrix<double> B(6,12);

		const double a1 = y(2)*z(1) - y(1)*z(2) + y(1)*z(3) - y(3)*z(1) - y(2)*z(3) + y(3)*z(2);
		const double a2 = y(0)*z(2) - y(2)*z(0) - y(0)*z(3) + y(3)*z(0) + y(2)*z(3) - y(3)*z(2);
		const double a3 = y(1)*z(0) - y(0)*z(1) + y(0)*z(3) - y(3)*z(0) - y(1)*z(3) + y(3)*z(1);
		const double a4 = y(0)*z(1) - y(1)*z(0) - y(0)*z(2) + y(2)*z(0) + y(1)*z(2) - y(2)*z(1);
		const double b1 = x(1)*z(2) - x(2)*z(1) - x(1)*z(3) + x(3)*z(1) + x(2)*z(3) - x(3)*z(2);
		const double b2 = x(2)*z(0) - x(0)*z(2) + x(0)*z(3) - x(3)*z(0) - x(2)*z(3) + x(3)*z(2);
		const double b3 = x(0)*z(1) - x(1)*z(0) - x(0)*z(3) + x(3)*z(0) + x(1)*z(3) - x(3)*z(1);
		const double b4 = x(1)*z(0) - x(0)*z(1) + x(0)*z(2) - x(2)*z(0) - x(1)*z(2) + x(2)*z(1);
		const double c1 = x(2)*y(1) - x(1)*y(2) + x(1)*y(3) - x(3)*y(1) - x(2)*y(3) + x(3)*y(2);
		const double c2 = x(0)*y(2) - x(2)*y(0) - x(0)*y(3) + x(3)*y(0) + x(2)*y(3) - x(3)*y(2);
		const double c3 = x(1)*y(0) - x(0)*y(1) + x(0)*y(3) - x(3)*y(0) - x(1)*y(3) + x(3)*y(1);
		const double c4 = x(0)*y(1) - x(1)*y(0) - x(0)*y(2) + x(2)*y(0) + x(1)*y(2) - x(2)*y(1);

		B(0,0)  = a1;
		B(0,3)  = a2;
		B(0,6)  = a3;
		B(0,9)  = a4;

		B(1,1)  = b1;
		B(1,4)  = b2;
		B(1,7)  = b3;
		B(1,10) = b4;

		B(2,2)  = c1;
		B(2,5)  = c2;
		B(2,8)  = c3;
		B(2,11) = c4;

		B(3,0)  = b1;
		B(3,1)  = a1;
		B(3,3)  = b2;
		B(3,4)  = a2;
		B(3,6)  = b3;
		B(3,7)  = a3;
		B(3,9)  = b4;
		B(3,10) = a4;

		B(4,1)  = c1;
		B(4,2)  = b1;
		B(4,4)  = c2;
		B(4,5)  = b2;
		B(4,7)  = c3;
		B(4,8)  = b3;
		B(4,10) = c4;
		B(4,11) = b4;

		B(5,0)  = c1;
		B(5,2)  = a1;
		B(5,3)  = c2;
		B(5,5)  = a2;
		B(5,6)  = c3;
		B(5,8)  = a3;
		B(5,9)  = c4;
		B(5,11) = a4;

		B *= 1.0/(detJ());
		
		return B;
	};

	dmatrix<double> tetra_element::Bt()
	{
		dmatrix<double> B(3,4);

		const double a1 = y(2)*z(1) - y(1)*z(2) + y(1)*z(3) - y(3)*z(1) - y(2)*z(3) + y(3)*z(2);
		const double a2 = y(0)*z(2) - y(2)*z(0) - y(0)*z(3) + y(3)*z(0) + y(2)*z(3) - y(3)*z(2);
		const double a3 = y(1)*z(0) - y(0)*z(1) + y(0)*z(3) - y(3)*z(0) - y(1)*z(3) + y(3)*z(1);
		const double a4 = y(0)*z(1) - y(1)*z(0) - y(0)*z(2) + y(2)*z(0) + y(1)*z(2) - y(2)*z(1);
		const double b1 = x(1)*z(2) - x(2)*z(1) - x(1)*z(3) + x(3)*z(1) + x(2)*z(3) - x(3)*z(2);
		const double b2 = x(2)*z(0) - x(0)*z(2) + x(0)*z(3) - x(3)*z(0) - x(2)*z(3) + x(3)*z(2);
		const double b3 = x(0)*z(1) - x(1)*z(0) - x(0)*z(3) + x(3)*z(0) + x(1)*z(3) - x(3)*z(1);
		const double b4 = x(1)*z(0) - x(0)*z(1) + x(0)*z(2) - x(2)*z(0) - x(1)*z(2) + x(2)*z(1);
		const double c1 = x(2)*y(1) - x(1)*y(2) + x(1)*y(3) - x(3)*y(1) - x(2)*y(3) + x(3)*y(2);
		const double c2 = x(0)*y(2) - x(2)*y(0) - x(0)*y(3) + x(3)*y(0) + x(2)*y(3) - x(3)*y(2);
		const double c3 = x(1)*y(0) - x(0)*y(1) + x(0)*y(3) - x(3)*y(0) - x(1)*y(3) + x(3)*y(1);
		const double c4 = x(0)*y(1) - x(1)*y(0) - x(0)*y(2) + x(2)*y(0) + x(1)*y(2) - x(2)*y(1);

		B(0,0)  = a1;
		B(0,1)  = a2;
		B(0,2)  = a3;
		B(0,3)  = a4;

		B(1,0)  = b1;
		B(1,1)  = b2;
		B(1,2)  = b3;
		B(1,3)  =  b4;

		B(2,0)  = c1;
		B(2,1)  = c2;
		B(2,2)  = c3;
		B(2,3)  = c4;

		B *= 1.0/(detJ());

		return B;
	};

	dmatrix<double> tetra_element::De()
	{
		dmatrix<double> D(6,6);

		D(0,0) = D(1,1) = D(2,2) = 1.0-nu;
		D(3,3) = D(4,4) = D(5,5) = 0.5-nu;
		D(0,1) = D(0,2) = nu;
		D(1,0) = D(1,2) = nu;
		D(2,0) = D(2,1) = nu;

		D *= E/((1.0+nu)*(1.0-2.0*nu));

		return D;
	};

	dmatrix<double> tetra_element::Dt()
	{
		dmatrix<double> D(3,3);
		
		D(0,0) = D(1,1) = D(2,2) = k;
		
		return D;
	};

	dmatrix<double> tetra_element::Ke()
	{
		dmatrix<double> b = Be();
		dmatrix<double> k = b.t()*De()*b;
		k *= V();

		//add in mixed boundary values
		k += Kbc;

		return k;
	};

	dmatrix<double> tetra_element::Kt()
	{
		dmatrix<double> b = Bt();
		dmatrix<double> k = b.t()*Dt()*b;
		k *= V();

		//add in mixed boundary values
		k += Kbc;

		return k;
	};



	dmatrix<double> tetra_element::K()
	{
		if(analysis_type == ELASTIC)
		{
			return Ke();
		} else if (analysis_type == THERMAL)
		{
			return Kt();
		} 

		return dmatrix<double>(0,0);
	}

	vector<double> tetra_element::F()
	{
		vector<double> f;

		f = Fbc;

		return f;
	};

	double tetra_element::detJ()
	{
		dmatrix<double> J(4,4);
		for(int n=0; n<4; n++)
		{
			J(0,n) = 1;
			for(int i=0; i<3; i++)
			{
				J(i+1,n) = nodes[n]->coord[i];
			}
		}

		return J.det();
	}


	//BC methods
	point tetra_element::get_center()
	{
		point p(0,0,0);

		for(int n=0; n<nnodes; n++)
		{
			for(int i=0; i<3; i++)
			{
				p[i] += nodes[n]->coord[i]/double(nnodes);
			}
		}

		return p;
	};
	vector<face> tetra_element::getfaces()
	{
		vector<face> fs;
		
		//the four faces
		{
			face f;
			f.push_back(nodes[0]);
			f.push_back(nodes[1]);
			f.push_back(nodes[2]);
			fs.push_back(f);
		}

		{
			face f;
			f.push_back(nodes[0]);
			f.push_back(nodes[3]);
			f.push_back(nodes[1]);
			fs.push_back(f);
		}

		{
			face f;
			f.push_back(nodes[1]);
			f.push_back(nodes[3]);
			f.push_back(nodes[2]);
			fs.push_back(f);
		}

		{
			face f;
			f.push_back(nodes[2]);
			f.push_back(nodes[3]);
			f.push_back(nodes[0]);
			fs.push_back(f);
		}

		return fs;
	}

	void tetra_element::assign_face_nbc(const face & f, const double val)
	{
		const int nn = 3;	//number of nodes in the face

		if(f.nodes.size() != nn)
		{
			//element must have a nn node face
			return;
		}

		double load = f.area()*val / double(nn);	//distributed nodal load
		point norm = f.normal(); //find the unit normal vector
		

		for(int n=0; n<nn; n++)
		{
			if(analysis_type == THERMAL)
			{
				//find which interal node cooresponds to the node of the face
				for(int i=0; i<nnodes; i++)
				{
					if(nodes[i]->id == f.nodes[n]->id)
					{
						//assign the load
						Fbc[i] += load;

						//break the interal for loop
						break;
					}
				}
			} else if (analysis_type == ELASTIC)
			{
				//find which interal node cooresponds to the node of the face
				for(int i=0; i<4; i++)
				{
					if(nodes[i]->id == f.nodes[n]->id)
					{
						//assign the load
						Fbc[i*3+0] += load*norm[0];
						Fbc[i*3+1] += load*norm[1];
						Fbc[i*3+2] += load*norm[2];

						//break the internal for loop
						break;
					}
				}
			}
		}
		
	};
	void tetra_element::assign_face_mbc(const face & f, const double h, const double phi)
	{
		const int nn = 3;	//number of nodes in the face

		if(f.nodes.size() != nn)
		{
			//element must have a nn node face
			return;
		}

		double hA = f.area()*h;

		//map the face node to the internal node
		vector<int> lid(nn, -1);	
		for(int n=0; n<nn; n++)
		{
			for(int i=0; i<nnodes; i++)
			{
				if(nodes[i]->id == f.nodes[n]->id)
				{
					//update the map
					lid[n] = i;

					//break the interal for loop
					break;
				}
			}
		}

		//assign the mbc
		for(int i=0; i<nn; i++)
		{
			if(analysis_type == THERMAL)
			{
				//update the F vector
				Fbc[lid[i]] += hA*phi/double(nn);
			
				//update the k matrix
				for(int j=0; j<nn; j++)
				{
					if(i == j)
					{
						Kbc(lid[i], lid[j]) += hA/6.0;
					} else {
						Kbc(lid[i], lid[j]) += hA/12.0;
					}
				}
			
			} else if (analysis_type == ELASTIC)
			{
				//Mixed face BC for ELASTIC analysis is not implemented
			}
		}
		
	};


	void tetra_element::assign_ele_therm(const double dt)
	{
		if(analysis_type == THERMAL)
		{
			//can't apply a thermal deformation load to a thermal analysis
			return;
		}

		//formulate the thermal strain
		dmatrix<double> eps(6,1);
		eps(0,0) = dt*alpha;
		eps(1,0) = dt*alpha;
		eps(2,0) = dt*alpha;

		dmatrix<double> Ft = Be().t()*De()*eps;
		Ft *= V();

		//push the dmatrix into the Fbc vector
		for(int i=0; i<nnodes*nodeDOF; i++)
		{
			Fbc[i] += Ft(i,0);
		}	
	}
	//Utility functions used to compactly access the nodal x or y or z
	__inline double tetra_element::x(const int n)
	{
		return nodes[n]->coord[0];
	};

	__inline double tetra_element::y(const int n)
	{
		return nodes[n]->coord[1];
	};

	__inline double tetra_element::z(const int n)
	{
		return nodes[n]->coord[2];
	};

};
