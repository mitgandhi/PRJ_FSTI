#include "CSlipperGap.h"
#include <omp.h>
#include <time.h>
#include <map>

#pragma once
#define OMP

extern void matlab(const Array<double,2>& data,const string file);

struct cell
{
	int nid[8];
	//int cid[1][1][1];

	cell()
	{
		for(int i=0; i<8; i++)
		{
			nid[i] = -1;
		}
	}

};

struct node
{
	double coord[3];
	int hitcnt;

	node()
	{
		hitcnt = 0;
		coord[0] = 0;
		coord[1] = 0;
		coord[2] = 0;
	}
};

struct octant
{
	int lat;	//S: -1, N: 1
	int lon; //W: -1, E: 1
	int alt; //B: -1, T: 1

	static const int LA = 100;
	static const int LO = 10;
	static const int AL = 1;

	octant()
	{
	};

	octant(int latitude, int longitude, int altitude) : lat(latitude), lon(longitude), alt(altitude)
	{
	};

	octant flip(int LaLoAl)
	{
		//create a copy of this octant
		octant tmp(lat,lon,alt);

		//Should be in the form:
		int Al = LaLoAl % 2;
		LaLoAl = (LaLoAl-Al)/10;
		int Lo = LaLoAl % 2;
		LaLoAl = (LaLoAl-Lo)/10;
		int La = LaLoAl % 2;

		//Flip whichever mask bits are set
		if(Al == 1)
		{
			tmp.alt *= -1;
		}
		if(Lo == 1)
		{
			tmp.lon *= -1;
		}
		if(La == 1)
		{
			tmp.lat *= -1;
		}

		return tmp;
	}

	octant(int n)
	{
		if(n == 0 || n == 1 || n == 4 || n == 5)
		{
			lat = -1;
		} else {
			lat = 1;
		}
		if(n == 0 || n == 3 || n == 4 || n == 7)
		{
			lon = -1;
		} else {
			lon = 1;
		}
		if(n == 0 || n == 1 || n == 2 || n == 3)
		{
			alt = -1;
		} else {
			alt = 1;
		}
	}

	int n()
	{
		
		const int S = -1;
		const int N = 1;
		const int W = -1;
		const int E = 1;
		const int B = -1;
		const int T = 1;

		if(alt == B && lat == S && lon == W)
		{
			return 0;
		}
		if(alt == B && lat == S && lon == E)
		{
			return 1;
		}
		if(alt == B && lat == N && lon == E)
		{
			return 2;
		}
		if(alt == B && lat == N && lon == W)
		{
			return 3;
		}
		if(alt == T && lat == S && lon == W)
		{
			return 4;
		}
		if(alt == T && lat == S && lon == E)
		{
			return 5;
		}
		if(alt == T && lat == N && lon == E)
		{
			return 6;
		}
		if(alt == T && lat == N && lon == W)
		{
			return 7;
		}

		//Something is wrong:
		return -1;
	}

};


/*
struct octant
{
	//ALTHOUGH THIS CLASS IS MUCH SMALLER THAN THE PREVIOUS ONE
	//THE ORDER OF THE NODES IS WRONG AND NEEDS TO SOMEHOW BE CHANGED

	static const char LA = 1;
	static const char LO = 4;
	static const char AL = 2;

	char c;

	octant(const int n) : c(n)
	{
	}

	//extract a specific bit and then convert to -1, 1
	int lat(void)
	{
		return 2*bool(c&LA) - 1;
	}
	int lon(void)
	{
		return 2*bool(c&LO) - 1;
	}
	int alt(void)
	{
		return 2*bool(c&AL) - 1;
	}

	char flip(const char f = 0)
	{
		return c^f;
	}

	//MOVE THIS TO THE VTK CLASS
	void setoctant(int m, int n, int q, int o, int nid)
	{
		//The n domain is connected
		n = (n+Fluid->N)%Fluid->N;

		if(m < 0 || m >= Fluid->M || n < 0 || n >= Fluid->N || q < 0 || q >= Fluid->Q)
		{
			//break if it is a boundary
			return;
		}

		mesh(m,n,q).nid[o] = nid;
	}
};
*/

class vtk
{
public:
	CSlipperGap * SlipperGap;

	caspar_input * gapinput;
	CSlipperGap::sFluid * Fluid;
	CSlipperGap::sSolid * slipper;
	CSlipperGap::sSolid * swashplate;
	CSlipperGap::sOperatingSlipperGap * operatingslippergap; 
	CSlipperGap::sGeometrySlipperGap * geometryslippergap; 
	
	Array<double,2> pnosqz;

	Array<cell,3> mesh;
	vector<node> nodes;

	int Q;

	vtk(CSlipperGap * S) : SlipperGap(S)
	{
		gapinput = SlipperGap->gapinput;
		Fluid = &(SlipperGap->Fluid);
		slipper = &(SlipperGap->slipper);
		swashplate = &(SlipperGap->swashplate);
		operatingslippergap = &(SlipperGap->operatingslippergap);
		geometryslippergap = &(SlipperGap->geometryslippergap);
	}
	
	void setoctant(int m, int n, int q, octant o, int nid, bool cheap = false)
	{
		//The n domain is connected
		n = (n+Fluid->N)%Fluid->N;

		int qmax;
		if(cheap)
		{
			qmax = 1;
		} else {
			qmax = Fluid->Q;
		}

		if(m < 0 || m >= Fluid->M || n < 0 || n >= Fluid->N || q < 0 || q >= qmax)
		{
			//break if it is a boundary
			return;
		}

		mesh(m,n,q).nid[o.n()] = nid;
	}
			

	void build(const bool GCS = false)
	{
		mesh.resize(Fluid->M, Fluid->N, Fluid->Q);
		Q = Fluid->Q;
		int nid=0;

		//Loop through all the fluid cells
		for(int m=0; m<Fluid->M; m++)
		{
			for(int n=0; n<Fluid->N; n++)
			{
				for(int q=0; q<Fluid->Q; q++)
				{
					//For each fluid cell, check each octant
					for(int nd=0; nd<8; nd++)
					{
						//Check if the octant has a node
						if(mesh(m,n,q).nid[nd] == -1)
						{
							//Since the octant doesn't have a node, 
							//we need to assign a node and update the fluid cells
							//which share this octant

							octant o(nd);	//the current octant

							//to find each neighbor cell, we add every combination of lat, lon, and alt of
							//the octant o, to our current cell
							//That neighbor cell shares the 'flipped' octant of our current one, o

							//setoctant() just skips out of bound cells
							setoctant(m, n, q, o, nid);
							setoctant(m+o.lat, n, q, o.flip(o.LA), nid);
							setoctant(m+o.lat, n+o.lon, q, o.flip(o.LA+o.LO), nid);
							setoctant(m+o.lat, n, q+o.alt, o.flip(o.LA+o.AL), nid);
							setoctant(m, n+o.lon, q, o.flip(o.LO), nid);
							setoctant(m, n+o.lon, q+o.alt, o.flip(o.LO+o.AL), nid);
							setoctant(m, n, q+o.alt, o.flip(o.AL), nid);
							setoctant(m+o.lat, n+o.lon, q+o.alt, o.flip(o.LA+o.LO+o.AL), nid);

							nid++;
						}
					}
				}
			}
		}

		//resize the nodes vector
		nodes.resize(nid);

		//Loop through all fluid cells
		for(int m=0; m<Fluid->M; m++)
		{
			for(int n=0; n<Fluid->N; n++)
			{
				for(int q=0; q<Fluid->Q; q++)
				{
					for(int nd=0; nd<8; nd++)
					{

						octant o(nd);

						double r;

						//Use 'half' control volumes at the raidal boundries
						if(m+o.lat < 0 || m+o.lat >= Fluid->M)
						{
							r = Fluid->r(m,n);
						} else {
							//Normal r						
							r = Fluid->r(m,n) + o.lat * 0.5 * Fluid->dr(m,n);
						}
						double t = Fluid->theta(m,n) + o.lon * 0.5 * Fluid->dtheta(m,n);
						
						//Create the LCS
						double Rx = r*sin(t);
						double Ry = r*cos(t);

						if(GCS)
						{
							const double phi = operatingslippergap->phi_rad;

							//Move LCS -> GCS
							//Rotate the LCS about its origin to the GCS alignment
							double GRx = Rx*cos(-phi) - Ry*sin(-phi);
							double GRy = Rx*sin(-phi) + Ry*cos(-phi);
							
							//Translate the LCS origin to the GCS origin
							Rx = GRx + geometryslippergap->rB*sin(phi);
							Ry = GRy + geometryslippergap->rB*cos(phi)/cos(operatingslippergap->beta_rad);
						} 
						
						nodes[mesh(m,n,q).nid[nd]].coord[0] = Rx;
						nodes[mesh(m,n,q).nid[nd]].coord[1] = Ry;
						nodes[mesh(m,n,q).nid[nd]].coord[2] += Fluid->z(m,n,q);
						nodes[mesh(m,n,q).nid[nd]].hitcnt++;
					}
				}
			}
		}

		//Average the nodes
		for(int n=0; n<nodes.size(); n++)
		{
			for(int i=2; i<3; i++)
			{
				nodes[n].coord[i] /= (double) nodes[n].hitcnt;
			}

			//Reset the hitcnt
			nodes[n].hitcnt = 0;
		}
	}

	void buildCheap(const bool GCS = false)
	{
		mesh.resize(Fluid->M, Fluid->N, 1);
		Q = 1;
		int nid=0;

		//Loop through all the fluid cells
		for(int m=0; m<Fluid->M; m++)
		{
			for(int n=0; n<Fluid->N; n++)
			{
				for(int q=0; q<1; q++)
				{
					//For each fluid cell, check each octant
					for(int nd=0; nd<8; nd++)
					{
						//Check if the octant has a node
						if(mesh(m,n,q).nid[nd] == -1)
						{
							//Since the octant doesn't have a node, 
							//we need to assign a node and update the fluid cells
							//which share this octant

							octant o(nd);	//the current octant

							//to find each neighbor cell, we add every combination of lat, lon, and alt of
							//the octant o, to our current cell
							//That neighbor cell shares the 'flipped' octant of our current one, o

							//setoctant() just skips out of bound cells
							setoctant(m, n, q, o, nid, true);
							setoctant(m+o.lat, n, q, o.flip(o.LA), nid, true);
							setoctant(m+o.lat, n+o.lon, q, o.flip(o.LA+o.LO), nid, true);
							setoctant(m+o.lat, n, q+o.alt, o.flip(o.LA+o.AL), nid, true);
							setoctant(m, n+o.lon, q, o.flip(o.LO), nid, true);
							setoctant(m, n+o.lon, q+o.alt, o.flip(o.LO+o.AL), nid, true);
							setoctant(m, n, q+o.alt, o.flip(o.AL), nid, true);
							setoctant(m+o.lat, n+o.lon, q+o.alt, o.flip(o.LA+o.LO+o.AL), nid, true);

							nid++;
						}
					}
				}
			}
		}

		//resize the nodes vector
		nodes.resize(nid);

		//Loop through all fluid cells
		for(int m=0; m<Fluid->M; m++)
		{
			for(int n=0; n<Fluid->N; n++)
			{
				for(int q=0; q<1; q++)
				{
					for(int nd=0; nd<8; nd++)
					{

						octant o(nd);
						
						double r;

						//Use 'half' control volumes at the raidal boundries
						if(m+o.lat < 0 || m+o.lat >= Fluid->M)
						{
							r = Fluid->r(m,n);
						} else {
							//Normal r						
							r = Fluid->r(m,n) + o.lat * 0.5 * Fluid->dr(m,n);
						}
						double t = Fluid->theta(m,n) + o.lon * 0.5 * Fluid->dtheta(m,n);
						
						//Create the LCS
						double Rx = r*sin(t);
						double Ry = r*cos(t);

						if(GCS)
						{
							const double phi = operatingslippergap->phi_rad;

							//Move LCS -> GCS
							//Rotate the LCS about its origin to the GCS alignment
							double GRx = Rx*cos(-phi) - Ry*sin(-phi);
							double GRy = Rx*sin(-phi) + Ry*cos(-phi);
							
							//Translate the LCS origin to the GCS origin
							Rx = GRx + geometryslippergap->rB*sin(phi);
							Ry = GRy + geometryslippergap->rB*cos(phi)/cos(operatingslippergap->beta_rad);
						} 
						
						if(nodes[mesh(m,n,q).nid[nd]].coord[0] == 0)
						{
							nodes[mesh(m,n,q).nid[nd]].coord[0] = Rx;
							nodes[mesh(m,n,q).nid[nd]].coord[1] = Ry;
						}
						
						//This should be averaged since they are face nodes, not centroidal
						double grv = 0;
						if(Fluid->hgroove(m,n) == 1)
						{
						//	grv = 10e-6;
						}
						nodes[mesh(m,n,q).nid[nd]].coord[2] += 0.5*(Fluid->h(m,n)+grv) + o.alt * 0.5 * (Fluid->h(m,n)+grv) + swashplate->ehd(m,n);
						nodes[mesh(m,n,q).nid[nd]].hitcnt++;
					}
				}
			}
		}

		//Average the nodes
		for(int n=0; n<nodes.size(); n++)
		{
			for(int i=2; i<3; i++)
			{
				nodes[n].coord[i] /= (double) nodes[n].hitcnt;
			}

			//Reset the hitcnt
			nodes[n].hitcnt = 0;
		}

	}


	//Utility functions
	string d2h(const double d)
	{
		union _DecHex
		{
			double d;
			unsigned long int h[2];
		} DH;	

		DH.d = d;
		ostringstream o("");
		o << hex << DH.h[0] << "\t" << DH.h[1];
		return o.str();
	}

	double h2d(const string s1, const string s2)
	{
		union _DecHex
		{
			double d;
			unsigned long int h[2];
		} DH;	

		DH.h[0] = strtoul(s1.c_str(), NULL, 16);
		DH.h[1] = strtoul(s2.c_str(), NULL, 16);
		
		return DH.d;
	}


	void writevtk(string fileName, double phi = 0, double rB = 0)
	{
		ofstream f(fileName.c_str(), ios::trunc);


		f << "# vtk DataFile Version 2.0" << '\n';
		f << "#" << d2h(phi) << "\t" << d2h(rB) << "\t" << d2h(operatingslippergap->beta_rad) << '\n';
		f << "ASCII" << '\n';
		f << "DATASET UNSTRUCTURED_GRID" << '\n';
		f << "POINTS " << nodes.size() << " double" << '\n';
		
		for(int n=0; n<nodes.size(); n++)
		{
			f << setprecision(10);
			for(int i=0; i<3; i++)
			{
				f << scientific << nodes[n].coord[i] << "\t";
			}
			f << '\n';
		}

		f << "CELLS " << mesh.size() << " " << mesh.size()*9 << '\n';
		for(int m=0; m<Fluid->M; m++)
		{
			for(int n=0; n<Fluid->N; n++)
			{
				for(int q=0; q<Q; q++)
				{
					f << "8\t";
					for(int i=0; i<8; i++)
					{
						f << mesh(m,n,q).nid[i] << "\t";
					}
					f << '\n';
				}
			}
		}

		f << "CELL_TYPES " << mesh.size() << '\n';
		for(int i=0; i<mesh.size(); i++)
		{
			f << "12" << '\n';
		}

		//Cell data		
		f << "CELL_DATA " << mesh.size() << '\n';
		
		//Write cell height as data
		f << "SCALARS height_m float 1" << '\n';
		f << "LOOKUP_TABLE hei" << '\n';
		for(int m=0; m<Fluid->M; m++)
		{
			for(int n=0; n<Fluid->N; n++)
			{
				for(int q=0; q<Q; q++)
				{
					f << Fluid->h(m,n) << '\n';
				}
			}
		}
		
		//Write cell pressure as data
		f << "SCALARS pressure_Pa float 1" << '\n';
		f << "LOOKUP_TABLE pre" << '\n';
		for(int m=0; m<Fluid->M; m++)
		{
			for(int n=0; n<Fluid->N; n++)
			{
				for(int q=0; q<Q; q++)
				{
					f << Fluid->p(m,n) << '\n';
				}
			}
		}
		
		/*
		for(int vn=0; vn<=10; vn++)
		{
			//Write reynolds source terms as data
			f << "SCALARS ReyVar" + n2s(vn) + " float 1" << '\n';
			f << "LOOKUP_TABLE reyvar" << '\n';
			for(int m=0; m<Fluid->M; m++)
			{
				for(int n=0; n<Fluid->N; n++)
				{
					for(int q=0; q<Q; q++)
					{
						f << Fluid->ReyVals(vn,m,n) << '\n';
					}
				}
			}
		}
		*/

		//Write cell pFullLoop as data
		f << "SCALARS pressureFullLoop_Pa float 1" << '\n';
		f << "LOOKUP_TABLE pFL" << '\n';
		for(int m=0; m<Fluid->M; m++)
		{
			for(int n=0; n<Fluid->N; n++)
			{
				for(int q=0; q<Q; q++)
				{
					f << Fluid->pFullLoop(m,n)  << '\n';
				}
			}
		}

		//Write cell dhdt as data
		f << "SCALARS squeezeVelocity_m/s float 1" << '\n';
		f << "LOOKUP_TABLE dht" << '\n';
		for(int m=0; m<Fluid->M; m++)
		{
			for(int n=0; n<Fluid->N; n++)
			{
				for(int q=0; q<Q; q++)
				{
					f << Fluid->dht(m,n)  << '\n';
				}
			}
		}
		
		//Write cell boundary as data
		f << "SCALARS boundaryType float 1" << '\n';
		f << "LOOKUP_TABLE bnd" << '\n';
		for(int m=0; m<Fluid->M; m++)
		{
			for(int n=0; n<Fluid->N; n++)
			{
				for(int q=0; q<Q; q++)
				{
					f << Fluid->boundary(m,n) << '\n';
				}
			}
		}

		/*
		//Reynolds number
		{
			//We actually have to create this value here
			Array<double,3> vrx(Fluid->M,Fluid->N,Fluid->Q);
			Array<double,3> vry(Fluid->M,Fluid->N,Fluid->Q);
			vrx = Fluid->vr * sin(Fluid->theta(tensor::i,tensor::j)) + Fluid->vtheta * cos(Fluid->theta(tensor::i,tensor::j));
			vry = Fluid->vr * cos(Fluid->theta(tensor::i,tensor::j)) - Fluid->vtheta * sin(Fluid->theta(tensor::i,tensor::j));

			Array<double,3> V(sqrt(vrx*vrx+vry*vry));
			Array<double,2> mu(mean(Fluid->oilviscosity,tensor::k));
			Array<double,2> rho(mean(Fluid->oildensity,tensor::k));
			
			f << "SCALARS Re float 1" << '\n';
			f << "LOOKUP_TABLE re" << '\n';
			for(int m=0; m<Fluid->M; m++)
			{
				for(int n=0; n<Fluid->N; n++)
				{
					for(int q=0; q<Q; q++)
					{
						f << (rho(m,n)*mean(V(m,n,Range::all()))*Fluid->h(m,n)/mu(m,n)) << '\n';
					}
				}
			}
		}
		*/

		//Write cell contact pressure as data
		f << "SCALARS contactPressure_Pa float 1" << '\n';
		f << "LOOKUP_TABLE pct" << '\n';
		for(int m=0; m<Fluid->M; m++)
		{
			for(int n=0; n<Fluid->N; n++)
			{
				for(int q=0; q<Q; q++)
				{
					f << Fluid->p_contact(m,n) << '\n';
				}
			}
		}		
		
		//Write cell ehd as data
		f << "SCALARS slipperPehd_m float 1" << '\n';
		f << "LOOKUP_TABLE slip" << '\n';
		for(int m=0; m<Fluid->M; m++)
		{
			for(int n=0; n<Fluid->N; n++)
			{
				for(int q=0; q<Q; q++)
				{
					f << slipper->ehd(m,n) << '\n';
				}
			}
		}

		f << "SCALARS swashplatePehd_m float 1" << '\n';
		f << "LOOKUP_TABLE swash" << '\n';
		for(int m=0; m<Fluid->M; m++)
		{
			for(int n=0; n<Fluid->N; n++)
			{
				for(int q=0; q<Q; q++)
				{
					f << swashplate->ehd(m,n) << '\n';
				}
			}
		}

		if(gapinput->options_slipper.general.SlipperThermoElastic)
		{
			f << "SCALARS slipperTehd_m float 1" << '\n';
			f << "LOOKUP_TABLE slipperT" << '\n';
			for(int m=0; m<Fluid->M; m++)
			{
				for(int n=0; n<Fluid->N; n++)
				{
					for(int q=0; q<Q; q++)
					{
						f << SlipperGap->t_slipper.pdeform(m,n) << '\n';
					}
				}
			}
		}

		if(gapinput->options_slipper.general.SwashplateThermoElastic)
		{
			f << "SCALARS swashplateTehd_m float 1" << '\n';
			f << "LOOKUP_TABLE slipperT" << '\n';
			for(int m=0; m<Fluid->M; m++)
			{
				for(int n=0; n<Fluid->N; n++)
				{
					for(int q=0; q<Q; q++)
					{
						f << SlipperGap->t_swashplate.pdeform(m,n) << '\n';
					}
				}
			}
		}

		f << "SCALARS Qflux_W/m^2 float 1" << '\n';
		f << "LOOKUP_TABLE qflux" << '\n';
		for(int m=0; m<Fluid->M; m++)
		{
			for(int n=0; n<Fluid->N; n++)
			{
				for(int q=0; q<Q; q++)
				{
					f << Fluid->Qflux(m,n) << '\n';
				}
			}
		}

		//Only write the fields for the full vtk, not the cheap version
		if(Q == Fluid->Q)
		{
			//Write cell temperature as data
			f << "SCALARS temperature_C float 1" << '\n';
			f << "LOOKUP_TABLE tem" << '\n';
			for(int m=0; m<Fluid->M; m++)
			{
				for(int n=0; n<Fluid->N; n++)
				{
					for(int q=0; q<Fluid->Q; q++)
					{
						f << Fluid->T(m,n,q) << '\n';
					}
				}
			}

			//Write cell visco as data
			f << "SCALARS viscosity_Pas float 1" << '\n';
			f << "LOOKUP_TABLE mu" << '\n';
			for(int m=0; m<Fluid->M; m++)
			{
				for(int n=0; n<Fluid->N; n++)
				{
					for(int q=0; q<Fluid->Q; q++)
					{
						f << Fluid->oilviscosity(m,n,q) << '\n';
					}
				}
			}

			//Write fluid velocity as data
			f << "VECTORS fluidVelocity_m/s float" << '\n';
			//f << "LOOKUP_TABLE tem" << '\n';
			{
				//We actually have to create this value here
				Array<double,3> vrx(Fluid->M,Fluid->N,Fluid->Q);
				Array<double,3> vry(Fluid->M,Fluid->N,Fluid->Q);
				vrx = Fluid->vr * sin(Fluid->theta(tensor::i,tensor::j)) + Fluid->vtheta * cos(Fluid->theta(tensor::i,tensor::j));
				vry = Fluid->vr * cos(Fluid->theta(tensor::i,tensor::j)) - Fluid->vtheta * sin(Fluid->theta(tensor::i,tensor::j));
				
				for(int m=0; m<Fluid->M; m++)
				{
					for(int n=0; n<Fluid->N; n++)
					{
						for(int q=0; q<Fluid->Q; q++)
						{
							f	<< vrx(m,n,q) << "\t"
								<< vry(m,n,q) << "\t"
								<< "0.0" << '\n';
						}
					}
				}
			}
		} else {
			//Write 2d mean temp data
			f << "SCALARS temperature_C float 1" << '\n';
			f << "LOOKUP_TABLE tem" << '\n';
			for(int m=0; m<Fluid->M; m++)
			{
				for(int n=0; n<Fluid->N; n++)
				{
					f << mean(Fluid->T(m,n,Range::all())) << '\n';
				}
			}

			{
				//We actually have to create this value here
				Array<double,3> vrx(Fluid->M,Fluid->N,Fluid->Q);
				Array<double,3> vry(Fluid->M,Fluid->N,Fluid->Q);
				vrx = Fluid->vr * sin(Fluid->theta(tensor::i,tensor::j)) + Fluid->vtheta * cos(Fluid->theta(tensor::i,tensor::j));
				vry = Fluid->vr * cos(Fluid->theta(tensor::i,tensor::j)) - Fluid->vtheta * sin(Fluid->theta(tensor::i,tensor::j));

				f << "VECTORS meanFluidVelocity_m/s float" << '\n';
				for(int m=0; m<Fluid->M; m++)
				{
					for(int n=0; n<Fluid->N; n++)
					{
						f << mean(vrx(m,n,Range::all())) << "\t" << 
							 mean(vry(m,n,Range::all())) << "\t" << 
							 "0.0" << '\n';
					}
				}
			}

			Array<double,2> mu(mean(Fluid->oilviscosity,tensor::k));
			f << "SCALARS viscosity_Pas float 1" << '\n';
			f << "LOOKUP_TABLE visco" << '\n';
			for(int m=0; m<Fluid->M; m++)
			{
				for(int n=0; n<Fluid->N; n++)
				{
					for(int q=0; q<Q; q++)
					{
						f << mu(m,n) << '\n';
					}
				}
			}
		}
		f.close();
		

		/*
		//binary
		fstream dat("slipper.bin",ios::out|ios::binary);
		for(int m=0; m<Fluid->M; m++)
		{
			for(int n=0; n<Fluid->N; n++)
			{
				dat.write((char*) &(Fluid->h(m,n)), sizeof(double));
			}
		}

		for(int m=0; m<Fluid->M; m++)
		{
			for(int n=0; n<Fluid->N; n++)
			{
				dat.write((char*) &(Fluid->p(m,n)), sizeof(double));
			}
		}

		for(int m=0; m<Fluid->M; m++)
		{
			for(int n=0; n<Fluid->N; n++)
			{
				for(int q=0; q<Fluid->Q; q++)
				{
					dat.write((char*) &(Fluid->T(m,n,q)), sizeof(double));
				}
			}
		}

		{
			//We actually have to create this value here
			Array<double,3> vrx(Fluid->M,Fluid->N,Fluid->Q);
			Array<double,3> vry(Fluid->M,Fluid->N,Fluid->Q);
			vrx = Fluid->vr * sin(Fluid->theta(tensor::i,tensor::j)) + Fluid->vtheta * cos(Fluid->theta(tensor::i,tensor::j));
			vry = Fluid->vr * cos(Fluid->theta(tensor::i,tensor::j)) - Fluid->vtheta * sin(Fluid->theta(tensor::i,tensor::j));
			
			for(int m=0; m<Fluid->M; m++)
			{
				for(int n=0; n<Fluid->N; n++)
				{
					for(int q=0; q<Fluid->Q; q++)
					{
						dat.write((char*) &(vrx(m,n,q)), sizeof(double));
						dat.write((char*) &(vry(m,n,q)), sizeof(double));
					}
				}
			}
		}
		dat.close();
		*/

	}


};
void CSlipperGap::writevtk(string fileName)
{
	//Update the 3d fluid height
	Fluid.z = Fluid.h(tensor::i,tensor::j)/(Fluid.Q-1)*tensor::k;


	//We need to create a nodal mesh based on the volume centroid of the fluid mesh
	vtk fluid(this);
	
	fluid.buildCheap(false);
	//fluid.build(false);

	//Solve for pnosqz and copy it into the vtk struct
	fluid.pnosqz.resize(Fluid.M,Fluid.N);
	fluid.pnosqz = 0;
/*
	Array<double,2> dht(Fluid.dht.copy());
	Array<double,2> p(Fluid.p.copy());
	Fluid.dht = 0;
	SlipperReynolds(1.0);
	fluid.pnosqz = Fluid.p;
	Fluid.dht = dht;
	Fluid.p = p;
*/
	fluid.writevtk(fileName, operatingslippergap.phi_rad, geometryslippergap.rB);


}
