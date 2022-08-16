#include "fem.h"

namespace CasparSlipperFEM
{
	void fem::stiffnessMatrix()
	{
		//Check the state of fem members to see if we can proceed with a valid analysis
		checkAnalysis();

		//Renumber the unconstrained DOF
		numberDOF();

		//Update the b load vector
		setLoads();

		//Resize the global stiffness matrix and load vector
		gmm::resize(K,uDOF,uDOF);
		gmm::clear(K);
		
		//Average abs K and nnz (number non-zero in K). Used for IR
		double avgK = 0;
		int nnz = 0;

		double lastpdisp = 0;	//used for GapLog'ing progress

		GapLog << "\tProgress: ";
		for(int e=0; e<elecnt; e++)
		{
			//used for GapLog'ing progress
			{
				double p = double(e)/double(elecnt);
				if(p  > lastpdisp + 0.05)
				{
					GapLog << ".";
					lastpdisp = p;
				}
			}

			//local element stiffness matrix
			dmatrix<double> Kloc = elements[e]->K();

			vector<double> Bloc = elements[e]->F();
			
			//Obtain the local2global DOF map
			vector<int> l2g = elements[e]->l2gDOF();
			
			//Assemble the local stiffness matrix into the global matrix	
			for(int i=0; i<Kloc.m; i++)
			{
				if(l2g[i] < 0)	//skip constrained nodes
				{
					continue;
				}

				//add load vector
				b[l2g[i]] += Bloc[i];

				for(int j=0; j<Kloc.n; j++)
				{
					if(l2g[j] >= 0)
					{
						K(l2g[i], l2g[j]) += Kloc(i,j);
						
						avgK += fabs(Kloc(i,j));
						nnz++;
					} else {
						//DOF is constrained so move to b vector
						b[l2g[i]] -= Kloc(i,j) * ebcs[elements[e]->niddof(j)];
					}
				}		
			}

		}

		//GapLog finished progress
		GapLog << ". finished!" << endl;
		
		//Set inertia relief
		apply_inertia_relief(avgK / double(nnz));
	};

	int fem::get_nodecnt()
	{
		return nodecnt;
	}

	int fem::get_elecnt()
	{
		return elecnt;
	}

	double fem::getele_temp(const int eid)
	{
		if(eid < 0 || eid > elecnt)
		{
			return 0;
		}

		double phi = 0;
		double cnt = 0;

		//get an average val for DOF 1
		for(int n=0; n<elements[eid]->nodes.size(); n++)
		{
			const int Gdof = elements[eid]->nodes[n]->DOF[0];
			if(Gdof >= 0)
			{
				phi += X[Gdof];
			} else {
				phi += ebcs[niddof(elements[eid]->nodes[n]->id,0)];
			}
			cnt++;
		}

		return phi/cnt;
	}

	double fem::getnode_temp(const int nid)
	{
		if(nid < 0 || nid > nodecnt || analysis_type != THERMAL)
		{
			return 0;
		}

		const int Gdof = nodes[nid].DOF[0];
		if(Gdof < 0)
		{
			return ebcs[niddof(nid, 0)];
		}

		return X[Gdof];
	}

	vector<double> fem::getnode_deform(const int nid)
	{
		if(nid < 0 || nid > nodecnt || analysis_type != ELASTIC)
		{
			return vector<double>(3,0);
		}

		vector<double> deform(3);

		for(int i=0; i<3; i++)
		{
			const int Gdof = nodes[nid].DOF[i];
			if(Gdof < 0)
			{
				deform[i] = ebcs[niddof(nid, i)];
			} else {
				deform[i] = X[Gdof];
			}
		}

		return deform;
	}


	vector<double> inverseI(const vector<double> A)
	{
		vector<double> B(9, 0);

		//det(A)
		double det = A[0]*( A[4]*A[8]-A[7]*A[5] )
						  - A[1]*( A[3]*A[8]-A[5]*A[6] )
						  + A[2]*( A[3]*A[7]-A[4]*A[6] );
		
		//Inverse of A
		double invdet = 1.0/det + 1.0e-12;
		B[0] =  ( A[4]*A[8] - A[7]*A[5] ) * invdet;
		B[3] = -( A[1]*A[8] - A[2]*A[7] ) * invdet;
		B[6] =  ( A[1]*A[5] - A[2]*A[4] ) * invdet;
		B[1] = -( A[3]*A[8] - A[5]*A[6] ) * invdet;
		B[4] =  ( A[0]*A[8] - A[2]*A[6] ) * invdet;
		B[7] = -( A[0]*A[5] - A[3]*A[2] ) * invdet;
		B[2] =  ( A[3]*A[7] - A[6]*A[4] ) * invdet;
		B[5] = -( A[0]*A[7] - A[6]*A[1] ) * invdet;
		B[8] =  ( A[0]*A[4] - A[3]*A[1] ) * invdet;

		return B;
		
	}

	void fem::apply_inertia_relief(const double avgK)
	{
		if(!inertia_relief)
		{
			//don't do anything if not using IR
			return;
		}

		GapLog << endl << "  * Inertia Relief ..." << endl;
		GapLog << "     * Constraints Removed ..." << endl;

		//Mass vectors
		vector<double> nodeMass(nodecnt, 0);
		vector<double> eleMass(elecnt, 0);

		//Body mass, COG
		double Mtot = 0;
		double xcg = 0, ycg = 0, zcg = 0;

		//Find the ele mass, total mass, nodal mass, and COF
		for(int i=0; i<elecnt; i++)
		{
			element * e = elements[i];

			//Find element mass
			eleMass[i] = e->V()*e->rho;

			//Set nodal mass
			const int nn = (int) e->nodes.size();
			for(int n=0; n<nn; n++)
			{
				nodeMass[e->nodes[n]->id] += eleMass[i] / double(nn);
			}

			//COG
			point p = e->get_center();
			xcg += eleMass[i]*p.x();
			ycg += eleMass[i]*p.y();
			zcg += eleMass[i]*p.z();

			Mtot += eleMass[i];
		}
		xcg /= Mtot;
		ycg /= Mtot;
		zcg /= Mtot;

		//Calculate the inertia tensor of each cell and body inertia tensor
		vector< vector<double> > I(elecnt);
		vector<double> Itot(9, 0);

		for(int i=0; i<elecnt; i++)
		{
			element * e = elements[i];

			//Find the cell centroid position vector from the COG
			point p = e->get_center();
			double rx = p.x()-xcg;
			double ry = p.y()-ycg;
			double rz = p.z()-zcg;

			//Find cell inertial tensor
			vector<double> Icell(9);

			//Ixx
			Icell[0] = eleMass[i]*(ry*ry + rz*rz);
			//Iyy
			Icell[4] = eleMass[i]*(rx*rx + rz*rz);
			//Izz
			Icell[8] = eleMass[i]*(rx*rx + ry*ry);

			//Ixy
			Icell[1] = Icell[3] = -eleMass[i]*rx*ry;
			//Ixz
			Icell[2] = Icell[6] = -eleMass[i]*rx*rz;
			//Iyz
			Icell[5] = Icell[7] = -eleMass[i]*ry*rz;

			I[i] = Icell;

			//The total body tensor
			for(int i=0; i<9; i++)
			{
				Itot[i] += Icell[i];
			}
		}
		
		//Inverse of total body inertial tensor
		vector<double> Itot_Inv(inverseI(Itot));

		//Calculate net force & torque about COG
		vector<double> F(3, 0);
		vector<double> T(3, 0);
		for(int i=0; i<nodes.size(); i++)
		{
			//loop through x,y,z
			for(int d=0; d<3; d++)
			{
				const int Gdof = nodes[i].DOF[d];		
				if(Gdof < 0) //because this is IR, all ebc should be removed. but just to be safe
				{
					continue;
				}
				const double f = b[Gdof];
				
				//force
				F[d] += f;
				
				//torque
				if(d == 0)
				{
					//Ty += Fx*Rz
					T[1] +=  f * (nodes[i].z() - zcg);
					//Tz += -Fx*Ry
					T[2] += -f * (nodes[i].y() - ycg);
				} else if (d == 1)
				{
					//Tx += -Fy*Rz
					T[0] +=  -f * (nodes[i].z() - zcg);
					//Tz += Fy*Rx
					T[2] += f * (nodes[i].x() - xcg);
				} else if (d == 2)
				{
					//Tx += Fz*Ry
					T[0] +=  f * (nodes[i].y() - ycg);
					//Ty += -Fz*Rx
					T[1] += -f * (nodes[i].x() - xcg);
				}
				
			}
		}

		GapLog << "     M: " << Mtot << endl;
		GapLog << "     COG = [ " << xcg << "\t" << ycg << "\t" << zcg << " ]" << endl;
		GapLog << "     dF = [ " << F[0] << "\t" << F[1] << "\t" << F[2] << " ]" << endl;
		GapLog << "     dM = [ " << T[0] << "\t" << T[1] << "\t" << T[2] << " ]" << endl;
		GapLog << "     I = [ " << Itot[0] << "\t" << Itot[1] << "\t" << Itot[2] << endl
			  << "       " << Itot[3] << "\t" << Itot[4] << "\t" << Itot[5] << endl
			  << "       " << Itot[6] << "\t" << Itot[7] << "\t" << Itot[8] << " ]" << endl;


		//Linear acceleration
		vector<double> A_lin(3,0);
		for(int i=0; i<3; i++)
		{
			A_lin[i] = - F[i] / Mtot;
		}

		//Angular acceleration
		vector<double> A_ang(3,0);
		for(int i=0; i<9; i++)
		{
			A_ang[i/3] += - Itot_Inv[i] * T[i%3];
		}

		GapLog << "     a = [ " << A_lin[0] << "\t" << A_lin[1] << "\t" << A_lin[2] << " ] m/s^2" << endl;
		GapLog << "     alpha = [ " << A_ang[0] << "\t" << A_ang[1] << "\t" << A_ang[2] << " ] rad/s^2" << endl;


		//Calculate the inertia load on each element
		for(int e=0; e<elecnt; e++)
		{
			//Element force
			vector<double> F_lin(3,0);
			for(int i=0; i<3; i++)
			{
				F_lin[i] = eleMass[e]*A_lin[i];
			}

			//Cell torque
			vector<double> T_cell(3,0);
			for(int i=0; i<9; i++)
			{
				T_cell[i/3] += I[e][i] * A_ang[i%3];
			}

			//Cell centroid position vector
			point p = elements[e]->get_center();
			double rx = p.x()-xcg;
			double ry = p.y()-ycg;
			double rz = p.z()-zcg;

			vector<double> F_ang(3,0);
			double norm_r = pow(rx*rx+ry*ry+rz*rz,0.5);
			F_ang[0] = 1.0/pow(norm_r,2.0) * (T_cell[1]*rz-T_cell[2]*ry);
			F_ang[1] = 1.0/pow(norm_r,2.0) * (T_cell[2]*rx-T_cell[0]*rz);
			F_ang[2] = 1.0/pow(norm_r,2.0) * (T_cell[0]*ry-T_cell[1]*rx);

			vector<double> F_tot(3,0);
			for(int i=0; i<3; i++)
			{
				F_tot[i] = (F_lin[i]+F_ang[i]) / double(elements[e]->nodes.size());
			}

			//Apply the inertial reaction force to the cell nodes
			for(int n=0; n<elements[e]->nodes.size(); n++)
			{
				for(int i=0; i<3; i++)
				{
					b[ elements[e]->nodes[n]->DOF[i] ] += F_tot[i];
				}
			}	
		}


		//Calculate net force & torque about COG now that the inertial reaction loads have been applied
		//Just for checking, not actually needed
		{
			vector<double> F(3, 0);
			vector<double> T(3, 0);
			for(int i=0; i<nodes.size(); i++)
			{
				//loop through x,y,z
				for(int d=0; d<3; d++)
				{
					const int Gdof = nodes[i].DOF[d];		
					if(Gdof < 0) //because this is IR, all ebc should be removed. but just to be safe
					{
						continue;
					}
					const double f = b[Gdof];
					
					//force
					F[d] += f;
					
					//torque
					if(d == 0)
					{
						//Ty += Fx*Rz
						T[1] +=  f * (nodes[i].z() - zcg);
						//Tz += -Fx*Ry
						T[2] += -f * (nodes[i].y() - ycg);
					} else if (d == 1)
					{
						//Tx += -Fy*Rz
						T[0] +=  -f * (nodes[i].z() - zcg);
						//Tz += Fy*Rx
						T[2] += f * (nodes[i].x() - xcg);
					} else if (d == 2)
					{
						//Tx += Fz*Ry
						T[0] +=  f * (nodes[i].y() - ycg);
						//Ty += -Fz*Rx
						T[1] += -f * (nodes[i].x() - xcg);
					}
				}
			}

			GapLog << "     * Inertial Loads Applied ..." << endl;
			GapLog << "     dF = [ " << F[0] << "\t" << F[1] << "\t" << F[2] << " ]" << endl;
			GapLog << "     dM = [ " << T[0] << "\t" << T[1] << "\t" << T[2] << " ]" << endl;
			GapLog << endl;
		}



		//Add the lambda const to the K matrix

		//find the avg nodal mass
		double avgNodeMass = 0;
		for(int n=0; n<nodecnt; n++)
		{
			avgNodeMass += nodeMass[n];
		}
		avgNodeMass /= double(avgNodeMass);

		//Scale factor = [abs mean stiffness matrix] / [mean node mass]
		double scale = avgK/avgNodeMass;
		for(int n=0; n<nodecnt; n++)
		{
			for(int d=0; d<3; d++)
			{
				const int Gdof = nodes[n].DOF[d];

				K(uDOF - 3 + d, Gdof) = nodeMass[n] * scale;
				K(Gdof, uDOF - 3 + d) = nodeMass[n] * scale;
			}
		}



	};
};
