#include "CSlipperGap.h"
#include <omp.h>
#include <time.h>
#pragma once
#define OMP

extern void matlab(const Array<double,2>& data,const string file);

int CSlipperGap::SlipperEnergy() {


	// fluid velocities at boundary faces
	Array<double,3> vr_s(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> vr_n(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> vtheta_w(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> vtheta_e(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> vr_b(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> vr_t(Fluid.M,Fluid.N,Fluid.Q); 
	Array<double,3> vtheta_b(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> vtheta_t(Fluid.M,Fluid.N,Fluid.Q);

	//Central difference method to find face velocities
	// radial velocity on north faces
	vr_n(0,all,all) = 0.0;
	vr_n(Fluid.M-1,all,all) = 0.0;
	vr_n(Range(1,Fluid.M-2),all,all) = 0.50*(Fluid.vr(Range(1,Fluid.M-2),all,all) + Fluid.vr(Range(2,Fluid.M-1),all,all));
	// radial velocity on south faces
	vr_s(0,all,all) = 0.0;
	vr_s(Fluid.M-1,all,all) = 0.0;
	vr_s(Range(1,Fluid.M-2),all,all) = 0.50*(Fluid.vr(Range(1,Fluid.M-2),all,all) + Fluid.vr(Range(0,Fluid.M-3),all,all));
	// circumferential velocity on east faces
	vtheta_e(all,Range(0,Fluid.N-2),all) = 0.50*(Fluid.vtheta(all,Range(0,Fluid.N-2),all) + Fluid.vtheta(all,Range(1,Fluid.N-1),all));
	vtheta_e(all,Fluid.N-1,all) = 0.50*(Fluid.vtheta(all,Fluid.N-1,all) + Fluid.vtheta(all,0,all));
	// circumferential velocity on west faces
	vtheta_w(all,0,all) = 0.50*(Fluid.vtheta(all,Fluid.N-1,all) + Fluid.vtheta(all,0,all));
	vtheta_w(all,Range(1,Fluid.N-1),all) = 0.50*(Fluid.vtheta(all,Range(0,Fluid.N-2),all) + Fluid.vtheta(all,Range(1,Fluid.N-1),all));
	// radial velocity on top faces
	vr_t(all,all,0) = 0.0;
	vr_t(all,all,Fluid.Q-1) = 0.0;
	vr_t(all,all,Range(1,Fluid.Q-2)) = 0.50*(Fluid.vr(all,all,Range(1,Fluid.Q-2)) + Fluid.vr(all,all,Range(2,Fluid.Q-1)));
	// radial velocity on bottom faces
	vr_b(all,all,0) = 0.0;
	vr_b(all,all,Fluid.Q-1) = 0.0;
	vr_b(all,all,Range(1,Fluid.Q-2)) = 0.50*(Fluid.vr(all,all,Range(1,Fluid.Q-2)) + Fluid.vr(all,all,Range(0,Fluid.Q-3)));
	// circumferential velocity on top faces
	vtheta_t(all,all,0) = 0.0;
	vtheta_t(all,all,Range(1,Fluid.Q-2)) = 0.50*(Fluid.vtheta(all,all,Range(1,Fluid.Q-2)) + Fluid.vtheta(all,all,Range(2,Fluid.Q-1)));
	vtheta_t(all,all,Fluid.Q-1) = 0.0;
	// circumferentia velocity on bottom faces
	vtheta_b(all,all,0) = 0.0;
	vtheta_b(all,all,Range(1,Fluid.Q-2)) = 0.50*(Fluid.vtheta(all,all,Range(1,Fluid.Q-2)) + Fluid.vtheta(all,all,Range(0,Fluid.Q-3)));
	vtheta_b(all,all,Fluid.Q-1) = 0.0;

	// mesh motion velocities
	Array<double,2> mesh_vr_s(Fluid.M,Fluid.N);
	Array<double,2> mesh_vr_n(Fluid.M,Fluid.N);
	Array<double,2> mesh_vtheta_w(Fluid.M,Fluid.N);
	Array<double,2> mesh_vtheta_e(Fluid.M,Fluid.N);

	// mesh motion radial velocity on north faces
	mesh_vr_n(0,all) = operatingslippergap.vgr(0,all);
	mesh_vr_n(Fluid.M-1,all) = operatingslippergap.vgr(Fluid.M-1,all);
	mesh_vr_n(Range(1,Fluid.M-2),all) = 0.50*(operatingslippergap.vgr(Range(1,Fluid.M-2),all) + operatingslippergap.vgr(Range(2,Fluid.M-1),all,all));
	// mesh motion radial velocity on south faces
	mesh_vr_s(0,all) = operatingslippergap.vgr(0,all);
	mesh_vr_s(Fluid.M-1,all) = operatingslippergap.vgr(Fluid.M-1,all);
	mesh_vr_s(Range(1,Fluid.M-2),all) = 0.50*(operatingslippergap.vgr(Range(1,Fluid.M-2),all) + operatingslippergap.vgr(Range(0,Fluid.M-3),all));
	// mesh motion circumferential velocity on east faces 
	mesh_vtheta_e(all,Range(0,Fluid.N-2)) = 0.50*(operatingslippergap.vgtheta(all,Range(0,Fluid.N-2)) + operatingslippergap.vgtheta(all,Range(1,Fluid.N-1)));
	mesh_vtheta_e(all,Fluid.N-1) = 0.50*(operatingslippergap.vgtheta(all,Fluid.N-1) + operatingslippergap.vgtheta(all,0));
	// mesh motion circumferential velocity on west faces
	mesh_vtheta_w(all,0) = 0.50*(operatingslippergap.vgtheta(all,Fluid.N-1) + operatingslippergap.vgtheta(all,0));
	mesh_vtheta_w(all,Range(1,Fluid.N-1)) = 0.50*(operatingslippergap.vgtheta(all,Range(0,Fluid.N-2)) + operatingslippergap.vgtheta(all,Range(1,Fluid.N-1)));

	// dz
	Array<double,2> dz(Fluid.M,Fluid.N);
	dz = Fluid.h/(Fluid.Q - 1);
	
	// ------------------------- coefficient due to diffusion ------------------------- //

	// oil thermal diffusivity
	//const double lambda = oilslippergap.oillambda;
	const double lambda = oil_properties->get_lambda();
	// oil heat capacity
	//const double cp = oilslippergap.oilC;
	const double cp = oil_properties->get_C();

	Array<double,2> DW(Fluid.M,Fluid.N);
	DW = (lambda*(Fluid.dr)*(dz))/(Fluid.r*(Fluid.dtheta));
	Array<double,2> DE(Fluid.M,Fluid.N);	
	DE = (lambda*(Fluid.dr)*(dz))/(Fluid.r*(Fluid.dtheta));
	Array<double,2> DS(Fluid.M,Fluid.N);
	DS = (lambda*(Fluid.r - Fluid.dr)*(Fluid.dtheta)*(dz))/(Fluid.dr);
	Array<double,2> DN(Fluid.M,Fluid.N);
	DN = (lambda*(Fluid.r + Fluid.dr)*(Fluid.dtheta)*(dz))/(Fluid.dr);
	Array<double,2> DB(Fluid.M,Fluid.N);
	DB = (lambda*(Fluid.r)*(Fluid.dtheta)*(Fluid.dr))/(dz);
	Array<double,2> DT(Fluid.M,Fluid.N);
	DT = (lambda*(Fluid.r)*(Fluid.dtheta)*(Fluid.dr))/(dz);

	// ------------------------- coefficients due to convection ------------------------- //

	firstIndex ii;
	secondIndex jj;
	thirdIndex kk;

	Array<double,3> FS(Fluid.M,Fluid.N,Fluid.Q);
	FS = Fluid.oildensity(ii,jj,kk)*cp*(Fluid.r(ii,jj) - Fluid.dr(ii,jj))*(Fluid.dtheta(ii,jj))*(dz(ii,jj))*(vr_s(ii,jj,kk) - mesh_vr_s(ii,jj));
	Array<double,3> FN(Fluid.M,Fluid.N,Fluid.Q);
	FN = Fluid.oildensity(ii,jj,kk)*cp*(Fluid.r(ii,jj) + Fluid.dr(ii,jj))*(Fluid.dtheta(ii,jj))*(dz(ii,jj))*(vr_n(ii,jj,kk) - mesh_vr_n(ii,jj));
	Array<double,3> FW(Fluid.M,Fluid.N,Fluid.Q);
	FW = Fluid.oildensity(ii,jj,kk)*cp*(dz(ii,jj))*(Fluid.dr(ii,jj))*(vtheta_w(ii,jj,kk) - mesh_vtheta_w(ii,jj));
	Array<double,3> FE(Fluid.M,Fluid.N,Fluid.Q);
	FE = Fluid.oildensity(ii,jj,kk)*cp*(dz(ii,jj))*(Fluid.dr(ii,jj))*(vtheta_e(ii,jj,kk) - mesh_vtheta_e(ii,jj));

	// ------------------------- coefficients for power law scheme ------------------------- //

	Array<double,3> FDW(Fluid.M,Fluid.N,Fluid.Q);
	FDW = blitz::pow(1.0 - 0.1*abs(FW(ii,jj,kk)/DW(ii,jj)),5.0);
	Array<double,3> FDE(Fluid.M,Fluid.N,Fluid.Q);
	FDE = blitz::pow(1.0 - 0.1*abs(FE(ii,jj,kk)/DE(ii,jj)),5.0);
	Array<double,3> FDS(Fluid.M,Fluid.N,Fluid.Q);
	FDS = blitz::pow(1.0 - 0.1*abs(FS(ii,jj,kk)/DS(ii,jj)),5.0);
	Array<double,3> FDN(Fluid.M,Fluid.N,Fluid.Q);
	FDN = blitz::pow(1.0 - 0.1*abs(FN(ii,jj,kk)/DN(ii,jj)),5.0);

	// ------------------------- coefficients of difference equation ------------------------- //

	Array<double,3> aW(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> aE(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> aS(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> aN(Fluid.M,Fluid.N,Fluid.Q);

	Array<double,3> F(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> FD(Fluid.M,Fluid.N,Fluid.Q);

	// west
	F = where(FW > 0.0, FW, 0.0);
	FD = where(FDW > 0.0, FDW, 0.0);
	aW = DW(ii,jj)*FD(ii,jj,kk) + F(ii,jj,kk);
	// east
	F = where(-1.0*FE > 0.0, -1.0*FE, 0.0);
	FD = where(FDE > 0.0, FDE, 0.0);
	aE = DE(ii,jj)*FD(ii,jj,kk) + F(ii,jj,kk);	
	// south
	F = where(FS > 0.0, FS, 0.0);
	FD = where(FDS > 0.0, FDS, 0.0);
	aS = DS(ii,jj)*FD(ii,jj,kk) + F(ii,jj,kk);
	// north
	F = where(-1.0*FN > 0.0, -1.0*FN, 0.0);
	FD = where(FDN > 0.0, FDN, 0.0);
	aN = DN(ii,jj)*FD(ii,jj,kk) + F(ii,jj,kk);

	Array<double,3> aB(Fluid.M,Fluid.N,Fluid.Q); // only diffusion on bottom 
	aB = DB(ii,jj,kk); // only diffusion on bottom 
	Array<double,3> aT(Fluid.M,Fluid.N,Fluid.Q); // only diffusion on top
	aT = DT(ii,jj,kk); // only diffusion on top
	
	Array<double,3> aP(aW + aE + aS + aN + aB + aT); // aP coefficient, equal to the sum of all others

	// source term (dissipation function in polar coordinates)
	Array<double,3> bb(Fluid.M,Fluid.N,Fluid.Q);


	bb = Fluid.oilviscosity(ii,jj,kk)*(
								(Fluid.r(ii,jj))*(Fluid.dtheta(ii,jj))*(Fluid.dr(ii,jj))*(dz(ii,jj)))*( // finite volume
								blitz::pow((vr_t(ii,jj,kk) - vr_b(ii,jj,kk))/(dz(ii,jj)),2.0) + // dur/dz
								blitz::pow((vtheta_t(ii,jj,kk) - vtheta_b(ii,jj,kk))/(dz(ii,jj)),2.0) + // dutheta/dz
								(4.0/3.0)*blitz::pow((Fluid.vr(ii,jj,kk)/Fluid.r(ii,jj)),2.0) + //4/3*(ur/r)^2
								blitz::pow((Fluid.vtheta(ii,jj,kk)/Fluid.r(ii,jj)),2.0) // (utheta/r)^2
							);
	
	// --------------------------------- initialization of solver --------------------------------- //

	double TB(0), TT(0), TW(0), TE(0), TS(0), TN(0), TP(0); // pressure on each element
	Array<double,3> TOld(Fluid.M,Fluid.N,Fluid.Q); // pressur of previous iteration
	Array<double,3> errors(Fluid.M,Fluid.N,Fluid.Q); // error for each element

	const double alpha = gapinput->options_slipper.numeric.AlphaEnergy; // overrelaxation corfficient
	
	int iterations(0); // iterations number
	
	// --------------------------------------  SOR loop  -------------------------------------- //
double tmp;
	do {

		TOld = Fluid.T; // old T field (from previous iteration)

		for (int k=1; k<Fluid.Q-1; k++) { // all the layers except first and last (wich are boundaries)
			// all the radius except first and last (wich are boundaries)
			for(int i=1; i<Fluid.M-1; i++) { 
				// all the circles 
				for(int j=0; j<Fluid.N; j++) { 

					if(Fluid.boundary(i,j) == -1)
					{
				
						// --------------------------- east and west --------------------------- //
						TE = (j != Fluid.N-1) ? Fluid.T(i,j+1,k) : Fluid.T(i,0,k);
						TW = (j != 0) ? Fluid.T(i,j-1,k) : Fluid.T(i,Fluid.N-1,k);
											
						// --------------------------- north and south --------------------------- //
						TN = Fluid.T(i+1,j,k);
						TS = Fluid.T(i-1,j,k);

						// --------------------------- top and bottom --------------------------- //
						TB = Fluid.T(i,j,k-1);
						TT = Fluid.T(i,j,k+1);

						// -------------------- new temperature value for element P -------------------- //
						tmp	=	
									Fluid.T(i,j,k) + 
									alpha*(
										(
											aW(i,j,k)*TW + aE(i,j,k)*TE + 
											aS(i,j,k)*TS + aN(i,j,k)*TN + 
											aB(i,j,k)*TB + aT(i,j,k)*TT + 
											bb(i,j,k) 
										)
										/aP(i,j,k) - Fluid.T(i,j,k)
									);	
						Fluid.T(i,j,k) = tmp;
					}
				}
			}
		}
	
		//0 flux at slip & swash
		//Fluid.T(Range::all(), Range::all(), 0) = Fluid.T(Range::all(), Range::all(), 1);
		//Fluid.T(Range::all(), Range::all(), Fluid.Q-1) = Fluid.T(Range::all(), Range::all(), Fluid.Q-2);

		errors = blitz::fabs(Fluid.T - TOld); // calculate the error at the end of iteration
		
		iterations++; // increase the iterator

	} while (max(errors) > 1e-1 && iterations < 1e3);

	if(max(errors) > 1e-1 && iterations >= 1e3)
	{
		GapLog.message("WARNING: GS Energy temperature loop failed to converge!");

		//to prevent the thermal from crashing if this contains NaN's
		SlipperInitializeTemperature();
		Fluid.Qflux = 0;
	}

	//														//
	//Calculate a Qflux : heat transfer into the solid bodies
	//There are a number of methods to choose from
	//														//
	/*
	//Sum of viscous source (send only 0.5 to the slipper & 0.5 to the swashplate)
	Fluid.Qflux = 0.5*sum(bb(Range::all(), Range::all(), Range(1,Fluid.Q-2)),tensor::k) / Fluid.dA;
	Fluid.Qflux(0, Range::all()) = 0;
	Fluid.Qflux(Fluid.M-1, Range::all()) = 0;
	*/
	//"Proper" heat transfer (gradient poor when used on a large mesh)  (send only 0.5 to the slipper & 0.5 to the swashplate)
	
	{
		////Calculate the flux gradient using cells 0 & 1
		//Range all = Range::all();
		//Fluid.Qflux =	(	aB(all,all,1)*(Fluid.T(all, all, 1)-Fluid.T(all, all, 0))
		//						+
		//						aT(all,all,Fluid.Q-2)*(Fluid.T(all, all, Fluid.Q-2)-Fluid.T(all, all, Fluid.Q-1))
		//					)/Fluid.dA*0.5;

		//Calculate the flux gradient using cells 1 & 2
		Fluid.Qflux =	(	aB(all,all,2)*(Fluid.T(all, all, 2)-Fluid.T(all, all, 1))
								+
							-aT(all,all,Fluid.Q-3)*(Fluid.T(all, all, Fluid.Q-2)-Fluid.T(all, all, Fluid.Q-3))
							)/Fluid.dA*0.5;
	}
	
	/*
	//Advanced qflux calculation methods
	{
		Range mall = Range(1,Fluid.M-2);	//this is actually the non-bounary m
		Range north = Range(2,Fluid.M-1);
		Range south  = Range(0,Fluid.M-3);
		
		Array<double, 3> TEast(Fluid.M,Fluid.N,Fluid.Q);
		TEast(all,Range(0,Fluid.N-2),all) = Fluid.T(all,Range(1,Fluid.N-1),all);
		TEast(all,Fluid.N-1,all) = Fluid.T(all,0,all);
		
		Array<double, 3> TWest(Fluid.M,Fluid.N,Fluid.Q);
		TWest(all,Range(1,Fluid.N-1),all) = Fluid.T(all,Range(0,Fluid.N-2),all);
		TWest(all,0,all) = Fluid.T(all,Fluid.N-1,all);
		
		Fluid.Qflux(0,all) = 0;
		Fluid.Qflux(Fluid.M-1,all) = 0;
		
		
		////Consider only convective heat transfer in the x & y directions
		//Fluid.Qflux(mall,all) = sum(
		//							bb(mall,all,all)+(
		//							-FN(mall,all,all)*where(FN < 0, Fluid.T(north,all,all), Fluid.T(mall,all,all)) +
		//							FS(mall,all,all)*where(FS > 0, Fluid.T(south,all,all), Fluid.T(mall,all,all)) +
		//							-FE(mall,all,all)*where(FE < 0, TEast, Fluid.T(mall,all,all)) +
		//							FW(mall,all,all)*where(FW > 0, TWest, Fluid.T(mall,all,all))
		//							),tensor::k)*0.5/Fluid.dA(mall,all);
		//

		//Consider both diffusive and convective (per the power law) in the x & y directions
		Fluid.Qflux(mall,all) =	sum( bb(mall,all,all) +
									 aN(mall,all,all) * Fluid.T(north,all,all) +
									 aS(mall,all,all) * Fluid.T(south,all,all) +
									 aE(mall,all,all) * TEast(mall,all,all) +
									 aW(mall,all,all) * TWest(mall,all,all) -
									 (aP(mall,all,all)-aB(mall,all,all)-aT(mall,all,all)) * Fluid.T(mall,all,all)
									, tensor::k)/Fluid.dA(mall,all)*0.5;
	}
	*/

	//limit the max fluxes
	Fluid.Qflux = where(Fluid.Qflux > 5e5, 5e5, Fluid.Qflux);
	Fluid.Qflux = where(Fluid.Qflux < -1e5, -1e5, Fluid.Qflux);

	Fluid.Qphisum = sum(bb);

	return iterations;
	
}
int CSlipperGap::SlipperEnergyCG() {
	

	// velocities at boundary faces
	Array<double,3> vr_s(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> vr_n(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> vtheta_w(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> vtheta_e(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> vr_b(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> vr_t(Fluid.M,Fluid.N,Fluid.Q); 
	Array<double,3> vtheta_b(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> vtheta_t(Fluid.M,Fluid.N,Fluid.Q);

	//Central difference method to find face velocities
	// radial velocity on north faces
	vr_n(0,all,all) = 0.0;
	vr_n(Fluid.M-1,all,all) = 0.0;
	vr_n(Range(1,Fluid.M-2),all,all) = 0.50*(Fluid.vr(Range(1,Fluid.M-2),all,all) + Fluid.vr(Range(2,Fluid.M-1),all,all));
	// radial velocity on south faces
	vr_s(0,all,all) = 0.0;
	vr_s(Fluid.M-1,all,all) = 0.0;
	vr_s(Range(1,Fluid.M-2),all,all) = 0.50*(Fluid.vr(Range(1,Fluid.M-2),all,all) + Fluid.vr(Range(0,Fluid.M-3),all,all));
	// circumferential velocity on east faces
	vtheta_e(all,Range(0,Fluid.N-2),all) = 0.50*(Fluid.vtheta(all,Range(0,Fluid.N-2),all) + Fluid.vtheta(all,Range(1,Fluid.N-1),all));
	vtheta_e(all,Fluid.N-1,all) = 0.50*(Fluid.vtheta(all,Fluid.N-1,all) + Fluid.vtheta(all,0,all));
	// circumferential velocity on west faces
	vtheta_w(all,0,all) = 0.50*(Fluid.vtheta(all,Fluid.N-1,all) + Fluid.vtheta(all,0,all));
	vtheta_w(all,Range(1,Fluid.N-1),all) = 0.50*(Fluid.vtheta(all,Range(0,Fluid.N-2),all) + Fluid.vtheta(all,Range(1,Fluid.N-1),all));
	// radial velocity on top faces
	vr_t(all,all,0) = 0.0;
	vr_t(all,all,Fluid.Q-1) = 0.0;
	vr_t(all,all,Range(1,Fluid.Q-2)) = 0.50*(Fluid.vr(all,all,Range(1,Fluid.Q-2)) + Fluid.vr(all,all,Range(2,Fluid.Q-1)));
	// radial velocity on bottom faces
	vr_b(all,all,0) = 0.0;
	vr_b(all,all,Fluid.Q-1) = 0.0;
	vr_b(all,all,Range(1,Fluid.Q-2)) = 0.50*(Fluid.vr(all,all,Range(1,Fluid.Q-2)) + Fluid.vr(all,all,Range(0,Fluid.Q-3)));
	// circumferential velocity on top faces
	vtheta_t(all,all,0) = 0.0;
	vtheta_t(all,all,Range(1,Fluid.Q-2)) = 0.50*(Fluid.vtheta(all,all,Range(1,Fluid.Q-2)) + Fluid.vtheta(all,all,Range(2,Fluid.Q-1)));
	vtheta_t(all,all,Fluid.Q-1) = 0.0;
	// circumferentia velocity on bottom faces
	vtheta_b(all,all,0) = 0.0;
	vtheta_b(all,all,Range(1,Fluid.Q-2)) = 0.50*(Fluid.vtheta(all,all,Range(1,Fluid.Q-2)) + Fluid.vtheta(all,all,Range(0,Fluid.Q-3)));
	vtheta_b(all,all,Fluid.Q-1) = 0.0;

	// mesh motion velocities
	Array<double,2> mesh_vr_s(Fluid.M,Fluid.N);
	Array<double,2> mesh_vr_n(Fluid.M,Fluid.N);
	Array<double,2> mesh_vtheta_w(Fluid.M,Fluid.N);
	Array<double,2> mesh_vtheta_e(Fluid.M,Fluid.N);

	// mesh motion radial velocity on north faces
	mesh_vr_n(0,all) = operatingslippergap.vgr(0,all);
	mesh_vr_n(Fluid.M-1,all) = operatingslippergap.vgr(Fluid.M-1,all);
	mesh_vr_n(Range(1,Fluid.M-2),all) = 0.50*(operatingslippergap.vgr(Range(1,Fluid.M-2),all) + operatingslippergap.vgr(Range(2,Fluid.M-1),all,all));
	// mesh motion radial velocity on south faces
	mesh_vr_s(0,all) = operatingslippergap.vgr(0,all);
	mesh_vr_s(Fluid.M-1,all) = operatingslippergap.vgr(Fluid.M-1,all);
	mesh_vr_s(Range(1,Fluid.M-2),all) = 0.50*(operatingslippergap.vgr(Range(1,Fluid.M-2),all) + operatingslippergap.vgr(Range(0,Fluid.M-3),all));
	// mesh motion circumferential velocity on east faces 
	mesh_vtheta_e(all,Range(0,Fluid.N-2)) = 0.50*(operatingslippergap.vgtheta(all,Range(0,Fluid.N-2)) + operatingslippergap.vgtheta(all,Range(1,Fluid.N-1)));
	mesh_vtheta_e(all,Fluid.N-1) = 0.50*(operatingslippergap.vgtheta(all,Fluid.N-1) + operatingslippergap.vgtheta(all,0));
	// mesh motion circumferential velocity on west faces
	mesh_vtheta_w(all,0) = 0.50*(operatingslippergap.vgtheta(all,Fluid.N-1) + operatingslippergap.vgtheta(all,0));
	mesh_vtheta_w(all,Range(1,Fluid.N-1)) = 0.50*(operatingslippergap.vgtheta(all,Range(0,Fluid.N-2)) + operatingslippergap.vgtheta(all,Range(1,Fluid.N-1)));

	// dz
	Array<double,2> dz(Fluid.M,Fluid.N);
	dz = Fluid.h/(Fluid.Q - 1);
	
	// ------------------------- coefficient due to diffusion ------------------------- //

	// oil thermal diffusivity
	//const double lambda = oilslippergap.oillambda;
	const double lambda = oil_properties->get_lambda();
	// oil heat capacity
	//const double cp = oilslippergap.oilC;
	const double cp = oil_properties->get_C();

	Array<double,2> DW(Fluid.M,Fluid.N);
	DW = (lambda*(Fluid.dr)*(dz))/(Fluid.r*(Fluid.dtheta));
	Array<double,2> DE(Fluid.M,Fluid.N);	
	DE = (lambda*(Fluid.dr)*(dz))/(Fluid.r*(Fluid.dtheta));
	Array<double,2> DS(Fluid.M,Fluid.N);
	DS = (lambda*(Fluid.r - Fluid.dr)*(Fluid.dtheta)*(dz))/(Fluid.dr);
	Array<double,2> DN(Fluid.M,Fluid.N);
	DN = (lambda*(Fluid.r + Fluid.dr)*(Fluid.dtheta)*(dz))/(Fluid.dr);
	Array<double,2> DB(Fluid.M,Fluid.N);
	DB = (lambda*(Fluid.r)*(Fluid.dtheta)*(Fluid.dr))/(dz);
	Array<double,2> DT(Fluid.M,Fluid.N);
	DT = (lambda*(Fluid.r)*(Fluid.dtheta)*(Fluid.dr))/(dz);

	// ------------------------- coefficients due to convection ------------------------- //

	firstIndex ii;
	secondIndex jj;
	thirdIndex kk;

	Array<double,3> FS(Fluid.M,Fluid.N,Fluid.Q);
	FS = Fluid.oildensity(ii,jj,kk)*cp*(Fluid.r(ii,jj) - Fluid.dr(ii,jj))*(Fluid.dtheta(ii,jj))*(dz(ii,jj))*(vr_s(ii,jj,kk) - mesh_vr_s(ii,jj));
	Array<double,3> FN(Fluid.M,Fluid.N,Fluid.Q);
	FN = Fluid.oildensity(ii,jj,kk)*cp*(Fluid.r(ii,jj) + Fluid.dr(ii,jj))*(Fluid.dtheta(ii,jj))*(dz(ii,jj))*(vr_n(ii,jj,kk) - mesh_vr_n(ii,jj));
	Array<double,3> FW(Fluid.M,Fluid.N,Fluid.Q);
	FW = Fluid.oildensity(ii,jj,kk)*cp*(dz(ii,jj))*(Fluid.dr(ii,jj))*(vtheta_w(ii,jj,kk) - mesh_vtheta_w(ii,jj));
	Array<double,3> FE(Fluid.M,Fluid.N,Fluid.Q);
	FE = Fluid.oildensity(ii,jj,kk)*cp*(dz(ii,jj))*(Fluid.dr(ii,jj))*(vtheta_e(ii,jj,kk) - mesh_vtheta_e(ii,jj));

	// ------------------------- coefficients for power law scheme ------------------------- //

	Array<double,3> FDW(Fluid.M,Fluid.N,Fluid.Q);
	FDW = blitz::pow(1.0 - 0.1*abs(FW(ii,jj,kk))/DW(ii,jj),5.0);
	Array<double,3> FDE(Fluid.M,Fluid.N,Fluid.Q);
	FDE = blitz::pow(1.0 - 0.1*abs(FE(ii,jj,kk))/DE(ii,jj),5.0);
	Array<double,3> FDS(Fluid.M,Fluid.N,Fluid.Q);
	FDS = blitz::pow(1.0 - 0.1*abs(FS(ii,jj,kk))/DS(ii,jj),5.0);
	Array<double,3> FDN(Fluid.M,Fluid.N,Fluid.Q);
	FDN = blitz::pow(1.0 - 0.1*abs(FN(ii,jj,kk))/DN(ii,jj),5.0);

	// ------------------------- coefficients of difference equation ------------------------- //

	Array<double,3> aW(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> aE(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> aS(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> aN(Fluid.M,Fluid.N,Fluid.Q);

	Array<double,3> F(Fluid.M,Fluid.N,Fluid.Q);
	Array<double,3> FD(Fluid.M,Fluid.N,Fluid.Q);

	// west
	F = where(FW > 0.0, FW, 0.0);
	FD = where(FDW > 0.0, FDW, 0.0);
	aW = DW(ii,jj)*FD(ii,jj,kk) + F(ii,jj,kk);
	// east
	F = where(-1.0*FE > 0.0, -1.0*FE, 0.0);
	FD = where(FDE > 0.0, FDE, 0.0);
	aE = DE(ii,jj)*FD(ii,jj,kk) + F(ii,jj,kk);	
	// south
	F = where(FS > 0.0, FS, 0.0);
	FD = where(FDS > 0.0, FDS, 0.0);
	aS = DS(ii,jj)*FD(ii,jj,kk) + F(ii,jj,kk);
	// north
	F = where(-1.0*FN > 0.0, -1.0*FN, 0.0);
	FD = where(FDN > 0.0, FDN, 0.0);
	aN = DN(ii,jj)*FD(ii,jj,kk) + F(ii,jj,kk);

	Array<double,3> aB(Fluid.M,Fluid.N,Fluid.Q); // only diffusion on bottom 
	aB = DB(ii,jj,kk); // only diffusion on bottom 
	Array<double,3> aT(Fluid.M,Fluid.N,Fluid.Q); // only diffusion on top
	aT = DT(ii,jj,kk); // only diffusion on top
	
	Array<double,3> aP(aW + aE + aS + aN + aB + aT); // aP coefficient, equal to the sum of all others

	// source term (dissipation function in polar coordinates)
	Array<double,3> bb(Fluid.M,Fluid.N,Fluid.Q);


	bb = Fluid.oilviscosity(ii,jj,kk)*(
								(Fluid.r(ii,jj))*(Fluid.dtheta(ii,jj))*(Fluid.dr(ii,jj))*(dz(ii,jj)))*( // finite volume
								blitz::pow((vr_t(ii,jj,kk) - vr_b(ii,jj,kk))/(dz(ii,jj)),2.0) + // dur/dz
								blitz::pow((vtheta_t(ii,jj,kk) - vtheta_b(ii,jj,kk))/(dz(ii,jj)),2.0) + // dutheta/dz
								(4.0/3.0)*blitz::pow((Fluid.vr(ii,jj,kk)/Fluid.r(ii,jj)),2.0) + //4/3*(ur/r)^2
								blitz::pow((Fluid.vtheta(ii,jj,kk)/Fluid.r(ii,jj)),2.0) // (utheta/r)^2
							);
	
	// --------------------------------- initialization of solver --------------------------------- //

	//We need to redefine a 3d boundary
	Array<int,3> boundary(Fluid.M,Fluid.N,Fluid.Q);
	boundary = where(Fluid.boundary(tensor::i,tensor::j) == -1, -1 , 0);
	//Set the top and bottom boundaries
	boundary(Range::all(), Range::all(), 0) = 0;
	boundary(Range::all(), Range::all(), Fluid.Q-1) = 0;
	
	// --------------------------------------  CG loop  -------------------------------------- //
	//Renumber the non-boundary cells with mid
	Array<int,3> mid(Fluid.M,Fluid.N,Fluid.Q);
	mid = -1;
	
	int cells = 0;
	{
		int id=0;
		for(int i=0; i<Fluid.M; i++)
		{
			for(int j=0; j<Fluid.N; j++)
			{
				for(int k=0; k<Fluid.Q; k++)
				{
					if(boundary(i,j,k) == -1)
					{
						mid(i,j,k) = id;
						id++;
						cells++;
					}
				}
			}
		}
	}

	//The A Matrix
	gmm::col_matrix<gmm::wsvector<double> > A(cells,cells);

	//The B & X  Vector
	vector<double> B(cells),X(cells);

	//Build the K array inside
	for(int i=0;i<Fluid.M;i++)
	{
		for(int j=0;j<Fluid.N;j++)	
		{
			for(int k=0;k<Fluid.Q;k++)	
			{
				if(boundary(i,j,k) == -1)
				{
					const int row = mid(i,j,k);
					
					//Copy the initial T to X
					X[row] = Fluid.T(i,j,k);

					int x,y;
					if (j==0)
					{
						x = Fluid.N - 1;
					} 
					else 
					{ 
						x = j - 1;
					}
					if (j==Fluid.N-1)
					{
						y = 0;
					} 
					else 
					{ 
						y = j + 1;
					}

					B[row] = bb(i,j,k);

					if(mid(i-1,j,k) >= 0)
					{
						A(row,mid(i-1,j,k)) = -aS(i,j,k);
					} else {
						//boundary
						B[row] += aS(i,j,k)*Fluid.T(i-1,j,k);
					}

					if(mid(i+1,j,k) >= 0)
					{
						A(row,mid(i+1,j,k)) = -aN(i,j,k);
					} else {
						//boundary
						B[row] += aN(i,j,k)*Fluid.T(i+1,j,k);
					}

					if(mid(i,y,k) >= 0)
					{
						A(row,mid(i,y,k)) = -aE(i,j,k);
					} else {
						//boundary
						B[row] += aE(i,j,k)*Fluid.T(i,y,k);
					}

					if(mid(i,x,k) >= 0)
					{
						A(row,mid(i,x,k)) = -aW(i,j,k);
					} else {
						//boundary
						B[row] += aW(i,j,k)*Fluid.T(i,x,k);
					}

					if(mid(i,j,k-1) >= 0)
					{
						A(row,mid(i,j,k-1)) = -aB(i,j,k);
					} else {
						//boundary
						B[row] += aB(i,j,k)*Fluid.T(i,j,k-1);
					}

					if(mid(i,j,k+1) >= 0)
					{
						A(row,mid(i,j,k+1)) = -aT(i,j,k);
					} else {
						//boundary
						B[row] += aT(i,j,k)*Fluid.T(i,j,k+1);
					}

					A(row,row) = aP(i,j,k);
				}
			}
		}
	}

	//Set up the proper sparse matrix
	gmm::csc_matrix<double> AS;
	gmm::copy(A,AS);

	//Set solver params
	gmm::iteration iter(1e-5);
	iter.set_maxiter(500);
	iter.set_noisy(0);

	//Precondition
	gmm::ildlt_precond< gmm::csc_matrix<double> > PR(AS);

	//Solve
	gmm::bicgstab(AS, X, B, PR, iter);

	if(iter.get_iteration() == iter.get_maxiter())
	{
		//the cheap ildlt precond didn't work, so let's try something a little stronger
		gmm::ildltt_precond< gmm::csc_matrix<double> > PR2(AS, 15, 1e-8);
		
		//reset the X
		for(int i=0;i<Fluid.M;i++)
		{
			for(int j=0;j<Fluid.N;j++)	
			{
				for(int k=0;k<Fluid.Q;k++)	
				{
					if(boundary(i,j,k) == -1)
					{
						const int row = mid(i,j,k);
						X[row] = Fluid.T(i,j,k);
					}
				}
			}
		}

		iter = gmm::iteration(1e-5);
		iter.set_maxiter(500);
		iter.set_noisy(0);

		gmm::bicgstab(AS, X, B, PR2, iter);

		if(iter.get_iteration() == iter.get_maxiter())
		{
			//still didn't work. throw a warning
			GapLog.message("\tInfo: CG Energy temperature loop failed to converge, using GS method instead.");

			//we're going to fall back to a GS method now
			return SlipperEnergy();
		}

	}

	//Copy the result back to the Fluid.T
	for(int i=0;i<Fluid.M;i++)
	{
		for(int j=0;j<Fluid.N;j++)	
		{
			for(int k=0;k<Fluid.Q;k++)	
			{
				if(boundary(i,j,k) == -1)
				{
					const int row = mid(i,j,k);
					Fluid.T(i,j,k) = X[row];
				}
			}
		}
	}

	//														//
	//Calculate a Qflux : heat transfer into the solid bodies
	//There are a number of methods to choose from
	//														//
	/*
	//Sum of viscous source (send only 0.5 to the slipper & 0.5 to the swashplate)
	Fluid.Qflux = 0.5*sum(bb(Range::all(), Range::all(), Range(1,Fluid.Q-2)),tensor::k) / Fluid.dA;
	Fluid.Qflux(0, Range::all()) = 0;
	Fluid.Qflux(Fluid.M-1, Range::all()) = 0;
	*/
	//"Proper" heat transfer (gradient poor when used on a large mesh)  (send only 0.5 to the slipper & 0.5 to the swashplate)
	
	{
		//Calculate the flux gradient using cells 0 & 1
		Range all = Range::all();
		Fluid.Qflux =	(	aB(all,all,1)*(Fluid.T(all, all, 1)-Fluid.T(all, all, 0))
								+
								aT(all,all,Fluid.Q-2)*(Fluid.T(all, all, Fluid.Q-2)-Fluid.T(all, all, Fluid.Q-1))
							)/Fluid.dA*0.5;

		//Calculate the flux gradient using cells 1 & 2
		Fluid.Qflux =	(	aB(all,all,2)*(Fluid.T(all, all, 2)-Fluid.T(all, all, 1))
								+
							-aT(all,all,Fluid.Q-3)*(Fluid.T(all, all, Fluid.Q-2)-Fluid.T(all, all, Fluid.Q-3))
							)/Fluid.dA*0.5;
	}
	
	/*
	//Advanced qflux calculation methods
	{
		Range mall = Range(1,Fluid.M-2);	//this is actually the non-bounary m
		Range north = Range(2,Fluid.M-1);
		Range south  = Range(0,Fluid.M-3);
		
		Array<double, 3> TEast(Fluid.M,Fluid.N,Fluid.Q);
		TEast(all,Range(0,Fluid.N-2),all) = Fluid.T(all,Range(1,Fluid.N-1),all);
		TEast(all,Fluid.N-1,all) = Fluid.T(all,0,all);
		
		Array<double, 3> TWest(Fluid.M,Fluid.N,Fluid.Q);
		TWest(all,Range(1,Fluid.N-1),all) = Fluid.T(all,Range(0,Fluid.N-2),all);
		TWest(all,0,all) = Fluid.T(all,Fluid.N-1,all);
		
		Fluid.Qflux(0,all) = 0;
		Fluid.Qflux(Fluid.M-1,all) = 0;
		
		
		////Consider only convective heat transfer in the x & y directions
		//Fluid.Qflux(mall,all) = sum(
		//							bb(mall,all,all)+(
		//							-FN(mall,all,all)*where(FN < 0, Fluid.T(north,all,all), Fluid.T(mall,all,all)) +
		//							FS(mall,all,all)*where(FS > 0, Fluid.T(south,all,all), Fluid.T(mall,all,all)) +
		//							-FE(mall,all,all)*where(FE < 0, TEast, Fluid.T(mall,all,all)) +
		//							FW(mall,all,all)*where(FW > 0, TWest, Fluid.T(mall,all,all))
		//							),tensor::k)*0.5/Fluid.dA(mall,all);
		//

		//Consider both diffusive and convective (per the power law) in the x & y directions
		Fluid.Qflux(mall,all) =	sum( bb(mall,all,all) +
									 aN(mall,all,all) * Fluid.T(north,all,all) +
									 aS(mall,all,all) * Fluid.T(south,all,all) +
									 aE(mall,all,all) * TEast(mall,all,all) +
									 aW(mall,all,all) * TWest(mall,all,all) -
									 (aP(mall,all,all)-aB(mall,all,all)-aT(mall,all,all)) * Fluid.T(mall,all,all)
									, tensor::k)/Fluid.dA(mall,all)*0.5;
	}
	*/
	
	return (int) iter.get_iteration();
}
