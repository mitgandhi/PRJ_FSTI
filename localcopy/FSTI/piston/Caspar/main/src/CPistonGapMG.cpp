#include "CPistonGap.h"
#include "logger.h"//debug dwm
#pragma once

extern class CPistonGap myPistonGap;
extern struct sGapResult myGapResult;
extern class CGapUtils myGapUtils;


//Set multigrid problem
void CPistonGap::PistonSetMG(void)
{
	Range all = Range::all();

	int nL,l;

	//finer mesh value equal to gap geometry
	levels[0].xyz = xyzf_gap;
	levels[0].phi = phi;
	levels[0].dx = dx;
	levels[0].dy = dy;
	levels[0].dphi = dphi;
	
	//lower levels
	l = 1;
	nL = levels[0].nL;
	while(l<nL)
	{
		//angular increment
		levels[l].dphi = 2.0 * PI / levels[l].N;
		//circumferential increment
		levels[l].dx = geometrypistongap.dK * PI / levels[l].N;
		//axial increment
		levels[l].dy = geometrypistongap.lvar / levels[l].M;
		//angular field
		for(int i=0;i<levels[l].N;i++)
		{
			levels[l].phi(Range(i*levels[l].M,(i+1)*levels[l].M-1)) = (i+0.5) * levels[l].dphi;
		};
		//field coordinates
		levels[l].xyz(all,0) = rK * cos(levels[l].phi);
		levels[l].xyz(all,1) = rK * sin(levels[l].phi);
		for(int i=0;i<levels[l].N;i++)
		{
			levels[l].xyz(Range(i*levels[l].M,(i+1)*levels[l].M-1),2) = (0.5 + tensor::i)*levels[l].dy;
		};
		l++;
	};

	//search neighbours for fine 2 coarse and coarse 2 fine interpolations
	int nn;
	double Ac,Af;
	Array<double,2> FaceDist;
	l = 0;
	nL = levels[0].nL;
	while(l<nL-1)
	{
		//coarse 2 fine - prolongation
		myGapUtils.SearchNeighbours(levels[l+1].xyz,levels[l].xyz,levels[l].FaceId_c2f,FaceDist,4);
		//fine 2 coarse - restriction
		Af = levels[l].dx * levels[l].dy;
		Ac = levels[l+1].dx * levels[l+1].dy;
		nn = (int) floor(Ac/Af);
		myGapUtils.SearchNeighbours(levels[l].xyz,levels[l+1].xyz,levels[l+1].FaceId_f2c,levels[l+1].FaceDist_f2c,nn);
		//new level
		l++;
	};

	//modify face neighbours for bilinear interpolation
	int MGInt = levels[0].MGInt;
	int id,nf;
	double xf,yf,xc,yc,dt,v1,v2,x,y,x1,y1,x2,y2,D;
	if(MGInt)
	{
		l=0;
		while(l<nL-1)
		{
			//number fo fine grid nodes
			nf = levels[l].FaceId_c2f.extent(0);
			//resize coefficient array
			levels[l].K_bilin.resize(nf,4);
			levels[l].K_bilin = 0.0;
			for(int i=0;i<nf;i++)
			{
				//face id
				id = levels[l].FaceId_c2f(i,0);
				//x coordinates fine mesh (length)
				xf = levels[l].xyz(i,2);
				//x coordinates coarse mesh (length)
				xc = levels[l+1].xyz(id,2);
				//y coordinate fine mesh (circumference)
				v1 = levels[l].xyz(i,0);
				v2 = levels[l].xyz(i,1);
				dt = atan2(v2,v1);
				if(dt<0){ dt+= 2*PI; };
				yf = rK * dt;
				//y coordinate coarse mesh (circumference)
				id = levels[l].FaceId_c2f(i,0);
				v1 = levels[l+1].xyz(id,0);
				v2 = levels[l+1].xyz(id,1);
				dt = atan2(v2,v1);
				if(dt<0){ dt+= 2*PI; };
				yc = rK * dt;
				//coarse face area and coordinates
				D = levels[l+1].dx * levels[l+1].dy;
				x1 = 0.0;	y1 = 0.0;
				x2 = levels[l+1].dy;	y2 = levels[l+1].dx;
				//check fine grid point position to define bilinear sequence on coarse grid
				Array<int,1> FaceId(4);
				FaceId=0;
				if(xf<=xc)
				{
					if(yf<=yc)
					{
						FaceId(0) = levels[l].FaceId_c2f(i,0) - levels[l+1].M - 1;
						FaceId(1) = levels[l].FaceId_c2f(i,0) - levels[l+1].M;
						FaceId(2) = levels[l].FaceId_c2f(i,0) - 1;
						FaceId(3) = levels[l].FaceId_c2f(i,0);
					}
					else
					{
						FaceId(0) = levels[l].FaceId_c2f(i,0) - 1;
						FaceId(1) = levels[l].FaceId_c2f(i,0);
						FaceId(2) = levels[l].FaceId_c2f(i,0) + levels[l+1].M - 1;
						FaceId(3) = levels[l].FaceId_c2f(i,0) + levels[l+1].M;
					}
					//correct boundary nodes to close circumference
					FaceId = where( FaceId<0, FaceId + levels[l+1].N*levels[l+1].M  , FaceId );
					FaceId = where( FaceId>=levels[l+1].N*levels[l+1].M, FaceId - levels[l+1].N*levels[l+1].M  , FaceId );
					//correct boundary nodes south boundary (DC side)
					if(xf<=0.5*levels[l+1].dy)
					{
						FaceId(0) = -1;		FaceId(2) = -1;
						D *= 0.5;
						x2 *= 0.5;
					}
				}
				else
				{
					if(yf<=yc)
					{
						FaceId(0) = levels[l].FaceId_c2f(i,0) - levels[l+1].M;
						FaceId(1) = levels[l].FaceId_c2f(i,0) - levels[l+1].M + 1;
						FaceId(2) = levels[l].FaceId_c2f(i,0);
						FaceId(3) = levels[l].FaceId_c2f(i,0) + 1;
					}
					else
					{
						FaceId(0) = levels[l].FaceId_c2f(i,0);
						FaceId(1) = levels[l].FaceId_c2f(i,0) + 1;
						FaceId(2) = levels[l].FaceId_c2f(i,0) + levels[l+1].M;
						FaceId(3) = levels[l].FaceId_c2f(i,0) + levels[l+1].M + 1;
					}
					//correct boundary nodes to close circumference
					FaceId = where( FaceId<0, FaceId + levels[l+1].N*levels[l+1].M  , FaceId );
					FaceId = where( FaceId>=levels[l+1].N*levels[l+1].M, FaceId - levels[l+1].N*levels[l+1].M  , FaceId );
					//correct boundary nodes north boundary (case side)
					if(xf>=(levels[l+1].M-0.5)*levels[l+1].dy)
					{
						FaceId(1) = -2;		FaceId(3) = -2;
						D *= 0.5;
						x2 *= 0.5;
					}
				}
				//assign face ids
				levels[l].FaceId_c2f(i,all) = FaceId(all);
				//x coordinate fine mesh (circumference)
				v1 = levels[l].xyz(i,0);
				v2 = levels[l].xyz(i,1);
				//angle
				dt = atan2(v2,v1);
				if(dt<0){ dt+= 2*PI; };
				yf = rK * dt;
				//x coordinate coarse mesh (circumference)
				if(FaceId(0)==-1)
				{
					v1 = levels[l+1].xyz(FaceId(1),0);
					v2 = levels[l+1].xyz(FaceId(1),1);
				}
				else
				{
					v1 = levels[l+1].xyz(FaceId(0),0);
					v2 = levels[l+1].xyz(FaceId(0),1);
				}
				//angle
				dt = atan2(v2,v1);
				if(dt<0){ dt+= 2*PI; };
				yc = rK * dt;
				//y coordinates (length)
				xf = levels[l].xyz(i,2);
				if(FaceId(0)==-1)
				{
					xc = 0.0;
				}
				else
				{
					xc = levels[l+1].xyz(FaceId(0),2);
				}
				//x&y fine relative locations
				x = xf - xc;
				if(yc>=yf)
				{
					yc -= 2.0*PI*rK;
				};
				y = yf - yc;
				//coefficients for bilinear interpolation
				levels[l].K_bilin(i,0) = (x2-x) * (y2-y) / D;
				levels[l].K_bilin(i,1) = (x-x1) * (y2-y) / D;
				levels[l].K_bilin(i,2) = (x2-x) * (y-y1) / D;
				levels[l].K_bilin(i,3) = (x-x1) * (y-y1) / D;
			};
			//new level
			l++;
		};
	};

};
//Initilaize Reynolds diffusive coefficients for MG levels
void CPistonGap::PistonReynoldsCalcCoefficientsMG(double dt)
{
	double vK,omega,speedK,pDC,pCase;
	int nL;
	Range all = Range::all( );

	//Geometry
	vK = operatingpistongap.vK;	
	omega = operatingpistongap.omega;		
	speedK = operatingpistongap.speedK;

	//Initializing the pressure 2D vector 
	pDC = operatingpistongap.pDC;
	pCase = operatingpistongap.pCase;

	//Fluid average viscosity
	mu = 0;
	rho2d.resize(M*N);
	rho2d = 0;
	for(int i=0;i<Q;i++)
	{
		mu(Range(0,N*M-1)) += oilviscosity(Range(i*N*M,(i+1)*N*M-1));
		rho2d(Range(0,N*M-1)) += oildensity(Range(i*N*M,(i+1)*N*M-1));
	};
	mu/=Q;
	rho2d /= Q;

	//Number of levels
	nL = levels[0].nL;
	//Fluid viscosity level 0
	levels[0].mu = mu;

	levels[0].rho = rho2d;
	//Film thickness level 0
	levels[0].h = h;
	//Assign first value to field
	levels[0].x = pnew;
	//Initial value of residual
	levels[0].Rx = 1.0;


	//Diffusive coefficients
	int l=0;
	while(l<nL)
	{
		//Reset variables
		levels[l].b = 0.0; levels[l].r = 0.0; levels[l].e = 0.0; levels[l].Rx = 1.0; 
		//Interpolate viscosity and film thickness to lower levels
		if(l!=nL-1)
		{
			//Viscosity
			levels[l+1].mu = PistonRestrictionMG(levels[l].mu,levels[l+1].FaceId_f2c,levels[l+1].FaceDist_f2c);
			levels[l+1].rho = PistonRestrictionMG(levels[l].rho,levels[l+1].FaceId_f2c,levels[l+1].FaceDist_f2c);
			//Film thickness
			levels[l+1].h = PistonRestrictionMG(levels[l].h,levels[l+1].FaceId_f2c,levels[l+1].FaceDist_f2c);
		}

		//-------------Diffusive coefficients calculations-------------//
		double mun,mus,mue,muw,hKn,hKs,hKe,hKw,hn,hs,he,hw,hp,defsqueezeK,defsqueezeB,expansion,rhon,rhos,rhoe,rhow,rhop;
		double nB,sB,xn,xs,xe,xw;
		nB=0.0;
		sB=0.0;
		
		//Grid size
		int N = levels[l].N;
		int M = levels[l].M;
		double dx = levels[l].dx;
		double dy = levels[l].dy;
		
		//Face diffusivity values
		for(int i=0;i<N*M;i++)
		{

			rhop = levels[l].rho(i);
			hp = levels[l].h(i);
			//North
			if(i%M==(M-1))
			{
				mun = levels[l].mu(i);
				hn = levels[l].h(i);
				rhon = levels[l].rho(i);
				if(l==0){
					hKn = hK(i);
				}
			}
			else
			{
				mun = 0.5*( levels[l].mu(i) + levels[l].mu(i+1) );
				hn = 0.5*( levels[l].h(i) + levels[l].h(i+1) );
				rhon = 0.5*( levels[l].rho(i) + levels[l].rho(i+1) );
				if(l==0){
					hKn = 0.5*( hK(i) + hK(i+1) );
				}
			}
			//South
			if(i%M==0)
			{
				mus = levels[l].mu(i);
				hs = levels[l].h(i);
				rhos = levels[l].rho(i);
				if(l==0){
					hKs = hK(i);
				}
			}
			else
			{
				mus = 0.5*( levels[l].mu(i) + levels[l].mu(i-1) );
				hs = 0.5*( levels[l].h(i) + levels[l].h(i-1) );
				rhos = 0.5*( levels[l].rho(i) + levels[l].rho(i-1) );
				if(l==0){
					hKs = 0.5*( hK(i) + hK(i-1) );
				}
			}
			//East
			if(i>=(N-1)*M)
			{
				mue = 0.5*( levels[l].mu(i) + levels[l].mu(i-(N-1)*M) );
				he = 0.5*( levels[l].h(i) + levels[l].h(i-(N-1)*M) );
				rhoe = 0.5*( levels[l].rho(i) + levels[l].rho(i-(N-1)*M) );
				if(l==0){
					hKe = 0.5*( hK(i) + hK(i-(N-1)*M) );
				}
			}
			else
			{
				mue = 0.5*( levels[l].mu(i) + levels[l].mu(i+M) );
				he = 0.5*( levels[l].h(i) + levels[l].h(i+M) );
				rhoe = 0.5*( levels[l].rho(i) + levels[l].rho(i+M) );
				if(l==0){
					hKe = 0.5*( hK(i) + hK(i+M) );
				}
			}
			//West
			if(i<M)
			{
				muw = 0.5*( levels[l].mu(i) + levels[l].mu(i+(N-1)*M) );
				hw = 0.5*( levels[l].h(i) + levels[l].h(i+(N-1)*M) );
				rhow = 0.5*( levels[l].rho(i) + levels[l].rho(i+(N-1)*M) );
				if(l==0){
					hKw = 0.5*( hK(i) + hK(i+(N-1)*M) );
				}
			}
			else
			{
				muw = 0.5*( levels[l].mu(i) + levels[l].mu(i-M) );
				hw = 0.5*( levels[l].h(i) + levels[l].h(i-M) );
				rhow = 0.5*( levels[l].rho(i) + levels[l].rho(i-M) );
				if(l==0){
					hKw = 0.5*( hK(i) + hK(i-M) );
				}
			}

			//North
			levels[l].an(i) = rhon * pow(hn,3.0) * dx / ( dy * 6.0 * mun );
			if(i%M==(M-1) && i>0){
				levels[l].an(i) *= 2.0;	};
			//South
			levels[l].as(i) = rhos * pow(hs,3.0) * dx / ( dy * 6.0 * mus );
			if(i%M==0){
				levels[l].as(i) *= 2.0;	};
			//East
			levels[l].ae(i) = rhoe * pow(he,3.0) * dy / ( dx * 6.0 * mue );
			//West
			levels[l].aw(i) = rhow * pow(hw,3.0) * dy / ( dx * 6.0 * muw );
			//Point
			levels[l].ap(i) = levels[l].an(i) + levels[l].as(i) + levels[l].ae(i) + levels[l].aw(i);
			//Source level 0
			if(l==0)
			{
				//Gradients sliding part
				double dhKx = (hKe - hKw)/dx;
				double dhKy = (hKn - hKs)/dy;
				//Squeeze due to deformation change
				defsqueezeK = (CPistonGap::defK_p_gap(i) - CPistonGap::defK_p_gap_squeeze(i))/dt;
				defsqueezeB = (CPistonGap::defB_p_gap(i) - CPistonGap::defB_p_gap_squeeze(i))/dt;
				expansion = (rho2d(i) - CPistonGap::density_expansion(i))/dt;
				//Source
				if(hp > geometrypistongap.hmin){
					levels[l].b(i) = ( - speedK * omega * rK * (rhoe*he - rhow*hw) / dx - vK * (rhon*hn - rhos*hs) / dy 
						+ 2.0 * (rhop * dht(i) + rhop * speedK * omega * rK * dhKx + rhop * vK * dhKy + rhop * defsqueezeK - rhop * defsqueezeB - hp * expansion ) ) * dx * dy;
				}
				else
					levels[l].b(i) =  ( - speedK * omega * rK * hp * (rhoe - rhow) / dx - vK * hp * (rhon - rhos) / dy 
						- 2.0 * ( hp * expansion ) ) * dx * dy;
				//Dirichlet
				if(i%M==(M-1))
				{
					levels[l].b(i) += levels[l].an(i)*pCase;
				};
				if(i%M==0)
				{
					levels[l].b(i) += levels[l].as(i)*pDC;
				};
				//Residual
				if(i%M==(M-1))
				{
					xn = nB;
				}
				else
				{
					xn = levels[l].x(i+1);
				}
				if(i%M==0)
				{
					xs = sB;
				}
				else
				{
					xs = levels[l].x(i-1);
				}
				if(i>=(N-1)*M)
				{
					xe = levels[l].x(i-(N-1)*M);
				}
				else
				{
					xe = levels[l].x(i+M);
				}
				if(i<M)
				{
					xw = levels[l].x(i+(N-1)*M);
				}
				else
				{
					xw = levels[l].x(i-M);
				}
				levels[l].r(i) = levels[l].b(i)  - ( levels[l].ap(i) * levels[l].x(i) - levels[l].an(i) * xn - levels[l].as(i) * xs - levels[l].ae(i) * xe - levels[l].aw(i) * xw  );
			};
		};

		l++;

	};
                                            
};
//Solve Reynolds using MG: V or W cycle
void CPistonGap::PistonReynoldsMG(void)
{
	int VW;

	//Cycle type
	VW = levels[0].VW;

	//Solve with V or W cycle
	if(VW)
	{
		PistonWcycleMG();
	}
	else
	{
		PistonVcycleMG();
	}

	//Assign and limit pressure
	pnew = levels[0].x;

	//Limit pressure
	if(PressureDeformation)
	{
		p = where(pnew>1.0e9,1.0e9,pnew);
	}
	else
	{
		p = where(pnew>1.0e9,1.0e9,pnew);
	};
	p = where(p<1.0e4,1.0e4,p);


	/*//Test MG
	if(VW)
	{
		fout.open("./outputs_piston/pW.dat");
		for(int i = 0; i < N ; i++)
		{
			for(int j = 0; j < M ; j++)
			{
				fout << scientific << pnew(j+i*M)*1.0e-5 << "\t" ;
			}
			fout << "\n";
		}
		fout.close();
		fout.clear();
		system("PAUSE");
	}
	else
	{
		fout.open("./outputs_piston/pV.dat");
		for(int i = 0; i < N ; i++)
		{
			for(int j = 0; j < M ; j++)
			{
				fout << scientific << pnew(j+i*M)*1.0e-5 << "\t" ;
			}
			fout << "\n";
		}
		fout.close();
		fout.clear();
		system("PAUSE");
	}*/

};
//GS SOR method in levels
void CPistonGap::PistonGaussSeidelMG(int l,int leg)
{
	double alpha,sB,nB,xn,xs,xe,xw;
	int iterations,v;

	alpha = 1.8;
	iterations = 0;
	int N = levels[l].N;
	int M = levels[l].M;
	int v1 = levels[l].v1;
	int v2 = levels[l].v2;

	//Boundaries Reynolds
	nB = 0.0;	sB = 0.0;

	//Number of sweeps
	if(leg)
	{
		v = (l+2) * v2;
		//Coarse grid correction
		levels[l].x += levels[l].e;
	}
	else
	{
		v = (l+2) * v1;
	}

	do
	{
		//SOR
		for(int i=0;i<N*M;i++)
		{
			//North
			if(i%M==(M-1) && i>0)
			{
				xn = nB;
			}
			else
			{
				xn = levels[l].x(i+1);
			}
			//South
			if(i%M==0)
			{
				xs = sB;
			}
			else
			{
				xs = levels[l].x(i-1);
			}
			//East
			if(i>=(N-1)*M)
			{
				xe = levels[l].x(i-(N-1)*M);
			}
			else
			{
				xe = levels[l].x(i+M);
			}
			//West
			if(i<M)
			{
				xw = levels[l].x(i+(N-1)*M);
			}
			else
			{
				xw = levels[l].x(i-M);
			}
			//Point value
			levels[l].x(i) +=  alpha * ( ( ( levels[l].an(i) * xn + levels[l].as(i) * xs + levels[l].ae(i) * xe + levels[l].aw(i) * xw + levels[l].b(i) ) / levels[l].ap(i) ) - levels[l].x(i) );
			if(l == 0)
				levels[0].x(i) = levels[0].x(i)<=1e4? 1e4:levels[0].x(i);
		}
		//Counter
		iterations++;
	}while(iterations < v);

	//Calculate residual
	if(leg==0 ||l==0)
	{
		for(int i=0;i<N*M;i++)
		{
			//North
			if(i%M==(M-1) && i>0)
			{
				xn = nB;
			}
			else
			{
				xn = levels[l].x(i+1);
			}
			//South
			if(i%M==0)
			{
				xs = sB;
			}
			else
			{
				xs = levels[l].x(i-1);
			}
			//East
			if(i>=(N-1)*M)
			{
				xe = levels[l].x(i-(N-1)*M);
			}
			else
			{
				xe = levels[l].x(i+M);
			}
			//West
			if(i<M)
			{
				xw = levels[l].x(i+(N-1)*M);
			}
			else
			{
				xw = levels[l].x(i-M);
			}
			//Residual
			levels[l].r(i) = levels[l].b(i)  - ( levels[l].ap(i) * levels[l].x(i) - levels[l].an(i) * xn - levels[l].as(i) * xs - levels[l].ae(i) * xe - levels[l].aw(i) * xw  );
		}
	};

};
//Restrict from finer to coarser level
Array<double,1> CPistonGap::PistonRestrictionMG(Array<double,1> field_f,Array<int,2> FaceId_f2c,Array<double,2> FaceDist_f2c)
{
	int nc = (int) FaceId_f2c.extent(0);
	int nn = (int) FaceId_f2c.extent(1);
	double dn,wi;


	Array<double,1> field_c(nc);
	field_c = 0.0;

	//weigthed average over neighbouring nodes distances 
	for(int i=0;i<nc;i++)
	{
		dn = 0.0;
		wi = 0.0;
		for(int j=0;j<nn;j++)
		{
			wi = 1.0 / FaceDist_f2c(i,j) ;
			field_c(i) += field_f(FaceId_f2c(i,j)) * wi;
			dn += wi ;
		}
		field_c(i) /= dn;
	}

	return field_c;


};
//Interpolate from coarser to finer level
Array<double,1> CPistonGap::PistonProlongationMG(Array<double,1> field_c,Array<int,2> FaceId_c2f)
{

	int nf = (int) FaceId_c2f.extent(0);
	Array<double,1> field_f(nf);
	field_f=0.0;

	//nearest neighbour interpolation
	for(int i=0;i<nf;i++)
	{
		field_f(i) = field_c(FaceId_c2f(i,0));
	}
		
	return field_f;

};
Array<double,1> CPistonGap::PistonProlongationBilinearMG(int l_c,Array<double,1> field_c,Array<int,2> FaceId_c2f)
{

	int nf,id11,id12,id21,id22;
	double nB,sB,Q11,Q12,Q21,Q22;

	//north and south boundaries
	nB = 0.0;	sB = 0.0;

	//fine field initialization
	nf = FaceId_c2f.extent(0);
	Array<double,1> field_f(nf);
	field_f=0.0;
	//Interpolation scheme
	for(int i=0;i<nf;i++)
	{
		//coarse field node ids
		id11 = FaceId_c2f(i,0);
		id21 = FaceId_c2f(i,1);
		id12 = FaceId_c2f(i,2);
		id22 = FaceId_c2f(i,3);
		//coarse field values
		if(id11==-1)
		{
			Q11 = sB;	Q12 = sB;
		}
		else
		{
			Q11 = field_c(id11);	Q12 = field_c(id12);
		}
		if(id21==-2)
		{
			Q21 = nB;	Q22 = nB;
		}
		else
		{
			Q21 = field_c(id21);	Q22 = field_c(id22);
		}
		//bilinear interpolation
		field_f(i) = Q11 * levels[l_c-1].K_bilin(i,0) + Q21 * levels[l_c-1].K_bilin(i,1) + Q12 * levels[l_c-1].K_bilin(i,2) + Q22 * levels[l_c-1].K_bilin(i,3);
	};

	return field_f;

};
//W cycle routine
void CPistonGap::PistonWcycleMG(void)
{
	int l,nL,iterations,MGInt,v1,v2,iter_max;
	double cost;
	pcon.resize(levels[0].x.size()); pcon = levels[0].x;

	//Number of levels
	nL = levels[0].nL;

	//Number of sweeps
	v1 = levels[0].v1;
	v2 = levels[0].v2;

	//Interpolation type
	MGInt = levels[0].MGInt;

	//Initialize
	levels[0].Rx = 1.0;
	iterations = 0;
	cost = 0.0;
	iter_max = 10000;

	/*//Test MG
	clock_t start,end,ticks;
	start = clock();*/
	do
	{
		/*//Test MG
		fout.open("./outputs_piston/RW.dat",ios::app);
		fout << levels[0].Rx << "\t" << cost << "\n";
		fout.close();*/

		//Level 0
		l=0;

		//W-cycle Down Leg
		do
		{
			//Solve on finer level
			PistonGaussSeidelMG(l,0);
			//Interpolate residual to coarser level
			if(l!=nL-1)
			{
				levels[l+1].b = PistonRestrictionMG(levels[l].r,levels[l+1].FaceId_f2c,levels[l+1].FaceDist_f2c);
			}
			//Increment level
			l++;
			/*//Test MG - Cost
			cost += v1 * (l+1) * pow(2.0,-2.0*l);*/
		}while(l<nL);

		//W-cycle Up Leg
		l=nL-1;
		do
		{
			//Interpolate error to finer level
			if(MGInt){
				levels[l-1].e = PistonProlongationBilinearMG(l,levels[l].x,levels[l-1].FaceId_c2f);
			}
			else{
				levels[l-1].e = PistonProlongationMG(levels[l].x,levels[l-1].FaceId_c2f);
			};
			//Solve on finer level
			PistonGaussSeidelMG(l-1,1);
			//Decrement level
			l--;
			/*//Test MG - Reset coarser field value
			levels[l].x = 0.0;
			//Test MG - Cost
			cost += v2 * ((l-1)+1) * pow(2.0,-2.0*(l-1));*/
		}while(l>=(int) nL/3);

		//W-cycle Down Leg
		l=(int) nL/3;
		do
		{
			//Solve on finer level
			PistonGaussSeidelMG(l,0);
			//Interpolate residual to coarser level
			if(l!=nL-1)
			{
				levels[l+1].b = PistonRestrictionMG(levels[l].r,levels[l+1].FaceId_f2c,levels[l+1].FaceDist_f2c);
			}
			//Increment level
			l++;
			/*//Test MG - Cost
			cost += v1 * (l+1) * pow(2.0,-2.0*l);*/
		}while(l<nL);

		//W-cycle Up Leg
		l=nL-1;
		do
		{
			//Interpolate error to finer level
			if(MGInt){
				levels[l-1].e = PistonProlongationBilinearMG(l,levels[l].x,levels[l-1].FaceId_c2f);
			}
			else{
				levels[l-1].e = PistonProlongationMG(levels[l].x,levels[l-1].FaceId_c2f);
			};
			//Solve on finer level
			PistonGaussSeidelMG(l-1,1);
			//Decrement level
			l--;
			/*//Test MG - Reset coarser field value
			levels[l].x = 0.0;
			//Test MG - Cost
			cost += v2 * ((l-1)+1) * pow(2.0,-2.0*(l-1));*/
		}while(l>=(int) 2*nL/3);

		//W-cycle Down Leg
		l=(int) 2*nL/3;
		do
		{
			//Solve on finer level
			PistonGaussSeidelMG(l,0);
			//Interpolate residual to coarser level
			if(l!=nL-1)
			{
				levels[l+1].b = PistonRestrictionMG(levels[l].r,levels[l+1].FaceId_f2c,levels[l+1].FaceDist_f2c);
			}
			//Increment level
			l++;
			/*//Test MG - Cost
			cost += v1 * (l+1) * pow(2.0,-2.0*l);*/
		}while(l<nL);

		//W-cycle Up Leg
		l=nL-1;
		do
		{
			//Interpolate error to finer level
			if(MGInt){
				levels[l-1].e = PistonProlongationBilinearMG(l,levels[l].x,levels[l-1].FaceId_c2f);
			}
			else{
				levels[l-1].e = PistonProlongationMG(levels[l].x,levels[l-1].FaceId_c2f);
			};
			//Solve on finer level
			PistonGaussSeidelMG(l-1,1);
			//Decrement level
			l--;
			/*//Test MG - Reset coarser field value
			levels[l].x = 0.0;
			//Test MG - Cost
			cost += v2 * ((l-1)+1) * pow(2.0,-2.0*(l-1));*/
		}while(l>=(int) nL/3);

		//W-cycle Down Leg
		l=(int) nL/3;
		do
		{
			//Solve on finer level
			PistonGaussSeidelMG(l,0);
			//Interpolate residual to coarser level
			if(l!=nL-1)
			{
				levels[l+1].b = PistonRestrictionMG(levels[l].r,levels[l+1].FaceId_f2c,levels[l+1].FaceDist_f2c);
			}
			//Increment level
			l++;
			/*//Test MG - Cost
			cost += v1 * (l+1) * pow(2.0,-2.0*l);*/
		}while(l<nL);

		//W-cycle Up Leg
		l=nL-1;
		do
		{
			//Interpolate error to finer level
			if(MGInt){
				levels[l-1].e = PistonProlongationBilinearMG(l,levels[l].x,levels[l-1].FaceId_c2f);
			}
			else{
				levels[l-1].e = PistonProlongationMG(levels[l].x,levels[l-1].FaceId_c2f);
			};
			//Solve on finer level
			PistonGaussSeidelMG(l-1,1);
			//Decrement level
			l--;
			/*//Test MG - Reset coarser field value
			levels[l].x = 0.0;
			//Test MG - Cost
			cost += v2 * ((l-1)+1) * pow(2.0,-2.0*(l-1));*/
		}while(l>=1);

		//Counter
		iterations++;

		//Scaled residual
		PistonReynoldsCalcResidualMG(0);

		//Residual in MPa
	}while( levels[0].Rx > (Rmin_R * 1e6) && iterations < iter_max );

	//Log << "MG Iterations: " << iterations << "\n";


	//Log file writing if max number of iterations is reached
	if(iterations>=iter_max)
	{
		fout.open ("./output/piston/matlab/PistonGapLoopLog.txt",ios::app);
		fout << "Max Number of Iterations Reached in Reynolds Equation Loop" << "\n";
		fout << "Shaft Angle: " << "\t" << operatingpistongap.phi_rad*180/PI << "\n";
		fout << "Residual: " << "\t" << levels[0].Rx << "\n";
		fout << " " << "\n";
		fout.close();
		fout.clear();
	}

	/*//Test MG
	end = clock();
	ticks = end - start;
	double time = ticks / (double) CLOCKS_PER_SEC;
	fout.open("./outputs_piston/nW.dat");
	fout << iterations << "\n";
	fout.close();
	fout.open("./outputs_piston/tW.dat");
	fout << time << "\n";
	fout.close();*/
                                                         
};
//V cycle routine
void CPistonGap::PistonVcycleMG(void)
{
	int l,nL,iterations,MGInt,v1,v2,iter_max;
	double cost;
	pcon.resize(levels[0].x.size()); pcon = levels[0].x;

	//Number of levels
	nL = levels[0].nL;

	//Number of sweeps
	v1 = levels[0].v1;
	v2 = levels[0].v2;

	//Interpolation type
	MGInt = levels[0].MGInt;

	//Initialize
	levels[0].Rx = 1.0;
	iterations = 0;
	cost = 1.0;
	iter_max = 10000;

	/*//Test MG
	clock_t start,end,ticks;
	start = clock();*/
	do
	{
		/*//Test MG
		fout.open("./outputs_piston/RV.dat",ios::app);
		fout << levels[0].Rx << "\t" << cost << "\n";
		fout.close();*/

		//V-cycle Down Leg
		l=0;
		do
		{
			//Solve on finer level
			PistonGaussSeidelMG(l,0);
			//Interpolate residual to coarser level
			if(l!=nL-1)
			{
				levels[l+1].b = PistonRestrictionMG(levels[l].r,levels[l+1].FaceId_f2c,levels[l+1].FaceDist_f2c);
			}
			//Increment level
			l++;
			/*//Test MG - Cost
			cost += v1 * (l+1) * pow(2.0,-2.0*l);*/
		}while(l<nL);

		//V-cycle Up Leg
		l=nL-1;
		do
		{
			//Interpolate error to finer level
			if(MGInt){
				levels[l-1].e = PistonProlongationBilinearMG(l,levels[l].x,levels[l-1].FaceId_c2f);
			}
			else{
				levels[l-1].e = PistonProlongationMG(levels[l].x,levels[l-1].FaceId_c2f);
			};
			//Solve on finer level
			PistonGaussSeidelMG(l-1,1);
			//Decrement level
			l--;
			/*//Test MG - Reset coarser field value
			levels[l].x = 0.0;
			//Test MG - Cost
			cost += v2 * ((l-1)+1) * pow(2.0,-2.0*(l-1));*/
		}while(l>=1);

		//Scaled residual
		PistonReynoldsCalcResidualMG(0);

		//Counter
		iterations++;

		//Residual in MPa
	}while( levels[0].Rx > (Rmin_R * 1e6) && iterations < iter_max );


	//Log file writing if max number of iterations is reached
	if(iterations>=iter_max)
	{
		fout.open("./output/piston/matlab/PistonGapLoopLog.txt",ios::app);
		fout << "Max Number of Iterations Reached in Reynolds Equation Loop" << "\n";
		fout << "Shaft Angle: " << "\t" << operatingpistongap.phi_rad*180/PI << "\n";
		fout << "Residual: " << "\t" << levels[0].Rx << "\n";
		fout << " " << "\n";
		fout.close();
		fout.clear();
	}


	/*//Test MG
	end = clock();
	ticks = end - start;
	double time = ticks / (double) CLOCKS_PER_SEC;

	fout.open("./outputs_piston/nV.dat");
	fout << iterations << "\n";
	fout.close();

	fout.open("./outputs_piston/tV.dat");
	fout << time << "\n";
	fout.close();*/


};
//Calculate Reynolds residual for MG solution
void CPistonGap::PistonReynoldsCalcResidualMG(int l)
{
	/*if(!pcon.size()){
		pcon.resize(levels[0].x.size());
		pcon = 0.0;
	}*/

	
	levels[l].Rx = max(fabs(pcon - levels[0].x));

	pcon = levels[0].x;
	

	//Log << levels[l].Rx << "\n";

	//Previously: dwm
	/*double n,d;
		
	//Calculate scaled residual for each level
	levels[l].Rx = 0.0;
	n = sum(fabs(levels[l].r));
	d = sum(fabs(levels[l].ap*levels[l].x));
	//Test MG
	n = sqrt(sum(pow(levels[l].r,2.0)));
	d = sqrt(sum(pow(levels[l].b,2.0)));
	levels[l].Rx = n / d;*/
	
};

