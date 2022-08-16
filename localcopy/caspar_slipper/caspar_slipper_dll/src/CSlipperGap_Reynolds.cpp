#include "CSlipperGap.h"
#include <omp.h>
#include <time.h>
#pragma once
#define OMP

extern void matlab(const Array<double,2>& data,const string file);

double CSlipperGap::presid()
{

	int M = Fluid.M;
	int N = Fluid.N;

	Array<double,2> mu(M,N);		mu = 0;			//Actual visco
	
	Array<double,2> hn(M,N);		hn = 0;			//Gap height, north face
	Array<double,2> hs(M,N);		hs = 0;			//Gap height, south face
	Array<double,2> he(M,N);		he = 0;			//Gap height, east face
	Array<double,2> hw(M,N);		hw = 0;			//Gap height, west face

	Array<double,2> htn(M,N);		htn = 0;			//Gap height top, north face
	Array<double,2> hts(M,N);		hts = 0;			//Gap height top, south face
	Array<double,2> hte(M,N);		hte = 0;			//Gap height top, east face
	Array<double,2> htw(M,N);		htw = 0;			//Gap height top, west face

	Array<double,2> hbn(M,N);		hbn = 0;			//Gap height bottom, north face
	Array<double,2> hbs(M,N);		hbs = 0;			//Gap height bottom, south face
	Array<double,2> hbe(M,N);		hbe = 0;			//Gap height bottom, east face
	Array<double,2> hbw(M,N);		hbw = 0;			//Gap height bottom, west face

	Array<double,2> mun(M,N);		mun = 0;			//Visco, north face
	Array<double,2> mus(M,N);		mus = 0;			//Visco, south face
	Array<double,2> mue(M,N);		mue = 0;			//Visco, east face
	Array<double,2> muw(M,N);		muw = 0;			//Visco, west face

	Array<double,2> an(M,N);		an = 0;			//Finite volume coefficent an
	Array<double,2> as(M,N);		as = 0;			//Finite volume coeffcient as
	Array<double,2> ae(M,N);		ae = 0;			//Finite volume coefficent ae
	Array<double,2> aw(M,N);		aw = 0;			//Finite volume coefficent aw
	Array<double,2> ap(M,N);		ap = 0;			//Finite volume coefficient ap
	Array<double,2> b(M,N);			b = 0;			//Reynolds source term
	
	//Linaer interpolate cell height to cell faces
	cell2face(Fluid.h, he, hw, hn, hs);

	//Average viscosity over gap height
	mu = mean(Fluid.oilviscosity,tensor::k);

	//Linaer interpolate cell visco to cell faces
	cell2face(mu, mue, muw, mun, mus);

	//Linaer interpolate top gap height to cell faces
	Array<double,2> hT (Fluid.h + swashplate.ehd);
	cell2face(hT, hte, htw, htn, hts);
	
	//Linaer interpolate bottom gap height to cell faces
	Array<double,2> hB (swashplate.ehd.copy());
	cell2face(hB, hbe, hbw, hbn, hbs);

	//
	Array<double,2> vr((operatingslippergap.tvr+operatingslippergap.bvr)/2.0);
	Array<double,2> vrn(M,N);
	Array<double,2> vrs(M,N);
	Array<double,2> vre(M,N);
	Array<double,2> vrw(M,N);
	cell2face(vr, vre, vrw, vrn, vrs);

	Array<double,2> vtheta((operatingslippergap.tvtheta+operatingslippergap.bvtheta)/2.0);
	Array<double,2> vthetan(M,N);
	Array<double,2> vthetas(M,N);
	Array<double,2> vthetae(M,N);
	Array<double,2> vthetaw(M,N);
	cell2face(vtheta, vthetae, vthetaw, vthetan, vthetas);

	//Spatial coefficients
	if (gapinput->options_slipper.general.reynolds_mu == 2)
	{
		//Full Mu
		ae = (pow(he,3.0)*Fluid.dr)/(mue*Fluid.r * (Fluid.dtheta+East(Fluid.dtheta))/2.0 );
		aw = (pow(hw,3.0)*Fluid.dr)/(muw*Fluid.r * (Fluid.dtheta+West(Fluid.dtheta))/2.0 );
		an = (pow(hn,3.0)*(Fluid.r+Fluid.dr/2.0)*Fluid.dtheta)/(mun * (Fluid.dr+North(Fluid.dr))/2.0 );
		as = (pow(hs,3.0)*(Fluid.r-Fluid.dr/2.0)*Fluid.dtheta)/(mus * (Fluid.dr+South(Fluid.dr))/2.0 );
		ap = an + as + ae + aw;
	} else {
		//Averaged mu
		ae = (pow(he,3.0)*Fluid.dr)/(mu*Fluid.r * (Fluid.dtheta+East(Fluid.dtheta))/2.0 );
		aw = (pow(hw,3.0)*Fluid.dr)/(mu*Fluid.r * (Fluid.dtheta+West(Fluid.dtheta))/2.0 );
		an = (pow(hn,3.0)*(Fluid.r+Fluid.dr/2.0)*Fluid.dtheta)/(mu * (Fluid.dr+North(Fluid.dr))/2.0 );
		as = (pow(hs,3.0)*(Fluid.r-Fluid.dr/2.0)*Fluid.dtheta)/(mu * (Fluid.dr+South(Fluid.dr))/2.0 );
		ap = an + as + ae + aw;
	}

	//Source Term
	b =	(
				-( (operatingslippergap.tvr+operatingslippergap.bvr)/2.0 * (hn-hs)/Fluid.dr )
				-( (operatingslippergap.tvtheta+operatingslippergap.bvtheta)/2.0 * (he-hw)/(Fluid.r*Fluid.dtheta) )
				-( (operatingslippergap.tvr+operatingslippergap.bvr)/2.0 * Fluid.h/Fluid.r)
				-( Fluid.h * (vrn-vrs)/Fluid.dr )
				-( Fluid.h * (vthetae-vthetaw)/(Fluid.r*Fluid.dtheta) )
				+( operatingslippergap.tvr * (htn-hts)/Fluid.dr )
				+( operatingslippergap.tvtheta * (hte-htw)/(Fluid.r*Fluid.dtheta) )
				-( operatingslippergap.bvr * (hbn-hbs)/Fluid.dr )
				-( operatingslippergap.bvtheta * (hbe-hbw)/(Fluid.r*Fluid.dtheta) )
				-( Fluid.dht )
			)*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	//b = where(Fluid.contact == 0, b, 0);

	if(curFullLoop <= 0 && gapinput->options_slipper.general.EHDsqueeze && gapinput->options_slipper.general.SlipperPressureDeformation)
	{
		//double IMval = 9.07295e-008/100e+5;
		//Array<double, 2> sqzfactor(IMval/gapinput->common.timestep*12*Fluid.r*Fluid.dr*Fluid.dtheta);
		Array<double, 2> sqzfactor(Fluid.newt_ehdsqzapprox(tensor::i)/gapinput->common.timestep*12*Fluid.r*Fluid.dr*Fluid.dtheta);
		ap += sqzfactor;
		b += Fluid.pFullLoop*sqzfactor;
	}

	//Calculate the initial residual L2 norm
	{
		double num = 0;
		double dom = 0;
		
		for(int i=1;i<M-1;i++)
		{
			for(int j=0;j<N;j++)	
			{
				
				if(Fluid.boundary(i,j) == -1 && 
					Fluid.p_uncut(i,j) == Fluid.p(i,j) &&
					Fluid.p_uncut(i,(N+j+1)%N) == Fluid.p(i,(N+j+1)%N) &&
					Fluid.p_uncut(i,(N+j-1)%N) == Fluid.p(i,(N+j-1)%N) &&
					Fluid.p_uncut(i+1,j) == Fluid.p(i+1,j) &&
					Fluid.p_uncut(i-1,j) == Fluid.p(i-1,j)
					)
				
				//if(Fluid.boundary(i,j) == -1)
				{
					double n = b(i,j) -	ap(i,j)*Fluid.p_uncut(i,j) + 
												ae(i,j)*Fluid.p_uncut(i,(N+j+1)%N) + 
												aw(i,j)*Fluid.p_uncut(i,(N+j-1)%N) + 
												an(i,j)*Fluid.p_uncut(i+1,j) + 
												as(i,j)*Fluid.p_uncut(i-1,j);
					num += n*n;
					dom += pow(ap(i,j)*Fluid.p_uncut(i,j),2.0);					
				}
			}
		}

		return pow(num,0.5)/pow(dom,0.5);		

	}

}
double CSlipperGap::presidPC()
{

	int M = Fluid.M;
	int N = Fluid.N;
	
	Array<double,2> mu(M,N);		mu = 0;			//Actual visco
	
	Array<double,2> hn(M,N);		hn = 0;			//Gap height, north face
	Array<double,2> hs(M,N);		hs = 0;			//Gap height, south face
	Array<double,2> he(M,N);		he = 0;			//Gap height, east face
	Array<double,2> hw(M,N);		hw = 0;			//Gap height, west face

	Array<double,2> htn(M,N);		htn = 0;			//Gap height top, north face
	Array<double,2> hts(M,N);		hts = 0;			//Gap height top, south face
	Array<double,2> hte(M,N);		hte = 0;			//Gap height top, east face
	Array<double,2> htw(M,N);		htw = 0;			//Gap height top, west face

	Array<double,2> hbn(M,N);		hbn = 0;			//Gap height bottom, north face
	Array<double,2> hbs(M,N);		hbs = 0;			//Gap height bottom, south face
	Array<double,2> hbe(M,N);		hbe = 0;			//Gap height bottom, east face
	Array<double,2> hbw(M,N);		hbw = 0;			//Gap height bottom, west face

	Array<double,2> mun(M,N);		mun = 0;			//Visco, north face
	Array<double,2> mus(M,N);		mus = 0;			//Visco, south face
	Array<double,2> mue(M,N);		mue = 0;			//Visco, east face
	Array<double,2> muw(M,N);		muw = 0;			//Visco, west face

	Array<double,2> an(M,N);		an = 0;			//Finite volume coefficent an
	Array<double,2> as(M,N);		as = 0;			//Finite volume coeffcient as
	Array<double,2> ae(M,N);		ae = 0;			//Finite volume coefficent ae
	Array<double,2> aw(M,N);		aw = 0;			//Finite volume coefficent aw
	Array<double,2> ap(M,N);		ap = 0;			//Finite volume coefficient ap
	Array<double,2> b(M,N);			b = 0;			//Reynolds source term
	
	//Linaer interpolate cell height to cell faces
	cell2face(Fluid.h, he, hw, hn, hs);

	//Average viscosity over gap height
	mu = mean(Fluid.oilviscosity,tensor::k);

	//Linaer interpolate cell visco to cell faces
	cell2face(mu, mue, muw, mun, mus);

	//Linaer interpolate top gap height to cell faces
	Array<double,2> hT (Fluid.h + swashplate.ehd);
	cell2face(hT, hte, htw, htn, hts);
	
	//Linaer interpolate bottom gap height to cell faces
	Array<double,2> hB (swashplate.ehd.copy());
	cell2face(hB, hbe, hbw, hbn, hbs);

	//Added phi flow factor (Patir and Cheng)
	const double sigma = gapinput->options_slipper.general.RoughnessRq*1.0e-6;

	//Spatial coefficients
	if (gapinput->options_slipper.general.reynolds_mu == 2)
	{
		//Full Mu
		ae = (pow(he,3.0)*Fluid.dr)/(mue*Fluid.r * (Fluid.dtheta+East(Fluid.dtheta))/2.0 );
		aw = (pow(hw,3.0)*Fluid.dr)/(muw*Fluid.r * (Fluid.dtheta+West(Fluid.dtheta))/2.0 );
		an = (pow(hn,3.0)*(Fluid.r+Fluid.dr/2.0)*Fluid.dtheta)/(mun * (Fluid.dr+North(Fluid.dr))/2.0 );
		as = (pow(hs,3.0)*(Fluid.r-Fluid.dr/2.0)*Fluid.dtheta)/(mus * (Fluid.dr+South(Fluid.dr))/2.0 );
		ap = an + as + ae + aw;
	} else {
		//Averaged mu
		ae = ( (1-0.9*exp(-0.56*(he/sigma))) * pow(he,3.0)*Fluid.dr)/(mu*Fluid.r * (Fluid.dtheta+East(Fluid.dtheta))/2.0 );
		aw = ( (1-0.9*exp(-0.56*(hw/sigma))) * pow(hw,3.0)*Fluid.dr)/(mu*Fluid.r * (Fluid.dtheta+West(Fluid.dtheta))/2.0 );
		an = ( (1-0.9*exp(-0.56*(hn/sigma))) * pow(hn,3.0)*(Fluid.r+Fluid.dr/2.0)*Fluid.dtheta)/(mu * (Fluid.dr+North(Fluid.dr))/2.0 );
		as = ( (1-0.9*exp(-0.56*(hs/sigma))) * pow(hs,3.0)*(Fluid.r-Fluid.dr/2.0)*Fluid.dtheta)/(mu * (Fluid.dr+South(Fluid.dr))/2.0 );
		ap = an + as + ae + aw;
	}

	//Source Term
	b =	(
				-( (operatingslippergap.tvr+operatingslippergap.bvr)/2.0 * (hn-hs)/Fluid.dr )
				-( (operatingslippergap.tvtheta+operatingslippergap.bvtheta)/2.0 * (he-hw)/(Fluid.r*Fluid.dtheta) )
				+( operatingslippergap.tvr * (htn-hts)/Fluid.dr )
				+( operatingslippergap.tvtheta * (hte-htw)/(Fluid.r*Fluid.dtheta) )
				-( operatingslippergap.bvr * (hbn-hbs)/Fluid.dr )
				-( operatingslippergap.bvtheta * (hbe-hbw)/(Fluid.r*Fluid.dtheta) )
				-( Fluid.dht )
			)*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	
	//Added contact factor (Wu and Zheng)
	{
		Array<double,2> H(Fluid.h/sigma);
		b *= where(H < 3, exp(-0.6912+0.782*H-0.304*H*H+0.0401*H*H*H), 1);
	}

	//b = where(Fluid.contact == 0, b, 0);

	//Calculate the initial residual L2 norm
	{
		double num = 0;
		double dom = 0;
		
		for(int i=1;i<M-1;i++)
		{
			for(int j=0;j<N;j++)	
			{
				
				if(Fluid.boundary(i,j) == -1 && 
					Fluid.p_uncut(i,j) == Fluid.p(i,j) &&
					Fluid.p_uncut(i,(N+j+1)%N) == Fluid.p(i,(N+j+1)%N) &&
					Fluid.p_uncut(i,(N+j-1)%N) == Fluid.p(i,(N+j-1)%N) &&
					Fluid.p_uncut(i+1,j) == Fluid.p(i+1,j) &&
					Fluid.p_uncut(i-1,j) == Fluid.p(i-1,j)
					)
				
				//if(Fluid.boundary(i,j) == -1)
				{
					double n = b(i,j) -	ap(i,j)*Fluid.p_uncut(i,j) + 
												ae(i,j)*Fluid.p_uncut(i,(N+j+1)%N) + 
												aw(i,j)*Fluid.p_uncut(i,(N+j-1)%N) + 
												an(i,j)*Fluid.p_uncut(i+1,j) + 
												as(i,j)*Fluid.p_uncut(i-1,j);
					num += n*n;
					dom += pow(ap(i,j)*Fluid.p_uncut(i,j),2.0);					
				}
			}
		}

		return pow(num,0.5)/pow(dom,0.5);		

	}

}
int CSlipperGap::SlipperReynolds(double alphaPold)		//SlipperGap Reynolds Equation
{
	using namespace blitz::tensor; 
	int M = Fluid.M;
	int N = Fluid.N;
	
	Array<double,2> mu(M,N);		mu = 0;			//Actual visco
	
	Array<double,2> hn(M,N);		hn = 0;			//Gap height, north face
	Array<double,2> hs(M,N);		hs = 0;			//Gap height, south face
	Array<double,2> he(M,N);		he = 0;			//Gap height, east face
	Array<double,2> hw(M,N);		hw = 0;			//Gap height, west face

	Array<double,2> htn(M,N);		htn = 0;			//Gap height top, north face
	Array<double,2> hts(M,N);		hts = 0;			//Gap height top, south face
	Array<double,2> hte(M,N);		hte = 0;			//Gap height top, east face
	Array<double,2> htw(M,N);		htw = 0;			//Gap height top, west face

	Array<double,2> hbn(M,N);		hbn = 0;			//Gap height bottom, north face
	Array<double,2> hbs(M,N);		hbs = 0;			//Gap height bottom, south face
	Array<double,2> hbe(M,N);		hbe = 0;			//Gap height bottom, east face
	Array<double,2> hbw(M,N);		hbw = 0;			//Gap height bottom, west face

	Array<double,2> mun(M,N);		mun = 0;			//Visco, north face
	Array<double,2> mus(M,N);		mus = 0;			//Visco, south face
	Array<double,2> mue(M,N);		mue = 0;			//Visco, east face
	Array<double,2> muw(M,N);		muw = 0;			//Visco, west face

	Array<double,2> an(M,N);		an = 0;			//Finite volume coefficent an
	Array<double,2> as(M,N);		as = 0;			//Finite volume coeffcient as
	Array<double,2> ae(M,N);		ae = 0;			//Finite volume coefficent ae
	Array<double,2> aw(M,N);		aw = 0;			//Finite volume coefficent aw
	Array<double,2> ap(M,N);		ap = 0;			//Finite volume coefficient ap
	Array<double,2> b(M,N);			b = 0;			//Reynolds source term
	
	//Linaer interpolate cell height to cell faces
	cell2face(Fluid.h, he, hw, hn, hs);

	//Average viscosity over gap height
	mu = mean(Fluid.oilviscosity,tensor::k);

	//Linaer interpolate cell visco to cell faces
	cell2face(mu, mue, muw, mun, mus);

	//Linaer interpolate top gap height to cell faces
	Array<double,2> hT (Fluid.h + swashplate.ehd);
	cell2face(hT, hte, htw, htn, hts);
	
	//Linaer interpolate bottom gap height to cell faces
	Array<double,2> hB (swashplate.ehd.copy());
	cell2face(hB, hbe, hbw, hbn, hbs);

	//
	Array<double,2> vr((operatingslippergap.tvr+operatingslippergap.bvr)/2.0);
	Array<double,2> vrn(M,N);
	Array<double,2> vrs(M,N);
	Array<double,2> vre(M,N);
	Array<double,2> vrw(M,N);
	cell2face(vr, vre, vrw, vrn, vrs);

	Array<double,2> vtheta((operatingslippergap.tvtheta+operatingslippergap.bvtheta)/2.0);
	Array<double,2> vthetan(M,N);
	Array<double,2> vthetas(M,N);
	Array<double,2> vthetae(M,N);
	Array<double,2> vthetaw(M,N);
	cell2face(vtheta, vthetae, vthetaw, vthetan, vthetas);


	//Spatial coefficients
	if (gapinput->options_slipper.general.reynolds_mu == 2)
	{
		//Full Mu
		ae = (pow(he,3.0)*Fluid.dr)/(mue*Fluid.r * (Fluid.dtheta+East(Fluid.dtheta))/2.0 );
		aw = (pow(hw,3.0)*Fluid.dr)/(muw*Fluid.r * (Fluid.dtheta+West(Fluid.dtheta))/2.0 );
		an = (pow(hn,3.0)*(Fluid.r+Fluid.dr/2.0)*Fluid.dtheta)/(mun * (Fluid.dr+North(Fluid.dr))/2.0 );
		as = (pow(hs,3.0)*(Fluid.r-Fluid.dr/2.0)*Fluid.dtheta)/(mus * (Fluid.dr+South(Fluid.dr))/2.0 );
		ap = an + as + ae + aw;
	} else {
		//Averaged mu
		ae = (pow(he,3.0)*Fluid.dr)/(mu*Fluid.r * (Fluid.dtheta+East(Fluid.dtheta))/2.0 );
		aw = (pow(hw,3.0)*Fluid.dr)/(mu*Fluid.r * (Fluid.dtheta+West(Fluid.dtheta))/2.0 );
		an = (pow(hn,3.0)*(Fluid.r+Fluid.dr/2.0)*Fluid.dtheta)/(mu * (Fluid.dr+North(Fluid.dr))/2.0 );
		as = (pow(hs,3.0)*(Fluid.r-Fluid.dr/2.0)*Fluid.dtheta)/(mu * (Fluid.dr+South(Fluid.dr))/2.0 );
		ap = an + as + ae + aw;
	}

	//Source Term
	b =	(
				-( (operatingslippergap.tvr+operatingslippergap.bvr)/2.0 * (hn-hs)/Fluid.dr )
				-( (operatingslippergap.tvtheta+operatingslippergap.bvtheta)/2.0 * (he-hw)/(Fluid.r*Fluid.dtheta) )
				-( (operatingslippergap.tvr+operatingslippergap.bvr)/2.0 * Fluid.h/Fluid.r)
				-( Fluid.h * (vrn-vrs)/Fluid.dr )
				-( Fluid.h * (vthetae-vthetaw)/(Fluid.r*Fluid.dtheta) )
				+( operatingslippergap.tvr * (htn-hts)/Fluid.dr )
				+( operatingslippergap.tvtheta * (hte-htw)/(Fluid.r*Fluid.dtheta) )
				-( operatingslippergap.bvr * (hbn-hbs)/Fluid.dr )
				-( operatingslippergap.bvtheta * (hbe-hbw)/(Fluid.r*Fluid.dtheta) )
				-( Fluid.dht )
			)*12*Fluid.r*Fluid.dr*Fluid.dtheta;

	if(curFullLoop <= 0 && gapinput->options_slipper.general.EHDsqueeze && gapinput->options_slipper.general.SlipperPressureDeformation)
	{
		//double IMval = 9.07295e-008/100e+5;
		//Array<double, 2> sqzfactor(IMval/gapinput->common.timestep*12*Fluid.r*Fluid.dr*Fluid.dtheta);
		Array<double, 2> sqzfactor(Fluid.newt_ehdsqzapprox(tensor::i)/gapinput->common.timestep*12*Fluid.r*Fluid.dr*Fluid.dtheta);
		ap += sqzfactor;
		b += Fluid.pFullLoop*sqzfactor;
	}

	//const double wearHeight = 2.0*Fluid.ContactHeight;
	//b = where(Fluid.h <= wearHeight, 0, b);
	//b = where(Fluid.contact == 0, b, 0);
	/*
	Fluid.ReyVals(0,all,all) = (-( (operatingslippergap.tvr+operatingslippergap.bvr)/2.0 * (hn-hs)/Fluid.dr ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(1,all,all) = (-( (operatingslippergap.tvtheta+operatingslippergap.bvtheta)/2.0 * (he-hw)/(Fluid.r*Fluid.dtheta) ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(2,all,all) = (( operatingslippergap.tvr * (htn-hts)/Fluid.dr ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;;
	Fluid.ReyVals(3,all,all) = ( operatingslippergap.tvtheta * (hte-htw)/(Fluid.r*Fluid.dtheta) )*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(4,all,all) = (-( operatingslippergap.bvr * (hbn-hbs)/Fluid.dr ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(5,all,all) = (-( operatingslippergap.bvtheta * (hbe-hbw)/(Fluid.r*Fluid.dtheta) ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(6,all,all) = (-( Fluid.dht ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(7,all,all) = ae;
	Fluid.ReyVals(8,all,all) = aw;
	Fluid.ReyVals(9,all,all) = an;
	Fluid.ReyVals(10,all,all) = as;
	
	Fluid.ReyVals(0,all,all) = where(Fluid.h <= wearHeight, 0, Fluid.ReyVals(0,all,all));
	Fluid.ReyVals(1,all,all) = where(Fluid.h <= wearHeight, 0, Fluid.ReyVals(1,all,all));
	Fluid.ReyVals(2,all,all) = where(Fluid.h <= wearHeight, 0, Fluid.ReyVals(2,all,all));
	Fluid.ReyVals(3,all,all) = where(Fluid.h <= wearHeight, 0, Fluid.ReyVals(3,all,all));
	Fluid.ReyVals(4,all,all) = where(Fluid.h <= wearHeight, 0, Fluid.ReyVals(4,all,all));
	Fluid.ReyVals(5,all,all) = where(Fluid.h <= wearHeight, 0, Fluid.ReyVals(5,all,all));
	Fluid.ReyVals(6,all,all) = where(Fluid.h <= wearHeight, 0, Fluid.ReyVals(6,all,all));
	*/
	/*
	Fluid.ReyVals(0,all,all) =	(-( (operatingslippergap.tvr+operatingslippergap.bvr)/2.0 * (hn-hs)/Fluid.dr ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(1,all,all) =	(-( (operatingslippergap.tvtheta+operatingslippergap.bvtheta)/2.0 * (he-hw)/(Fluid.r*Fluid.dtheta) ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(2,all,all) =	( operatingslippergap.tvr * (htn-hts)/Fluid.dr )*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(3,all,all) =	( operatingslippergap.tvtheta * (hte-htw)/(Fluid.r*Fluid.dtheta) )*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(4,all,all) =	(-( Fluid.dht ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;

	Fluid.ReyVals(5,all,all) =	hn-hs;
	Fluid.ReyVals(6,all,all) =	he-hw;
	Fluid.ReyVals(7,all,all) =	htn-hts;
	Fluid.ReyVals(8,all,all) =	hte-htw;

	Fluid.ReyVals(9,all,all) =	operatingslippergap.tvr;
	Fluid.ReyVals(10,all,all) =	operatingslippergap.tvtheta;
	*/

	//Fluid.ReyVals(all,all,all) = where(Fluid.contact(tensor::j,tensor::k) == 0, Fluid.ReyVals(tensor::i,tensor::j,tensor::k), 0);

	//Renumber the non-boundary cells with mid
	Array<int,2> mid(Fluid.M,Fluid.N);
	mid = -1;
	int cells = 0;
	{
		int id=0;

		for(int i=0; i<Fluid.M; i++)
		{
			for(int j=0; j<Fluid.N; j++)
			{
				if(Fluid.boundary(i,j) == -1)
				{
					mid(i,j) = id;
					id++;
					cells++;
				}
			}
		}
	}

	//The A Matrix
	gmm::col_matrix<gmm::wsvector<double> > A(cells,cells);

	//The B & X  Vector
	vector<double> B(cells),X(cells);

	//Copy the initial P to X
	for(int i=0;i<M-2;i++)
	{
		for(int j=0;j<N;j++)	
		{
			if(Fluid.boundary(i,j) == -1)
			{
				const int row = mid(i,j);
				X[row] = Fluid.p_uncut(i,j);
			}
		}
	}

	//Build the K array inside
	for(int i=0;i<M;i++)
	{
		for(int j=0;j<N;j++)	
		{
			if(Fluid.boundary(i,j) == -1)
			{
				int x,y;

				if (j==0)
				{
					x = N - 1;
				} 
				else 
				{ 
					x = j - 1;
				}
				if (j==N-1)
				{
					y = 0;
				} 
				else 
				{ 
					y = j + 1;
				}

				/*
				an(i,j) * pnew(i+1,j)
				as(i,j) * pnew(i-1,j) 
				ae(i,j) * pnew(i,y) 
				aw(i,j) * pnew(i,x) 
				b(i,j)
				ap(i,j)
				*/

				const int row = mid(i,j);

				B[row] = b(i,j);

				if(mid(i-1,j) >= 0)
				{
					A(row,mid(i-1,j)) = -as(i,j);
				} else {
					//boundary
					B[row] += as(i,j)*Fluid.p_uncut(i-1,j);
				}

				if(mid(i+1,j) >= 0)
				{
					A(row,mid(i+1,j)) = -an(i,j);
				} else {
					//boundary
					B[row] += an(i,j)*Fluid.p_uncut(i+1,j);
				}

				if(mid(i,y) >= 0)
				{
					A(row,mid(i,y)) = -ae(i,j);
				} else {
					//boundary
					B[row] += ae(i,j)*Fluid.p_uncut(i,y);
				}

				if(mid(i,x) >= 0)
				{
					A(row,mid(i,x)) = -aw(i,j);
				} else {
					//boundary
					B[row] += aw(i,j)*Fluid.p_uncut(i,x);
				}

				//A(row,row) = ap(i,j);

				A(row,row) = ap(i,j)/alphaPold;
				B[row] += (1-alphaPold)*ap(i,j)*Fluid.p(i,j)/alphaPold;
			}
		}
	}

	//Set up the proper sparse matrix
	gmm::csc_matrix<double> AS;
	gmm::copy(A,AS);

	//Set solver params
	gmm::iteration iter(1e-9);
	iter.set_maxiter(5000);
	iter.set_noisy(0);
	
	//Precondition
	gmm::ildlt_precond< gmm::csc_matrix<double> > PR(AS);
	//gmm::diagonal_precond< gmm::csc_matrix<double> > PR(AS);

	//Solve
	gmm::bicgstab(AS, X, B, PR, iter);

	if(iter.get_iteration() == iter.get_maxiter())
	{
		GapLog.message("WARNING: Reynolds pressure loop failed to converge!");
	}

	//Copy the result back to the Fluid.p
	for(int i=0;i<M;i++)
	{
		for(int j=0;j<N;j++)	
		{
			if(Fluid.boundary(i,j) == -1)
			{
				const int row = mid(i,j);
				Fluid.p_uncut(i,j) = X[row];
			}
		}
	}

	/*
	Fluid.p = where(Fluid.p_uncut<0.1e+5,Fluid.p - alphaPold * (Fluid.p - 0.1e+5),			//low pressure cut relax
						where(Fluid.p_uncut>1000.0e+5,Fluid.p - alphaPold * (Fluid.p - 1000.0e+5),		//high pressure cut relax
							Fluid.p_uncut));		//no extra relax needed
	*/
	Fluid.p = Fluid.p_uncut;

	//Just to make certain
	Fluid.p = where(Fluid.p<0.8e+5, 0.8e+5 ,Fluid.p);
	Fluid.p = where(Fluid.p>2000.0e+5,2000.0e+5,Fluid.p);

	return (int) iter.get_iteration();
}
int CSlipperGap::SlipperReynoldsPC(double alphaPold)		//SlipperGap Reynolds Equation
{
	using namespace blitz::tensor; 
	int M = Fluid.M;
	int N = Fluid.N;
	
	Array<double,2> mu(M,N);		mu = 0;			//Actual visco
	
	Array<double,2> hn(M,N);		hn = 0;			//Gap height, north face
	Array<double,2> hs(M,N);		hs = 0;			//Gap height, south face
	Array<double,2> he(M,N);		he = 0;			//Gap height, east face
	Array<double,2> hw(M,N);		hw = 0;			//Gap height, west face

	Array<double,2> htn(M,N);		htn = 0;			//Gap height top, north face
	Array<double,2> hts(M,N);		hts = 0;			//Gap height top, south face
	Array<double,2> hte(M,N);		hte = 0;			//Gap height top, east face
	Array<double,2> htw(M,N);		htw = 0;			//Gap height top, west face

	Array<double,2> hbn(M,N);		hbn = 0;			//Gap height bottom, north face
	Array<double,2> hbs(M,N);		hbs = 0;			//Gap height bottom, south face
	Array<double,2> hbe(M,N);		hbe = 0;			//Gap height bottom, east face
	Array<double,2> hbw(M,N);		hbw = 0;			//Gap height bottom, west face

	Array<double,2> mun(M,N);		mun = 0;			//Visco, north face
	Array<double,2> mus(M,N);		mus = 0;			//Visco, south face
	Array<double,2> mue(M,N);		mue = 0;			//Visco, east face
	Array<double,2> muw(M,N);		muw = 0;			//Visco, west face

	Array<double,2> an(M,N);		an = 0;			//Finite volume coefficent an
	Array<double,2> as(M,N);		as = 0;			//Finite volume coeffcient as
	Array<double,2> ae(M,N);		ae = 0;			//Finite volume coefficent ae
	Array<double,2> aw(M,N);		aw = 0;			//Finite volume coefficent aw
	Array<double,2> ap(M,N);		ap = 0;			//Finite volume coefficient ap
	Array<double,2> b(M,N);			b = 0;			//Reynolds source term
	
	//Linaer interpolate cell height to cell faces
	cell2face(Fluid.h, he, hw, hn, hs);

	//Average viscosity over gap height
	mu = mean(Fluid.oilviscosity,tensor::k);

	//Linaer interpolate cell visco to cell faces
	cell2face(mu, mue, muw, mun, mus);

	//Linaer interpolate top gap height to cell faces
	Array<double,2> hT (Fluid.h + swashplate.ehd);
	cell2face(hT, hte, htw, htn, hts);
	
	//Linaer interpolate bottom gap height to cell faces
	Array<double,2> hB (swashplate.ehd.copy());
	cell2face(hB, hbe, hbw, hbn, hbs);
	
	//Added phi flow factor (Patir and Cheng)
	const double sigma = gapinput->options_slipper.general.RoughnessRq*1.0e-6;

	//Spatial coefficients
	if (gapinput->options_slipper.general.reynolds_mu == 2)
	{
		//Full Mu
		ae = (pow(he,3.0)*Fluid.dr)/(mue*Fluid.r * (Fluid.dtheta+East(Fluid.dtheta))/2.0 );
		aw = (pow(hw,3.0)*Fluid.dr)/(muw*Fluid.r * (Fluid.dtheta+West(Fluid.dtheta))/2.0 );
		an = (pow(hn,3.0)*(Fluid.r+Fluid.dr/2.0)*Fluid.dtheta)/(mun * (Fluid.dr+North(Fluid.dr))/2.0 );
		as = (pow(hs,3.0)*(Fluid.r-Fluid.dr/2.0)*Fluid.dtheta)/(mus * (Fluid.dr+South(Fluid.dr))/2.0 );
		ap = an + as + ae + aw;
	} else {
		//Averaged mu
		ae = ( (1-0.9*exp(-0.56*(he/sigma))) * pow(he,3.0)*Fluid.dr)/(mu*Fluid.r * (Fluid.dtheta+East(Fluid.dtheta))/2.0 );
		aw = ( (1-0.9*exp(-0.56*(hw/sigma))) * pow(hw,3.0)*Fluid.dr)/(mu*Fluid.r * (Fluid.dtheta+West(Fluid.dtheta))/2.0 );
		an = ( (1-0.9*exp(-0.56*(hn/sigma))) * pow(hn,3.0)*(Fluid.r+Fluid.dr/2.0)*Fluid.dtheta)/(mu * (Fluid.dr+North(Fluid.dr))/2.0 );
		as = ( (1-0.9*exp(-0.56*(hs/sigma))) * pow(hs,3.0)*(Fluid.r-Fluid.dr/2.0)*Fluid.dtheta)/(mu * (Fluid.dr+South(Fluid.dr))/2.0 );
		ap = an + as + ae + aw;
	}

	//Source Term
	b =	(
				-( (operatingslippergap.tvr+operatingslippergap.bvr)/2.0 * (hn-hs)/Fluid.dr )
				-( (operatingslippergap.tvtheta+operatingslippergap.bvtheta)/2.0 * (he-hw)/(Fluid.r*Fluid.dtheta) )
				+( operatingslippergap.tvr * (htn-hts)/Fluid.dr )
				+( operatingslippergap.tvtheta * (hte-htw)/(Fluid.r*Fluid.dtheta) )
				-( operatingslippergap.bvr * (hbn-hbs)/Fluid.dr )
				-( operatingslippergap.bvtheta * (hbe-hbw)/(Fluid.r*Fluid.dtheta) )
				-( Fluid.dht )
			)*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	
	//Added contact factor (Wu and Zheng)
	{
		Array<double,2> H(Fluid.h/sigma);
		b *= where(H < 3, exp(-0.6912+0.782*H-0.304*H*H+0.0401*H*H*H), 1);
	}

	//b = where(Fluid.contact == 0, b, 0);

	/*
	Fluid.ReyVals(0,all,all) =	(-( (operatingslippergap.tvr+operatingslippergap.bvr)/2.0 * (hn-hs)/Fluid.dr ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(1,all,all) =	(-( (operatingslippergap.tvtheta+operatingslippergap.bvtheta)/2.0 * (he-hw)/(Fluid.r*Fluid.dtheta) ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(2,all,all) =	( operatingslippergap.tvr * (htn-hts)/Fluid.dr )*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(3,all,all) =	( operatingslippergap.tvtheta * (hte-htw)/(Fluid.r*Fluid.dtheta) )*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(4,all,all) =	(-( Fluid.dht ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;

	Fluid.ReyVals(5,all,all) =	hn-hs;
	Fluid.ReyVals(6,all,all) =	he-hw;
	Fluid.ReyVals(7,all,all) =	htn-hts;
	Fluid.ReyVals(8,all,all) =	hte-htw;

	Fluid.ReyVals(9,all,all) =	operatingslippergap.tvr;
	Fluid.ReyVals(10,all,all) =	operatingslippergap.tvtheta;
	*/

	//Fluid.ReyVals(all,all,all) = where(Fluid.contact(tensor::j,tensor::k) == 0, Fluid.ReyVals(tensor::i,tensor::j,tensor::k), 0);

	//Renumber the non-boundary cells with mid
	Array<int,2> mid(Fluid.M,Fluid.N);
	mid = -1;
	int cells = 0;
	{
		int id=0;

		for(int i=0; i<Fluid.M; i++)
		{
			for(int j=0; j<Fluid.N; j++)
			{
				if(Fluid.boundary(i,j) == -1)
				{
					mid(i,j) = id;
					id++;
					cells++;
				}
			}
		}
	}

	//The A Matrix
	gmm::col_matrix<gmm::wsvector<double> > A(cells,cells);

	//The B & X  Vector
	vector<double> B(cells),X(cells);

	//Copy the initial P to X
	for(int i=0;i<M-2;i++)
	{
		for(int j=0;j<N;j++)	
		{
			if(Fluid.boundary(i,j) == -1)
			{
				const int row = mid(i,j);
				X[row] = Fluid.p_uncut(i,j);
			}
		}
	}

	//Build the K array inside
	for(int i=0;i<M;i++)
	{
		for(int j=0;j<N;j++)	
		{
			if(Fluid.boundary(i,j) == -1)
			{
				int x,y;

				if (j==0)
				{
					x = N - 1;
				} 
				else 
				{ 
					x = j - 1;
				}
				if (j==N-1)
				{
					y = 0;
				} 
				else 
				{ 
					y = j + 1;
				}

				/*
				an(i,j) * pnew(i+1,j)
				as(i,j) * pnew(i-1,j) 
				ae(i,j) * pnew(i,y) 
				aw(i,j) * pnew(i,x) 
				b(i,j)
				ap(i,j)
				*/

				const int row = mid(i,j);

				B[row] = b(i,j);

				if(mid(i-1,j) >= 0)
				{
					A(row,mid(i-1,j)) = -as(i,j);
				} else {
					//boundary
					B[row] += as(i,j)*Fluid.p_uncut(i-1,j);
				}

				if(mid(i+1,j) >= 0)
				{
					A(row,mid(i+1,j)) = -an(i,j);
				} else {
					//boundary
					B[row] += an(i,j)*Fluid.p_uncut(i+1,j);
				}

				if(mid(i,y) >= 0)
				{
					A(row,mid(i,y)) = -ae(i,j);
				} else {
					//boundary
					B[row] += ae(i,j)*Fluid.p_uncut(i,y);
				}

				if(mid(i,x) >= 0)
				{
					A(row,mid(i,x)) = -aw(i,j);
				} else {
					//boundary
					B[row] += aw(i,j)*Fluid.p_uncut(i,x);
				}

				//A(row,row) = ap(i,j);

				A(row,row) = ap(i,j)/alphaPold;
				B[row] += (1-alphaPold)*ap(i,j)*Fluid.p(i,j)/alphaPold;
			}
		}
	}

	//Set up the proper sparse matrix
	gmm::csc_matrix<double> AS;
	gmm::copy(A,AS);

	//Set solver params
	gmm::iteration iter(1e-9);
	iter.set_maxiter(5000);
	iter.set_noisy(0);
	
	//Precondition
	gmm::ildlt_precond< gmm::csc_matrix<double> > PR(AS);
	//gmm::diagonal_precond< gmm::csc_matrix<double> > PR(AS);

	//Solve
	gmm::bicgstab(AS, X, B, PR, iter);

	if(iter.get_iteration() == iter.get_maxiter())
	{
		GapLog.message("WARNING: Reynolds pressure loop failed to converge!");
	}

	//Copy the result back to the Fluid.p
	for(int i=0;i<M;i++)
	{
		for(int j=0;j<N;j++)	
		{
			if(Fluid.boundary(i,j) == -1)
			{
				const int row = mid(i,j);
				Fluid.p_uncut(i,j) = X[row];
			}
		}
	}

	/*
	Fluid.p = where(Fluid.p_uncut<0.1e+5,Fluid.p - alphaPold * (Fluid.p - 0.1e+5),			//low pressure cut relax
						where(Fluid.p_uncut>1000.0e+5,Fluid.p - alphaPold * (Fluid.p - 1000.0e+5),		//high pressure cut relax
							Fluid.p_uncut));		//no extra relax needed
	*/
	Fluid.p = Fluid.p_uncut;

	//Just to make certain
	Fluid.p = where(Fluid.p<0.2e+5, 0.2e+5 ,Fluid.p);
	Fluid.p = where(Fluid.p>800.0e+5,800.0e+5,Fluid.p);

	return (int) iter.get_iteration();
}
int CSlipperGap::SlipperReynoldsFkSqz(double alphaPold)		//SlipperGap Reynolds Equation
{
	using namespace blitz::tensor; 
	int M = Fluid.M;
	int N = Fluid.N;
	
	Array<double,2> mu(M,N);		mu = 0;			//Actual visco
	
	Array<double,2> hn(M,N);		hn = 0;			//Gap height, north face
	Array<double,2> hs(M,N);		hs = 0;			//Gap height, south face
	Array<double,2> he(M,N);		he = 0;			//Gap height, east face
	Array<double,2> hw(M,N);		hw = 0;			//Gap height, west face

	Array<double,2> htn(M,N);		htn = 0;			//Gap height top, north face
	Array<double,2> hts(M,N);		hts = 0;			//Gap height top, south face
	Array<double,2> hte(M,N);		hte = 0;			//Gap height top, east face
	Array<double,2> htw(M,N);		htw = 0;			//Gap height top, west face

	Array<double,2> hbn(M,N);		hbn = 0;			//Gap height bottom, north face
	Array<double,2> hbs(M,N);		hbs = 0;			//Gap height bottom, south face
	Array<double,2> hbe(M,N);		hbe = 0;			//Gap height bottom, east face
	Array<double,2> hbw(M,N);		hbw = 0;			//Gap height bottom, west face

	Array<double,2> mun(M,N);		mun = 0;			//Visco, north face
	Array<double,2> mus(M,N);		mus = 0;			//Visco, south face
	Array<double,2> mue(M,N);		mue = 0;			//Visco, east face
	Array<double,2> muw(M,N);		muw = 0;			//Visco, west face

	Array<double,2> an(M,N);		an = 0;			//Finite volume coefficent an
	Array<double,2> as(M,N);		as = 0;			//Finite volume coeffcient as
	Array<double,2> ae(M,N);		ae = 0;			//Finite volume coefficent ae
	Array<double,2> aw(M,N);		aw = 0;			//Finite volume coefficent aw
	Array<double,2> ap(M,N);		ap = 0;			//Finite volume coefficient ap
	Array<double,2> b(M,N);			b = 0;			//Reynolds source term
	
	//Linaer interpolate cell height to cell faces
	cell2face(Fluid.h, he, hw, hn, hs);

	//Average viscosity over gap height
	mu = mean(Fluid.oilviscosity,tensor::k);

	//Linaer interpolate cell visco to cell faces
	cell2face(mu, mue, muw, mun, mus);

	//Linaer interpolate top gap height to cell faces
	Array<double,2> hT (Fluid.h + swashplate.ehd);
	cell2face(hT, hte, htw, htn, hts);
	
	//Linaer interpolate bottom gap height to cell faces
	Array<double,2> hB (swashplate.ehd.copy());
	cell2face(hB, hbe, hbw, hbn, hbs);

	//Spatial coefficients
	if (gapinput->options_slipper.general.reynolds_mu == 2)
	{
		//Full Mu
		ae = (pow(he,3.0)*Fluid.dr)/(mue*Fluid.r * (Fluid.dtheta+East(Fluid.dtheta))/2.0 );
		aw = (pow(hw,3.0)*Fluid.dr)/(muw*Fluid.r * (Fluid.dtheta+West(Fluid.dtheta))/2.0 );
		an = (pow(hn,3.0)*(Fluid.r+Fluid.dr/2.0)*Fluid.dtheta)/(mun * (Fluid.dr+North(Fluid.dr))/2.0 );
		as = (pow(hs,3.0)*(Fluid.r-Fluid.dr/2.0)*Fluid.dtheta)/(mus * (Fluid.dr+South(Fluid.dr))/2.0 );
		ap = an + as + ae + aw;
	} else {
		//Averaged mu
		ae = (pow(he,3.0)*Fluid.dr)/(mu*Fluid.r * (Fluid.dtheta+East(Fluid.dtheta))/2.0 );
		aw = (pow(hw,3.0)*Fluid.dr)/(mu*Fluid.r * (Fluid.dtheta+West(Fluid.dtheta))/2.0 );
		an = (pow(hn,3.0)*(Fluid.r+Fluid.dr/2.0)*Fluid.dtheta)/(mu * (Fluid.dr+North(Fluid.dr))/2.0 );
		as = (pow(hs,3.0)*(Fluid.r-Fluid.dr/2.0)*Fluid.dtheta)/(mu * (Fluid.dr+South(Fluid.dr))/2.0 );
		ap = an + as + ae + aw;
	}

	//We need to construct the modified dht term
	double IMslip = 0;
	double IMswash = 0;
	const double TimeStep = gapinput->common.timestep;
	if(gapinput->options_slipper.general.SlipperPressureDeformation == 1)
	{
		IMslip = 2.0*slipper.avgIM / 100e+5;
	}
	if(gapinput->options_slipper.general.SwashplatePressureDeformation == 1)
	{
		IMswash = 2.0*swashplate.avgIM / 100e+5;
	}


	//Source Term
	b =	(
				-( (operatingslippergap.tvr+operatingslippergap.bvr)/2.0 * (hn-hs)/Fluid.dr )
				-( (operatingslippergap.tvtheta+operatingslippergap.bvtheta)/2.0 * (he-hw)/(Fluid.r*Fluid.dtheta) )
				+( operatingslippergap.tvr * (htn-hts)/Fluid.dr )
				+( operatingslippergap.tvtheta * (hte-htw)/(Fluid.r*Fluid.dtheta) )
				-( operatingslippergap.bvr * (hbn-hbs)/Fluid.dr )
				-( operatingslippergap.bvtheta * (hbe-hbw)/(Fluid.r*Fluid.dtheta) )
				-( Fluid.dht - (IMslip-IMswash)/TimeStep*Fluid.p )
			)*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	b = where(Fluid.contact == 0, b, 0);

	ap -= -( (IMslip-IMswash)/TimeStep )*12*Fluid.r*Fluid.dr*Fluid.dtheta;

	/*
	Fluid.ReyVals(0,all,all) =	(-( (operatingslippergap.tvr+operatingslippergap.bvr)/2.0 * (hn-hs)/Fluid.dr ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(1,all,all) =	(-( (operatingslippergap.tvtheta+operatingslippergap.bvtheta)/2.0 * (he-hw)/(Fluid.r*Fluid.dtheta) ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(2,all,all) =	( operatingslippergap.tvr * (htn-hts)/Fluid.dr )*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(3,all,all) =	( operatingslippergap.tvtheta * (hte-htw)/(Fluid.r*Fluid.dtheta) )*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(4,all,all) =	(-( Fluid.dht ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;

	Fluid.ReyVals(5,all,all) =	hn-hs;
	Fluid.ReyVals(6,all,all) =	he-hw;
	Fluid.ReyVals(7,all,all) =	htn-hts;
	Fluid.ReyVals(8,all,all) =	hte-htw;

	Fluid.ReyVals(9,all,all) =	operatingslippergap.tvr;
	Fluid.ReyVals(10,all,all) =	operatingslippergap.tvtheta;
	*/

	//Fluid.ReyVals(all,all,all) = where(Fluid.contact(tensor::j,tensor::k) == 0, Fluid.ReyVals(tensor::i,tensor::j,tensor::k), 0);

	//Renumber the non-boundary cells with mid
	Array<int,2> mid(Fluid.M,Fluid.N);
	mid = -1;
	int cells = 0;
	{
		int id=0;

		for(int i=0; i<Fluid.M; i++)
		{
			for(int j=0; j<Fluid.N; j++)
			{
				if(Fluid.boundary(i,j) == -1)
				{
					mid(i,j) = id;
					id++;
					cells++;
				}
			}
		}
	}

	//The A Matrix
	gmm::col_matrix<gmm::wsvector<double> > A(cells,cells);

	//The B & X  Vector
	vector<double> B(cells),X(cells);

	//Copy the initial P to X
	for(int i=0;i<M-2;i++)
	{
		for(int j=0;j<N;j++)	
		{
			if(Fluid.boundary(i,j) == -1)
			{
				const int row = mid(i,j);
				X[row] = Fluid.p_uncut(i,j);
			}
		}
	}

	//Build the K array inside
	for(int i=0;i<M;i++)
	{
		for(int j=0;j<N;j++)	
		{
			if(Fluid.boundary(i,j) == -1)
			{
				int x,y;

				if (j==0)
				{
					x = N - 1;
				} 
				else 
				{ 
					x = j - 1;
				}
				if (j==N-1)
				{
					y = 0;
				} 
				else 
				{ 
					y = j + 1;
				}

				/*
				an(i,j) * pnew(i+1,j)
				as(i,j) * pnew(i-1,j) 
				ae(i,j) * pnew(i,y) 
				aw(i,j) * pnew(i,x) 
				b(i,j)
				ap(i,j)
				*/

				const int row = mid(i,j);

				B[row] = b(i,j);

				if(mid(i-1,j) >= 0)
				{
					A(row,mid(i-1,j)) = -as(i,j);
				} else {
					//boundary
					B[row] += as(i,j)*Fluid.p_uncut(i-1,j);
				}

				if(mid(i+1,j) >= 0)
				{
					A(row,mid(i+1,j)) = -an(i,j);
				} else {
					//boundary
					B[row] += an(i,j)*Fluid.p_uncut(i+1,j);
				}

				if(mid(i,y) >= 0)
				{
					A(row,mid(i,y)) = -ae(i,j);
				} else {
					//boundary
					B[row] += ae(i,j)*Fluid.p_uncut(i,y);
				}

				if(mid(i,x) >= 0)
				{
					A(row,mid(i,x)) = -aw(i,j);
				} else {
					//boundary
					B[row] += aw(i,j)*Fluid.p_uncut(i,x);
				}

				//A(row,row) = ap(i,j);

				A(row,row) = ap(i,j)/alphaPold;
				B[row] += (1-alphaPold)*ap(i,j)*Fluid.p(i,j)/alphaPold;
			}
		}
	}

	//Set up the proper sparse matrix
	gmm::csc_matrix<double> AS;
	gmm::copy(A,AS);

	//Set solver params
	gmm::iteration iter(1e-9);
	iter.set_maxiter(5000);
	iter.set_noisy(0);
	
	//Precondition
	gmm::ildlt_precond< gmm::csc_matrix<double> > PR(AS);
	//gmm::diagonal_precond< gmm::csc_matrix<double> > PR(AS);

	//Solve
	gmm::bicgstab(AS, X, B, PR, iter);

	if(iter.get_iteration() == iter.get_maxiter())
	{
		GapLog.message("WARNING: Reynolds pressure loop failed to converge!");
	}

	//Copy the result back to the Fluid.p
	for(int i=0;i<M;i++)
	{
		for(int j=0;j<N;j++)	
		{
			if(Fluid.boundary(i,j) == -1)
			{
				const int row = mid(i,j);
				Fluid.p_uncut(i,j) = X[row];
			}
		}
	}

	/*
	Fluid.p = where(Fluid.p_uncut<0.1e+5,Fluid.p - alphaPold * (Fluid.p - 0.1e+5),			//low pressure cut relax
						where(Fluid.p_uncut>1000.0e+5,Fluid.p - alphaPold * (Fluid.p - 1000.0e+5),		//high pressure cut relax
							Fluid.p_uncut));		//no extra relax needed
	*/
	Fluid.p = Fluid.p_uncut;

	//Just to make certain
	Fluid.p = where(Fluid.p<0.2e+5, 0.2e+5 ,Fluid.p);
	Fluid.p = where(Fluid.p>1000.0e+5,1000.0e+5,Fluid.p);

	return (int) iter.get_iteration();
}
int CSlipperGap::SlipperReynoldsCut(double alphaPold)		//SlipperGap Reynolds Equation
{
	using namespace blitz::tensor; 
	int M = Fluid.M;
	int N = Fluid.N;
	
	Array<double,2> mu(M,N);		mu = 0;			//Actual visco
	
	Array<double,2> hn(M,N);		hn = 0;			//Gap height, north face
	Array<double,2> hs(M,N);		hs = 0;			//Gap height, south face
	Array<double,2> he(M,N);		he = 0;			//Gap height, east face
	Array<double,2> hw(M,N);		hw = 0;			//Gap height, west face

	Array<double,2> htn(M,N);		htn = 0;			//Gap height top, north face
	Array<double,2> hts(M,N);		hts = 0;			//Gap height top, south face
	Array<double,2> hte(M,N);		hte = 0;			//Gap height top, east face
	Array<double,2> htw(M,N);		htw = 0;			//Gap height top, west face

	Array<double,2> hbn(M,N);		hbn = 0;			//Gap height bottom, north face
	Array<double,2> hbs(M,N);		hbs = 0;			//Gap height bottom, south face
	Array<double,2> hbe(M,N);		hbe = 0;			//Gap height bottom, east face
	Array<double,2> hbw(M,N);		hbw = 0;			//Gap height bottom, west face

	Array<double,2> mun(M,N);		mun = 0;			//Visco, north face
	Array<double,2> mus(M,N);		mus = 0;			//Visco, south face
	Array<double,2> mue(M,N);		mue = 0;			//Visco, east face
	Array<double,2> muw(M,N);		muw = 0;			//Visco, west face

	Array<double,2> an(M,N);		an = 0;			//Finite volume coefficent an
	Array<double,2> as(M,N);		as = 0;			//Finite volume coeffcient as
	Array<double,2> ae(M,N);		ae = 0;			//Finite volume coefficent ae
	Array<double,2> aw(M,N);		aw = 0;			//Finite volume coefficent aw
	Array<double,2> ap(M,N);		ap = 0;			//Finite volume coefficient ap
	Array<double,2> b(M,N);			b = 0;			//Reynolds source term
	
	//Linaer interpolate cell height to cell faces
	cell2face(Fluid.h, he, hw, hn, hs);

	//Average viscosity over gap height
	mu = mean(Fluid.oilviscosity,tensor::k);

	//Linaer interpolate cell visco to cell faces
	cell2face(mu, mue, muw, mun, mus);

	//Linaer interpolate top gap height to cell faces
	Array<double,2> hT (Fluid.h + swashplate.ehd);
	cell2face(hT, hte, htw, htn, hts);
	
	//Linaer interpolate bottom gap height to cell faces
	Array<double,2> hB (swashplate.ehd.copy());
	cell2face(hB, hbe, hbw, hbn, hbs);

	//Spatial coefficients
	if (gapinput->options_slipper.general.reynolds_mu == 2)
	{
		//Full Mu
		ae = (pow(he,3.0)*Fluid.dr)/(mue*Fluid.r * (Fluid.dtheta+East(Fluid.dtheta))/2.0 );
		aw = (pow(hw,3.0)*Fluid.dr)/(muw*Fluid.r * (Fluid.dtheta+West(Fluid.dtheta))/2.0 );
		an = (pow(hn,3.0)*(Fluid.r+Fluid.dr/2.0)*Fluid.dtheta)/(mun * (Fluid.dr+North(Fluid.dr))/2.0 );
		as = (pow(hs,3.0)*(Fluid.r-Fluid.dr/2.0)*Fluid.dtheta)/(mus * (Fluid.dr+South(Fluid.dr))/2.0 );
		ap = an + as + ae + aw;
	} else {
		//Averaged mu
		ae = (pow(he,3.0)*Fluid.dr)/(mu*Fluid.r * (Fluid.dtheta+East(Fluid.dtheta))/2.0 );
		aw = (pow(hw,3.0)*Fluid.dr)/(mu*Fluid.r * (Fluid.dtheta+West(Fluid.dtheta))/2.0 );
		an = (pow(hn,3.0)*(Fluid.r+Fluid.dr/2.0)*Fluid.dtheta)/(mu * (Fluid.dr+North(Fluid.dr))/2.0 );
		as = (pow(hs,3.0)*(Fluid.r-Fluid.dr/2.0)*Fluid.dtheta)/(mu * (Fluid.dr+South(Fluid.dr))/2.0 );
		ap = an + as + ae + aw;
	}

	//Source Term
	b =	(
				-( (operatingslippergap.tvr+operatingslippergap.bvr)/2.0 * (hn-hs)/Fluid.dr )
				-( (operatingslippergap.tvtheta+operatingslippergap.bvtheta)/2.0 * (he-hw)/(Fluid.r*Fluid.dtheta) )
				+( operatingslippergap.tvr * (htn-hts)/Fluid.dr )
				+( operatingslippergap.tvtheta * (hte-htw)/(Fluid.r*Fluid.dtheta) )
				-( operatingslippergap.bvr * (hbn-hbs)/Fluid.dr )
				-( operatingslippergap.bvtheta * (hbe-hbw)/(Fluid.r*Fluid.dtheta) )
				-( Fluid.dht )
			)*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	b = where(Fluid.contact == 0, b, 0);

	/*
	Fluid.ReyVals(0,all,all) =	(-( (operatingslippergap.tvr+operatingslippergap.bvr)/2.0 * (hn-hs)/Fluid.dr ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(1,all,all) =	(-( (operatingslippergap.tvtheta+operatingslippergap.bvtheta)/2.0 * (he-hw)/(Fluid.r*Fluid.dtheta) ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(2,all,all) =	( operatingslippergap.tvr * (htn-hts)/Fluid.dr )*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(3,all,all) =	( operatingslippergap.tvtheta * (hte-htw)/(Fluid.r*Fluid.dtheta) )*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(4,all,all) =	(-( Fluid.dht ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;

	Fluid.ReyVals(5,all,all) =	hn-hs;
	Fluid.ReyVals(6,all,all) =	he-hw;
	Fluid.ReyVals(7,all,all) =	htn-hts;
	Fluid.ReyVals(8,all,all) =	hte-htw;

	Fluid.ReyVals(9,all,all) =	operatingslippergap.tvr;
	Fluid.ReyVals(10,all,all) =	operatingslippergap.tvtheta;
	*/

	//Remember orignal pressure boundry
	Array<double,2> bndry(Fluid.boundary.copy());

	//Renumber the non-boundary cells with mid
	Array<int,2> mid(Fluid.M,Fluid.N);
	mid = -1;
	int cells = 0;
	{
		int id=0;

		for(int i=0; i<Fluid.M; i++)
		{
			for(int j=0; j<Fluid.N; j++)
			{
				if(Fluid.contact(i,j) == 1)
				{
					Fluid.boundary(i,j) = Fluid.p(i,j);
				}

				if(Fluid.boundary(i,j) == -1)
				{
					mid(i,j) = id;
					id++;
					cells++;
				}
			}
		}
	}

	//The A Matrix
	gmm::col_matrix<gmm::wsvector<double> > A(cells,cells);

	//The B & X  Vector
	vector<double> B(cells),X(cells);

	//Copy the initial P to X
	for(int i=0;i<M-2;i++)
	{
		for(int j=0;j<N;j++)	
		{
			if(Fluid.boundary(i,j) == -1)
			{
				const int row = mid(i,j);
				X[row] = Fluid.p(i,j);
			}
		}
	}

	//Build the K array inside
	for(int i=0;i<M;i++)
	{
		for(int j=0;j<N;j++)	
		{
			if(Fluid.boundary(i,j) == -1)
			{
				int x,y;

				if (j==0)
				{
					x = N - 1;
				} 
				else 
				{ 
					x = j - 1;
				}
				if (j==N-1)
				{
					y = 0;
				} 
				else 
				{ 
					y = j + 1;
				}

				/*
				an(i,j) * pnew(i+1,j)
				as(i,j) * pnew(i-1,j) 
				ae(i,j) * pnew(i,y) 
				aw(i,j) * pnew(i,x) 
				b(i,j)
				ap(i,j)
				*/

				const int row = mid(i,j);

				B[row] = b(i,j);

				if(mid(i-1,j) >= 0)
				{
					A(row,mid(i-1,j)) = -as(i,j);
				} else {
					//boundary
					B[row] += as(i,j)*Fluid.p(i-1,j);
				}

				if(mid(i+1,j) >= 0)
				{
					A(row,mid(i+1,j)) = -an(i,j);
				} else {
					//boundary
					B[row] += an(i,j)*Fluid.p(i+1,j);
				}

				if(mid(i,y) >= 0)
				{
					A(row,mid(i,y)) = -ae(i,j);
				} else {
					//boundary
					B[row] += ae(i,j)*Fluid.p(i,y);
				}

				if(mid(i,x) >= 0)
				{
					A(row,mid(i,x)) = -aw(i,j);
				} else {
					//boundary
					B[row] += aw(i,j)*Fluid.p(i,x);
				}

				//A(row,row) = ap(i,j);

				A(row,row) = ap(i,j)/alphaPold;
				B[row] += (1-alphaPold)*ap(i,j)*Fluid.p(i,j)/alphaPold;
			}
		}
	}

	//Set up the proper sparse matrix
	gmm::csc_matrix<double> AS;
	gmm::copy(A,AS);

	//Set solver params
	gmm::iteration iter(1e-9);
	iter.set_maxiter(5000);
	iter.set_noisy(0);
	
	//Precondition
	gmm::ildlt_precond< gmm::csc_matrix<double> > PR(AS);
	//gmm::diagonal_precond< gmm::csc_matrix<double> > PR(AS);

	//Solve
	gmm::bicgstab(AS, X, B, PR, iter);

	if(iter.get_iteration() == iter.get_maxiter())
	{
		GapLog.message("WARNING: Reynolds pressure loop failed to converge!");
	}

	//Copy the result back to the Fluid.p
	for(int i=0;i<M;i++)
	{
		for(int j=0;j<N;j++)	
		{
			if(Fluid.boundary(i,j) == -1)
			{
				const int row = mid(i,j);
				Fluid.p(i,j) = X[row];
				Fluid.p_uncut(i,j) = X[row];
			}
		}
	}

	//Fluid.p = Fluid.p;

	//Lower pressure limit
	Fluid.p = where(Fluid.p<0.2e+5,0.2e+5,Fluid.p);

	//Higher pressure limit
	//This is only used in the Visco / Density models
	Fluid.p = where(Fluid.p>10000.0e+5,10000.0e+5,Fluid.p);

	//Restore the old boundary
	Fluid.boundary = bndry;

	return (int) iter.get_iteration();
}
int CSlipperGap::SlipperReynoldsGS(double alphaPold)		//SlipperGap Reynolds Equation
{
	using namespace blitz::tensor; 
	int M = Fluid.M;
	int N = Fluid.N;
	
	Array<double,2> mu(M,N);		mu = 0;			//Actual visco
	
	Array<double,2> hn(M,N);		hn = 0;			//Gap height, north face
	Array<double,2> hs(M,N);		hs = 0;			//Gap height, south face
	Array<double,2> he(M,N);		he = 0;			//Gap height, east face
	Array<double,2> hw(M,N);		hw = 0;			//Gap height, west face

	Array<double,2> htn(M,N);		htn = 0;			//Gap height top, north face
	Array<double,2> hts(M,N);		hts = 0;			//Gap height top, south face
	Array<double,2> hte(M,N);		hte = 0;			//Gap height top, east face
	Array<double,2> htw(M,N);		htw = 0;			//Gap height top, west face

	Array<double,2> hbn(M,N);		hbn = 0;			//Gap height bottom, north face
	Array<double,2> hbs(M,N);		hbs = 0;			//Gap height bottom, south face
	Array<double,2> hbe(M,N);		hbe = 0;			//Gap height bottom, east face
	Array<double,2> hbw(M,N);		hbw = 0;			//Gap height bottom, west face

	Array<double,2> mun(M,N);		mun = 0;			//Visco, north face
	Array<double,2> mus(M,N);		mus = 0;			//Visco, south face
	Array<double,2> mue(M,N);		mue = 0;			//Visco, east face
	Array<double,2> muw(M,N);		muw = 0;			//Visco, west face

	Array<double,2> an(M,N);		an = 0;			//Finite volume coefficent an
	Array<double,2> as(M,N);		as = 0;			//Finite volume coeffcient as
	Array<double,2> ae(M,N);		ae = 0;			//Finite volume coefficent ae
	Array<double,2> aw(M,N);		aw = 0;			//Finite volume coefficent aw
	Array<double,2> ap(M,N);		ap = 0;			//Finite volume coefficient ap
	Array<double,2> b(M,N);			b = 0;			//Reynolds source term
	
	//Linaer interpolate cell height to cell faces
	cell2face(Fluid.h, he, hw, hn, hs);

	//Average viscosity over gap height
	mu = mean(Fluid.oilviscosity,tensor::k);

	//Linaer interpolate cell visco to cell faces
	cell2face(mu, mue, muw, mun, mus);

	//Linaer interpolate top gap height to cell faces
	Array<double,2> hT (Fluid.h + swashplate.ehd);
	cell2face(hT, hte, htw, htn, hts);
	
	//Linaer interpolate bottom gap height to cell faces
	Array<double,2> hB (swashplate.ehd.copy());
	cell2face(hB, hbe, hbw, hbn, hbs);

	//
	Array<double,2> vr((operatingslippergap.tvr+operatingslippergap.bvr)/2.0);
	Array<double,2> vrn(M,N);
	Array<double,2> vrs(M,N);
	Array<double,2> vre(M,N);
	Array<double,2> vrw(M,N);
	cell2face(vr, vre, vrw, vrn, vrs);

	Array<double,2> vtheta((operatingslippergap.tvtheta+operatingslippergap.bvtheta)/2.0);
	Array<double,2> vthetan(M,N);
	Array<double,2> vthetas(M,N);
	Array<double,2> vthetae(M,N);
	Array<double,2> vthetaw(M,N);
	cell2face(vtheta, vthetae, vthetaw, vthetan, vthetas);


	//Spatial coefficients
	if (gapinput->options_slipper.general.reynolds_mu == 2)
	{
		//Full Mu
		ae = (pow(he,3.0)*Fluid.dr)/(mue*Fluid.r * (Fluid.dtheta+East(Fluid.dtheta))/2.0 );
		aw = (pow(hw,3.0)*Fluid.dr)/(muw*Fluid.r * (Fluid.dtheta+West(Fluid.dtheta))/2.0 );
		an = (pow(hn,3.0)*(Fluid.r+Fluid.dr/2.0)*Fluid.dtheta)/(mun * (Fluid.dr+North(Fluid.dr))/2.0 );
		as = (pow(hs,3.0)*(Fluid.r-Fluid.dr/2.0)*Fluid.dtheta)/(mus * (Fluid.dr+South(Fluid.dr))/2.0 );
		ap = an + as + ae + aw;
	} else {
		//Averaged mu
		ae = (pow(he,3.0)*Fluid.dr)/(mu*Fluid.r * (Fluid.dtheta+East(Fluid.dtheta))/2.0 );
		aw = (pow(hw,3.0)*Fluid.dr)/(mu*Fluid.r * (Fluid.dtheta+West(Fluid.dtheta))/2.0 );
		an = (pow(hn,3.0)*(Fluid.r+Fluid.dr/2.0)*Fluid.dtheta)/(mu * (Fluid.dr+North(Fluid.dr))/2.0 );
		as = (pow(hs,3.0)*(Fluid.r-Fluid.dr/2.0)*Fluid.dtheta)/(mu * (Fluid.dr+South(Fluid.dr))/2.0 );
		ap = an + as + ae + aw;
	}

	//Source Term
	b =	(
				-( (operatingslippergap.tvr+operatingslippergap.bvr)/2.0 * (hn-hs)/Fluid.dr )
				-( (operatingslippergap.tvtheta+operatingslippergap.bvtheta)/2.0 * (he-hw)/(Fluid.r*Fluid.dtheta) )
				-( (operatingslippergap.tvr+operatingslippergap.bvr)/2.0 * Fluid.h/Fluid.r)
				-( Fluid.h * (vrn-vrs)/Fluid.dr )
				-( Fluid.h * (vthetae-vthetaw)/(Fluid.r*Fluid.dtheta) )
				+( operatingslippergap.tvr * (htn-hts)/Fluid.dr )
				+( operatingslippergap.tvtheta * (hte-htw)/(Fluid.r*Fluid.dtheta) )
				-( operatingslippergap.bvr * (hbn-hbs)/Fluid.dr )
				-( operatingslippergap.bvtheta * (hbe-hbw)/(Fluid.r*Fluid.dtheta) )
				-( Fluid.dht )
			)*12*Fluid.r*Fluid.dr*Fluid.dtheta;

	if(curFullLoop <= 0 && gapinput->options_slipper.general.EHDsqueeze && gapinput->options_slipper.general.SlipperPressureDeformation)
	{
		//double IMval = 9.07295e-008/100e+5;
		//Array<double, 2> sqzfactor(IMval/gapinput->common.timestep*12*Fluid.r*Fluid.dr*Fluid.dtheta);
		Array<double, 2> sqzfactor(Fluid.newt_ehdsqzapprox(tensor::i)/gapinput->common.timestep*12*Fluid.r*Fluid.dr*Fluid.dtheta);
		ap += sqzfactor;
		b += Fluid.pFullLoop*sqzfactor;
	}

	//const double wearHeight = 2.0*Fluid.ContactHeight;
	//b = where(Fluid.h <= wearHeight, 0, b);
	//b = where(Fluid.contact == 0, b, 0);
	/*
	Fluid.ReyVals(0,all,all) = (-( (operatingslippergap.tvr+operatingslippergap.bvr)/2.0 * (hn-hs)/Fluid.dr ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(1,all,all) = (-( (operatingslippergap.tvtheta+operatingslippergap.bvtheta)/2.0 * (he-hw)/(Fluid.r*Fluid.dtheta) ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(2,all,all) = (( operatingslippergap.tvr * (htn-hts)/Fluid.dr ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;;
	Fluid.ReyVals(3,all,all) = ( operatingslippergap.tvtheta * (hte-htw)/(Fluid.r*Fluid.dtheta) )*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(4,all,all) = (-( operatingslippergap.bvr * (hbn-hbs)/Fluid.dr ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(5,all,all) = (-( operatingslippergap.bvtheta * (hbe-hbw)/(Fluid.r*Fluid.dtheta) ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(6,all,all) = (-( Fluid.dht ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(7,all,all) = ae;
	Fluid.ReyVals(8,all,all) = aw;
	Fluid.ReyVals(9,all,all) = an;
	Fluid.ReyVals(10,all,all) = as;
	
	Fluid.ReyVals(0,all,all) = where(Fluid.h <= wearHeight, 0, Fluid.ReyVals(0,all,all));
	Fluid.ReyVals(1,all,all) = where(Fluid.h <= wearHeight, 0, Fluid.ReyVals(1,all,all));
	Fluid.ReyVals(2,all,all) = where(Fluid.h <= wearHeight, 0, Fluid.ReyVals(2,all,all));
	Fluid.ReyVals(3,all,all) = where(Fluid.h <= wearHeight, 0, Fluid.ReyVals(3,all,all));
	Fluid.ReyVals(4,all,all) = where(Fluid.h <= wearHeight, 0, Fluid.ReyVals(4,all,all));
	Fluid.ReyVals(5,all,all) = where(Fluid.h <= wearHeight, 0, Fluid.ReyVals(5,all,all));
	Fluid.ReyVals(6,all,all) = where(Fluid.h <= wearHeight, 0, Fluid.ReyVals(6,all,all));
	*/
	/*
	Fluid.ReyVals(0,all,all) =	(-( (operatingslippergap.tvr+operatingslippergap.bvr)/2.0 * (hn-hs)/Fluid.dr ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(1,all,all) =	(-( (operatingslippergap.tvtheta+operatingslippergap.bvtheta)/2.0 * (he-hw)/(Fluid.r*Fluid.dtheta) ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(2,all,all) =	( operatingslippergap.tvr * (htn-hts)/Fluid.dr )*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(3,all,all) =	( operatingslippergap.tvtheta * (hte-htw)/(Fluid.r*Fluid.dtheta) )*12*Fluid.r*Fluid.dr*Fluid.dtheta;
	Fluid.ReyVals(4,all,all) =	(-( Fluid.dht ))*12*Fluid.r*Fluid.dr*Fluid.dtheta;

	Fluid.ReyVals(5,all,all) =	hn-hs;
	Fluid.ReyVals(6,all,all) =	he-hw;
	Fluid.ReyVals(7,all,all) =	htn-hts;
	Fluid.ReyVals(8,all,all) =	hte-htw;

	Fluid.ReyVals(9,all,all) =	operatingslippergap.tvr;
	Fluid.ReyVals(10,all,all) =	operatingslippergap.tvtheta;
	*/
	//Fluid.ReyVals(all,all,all) = where(Fluid.contact(tensor::j,tensor::k) == 0, Fluid.ReyVals(tensor::i,tensor::j,tensor::k), 0);

	//Renumber the non-boundary cells with mid
	Array<int,2> mid(Fluid.M,Fluid.N);
	mid = -1;
	int cells = 0;
	{
		int id=0;

		for(int i=0; i<Fluid.M; i++)
		{
			for(int j=0; j<Fluid.N; j++)
			{
				if(Fluid.boundary(i,j) == -1)
				{
					mid(i,j) = id;
					id++;
					cells++;
				}
			}
		}
	}

	double alpha = gapinput->options_slipper.numeric.AlphaReynolds;
	int iterations = 0;
	int x,y;
	double error,errornew;
	Array<double, 2> pnew(Fluid.p.copy());
	Array<double, 2> pold(Fluid.p.copy());

	//Gauss seidel
	do
	{
		error = 0.0;
		errornew = 0.0;

		//SOR
		for(int i=1;i<M-1;i++)
		{
			for(int j=0;j<N;j++)	
			{
				if(Fluid.boundary(i,j) == -1)
				{
					if (j==0)
					{
						x = N - 1;
					} 
					else 
					{ 
						x = j - 1;
					}
					if (j==N-1)
					{
						y = 0;
					} 
					else 
					{ 
						y = j + 1;
					}

					/*
					pnew(i,j) = pnew(i,j) + alpha * (  
									(	an(i,j)*pnew(i+1,j) +
										as(i,j)*pnew(i-1,j) + 
										ae(i,j)*pnew(i,y) + 
										aw(i,j)*pnew(i,x) + b(i,j) ) / ap(i,j) 
									- pnew(i,j) );
					*/
					
					
					pnew(i,j) = pnew(i,j) + alpha * (  
									(	an(i,j)*pnew(i+1,j) +
										as(i,j)*pnew(i-1,j) + 
										ae(i,j)*pnew(i,y) + 
										aw(i,j)*pnew(i,x) + b(i,j) + (1-alphaPold)*ap(i,j)*Fluid.p(i,j)/alphaPold ) * alphaPold / ap(i,j) 
									- pnew(i,j) );
					
					
					Fluid.p_uncut(i,j) = pnew(i,j);
 					
					//Lower pressure limit
					if (pnew(i,j) < 0.8e+5)
					{
						pnew(i,j) = 0.8e+5;
					}
					//Higher pressure limit
					if (pnew(i,j) > 2.0e8)
					{
						pnew(i,j) = 2.0e8;
					}

					//Error
					errornew = fabs(pnew(i,j) - pold(i,j));
					pold(i,j) = pnew(i,j);
					if (errornew > error)
					{
						error = errornew;
					}
				}
			}
		}

		//Counter
		iterations++;
	
	}while(error > 0.01 && iterations < 1e+5);
	
	if(error > 0.01 && iterations >= 1e+5 && gapinput->options_slipper.general.DebugMode != 1)
	//if(error > 0.01 && iterations >= 1e4 && GapInput.SlipperOptions.SlipperDebugMode != 1)
	{
		GapLog.message("WARNING: Reynolds pressure loop failed to converge!");
	}

	//Assign new pressure to p field
	Fluid.p = pnew;

	return iterations;

}
void CSlipperGap::cell2face(Array<double,2> & p, Array<double,2> & pe, Array<double,2> & pw, Array<double,2> & pn, Array<double,2> & ps)
{
	const int M = Fluid.M;
	const int N = Fluid.N;

	//pn = (p*drp+pn*drn)/(drp+drn)
	pn(Range(0,M-2),all) = 
		(p(Range(0,M-2),all)*Fluid.dr(Range(0,M-2),all) +					//(p*drp +
		p(Range(1,M-1),all)*Fluid.dr(Range(1,M-1),all))						//pn*drn)
		/(Fluid.dr(Range(0,M-2),all)+Fluid.dr(Range(1,M-1),all));		// /(drp+drn);
	//"upwind" the boundaries
	pn(Range(0,M-2),all) = where(Fluid.boundary(Range(1,M-1),all) == -1, pn(Range(0,M-2),all), p(Range(0,M-2),all));
	pn(Range(0,M-2),all) = where(Fluid.hgroove(Range(1,M-1),all) == 0, pn(Range(0,M-2),all), p(Range(0,M-2),all));
	pn(M-1,all) = p(M-1,all);
	

	//ps = (p*drp+ps*drs)/(drp+drs)
	ps(Range(1,M-1),all) = 
		(p(Range(1,M-1),all)*Fluid.dr(Range(1,M-1),all) +					//(p*drp +
		p(Range(0,M-2),all)*Fluid.dr(Range(0,M-2),all))						//ps*drs)
		/(Fluid.dr(Range(1,M-1),all)+Fluid.dr(Range(0,M-2),all));		// /(drp+drs);
	//"upwind" the boundarie
	ps(Range(1,M-1),all) = where(Fluid.boundary(Range(0,M-2),all) == -1, ps(Range(1,M-1),all), p(Range(1,M-1),all));
	ps(Range(1,M-1),all) = where(Fluid.hgroove(Range(0,M-2),all) == 0, ps(Range(1,M-1),all), p(Range(1,M-1),all));
	ps(0,all) = p(0,all);
	
	//pe = (p*dtp+pe*dte)/(dtp+dte)
	pe(all,Range(0,N-2)) =
		(p(all,Range(0,N-2))*Fluid.dtheta(all,Range(0,N-2)) +						//(p*dtp +
		p(all,Range(1,N-1))*Fluid.dtheta(all,Range(1,N-1)))						//pe*dte)
		/(Fluid.dtheta(all,Range(0,N-2))+Fluid.dtheta(all,Range(1,N-1)));		// /(dtp+dte);
	pe(all,Range(0,N-2)) = where(Fluid.boundary(all,Range(1,N-1)) == -1, pe(all,Range(0,N-2)), p(all,Range(0,N-2)));
	pe(all,Range(0,N-2)) = where(Fluid.hgroove(all,Range(1,N-1)) == 0, pe(all,Range(0,N-2)), p(all,Range(0,N-2)));
	
	pe(all,N-1) =
		(p(all,N-1)*Fluid.dtheta(all,N-1) +												//(p*dtp +
		p(all,0)*Fluid.dtheta(all,0))														//pe*dte)
		/(Fluid.dtheta(all,N-1)+Fluid.dtheta(all,0));								// /(dtp+dte);
	pe(all,N-1) = where(Fluid.boundary(all,0) == -1, pe(all,N-1), p(all,N-1));
	pe(all,N-1) = where(Fluid.hgroove(all,0) == 0, pe(all,N-1), p(all,N-1));

	//pw = (p*dtp+pw*dtw)/(dtp+dtw)
	pw(all,Range(1,N-1)) =
		(p(all,Range(1,N-1))*Fluid.dtheta(all,Range(1,N-1)) +						//(p*dtp +
		p(all,Range(0,N-2))*Fluid.dtheta(all,Range(0,N-2)))						//pw*dtw)
		/(Fluid.dtheta(all,Range(1,N-1))+Fluid.dtheta(all,Range(0,N-2)));		// /(dtp+dtw);
	pw(all,Range(1,N-1)) = where(Fluid.boundary(all,Range(0,N-2)) == -1, pw(all,Range(1,N-1)), p(all,Range(1,N-1)));
	pw(all,Range(1,N-1)) = where(Fluid.hgroove(all,Range(0,N-2)) == 0, pw(all,Range(1,N-1)), p(all,Range(1,N-1)));

	pw(all,0) =
		(p(all,0)*Fluid.dtheta(all,0) +													//(p*dtp +
		p(all,N-1)*Fluid.dtheta(all,N-1))												//pw*dtw)
		/(Fluid.dtheta(all,0)+Fluid.dtheta(all,N-1));								// /(dtp+dtw);
	pw(all,0) = where(Fluid.boundary(all,N-1) == -1, pw(all,0), p(all,0));
	pw(all,0) = where(Fluid.hgroove(all,N-1) == 0, pw(all,0), p(all,0));
}
Array<double,2> CSlipperGap::North(Array<double,2> & phi)
{
	blitz::Range all = Range::all();
	Array<double,2> n(Fluid.M, Fluid.N);
	n(Range(0,Fluid.M-2),all) = phi(Range(1,Fluid.M-1),all);

	//upwind the boundary
	n(Fluid.M-1,all) = phi(Fluid.M-1,all);
	
	return n;
}
Array<double,2> CSlipperGap::South(Array<double,2> & phi)
{
	blitz::Range all = Range::all();
	Array<double,2> s(Fluid.M, Fluid.N);
	s(Range(1,Fluid.M-1),all) = phi(Range(0,Fluid.M-2),all);

	//upwind the boundary
	s(0,all) = phi(0,all);
	
	return s;
}
Array<double,2> CSlipperGap::East(Array<double,2> & phi)
{
	blitz::Range all = Range::all();
	Array<double,2> e(Fluid.M, Fluid.N);
	e(all,Range(0,Fluid.N-2)) = phi(all,Range(1,Fluid.N-1));

	//"connect" the boundary
	e(all,Fluid.N-1) = phi(all,0);
	
	return e;
}
Array<double,2> CSlipperGap::West(Array<double,2> & phi)
{
	blitz::Range all = Range::all();
	Array<double,2> w(Fluid.M, Fluid.N);
	w(all,Range(1,Fluid.N-1)) = phi(all,Range(0,Fluid.N-2));

	//"connect" the boundary
	w(all,0) = phi(all,Fluid.N-1);
	
	return w;
}
