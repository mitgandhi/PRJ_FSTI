#include "CSlipperGap.h"
#include <omp.h>
#include <time.h>
#include <deque>
#pragma once
#define OMP

extern void matlab(const Array<double,2>& data,const string file);

void CSlipperGap::testDeform(void)
{
	GapLog << endl << "Fluid.newt_ehdsqzapprox = " << endl;
	for(int i=0; i<Fluid.M; i++)
	{
		GapLog << Fluid.newt_ehdsqzapprox(i) << endl;
	}
	exit(1);



	//we want to test simple deformations
	for(int i=1; i<28; i++)
	{
		operatingslippergap.pG = 0;
		operatingslippergap.pcase = 0;
		forcesslippergap.psocket = 0;
		Fluid.p = 0;

		Fluid.p(Range(i-1,i+1),Range(59,61)) = 100e5;
		SlipperCalcEHD(1);

		matlab(Fluid.p, "p_pm1_"+n2s(i)+".txt");
		matlab(slipper.ehd, "slip_pm1_"+n2s(i)+".txt");
	}


	exit(1);

	swashplate.ehd = 0;
	operatingslippergap.pG = 300e5;
	SlipperPressureBounds();

	int h;

	h = 2;
	Fluid.h = double(h)*1e-6;
	SlipperReynolds(1.0);
	SlipperCalcFluidV();
	SlipperCalcDensity();
	SlipperEnergy();
	writevtk("./"+n2s(h)+".vtk");

	h = 14;
	Fluid.h = double(h)*1e-6;
	SlipperReynolds(1.0);
	SlipperCalcFluidV();
	SlipperCalcDensity();
	SlipperEnergy();
	writevtk("./"+n2s(h)+".vtk");
	
	exit(1);


	/*
	SlipperCalch();
	SlipperCalcdht();
	for(double pg=0.1e+5; pg<5e+5; pg += 0.1e+5)
	{
		operatingslippergap.pG = pg;
		SlipperPressureBounds();
		for(int i=0; i<10; i++)
		{
			SlipperReynolds(0.4);
			SlipperEnergy();
			SlipperCalcViscosity(1.0);
			SlipperCalcEHD(1.0);
		}
		SlipperCalcFluidV();
		cout << pg << "\t" << operatingslippergap.QSG << endl;
	}

	exit(1111);
	*/

/*
	double Rg = geometryslippergap.routG;
	Array<double,2> A(3,3);
	TinyVector<double,3> b;
	TinyVector<double,3> x;

	//Create geometry matrix
	A = 1,1,1,Rg,-.5*Rg,-.5*Rg,0,-sqrt(3.0)/2*Rg,sqrt(3.0)/2*Rg;


		Array<double,2> hT (Fluid.h + swashplate.ehd);
		Array<double,2> contact(hT - geometryslippergap.hmaxG);
		contact = where(contact < 0, 0, contact);
		contact(Range(0,Fluid.M-2),Range::all()) = 0;
		
		contact = 0;
		contact(Fluid.M-1,0) = 1e-6;
		contact(Fluid.M-1, Fluid.N-1) = 1e-6;
		//contact(Fluid.M-1, 1) = 1e-6;

		Array<double,2> contactp ( contact*geometryslippergap.Fslipper / (double)Fluid.N );


		double Fz = sum(contactp);
		double Mx = sum(contactp*Fluid.Ly);
		double My = sum(-contactp*Fluid.Lx);

		//b = Fz,Mx,My;
		b(0) = Fz;
		b(1) = Mx;
		b(2) = My;

		cout << A << endl;
		
		for(int i=0; i<3; i++)
		{
			cout << b(i) << "\t";
		}
		cout << endl;

		x = Solve3(A,b);

		for(int i=0; i<3; i++)
		{
			cout << x(i) << "\t";
		}


		//for debugging
		
		matlab(contactp,"./contactp.txt");




	exit(1987);
*/

	/*
	operatingslippergap.pDC = 100e+5;
	operatingslippergap.pG = operatingslippergap.pDC;
	operatingslippergap.pcase = 1e+5;

	double minh = 5e-6;
	double dh = 5e-6;
	double slope = dh/(max(Fluid.Lx)-min(Fluid.Lx));
	Fluid.h = slope*(Fluid.Lx-min(Fluid.Lx))+minh;
	
	slipper.ehd = 0;
	swashplate.ehd = 0;
	Fluid.dht = 0;

	Fluid.oilviscosity = 0.02;
	SlipperPressureBounds();
	SlipperReynolds(1.0);

	writevtk("./fluid.vtk");

	matlab(Fluid.p, "p.txt");
	matlab(operatingslippergap.vgx, "vgx.txt");
	matlab(operatingslippergap.vgy, "vgy.txt");


	exit(1987);
	*/
	/*
	//THIS SECTION IS USED TO GENERATE A FILEPRESSURE.TXT FILE THAT CAN BE USED WITH INFLUGEN
	//IT GENERATES A HYDROSTATIC PRESSURE FIELD AND THE ASSOCIATED POCKET / SOCKET PRESSURES

		operatingslippergap.pDC = 300e+5;
		operatingslippergap.pG = operatingslippergap.pDC;
		operatingslippergap.pcase = 1e+5;

		Fluid.h = 5e-6;
		Fluid.dht = 0;
		slipper.ehd = 0;
		swashplate.ehd = 0;

		GetFSK();	//update the socket pressure
		
		Fluid.oilviscosity = 0.02;
		SlipperPressureBounds();
		SlipperReynolds(1.0);

		ofstream fp("filepressure.txt");
		fp << "p=" << operatingslippergap.pG << endl;
		fp << "s=" << forcesslippergap.psocket << endl;
		for(int i=0; i<slipper.facePressure.size(); i++)
		{
			fp << slipper.facePressure[i] << endl;
		}
		fp.close();
		writevtk("./fluid.vtk");
		exit(1988);
	//
	//
	*/

/*
	fstream fp("./fp.bin", ios::out | ios::binary);
	for(int i=0; i<swashplate.facePressure.size(); i++)
	{
		fp.write((char*) &(swashplate.facePressure[i]), sizeof(double));
	}
	fp.close();

	fstream nd("./nd.bin", ios::out | ios::binary);
	for(int i=0; i<swashplate.nodeDeformation.size(); i++)
	{
		nd.write((char*) &(swashplate.nodeDeformation[i]), sizeof(double));
	}
	nd.close();
	*/

	ofstream fp("fp.txt");
	for(int i=0; i<swashplate.facePressure.size(); i++)
	{
		fp << setprecision(20) << swashplate.facePressure[i] << endl;
	}
	fp.close();

	ofstream nd("nd.txt");
	for(int i=0; i<swashplate.nodeDeformation.size(); i++)
	{
		nd << setprecision(20) <<  swashplate.nodeDeformation[i] << endl;
	}
	nd.close();

	exit(1);

	swashplate.writevtk("0.vtk");
	//SlipperReynolds(1.0);
	Fluid.p = 200e+5;
	SwashplateCalcEHD(1.0);
	swashplate.writevtk("1.vtk");

	operatingslippergap.phi_deg = 10.0; //[deg] between 0 and 360
	operatingslippergap.phi_rad = operatingslippergap.phi_deg*PI/180; //[rad] between 0 and 2PI
	SlipperSetCoordSys();
	{
		//Because the phi angle has changed we need to update which slippers are in the FullGroup

		int steps = int(floor(360.0/gapinput->common.phistep_deg));
		int step = int(floor(operatingslippergap.phi_deg/gapinput->common.phistep_deg+0.5)) % steps;

		for(int np = 0, cnt = 0; np<geometryslippergap.npistons; np++)
		{
			int s = (step + np*steps/geometryslippergap.npistons) % steps;
			for(int i=0; i<Fluid.M; i++)
			{
				for(int j=0; j<Fluid.N; j++, cnt++)
				{
					FullGroup.KDslip.points[cnt][0] = FullGroup.Gx(i,j,s);
					FullGroup.KDslip.points[cnt][1] = FullGroup.Gy(i,j,s);
				}
			}
		}

		if(FullGroup.KDslip.kdtree == NULL)
		{
			FullGroup.KDslip.kdtree = new ANNkd_tree(FullGroup.KDslip.points, Fluid.M*Fluid.N*geometryslippergap.npistons, FullGroup.KDslip.dim);
		} else {
			delete FullGroup.KDslip.kdtree;
			FullGroup.KDslip.kdtree = new ANNkd_tree(FullGroup.KDslip.points, Fluid.M*Fluid.N*geometryslippergap.npistons, FullGroup.KDslip.dim);
		}
	}
	SwashplateInitEHD(1.0);

	swashplate.writevtk("2.vtk");

	exit(0);


}

void CSlipperGap::testDeform(vector<double> &xg,vector<double> &vg)
{
	double vorg = vg[1];
	for(double dv=-0.004; dv<0.004; dv+=0.0005)
	{
		vg[1] = vorg+dv;

		for(int i = 0; i < 3; i++)
		{
			control_point_velocity[i] = vg[i];
			xg[i] = pold[i]+gapinput->common.timestep*(vold[i]+vg[i])/2.0;
		}
		SlipperCalcdht();
		SlipperCalch(xg);
		Fluid.dht = 0;
		SlipperPressureBounds();
		
		SlipperReynolds(1.0);
		SlipperCalcFluidV();
		SlipperCalcDensity();
		SlipperEnergy();
		SlipperpG(1);
		SlipperReynolds(1.0);
		SlipperCalcFluidForces();
		SlipperCalcExternalForces();
		vector<double> dF(3,0);
		SlipperCalcdF(dF, xg);
	
		ofstream f("lt.txt", ios::app);
		f << vg[0] << "\t" << vg[1] << "\t" << vg[2] << "\t" << dF[0] << "\t" << dF[1] << "\t" << dF[2] << "\t" <<
				GapResult->SlipperResult.dFcomp[0] << "\t" << GapResult->SlipperResult.dFcomp[1] << "\t" << GapResult->SlipperResult.dFcomp[2] << endl;

		f.close();

		writevtk("./vtk/"+n2s(vg[0])+"$"+n2s(vg[1])+"$"+n2s(vg[2])+".vtk");
	}
	
	exit(1);

	Fluid.dht = 0;

	//we need to setup a stable ehd field
	{
		Fluid.h = 5e-6;
		SlipperPressureBounds();
		SlipperCalcEHD(1.0);
	}
	

	double lift = -6.0e-6;
	double maxtilt = 2.0e-6;
	
	for(double tilt = -maxtilt; tilt <= maxtilt; tilt += 1.0e-6)
	{
		xg[0] = lift;
		xg[1] = lift+tilt;
		xg[2] = lift-tilt;


		SlipperCalch(xg);

		GetFSK();

		SlipperReynolds(1.0);
		SlipperCalcFluidV();
		SlipperCalcDensity();
		SlipperEnergy();
		SlipperpG(1);
		SlipperReynolds(1.0);
		SlipperCalcFluidForces();
		SlipperCalcExternalForces();
		vector<double> dF(3,0);
		SlipperCalcdF(dF, xg);

		ofstream f("lt.txt", ios::app);
		f << lift << "\t" << tilt << "\t" << dF[0] << "\t" << dF[1] << "\t" << dF[2] << "\t" <<
				GapResult->SlipperResult.dFcomp[0] << "\t" << GapResult->SlipperResult.dFcomp[1] << "\t" << GapResult->SlipperResult.dFcomp[2] << endl;

		f.close();

		writevtk("fluid$"+n2s(lift)+"$"+n2s(tilt)+".vtk");

	}



	exit(1987);


	/*
	xg[0] = 5e-6;
	xg[1] = 5e-6;
	xg[2] = 5e-6;
*/

	//vg[0] = 0;
	//vg[1] = 0;
	//vg[2] = 0;


}

void CSlipperGap::debug_makeEHDfluidpressurefile()
{
	swashplate.ehd = 0;
	operatingslippergap.pG = 300e5;
	operatingslippergap.pcase = 1e5;
	SlipperPressureBounds();

	operatingslippergap.pDC = operatingslippergap.pG;
	GetFSK();

	int h;

	h = 2;
	Fluid.h = double(h)*1e-6;

	if(gapinput->options_slipper.general.EnableSlipperMacro == 1)
	{
		Fluid.h += geometryslippergap.SlipperMacro;
	}

	SlipperReynolds(1.0);
	//SlipperCalcFluidV();
	//SlipperCalcDensity();
	//SlipperEnergy();

	//SlipperCalcEHD(1.0);

	SlipperCalcFluidForces();

	cout << "forcesslippergap.FfGz = " << forcesslippergap.FfGz << endl;
	cout << "forcesslippergap.MfGx = " << forcesslippergap.MfGx << endl;
	cout << "forcesslippergap.MfGy = " << forcesslippergap.MfGy << endl;
	cout << "forcesslippergap.FSK = " << forcesslippergap.FSK << endl;
	
	
	writevtk("./"+n2s(h)+".vtk");

	
	exit(1987);
}


void matlab(const Array<double,2>& data,const string file)
{
	ofstream f;
	f.open(file.c_str());
	for(int i=0;i<data.extent(0);i++)
	{
		for(int j=0;j<data.extent(1);j++)
		{
			f << setprecision (10) << data(i,j)  << " ";
		}
		f << endl;
	}
	f.close();

}




void CSlipperGap::testGap(void)
{
	double resid = 0;

	double alphaPold = 1.00;
	bool FullLoop = true;
	int slippercnt=0;
	int MaxpdequeSize = 1;

	deque< Array<double,2> > pdeque;
	pdeque.push_front(Fluid.p_uncut.copy());
	
	Array<double,2> pold(Fluid.M,Fluid.N);
	pold = Fluid.p;

	int goodctr = 0;
	int resetcnt = 0;
	do
	{
		double oldresid = resid;
		resid = presid();

		if (sum(where(Fluid.h<Fluid.ContactHeight,1,0)) > 0)
		{
			goodctr = 0;
			resetcnt = slippercnt;

			//Used to 'reset' convergence
			operatingslippergap.pG = operatingslippergap.pDC;
			SlipperPressureBounds();
			SlipperInitEHD();
			SlipperInitializeTemperature();
			SlipperCalcViscosity(1.0);
			pdeque.clear();
			pdeque.push_front(Fluid.p_uncut.copy());        

			
			if(alphaPold < 0.1)
			{
				alphaPold *= 0.9;
			} else {
				alphaPold -= 0.05;
			}
			
			if(alphaPold < 0.01)
			{
				alphaPold = 0.01;
			}

			if(alphaPold < 0.74)
			{
				MaxpdequeSize++;
			}
			
		} else if (slippercnt > 0 && oldresid < resid)
		{
			goodctr = 0;
			
			alphaPold -= 0.05;
			if(alphaPold < 0.1)
			{
				alphaPold = 0.1;
			}

		}

		if(goodctr > 25)
		{
			alphaPold += 0.025;
			if(alphaPold >= 1.0)
			{
				alphaPold = 1.0;
			}	
			goodctr = 0;
		}
		goodctr++;

		/*
		//"Cut the gap" if needed
		Fluid.contact = where(Fluid.h<Fluid.ContactHeight,1,0);
		if(sum(Fluid.contact) > 0)
		{
			GapLog.message("Contact: " + n2s(sum(Fluid.contact)) + " Min h: " + n2s(min(Fluid.h)));
		}
		Fluid.h = where(Fluid.h<Fluid.ContactHeight,Fluid.ContactHeight,Fluid.h);
		*/
		
		//Solve Reynolds
		int nRiterations = SlipperReynolds(alphaPold);

	
		//Calculate contact pressure
		//SlipperCalcContactPressure();

		if(slippercnt-resetcnt < 10 || (slippercnt-resetcnt) % 10 == 0)
		{
			//Update fluid velocities
			SlipperCalcFluidV();
		}

		if(FullLoop)
		{
			if(slippercnt-resetcnt < 10 || (slippercnt-resetcnt) % 10 == 0)
			{
				//Calculate fluid density
				SlipperCalcDensity();

				int nEiterations;
				if(gapinput->options_slipper.general.CalcEnergy)
				{
					//Solve Energy
					nEiterations = SlipperEnergy();
				}

				//Calculate viscosity
				SlipperCalcViscosity(1.0);
			}

			//Calculate pressure deformation
			if(gapinput->options_slipper.general.SlipperPressureDeformation)
			{
				SlipperCalcEHD(1.0);

			}

			if(slippercnt-resetcnt < 10 || (slippercnt-resetcnt) % 10 == 0)
			{
				//Update pocket pressure
				SlipperpG(0.5);
			}

			double perr1 = sqrt(sum(pow((pold - Fluid.p)/Fluid.p,2.0)));
			pold = Fluid.p;

			if(gapinput->options_slipper.general.DebugMode == 1)
			{
				GapLog.message(
					"\tpG: " + n2s(operatingslippergap.pG) + 
					"\tm(p): " + n2s(mean(Fluid.p)) +
//					"\tm(T): " + n2s(mean(Fluid.T)) +
					"\tm(h): " + n2s(mean(Fluid.h)) +
//					"\tm(mu): " + n2s(mean(Fluid.oilviscosity)) +
//					"\tn(p): " + n2s(nRiterations) +
//					"\tn(T): " + n2s(nEiterations) +
					"\tresid: " + n2s(resid) +
					"\talp: " + n2s(alphaPold) +
					"\tperr: " + n2s(perr1)
					);

				cout <<
//					"\tpG: " + n2s(operatingslippergap.pG) + 
//					"\tm(p): " + n2s(mean(Fluid.p)) +
//					"\tm(T): " + n2s(mean(Fluid.T)) +
//					"\tm(h): " + n2s(mean(Fluid.h)) +
//					"\tm(mu): " + n2s(mean(Fluid.oilviscosity)) +
//					"\tn(p): " + n2s(nRiterations) +
//					"\tn(T): " + n2s(nEiterations) +
					"\tresid: " + n2s(resid) +
					"\talp: " + n2s(alphaPold) +
					"\tqsz: " + n2s(MaxpdequeSize) +
					"\tperr: " + n2s(perr1)
					<< endl;
			}

			//matlab(Fluid.p,"./dense/p" + n2s(slippercnt) + ".txt");
			//matlab(Fluid.h,"./dense/h" + n2s(slippercnt) + ".txt");
		}

		//Compute the relative error between iterations
		//perr = max(abs(Fluid.p - pold));
		
		//Increment fixed point itteration loop counter
		slippercnt++;

	} while (resid > 1e-5 && slippercnt < 1000 && FullLoop);	//Convergence criteria

	matlab(Fluid.h,"h.txt");
	matlab(Fluid.p,"p.txt");

}

void CSlipperGap::alphaP(const double init_residual, const double resid, double & alphaPold)
{
	/*
	if(resid <= init_residual)
	{
	
		alphaPold *= 1.01;
		if(alphaPold > 0.95)
		{
			alphaPold = 0.95;
		}
		
	} else {
		alphaPold -= 0.05;
		if(alphaPold < 0.1)
		{
			alphaPold = 0.1;
		}
	}
	*/
	

}
