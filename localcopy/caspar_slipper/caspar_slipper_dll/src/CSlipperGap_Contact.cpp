#include "CSlipperGap.h"
#include <time.h>
#include <deque>
#include "QuadProg++.h"
#pragma once
#define OMP

extern void matlab(const Array<double,2>& data,const string file);

bool CSlipperGap::SlipperCalcContactPressure(const double alphaEHD, const int FullLoop, const double alphaContactP)
{
	//Saturation and relaxation parameters
	double AlphaContact = gapinput->options_slipper.numeric.AlphaContact;  //Reduce contact pressure by this factor
	double contact_p_lim = gapinput->options_slipper.numeric.contact_p_lim;  //Limit contact pressure to this value
	double contact_def_lim = gapinput->options_slipper.numeric.contact_def_lim;  //Limit contact deformation to this value 


	//at twice the contact height, set the p_contact_old to = 0
	const double ContactHeight = 2.0*Fluid.ContactHeight;
	Fluid.p_contact = where(Fluid.h > 2.0*ContactHeight, 0, Fluid.p_contact);
	
	Array<double, 2> p_contact_old(Fluid.p_contact.copy());
	
	//initial checking
	{
		double sumcontact = sum(Fluid.contact);
		if(sumcontact > 250)
		{
			//this is a _large_ contact problem and the QuadProg library used takes forever even though Matlab does it quick.
			//thus it is _possible_, but for now, fuck it.
			// - correction to the above comment, I think the major time is spent building the IMlib needed for a large contact
			
			//found by using the below method to find the average of a small (8 cell) patch on the H1P130
			double Kstiff = 3.6943e+13; //units = Pa/m
			Fluid.p_contact = Fluid.h_contact * Kstiff;

			/*
			//limit pressure
			Fluid.p_contact = where(Fluid.p_contact > 800e5, 800e5, Fluid.p_contact);

			//don't relax the problem
			Fluid.p += Fluid.p_contact;
			*/

			//limit pressure
			Fluid.p_contact = where(Fluid.p_contact > contact_p_lim, contact_p_lim, Fluid.p_contact);

			//don't relax the problem
			Fluid.p += AlphaContact * Fluid.p_contact;

			Fluid.h = where(Fluid.h < ContactHeight , ContactHeight , Fluid.h);

			return true;
		}

		if(sumcontact == 0)
		{
			//no contact so return

			//but we should relax somewhat to prevent contact chattering.. i don't know if this has an impact
			Fluid.p_contact *= alphaContactP;

			return true;
		}
	}
	

	//this should go to the CSlipperGap class.. cheating for now
	static map<pair<int, int>, Array<double,2> > IMlib;

	//The latest way

	//save the pressure and the slipper ehd
	Array<double, 2> p(Fluid.p.shape());
	p = Fluid.p;
	Array<double, 2> slip_ehd(slipper.ehd.shape());
	slip_ehd = slipper.ehd;

	//build a vector of the cells in contact	
	vector<pair<int, int> > contact;

	for(int m=0; m<Fluid.M; m++)
	{
		for(int n=0; n<Fluid.N; n++)
		{
			if(Fluid.contact(m,n) == 1)
			{
				pair<int, int> cell(m,n);

				if(IMlib.count(cell) == 0)
				{
					//we need to add it to the library
					
					//try to use the actual IM's to determine the approprate value
					//this will fail if the IM mesh is too course wrt the fluid mesh
					//due to the simplified interpolation method currently in use
					/*
					Fluid.p = 0;
					Fluid.p(m,n) = 100e5;

					Fluid2Slipper();
					slipper.calcDeform(0, 0, 0);
					Slipper2Fluid(1);
					*/

					//use an approximated method that simply guesses
					//at what the value of penetration should be
					IMlib[cell].resize(slipper.ehd.shape());
					slipper.ehd = 0.1e-9;
					
					for(int i=(m-2<0?0:m-2); i<=(m+2>=Fluid.M?Fluid.M-1:m+2); i++)
					{
						for(int j=Fluid.N+n-2; j<=Fluid.N+n+2; j++)
						{
							slipper.ehd(i, j%Fluid.N) = 1e-9;
						}
					}
					slipper.ehd(m,n) = 3e-9;

					IMlib[cell] = slipper.ehd;
				}

				contact.push_back(cell);
			}
		}
	}

	int contactpts = (int) contact.size();
	
	//build the K matrix / b vector
	Array<double, 2> K(contactpts, contactpts);
	K = 0;

	vector<double> b(contactpts, 0);

	for(int i=0; i<contactpts; i++)
	{
		for(int j=0; j<contactpts; j++)
		{
			K(i,j) = IMlib[contact[i]](contact[j].first, contact[j].second);
		}

		b[i] = Fluid.h_contact(contact[i].first, contact[i].second);
	}

	/*
	matlab(K, "./k.txt");
	//write the b vector
	ofstream bf("./b.txt");
	for(int i=0; i<b.size(); i++)
	{
		bf << b[i] << endl;
	}
	bf.close();
	*/
	
	//setup the quad programming
	QuadProgPP::Matrix<double> G, CE, CI;
	QuadProgPP::Vector<double> g0, ce0, ci0, x;
		
	G.resize(contactpts, contactpts);
	for (int i = 0; i < contactpts; i++)	
	{
		for (int j = 0; j < contactpts; j++)
		{
			if(i == j)
			{
				G[i][j] = 1;
			} else {
				G[i][j] = 0;
			}
		}
	}

	g0.resize(contactpts);
	for (int i = 0; i < contactpts; i++)
	{
		g0[i] = 0;
	}

	CE.resize(contactpts, 0);
	ce0.resize(0);

	CI.resize(contactpts, 2*contactpts);

	//the compliance matrix
	for (int i = 0; i < contactpts; i++)
	{
		for (int j = 0; j < contactpts; j++)
		{
			CI[i][j] = K(j,i)*1e4;
		}
	}

	//the p>0 constraint
	for (int i = 0; i < contactpts; i++)
	{
		for (int j = contactpts; j < 2*contactpts; j++)
		{
			if(i == j-contactpts)
			{
				CI[i][j] = 1;
			} else {
				CI[i][j] = 0;
			}
		}
	}

	ci0.resize(2*contactpts);
	for (int i = 0; i < contactpts; i++)
	{
		ci0[i] = -b[i]*1e6;
	}
	for (int i = contactpts; i < 2*contactpts; i++)
	{
		ci0[i] = 0;
	}

	//x is the contact pressure in bar
	x.resize(contactpts);
	QuadProgPP::solve_quadprog(G, g0, CE, ce0, CI, ci0, x);

	Fluid.p_contact = 0;
	for(int i=0; i<contactpts; i++)
	{
		Fluid.p_contact(contact[i].first, contact[i].second) = x[i]*1e5;
	}

	/*
	//limit pressure
	Fluid.p_contact = where(Fluid.p_contact > 800e5, 800e5, Fluid.p_contact);
	*/
	//limit pressure
	Fluid.p_contact = where(Fluid.p_contact > contact_p_lim, contact_p_lim, Fluid.p_contact);

	//restore the pressure / slipper ehd
	Fluid.p = p;
	slipper.ehd = slip_ehd;

	//relax the contact pressure where there was previous contact
	Fluid.p_contact = where(p_contact_old > 0, p_contact_old + alphaContactP*(Fluid.p_contact-p_contact_old), Fluid.p_contact);

	//Add in the contact pressure
	Fluid.p += AlphaContact * Fluid.p_contact;

	
	//Calculate the contact deformation
	if(FullLoop>0)
	{
		//option 1 - best but expensive
		//SlipperCalcEHD(alphaEHD);

		//option 2 - cheap but a flat contact region
		//Fluid.h = where(Fluid.h < ContactHeight , ContactHeight , Fluid.h);
		
		
		//option 3 - hybrid.
		for (int i = 0; i < contactpts; i++)
		{
			double contact_deform = 0;
			for (int j = 0; j < contactpts; j++)
			{
				 contact_deform += K(i,j)*Fluid.p_contact(contact[j].first, contact[j].second)/100e5;
			}

			/*
			//limit max contact deform
			if(contact_deform > 5e-6)
			{
				contact_deform = 5e-6;
			}
			*/

			//limit max contact deform
			if(contact_deform > contact_def_lim)
			{
				contact_deform = contact_def_lim;
			}

			slipper.ehd(contact[i].first, contact[i].second) += contact_deform;
		}

		SlipperCalch();
		
	}

	//always ensure the film is above contact thickness
	Fluid.h = where(Fluid.h < ContactHeight , ContactHeight , Fluid.h);
	
	return true;
}
