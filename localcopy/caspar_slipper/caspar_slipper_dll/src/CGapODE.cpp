#include "CGapODE.h"
#ifdef _WIN32
#include <ppl.h>
#endif
#pragma once

extern class CGapLog GapLog;

#ifdef _WIN32
//extern variables used in coupled caspar
extern Concurrency::event * coupled_caspar_local_event;
extern Concurrency::event * coupled_caspar_global_event;
extern bool coupled_caspar_simulation;
#endif

CGapODE::CGapODE(caspar_input * gapinputs) : gapinput(gapinputs), NewtonIteration(gapinputs)
{
	GapLog.message("ODE -> Constructing object ODE... ");

	//Number of states for the solver
	n = 3;

};
CGapODE::~CGapODE(void)
{

	GapLog.message("ODE -> Object ODE destructed successfully!");

};

void CGapODE::ODEmain(const bool Resume)
{	
	gapinput->common.resumed = Resume;
	ODESetupSolver(Resume);
	ODEIntegrate();
	ODECleanSolver();
};
void CGapODE::ODESetupSolver(const bool Resume)
{
	
	p.resize(n);
	v.resize(n);
	pold.resize(n);
	vold.resize(n);
	
	//Initialize position / height
	
	p[0] = gapinput->options_slipper.position.hG1;
	p[1] = gapinput->options_slipper.position.hG2;
	p[2] = gapinput->options_slipper.position.hG3;
	v[0] = 0.0;
	v[1] = 0.0;
	v[2] = 0.0;

	//Old is same as initial
	pold = p;
	vold = v;

	//Initial and value for time
	timeStart = 0;	//now hard coded to start at 0

	//Set the present time
	time = timeStart;
	updateTime(time);

	//Resume mode
	if(Resume)
	{
		GapLog.message("CGapODE::ODESetupSolver -> Resuming from resume.bin");

		//Read the resume file
		if(ODEResumeRead())
		{	
			GapLog.message("CGapODE::ODESetupSolver -> Resume time: " + n2s(time));
		}
		updateTime(time);
	}

};
void CGapODE::ODECleanSolver(void)
{
};
void CGapODE::ODEIntegrate(void)
{
	while(gapinput->common.revolution_index < gapinput->lubrication_module.n_lubrication_revolutions) //main program loop
	{

		GapLog << "\nTime = " << time << " s " << "- Angle = " << gapinput->common.phi_deg << " deg" << endl ;

		//Is this the start of a revolution?	

		//Run the newton iteration
		NewtonIteration.NewtonCalcNewtonIteration(p,v,pold,vold);

		//Handle Result Output
		//Resume file update to dirty plots
		ODEResumeWrite(true);

		NewtonIteration.GapResult.plot1d();
		NewtonIteration.GapResult.plot2d(NewtonIteration.SlipperGap);

		//Check for bad gap heights
		//find abs max p
		double mx = abs(p[0]);
		for(unsigned int i=1;i<p.size();i++)
		{
			if(abs(p[i]) > mx)
			{
				mx = abs(p[i]);
			}
		}
		
		if(mx > 1e-3)
		{
			GapLog << "ERROR: Strange gap heights > 1e-3 [m]. Aborting GapModule" << endl;
			return;
		}

		//Advance time
		time += gapinput->common.timestep;
		updateTime(time);

		//Update the old gap heights / velocities
		pold = p;
		vold = v;

		if(gapinput->options_slipper.general.Explicit != 0)
		{
			//Explicit ODE option
			//Integrate the velocity using a first order euler scheme
			for(unsigned int i=0;i<p.size();i++)
			{
				p[i] += v[i]*gapinput->common.timestep;
			}
		}

		//Resume file update to clean plots and new time / p / v
		ODEResumeWrite(false);
		

		#ifdef _WIN32
		//Handle coupled caspar blocking events
		if(coupled_caspar_simulation)
		{
			if(gapinput->common.first_rev_step)
			{
				//this is a new revolution
				coupled_caspar_local_event->set();			//signal that the slipper has finished a revolution
				coupled_caspar_global_event->wait();		//wait for the all interfaces to finish
				coupled_caspar_global_event->reset();		//reset the global event flag
			}

		}
		#endif
		
	} //end of main program loop

}
void CGapODE::ODEResumeWrite(const bool DirtyPlot)
{
	//THE RESUME FILE IS WRITTEN IN BINARY FORMAT
	//THE FORMAT IS AS FOLLOWS:
	/*
		char(2)
		double time
		int n
		vector<double> pold(n)
		vector<double> vold(n)
		char DirtyPlot
		char(3)
	*/

	//NOTE THAT time IS THE TIME THAT HAS NOT YET BEEN SOLVED FOR
	//AND WHERE THE SIMULATION SHOULD RESUME

	fstream f("./output/slipper/resume.bin", ios::out|ios::binary);

	//open byte
	f.put(char(2));
	
	//time
	f.write((char*) &time, sizeof(time));

	//number of variables
	f.write((char*) &n, sizeof(n));

	//position
	for(int i=0;i<n;i++)
	{
		f.write((char*) &(p[i]), sizeof(double));
		f.write((char*) &(pold[i]), sizeof(double));
	}

	//velocity
	for(int i=0;i<n;i++)
	{
		f.write((char*) &(v[i]), sizeof(double));
		f.write((char*) &(vold[i]), sizeof(double));
	}

	//The Slipper Fluid alphaold
	f.write((char*) &(NewtonIteration.SlipperGap->Fluid.alphaPold), sizeof(double));
	f.write((char*) &(NewtonIteration.SlipperGap->operatingslippergap.pG), sizeof(double));
	
	//Slipper deformation
	for(int i=0; i<NewtonIteration.SlipperGap->Fluid.M; i++)
	{
		for(int j=0; j<NewtonIteration.SlipperGap->Fluid.N; j++)
		{
			f.write((char*) &(NewtonIteration.SlipperGap->Fluid.pold(i,j)), sizeof(double));

			if(gapinput->options_slipper.general.SlipperPressureDeformation == 1)
			{
				f.write((char*) &(NewtonIteration.SlipperGap->slipper.ehd(i,j)), sizeof(double));
				f.write((char*) &(NewtonIteration.SlipperGap->Fluid.ehdsqzOld(i,j)), sizeof(double));

				if(gapinput->options_slipper.general.SwashplatePressureDeformation == 1)
				{
					f.write((char*) &(NewtonIteration.SlipperGap->swashplate.ehd(i,j)), sizeof(double));
				}
			}
		}
	}

	//Slipper thermo elastic analysis
	if(gapinput->options_slipper.general.SlipperThermoElastic == 1)
	{
		for(int i=0; i<NewtonIteration.SlipperGap->t_slipper.nodecnt; i++)
		{
			f.write((char*) (&(NewtonIteration.SlipperGap->t_slipper.nodeDeformation[i])), sizeof(double));
			f.write((char*) (&(NewtonIteration.SlipperGap->t_slipper.nodeTemperature[i])), sizeof(double));
		}
		for(int i=0; i<NewtonIteration.SlipperGap->t_slipper.facecnt; i++)
		{
			f.write((char*) (&(NewtonIteration.SlipperGap->t_slipper.faceFlux[i])), sizeof(double));
			f.write((char*) (&(NewtonIteration.SlipperGap->t_slipper.faceFlux_old[i])), sizeof(double));
		}
		for(size_t i=0; i<NewtonIteration.SlipperGap->t_slipper.nodeTemperature_old_FullMesh.size(); i++)
		{
			f.write((char*) (&(NewtonIteration.SlipperGap->t_slipper.nodeTemperature_old_FullMesh[i])), sizeof(double));
		}
		f.write((char*) (&(NewtonIteration.SlipperGap->t_slipper.faceFluxcnt)), sizeof(double));
		f.write((char*) (&(NewtonIteration.SlipperGap->t_slipper.last_run_step)), sizeof(int));
	}

	//Swashplate thermo elastic analysis
	if(gapinput->options_slipper.general.SwashplateThermoElastic)
	{
		for(int i=0; i<NewtonIteration.SlipperGap->t_swashplate.nodecnt; i++)
		{
			f.write((char*) (&(NewtonIteration.SlipperGap->t_swashplate.nodeDeformation[i])), sizeof(double));
			f.write((char*) (&(NewtonIteration.SlipperGap->t_swashplate.nodeTemperature[i])), sizeof(double));
		}
		for(int i=0; i<NewtonIteration.SlipperGap->t_swashplate.facecnt; i++)
		{
			f.write((char*) (&(NewtonIteration.SlipperGap->t_swashplate.faceFlux[i])), sizeof(double));
			f.write((char*) (&(NewtonIteration.SlipperGap->t_swashplate.faceFlux_old[i])), sizeof(double));
		}
		for(size_t i=0; i<NewtonIteration.SlipperGap->t_swashplate.nodeTemperature_old_FullMesh.size(); i++)
		{
			f.write((char*) (&(NewtonIteration.SlipperGap->t_swashplate.nodeTemperature_old_FullMesh[i])), sizeof(double));
		}
		f.write((char*) (&(NewtonIteration.SlipperGap->t_swashplate.faceFluxcnt)), sizeof(double));
		f.write((char*) (&(NewtonIteration.SlipperGap->t_swashplate.last_run_step)), sizeof(int));
	}

	if(gapinput->options_slipper.general.SwashplatePressureDeformation == 1)
	{
		//FullGroup presssures - Needed only for swashplate deformation	
		
		for(int i=0; i<NewtonIteration.SlipperGap->FullGroup.p.extent(0); i++)
		{
			for(int j=0; j<NewtonIteration.SlipperGap->FullGroup.p.extent(1); j++)
			{
				for(int k=0; k<NewtonIteration.SlipperGap->FullGroup.p.extent(2); k++)
				{
					f.write((char*) &(NewtonIteration.SlipperGap->FullGroup.p(i,j,k)), sizeof(double));
				}
			}
		}

		//Swashplate deformation
		for(int i=0; i<NewtonIteration.SlipperGap->swashplate.nodecnt; i++)
		{
			f.write((char*) &(NewtonIteration.SlipperGap->swashplate.nodeDeformation[i]), sizeof(double));
		}
	}

	//dirty plots
	if(DirtyPlot)
	{
		f.put(char(121));	// y
	} else {
		f.put(char(110));	// n
	}
	
	//close byte
	f.put(char(3)); //close byte

	f.close();

	if(gapinput->options_slipper.general.DebugMode == 1 && !DirtyPlot)
	{
		//Only write the file when we finally have a clean plot

		//Copy the resume file to the Resume dir
		#if OS == WINDOWS
			//system("IF NOT EXIST .\\Resume ( mkdir .\\Resume )");
			ifstream ifs("./output/slipper/resume.bin", std::ios::binary);
			ofstream ofs(string("./output/slipper/resume/" + n2s(gapinput->common.phi_deg) + ".bin").c_str(), std::ios::binary);
			ofs << ifs.rdbuf();
			ofs.close();
			ifs.close();
		#endif
	}

}
bool CGapODE::ODEResumeRead(void)
{
	//THE RESUME FILE IS WRITTEN IN BINARY FORMAT
	//THE FORMAT IS AS FOLLOWS:
	/*
		char(2)
		double time
		int n
		vector<double> pold(n)
		vector<double> vold(n)
		char DirtyPlot
		char(3)
	*/

	//NOTE THAT time IS THE TIME THAT HAS NOT YET BEEN SOLVED FOR
	//AND WHERE THE SIMULATION SHOULD RESUME

	fstream f("./output/slipper/resume.bin", ios::in|ios::binary);
	if (!f.is_open()) 	{
		GapLog << "CGapODE::ODEResumeRead -> Unable to open resume.bin! Starting simulation as normal." << endl;
		return false;
	}

	char c;
	f.get(c);
	if (char(2) != c)
	{
		f.close();
		GapLog << "CGapODE::ODEResumeRead -> Error in open byte of resume.bin!" << endl;
		GapLog << CGapLog::error(1) << "Error in resume.bin file! Aborting." << endl;
	}

	f.read(reinterpret_cast<char*> (&time), sizeof(time));
	int nResume;
	f.read(reinterpret_cast<char*> (&nResume), sizeof(nResume));

	if(n != nResume)
	{
		f.close();
		GapLog.message("CGapODE::ODEResumeRead -> n does not match General.gen in resume.bin!");
		GapLog << CGapLog::error(1) << "Error in resume.bin file! Aborting." << endl;
	}
	
	for(int i=0;i<n;i++)
	{
		f.read(reinterpret_cast<char*> (&p[i]), sizeof(double));
		f.read(reinterpret_cast<char*> (&pold[i]), sizeof(double));
	}
	for(int i=0;i<n;i++)
	{
		f.read(reinterpret_cast<char*> (&v[i]), sizeof(double));
		f.read(reinterpret_cast<char*> (&vold[i]), sizeof(double));
	}

	//The Slipper Fluid alphaold
	f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->Fluid.alphaPold)), sizeof(double));
	f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->operatingslippergap.pG)), sizeof(double));

	//Slipper deformation
	for(int i=0; i<NewtonIteration.SlipperGap->Fluid.M; i++)
	{
		for(int j=0; j<NewtonIteration.SlipperGap->Fluid.N; j++)
		{
			f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->Fluid.pold(i,j))), sizeof(double));
				
			if(gapinput->options_slipper.general.SlipperPressureDeformation == 1)
			{
				f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->slipper.ehd(i,j))), sizeof(double));
				f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->Fluid.ehdsqzOld(i,j))), sizeof(double));
				
				if(gapinput->options_slipper.general.SwashplatePressureDeformation == 1)
				{
					f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->swashplate.ehd(i,j))), sizeof(double));
				}
			}
		}
	}

	//Slipper thermo elastic analysis
	if(gapinput->options_slipper.general.SlipperThermoElastic == 1)
	{
		for(int i=0; i<NewtonIteration.SlipperGap->t_slipper.nodecnt; i++)
		{
			f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->t_slipper.nodeDeformation[i])), sizeof(double));
			f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->t_slipper.nodeTemperature[i])), sizeof(double));
		}
		for(int i=0; i<NewtonIteration.SlipperGap->t_slipper.facecnt; i++)
		{
			f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->t_slipper.faceFlux[i])), sizeof(double));
			f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->t_slipper.faceFlux_old[i])), sizeof(double));
		}
		for(size_t i=0; i<NewtonIteration.SlipperGap->t_slipper.nodeTemperature_old_FullMesh.size(); i++)
		{
			f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->t_slipper.nodeTemperature_old_FullMesh[i])), sizeof(double));
		}
		f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->t_slipper.faceFluxcnt)), sizeof(double));
		f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->t_slipper.last_run_step)), sizeof(int));
	}

	//Swashplate thermo elastic analysis
	if(gapinput->options_slipper.general.SwashplateThermoElastic)
	{
		for(int i=0; i<NewtonIteration.SlipperGap->t_swashplate.nodecnt; i++)
		{
			f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->t_swashplate.nodeDeformation[i])), sizeof(double));
			f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->t_swashplate.nodeTemperature[i])), sizeof(double));
		}
		for(int i=0; i<NewtonIteration.SlipperGap->t_swashplate.facecnt; i++)
		{
			f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->t_swashplate.faceFlux[i])), sizeof(double));
			f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->t_swashplate.faceFlux_old[i])), sizeof(double));
		}
		for(size_t i=0; i<NewtonIteration.SlipperGap->t_swashplate.nodeTemperature_old_FullMesh.size(); i++)
		{
			f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->t_swashplate.nodeTemperature_old_FullMesh[i])), sizeof(double));
		}
		f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->t_swashplate.faceFluxcnt)), sizeof(double));
		f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->t_swashplate.last_run_step)), sizeof(int));
	}
	
	if(gapinput->options_slipper.general.SwashplatePressureDeformation == 1)
	{
		//FullGroup presssures - Needed only for swashplate deformation	
		
		for(int i=0; i<NewtonIteration.SlipperGap->FullGroup.p.extent(0); i++)
		{
			for(int j=0; j<NewtonIteration.SlipperGap->FullGroup.p.extent(1); j++)
			{
				for(int k=0; k<NewtonIteration.SlipperGap->FullGroup.p.extent(2); k++)
				{
					f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->FullGroup.p(i,j,k))), sizeof(double));
				}
			}
		}

		//Swashplate deformation
		for(int i=0; i<NewtonIteration.SlipperGap->swashplate.nodecnt; i++)
		{
			f.read(reinterpret_cast<char*> (&(NewtonIteration.SlipperGap->swashplate.nodeDeformation[i])), sizeof(double));
		}
	}
	
	f.get(c);
	if (char(110) != c)
	{
		f.close();
		GapLog.message("CGapODE::ODEResumeRead -> A DirtyPlot status was reported in resume.bin!");
		GapLog << "The resume.bin file indicates that the simulation was interrupted in the middle " <<
			"of writing output files." << endl << "Resuming now requires advanced techniques." << endl;
		GapLog << CGapLog::error(1) << "Aborting." << endl;
	}

	f.get(c);
	if (char(3) != c)
	{
		f.close();
		GapLog.message("CGapODE::ODEResumeRead -> Error in close byte of resume.bin!");
		GapLog << CGapLog::error(1) << "Error in resume.bin file! Aborting." << endl;
	}

	f.close();

	return true;
}
void CGapODE::updateTime(const double time_s)
{
	//reset the step flags
	gapinput->common.first_rev_step = false;
	gapinput->common.last_rev_step = false;
	gapinput->common.first_step = false;
	gapinput->common.last_step = false;
	
	//update the continous time / phi
	gapinput->common.time = time_s;
	gapinput->common.phi_rad = gapinput->common.time / gapinput->common.rev_period * 2.0 * PI;
	gapinput->common.phi_deg = gapinput->common.time / gapinput->common.rev_period * 360.0;
		
	//update the rev phi
	gapinput->common.phi_rev_deg = fmod(gapinput->common.phi_deg, 360.0);
	gapinput->common.phi_rev_rad = fmod(gapinput->common.phi_rad, 2.0*PI);

	//first revolution step test
	if(gapinput->common.phi_rev_deg <= gapinput->common.phi_deg_tol || 
	   gapinput->common.phi_rev_deg >= (360.0 - gapinput->common.phi_deg_tol))
	{
		gapinput->common.phi_rev_deg = 0;
		gapinput->common.phi_rev_rad = 0;
		gapinput->common.first_rev_step = true;
		gapinput->common.revolution_index++;
	}

	//last revolution step test
	if(gapinput->common.phi_rev_deg + gapinput->common.phistep_deg >= (360.0 - gapinput->common.phi_deg_tol))
	{
		gapinput->common.last_rev_step = true;
	}

	//first step
	if(gapinput->common.time == timeStart)
	{
		gapinput->common.first_step = true;
		gapinput->common.revolution_index = 0;
	}

	//last step
	if(gapinput->common.time + gapinput->common.timestep > 
		gapinput->common.rev_period * (double) gapinput->lubrication_module.n_lubrication_revolutions
		)
	{
		gapinput->common.last_step = true;
	}
	


}
