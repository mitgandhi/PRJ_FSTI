#include "CSlipperGap.h"
#include "CGapResult.h"
#include <stdio.h>
#pragma once

void CGapResult::make_summed(vector<double> &sums)
{
	//init the sums vector
	sums.clear();
	sums.resize(3, 0);
	
	//a vector that will contain the data
	vector<vector<double> > txt_data;
	
	//load the .txt file
	ifstream txt("./output/slipper/slipper.txt");

	if(!txt.is_open())
	{
		//unable to open
		return;
	}

	while(!txt.eof())
	{
		string l;
		getline(txt, l);
		if(l.find("%") != string::npos)
		{
			//technically the find isn't correct. % should be the first non-white char

			//this is a commented line, so skip 
			continue;
		}

		vector<string> fields;
		CSlipperGap::Tokenize(l, fields, "\t ");

		vector<double> dfields;
		for(size_t i=0; i<fields.size(); i++)
		{
			dfields.push_back(atof(fields[i].c_str()));
		}

		if(dfields.size() > 0)
		{
			txt_data.push_back(dfields);
		}

	}

	txt.close();
	
	if(txt_data.size() == 0)
	{
		//don't bother if empty
		return;
	}

	//done loading the txt, now we need to add in the value from this revolution
	//only bother with the values we'll actually use
	{
		vector<double> tmp(txt_data[0].size());
		tmp[12] = SlipperResult.QSG;
		tmp[13] = SlipperResult.TorqueLoss;
		tmp[14] = SlipperResult.Ploss;
		txt_data.push_back(tmp);
	}

	//offset vector

	vector<int> offset(gapinput->operating_conditions.npistons);
	
	double spacing = double(gapinput->common.rev_steps) / double(gapinput->operating_conditions.npistons);
	for(size_t i=0; i<gapinput->operating_conditions.npistons; i++)
	{
		offset[i] = (int) floor(spacing * double(i));
	}

	//and now compute the three sums
	{
		size_t row = txt_data.size()-1;
		vector<size_t> cols(3);
		cols[0] = 12;
		cols[1] = 13;
		cols[2] = 14;

		if(row >= offset[gapinput->operating_conditions.npistons-1])
		{
			for(size_t i=0; i<gapinput->operating_conditions.npistons; i++)
			{
				for(size_t j=0; j<3; j++)
				{
					sums[j] += txt_data[row-offset[i]][cols[j]];
				}
			}
		} 
	}
	
}

void CGapResult::plot1d(void)
{


	if(gapinput->lubrication_module.solve_slipper)
	{
	
		//we first need to compute 'summed' values from this timestep
		vector<double> summed_vals;
		make_summed(summed_vals);
		

		//The revised slipper.caspar output file
		{
			//should we create a header?
			ifstream f("./output/slipper/slipper.txt");
			if(f.is_open())
			{
				//file exists, don't need to create a header
				f.close();
			} else {
				//file doesn't exist, so try to create it and write a header
				ofstream of("./output/slipper/slipper.txt");
				if(of.is_open())
				{
					of <<
					"%t" << "\t" <<
					"rev" << "\t" <<
					"phi" << "\t" <<
					"h1" << "\t" <<
					"h2" << "\t" <<
					"h3" << "\t" <<
					"dhdt1" << "\t" <<
					"dhdt2" << "\t" <<
					"dhdt3" << "\t" <<
					"minheight" << "\t" <<
					"meanheight" << "\t" <<
					"maxheight" << "\t" <<
					"QSG" << "\t" <<
					"TorqueLoss" << "\t" <<
					"Powerloss" << "\t" <<
					"QSG_total" << "\t" <<
					"TorqueLoss_total" << "\t" <<
					"Powerloss_total" << "\t" <<
					"FfGz" << "\t" <<
					"MfGx" << "\t" <<
					"MfGy" << "\t" <<
					"M_FTGx" << "\t" <<
					"M_FTGy" << "\t" <<
					"F_FTGx" << "\t" <<
					"F_FTGy" << "\t" <<
					"FTG" << "\t" <<
					"FTK" << "\t" <<
					"FSK" << "\t" <<
					"MGx_centrifugal" << "\t" <<
					"F_fluid1" << "\t" <<
					"F_fluid2" << "\t" <<
					"F_fluid3" << "\t" <<
					"F_external1" << "\t" <<
					"F_external2" << "\t" <<
					"F_external3" << "\t" <<
					"F_holder1" << "\t" <<
					"F_holder2" << "\t" <<
					"F_holder3" << "\t" <<
					"dF1" << "\t" <<
					"dF2" << "\t" <<
					"dF3" << "\t" <<
					"dFz" << "\t" <<
					"dFMx" << "\t" <<
					"dFMy" << "\t" <<
					"pHP" << "\t" <<
					"pLP" << "\t" <<
					"pDC" << "\t" <<
					"pG" << "\t" <<
					"real_h1" << "\t" <<
					"real_h2" << "\t" <<
					"real_h3" << "\t" <<
					"Q_S_pois" << "\t" <<
					"Q_S_couette" << "\t" <<
					"M_TJx" << "\t" <<
					"M_TJy" << "\t" <<
					"avg_p_contact" << "\t" <<
					"visc_sock" << "\t" <<
					"tilt_speed" << "\t" <<
					"psocket" << "\t" <<
					"mu_coef" << "\t" <<
					"F_fric" << "\t" <<
					endl;

					of.close();
				}
			}


		}

		//prepare to output the _actual_ fluid film thickness in the nearest "cell" to the three control points
		//this should probably be done in the CSlipperGap class, as the code has been structured but oh well
		double actualFluidH1, actualFluidH2, actualFluidH3;
		{
			const int M = SlipperResult.h_slipper.extent(0)-1;
			const int N = SlipperResult.h_slipper.extent(1)-1;
			const int dN = N/3;

			actualFluidH1 = SlipperResult.h_slipper(M,0);
			actualFluidH2 = SlipperResult.h_slipper(M,dN);
			actualFluidH3 = SlipperResult.h_slipper(M,2*dN);
		}

		ofstream foutSlipper("./output/slipper/slipper.txt",ios::app);
		foutSlipper << scientific << 
			time << "\t" << scientific <<										// 1
			gapinput->common.phi_deg/360.0 << "\t" << scientific <<				// 2
			gapinput->common.phi_rev_deg << "\t" << scientific <<				// 3
			gappositions[0] << "\t" << scientific <<							// 4
			gappositions[1] << "\t" << scientific <<							// 5
			gappositions[2]	<< "\t" << scientific <<							// 6
			gapvelocities[0] << "\t" << scientific <<							// 7
			gapvelocities[1] << "\t" << scientific <<							// 8
			gapvelocities[2] << "\t" << scientific <<							// 9
			SlipperResult.minh << "\t" << scientific <<							// 10
			SlipperResult.meanh << "\t" << scientific <<						// 11
			SlipperResult.maxh << "\t" << scientific <<							// 12
			SlipperResult.QSG << "\t" << scientific <<							// 13
			SlipperResult.TorqueLoss << "\t" << scientific <<					// 14
			SlipperResult.Ploss << "\t" << scientific << 						// 15
			summed_vals[0] << "\t" << scientific << 							// 16
			summed_vals[1] << "\t" << scientific << 							// 17
			summed_vals[2] << "\t" << scientific << 							// 18
			SlipperResult.FfGz << "\t" << scientific <<							// 19
			SlipperResult.MfGx << "\t" << scientific <<							// 20
			SlipperResult.MfGy << "\t" << scientific <<							// 21
			SlipperResult.M_FTGx << "\t" << scientific <<						// 22
			SlipperResult.M_FTGy << "\t" << scientific <<						// 23
			SlipperResult.F_FTGx << "\t" << scientific <<						// 24
			SlipperResult.F_FTGy << "\t" << scientific <<						// 25
			SlipperResult.FTG << "\t" << scientific <<							// 26
			SlipperResult.FTK << "\t" << scientific <<							// 27
			SlipperResult.FSK << "\t" << scientific <<							// 28
			SlipperResult.MGx_centrifugal << "\t" << scientific <<				// 29
			SlipperResult.F_fluid[0] << "\t" << scientific <<					// 30
			SlipperResult.F_fluid[1] << "\t" << scientific <<					// 31
			SlipperResult.F_fluid[2] << "\t" << scientific <<					// 32
			SlipperResult.F_external[0] << "\t" << scientific <<				// 33
			SlipperResult.F_external[1] << "\t" << scientific <<				// 34
			SlipperResult.F_external[2] << "\t" << scientific <<				// 35
			SlipperResult.F_contact[0] << "\t" << scientific <<					// 36
			SlipperResult.F_contact[1] << "\t" << scientific <<					// 37
			SlipperResult.F_contact[2] << "\t" << scientific <<					// 38
			SlipperResult.dF[0] << "\t" << scientific <<						// 39
			SlipperResult.dF[1] << "\t" << scientific <<						// 40
			SlipperResult.dF[2] << "\t" << scientific <<						// 41
			SlipperResult.dFcomp[0] << "\t" << scientific <<					// 42
			SlipperResult.dFcomp[1] << "\t" << scientific <<					// 43
			SlipperResult.dFcomp[2] << "\t" << scientific <<					// 44
			pHP << "\t" << scientific <<										// 45
			pLP << "\t" << scientific <<										// 46
			pDC << "\t" << scientific <<										// 47
			SlipperResult.pG << "\t" << scientific << 							// 48
			actualFluidH1 << "\t" << scientific << 								// 49
			actualFluidH2 << "\t" << scientific << 								// 50
			actualFluidH3 << "\t" << scientific << 								// 51
			SlipperResult.Q_S_pois << "\t" << scientific <<						// 52
			SlipperResult.Q_S_couette << "\t" << scientific <<					// 53
			SlipperResult.M_TJx << "\t" << scientific <<						// 54
			SlipperResult.M_TJy << "\t" << scientific <<						// 55
			SlipperResult.avg_p_contact <<  "\t" << scientific <<				// 56
			SlipperResult.visc_sock << "\t" << scientific <<					// 57
			SlipperResult.tilt_speed << "\t" << scientific <<					// 58
			SlipperResult.psocket << "\t" << scientific <<						// 59
			SlipperResult.mu_coef << "\t" << scientific <<						// 60
			SlipperResult.F_fric << "\t" << scientific <<						// 61
		endl;
		foutSlipper.close();

		ofstream foutSlipperForce("./output/slipper/ftg.txt",ios::app);
		foutSlipperForce << scientific << 
			time << "\t" << scientific <<										// 1
			phi << "\t" << scientific <<										// 2
			SlipperResult.FTG << "\t" << scientific <<							// 3
			endl;
		foutSlipperForce.close();


	}
}
void CGapResult::plot2d(CSlipperGap * SlipperGap)
{

	if(gapinput->lubrication_module.solve_slipper)
	{
		/*
		//slipper ehd squeeze pressure
		ofstream foutpGehd("./Outputs/Slipper_Pehd.pG",ios::app);
		foutpGehd << "%" << scientific << phi << endl;
		for(int i=0;i<SlipperResult.P_slipper.extent(0);i++)
		{
			for(int j=0;j<SlipperResult.P_slipper.extent(1);j++)
			{
				foutpGehd << SlipperResult.P_ehdsq(i,j)  << " ";
			}
			foutpGehd << endl;
		}
		foutpGehd.close();
		*/

		//slipper velocity (normally only needed for debugging)
		/*
		ofstream foutvr("./Outputs/Slipper_velocity.vr",ios::app);
		foutvr << "%" << scientific << phi << endl;
		for(int i=0;i<SlipperResult.vr_slipper.extent(0);i++)
		{
			for(int j=0;j<SlipperResult.vr_slipper.extent(1);j++)
			{
				foutvr << SlipperResult.vr_slipper(i,j,SlipperResult.vtheta_slipper.extent(2)-1)  << " ";
			}
			foutvr << endl;
		}
		foutvr.close();
		ofstream foutvt("./Outputs/Slipper_velocity.vt",ios::app);
		foutvt << "%" << scientific << phi << endl;
		for(int i=0;i<SlipperResult.vtheta_slipper.extent(0);i++)
		{
			for(int j=0;j<SlipperResult.vtheta_slipper.extent(1);j++)
			{
				foutvt << SlipperResult.vtheta_slipper(i,j,SlipperResult.vtheta_slipper.extent(2)-1)  << " ";
			}
			foutvt << endl;
		}
		foutvt.close();
		*/

		/*
		//Normally only needed for debugging
		//slipper velocity in rect coords
		ofstream foutvgx("./Outputs/Slipper_vgx.vg",ios::app);
		foutvgx << "%" << scientific << phi << endl;
		for(int i=0;i<SlipperResult.vgx.extent(0);i++)
		{
			for(int j=0;j<SlipperResult.vgx.extent(1);j++)
			{
				foutvgx << SlipperResult.vgx(i,j)  << " ";
			}
			foutvgx << endl;
		}
		foutvgx.close();

		ofstream foutvgy("./Outputs/Slipper_vgy.vg",ios::app);
		foutvgy << "%" << scientific << phi << endl;
		for(int i=0;i<SlipperResult.vgy.extent(0);i++)
		{
			for(int j=0;j<SlipperResult.vgy.extent(1);j++)
			{
				foutvgy << SlipperResult.vgy(i,j)  << " ";
			}
			foutvgy << endl;
		}
		foutvgy.close();
		*/
		/*
		//slipper pressure
		ofstream foutpG((txt_dir+"Slipper_Pressure.pG").c_str(),ios::app);
		foutpG << "%" << scientific << phi << endl;
		for(int i=0;i<SlipperResult.P_slipper.extent(0);i++)
		{
			for(int j=0;j<SlipperResult.P_slipper.extent(1);j++)
			{
				foutpG << SlipperResult.P_slipper(i,j)  << " ";
			}
			foutpG << '\n';
		}
		foutpG.close();

		if (gapinput->options_slipper.general.SlipperPressureDeformation)
		{
			//slipper ehd
			ofstream foutEHD((txt_dir+"Slipper_EHD.hG").c_str(),ios::app);
			foutEHD << "%" << scientific << phi << endl;
			for(int i=0;i<SlipperResult.slipperEHD.extent(0);i++)
			{
				for(int j=0;j<SlipperResult.slipperEHD.extent(1);j++)
				{
					foutEHD << SlipperResult.slipperEHD(i,j)  << " ";
				}
				foutEHD << '\n';
			}
			foutEHD.close();
		}

		if (gapinput->options_slipper.general.SwashplatePressureDeformation)
		{
			//slipper ehd
			ofstream foutEHD((txt_dir+"Swashplate_EHD.hG").c_str(),ios::app);
			foutEHD << "%" << scientific << phi << endl;
			for(int i=0;i<SlipperResult.swashplateEHD.extent(0);i++)
			{
				for(int j=0;j<SlipperResult.swashplateEHD.extent(1);j++)
				{
					foutEHD << SlipperResult.swashplateEHD(i,j)  << " ";
				}
				foutEHD << '\n';
			}
			foutEHD.close();
		}
		*/

		//slipper gap height
		ofstream fouthG((txt_dir+"Slipper_height.hG").c_str(),ios::app);
		fouthG << "%" << scientific << phi << endl;
		for(int i=0;i<SlipperResult.h_slipper.extent(0);i++)
		{
			for(int j=0;j<SlipperResult.h_slipper.extent(1);j++)
			{
				fouthG << SlipperResult.h_slipper(i,j)  << " ";
			}
			fouthG << '\n';
		}
		fouthG.close();


		SlipperGap->writevtk(vtk_dir+"fluid." + n2s(phi) + ".vtk");
		if(phi - 721 >= 0)
		{
			//remove vtk's older than two revolutions
			remove((vtk_dir+"fluid." + n2s(int(floor(phi-721.0+0.5))) + ".vtk").c_str());
		}
		if(gapinput->options_slipper.general.SwashplatePressureDeformation)
		{
			if(gapinput->options_slipper.general.SwashplateThermoElastic)
			{
				vector<double> qflux_tmp(SlipperGap->t_swashplate.faceFlux);
				for(int i=0; i<qflux_tmp.size(); i++)
				{
					//better to use revsteps for the avg which should be equal to SlipperGap->t_swashplate.faceFluxcnt
					//at the end of the revolution when the proper qflux will be calculated
					qflux_tmp[i] = double(gapinput->operating_conditions.npistons) * qflux_tmp[i] / double(gapinput->common.rev_steps);
				}
				SlipperGap->swashplate.writevtk(vtk_dir+"swashplate." + n2s(phi) + ".vtk", qflux_tmp);
			} else {
				SlipperGap->swashplate.writevtk(vtk_dir+"swashplate." + n2s(phi) + ".vtk");
			}
			if(phi - 721 >= 0)
			{
				//remove vtk's older than two revolutions
				remove((vtk_dir+"swashplate." + n2s(int(floor(phi-721.0+0.5))) + ".vtk").c_str());
			}
		}

	}

	
}
CGapResult::CGapResult(caspar_input * gapinputs) : gapinput(gapinputs)
{
	//define the dir's
	txt_dir = "./output/slipper/matlab/";
	vtk_dir = "./output/slipper/vtk/";
};
