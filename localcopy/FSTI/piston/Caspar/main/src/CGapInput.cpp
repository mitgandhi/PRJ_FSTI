#include "CGapInput.h"
#include "logger.h"
#include "../../caspar_input/input.h"
#include <iomanip>
#include <iostream>
#pragma once

extern class input myinput;

CGapInput::CGapInput(void)
{

}
CGapInput::~CGapInput(void)	//Destructor
{

}; 
/*void CGapInput::readGeneral(void) {

	// open the input file
	ifstream fgeneral("./inputs_piston/General.gen");
	
	if (!fgeneral.is_open()) 	{
		Log << "Unable to open General.gen! - Generate General.gen and restart..." << "\n";
		system("PAUSE");
		exit(1);
	}

	string tmp;
	string line;

	while(fgeneral) {
		getline(fgeneral,line);
		istringstream iss (line,istringstream::in);
		while(iss) {
			iss >> tmp;
			
			// check if there is a comment
			size_t comment  = tmp.find("//");
			if (comment!=string::npos)
			{
				break;
			}
			if (tmp == "mode"){
				iss >> General.mode;
				break;
			}
			if (tmp == "npistons"){
				iss >> General.npistons;
				break;
			}
			if (tmp == "nrevolutions"){
				iss >> General.nrevolutions;
				break;
			}
			if (tmp == "speed"){
				iss >> General.speed;
				break;
			}
			if (tmp == "dB"){
				iss >> General.dB;
				break;
			}
			if (tmp == "dK"){
				iss >> General.dK;
				break;
			}
			if (tmp == "beta"){
				iss >> General.beta;
				break;
			}
			if (tmp == "gamma"){
				iss >> General.gamma;
				break;
			}
			if (tmp == "betamax"){
				iss >> General.betamax;
				break;
			}
			if (tmp == "pHP"){
				iss >> General.pHP;
				break;
			}
			if (tmp == "pLP"){
				iss >> General.pLP;
				break;
			}
			if (tmp == "pCase"){
				iss >> General.pCase;
				break;
			}
			if (tmp == "THP"){
				iss >> General.THP;
				break;
			}
			if (tmp == "TLP"){
				iss >> General.TLP;
				break;
			}
			if (tmp == "TCase"){
				iss >> General.TCase;
				break;
			}
			if (tmp == "oiltype"){
				iss >> General.oiltype;
				break;
			}
			if (tmp == "oildensity"){
				iss >> General.oildensity;
				break;
			}
			if (tmp == "oilbetaP"){
				iss >> General.oilbetaP;
				break;
			}
			if (tmp == "oilbetaT"){
				iss >> General.oilbetaT;
				break;
			}
			if (tmp == "oilK"){
				iss >> General.oilK;
				break;
			}
			if (tmp == "oilbetaKP"){
				iss >> General.oilbetaKP;
				break;
			}
			if (tmp == "oilbetaKT"){
				iss >> General.oilbetaKT;
				break;
			}
			if (tmp == "oilviscosity"){
				iss >> General.oilviscosity;
				break;
			}
			if (tmp == "oilW"){
				iss >> General.oilW;
				break;
			}
			if (tmp == "oilTc1"){
				iss >> General.oilTc1;
				break;
			}
			if (tmp == "oilPc1"){
				iss >> General.oilPc1;
				break;
			}
			if (tmp == "oilTc2"){
				iss >> General.oilTc2;
				break;
			}
			if (tmp == "oilPc2"){
				iss >> General.oilPc2;
				break;
			}
			if (tmp == "oillambda"){
				iss >> General.oillambda;
				break;
			}
			if (tmp == "oilC"){
				iss >> General.oilC;
				break;
			}
			if (tmp == "alpha1"){
				iss >> General.alpha1;
				break;
			}
			if (tmp == "alpha2"){
				iss >> General.alpha2;
				break;
			}
			if (tmp == "alpha3"){
				iss >> General.alpha3;
				break;
			}
		}
	}

	fgeneral.close();

}
//removed dwm
*/

/*void CGapInput::readGeometry(void) {

	ifstream fgeometry("./inputs_piston/Geometry.geo");
	if (!fgeometry.is_open()) 	{
		Log << "Unable to open Geometry.geo! - Generate Geometry.geo and restart..." << "\n";
		system("PAUSE");
		exit(1);
	}

	string tmp;
	string line;

	while(fgeometry) {
		getline(fgeometry,line);
		istringstream iss (line,istringstream::in);
		while(iss) {
			iss >> tmp;
			// check if there is a comment
			size_t comment  = tmp.find("//");
			if (comment!=string::npos)
			{
				break;
			}
			if (tmp == "dZ"){
				iss >> Geometry.dZ;
				break;
			}
			if (tmp == "lK"){
				iss >> Geometry.lK;
				break;
			}
			if (tmp == "lF"){
				iss >> Geometry.lF;
				break;
			}
			if (tmp == "le"){
				iss >> Geometry.le;
				break;
			}
			if (tmp == "lZ0"){
				iss >> Geometry.lZ0;
				break;
			}
			if (tmp == "dDK"){
				iss >> Geometry.dDK;
				break;
			}
			if (tmp == "doutG"){
				iss >> Geometry.doutG;
				break;
			}
			if (tmp == "dinG"){
				iss >> Geometry.dinG;
				break;
			}
			if (tmp == "lSK"){
				iss >> Geometry.lSK;
				break;
			}
			if (tmp == "mK"){
				iss >> Geometry.mK;
				break;
			}
			if (tmp == "lKG"){
				iss >> Geometry.lKG;
				break;
			}
			if (tmp == "lch"){
				iss >> Geometry.lch;
				break;
			}
			if (tmp == "rK_red"){
				iss >> Geometry.rK_red;
				break;
			}
			if (tmp == "lK_hs"){
				iss >> Geometry.lK_hs;
				break;
			}
			if (tmp == "hmin"){
				iss >> Geometry.hmin;
				break;
			}
			if (tmp == "lengthB"){
				iss >> Geometry.lengthB;
				break;
			}
			if (tmp == "lengthcanalB"){
				iss >> Geometry.lengthcanalB;
				break;
			}
			if (tmp == "speedK"){
				iss >> Geometry.speedK;
				break;
			}
			if (tmp == "PistonMacroGeometry"){
				iss >> Geometry.PistonMacroGeometry;
				break;
			}
			if (tmp == "CylinderMacroGeometry"){
				iss >> Geometry.CylinderMacroGeometry;
				break;
			}
		}
	}

	fgeometry.close();


}*/

/*void CGapInput::readBoundary(void) {

	ifstream fboundary("./inputs_piston/boundary.bpd");
	if (!fboundary.is_open()) 	{
		Log << "Unable to open Boundary.bpd! - Generate Boundary.bpd and restart..." << "\n";
		system("PAUSE");
		exit(1);
	}

	string tmp;
	string line;

	while(fboundary) {
		getline(fboundary,line);
		istringstream iss (line,istringstream::in);
		while(iss) {
			iss >> tmp;
			
			//Check if there is a comment
			size_t comment  = tmp.find("//");
			if (comment!=string::npos)
			{
				break;
			}

			//Temperature boundaries
			if (tmp == "Tmax"){
				iss >> Boundary.Tmax;
				break;
			}
			if (tmp == "AlphaDC"){
			iss >> Boundary.AlphaDC;
			break;
			}
			if (tmp == "AlphaCase"){
				iss >> Boundary.AlphaCase;
				break;
			}

			//Materials
			if (tmp == "EmodK"){
				iss >> Boundary.EmodK;
				break;
			}
			if (tmp == "EmodB"){
				iss >> Boundary.EmodB;
				break;
			}
			if (tmp == "vK"){
				iss >> Boundary.vK;
				break;
			}
			if (tmp == "vB"){
				iss >> Boundary.vB;
				break;
			}

			//Initial positions
			if (tmp == "PistonInPos"){
				iss >> Boundary.xA >> Boundary.yA >> Boundary.xB >> Boundary.yB;
				break;
			}

			//ODE solver paramters
			if (tmp == "Simalphastep"){
				iss >> Boundary.Simalphastep;
				break;
			}
			if (tmp == "Simalphaplot"){
				iss >> Boundary.Simalphaplot;
				break;
			}
		
			//Convergence loops parameters
			if (tmp == "AlphaP"){
				iss >> Boundary.AlphaP;
				break;
			}
			if (tmp == "AlphaDef"){
				iss >> Boundary.AlphaDef;
				break;
			}
			if (tmp == "AlphaMu"){
				iss >> Boundary.AlphaMu;
				break;
			}
			if (tmp == "AlphaTh"){
				iss >> Boundary.AlphaTh;
				break;
			}
			if (tmp == "Rmin_h"){
				iss >> Boundary.Rmin_h;
				break;
			}
			if (tmp == "Rmin_p"){
				iss >> Boundary.Rmin_p;
				break;
			}
			if (tmp == "nmax"){
				iss >> Boundary.nmax;
				break;
			}
			if (tmp == "penCells"){
				iss >> Boundary.penCells;
				break;
			}
			if (tmp == "Rmin_R"){
				iss >> Boundary.Rmin_R;
				break;
			}
			if (tmp == "Rmin_E"){
				iss >> Boundary.Rmin_E;
				break;
			}

			//Newton loop parameters
			if (tmp == "epsilonK"){
				iss >> Boundary.epsilonK;
				break;
			}
			if (tmp == "jmax"){
				iss >> Boundary.jmax;
				break;
			}
			if (tmp == "kmax"){
				iss >> Boundary.kmax;
				break;
			}
			if (tmp == "delta_v"){
				iss >> Boundary.delta_v;
				break;
			}
			if (tmp == "max_v"){
				iss >> Boundary.p_max;
				break;
			}
		}
	
	}

	fboundary.close();


	//Time step exceptions
	if(Boundary.Simalphastep<0.1)
	{
		Log << "\n";
		Log << "\n";
		Log << "Simulation time step is lower than 0.1 deg! Increase time step and restart..." << "\t";
		Log << "\n";
		system("PAUSE");
	};
	if(Boundary.Simalphastep>0.5)
	{
		Log << "\n";
		Log << "\n";
		Log << "Simulation time step is higher than 0.5 deg! Decrease time step and restart..." << "\t";
		Log << "\n";
		system("PAUSE");
	};


}*/


/*void CGapInput::readGrids(void) {

	ifstream fgrid("./inputs_piston/Grids.grs");
	if (!fgrid.is_open()) 	{
		Log << "Unable to open Grids.grs! - Generate Grids.grs and restart..." << "\n";
		system("PAUSE");
		exit(1);
	}

	string tmp;
	string line;
	vector<int> nMeshMG;
	nMeshMG.resize(2);

	while(fgrid)
	{
		getline(fgrid,line);
		istringstream outiss(line,istringstream::in);
		while(outiss)
		{
			outiss >> tmp;
			
			// check if there is a comment
			size_t comment  = tmp.find("//");
			if (comment!=string::npos)
			{
				break;
			}
			//read gauss-seidel properties
			if(tmp == "nMeshGS")
			{
				getline(fgrid,line);
				istringstream iniss(line);
				while(iniss.str() != "}")
				{
					while(iniss)
					{
						iniss >> tmp;
						if(tmp == "N")
						{
							iniss >> PistonGapGrid.N;
							break;
						}
						if(tmp == "M")
						{
							iniss >> PistonGapGrid.M;
							break;
						}
						if(tmp == "Q")
						{
							iniss >> PistonGapGrid.Q;
							break;
						}
					}
					iniss.clear();
					getline(fgrid,line);
					iniss.str(line);
				}
			}
			//read mutigrid properties
			if(PistonOptions.ReynoldsMultiGrid)
			{
				if(tmp == "nMeshMG")
				{
					getline(fgrid,line);
					istringstream iniss(line);
					while(iniss.str() != "}")
					{
						while(iniss)
						{
							iniss >> nMeshMG[0] >> nMeshMG[1];
							PistonGapMultiGrid.N.push_back(nMeshMG[0]);
							PistonGapMultiGrid.M.push_back(nMeshMG[1]);
							break;
						}
						iniss.clear();
						getline(fgrid,line);
						iniss.str(line);
					}
				}
				if (tmp == "Q")
				{
					outiss >> PistonGapMultiGrid.Q;
					break;
				}
				if (tmp == "VW")
				{
					outiss >> PistonGapMultiGrid.VW;
					break;
				}
				if (tmp == "MGInt")
				{
					outiss >> PistonGapMultiGrid.MGInt;
					break;
				}
				if (tmp == "v1")
				{
					outiss >> PistonGapMultiGrid.v1;
					break;
				}
				if (tmp == "v2")
				{
					outiss >> PistonGapMultiGrid.v2;
					break;
				}
			};
		};
	};
	//assign multigrid levels
	if(PistonOptions.ReynoldsMultiGrid)
	{
		PistonGapMultiGrid.nL = (int) PistonGapMultiGrid.N.size();
	}

	fgrid.close();

}*/

void CGapInput::readpFile(void)
{

	char str[256];
	char cbuffer;
	double ttemp;
	double pDCtemp;
	double pHPtemp;
	double pLPtemp;

	//Polynomial piston
	ifstream fpFile(myinput.data.operating_conditions.pModuleFile);

	if (!fpFile.is_open()) 	{
	Log << "Unable to open pFile.dat! - Generate pFile.dat and restart..." << endl;
	//system("PAUSE");
	exit(1);
	}
	
	while(!fpFile.eof()) 
	{
		if( fpFile.get(cbuffer) && (cbuffer =='/'))
		{
			fpFile.getline(str,sizeof(str));continue;
		}
		else
		{
			fpFile.putback(cbuffer);
		}
			
		fpFile >> ttemp >> pDCtemp >> pHPtemp >> pLPtemp;
		pFile.time.push_back(ttemp);
		pFile.pDC.push_back(pDCtemp);
		pFile.pHP.push_back(pHPtemp);
		pFile.pLP.push_back(pLPtemp);
	}
	fpFile.close();

}
/*void CGapInput::readOptions(void) {

	ifstream foptions("./inputs_piston/Options.opt");
	if (!foptions.is_open()) 	{
		Log << "Unable to open Options.opt! - Generate Options.opt and restart..." << "\n";
		system("PAUSE");
		exit(1);
	}

	string tmp;
	string line;
	size_t comment;
	
	while(foptions) {
		getline(foptions,line);
		istringstream outiss (line,istringstream::in);
		outiss >> tmp;
		// piston gap mesh
		if(tmp == "PistonOptions") {
			getline(foptions,line);
			istringstream iniss(line);
			while(iniss.str() != "}") {
				while(iniss) {
					iniss >> tmp;
					comment  = tmp.find("//"); // check if there is a comment
					if (comment!=string::npos) {
						break;
					}

					if(tmp == "ReadpFile"){
						iniss >> PistonOptions.ReadpFile;
						break;
					}
					if(tmp == "ReynoldsMultiGrid"){
						iniss >> PistonOptions.ReynoldsMultiGrid;
						break;
					}
					if(tmp == "EnergyEquation"){
						iniss >> PistonOptions.EnergyEquation;
						break;
					}
					if(tmp == "HeatTransfer"){
						iniss >> PistonOptions.HeatTransfer;
						break;
					}
					if(tmp == "PressureDeformation"){
						iniss >> PistonOptions.PressureDeformation;
						break;
					}
					if(tmp == "PressureDeformationOMP"){
						iniss >> PistonOptions.PressureDeformationOMP;
						break;
					}
					if(tmp == "ThermalDeformation"){
						iniss >> PistonOptions.ThermalDeformation;
						break;
					}
					if(tmp == "EHDTestRig"){
						iniss >> PistonOptions.EHDTestRig;
						break;
					}
					if(tmp == "TriboTestRig"){
						iniss >> PistonOptions.TriboTestRig;
						break;
					}
				}
				iniss.clear();	// this is essential!!
				getline(foptions,line);
				iniss.str(line);
			}
		}
	}
	
	foptions.close();

}*/

void CGapInput::readMacroGeometryPiston(void)
{
	double polytemp;
	double dtemp;
	double ltemp;
	string line;
	size_t comment;
	size_t eof;
	istringstream iniss;
	

	//Polynomial piston
	if(myinput.data.options_piston.general.McrK == 5 || myinput.data.options_piston.general.McrK == 1)
	{	
		ifstream fmacrogeometry(myinput.data.options_piston.general.McrK_file);
		if (!fmacrogeometry.is_open()) 	{
			Log << "Unable to open " << myinput.data.options_piston.general.McrK_file << "! - Generate and restart..." << endl;
			//system("PAUSE");
			exit(1);
		}
		//Polynomial piston
		if( myinput.data.options_piston.general.McrK == 5)
		{	
			while(fmacrogeometry)
			{
				iniss.clear();
				getline(fmacrogeometry,line);
				iniss.str(line);
				eof = line.find("*");
				comment = line.find("//");
				if(eof!=string::npos)
				{ 
					break; 
				}
				else if(comment!=string::npos)
				{
					continue;
				}
				else
				{
					iniss >> polytemp;
					Geometry.polygap_coeff.push_back(polytemp);
				};
			};
			fmacrogeometry.clear();
			fmacrogeometry.close();
		}
		//Stepwise piston
		else if( myinput.data.options_piston.general.McrK == 1)
		{

			while(fmacrogeometry)
			{
				iniss.clear();
				getline(fmacrogeometry,line);
				iniss.str(line);
				eof = line.find("*");
				comment = line.find("//");
				if(eof!=string::npos)
				{ 
					break; 
				}
				else if(comment!=string::npos)
				{
					continue;
				}
				else
				{
					iniss >> ltemp >> dtemp;
					Geometry.stepwisegap_l_K.push_back(ltemp*1e-3);
					Geometry.stepwisegap_d_K.push_back(dtemp*1e-6);
				};
			};
			fmacrogeometry.clear();
			fmacrogeometry.close();
		};
	}
	//2D piston
	if(myinput.data.options_piston.general.McrK == 2){

		ifstream fmacrogeometry(myinput.data.options_piston.general.McrK_file);
		if (!fmacrogeometry.is_open()) 	{
			Log << "Unable to open " << myinput.data.options_piston.general.McrK_file << "! - Generate and restart..." << endl;
			//system("PAUSE");
			exit(1);
		}

		while(fmacrogeometry)
		{
			
			iniss.clear();
			getline(fmacrogeometry,line);
			iniss.str(line);
			eof = line.find("*");
			comment = line.find("//");
			if(eof!=string::npos)
			{ 
				break; 
			}
			else if(comment!=string::npos)
			{
				continue;
			}
			else if(line.find("Circumference")!=string::npos)
			{
				getline(fmacrogeometry,line);
				iniss.str(line);
				iniss >> dtemp;
				Geometry.PistonCirc = dtemp;
				//Log << dtemp << "\n";
			}
			else if(line.find("Axial")!=string::npos)
			{
				getline(fmacrogeometry,line);
				iniss.str(line);
				iniss >> dtemp;
				Geometry.PistonAx = dtemp;
				//Log << dtemp << "\n";
			}
			else
			{
				while(iniss){
					iniss >> dtemp;
					Geometry.pistonsurface.push_back(dtemp*1e-6);//m
					//Log << dtemp;
				};
				//Log << "\n";
			};
		};
		fmacrogeometry.clear();
		fmacrogeometry.close();
	};

}

void CGapInput::readMacroGeometryCylinder(void)
{
	double dtemp;
	double ltemp;
	string line;
	size_t comment;
	size_t eof;
	istringstream iniss;
	

	//Stepwise bushing
	if( myinput.data.options_piston.general.McrB == 1)
	{	
		ifstream fmacrogeometry(myinput.data.options_piston.general.McrB_file);
		if (!fmacrogeometry.is_open()) 	{
			Log << "Unable to open " << myinput.data.options_piston.general.McrB_file << "! - Generate and restart..." << endl;
			//system("PAUSE");
			exit(1);
		}
		//Stepwise bushing
		while(fmacrogeometry)
		{
			iniss.clear();
			getline(fmacrogeometry,line);
			iniss.str(line);
			eof = line.find("*");
			comment = line.find("//");
			if(eof!=string::npos)
			{ 
				break; 
			}
			else if(comment!=string::npos)
			{
				continue;
			}
			else
			{
				iniss >> ltemp >> dtemp;
				Geometry.stepwisegap_l_B.push_back(ltemp*1e-3);//m
				Geometry.stepwisegap_d_B.push_back(dtemp*1e-6);
			};
		};
		fmacrogeometry.clear();
		fmacrogeometry.close();
	};

	//2D bushing
	if(myinput.data.options_piston.general.McrB == 2){

		ifstream fmacrogeometry(myinput.data.options_piston.general.McrB_file);
		if (!fmacrogeometry.is_open()) 	{
			Log << "Unable to open " << myinput.data.options_piston.general.McrB_file << "! - Generate and restart..." << endl;
			//system("PAUSE");
			exit(1);
		}

		while(fmacrogeometry)
		{
			
			iniss.clear();
			getline(fmacrogeometry,line);
			iniss.str(line);
			eof = line.find("*");
			comment = line.find("//");
			if(eof!=string::npos)
			{ 
				break; 
			}
			else if(comment!=string::npos)
			{
				continue;
			}
			else if(line.find("Circumference")!=string::npos)
			{
				getline(fmacrogeometry,line);
				iniss.str(line);
				iniss >> dtemp;
				Geometry.BushingCirc = dtemp;
				//Log << dtemp << "\n";
			}
			else if(line.find("Axial")!=string::npos)
			{
				getline(fmacrogeometry,line);
				iniss.str(line);
				iniss >> dtemp;
				Geometry.BushingAx = dtemp;
				//Log << dtemp << "\n";
			}
			else
			{
				while(iniss){
					iniss >> dtemp;
					Geometry.bushingsurface.push_back(dtemp*1e-6);//m
					//Log << dtemp;
				};
				//Log << "\n";
			};
		};
		fmacrogeometry.clear();
		fmacrogeometry.close();
	};

}

void CGapInput::readBodySurfacexyzPressure()
{

	int size,i;
	ifstream fin;
	string line;
	size_t eof;
	istringstream iniss;
	string input;
	string appended;

	appended = "/IM_bin/";
	input = myinput.data.options_piston.general.IM_piston_path + "/im_piston_gap.bin";
	fin.open(input.c_str());
	if(fin.is_open())
		appended = "/";
	fin.clear();
	fin.close();
	
	//read piston surface faces coordinates
	input = myinput.data.options_piston.general.IM_piston_path + appended + "xyz_faces.dat";
	fin.open(input.c_str());
	if (!fin.is_open()) 	{
		Log << "Unable to open pressure piston xyz_faces.dat! - Perhaps the IM file path is wrong..." << endl;
		//system("PAUSE");
		exit(1);
	}
	getline(fin,line);
	iniss.clear();
	iniss.str(line);
	iniss >> size;
	//size coordinates array
	xyzfK_p.resize(size,3);	xyzfK_p=0.0;
	i=0;
	while(fin)
	{
		iniss.clear();
		getline(fin,line);
		iniss.str(line);
		eof = line.find("*");
		if(eof!=string::npos)
		{ 
			break; 
		}
		else
		{
			iniss >> xyzfK_p(i,0) >> xyzfK_p(i,1) >> xyzfK_p(i,2);
		};
		i++;
	};
	fin.clear();
	fin.close();

	//read piston surface nodes coordinates
	input = myinput.data.options_piston.general.IM_piston_path + appended + "xyz_nodes.dat";
	fin.open(input.c_str());
	if (!fin.is_open()) 	{
		Log << "Unable to open pressure piston xyz_nodes.dat! - Perhaps the IM file path is wrong..." << endl;
		Log << input << "\n";
		//system("PAUSE");
		exit(1);
	}
	getline(fin,line);
	iniss.clear();
	iniss.str(line);
	iniss >> size;
	//size coordinates array
	xyznK_p.resize(size,3);	xyznK_p=0.0;
	i=0;
	while(fin)
	{
		iniss.clear();
		getline(fin,line);
		iniss.str(line);
		eof = line.find("*");
		if(eof!=string::npos)
		{ 
			break; 
		}
		else
		{
			iniss >> xyznK_p(i,0) >> xyznK_p(i,1) >> xyznK_p(i,2);
		};
		i++;
	};
	fin.clear();
	fin.close();

	appended = "/IM_bin/";
	input = myinput.data.options_piston.general.IM_bushing_path + "/im_bushing_gap.bin";
	fin.open(input.c_str());
	if(fin.is_open())
		appended = "/";
	fin.clear();
	fin.close();

	//read cylinder surface faces coordinates
	input = myinput.data.options_piston.general.IM_bushing_path + appended + "xyz_faces.dat";
	fin.open(input.c_str());
	if (!fin.is_open()) 	{
		Log << "Unable to open pressure cylinder xyz_faces.dat! - Perhaps the IM file path is wrong..." << "\n";
		Log << input << endl;
		//system("PAUSE");
		exit(1);
	}
	getline(fin,line);
	iniss.clear();
	iniss.str(line);
	iniss >> size;
	//size coordinates array
	xyzfB_p.resize(size,3);	xyzfB_p=0.0;
	i=0;
	while(fin)
	{
		iniss.clear();
		getline(fin,line);
		iniss.str(line);
		eof = line.find("*");
		if(eof!=string::npos)
		{ 
			break; 
		}
		else
		{
			iniss >> xyzfB_p(i,0) >> xyzfB_p(i,1) >> xyzfB_p(i,2);
		};
		i++;
	};
	fin.clear();
	fin.close();

	//read cylinder surface nodes coordinates
	input = " ";
	input = myinput.data.options_piston.general.IM_bushing_path + appended + "xyz_nodes.dat";
	fin.open(input.c_str());
	if (!fin.is_open()) 	{
		Log << "Unable to open pressure cylinder xyz_nodes.dat! - Perhaps the IM file path is wrong..." << endl;
		//system("PAUSE");
		exit(1);
	}
	getline(fin,line);
	iniss.clear();
	iniss.str(line);
	iniss >> size;
	//size coordinates array
	xyznB_p.resize(size,3);	xyznB_p=0.0;
	i=0;
	while(fin)
	{
		iniss.clear();
		getline(fin,line);
		iniss.str(line);
		eof = line.find("*");
		if(eof!=string::npos)
		{ 
			break; 
		}
		else
		{
			iniss >> xyznB_p(i,0) >> xyznB_p(i,1) >> xyznB_p(i,2);
		};
		i++;
	};
	fin.clear();
	fin.close();

};
void CGapInput::readBodySurfacexyzThermal(string body,string body_path_th)
{

	int size,i;
	ifstream fin;
	string line;
	size_t eof;
	istringstream iniss;
	string input;

	//piston read coordinates
	if(body=="piston")
	{
		//read piston surface faces coordinates
		input = "./temp/piston/piston/xyz_faces.dat";// body_path_th + "/" + "xyz_faces.dat";
		fin.open(input.c_str());
		if (!fin.is_open()) 	{
			Log << "Unable to open thermal piston xyz_faces.dat! - Generate xyz_faces.dat and restart..." << endl;
			//system("PAUSE");
			exit(1);
		}
		getline(fin,line);
		iniss.str(line);
		iniss >> size;
		//size coordinates array
		xyzfK_th.resize(size,3); xyzfK_th=0.0;
		i=0;
		while(fin)
		{
			iniss.clear();
			getline(fin,line);
			iniss.str(line);
			eof = line.find("*");
			if(eof!=string::npos)
			{ 
				break; 
			}
			else
			{
				iniss >> xyzfK_th(i,0) >> xyzfK_th(i,1) >> xyzfK_th(i,2);
			};
			i++;
		};
		fin.clear();
		fin.close();

		//read piston surface nodes coordinates
		input = "./temp/piston/piston/xyz_nodes.dat"; //body_path_th + "/" + "xyz_nodes.dat";
		fin.open(input.c_str());
		if (!fin.is_open()) 	{
			Log << "Unable to open thermal piston xyz_nodes.dat! - Generate xyz_nodes.dat and restart..." << endl;
			//system("PAUSE");
			exit(1);
		}
		getline(fin,line);
		iniss.str(line);
		iniss >> size;
		//size coordinates array
		xyznK_th.resize(size,3); xyznK_th=0.0;
		i=0;
		while(fin)
		{
			iniss.clear();
			getline(fin,line);
			iniss.str(line);
			eof = line.find("*");
			if(eof!=string::npos)
			{ 
				break; 
			}
			else
			{
				iniss >> xyznK_th(i,0) >> xyznK_th(i,1) >> xyznK_th(i,2);
			};
			i++;
		};
		fin.clear();
		fin.close();
	};


	//cylinder read coordinates
	if(body=="cylinder")
	{
		//read cylinder surface faces coordinates
		input = "./temp/piston/bushing/xyz_faces.dat"; //body_path_th + "/" + "xyz_faces.dat";
		fin.open(input.c_str());
		if (!fin.is_open()) 	{
			Log << "Unable to open thermal cylinder xyz_faces.dat! - Generate xyz_faces.dat and restart..." << endl;
			//system("PAUSE");
			exit(1);
		}
		getline(fin,line);
		iniss.str(line);
		iniss >> size;
		//size coordinates array
		xyzfB_th.resize(size,3); xyzfB_th=0.0;
		i=0;
		while(fin)
		{
			iniss.clear();
			getline(fin,line);
			iniss.str(line);
			eof = line.find("*");
			if(eof!=string::npos)
			{ 
				break; 
			}
			else
			{
				iniss >> xyzfB_th(i,0) >> xyzfB_th(i,1) >> xyzfB_th(i,2);
			};
			//xyzfB_th(i,2) -= meshshift;
			i++;
		};
		fin.clear();
		fin.close();

		//read cylinder surface nodes coordinates
		input = "./temp/piston/bushing/xyz_nodes.dat";// body_path_th + "/" + "xyz_nodes.dat";
		fin.open(input.c_str());
		if (!fin.is_open()) 	{
			Log << "Unable to open thermal cylinder xyz_nodes.dat! - Generate xyz_nodes.dat and restart..." << endl;
			//system("PAUSE");
			exit(1);
		}
		getline(fin,line);
		iniss.str(line);
		iniss >> size;
		//size coordinates array
		xyznB_th.resize(size,3); xyznB_th=0.0;
		i=0;
		while(fin)
		{
			iniss.clear();
			getline(fin,line);
			iniss.str(line);
			eof = line.find("*");
			if(eof!=string::npos)
			{ 
				break; 
			}
			else
			{
				iniss >> xyznB_th(i,0) >> xyznB_th(i,1) >> xyznB_th(i,2);
			};
			//xyznB_th(i,2) -= meshshift;
			i++;
		};
		fin.clear();
		fin.close();
	}

};
void CGapInput::readInfluenceMatricesPistonCylinder(string IM_piston_path,string IM_cylinder_path)
{
	//read influence matrices piston body
	int nFaces_gap,nNodes_gap,nIM;
	char node[50];
	ifstream fin;
	string number;
	string im_name;
	bool oldinflugen;
	short dcimflag = 0;

	//----------------read piston influence matrices in 1-D vector----------------//
	nFaces_gap = (int) xyzfK_p.extent(0);
	nNodes_gap = (int) xyznK_p.extent(0);
	//number of matrices: gap_faces + DC + CASE + socket_0 + socket_1
	nIM = nFaces_gap + 4;//4;
	//size influence matrices total container
	IM_piston.resize(nNodes_gap*nIM);
	IM_piston=0.0;

	oldinflugen = true;
	im_name = IM_piston_path + "/im_piston_gap.bin";
	fin.open(im_name.c_str(),std::ios::binary | std::ios::in);
	if(fin.is_open()){
		oldinflugen = false;
		//temporary array to contain matrix
		double* IM = new double[nNodes_gap*nFaces_gap];
		fin.read(reinterpret_cast<char*>(IM), nNodes_gap*nFaces_gap*sizeof(double));
		for(int j=0;j<nNodes_gap*nFaces_gap;j++)
		{
			IM_piston(j) = IM[j];
		}
		fin.clear();
		fin.close();
		//DC IM
		im_name = IM_piston_path + "/dc.bin";
		fin.open(im_name.c_str(),std::ios::binary | std::ios::in);
		if(!fin.is_open()){
			Log << "Cannot find Piston IM dc.bin. Perhaps the IM file path is wrong..." << endl;
			exit(1);
		}
		fin.read(reinterpret_cast<char*>(IM), nNodes_gap*nFaces_gap*sizeof(double));
		for(int j=0;j<nNodes_gap;j++)
		{
			IM_piston(j+nNodes_gap*nFaces_gap) = IM[j];
		}
		fin.clear();
		fin.close();
		//Case IM
		im_name = IM_piston_path + "/case.bin";
		fin.open(im_name.c_str(),std::ios::binary | std::ios::in);
		if(!fin.is_open()){
			Log << "Cannot find Piston IM case.bin. Perhaps the IM file path is wrong..." << endl;
			exit(1);
		}
		fin.read(reinterpret_cast<char*>(IM), nNodes_gap*nFaces_gap*sizeof(double));
		for(int j=0;j<nNodes_gap;j++)
		{
			IM_piston(j+nNodes_gap*(nFaces_gap+1)) = 0.01 * IM[j];
		}
		fin.clear();
		fin.close();
		//Socket0 IM
		im_name = IM_piston_path + "/socket0.bin";
		fin.open(im_name.c_str(),std::ios::binary | std::ios::in);
		if(!fin.is_open()){
			Log << "Cannot find Piston IM socket0.bin. Perhaps the IM file path is wrong..." << endl;
			//exit(1);
		}
		else
		{
			fin.read(reinterpret_cast<char*>(IM), nNodes_gap*nFaces_gap*sizeof(double));
			for(int j=0;j<nNodes_gap;j++)
			{
				IM_piston(j+nNodes_gap*(nFaces_gap+2)) = IM[j];
			}
		}
		fin.clear();
		fin.close();
		//Socket1 IM
		im_name = IM_piston_path + "/socket1.bin";
		fin.open(im_name.c_str(),std::ios::binary | std::ios::in);
		if(!fin.is_open()){
			Log << "Cannot find Piston IM socket1.bin. Perhaps the IM file path is wrong..." << endl;
			//exit(1);
		}
		else
		{
			fin.read(reinterpret_cast<char*>(IM), nNodes_gap*nFaces_gap*sizeof(double));
			for(int j=0;j<nNodes_gap;j++)
			{
				IM_piston(j+nNodes_gap*(nFaces_gap+3)) = IM[j];
			}
		}
		delete [] IM;
	}
	fin.clear();
	fin.close();

	if(oldinflugen){
		//reading loop over number of matrices
		for(int i=0;i<nIM;i++)
		{
			if(i<nIM-4)//4)
			{
				//IM name
				_itoa_s(i+1,node,10);
				number = node;
			}
			else if(i==(nIM-4))//4))
			{
				//IM name
				number = "DC";
			}
			else if(i==(nIM-3))//3))
			{
				//IM name
				number = "CASE";
			}
			else if(i==(nIM-2))
			{
				number = "socket_0";
			}
			else
			{
				number = "socket_1";
			};
			//temporary array to contain matrix
			double* IM = new double[nNodes_gap];
			//read binary form
			im_name = IM_piston_path + "/IM_bin/" + number + ".bin";
			fin.open(im_name.c_str(),std::ios::binary | std::ios::in);
			if (!fin.is_open())
			{
				Log << "Unable to open piston IM " + number + ".bin! - Perhaps the IM file path is wrong..." << endl;
				//system("PAUSE");
				exit(1);
			}
			fin.read(reinterpret_cast<char*>(IM), nNodes_gap*sizeof(double));
			//assign matrix to matrices container
			for(int j=0;j<nNodes_gap;j++)
			{
				IM_piston(i*nNodes_gap+j) =  IM[j];
			}

			//clear
			delete [] IM;
			fin.clear();
			fin.close();
		}
	}

	//-------------------read cylinder influence matrices in 1-D vector-----------------//
	nFaces_gap = (int) xyzfB_p.extent(0);
	nNodes_gap = (int) xyznB_p.extent(0);
	//number of matrices: gap_faces + (npistons-1)*gap + DC*npistons + CASE
	nIM = nFaces_gap + 2*myinput.data.operating_conditions.npistons;
	//size influence matrices total container
	IM_cylinder.resize(nNodes_gap*nIM);
	IM_cylinder=0.0;

	oldinflugen = true;
	im_name = IM_cylinder_path + "/im_bushing_gap.bin";
	fin.open(im_name.c_str(),std::ios::binary | std::ios::in);
	if(fin.is_open()){
		oldinflugen = false;
		//temporary array to contain matrix
		double* IM = new double[nNodes_gap*nFaces_gap];
		fin.read(reinterpret_cast<char*>(IM), nNodes_gap*nFaces_gap*sizeof(double));
		for(int j=0;j<nNodes_gap*nFaces_gap;j++)
		{
			IM_cylinder(j) = IM[j];
			//Log << IM_cylinder(j) << "\n";
		}
		fin.clear();
		fin.close();
		//other gap IMs
		for(short dcim = 2;dcim <= myinput.data.operating_conditions.npistons;dcim++){
			char temp[2];
			itoa(dcim,temp,10);
			im_name = IM_cylinder_path + "/cyl" + temp + ".bin";
			fin.open(im_name.c_str(),std::ios::binary | std::ios::in);
			if(!fin.is_open()){
				Log << "Cannot find Cylinder IM cyl" << dcim << ".bin. Perhaps the IM file path is wrong..." << "\n";
			}
			else{
				fin.read(reinterpret_cast<char*>(IM), nNodes_gap*nFaces_gap*sizeof(double));
				for(int j=0;j<nNodes_gap;j++)
				{
					IM_cylinder(j+nNodes_gap*(nFaces_gap+dcim-2)) = IM[j];
				}
			}
			fin.clear();
			fin.close();
		}
		//DC IM
		im_name = IM_cylinder_path + "/dc.bin";
		fin.open(im_name.c_str(),std::ios::binary | std::ios::in);
		if(!fin.is_open()){
			Log << "Cannot find Cylinder IM dc.bin. Perhaps the IM file path is wrong..." << "\n";
		}
		else{
			dcimflag++;
			fin.read(reinterpret_cast<char*>(IM), nNodes_gap*nFaces_gap*sizeof(double));
			for(int j=0;j<nNodes_gap;j++)
			{
				IM_cylinder(j+nNodes_gap*(nFaces_gap+myinput.data.operating_conditions.npistons-2)) = IM[j];
			}
		}
		fin.clear();
		fin.close();
		for(short dcim = 1;dcim <= myinput.data.operating_conditions.npistons;dcim++){
			char temp[2];
			itoa(dcim,temp,10);
			im_name = IM_cylinder_path + "/dc" + temp + ".bin";
			fin.open(im_name.c_str(),std::ios::binary | std::ios::in);
			if(!fin.is_open()){
				Log << "Cannot find Cylinder IM dc" << dcim << ".bin. Perhaps the IM file path is wrong..." << "\n";
			}
			else{
				if(dcim == 1) dcimflag++;
				fin.read(reinterpret_cast<char*>(IM), nNodes_gap*nFaces_gap*sizeof(double));
				for(int j=0;j<nNodes_gap;j++)
				{
					IM_cylinder(j+nNodes_gap*(nFaces_gap+myinput.data.operating_conditions.npistons+dcim-3)) = IM[j];
				}
			}
			fin.clear();
			fin.close();
		}
		//Case IM
		im_name = IM_cylinder_path + "/case.bin";
		fin.open(im_name.c_str(),std::ios::binary | std::ios::in);
		if(!fin.is_open()){
			Log << "Cannot find Cylinder IM case.bin. Perhaps the IM file path is wrong..." << endl;
			exit(1);
		}
		fin.read(reinterpret_cast<char*>(IM), nNodes_gap*nFaces_gap*sizeof(double));
		for(int j=0;j<nNodes_gap;j++)
		{ 
			IM_cylinder(j+nNodes_gap*(nFaces_gap+2*myinput.data.operating_conditions.npistons-2)) = 0.01 * IM[j];
		}
		delete [] IM;
	}
	fin.clear();
	fin.close();

	if(oldinflugen){

		//reading loop over number of matrices
		for(int i=0;i<nIM-2*myinput.data.operating_conditions.npistons;i++)
		{
			//if(i<nIM-2)
			//{
				//IM name
				_itoa_s(i+1,node,10);
				number = node;
			//}
			/*else if(i==(nIM-2))
			{
				//IM name
				number = "DC";
			}
			else
			{
				//IM name
				number = "CASE";
			};*/
			//temporary array to contain matrix
			double* IM = new double[nNodes_gap];
			//read binary form
			im_name = IM_cylinder_path + "/IM_bin/" + number + ".bin";
			fin.open(im_name.c_str(),std::ios::binary | std::ios::in);
			if (!fin.is_open())
			{
				Log << "Unable to open cylinder IM " + number + ".bin! - Perhaps the IM file path is wrong..." << endl;
				//system("PAUSE");
				exit(1);
			}
			fin.read(reinterpret_cast<char*>(IM), nNodes_gap*sizeof(double));
			//assign matrix to matrices container
			for(int j=0;j<nNodes_gap;j++)
			{
				IM_cylinder(i*nNodes_gap+j) = IM[j];
			}
			//clear
			delete [] IM;
			fin.clear();
			fin.close();
		}

		//other gap ims
		for(int i = 2;i<=myinput.data.operating_conditions.npistons;i++){
			double* IM = new double[nNodes_gap];
			//read binary form
			char numb[2];
			itoa(i,numb,10);
			im_name = IM_cylinder_path + "/IM_bin/Cyl" + numb + ".bin";
			fin.open(im_name.c_str(),std::ios::binary | std::ios::in);
			if (!fin.is_open())
			{
				Log << "Unable to open cylinder Cyl" << i << ".bin! - Perhaps the IM file path is wrong..." << "\n";
				//system("PAUSE");
			}
			fin.read(reinterpret_cast<char*>(IM), nNodes_gap*sizeof(double));
			//assign matrix to matrices container
			for(int j=0;j<nNodes_gap;j++)
			{
				IM_cylinder(j+nNodes_gap*(nFaces_gap+i-2)) = IM[j];
			}
			//clear
			delete [] IM;
			fin.clear();
			fin.close();
		}

		//dc matrices
		double* IM = new double[nNodes_gap];
		//read binary form
		im_name = IM_cylinder_path + "/IM_bin/DC.bin";
		fin.open(im_name.c_str(),std::ios::binary | std::ios::in);
		if (!fin.is_open())
		{
			Log << "Unable to open cylinder DC.bin! - Perhaps the IM file path is wrong..." << "\n";
			//system("PAUSE");
		}
		else{
			dcimflag++;
			fin.read(reinterpret_cast<char*>(IM), nNodes_gap*sizeof(double));
			//assign matrix to matrices container
			for(int j=0;j<nNodes_gap;j++)
			{
				IM_cylinder(j+nNodes_gap*(nFaces_gap+myinput.data.operating_conditions.npistons-2)) = IM[j];
			}
		}
		//clear
		//delete [] IM;
		fin.clear();
		fin.close();
		for(int i = 1;i<=myinput.data.operating_conditions.npistons;i++){
			double* IM = new double[nNodes_gap];
			//read binary form
			char numb[2];
			itoa(i,numb,10);
			im_name = IM_cylinder_path + "/IM_bin/dc" + numb + ".bin";
			fin.open(im_name.c_str(),std::ios::binary | std::ios::in);
			if (!fin.is_open())
			{
				Log << "Unable to open cylinder dc" << i << ".bin! - Perhaps the IM file path is wrong..." << "\n";
				//system("PAUSE");
			}
			else{
				if(i == 1) dcimflag++;
				fin.read(reinterpret_cast<char*>(IM), nNodes_gap*sizeof(double));
				//assign matrix to matrices container
				for(int j=0;j<nNodes_gap;j++)
				{
					IM_cylinder(j+nNodes_gap*(nFaces_gap+myinput.data.operating_conditions.npistons+i-3)) = IM[j];
				}
			}
			//clear
			delete [] IM;
			fin.clear();
			fin.close();
		}

		//case matrix
		//double* IM = new double[nNodes_gap];
		//read binary form
		im_name = IM_cylinder_path + "/IM_bin/CASE.bin";
		fin.open(im_name.c_str(),std::ios::binary | std::ios::in);
		if (!fin.is_open())
		{
			Log << "Unable to open cylinder CASE.bin! - Perhaps the IM file path is wrong..." << endl;
			//system("PAUSE");
			exit(1);
		}
		fin.read(reinterpret_cast<char*>(IM), nNodes_gap*sizeof(double));
		//assign matrix to matrices container
		for(int j=0;j<nNodes_gap;j++)
		{
			IM_cylinder(j+nNodes_gap*(nFaces_gap+2*myinput.data.operating_conditions.npistons-2)) = IM[j];
		}
		//clear
		delete [] IM;
		fin.clear();
		fin.close();

	}
	//check to make sure one and only one dc1 IM has been read.
	if(dcimflag != 1){
		Log << "Found " << dcimflag << " IMs for DC 1. Check that you have exactly one IM for DC 1!" << endl;
		exit(1);
	}
};
//dwm Read FTG.txt file at the beginning of each revolution.
void CGapInput::readFTG(void)
{
	//clear previous contents
	FTGFile.FTG.clear();
	//FTGFile.phi.clear();
	FTGFile.time.clear();

	char str[256];
	char cbuffer;
	double ttemp;
	double phitemp;
	double FTGtemp;

	ifstream fFTGFile("./output/slipper/ftg.txt");

	//if (!fFTGFile.is_open()) 	{
	//cout << "Unable to open pFile.dat! - Generate pFile.dat and restart..." << "\n";
	//system("PAUSE");
	//exit(1);
	//	FTGFile.time.push_back(0);
	//	FTGFile.phi.push_back(0);
	//	FTGFile.FTG.push_back(0);
	//}
	if (fFTGFile.is_open()){
		while(!fFTGFile.eof()) 
		{
			if( fFTGFile.get(cbuffer) && (cbuffer =='/'))
			{
				fFTGFile.getline(str,sizeof(str));continue;
			}
			else
			{
				fFTGFile.putback(cbuffer);
			}
			
			fFTGFile >> ttemp >> phitemp >> FTGtemp;
			FTGFile.time.push_back(ttemp);
			//FTGFile.phi.push_back(phitemp);
			FTGFile.FTG.push_back(FTGtemp);
		}
	fFTGFile.close();
	};
};

void CGapInput::readVPFlux(void)
{
	//clear previous contents
	VPFluxFile.qvp = 0.0;
	
	char str[256];
	char cbuffer;
	double fluxtemp;

	ifstream fVPFluxFile("./output/block/block_flux.txt");

	if (fVPFluxFile.is_open()){
		while(!fVPFluxFile.eof()) 
		{
			if( fVPFluxFile.get(cbuffer) && (cbuffer =='/'))
			{
				fVPFluxFile.getline(str,sizeof(str));continue;
			}
			else
			{
				fVPFluxFile.putback(cbuffer);
			}
			
			fVPFluxFile >> fluxtemp;
			VPFluxFile.qvp = fluxtemp;
		}
	fVPFluxFile.close();
	//Log << "Reading VP Flux: " << VPFluxFile.qvp << "\n";
	};
};