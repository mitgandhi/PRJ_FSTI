#include "CSlipperGap.h"
//#include <omp.h>  //Replaced this line of code with the following line because there were conflicting versions of omp.h after I switched to the Intel compiler. 
//One version of omp.h from Visual Studio. Another version of omp.h from Intel XE. 
//Visual Studio was including the Intel XE version, but the slipper code was written using the Visual Studio version syntax. 
#include <c:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\include\omp.h>
#include <time.h>
#include <map>
#include <set>

#pragma once
#define OMP

extern void matlab(const Array<double,2>& data,const string file);

void CSlipperGap::SlipperCalch(vector<double> &xg)
{
	double xg1,xg2,xg3;
	double R = geometryslippergap.routG;
	
	xg1 = xg[0];
	xg2 = xg[1];
	xg3 = xg[2];

	Fluid.hrigid = Fluid.r*sin(Fluid.theta)*sqrt(1.0/3)/R*(xg2-xg3)+Fluid.r*cos(Fluid.theta)/(3*R)*(2*xg1-xg2-xg3)+(xg1+xg2+xg3)/3;

	//Macro geometry
	if(gapinput->options_slipper.general.EnableSlipperMacro == 1)
	{
		Fluid.hrigid += geometryslippergap.SlipperMacro;
	}

	//thermal deformation (goes into hrigid)
	if(gapinput->options_slipper.general.SlipperThermoElastic == 1)
	{
		Fluid.hrigid += t_slipper.pdeform;
	}
	if(gapinput->options_slipper.general.SwashplateThermoElastic)
	{
		//Currently swashplate thermoelastic analyis is used for the thermal boundary only
		//Fluid.hrigid -= t_swashplate.pdeform;
	}

	//Groove geometry
	Fluid.hrigid += Fluid.hgroove*(5e-6);
	
	if(gapinput->options_slipper.general.SlipperPressureDeformation == 1)
	{
		Fluid.h = Fluid.hrigid + slipper.ehd - swashplate.ehd;		
	} else {
		Fluid.h = Fluid.hrigid;
	}

	Fluid.z = Fluid.h(tensor::i,tensor::j)/(Fluid.Q-1)*tensor::k;
}
void CSlipperGap::SlipperCalch(void)
{
	SlipperCalch(control_point_position);
}
void CSlipperGap::SlipperCalcdht(void)
{
	double vg1,vg2,vg3;
	double R = geometryslippergap.routG;
	
	vg1 = control_point_velocity[0];
	vg2 = control_point_velocity[1];
	vg3 = control_point_velocity[2];

	Fluid.dht = Fluid.r*sin(Fluid.theta)*sqrt(1.0/3)/R*(vg2-vg3)+Fluid.r*cos(Fluid.theta)/(3*R)*(2*vg1-vg2-vg3)+(vg1+vg2+vg3)/3;

	//Fluid.ReyVals(0,all,all) = Fluid.dht;

	if(gapinput->options_slipper.general.EHDsqueeze == 1)
	{
		if(sum(slipper.oldehd) != 0)
		{
			
			Array<double, 2> ehdsqz(Fluid.M,Fluid.N);
			ehdsqz = ( (slipper.ehd-slipper.oldehd) - (swashplate.ehd-swashplate.oldehd) )/gapinput->common.timestep;
			//ehdsqz = ( (slipper.ehd-slipper.oldehd) )/gapinput->common.timestep;

			//Fluid.ReyVals(1,all,all) = (slipper.ehd-slipper.oldehd)/GapInput.General.TimeStep;
			//Fluid.ReyVals(2,all,all) = (- (swashplate.ehd-swashplate.oldehd))/GapInput.General.TimeStep;

			/*
			if(alpha != 1 && sum(Fluid.ehdsqzOld) != 0)
			{
				Array<double, 2> dehdsqz((ehdsqz - Fluid.ehdsqzOld) / Fluid.ehdsqzOld);
				
				cout << "eo: " << sum(Fluid.ehdsqzOld) << "\t";
				cout << "e1: " << sum(ehdsqz) << "\t";
				
				cout << "m: " << max(abs(dehdsqz)) << "\t";
				double lim = 1.0;
				if(alpha != 1 && max(abs(dehdsqz)) > lim)
				{
					double a = lim/max(fabs(dehdsqz));
					
					cout << "a: " << n2s(a) << "\t";
					ehdsqz = a*(ehdsqz - Fluid.ehdsqzOld) + Fluid.ehdsqzOld;
					cout << "e2: " << sum(ehdsqz) << endl;
				}
			}
			*/

			
			//ehdsqz = Fluid.ehdsqzOld - alpha * ( Fluid.ehdsqzOld - ehdsqz );
			

			Fluid.dht += ehdsqz;
			//Fluid.ehdsqzOld = ehdsqz;
		}
	}
}
void CSlipperGap::sSolid::loadIM()
{
	//test to see if we should still load the old ig1 format
	{
		ifstream n((path + "nodes." + suffix).c_str());
		if(n.is_open())
		{
			n.close();
			loadIM_ig1();
			return;
		}
	}

	//load the new ig2 format
	
	//Will be used to renumber the nodes in faces
	map<int,int> NodeNumbers;

	ifstream n((path + "nodes.txt").c_str());
	{
		nodes.clear();

		string line;
		int id = 0;
		while(getline(n,line)) 
		{
			vector<string> data;
			
			Tokenize(line,data,"\t");
			if(data.size() == 3)
			{
				nodes.push_back(node(
											atof(data[1].c_str()),
											atof(data[2].c_str())
											));
				NodeNumbers[atoi(data[0].c_str())] = id;
				id++;
			}
		}
		
	}
	n.close();
	nodecnt = (int) nodes.size();

	if(nodecnt < 1)
	{
		GapLog->message("ERROR: The nodes file " + path + "nodes.txt does not seem to exist or is empty!");
		exit(1);
	}

	ifstream f((path + "faces.txt").c_str());
	{
		faces.clear();

		string line;
		while(getline(f,line)) 
		{
			vector<string> data;
			Tokenize(line,data,"\t");
			if(data.size() == 5)
			{
				faces.push_back(face(
											NodeNumbers.find(atoi(data[0].c_str()))->second,
											NodeNumbers.find(atoi(data[1].c_str()))->second,
											NodeNumbers.find(atoi(data[2].c_str()))->second,
											atof(data[3].c_str()),
											atof(data[4].c_str())
											));
			}
		}
		
	}
	f.close();
	facecnt = (int) faces.size();

	if(facecnt < 1)
	{
		GapLog->message("ERROR: The nodes file " + path + "faces.txt does not seem to exist or is empty!");
		exit(1);
	}

	NodeNumbers.clear();

	//Allocate memory

	//For the influence matricies
	IM = new double * [facecnt];
	for(int i=0; i<facecnt; i++)
	{
		IM[i] = new double [nodecnt];
	}
	IMpocket = new double [nodecnt];
	IMcase = new double [nodecnt];
	IMsocket = new double [nodecnt];

	//Resize other sSolid members
	facePressure.resize(facecnt);
	nodeDeformation.resize(nodecnt);
	
	//Load the gap IM

	//First try to see if there is a merged gap.suffix file and if so load that
	fstream gapim((path + "im_" + suffix + "_gap.bin").c_str(), ios::in|ios::binary);
	if(gapim.is_open())
	{
		for(int i=0; i<facecnt; i++)
		{
			gapim.read(reinterpret_cast<char*> (&IM[i][0]), nodecnt*sizeof(double));
		}
	} else {	
		GapLog->message("ERROR! Unable to open: " + path + "im_" + suffix + "_gap.bin");
		exit(1);
	}
	gapim.close();

	//Load pocket IM
	{
		string file = path + "im_" + suffix + "_pocket.bin";
		fstream m(file.c_str(), ios::in|ios::binary);
		if(!m.is_open())
		{
			GapLog->message("SlipperGap -> INFO: " + suffix + " pocket IM failed to load. This is okay if there is no pocket IM.");

			//Zero the IM
			for(int i=0; i<nodecnt; i++)
			{
				IMpocket[i] = 0;
			}
		} else {
			m.read(reinterpret_cast<char*> (&IMpocket[0]), nodecnt*sizeof(double));
		}

		m.close();
	}


	//Load case IM
	{

		string file = path + "im_" + suffix + "_case.bin";
		fstream m(file.c_str(), ios::in|ios::binary);
		if(!m.is_open())
		{
			GapLog->message("SlipperGap -> INFO: " + suffix + " case IM failed to load. This is okay if there is no case IM.");
			
			//Zero the IM
			for(int i=0; i<nodecnt; i++)
			{
				IMcase[i] = 0;
			}

		} else {
			m.read(reinterpret_cast<char*> (&IMcase[0]), nodecnt*sizeof(double));
		}

		m.close();
	}

	//Load socket IM
	{

		string file = path + "im_" + suffix + "_socket.bin";
		fstream m(file.c_str(), ios::in|ios::binary);
		if(!m.is_open())
		{
			GapLog->message("SlipperGap -> INFO: " + suffix + " socket IM failed to load. This is okay if there is no socket IM.");
			
			//Zero the IM
			for(int i=0; i<nodecnt; i++)
			{
				IMsocket[i] = 0;
			}

		} else {
			m.read(reinterpret_cast<char*> (&IMsocket[0]), nodecnt*sizeof(double));
		}

		m.close();
	}

	//ReSize the faces KD tree
	KDfaces.points = annAllocPts(facecnt,KDfaces.dim);

	//Load the points
	for(int f=0; f<facecnt; f++)
	{
		KDfaces.points[f][0] = faces[f].x;
		KDfaces.points[f][1] = faces[f].y;
	}

	//Build the tree
	KDfaces.kdtree = new ANNkd_tree(KDfaces.points, facecnt, KDfaces.dim);

	//Find avg IM
	avgIM = 0;
	int avgCnt = 0;
	for(int i=0; i<facecnt; i++)
	{
		for(int j=0; j<3; j++)
		{
			avgIM += IM[i][faces[i].nodes[j]];
			avgCnt++;
		}
	}
	avgIM /= (double) avgCnt;

	GapLog->message("Average IM stiffness of " + suffix + " = " + n2s(avgIM));

}
void CSlipperGap::sSolid::loadIM_ig1()
{
	//test to see if we should update the path to append /output/ for legacy IMgen code
	{
		ifstream n((path + "nodes." + suffix).c_str());
		if(n.is_open())
		{
			n.close();
		} else {
			ifstream n2((path + "/output/nodes." + suffix).c_str());
			if(n2.is_open())
			{
				n2.close();
				path += "/output/";
			}
		}
	}

	//Will be used to renumber the nodes in faces
	map<int,int> NodeNumbers;

	ifstream n((path + "nodes." + suffix).c_str());
	{
		nodes.clear();

		string line;
		int id = 0;
		while(getline(n,line)) 
		{
			vector<string> data;
			
			Tokenize(line,data,"\t");
			if(data.size() == 3)
			{
				nodes.push_back(node(
											atof(data[1].c_str()),
											atof(data[2].c_str())
											));
				NodeNumbers[atoi(data[0].c_str())] = id;
				id++;
			}
		}
		
	}
	n.close();
	nodecnt = (int) nodes.size();

	if(nodecnt < 1)
	{
		GapLog->message("ERROR: The nodes file " + path + "nodes." + suffix + " does not seem to exist or is empty!");
		exit(1);
	}

	ifstream f((path + "faces." + suffix).c_str());
	{
		faces.clear();

		string line;
		while(getline(f,line)) 
		{
			vector<string> data;
			Tokenize(line,data,"\t");
			if(data.size() == 5)
			{
				faces.push_back(face(
											NodeNumbers.find(atoi(data[0].c_str()))->second,
											NodeNumbers.find(atoi(data[1].c_str()))->second,
											NodeNumbers.find(atoi(data[2].c_str()))->second,
											atof(data[3].c_str()),
											atof(data[4].c_str())
											));
			}
		}
		
	}
	f.close();
	facecnt = (int) faces.size();

	if(facecnt < 1)
	{
		GapLog->message("ERROR: The nodes file " + path + "faces." + suffix + " does not seem to exist or is empty!");
		exit(1);
	}

	NodeNumbers.clear();

	//Allocate memory

	//For the influence matricies
	IM = new double * [facecnt];
	for(int i=0; i<facecnt; i++)
	{
		IM[i] = new double [nodecnt];
	}
	IMpocket = new double [nodecnt];
	IMcase = new double [nodecnt];
	IMsocket = new double [nodecnt];

	//Resize other sSolid members
	facePressure.resize(facecnt);
	nodeDeformation.resize(nodecnt);
	
	//Load the gap IM

	//First try to see if there is a merged gap.suffix file and if so load that
	fstream gapim((path + "gap." + suffix).c_str(), ios::in|ios::binary);
	if(gapim.is_open())
	{
		for(int i=0; i<facecnt; i++)
		{
			char c;

			gapim.get(c);
			if (char(2) != c)
			{
				GapLog->message("Error in open byte of " + n2s(i) + " gap IM! Aborting.");
				exit(1);
			}

			gapim.read(reinterpret_cast<char*> (&IM[i][0]), nodecnt*sizeof(double));

			gapim.get(c);
			if (char(3) != c)
			{
				GapLog->message("Error in close byte of " + n2s(i) + " gap IM! Aborting.");
				exit(1);
			}
		}
	} else {
		//couldn't read the gap.suffix file so try and load individual gap IM's

		for(int i=0; i<facecnt; i++)
		{
			string file = path + n2s(i) + "." + suffix;
			fstream m(file.c_str(), ios::in|ios::binary);
			if(!m.is_open())
			{
				GapLog->message("ERROR! Unable to open: " + file);
				exit(1);
			}

			char c;

			m.get(c);
			if (char(2) != c)
			{
				GapLog->message("Error in open byte of " + file + " file! Aborting.");
				exit(1);
			}

			m.read(reinterpret_cast<char*> (&IM[i][0]), nodecnt*sizeof(double));

			m.get(c);
			if (char(3) != c)
			{
				GapLog->message("Error in close byte of " + file + " file! Aborting.");
				exit(1);
			}

			m.close();
		}
	}
	gapim.close();

	//Load pocket IM
	{

		string file = path + "pocket." + suffix;
		fstream m(file.c_str(), ios::in|ios::binary);
		if(!m.is_open())
		{
			GapLog->message("SlipperGap -> INFO: " + suffix + " pocket IM failed to load. This is okay if there is no pocket IM.");

			//Zero the IM
			for(int i=0; i<nodecnt; i++)
			{
				IMpocket[i] = 0;
			}
		} else {
			char c;
			m.get(c);
			if (char(2) != c)
			{
				*GapLog << "Error in open byte of " << file << " file! Aborting." << endl;
				exit(1);
			}

			m.read(reinterpret_cast<char*> (&IMpocket[0]), nodecnt*sizeof(double));

			m.get(c);
			if (char(3) != c)
			{
				*GapLog << "Error in close byte of " << file << " file! Aborting." << endl;
				exit(1);
			}
		}

		m.close();
	}


	//Load case IM
	{

		string file = path + "case." + suffix;
		fstream m(file.c_str(), ios::in|ios::binary);
		if(!m.is_open())
		{
			GapLog->message("SlipperGap -> INFO: " + suffix + " case IM failed to load. This is okay if there is no case IM.");
			
			//Zero the IM
			for(int i=0; i<nodecnt; i++)
			{
				IMcase[i] = 0;
			}

		} else {
			char c;
			m.get(c);
			if (char(2) != c)
			{
				*GapLog << "Error in open byte of " << file << " file! Aborting." << endl;
				exit(1);
			}

			m.read(reinterpret_cast<char*> (&IMcase[0]), nodecnt*sizeof(double));

			m.get(c);
			if (char(3) != c)
			{
				*GapLog << "Error in close byte of " << file << " file! Aborting." << endl;
				exit(1);
			}
		}

		m.close();
	}

	//Load socket IM
	{

		string file = path + "socket." + suffix;
		fstream m(file.c_str(), ios::in|ios::binary);
		if(!m.is_open())
		{
			GapLog->message("SlipperGap -> INFO: " + suffix + " socket IM failed to load. This is okay if there is no socket IM.");
			
			//Zero the IM
			for(int i=0; i<nodecnt; i++)
			{
				IMsocket[i] = 0;
			}

		} else {
			char c;
			m.get(c);
			if (char(2) != c)
			{
				*GapLog << "Error in open byte of " << file << " file! Aborting." << endl;
				exit(1);
			}

			m.read(reinterpret_cast<char*> (&IMsocket[0]), nodecnt*sizeof(double));

			m.get(c);
			if (char(3) != c)
			{
				*GapLog << "Error in close byte of " << file << " file! Aborting." << endl;
				exit(1);
			}
		}

		m.close();
	}

	//ReSize the faces KD tree
	KDfaces.points = annAllocPts(facecnt,KDfaces.dim);

	//Load the points
	for(int f=0; f<facecnt; f++)
	{
		KDfaces.points[f][0] = faces[f].x;
		KDfaces.points[f][1] = faces[f].y;
	}

	//Build the tree
	KDfaces.kdtree = new ANNkd_tree(KDfaces.points, facecnt, KDfaces.dim);

	//Find avg IM
	avgIM = 0;
	int avgCnt = 0;
	for(int i=0; i<facecnt; i++)
	{
		for(int j=0; j<3; j++)
		{
			avgIM += IM[i][faces[i].nodes[j]];
			avgCnt++;
		}
	}
	avgIM /= (double) avgCnt;

	GapLog->message("Average IM stiffness of " + suffix + " = " + n2s(avgIM));

	//Modify the gap IM to account for high local stiffness
	if(suffix.compare("slipper") == 0 && false)
	{
		*GapLog << "Modifying the local gap stiffness for the " << suffix << " IM ... ";
		int modcnt = 0;

		//First build a map of node -> element to find neighbor elements
		vector< vector<int> > n2e(nodecnt);
		
		//populate n2e with empty vectors
		for(int n=0; n<nodecnt; n++)
		{
			vector<int> eleids(0);
			n2e[n] = eleids;
		}
		
		//assign n2e
		for(int f=0; f<faces.size(); f++)
		{
			for(int n=0; n<3; n++)
			{
				n2e[ faces[f].nodes[n] ].push_back(f);
			}
		}

		//update each IM
		for(int f=0; f<faces.size(); f++)
		{
			//keep track of which nodes have already been modified
			set<int> modnode;

			//modify the primary nodes
			for(int n=0; n<3; n++)
			{
				const int nid = faces[f].nodes[n];
				modnode.insert(nid);
				IM[f][nid] *= 2.25;
				modcnt++;
			}

			//modify the secondary nodes
			for(int n=0; n<3; n++)
			{
				//the primary node id
				const int pnid = faces[f].nodes[n];
				
				//loop through all the attached elements
				for(int e=0; e<n2e[pnid].size(); e++)
				{
					//the attached element id
					const int eid = n2e[pnid][e];

					//for each attached element, loop through its nodes
					for(int sn=0; sn<3; sn++)
					{
						//the secondary node id
						const int snid = faces[eid].nodes[n];

						//check if hasn't been modified
						if(modnode.count(snid) == 0)
						{
							IM[f][snid] *= 1.05;
							modnode.insert(snid);
							modcnt++;
						}
					}
				}
			}
			
		}

		*GapLog << "done! (modified " << modcnt << " values)" << endl;
	}

}
void CSlipperGap::sSolid::calcDeform(double pG, double pcase, double psocket)
{
	//The simplest way possible to start
	
	
	//Pocket / case

	//Scale the pressure
	pG = pG/100e+5;
	pcase = pcase/100e+5;
	psocket = psocket/100e+5;
	for(int n=0; n<nodecnt; n++)
	{
		nodeDeformation[n] = IMpocket[n]*pG+IMcase[n]*pcase+IMsocket[n]*psocket;
	}
	
	/*

	//Single Thread Version
	//Note 2x speedup can be gained by using a simple nD local var
	for(int f=0; f<facecnt; f++)
	{
		const double scale = facePressure[f]/100e+5;
		
		for(int n=0; n<nodecnt; n++)
		{
			nodeDeformation[n] += IM[f][n]*scale;	
		}
	}
	*/

	
	//OpenMP Version
	int NT = 2;
	omp_set_num_threads(NT);
	double ** nD = new double * [NT];
	for(int t=0; t<NT; t++)
	{
		nD[t] = new double [nodecnt];
		for(int n=0; n<nodecnt; n++)
		{
			nD[t][n] = 0;
		}
	}

	#pragma omp parallel
	{
		#pragma omp for
		for(int f=0; f<facecnt; f++)
		{
			const int t = omp_get_thread_num();
			const double scale = facePressure[f]/100e+5;
			
			for(int n=0; n<nodecnt; n++)
			{
				nD[t][n] += IM[f][n]*scale;	
			}
		}

	}

	for(int n=0; n<nodecnt; n++)
	{
		for(int t=0; t<NT; t++)
		{
			nodeDeformation[n] += nD[t][n];
		}
	}
	
	for(int t=0; t<NT; t++)
	{
		delete [] nD[t];
	}
	delete [] nD;

}
void CSlipperGap::sSolid::writevtk(const string file)
{
		ofstream f(file.c_str());

		f << "# vtk DataFile Version 2.0" << endl;
		f << "# comment" << endl;
		f << "ASCII" << endl;
		f << "DATASET UNSTRUCTURED_GRID" << endl;
		f << "POINTS " << nodecnt << " double" << endl;
		
		for(int n=0; n<nodecnt; n++)
		{
			f << setprecision(10);
			f << scientific	<< nodes[n].x << "\t"
									<< nodes[n].y << "\t"
									<< "0.0";
			f << endl;
		}

		f << "CELLS " << facecnt << " " << facecnt * 4 << endl;
		for(int c=0; c<facecnt; c++)
		{
			f << "3\t";
			for(int i=0; i<3; i++)
			{
				f << faces[c].nodes[i] << "\t";
			}
			f << endl;
			
		}

		f << "CELL_TYPES " << facecnt << endl;
		for(int i=0; i<facecnt; i++)
		{
			f << "5" << endl;
		}

		f << "CELL_DATA " << facecnt << endl;
		f << "SCALARS pressure float 1" << endl;
		f << "LOOKUP_TABLE press" << endl;
		for(int i=0; i<facecnt; i++)
		{
			f << scientific << facePressure[i] << endl;
		}

		 
		f << "POINT_DATA " << nodecnt << endl;
		f << "SCALARS ehd float 1" << endl;
		f << "LOOKUP_TABLE ehd" << endl;
		for(int i=0; i<nodecnt; i++)
		{
			f << scientific << nodeDeformation[i] << endl;
		}
		
		f.close();
}
void CSlipperGap::sSolid::writevtk(const string file, const vector<double> &qflux)
{
		ofstream f(file.c_str());

		f << "# vtk DataFile Version 2.0" << endl;
		f << "# comment" << endl;
		f << "ASCII" << endl;
		f << "DATASET UNSTRUCTURED_GRID" << endl;
		f << "POINTS " << nodecnt << " double" << endl;
		
		for(int n=0; n<nodecnt; n++)
		{
			f << setprecision(10);
			f << scientific	<< nodes[n].x << "\t"
									<< nodes[n].y << "\t"
									<< "0.0";
			f << endl;
		}

		f << "CELLS " << facecnt << " " << facecnt * 4 << endl;
		for(int c=0; c<facecnt; c++)
		{
			f << "3\t";
			for(int i=0; i<3; i++)
			{
				f << faces[c].nodes[i] << "\t";
			}
			f << endl;
			
		}

		f << "CELL_TYPES " << facecnt << endl;
		for(int i=0; i<facecnt; i++)
		{
			f << "5" << endl;
		}

		f << "CELL_DATA " << facecnt << endl;
		f << "SCALARS pressure float 1" << endl;
		f << "LOOKUP_TABLE press" << endl;
		for(int i=0; i<facecnt; i++)
		{
			f << scientific << facePressure[i] << endl;
		}

		if(facecnt == qflux.size())
		{
			f << "SCALARS q_flux float 1" << endl;
			f << "LOOKUP_TABLE flux" << endl;
			for(int i=0; i<facecnt; i++)
			{
				f << scientific << qflux[i] << endl;
			}
		}
		 
		f << "POINT_DATA " << nodecnt << endl;
		f << "SCALARS ehd float 1" << endl;
		f << "LOOKUP_TABLE ehd" << endl;
		for(int i=0; i<nodecnt; i++)
		{
			f << scientific << nodeDeformation[i] << endl;
		}
		
		f.close();
}
void CSlipperGap::Fluid2Slipper()
{
	ANNpoint qPt = annAllocPt(2);
	//Interpolation points
	int k = 1;
	ANNidxArray	  nnIdx = new ANNidx[k];
	ANNdistArray  dists = new ANNdist[k];

	for(int f=0; f<slipper.facecnt; f++)
	{

		qPt[0] = slipper.faces[f].x;
		qPt[1] = slipper.faces[f].y;

		Fluid.KDslip.kdtree->annkSearch
		(
			qPt,
			k,
			nnIdx,
			dists,
			Fluid.KDslip.eps
		);	

		//The closest fluid volume centroid
		int c1 = nnIdx[0];
		int m1 = nnIdx[0]/Fluid.N;
		int n1 = nnIdx[0]%Fluid.N;

		//We need to find two other points to use
		int m2, n2;
		int m3, n3;

		double x = qPt[0];
		double y = qPt[1];
		double x1 = Fluid.LRx(m1,n1);
		double y1 = Fluid.LRy(m1,n1);

		//to polar!
		double r = sqrt(x*x+y*y);
		double t = atan2(y,x);
		double r1 = sqrt(x1*x1+y1*y1);
		double t1 = atan2(y1,x1);

		if(r>r1)
		{
			m2 = m1+1;
			//Special checks to use the 2nd inner point if at a boundary
			m2 = m2<0 ? 1 : m2;
			m2 = m2>=Fluid.M ? Fluid.M-2 : m2;

			n2 = n1;
		} else {
			m2 = m1-1;
			//Special checks to use the 2nd inner point if at a boundary
			m2 = m2<0 ? 1 : m2;
			m2 = m2>=Fluid.M ? Fluid.M-2 : m2;

			n2 = n1;
		}
		
		double x2 = Fluid.LRx(m2,n2);
		double y2 = Fluid.LRy(m2,n2);

		if(t>t1)
		{
			m3 = m1;
			n3 = n1+1;

			//Special trick if n is at a boundary since n is "wraped"
			n3 = (n3+Fluid.N)%Fluid.N;
		} else {
			m3 = m1;			
			n3 = n1-1;

			//Special trick if n is at a boundary since n is "wraped"
			n3 = (n3+Fluid.N)%Fluid.N;
		}

		double x3 = Fluid.LRx(m3,n3);
		double y3 = Fluid.LRy(m3,n3);

		//Barycentric Interpolation
		//http://en.wikipedia.org/wiki/Barycentric_coordinates_(mathematics)

		const double detT = (y2-y3)*(x1-x3)+(x3-x2)*(y1-y3);
		const double lam1 = ((y2-y3)*(x-x3)+(x3-x2)*(y-y3))/detT;
		const double lam2 = ((y3-y1)*(x-x3)+(x1-x3)*(y-y3))/detT;
		const double lam3 = 1.0 - lam1 - lam2;
		
		if(detT == 0)
		{
			//Two of the points are the same
			//This really should not happen with the current scheme

			//If it does an error needs to be thrown as interpolation will fail with NaN
		}

		if(lam1 < 0 || lam2 < 0 || lam3 < 0)
		{
			//The query point is not bounded by the interpolation points
			
			const double lamS = fabs(lam1) + fabs(lam2) + fabs(lam3);
		
			//This is only a big deal if lamS is significantly greater than 1.
			//A warning should be thrown if it is greater than 2-3??

			if(lamS > 3.0)
			{
				//cout << "WARNING: Forcing F2S f: " << f << " with lamS: " << lamS << endl;
				GapLog.message("WARNING: There is a mismatch between the fluid and slipper pressure deformation mesh. Forcing F2S interpolation at f=" + n2s(f) + ", lamS=" + n2s(lamS));
				// To fix this error, make sure that the mesh dimensions agree with doutG, dinG, and dDG in geometry.txt.
			}
		}

		//interp
		slipper.facePressure[f] =	lam1*Fluid.p(m1,n1) + 
											lam2*Fluid.p(m2,n2) + 
											lam3*Fluid.p(m3,n3);

		//limit check
		slipper.facePressure[f] = slipper.facePressure[f] < 1.0 ? 1.0 : slipper.facePressure[f];
	}

	annDeallocPt(qPt);
	delete [] nnIdx;
	delete [] dists;

}
void CSlipperGap::Fluid2Swashplate()
{
	ANNpoint qPt = annAllocPt(2);
	//Interpolation points
	int k = 1;
	ANNidxArray	  nnIdx = new ANNidx[k];
	ANNdistArray  dists = new ANNdist[k];

	swashplate.quickF2S.clear();

	for(int f=0; f<swashplate.facecnt; f++)
	{

		qPt[0] = swashplate.faces[f].x;
		qPt[1] = swashplate.faces[f].y;

		FullGroup.KDslip.kdtree->annkSearch
		(
			qPt,
			k,
			nnIdx,
			dists,
			FullGroup.KDslip.eps
		);	

		//The closest fluid volume centroid
		int steps = gapinput->common.rev_steps;
		int step = int(floor((gapinput->common.phi_rev_deg+gapinput->common.phi_deg_tol) / gapinput->common.phistep_deg)) % steps;
		
		int c1 = nnIdx[0];
		int s1 = nnIdx[0]/(Fluid.M*Fluid.N);
		int d1 = (step + s1*steps/geometryslippergap.npistons) % steps;
		int m1 = (nnIdx[0]%(Fluid.M*Fluid.N))/Fluid.N;
		int n1 = (nnIdx[0]%(Fluid.M*Fluid.N))%Fluid.N;

		//interp
		swashplate.facePressure[f] =	FullGroup.p(m1,n1,d1);
		
		if(d1 == step)
		{
			//Store for quick interp
			vector<int> tmp(3);
			tmp[0] = m1;
			tmp[1] = n1;
			tmp[2] = f;
			swashplate.quickF2S.push_back(tmp);
		}
	}

	annDeallocPt(qPt);
	delete [] nnIdx;
	delete [] dists;

}
void CSlipperGap::Swashplate2Fluid(const double alpha)
{
	//how many pts should be tried for interpolation boundedness
	int maxBoundTry = 10;
	ANNpoint qPt = annAllocPt(2);

	for(int m=0; m<Fluid.M; m++)
	{
		for(int n=0; n<Fluid.N; n++)
		{

			qPt[0] = Fluid.Gx(m,n);
			qPt[1] = Fluid.Gy(m,n);

			//Used to find the "best" pts in case of extrapolation
			vector<double> lamS;

			//Interpolation points
			for(int k = 1; k<=maxBoundTry+1; k++)
			{
				bool force = false;
				if(k == maxBoundTry+1)
				{
					//The closest `maxBoundTry` points didn't bound the query pt
					//so just force the min(sum(abs(lam))) one and extrapolate

					force = true;

					//Use whatever k value had the lowest lamS
					
					double minlamS = lamS[0];
					k = 1;
					for(unsigned int i=1; i<lamS.size(); i++)
					{
						if(lamS[i] < minlamS)
						{
							k = i+1;
							minlamS = lamS[i];
						}
					}
					
				}

				ANNidxArray	  nnIdx = new ANNidx[k];
				ANNdistArray  dists = new ANNdist[k];

				swashplate.KDfaces.kdtree->annkSearch
				(
					qPt,
					k,
					nnIdx,
					dists,
					swashplate.KDfaces.eps
				);	


				//Barycentric Interpolation
				//http://en.wikipedia.org/wiki/Barycentric_coordinates_(mathematics)

				const double x = qPt[0];
				const double y = qPt[1];
				const double x1 = swashplate.nodes[swashplate.faces[nnIdx[k-1]].nodes[0]].x;
				const double y1 = swashplate.nodes[swashplate.faces[nnIdx[k-1]].nodes[0]].y;
				const double x2 = swashplate.nodes[swashplate.faces[nnIdx[k-1]].nodes[1]].x;
				const double y2 = swashplate.nodes[swashplate.faces[nnIdx[k-1]].nodes[1]].y;
				const double x3 = swashplate.nodes[swashplate.faces[nnIdx[k-1]].nodes[2]].x;
				const double y3 = swashplate.nodes[swashplate.faces[nnIdx[k-1]].nodes[2]].y;

				const double detT = (y2-y3)*(x1-x3)+(x3-x2)*(y1-y3);
				const double lam1 = ((y2-y3)*(x-x3)+(x3-x2)*(y-y3))/detT;
				const double lam2 = ((y3-y1)*(x-x3)+(x1-x3)*(y-y3))/detT;
				const double lam3 = 1.0 - lam1 - lam2;

				if(force)
				{
					const double lamS = fabs(lam1) + fabs(lam2) + fabs(lam3);
					
					//This is only a big deal if lamS is significantly greater than 1.
					//A warning should be thrown if it is greater than 2-3??
					
					if(lamS > 3.0)
					{
						//Only throw an error message for non-boundary volumes
						if(Fluid.boundary(m,n) == -1)
						{
							//cout << "WARNING: Forcing S2F m: " << m << " n: " << n << " with lamS: " << lamS << endl;
							GapLog.message("WARNING: There is a mismatch between the fluid and swashplate pressure deformation mesh. Forcing S2F interpolation at m=" + n2s(m) + ", n=" + n2s(n) + ", lamS=" + n2s(lamS));
						}
					}
				}

				if( (lam1 < 0 || lam2 < 0 || lam3 < 0 || detT == 0) && !force )
				{
					lamS.push_back(fabs(lam1) + fabs(lam2) + fabs(lam3));
					
					delete [] nnIdx;
					delete [] dists;
				
					continue;
				}

				//Deformation
				double newehd =
								lam1*swashplate.nodeDeformation[swashplate.faces[nnIdx[k-1]].nodes[0]] + 
								lam2*swashplate.nodeDeformation[swashplate.faces[nnIdx[k-1]].nodes[1]] + 
								lam3*swashplate.nodeDeformation[swashplate.faces[nnIdx[k-1]].nodes[2]];

				//Add in the fake ehd
				//swashplate.fakeehd(m,n) = -Fluid.p(m,n) * 0.002 / 300e+9;
				//newehd += swashplate.fakeehd(m,n);

				//damp the new deformation
				swashplate.ehd(m,n) = swashplate.oldehd(m,n) + alpha*(newehd - swashplate.oldehd(m,n));

				//sucuessful so
				delete [] nnIdx;
				delete [] dists;
				break;
			}

		}
	}

	//clean up the ann pt
	annDeallocPt(qPt);
}
void CSlipperGap::Slipper2Fluid(const double alpha)
{
	//how many pts should be tried for interpolation boundedness
	int maxBoundTry = 10;
	ANNpoint qPt = annAllocPt(2);

	for(int m=0; m<Fluid.M; m++)
	{
		for(int n=0; n<Fluid.N; n++)
		{

			qPt[0] = Fluid.LRx(m,n);
			qPt[1] = Fluid.LRy(m,n);

			//Used to find the "best" pts in case of extrapolation
			vector<double> lamS;

			//Interpolation points
			for(int k = 1; k<=maxBoundTry+1; k++)
			{
				bool force = false;
				if(k == maxBoundTry+1)
				{
					//The closest `maxBoundTry` points didn't bound the query pt
					//so just force the min(sum(abs(lam))) one and extrapolate

					force = true;

					//Use whatever k value had the lowest lamS
					
					double minlamS = lamS[0];
					k = 1;
					for(unsigned int i=1; i<lamS.size(); i++)
					{
						if(lamS[i] < minlamS)
						{
							k = i+1;
							minlamS = lamS[i];
						}
					}
					
				}

				ANNidxArray	  nnIdx = new ANNidx[k];
				ANNdistArray  dists = new ANNdist[k];

				slipper.KDfaces.kdtree->annkSearch
				(
					qPt,
					k,
					nnIdx,
					dists,
					slipper.KDfaces.eps
				);	


				//Barycentric Interpolation
				//http://en.wikipedia.org/wiki/Barycentric_coordinates_(mathematics)

				const double x = qPt[0];
				const double y = qPt[1];
				const double x1 = slipper.nodes[slipper.faces[nnIdx[k-1]].nodes[0]].x;
				const double y1 = slipper.nodes[slipper.faces[nnIdx[k-1]].nodes[0]].y;
				const double x2 = slipper.nodes[slipper.faces[nnIdx[k-1]].nodes[1]].x;
				const double y2 = slipper.nodes[slipper.faces[nnIdx[k-1]].nodes[1]].y;
				const double x3 = slipper.nodes[slipper.faces[nnIdx[k-1]].nodes[2]].x;
				const double y3 = slipper.nodes[slipper.faces[nnIdx[k-1]].nodes[2]].y;

				const double detT = (y2-y3)*(x1-x3)+(x3-x2)*(y1-y3);
				const double lam1 = ((y2-y3)*(x-x3)+(x3-x2)*(y-y3))/detT;
				const double lam2 = ((y3-y1)*(x-x3)+(x1-x3)*(y-y3))/detT;
				const double lam3 = 1.0 - lam1 - lam2;

				if(force)
				{
					const double lamS = fabs(lam1) + fabs(lam2) + fabs(lam3);
					
					//This is only a big deal if lamS is significantly greater than 1.
					//A warning should be thrown if it is greater than 2-3??
					
					if(lamS > 3.0)
					{
						//Only throw an error message for non-boundary volumes
						if(Fluid.boundary(m,n) == -1)
						{
							//cout << "WARNING: Forcing S2F m: " << m << " n: " << n << " with lamS: " << lamS << endl;
							GapLog.message("WARNING: There is a mismatch between the fluid and slipper pressure deformation mesh. Forcing S2F interpolation at m=" + n2s(m) + ", n=" + n2s(n) + ", lamS=" + n2s(lamS));
						}
					}
				}

				if( (lam1 < 0 || lam2 < 0 || lam3 < 0 || detT == 0) && !force )
				{
					lamS.push_back(fabs(lam1) + fabs(lam2) + fabs(lam3));
					
					delete [] nnIdx;
					delete [] dists;
				
					continue;
				}

				//Deformation
				double newehd =
								lam1*slipper.nodeDeformation[slipper.faces[nnIdx[k-1]].nodes[0]] + 
								lam2*slipper.nodeDeformation[slipper.faces[nnIdx[k-1]].nodes[1]] + 
								lam3*slipper.nodeDeformation[slipper.faces[nnIdx[k-1]].nodes[2]];

				//damp the new deformation
				slipper.ehd(m,n) = slipper.oldehd(m,n) + alpha*(newehd - slipper.oldehd(m,n));

				//sucuessful so
				delete [] nnIdx;
				delete [] dists;
				break;
			}

		}
	}

	//clean up the ann pt
	annDeallocPt(qPt);
}
void CSlipperGap::LoadSlipperIM(void)
{

	slipper.path = gapinput->options_slipper.general.IMpath + "/slipper/";
	slipper.suffix = "slipper";
	slipper.M = Fluid.M;
	slipper.N = Fluid.N;
	slipper.loadIM();

}
void CSlipperGap::LoadSwashplateIM(void)
{

	swashplate.path = gapinput->options_slipper.general.IMpath + "/swashplate/";
	swashplate.suffix = "swashplate";
	swashplate.M = Fluid.M;
	swashplate.N = Fluid.N;
	swashplate.loadIM();

}
void CSlipperGap::SlipperVTK(const string file)
{
		ofstream f(file.c_str());


		f << "# vtk DataFile Version 2.0" << endl;
		f << "# comment" << endl;
		f << "ASCII" << endl;
		f << "DATASET UNSTRUCTURED_GRID" << endl;
		f << "POINTS " << slipper.nodecnt << " double" << endl;
		
		for(int n=0; n<slipper.nodecnt; n++)
		{
			f << setprecision(10);
			f << scientific	<< slipper.nodes[n].x << "\t"
									<< slipper.nodes[n].y << "\t"
									<< "0.0";
			f << endl;
		}

		f << "CELLS " << slipper.facecnt << " " << slipper.facecnt * 4 << endl;
		for(int c=0; c<slipper.facecnt; c++)
		{
			f << "3\t";
			for(int i=0; i<3; i++)
			{
				f << slipper.faces[c].nodes[i] << "\t";
			}
			f << endl;
			
		}

		f << "CELL_TYPES " << slipper.facecnt << endl;
		for(int i=0; i<slipper.facecnt; i++)
		{
			f << "5" << endl;
		}

		f << "CELL_DATA " << slipper.facecnt << endl;
		f << "SCALARS pressure float 1" << endl;
		f << "LOOKUP_TABLE press" << endl;
		for(int i=0; i<slipper.facecnt; i++)
		{
			f << scientific << slipper.facePressure[i] << endl;
		}

		 
		f << "POINT_DATA " << slipper.nodecnt << endl;
		f << "SCALARS ehd float 1" << endl;
		f << "LOOKUP_TABLE ehd" << endl;
		for(int i=0; i<slipper.nodecnt; i++)
		{
			f << scientific << slipper.nodeDeformation[i] << endl;
		}
		


		f.close();
}
void CSlipperGap::SlipperInitEHD(void)
{
		if(gapinput->options_slipper.general.SlipperPressureDeformation == 1)
		{
			//Init the solid deform using init pressure field
			SlipperCalcEHD(1.0);
		}
}
void CSlipperGap::SlipperCalcEHD(const double alpha)
{
	Fluid2Slipper();
	slipper.calcDeform(operatingslippergap.pG, operatingslippergap.pcase, forcesslippergap.psocket);
	Slipper2Fluid(alpha);
	/*
	for(int k=0; k<2; k++)
	{
		//Smooth the ehd
		for(int j=Fluid.N-1; j>=0; j--)
		{
			for(int i=0; i<Fluid.M; i++)
			{
				double div = 4.0;
				double avg = 2*slipper.ehd(i,j);
				
				if(i < Fluid.M-2)
				{
					avg += slipper.ehd(i+1,j);
					div += 1.0;
				}

				if(i > 0)
				{
					avg += slipper.ehd(i-1,j);
					div += 1.0;
				} 

				avg += slipper.ehd(i,(Fluid.N+j+1)%Fluid.N);
				avg += slipper.ehd(i,(Fluid.N+j-1)%Fluid.N);

				slipper.ehd(i,j) = avg/div;
			}

			for(int i=Fluid.M-1; i>=0; i--)
			{
				double div = 4.0;
				double avg = 2*slipper.ehd(i,j);
				
				if(i < Fluid.M-2)
				{
					avg += slipper.ehd(i+1,j);
					div += 1.0;
				}

				if(i > 0)
				{
					avg += slipper.ehd(i-1,j);
					div += 1.0;
				} 

				avg += slipper.ehd(i,(Fluid.N+j+1)%Fluid.N);
				avg += slipper.ehd(i,(Fluid.N+j-1)%Fluid.N);

				slipper.ehd(i,j) = avg/div;
			}
		}
	}
	*/
	SlipperCalch();

}
void CSlipperGap::SwashplateInitEHD(const double alpha)
{
	//Update the FullGroup pressure for this degree
	int steps = gapinput->common.rev_steps;
	int step = int(floor((gapinput->common.phi_rev_deg+gapinput->common.phi_deg_tol) / gapinput->common.phistep_deg)) % steps;
	FullGroup.p(all,all,step) = Fluid.p;

	Fluid2Swashplate();
	swashplate.calcDeform(0, operatingslippergap.pcase, 0);
	Swashplate2Fluid(alpha);
	/*
	for(int k=0; k<2; k++)
	{
		//Smooth the ehd
		for(int j=Fluid.N-1; j>=0; j--)
		{
			for(int i=0; i<Fluid.M; i++)
			{
				double div = 4.0;
				double avg = 2*swashplate.ehd(i,j);
				
				if(i < Fluid.M-2)
				{
					avg += swashplate.ehd(i+1,j);
					div += 1.0;
				}

				if(i > 0)
				{
					avg += swashplate.ehd(i-1,j);
					div += 1.0;
				} 

				avg += swashplate.ehd(i,(Fluid.N+j+1)%Fluid.N);
				avg += swashplate.ehd(i,(Fluid.N+j-1)%Fluid.N);

				swashplate.ehd(i,j) = avg/div;
			}

			for(int i=Fluid.M-1; i>=0; i--)
			{
				double div = 4.0;
				double avg = 2*swashplate.ehd(i,j);
				
				if(i < Fluid.M-2)
				{
					avg += swashplate.ehd(i+1,j);
					div += 1.0;
				}

				if(i > 0)
				{
					avg += swashplate.ehd(i-1,j);
					div += 1.0;
				} 

				avg += swashplate.ehd(i,(Fluid.N+j+1)%Fluid.N);
				avg += swashplate.ehd(i,(Fluid.N+j-1)%Fluid.N);

				swashplate.ehd(i,j) = avg/div;
			}
		}
	}
	*/
	SlipperCalch();
}
void CSlipperGap::SwashplateCalcEHD(const double alpha)
{
	//Update the FullGroup pressure for this degree
	int steps = gapinput->common.rev_steps;
	int step = int(floor((gapinput->common.phi_rev_deg+gapinput->common.phi_deg_tol) / gapinput->common.phistep_deg)) % steps;
	FullGroup.p(all,all,step) = Fluid.p;

	//This should be much quicker than calling Fluid2Swashplate and swashplate.calcDeform even without OpenMP
	for(int i=0; i<swashplate.quickF2S.size(); i++)
	{
		const int m = swashplate.quickF2S[i][0];
		const int n = swashplate.quickF2S[i][1];
		const int f = swashplate.quickF2S[i][2];

		//Update the swashplate nodel deformation
		const double scale = (Fluid.p(m,n) - swashplate.facePressure[f])/100e+5;
		for(int n=0; n<swashplate.nodecnt; n++)
		{
			swashplate.nodeDeformation[n] += swashplate.IM[f][n]*scale;	
		}

		//Update the face pressure
		swashplate.facePressure[f] = Fluid.p(m,n);
	}

	Swashplate2Fluid(alpha);
	/*
	for(int k=0; k<2; k++)
	{
		//Smooth the ehd
		for(int j=Fluid.N-1; j>=0; j--)
		{
			for(int i=0; i<Fluid.M; i++)
			{
				double div = 4.0;
				double avg = 2*swashplate.ehd(i,j);
				
				if(i < Fluid.M-2)
				{
					avg += swashplate.ehd(i+1,j);
					div += 1.0;
				}

				if(i > 0)
				{
					avg += swashplate.ehd(i-1,j);
					div += 1.0;
				} 

				avg += swashplate.ehd(i,(Fluid.N+j+1)%Fluid.N);
				avg += swashplate.ehd(i,(Fluid.N+j-1)%Fluid.N);

				swashplate.ehd(i,j) = avg/div;
			}

			for(int i=Fluid.M-1; i>=0; i--)
			{
				double div = 4.0;
				double avg = 2*swashplate.ehd(i,j);
				
				if(i < Fluid.M-2)
				{
					avg += swashplate.ehd(i+1,j);
					div += 1.0;
				}

				if(i > 0)
				{
					avg += swashplate.ehd(i-1,j);
					div += 1.0;
				} 

				avg += swashplate.ehd(i,(Fluid.N+j+1)%Fluid.N);
				avg += swashplate.ehd(i,(Fluid.N+j-1)%Fluid.N);

				swashplate.ehd(i,j) = avg/div;
			}
		}
	}
	*/
	SlipperCalch();
}
void CSlipperGap::SwashplateInitPressure(vector<double> & pDC)
{
	int steps = gapinput->common.rev_steps;
	
	//Resize the arrays
	FullGroup.Gx.resize(Fluid.M,Fluid.N,steps);
	FullGroup.Gy.resize(Fluid.M,Fluid.N,steps);
	FullGroup.p.resize(Fluid.M,Fluid.N,steps);

	//Populate the arrays
	for(int i=0; i<steps; i++)
	{
		operatingslippergap.phi_deg = i*gapinput->common.phistep_deg;
		operatingslippergap.phi_rad = operatingslippergap.phi_deg*PI/180.0; //[rad] between 0 and 2PI

		//Initialize pG to pDC
		operatingslippergap.pG = pDC[i];
		
		//Update the coordinate systems
		SlipperSetCoordSys();

		//Define the pressure boundaries (note that the KD tree is required for theta bounds)	
		SlipperDefinePressureBounds();

		//Pressure bounds need to be defined first
		//This is needed to initialize temperature (and thus visco) and deformation below
		SlipperPressureBounds();

		FullGroup.Gx(all,all,i) = Fluid.Gx;
		FullGroup.Gy(all,all,i) = Fluid.Gy;
		
		//if the simulation is being resumed, this will be overwritten. must do it here anyways because
		//gapinput->common.resumed is not set yet!
		FullGroup.p(all,all,i) = Fluid.p;
	}


}
