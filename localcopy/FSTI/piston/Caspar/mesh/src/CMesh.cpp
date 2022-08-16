#include "CMesh.h"
#include "../../main/src/logger.h"
#include <iomanip>
#include "../../caspar_input/input.h"
#pragma once

extern class input myinput;

CMesh::CMesh()
{
	
};
CMesh::~CMesh()
{
	
	
};
void CMesh::ReadABAQUSTeth(string mesh_path)
{


	Log << "\n";
	Log << "Reading " + mesh_path + " Teth Mesh ..." << "\n";

	//Check to see if gap_block is defined... dwm


	//read boundaries and materials file
	ReadBoundaries();
	ReadMaterials();


	//size boundary vectors according to specified boundaries
 	if(body)
		qbABQS.resize(qbi+myinput.data.thermal.piston.neumann_bc.size());
	else
		qbABQS.resize(qbi+myinput.data.thermal.block.neumann_bc.size());
	phibABQS.resize((int) phibi.size());
	mixbABQS.resize((int) hbi.size());

	
	//open the input file
	ifstream fabqs;
	ifstream fconn;
	string line;


	//nodes coordinates
	vector<Array<double,1>> xyzABQS;
	//resize constraint arrays
	cxyz_th.resize(3);
	//initialize cells connectivity collectors
	int clr=0;
	vector<vector<Array<int,1>>> connABQS;
	connABQS.resize(myinput.data.materials.size());
	//number of nodes per face
	int nf=3;
	//initialize total boundary faces,nodes
	nf_bd_user=0;
	nn_bd_user=0;
	
	//read coordinates
	string input;
	input = mesh_path;// + "/" + mesh_name;
	fabqs.open(input.c_str());
	getline(fabqs,line);

	//check to see if file is open
	if(!fabqs.good()){
		Log << "Error when opening " << mesh_path << " mesh." << "\n";
		Log << "Perhaps the mesh was not found." << endl;
		exit(1);
	}

	//setup an array for each input file section to see if each input is used.
	Array<bool,1> materialcheck;
	Array<bool,1> mixbcheck;
	Array<bool,1> phibcheck;
	Array<bool,1> constraintcheck;
	Array<bool,1> qbcheck;
	int gapcheck = 0;
	if(body){
		materialcheck.resize(myinput.data.thermal.piston.materials.size());
		mixbcheck.resize(myinput.data.thermal.piston.mixed_bc.size());
		phibcheck.resize(myinput.data.thermal.piston.dirichlet_bc.size());
		constraintcheck.resize(myinput.data.thermal.piston.constraints.size());
		qbcheck.resize(myinput.data.thermal.piston.neumann_bc.size());
	}
	else{
		materialcheck.resize(myinput.data.thermal.block.materials.size());
		mixbcheck.resize(myinput.data.thermal.block.mixed_bc.size());
		phibcheck.resize(myinput.data.thermal.block.dirichlet_bc.size());
		constraintcheck.resize(myinput.data.thermal.block.constraints.size());
		qbcheck.resize(myinput.data.thermal.block.neumann_bc.size());
	}
	materialcheck = 1;
	mixbcheck = 1;
	phibcheck = 1;
	constraintcheck = 1;
	qbcheck = 1;

	//read loop
	Log << "Reading Mesh" << "\n";
	//Log << line << "\n";
	while(!fabqs.eof())
	{
		//----------------read node coordinates
		if(line=="*NODE")
		{
			Log << "Reading Node Coordinates" << "\n";
			int i=0;
			double c;
			while(!fabqs.eof())
			{
				getline(fabqs,line);
				size_t found = line.find(",");
				size_t comment = line.find("*");
				//break reading if comment
				if(comment!=string::npos){
					break;
				};
				//skip comma in the line
				while(found!=string::npos)
				{
					line[found]=' ';
					found=line.find(",",found+1);
				};
				istringstream iss(line,istringstream::in);
				//[index x y z ] array
				Array<double,1> nodes;
				//push back every node
				xyzABQS.push_back(nodes);
				xyzABQS[i].resize(4);
				xyzABQS[i]=0.0;
				//assign index and coordinates
				for(int j=0;j<4;j++)
				{
					iss >> c;
					xyzABQS[i](j) = c;
				};
				i++;
			};
		}
		//-------------------read node connectivity
		else if(line.find("TYPE=C3D4")!=string::npos)
		{
			//Find the name of the mesh
			std::string name = line.substr(line.find("ELSET=") + 6);
			Log << "Matching Mesh " << name << " to Material Properties Definition..." << "\n";
			Log << "\n";
			int clr; //mesh index
			bool found = 0;
			if(body){
				for(int i = 0;i<myinput.data.thermal.piston.materials.size();i++){
					if(name.compare(myinput.data.thermal.piston.materials[i].first) == 0){
						clr = i;
						found = 1;
						materialcheck(i) = 0;//this has been used.
					}
				}
			}
			else{
				for(int i = 0; i < myinput.data.thermal.block.materials.size();i++){
					if(name.compare(myinput.data.thermal.block.materials[i].first ) == 0){
						clr = i;
						found = 1;
						materialcheck(i) = 0;//this has been used.
					}
				}
			}
			if(!found)
			{
				Log << "Error: Material for element set \"" << name << "\" not defined in input file!" << endl;
				exit(1);
			}

			int temp = clr + 1;
			Log << "Success!  " << name << " matched with mesh " << temp << "\n";
			Log << "\n";
			//read mesh
			int i=0;
			int id;
			while(!fabqs.eof())
			{
				getline(fabqs,line);
				size_t found = line.find(",");
				size_t comment = line.find("*");
				//break reading if comment
				if(comment!=string::npos){
					break;
				};
				//skip comma in the line
				while (found!=string::npos)
				{
					line[found]=' ';
					found=line.find(",",found+1);
				};
				istringstream iss(line,istringstream::in);
				//[index n1 n2 n3 n4] array
				Array<int,1> ids;
				//Log << "clr: " << clr << " i: " << i << "\n";
				//push back every element in container
				//Log << "Pushing Back..." << "\n";
				connABQS[clr].push_back(ids);
				//Log << "Resizing..." << "\n";
				connABQS[clr][i].resize(5);
				//Log << "Setting to Zero..." << "\n";
				connABQS[clr][i]=0;
				//assign element connectivity
				for(int jj=0;jj<5;jj++)
				{
					iss >> id;
					connABQS[clr][i](jj) = id-1;
				};
				i++;
			};
		}
		//-----------------read thermal node sets (xyz direction fem constraints and boundary conditions)
		else if(line.find("*NSET")!=string::npos)//check if nodeset
		{
			Log << "Matching Node Set... " << "\n";
			Log << "\n";
			std::string name = line.substr(line.find("NSET=") + 5);
			bool cons = 0;
			int cons_k; //index of set name in list of set names 
			bool phib = 0;
			int phib_k;
			bool mixb = 0;
			int mixb_k;
			bool qb = 0;
			int qb_k;
			bool found = 0;
			Array<bool,1> run;
			run.resize(3);
			run = 0;
			if(body){
				for(int i = 0;i<myinput.data.thermal.piston.constraints.size();i++){
					if(name.compare(myinput.data.thermal.piston.constraints[i].setname) == 0){
						cons = 1;
						cons_k = i;
						found = 1;
						constraintcheck(i) = 0;//this has been used.
						Log << "*********************************************" << "\n";
						Log << "Node Set: " << name << "\n";
						Log << "Type: Constraint" << "\n";
						Log << "Directions: ";
						if(myinput.data.thermal.piston.constraints[i].x){
							Log << "x, ";
							run(0) = 1;
						};
						if(myinput.data.thermal.piston.constraints[i].y){
							Log << "y, ";
							run(1) = 1;
						};
						if(myinput.data.thermal.piston.constraints[i].z){
							Log << "z";
							run(2) = 1;
						};
						Log << "\n";
						Log << "\n";
					};
				};
				for(int i = 0;i<myinput.data.thermal.piston.dirichlet_bc.size();i++){
					if(name.compare(myinput.data.thermal.piston.dirichlet_bc[i].setname) == 0){
						phib = 1;
						phib_k = i;
						found = 1;
						phibcheck(i) = 0;//this has been used.
						Log << "*********************************************" << "\n";
						Log << "Node Set: " << name << "\n";
						Log << "Type: Dirichlet" << "\n";
						Log << "Tp: " << myinput.data.thermal.piston.dirichlet_bc[i].Tp << "\n";
						Log << "\n";
					};
				};
				for(int i = 0;i<myinput.data.thermal.piston.mixed_bc.size();i++){
					if(name.compare(myinput.data.thermal.piston.mixed_bc[i].setname) == 0){
						mixb = 1;
						mixb_k = i;
						found = 1;
						mixbcheck(i) = 0;//this has been used.
						Log << "*********************************************" << "\n";
						Log << "Node Set: " << name << "\n";
						Log << "Type: Mixed" << "\n";
						Log << "T_inf: " << myinput.data.thermal.piston.mixed_bc[i].Tinf << "\n";
						Log << "h: " << myinput.data.thermal.piston.mixed_bc[i].h << "\n";
						Log << "\n";
					};
				};
				for(int i = 0;i<myinput.data.thermal.piston.neumann_bc.size();i++){
					if(name.compare(myinput.data.thermal.piston.neumann_bc[i].setname) == 0){
						qb = 1;
						qb_k = i + qbi;
						found = 1;
						qbcheck(i) = 0;//this has been used.
						Log << "*********************************************" << "\n";
						Log << "Node Set: " << name << "\n";
						Log << "Type: Neumann" << "\n";
						Log << "q: " << myinput.data.thermal.piston.neumann_bc[i].q[0] << "\n";
						Log << "\n";
					};
				};
			}
			else{
				for(int i = 0; i < myinput.data.thermal.block.constraints.size();i++){
					if(name.compare(myinput.data.thermal.block.constraints[i].setname ) == 0){
						cons = 1;
						cons_k = i;
						found = 1;
						constraintcheck(i) = 0;//this has been used.
						Log << "*********************************************" << "\n";
						Log << "Node Set: " << name << "\n";
						Log << "Type: Constraint" << "\n";
						Log << "Directions: ";
						if(myinput.data.thermal.block.constraints[i].x){
							Log << "x, ";
							run(0) = 1;
						};
						if(myinput.data.thermal.block.constraints[i].y){
							Log << "y, ";
							run(1) = 1;
						};
						if(myinput.data.thermal.block.constraints[i].z){
							Log << "z";
							run(2) = 1;
						};
						Log << "\n";
						Log << "\n";
					};
				};
				for(int i = 0;i<myinput.data.thermal.block.dirichlet_bc.size();i++){
					if(name.compare(myinput.data.thermal.block.dirichlet_bc[i].setname) == 0){
						phib = 1;
						phib_k = i;
						found = 1;
						phibcheck(i) = 0;//this has been used.
						Log << "*********************************************" << "\n";
						Log << "Node Set: " << name << "\n";
						Log << "Type: Dirichlet" << "\n";
						Log << "Tp: " << myinput.data.thermal.block.dirichlet_bc[i].Tp << "\n";
						Log << "\n";

					};
				};
				for(int i = 0;i<myinput.data.thermal.block.mixed_bc.size();i++){
					if(name.compare(myinput.data.thermal.block.mixed_bc[i].setname) == 0){
						mixb = 1;
						mixb_k = i;
						found = 1;
						mixbcheck(i) = 0;//this has been used.
						Log << "*********************************************" << "\n";
						Log << "Node Set: " << name << "\n";
						Log << "Type: Mixed" << "\n";
						Log << "T_inf: " << myinput.data.thermal.block.mixed_bc[i].Tinf << "\n";
						Log << "h: " << myinput.data.thermal.block.mixed_bc[i].h << "\n";
						Log << "\n";
					};
				};
				for(int i = 0;i<myinput.data.thermal.block.neumann_bc.size();i++){
					if(name.compare(myinput.data.thermal.block.neumann_bc[i].setname) == 0){
						qb = 1;
						qb_k = i + qbi;
						found = 1;
						qbcheck(i) = 0;//this has been used.
						Log << "*********************************************" << "\n";
						Log << "Node Set: " << name << "\n";
						Log << "Type: Neumann" << qb_k << "\n";
						Log << "q: " << myinput.data.thermal.block.neumann_bc[i].q[0] << "\n";
						Log << "\n";
					};
				};
			};
			if(!found && (name.find("gap_piston_")==string::npos))
			{
				Log << "Node set " << name << " not defined in input file." << "\n";
				getline(fabqs,line);
			};
			if(cons)
			{
				//int j;
				/*/choose degree of freedom 
				if(myinput.data.thermal.piston.constraints[k].x)
				{
					j=0;
				}
				if(line.find("cy_th")!=string::npos)
				{
					j=1;
				}
				if(line.find("cz_th")!=string::npos)
				{
					j=2;
				}*/
				//assign nodes
				int id;
				while(!fabqs.eof())
				{
					getline(fabqs,line);
					size_t found = line.find(",");
					size_t comment = line.find("*");
					//break if comment
					if(comment!=string::npos){
						break;
					};
					//skipe if comma in line
					while(found!=string::npos)
					{
						line[found]=' ';
						found=line.find(",",found+1);
					};
					

					for(int j = 0; j<3 ; j++){
						//Log << "Loop " << j << " Run = " << run(j) << "\n";
						if(run(j)){
							//push back nodes to current j degree of freedom
							istringstream iss(line,istringstream::in);
							while(iss)
							{
								iss >> id;
								cxyz_th(j).push_back(id-1);
								//int temp = id - 1 ;
								//Log << "Assigning Constraint " << j << " = " << temp << "\n";
							};
							//check last node is not read twice
							if((int) cxyz_th(j).size()>1 && cxyz_th(j)[(int) cxyz_th(j).size()-1]==cxyz_th(j)[(int) cxyz_th(j).size()-2])
							{
								cxyz_th(j).pop_back();
							};
						};
					};
				};
				
			}
			else if(phib){
				//assign nodes
				int m = phib_k;
				int id;
				while(!fabqs.eof())
				{
					getline(fabqs,line);
					size_t found = line.find(",");
					size_t comment = line.find("*");
					if(comment!=string::npos){
						break;
					};
					while(found!=string::npos)
					{
						line[found]=' ';
						found=line.find(",",found+1);
					};
					istringstream iss(line,istringstream::in);
					while(iss)
					{
						iss >> id;
						phibABQS[m].push_back(id-1);
						nn_bd_user++;
					};
					if((int) phibABQS[m].size()>1 && phibABQS[m][(int) phibABQS[m].size()-1]==phibABQS[m][(int) phibABQS[m].size()-2])
					{
						phibABQS[m].pop_back();
						nn_bd_user--;
					};
				};
			}

			else if(mixb)
			{
				int n = mixb_k;
				//assign nodes
				int id;
				while(!fabqs.eof())
				{
					getline(fabqs,line);
					size_t found = line.find(",");
					size_t comment = line.find("*");
					if(comment!=string::npos){
						break;
					};
					while(found!=string::npos)
					{
						line[found]=' ';
						found=line.find(",",found+1);
					};
					istringstream iss(line,istringstream::in);
					while(iss)
					{
						iss >> id;
						mixbABQS[n].push_back(id-1);
						nn_bd_user++;
					};
					if((int) mixbABQS[n].size()>1 && mixbABQS[n][(int) mixbABQS[n].size()-1]==mixbABQS[n][(int) mixbABQS[n].size()-2])
					{
						mixbABQS[n].pop_back();
						nn_bd_user--;
					};
				};
			}
			else if(qb)
			{
				int n = qb_k;
				//if(qbABQS.size()>=n){
					//qbABQS.resize(n);
				//	vector<int> temp;
				//	qbABQS.push_back(temp);
				//}; removed dwm
				//assign nodes
				int id;
				while(!fabqs.eof())
				{
					getline(fabqs,line);
					size_t found = line.find(",");
					size_t comment = line.find("*");
					if(comment!=string::npos){
						break;
					};
					while(found!=string::npos)
					{
						line[found]=' ';
						found=line.find(",",found+1);
					};
					istringstream iss(line,istringstream::in);
					while(iss)
					{
						iss >> id;
						qbABQS[n].push_back(id-1);
						nn_bd_user++;
					};
					if((int) qbABQS[n].size()>1 && qbABQS[n][(int) qbABQS[n].size()-1]==qbABQS[n][(int) qbABQS[n].size()-2])
					{
						qbABQS[n].pop_back();
						nn_bd_user--;
					};
				};
			}

			else if(name.find("gap_piston_")!=string::npos){//check if a qb
				Log << "*********************************************" << "\n";
				Log << "Node Set: " << name << "\n";
				Log << "Type: Gap" << "\n";
				Log << "\n";
				gapcheck += 1;
				//get boundary number
				size_t found = line.find_first_of("0123456789");
				string idx;
				while(found!=string::npos)
				{
					idx += line[found];
					found = line.find_first_of("0123456789",found+1);
				};
				//covert string to integer
				const char *idx_c;
				idx_c = idx.c_str();
				int l = atoi(idx_c)-1;
				//check that there aren't too many gap interfaces defined
				if(l >= qbi)
				{
					Log << "Error: Too many gap_piston interfaces defined in mesh!  Expected " << qbi << ", found " << l << "." << endl;
					exit(1);
				};
				//assign nodes
				int id;
				while(!fabqs.eof())
				{
					getline(fabqs,line);
					size_t found = line.find(",");
					size_t comment = line.find("*");
					if(comment!=string::npos){
						break;
					};
					while(found!=string::npos)
					{
						line[found]=' ';
						found=line.find(",",found+1);
					};
					istringstream iss(line,istringstream::in);
					while(iss)
					{
						iss >> id;
						qbABQS[l].push_back(id-1);
						nn_bd_user++;
					};
					if((int) qbABQS[l].size()>1 && qbABQS[l][(int) qbABQS[l].size()-1]==qbABQS[l][(int) qbABQS[l].size()-2])
					{
						qbABQS[l].pop_back();
						nn_bd_user--;
					};
				};
			};

			
			//else{
			//	getline(fabqs,line);
			//};

			Log << "Done!" << "\n";
			Log << "\n";
		}

		//----------------------read face sets
		else if(line.find("TYPE=S3")!=string::npos){//check if face set
			Log << "Matching Face Set..." << "\n";
			Log << "\n";
			std::string name = line.substr(line.find("ELSET=") + 6);
			bool phibb = 0;
			int phib_k;
			bool mixbb = 0;
			int mixb_k;
			bool qbb = 0;
			int qb_k;
			bool found = 0;
			if(body){
				for(int i = 0;i<myinput.data.thermal.piston.dirichlet_bc.size();i++){
					if(name.compare(myinput.data.thermal.piston.dirichlet_bc[i].setname) == 0){
						phibb = 1;
						phib_k = i;
						found = 1;
						phibcheck(i) = 0;
						Log << "*********************************************" << "\n";
						Log << "Face Set: " << name << "\n";
						Log << "Type: Dirichlet" << "\n";
						Log << "Tp: " << myinput.data.thermal.piston.dirichlet_bc[i].Tp << "\n";
						Log << "\n";
					};
				};
				for(int i = 0;i<myinput.data.thermal.piston.mixed_bc.size();i++){
					if(name.compare(myinput.data.thermal.piston.mixed_bc[i].setname) == 0){
						mixbb = 1;
						mixb_k = i;
						found = 1;
						mixbcheck(i) = 0;
						Log << "*********************************************" << "\n";
						Log << "Face Set: " << name << "\n";
						Log << "Type: Mixed" << "\n";
						Log << "Tinf: " << myinput.data.thermal.piston.mixed_bc[i].Tinf << "\n";
						Log << "h: " << myinput.data.thermal.piston.mixed_bc[i].h << "\n";
						Log << "\n";
					};
				};
				for(int i = 0;i<myinput.data.thermal.piston.neumann_bc.size();i++){
					if(name.compare(myinput.data.thermal.piston.neumann_bc[i].setname) == 0){
						qbb = 1;
						qb_k = i + qbi;
						found = 1;
						qbcheck(i) = 0;
						Log << "*********************************************" << "\n";
						Log << "Face Set: " << name << "\n";
						Log << "Type: Neumann" << "\n";
						Log << "q: " << myinput.data.thermal.piston.neumann_bc[i].q[0] << "\n";
						Log << "\n";
					};
				};
			}
			else{
				for(int i = 0;i<myinput.data.thermal.block.dirichlet_bc.size();i++){
					if(name.compare(myinput.data.thermal.block.dirichlet_bc[i].setname) == 0){
						phibb = 1;
						phib_k = i;
						found = 1;
						phibcheck(i) = 0;
						Log << "*********************************************" << "\n";
						Log << "Face Set: " << name << "\n";
						Log << "Type: Dirichlet" << "\n";
						Log << "Tp: " << myinput.data.thermal.block.dirichlet_bc[i].Tp << "\n";
						Log << "\n";
					};
				};
				for(int i = 0;i<myinput.data.thermal.block.mixed_bc.size();i++){
					if(name.compare(myinput.data.thermal.block.mixed_bc[i].setname) == 0){
						mixbb = 1;
						mixb_k = i;
						found = 1;
						mixbcheck(i) = 0;
						Log << "*********************************************" << "\n";
						Log << "Face Set: " << name << "\n";
						Log << "Type: Mixed" << "\n";
						Log << "Tinf: " << myinput.data.thermal.block.mixed_bc[i].Tinf << "\n";
						Log << "h: " << myinput.data.thermal.block.mixed_bc[i].h << "\n";
						Log << "\n";
					};
				};
				for(int i = 0;i<myinput.data.thermal.block.neumann_bc.size();i++){
					if(name.compare(myinput.data.thermal.block.neumann_bc[i].setname) == 0){
						qbb = 1;
						qb_k = i + qbi;
						found = 1;
						qbcheck(i) = 0;
						Log << "*********************************************" << "\n";
						Log << "Face Set: " << name << "\n";
						Log << "Type: Neumann" << "\n";
						Log << "q: " << myinput.data.thermal.block.neumann_bc[i].q[0] << "\n";
						Log << "\n";
					};
				};
			};
			if(!found && (name.find("gap_piston_")==string::npos)){
				Log << "Face set " << name << " is not defined in input file." << "\n";
				getline(fabqs,line);
			};
			if(phibb){
				int m = phib_k;
				//assign nodes
				int id;
				vector<int> phib(nf+1);
				while(!fabqs.eof())
				{
					getline(fabqs,line);
					size_t found = line.find(",");
					size_t comment = line.find("*");
					//break if comment
					if(comment!=string::npos){
						break;
					};
					//skip if comma in line
					while(found!=string::npos)
					{
						line[found]=' ';
						found=line.find(",",found+1);
					};
					//assign [index n1 n2 n3] face array to temporary vector
					istringstream iss(line,istringstream::in);
					for(int jj=0;jj<nf+1;jj++)
					{
						iss >> id;
						phib[jj] = id-1;
						//new boundary face
						if(jj==0)
						{
							nf_bd_user++;
						};
					};
					//push back to global m phib boundary container
					for(int jj=0;jj<nf;jj++)
					{
						phibABQS[m].push_back(phib[jj+1]);
						//eliminate duplicate nodes
						for(int ii=0;ii<phibABQS[m].size()-1;ii++)
						{
							if(phibABQS[m][ii]==phib[jj+1])
							{
								phibABQS[m].pop_back();
								break;
							};
						};
					};
				};
			}
			else if(mixbb){
				int n = mixb_k;
				//assign nodes
				int id;
				vector<int> mixb(nf+1);
				while(!fabqs.eof())
				{
					getline(fabqs,line);
					size_t found = line.find(",");
					size_t comment = line.find("*");
					//break if comment
					if(comment!=string::npos){
						break;
					};
					//skip if comma in line
					while(found!=string::npos)
					{
						line[found]=' ';
						found=line.find(",",found+1);
					};
					//assign [index n1 n2 n3] face array to temporary vector
					istringstream iss(line,istringstream::in);
					for(int jj=0;jj<nf+1;jj++)
					{
						iss >> id;
						mixb[jj] = id-1;
						//new boundary face
						if(jj==0)
						{
							nf_bd_user++;
						};
					};
					//push back to global n mixb boundary container
					for(int jj=0;jj<nf;jj++)
					{
						mixbABQS[n].push_back(mixb[jj+1]);
						//eliminate duplicate nodes
						for(int ii=0;ii<mixbABQS[n].size()-1;ii++)
						{
							if(mixbABQS[n][ii]==mixb[jj+1])
							{
								mixbABQS[n].pop_back();
								break;
							};
						};
					};
				};
			}
			else if(qbb){
				int n = qb_k;
				//if(qbABQS.size()>=n)
				//	qbABQS.resize(n-1); removed dwm
				//assign nodes
				int id;
				vector<int> qb(nf+1);
				while(!fabqs.eof())
				{
					getline(fabqs,line);
					size_t found = line.find(",");
					size_t comment = line.find("*");
					//break if comment
					if(comment!=string::npos){
						break;
					};
					//skip if comma in line
					while(found!=string::npos)
					{
						line[found]=' ';
						found=line.find(",",found+1);
					};
					//assign [index n1 n2 n3] face array to temporary vector
					istringstream iss(line,istringstream::in);
					for(int jj=0;jj<nf+1;jj++)
					{
						iss >> id;
						qb[jj] = id-1;
						//new boundary face
						if(jj==0)
						{
							nf_bd_user++;
						};
					};
					//push back to global n mixb boundary container
					for(int jj=0;jj<nf;jj++)
					{
						qbABQS[n].push_back(qb[jj+1]);
						//eliminate duplicate nodes
						for(int ii=0;ii<qbABQS[n].size()-1;ii++)
						{
							if(qbABQS[n][ii]==qb[jj+1])
							{
								qbABQS[n].pop_back();
								break;
							};
						};
					};
				};
			}
			else if(line.find("gap_piston_")!=string::npos)//check if a qb
			{
				Log << "*********************************************" << "\n";
				Log << "Face Set: " << name << "\n";
				Log << "Type: Gap" << "\n";
				Log << "\n";
				gapcheck += 1;
				//get boundary number
				size_t found = line.find_first_of("0123456789");
				string idx;
				while(found!=string::npos)
				{
					idx += line[found];
					found = line.find_first_of("0123456789",found+1);
				};
				idx.erase(0,1);
				//covert string to integer
				const char *idx_c;
				idx_c = idx.c_str();
				int l = atoi(idx_c)-1;
				//check that there aren't too many gap interfaces.
				if(l >= qbi)
				{
					Log << "Error: Too many gap interfaces defined in mesh!  Expected " << qbi << ", found " << l << "." << endl;
					exit(1);
				};
				//assign nodes
				int id;
				vector<int> qb(nf+1);
				while(!fabqs.eof())
				{
					getline(fabqs,line);
					size_t found = line.find(",");
					size_t comment = line.find("*");
					//break if comment
					if(comment!=string::npos){
						break;
					};
					//skip if comma in line
					while(found!=string::npos)
					{
						line[found]=' ';
						found=line.find(",",found+1);
					};
					//assign [index n1 n2 n3] face array to temporary vector
					istringstream iss(line,istringstream::in);
					for(int jj=0;jj<nf+1;jj++)
					{
						iss >> id;
						qb[jj] = id-1;
						//new boundary face
						if(jj==0)
						{
							nf_bd_user++;
						};
					};
					//push back to global l qb boundary container
					for(int jj=0;jj<nf;jj++)
					{
						qbABQS[l].push_back(qb[jj+1]);
						//eliminate duplicate nodes
						for(int ii=0;ii<qbABQS[l].size()-1;ii++)
						{
							if(qbABQS[l][ii]==qb[jj+1])
							{
								qbABQS[l].pop_back();
								break;
							};
						};
					};
				};
				Log << "Done!" << "\n";
				Log << "\n";
			};

		//--------------------read rest
			
		}
		else
		{
			getline(fabqs,line);	
		};
	};
	fabqs.close();
	Log << "Done Reading Mesh" << "\n";


	//grid sizes
	nNodes = (int) xyzABQS.size();
	nCells = 0;
	for(int j=0;j<(int)connABQS.size();j++)
	{
		nCells += (int) connABQS[j].size();
	}


	//assign inputs to coordinate field
	Log << "Geometry ..." << " ";
	xyz.resize(nNodes,3);
	for(int i=0;i<nNodes;i++)
	{
		for(int j=0;j<3;j++)
		{
			xyz(i,j) = xyzABQS[i](j+1) * 1.0e-3;
		}
	};
	Log << "done!" << "\n";


	//connectivity
	conn.resize(nCells,4);

	//diffusivity and material properties
	Tf.resize(nCells);
	E.resize(nCells);
	v.resize(nCells);
	alpha.resize(nCells);
	rho.resize(nCells);


	//assign inputs to connectivity, diffusivity and materials fields
	int nC=0;
	int nCi=0;
	Log << "Connectivity and materials ..." << " ";
	for(int i=0;i<(int)connABQS.size();i++)
	{
		nCi=nC;
		for(int j=0;j<(int)connABQS[i].size();j++)
		{
			for(int k=0;k<4;k++)
			{
				conn(nCi+j,k) = connABQS[i][j](k+1);
			}
			//diffusivity and material properties
			Tf(nCi+j)=Tfi[i];
			E(nCi+j)=Ei[i];
			v(nCi+j)=vi[i];
			alpha(nCi+j)=alphai[i];
			rho(nCi+j)=rhoi[i];
			nC++;
		};
	}
	Log << "done!" << "\n";


	//clear mesh inputs
	Tfi.clear();		Ei.clear();		vi.clear();		alphai.clear();

	Log << "Writing FEM inputs ..." << " ";


	//Output geometry
	fout.open("./temp/piston/input/geometry.fem");
	for(int i = 0; i < nNodes ; i++)
	{
		fout << scientific << xyz(i,0) << "\t" << xyz(i,1) << "\t" << xyz(i,2) << "\n";
	}
	fout.close();
	fout.clear();


	//Output connectivity
	fout.open("./temp/piston/input/connectivity.fem");
	for(int i = 0; i < nCells ; i++)
	{
		fout << scientific << conn(i,0) << "\t" << conn(i,1) << "\t" << conn(i,2) << "\t" << conn(i,3) << "\n";
	}
	fout.close();
	fout.clear();


	//Output materials
	fout.open("./temp/piston/input/materials.fem");
	for(int i = 0; i < nCells ; i++)
	{
		fout << scientific << E(i) << "\t" << v(i) <<  "\n";
	}
	fout.close();
	fout.clear();

	Log << "Done!" << "\n";

	Log << "Checking Mesh..." << "\n";

	if(gapcheck != qbi){
		Log << "Error: Expected " << qbi << " gap boundaries.  Found " << gapcheck << " defined in mesh!" << endl;
		exit(1);
	};
	for(int i = 0;i<phibcheck.size();i++){
		if(phibcheck(i)){
			if(body)
				Log << "Error, dirichlet boundary set " << myinput.data.thermal.piston.dirichlet_bc[i].setname << " not found in mesh!" << "\n";
			else
				Log << "Error, dirichlet boundary set " << myinput.data.thermal.block.dirichlet_bc[i].setname << " not found in mesh!" << endl;
			exit(1);
		};
	};
	for(int i = 0;i<materialcheck.size();i++){
		if(materialcheck(i)){
			if(body)
				Log << "Error, material set " << myinput.data.thermal.piston.materials[i].first << " not found in mesh!" << "\n";
			else
				Log << "Error, material set " << myinput.data.thermal.block.materials[i].first << " not found in mesh!" << endl;
			exit(1);
		};
	};
	for(int i = 0;i<constraintcheck.size();i++){
		if(constraintcheck(i)){
			if(body)
				Log << "Error, constraint node set " << myinput.data.thermal.piston.constraints[i].setname << " not found in mesh!" << "\n";
			else
				Log << "Error, constraint node set " << myinput.data.thermal.block.constraints[i].setname << " not found in mesh!" << endl;
			exit(1);
		};
	};
	for(int i = 0;i<mixbcheck.size();i++){
		if(mixbcheck(i)){
			if(body)
				Log << "Error, mixed boundary set " << myinput.data.thermal.piston.mixed_bc[i].setname << " not found in mesh!" << "\n";
			else
				Log << "Error, mixed boundary set " << myinput.data.thermal.block.mixed_bc[i].setname << " not found in mesh!" << endl;
			exit(1);
		};
	};
	for(int i = 0;i<qbcheck.size();i++){
		if(qbcheck(i)){
			if(body)
				Log << "Error, neumann boundary set " << myinput.data.thermal.piston.neumann_bc[i].setname << " not found in mesh!" << "\n";
			else
				Log << "Error, neumann boundary set " << myinput.data.thermal.block.neumann_bc[i].setname << " not found in mesh!" << endl;
			exit(1);
		};
	};
	if(body && !myinput.data.thermal.piston.IR){
		for(int i = 0;i<3;i++){
			if(!cxyz_th(i).size()){
				Log << "Error, piston not fully constrained!" << endl;
				exit(1);
			};
		};
	}
	else if(!body && !myinput.data.thermal.block.IR){
		for(int i = 0;i<3;i++){
			if(!cxyz_th(i).size()){
				Log << "Error, block not fully constrained!" << i << endl;
				exit(i);
			};
		};
	};


	Log << "Teth Mesh Read Successfully!" << "\n";
	Log << "\n";



};
void CMesh::ReadMaterials(void)
{
	if(body == 1){//piston
		Ei.resize(myinput.data.thermal.piston.materials.size());
		vi.resize(myinput.data.thermal.piston.materials.size());
		alphai.resize(myinput.data.thermal.piston.materials.size());
		Tfi.resize(myinput.data.thermal.piston.materials.size());
		rhoi.resize(myinput.data.thermal.piston.materials.size());
		for( int i = 0 ; i < myinput.data.thermal.piston.materials.size() ; i++ ){
			Log << "Finding Material Properties for " << myinput.data.thermal.piston.materials[i].first << "." << "\n";
			Log << "\n";
			bool found = 0;
			for( int j = 0 ; j < myinput.data.materials.size() ; j++ ){
				//Log << "Comparing: " << myinput.data.materials[j].name << " to " << myinput.data.thermal.piston.materials[i].second << "." << "\n";
				//int test = myinput.data.thermal.piston.materials[i].second.compare(myinput.data.materials[j].name);
				//Log << test << "\n";
				if(myinput.data.thermal.piston.materials[i].second.compare(myinput.data.materials[j].name) == 0){
					found = 1;
					Ei[i] = myinput.data.materials[j].E;
					vi[i] = myinput.data.materials[j].nu;
					alphai[i] = myinput.data.materials[j].alpha;
					Tfi[i] = myinput.data.materials[j].lambda;
					rhoi[i] = myinput.data.materials[j].rho;
					Log << "\n";
					int temp = i+1;
					Log << "Material Properties for " << myinput.data.materials[j].name << ", mesh " << temp << "\n";
					Log << "E = " << myinput.data.materials[j].E << "\t[Pa]" << "\n";
					Log << "v = " << myinput.data.materials[j].nu << "\t[-]" << "\n";
					Log << "alpha = " << myinput.data.materials[j].alpha << "\t[-]" << "\n";
					Log << "Tf = " << myinput.data.materials[j].lambda << "\t[W/mK]" << "\n";
					Log << "rho = " << myinput.data.materials[j].rho << "\t[kg/m3]" << "\n";
					Log << "\n";
				};
			};
			if(!found){
				Log << "Error: Material properties for " << myinput.data.thermal.piston.materials[i].first << " not found!" << endl;
				exit(1);
			}
		};
	}
	else if(body == 0){//bushing
		Ei.resize(myinput.data.thermal.block.materials.size());
		vi.resize(myinput.data.thermal.block.materials.size());
		alphai.resize(myinput.data.thermal.block.materials.size());
		Tfi.resize(myinput.data.thermal.block.materials.size());
		rhoi.resize(myinput.data.thermal.block.materials.size());
		for( int i = 0 ; i < myinput.data.thermal.block.materials.size() ; i++ ){
			Log << "Finding Material Properties for " << myinput.data.thermal.block.materials[i].first << "." << "\n";
			bool found = 0;
			for( int j = 0 ; j < myinput.data.materials.size() ; j++ ){
				if(myinput.data.thermal.block.materials[i].second.compare(myinput.data.materials[j].name) == 0){
					found = 1;
					Ei[i] = myinput.data.materials[j].E;
					vi[i] = myinput.data.materials[j].nu;
					alphai[i] = myinput.data.materials[j].alpha;
					Tfi[i] = myinput.data.materials[j].lambda;
					rhoi[i] = myinput.data.materials[j].rho;
					Log << "\n";
					int temp = i + 1;
					Log << "Material Properties for " << myinput.data.materials[j].name << ", mesh " << temp << "\n";
					Log << "E = " << myinput.data.materials[j].E << "\t[Pa]" << "\n";
					Log << "v = " << myinput.data.materials[j].nu << "\t[-]" << "\n";
					Log << "alpha = " << myinput.data.materials[j].alpha << "\t[-]" << "\n";
					Log << "Tf = " << myinput.data.materials[j].lambda << "\t[W/mK]" << "\n";
					Log << "rho = " << myinput.data.materials[j].rho << "\t[kg/m3]" << "\n";
					Log << "\n";
				};
			};
			if(!found){
				Log << "Error: Material properties for " << myinput.data.thermal.block.materials[i].first << " not found!" << endl;
				exit(1);
			}
		};
	};


	/*for(int i = 0;i<myinput.data.materials.size();i++){
		namei.push_back(myinput.data.materials[i].name);
		Ei.push_back(myinput.data.materials[i].E);
		vi.push_back(myinput.data.materials[i].nu);
		alphai.push_back(myinput.data.materials[i].alpha);
		Tfi.push_back(myinput.data.materials[i].lambda);
		rhoi.push_back(myinput.data.materials[i].rho);
	};

	first try dwm*/

	// open the input file
	/*string input;
	input = mesh_path + "/materials.mtr";
	ifstream fmaterials(input.c_str());
	
	if (!fmaterials.is_open()) 	{
		Log << "Unable to open materials.mtr! - Generate materials.mtr and restart..." << "\n";
		system("PAUSE");
		exit(1);
	}

	string tmp;
	string line;
	size_t comment;

	double E,v,alpha,Tf,rho;
	while(fmaterials)
	{
		getline(fmaterials,line);
		istringstream outiss (line,istringstream::in);
		outiss >> tmp;
		if(tmp == "Body") {
			getline(fmaterials,line);
			istringstream iniss(line);
			while(iniss.str() != "}") {
				while(iniss) {
					iniss >> tmp;
					comment  = tmp.find("//"); // check if there is a comment
					if (comment!=string::npos) {
						break;
					}
					if(tmp == "Emod"){
						iniss >> E;
						Ei.push_back(E);
						break;
					}
					if(tmp == "v"){
						iniss >> v;
						vi.push_back(v);
						break;
					}
					if(tmp == "alpha"){
						iniss >> alpha;
						alphai.push_back(alpha);
						break;
					}
					if(tmp == "lambda"){
						iniss >> Tf;
						Tfi.push_back(Tf);
						break;
					}
					if(tmp == "rho"){
						iniss >> rho;
						rhoi.push_back(rho);
						break;
					}
				}
				iniss.clear();	// this is essential!!
				getline(fmaterials,line);
				iniss.str(line);
			}
		};
		
	};*/

};
void CMesh::ReadBoundaries(void)
{
	if(body == 1)//piston
		qbi = 1; // only one gap boundary
	else if(body == 0){//bushing
		if(myinput.data.options_piston.general.EHDTestRig || myinput.data.options_piston.general.TriboTestRig)//test rigs with only one bushing simulated
			qbi = 1;
		else //simulate all bushings
			qbi = myinput.data.operating_conditions.npistons;
	};
	
	//dirichlet boundaries
	phibi.resize(0);//test dwm
	if(body == 1 && myinput.data.thermal.piston.dirichlet_bc.size()){//piston
		//phibi.resize(myinput.data.thermal.piston.dirichlet_bc.size();
		for(int i = 0;i< myinput.data.thermal.piston.dirichlet_bc.size(); i++)
		{
			phibi.push_back(myinput.data.thermal.piston.dirichlet_bc[i].Tp);
		};
	}
	else if(body == 0 && myinput.data.thermal.block.dirichlet_bc.size()){//bushing
		for(int i = 0;i < myinput.data.thermal.block.dirichlet_bc.size(); i++)
		{
			phibi.push_back(myinput.data.thermal.block.dirichlet_bc[i].Tp);
		};
	};

	//mixed boundaries
	hbi.resize(0);//test dwm
	phibinfi.resize(0);//test dwm
	if(body == 1 && myinput.data.thermal.piston.mixed_bc.size()){//piston
		for(int i=0;i<myinput.data.thermal.piston.mixed_bc.size();i++){
			hbi.push_back(myinput.data.thermal.piston.mixed_bc[i].h);
			phibinfi.push_back(myinput.data.thermal.piston.mixed_bc[i].Tinf);
		};
	}
	else if(body == 0 && myinput.data.thermal.block.mixed_bc.size()){//bushing
		for(int i=0;i<myinput.data.thermal.block.mixed_bc.size();i++){
			hbi.push_back(myinput.data.thermal.block.mixed_bc[i].h);
			phibinfi.push_back(myinput.data.thermal.block.mixed_bc[i].Tinf);
		};
	};


};
/*void CMesh::ReadMeshMain(void)
{
	// open the input file
	string input;
	input = "./main.dat";
	ifstream fmain(input.c_str());
	
	if (!fmain.is_open()) 	{
		Log << "Unable to open main.dat! - Generate main.dat and restart..." << "\n";
		system("PAUSE");
		exit(1);
	}

	string tmp;
	string line;
	size_t comment;
	
	while(fmain)
	{
		getline(fmain,line);
		istringstream outiss (line,istringstream::in);
		outiss >> tmp;
		// piston thermal
		if(tmp == "PistonThermal") {
			getline(fmain,line);
			istringstream iniss(line);
			while(iniss.str() != "}") {
				while(iniss) {
					iniss >> tmp;
					comment  = tmp.find("//"); // check if there is a comment
					if (comment!=string::npos) {
						break;
					}
					if(tmp == "Path"){
						iniss >> piston_path_th;
						break;
					}
					if(tmp == "Body"){
						iniss >> piston_name;
						break;
					}
					if(tmp == "Type"){
						iniss >> piston_type;
						break;
					}
					if(tmp == "IR"){
						iniss >> IR_p;
						break;
					}
					if(tmp == "d_tol"){
						iniss >> d_tol_p;
						break;
					}
				}
				iniss.clear();	// this is essential!!
				getline(fmain,line);
				iniss.str(line);
			}
		};
		
		// cylinder thermal
		if(tmp == "CylinderThermal") {
			getline(fmain,line);
			istringstream iniss(line);
			while(iniss.str() != "}") {
				while(iniss) {
					iniss >> tmp;
					comment  = tmp.find("//"); // check if there is a comment
					if (comment!=string::npos) {
						break;
					}
					if(tmp == "Path"){
						iniss >> cylinder_path_th;
						break;
					}
					if(tmp == "Body"){
						iniss >> cylinder_name;
						break;
					}
					if(tmp == "Type"){
						iniss >> cylinder_type;
						break;
					}
					if(tmp == "IR"){
						iniss >> IR_c;
						break;
					}
					if(tmp == "d_tol"){
						iniss >> d_tol_c;
						break;
					}
				}
				iniss.clear();	// this is essential!!
				getline(fmain,line);
				iniss.str(line);
			}
		};
		// piston pressure
		if(tmp == "PistonPressure") {
			getline(fmain,line);
			istringstream iniss(line);
			while(iniss.str() != "}") {
				while(iniss) {
					iniss >> tmp;
					comment  = tmp.find("//"); // check if there is a comment
					if (comment!=string::npos) {
						break;
					}
					if(tmp == "PathIM"){
						iniss >> IM_piston_path;
						break;
					}
				}
				iniss.clear();	// this is essential!!
				getline(fmain,line);
				iniss.str(line);
			}
		};
		// cylinder thermal
		if(tmp == "CylinderPressure") {
			getline(fmain,line);
			istringstream iniss(line);
			while(iniss.str() != "}") {
				while(iniss) {
					iniss >> tmp;
					comment  = tmp.find("//"); // check if there is a comment
					if (comment!=string::npos) {
						break;
					}
					if(tmp == "PathIM"){
						iniss >> IM_cylinder_path;
						break;
					}
				}
				iniss.clear();	// this is essential!!
				getline(fmain,line);
				iniss.str(line);
			}
		};
	};

};*/



