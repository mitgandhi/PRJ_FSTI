#include "fem.h"

namespace CasparSlipperFEM
{
	//Utility functions used for text reading
	vector<string> Tokenize(const string& str,const string& delimiters)
	{
		vector<string> tokens;

		 // Skip delimiters at beginning.
		 string::size_type lastPos = str.find_first_not_of(delimiters, 0);
		 // Find first "non-delimiter".
		 string::size_type pos = str.find_first_of(delimiters, lastPos);

		while (string::npos != pos || string::npos != lastPos)
		 {
			  // Found a token, add it to the vector.
			  tokens.push_back(str.substr(lastPos, pos - lastPos));
			  // Skip delimiters.  Note the "not_of"
			  lastPos = str.find_first_not_of(delimiters, pos);
			  // Find next "non-delimiter"
			  pos = str.find_first_of(delimiters, lastPos);
		 }

		return tokens;
	};

	double s2d(const string s)
	{
		istringstream iss (s);
		double n;
		iss >> n;
		return n;
	};
	int s2i(const string s)
	{
		istringstream iss (s);
		int n;
		iss >> n;
		return n;
	};

	//the primary fem input methods
			
	struct section //used in this type of option file reader
	{
		string name;
		vector<istringstream*> fields;
	};

	void fem::readoptions(string file)
	{
		//clear any current options
		options.materials.clear();
		options.thermal_boundries.clear();
		options.elastic_constraints.clear();

		ifstream in(file.c_str());
		if (!in.is_open()) 	
		{
			GapLog << "input::read(): Unable to open " << file << " input file!" << endl;
			exit(1);
		}

		//first load the materials from the materials.txt file
		//also create a map of mat name 2 material id
		map<string, size_t> matname2id;
		if(caspar_gap_input != NULL)
		{
			for(size_t i=0; i<caspar_gap_input->materials.size(); i++)
			{
				femoptions::material mat;

				mat.name = caspar_gap_input->materials[i].name;
				mat.alpha = caspar_gap_input->materials[i].alpha;
				mat.E = caspar_gap_input->materials[i].E;
				mat.lambda = caspar_gap_input->materials[i].lambda;
				mat.nu = caspar_gap_input->materials[i].nu;
				mat.rho = caspar_gap_input->materials[i].rho;
				
				options.materials.push_back(mat);

				matname2id[mat.name] = options.materials.size()-1;
			}
		}

		string tmp;
		string line;
		double val = 0;
		vector<section* > sections;

		//read the file
		while(getline(in,line)) 
		{
			if(line.size() > 0)
			{
				istringstream iss (line,istringstream::in);
				iss >> tmp;
				size_t comment  = tmp.find("//"); // check if there is a comment
				// if is not a comment read
				if (comment == string::npos) 
				{
					// new section
					if(tmp.compare("section") == 0)
					{
						section* s = new section;
						sections.push_back(s);
						iss >> tmp;	
						if(tmp.compare("section") != 0)
						{
							s -> name = tmp;
							
							while(tmp.compare("endsection") != 0)
							{
								getline(in,line);
								if(line.size() > 0)
								{
									tmp.resize(0);
									istringstream iss(line,istringstream::in);
									iss >> tmp;
									if(tmp.size() > 0)
									{
										s -> fields.push_back(new istringstream (line,istringstream::in));
									}
								}
							}
						}
						else
						{
							GapLog << "Missing section name!" << endl;
							exit(1);
						}
					}
				}
			}
		}

		//process the sections
		string fieldname;
		string fieldvalue;

		bool read_general = false;
		bool read_FEM_solver = false;

		for(unsigned int i=0; i<sections.size(); i++)
		{
			// read section general
			if((sections[i] -> name).compare("general") == 0)
			{
				read_general = true;
				int validfields = 0;
				for(unsigned int j=0; j<sections[i] -> fields.size(); j++)
				{
					*(sections[i] -> fields[j]) >> fieldname;
					if(fieldname.compare("meshFile") == 0)
					{
						*(sections[i] -> fields[j]) >> options.meshfile;
						validfields++;
					}
					if(fieldname.compare("IR") == 0)
					{
						*(sections[i] -> fields[j]) >> options.inrel;
						validfields++;
					}
					if(fieldname.compare("maxIters") == 0)
					{
						*(sections[i] -> fields[j]) >> options.maxIters; 
						validfields++;
					}
					if(fieldname.compare("tolerance") == 0)
					{
						*(sections[i] -> fields[j]) >> options.tolerance; 
						validfields++;
					}
				}
				if(validfields != 4)
				{	
					GapLog << "Some field in section general is missing."
							 << " Found " << validfields << " expected " << 4 << endl;
				}
			}
			
			// read new material
			if((sections[i] -> name).compare("materials") == 0)
			{		
				for(unsigned int j=0; j<sections[i] -> fields.size(); j++)
				{
					string s;
					getline(*(sections[i] -> fields[j]), s);
					vector<string> line = Tokenize(s, "\t ");

					if(line.size() == 2)
					{
						if(matname2id.count(line[1]) > 0)
						{
							options.materials[matname2id[line[1]]].sets.push_back(line[0]);
						}
					}
					
				}
			}

			// read new boundary
			if((sections[i] -> name).compare("boundary") == 0)
			{
				femoptions::thermal_boundary boundary;

				int validfields = 0;
				for(unsigned int j=0; j<sections[i] -> fields.size(); j++)
				{
					*(sections[i] -> fields[j]) >> fieldname;
					if(fieldname.compare("type") == 0)
					{
						validfields++;
						string bc_type;
						*(sections[i] -> fields[j]) >> bc_type;
						if(bc_type.compare("dirichlet") == 0)
						{
							boundary.type = boundary.DIRICHLET;
						} else if (bc_type.compare("neumann") == 0)
						{
							boundary.type = boundary.NEUMANN;
						} else if (bc_type.compare("mixed") == 0)
						{
							boundary.type = boundary.ROBBINS;
						} else {
							GapLog << "Wrong boundary definition: specify dirichlet, " 
									 << "neumann or mixed!" << endl;
							exit(1);
						}
					}
					else if(fieldname.compare("set") == 0)
						*(sections[i] -> fields[j]) >> boundary.set, validfields++;
					else if(fieldname.compare("Tp") == 0)
						*(sections[i] -> fields[j]) >> boundary.Tp, validfields++;
					else if(fieldname.compare("Tinf") == 0)
						*(sections[i] -> fields[j]) >> boundary.Tinf, validfields++;
					else if(fieldname.compare("h") == 0)
						*(sections[i] -> fields[j]) >> boundary.h, validfields++;
					else if(fieldname.compare("q") == 0)
						*(sections[i] -> fields[j]) >> boundary.q, validfields++;
					else if(fieldname.compare("//") != 0 && fieldname.compare("endsection") != 0)
					{
						GapLog << "Wrong keyword (" << fieldname << ") in " 
							 << sections[i] -> name << " definition!" << endl;
						exit(1);
					}
				}

				options.thermal_boundries.push_back(boundary);
			}

			// read new constraint
			if((sections[i] -> name).compare("constraint") == 0)
			{
				femoptions::elastic_constraint constraint;

				int validfields = 0;
				for(unsigned int j=0; j<sections[i] -> fields.size(); j++)
				{
					*(sections[i] -> fields[j]) >> fieldname;
					if(fieldname.compare("set") == 0)
					{	
						*(sections[i] -> fields[j]) >> constraint.set, validfields++;
					}
					else if(fieldname.compare("x") == 0)
					{
						*(sections[i] -> fields[j]) >> constraint.x, validfields++;
					}
					else if(fieldname.compare("y") == 0)
					{
						*(sections[i] -> fields[j]) >> constraint.y, validfields++;
					}
					else if(fieldname.compare("z") == 0)
					{
						*(sections[i] -> fields[j]) >> constraint.z, validfields++;
					}
					else if(fieldname.compare("//") != 0 && fieldname.compare("endsection") != 0)
					{
						GapLog << "Wrong keyword (" << fieldname << ") in " 
							 << sections[i] -> name << " definition!" << endl;
						exit(1);
					}
				}
				
				if(constraint.set.length() > 0)
				{
					//normal constraint
					options.elastic_constraints.push_back(constraint);
				} 
				
				//either empty section or inertia relief so do nothing further

			}

		}	//done processing the sections

		//check for mandatory sections

		if(read_general = false)
		{
			GapLog << endl << "In input file " << file
					 << " the general section ismissing!" << endl;
			exit(1);
		}
		if(read_FEM_solver = false)
		{
			GapLog << endl << "In input file " << file
					 << " the FEM_solver section ismissing!" << endl;
			exit(1);
		}

	};
	void fem::loadinp()
{
	string inpfile = options.meshfile;
	const double scale = options.scalefactor;

	//Unless node / element renumbering is implemented (not difficult) these must be cleared
	nodes.clear();
	nodecnt = 0;
	clearElements();
	facesets.clear();
	nodesets.clear();

	ifstream inp(inpfile.c_str());
	if(!inp.is_open())
	{
		GapLog << "Error! Unable to open: " << inpfile << endl;
		exit(1);
	}

	while(!inp.eof())
	{
		string l;
		getline(inp, l);;

		if(l.find("*NODE") == 0)
		{
			GapLog << "\tReading nodes... ";
			nodecnt = 0;
			while(true)
			{
				streamoff pos = inp.tellg();
				getline(inp, l);
				if(l.find("*") != string::npos)
				{
					inp.seekg(pos);
					break;	//finished with the node section
				}
				vector<string> line = Tokenize(l, ",");
				if( s2i(line[0])-1 != nodecnt)
				{
					GapLog << "Error! Node IDs are not numbered correctly. Problem at nid = " << line[0] << endl;
					exit(1);
				}

				nodes.push_back(node(scale*s2d(line[1]), scale*s2d(line[2]), scale*s2d(line[3])));		
				nodecnt++;
			}
			nodecnt = (int) nodes.size();

			//update the node id's
			for(int n=0; n<nodecnt; n++)
			{
				nodes[n].id = n;
			}

			GapLog << "read " << nodecnt << " nodes." << endl;

			//restart the main file loop
			continue;
		}

		if(l.find("*ELEMENT,TYPE=S3,ELSET=") == 0)
		{
			//get the face set name
			string name = l.substr(23);
			GapLog << "\tReading face set \"" << name << "\" ... ";

			while(true)
			{
				streamoff pos = inp.tellg();
				getline(inp, l);
				if(l.find("*") != string::npos)
				{
					inp.seekg(pos);
					break;	//finished with the node section
				}
				vector<string> line = Tokenize(l, ",");
				
				face Face;
				for(unsigned int i=1; i<line.size(); i++)
				{
					int nid = s2i(line[i])-1;
					if(nid < 0 || nid >= nodecnt)
					{
						GapLog << "Error! Problem with nid = " << nid << " of face = " << line[0] << " in faceset = " << name << endl;
						exit(1);
					}
					Face.push_back(&nodes[nid]);
				}

				facesets[name].push_back(Face);
			}
			GapLog << "read " << facesets[name].size() << " faces." << endl;
			
			//restart the main file loop
			continue;
		}

		if(l.find("*ELEMENT,TYPE=C3D4,ELSET=") == 0)
		{
			//get the element set name
			string name = l.substr(25);
			GapLog << "\tReading element set \"" << name << "\" ... ";
			
			int start_elecnt = elecnt;

			while(true)
			{
				streamoff pos = inp.tellg();
				getline(inp, l);
				if(l.find("*") != string::npos)
				{
					inp.seekg(pos);
					break;	//finished with the node section
				}
				vector<string> line = Tokenize(l, ",");
				if(line.size() != 5)
				{
					GapLog << "Error! Problem with element = " << line[0] << " in element set = " << name << endl;
					exit(1);
				}
				
				tetra_element * e = new tetra_element();

				for(int i=1; i<5; i++)
				{
					int nid = s2i(line[i])-1;
					if(nid < 0 || nid >= nodecnt)
					{
						GapLog << "Error! Problem with nid = " << nid << " in element = " << line[0] << " of element set = " << name << endl;
						exit(1);
					}
					e->nodes[i-1] = &nodes[nid];
				}

				e->assign_analysis(analysis_type);

				elements.push_back(e);
				elementsets[name].push_back(elecnt);
				elecnt++;
			}
			elecnt = (int) elements.size();

			GapLog << "read " << elecnt-start_elecnt << " elements." << endl;

			//restart the main file loop
			continue;
		}

		if(l.find("*NSET, NSET=") == 0)
		{
			//get the face set name
			string name = l.substr(12);
			GapLog << "\tReading node set \"" << name << "\" ... ";

			while(true)
			{
				streamoff pos = inp.tellg();
				getline(inp, l);
				if(l.find("*") != string::npos)
				{
					inp.seekg(pos);
					break;	//finished with the node section
				}
				vector<string> line = Tokenize(l, ",");
				
				for(unsigned int i=0; i<line.size(); i++)
				{
					int nid = s2i(line[i])-1;
					if(nid < 0 || nid >= nodecnt)
					{
						GapLog << "Error! Problem with nid = " << nid << " in node set = " << name << endl;
						exit(1);
					}
					nodesets[name].push_back(nid);
				}

			}
			GapLog << "read " << nodesets[name].size() << " nodes." << endl;
			
			//restart the main file loop
			continue;
		}

	}

	inp.close();

	GapLog << "\tBuilding element faces ... ";
	
	//update the elementfaces 
	for(int e=0; e<elecnt; e++)
	{
		vector<face> fs = elements[e]->getfaces();
		for(int f=0; f<fs.size(); f++)
		{
			elementfaces[fs[f]].push_back(e);
		}
	}

	//create a local vector of all EXTERNAL faces in the mesh
	vector<face> meshfaces;
	for(map<face, vector<int> >::iterator f = elementfaces.begin(); f != elementfaces.end(); f++)
	{
		if(f->second.size() == 1)
		{
			meshfaces.push_back(f->first);
		}
	}

	GapLog << "found " << meshfaces.size() << " external and " << elementfaces.size() << " total faces." << endl;

	//Build facesets from each node set
	for(map<string, nodeset>::iterator it=nodesets.begin(); it != nodesets.end(); it++)
	{
		//it->first - set name
		//it->second - the node set
		
		string name = it->first;

		//check so we don't clobber an existing faceset
		if(facesets.count(name) != 0)
		{
			error("face set = " + name + " is named the same as an existing nodeset and must be renamed!");
		}
		
		//create a c++ set that will allow for a fast lookup
		set<int> ns;
		for(int n=0; n<it->second.size(); n++)
		{
			ns.insert(it->second[n]);
		}

		//now we need to loop through all of the meshfaces and test for each face if
		//all the nodes are present in ns
		for(vector<face>::iterator f=meshfaces.begin(); f != meshfaces.end(); f++)
		{
			bool faceINnset = true;	//assume all the nodes are present

			for(int n=0; n<f->nodes.size(); n++)
			{
				if(ns.count(f->nodes[n]->id) == 0)
				{
					faceINnset = false;	//this node isn't present so set to false
					break;	//break the inner for loop
				}
			}

			if(faceINnset)
			{
				//add the face to the faceset
				facesets[name].push_back(*f);
			}
		}

		GapLog << "\tCreated face set \"" << name << "\" from the node set with " << facesets[name].size() << " faces." << endl;
	}

	//the final step is to remove any newly created facesets that might 'clobber' faces in the gap faceset

	//first just check that gap exists
	if(!checkFaceset("gap"))
	{
		error("a node or face set named gap MUST exist!");
	}

	//let's create a c++ set for the gap faceset for quick lookup
	set<face> gapfaceset;
	for(int f=0; f<facesets["gap"].size(); f++)
	{
		gapfaceset.insert(facesets["gap"][f]);
	}

	//now loop through all NODE created facesets (except "gap") and remove any duplicate faces
	//if a faceset is actually defined in the input, we will not correct it
	for(map<string, nodeset>::iterator it=nodesets.begin(); it != nodesets.end(); it++)
	{
		string name = it->first;

		//obviously we need to skip the gap faceset
		if(name.compare("gap") == 0)
		{
			continue;
		}
		
		//let's just keep a count of how many faces we delete for logging purposes
		int del = 0;

		for(int f = 0; f<facesets[name].size(); f++)
		{
			if(gapfaceset.count(facesets[name][f]) > 0)
			{
				facesets[name].erase(facesets[name].begin() + f);
				f--;	//becuase we removed the current element, update the counter
				del++;
			}
		}

		if(del > 0)
		{
			GapLog << "\tRemoving " << del << " face(s) from the " << name << " face set because they duplicate the gap face set." << endl;
		}
	}
};

	void fem::writeK(const string filename)
	{
		//This will write the K matrix in ijv format to a text file
		ofstream kf(filename.c_str());
		if(!kf.is_open())
		{
			return;
		}

		for(int j=0; j<K.ncols(); j++)
		{
			gmm::wsvector_iterator<double> it;
			for(it = K[j].begin(); it != K[j].end(); it++)
			{
				kf << it.index() << "\t" << j << "\t" << (*it) << endl;
			}
		}
		
		kf.close();
	}

	void fem::writeB(const string filename)
	{
		//This will write the b vector in dense format to a text file
		ofstream bf(filename.c_str());
		if(!bf.is_open())
		{
			return;
		}

		for(int i=0; i<b.size(); i++)
		{
			bf << b[i] << endl;
		}
		
		bf.close();
	}

	void fem::writeVTK(const string filename, const vector<double> & nodetmp, const vector<double> & dirchlet,const vector<double> & neumann,const vector<double> & mixed_h,const vector<double> & mixed_phi)
{
	//define the special binary vtk writing class
	//this vtkwriter class is special because it will write double / float values in binary **Big-Endian** format
	class vtkwriter : public ofstream
	{
	private:
		union fic	//this union will make it easy to reverse the float bytes
		{
			float f;
			int i;
			char b[4];
		};
		inline void bwrite(const fic &littleE)
		{
			fic bigE;
			bigE.b[0] = littleE.b[3];
			bigE.b[1] = littleE.b[2];
			bigE.b[2] = littleE.b[1];
			bigE.b[3] = littleE.b[0];
	
			write(bigE.b, sizeof(fic));
		}

	public:
		//overload the << operator for double, float, or int
		vtkwriter& operator<<(const double& d)
		{
			fic v;
			v.f = static_cast<float>(d);
			bwrite(v);
			return *this;
		}
		vtkwriter& operator<<(const float& f)
		{
			fic v;
			v.f = f;
			bwrite(v);
			return *this;
		}
		vtkwriter& operator<<(const int& i)
		{
			fic v;
			v.i = i;
			bwrite(v);
			return *this;
		}
	};

	# define VTK_VOL_ELM_TYPE 10
	# define VTK_SURF_ELM_TYPE 5

	GapLog << "Writing VTK output ... ";

	//ofstream vtk;
	vtkwriter vtk;
	vtk.open(filename.c_str(), ios::binary);

	if (!vtk.is_open()) 
	{
		GapLog << "Error opening " << filename << "\n";
		system("pause");
		exit(1);
	}

	// write VTK header
	vtk << 
		"# vtk DataFile Version 3.0" << "\n" <<
		"vtk output" << "\n" <<
		"BINARY" << "\n" <<
		"DATASET UNSTRUCTURED_GRID" << "\n" << 
		"POINTS " << nodecnt << " float" << "\n";
		
	// ----------------------------- mesh nodes -------------------------------//
	
	for(unsigned int i = 0; i < nodes.size(); i++) 
	{
		vtk << nodes[i][0];
		vtk << nodes[i][1];
		vtk << nodes[i][2]; 
	}

	// ------------------------ elements definition -------------------------- //
		
	//NOTE I HAVE FORCED ELEMENTE NODE CNT = 4 for now. -> need to be more general!
	const int elenodecnt = 4;
	vtk << "\nCELLS " << elecnt << "\t" << (1 + elenodecnt)*elecnt << "\n";

	for(int i = 0; i < elecnt; i++) {
		vtk << elenodecnt;
		for(unsigned int j = 0; j < elements[i]->nodes.size(); j++)
			vtk << elements[i]->nodes[j]->id; // index start from 0
	}
		
	// --------------------------- elements type ----------------------------- //

	vtk << "\nCELL_TYPES " << elecnt << "\n";
				
	for (int i = 0; i < elecnt; i++) 
	{
		vtk << int(VTK_VOL_ELM_TYPE);
	}

	vtk << "\nPOINT_DATA " << nodecnt << "\n";

	if(analysis_type == ELASTIC)
	{
		vtk << "\nVECTORS displacement float" << "\n";			
		for(int i=0; i<nodecnt; i++)
		{
			if(nodes[i].DOF.size() < 3)
			{
				error("Invalid nodal DOF for elastic analysis!");
			}

			for(int j=0; j<3; j++)
			{
				const int Gdof = nodes[i].DOF[j];
				if(Gdof >= 0)
				{
					vtk << X[Gdof];
				} else {
					vtk << ebcs[niddof(i,j)];
				}
			}
		}
		
		//Applied Loads
		vtk << "\nVECTORS thermal_load float" << "\n";
		for(int i=0; i<nodecnt; i++)
		{
			if(nodes[i].DOF.size() < 3)
			{
				error("Invalid nodal DOF for elastic analysis!");
			}

			for(int j=0; j<3; j++)
			{
				const int Gdof = nodes[i].DOF[j];
				if(Gdof >= 0)
				{
					vtk << b[Gdof];
				} else {
					vtk << ebcs[niddof(i,j)];
				}
			}
		}

		//temperature
		if(nodetmp.size() == nodes.size())
		{
			vtk << "\nSCALARS temperature float 1" << "\n";
			vtk << "\nLOOKUP_TABLE temperature" << "\n";
			for(int n=0; n<nodetmp.size(); n++)
			{
				vtk << nodetmp[n];
			}
		}

		//thermal boundries
		if(dirchlet.size() == nodecnt)
		{
			vtk << "\nSCALARS boundary_dirchlet float 1" << "\n";
			vtk << "\nLOOKUP_TABLE dirchlet" << "\n";
			for(int n=0; n<dirchlet.size(); n++)
			{
				vtk << dirchlet[n];
			}
		}

		if(neumann.size() == nodecnt)
		{
			vtk << "\nSCALARS boundary_neumann float 1" << "\n";
			vtk << "\nLOOKUP_TABLE neumann" << "\n";
			for(int n=0; n<neumann.size(); n++)
			{
				vtk << neumann[n];
			}
		}

		if(mixed_h.size() == nodecnt)
		{
			vtk << "\nSCALARS boundary_mixed_h float 1" << "\n";
			vtk << "\nLOOKUP_TABLE mixed_h" << "\n";
			for(int n=0; n<mixed_h.size(); n++)
			{
				vtk << mixed_h[n];
			}
		}

		if(mixed_phi.size() == nodecnt)
		{
			vtk << "\nSCALARS boundary_mixed_phi float 1" << "\n";
			vtk << "\nLOOKUP_TABLE mixed_phi" << "\n";
			for(int n=0; n<mixed_phi.size(); n++)
			{
				vtk << mixed_phi[n];
			}
		}

	}

	if(analysis_type == THERMAL)
	{
		vtk << "\nSCALARS temperature float 1" << "\n";
		vtk << "\nLOOKUP_TABLE temperature" << "\n";
		for(int i=0; i<nodecnt; i++)
		{
			if(nodes[i].DOF.size() != 1)
			{
				error("Invalid nodal DOF for thermal analysis!");
			}

			const int Gdof = nodes[i].DOF[0];
			if(Gdof >= 0)
			{
				vtk << X[Gdof];
			} else {
				vtk << ebcs[niddof(i,0)];
			}
		}
	}

	vtk.close();
	vtk.clear();

	GapLog << "done." << "\n";
	}

	void fem::writeFacesetVTK(const string filename, const string name)
	{
		if(!checkFaceset(name))
		{
			//Faceset doesn't exist so return
			return;
		}

		faceset * fs = &facesets[name];	//a pointer to our faceset

		int vtkelecnt = (int) fs->size();

		//we need to build a map of nodes in this faceset and assign a vtk node id
		map<node *, int> vtknodesmap;

		//we also need a vector of the nodes in this faceset ordered according to the local vtk node id
		vector<node *> vtknodes;
	
		for(faceset::iterator f=fs->begin(); f != fs->end(); f++)
		{
			for(vector<node*>::iterator n=f->nodes.begin(); n != f->nodes.end(); n++)
			{
				if(vtknodesmap.count((*n)) == 0)
				{
					//we need to insert the node
					vtknodesmap[*n] = (int) vtknodes.size();
					vtknodes.push_back(*n);
				}
			}
		}

		int vtknodecnt = (int) vtknodes.size();

		//now we can write the vtk file
		ofstream vtk(filename.c_str());
		if (!vtk.is_open()) 
		{
			error("Error opening " + filename + "!");
		}

		// write VTK header
		vtk << 
			"# vtk DataFile Version 2.0" << endl <<
			"vtk output" << endl <<
			"ASCII" << endl <<
			"DATASET UNSTRUCTURED_GRID" << endl << 
			"POINTS " << vtknodecnt << " double" << endl;
	
		// ----------------------------- mesh nodes -------------------------------//
	
	
		for(unsigned int i = 0; i < vtknodes.size(); i++) 
		{
			vtk << vtknodes[i]->x() << "\t" 
					<< vtknodes[i]->y() << "\t" 
					<< vtknodes[i]->z() << endl; 
		}
	
		vtk << endl;

		// ------------------------ elements definition -------------------------- //
	
		//NOTE I HAVE FORCED ELEMENTE NODE CNT = 3 for now. -> need to be more general!
		const int elenodecnt = 3;
		vtk << "CELLS " << vtkelecnt << "\t" 
			<< (1 + elenodecnt)*vtkelecnt << endl;

		for(faceset::iterator f=fs->begin(); f != fs->end(); f++)
		{
			vtk << elenodecnt << "\t";
			for(vector<node*>::iterator n=f->nodes.begin(); n != f->nodes.end(); n++)
			{
				vtk << vtknodesmap[*n] << "\t";
			}
			vtk << endl;
		}

		vtk << endl;

		// --------------------------- elements type ----------------------------- //

		vtk << "CELL_TYPES " << vtkelecnt << endl;
			
		for (int i = 0; i < vtkelecnt; i++) 
		{
			vtk << VTK_SURF_ELM_TYPE << endl;
		}

		vtk.close();

}
};
