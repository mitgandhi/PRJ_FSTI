# include "./mesh.h"
# include <string>
# include <fstream>
# include <sstream>
# include <set>
# include <algorithm>
# include "../FSTI_Block_dll/log.h"

using namespace std;
extern class gaplog Log;

# define VTK_TETRA 10
# define VTK_TRI 5

// ---------------------------------------------------------------------- //
void tetra_mesh::read_abaqus(const char* filename, double scale_factor)
{
	
	string line;

	string inputfile;

	ifstream infile(filename);
	if (!infile.is_open()) 	
	{
		Log << "\nMesh file " << filename << " not found!" << gaplog::endl; 
		exit(1);
	}

	Log << gaplog::endl << "  * Reading mesh " << filename  << " in Abaqus format" 
			 << gaplog::endl << gaplog::endl;

	map<int,point> _nodes;
	
	while(getline(infile,line))
	{

		// reading list of nodes
		if(line.compare("*NODE") == 0) 
		{
			Log << "     - Reading nodes ... ";
			while(true) 
			{
				// store the pointer location before reading
				int pos = infile.tellg(); 
				getline(infile,line);
				if(line.find("*") != string::npos)
				{
					infile.seekg(pos); // restore the pointer location
					break;
				}
				else 
				{
					istringstream iss(line.c_str());
					string substr;
					getline(iss,substr,',');
					int nidx;
					istringstream(substr) >> nidx;
					
					double x,y,z;
					getline(iss,substr,',');
					istringstream(substr) >> x;
					getline(iss,substr,',');
					istringstream(substr) >> y;
					getline(iss,substr,',');
					istringstream(substr) >> z;

					// scale the nodes
					x *= scale_factor;
					y *= scale_factor;
					z *= scale_factor;
					
					// add the node
					_nodes[nidx] = point(x,y,z);
				}				
			}
			
			Log << "done! (Found " << _nodes.size() << " nodes)" << gaplog::endl;

			// reorder the nodes
			int c = 0;
			nodes.resize(_nodes.size());
			for(map<int,point>::const_iterator it = _nodes.begin(); it!=_nodes.end(); it++)
			{
				o2n[it->first] = c;
				n2o[c] = it->first;
				nodes[c++] = it->second;
			}
		}

		// reading list of elements
		else if(line.find("*ELEMENT") != string::npos)
		{
			// new volume set
			if(line.find("TYPE=C3D4") != string::npos)
			{
				string setname;
				// extract the name
				string keyword("ELSET="); 
				size_t f = line.find(keyword);
				if(f != std::string::npos)
					setname = line.substr(f + keyword.size(),line.size() - 1);

				elm_sets[setname] = vector<int>(0);

				Log << "     - Reading element set " << setname << " ... ";

				while(true) 
				{
					int pos = infile.tellg(); // store the pointer location
					getline(infile,line);
					if (line.find("*") != string::npos)
					{
						infile.seekg(pos);
						break;
					}
					else
					{
						// create the new element
						tetra thiselement;
						
						istringstream iss(line.c_str());
						
						string substr;
						getline(iss,substr,',');
						int elmid;
						istringstream(substr) >> elmid;	// just throw it away
						
						for(unsigned int h=0; h<4; h++)
						{
							int nid;
							getline(iss,substr,',');
							istringstream(substr) >> nid;
							thiselement.nodes[h] = o2n[nid];
						}
						
						// assign the element set id containint this element
						thiselement.elm_set = setname;
						thiselement.id = elements.size();
						elements.push_back(thiselement);
						// this id is only referred to the volume elements, 
						// we don't use the id specified in the mesh, because 
						// those ids are defined counting also th_boundary faces 
						// elements" type S3
						elm_sets[setname].push_back(elements.size()-1);
					}
				}
				
				Log << "done! (Found " << elm_sets[setname].size() << " elements)." << gaplog::endl;
			}
			
		}
		// new node set
		else if(line.find("*NSET") != string::npos)
		{
			// extract the name
			string setname;
			string keyword("NSET="); 
			size_t f = line.find(keyword);
			if(f != std::string::npos)
				setname = line.substr(f + keyword.size(),line.size() - 1);
			//else
			//	Log << line << gaplog::endl;

			Log << "     - Reading node set " << setname << " ... ";

			node_sets[setname] = vector<int>(0);
			
			while(true) 
			{
				int pos = infile.tellg(); // store the pointer location
				getline(infile,line);
				if(line.find("*") != string::npos)
				{
					infile.seekg(pos);
					break;
				}
				else
				{
					istringstream iss(line.c_str());
					string substr;
					int cnt = 0;
					while(iss.good())
					{
						int nidx;
						getline(iss,substr,',');
						if(substr.size() > 0)
						{
							istringstream(substr) >> nidx;
							node_sets[setname].push_back(o2n[nidx]); 
							cnt++;
						}
					}
				}
			}

			Log << "done! (Found " << node_sets[setname].size() << " nodes)." << gaplog::endl;
		}
		
	}
		
	infile.close();

	Log << gaplog::endl;

}
// ---------------------------------------------------------------------- //
void tetra_mesh::check_free_nodes()
{
	Log << "  * Checking for free nodes ... ";

	map<int,bool> usednodes;
	for(int i=0; i<elements.size(); i++)
	{
		for(int j=0; j<4; j++)
			usednodes[elements[i].nodes[j]] = true;
	}

	int freenodes = nodes.size() - usednodes.size();

	if(freenodes > 0)
	{  
		Log << "\nIn the mesh there are " << freenodes << " free nodes!" 
				<< " Please fix this problem!" << gaplog::endl;
		exit(1);
	}

	Log << "ok: no free nodes found!" << gaplog::endl;
}
// ---------------------------------------------------------------------- //
void tetra_mesh::define_elements()
{
	for(unsigned int i=0; i<elements.size(); i++)
	{
		tetra& t = elements[i];
		
		t.id = i;

		const point V0 = nodes[t.nodes[0]];
		const point V1 = nodes[t.nodes[1]];
		const point V2 = nodes[t.nodes[2]];
		const point V3 = nodes[t.nodes[3]];

		point a =  V0 - V3;
		point b =  V1 - V3;
		point c =  V2 - V3;
		
		t.center = (V0 + V1 + V2 + V3)/4.0;
		t.volume =  fabs(a*(b^c))/6.0;	
		t.belm = false;

	}
}
// ---------------------------------------------------------------------- //
void tetra_mesh::define_faces()
{
	multimap<vector<int>, int> allfaces;
	map<vector<int>, int> _faces;
	
	int perm[4][3];
	perm[0][0] = 0, perm[0][1] = 1, perm[0][2] = 2;
	perm[1][0] = 0, perm[1][1] = 2, perm[1][2] = 3;
	perm[2][0] = 0, perm[2][1] = 3, perm[2][2] = 1;
	perm[3][0] = 1, perm[3][1] = 3, perm[3][2] = 2;

	Log << "  * Defining boundary faces ... ";

	// loop to all the elements
	for(int i=0, id=0, f=0; i<elements.size(); i++)
	{
		// loop to all the faces in the thermalElement
		for(int j=0; j<4; j++, id++)
		{
			vector<int> tmp(3);
			for(int h=0; h<3; h++)
			{
				// define face nodes
				tmp[h] = elements[i].nodes[perm[j][h]];
			}
			
			sort(tmp.begin(),tmp.end());
			
			allfaces.insert(pair<vector<int>, int>(tmp, id));
			
			if(_faces.find(tmp) == _faces.end())
			{
				_faces[tmp] = f++;
			}
		}

	}
	
	faces.resize(0);
		
	double A = 0;

	// fill the faces and b_faces containers
	int fid = 0;
	map<vector<int>, int>::const_iterator itm;
	for(itm=_faces.begin(); itm!=_faces.end(); itm++)
	{
		
		// define shared chells
		vector<int> sharedcells(2,-1);
		pair<multimap<vector<int>,int>::iterator,multimap<vector<int>,int>::iterator> rep;
		rep = allfaces.equal_range(itm->first);	// get repetitions

		// set shared cells index can be:
		//   one index -> th_boundary face
		//   two indices -> internal face
		int shf = 0;
		for (multimap<vector<int>,int>::iterator itmm=rep.first; itmm!=rep.second; ++itmm)
			sharedcells[shf++] = static_cast<int>(itmm->second)/4.0;

		// Neighbors of cell
		int c0 = sharedcells[0];
		int c1 = sharedcells[1];
		
		// Boundary cell
		if(c1 == -1)
		{
			// create a new face
			tri face;
			face.nodes[0] = itm->first[0];
			face.nodes[1] = itm->first[1];
			face.nodes[2] = itm->first[2];
			
			face.elm = c0;																	// label the associated element
			elements[c0].belm = true;												// label as th_boundary element
			elements[c0].ext_faces.push_back(faces.size());	// add the ext face id

			// determine the elm_fid
			map<int,int> elm_nid;
			for(int h=0; h<4; h++)
				elm_nid[elements[c0].nodes[h]] = h;

			int sum = elm_nid[face.nodes[0]] + elm_nid[face.nodes[1]] + elm_nid[face.nodes[2]];
			
			if(sum == 3)			face.elm_fid = 3;		// 0 + 1 + 2 -> 3
			else if(sum == 5) face.elm_fid = 1;		// 0 + 2 + 3 -> 1
			else if(sum == 6) face.elm_fid = 0;		// 1 + 2 + 3 -> 0
			else if(sum == 4) face.elm_fid = 2;		// 1 + 3 + 0 -> 2
				
			// define center
			face.center = (nodes[face.nodes[0]] + nodes[face.nodes[1]] + nodes[face.nodes[2]])/3.0;
			// define area and normal
			point r1 = nodes[face.nodes[1]] - nodes[face.nodes[0]];
			point r2 = nodes[face.nodes[2]] - nodes[face.nodes[0]];
			face.normal = r1^r2;
			face.area = 0.5*face.normal.mag();
			face.normal = face.normal/face.normal.mag();
			// check the orientation
			point r = face.center - elements[c0].center;
			r = r/r.mag();
			if(face.normal*r < 0)
				face.normal = -1.0*face.normal;
			face.id = fid++;

			A += face.area;

			// push back
			faces.push_back(face);


		};
	}

	Log << "done! Found " << faces.size() 
			 << " boundary faces. Total area = " << A << gaplog::endl;
}
// ---------------------------------------------------------------------- //
void tetra_mesh::ns2fs()
{
	Log << "  * Defining face sets ... " << gaplog::endl << gaplog::endl;

	int total_faces = 0;
	vector<bool> included(faces.size(), false);	// list of included faces
	
	for(geom_set::const_iterator it = node_sets.begin(); it != node_sets.end(); it++)
	{
		map<int,int> node_map;
		for(unsigned int i=0; i<it->second.size(); i++)
			node_map[it->second[i]] = 1;

		// new face set
		face_sets[it->first] = vector<int>(0);
		
		for(int f=0; f<faces.size(); f++)
		{
			int found = 0;
			for(int j=0; j<3; j++)
			{
				if(node_map[faces[f].nodes[j]] == 1)
					found++;
			}
			if(found == 3 && !included[f]) // check if the face was not already included
			{
				// push back this face
				face_sets[it->first].push_back(f);
				// update the face_sets vector
				faces[f].face_set = it->first;
				// update the list of included faces
				included[f] = true;
			}
		}

		Log << "    - Defined face set of size " << face_sets[it->first].size()
				 << " associated with the node set " << it->first
				 << gaplog::endl;
		//Log << "Checking... " <<face_sets[it->first]<< gaplog::endl; 
		//std::cout<<"a is of type: "<<typeid(it->first).name()<<std::endl; // Output 'a is of type int'
		//write_fset_vtk( it->first );
		//write_fset_vtk(it->first.c_str);
		total_faces += face_sets[it->first].size();
	}
	
	Log << "\n  Defined " << total_faces << " faces" << gaplog::endl;
	
	if(total_faces != faces.size())
	{
		face_sets["undefined"];
		for(unsigned int f=0; f<faces.size(); f++)
		{
			if(!included[f])
				face_sets["undefined"].push_back(f);
		}

		Log << "\n  Found " << faces.size() - total_faces << " undefined faces.\n";
		Log << "\n  Writing undefined faces set.\n";
		write_fset_vtk("undefined");
		

	}
}
// ---------------------------------------------------------------------- //
void tetra_mesh::build(const char* filename, double scale)
{
	read_abaqus(filename, scale); // mm to m
	check_free_nodes();
	define_elements();
	define_faces();
	ns2fs();
}
// ---------------------------------------------------------------------- //
void tetra_mesh::write_fset_vtk(const char* fset) const
{

	geom_set::const_iterator it = face_sets.find(fset);
	if(it == face_sets.end())
	{
		Log << "\nCould not find face set " << fset << "\n\n";
		exit(1);
	}
	const vector<int>& fs = it->second;

	set<int> _nodes;
	for(unsigned int i=0; i<fs.size(); i++)
	{
		for(int j=0; j<3; j++)
			_nodes.insert(faces[fs[i]].nodes[j]);
	}

	map<int,int> g2l;
	set<int>::const_iterator itt;
	int c=0;
	for(itt=_nodes.begin(); itt!=_nodes.end(); itt++)
		g2l[*itt] = c++;

	ostringstream oss;
	oss << fset << ".vtk";

	ofstream vtk(oss.str().c_str());

	// write VTK header
	vtk << 
		"# vtk DataFile Version 2.0" << endl <<
		"vtk output" << endl <<
		"ASCII" << endl <<
		"DATASET UNSTRUCTURED_GRID" << endl << 
		"POINTS " << _nodes.size() << " float" << endl;
	
	// ----------------------------- mesh nodes -------------------------------//
	
	for(itt=_nodes.begin(); itt!=_nodes.end(); itt++)
	{
		vtk << nodes[*itt].x() << "\t" 
				<< nodes[*itt].y() << "\t" 
				<< nodes[*itt].z() << endl; 
	}

	vtk << endl;

	// ------------------------ elements definition -------------------------- //

	vtk << "CELLS " << fs.size() << "\t" 
			<< (1 + 3)*fs.size() << endl;

	for(unsigned int i=0; i<fs.size(); i++) 
	{
		vtk << 3 << "\t";
		for(unsigned int j=0; j<3; j++)
			vtk << g2l[faces[fs[i]].nodes[j]] << "\t";
		vtk << endl;
	}
	
	vtk << endl;

	// --------------------------- elements type ----------------------------- //

	vtk << "CELL_TYPES " << fs.size() << endl;
			
	for(int i=0; i<fs.size(); i++) 
	{
		vtk << 5 << endl;
	}

}
// ---------------------------------------------------------------------- //
// ---------------------------------------------------------------------- //
void tetra_sub_mesh::build(std::vector<std::string> sets_list)
{

	// update the list of sets
	sets = sets_list;

	// define the list of elements
	for(unsigned int i=0; i<msh.elements.size(); i++)
	{
		for(unsigned int j=0; j<sets.size(); j++)
		{
			if(msh.elements[i].elm_set.compare(sets[j]) == 0)
			{
				geid.push_back(i);				
				break;
			}
		}
	}

	// define the local list of nodes
	
	set<int> nds;
	
	for(unsigned int i=0; i<geid.size(); i++)
	{
		for(int j=0; j<4; j++)
			nds.insert(msh.elements[geid[i]].nodes[j]);
	}

	for(set<int>::const_iterator it=nds.begin(); it!=nds.end(); it++)
	{
		g2l_nn[*it] = gnid.size();
		gnid.push_back(*it);
	}

}
// ---------------------------------------------------------------------- //
void tetra_sub_mesh::write_vtk(const char* fname)
{
	ofstream vtk(fname);

	// write VTK header
	vtk << 
		"# vtk DataFile Version 2.0" << endl <<
		"vtk output" << endl <<
		"ASCII" << endl <<
		"DATASET UNSTRUCTURED_GRID" << endl << 
		"POINTS " << nn() << " float" << endl;

	// ----------------------------- mesh nodes -------------------------------//
	
	for(unsigned int i=0; i<nn(); i++) 
	{
		vtk << get_node(i).x() << "\t" 
				<< get_node(i).y() << "\t" 
				<< get_node(i).z() << endl; 
	}

	vtk << endl;

	// ------------------------ elements definition -------------------------- //

	vtk << "CELLS " << ne() << "\t" << (1 + 4)*ne() << endl;

	for(int i=0; i<ne(); i++) 
	{

		vtk << 4 << "\t";
		for(unsigned int j=0; j<4; j++)
			vtk << g2l_nn[get_elm(i).nodes[j]] << "\t"; 
		vtk << endl;
	}
	
	vtk << endl;

	// --------------------------- elements type ----------------------------- //

	vtk << "CELL_TYPES " << ne() << endl;
			
	for(int i=0; i<ne(); i++) 
	{
		vtk << VTK_TETRA << endl;
	}

}
// ---------------------------------------------------------------------- //
