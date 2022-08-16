#include "CasparSlipperFEM.h"

namespace CasparSlipperFEM
{

	void thermoelastic::run_analysis(const string heatfluxsetname, const double alpha, const vector<double> & old_temp)
	{
		fem body;
		
		//set the gapinput object
		body.caspar_gap_input = caspar_gap_input;

		//load the options file
		GapLog << "Loading the options files ... ";
		body.readoptions(option_file);
		GapLog << "done." << endl << endl;

		GapLog << "Loading \"" << body.options.meshfile << "\"... " << endl;
		body.loadinp();
		GapLog << endl;

		//assing mats
		body.assign_mats();

		//thermal analysis
		body.assign_analysis(THERMAL);

		GapLog << "Applying thermal BCs ... ";
		
		//temporarily remove faces from the face set heatfluxsetname that have a
		//non-zero heat_flux value to prevent the default constraint from being applied
		fem::faceset heatfluxset_backup(body.facesets[heatfluxsetname]);

		for(int i=(int) heat_flux.size()-1; i>=0; i--)
		{
			if(fabs(heat_flux[i]) > 1.0e-3)	//a 'zero' tolerance
			{
				body.facesets[heatfluxsetname].erase(body.facesets[heatfluxsetname].begin()+i);
			}
		}

		//apply constraints
		body.assign_bcs();

		//save the thermal boundaries for the vtk
		vector<double> dirchlet;
		vector<double> neumann;
		vector<double> mixed_h;
		vector<double> mixed_phi;
		body.get_thermal_bcs(dirchlet, neumann, mixed_h, mixed_phi, heatfluxsetname, heatfluxset_backup, heat_flux);
		
		//reset the heatfluxset faceset
		body.facesets[heatfluxsetname] = heatfluxset_backup;

		//apply loads
		body.apply_faceset_nbc(heatfluxsetname, heat_flux);
		GapLog << "done." << endl << endl;
		
		//build the stiffness matrix
		GapLog << "Assembling the thermal stiffness matrix... " << endl;
		body.stiffnessMatrix();
		GapLog << endl;

		//solve the system
		GapLog << "Solving the thermal linear system..." << endl;
		body.solve();
		GapLog << endl;

		//obtain avg element temps
		vector<double> eleTemp(body.get_elecnt(), 0);
		for(int e=0; e<body.get_elecnt(); e++)
		{
			eleTemp[e] = body.getele_temp(e);
		}

		//also save the nodal temps
		node_temp.resize(body.get_nodecnt(), 0);
		for(int n=0; n<body.get_nodecnt(); n++)
		{
			if(old_temp.size() == node_temp.size())
			{
				//relax the new node temperature
				node_temp[n] = old_temp[n] + alpha*(body.getnode_temp(n)-old_temp[n]);
			} else	{
				node_temp[n] = body.getnode_temp(n);
			}
		}

		//Reset the slipper for elastic
		GapLog << "Resetting the analysis ... ";
		body.reset();
		GapLog << "done." << endl << endl;

		body.assign_analysis(ELASTIC);

		GapLog << "Applying elastic BCs ... ";
		if(body.options.inrel == 1)
		{
			//set inertia relief
			body.inertia_relief = true;
		} else {
			body.assign_bcs();
		}

		//apply loads
		for(int e = 0; e<eleTemp.size(); e++)
		{
			body.apply_ele_therm(e, eleTemp[e]-body.options.reftemp);
		}
		GapLog << "done." << endl << endl;

		//build the stiffness matrix
		GapLog << "Assembling the elastic stiffness matrix... " << endl;
		body.stiffnessMatrix();
		GapLog << endl;

		//solve the system
		GapLog << "Solving the elastic linear system..." << endl;
		body.solve();
		GapLog << endl;

		//write output
		body.writeVTK(vtk_file, node_temp, dirchlet, neumann, mixed_h, mixed_phi);
		GapLog << endl;

		//save the nodal displacements
		node_deform.resize(body.get_nodecnt());
		for(int n=0; n<body.get_nodecnt(); n++)
		{
			node_deform[n] = body.getnode_deform(n);
		}

	};


	void thermoelastic::get_faceset(fem * body, const string name, vector<int> & nid, vector<vector<double> > & nodexyz, vector<vector<int> > & setfaces)
	{
		//explalination of vars

		//vector of global nids for this set to be used by future calls
		//vector<int> nid;

		//vector of set nodexyz
		//vector<vector<double> > nodexyz;

		//vector of vector of set face nid's
		//vector<vector<int> > setfaces;
		
		if(!body->checkFaceset(name))
		{
			return;
		}

		//get a pointer to the faceset
		fem::faceset * fs = & body->facesets[name];

		//build a map of the global nid to local set nid - used internally only
		map<int, int> glob2loc_nid;

		//populate the nid vector
		for(int i=0; i<fs->size(); i++)
		{
			vector<int> face_nid;

			face * f = & (*fs)[i];
			for(int j=0; j<f->nodes.size(); j++)
			{
				//check if this node exists in the local set
				if(glob2loc_nid.count(f->nodes[j]->id) == 0)
				{
					//add the node to the set map
					glob2loc_nid[f->nodes[j]->id] = (int) nid.size();

					//push the global id to the nid vector
					nid.push_back(f->nodes[j]->id);

					//push the xyz of the node
					vector<double> xyz(3);
					xyz[0] = f->nodes[j]->x();
					xyz[1] = f->nodes[j]->y();
					xyz[2] = f->nodes[j]->z();
					nodexyz.push_back(xyz);
				}

				//push the local nid into the face vector
				face_nid.push_back(glob2loc_nid[f->nodes[j]->id]);
			}

			//push the face into setface
			setfaces.push_back(face_nid);

		}


	};

	unsigned int thermoelastic::get_nodecnt(fem * body)
	{
		return body->nodecnt;
	}

	fem * thermoelastic::open()
	{
		fem * solid = new fem;

		//set the gapinput object
		solid->caspar_gap_input = caspar_gap_input;

		//load the options file
		GapLog << "Loading the options files ... ";
		solid->readoptions(option_file);
		GapLog << "done." << endl << endl;

		GapLog << "Loading \"" << solid->options.meshfile << "\"... " << endl;
		solid->loadinp();
		GapLog << endl;

		return solid;
	};

	void thermoelastic::close(fem * body)
	{
		delete body;
	};

}
