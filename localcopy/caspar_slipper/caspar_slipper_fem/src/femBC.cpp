#include "fem.h"

namespace CasparSlipperFEM
{
	void fem::apply_ebc(const string ns, const int dof, const double val)
	{
		if(!checkNodeset(ns))
		{
			//Nodeset doesn't exist so return
			return;
		}
		for(int i=0; i<nodesets[ns].size(); i++)
		{
			ebcs[niddof (nodesets[ns][i], dof)] = val;
		}
	}

	void fem::apply_faceset_nbc(const string fs, const double val)
	{
		if(!checkFaceset(fs))
		{
			//Faceset doesn't exist so return
			return;
		}

		for(int i=0; i<facesets[fs].size(); i++)
		{
			face * f = &facesets[fs][i];

			//get the elements that have this face
			vector<int> e = elementfaces[(*f)];
			if(e.size() > 1)
			{
				warning("Loading multiple elements for a face load. Is this what you want?");
			}

			if(e.size() == 0)
			{
				warning("Not loading any elements for a face load. Is this what you want?");
			}

			for(int j=0; j<e.size(); j++)
			{
				elements[e[j]]->assign_face_nbc(*f, val);
			}
		}
	};


	void fem::apply_faceset_nbc(const string fs, vector<double> vals)
	{
		if(!checkFaceset(fs))
		{
			//Faceset doesn't exist so return
			return;
		}

		if(facesets[fs].size() != vals.size())
		{
			//Vals don't match faceset size so impossible to load
			warning("Faceset " + n2s(fs) + " size doesn't match natural bc loading vector size!");
			return;
		}

		for(int i=0; i<facesets[fs].size(); i++)
		{
			face * f = &facesets[fs][i];

			//get the elements that have this face
			vector<int> e = elementfaces[(*f)];
			if(e.size() > 1)
			{
				warning("Loading multiple elements for a face load. Is this what you want?");
			}

			if(e.size() == 0)
			{
				warning("Not loading any elements for a face load. Is this what you want?");
			}

			for(int j=0; j<e.size(); j++)
			{
				elements[e[j]]->assign_face_nbc(*f, vals[i]);
			}
		}
	};

	void fem::apply_faceset_mbc(const string fs, const double h, const double phi)
	{
		if(!checkFaceset(fs))
		{
			//Faceset doesn't exist so return
			return;
		}

		for(int i=0; i<facesets[fs].size(); i++)
		{
			face * f = &facesets[fs][i];

			//get the elements that have this face
			vector<int> e = elementfaces[(*f)];
			if(e.size() > 1)
			{
				warning("Loading multiple elements for a face load. Is this what you want?");
			}

			if(e.size() == 0)
			{
				warning("Not loading any elements for a face load. Is this what you want?");
			}

			for(int j=0; j<e.size(); j++)
			{
				elements[e[j]]->assign_face_mbc(*f, h, phi);
			}
		}
	};

	void fem::apply_ele_therm(const int eid, const double t)
	{
		if(eid < 0 || eid > elecnt)
		{
			return;
		}
		elements[eid]->assign_ele_therm(t);
	};


	void fem::assign_E(const string elset, const double E)
	{
		if(!checkElementset(elset))
		{
			//Elementset doesn't exist so return
			return;
		}
		elementset * solid = & elementsets[elset];
		for(int i=0; i<solid->size(); i++)
		{
			elements[(*solid)[i]]->E = E;
		}
	};

	void fem::assign_nu(const string elset, const double nu)
	{
		if(!checkElementset(elset))
		{
			//Elementset doesn't exist so return
			return;
		}
		elementset * solid = & elementsets[elset];
		for(int i=0; i<solid->size(); i++)
		{
			elements[(*solid)[i]]->nu = nu;
		}
	};

	void fem::assign_rho(const string elset, const double rho)
	{
		if(!checkElementset(elset))
		{
			//Elementset doesn't exist so return
			return;
		}
		elementset * solid = & elementsets[elset];
		for(int i=0; i<solid->size(); i++)
		{
			elements[(*solid)[i]]->rho = rho;
		}
	};

	void fem::assign_k(const string elset, const double k)
	{
		if(!checkElementset(elset))
		{
			//Elementset doesn't exist so return
			return;
		}
		elementset * solid = & elementsets[elset];
		for(int i=0; i<solid->size(); i++)
		{
			elements[(*solid)[i]]->k = k;
		}
	};

	void fem::assign_alpha(const string elset, const double alpha)
	{
		if(!checkElementset(elset))
		{
			//Elementset doesn't exist so return
			return;
		}
		elementset * solid = & elementsets[elset];
		for(int i=0; i<solid->size(); i++)
		{
			elements[(*solid)[i]]->alpha = alpha;
		}
	};

	void fem::assign_analysis(analysis_types type)
	{
		analysis_type = type;

		for(int e=0; e<elecnt; e++)
		{
			elements[e]->assign_analysis(type);
		}
	};


	void fem::assign_mats()
	{
		for(int i=0; i<options.materials.size(); i++)
		{
			for(int j=0; j<options.materials[i].sets.size(); j++)
			{
				const string setname = options.materials[i].sets[j];

				assign_k(setname, options.materials[i].lambda);
				assign_E(setname, options.materials[i].E);
				assign_nu(setname, options.materials[i].nu);
				assign_alpha(setname, options.materials[i].alpha);
				assign_rho(setname, options.materials[i].rho);
			}
		}
	}

	void fem::assign_bcs()
	{
		//switch case depending on analysis type
		if(analysis_type == THERMAL)
		{
			for(int i=0; i<options.thermal_boundries.size(); i++)
			{
				femoptions::thermal_boundary * b = & options.thermal_boundries[i];

				if(b->type == b->DIRICHLET)
				{
					//we will assign the same node multiple times, but the ebcs map will forgive this

					//we will need to find the nodes in this face set
					faceset * fs = & facesets[b->set];
					for(int i=0; i<fs->size(); i++)	//for each face
					{
						for(int j=0; j<(*fs)[i].nodes.size(); j++) //for each face node
						{
							//dof = 0 because this is a thermal analysis
							ebcs[niddof((*fs)[i].nodes[j]->id, 0)] = b->Tp;
						}
					}
					
				} else if (b->type == b->NEUMANN)
				{
					apply_faceset_nbc(b->set, b->q);
				} else if (b->type == b->ROBBINS)
				{
					apply_faceset_mbc(b->set, b->h, b->Tinf);
				}

			}

		} else if (analysis_type == ELASTIC)
		{
			//loop through all the const's
			for(int i=0; i<options.elastic_constraints.size(); i++)
			{
				femoptions::elastic_constraint * c = & options.elastic_constraints[i];		
			
				//apply_ebc(c->set, 0, c->x);
				//apply_ebc(c->set, 1, c->y);
				//apply_ebc(c->set, 2, c->z);

				//the above method is dangerous because it conflicts with other interface standards
				//ie. to activate a constraint above you would set it = 0 for fixed
				//setting to 1 would specify a 1 meter displacement

				//instead here 1 turns a fixed constraint on
				if(c->x == 1)
				{
					apply_ebc(c->set, 0, 0);
				}

				if(c->y == 1)
				{
					apply_ebc(c->set, 1, 0);
				}

				if(c->z == 1)
				{
					apply_ebc(c->set, 2, 0);
				}

			}
		}
	}

	void fem::get_thermal_bcs(vector<double> &dirchlet, vector<double> &neumann, vector<double> &mixed_h, vector<double> &mixed_phi, const string heatfluxsetname, const fem::faceset &fullheatfluxset, const vector<double> &heat_flux_vals)
	{
		//only an option in THERMAL mode
		if(analysis_type == ELASTIC)
		{
			return;
		}

		dirchlet.clear();
		neumann.clear();
		mixed_h.clear();
		mixed_phi.clear();
		
		double NaN = std::numeric_limits<double>::quiet_NaN ();

		dirchlet.resize(nodecnt, NaN);
		neumann.resize(nodecnt, NaN);
		mixed_h.resize(nodecnt, NaN);
		mixed_phi.resize(nodecnt, NaN);

		vector<int> nodevalctr_neumann(nodecnt, 0);
		vector<int> nodevalctr_robbins(nodecnt, 0);

		// go through all the bc's
		for(int i=0; i<options.thermal_boundries.size(); i++)
		{
			femoptions::thermal_boundary * b = & options.thermal_boundries[i];

			if(b->type == b->DIRICHLET)
			{
				//we will assign the same node multiple times, but the ebcs map will forgive this

				//we will need to find the nodes in this face set
				faceset * fs = & facesets[b->set];
				for(int i=0; i<fs->size(); i++)	//for each face
				{
					for(int j=0; j<(*fs)[i].nodes.size(); j++) //for each face node
					{
						dirchlet[(*fs)[i].nodes[j]->id] = b->Tp;
					}
				}
					
			} else if (b->type == b->NEUMANN)
			{
				for(int i=0; i<facesets[b->set].size(); i++)
				{
					face * f = &facesets[b->set][i];

					for(int j=0; j<f->nodes.size(); j++)
					{
						//we want to push the value to the nodes
						const int nid = f->nodes[j]->id;
						if(_isnan(neumann[nid]))
						{
							neumann[nid] = 0;
						}
						neumann[nid] += b->q;
						nodevalctr_neumann[nid]++;
					}
				}

			} else if (b->type == b->ROBBINS)
			{
				

				for(int i=0; i<facesets[b->set].size(); i++)
				{
					face * f = &facesets[b->set][i];
															
					for(int j=0; j<f->nodes.size(); j++)
					{
						//we want to push the value to the nodes
						const int nid = f->nodes[j]->id;
						if(_isnan(mixed_h[nid]))
						{
							mixed_h[nid] = 0;
							mixed_phi[nid] = 0;
						}
						mixed_h[nid] += b->h;
						mixed_phi[nid] += b->Tinf;
						nodevalctr_robbins[nid]++;	
					}
				}
			}
		}

		//apply the fullheatflux set
		//the normal heatfluxset has currently been reduced to exclude faces that are having a flux applied from the gap
		for(int i=0; i<fullheatfluxset.size(); i++)
		{
			//only include faces that have a non-zero flux value
			if(fabs(heat_flux_vals[i]) > 1.0e-3)	//a 'zero' tolerance
			{
				const face * f = &fullheatfluxset[i];

				for(int j=0; j<f->nodes.size(); j++)
				{
					//we want to push the value to the nodes
					const int nid = f->nodes[j]->id;

					if(_isnan(neumann[nid]))
					{
						neumann[nid] = 0;
					}
					neumann[nid] += heat_flux_vals[i];
					nodevalctr_neumann[nid]++;				
				}
			}
		}
		
		//take the average
		for(int i=0; i<nodecnt; i++)
		{
			if(nodevalctr_neumann[i] > 0)
			{
				neumann[i] /= (double) nodevalctr_neumann[i];
			}

			if(nodevalctr_robbins[i] > 0)
			{
				mixed_h[i] /= (double) nodevalctr_robbins[i];
				mixed_phi[i] /= (double) nodevalctr_robbins[i];
			}
		}

}
};
