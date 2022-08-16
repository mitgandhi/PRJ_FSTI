#include "fem.h"

namespace CasparSlipperFEM
{
	void fem::numberDOF()
	{
		//define nodal degrees of freedom
		int nDOF;
		if(analysis_type == ELASTIC)
		{
			nDOF = 3;
		} else if (analysis_type == THERMAL)
		{
			nDOF = 1;
		}

		//Clear ebcs for inertia relief
		if(inertia_relief)
		{
			ebcs.clear();
		}

		//Define the global DOF for unconstrained nodes
		uDOF = 0;
		DOF = 0;
		for(int n=0; n<nodecnt; n++)
		{
			nodes[n].DOF.resize(nDOF);
			for(int d=0; d<nDOF; d++)
			{
				if(ebcs.count(niddof(n,d)) == 0)	//check to see if the nodal DOF is constrained
				{
					//check to see if there is a RBE connecting nodal DOF
					if(rbe.count(niddof(n,d)) > 0)
					{
						//RBE will have two niddof that become constrained together.
						//we will default the lower node id to be the independent node with the higher id, the dependent id

						niddof * rbenode2 = &rbe[niddof(n,d)];
						if(n < rbenode2->first)
						{
							//this node is the independent node, treat as a normal node
							nodes[n].DOF[d] = uDOF;
							uDOF++;
						} else {
							//this is the dependent node, so assign the global DOF to be the same as the independent node
							nodes[n].DOF[d] = nodes[rbenode2->first].DOF[rbenode2->second];
						}
					} else {
						//normal node, not essential constrained, no RBE
						nodes[n].DOF[d] = uDOF;
						uDOF++;
					}
				} else {
					nodes[n].DOF[d] = -1;
				}
				DOF++;
			}
		}

		//Add 3 DOF for lambda if using IR
		if(inertia_relief)
		{
			uDOF += 3;
		}

		//resize the x, X vector
		X.resize(uDOF, 0);

		//Estimate the RAM based on element count & element type
		double ram = 0;
		for(int e=0; e<elements.size(); e++)
		{
			if(elements[e]->element_type == element::TETRA)
			{
				//tetra's have 4 nodes:
				ram += (nDOF*nDOF*4*4) * 14; //DOF in element stiffness matrix * 14 byte estimate
			}
		}
		
		GapLog << "\tNumber of DOF = " << uDOF << endl;
		//GapLog << "\tEstimated RAM for K matrix = " << ram/1024/1024 << " MB (Could be very wrong)" << endl;
	}

	void fem::setLoads()
	{
		//resize the load vector
		b.resize(uDOF, 0);
	};

	void fem::clearElements()
	{
		for(unsigned int i=0; i<elements.size(); i++)
		{
			delete elements[i];
		}
		elements.clear();
		elecnt = 0;
	};


	void fem::checkAnalysis()
	{
		//This method checks the current state of the analysis to see if we can proceed 
		//to building K and solve

		if(analysis_type == NONE)
		{
			GapLog << endl << endl;
			warning("Analysis type was specified as none. Is this what you want?");
			GapLog << "Finished processing input, exiting." << endl;
			exit(1);
		}

		if(analysis_type == THERMAL)
		{
			//inertia relief applies only to elastic analysis, so force false
			inertia_relief = false;
		}

	};

	bool fem::checkNodeset(const string ns)
	{
		if(nodesets.count(ns) == 0)
		{
			warning("node set = " + n2s(ns) + " does not exist!");
			return false;
		}
		return true;
	};

	bool fem::checkFaceset(const string fs)
	{
		if(facesets.count(fs) == 0)
		{
			warning("face set = " + n2s(fs) + " does not exist!");
			return false;
		}
		return true;
	};

	bool fem::checkElementset(const string es)
	{
		if(elementsets.count(es) == 0)
		{
			warning("element set = " + n2s(es) + " does not exist!");
			return false;
		}
		return true;
	};

	void fem::warning(const string w)
	{
			GapLog << "#----------------------------------------------#" << endl;
			GapLog << "# WARNING: " << w << endl;
			GapLog << "#----------------------------------------------#" << endl;
	};

	void fem::error(const string e)
	{
			GapLog << "#----------------------------------------------#" << endl;
			GapLog << "# ERROR: " << e << endl;
			GapLog << "#----------------------------------------------#" << endl;
			exit(1);
	};

	fem::fem()
	{
		nodecnt = (int) nodes.size();
		elecnt = (int) elements.size();

		//Set default variables
		reset();
	}

	fem::~fem()
	{
		clearElements();
	};

	void fem::reset()
	{
		//used to reset the analysis but DOES NOT clear input geometry
		analysis_type = NONE;

		//Default IR to false
		inertia_relief = false;

		//Remove ebc's
		ebcs.clear();

		//Clear the main var's
		b.clear();
		X.clear();
		gmm::clear(K);
		uDOF = 0;
		DOF = 0;

		//clear node dof's
		for(int n = 0; n<nodecnt; n++)
		{
			nodes[n].DOF.clear();
		} 

		//clear element analysis type
		for(int e = 0; e<elecnt; e++)
		{
			elements[e]->assign_analysis(NONE);
		} 

	};

};
