# include <valarray>
# include <vector>
# include <iostream>
# include <string>
# include <fstream>
# include <algorithm>
# include "../../main/src/logger.h"
# include "femFunctions.h"
# include <iomanip>
#include "../../caspar_input/input.h"

using namespace std;

extern class input myinput;

void fem::loads (const char* path, const char* inputType) {

	int elmType;
	int intOp; 
	int nDof; 
	int nElements;
	int nNodes;
	int nc;
	int nst; 
	int nAllDof;
	int nnc;

	// ------------------------------ read general file ------------------------------ //
	
	readGeneral(path,elmType,intOp,nDof,nElements,nNodes,nc,nst);
	nAllDof = nNodes * nDof; 
	nnc = nAllDof - nc;		
	
	// ------------------------------ read loads file ------------------------------ //
	
	vector<int> iForces;
	vector<double> forces;	
	readLoads(path,nst,iForces,forces,inputType);

	// ------------------------- read constrains ------------------------- //

	vector<int> iConstrains;
	readConstraints(path,nc,iConstrains,inputType);

	// ---------------- check is there are some loads that are also constrained ------------ //
	
	Log << "Setting loads ... ";
	// remove loaded nodes if they are also constrined
	int ix;
	double temp;
	vector<double> correctLoads;
	for(int i=0;i<nst;i++)
	{
		ix = iForces[i];
		temp = forces[i];
		correctLoads.push_back(temp);
		for(int j=0;j<nc;j++)
		{
			if(ix==iConstrains[j])
			{
				correctLoads.pop_back();
			}
		}
	}
	nst=(int) correctLoads.size();


	ofstream outLoads(string(string(path) + string("/temp/piston/solver/loads.bin")).c_str(), std::ios::binary | std::ios::out);
	outLoads.write(reinterpret_cast<char*>(&nnc),sizeof(int));
	for(int i=0; i < nnc ; i++) {
		double v = correctLoads[i];
		outLoads.write(reinterpret_cast<char*>(&v),sizeof(double));
	}
	
	outLoads.close();
	
	Log << "done!" << "\n";


	/*sort(iConstrains.begin(),iConstrains.end());
	int iCons_min = iConstrains[0];
	int iCons_max = iConstrains[iConstrains.size() - 1];

	// remove loaded nodes if they are also constrined
	for(unsigned int i=0; i<iForces.size(); i++)
	{
		if(iForces[i] >= iCons_min && iForces[i] <= iCons_max)
		{
			vector<int>::iterator it = find(iConstrains.begin(),iConstrains.end(),iForces[i]);
			if(it != iConstrains.end()) {
				iForces.erase(iForces.begin() + i);
				forces.erase(forces.begin() + i);
				nst--;
			}
		}
	}

	// sort the array
	vector<double> correctLoads(nnc,0);
 
	for(int i_load = 0; i_load < nst; i_load++) {
		// force index is bigger than the maximum index of constrained dof
		if(iForces[i_load] > iCons_max) {
			int i_real = iForces[i_load] - iConstrains.size();
			correctLoads[i_real] = forces[i_load];
		}
		// force index is in-between in the list of constrained dof (is worth to search)
		else if(iForces[i_load] >= iCons_min && iForces[i_load] <= iCons_max) {
			int j = 0;
			while(iConstrains[j] <= iForces[i_load]) {
				j++;
			}
			int i_real = iForces[i_load] - j;
			correctLoads[i_real] = forces[i_load];
		}
		// iForces[i_load] < iCons_min
		else { 
			int i_real = iForces[i_load];
			correctLoads[i_real] = forces[i_load];
		}	
	}*/

}