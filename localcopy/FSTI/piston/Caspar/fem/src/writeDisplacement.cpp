# include "../../main/src/logger.h"
# include <fstream>
# include <vector>
# include <algorithm>
# include <iostream>
# include <string>
# include <iomanip>
# include "femFunctions.h"

using namespace std;

void fem::writeDisplacement(const char* path, const char* inputType) {

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

	// ------------------------- read constrains ------------------------- //

	vector<int> iConstrains;
	readConstraints(path,nc,iConstrains,inputType);
	
	int min_ic = *min_element(iConstrains.begin(),iConstrains.end());
	int max_ic = *max_element(iConstrains.begin(),iConstrains.end());

	Log << "Writing the displacement ... ";

	std::ifstream in(string(string(path) + string("/temp/piston/solver/x_tmp.bin")).c_str(), ios::in | ios::binary);
	if (!in.is_open()) {
		Log << "Error opening result file from solver!\n" << endl;
		//system("pause");
		exit(1);
	}
	
	double* displ = new double[nnc];
	in.read(reinterpret_cast<char*>(displ), nnc*sizeof(double));
	in.close();

	double* allDispl = new double[nAllDof];

	// riassemble the displacement vector
	for(int i = 0, k = 0; i < nAllDof; i += nDof) {
		for(int j = 0; j < nDof; j++) {
			if( (i + j) >= min_ic && (i + j) <= max_ic ) {
				vector<int>::iterator it;
				it = find(iConstrains.begin(), iConstrains.end(), i + j);
				if(it == iConstrains.end()) { // the displacement is not constrained
					allDispl[i+j] = displ[k]; 
					k++;
				}
				else {
					allDispl[i+j] = 0.0;
				}
			}
			else {
				allDispl[i+j] = displ[k]; // the displacement is not constrained
				k++;
			}
		}
	}
	
	/*string displfile = string(string(path) + string("/output/displacement.fem"));
	ofstream out(displfile.c_str(), ios::binary | ios::out);

	if (!out.is_open()) 
	{
		cout << "Error opening " << displfile << "\n";
		system("pause");
		exit(1);
	}
	out.write(reinterpret_cast<char*>(&nAllDof),sizeof(int));
	out.write(reinterpret_cast<char*>(allDispl),nAllDof*sizeof(double));
	
	out.close();
	*/

	string displfile = string(string(path) + string("/temp/piston/output/displacement.fem"));
	ofstream fout(displfile.c_str());
	for(int i=0; i<nNodes; i++)
		fout << allDispl[3*i] << "\t" << allDispl[3*i + 1] << "\t" << allDispl[3*i + 2] << "\n";
	fout.close();
	fout.clear();

	delete[] displ;
	delete[] allDispl;

	Log << "done!" << "\n";

}