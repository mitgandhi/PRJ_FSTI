# include <vector>
# include <algorithm>
# include <iomanip>
# include <fstream>
# include <string>
# include "../../main/src/logger.h"
# include "femElement.h"
# include "femFunctions.h"

using namespace std;


void fem::stiffness(const char* path, const char* inputType) {

	int elmType;
	int intOp; 
	int nDof; 
	int nElements;
	int nNodes;
	int nc;
	int nst; 
	int nAllDof;
	int nnc;	

	// ------------------------------ setting up input ------------------------------ //
	
	readGeneral(path,elmType,intOp,nDof,nElements,nNodes,nc,nst);

	int elementNodes = 0;
	
	if(elmType == BRICK)
		elementNodes = 8;
	else if(elmType == TETRA)
		elementNodes = 4;

	nAllDof = nNodes * nDof;		// Number of all DOF -> number of row or column of global matrix
	nnc = nAllDof - nc;		
	vector<double> x,y,z;			// arrays for nodes coordinate
	vector<int> connect;			// array for connections
	vector<int> iConstrains;		// array for global index of constrained nodes
	vector<double> young;			// array for young modulus
	vector<double> poisson;			// array for poisson ratio

	// ------------------------- setting up coordinates vectors ------------------------- //
	
	readGeometry(path,nNodes,x,y,z,inputType);
	
	// ------------------------- setting up connectivity vector ------------------------- //
	
	readConnectivity(path,elmType,nElements,connect,inputType);

	// ------------------------------ setting up material ------------------------------ //

	readMaterials(path,nElements,young,poisson,inputType);
	
	// ------------------------------ setting up constrains ----------------------------- //

	readConstraints(path,nc,iConstrains,inputType);
	sort(iConstrains.begin(),iConstrains.end());
	int iConstr_min = iConstrains[0];
	int iConstr_max = iConstrains[(int) iConstrains.size() - 1];
	

	// ---------------------------------------------------------------------------------- //

	// the global stiffness matrix is assembled element by element
	// the stiffness matrix is going to be stored in compress row format:
	// rowIdx	colIdx	value
	// because of the stiffness matrix is symmetric, only the lower triangular part will be stored 

	// 2D dimension vector, for each row are stored the corresponding column indexes
	
	vector<vector<int>> columns(nnc);
	vector<vector<double>> val(nnc);
		

	// resize to the actual dimension of the stiffness matrix:
	//	"nAllDof" is the total number of the degrees of freedom, the stiffness matrix would have
	//	this dimension if no constrains are specified. So if we have "nc" dof constrained, the dimension
	//	of the matrix will be nnc = nAllDof - nc
	
	
	int nnz = 0;	// number of non zero value in stiffness matrix

	// ------------------ this is for the operation progress information ------------------ //
	double singleStep = 1.0; // percentage
	double nextStep = singleStep; // percentage

	//streamsize default_value = Log.precision(); dwm logger
	//cout.precision(1); dwm logger
	double tempdouble = 0.00;
	Log << "\nAssembling the stiffness matrix " << fixed << tempdouble << " %";

	// ------------------------------------------------------------------------------------ //
	
	
	for (int i_element = 0; i_element < nElements; i_element++) {

		// ------------------ gives information on operation progress
		double current = (double) i_element/nElements;
		if(current > nextStep/100) {
			nextStep += singleStep;
			tempdouble = 100*current;
			Log << "\rAssembling the stiffness matrix " << fixed << tempdouble << " %";
		}
		
		// ----------------- calculate the stiffness matrix for one element ----------------- //


		double* xEle = new double[elementNodes]; 
		double* yEle = new double[elementNodes]; 
		double* zEle = new double[elementNodes]; 
		int* cEle = new int[elementNodes];

				
		for(int i_node = 0; i_node < elementNodes; i_node++) {
			xEle[i_node] = x[connect[i_element*elementNodes + i_node]];
			yEle[i_node] = y[connect[i_element*elementNodes + i_node]];
			zEle[i_node] = z[connect[i_element*elementNodes + i_node]];
			cEle[i_node] = connect[i_element*elementNodes + i_node];
		}

		femElement* element = 0;

					
		if(elmType == BRICK) {
			element = new brick(intOp, young[i_element], poisson[i_element], xEle, yEle, zEle, cEle);
		}
		else if(elmType == TETRA) {
			element = new tetra(young[i_element], poisson[i_element], xEle, yEle, zEle, cEle);
		}
		else {
			Log << "Wrong element type, use BRICK or TETRA!" << endl;
			exit(1);
		}

		element -> calcStiff();	

		
		
		// store the global indexes corresponding to the nodes in the element considering constrains!
		int* globalIdx = new int[nDof*elementNodes];
		memcpy(globalIdx,element -> globalNumbering,nDof*elementNodes*sizeof(int));

		// array to store information about how much we have to shift the global numbering
		// considering the dof constrained before each global index
		int* globalShift = new int[nDof*elementNodes];
		memset(globalShift,0,nDof*elementNodes*sizeof(int));

		for(int i=0; i<nDof*elementNodes; i++) {
			if(globalIdx[i] >= iConstr_min && globalIdx[i] <= iConstr_max) {
				vector<int>::iterator it = find (iConstrains.begin(), iConstrains.end(), globalIdx[i]);
				if(it != iConstrains.end()) {
					globalIdx[i] = -1;
				}
			}
		}

		for(int i=0; i<nDof*elementNodes; i++) {
			if(globalIdx[i] < iConstr_min) {
				globalShift[i] = 0;
			}
			else if(globalIdx[i] >= iConstr_min && globalIdx[i] <= iConstr_max) {
				int c = 0;
				while(iConstrains[c] <= globalIdx[i])
					c++;
				globalShift[i] = c;
			}
			else if(globalIdx[i] > iConstr_max) {
				globalShift[i] = (int) iConstrains.size();
			}
		}		

		int min = 0, max = 0;

		// loop for all the local nodes in the element
		for(int j_local = 0; j_local < nDof*elementNodes; j_local++) {
			// check if the column must be escluded because of constraint
			if(globalIdx[j_local] > -1) { // the column doesn't have to be excluded	
				// index in row array considering excluded rows						
				int thisColumn = globalIdx[j_local] - globalShift[j_local]; 
				for(int i_local = 0; i_local < nDof*elementNodes; i_local++) {
					// store just the lower triangular part
					//if(globalIdx[i_local] >= globalIdx[j_local]) {  
						// check if the row must be escluded because of constrain
						if(globalIdx[i_local] > -1) { // the row doesn't have to be excluded	
							int thisRow = globalIdx[i_local] - globalShift[i_local];
							// look for thisColumns in columns[thisColumn] vector
							vector<int>::iterator it = 
								find (columns[thisColumn].begin(), columns[thisColumn].end(),thisRow);			
							if(it == columns[thisColumn].end()) { // the row is not present
								nnz++;	// we are going to store a new element
								columns[thisColumn].push_back(thisRow);	// store a new index
								val[thisColumn].push_back(element -> stiffMat[i_local][j_local]); // store a new value
							}
							else { // the column is already present
								int idx = int(it - columns[thisColumn].begin());	// get the index of the column
								val[thisColumn][idx] += element -> stiffMat[i_local][j_local]; // add the value
							}
						}
					//}
				} // end of "for(int j_local = 0; j_local<24; j_local++)"

			} // end of "if (it == myvector.end())"
		
		} // end of "for(int i_local = 0; i_local<24; i_local++)" 
		
		
		delete element;
		delete [] globalIdx;
		delete [] globalShift;
				
		delete [] xEle;
		delete [] yEle;
		delete [] zEle;
		delete [] cEle;
			
	}

	tempdouble = 100.0;
	Log << "\rAssembling the stiffness matrix " << fixed << tempdouble << " %" << "\n";
	
	// --------------------------------------------------------------------------------------------------- //

	Log << "Writing the stiffness matrix file ... " ;

	ofstream out(string(string(path) + string("/temp/piston/solver/stiffnessMatrix.bin")).c_str(), std::ios::binary | std::ios::out);
	out.write(reinterpret_cast<char*>(&nnc),sizeof(int));
	out.write(reinterpret_cast<char*>(&nnz),sizeof(int));
	for(unsigned int i=0; i<(int) columns.size(); i++)
	{
		for(unsigned int j=0; j<(int) columns[i].size(); j++)
		{
			int idx = i;
			out.write(reinterpret_cast<char*>(&idx),sizeof(int));
		}
	}
	for(unsigned int i=0; i<(int) columns.size(); i++)
	{
		for(unsigned int j=0; j<(int) columns[i].size(); j++)
		{
			int idx = columns[i][j];
			out.write(reinterpret_cast<char*>(&idx),sizeof(int));
		}
	}
	for(unsigned int i=0; i<(int) columns.size(); i++)
	{
		for(unsigned int j=0; j<(int) columns[i].size(); j++)
		{
			double v = val[i][j];
			out.write(reinterpret_cast<char*>(&v),sizeof(double));
		}
	}

	out.close();

	Log << "done!" << "\n";

	//cout.precision(default_value); dwm logger



	/*ofstream out(string(string(path) + string("/solver/stiffnessMatrix.bin")).c_str(), std::ios::binary | std::ios::out);
		
	// number of rows 
	out.write(reinterpret_cast<char*>(&nnc),sizeof(int)); 
	// number of columns
	out.write(reinterpret_cast<char*>(&nnc),sizeof(int));
	// flags
	int flags = TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER;
	out.write(reinterpret_cast<char*>(&flags),sizeof(int));

	// write columns pointer 
	int col_ptr = 0;
	out.write(reinterpret_cast<char*>(&col_ptr),sizeof(int));
	for(int i = 0; i < nnc; i++) {
		col_ptr += columns[i].size();
		out.write(reinterpret_cast<char*>(&col_ptr),sizeof(int));
	}
	
	// write row indexes 
	for(int i = 0; i < nnc; i++) {
		vector<int> ord(columns[i]);
		sort(ord.begin(),ord.end());
		for(unsigned int j=0; j < ord.size(); j++) {
			int idx = ord[j];
			out.write(reinterpret_cast<char*>(&idx),sizeof(int));
		}
	}

	// write matrix entries 
	for(int i = 0; i < nnc; i++) {
		vector<int> ord(columns[i]);
		sort(ord.begin(),ord.end());
		for(unsigned int j=0; j < ord.size(); j++) {
			vector<int>::iterator it;
			it = find(columns[i].begin(),columns[i].end(),ord[j]);
			int idx = it - columns[i].begin();
			double value = val[i][idx];
			out.write(reinterpret_cast<char*>(&value),sizeof(double));
		}
	}
		
	out.close();*/


}

