# include <cstdlib>
# include <cstdio>
# include <fstream>
# include <iostream>
# include <string>
# include <vector>
# include "gmm/gmm_kernel.h"
# include "gmm/gmm_iter_solvers.h"
# include "gmm/gmm_inoutput.h"
# include "femFunctions.h" // namespace of fem functions
# include "../../main/src/logger.h"
# include <iomanip>


// solving the system using the ICCG solver

int fem::solve(const char* path, const int iterMax, const double convergeTol) {

	Log << "Reading the stiffness matrix ... ";

	std::ifstream in(std::string(std::string(path) + std::string("/temp/piston/solver/stiffnessMatrix.bin")).c_str(), std::ios::binary | std::ios::in);
	int size;
	in.read(reinterpret_cast<char*>(&size), sizeof(int));
	int nnz;
	in.read(reinterpret_cast<char*>(&nnz), sizeof(int));
	int* row = new int[nnz];
	int* col = new int[nnz];
	double* val = new double[nnz];
	in.read(reinterpret_cast<char*>(row), nnz*sizeof(int));
	in.read(reinterpret_cast<char*>(col), nnz*sizeof(int));
	in.read(reinterpret_cast<char*>(val), nnz*sizeof(double));
	in.close();

	gmm::col_matrix<gmm::wsvector<double>> _K(size,size);
	for(int i=0; i<nnz; i++)
	{
		_K(row[i],col[i]) = val[i];
	}
	gmm::csc_matrix<double> K;
	gmm::copy(_K,K);

	delete [] row;
	delete [] col;
	delete [] val;

	Log << "done!" << "\n";
		
	// --------------------------- read the load vector --------------------------- //

	std::vector<double> b(size);
	
	Log << "Reading the load vector ... ";
	std::ifstream bIn (std::string(std::string(path) + std::string("/temp/piston/solver/loads.bin")).c_str(), std::ios::binary | std::ios::out);
	if (!bIn.is_open()) {
		std::cerr << "Could not read the loads vector!" << endl;
		exit(1);
	}
	int bsize;
	bIn.read(reinterpret_cast<char*>(&bsize), sizeof(int));
	if (bsize != size) {
		
		Log << "b size different from A size! A size is " << size << " b size is " << bsize << endl;
		exit(1);
	}

	double* _b = new double[size];
	bIn.read(reinterpret_cast<char*>(_b), size*sizeof(double));

	for(int i=0; i<size; i++)
	{
		b[i] = _b[i];
	}

	delete [] _b;

	Log << "done!" << "\n";

	// --------------------------- solve the system --------------------------- //

	std::vector<double> x(size,0);

	clock_t begin;
	clock_t end;

	gmm::iteration iter(convergeTol);
	iter.set_maxiter(iterMax);
	iter.set_noisy(1);

	begin = clock();
	// symmetric ilu preconditioner
	gmm::ildlt_precond<gmm::csc_matrix<double>>P(K);
	gmm::identity_matrix PS;

	Log << "Solving the linear system ... ";
	gmm::cg(K, x, b, PS, P, iter); // Conjugate gradient

	end = clock();

	if(iter.converged())
	{
		int tempint = iter.get_iteration();
		Log << "done in " << tempint << " iterations!" << "\n";
		double tempdouble = (end - begin)/1000;
		tempint = iter.get_res();
		Log << "Final relative residual " << tempint 
				  << " Solution time: " << tempdouble << "[s]" << "\n";
	}
	else
	{
		Log << "problem in solving the linear system!" << endl;
		exit(1);
	}
	
	Log << "Writing the solution vector ... ";
	std::ofstream out(std::string(std::string(path) + std::string("/temp/piston/solver/x_tmp.bin")).c_str(), std::ios::binary | std::ios::out);
	for(int i=0; i<size; i++)
	{
		double val = x[i];
		out.write(reinterpret_cast<char*>(&val),sizeof(double));
	}
	out.close();
	Log << "done!" << "\n";

	return 0;
	
}
