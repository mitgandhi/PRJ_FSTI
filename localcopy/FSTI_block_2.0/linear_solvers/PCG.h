# ifndef __PCG__
# define __PCG__

// # define GMM_USES_LAPACK	// use the atals library

# include <gmm/gmm_kernel.h>
# include <gmm/gmm_iter_solvers.h>
# include "preconditioners.h"
# include <vector>

class PCG
{

	// system size
	int size;

	// const pointer to the sparse matrix
	const gmm::col_matrix<gmm::wsvector<double>>* A;
	// const pointer to the rhs
	const std::vector<double>* b;

	// sparse matrix effectively used to solve the linear system
	gmm::csc_matrix<double>* M;
	
	// precondition matrices
	gmm::ildlt_precond<gmm::csc_matrix<double>>* P_ildlt;
	gmm::ildltt_precond<gmm::csc_matrix<double>>* P_ildltt;	
	gmm::diagonal_precond<gmm::csc_matrix<double>>* P_diag;
	gmm::mr_approx_inverse_precond<gmm::csc_matrix<double>>* P_inv;
	
	// identity matrix
	gmm::identity_matrix PS;

public:

	// solution vector 
	std::vector<double> x;

	// constructor
	PCG
	(
		const gmm::col_matrix<gmm::wsvector<double>>* A,
		const std::vector<double>* b
	);
	
	// destructor
	~PCG();

	// update the system size and M matrix
	void update();
	
	// calculate the preconditioner
	void preconditioner(int fill_in = 1, double threshold = 1e-6);

	void initialize_x();

	// clear the M and P
	void clear();

	// solve the system
	int solve();

	// get normalized residual ||b - Ax||/||b||
	double residual();

	// write the matrix
	void write();

	// solver settings
	struct
	{
		double tol;					// solver tolerance
		int itermax;				// max iterations
		bool quiet;					// quiet mode
		preconditioners P;	// preconditioner type	
	} system;
	
};

# endif