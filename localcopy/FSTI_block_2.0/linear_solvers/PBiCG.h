# ifndef __PBiCG__
# define __PBiCG__

// # define GMM_USES_LAPACK	// use the atals library

# include <gmm/gmm_kernel.h>
# include <gmm/gmm_iter_solvers.h>
# include "./preconditioners.h"


class PBiCG
{

	// system size
	int size;

	// const pointer to the sparse matrix
	const gmm::col_matrix<gmm::wsvector<double>>* A;
	// const pointer to the rhs
	const std::vector<double>* b;

	// sparse matrix used to effectively sove the system
	gmm::csc_matrix<double>* M;
	
	// precondition matrices
	gmm::ilu_precond<gmm::csc_matrix<double>>* P_ilu;
	gmm::ilut_precond<gmm::csc_matrix<double>>* P_ilut;	
	gmm::diagonal_precond<gmm::csc_matrix<double>>* P_diag;
	
	// identity matrix
	gmm::identity_matrix PS;

public:

	// solution vector 
	std::vector<double> x;

	// constructor
	PBiCG
	(
		const gmm::col_matrix<gmm::wsvector<double>>* A,
		const std::vector<double>* b
	);
	
	// destructor
	~PBiCG();

	// update the system size and M matrix
	void update();
	
	// calculate the preconditioner
	void preconditioner();

	// initialize the solution vector to zero
	void initialize_x();

	// clear M and P
	void clear();

	// solve the system
	int solve();

	void write();

	// get normalized residual ||b - Ax||/||b||
	double residual();

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