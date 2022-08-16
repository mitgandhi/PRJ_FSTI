// ------------------------------------------------------------------------- //
// -------------------- class to solve the energy equation ----------------- //
// ------------------------------------------------------------------------- //

# ifndef __energy__
# define __energy__

# include "./scalar_field.h"
# include "./block_gap_main.h"
# include "../linear_solvers/PBiCG.h"

class Energy
{

	friend class cblock_gap_main;

	cblock_gap_main&	gap;	// reference to the main gap structure
	
	double TLP;							// low pressure line temperature
	double THP;							// high pressure line temperature
	double Tcase;						// case temperature

	// -------------------------------- fields ------------------------------- //

	const vector_field& V;		// fluid velocity
	const vector_field& Vc;		// fluid velocity (couette)
	const scalar_field& mu;		// fluid viscosity
	const scalar_field& hcb;	// film thickness (block)
	const scalar_field& hvp;	// film thickness (valve plate)
	scalar_field& phid;				// dissipation field
	scalar_field& T;					// fluid temperature
	
	// ---------------------------- linear system ---------------------------- //

	int size;																		// linear system size
	gmm::col_matrix<gmm::wsvector<double>> A;		// fvm matrix
	std::vector<double> b;											// known term
	PBiCG solver;																// linear solver

	// ---------------------------- private functions ------------------------ //

	// calculate the diffusive term at the boundary faces
	std::vector<double> getd(int i, int j, int k);
	// populate the fvm matrix, and calculate the b term
	void discretize();	

public:

	// ----------------------------- public functions ------------------------ //

	Energy(cblock_gap_main& gap_main);		// constructor
	~Energy();													// destructor
	void update_dimension();						// update the linear system dimension
	int solve();												// solve the energy equation 
	
};

# endif
