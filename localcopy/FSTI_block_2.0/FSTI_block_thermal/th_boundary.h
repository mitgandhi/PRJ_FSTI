# ifndef __th_boundary__
# define __th_boundary__

# include <vector>
# include <string>

enum th_bc_type
{
	dirichlet,
	mixed,
	neumann
};

struct th_boundary
{
	
	std::string setname;	// name of the set where the bc is applied
	th_bc_type bctype;			// type of th_boundary condition

	// --------- Dirichlet ---------- //
	// dirichlet th_boundary is a constraint, which can be
	// a specified node position or temperature
	
	struct 
	{
		// for thermal use x = y = z = -1
		int x;
		int y;
		int z;
		double Ts;	
	} dirichlet;

	// --------- Neumann ---------- //
	// heat flux, uniform (size = 1) or non-uniform
	struct
	{
		std::vector<double> q;
	} neumann;

	// --------- Mixed ---------- //
	// convection coeff and surrounding T
	struct
	{
		double h;
		double Tinf;
	} mixed;
};


# endif