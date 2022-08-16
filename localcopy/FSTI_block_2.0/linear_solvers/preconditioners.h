# ifndef __preconditioners__
# define __preconditioners__

enum preconditioners
{
	IDENTITY,	// no preconditioner, identity matrix
	ILU,			// incomplete level 0 LU decomposition
	ILUT,			// incomplete k-epsilon LU decomposition
	ILDLT,		// incomplete level 0 LDLT decomposition (symmetric matrices)
	ILDLTT,		// incomplete k-epsilon LDLT decomposition (symmetric matrices)
	INV,			// approximate inverse preconditioner
	JACOBI		// simple diagonal preconditioner
};

# endif