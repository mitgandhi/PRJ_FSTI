# include "./PBiCG.h"
# include <gmm/gmm_inoutput.h>
# include <vector>

PBiCG::PBiCG
(
	const gmm::col_matrix<gmm::wsvector<double>>* _A,
	const std::vector<double>* _b
) 
: 
	size(_A->ncols()),
	A(_A),
	b(_b),
	x(size,0.0),
	P_ilu(0),
	P_ilut(0),
	P_diag(0),
	M(0)
{
	system.tol = 1e-6;
	system.itermax = 3000;
	system.quiet = true;
	system.P = ILU;

	// populate M
	M = new gmm::csc_matrix<double> (size, size);
	gmm::copy(*A,*M);
}
// ------------------------------------------------------------------------- //
PBiCG::~PBiCG()
{
	if(P_ilu != 0)
		delete P_ilu;
	if(P_ilut != 0)
		delete P_ilut;
	if(P_diag != 0)
		delete P_diag;
	if(M != 0)
		delete M;
}
// ------------------------------------------------------------------------- //
void PBiCG::update()
{
	int oldsize = size;
	size = A->ncols();
	
	if(size != oldsize)
		x.resize(size,0.0);

	if(M != 0)
	{
		delete M;
		M = 0;
	}

	M = new gmm::csc_matrix<double> (size, size);
	gmm::copy((*A), (*M));

	// delete the preconditioner
	if(P_ilu != 0)
	{
		delete P_ilu;
		P_ilu = 0;
	}
	if(P_ilut != 0)
	{
		delete P_ilut;
		P_ilut = 0;
	}
	if(P_diag != 0)
	{
		delete P_diag;
		P_diag = 0;
	}
}
// ------------------------------------------------------------------------- //
void PBiCG::preconditioner()
{
	// ILU preconditioner
	if(system.P == ILU)
	{
		if(P_ilu != 0)
		{
			delete P_ilu;
			P_ilu = 0;
		}
		P_ilu = new gmm::ilu_precond<gmm::csc_matrix<double>>(*M);
	}
	// ILUT preconditioner
	if(system.P == ILUT)
	{
		if(P_ilut != 0)
		{
			delete P_ilut;
			P_ilut = 0;
		}
		P_ilut = new gmm::ilut_precond<gmm::csc_matrix<double>>(*M,10, 1e-12);
	}
	// JACOBI preconditioner
	if(system.P == JACOBI)
	{
		if(P_diag != 0)
		{
			delete P_diag;
			P_diag = 0;
		}
		P_diag = new gmm::diagonal_precond<gmm::csc_matrix<double>>(*M);
	}
}
// ------------------------------------------------------------------------- //
void PBiCG::clear()
{
	if(M != 0)
	{
		delete M;
		M = 0;
	}
	// delete the preconditioner
	if(P_ilu != 0)
	{
		delete P_ilu;
		P_ilu = 0;
	}
	if(P_ilut != 0)
	{
		delete P_ilut;
		P_ilut = 0;
	}
	if(P_diag != 0)
	{
		delete P_diag;
		P_diag = 0;
	}
}
// ------------------------------------------------------------------------- //
void PBiCG::initialize_x()
{
	gmm::clear(x);
}
// ------------------------------------------------------------------------- //
int PBiCG::solve()
{

	gmm::iteration iter(system.tol);
	iter.set_maxiter(system.itermax);
	
	if(system.quiet)
		iter.set_noisy(0);
	else
	{
		iter.set_noisy(1);
		// set cout to scientific
		cout.setf(std::ios::scientific,std::ios::floatfield);  
	}

	// Biconjugate Grasient stabilized , with ILU preconditioner
	if(system.P == ILU && P_ilu != 0)
	{
		gmm::bicgstab((*M), x, (*b), (*P_ilu), iter);	
	}
	// Biconjugate Grasient stabilized , with ILUT preconditioner
	else if(system.P == ILUT && P_ilut != 0)
	{
		gmm::bicgstab((*M), x, (*b), (*P_ilut), iter);	
	}
	// Biconjugate Grasient stabilized , with JACOBI preconditioner
	else if(system.P == JACOBI && P_diag != 0)
	{
		gmm::bicgstab((*M), x, (*b), (*P_diag), iter);	
	}
	else
	{
		cout << "\nPBiCG::solve() Warning: solving with no PC." << endl;
		gmm::bicgstab((*M), x, (*b), PS, iter);	
	}

	if(!system.quiet)
	{
		// reset cout to fixed
		cout.setf(std::ios::fixed, std::ios::floatfield);  
	}
	
	return iter.get_iteration();
}
// ------------------------------------------------------------------------- //
double PBiCG::residual()
{

	std::vector<double> res(size);
	std::vector<double> Ax(size);

	gmm::mult(*M, x, Ax);	// M1 * V2 --> V1
	gmm::add((*b), gmm::scaled(Ax, -1.0), res); // V1 - 2.0 * V2 --> V2
	double num = 0;
	double den = 0;
	for(int i = 0; i < size; i++)
		num += pow(res[i],2.0), den += pow((*b)[i],2.0);

	if(den != 0)
		return sqrt(num)/sqrt(den);
	else
		return sqrt(num);
}
// ------------------------------------------------------------------------- //
void PBiCG::write()
{
	MatrixMarket_save("./M.txt", *M);
}
// ------------------------------------------------------------------------- //