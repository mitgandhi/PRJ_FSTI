# include "PCG.h"
# include <vector>
# include <gmm/gmm_inoutput.h>

// ------------------------------------------------------------------------- //
PCG::PCG
(
	const gmm::col_matrix<gmm::wsvector<double>>* _A,
	const std::vector<double>* _b
) 
: 
	size(_A->ncols()),
	A(_A),
	b(_b),
	x(size,0.0),
	P_ildlt(0),
	P_ildltt(0),
	P_diag(0),
	M(0)
{
	system.tol = 1e-6;
	system.itermax = 3000;
	system.quiet = true;
	system.P = IDENTITY;

	// populate M
	M = new gmm::csc_matrix<double> (size, size);
	gmm::copy(*A,*M);
}
// ------------------------------------------------------------------------- //
PCG::~PCG()
{
	if(P_ildlt != 0)
		delete P_ildlt;
	if(P_ildltt != 0)
		delete P_ildltt;
	if(P_diag != 0)
		delete P_diag;
	if(M != 0)
		delete M;
}
// ------------------------------------------------------------------------- //
void PCG::update()
{
	int oldsize = size;
	size = A->ncols();
	
	if(size != oldsize)
		x.resize(size,0.0);

	if(M != 0)
		delete M;
	
	M = new gmm::csc_matrix<double> (size, size);
	gmm::copy((*A), (*M));

}
// ------------------------------------------------------------------------- //
void PCG::preconditioner(int fill_in, double threshold)
{
	
	// ILDLT preconditioner
	if(system.P == ILDLT)
	{
		//cout << "  * Calculating ILDLT preconditioner ... ";
		if(P_ildlt != 0)
		{
			delete P_ildlt;
			P_ildlt = 0;
		}
		P_ildlt = new gmm::ildlt_precond<gmm::csc_matrix<double>>(*M);
		//cout << "done!" << endl;
	}
	// ILDLTT preconditioner
	if(system.P == ILDLTT)
	{
		//cout << "  * Calculating ILDLTT preconditioner ... ";
		if(P_ildltt != 0)
		{
			delete P_ildltt;
			P_ildltt = 0;
		}
		P_ildltt = new gmm::ildltt_precond<gmm::csc_matrix<double>>(*M,fill_in,threshold);
		//cout << "done!" << endl;
	}
	// INV preconditioner
	if(system.P == INV)
	{
		//cout << "  * Calculating APPROX_INV preconditioner ... ";
		if(P_inv != 0)
		{
			delete P_inv;
			P_inv = 0;
		}
		P_inv = new gmm::mr_approx_inverse_precond<gmm::csc_matrix<double>>(*M,fill_in,threshold);
		//cout << "done!" << endl;
	}
	// JACOBI preconditioner
	if(system.P == JACOBI)
	{
		//cout << "  * Calculating JACOBI preconditioner ... ";
		if(P_diag != 0)
		{
			delete P_diag;
			P_diag = 0;
		}
		P_diag = new gmm::diagonal_precond<gmm::csc_matrix<double>>(*M);
		//cout << "done!" << endl;
	}
}
// ------------------------------------------------------------------------- //
void PCG::clear()
{
	if(M != 0)
	{
		delete M;
		M = 0;
	}
	// delete the preconditioner
	if(P_ildlt != 0)
	{
		delete P_ildlt;
		P_ildlt = 0;
	}
	if(P_ildltt != 0)
	{
		delete P_ildltt;
		P_ildltt = 0;
	}
	if(P_diag != 0)
	{
		delete P_diag;
		P_diag = 0;
	}
}
// ------------------------------------------------------------------------- //
void PCG::initialize_x()
{
	gmm::clear(x);
}
// ------------------------------------------------------------------------- //
int PCG::solve()
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


	// Conjugate Gradient , with ILDLT preconditioner
	if(system.P == ILDLT && P_ildlt != 0)
	{
		gmm::cg((*M), x, (*b), PS, (*P_ildlt), iter);
	}
	// Conjugate Gradient , with ILDLT preconditioner
	else if(system.P == ILDLTT && P_ildltt != 0)
	{
		// solve with preconditioner
		gmm::cg((*M), x, (*b), PS, (*P_ildltt), iter);
	}
	// Conjugate Gradient , with INV preconditioner
	else if(system.P == INV && P_inv != 0)
	{
		gmm::cg((*M), x, (*b), PS, (*P_inv), iter);
	}
	// Conjugate Gradient , with JACOBI preconditioner
	else if(system.P == JACOBI && P_diag != 0)
	{
		gmm::cg((*M), x, (*b), PS, (*P_diag), iter);
	}
	else
	{
		cout << "\nPCG::solve() Warning: solving with no PC." << endl;
		gmm::cg((*M), x, (*b), PS, PS, iter);
	}

	if(!system.quiet)
	{
		// reset cout to fixed
		cout.setf(std::ios::fixed,std::ios::floatfield);  
	}
	
	return iter.get_iteration();
}
// ------------------------------------------------------------------------- //
double PCG::residual()
{

	std::vector<double> res(size);
	std::vector<double> Ax(size);
	gmm::mult((*M), x, Ax);	// M1 * V2 --> V1
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
void PCG::write()
{
	MatrixMarket_save("./M.txt", *M);
}
// ------------------------------------------------------------------------- //