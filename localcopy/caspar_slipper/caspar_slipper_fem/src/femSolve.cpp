#include "fem.h"

namespace CasparSlipperFEM
{
	bool fem::solve(int SolSequence)
	{
		for(size_t i=0; i<X.size(); i++)
		{
			X[i] = 0;
		}
		
		//Set up the proper sparse matrix
		gmm::csc_matrix<double> KS;
		gmm::copy(K,KS);

		//Set solver params
		gmm::iteration iter(options.tolerance);
		iter.set_maxiter(options.maxIters);
		iter.set_noisy(0);
				
		//Precondition
		if(SolSequence == 0)
		{
			gmm::ildlt_precond< gmm::csc_matrix<double> > PR(KS);
			gmm::cg(KS, X, b, PR, iter);
		} else if (SolSequence == 1)
		{
			gmm::diagonal_precond< gmm::csc_matrix<double> > PR(KS);
			gmm::cg(KS, X, b, PR, iter);
		} else if (SolSequence == 2)
		{
			gmm::ildlt_precond< gmm::csc_matrix<double> > PR(KS);
			gmm::bicgstab(KS, X, b, PR, iter);
		}

		if(iter.get_iteration() == iter.get_maxiter())
		{
			warning("Linear solver failed to converge in "+n2s(iter.get_maxiter())+" iterations!");

			SolSequence++;

			if(SolSequence < 3)
			{
				cout << endl << "   Relaunching the linear solver with different options .. (SolSequence = " << SolSequence << " )" << endl;
				cout << endl;
				return solve(SolSequence);
			}
		} else {
			cout << "\tSolved linear system in " << iter.get_iteration() << " iterations." << endl;
		}

		if(inertia_relief)
		{
			cout << "\tInertia Relief lambda = [ " << X[X.size()-3] << "\t"
																<< X[X.size()-2] << "\t"
																<< X[X.size()-1] << " ]" << endl;
		}

		if(iter.get_iteration() == iter.get_maxiter())
		{
			return false;
		}
		return true;

	};

};
