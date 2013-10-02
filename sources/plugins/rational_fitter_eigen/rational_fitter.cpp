#include "rational_fitter.h"

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <core/common.h>

ALTA_DLL_EXPORT fitter* provide_fitter()
{
	return new rational_fitter_eigen();
}

rational_fitter_eigen::rational_fitter_eigen() 
{
}
rational_fitter_eigen::~rational_fitter_eigen() 
{
}

bool rational_fitter_eigen::fit_data(const data* dat, function* fit, const arguments &args)
{
	rational_function* r = dynamic_cast<rational_function*>(fit) ;
	const vertical_segment* d = dynamic_cast<const vertical_segment*>(dat) ;
	if(r == NULL || d == NULL)
	{
		std::cerr << "<<ERROR>> not passing the correct class to the fitter" << std::endl ;
		return false ;
	}

	// I need to set the dimension of the resulting function to be equal
	// to the dimension of my fitting problem
	r->setDimX(d->dimX()) ;
	r->setDimY(d->dimY()) ;
	r->setMin(d->min()) ;
	r->setMax(d->max()) ;

	std::cout << "<<INFO>> np =" << _np << "& nq =" << _nq  << std::endl ;
	r->setSize(_np, _nq);

    timer time ;
	time.start() ;


	if(fit_data(d, _np, _nq, r))
	{
        time.stop();
		std::cout << "<<INFO>> got a fit" << std::endl ;
        std::cout << "<<INFO>> it took " << time << std::endl ;

		return true ;
	}

	std::cout << "<<INFO>> fit failed\r"  ;
	std::cout.flush() ;

	return false ;
}

void rational_fitter_eigen::set_parameters(const arguments& args)
{
	_np = args.get_float("np", 10) ;
	_nq = args.get_float("nq", 10) ;
}
		
bool rational_fitter_eigen::fit_data(const vertical_segment* d, int np, int nq, rational_function* r) 
{
    // For each output dimension (color channel for BRDFs) perform
    // a separate fit on the y-1D rational function.
    for(int j=0; j<d->dimY(); ++j)
    {
        rational_function_1d* rs = r->get(j);
        if(!fit_data(d, np, nq, j, rs))
        {
            return false ;
        }
    }

    return true ;
}

// dat is the data object, it contains all the points to fit
// np and nq are the degree of the RP to fit to the data
// y is the dimension to fit on the y-data (e.g. R, G or B for RGB signals)
// the function return a ration BRDF function and a boolean
bool rational_fitter_eigen::fit_data(const vertical_segment* d, int np, int nq, int ny, rational_function_1d* r)
{
	// Each constraint (fitting interval or point
	// add another dimension to the constraint
	// matrix
	Eigen::MatrixXd CI(np+nq, d->size()) ;
	for(int i=0; i<d->size(); ++i)	
	{		
		// A row of the constraint matrix has this 
		// form: [p_{0}(x_i), .., p_{np}(x_i), -f(x_i) q_{0}(x_i), .., -f(x_i) q_{nq}(x_i)]
		for(int j=0; j<np+nq; ++j)
		{
			// Filling the p part
			if(j<np)
			{
				const double pi = r->p(d->get(i), j) ;
				CI(j, i) =  pi ;
			}
			// Filling the q part
			else
			{
				vec yl, yu ; 
				d->get(i, yl, yu) ;

				const double y  = 0.5*(yl[ny] + yu[ny]) ;
				const double qi = r->q(d->get(i), j-np) ;
				CI(j, i) = -y * qi ;
			}
		}
	}

//	std::cout << CI << std::endl << std::endl ;

	// Perform the Eigen decomposition of CI CI'
	Eigen::MatrixXd M(np+nq, np+nq) ;
	M.noalias() = CI * CI.transpose();
  // Faster alternative for large np+nq (only the lower triangular half gets computed):
  // M.setZero();
  // M.selfadjointView<Lower>().rankUpdate(CI.transpose());
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(M) ;
/*
	std::cout << M.size() << std::endl << std::endl ;
	std::cout << M << std::endl << std::endl ;
*/
	if(solver.info() == Eigen::Success)
	{
/*
		std::cout << solver.eigenvalues() << std::endl ;
*/
		// Calculate the minimum eigen value and
		// its position
#ifdef NOT_WORKING
		int    min_id = 0;
		double min_val = solver.eigenvalues().minCoeff(&min_id);
#else
		int    min_id = 0 ;
		double min_val = std::numeric_limits<double>::max() ;
		for(int i=0; i<solver.eigenvalues().size(); ++i)
		{
			double value = solver.eigenvalues()[i] ;
			if(value >= 0 && value < min_val)
			{
				min_id  = i ;
				min_val = value ;
			}
		}
		
		if(min_val == std::numeric_limits<double>::max())
		{
    		return false;
		}
#endif
		// Recopy the vector d
		vec p(np), q(nq);
		Eigen::VectorXd::Map(&p[0], np) = solver.eigenvectors().col(min_id).head(np);
		Eigen::VectorXd::Map(&q[0], nq) = solver.eigenvectors().col(min_id).tail(nq);
		
		r->update(p, q) ;
		std::cout << "<<INFO>> got solution " << *r << std::endl ;
		return true;
	}
	else
	{
		return false ;
	}

}
