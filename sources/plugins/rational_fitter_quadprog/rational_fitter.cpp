#include "rational_fitter.h"

#include <Eigen/SVD>
#include <Array.hh>
#include <QuadProg++.hh>

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#ifdef WIN32
#define isnan(X) ((X != X))
#endif

#include <core/common.h>

ALTA_DLL_EXPORT fitter* provide_fitter()
{
	return new rational_fitter_quadprog();
}

rational_fitter_quadprog::rational_fitter_quadprog() : _boundary(1.0)
{
}
rational_fitter_quadprog::~rational_fitter_quadprog() 
{
}

bool rational_fitter_quadprog::fit_data(const data* dat, function* fit, const arguments &args)
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

	std::cout << "<<INFO>> np in  [" << _min_np << ", " << _max_np 
	          << "] & nq in [" << _min_nq << ", " << _max_nq << "]" << std::endl ;

	int temp_np = _min_np, temp_nq = _min_nq ;
	while(temp_np <= _max_np || temp_nq <= _max_nq)
	{
        timer time ;
		time.start() ;
		
		r->setSize(temp_np, temp_nq);
		if(fit_data(d, temp_np, temp_nq, r))
		{
            time.stop() ;
            std::cout << "<<INFO>> got a fit using np = " << temp_np << " & nq =  " << temp_nq << "      " << std::endl ;
            std::cout << "<<INFO>> it took " << time << std::endl ;

			return true ;
		}

		
		std::cout << "<<INFO>> fit using np = " << temp_np << " & nq =  " << temp_nq << " failed" << std::endl  ;
		std::cout.flush() ;

		if(temp_np == _max_np && temp_nq == _max_nq)
		{
			return false;
		}

      if(temp_np < _max_np)
		{
			++temp_np ;
		}
      if(temp_nq < _max_nq)
		{
			++temp_nq ;
		}
	}
	return false ;
}

void rational_fitter_quadprog::set_parameters(const arguments& args)
{
	_max_np = args.get_float("np", 10) ;
	_max_nq = args.get_float("nq", 10) ;
	_min_np = args.get_float("min-np", _max_np) ;
	_min_nq = args.get_float("min-nq", _max_nq) ;

	_max_np = std::max<int>(_max_np, _min_np);
	_max_nq = std::max<int>(_max_nq, _min_nq);

	_boundary = args.get_float("boundary-constraint", 1.0f);
}
		

bool rational_fitter_quadprog::fit_data(const vertical_segment* d, int np, int nq, rational_function* r)
{
    // For each output dimension (color channel for BRDFs) perform
    // a separate fit on the y-1D rational function.
	for(int j=0; j<d->dimY(); ++j)
	{
		rational_function_1d* rs = r->get(j);
		rs->resize(np, nq);

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
bool rational_fitter_quadprog::fit_data(const vertical_segment* d, int np, int nq, int ny, rational_function_1d* r) 
{
	// Size of the problem
	const int N = np+nq ;
	const int M = d->size() ;

	// Matrices of the problem
	QuadProgPP::Matrix<double> G (0.0, N, N) ;
	QuadProgPP::Vector<double> g (0.0, N) ;
	QuadProgPP::Matrix<double> CI(0.0, N, 2*M) ;
	QuadProgPP::Vector<double> ci(0.0, 2*M) ;
	QuadProgPP::Matrix<double> CE(0.0, N, 0) ;
	QuadProgPP::Vector<double> ce(0.0, 0) ;

	Eigen::MatrixXd eCI(2*M, N) ;

	// Select the size of the result vector to
	// be equal to the dimension of p + q
	for(int i=0; i<N; ++i)
	{
		G[i][i] = 1.0 ; 
	}
	
	// Each constraint (fitting interval or point
	// add another dimension to the constraint
	// matrix
	for(int i=0; i<M; ++i)	
	{		
		// Norm of the row vector
		double a0_norm = 0.0 ;
		double a1_norm = 0.0 ;

        vec xi = d->get(i) ;

		// A row of the constraint matrix has this 
		// form: [p_{0}(x_i), .., p_{np}(x_i), -f(x_i) q_{0}(x_i), .., -f(x_i) q_{nq}(x_i)]
		// For the lower constraint and negated for 
		// the upper constraint
		for(int j=0; j<N; ++j)
		{
			// Filling the p part
			if(j<np)
			{
				const double pi = r->p(xi, j);

				CI[j][2*i+0] =  pi;
				CI[j][2*i+1] = -pi;

				// Updating Eigen matrix
				eCI(2*i+0, j) = CI[j][2*i+0];
				eCI(2*i+1, j) = CI[j][2*i+1];
			}
			// Filling the q part
			else
			{
				vec yl, yu ; 
				d->get(i, yl, yu) ;

				const double qi = r->q(xi, j-np);

				CI[j][2*i+0] = -yu[ny] * qi;
				CI[j][2*i+1] =  yl[ny] * qi;
				
				// Updating Eigen matrix
				eCI(2*i+0, j) = CI[j][2*i+0];
				eCI(2*i+1, j) = CI[j][2*i+1];
			}

			// Update the norm of the row
			a0_norm += CI[j][2*i+0]*CI[j][2*i+0];
			a1_norm += CI[j][2*i+1]*CI[j][2*i+1];
		}
	
		// Set the c vector, will later be updated using the
		// delta parameter.
		ci[2*i+0] = -sqrt(a0_norm) ;
		ci[2*i+1] = -sqrt(a1_norm) ;
	}
#ifdef DEBUG
	std::cout << "CI = [" ;
	for(int j=0; j<2*M; ++j)
	{
		for(int i=0; i<N; ++i)
		{
			std::cout << CI[i][j] ;
			if(i != N-1) std::cout << ", ";
		}
		if(j != 2*M-1)
			std::cout << ";" << std::endl; 
		else
			std::cout << "]" << std::endl ;
	}
#endif
	
	// Update the ci column with the delta parameter
	// (See Celis et al. 2007 p.12)
	Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::HouseholderQRPreconditioner> svd(eCI,  Eigen::ComputeThinU | Eigen::ComputeThinV);
	const double sigma_m = svd.singularValues()(std::min(2*M, N)-1) ;
	const double sigma_M = svd.singularValues()(0) ;

#ifdef DEBUG
	std::cout << "<<DEBUG>> SVD = [ " ;
	for(int i=0; i<std::min(2*M, N); ++i)
	{
		std::cout << svd.singularValues()(i) << ", " ;
	}
	std::cout << " ]" << std::endl ;
#endif
	
	double delta = sigma_m / sigma_M ;

	if(isnan(delta) || (std::abs(delta) == std::numeric_limits<double>::infinity()))
	{
		std::cerr << "<<ERROR>> delta factor is NaN of Inf" << std::endl ;
		return false ;
	}
   else if(delta <= 0.0)
	{
		delta = 1.0 ;
	}

#ifndef DEBUG
	std::cout << "<<DEBUG>> delta factor: " << sigma_m << " / " << sigma_M << " = " << delta << std::endl ;
#endif
	for(int i=0; i<2*M; ++i)	
	{		
		ci[i] = ci[i] * delta ; 
#ifdef DEBUG
        std::cout << i << "\t" << -ci[i] << std::endl ;
#endif
	}
#ifdef DEBUG
	std::cout << std::endl << std::endl ;

	std::cout << eCI << std::endl << std::endl ;
#endif
	// Compute the solution
	QuadProgPP::Vector<double> x;
	double cost = QuadProgPP::solve_quadprog(G, g, CE, ce, CI, ci, x);


	bool solves_qp = !(cost == std::numeric_limits<double>::infinity());
	for(int i=0; i<np+nq; ++i)
	{
		const double v = x[i];
		solves_qp = solves_qp && !isnan(v) && (v != std::numeric_limits<double>::infinity()) ;
	}

	if(solves_qp)
	{
		// Recopy the vector d
		vec p(np), q(nq);
		double norm = 0.0 ;
		for(int i=0; i<N; ++i)
		{
			const double v = x[i];
			norm += v*v ;
			if(i < np)
			{
				p[i] = v ;
			}
			else
			{
				q[i - np] = v ;
			}
		}

		r->update(p, q);
#ifdef DEBUG
		std::cout << "<<INFO>> got solution " << *r << std::endl ;
#endif
		return norm > 0.0;
	}
	else

	{
		return false; 
	}
}
