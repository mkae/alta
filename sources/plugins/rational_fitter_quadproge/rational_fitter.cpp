#include "rational_fitter.h"

#include <Eigen/Dense>
#include <Eigen/SVD>
#include "eiquadprog.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <QTime>


data* rational_fitter_quadprog::provide_data() const
{
	return new vertical_segment() ;
}

function* rational_fitter_quadprog::provide_function() const 
{
	return new rational_function() ;
}

rational_fitter_quadprog::rational_fitter_quadprog() : QObject()
{
}
rational_fitter_quadprog::~rational_fitter_quadprog() 
{
}

bool rational_fitter_quadprog::fit_data(const data* dat, function* fit)
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
		QTime time ;
		time.start() ;
		
		if(fit_data(d, temp_np, temp_nq, r))
		{
			int msec = time.elapsed() ;
			int sec  = (msec / 1000) % 60 ;
			int min  = (msec / 60000) % 60 ;
			int hour = (msec / 3600000) ;
			std::cout << "<<INFO>> got a fit using np = " << temp_np << " & nq =  " << temp_nq << "      " << std::endl ;
			std::cout << "<<INFO>> it took " << hour << "h " << min << "m " << sec << "s" << std::endl ;

			return true ;
		}

		std::cout << "<<INFO>> fit using np = " << temp_np << " & nq =  " << temp_nq << " failed\r"  ;
		std::cout.flush() ;

		if(temp_np <= _max_np)
		{
			++temp_np ;
		}
		if(temp_nq <= _max_nq)
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
}
		

bool rational_fitter_quadprog::fit_data(const vertical_segment* d, int np, int nq, rational_function* r) 
{

	// Multidimensional coefficients
	std::vector<double> Pn ; Pn.reserve(d->dimY()*np) ;
	std::vector<double> Qn ; Qn.reserve(d->dimY()*nq) ;

	for(int j=0; j<d->dimY(); ++j)
	{
		if(!fit_data(d, np, nq, j, r))
			return false ;

		for(int i=0; i<np; ++i) { Pn.push_back(r->getP(i)) ; }
		for(int i=0; i<nq; ++i) { Qn.push_back(r->getQ(i)) ; }
	}

	r->update(Pn, Qn) ;
	return true ;
}

// dat is the data object, it contains all the points to fit
// np and nq are the degree of the RP to fit to the data
// y is the dimension to fit on the y-data (e.g. R, G or B for RGB signals)
// the function return a ration BRDF function and a boolean
bool rational_fitter_quadprog::fit_data(const vertical_segment* dat, int np, int nq, int ny, rational_function* rf) 
{
	rational_function* r = dynamic_cast<rational_function*>(rf) ;
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

	// Get the maximum value in data to scale the input parameter space
	// so that it reduces the values of the polynomial
	vec dmax = d->max() ;

	// Matrices of the problem
	Eigen::MatrixXd G (np+nq, np+nq) ;
	Eigen::VectorXd g (np+nq) ;
	Eigen::MatrixXd CI(np+nq, 2*d->size()) ;
	Eigen::VectorXd ci(2*d->size()) ;
	Eigen::MatrixXd CE(np+nq, 2*d->size()) ;
	Eigen::VectorXd ce(2*d->size()) ;

	// Select the size of the result vector to
	// be equal to the dimension of p + q
	for(int i=0; i<np+nq; ++i)
	{
		for(int j=0; j<np+nq; ++j)
		{
			if(i == j)
			{
				G(i,j) = 1.0 ; 
			}
			else
			{
				G(i,j) = 0.0 ;
			}
		}
		g(i) = 0.0 ;
	}
	
	// Each constraint (fitting interval or point
	// add another dimension to the constraint
	// matrix
	for(int i=0; i<d->size(); ++i)	
	{		
		// Norm of the row vector
		double a0_norm = 0.0 ;
		double a1_norm = 0.0 ;

		vec xi = d->get(i) ;
/*		for(int k=0; k<d->dimX(); ++k)
		{
			xi[k] /= dmax[k] ;
		}
*/
		// A row of the constraint matrix has this 
		// form: [p_{0}(x_i), .., p_{np}(x_i), -f(x_i) q_{0}(x_i), .., -f(x_i) q_{nq}(x_i)]
		// For the lower constraint and negated for 
		// the upper constraint
		for(int j=0; j<np+nq; ++j)
		{
			CE(j, i) = 0.0 ;
			CE(j, i+d->size()) = 0.0 ;

			// Filling the p part
			if(j<np)
			{
				const double pi = r->p(xi, j) ;
				a0_norm += pi*pi ;
				a1_norm += pi*pi ;
				CI(i,j) =  pi ;
				CI(i+d->size(),j) = -pi ;
			}
			// Filling the q part
			else
			{
				vec yl, yu ; 
				d->get(i, yl, yu) ;

				const double qi = r->q(xi, j-np) ;
				a0_norm += qi*qi * (yl[ny]*yl[ny]) ;
				CI(j, i) = -yl[ny] * qi ;
				
				a1_norm += qi*qi * (yu[ny]*yu[ny]) ;
				CI(j, i+d->size()) = yu[ny] * qi ;
			}
		}
	
		// Set the c vector, will later be updated using the
		// delta parameter.
		ci(i)           = 0.0 ; // sqrt(a0_norm) ;
		ci(i+d->size()) = 0.0 ; // sqrt(a1_norm) ;
		ce(i)           = 0.0 ;
		ce(i+d->size()) = 0.0 ;
	
		CI.row(i)           /= sqrt(a0_norm) ;
		CI.row(i+d->size()) /= sqrt(a1_norm) ;
	}
#ifdef DEBUG
	std::cout << "CI = [" ;
	for(int j=0; j<d->size()*2; ++j)
	{
		for(int i=0; i<np+nq; ++i)
		{
			std::cout << CI(i,j) ;
			if(i != np+nq-1) std::cout << ", ";
		}
		if(j != d->size()*2-1)
			std::cout << ";" << std::endl; 
		else
			std::cout << "]" << std::endl ;
	}
#endif
	// Update the ci column with the delta parameter
	// (See Celis et al. 2007 p.12)
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(CI);
	const double sigma_m = svd.singularValues()(std::min(2*d->size(), np+nq)-1) ;
	const double sigma_M = svd.singularValues()(0) ;

#ifdef DEBUG
	std::cout << "<<DEBUG>> SVD = [ " ;
	for(int i=0; i<std::min(2*d->size(), np+nq); ++i)
	{
		std::cout << svd.singularValues()(i) << ", " ;
	}
	std::cout << " ]" << std::endl ;
#endif
	
	double delta = sigma_m / sigma_M ;
	if(std::isnan(delta) || (std::abs(delta) == std::numeric_limits<double>::infinity()))
	{
#ifdef DEBUG
		std::cerr << "<<ERROR>> delta factor is NaN of Inf" << std::endl ;
#endif
		return false ;
	}
	else if(delta == 0.0)
	{
		delta = 1.0 ;
	}


#ifndef DEBUG
	std::cout << "<<DEBUG>> delta factor: " << sigma_m << " / " << sigma_M << " = " << delta << std::endl ;
#endif
	for(int i=0; i<2*d->size(); ++i)	
	{
//		CI.row(i) /= ci(i) ;
		ci(i) = -delta ; 
	}

	// Compute the solution
	Eigen::VectorXd x(np+nq);
	double cost = solve_quadprog(G, g, CE, ce, CI, ci, x);


	bool solves_qp = !(cost == std::numeric_limits<double>::infinity());
	for(int i=0; i<np+nq; ++i)
	{
		const double v = x[i];
		solves_qp = solves_qp && !std::isnan(v) && (v != std::numeric_limits<double>::infinity()) ;
	}

	if(solves_qp)
	{
		// Recopy the vector d
		std::vector<double> p, q;
		double norm = 0.0 ;
		for(int i=0; i<np+nq; ++i)
		{
			const double v = x[i];
			norm += v*v ;
			if(i < np)
			{
				p.push_back(v) ;
			}
			else
			{
				q.push_back(v) ;
			}
		}

		r->update(p, q);
		std::cout << "<<INFO>> got solution " << *r << std::endl ;
		
		return norm > 0.0;
	}
	else

	{
		return false; 
	}
}

Q_EXPORT_PLUGIN2(rational_fitter_quadprog, rational_fitter_quadprog)
