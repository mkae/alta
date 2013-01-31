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

#include <rational_function.h>
#include <rational_data.h>

rational_fitter_quadprog::rational_fitter_quadprog() : QObject()
{
}
rational_fitter_quadprog::~rational_fitter_quadprog() 
{
}

bool rational_fitter_quadprog::fit_data(const data* d, function* fit)
{
	std::cout << "<<INFO>> np in  [" << _min_np << ", " << _max_np << "] & nq in [" << _min_nq << ", " << _max_nq << "]" << std::endl ;
	int temp_np = _min_np, temp_nq = _min_nq ;
	while(temp_np <= _max_np || temp_nq <= _max_nq)
	{
		if(fit_data(d, temp_np, temp_nq, fit))
		{
			return true ;
		}

		std::cout << "<<INFO>> fitt using np = " << temp_np << " & nq =  " << temp_nq << " failed\r"  ;
		std::cout.flush() ;

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
}
		

bool rational_fitter_quadprog::fit_data(const data* dat, int np, int nq, function* rf) 
{
	// by default, we have a nonnegative QP with Ax - b >= 0
//	Program qp (CGAL::LARGER, false, 0, false, 0) ; 

	rational_function* r = dynamic_cast<rational_function*>(rf) ;
	const rational_data* d = dynamic_cast<const rational_data*>(dat) ;
	if(r == NULL || d == NULL)
	{
		std::cerr << "<<ERROR>> not passing the correct class to the fitter" << std::endl ;
		return false ;
	}

	// I need to set the dimension of the resulting function to be equal
	// to the dimension of my fitting problem
	r->setDimX(d->dimX()) ;
	r->setDimY(d->dimY()) ;

	// Matrices of the problem
	QuadProgPP::Matrix<double> G (0.0, np+nq, np+nq) ;
	QuadProgPP::Vector<double> g (0.0, np+nq) ;
	QuadProgPP::Matrix<double> CI(0.0, np+nq, 2*d->size()) ;
	QuadProgPP::Vector<double> ci(0.0, 2*d->size()) ;
	QuadProgPP::Matrix<double> CE(0.0, np+nq, 2*d->size()) ;
	QuadProgPP::Vector<double> ce(0.0, 2*d->size()) ;

	// Select the size of the result vector to
	// be equal to the dimension of p + q
	for(int i=0; i<np+nq; ++i)
	{
		G[i][i] = 1.0 ; 
	}
	
	// Each constraint (fitting interval or point
	// add another dimension to the constraint
	// matrix
	for(int i=0; i<d->size(); ++i)	
	{		
		// Norm of the row vector
		double a0_norm = 0.0 ;
		double a1_norm = 0.0 ;

		// A row of the constraint matrix has this 
		// form: [p_{0}(x_i), .., p_{np}(x_i), -f(x_i) q_{0}(x_i), .., -f(x_i) q_{nq}(x_i)]
		// For the lower constraint and negated for 
		// the upper constraint
		for(int j=0; j<np+nq; ++j)
		{
			// Filling the p part
			if(j<np)
			{
				const double pi = r->p(d->get(i), j) ;
				a0_norm += pi*pi ;
				a1_norm += pi*pi ;
				CI[j][i] =  pi ;
				CI[j][i+d->size()] = -pi ;
			}
			// Filling the q part
			else
			{
				vec yl, yu ; 
				d->get(i, yl, yu) ;

				const double qi = r->q(d->get(i), j-np) ;
				a0_norm += qi*qi * (yl[0]*yl[0]) ;
				CI[j][i] = -yl[0] * qi ;
				
				a1_norm += qi*qi * (yu[0]*yu[0]) ;
				CI[j][i+d->size()] = yu[0] * qi ;
			}
		}
	
		// Set the c vector, will later be updated using the
		// delta parameter.
		ci[i] = -sqrt(a0_norm) ;
		ci[i+d->size()] = -sqrt(a1_norm) ;
	}
#ifdef DEBUG
	for(int j=0; j<d->size()*2; ++j)
	{
		for(int i=0; i<np+nq; ++i)
			std::cout << CI[i][j] << "\t";
		std::cout << std::endl; 
	}
	std::cout << std::endl ;
#endif
	// Update the ci column with the delta parameter
	// (See Celis et al. 2007 p.12)
/*
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(CI);
	const double sigma_m = svd.singularValues()(std::min(2*d->size(), np+nq)-1) ;
	const double sigma_M = svd.singularValues()(0) ;
*/
#ifdef DEBUG
/*
	std::cout << "<<DEBUG>> SVD = [ " ;
	for(int i=0; i<std::min(2*d->size(), np+nq); ++i)
	{
		std::cout << svd.singularValues()(i) << ", " ;
	}
	std::cout << " ]" << std::endl ;
*/
#endif
/*	
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


#ifdef DEBUG
	std::cout << "<<DEBUG>> delta factor: " << sigma_m << " / " << sigma_M << " = " << delta << std::endl ;
#endif
	for(int i=0; i<2*d->size(); ++i)	
	{		
//		qp.set_b(i, ET(delta * ci(i))) ;
	}
*/
#ifdef DEBUG
	// Export some informations on the problem to solve
//	std::cout << "<<DEBUG>> " << qp.get_n() << " variables" << std::endl ;
//	std::cout << "<<DEBUG>> " << qp.get_m() << " constraints" << std::endl ;
#endif

/*
	if(qp.is_linear() && !qp.is_nonnegative())
	{
		std::cerr << "<<ERROR>> the problem should not be linear or negative!" << std::endl ;
		return false ;
	}
*/
#ifdef EXTREM_DEBUG
/*
	// Iterate over the rows
	std::cout << std::endl ;
	std::cout << "A = [" ;
	for(int j=0; j<qp.get_m(); ++j)
	{
		if(j == 0) std::cout << "   " ;

		// Iterate over the columns
		for(int i=0; i<qp.get_n(); ++i)
		{
			if(i > 0) std::cout << ",\t" ;
			std::cout << *(*(qp.get_a()+i) +j) ;
		}
		std::cout << ";" ;
		if(j < qp.get_n()-1) std::cout << std::endl ;
	}
	std::cout << "]" << std::endl << std::endl ;
	
	std::cout << std::endl ;
	std::cout << "D = [" ;
	for(int j=0; j<np+nq; ++j)
	{
		if(j == 0) std::cout << "   " ;

		// Iterate over the columns
		for(int i=0; i<np+nq; ++i)
		{
			if(i > 0) std::cout << ",\t" ;
			std::cout << *(*(qp.get_d()+i) +j) ;
		}
		std::cout << ";" ;
	}
	std::cout << "]" << std::endl << std::endl ;
*/
#endif

	// Compute the solution
	QuadProgPP::Vector<double> x;
	double cost = QuadProgPP::solve_quadprog(G, g, CE, ce, CI, ci, x);


	bool solves_qp = !cost == std::numeric_limits<double>::infinity();
	for(int i=0; i<np+nq; ++i)
	{
		const double v = x[i];
		solves_qp = solves_qp && !std::isnan(v) && (v != std::numeric_limits<double>::infinity()) ;
	}

	if(solves_qp)
	{
		// Recopy the vector d
		std::vector<double> p, q;
		for(int i=0; i<np+nq; ++i)
		{
			const double v = x[i];

			if(i < np)
			{
				p.push_back(v) ;
			}
			else
			{
				q.push_back(v) ;
			}
		}

		if(r != NULL)
		{
			delete r ;
		}
		r = new rational_function(p, q);
			
		std::cout << "<<INFO>> got solution " << *r << std::endl ;
		return true;
	}
	else

	{
		return false; 
	}
}

Q_EXPORT_PLUGIN2(rational_fitter_quadprog, rational_fitter_quadprog)
