#include "rational_1d_fitter_cgal.h"

#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <CGAL/MP_Float.h>
#include <Eigen/SVD>

#include <boost/regex.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>


typedef CGAL::MP_Float ET ;
typedef CGAL::Quadratic_program<ET> Program ;
typedef CGAL::Quadratic_program_solution<ET> Solution ;
typedef CGAL::Quadratic_program_options Options ;
		

bool rational_1d_fitter_cgal::fit_data(const rational_1d_data* data, int np, int nq, rational_1d*& r) 
{
	// by default, we have a nonnegative QP with Ax - b >= 0
	Program qp (CGAL::LARGER, false, 0, false, 0) ; 

	// Select the size of the result vector to
	// be equal to the dimension of p + q
	for(int i=0; i<np+nq; ++i)
	{
		qp.set_d(i, i, 1.0) ; 
	}
	
	// Each constraint (fitting interval or point
	// add another dimension to the constraint
	// matrix
	Eigen::MatrixXd CI(np+nq, 2*data->size()) ;
	Eigen::VectorXd ci(2*data->size()) ;
	for(int i=0; i<data->size(); ++i)	
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
				const double pi = r->p((*data)[i][0], j) ;
				a0_norm += pi*pi ;
				a1_norm += pi*pi ;
				qp.set_a(j, i,  ET(pi)) ;
				qp.set_a(j, i+data->size(), ET(-pi)) ;
				CI(j, i) =  pi ;
				CI(j, i+data->size()) = -pi ;
			}
			// Filling the q part
			else
			{
				const double qi = r->q((*data)[i][0], j-np) ;
				a0_norm += qi*qi * ((*data)[i][1]*(*data)[i][1]) ;
				qp.set_a(j, i, ET(-(*data)[i][1] * qi)) ;
				CI(j, i) = -(*data)[i][1] * qi ;
				
				a1_norm += qi*qi * ((*data)[i][2]*(*data)[i][2]) ;
				qp.set_a(j, i+data->size(),  ET((*data)[i][2] * qi)) ;
				CI(j, i+data->size()) = (*data)[i][2] * qi ;
			}
		}
	
		// Set the c vector, will later be updated using the
		// delta parameter.
		ci(i) = sqrt(a0_norm) ;
		ci(i+data->size()) = sqrt(a1_norm) ;
	}

	// Update the ci column with the delta parameter
	// (See Celis et al. 2007 p.12)
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(CI);
	const double sigma_m = svd.singularValues()(std::min(2*data->size(), np+nq)-1) ;
	const double sigma_M = svd.singularValues()(0) ;

#ifdef DEBUG
	std::cout << "<<DEBUG>> SVD = [ " ;
	for(int i=0; i<std::min(2*data->size(), np+nq); ++i)
	{
		std::cout << svd.singularValues()(i) << ", " ;
	}
	std::cout << " ]" << std::endl ;
#endif
	
	const double delta = sigma_m / sigma_M ;
	if(isnan(delta) || (abs(delta) == std::numeric_limits<double>::infinity()))
	{
#ifdef DEBUG
		std::cerr << "<<ERROR>> delta factor is NaN of Inf" << std::endl ;
#endif
		return false ;
	}

#ifdef DEBUG
	std::cout << "<<DEBUG>> delta factor: " << sigma_m << " / " << sigma_M << " = " << delta << std::endl ;
#endif
	for(int i=0; i<2*data->size(); ++i)	
	{		
		qp.set_b(i, ET(delta * ci(i))) ;
	}

#ifdef DEBUG
	// Export some informations on the problem to solve
	std::cout << "<<DEBUG>> " << qp.get_n() << " variables" << std::endl ;
	std::cout << "<<DEBUG>> " << qp.get_m() << " constraints" << std::endl ;
#endif

	if(qp.is_linear() && !qp.is_nonnegative())
	{
		std::cerr << "<<ERROR>> the problem should not be linear or negative!" << std::endl ;
		return false ;
	}

#ifdef EXTREM_DEBUG
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
#endif

	// solve the program, using ET as the exact type
	Options  o ;
	o.set_auto_validation(true) ;
	Solution s = CGAL::solve_quadratic_program(qp, ET(), o) ;

#ifdef DEBUG
	if(s.is_infeasible())
	{
		std::cout << "<<DEBUG>> the current program is infeasible" << std::endl ;
	}
#endif

	bool solves_qp = !s.is_infeasible() && s.solves_quadratic_program(qp)  ;
	for(int i=0; i<np+nq; ++i)
	{
		const double v = (double)CGAL::to_double(*(s.variable_numerators_begin()+i)) ;
		solves_qp = solves_qp && !isnan(v) && (v != std::numeric_limits<double>::infinity()) ;
	}

	if(solves_qp)
	{
		// Recopy the vector data
		std::vector<double> p, q;
		for(int i=0; i<np+nq; ++i)
		{
			const double v = CGAL::to_double(*(s.variable_numerators_begin()+i)) ;

			if(i < np)
			{
				p.push_back(v) ;
			}
			else
			{
				q.push_back(v) ;
			}
		}

		if(r != nullptr)
		{
			delete r ;
		}
		r = new rational_1d(p, q);
		return true;
	}
	else
	{
		return false; 
	}
}
