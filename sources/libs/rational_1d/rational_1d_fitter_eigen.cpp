#include "rational_1d_fitter_eigen.h"

#include <Eigen/Dense>
#include <Eigen/SVD>
#include "eiquadprog.hpp"

#include <iostream>
#include <cfenv>

bool rational_1d_fitter_eigen::fit_data(const rational_1d_data& data, int np, int nq, rational_1d& r) 
{
	// Select the size of the result vector to be equal to the dimension 
	// of p + q.
	//
	// Note: since the QP has the following form in the Eigen code:
	//     f(x) = 0.5 * x' G x
	// I have to multiply my matrix by a factor 2
	Eigen::MatrixXd G(np+nq, np+nq) ;
	for(int i=0; i<np+nq; ++i)
	{
		G(i,i) = 2.0 ; 
	}
	
	// Each constraint (fitting interval or point
	// add another dimension to the constraint
	// matrix
	Eigen::MatrixXd CI(np+nq, 2*data.size()) ;
	Eigen::VectorXd ci(2*data.size()) ;
	for(int i=0; i<data.size(); ++i)	
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
				const double pi = r.p(data[i][0], j) ;
				a0_norm += pi*pi ;
				a1_norm += pi*pi ;
				CI(j, i) =   pi ;
				CI(j, i+data.size()) = - pi ;
			}
			// Filling the q part
			else
			{
				const double qi = r.q(data[i][0], j-np) ;
				a0_norm += qi*qi * (data[i][1]*data[i][1]) ;
				CI(j, i) = - data[i][1] * qi ;
				
				a1_norm += qi*qi * (data[i][2]*data[i][2]) ;
				CI(j, i+data.size()) =   data[i][2] * qi ;
			}
		}
	
		// Set the c vector, will later be updated using the
		// delta parameter.
		ci(i) = -sqrt(a0_norm) ;
		ci(i+data.size()) = -sqrt(a1_norm) ;
	}
	
	// Update the ci column with the delta parameter
	// (See Celis et al. 2007 p.12)
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(CI);
	const double sigma_m = svd.singularValues()(std::min(2*data.size(), np+nq)-1) ;
	const double sigma_M = svd.singularValues()(0) ;

#ifdef DEBUG
	std::cout << "<<DEBUG>> SVD = [ " ;
	for(int i=0; i<std::min(2*data.size(), np+nq); ++i)
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
	std::cout << "<<DEBUG>> delta factor: " << sigma_M << " / " << sigma_m << " = " << delta << std::endl ;
#endif
	for(int i=0; i<2*data.size(); ++i)	
	{		
		ci(i) *= delta ;
	}
	
	// solve the program, using ET as the exact type
	Eigen::VectorXd g0(np+nq) ;
	Eigen::MatrixXd CE(np+nq, 2*data.size()) ;
	Eigen::VectorXd ce(2*data.size()) ;
	Eigen::VectorXd x(np+nq) ;

	double cost = solve_quadprog(G, g0, CE, ce, CI, ci, x) ;
	bool solves_qp = abs(cost) != std::numeric_limits<double>::infinity() /*&& cost > 0*/ ;
	
	// Check the data
	for(int i=0; i<np+nq; ++i)
	{
		solves_qp = solves_qp && !isnan(x(i)) && (x(i) != std::numeric_limits<double>::infinity()) ;
	}

#ifdef DEBUG
	std::cout << "<<DEBUG>> cost: " << cost << std::endl ;
#endif

	if(solves_qp)
	{
		// Recopy the vector data
		std::vector<double> p, q;
		for(int i=0; i<np+nq; ++i)
		{
			const double v = x(i) ;

			if(i < np)
			{
				p.push_back(v) ;
			}
			else
			{
				q.push_back(v) ;
			}
		}

		r = rational_1d(p, q);
		return true ;
	}
	else
	{
		return false; 
	}
}

