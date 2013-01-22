#include "rational_1d_fitter_eigen.h"

#include <Eigen/Dense>
#include <Eigen/SVD>
#include "eiquadprog.hpp"

#include <iostream>
#include <cfenv>

// Fitting a data
bool rational_1d_fitter_eigen::fit_data(const rational_1d_data& data, rational_1d& fit)
{
	return fit_data(data, 10, 10, fit) ;
}

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
		float a0_norm = 0.0f ;
		float a1_norm = 0.0f ;

		// A row of the constraint matrix has this 
		// form: [p_{0}(x_i), .., p_{np}(x_i), -f(x_i) q_{0}(x_i), .., -f(x_i) q_{nq}(x_i)]
		// For the lower constraint and negated for 
		// the upper constraint
		for(int j=0; j<np+nq; ++j)
		{
			// Filling the p part
			if(j<np)
			{
				const float pi = r.p(data[i][0], j) ;
				a0_norm += pi*pi ;
				a1_norm += pi*pi ;
				CI(j, 2*i+0) =  pi ;
				CI(j, 2*i+1) = -pi ;
			}
			// Filling the q part
			else
			{
				const float qi = r.q(data[i][0], j-np) ;
				a0_norm += qi*qi * (data[i][1]*data[i][1]) ;
				CI(j, 2*i+0) = -data[i][1] * qi ;
				
				a1_norm += qi*qi * (data[i][2]*data[i][2]) ;
				CI(j, 2*i+1) = data[i][2] * qi ;
			}
		}
	
		// Set the c vector, will later be updated using the
		// delta parameter.
		ci(2*i+0) = -sqrt(a0_norm) ;
		ci(2*i+1) = -sqrt(a1_norm) ;
	}
	
	// Update the ci column with the delta parameter
	// (See Celis et al. 2007 p.12)
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(CI);
	const double sigma_m = svd.singularValues()(std::min(2*data.size(), np+nq)-1) ;
	const double sigma_M = svd.singularValues()(0) ;
	const double delta = sigma_m / sigma_M ;
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
	bool prog_solved = cost != std::numeric_limits<double>::infinity() ;

#ifdef DEBUG
	std::cout << "<<DEBUG>> cost: " << cost << std::endl ;
#endif

	if(prog_solved)
	{
		// Recopy the vector data
		std::vector<float> p, q;
		for(int i=0; i<np+nq; ++i)
		{
			const float v = x(i) ;

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

