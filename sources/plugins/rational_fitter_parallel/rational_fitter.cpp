#include "rational_fitter.h"

#include <core/plugins_manager.h>

#include <Eigen/SVD>
#include <Array.hh>
#include <QuadProg++.hh>

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <QTime>


data* rational_fitter_parallel::provide_data() const
{
	return new vertical_segment() ;
}

function* rational_fitter_parallel::provide_function() const 
{
	return new rational_function() ;
}

rational_fitter_parallel::rational_fitter_parallel() : QObject()
{
}
rational_fitter_parallel::~rational_fitter_parallel() 
{
}

bool rational_fitter_parallel::fit_data(const data* dat, function* fit)
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

	std::cout << "<<INFO>> N in  [" << _min_np << ", " << _max_np << "]"  << std::endl ;

	for(int i=_min_np; i<=_max_np; ++i)
	{
		std::cout << "<<INFO>> fit using np+nq = " << i << "\r"  ;
		std::cout.flush() ;
		QTime time ;
		time.start() ;

		// Allocate enough processor to run fully in parallel, but account for
		// the fact that each thread will require its own memory.
		size_t need_memory = ((3*i+1)*d->size()*2 + 2*i + 2*d->dimY()*d->size()) * sizeof(double); 
		size_t avai_memory = plugins_manager::get_system_memory();
		int nb_cores = (60 * avai_memory) / (100 * need_memory) ;
		nb_cores = std::min<int>(nb_cores, omp_get_num_procs());

		if(nb_cores == 0)
		{
			std::cerr << "<<ERROR>> not enough memory to perform the fit" << std::endl ;
#ifdef DEBUG
			std::cout << "<<DEBUG>> " << need_memory / (1024*1024) << "MB required / " 
			          << avai_memory / (1024*1024) << "MB available" << std::endl;
#endif
			return false;
		}
#ifdef DEBUG
		std::cout << "<<DEBUG>> will use " << nb_cores << " threads to compute the quadratic programs" << std::endl ;
#endif
		omp_set_num_threads(nb_cores) ;

		double min_delta = std::numeric_limits<double>::max();
		int nb_sol_found = 0;
		int np, nq ;
		#pragma omp parallel for
		for(int j=1; j<i; ++j)
		{
			int temp_np = i - j;
			int temp_nq = j;

			vec p(temp_np*r->dimY()), q(temp_nq*r->dimY());

			double delta;
			bool is_fitted = fit_data(d, temp_np, temp_nq, r, p, q, delta);
			if(is_fitted)
			{
				#pragma omp critical
				{
					++nb_sol_found ;
					if(delta < min_delta)
					{
						min_delta = delta ;
#ifdef DEBUG
						std::cout << p << std::endl ;
						std::cout << q << std::endl ;
#endif
						r->update(p, q);
						np = p.size() / d->dimY();
						nq = q.size() / d->dimY();
					}
				}
			}
		}
				
		if(min_delta < std::numeric_limits<double>::max())
		{
			int msec = time.elapsed() ;
			int sec  = (msec / 1000) % 60 ;
			int min  = (msec / 60000) % 60 ;
			int hour = (msec / 3600000) ;
			std::cout << "<<INFO>> got a fit using N = " << i << std::endl ;
			std::cout << "<<INFO>> got a fit using {np, nq} = " << np << ", " << nq << std::endl ;
			std::cout << "<<INFO>> it took " << hour << "h " << min << "m " << sec << "s" << std::endl ;
			std::cout << "<<INFO>> I got " << nb_sol_found << " solutions to the QP" << std::endl ;
			return true ;
		}
	}

	return false ;
}

void rational_fitter_parallel::set_parameters(const arguments& args)
{
	_max_np = args.get_float("np", 10) ;
	_max_nq = args.get_float("nq", 10) ;
	_min_np = args.get_float("min-np", _max_np) ;
	_min_nq = args.get_float("min-nq", _max_nq) ;
}
		

bool rational_fitter_parallel::fit_data(const vertical_segment* d, int np, int nq, rational_function* r, vec& P, vec& Q, double& delta) 
{
	for(int j=0; j<d->dimY(); ++j)
	{
		vec p(np), q(nq);
		if(!fit_data(d, np, nq, j, r, p, q, delta))
			return false ;

		for(int i=0; i<np; ++i) { P[j*np + i] = p[i] ; }
		for(int i=0; i<nq; ++i) { Q[j*nq + i] = q[i] ; }
	}

	return true ;
}

// dat is the data object, it contains all the points to fit
// np and nq are the degree of the RP to fit to the data
// y is the dimension to fit on the y-data (e.g. R, G or B for RGB signals)
// the function return a ration BRDF function and a boolean
bool rational_fitter_parallel::fit_data(const vertical_segment* dat, int np, int nq, int ny, rational_function* rf, vec& p, vec& q, double& delta) 
{
	rational_function* r = dynamic_cast<rational_function*>(rf) ;
	const vertical_segment* d = dynamic_cast<const vertical_segment*>(dat) ;
	if(r == NULL || d == NULL)
	{
		std::cerr << "<<ERROR>> not passing the correct class to the fitter" << std::endl ;
		return false ;
	}

	// Get the maximum value in data to scale the input parameter space
	// so that it reduces the values of the polynomial
	vec dmax = d->max() ;

	// Matrices of the problem
	QuadProgPP::Matrix<double> G (0.0, np+nq, np+nq) ;
	QuadProgPP::Vector<double> g (0.0, np+nq) ;
	QuadProgPP::Matrix<double> CI(0.0, np+nq, 2*d->size()) ;
	QuadProgPP::Vector<double> ci(0.0, 2*d->size()) ;
	QuadProgPP::Matrix<double> CE(0.0, np+nq, 0) ;
	QuadProgPP::Vector<double> ce(0.0, 0) ;

	Eigen::MatrixXd eCI(np+nq, 2*d->size()) ;

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

		vec xi = d->get(i) ;

		// A row of the constraint matrix has this 
		// form: [p_{0}(x_i), .., p_{np}(x_i), -f(x_i) q_{0}(x_i), .., -f(x_i) q_{nq}(x_i)]
		// For the lower constraint and negated for 
		// the upper constraint
		for(int j=0; j<np+nq; ++j)
		{
			// Filling the p part
			if(j<np)
			{
				const double pi = r->p(xi, j) ;
				a0_norm += pi*pi ;
				a1_norm += pi*pi ;
				CI[j][i] =  pi ;
				CI[j][i+d->size()] = -pi ;

				// Updating Eigen matrix
				eCI(j,i) = pi ;
				eCI(j,i+d->size()) = -pi ;
			}
			// Filling the q part
			else
			{
				vec yl, yu ; 
				d->get(i, yl, yu) ;

				const double qi = r->q(xi, j-np) ;
				a0_norm += qi*qi * (yl[ny]*yl[ny]) ;
				CI[j][i] = -yl[ny] * qi ;
				
				a1_norm += qi*qi * (yu[ny]*yu[ny]) ;
				CI[j][i+d->size()] = yu[ny] * qi ;
				
				// Updating Eigen matrix
				eCI(j,i) = -yl[ny] * qi ;
				eCI(j,i+d->size()) = yu[ny] * qi ;
			}
		}
	
		// Set the c vector, will later be updated using the
		// delta parameter.
		ci[i] = -sqrt(a0_norm) ;
		ci[i+d->size()] = -sqrt(a1_norm) ;
	}
#ifdef DEBUG
	std::cout << "CI = [" ;
	for(int j=0; j<d->size()*2; ++j)
	{
		for(int i=0; i<np+nq; ++i)
		{
			std::cout << CI[i][j] ;
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
	Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::HouseholderQRPreconditioner> svd(eCI);
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
	
	delta = sigma_M / sigma_m ;
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
		ci[i] = ci[i] / delta ; 
#ifdef DEBUG
		std::cout << ci[i] << "\t" ;
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
		solves_qp = solves_qp && !std::isnan(v) && (v != std::numeric_limits<double>::infinity()) ;
	}

	if(solves_qp)
	{
		// Recopy the vector d
		double norm = 0.0 ;
		for(int i=0; i<np+nq; ++i)
		{
			const double v = x[i];
			norm += v*v ;
			if(i < np)
			{
				p[i] = v ;
			}
			else
			{
				q[i-np] = v ;
			}
		}

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

Q_EXPORT_PLUGIN2(rational_fitter_parallel, rational_fitter_parallel)
