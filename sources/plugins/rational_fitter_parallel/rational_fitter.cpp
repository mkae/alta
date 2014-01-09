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
#include <string>

#include "quadratic_program.h"

ALTA_DLL_EXPORT fitter* provide_fitter()
{
	return new rational_fitter_parallel();
}

rational_fitter_parallel::rational_fitter_parallel()
{
}
rational_fitter_parallel::~rational_fitter_parallel() 
{
}

bool rational_fitter_parallel::fit_data(const data* dat, function* fit, const arguments &args)
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

	const int _min_np = args.get_int("min-np", 10);
	const int _max_np = args.get_int("np", _min_np);
	std::cout << "<<INFO>> N in  [" << _min_np << ", " << _max_np << "]"  << std::endl ;

    int k = 1;
    for(int i=_min_np; i<=_max_np; ++i)
	{
		std::cout << "<<INFO>> fit using np+nq = " << i << std::endl ;
		std::cout.flush() ;
		timer time ;
		time.start() ;

#ifdef NOT_WORKING
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
			//return false;
			nb_cores = 1;
		}
#endif
		int nb_cores = args.get_int("nbcores", omp_get_num_procs());

#ifdef DEBUG
		std::cout << "<<DEBUG>> will use " << nb_cores << " threads to compute the quadratic programs" << std::endl ;
#endif

		omp_set_num_threads(nb_cores) ;
		
		std::vector<rational_function*> rs;
		for(int j=0; j<nb_cores; ++j)
		{
			rational_function* rj = dynamic_cast<rational_function*>(plugins_manager::get_function(args));

			rj->setDimX(d->dimX()) ;
			rj->setDimY(d->dimY()) ;
			rj->setMin(d->min()) ;
			rj->setMax(d->max()) ;

			if(rj == NULL)
			{
				std::cerr << "<<ERROR>> unable to obtain a rational function from the plugins manager" << std::endl;
				return false;
			}
			rs.push_back(rj);
		}

		// Solution for the case of optimizing the L2 norm
		double min_l2 = std::numeric_limits<double>::max();
		rational_function* min_l2_fun = NULL;

		// Solution for the case of optimizing the LINF norm
		double min_linf = std::numeric_limits<double>::max();
		rational_function* min_linf_fun = NULL;

		double min_delta  = std::numeric_limits<double>::max();
		double mean_delta = 0.0;
		int nb_sol_found  = 0;
		int nb_sol_tested = 0;

		#pragma omp parallel for
		for(int j=1; j<i; ++j)
		{
			int temp_np = i - j;
			int temp_nq = j;

			vec p(temp_np*r->dimY()), q(temp_nq*r->dimY());
			rs[omp_get_thread_num()]->setSize(temp_np, temp_nq);

			double delta, linf_dist, l2_dist;
			bool is_fitted = fit_data(d, temp_np, temp_nq, rs[omp_get_thread_num()], args, p, q, delta, linf_dist, l2_dist);
			if(is_fitted)
			{
				#pragma omp critical
				{
					++nb_sol_found ;
               mean_delta += delta ;

               std::cout << "<<INFO>> found a solution with np=" << temp_np << ", nq = " << temp_nq << std::endl;
               std::cout << "<<INFO>> Linf error = " << linf_dist << std::endl;
               std::cout << "<<INFO>> L2   error = " << l2_dist << std::endl;
               std::cout << "<<INFO>>      delta = " << delta << std::endl;
               std::cout << std::endl;

               // Get the solution with the minimum delta
					if(delta < min_delta)
					{
						min_delta   = delta ;
						r->setSize(temp_np, temp_nq);
						for(int y=0; y<r->dimY(); ++y)
						{
							r->update(y, rs[omp_get_thread_num()]->get(y));
						}
					}

					// Get the solution with the minumum L2 norm
					if(l2_dist < min_l2)
					{
						min_l2 = l2_dist;
						min_l2_fun = rs[omp_get_thread_num()];
					}

					// Get the solution with the minumum LINF norm
					if(linf_dist < min_linf)
					{
						min_linf = linf_dist;
						min_linf_fun = rs[omp_get_thread_num()];
					}
				}
			}

			#pragma omp critical
			{
				// Update the solution
				nb_sol_tested++;
				std::cout << "<<DEBUG>> nb solutions tested: " << nb_sol_tested << " / " << i << "\r";
				std::cout.flush();
			}
		}

/*
		// Clean memory
		for(int j=0; j<nb_cores; ++j)
		{
			delete rs[j];
		}
*/
		std::cout << "                                                                           \r";
		if(min_l2 < std::numeric_limits<double>::max())
		{
			std::cout << "<<INFO>>  min L2 = " << min_l2 << std::endl;
			std::cout << *min_l2_fun << std::endl;
			std::cout << std::endl;
		}
//		delete min_l2_fun;

		if(min_linf < std::numeric_limits<double>::max())
		{
			std::cout << "<<INFO>>  min Linf = " << min_linf << std::endl;
			std::cout << *min_linf_fun << std::endl;
			std::cout << std::endl;
		}
//		delete min_linf_fun;

		if(min_delta < std::numeric_limits<double>::max())
		{
			std::cout << "<<INFO>> mean delta = " << mean_delta/nb_sol_found << std::endl;
			std::cout << "<<INFO>>  min delta = " << min_delta << std::endl;
			std::cout << *r<< std::endl;

			time.stop();
			std::cout << "<<INFO>> got a fit using N = " << i << std::endl ;
			std::cout << "<<INFO>> it took " << time << std::endl ;
			std::cout << "<<INFO>> I got " << nb_sol_found << " solutions to the QP" << std::endl ;
			return true ;
		}
	}

	return false ;
}

void rational_fitter_parallel::set_parameters(const arguments& args)
{
}


bool rational_fitter_parallel::fit_data(const vertical_segment* d, int np, int nq, 
        rational_function* r, const arguments &args,
        vec& P, vec& Q, double& delta, double& linf_dist, double& l2_dist)
{
	// Fit the different outptu dimension independantly
	for(int j=0; j<d->dimY(); ++j)
	{
		vec p(np), q(nq);
		rational_function_1d* rf = r->get(j);
		rf->resize(np, nq);
		if(!fit_data(d, np, nq, j, rf, p, q, delta))
		{
			return false ;
		}

		rf->update(p, q);
	}

	linf_dist = r->Linf_distance(d);
	l2_dist   = r->L2_distance(d);

	return true ;
}

// dat is the data object, it contains all the points to fit
// np and nq are the degree of the RP to fit to the data
// y is the dimension to fit on the y-data (e.g. R, G or B for RGB signals)
// the function return a ration BRDF function and a boolean
bool rational_fitter_parallel::fit_data(const vertical_segment* d, int np, int nq, int ny,
                                        rational_function_1d* r, vec& p, vec& q, double& delta)
{
	const int m = d->size(); // 2*m = number of constraints
	const int n = np+nq;     // n = np+nq

    quadratic_program qp(np, nq);


    // Starting with only a 100 vertical segments (meaning a 200xn matrix)
    const int di = (m-1) / 99;
    for(int i=0; i<m; i+=di)
	{

		// Create two vector of constraints
		vec c1(n), c2(n);
        get_constraint(i, np, nq, ny, d, r, c1, c2);

		qp.add_constraints(c1);
		qp.add_constraints(c2);
	}

	while(qp.nb_constraints() < 2*m)
	{
#ifdef DEBUG
        std::cout << "<<DEBUG>> thread " << omp_get_thread_num() << ", number of intervals tested = " << qp.nb_constraints()/2 << std::endl ;
#endif
		QuadProgPP::Vector<double> x(n);
		bool solves_qp = qp.solve_program(x, delta, p, q);
		r->update(p, q);

		if(solves_qp)
		{
			if(qp.test_constraints(ny, r, d))
			{
#ifdef DEBUG
				std::cout << "<<INFO>> got solution " << *r << std::endl ;
#endif
				/*
					int current = 0, i=0;
					while(i < 100 && current < m)
					{

					int next = quadratic_program::next_unmatching_constraint(current, ny, );

				// Create two vector of constraints
				vec c1(n), c2(n);
				get_constraint(next, np, nq, ny, d, r, c1, c2);

				qp.add_constraints(c1);
				qp.add_constraints(c2);

				++i;
				current = next;
				}
				*/
				return true;
            }
        }
        else
		{
#ifdef DEBUG
			std::cout << "<<DEBUG>> not enough coefficients" << std::endl;
#endif
			return false;
		}
	}


	return false;
}

void rational_fitter_parallel::get_constraint(int i, int np, int nq, int ny, 
		                                        const vertical_segment* data, 
															 const rational_function_1d* func,
															 vec& cu, vec& cl)
{
	const vec xi = data->get(i) ;
	cu.resize(np+nq);
	cl.resize(np+nq);

	// Create two vector of constraints
	for(int j=0; j<np+nq; ++j)
	{
		// Filling the p part
		if(j<np)
		{
			const double pi = func->p(xi, j) ;
			cu[j] =  pi ;
			cl[j] = -pi ;

		}
		// Filling the q part
		else
		{
			vec yl, yu ;
			data->get(i, yl, yu) ;
			const double qi = func->q(xi, j-np) ;

			cu[j] = -yu[ny] * qi ;
			cl[j] =  yl[ny] * qi ;
		}
	}
}
