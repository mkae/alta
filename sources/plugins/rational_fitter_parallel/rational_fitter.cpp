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
#include <list>
#include <omp.h>


#include "quadratic_program.h"

ALTA_DLL_EXPORT fitter* provide_fitter()
{
	return new rational_fitter_parallel();
}

rational_fitter_parallel::rational_fitter_parallel() : nb_starting_points(100)
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

	const int nb_starting_points = args.get_int("nb-starting-points", 100);
	std::cout << "<<INFO>> number of data point used in start: " << nb_starting_points << std::endl;

	const int step = args.get_int("np-step", 1);

    for(int i=_min_np; i<=_max_np; i+=step)
	{
		std::cout << "<<INFO>> fit using np+nq = " << i << std::endl ;
		std::cout.flush() ;
		timer time ;
		time.start() ;

        int nb_cores = args.get_int("nb-cores", omp_get_num_procs());
#ifdef DEBUG
		std::cout << "<<DEBUG>> will use " << nb_cores << " threads to compute the quadratic programs" << std::endl ;
#endif

		omp_set_num_threads(nb_cores) ;


		double min_delta  = std::numeric_limits<double>::max();
		double mean_delta = 0.0;
		int nb_sol_found  = 0;
		int nb_sol_tested = 0;

        #pragma omp parallel for shared(args, nb_sol_found, nb_sol_tested, min_delta, mean_delta), schedule(dynamic,1)
        for(int j=1; j<i; ++j)
        {
            // Compute the number of coefficients in the numerator and in the denominator
            // from the current number of coefficients i and the current index in the
            // loop j.
            int temp_np = i - j;
            int temp_nq = j;

            vec p(temp_np*r->dimY()), q(temp_nq*r->dimY());

            // Allocate a rational function and set it to the correct size, dimensions
            // and parametrizations.
				rational_function* rk = NULL;
            #pragma omp critical (args)
            {
					rk = dynamic_cast<rational_function*>(plugins_manager::get_function(args));
				}
				if(rk == NULL)
            {
                std::cerr << "<<ERROR>> unable to obtain a rational function from the plugins manager" << std::endl;
                throw;
            }
            rk->setParametrization(r->input_parametrization());
            rk->setParametrization(r->output_parametrization());
            rk->setDimX(r->dimX()) ;
            rk->setDimY(r->dimY()) ;
            rk->setMin(r->min()) ;
            rk->setMax(r->max()) ;

            // Set the rational function size
            rk->setSize(temp_np, temp_nq);

            double delta, linf_dist, l2_dist;
            bool is_fitted = fit_data(d, temp_np, temp_nq, rk, args, p, q, delta, linf_dist, l2_dist);
            if(is_fitted)
            {
                #pragma omp critical (nb_sol_found)
                {
                    ++nb_sol_found ;
                    mean_delta += delta ;

                    std::cout << "<<INFO>> found a solution with np=" << temp_np << ", nq = " << temp_nq << std::endl;
                    std::cout << "<<INFO>> Linf error = " << linf_dist << std::endl;
                    std::cout << "<<INFO>> L2   error = " << l2_dist << std::endl;
                    std::cout << "<<INFO>>      delta = " << delta << std::endl;
                    std::cout << std::endl;

                    // Get the solution with the minimum delta, and update the main
                    // rational function r.
                    if(delta < min_delta)
                    {
                        min_delta = delta ;
                        r->setSize(temp_np, temp_nq);
                        for(int y=0; y<r->dimY(); ++y)
                        {
                            r->update(y, rk->get(y));
                        }
                    }
                }
            }

				if(rk != NULL)
	            delete rk; // memory clean

            #pragma omp critical (nb_sol_tested)
            {
                // Update the solution
                nb_sol_tested++;

                std::cout << "<<DEBUG>> nb solutions tested: " << nb_sol_tested << " / " << i << "\r";
                std::cout.flush();
            }
        }

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
	// Fit the different output dimension independantly
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
// the function returns a rational BRDF function and a boolean
bool rational_fitter_parallel::fit_data(const vertical_segment* d, int np, int nq, int ny,
                                        rational_function_1d* r, vec& p, vec& q, double& delta)
{
	const int m = d->size(); // 2*m = number of constraints
	const int n = np+nq;     // n = np+nq

    quadratic_program qp(np, nq);

    // Starting with only a nb_starting_points vertical segments
    std::list<unsigned int> training_set;
    const int di = std::max((m-1) / (nb_starting_points-1), 1);
    for(int i=0; i<m; ++i)
	{
        if(i % di == 0)
        {
            // Create two vector of constraints
            vec c1(n), c2(n);
            get_constraint(i, np, nq, ny, d, r, c1, c2);

            qp.add_constraints(c1);
            qp.add_constraints(c2);

        }
        else
        {
            training_set.push_back(i);
        }
	}
    qp.set_training_set(training_set);

    do
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
    } while(qp.nb_constraints() < 2*m);

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

	// Create two vectors of constraints
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
