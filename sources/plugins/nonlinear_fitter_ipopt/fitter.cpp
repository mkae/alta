#include "fitter.h"

#include <Eigen/Dense>

#include <coin/IpTNLP.hpp>
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <core/common.h>

ALTA_DLL_EXPORT fitter* provide_fitter()
{
    return new nonlinear_fitter_ipopt();
}

class altaNLP : public Ipopt::TNLP
{
	public:
		altaNLP(nonlinear_function* f, const data* d) : TNLP(), _f(f), _d(d)
		{
		}

		// Return the size of the NL problem
		virtual bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
		                          Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style)
		{
			// Size of the problem
			n = _f->nbParameters();
			m = 0;
			
			// Nonzeros in the Jacobian and Hessian
			nnz_jac_g = 0;
			nnz_h_lag = 0;

			// C++ indexing of arrays
			index_style = C_STYLE;

			return true;
		}

		virtual bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
		                             Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u)
		{
			// Check the size of the problem
			assert(n == _f->nbParameters());
			assert(m == 0);

			// min and Max values for the parameters
			vec p_min = _f->getParametersMin();
			vec p_max = _f->getParametersMax();

			for(int i=0; i<n; ++i)
			{
				x_l[i] = p_min[i];
				x_u[i] = p_max[i];
			}

			return true;
		}

		virtual bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
		                                bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
		                                Ipopt::Index m, bool init_lambda,
		                                Ipopt::Number* lambda)
		{
			// Check the input
			assert(n == _f->nbParameters());
			assert(init_x == true);

			vec p = _f->parameters();
			for(int i=0; i<n; ++i)
			{
				x[i] = p[i];
			}

			return true;
		}
		

		virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, 
		                    Ipopt::Number& obj_value)
		{
			// Update the parameters vector
			vec _p(n);
			for(int i=0; i<n; ++i) { _p[i] = x[i]; }
			_f->setParameters(_p);

			obj_value = 0.0;
			for(int s=0; s<_d->size(); ++s)
			{
				vec _x  = _d->get(s);

				// Extract the objective from the current vector
				vec _di = vec(_f->dimY());
				for(int i=0; i<_f->dimY(); ++i)
				{
					_di[i] = _x[_f->dimX() + i];
				}

				// Convert the sample point into the function space
				vec x(_f->dimX());
				params::convert(&_x[0], _d->input_parametrization(), _f->input_parametrization(), &x[0]);

				// Compute the difference vector and add its
				// components to the obj_value
				vec _y = _di - (*_f)(x);
				for(int i=0; i<_f->dimY(); ++i)
				{
					obj_value += pow(_y[i], 2);
				}
			}

			return true;
		}
		
		virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, 
		                         bool new_x, Ipopt::Number* grad_f)
		{
			// Update the parameters vector
			vec _p(n);
			for(int i=0; i<n; ++i) { _p[i] = x[i]; }
			_f->setParameters(_p);

			// Clean the value
			for(int i=0; i<n; ++i) { grad_f[i] = 0.0; }

			// Add all the gradients
			for(int s=0; s<_d->size(); ++s)
			{
				vec _x  = _d->get(s);
				
				// Extract the objective from the current vector
				vec _di = vec(_f->dimY());
				for(int i=0; i<_f->dimY(); ++i)
				{
					_di[i] = _x[_f->dimX() + i];
				}
				
				// Convert the sample point into the function space
				vec x(_f->dimX());
				params::convert(&_x[0], _d->input_parametrization(), _f->input_parametrization(), &x[0]);

				// Compute the difference vector and add its
				// components to the obj_value
				vec _y = (*_f)(x) - _di;

				// Get the jacobian of the function at position x_i for the current
				// set of parameters (set prior to function call)
				vec _jac = _f->parametersJacobian(x);

				// Fill the columns of the matrix
				for(int j=0; j<_f->nbParameters(); ++j)
				{
					// For each output channel, update the subpart of the
					// vector row
					for(int i=0; i<_f->dimY(); ++i)
					{
						grad_f[j] += 2 * _y[i] * _jac[i*_f->nbParameters() + j];
					}
				}
			}

			return true;
		}
		
		virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, 
		                    Ipopt::Index m, Ipopt::Number* g)
		{
			return true;
		}
		
		virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                              Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, 
										Ipopt::Index *jCol, Ipopt::Number* values)
		{
			return true;
		}

		/*
		virtual bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                          Ipopt::Number obj_factor, Ipopt::Index m, 
								  const Ipopt::Number* lambda, bool new_lambda, 
								  Ipopt::Index nele_hess, Ipopt::Index* iRow,
                          Ipopt::Index* jCol, Ipopt::Number* values)
		{
			return false;
		}
		*/

		virtual void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, 
		                               const Ipopt::Number* x, const Ipopt::Number* z_L, 
												 const Ipopt::Number* z_U, Ipopt::Index m, 
												 const Ipopt::Number* g, const Ipopt::Number* lambda,
                                     Ipopt::Number obj_value, const Ipopt::IpoptData* ip_data,
		                               Ipopt::IpoptCalculatedQuantities* ip_cq)
		{
			// Update the parameters vector
			vec _p(n);
			for(int i=0; i<n; ++i) { _p[i] = x[i]; }
			_f->setParameters(_p);
		}

	protected:

		const data* _d;
		nonlinear_function* _f;
};

nonlinear_fitter_ipopt::nonlinear_fitter_ipopt() 
{
}
nonlinear_fitter_ipopt::~nonlinear_fitter_ipopt() 
{
}

bool nonlinear_fitter_ipopt::fit_data(const data* d, function* fit, const arguments &args)
{
	// I need to set the dimension of the resulting function to be equal
	// to the dimension of my fitting problem
	fit->setDimX(d->dimX()) ;
	fit->setDimY(d->dimY()) ;
	fit->setMin(d->min()) ;
	fit->setMax(d->max()) ;

	// Convert the function and bootstrap it with the data
	if(dynamic_cast<nonlinear_function*>(fit) == NULL)
	{
		std::cerr << "<<ERROR>> the function is not a non-linear function" << std::endl;
		return false;
	}
	nonlinear_function* nf = dynamic_cast<nonlinear_function*>(fit);
	nf->bootstrap(d, args);

#ifndef DEBUG
	std::cout << "<<DEBUG>> number of parameters: " << nf->nbParameters() << std::endl;
#endif
	if(nf->nbParameters() == 0)
	{
		return true;
	}

	/* the following starting values provide a rough fit. */
	vec p = nf->parameters();

	// Create the problem
	Ipopt::SmartPtr<Ipopt::TNLP> nlp = new altaNLP(nf, d);
	Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();


	Ipopt::ApplicationReturnStatus status;
	status = app->Initialize();
	if(status != Ipopt::Solve_Succeeded) 
	{
		std::cout << "<<ERROR>> unable to create Ipopt solver" << std::endl;
		return false;
	}

	// Static parameters for the Solver
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");
#ifdef DEBUG
	app->Options()->SetIntegerValue("print_level", 5);
#else
	app->Options()->SetIntegerValue("print_level", 0);
#endif

	// Solver's options
	app->Options()->SetIntegerValue("max_iter", args.get_int("ipopt-max-iter", 10));
	app->Options()->SetStringValue("linear_solver", args.get_string("ipopt-solver", "mumps"));
	

	// Solves the NL problem
	status = app->OptimizeTNLP(nlp);
	if(status == Ipopt::Solve_Succeeded)
	{
#ifdef DEBUG
		Ipopt::Index iter_count = app->Statistics()->IterationCount();
		std::cout << "<<DEBUG>> number of iterations: " << iter_count << std::endl;
#endif
			std::cout << "<<INFO>> found parameters: " << nf->parameters() << std::endl;

		return true;
	}
	else if(status == Ipopt::Maximum_Iterations_Exceeded)
	{
		std::cout << "<<INFO>> the maximum number of iteration has been reached" << std::endl;
#ifdef DEBUG
		Ipopt::Index iter_count = app->Statistics()->IterationCount();
		std::cout << "<<DEBUG>> number of iterations: " << iter_count << std::endl;
#endif
		std::cout << "<<INFO>> found parameters: " << nf->parameters() << std::endl;

		return true;
	}
	else
	{
		return false;
	}
}

void nonlinear_fitter_ipopt::set_parameters(const arguments& args)
{
}
