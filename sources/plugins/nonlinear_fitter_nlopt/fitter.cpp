#include "fitter.h"

#include <nlopt.h>

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <QTime>

#include <core/common.h>

ALTA_DLL_EXPORT fitter* provide_fitter()
{
	return new nonlinear_fitter_nlopt();
}

void print_nlopt_error(nlopt_result res, const std::string& string)
{
	if(res == NLOPT_FAILURE)
	{
		std::cerr << "<<ERROR>> generic failure for \"" << string << "\"" << std::endl;
	}
	else if(res == NLOPT_INVALID_ARGS)
	{
		std::cerr << "<<ERROR>> invalid arguments for \"" << string << "\"" << std::endl;
	}
	else if(res == NLOPT_OUT_OF_MEMORY)
	{
		std::cerr << "<<ERROR>> not enough memory for \"" << string << "\"" << std::endl;
	}
}

// The parameter of the function _f should be set prior to this function
// call. If not it will produce undesirable results.
void df(double* fjac, const nonlinear_function* f, const data* d)
{
	// Clean memory
	memset(fjac, 0.0, f->nbParameters()*sizeof(double));

	// Each constraint is of the form data point * color channel
	for(int s=0; s<d->size(); ++s)
	{
		vec xi = d->get(s);

		// Get the jacobian of the function at position x_i for the current
		// set of parameters (set prior to function call)
		vec _jac = f->parametersJacobian(xi);

		vec _di = vec(f->dimY());
		for(int i=0; i<f->dimY(); ++i)
		{
			_di[i] = xi[f->dimX() + i];
		}

		// Should add the resulting vector completely
		vec _y = (*f)(xi) - _di;

		// For each output channel, update the subpart of the
		// vector row
		for(int i=0; i<f->dimY(); ++i)
		{
			// Fill the columns of the matrix
			for(int j=0; j<f->nbParameters(); ++j)
			{
				fjac[j] = 2 * _y[i] * _jac[i*f->nbParameters() + j];
			}
		}
	}
}

double f(unsigned n, const double* x, double* dy, void* dat)
{
	nonlinear_function* _f = (nonlinear_function*)(((void**)dat)[0]);
	const data* _d = (const data*)(((void**)dat)[1]);

	// Update the parameters vector
	vec _p(_f->nbParameters());
	for(int i=0; i<_f->nbParameters(); ++i) { _p[i] = x[i]; }
	_f->setParameters(_p);

	// Create the result vector
	double y = 0.0;

	// Each constraint is of the form data point * color channel
	for(int s=0; s<_d->size(); ++s)
	{
		vec xi = _d->get(s);
		vec _di = vec(_f->dimY());
		for(int i=0; i<_f->dimY(); ++i)
		{
			_di[i] = xi[_f->dimX() + i];
		}

		// Should add the resulting vector completely
		vec _y = (*_f)(xi) - _di;
		for(int i=0; i<_f->dimY(); ++i)
		{
			y += pow(_y[i], 2);
		}
	}

	if(dy != NULL)
	{
		df(dy, _f, _d);
	}

	return y;
}

nonlinear_fitter_nlopt::nonlinear_fitter_nlopt() 
{
}
nonlinear_fitter_nlopt::~nonlinear_fitter_nlopt() 
{
}

bool nonlinear_fitter_nlopt::fit_data(const data* d, function* fit, const arguments &args)
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

	// the following starting values provide a rough fit is the bootstrap flag is
	// enabled
	vec p = nf->parameters();

	// Check the optimizer to use 
	nlopt_algorithm algorithm;
	nlopt_result res;
	std::string optimizerName = args["nlopt-optimizer"];
	if(optimizerName == "COBYLA")
	{
		algorithm = NLOPT_LN_COBYLA;
	}
	else if(optimizerName == "BOBYQA")
	{
		algorithm = NLOPT_LN_BOBYQA;
	}
	else if(optimizerName == "NEWUOA")
	{
		algorithm = NLOPT_LN_NEWUOA_BOUND;
	}
	else if(optimizerName == "PRAXIS")
	{
		algorithm = NLOPT_LN_PRAXIS;
	}
    else if(optimizerName == "Nelder-Mead")
	{
		algorithm = NLOPT_LN_NELDERMEAD;
	}
	else if(optimizerName == "Sbplx")
	{
		algorithm = NLOPT_LN_SBPLX;
	}
	else
	{
		// The default option one can find in the NLOpt documentation
		algorithm = NLOPT_LD_MMA; 
	}

	// Create the optimizer
	nlopt_opt opt = nlopt_create(algorithm, nf->nbParameters());
	if(opt == NULL)
	{
		std::cerr << "<<ERROR>> unable to create the optimizer" << std::endl;
		return false;
	}


	// Set the bounds to the parameters
	vec p_min = nf->getParametersMin();
	vec p_max = nf->getParametersMax();
	nlopt_set_lower_bounds(opt, &p_min[0]);
	nlopt_set_upper_bounds(opt, &p_max[0]);

	// Set the stopping criterion to a maximum number of evaluation
	if(args.is_defined("nlop-max-num-iterations"))
	{
		res = nlopt_set_maxeval(opt, args.get_int("nlop-max-num-iterations", 10));
		if(res < 0) { print_nlopt_error(res, "nlopt_set_maxeval"); }
	}
	nlopt_set_xtol_rel(opt, 1e-4);

	// Create the problem
	void* dat[2];
	dat[0] = (void*)nf;
	dat[1] = (void*)d;
	res = nlopt_set_min_objective(opt, f, dat);
	if(res < 0)
	{
		print_nlopt_error(res, "nlopt_set_min_objective");
		return false;
	}


	// Launch minimization
	double f_end;
	res = nlopt_optimize(opt, &p[0], &f_end);

	if(res > 0)
	{
		std::cout << "<<INFO>> optimized distance: " << f_end << std::endl;
		std::cout << "<<INFO>> found parameters: " << p << std::endl;
		nf->setParameters(p);
	}
	else
	{
		print_nlopt_error(res, "nlopt_optimize");
	}

	nlopt_destroy(opt);
	return res > 0;
}

void nonlinear_fitter_nlopt::set_parameters(const arguments& args)
{
}
