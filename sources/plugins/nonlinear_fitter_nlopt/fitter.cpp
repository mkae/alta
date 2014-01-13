#include "fitter.h"

#include <nlopt.h>

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

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
		// Get the data sample and extract the value part
		vec xi = d->get(s);
		vec _di = vec(d->dimY());
		for(int i=0; i<d->dimY(); ++i)
		{
			_di[i] = xi[d->dimX() + i];
		}
		
		// Convert the sample point into the function space
		vec x(f->dimX());
		params::convert(&xi[0], d->input_parametrization(), f->input_parametrization(), &x[0]);

		// Get the jacobian of the function at position x_i for the current
		// set of parameters (set prior to function call)
		vec _jac = f->parametersJacobian(x);

		// Should add the resulting vector completely
		vec _y = (*f)(x) - _di;

		// Fill the columns of the matrix
		for(int j=0; j<f->nbParameters(); ++j)
		{
			// For each output channel, update the subpart of the
			// vector row
			for(int i=0; i<f->dimY(); ++i)
			{
				fjac[j] += 2 * _y[i] * _jac[i*f->nbParameters() + j];
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
		// Get the data sample and extract the value part
		vec xi = _d->get(s);
		vec _di = vec(_d->dimY());
		for(int i=0; i<_d->dimY(); ++i)
		{
			_di[i] = xi[_d->dimX() + i];
		}

		// Convert the sample point into the function space
		vec x(_f->dimX());
		params::convert(&xi[0], _d->input_parametrization(), _f->input_parametrization(), &x[0]);

		// Should add the resulting vector completely
		vec _y = (*_f)(x) - _di;
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
	if(optimizerName == "local-COBYLA")
	{
		algorithm = NLOPT_LN_COBYLA;
	}
	else if(optimizerName == "local-BOBYQA")
	{
		algorithm = NLOPT_LN_BOBYQA;
	}
	else if(optimizerName == "local-NEWUOA")
	{
		algorithm = NLOPT_LN_NEWUOA_BOUND;
	}
	else if(optimizerName == "local-PRAXIS")
	{
		algorithm = NLOPT_LN_PRAXIS;
	}
    else if(optimizerName == "local-Nelder-Mead")
	{
		algorithm = NLOPT_LN_NELDERMEAD;
	}
	else if(optimizerName == "local-spblx")
	{
		algorithm = NLOPT_LN_SBPLX;
	}
	else if(optimizerName == "global-crs")
	{
		algorithm = NLOPT_GN_CRS2_LM;
	}
	else if(optimizerName == "global-stogo")
	{
		algorithm = NLOPT_GD_STOGO;
	}
	else if(optimizerName == "global-isres")
	{
		algorithm = NLOPT_GN_ISRES;
	}
	else if(optimizerName == "local-mma")
	{
		algorithm = NLOPT_LD_MMA; 
	}
	else if(optimizerName == "local-sqp")
	{
		algorithm = NLOPT_LD_SLSQP;
	}
	else
	{
		algorithm = NLOPT_LD_SLSQP;
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
