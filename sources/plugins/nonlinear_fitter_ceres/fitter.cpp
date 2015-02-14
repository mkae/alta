/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "fitter.h"

#include <Eigen/Dense>
#include <ceres/ceres.h>

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <core/common.h>

ALTA_DLL_EXPORT fitter* provide_fitter()
{
    return new nonlinear_fitter_ceres();
}

class CeresFunctor : public ceres::CostFunction
{
	public:
		CeresFunctor(const ptr<nonlinear_function>& f, const vec& xi, const arguments& args) : _f(f), _xi(xi)
		{
			set_num_residuals(f->dimY());
			mutable_parameter_block_sizes()->push_back(f->nbParameters());

			_log_fit = args.is_defined("log-fit");
		}

		virtual bool Evaluate(double const* const* x, double* y, double** dy) const
		{
			// Check that the parameters used are within the bounds defined
			// by the function
			vec p_min = _f->getParametersMin();
			vec p_max = _f->getParametersMax();
			for(int i=0; i<_f->nbParameters(); ++i)
			{
				if(x[0][i] < p_min[i] || x[0][i] > p_max[i])
				{
					return false;
				}
			}

			// Update the parameters vector
			vec _p(_f->nbParameters());
			for(int i=0; i<_f->nbParameters(); ++i) { _p[i] = x[0][i]; }
			_f->setParameters(_p);

			vec _di = vec(_f->dimY());
			for(int i=0; i<_f->dimY(); ++i)
			{
				_di[i] = _xi[_f->dimX() + i];
			}

			const vec _yi = _f->value(_xi);
			const vec _y = _di - _yi;
			for(int i=0; i<_f->dimY(); ++i)
			{
				y[i] = (_log_fit) ? log(1.0 + _di[i]) - log(1.0 + _yi[i]) : _y[i] ;
			}

			if(dy != NULL)
			{
				df(_di, dy);
			}

			return true;
		}

		// The parameter of the function _f should be set prior to this function
		// call. If not it will produce undesirable results.
		virtual void df(const vec& di, double ** fjac) const
		{
			// Get the jacobian of the function at position x_i for the current
			// set of parameters (set prior to function call)
			vec _jac = _f->parametersJacobian(_xi);

			// For each output channel, update the subpart of the
			// vector row
			for(int i=0; i<_f->dimY(); ++i)
			{
				// Fill the columns of the matrix
				for(int j=0; j<_f->nbParameters(); ++j)
				{
                    fjac[0][i*_f->nbParameters() + j] = - ((_log_fit) ? _jac[i*_f->nbParameters() + j]/(1.0 + di[i]) : _jac[i*_f->nbParameters() + j]);
				}
         }
		}

	protected:

		// Data point and function to optimize
		const ptr<nonlinear_function>& _f;
		const vec _xi;

		// Arguments of the fitting procedure
		bool _log_fit;
};

nonlinear_fitter_ceres::nonlinear_fitter_ceres()
{
}
nonlinear_fitter_ceres::~nonlinear_fitter_ceres()
{
}

bool nonlinear_fitter_ceres::fit_data(const ptr<data>& d, ptr<function>& fit, const arguments &args)
{
    // I need to set the dimension of the resulting function to be equal
    // to the dimension of my fitting problem
    fit->setDimX(d->dimX()) ;
    fit->setDimY(d->dimY()) ;
    fit->setMin(d->min()) ;
    fit->setMax(d->max()) ;

	 // Convert the function and bootstrap it with the data
   ptr<nonlinear_function> nf = dynamic_pointer_cast<nonlinear_function>(fit);
    if(!nf)
    {
        std::cerr << "<<ERROR>> the function is not a non-linear function" << std::endl;
        return false;
    }


#ifndef DEBUG
	 std::cout << "<<DEBUG>> number of parameters: " << nf->nbParameters() << std::endl;
#endif
	 if(nf->nbParameters() == 0)
	 {
		 return true;
	 }

	 /* Bootstrap the function */
	 nf->bootstrap(d, args);

	 /* the following starting values provide a rough fit. */
	 vec p = nf->parameters();

	 std::cout << "<<DEBUG>> Starting vector: " << p << std::endl;
	 std::cout << "<<DEBUG>> Final vector should be between " << nf->getParametersMin() << " and " << nf->getParametersMax() << std::endl;

	 // Create the problem
	 ceres::Problem problem;
	 for(int i=0; i<d->size(); ++i)
	 {
		 vec xi = d->get(i);
		 vec xf(nf->dimX() + nf->dimY());

		 // Convert the sample to be in the parametrizatio of the function
		 params::convert(&xi[0], d->input_parametrization(), nf->input_parametrization(), &xf[0]);
		 for(int k=0; k<nf->dimY(); ++k)
		 {
			 xf[nf->dimX() + k] = xi[d->dimX() + k];
		 }

		 problem.AddResidualBlock(new CeresFunctor(nf, xf, args), NULL, &p[0]);
	 }

	 // Solves the NL problem
	 ceres::Solver::Summary summary;
	 ceres::Solve(options, &problem, &summary);


#ifdef DEBUG
	 std::cout << summary.BriefReport() << std::endl;
#endif
	 std::cout << "<<INFO>> found parameters: " << p << std::endl;

	 nf->setParameters(p);
	 return true;
}

void nonlinear_fitter_ceres::set_parameters(const arguments& args)
{
    std::string const MAX_NUM_ITER_OPTION = "ceres-max-num-iterations";
    if(args.is_defined(MAX_NUM_ITER_OPTION))
    {
      options.max_num_iterations = args.get_int(MAX_NUM_ITER_OPTION, 50); // Default value = 50
    }

    if(args["ceres-factorization"] == "normal-cholesky")
    {
      options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
    }
    else
    {
      options.linear_solver_type = ceres::DENSE_QR;
    }

    if(args.is_defined("ceres-debug"))
    {
      options.minimizer_progress_to_stdout = true; // Default value = false;
    }

    std::string const FUNCTION_TOLERANCE_OPTION = "ceres-function-tolerance";
    if( args.is_defined(FUNCTION_TOLERANCE_OPTION) )
    {
      options.function_tolerance = args.get_double(FUNCTION_TOLERANCE_OPTION, 1e-6); //Default value = 1e-6
    }

    std::string const GRADIENT_TOLERANCE_OPTION = "ceres-gradient-tolerance";
    if( args.is_defined(GRADIENT_TOLERANCE_OPTION))
    {
      options.function_tolerance = args.get_double(GRADIENT_TOLERANCE_OPTION, 1e-10); //Default value = 1e-10
    }

    std::string const PARAMETER_TOLERANCE_OPTION = "ceres-parameter-tolerance";
    if( args.is_defined(PARAMETER_TOLERANCE_OPTION))
    {
      options.function_tolerance = args.get_double(PARAMETER_TOLERANCE_OPTION, 1e-8); //Default value = 1e-8
    }


}
