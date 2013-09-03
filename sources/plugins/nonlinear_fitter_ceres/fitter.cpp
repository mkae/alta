#include "fitter.h"

#include <Eigen/Dense>
#include <ceres/ceres.h>

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
    return new nonlinear_fitter_ceres();
}

class CeresFunctor : public ceres::CostFunction
{
	public:
		CeresFunctor(nonlinear_function* f, const vec xi) : _f(f), _xi(xi)
		{
			set_num_residuals(f->dimY());
			mutable_parameter_block_sizes()->push_back(f->nbParameters());
		}

		virtual bool Evaluate(double const* const* x, double* y, double** dy) const 
		{
			// Update the parameters vector
			vec _p(_f->nbParameters());
			for(int i=0; i<_f->nbParameters(); ++i) { _p[i] = x[0][i]; }
			_f->setParameters(_p);

			vec _di = vec(_f->dimY());
			for(int i=0; i<_f->dimY(); ++i)
			{
				_di[i] = _xi[_f->dimX() + i];
			}

			// Should add the resulting vector completely
			vec _y = _di - (*_f)(_xi);
			for(int i=0; i<_f->dimY(); ++i)
			{
				y[i] = _y[i];
			}

			if(dy != NULL)
			{
				df(dy);
			}

			return true;
		}

		// The parameter of the function _f should be set prior to this function
		// call. If not it will produce undesirable results.
		virtual void df(double ** fjac) const
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
					fjac[0][i*_f->nbParameters() + j] = -_jac[i*_f->nbParameters() + j];
				}
			}
		}

	protected:

		const vec _xi;
		nonlinear_function* _f;
};

nonlinear_fitter_ceres::nonlinear_fitter_ceres() 
{
}
nonlinear_fitter_ceres::~nonlinear_fitter_ceres() 
{
}

bool nonlinear_fitter_ceres::fit_data(const data* d, function* fit, const arguments &args)
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
	 ceres::Problem problem;
	 for(int i=0; i<d->size(); ++i)
	 {
		 vec xi = d->get(i);
		 problem.AddResidualBlock(new CeresFunctor(nf, xi), NULL, &p[0]);
	 }

	 // Solver's options
	 ceres::Solver::Options options;
	 if(args.is_defined("ceres-max-num-iterations"))
	 {
		 options.max_num_iterations = args.get_int("ceres-max-num-iterations", 50); // Default value = 50
	 }
	 options.linear_solver_type = ceres::DENSE_QR;
	 if(args.is_defined("ceres-debug"))
	 {
		 options.minimizer_progress_to_stdout = true; // Default value = false;
	 }


	 // Solves the NL problem
	 ceres::Solver::Summary summary;
	 ceres::Solve(options, &problem, &summary);

	 std::cout << summary.BriefReport() << std::endl;
    std::cout << "<<INFO>> found parameters: " << p << std::endl;
    nf->setParameters(p);
    return true;
}

void nonlinear_fitter_ceres::set_parameters(const arguments& args)
{
}
