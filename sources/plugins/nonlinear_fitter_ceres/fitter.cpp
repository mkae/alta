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
		CeresFunctor(nonlinear_function* f, const vec xi) : _f(f), _xi(xi)
		{
			set_num_residuals(f->dimY());
			mutable_parameter_block_sizes()->push_back(f->nbParameters());
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

/*! \brief A CostFunction to fit one color channel at a time.
 *
 *  \details This forces the functions to have a full formulation. All the
 *  parameters need to have a value per color channel.
 */
class ColorChannelCost : public ceres::CostFunction
{
    public:
    ColorChannelCost(nonlinear_function* f, const vec xi, int j) : _f(f), _xi(xi), _j(j)
        {
            // The function should be set 1D prior to this call
            assert(f->dimY() == 1);

            set_num_residuals(1);
            mutable_parameter_block_sizes()->push_back(f->nbParameters());
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

            double _di = _xi[_f->dimX() + _j];

            // Should add the resulting vector completely
            y[0] = _di - (*_f)(_xi)[0];

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

            // Fill the columns of the matrix
            for(int j=0; j<_f->nbParameters(); ++j)
            {
                fjac[0][j] = -_jac[j];
            }
        }

    protected:

        const vec _xi;
        int _j;
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
    fit->setMin(d->getMin()) ;
    fit->setMax(d->getMax()) ;

	 // Convert the function and bootstrap it with the data
    if(dynamic_cast<nonlinear_function*>(fit) == NULL)
    {
        std::cerr << "<<ERROR>> the function is not a non-linear function" << std::endl;
        return false;
    }
    nonlinear_function* nf = dynamic_cast<nonlinear_function*>(fit);


#ifndef DEBUG
	 std::cout << "<<DEBUG>> number of parameters: " << nf->nbParameters() << std::endl;
#endif
	 if(nf->nbParameters() == 0)
	 {
		 return true;
	 }
#ifdef FIT_CHANNELS
     if(args.is_defined("ceres-channels"))
     {
         std::cout << "<<WARNING>> will fit the output dimensions separately" << std::endl;
         std::cout << "<<WARNING>> make sur the function is separable." << std::endl;
         return fit_channel(d, nf, args);
     }
     else
#endif
     {
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
             problem.AddResidualBlock(new CeresFunctor(nf, xi), NULL, &p[0]);
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
}

#ifdef FIT_CHANNELS
bool nonlinear_fitter_ceres::fit_channel(const data* d, nonlinear_function* nf,
                                         const arguments& args)
{
    /* the following starting values provide a rough fit. */
    vec p = nf->parameters();

    // Convert the current function to monochromatic
    nf->setDimY(1);


    for(int c=0; c<d->dimY(); ++c)
    {
        nf->bootstrap(d, args);

        // Temp parameter vector
        vec temp_p = nf->parameters();

        // Create the problem
        ceres::Problem problem;
        for(int i=0; i<d->size(); ++i)
        {
            vec xi = d->get(i);
            problem.AddResidualBlock(new ColorChannelCost(nf, xi, c), NULL, &temp_p[0]);
        }

        // Solves the NL problem
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);

#ifdef DEBUG
        std::cout << summary.BriefReport() << std::endl;
#endif

        // Update the resulting parameter vector
        for(int i=0; i<nf->nbParameters(); ++i)
        {
            p[c*nf->nbParameters() + i] = temp_p[i];
        }
    }

    std::cout << "<<INFO>> found parameters: " << p << std::endl;
    nf->setDimY(d->dimY());
    nf->setParameters(p);

    return true;
}
#endif


void nonlinear_fitter_ceres::set_parameters(const arguments& args)
{
    if(args.is_defined("ceres-max-num-iterations"))
    {
        options.max_num_iterations = args.get_int("ceres-max-num-iterations", 50); // Default value = 50
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
}
