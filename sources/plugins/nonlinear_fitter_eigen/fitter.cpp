#include "fitter.h"

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <unsupported/Eigen/LevenbergMarquardt>

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <core/common.h>

ALTA_DLL_EXPORT fitter* provide_fitter()
{
    return new nonlinear_fitter_eigen();
}

struct EigenFunctor: Eigen::DenseFunctor<double>
{
	EigenFunctor(const ptr<nonlinear_function>& f, const ptr<data> d, bool use_cosine) :
		Eigen::DenseFunctor<double>(f->nbParameters(), d->dimY()*d->size()), _f(f), _d(d), _cosine(use_cosine)
	{
#ifndef DEBUG
		std::cout << "<<DEBUG>> constructing an EigenFunctor for n=" << inputs() << " parameters and m=" << values() << " points" << std::endl ;
#endif
	}

	int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& y) const
	{
#ifdef DEBUG
		std::cout << "parameters:" << std::endl << x << std::endl << std::endl ;
#endif

		// Check that the parameters used are within the bounds defined
		// by the function
		vec p_min = _f->getParametersMin();
		vec p_max = _f->getParametersMax();
		for(int i=0; i<_f->nbParameters(); ++i)
		{
			if(x[i] < p_min[i] || x[i] > p_max[i])
			{
				return -1;
			}
		}

		// Update the parameters vector
		vec _p(inputs());
		for(int i=0; i<inputs(); ++i) { _p[i] = x(i); }
		_f->setParameters(_p);

		const int nx = _f->dimX();
		const int ny = _f->dimY();

		for(int s=0; s<_d->size(); ++s)
		{
			vec _x  = _d->get(s);
		
			// Convert the sample point into the function space
			vec x(nx);
			params::convert(&_x[0], _d->input_parametrization(), _f->input_parametrization(), &x[0]);

			// Compute the cosine factor. Only update the constant if the flag
			// is set in the object.
			double cos = 1.0;
			if(_cosine)
			{
				double cart[6]; params::convert(&_x[0], _d->input_parametrization(), params::CARTESIAN, cart);
				cos = cart[5];
			}

			vec _di = vec(ny);
			for(int i=0; i<ny; ++i)
				_di[i] = _x[_d->dimX() + i];

			// Should add the resulting vector completely
			vec _y = _di - cos*_f->value(x);
			for(int i=0; i<ny; ++i)
				y(i*_d->size() + s) = _y[i];

		}
#ifdef DEBUG
		std::cout << "diff vector:" << std::endl << y << std::endl << std::endl ;
#endif

		return 0;
	}

	int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) const
	{
		// Update the paramters
		vec _p(inputs());
		for(int i=0; i<inputs(); ++i) { _p[i] = x(i); }
		_f->setParameters(_p);

		// For each element to fit, fill the rows of the matrix
		for(int s=0; s<_d->size(); ++s)
		{
			// Get the position
			vec xi = _d->get(s);
			
			// Convert the sample point into the function space
			vec x(_f->dimX());
			params::convert(&xi[0], _d->input_parametrization(), _f->input_parametrization(), &x[0]);

			// Compute the cosine factor. Only update the constant if the flag
			// is set in the object.
			double cos = 1.0;
			if(_cosine)
			{
				double cart[6]; params::convert(&xi[0], _d->input_parametrization(), params::CARTESIAN, cart);
				cos = cart[5];
			}

			// Get the associated jacobian
			vec _jac = _f->parametersJacobian(x);

			// Fill the columns of the matrix
			for(int j=0; j<_f->nbParameters(); ++j)
			{
				// For each output channel, update the subpart of the
				// vector row
				for(int i=0; i<_f->dimY(); ++i)
				{
					fjac(i*_d->size() + s, j) = - cos * _jac[i*_f->nbParameters() + j];
				}
			}
		}
		return 0;
	}


	const ptr<nonlinear_function>& _f;
	const ptr<data> _d;

	// Flags
	bool _cosine;
};

struct CompoundFunctor: Eigen::DenseFunctor<double>
{
	CompoundFunctor(compound_function* f, int index, const ptr<data> d, bool use_cosine) :
		Eigen::DenseFunctor<double>((*f)[index]->nbParameters(), d->dimY()*d->size()), _f(f), _d(d), _index(index), _cosine(use_cosine)
	{
#ifndef DEBUG
		std::cout << "<<DEBUG>> constructing an EigenFunctor for n=" << inputs() << " parameters and m=" << values() << " points" << std::endl ;
#endif
	}

	int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& y) const
	{
#ifdef DEBUG
		std::cout << "parameters:" << std::endl << x << std::endl << std::endl ;
#endif

		// Update the parameters vector
		vec _p(inputs());
		for(int i=0; i<inputs(); ++i) { _p[i] = x(i); }
		nonlinear_function* f = (*_f)[_index];
		f->setParameters(_p);

		for(int s=0; s<_d->size(); ++s)
		{
			vec _x  = _d->get(s);

			// Compute the cosine factor. Only update the constant if the flag
			// is set in the object.
			double cos = 1.0;
			if(_cosine)
			{
				double cart[6]; params::convert(&_x[0], _d->input_parametrization(), params::CARTESIAN, cart);
				cos = cart[5];
			}

			vec _di = vec(f->dimY());
			for(int i=0; i<f->dimY(); ++i)
				_di[i] = _x[f->dimX() + i];

			// Compute the value of the preceding functions
			vec _fy = vec::Zero(f->dimY());
			for(int i=0; i<_index+1; ++i)
			{
                const nonlinear_function* f = (*_f)[i];
                if(f->input_parametrization() != _d->input_parametrization())
                {
                    vec x(f->dimX());
                    params::convert(&_x[0], _d->input_parametrization(), f->input_parametrization(), &x[0]);

                    _fy += (*f)(x);
                }
                else
                {
                    _fy += (*f)(_x);
                }
			}

			// Should add the resulting vector completely
			vec _y = _di - cos*_fy;
			for(int i=0; i<f->dimY(); ++i)
				y(i*_d->size() + s) = _y[i];

		}
#ifdef DEBUG
		std::cout << "diff vector:" << std::endl << y << std::endl << std::endl ;
#endif

		return 0;
	}

	int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) const
	{
		// Update the paramters
		vec _p(inputs());
		for(int i=0; i<inputs(); ++i) { _p[i] = x(i); }

		nonlinear_function* f = (*_f)[_index];
		f->setParameters(_p);

		// For each element to fit, fill the rows of the matrix
		for(int s=0; s<_d->size(); ++s)
		{
			// Get the position
			vec xi = _d->get(s);

			// Compute the cosine factor. Only update the constant if the flag
			// is set in the object.
			double cos = 1.0;
			if(_cosine)
			{
				double cart[6]; params::convert(&xi[0], _d->input_parametrization(), params::CARTESIAN, cart);
				cos = cart[5];
			}

			// Get the associated jacobian
            vec _jac;
            if(f->input_parametrization() != _d->input_parametrization())
            {
                vec x(f->dimX());
                params::convert(&xi[0], _d->input_parametrization(), f->input_parametrization(), &x[0]);

                _jac = f->parametersJacobian(x);
            }
            else
            {
                _jac = f->parametersJacobian(xi);
            }


			// Fill the columns of the matrix
			for(int j=0; j<f->nbParameters(); ++j)
			{
				// For each output channel, update the subpart of the
				// vector row
				for(int i=0; i<_f->dimY(); ++i)
				{
					fjac(i*_d->size() + s, j) = - cos * _jac[i*f->nbParameters() + j];
				}
			}
		}
		return 0;
	}


	compound_function* _f;
	const ptr<data> _d;

	// Flags
	bool _cosine;
	int _index;
};

nonlinear_fitter_eigen::nonlinear_fitter_eigen() 
{
}
nonlinear_fitter_eigen::~nonlinear_fitter_eigen() 
{
}

bool nonlinear_fitter_eigen::fit_data(const ptr<data>& d, ptr<function>& fit, const arguments &args)
{
    // I need to set the dimension of the resulting function to be equal
    // to the dimension of my fitting problem
    fit->setDimX(d->dimX()) ;
    fit->setDimY(d->dimY()) ;
    fit->setMin(d->min()) ;
    fit->setMax(d->max()) ;

	 // Convert the function and bootstrap it with the data
    if(!dynamic_pointer_cast<nonlinear_function>(fit))
    {
        std::cerr << "<<ERROR>> the function is not a non-linear function" << std::endl;
        return false;
    }
    ptr<nonlinear_function> nf = dynamic_pointer_cast<nonlinear_function>(fit);
	 nf->bootstrap(d, args);

#ifndef DEBUG
	 std::cout << "<<DEBUG>> number of parameters: " << nf->nbParameters() << std::endl;
#endif
	 if(nf->nbParameters() == 0)
	 {
		 return true;
	 }

    /* the following starting values provide a rough fit. */
    vec nf_x = nf->parameters();

    int info;
#ifdef OLD
	 if(args.is_defined("fit-compound"))
	 {
		 ptr<compound_function> compound = dynamic_pointer_cast<compound_function>(nf);
		 if(!compound_function)
		 {
			 std::cerr << "<<ERROR>> you should use --fit-compound with a compound function" << std::endl;
			 return false;
		 }

		 // Register how many parameter are already fitted
		 int already_fit = 0;

		 // Get the i-th function of the compound
		 for(int index=0; index<compound->size(); ++index)
		 {
			 nonlinear_function* f = (*compound)[index];
			 if(f->nbParameters() == 0)
				 continue;

			 Eigen::VectorXd x(f->nbParameters());
			 for(int i=0; i<f->nbParameters(); ++i)
			 {
				 x[i] = nf_x[i+already_fit];
			 }

			 CompoundFunctor functor(compound, index, d, args.is_defined("fit-with-cosine"));
			 Eigen::LevenbergMarquardt<CompoundFunctor> lm(functor);

			 info = lm.minimize(x);

			 if(info == Eigen::LevenbergMarquardtSpace::ImproperInputParameters)
			 {
				 std::cerr << "<<ERROR>> incorrect parameters" << std::endl;
				 return false;
			 }
			else if(info == Eigen::LevenbergMarquardtSpace::UserAsked)
			 {
				 std::cerr << "<<ERROR>> the search is using improper parameters: stopping" << std::endl;
				 return false;
			 }

			 // Update the vector of parameters
			 for(int i=0; i<f->nbParameters(); ++i)
			 {
				 nf_x[i+already_fit] = x[i];
			 }

			 // Update the number of already fitted parameters
			 already_fit += f->nbParameters();

#ifndef DEBUG
    std::cout << "<<DEBUG>> function " << index+1 << " using " << lm.iterations() << " iterations" << std::endl;
#endif
		 }
	 }
	 else
#endif
	 {
		 Eigen::VectorXd x(nf->nbParameters());
		 for(int i=0; i<nf->nbParameters(); ++i)
		 {
			 x[i] = nf_x[i];
		 }

		 EigenFunctor functor(nf, d, args.is_defined("fit-with-cosine"));
		 Eigen::LevenbergMarquardt<EigenFunctor> lm(functor);

		 info = lm.minimize(x);

		 if(info == Eigen::LevenbergMarquardtSpace::ImproperInputParameters)
		 {
			 std::cerr << "<<ERROR>> incorrect parameters" << std::endl;
			 return false;
		 }
		 else if(info == Eigen::LevenbergMarquardtSpace::UserAsked)
		 {
			 std::cerr << "<<ERROR>> the search is using improper parameters: stopping" << std::endl;
			 return false;
		 }

		 for(int i=0; i<nf->nbParameters(); ++i)
		 {
			 nf_x[i] = x(i);
		 }
#ifndef DEBUG
    std::cout << "<<DEBUG>> using " << lm.iterations() << " iterations" << std::endl;
#endif
	 }

    std::cout << "<<INFO>> found parameters: " << nf_x << std::endl;
    nf->setParameters(nf_x);
    return true;

}

void nonlinear_fitter_eigen::set_parameters(const arguments& args)
{
}
