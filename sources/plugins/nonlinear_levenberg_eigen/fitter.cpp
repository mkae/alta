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

#include <QTime>

#include <core/common.h>

ALTA_DLL_EXPORT fitter* provide_fitter()
{
    return new nonlinear_fitter_eigen();
}

struct EigenFunctor: Eigen::DenseFunctor<double>
{
    EigenFunctor(nonlinear_function* f, const data* d) : 
		 DenseFunctor<double>(f->nbParameters(), d->dimY()*d->size()), _f(f), _d(d)
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
		_f->setParameters(_p);

		for(int s=0; s<_d->size(); ++s)
		{
			vec _x  = _d->get(s);

			vec _di = vec(_f->dimY());
			for(int i=0; i<_f->dimY(); ++i)
				_di[i] = _x[_f->dimX() + i];

			// Should add the resulting vector completely
			vec _y = _di - (*_f)(_x);
			for(int i=0; i<_f->dimY(); ++i)
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

			// Get the associated jacobian
			vec _jac = _f->parametersJacobian(xi);

			// Fill the columns of the matrix
#ifdef DEBUG
			Eigen::MatrixXd temp (_f->dimY(), _f->nbParameters());
#endif
			for(int j=0; j<_f->nbParameters(); ++j)
			{
				// For each output channel, update the subpart of the
				// vector row
				for(int i=0; i<_f->dimY(); ++i)
				{
					fjac(i*_d->size() + s, j) = -_jac[i*_f->nbParameters() + j];
#ifdef DEBUG
					temp(i, j) = _jac[i*_f->nbParameters() + j];
#endif
				}
			}

#ifdef DEBUG
			std::cout << temp << std::endl << std::endl ;
#endif
		}
#ifdef DEBUG
			std::cout << "jacobian :" << std::endl << fjac << std::endl << std::endl;
#endif
		return 0;
	}


	nonlinear_function* _f;
    const data* _d;
};

nonlinear_fitter_eigen::nonlinear_fitter_eigen() : QObject()
{
}
nonlinear_fitter_eigen::~nonlinear_fitter_eigen() 
{
}

bool nonlinear_fitter_eigen::fit_data(const data* d, function* fit, const arguments &args)
{
    // I need to set the dimension of the resulting function to be equal
    // to the dimension of my fitting problem
    fit->setDimX(d->dimX()) ;
    fit->setDimY(d->dimY()) ;
    fit->setMin(d->min()) ;
    fit->setMax(d->max()) ;

	 // Convert the function and boostrap it with the data
    if(dynamic_cast<nonlinear_function*>(fit) == NULL)
    {
        std::cerr << "<<ERROR>> the function is not a non-linear function" << std::endl;
        return false;
    }
    nonlinear_function* nf = dynamic_cast<nonlinear_function*>(fit);
	 nf->boostrap(d, args);

#ifndef DEBUG
	 std::cout << "<<DEBUG>> number of parameters: " << nf->nbParameters() << std::endl;
#endif

    /* the following starting values provide a rough fit. */
    vec nf_x = nf->parameters();
    int info;
    Eigen::VectorXd x(nf->nbParameters());
    for(int i=0; i<nf->nbParameters(); ++i)
    {
        x[i] = nf_x[i];
    }


    EigenFunctor functor(nf, d);
    Eigen::LevenbergMarquardt<EigenFunctor> lm(functor);

	 info = lm.minimize(x);

    if(info == Eigen::LevenbergMarquardtSpace::ImproperInputParameters)
    {
        std::cerr << "<<ERROR>> incorrect parameters" << std::endl;
        return false;
    }


    vec p(nf->nbParameters());

    for(int i=0; i<p.size(); ++i)
    {
        p[i] = x(i);
    }
    std::cout << "<<INFO>> found parameters: " << p << std::endl;
#ifndef DEBUG
    std::cout << "<<DEBUG>> using " << lm.iterations() << " iterations" << std::endl;
#endif
    nf->setParameters(p);
    return true;

}

void nonlinear_fitter_eigen::set_parameters(const arguments& args)
{
}
		
data*     nonlinear_fitter_eigen::provide_data() const { return NULL; }
function* nonlinear_fitter_eigen::provide_function() const { return NULL; }

Q_EXPORT_PLUGIN2(nonlinear_fitter_eigen, nonlinear_fitter_eigen)
