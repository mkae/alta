#include "fitter.h"

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <unsupported/Eigen/NonLinearOptimization>

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <QTime>

fitter* provide_fitter()
{
    return new nonlinear_fitter_eigen();
}

struct EigenFunctor
{
    EigenFunctor(nonlinear_function* f, const data* d) : _f(f)
	{
	}

    inline int inputs() const { return _f->nbParameters(); }
    inline int values() const { return _f->dimY(); }

	int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& y) const
	{
        // Update the parameters vector
        vec _p(inputs());
        for(int i=0; i<inputs(); ++i)
        {
            _p[i] = x(i);
        }

        vec _x  = _d->get(0);
        vec _di = vec(_f->dimY());
        for(int i=0; i<_f->dimY(); ++i)
            _di[i] = _x[_f->dimX() + i];

        vec _y = (*_f)(_x) - _di;
		for(int i=0; i<values(); ++i) y(i) = _y[i];

        std::cout << "f(x) = " << y << std::endl;
		return 0;
	}

	int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) const
	{
		vec _x(inputs());
		for(int i=0; i<inputs(); ++i) _x[i] = x(i);

		vec _jac = _f->parametersJacobian(_x);
		for(int j=0; j<_f->nbParameters(); ++j)
			for(int i=0; i<values(); ++i)
			{
				fjac(i,j) = _jac[values()*j + i];
			}

        std::cout << "J = " <<  fjac << std::endl;

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

    if(dynamic_cast<nonlinear_function*>(fit) == NULL)
    {
        std::cerr << "<<ERROR>> the function is not a non-linear function" << std::endl;
        return false;
    }
    nonlinear_function* nf = dynamic_cast<nonlinear_function*>(fit);

    /* the following starting values provide a rough fit. */
    int info;
    Eigen::VectorXd x;
    x.setConstant(nf->nbParameters(), 1.);


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
    nf->setParameters(p);
    return true;

}

void nonlinear_fitter_eigen::set_parameters(const arguments& args)
{
}
		
data*     nonlinear_fitter_eigen::provide_data() const { return NULL; }
function* nonlinear_fitter_eigen::provide_function() const { return NULL; }

Q_EXPORT_PLUGIN2(nonlinear_fitter_eigen, nonlinear_fitter_eigen)
