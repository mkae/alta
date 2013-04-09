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

struct EigenFunctor
{
	EigenFunctor(nonlinear_function* f) : _f(f)
	{
	}

	inline int inputs() const { return _f->dimX(); }
	inline int values() const { return _f->dimY(); }

	int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& y) const
	{
		vec _x(inputs());
		for(int i=0; i<inputs(); ++i) _x[i] = x(i);

		vec _y = (*_f)(_x);
		for(int i=0; i<values(); ++i) y(i) = _y[i];

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

		return 0;
	}


	nonlinear_function* _f;
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
		std::cerr << "<<ERROR>> the function is non a non-linear function" << std::endl;
		return false;
	}
	nonlinear_function* nf = dynamic_cast<nonlinear_function*>(fit);

	/* the following starting values provide a rough fit. */
	int info;
	Eigen::VectorXd x;
	x.setConstant(nf->parameters().size(), 1.);

	EigenFunctor functor(nf);
	Eigen::LevenbergMarquardt<EigenFunctor> lm(functor);
	info = lm.minimize(x);

	vec p(nf->parameters().size());
	for(int i=0; i<p.size(); ++i)
	{
		p[i] = x(i);
	}
	nf->setParameters(p);

	return info == 1 ;
}

void nonlinear_fitter_eigen::set_parameters(const arguments& args)
{
}
		
data*     nonlinear_fitter_eigen::provide_data() const { return NULL; }
function* nonlinear_fitter_eigen::provide_function() const { return NULL; }

Q_EXPORT_PLUGIN2(nonlinear_fitter_eigen, nonlinear_fitter_eigen)
