#include "rational_function.h"

#include <core/common.h>
#include <core/rational_function.h>

ALTA_DLL_EXPORT function* provide_function()
{
    return new rational_function_legendre();
}

rational_function_legendre_1d::rational_function_legendre_1d()
{
}

rational_function_legendre_1d::rational_function_legendre_1d(int np, int nq, params::input param) :
    rational_function_1d(np, nq)
{
	setParametrization(param);
}

rational_function_legendre_1d::rational_function_legendre_1d(const vec& a, 
                                                             const vec& b) :
    rational_function_1d(a, b)
{
}

double rational_function_legendre_1d::legendre(double x, int i) const
{
	if(i == 0)
	{
		return 1;
	}
	else if(i == 1)
	{
		return x;
	}
	else
	{
		return ((2*i-1)*x*legendre(x, i-1) - (i-1)*legendre(x, i-2)) / (double)i ;
	}
}

// Get the p_i and q_j function
double rational_function_legendre_1d::p(const vec& x, int i) const
{
	std::vector<int> deg = index2degree(i);
	double res = 1.0;
	for(int k=0; k<dimX(); ++k)
	{
		res *= legendre(2.0*((x[k] - _min[k]) / (_max[k]-_min[k]) - 0.5), deg[k]);
	}

	// Apply cosine factor to the result if a parametrization was
	// set. I apply both cos(theta_l) and cos(theta_v).
	if(input_parametrization() != params::UNKNOWN_INPUT)
	{
		double cart[6];
		params::convert(&x[0], input_parametrization(), params::CARTESIAN, cart);

		res *= cart[2]*cart[5];
	}

	return res ;
}
double rational_function_legendre_1d::q(const vec& x, int i) const 
{
	std::vector<int> deg = index2degree(i);
	double res = 1.0;
	for(int k=0; k<dimX(); ++k)
	{
		res *= legendre(2.0*((x[k] - _min[k]) / (_max[k]-_min[k]) - 0.5), deg[k]);
	}
	
	return res ;
}


rational_function_legendre::rational_function_legendre()
{
}

rational_function_legendre::~rational_function_legendre()
{
}

