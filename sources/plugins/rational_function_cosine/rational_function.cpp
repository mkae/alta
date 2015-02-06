/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2014 CNRS
   Copyright (C) 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

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

rational_function_legendre_1d::rational_function_legendre_1d(int nX, int np, int nq, params::input param) :
    rational_function_1d(nX, np, nq)
{
	setParametrization(param);
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
	else if (i == 2)
	{
		return (3*x*x - 1) * 0.5;
	}
	else if (i == 3)
	{
		return ( (5*x*x - 3) * x ) * 0.5;
	}
	else if (i == 4)
	{
		return ((7*x*x - 6) * 5*x*x + 3) * 0.125;
	}
	else if (i== 5)
	{
		return  (((63*x*x - 70)*x*x + 15)*x)*0.125; 
	}
	else if (i== 6)
	{
		return (((231*x*x-315)*x*x + 105)*x*x - 5) * 0.0625;
	}
	else if (i==7)
	{
		return ((((429*x*x-693)*x*x + 315)*x*x - 35)*x) * 0.0625;
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

