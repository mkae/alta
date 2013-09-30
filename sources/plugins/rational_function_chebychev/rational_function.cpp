#include "rational_function.h"

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <core/common.h>

ALTA_DLL_EXPORT function* provide_function()
{
    return new rational_function_chebychev();
}

rational_function_chebychev::rational_function_chebychev() : rational_function() 
{
}

rational_function_chebychev::~rational_function_chebychev()
{
}


rational_function_chebychev_1d::rational_function_chebychev_1d()
{
}

rational_function_chebychev_1d::rational_function_chebychev_1d(int np, int nq) :
    rational_function_1d(np, nq)
{
}

rational_function_chebychev_1d::rational_function_chebychev_1d(const vec& a,
                                                               const vec& b) :
    rational_function_1d(a, b)
{
}
		

double chebychev(double x, int i)
{
#ifdef RECURSIVE_FORM
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
		return 2.0*x*chebychev(x, i-1) - chebychev(x, i-2) ;
	}
#else
    return cos(i * acos(x));
#endif
}

// Get the p_i and q_j function
double rational_function_chebychev_1d::p(const vec& x, int i) const
{
	std::vector<int> deg = index2degree(i);
	double res = 1.0;
	for(int k=0; k<dimX(); ++k)
	{
		double xk = 2.0*((x[k] - _min[k]) / (_max[k]-_min[k]) - 0.5);
		res *= chebychev(xk, deg[k]);
	}

	return res ;
}
double rational_function_chebychev_1d::q(const vec& x, int i) const
{
	std::vector<int> deg = index2degree(i);
	double res = 1.0; 
	for(int k=0; k<dimX(); ++k)
	{
		double xk = 2.0*((x[k] - _min[k]) / (_max[k]-_min[k]) - 0.5);
		res *= chebychev(xk, deg[k]);
	}

	return res ;
}
		
rational_function_1d* rational_function_chebychev::get(int i)
{
	// Check for consistency in the index of color channel
	if(i < _nY)
	{
		if(rs[i] == NULL)
		{
			rs[i] = new rational_function_chebychev_1d(np, nq);
			rs[i]->setDimX(dimX());
			rs[i]->setDimY(dimY());
			rs[i]->setMin(getMin()) ;
			rs[i]->setMax(getMax()) ;
		}
		return rs[i];
	}
	else
	{
		std::cout << "<<ERROR>> tried to access out of bound 1D RF" << std::endl;
		return NULL;
	}
}
