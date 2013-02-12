#include "rational_function.h"

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

rational_function_chebychev::rational_function_chebychev() : rational_function() 
{
}

rational_function_chebychev::rational_function_chebychev(const std::vector<double>& a,
                         const std::vector<double>& b) : rational_function(a, b)
{
}
rational_function_chebychev::~rational_function_chebychev()
{
}
		

double chebychev(double x, int i)
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
		return 2.0*x*chebychev(x, i-1) - chebychev(x, i-2) ;
	}
}

// Get the p_i and q_j function
double rational_function_chebychev::p(const vec& x, int i) const
{
	std::vector<int> deg = index2degree(i);
	double res = 1.0;
	for(int k=0; k<dimX(); ++k)
	{
		res *= chebychev(x[k], deg[k]);
	}

	return res ;
}
double rational_function_chebychev::q(const vec& x, int i) const 
{
	std::vector<int> deg = index2degree(i);
	double res = 1.0; 
	for(int k=0; k<dimX(); ++k)
	{
		res *= chebychev(x[k], deg[k]);
	}

	return res ;
}


Q_EXPORT_PLUGIN2(rational_function_chebychev, rational_function_chebychev)
