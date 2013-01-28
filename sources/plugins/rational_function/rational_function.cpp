#include "rational_function.h"

#include <boost/regex.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>


rational_function::rational_function() 
{
}

rational_function::rational_function(const std::vector<double>& a,
                         const std::vector<double>& b) :
	a(a), b(b)
{
}
rational_function::~rational_function()
{
}

// Overload the function operator
vec rational_function::operator()(const vec& x) const 
{
	vec res ;
	res.resize(dimY()) ;

	for(int k=0; k<_nY; ++k)
	{
		double p = 0.0f ;
		double q = 0.0f ;
		for(int l=0; l<_nX; ++l)
		{

			for(int i=a.size()-1; i>=0; --i)
			{
				p = x[l]*p + a[i] ;
			}

			for(int i=b.size()-1; i>=0; --i)
			{
				q = x[l]*q + b[i] ;
			}
		}
		res[k] = p/q ;
	}
	return res ;
}


// Get the p_i and q_j function
double rational_function::p(const vec& x, int i) const
{
	std::vector<int> deg ; deg.assign(dimY(), 0) ;
	
	double res = 1.0 ;

	if(i == 0)
		return res ;

	int temp_i = i ;
	int temp_c ;
	while(temp_i != 0)
	{
		temp_c = (temp_i-1) % dimX() ;
		temp_i = (temp_i - temp_c) / dimX() ;

		deg[temp_c] += 1 ;
	}

	for(int k=0; k<dimX(); ++k)
	{
		res *= pow(x[k], deg[k]) ;
	}

	return res ;
}
double rational_function::q(const vec& x, int i) const 
{
	std::vector<int> deg ; deg.assign(dimY(), 0) ;
	
	double res = 1.0 ;

	if(i == 0)
		return res ;

	int temp_i = i ;
	int temp_c ;
	while(temp_i != 0)
	{
		temp_c = (temp_i-1) % dimX() ;
		temp_i = (temp_i - temp_c) / dimX() ;

		deg[temp_c] += 1 ;
	}

	for(int k=0; k<dimX(); ++k)
	{
		res *= pow(x[k], deg[k]) ;
	}

	return res ;
}

// IO function to text files
void rational_function::load(const std::string& filename)
{
}
void rational_function::save() const
{
}

std::ostream& operator<< (std::ostream& out, const rational_function& r) 
{
	std::cout << "p = [" ;
	for(int i=0; i<r.a.size(); ++i)
	{
		if(i != 0)
		{
			std::cout << ", " ;
		}
		std::cout << r.a[i] ;
	}
	std::cout << "]" << std::endl ;

	std::cout << "q = [" ;
	for(int i=0; i<r.b.size(); ++i)
	{
		if(i != 0)
		{
			std::cout << ", " ;
		}
		std::cout << r.b[i] ;
	}
	std::cout << "]" << std::endl ;

}


Q_EXPORT_PLUGIN2(rational_function, rational_function)
