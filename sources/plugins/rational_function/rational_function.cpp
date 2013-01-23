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
double rational_function::operator()(double x) const 
{
	double p = 0.0f ;
	double q = 0.0f ;

	for(int i=a.size()-1; i>=0; --i)
	{
		p = x*p + a[i] ;
	}

	for(int i=b.size()-1; i>=0; --i)
	{
		q = x*q + b[i] ;
	}

	return p/q ;
}
		
// Get the p_i and q_j function
double rational_function::p(double x, int i) const
{
	return pow(x, i) ;
}
double rational_function::q(double x, int j) const 
{
	return pow(x, j) ;
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
