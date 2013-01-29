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
vec rational_function::value(const vec& x) const 
{
	vec res ;
	res.reserve(_nY) ;

	for(int k=0; k<_nY; ++k)
	{
		double p = 0.0f ;
		double q = 0.0f ;
		
		for(int i=a.size()-1; i>=0; --i)
		{
			p += a[i]*this->p(x, i) ;
		}

		for(int i=b.size()-1; i>=0; --i)
		{
			q += b[i]*this->q(x, i) ;
		}

		res[k] = p/q ;
	}
	return res ;
}


// Get the p_i and q_j function
double rational_function::p(const vec& x, int i) const
{
	if(dimX() == 1)
	{
		return pow(x[0], i) ;
	}
	
	std::vector<int> deg ; deg.assign(dimX(), 0) ;
	
	double res = 1.0 ;

	if(i == 0)
		return res ;

	int temp_i = i ;
	int temp_c ;
	while(temp_i > 1)
	{
		temp_c = temp_i % dimX() ;
		temp_i = (temp_i - temp_c) / dimX() ;

		deg[temp_c] += 1 ;
	}
	deg[0] += temp_i ;

	for(int k=0; k<dimX(); ++k)
	{
		res *= pow(x[k], deg[k]) ;
	}

	return res ;
}
double rational_function::q(const vec& x, int i) const 
{
	if(dimX() == 1)
	{
		return pow(x[0], i) ;
	}

	std::vector<int> deg ; deg.assign(dimX(), 0) ;
	
	double res = 1.0 ;

	if(i == 0)
		return res ;

	int temp_i = i ;
	int temp_c ;
	while(temp_i > 1)
	{
		temp_c = temp_i % dimX() ;
		temp_i = (temp_i - temp_c) / dimX() ;

		deg[temp_c] += 1 ;
	}
	deg[0] += temp_i ;

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
void rational_function::save(const std::string& filename, const arguments& args) const
{
	vec min, max ;
	min.assign(_nX, args.get_float("min", 0.0f)) ;
	max.assign(_nX, args.get_float("max", 1.5f)) ;

	int nb_samples = args.get_int("nb_samples", 100) ;
	double dt = (max[0]-min[0]) / nb_samples ;

	std::ofstream file(filename.c_str(), std::ios_base::trunc);
	for(double x=min[0]; x<=max[0]; x+=dt)
	{
		vec vx ; for(int i=0;i<_nX; ++i) { vx.push_back(x) ; }
		file << x << "\t" << value(vx)[0] << std::endl ;
		std::cout << x << "\t" << value(vx)[0] << std::endl ;
	}
}

std::ostream& operator<< (std::ostream& out, const rational_function& r) 
{
	std::cout << "p = [" ;
	for(unsigned int i=0; i<r.a.size(); ++i)
	{
		if(i != 0)
		{
			std::cout << ", " ;
		}
		std::cout << r.a[i] ;
	}
	std::cout << "]" << std::endl ;

	std::cout << "q = [" ;
	for(unsigned int i=0; i<r.b.size(); ++i)
	{
		if(i != 0)
		{
			std::cout << ", " ;
		}
		std::cout << r.b[i] ;
	}
	std::cout << "]" << std::endl ;

	return out ;
}


Q_EXPORT_PLUGIN2(rational_function, rational_function)
