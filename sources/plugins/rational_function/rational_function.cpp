#include "rational_function.h"

#include <boost/regex.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>


rational_function::rational_function() : a(), b()
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
		
void rational_function::update(const std::vector<double>& in_a,
                               const std::vector<double>& in_b)
{
	a.reserve(in_a.size()) ;
	b.reserve(in_b.size()) ;
	a = in_a ;
	b = in_b ;
}

// Overload the function operator
vec rational_function::value(const vec& x) const 
{
	vec res(_nY) ;

	const int np = a.size() / _nY ;
	const int nq = b.size() / _nY ;

	for(int k=0; k<_nY; ++k)
	{
		double p = 0.0f ;
		double q = 0.0f ;
		
		for(unsigned int i=0; i<np; ++i)
		{
			p += a[i*_nY + k]*this->p(x, i) ;
		}

		for(unsigned int i=0; i<nq; ++i)
		{
			q += b[i*_nY + k]*this->q(x, i) ;
		}

		res[k] = p/q ;
	}
	return res ;
}

std::vector<int> rational_function::index2degree(int i) const
{
	std::vector<int> deg ; deg.assign(dimX(), 0) ;
	if(dimX() == 1)
	{
		deg[0] = i ;
		return deg ;
	}
	
	int temp_i = i-1 ;
	int temp_c ;
	while(temp_i > 0)
	{
		temp_c = temp_i % dimX() ;
		temp_i = temp_i / dimX() ;

		deg[temp_c] += 1 ;
	}
//	deg[0] += temp_i ;
	return deg ;
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
	save_rational_function(filename) ;
}

void rational_function::save_gnuplot(const std::string& filename, const data* d, const arguments& args) const 
{
	std::ofstream file(filename.c_str(), std::ios_base::trunc);
	for(int i=0; i<d->size(); ++i)
	{
		vec v = d->get(i) ;
		vec y1 ; y1.assign(d->dimY(), 0.0) ;
		for(int k=0; k<d->dimY(); ++k) { y1[k] = v[d->dimX() + k] ; }

		vec y2 = value(v) ;
		for(int u=0; u<d->dimX(); ++u)
			file << v[u] << "\t" ;

		for(int u=0; u<d->dimY(); ++u)
			file << y2[u] << "\t" ;

		file << std::endl ;
	}
}

void rational_function::save_rational_function(const std::string& filename) const 
{
	std::ofstream file(filename.c_str(), std::ios_base::trunc);
	file << "#DIM " << _nX << " " << _nY << std::endl ;
	file << "#NP " << a.size() << std::endl ;
	file << "#NQ " << b.size() << std::endl ;
	file << "#BASIS poly" << std::endl ;

	for(unsigned int i=0; i<a.size(); ++i)
	{
		std::vector<int> index = index2degree(i) ;
		for(unsigned int j=0; j<index.size(); ++j)
		{
			file << index[j] << "\t" ;
		}
		file << a[i] << std::endl ;
	}
	
	for(unsigned int i=0; i<b.size(); ++i)
	{
		std::vector<int> index = index2degree(i) ;
		for(unsigned int j=0; j<index.size(); ++j)
		{
			file << index[j] << "\t" ;
		}
		file << b[i] << std::endl ;
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