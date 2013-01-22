#include "rational_1d_fitter.h"

#include <boost/regex.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>

rational_1d::rational_1d() 
{
}

rational_1d::rational_1d(const std::vector<float>& a,
                         const std::vector<float>& b) :
	a(a), b(b)
{
}
rational_1d::~rational_1d()
{
}

// Overload the function operator
float rational_1d::operator()(float x) const 
{
	float p = 0.0f ;
	float q = 0.0f ;

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
float rational_1d::p(float x, int i) const
{
	return pow(x, i) ;
}
float rational_1d::q(float x, int j) const 
{
	return pow(x, j) ;
}

// IO function to text files
void rational_1d::load(const std::string& filename)
{
}
void rational_1d::save() const
{
}

std::ostream& operator<< (std::ostream& out, const rational_1d& r) 
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

void rational_1d_data::load(const std::string& filename) 
{
	load(filename, -std::numeric_limits<float>::max(), std::numeric_limits<float>::max()) ;

}
		
void rational_1d_data::load(const std::string& filename, float min, float max) 
{
	std::ifstream file(filename) ;
	_min =  std::numeric_limits<float>::max() ;
	_max = -std::numeric_limits<float>::max() ;

	if(!file.is_open())
	{
		std::cerr << "<<ERROR>> unable to open file \"" << filename << "\"" << std::endl ;
	}

	// N-Floats regexp
	boost::regex e ("^([0-9]*\.?[0-9]+[\\t ]?)+");

	float x, y, dy ;
	while(file.good())
	{
		std::string line ;
		std::getline(file, line) ;
		std::stringstream linestream(line) ;
		
		// Discard incorrect lines
		if(!boost::regex_match(line,e))
		{
			continue ;
		}

		linestream >> x >> y ;
		if(linestream.good()) {
			linestream >> dy ;
		} else {
			// TODO Specify the delta in case
			dy = 0.001f ;
		}

		if(x <= max && x >= min)
		{
			std::vector<float> v ;
			v.push_back(x) ;
			v.push_back(y-dy) ;
			v.push_back(y+dy) ;
			_data.push_back(v) ;

			// Update min and max
			_min = std::min(_min, x) ;
			_max = std::max(_max, x) ;
		}
	}

	// Sort the vector
	std::sort(_data.begin(), _data.end(), [](const std::vector<float>& a, const std::vector<float>& b){return (a[0]<b[0]);});

	std::cout << "<<INFO>> loaded file \"" << filename << "\"" << std::endl ;
	std::cout << "<<INFO>> data inside [" << _min << ", " << _max << "]" << std::endl ;
}

bool rational_1d_data::get(int i, float& x, float& yl, float& yu) const
{
	if(i >= (int)_data.size())
	{
		return false ;
	}

	x  = _data[i][0] ;
	yl = _data[i][1] ;
	yu = _data[i][2] ;

	return true ;
}

const std::vector<float>& rational_1d_data::operator[](int i) const
{
	return _data[i] ;
}

int rational_1d_data::size() const
{
	return _data.size() ;
}

float rational_1d_data::min() const 
{
	return _min ;
}

float rational_1d_data::max() const 
{
	return _max ;
}
