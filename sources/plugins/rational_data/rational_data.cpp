#include "rational_data.h"

#include <boost/regex.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

void rational_data::load(const std::string& filename) 
{
	load(filename, -std::numeric_limits<double>::max(), std::numeric_limits<double>::max()) ;

}
		
void rational_data::load(const std::string& filename, double min, double max) 
{
	std::ifstream file(filename) ;
	_min =  std::numeric_limits<double>::max() ;
	_max = -std::numeric_limits<double>::max() ;

	if(!file.is_open())
	{
		std::cerr << "<<ERROR>> unable to open file \"" << filename << "\"" << std::endl ;
	}

	// N-Floats regexp
	boost::regex e ("^([0-9]*\.?[0-9]+[\\t ]?)+");

	double x, y, dy ;
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
			dy = 0.01f ;
		}

		if(x <= max && x >= min)
		{
			std::vector<double> v ;
			v.push_back(x) ;
			v.push_back(y-dy) ;
			v.push_back(y+dy) ;
			_data.push_back(v) ;

			// Update min and max
			_min = std::min(_min, x) ;
			_max = std::max(_max, x) ;
		}
	}

	//TODO Test for small data
/*	std::vector<std::vector<double> > temp ;
	std::ofstream temp_out("temp.gnuplot", std::ios_base::trunc) ;
	for(int i=0; i<20; ++i)
	{
		int k = (i * _data.size()) / 20 ;
		temp_out << _data[k][0] << "\t" << _data[k][1] << "\t" << 0.5*(_data[k][2] - _data[k][1]) << std::endl ;
		temp.push_back(_data[k]) ;
	}
	_data = temp ;
*/
	// Sort the vector
	std::sort(_data.begin(), _data.end(), [](const std::vector<double>& a, const std::vector<double>& b){return (a[0]<b[0]);});

	std::cout << "<<INFO>> loaded file \"" << filename << "\"" << std::endl ;
	std::cout << "<<INFO>> data inside [" << _min << ", " << _max << "]" << std::endl ;
	std::cout << "<<INFO>> " << _data.size() << " elements to fit" << std::endl ;
}

bool rational_data::get(int i, double& x, double& yl, double& yu) const
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
		
int rational_data::input_dimension() const 
{
	return 1 ;
}
int rational_data::output_dimension() const
{
	return 1 ;
}

const std::vector<double>& rational_data::operator[](int i) const
{
	return _data[i] ;
}
const std::vector<double>& rational_data::get(int i) const 
{
	return _data[i] ;
}

int rational_data::size() const
{
	return _data.size() ;
}

double rational_data::min() const 
{
	return _min ;
}

double rational_data::max() const 
{
	return _max ;
}

Q_EXPORT_PLUGIN2(rational_data, rational_data)
