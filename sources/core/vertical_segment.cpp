#include "vertical_segment.h"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>
#include <cassert>

void vertical_segment::load(const std::string& filename) 
{
	arguments args ;
	load(filename, args) ;
}
		
void vertical_segment::load(const std::string& filename, const arguments& args) 
{
	std::ifstream file(filename.c_str()) ;
	if(!file.is_open())
	{
		std::cerr << "<<ERROR>> unable to open file \"" << filename << "\"" << std::endl ;
		throw ;
	}

	vec min, max ;

	_nX = 0 ; _nY = 0 ;
	std::vector<int> vs ; int current_vs = 0 ;

	double x, y, dy ;
	while(file.good())
	{
		std::string line ;
		std::getline(file, line) ;
		std::stringstream linestream(line) ;
		
		// Discard incorrect lines
		if(linestream.peek() == '#')
		{
			linestream.ignore(1) ;

			std::string comment ;
			linestream >> comment ;

			if(comment == std::string("DIM"))
			{
				linestream >> _nX >> _nY ;

				vs.assign(dimY(), 0) ;
				for(int k=0; k<dimY(); ++k)
				{
					vs[k] = 0 ;
				}

				_min.resize(dimX()) ;
				_max.resize(dimX()) ;
				
				min = args.get_vec("min", _nX, -std::numeric_limits<double>::max()) ;
				max = args.get_vec("max", _nX,  std::numeric_limits<double>::max()) ;

				for(int k=0; k<dimX(); ++k)
				{
					_min[k] =  std::numeric_limits<double>::max() ;
					_max[k] = -std::numeric_limits<double>::max() ;
				}
			}
			else if(comment == std::string("VS"))
			{
				int t ;
				linestream >> t ;
				vs[current_vs] = t ; ++current_vs ;
			}
			continue ;
		} 
		else if(line.empty())
		{
			continue ;
		}

		vec v(dimX() + 3*dimY()) ;
		for(int i=0; i<dimX(); ++i)
			linestream >> v[i] ;

		for(int i=0; i<dimY(); ++i)
			linestream >> v[dimX() + i] ;

		for(int i=0; i<dimY(); ++i)
		{
			// TODO, the firts case does not account for the
			// dimension of the ouput vector
			if(vs[i] == 2) 
			{
				linestream >> v[dimX() + dimY()+i] ;
				linestream >> v[dimX() + 2*dimY()+i] ;
			} 
			else if(vs[i] == 1)
			{
				double dt ;
				linestream >> dt ;
				v[dimX() + dimY()+i] = v[dimX() + i] * (1.0 - dt) ;
				v[dimX() + 2*dimY()+i] = v[dimX() + i] * (1.0 + dt) ;
			}
			else 
			{
				// TODO Specify the delta in case
				// Handle multiple dim
				double dt = args.get_float("dt", 0.1);
				v[dimX() +   dimY()+i] = v[dimX() + i] * (1.0 - dt) ;
				v[dimX() + 2*dimY()+i] = v[dimX() + i] * (1.0 + dt) ;
			}
		}
		
		// If data is not in the interval of fit
		// TODO: Update to more dims
		bool is_in = true ;
		for(int i=0; i<dimX(); ++i)
		{
			if(v[i] < min[i] || v[i] > max[i])
			{
				is_in = false ;
			}
		}
		if(!is_in)
		{
			continue ;
		}

		_data.push_back(v) ;

		// Update min and max
		for(int k=0; k<dimX(); ++k)
		{
			_min[k] = std::min(_min[k], v[k]) ;
			_max[k] = std::max(_max[k], v[k]) ;
		}
	}

	std::cout << "<<INFO>> loaded file \"" << filename << "\"" << std::endl ;
	std::cout << "<<INFO>> data inside " << _min << " ... " << _max << std::endl ;
	std::cout << "<<INFO>> loading data file of R^" << dimX() << " -> R^" << dimY() << std::endl ;
	std::cout << "<<INFO>> " << _data.size() << " elements to fit" << std::endl ;
}

void vertical_segment::get(int i, vec& x, vec& yl, vec& yu) const
{
#ifdef DEBUG
    assert(i >= 0 && i < _data.size());
#endif
    x.resize(dimX()); yl.resize(dimY()) ; yu.resize(dimY()) ;
    for(int j=0; j<dimX(); ++j)
    {
        x[j] = _data[i][j] ;
    }
    for(int j=0; j<dimY(); ++j)
    {
        yl[j] = _data[i][dimX() + dimY() + j] ;
        yu[j] = _data[i][dimX() + 2*dimY() + j] ;
    }

}
		
void vertical_segment::get(int i, vec& yl, vec& yu) const
{
	yl.resize(dimY()) ; yu.resize(dimY()) ;
	for(int j=0; j<dimY(); ++j)	
	{
		yl[j] = _data[i][dimX() + dimY() + j] ;
		yu[j] = _data[i][dimX() + 2*dimY() + j] ;
	}
}
		
vec vertical_segment::operator[](int i) const
{
	return _data[i] ;
}
vec vertical_segment::get(int i) const 
{
	return _data[i] ;
}

int vertical_segment::size() const
{
	return _data.size() ;
}

vec vertical_segment::min() const 
{
	return _min ;
}

vec vertical_segment::max() const 
{
	return _max ;
}
