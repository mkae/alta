/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "vertical_segment.h"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>
#include <cassert>

//#define RELATIVE_ERROR

vertical_segment::vertical_segment(unsigned int dim_X, 
                                   unsigned int dim_Y, 
                                   unsigned int size)
	: data(dim_X, dim_Y)
{
	initializeToZero( size );
}

vertical_segment::vertical_segment( params::input in_param, 
                                   	params::output out_param,
                                   	unsigned int size )
	: data( in_param, out_param)
{
	initializeToZero( size );
}

void 
vertical_segment::initializeToZero( unsigned int number_of_data_elements )
{
	_data.clear();

	unsigned int const  size_of_one_element = dimX() + dimY();
	
	for( unsigned int i=0; i < number_of_data_elements; i++)
	{
		vec initial_element = vec::Zero( size_of_one_element );
		_data.push_back( initial_element );
	}
}



void vertical_segment::load(const std::string& filename) 
{
	arguments args ;
	load(filename, args) ;
}

void vertical_segment::load(const std::string& filename, const arguments& args) 
{
	std::ifstream file;

	// Raise an exception when 'open' fails.
	file.exceptions (std::ios::failbit);
	file.open(filename.c_str());
	file.exceptions (std::ios::goodbit);

	load_data_from_text(file, *this, args);

	file.close();
}

void vertical_segment::load_data_from_text(std::istream& input,
																					 vertical_segment& result,
																					 const arguments& args)
{
	vec min, max ;
	vec ymin, ymax;

	result._nX = 0 ; result._nY = 0 ;
	std::vector<int> vs ; int current_vs = 0 ;
	while(input.good())
	{
		std::string line ;
		std::getline(input, line) ;
		std::stringstream linestream(line) ;

		// Discard incorrect lines
		if(linestream.peek() == '#')
		{
			linestream.ignore(1) ;

			std::string comment ;
			linestream >> comment ;

			if(comment == std::string("DIM"))
			{
				linestream >> result._nX >> result._nY ;

				vs.reserve(result.dimY()) ;
				for(int k=0; k<result.dimY(); ++k)
				{
					vs[k] = 0 ;
				}

				result._min.resize(result.dimX()) ;
				result._max.resize(result.dimX()) ;

				min = args.get_vec("min", result._nX, -std::numeric_limits<float>::max()) ;
				max = args.get_vec("max", result._nX,  std::numeric_limits<float>::max()) ;
#ifdef DEBUG
				std::cout << "<<DEBUG>> data will remove outside of " << min << " -> " << max << " x-interval" << std::endl;
#endif

				ymin = args.get_vec("ymin", result._nY, -std::numeric_limits<float>::max()) ;
				ymax = args.get_vec("ymax", result._nY,  std::numeric_limits<float>::max()) ;
#ifdef DEBUG
				std::cout << "<<DEBUG>> data will remove outside of " << ymin << " -> " << ymax << " y-interval" << std::endl;
#endif
				
				for(int k=0; k<result.dimX(); ++k)
				{
					result._min[k] =  std::numeric_limits<double>::max() ;
					result._max[k] = -std::numeric_limits<double>::max() ;
				}
			}
			else if(comment == std::string("VS"))
			{
				int t ;
				linestream >> t ;
				vs[current_vs] = t ; ++current_vs ;
			}
			else if(comment == std::string("PARAM_IN"))
			{
				std::string param;
				linestream >> param;
				result._in_param = params::parse_input(param);
			}
			else if(comment == std::string("PARAM_OUT"))
			{
				std::string param;
				linestream >> param;
				result._out_param = params::parse_output(param);
			}
			continue ;
		} 
		else if(line.empty())
		{
			continue ;
		}
		else
		{
			// Read the data point x and y coordinates
			vec v = vec::Zero(result.dimX() + 3*result.dimY()) ;
			for(int i=0; i<result.dimX()+result.dimY(); ++i) 
			{
				linestream >> v[i] ;
			}

			// If data is not in the interval of fit
			bool is_in = true ;
			for(int i=0; i<result.dimX(); ++i)
			{
				if(v[i] < min[i] || v[i] > max[i])
				{
					is_in = false ;
				}
			}
			for(int i=0; i<result.dimY(); ++i)
			{
				if(v[result.dimX()+i] < ymin[i] || v[result.dimX()+i] > ymax[i])
				{
					is_in = false ;
				}
			}
			if(!is_in)
			{
				continue ;
			}

//			/*
			// Correction of the data by 1/cosine(theta_L)
			double factor = 1.0;
			if(args.is_defined("data-correct-cosine"))
			{
				double cart[6];
				params::convert(&v[0], result.input_parametrization(), params::CARTESIAN, cart);
				if(cart[5] > 0.0 && cart[2] > 0.0)
				{
					factor = 1.0/cart[5]*cart[2];
					for(int i=0; i<result.dimY(); ++i) 
					{
						v[i + result.dimX()] /= factor;
					}
				}
				else
				{
					continue;
				}
			}
			// End of correction
//			*/

			// Check if the data containt a vertical segment around the mean
			// value.
			for(int i=0; i<result.dimY(); ++i)
			{
				double min_dt = 0.0;
				double max_dt = 0.0;


				if(vs[i] == 2)
				{
					linestream >> min_dt ;
					linestream >> max_dt ;
					min_dt = min_dt-v[result.dimX()+i];
					max_dt = max_dt-v[result.dimX()+i];
				}
				else if(vs[i] == 1)
				{
					double dt ;
					linestream >> dt ;
					min_dt = -dt;
					max_dt =  dt;
				}
				else
				{
					double dt = args.get_float("dt", 0.1f);
					min_dt = -dt;
					max_dt =  dt;
				}

				if(args.is_defined("dt-relative"))
				{
               v[result.dimX() +   result.dimY()+i] = v[result.dimX() + i] * (1.0 + min_dt) ;
					v[result.dimX() + 2*result.dimY()+i] = v[result.dimX() + i] * (1.0 + max_dt) ;
				}
				else if(args.is_defined("dt-max"))
				{
               v[result.dimX() +   result.dimY()+i] = v[result.dimX() + i] + std::max(v[result.dimX() + i] * min_dt, min_dt);
					v[result.dimX() + 2*result.dimY()+i] = v[result.dimX() + i] + std::max(v[result.dimX() + i] * max_dt, max_dt);
				}
				else
				{
					v[result.dimX() +   result.dimY()+i] = v[result.dimX() + i] + min_dt ;
					v[result.dimX() + 2*result.dimY()+i] = v[result.dimX() + i] + max_dt ;
				}

				// You can enforce the vertical segment to stay in the positive
				// region using the --data-positive command line argument. Note
				// that the data point is also clamped to zero if negative.
				if(args.is_defined("dt-positive"))
				{
					v[result.dimX() +          i] = std::max(v[result.dimX() +          i], 0.0);
					v[result.dimX() +   result.dimY()+i] = std::max(v[result.dimX() +   result.dimY()+i], 0.0);
					v[result.dimX() + 2*result.dimY()+i] = std::max(v[result.dimX() + 2*result.dimY()+i], 0.0);
				}

#ifdef DEBUG
                std::cout << "<<DEBUG>> vs = [" << v[result.dimX() +   result.dimY()+i] << ", " << v[result.dimX() + 2*result.dimY()+i] << "]" << std::endl;
#endif
			}

			result._data.push_back(v) ;

			// Update min and max
			for(int k=0; k<result.dimX(); ++k)
			{
				result._min[k] = std::min(result._min[k], v[k]) ;
				result._max[k] = std::max(result._max[k], v[k]) ;
			}
		}
	}
			
	if(args.is_defined("data-correct-cosine"))
		result.save("/tmp/data-corrected.dat");

	std::cout << "<<INFO>> loaded input stream" << std::endl ;
	std::cout << "<<INFO>> data inside " << result._min << " ... " << result._max << std::endl ;
	std::cout << "<<INFO>> loading data input of R^" << result.dimX() << " -> R^" << result.dimY() << std::endl ;
	std::cout << "<<INFO>> " << result._data.size() << " elements to fit" << std::endl ;
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
	//SLOW !!!!! and useless 
	// call _data[i]
    const int n = dimX() + dimY();
    vec res(n);
    
    for(int k=0; k<n; ++k) 
    { 
			res[k] = _data[i][k]; 
    }    
    return res ;
}

//! \todo Check the vertical segment size and if the data
//! is not already present.
void vertical_segment::set(const vec& x)
{
//	assert(x.size() == dimX() + dimY() || x.size() == dimX() + 3*dimY());
	_data.push_back(x);
}

void vertical_segment::set(int i, const vec& x)
{
//	assert(x.size() == dimX() + dimY() || x.size() == dimX() + 3*dimY());
	_data[i] = x;
}

int vertical_segment::size() const
{
	return _data.size() ;
}
