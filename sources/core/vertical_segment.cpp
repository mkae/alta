/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "vertical_segment.h"
#include "data_storage.h"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
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

	header header(file);

	if (header["FORMAT"].string() == "binary")
			load_data_from_binary(file, header, *this);
	else if (header["FORMAT"].string() == "text")
			load_data_from_text(file, header, *this, args);
	else throw;															// FIXME: Throw a usable exception.

	file.close();
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
