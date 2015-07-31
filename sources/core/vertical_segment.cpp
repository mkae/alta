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
	: data(dim_X, dim_Y), _is_absolute(true), _dt(0.1) 
{
	initializeToZero( size );
}

vertical_segment::vertical_segment( params::input in_param, 
                                   	params::output out_param,
                                   	unsigned int size )
	: data( in_param, out_param), _is_absolute(true), _dt(0.1)
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

	arguments header = arguments::parse_header(file);

	// Default behaviour: parsing a file as TEXT file. Send a message error in case
	// the user did not set it.
	if(! header.is_defined("FORMAT")) {
		std::cerr << "<<DEBUG>> The file format is undefined, assuming TEXT" << std::endl;
	}

	if (header["FORMAT"] == "binary") {

#ifdef _WIN32
      // On MS Windows strange things happen when the ifstream is not opened in
      // binary mode.  To avoid having badbit set after a few entries read, we
      // close and re-open the file with correct flags.
      file.close();
      file.open(filename.c_str(), std::ifstream::binary);
      header = arguments::parse_header(file);
#endif

		load_data_from_binary(file, header, *this);
	} else {
		load_data_from_text(file, header, *this, args);
	}

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
   // Check if the input data 'x' has the size of a vertical segment (i.e. dimX+3*dimY),
   // if it only has the size of a single value, then create the segment.
   if(x.size() == dimX() + 3*dimY()) {
      _data.push_back(x);

   } else if(x.size() == dimX() + dimY()) {
      vec y(dimX() + 3*dimY());
      _data.push_back(vs(x));

   } else {
      std::cerr << "<<ERROR>> Passing an incorrect element to vertical_segment::set" << std::endl;
      throw;
   }
}

void vertical_segment::set(int i, const vec& x)
{
   // Check if the input data 'x' has the size of a vertical segment (i.e. dimX+3*dimY),
   // if it only has the size of a single value, then create the segment.
   if(x.size() == dimX() + 3*dimY()) {
      _data[i] = x;

   } else if(x.size() == dimX() + dimY()) {
      _data[i] = vs(x);

   } else {
      std::cerr << "<<ERROR>> Passing an incorrect element to vertical_segment::set" << std::endl;
      throw;
   }
}

int vertical_segment::size() const
{
	return _data.size() ;
}

vec vertical_segment::vs(const vec& x) const {
   vec y(dimX() + 3*dimY());
   for(unsigned int i=0; i<dimX()+dimY(); ++i) {
      y[i] = x[i];
   }
   
   for(unsigned int i=0; i<dimY(); ++i) {
      const double val = x[i + dimX()];
      y[i + dimX()+1*dimY()] = val - ((_is_absolute) ? _dt : _dt*val);
      y[i + dimX()+2*dimY()] = val + ((_is_absolute) ? _dt : _dt*val);
   }

   return y;
}
