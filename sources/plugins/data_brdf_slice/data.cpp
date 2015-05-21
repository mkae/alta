/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "data.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "EXR_IO.h"

data_brdf_slice::data_brdf_slice()
{
	// Allocate data
	width = 512; height = 512;
	slice = 1;
	_data = new double[3*width*height*slice];

    // Set the input and output parametrization
    _in_param  = params::STARK_2D;
    _out_param = params::RGB_COLOR;
    _nX = 2;
    _nY = 3;
}

data_brdf_slice::~data_brdf_slice()
{
	delete[] _data;
}

// Load data from a file
void data_brdf_slice::load(const std::string& filename) 
{
	delete[] _data;
	t_EXR_IO<double>::LoadEXR(filename.c_str(), width, height, _data);
}
void data_brdf_slice::load(const std::string& filename, const arguments&)
{
	load(filename);
}

void data_brdf_slice::save(const std::string& filename) const 
{
	if(!t_EXR_IO<double>::SaveEXR(filename.c_str(), width, slice*height, _data))
	{
		std::cerr << "<<ERROR>> unable to save image file" << std::endl;
	}
}

// Acces to data
vec data_brdf_slice::get(int id) const 
{
	vec res(5) ;
	const int i = id % width;
	const int k = id / (width*height);
	const int j = (id - k*width*height) / width;

	res[0] = (i+0.5) / double(width);
	res[1] = (j+0.5) / double(height);
	//res[2] = M_PI*(k+0.5) / double(slice);

	res[2] = _data[3*id + 0];
	res[3] = _data[3*id + 1];
	res[4] = _data[3*id + 2];

	return res ;
}
vec data_brdf_slice::operator[](int i) const 
{
	return get(i) ;
}

//! \todo Test this function
void data_brdf_slice::set(const vec& x)
{
	assert(x.size() == 5);
	assert(x[0] <= 1.0/*0.5*M_PI*/ && x[0] >= 0.0);
	assert(x[1] <= 1.0/*0.5*M_PI*/ && x[1] >= 0.0);

	const int i  = floor(x[0] * width  / /*(0.5*M_PI)*/ 1.0);
	const int j  = floor(x[1] * height / /*(0.5*M_PI)*/ 1.0);
	const int k  = 0; 
	//const int k  = floor(x[2] * slice  / (M_PI));
	const int id = i + j*width + k*width*height;

	_data[3*id + 0] = x[2];
	_data[3*id + 1] = x[3];
	_data[3*id + 2] = x[4];
}
void data_brdf_slice::set(int id, const vec& x)
{
	assert(x.size() == 3);

	_data[3*id + 0] = x[0];
	_data[3*id + 1] = x[1];
	_data[3*id + 2] = x[2];
}

vec data_brdf_slice::value(const vec& x) const
{
	assert(x[0] <= /*0.5*M_PI*/ 1.0 && x[0] >= 0.0);
	assert(x[1] <= /*0.5*M_PI*/ 1.0 && x[1] >= 0.0);

	const int i  = floor(x[0] * width  / /*(0.5*M_PI)*/ 1.0);
	const int j  = floor(x[1] * height / /*(0.5*M_PI)*/ 1.0);
	const int k  = 1; 
	//const int k  = floor(x[2] * slice  / (M_PI));
	const int id = (i + j*width)*k;

	if(i < 0 || i >= width)  { std::cerr << "<<ERROR>> out of bounds: " << x << std::endl; }
	if(j < 0 || j >= height) { std::cerr << "<<ERROR>> out of bounds: " << x << std::endl; }

	vec res(3);
	res[0] = _data[3*id + 0];
	res[1] = _data[3*id + 1];
	res[2] = _data[3*id + 2];
   return res;
}

// Get data size, e.g. the number of samples to fit
int data_brdf_slice::size() const 
{
	return width*height*slice;
}

// Get min and max input space values
vec data_brdf_slice::min() const 
{
	vec res(2);
	res[0] = 0.0 ;
	res[1] = 0.0 ;
	//res[2] = 0.0 ;
	return res ;
}
vec data_brdf_slice::max() const
{
	vec res(2);
	res[0] = M_PI / 2 ;
	res[1] = M_PI / 2 ;
	//res[2] = M_PI;
	return res ;
}

int data_brdf_slice::dimX() const 
{ 
	return 2 ; 
}
int data_brdf_slice::dimY() const 
{ 
	return 3;
}

ALTA_DLL_EXPORT data* provide_data()
{
    return new data_brdf_slice();
}


