/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "data.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include <core/vertical_segment.h>

mxArray *X, *Y, *x, *y;

#define BUFFER_SIZE 10000
char* output = new char[BUFFER_SIZE+1];

data_interpolant::data_interpolant()
{
	_data = new vertical_segment();

	// Create matlab engine
#ifdef WIN32
    if (!(ep = engOpen(NULL)))
#else
    if (!(ep = engOpen("matlab -nosplash")))
#endif
	{
		std::cerr << "<ERROR>> can't start MATLAB engine" << std::endl ;
	}

	output[BUFFER_SIZE] = '\0';
	engOutputBuffer(ep, output, BUFFER_SIZE) ;
}

data_interpolant::~data_interpolant()
{
	delete _data;
	delete[] output;

	mxDestroyArray(x);
	mxDestroyArray(y);
	mxDestroyArray(X);
	mxDestroyArray(Y);
	engClose(ep); 	
}

// Load data from a file
void data_interpolant::load(const std::string& filename) 
{
	// Load the data
	_data->load(filename);

	// Copy the informations
	setDimX(_data->dimX());
	setDimY(_data->dimY());
	setMin(_data->min());
	setMax(_data->max());

	Eigen::MatrixXd eX(_data->size(), dimX()*dimY());
	Eigen::MatrixXd eY(_data->size(), dimY());
	for(int i=0; i<_data->size(); ++i)
	{
		vec x = _data->get(i);

		// For a point extract all its coordinates
		for(int j=0; j<dimX(); ++j)
		{
			eX(i, j) =  x[j];
		}

		// Extract all its values
		for(int k=0; k<dimY(); ++k)
		{
			eY(i, k) = x[dimX() + k];
		}
	}

	// Create matrices
	X = mxCreateDoubleMatrix(_data->size(), dimX(), mxREAL);
	memcpy((void *)mxGetPr(X), (void *) eX.data(),  _data->size()*dimY()*dimX()*sizeof(double));
	engPutVariable(ep, "X", X);
	
	Y = mxCreateDoubleMatrix(_data->size(), dimY(), mxREAL);
	memcpy((void *)mxGetPr(Y), (void *) eY.data(),  _data->size()*dimY()*sizeof(double));
	engPutVariable(ep, "Y", Y);
	
	x = mxCreateDoubleMatrix(1, dimX(), mxREAL);
}
void data_interpolant::load(const std::string& filename, const arguments&)
{
	load(filename);
}

void data_interpolant::save(const std::string& filename) const 
{
}

// Acces to data
vec data_interpolant::get(int id) const 
{
	vec res(dimX() + dimY()) ;
	return res ;
}
vec data_interpolant::operator[](int i) const 
{
	return get(i) ;
}

//! \todo Test this function
void data_interpolant::set(vec x)
{
	assert(x.size() == dimX());
}

vec data_interpolant::value(vec, vec) const
{
	vec res(dimY());
	std::cerr << "<<ERROR>> Deprecated function: " << __func__ << std::endl;
	return res;
}
vec data_interpolant::value(vec ax) const
{
	vec res = vec::Zero(dimY());
	
	// Copy the input vector to matlab code
	memcpy((void *)mxGetPr(x), (void *)&ax[0],  dimX()*sizeof(double));
	engPutVariable(ep, "x", x);

	std::stringstream cmd;
	
	// Iterate over the output dimension
	for(int i=0; i<dimY(); ++i)
	{
		// Evaluate the matlab routine
		if(dimX() == 2)
		{
			cmd << "y(" << i+1 << ") = griddata(X(:,1), X(:,2), Y(:," << i+1 << "), x(1), x(2), 'cubic');";
		}
		else if(dimX() == 3)
		{
			cmd << "y(" << i+1 << ") = griddata(X(:,1), X(:,2), X(:,3), Y(:," << i+1 << "), x(1), x(2), x(3), 'cubic');";
		}
		else
		{
			cmd << "y(" << i+1 << ") = griddatan(X(:, 1:" << dimX() <<"), Y(:, " << i+1 << "), x, 'linear');";
		}
		engEvalString(ep, cmd.str().c_str());
#ifdef DEBUG
		std::cout << output;
#endif
	}

	// Get results and copy it
	y = engGetVariable(ep, "y") ;
	double* y_val = (double*)mxGetData(y);
	memcpy(&res[0], y_val, dimY()*sizeof(double));

	// Fail safe: if the query point is outside of the convex hull of the
	// data, rerun the algorithm using a nearest option.
	if(isnan(res[0]))
	{
		for(int i=0; i<dimY(); ++i)
		{
			cmd.str(std::string());
			cmd << "y("<< i+1 <<") = griddatan(X(:, 1:" << dimX() <<"), Y(:, " << i+1 << "), x, 'nearest');";
			engEvalString(ep, cmd.str().c_str());
		}
		// Get results and copy it
		y = engGetVariable(ep, "y") ;
		double* y_val = (double*)mxGetData(y);
		memcpy(&res[0], y_val, dimY()*sizeof(double));
	}

   return res;
}

// Get data size, e.g. the number of samples to fit
int data_interpolant::size() const 
{
	assert(_data != NULL);
	return _data->size();
}

ALTA_DLL_EXPORT data* provide_data()
{
    return new data_interpolant();
}


