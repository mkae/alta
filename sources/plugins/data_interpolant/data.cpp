/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015 Inria

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

//#define USE_DELAUNAY
#ifdef USE_DELAUNAY
#include <CGAL/Cartesian_d.h>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/leda_integer.h>
#include <CGAL/Delaunay_d.h>


typedef CGAL::Homogeneous_d<double> Kernel;
typedef CGAL::Delaunay_d<Kernel> Delaunay_d;
typedef Delaunay_d::Point_d Point;
typedef Delaunay_d::Simplex_handle S_handle;
typedef Delaunay_d::Point_const_iterator P_iter;

Delaunay_d* D;
int dD;
#else
#include <flann/flann.hpp>
#endif

data_interpolant::data_interpolant() : _data(new vertical_segment())
{
	_knn = 3;
}

data_interpolant::~data_interpolant()
{
#ifndef USE_DELAUNAY
	if(_kdtree != NULL)
		delete _kdtree;
#else
	if(D != NULL)
		delete D;
#endif
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
	setParametrization(_data->input_parametrization());
	setParametrization(_data->output_parametrization());

#ifdef USE_DELAUNAY
	dD = dimX()+dimY();
	D  = new Delaunay_d(dD);
	for(int i=0; i<_data->size(); ++i)
	{
		vec x = _data->get(i);

		Point pt(dD, &x[0], &x[dD]);
		D->insert(pt);
	}

	std::cout << "<<DEBUG>> number of points in the Delaunay triangulation: " << D->all_points().size() << std::endl;
	std::cout << "<<DEBUG>> number of points in input: " << _data->size() << std::endl;
#else
	// Update the KDtreee by inserting all points
	double* _d = new double[dimX()*_data->size()];
	flann::Matrix<double> pts(_d, _data->size(), dimX());
	for(int i=0; i<_data->size(); ++i)
	{
		vec x = _data->get(i);
		memcpy(pts[i], &x[0], dimX()*sizeof(double)); 
	}
	_kdtree = new flann::Index< flann::L2<double> >(pts, flann::KDTreeIndexParams(4));
	_kdtree->buildIndex();
#endif
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

//! \todo Test this function
void data_interpolant::set(const vec& x)
{
	NOT_IMPLEMENTED();
}
void data_interpolant::set(int i, const vec& x)
{
	NOT_IMPLEMENTED();	
}

vec data_interpolant::value(const vec& x) const
{
	vec res = vec::Zero(dimY());

#ifndef USE_DELAUNAY
	// Query point
	vec xc(x);
	flann::Matrix<double> pts(&xc[0], 1, dimX());
	std::vector< std::vector<int> >    indices;
	std::vector< std::vector<double> > dists;

	_kdtree->knnSearch(pts, indices, dists, _knn, flann::SearchParams());

	// Interpolate the value using the indices
	double cum_dist = 0.0;
	for(int i=0; i<indices[0].size(); ++i)
	{
		int indice = indices[0][i];
		vec y = _data->get(indice);

		for(int j=0; j<dimY(); ++j)
		{
			res[j] += dists[0][i] * y[dimX() + j];
		}
		cum_dist += dists[0][i];
	}
	if(cum_dist > 0.0)
	{
		for(int j=0; j<dimY(); ++j)
		{
			res[j] /= cum_dist;
		}
	}
#else

	Point pt_x(dD);
	for(int j=0; j<dimX(); ++j) { pt_x[j] = x[j]; }
	S_handle simplex = D->locate(pt_x);

	if(simplex == Delaunay_d::Simplex_handle())
	{
		return res;
	}

	// Interpolate the value of the vertices of the Simplex to estimate
	// the value of x
	double cum_dist = 0.0;
	for(int j=0; j<=dD; ++j)
	{
		Point y = D->point_of_simplex(simplex, j);

		// Compute the distance between y and x
		double dist = 0.0;
		for(int i=0; i<dimX(); ++i) { dist += pow(y.homogeneous(i)-x[i], 2); }
		dist = sqrt(dist);

		// Interpolate the data
		cum_dist += dist;
		for(int i=0; i<dimY(); ++i)
		{
			res[i] += dist * y.homogeneous(dimX() + i);
		}

	}

	if(cum_dist > 0.0)
	{
		for(int j=0; j<dimY(); ++j)
		{
			res[j] /= cum_dist;
		}
	}
#endif

   return res;
}

// Get data size, e.g. the number of samples to fit
int data_interpolant::size() const 
{
	assert(_data);
	return _data->size();
}

ALTA_DLL_EXPORT data* provide_data(const arguments&)
{
    return new data_interpolant();
}


