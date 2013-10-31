#include "data.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include <core/vertical_segment.h>

data_interpolant::data_interpolant()
{
	_kdtree = new flann::Index< flann::L2<double> >(flann::KDTreeIndexParams(4));
	_data = new vertical_segment();

	_knn = 10;
}

data_interpolant::~data_interpolant()
{
	delete _data;
	delete _kdtree;
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

	std::cout << "   " << _kdtree->veclen() << std::endl;

	// Update the KDtreee by inserting all points
	for(int i=0; i<_data->size(); ++i)
	{
		vec x = _data->get(i);
		flann::Matrix<double> pts(&x[0], dimX(), 1);

		_kdtree->addPoints(pts);
	}
	_kdtree->buildIndex();
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
vec data_interpolant::value(vec x) const
{
	vec res = vec::Zero(dimY());

	// Query point
	flann::Matrix<double> pts(&x[0], dimX(), 1);
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
	for(int j=0; j<dimY(); ++j)
	{
		res[j] /= cum_dist;
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


