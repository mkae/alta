#include "data.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "EXR_IO.h"

data_brdf_slice::data_brdf_slice()
{
	// Allocate data
	width = 512; height = 512;
	_data = new double[3*width*height];

    // Set the input and output parametrization
    _in_param  = params::RUSIN_TH_TD;
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
	if(!t_EXR_IO<double>::SaveEXR(filename.c_str(), width, height, _data))
	{
		std::cerr << "<<ERROR>> unable to save image file" << std::endl;
	}
}

// Acces to data
vec data_brdf_slice::get(int id) const 
{
	vec res(5) ;
	const int i = id % width;
	const int j = id / width;

	res[0] = 0.5*M_PI*(i+0.5) / double(width);
	res[1] = 0.5*M_PI*(j+0.5) / double(height);

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
void data_brdf_slice::set(vec x)
{
	assert(x.size() == 5);
	assert(x[0] <= 0.5*M_PI && x[0] >= 0.0);
	assert(x[1] <= 0.5*M_PI && x[1] >= 0.0);

	const int i  = floor(x[0] * width / (0.5*M_PI));
	const int j  = floor(x[1] * height / (0.5*M_PI));

	const int id = i + j*width;

	_data[3*id + 0] = x[2];
	_data[3*id + 1] = x[3];
	_data[3*id + 2] = x[4];
}

vec data_brdf_slice::value(vec, vec) const
{
	vec res(3);
	return res;
}
vec data_brdf_slice::value(vec x) const
{
	assert(x[0] <= 0.5*M_PI && x[0] >= 0.0);
	assert(x[1] <= 0.5*M_PI && x[1] >= 0.0);

	const int i  = floor(x[0] * width / (0.5*M_PI));
	const int j  = floor(x[1] * height / (0.5*M_PI));
	const int id = i + j*width;

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
	return width*height;
}

// Get min and max input space values
vec data_brdf_slice::min() const 
{
	vec res(2);
	res[0] = 0.0 ;
	res[1] = 0.0 ;
	return res ;
}
vec data_brdf_slice::max() const
{
	vec res(2);
	res[0] = M_PI / 2 ;
	res[1] = M_PI / 2 ;
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


