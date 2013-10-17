#include "data.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

data_brdf_slice::data_brdf_slice()
{
	// Allocate data
	_data = new fipImage();
	//_data->setSize(FIT_FLOAT, 512, 512, 16);

    // Set the input and output parametrization
    _in_param  = params::RUSIN_TH_TD;
    _out_param = params::RGB_COLOR;
    _nX = 2;
    _nY = 3;
}

data_brdf_slice::~data_brdf_slice()
{
	delete _data;
}

// Load data from a file
void data_brdf_slice::load(const std::string& filename) 
{
	_data->load(filename.c_str());
	_data->convertTo32Bits();
	width  = _data->getWidth();
	height = _data->getHeight();
}
void data_brdf_slice::load(const std::string& filename, const arguments&)
{
	load(filename);
}

void data_brdf_slice::save(const std::string& filename) const 
{
	if(!_data->save(filename.c_str()))
	{
		std::cerr << "<<ERROR>> unable to save image file" << std::endl;
	}
}

// Acces to data
vec data_brdf_slice::get(int id) const 
{
	vec res(3) ;
	int i = id % width;
	int j = id / width;

	RGBQUAD pixel;
	_data->getPixelColor(i, j, &pixel);

	res[0] = pixel.rgbRed / 255.0;
	res[1] = pixel.rgbGreen / 255.0;
	res[2] = pixel.rgbBlue / 255.0;

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

	int i = floor(x[0] * width / (0.5*M_PI));
	int j = floor(x[1] * height / (0.5*M_PI));

	RGBQUAD* pixel;
	_data->getPixelColor(i, j, pixel);
	pixel->rgbRed   = x[2];
	pixel->rgbGreen = x[3];
	pixel->rgbBlue  = x[4];
}

vec data_brdf_slice::value(vec, vec) const
{
	vec res(3);
	return res;
}
vec data_brdf_slice::value(vec x) const
{
	int i = floor(x[0] * width / (0.5*M_PI));
	int j = floor(x[1] * height / (0.5*M_PI));

	if(i < 0 || i >= width)  { std::cerr << "<<ERROR>> out of bounds: " << x << std::endl; }
	if(j < 0 || j >= height) { std::cerr << "<<ERROR>> out of bounds: " << x << std::endl; }

   return get(i + j*width);
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


