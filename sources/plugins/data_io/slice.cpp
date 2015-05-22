/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include <core/data.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "EXR_IO.h"

/*! \ingroup datas
 *  \brief Data interface for the BRDF slice file format.
 *  Plugin name: *data_brdf_slice*
 *
 *  \details 
 *
 *  This BRDF data format implements the 2D representation of a BRDF presented
 *  by Brent Burley in BRDF Explorer. It only stores a slice of the BRDF for
 *  different 2D parametrizations. This plugin stores the BRDF data into an
 *  EXR file.
 *
 *  It is possible to select the parametrization using the --param NAME 
 *  argument when loading the BRDF. The default parametrization is
 *  \ref STARK_2D
 *
 *  \author Laurent Belcour <laurent.belcour@gmail.com>
 *
 */
class BrdfSlice : public data {
	public:

		int width, height, slice;
		double* _data;

		BrdfSlice(const arguments& args) : data()
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
			
			// Allow to load a different parametrization depending on the 
			// parameters provided.
			if(args.is_defined("param")) {
				params::input param = params::parse_input(args["param"]);
				if(params::dimension(param) == 2) {
					std::cout << "<<INFO>> Specified param \"" << args["param"] << "\"" << std::endl;
					this->setParametrization(param);
				} else {
					std::cout << "<<ERROR>> Invalid specified param \"" << args["param"] << "\"" << std::endl;
					std::cout << "<<ERROR>> Must have 2D input dimension" << std::endl;
				}
			}
		}

		~BrdfSlice()
		{
			delete[] _data;
		}

		// Load data from a file
		void load(const std::string& filename) 
		{
			delete[] _data;
			t_EXR_IO<double>::LoadEXR(filename.c_str(), width, height, _data);
		}
		void load(const std::string& filename, const arguments& args)
		{
			load(filename);

		}

		void save(const std::string& filename) const 
		{
			if(!t_EXR_IO<double>::SaveEXR(filename.c_str(), width, slice*height, _data))
			{
				std::cerr << "<<ERROR>> unable to save image file" << std::endl;
			}
		}

		// Acces to data
		vec get(int id) const 
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
		vec operator[](int i) const 
		{
			return get(i) ;
		}

		//! \todo Test this function
		void set(const vec& x)
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
		void set(int id, const vec& x)
		{
			assert(x.size() == 3);

			_data[3*id + 0] = x[0];
			_data[3*id + 1] = x[1];
			_data[3*id + 2] = x[2];
		}

		vec value(const vec& x) const
		{
			// Safeguard. We can use either asserting or returning zero in case
			// the query values are not within reach.
			if(x[0] > 1.0 || x[0] < 0.0 || x[1] > 1.0 || x[1] < 0.0 ||
				isnan(x[0]) || isnan(x[1])) {
					vec res(3);
					return res;
			}
			/*
			assert(x[0] <= 1.0 && x[0] >= 0.0);
			assert(x[1] <= 1.0 && x[1] >= 0.0);
			*/

			const int i  = floor(x[0] * width  / 1.0);
			const int j  = floor(x[1] * height / 1.0);
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
		int size() const 
		{
			return width*height*slice;
		}

		// Get min and max input space values
		vec min() const 
		{
			vec res(2);
			res[0] = 0.0 ;
			res[1] = 0.0 ;
			//res[2] = 0.0 ;
			return res ;
		}
		vec max() const
		{
			vec res(2);
			res[0] = M_PI / 2 ;
			res[1] = M_PI / 2 ;
			//res[2] = M_PI;
			return res ;
		}

		int dimX() const 
		{ 
			return 2 ; 
		}
		int dimY() const 
		{ 
			return 3;
		}
};

ALTA_DLL_EXPORT data* provide_data(const arguments& args)
{
    return new BrdfSlice(args);
}


