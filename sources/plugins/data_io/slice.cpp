/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015 Inria

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
 *  \class data_brdf_slice
 *  \brief Data interface for the BRDF slice file format.
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
 *  \ref params::STARK_2D "STARK_2D"
 *
 *  \author Laurent Belcour <laurent.belcour@umontreal.ca>
 *
 */
class BrdfSlice : public data {
	public:

		int width, height, slice;
		vec _max, _min;
		double _phi;
		double* _data;
		bool _reverse;

		BrdfSlice(const arguments& args) : data()
		{
			// Allocate data
			width = 512; height = 512;
			slice = 1;
			_data = new double[3*width*height*slice];
			_phi = 0.5*M_PI;
	
			// Set the input and output parametrization
			_in_param  = params::STARK_2D;
			_out_param = params::RGB_COLOR;
			_nX = 2;
			_nY = 3;

			// Allow to load a different parametrization depending on the 
			// parameters provided.
			if(args.is_defined("param")) {
				params::input param = params::parse_input(args["param"]);

				// The param is a 2D param
				if(params::dimension(param) == 2) {
					std::cout << "<<INFO>> Specified param \"" << args["param"] << "\"" << std::endl;
					this->setParametrization(param);

				// The oaram is a 3D param
				} else if(params::dimension(param) == 3) {
					std::cout << "<<INFO>> Specified param \"" << args["param"] << "\"" << std::endl;
					this->setParametrization(param);
					this->_phi = (M_PI / 180.0) * args.get_float("phi", 90);
					_nX = 3;


				} else {
					std::cout << "<<ERROR>> Invalid specified param \"" << args["param"] << "\"" << std::endl;
					std::cout << "<<ERROR>> Must have 2D input dimension" << std::endl;
				}
			}

			// Is the position of the slice componnent (third coordinate)
			// reversed? This ensure that some params can be displayed.
			_reverse = _in_param == params::ISOTROPIC_TL_TV_PROJ_DPHI ||
			           _in_param == params::SCHLICK_TL_TK_PROJ_DPHI   ||
						  _in_param == params::RETRO_TL_TVL_PROJ_DPHI;
			
			// Update the domain
			_max = max();
			_min = min();
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
			vec res(_nX+_nY) ;
			const int i = id % width;
			const int k = id / (width*height);
			const int j = (id - k*width*height) / width;

			res[0] = (i+0.5) * (_max[0]-_min[0]) / double(width)  + _min[0];
			res[1] = (j+0.5) * (_max[1]-_min[1]) / double(height) + _min[1];
			if(_nX == 3) {
				res[2] = _phi;
			}

			// Reverse the first part of the vector
			if(_reverse) {
				res.segment(0, _nX).reverseInPlace();
			}

			res[_nX+0] = _data[3*id + 0];
			res[_nX+1] = _data[3*id + 1];
			res[_nX+2] = _data[3*id + 2];

			return res ;
		}

		//! \todo Test this function
		void set(const vec& x)
		{
			// Copy vector is required
			vec _x = x;

			// Reverse the first part of the vector
			if(_reverse) {
				_x.segment(0, _nX).reverseInPlace();
			}

			assert(_x.size() == _nX+_nY);
			assert(_x[0] <= _max[0] && _x[0] >= _min[0]);
			assert(_x[1] <= _max[1] && _x[1] >= _min[1]);

			const int i  = floor((_x[0]-_min[0]) * width  / (_max[0] - _min[0]));
			const int j  = floor((_x[1]-_min[1]) * height / (_max[1] - _min[1]));
			const int k  = 0; 
			//const int k  = floor(x[2] * slice  / (M_PI));
			const int id = i + j*width + k*width*height;

			_data[3*id + 0] = _x[_nX+0];
			_data[3*id + 1] = _x[_nX+1];
			_data[3*id + 2] = _x[_nX+2];
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
			// Copy vector is required
			vec _x = x;

			// Reverse the first part of the vector
			if(_reverse) {
				_x.segment(0, _nX).reverseInPlace();
			}

			// Safeguard. We can use either asserting or returning zero in case
			// the query values are not within reach.
			if(_x[0] > _max[0] || x[0] < _min[0] || x[1] > _max[1] || x[1] < _min[1] ||
				isnan(_x[0]) || isnan(x[1])) {
					vec res(3);
					return res;
			}
			/*
			assert(_x[0] <= 1.0 && x[0] >= 0.0);
			assert(_x[1] <= 1.0 && x[1] >= 0.0);
			*/
			const int i  = floor((_x[0]-_min[0]) * width  / (_max[0] - _min[0]));
			const int j  = floor((_x[1]-_min[1]) * height / (_max[1] - _min[1]));
			const int k  = 1; 
			//const int k  = floor(_x[2] * slice  / (M_PI));
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
			vec res(_nX);
			
			// First fill the third dimension. It can be overwritten by the next
			// part when the parametrization is reversed (projected ones).
			if(_nX == 3) {
				res[2] = 0.0 ;
			}

			// Fill the first two dimension unless the parametrization is a
			// projected one then it will fill the three components.
			if(_in_param == params::ISOTROPIC_TL_TV_PROJ_DPHI ||
				_in_param == params::ISOTROPIC_TV_TL_DPHI) {	
				res[0] = -0.5*M_PI ;
				res[1] = -0.5*M_PI ;
			} else if(_in_param == params::ISOTROPIC_TL_TV_PROJ_DPHI ||
					    _in_param == params::SCHLICK_TL_TK_PROJ_DPHI   ||
						 _in_param == params::RETRO_TL_TVL_PROJ_DPHI) {
				res[0] = -0.5*M_PI ;
				res[1] = -0.5*M_PI ;
				res[2] = -0.5*M_PI ;
			} else {
				res[0] = 0.0 ;
				res[1] = 0.0 ;
			}
			return res ;
		}
		vec max() const
		{
			vec res(_nX);

			// First fill the third dimension. It can be overwritten by the next
			// part when the parametrization is reversed (projected ones).
			if(_nX == 3) {
				res[2] = 2.0*M_PI;
			}

			// Fill the first two dimension unless the parametrization is a
			// projected one then it will fill the three components.
			if(_in_param == params::RUSIN_TH_TD_PD || 
				_in_param == params::ISOTROPIC_TV_TL_DPHI) {	
				res[0] = 0.5*M_PI ;
				res[1] = 0.5*M_PI ;
			} else if(_in_param == params::ISOTROPIC_TL_TV_PROJ_DPHI ||
					    _in_param == params::SCHLICK_TL_TK_PROJ_DPHI   ||
						 _in_param == params::RETRO_TL_TVL_PROJ_DPHI) {
				res[0] = 0.5*M_PI ;
				res[1] = 0.5*M_PI ;
				res[2] = 0.5*M_PI ;
			} else {
				res[0] = 1.0 ;
				res[1] = 1.0 ;
			}
			return res ;
		}

		int dimX() const 
		{ 
			return _nX ; 
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


