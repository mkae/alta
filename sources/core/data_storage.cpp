/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "data.h"
#include "data_storage.h"
#include "vertical_segment.h"

#include <iostream>
#include <limits>
#include <iomanip>
#include <cassert>

#ifdef __GLIBC__
# include <endian.h>
#endif

void vertical_segment::load_data_from_text(std::istream& input,
																					 const arguments& header,
																					 vertical_segment& result,
																					 const arguments& args)
{
	vec min, max ;
	vec ymin, ymax;

	result._nX = 0 ; result._nY = 0 ;

	if(! header.is_defined("DIM")) {
		std::cerr << "<<ERROR>> Undefined dimensions ! ";
		std::cerr << "Please add DIM [int] [int] into the file header." << std::endl;
		throw;
	}

	params::input  in_param  = params::parse_input(header.get_string("PARAM_IN", "UNKNOWN_INPUT"));
	params::output out_param = params::parse_output(header.get_string("PARAM_OUT", "UNKNOWN_OUTPUT"));
	result.setParametrizations(in_param, out_param);

	std::pair<int, int> dim = header.get_pair<int>("DIM");
	result.setDimX(dim.first);
	result.setDimY(dim.second);

	min = args.get_vec("min", result._nX, -std::numeric_limits<float>::max()) ;
	max = args.get_vec("max", result._nX,  std::numeric_limits<float>::max()) ;
#ifdef DEBUG
	std::cout << "<<DEBUG>> data will remove outside of " << min << " -> " << max << " x-interval" << std::endl;
#endif

	ymin = args.get_vec("ymin", result._nY, -std::numeric_limits<float>::max()) ;
	ymax = args.get_vec("ymax", result._nY,  std::numeric_limits<float>::max()) ;
#ifdef DEBUG
	std::cout << "<<DEBUG>> data will remove outside of " << ymin << " -> " << ymax << " y-interval" << std::endl;
#endif

	for(int k=0; k<result.dimX(); ++k)
	{
			result._min[k] =  std::numeric_limits<double>::max() ;
			result._max[k] = -std::numeric_limits<double>::max() ;
	}

	int vs_value = header.get_int("VS");

	// Now read the body.
	while(input.good())
	{
		std::string line ;
		std::getline(input, line) ;
		std::stringstream linestream(line) ;

		// Discard comments and empty lines.
		if(line.empty() || linestream.peek() == '#')
		{
			continue ;
		}
		else
		{
			// Read the data point x and y coordinates
			vec v = vec::Zero(result.dimX() + 3*result.dimY()) ;
			for(int i=0; i<result.dimX()+result.dimY(); ++i) 
			{
				linestream >> v[i] ;
			}

			// If data is not in the interval of fit
			bool is_in = true ;
			for(int i=0; i<result.dimX(); ++i)
			{
				if(v[i] < min[i] || v[i] > max[i])
				{
					is_in = false ;
				}
			}
			for(int i=0; i<result.dimY(); ++i)
			{
				if(v[result.dimX()+i] < ymin[i] || v[result.dimX()+i] > ymax[i])
				{
					is_in = false ;
				}
			}
			if(!is_in)
			{
				continue ;
			}

//			/*
			// Correction of the data by 1/cosine(theta_L)
			double factor = 1.0;
			if(args.is_defined("data-correct-cosine"))
			{
				double cart[6];
				params::convert(&v[0], result.input_parametrization(), params::CARTESIAN, cart);
				if(cart[5] > 0.0 && cart[2] > 0.0)
				{
					factor = 1.0/cart[5]*cart[2];
					for(int i=0; i<result.dimY(); ++i) 
					{
						v[i + result.dimX()] /= factor;
					}
				}
				else
				{
					continue;
				}
			}
			// End of correction
//			*/

			// Check if the data containt a vertical segment around the mean
			// value.
			for(int i=0; i<result.dimY(); ++i)
			{
				double min_dt = 0.0;
				double max_dt = 0.0;


				if(i == 0 && vs_value == 2)
				{
					linestream >> min_dt ;
					linestream >> max_dt ;
					min_dt = min_dt-v[result.dimX()+i];
					max_dt = max_dt-v[result.dimX()+i];
				}
				else if(i == 0 && vs_value == 1)
				{
					double dt ;
					linestream >> dt ;
					min_dt = -dt;
					max_dt =  dt;
				}
				else
				{
					double dt = args.get_float("dt", 0.1f);
					min_dt = -dt;
					max_dt =  dt;
				}

				if(args.is_defined("dt-relative"))
				{
               v[result.dimX() +   result.dimY()+i] = v[result.dimX() + i] * (1.0 + min_dt) ;
					v[result.dimX() + 2*result.dimY()+i] = v[result.dimX() + i] * (1.0 + max_dt) ;
				}
				else if(args.is_defined("dt-max"))
				{
               v[result.dimX() +   result.dimY()+i] = v[result.dimX() + i] + std::max(v[result.dimX() + i] * min_dt, min_dt);
					v[result.dimX() + 2*result.dimY()+i] = v[result.dimX() + i] + std::max(v[result.dimX() + i] * max_dt, max_dt);
				}
				else
				{
					v[result.dimX() +   result.dimY()+i] = v[result.dimX() + i] + min_dt ;
					v[result.dimX() + 2*result.dimY()+i] = v[result.dimX() + i] + max_dt ;
				}

				// You can enforce the vertical segment to stay in the positive
				// region using the --data-positive command line argument. Note
				// that the data point is also clamped to zero if negative.
				if(args.is_defined("dt-positive"))
				{
					v[result.dimX() +          i] = std::max(v[result.dimX() +          i], 0.0);
					v[result.dimX() +   result.dimY()+i] = std::max(v[result.dimX() +   result.dimY()+i], 0.0);
					v[result.dimX() + 2*result.dimY()+i] = std::max(v[result.dimX() + 2*result.dimY()+i], 0.0);
				}

#ifdef DEBUG
                std::cout << "<<DEBUG>> vs = [" << v[result.dimX() +   result.dimY()+i] << ", " << v[result.dimX() + 2*result.dimY()+i] << "]" << std::endl;
#endif
			}

			result._data.push_back(v) ;

			// Update min and max
			for(int k=0; k<result.dimX(); ++k)
			{
				result._min[k] = std::min(result._min[k], v[k]) ;
				result._max[k] = std::max(result._max[k], v[k]) ;
			}
		}
	}

	if(args.is_defined("data-correct-cosine"))
		result.save("/tmp/data-corrected.dat");

	std::cout << "<<INFO>> loaded input stream" << std::endl ;
	std::cout << "<<INFO>> data inside " << result._min << " ... " << result._max << std::endl ;
	std::cout << "<<INFO>> loading data input of R^" << result.dimX() << " -> R^" << result.dimY() << std::endl ;
	std::cout << "<<INFO>> " << result._data.size() << " elements to fit" << std::endl ;
}

void save_data_as_text(std::ostream& out, const data &data)
{
		out << "#DIM " << data.dimX() << " " << data.dimY() << std::endl;
		out << "#PARAM_IN  " << params::get_name(data.input_parametrization())  << std::endl;
		out << "#PARAM_OUT " << params::get_name(data.output_parametrization()) << std::endl;
		for(int i=0; i < data.size(); ++i)
		{
				vec x = data.get(i);
				for(int j=0; j< data.dimX() + data.dimY(); ++j)
				{
            out << std::setprecision(std::numeric_limits<double>::digits10)
                << x[j] << "\t";
				}
				out << std::endl;
		}
}

void save_data_as_binary(std::ostream &out, const data& data)
{
		out << "#DIM " << data.dimX() << " " << data.dimY() << std::endl;
		out << "#PARAM_IN  " << params::get_name(data.input_parametrization())  << std::endl;
		out << "#PARAM_OUT " << params::get_name(data.output_parametrization()) << std::endl;
		out << "#FORMAT binary" << std::endl;
		out << "#VERSION 0" << std::endl;
		out << "#PRECISION ieee754-double" << std::endl;
		out << "#SAMPLE_COUNT " << data.size() << std::endl;

		// FIXME: Note: on non-glibc systems, both macros may be undefined, so
		// the conditional is equivalent to "#if 0 == 0", which is usually what
		// we want.
#if __BYTE_ORDER == __LITTLE_ENDIAN
		out << "#ENDIAN little" << std::endl;
#else
		out << "#ENDIAN big" << std::endl;
#endif

		out << "#BEGIN_STREAM" << std::endl;

		for(int i=0; i < data.size(); ++i)
		{
				vec sample = data.get(i);
				const double *numbers = sample.data();

				assert(sample.size() == data.dimX() + data.dimY());
				out.write((const char *)numbers, sample.size() * sizeof(*numbers));
		}

		out << std::endl << "#END_STREAM" << std::endl;
}

void load_data_from_binary(std::istream& in, const arguments& header, data& data)
{
		// FIXME: For now we make a number of assumptions.
		assert(header["FORMAT"] == "binary");
		assert(header.get_int("VERSION") == 0);
		assert(header["PRECISION"] == "ieee754-double");
#if __BYTE_ORDER == __LITTLE_ENDIAN
		assert(header["ENDIAN"] == "little");
#else
		assert(header["ENDIAN"] == "big");
#endif

		data.setParametrizations(params::parse_input(header["PARAM_IN"]),
										 params::parse_output(header["PARAM_OUT"]));

		std::pair<int, int> dim = header.get_pair<int>("DIM");
		data.setDimX(dim.first);
		data.setDimY(dim.second);
      assert(data.dimX() > 0 && data.dimY() > 0);

		in.exceptions(std::ios_base::failbit);

		int sample_count = header.get_int("SAMPLE_COUNT");
      if(sample_count <= 0) {
         std::cerr << "<<ERROR>> Uncorrect or not samples count in the header, please check \'SAMPLE_COUNT\'" << std::endl;
      }

		// Initialize the mininum and maximum values.
		vec min(dim.first), max(dim.first);

      for(int k = 0; k < dim.first; k++)
      {
         min[k] =  std::numeric_limits<double>::max();
         max[k] = -std::numeric_limits<double>::max();
      }

		// TODO: Arrage to use mmap and make it zero-alloc and zero-copy.
      for (int i = 0; i < sample_count; i++)
      {
         vec row = vec::Zero(data.dimX() + data.dimY());
         std::streamsize expected = row.size() * sizeof(double);

         for (std::streamsize total = 0;
               total < expected && !in.eof();
               total += in.gcount())
         {
            char* ptr = (char*)row.data();
            in.read(ptr + total, expected - total);
         }

         data.set(row);

         // Update min and max.
         for(int k = 0; k < dim.first; k++)
         {
            min[k] = std::min(min[k], row[k]);
            max[k] = std::max(max[k], row[k]);
         }
      }

		// TODO: Factorize this in data::set.
		data.setMin(min);
		data.setMax(max);
}
