/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015, 2016 Inria

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

//using namespace alta;

alta::data* alta::load_data_from_text(std::istream& input,
                                      const alta::arguments& header)
{
  // FIXME: Eventually reinstate support for extra arguments when loading a
  // file.
  static alta::arguments args;

  vec min, max ;
  vec ymin, ymax;

  if(! header.is_defined("DIM")) {
    std::cerr << "<<ERROR>> Undefined dimensions ! ";
    std::cerr << "Please add DIM [int] [int] into the file header." << std::endl;
    throw;
  }

  params::input  in_param  = params::parse_input(header.get_string("PARAM_IN", "UNKNOWN_INPUT"));
  params::output out_param = params::parse_output(header.get_string("PARAM_OUT", "UNKNOWN_OUTPUT"));

  std::pair<int, int> dim = header.get_pair<int>("DIM");

  min = args.get_vec("min", dim.first, -std::numeric_limits<float>::max()) ;
  max = args.get_vec("max", dim.first,  std::numeric_limits<float>::max()) ;
#ifdef DEBUG
  std::cout << "<<DEBUG>> data will remove outside of " << min << " -> " << max << " x-interval" << std::endl;
#endif

  ymin = args.get_vec("ymin", dim.second, -std::numeric_limits<float>::max()) ;
  ymax = args.get_vec("ymax", dim.second,  std::numeric_limits<float>::max()) ;
#ifdef DEBUG
  std::cout << "<<DEBUG>> data will remove outside of " << ymin << " -> " << ymax << " y-interval" << std::endl;
#endif

  int vs_value = header.get_int("VS");

  std::vector<vec> content;

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
      vec v = vec::Zero(dim.first + 3*dim.second) ;
      for(int i=0; i<dim.first+dim.second; ++i)
      {
        linestream >> v[i] ;
      }

      // If data is not in the interval of fit
      bool is_in = true ;
      for(int i=0; i<dim.first; ++i)
      {
        if(v[i] < min[i] || v[i] > max[i])
        {
          is_in = false ;
        }
      }
      for(int i=0; i<dim.second; ++i)
      {
        if(v[dim.first+i] < ymin[i] || v[dim.first+i] > ymax[i])
        {
          is_in = false ;
        }
      }
      if(!is_in)
      {
        continue ;
      }

//      /*
      // Correction of the data by 1/cosine(theta_L)
      double factor = 1.0;
      if(args.is_defined("data-correct-cosine"))
      {
        double cart[6];
        params::convert(&v[0], in_param, params::CARTESIAN, cart);
        if(cart[5] > 0.0 && cart[2] > 0.0)
        {
          factor = 1.0/cart[5]*cart[2];
          for(int i=0; i<dim.second; ++i)
          {
            v[i + dim.first] /= factor;
          }
        }
        else
        {
          continue;
        }
      }
      // End of correction
//      */

      // Check if the data containt a vertical segment around the mean
      // value.
      for(int i=0; i<dim.second; ++i)
      {
        double min_dt = 0.0;
        double max_dt = 0.0;


        if(i == 0 && vs_value == 2)
        {
          linestream >> min_dt ;
          linestream >> max_dt ;
          min_dt = min_dt-v[dim.first+i];
          max_dt = max_dt-v[dim.first+i];
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
               v[dim.first +   dim.second+i] = v[dim.first + i] * (1.0 + min_dt) ;
          v[dim.first + 2*dim.second+i] = v[dim.first + i] * (1.0 + max_dt) ;
        }
        else if(args.is_defined("dt-max"))
        {
               v[dim.first +   dim.second+i] = v[dim.first + i] + std::max(v[dim.first + i] * min_dt, min_dt);
          v[dim.first + 2*dim.second+i] = v[dim.first + i] + std::max(v[dim.first + i] * max_dt, max_dt);
        }
        else
        {
          v[dim.first +   dim.second+i] = v[dim.first + i] + min_dt ;
          v[dim.first + 2*dim.second+i] = v[dim.first + i] + max_dt ;
        }

        // You can enforce the vertical segment to stay in the positive
        // region using the --data-positive command line argument. Note
        // that the data point is also clamped to zero if negative.
        if(args.is_defined("dt-positive"))
        {
          v[dim.first +          i] = std::max(v[dim.first +          i], 0.0);
          v[dim.first +   dim.second+i] = std::max(v[dim.first +   dim.second+i], 0.0);
          v[dim.first + 2*dim.second+i] = std::max(v[dim.first + 2*dim.second+i], 0.0);
        }

#ifdef DEBUG
                std::cout << "<<DEBUG>> vs = [" << v[dim.first +   dim.second+i] << ", " << v[dim.first + 2*dim.second+i] << "]" << std::endl;
#endif
      }

      content.push_back(std::move(v));
    }
  }

  std::cout << "<<INFO>> loaded input stream" << std::endl ;
  std::cout << "<<INFO>> loading data input of R^"
            << dim.first
            << " -> R^" << dim.second << std::endl ;
  std::cout << "<<INFO>> " << content.size() << " elements to fit" << std::endl ;

  parameters param(dim.first, dim.second, in_param, out_param);
  data* result = new vertical_segment(param, std::move(content));
  if(args.is_defined("data-correct-cosine"))
      result->save("/tmp/data-corrected.dat");

  return result;
}

void alta::save_data_as_text(std::ostream& out, const alta::data &data)
{
        using namespace alta;

    out << "#DIM " << data.parametrization().dimX() << " " << data.parametrization().dimY() << std::endl;
    out << "#PARAM_IN  "
        << params::get_name(data.parametrization().input_parametrization())
        << std::endl;
    out << "#PARAM_OUT "
        << params::get_name(data.parametrization().output_parametrization())
        << std::endl;

    for(int i=0; i < data.size(); ++i)
    {
        vec x = data.get(i);
        for(int j=0; j< data.parametrization().dimX() + data.parametrization().dimY(); ++j)
        {
                    out << std::setprecision(std::numeric_limits<double>::digits10)
                        << x[j] << "\t";
        }
        out << std::endl;
    }
}

void alta::save_data_as_binary(std::ostream &out, const alta::data& data)
{
        using namespace alta;

    out << "#DIM " << data.parametrization().dimX() << " " << data.parametrization().dimY() << std::endl;
    out << "#PARAM_IN  "
        << params::get_name(data.parametrization().input_parametrization())
        << std::endl;
    out << "#PARAM_OUT "
        << params::get_name(data.parametrization().output_parametrization())
        << std::endl;
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

        assert(sample.size() == data.parametrization().dimX() + data.parametrization().dimY());
        out.write((const char *)numbers, sample.size() * sizeof(*numbers));
    }

    out << std::endl << "#END_STREAM" << std::endl;
}

alta::data* alta::load_data_from_binary(std::istream& in, const alta::arguments& header)
{
    using namespace alta;

    // FIXME: For now we make a number of assumptions.
    assert(header["FORMAT"] == "binary");
    assert(header.get_int("VERSION") == 0);
    assert(header["PRECISION"] == "ieee754-double");
#if __BYTE_ORDER == __LITTLE_ENDIAN
    assert(header["ENDIAN"] == "little");
#else
    assert(header["ENDIAN"] == "big");
#endif

    std::pair<int, int> dim = header.get_pair<int>("DIM");
    assert(dim.first > 0 && dim.second > 0);

    in.exceptions(std::ios_base::failbit);

    int sample_count = header.get_int("SAMPLE_COUNT");
      if(sample_count <= 0) {
         std::cerr << "<<ERROR>> Uncorrect or not samples count in the header, please check \'SAMPLE_COUNT\'" << std::endl;
      }

      // TODO: Arrange to use mmap and make it zero-alloc and zero-copy.
      std::vector<vec> content(sample_count);
      for (int i = 0; i < sample_count; i++)
      {
         vec row = vec::Zero(dim.first + dim.second);
         std::streamsize expected = row.size() * sizeof(double);

         for (std::streamsize total = 0;
               total < expected && !in.eof();
               total += in.gcount())
         {
            char* ptr = (char*)row.data();
            in.read(ptr + total, expected - total);
         }

         content[i] = row;
      }

    parameters param(dim.first, dim.second,
                     params::parse_input(header["PARAM_IN"]),
                     params::parse_output(header["PARAM_OUT"]));

    return new alta::vertical_segment(param, std::move(content));
}
