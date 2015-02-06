/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include <core/args.h>
#include <core/data.h>
#include <core/function.h>
#include <core/vertical_segment.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>

int main(int argc, char** argv)
{
	arguments args(argc, argv);

	if(args.is_defined("help")) 
	{
		std::cout << "Usage: brdf2gnuplot [options] --input brdf.file --output gnuplot.file" << std::endl ;
		std::cout << "   ->  input, and output parameters are mandatory parameters" << std::endl ;
		std::cout << std::endl;
		std::cout << "Options:" << std::endl;
		std::cout << "  --data [filename]      data plugin module used to define the abcissas for the." << std::endl;
		std::cout << "                         export, or to load the --in-data file.." << std::endl;
		std::cout << "  --data-file [filename] produce a data by using the abcissas of the [filename]" << std::endl;
		std::cout << "                         data file." << std::endl;
		std::cout << "  --polar-plot           produce a polar data by sampling regularly the elevation" << std::endl;
		std::cout << "                         angle. Use parameter --inc-angle to define the incoming" << std::endl;
		std::cout << "                         angle. Use parameter --samples to define the number of" << std::endl;
		std::cout << "                         samples on this domain." << std::endl;
		std::cout << "  --inc-angle [float]    set the incoming light elevation in radian for the polar" << std::endl;
		std::cout << "                         plot export." << std::endl;
		std::cout << "  --samples [int]        set the points used for the polar plot export." << std::endl;
		std::cout << "  --cos-plot             export the BRDF*cosine instead of the BRDF alone." << std::endl;
		return 0;
	}
	
	if(! args.is_defined("input")) 
	{
		std::cerr << "<<ERROR>> the input filename is not defined" << std::endl;
		return 1;
	}
	if(! args.is_defined("output")) 
	{
		std::cerr << "<<ERROR>> the output filename is not defined" << std::endl;
		return 1;
	}


	// Load a function file
	function* f = plugins_manager::get_function(args["input"]) ;
	if(f == NULL)
	{
		return 1;
	}

	// Create output file
	std::ofstream file(args["output"].c_str(), std::ios_base::trunc);
	file.precision(10);

	// Should I export a BRDF or BRDF*cos ?
	// Cannot export a cosine term if no parametrization is defined for the
	// input function
	const bool cos_plot = args.is_defined("cos-plot") && f->input_parametrization() != params::UNKNOWN_INPUT;


	// Load a data file
	if(args.is_defined("data") || args.is_defined("data-file"))
	{
		ptr<data> d = plugins_manager::get_data(args["data"]);

		// Load data file if the plugin manager created a plugin object.
		if(d)
		{
			d->load(args["data-file"]);
		}
		else
		{
			std::cerr << "<<ERROR>> unable to load the data plugin" << std::endl;
			return 1;
		}

		// Print the distance to the data to check if it correspond to the value
		// computed prior.
		double L2   = f->L2_distance(d);
		double Linf = f->Linf_distance(d);
		std::cout << "<<INFO>> L2   distance to data = " << L2   << std::endl;
		std::cout << "<<INFO>> Linf distance to data = " << Linf << std::endl;

		// Check the kind of plot to do
		bool plot_error = false ;
		bool linear_plot = false;
		if(args.is_defined("error"))
		{
			std::cout << "<<INFO>> Exporting an error plot" << std::endl;
			plot_error = true ;
		}
		else if(args.is_defined("linear_error"))
		{
			std::cout << "<<INFO>> Exporting linear error plot" << std::endl;
			linear_plot = true;
		}

		for(int i=0; i<d->size(); ++i)
		{
			vec v = d->get(i) ;
			vec x(f->dimX());

			// Convert the data to the function's input space.
			if(f->input_parametrization() == params::UNKNOWN_INPUT)
			{
				memcpy(&x[0], &v[0], f->dimX()*sizeof(double));
			}
			else
			{
				params::convert(&v[0], d->input_parametrization(), f->input_parametrization(), &x[0]);
			}

			// Evaluate the function. I can add the cosine term to the BRDF
			// value.
			double costerm = 1.0;
			if(cos_plot) 
			{
				double cart[6];
				params::convert(&x[0], f->input_parametrization(), params::CARTESIAN, cart);
				costerm = cart[5]*cart[2];
			}
			vec y2 = costerm * f->value(x) ;


			if(!linear_plot)
			{
				for(int u=0; u<d->dimX(); ++u)
					file << v[u] << "\t" ;
			}
			else
			{
				file << i << "\t" ;
			}

			for(int u=0; u<d->dimY(); ++u)
			{
				if(plot_error)
				{
					file << (v[d->dimX() + u] - y2[u]) << "\t" ;
				}
				else if(linear_plot)
				{
					file << (v[d->dimX() + u] - y2[u])/v[d->dimX()+u] << "\t" ;
				}
				else
				{
					file << y2[u] << "\t" ;
				}
			}

			file << std::endl ;
		}
	}
	else if(args.is_defined("polar-plot"))
	{
		vec spherical(4);
		spherical[0] = args.get_float("inc-angle", 0.0);
		spherical[1] = M_PI;

		const int N = args.get_int("samples", 100) / 2;

		// Plot retro direction
		for(int i=0; i<N; ++i)
		{
			spherical[2] = 0.5 * M_PI * double(i) / double(N);
			spherical[3] = M_PI;

			vec x(f->dimX());
			params::convert(&spherical[0], params::SPHERICAL_TL_PL_TV_PV, f->input_parametrization(), &x[0]);

			vec y = f->value(x);
			file << -spherical[2] << "\t";
			for(int k=0; k<f->dimY(); ++k) { file << y[k] << "\t"; }
			file << std::endl;
		}

		// Plot forward direction
		for(int i=0; i<N; ++i)
		{
			spherical[2] = 0.5 * M_PI * double(i) / double(N);
			spherical[3] = 0.0;

			vec x(f->dimX());
			params::convert(&spherical[0], params::SPHERICAL_TL_PL_TV_PV, f->input_parametrization(), &x[0]);

			vec y = f->value(x);
			file << spherical[2] << "\t";
			for(int k=0; k<f->dimY(); ++k) { file << y[k] << "\t"; }
			file << std::endl;
		}

	}

	file.close();
	return 0 ;
}
