/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

/*! \package brdf2data
 *  \ingroup commands
 *  \brief
 *  This command allows to convert a \ref function object to a \ref data 
 *  object. And to save the \ref data object in a file specified by the
 *  \ref data plugin.
 *  \details
 */
#include <core/args.h>
#include <core/data.h>
#include <core/vertical_segment.h>
#include <core/params.h>
#include <core/function.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <cmath>

using namespace alta;

int main(int argc, char** argv)
{
	arguments args(argc, argv) ;

	if(args.is_defined("help")) {
		std::cout << "Usage: brdf2data --input brdf.file --output data.file [--data exporter.so --data-file data.file]" << std::endl ;
		std::cout << "Convert a function object to a data object."<< std::endl ;
		std::cout << std::endl;
		std::cout << "Mandatory arguments:" << std::endl;
		std::cout << "  --input     [filename]" << std::endl;
		std::cout << "  --output    [filename]" << std::endl;
		std::cout << std::endl;
		std::cout << "Optional arguments:" << std::endl;
		std::cout << "  --data      [filename] Name of the data plugin used to save the output" << std::endl ;
		std::cout << "                         data file. If no plugin is defined, the data file" << std::endl ;
		std::cout << "                         will be load using ALTA format but require an" << std::endl;
		std::cout << "                         ALTA file as template." << std::endl ;
		std::cout << "  --data-file [filename] ALTA data file used as a template if no data" << std::endl ;
		std::cout << "                         plugin is specified to export data." << std::endl ;
		return 0 ;
	}

	if(! args.is_defined("input")) {
		std::cerr << "<<ERROR>> the input filename is not defined" << std::endl ;
		return 1 ;
	}
	if(! args.is_defined("output")) {
		std::cerr << "<<ERROR>> the output filename is not defined" << std::endl ;
		return 1 ;
	}
	if(! args.is_defined("data") && ! args.is_defined("data-file")) {
		std::cerr << "<<ERROR>> the data exporter is not defined" << std::endl ;
		return 1 ;
	}
	
	// Get the associated data object and load the file is any
	ptr<data> d;

  try
  {
      d = plugins_manager::load_data(args["data-file"], args["data"], args);
	}
  CATCH_FILE_IO_ERROR(args["data-file"]);

    // Get the output object. In the case where it is not a VS file, we use
    // the load object.
    ptr<data> d_out = plugins_manager::get_data(args["data"], args);
    if(dynamic_pointer_cast<vertical_segment>(d) != NULL)
    {
        parameters p(d->parametrization().dimX(),
                     d->parametrization().dimY(),
                     d->parametrization().input_parametrization(),
                     d->parametrization().output_parametrization());
        d_out->setParametrization(p);
    }

	// Get the function file
	function* f = NULL;
	f = plugins_manager::load_function(args["input"]);
	if(f == NULL)
	{
		std::cerr << "<<ERROR>> cannot open the function file" << std::endl;
		return 1;
	}

	if(d && f != NULL)
	{
		// Is the output data file already allocated and has the same size
		// than the training data ?
		const bool out_filled = d->size() == d_out->size();
		const bool output_dif = args.is_defined("export-diff");

		vec temp(f->parametrization().dimX());
		for(int i=0; i<d->size(); ++i)
		{
        // Copy the input vector
        vec x = d->get(i);
        // Convert the data to the function's input space.
        if(f->parametrization().input_parametrization() == params::UNKNOWN_INPUT)
        {
            memcpy(&temp[0], &x[0], f->parametrization().dimX()*sizeof(double));
        }
        else
        {
            params::convert(&x[0],
                            d->parametrization().input_parametrization(),
                            f->parametrization().input_parametrization(),
                            &temp[0]);
        }
        vec y = f->value(temp);

        for(int j=0; j<d->parametrization().dimY(); ++j) {
            x[d->parametrization().dimX() + j] =
                (output_dif) ? x[d->parametrization().dimX() + j] - y[j] : y[j];
        }

        // If the output data is already allocated and has the same size
        // than the training data, we do simple copy of the index elements.
        if(out_filled) {
            d_out->set(i, y);
        } else {
            d_out->set(x);
        }
		}	

        d_out->save(args["output"]);
	}	
	else
	{
		std::cerr << "<<ERROR>> cannot import function or export data" << std::endl ;
	}

	return 0 ;
}
