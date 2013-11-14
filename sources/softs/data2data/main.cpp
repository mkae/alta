/*! \package data2data
 *  \ingroup commands
 *  \brief
 *  This command transform a \ref data object into another data object.
 *  \details
 *  This command is useful to change the parametrization of a data file,
 *  or to perform interpolation of sample values (i.e. to fill gaps).
 *
 *  <h3>Usage</h3>
 *  \verbatim	
       data2data --input data.file --output data.file --out-data exporter.so --in-data importer.so
    \endverbatim
 *
 *  <h3>Parameters</h3>
 *  <ul>
 *		<li><b>\-\-input <i>filename</i></b> data file to be loaded. The
 *		loading plugin is defined using the option <b>\-\-in-data <i>
 *		filename</i></b>.
 *		</li>
 *		<li><b>\-\-in-data <i>filename</i></b></li> specify the data plugin
 *		used to load the input file. If this option is not specified, the
 *		loading plugin will be a \ref vertical_segment plugin. \note If
 *		the input plugin is not interpolating, like \ref vertical_segment,
 *		you can only use the reparametrization tool.
 *		<li><b>\-\-input <i>filename</i></b> the resulting data file.
 *		</li>
 *  </ul>
 */
#include <core/args.h>
#include <core/data.h>
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

int main(int argc, char** argv)
{
	arguments args(argc, argv) ;

	if(args.is_defined("help")) {
		std::cout << "<<HELP>> data2data --input data.file --output data.file --out-data exporter.so --in-data importer.so" << std::endl ;
		std::cout << " - input, output, out-data are mandatory parameters" << std::endl ;
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
	/*
	if(! args.is_defined("out-data")) {
		std::cerr << "<<ERROR>> the data exporter is not defined" << std::endl ;
		return 1 ;
	}
	*/
	
	
	// Import data
	data* d_in = NULL ;
	d_in = plugins_manager::get_data(args["in-data"]) ;
	d_in->load(args["input"]);

	data* d_out = NULL;
	d_out = plugins_manager::get_data(args["out-data"]) ;

	if(d_out->dimY() != d_in->dimY())
	{
		std::cerr << "<<ERROR>> data types have incompatible output dimensions" << std::endl;
	}

	if(d_in != NULL && d_out != NULL)
	{
		vec temp(d_in->dimX());
		for(int i=0; i<d_out->size(); ++i)
		{
			// Copy the input vector
			vec x = d_out->get(i);
			params::convert(&x[0], d_out->parametrization(), d_in->parametrization(), &temp[0]);
			
			vec y = d_in->value(temp);

			for(int j=0; j<d_in->dimY(); ++j)
			{
				x[d_out->dimX() + j] = y[j];
			}

			d_out->set(x);
		}	

		d_out->save(args["output"]);
		return 0 ;
	}	
	else
	{
		std::cerr << "<<ERROR>> cannot import or export data" << std::endl ;
		return 1;
	}

}
