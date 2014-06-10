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
 *		<li><b>\-\-in-data <i>filename</i></b> specify the data plugin
 *		used to load the input file. If this option is not specified, the
 *		loading plugin will be a \ref vertical_segment plugin. \note If
 *		the input plugin is not interpolating, like \ref vertical_segment,
 *		you can only use the reparametrization tool.</li>
 *		<li><b>\-\-output <i>filename</i></b> the resulting data file.
 *		</li>
 *		<li><b>\-\-out-data <i>filename</i></b></li> specify the data plugin
 *		used to export the data. This parameter is optional. If not defined,
 *		the output format will be ALTA's \ref format.
 *		</li>
 *  </ul>
 */
#include <core/args.h>
#include <core/data.h>
#include <core/params.h>
#include <core/function.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>
#include <core/vertical_segment.h>

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
		std::cout << "Usage: data2data [options] --input data.file --output data.file" << std::endl ;
		std::cout << "Convert a data object to another data object."<< std::endl ;
		std::cout << std::endl;
		std::cout << "Mandatory arguments:" << std::endl;
		std::cout << "  --input    [filename]" << std::endl;
		std::cout << "  --output   [filename]" << std::endl;
		std::cout << std::endl;
		std::cout << "Optional arguments:" << std::endl;
		std::cout << "  --out-data [filename]  Name of the plugin used to save the outputed data file" << std::endl ;
		std::cout << "                         If none is provided, the exporter will export in ALTA" << std::endl ;
		std::cout << "                         by default." << std::endl ;
		std::cout << "  --in-data  [filename]  Name of the plugin used to load the input data file" << std::endl ;
		std::cout << "                         If none is provided, the exporter will export in ALTA" << std::endl ;
		std::cout << "                         by default." << std::endl ;
		std::cout << "  --param    [NAME]      Name of the parametrization used to output data when" << std::endl;
		std::cout << "                         No output data plugin is specified. Please see " << std::endl;
		std::cout << "                         --help-params for the list of available " << std::endl ;
		std::cout << "                         parametrizations." << std::endl ;
		std::cout << "  --data-correct-cosine  Divide the value of the data points by the product of" << std::endl;
		std::cout << "                         the light and view vector dot product with the normal." << std::endl ;
		std::cout << std::endl;
		std::cout << "Helps:" << std::endl;
		std::cout << "  --help                 Display this help." << std::endl;
		std::cout << "  --help-params          Display the available input parametrizations." << std::endl;
		return 0 ;
	}
	if(args.is_defined("help-params")) {
		params::print_input_params();
		return 0;
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
	ptr<data> d_in = plugins_manager::get_data(args["in-data"]) ;
	d_in->load(args["input"], args);

	ptr<data> d_out = plugins_manager::get_data(args["out-data"]) ;
	
	if(!d_in && !d_out)
	{
		std::cerr << "<<ERROR>> cannot import or export data" << std::endl ;
		return 1;
	}

	if(dynamic_pointer_cast<vertical_segment>(d_out))
	{
		params::input param = params::parse_input(args["param"]);
        if(param == params::UNKNOWN_INPUT && d_in->input_parametrization() != params::UNKNOWN_INPUT)
		{
            std::cout << "<<DEBUG>> using the input file input param for the output file." << std::endl;
            param = d_in->input_parametrization();
		}
        else if(param == params::UNKNOWN_INPUT)
        {
            std::cerr << "<<ERROR>> not parametrization defined for input and output files." << std::endl;
            return -1;
        }

		d_out->setParametrization(param);
		d_out->setDimX(params::dimension(param));
		d_out->setDimY(d_in->dimY());

		std::cout << "<<INFO>> output DIM = " << d_out->dimX() << ", " << d_out->dimY() << std::endl;

		vec temp(d_out->dimX() + d_out->dimY());
		for(int i=0; i<d_in->size(); ++i)
		{
			// Copy the input vector
			vec x = d_in->get(i);
			params::convert(&x[0], d_in->parametrization(), d_out->parametrization(), &temp[0]);
			
			for(int j=0; j<d_in->dimY(); ++j)
			{
				temp[d_out->dimX() + j] = x[d_in->dimX() + j];
			}

			d_out->set(temp);

		}	
	}
	else
	{
		if(d_out->dimY() != d_in->dimY())
		{
			std::cerr << "<<ERROR>> data types have incompatible output dimensions" << std::endl;
		}

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
	}	
	d_out->save(args["output"]);
	return 0 ;
}
