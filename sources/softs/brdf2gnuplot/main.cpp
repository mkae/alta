#include <core/args.h>
#include <core/data.h>
#include <core/function.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>

#include <QCoreApplication>

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>

int main(int argc, char** argv)
{
	QCoreApplication app(argc, argv, false);
	arguments args(argc, argv) ;

	plugins_manager manager(args) ;

	if(args.is_defined("help")) 
	{
		std::cout << "<<HELP>> brdf2gnuplot --input brdf.file --output gnuplot.file --func function.lib --data data.file" << std::endl ;
		std::cout << " - input, output and data are mandatory parameters" << std::endl ;
	}

	if(! args.is_defined("input")) 
	{
		std::cerr << "<<ERROR>> the input filename is not defined" << std::endl ;
		return 1 ;
	}
	if(! args.is_defined("output")) 
	{
		std::cerr << "<<ERROR>> the output filename is not defined" << std::endl ;
		return 1 ;
	}

	// Load a function file
	function* f = manager.get_function(args["input"]) ;

	// Load a data file
	data* d = NULL ;
	if(args.is_defined("data"))
	{
		std::cout << "<<INFO>> Using data \"" << args["data"] << "\"" << std::endl ;
		d = manager.get_data() ;
		d->load(args["data"]) ;
	}

    /*
    // Print the distance to the data to check if it correspond to the value
    // computed prior.
    double L2   = f->L2_distance(d);
    double Linf = f->Linf_distance(d);
    std::cout << "<<INFO>> L2   distance to data = " << L2   << std::endl;
    std::cout << "<<INFO>> Linf distance to data = " << Linf << std::endl;
    */

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

	// Create output file
	std::ofstream file(args["output"].c_str(), std::ios_base::trunc);

	if(d != NULL)
	{
		for(int i=0; i<d->size(); ++i)
		{
			vec v = d->get(i) ;

            vec x(f->dimX());
            params::convert(&v[0], d->input_parametrization(), f->input_parametrization(), &x[0]);

            vec y2 = f->value(x) ;
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
					file << (v[d->dimX() + u] - y2[u])/v[d->dimX()+u] << "\t" ;
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
	else
	{
		std::cerr << "<<ERROR>> --data is not defined" << std::endl ;
	}

	file.close();
	return 0 ;
}
