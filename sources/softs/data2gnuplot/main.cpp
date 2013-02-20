#include <core/args.h>
#include <core/data.h>
#include <core/function.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>

#include <QApplication>

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <cmath>

int main(int argc, char** argv)
{
	QApplication app(argc, argv, false);
	arguments args(argc, argv) ;

	plugins_manager manager(args) ;

	if(args.is_defined("help")) {
		std::cout << "<<HELP>> data2gnuplot --input data.file --output gnuplot.file" << std::endl ;
		std::cout << " - input and output are mandatory parameters" << std::endl ;
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

	data* d = NULL ;
	d = manager.get_data(args["loader"]) ;
	d->load(args["input"]);

	// Create output file
	std::ofstream file(args["output"].c_str(), std::ios_base::trunc);

	if(d != NULL)
	{
		std::cout << "<<INFO>> will export " << d->size() << " elements" << std::endl ;
	
		vec in(3), out(3) ;
		in[0] = 0.0;
		in[1] = 0.0;
		in[2] = 1.0;
		for(int i=0; i<90; ++i)
			for(int j=0; j<90; ++j)
			{
				double phi   = i * M_PI / 89 ;
				double theta = j * M_PI / 89 * 0.5 ;
				out[0] = cos(phi)*sin(theta);
				out[1] = sin(phi)*sin(theta);
				out[2] = cos(theta);
				vec v = d->value(in, out) ;

				file << phi << "\t" << theta << "\t" ;
				for(int u=0; u<d->dimY(); ++u)
					file << v[u] << "\t" ;
			
				file << std::endl ;
			}
	}	
	else
	{
		std::cerr << "<<ERROR>> cannot export data" << std::endl ;
	}

	return 0 ;
}
