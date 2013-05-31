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
		std::cout << "<<HELP>> data2diff --input data.file --output gnuplot.file --data loader.so" << std::endl ;
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

	// Import data
	data* d = NULL ;
	d = manager.get_data(args["data"]) ;
	d->load(args["input"]);

	// Create output file
	std::ofstream file(args["output"].c_str(), std::ios_base::trunc);

	if(d != NULL)
	{
		vec L(3), V(3), temp(6);
		for(int i=0; i<d->size(); ++i)
		{
			// Copy the input vector
			vec x = d->get(i);
			vec dx = x;
		
			// Convert input to CARTESIAN
			params::convert(&x[0], d->parametrization(), params::CARTESIAN, &temp[0]);
			L[0] = temp[0]; L[1] = temp[1]; L[2] = temp[2];
			V[0] = temp[3]; V[1] = temp[4]; V[2] = temp[5];
			vec y1 = d->value(L, V);
			
			// Convert perturbed input to CARTESIAN
			dx[0] += 0.01;
			params::convert(&dx[0], d->parametrization(), params::CARTESIAN, &temp[0]);
			L[0] = temp[0]; L[1] = temp[1]; L[2] = temp[2];
			V[0] = temp[3]; V[1] = temp[4]; V[2] = temp[5];
			vec y2 = d->value(L, V);


			// Print the input vector	
			for(int j=0; j<d->dimX(); ++j)
				file << x[j] << "\t";

			// And the diff values
			for(int j=0; j<d->dimY(); ++j)
				file << 100.0 * (y1[j]-y2[j]) << "\t";

			file << std::endl;
		}	

		file.close();
	}	
	else
	{
		std::cerr << "<<ERROR>> cannot export data" << std::endl ;
	}

	return 0 ;
}
