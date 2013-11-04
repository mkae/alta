#include <core/args.h>
#include <core/data.h>
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
		std::cout << "<<HELP>> data2gnuplot --input data.file --output gnuplot.file --data loader.so" << std::endl ;
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
    d = plugins_manager::get_data(args["data"]) ;
	d->load(args["input"]);

	// Create output file
	std::ofstream file(args["output"].c_str(), std::ios_base::trunc);

	if(d != NULL)
	{
		std::cout << "<<INFO>> will export " << d->size() << " elements" << std::endl ;

		double theta_in = (double)args.get_float("theta", 0.0f);
		double phi_in   = (double)args.get_float("phi", 0.0f);
		vec in(3), out(3) ;
		in[0] = cos(phi_in)*sin(theta_in);
		in[1] = sin(phi_in)*sin(theta_in);
		in[2] = cos(theta_in);

		vec _min = d->min();
		vec _max = d->max();

		vec x(d->dimX());
		x = _min + 0.5*(_max - _min);

		const int N = 100;
		for(int i=0; i<N; ++i)
			for(int j=0; j<N; ++j)
			{
				x[0] = _min[0] + (_max[0] - _min[0]) * double(i) / double(N-1);
				x[1] = _min[1] + (_max[1] - _min[1]) * double(j) / double(N-1);
				vec v = d->value(x) ;

				for(int u=0; u<d->dimX(); ++u)
				{
					file << x[u] << "\t";
				}
				for(int u=0; u<d->dimY(); ++u)
				{
					file << v[u] << "\t" ;
				}

				file << std::endl ;
			}
	}	
	else
	{
		std::cerr << "<<ERROR>> cannot export data" << std::endl ;
	}

	return 0 ;
}
