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
		std::cout << "<<HELP>> data2gnuplot --input data.file --output gnuplot.file --data loader.so --slice [0, 1, 2]" << std::endl ;
		std::cout << " - input and output are mandatory parameters" << std::endl ;
		std::cout << " - if no data parameter is set, the data will be loaded using vertical segments plugin" << std::endl ;
		std::cout << " - the slice parameter allows to select the dimensions to output" << std::endl ;
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

	data* d = plugins_manager::get_data(args["data"]) ;
	d->load(args["input"], args);

	// Create output file
	std::ofstream file(args["output"].c_str(), std::ios_base::trunc);

	// Slices of the data to be extracted
	std::vector<int> slice = args.get_vec<int>("slice");
	assert(slice.size() <= d->dimX());
	if(slice.size() == 0)
	{
		slice.resize(d->dimX());
		for(int u=0; u<d->dimX(); ++u)
		{
			slice[u] = u;
		}
	}
	else
	{
		for(int u=0; u<slice.size(); ++u)
		{
			assert(slice[u] < d->dimX());
		}
	}

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

		const int N = 10000;
		const int n = int(pow(N, 1.0/slice.size()));

		for(int i=0; i<N; ++i)
		{
			// Set the point of evaluation and get the value
			for(int u=0; u<d->dimX(); ++u)
			{
				int j = (i / int(pow(n, u))) % n;
				int v = slice[u];
				x[v] = _min[v] + (_max[v] - _min[v]) * double(j) / double(n);
			}
			vec v = d->value(x) ;

			// Copy the value in the file
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
