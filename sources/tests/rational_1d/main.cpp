#include <rational_1d_fitter_cgal.h>
#include <rational_1d_fitter_eigen.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <limits>

#include <core/args.h>

int main(int argc, char** argv)
{
	arguments args(argc, argv) ;
	if(args.is_defined("help")) {
		std::cout << argv[0] << " --np <int> --nq <int> --input <filename> --output <filename>" << std::endl ;
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

	// Load the data
	rational_1d_data data ;
	if(args.is_defined("min") && args.is_defined("max"))
	{
		data.load(args["input"], args.get_float("min", 0.0f), args.get_float("max", 1.0f));
	}
	else if(args.is_defined("min") && !args.is_defined("max"))
	{
		data.load(args["input"], args.get_float("min", 0.0f), std::numeric_limits<double>::max());
	}
	else if(args.is_defined("min") && !args.is_defined("max"))
	{
		data.load(args["input"], -std::numeric_limits<double>::max(), args.get_float("min", 0.0f));
	}
	else
	{
		data.load(args["input"]);
	}

	// Fitting call
	rational_1d_fitter* fitter ;
	if(args.is_defined("algorithm") && args["algorithm"] == std::string("eigen"))
	{
		std::cout << "<<INFO>> using Eigen method" << std::endl ;
		fitter = new rational_1d_fitter_eigen() ;
	}
	else
	{
		std::cout << "<<INFO>> using CGAL method" << std::endl ;
		fitter = new rational_1d_fitter_cgal() ;
	}
	fitter->set_parameters(args.get_int("min-np", 10), args.get_int("np", 10), args.get_int("min-nq", 10), args.get_int("nq", 10)) ;
	
	rational_1d r ;
	bool is_fitted = fitter->fit_data(data, r) ;

	// Display the result
	if(is_fitted)
	{
		std::cout << r << std::endl ;

		//*
		std::ofstream file(args["output"], std::ios_base::trunc);
		const double dt = (data.max() - data.min()) / 100.0f ;
		for(double x=data.min(); x<=data.max(); x+=dt)
		{
			file << x << "\t" << r(x) << std::endl ;
		}
		//*/
	}
	else
	{
		std::cout << "<<ERROR>> unable to fit the data" << std::endl ;
	}

	// Clean data
	delete fitter ;

	return 0 ;
}
