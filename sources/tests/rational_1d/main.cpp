#include <rational_1d_fitter.h>

#include <vector>
#include <iostream>
#include <fstream>

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

	rational_1d_data data ;
	data.load(args["input"]);
	rational_1d_fitter fitter ;
	rational_1d r = fitter.fit_data(data, args.get_int("np", 10), args.get_int("nq", 10)) ;

	std::cout << r << std::endl ;

//*
	std::ofstream file(args["output"], std::ios_base::trunc);
	const float dt = (data.max() - data.min()) / 100.0f ;
	for(float x=data.min(); x<=data.max(); x+=dt)
	{
		file << x << "\t" << r(x) << std::endl ;
	}
//*/
	return 0 ;
}
