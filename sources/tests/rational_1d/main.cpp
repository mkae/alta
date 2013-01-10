#include <rational_1d_fitter.h>

#include <vector>
#include <iostream>
#include <fstream>

int main(int argc, char** argv)
{

	rational_1d_data data ;
	data.load(std::string("input.gnuplot"));
	rational_1d_fitter fitter ;
	rational_1d r = fitter.fit_data(data) ;

//*
	std::ofstream file("output.gnuplot", std::ios_base::trunc);
	for(float x=0.0f; x<=1.09f; x+=0.01f)
	{
		file << x << "\t" << r(x) << std::endl ;
	}
//*/
	return 0 ;
}
