#include <iostream>
#include <fstream>
#include <cmath>

int main(int argc, int argv)
{
	std::ofstream f("input.gnuplot") ;

	for(int i=0; i<100; ++i)
	{
		const float x = i / (float)99.0f ;
		const float y = 100.0f * exp(-100.0f * x*x) * x*x ;
			
		f << x << "\t" << y << "\t" << 0.01f << std::endl ;
	}

	return 0 ;
}
