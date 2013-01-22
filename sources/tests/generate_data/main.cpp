#include <iostream>
#include <fstream>
#include <cmath>

int main(int argc, int argv)
{
	std::ofstream f("input.gnuplot") ;

	for(int i=0; i<100; ++i)
	{
		const float x = i / (float)100.0f ;
		const float y = 100.0f * exp(-10.0 * x*x) * x*x - 0.01 *x*x*x ;
			
		f << x << "\t" << y << "\t" << 0.1f << std::endl ;
	}

	return 0 ;
}
