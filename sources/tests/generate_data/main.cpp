#include <iostream>
#include <fstream>
#include <cmath>

int main(int argc, int argv)
{
	std::ofstream f("input.gnuplot") ;

	const int k = 2 ;
	if(k == 1)
	{
		f << "#DIM 1 1" << std::endl ;
		for(int i=0; i<100; ++i)
		{
			const float x = i / (float)100.0f ;
			const float y = 100.0f * exp(-10.0 * x*x) * x*x - 0.01 *x*x*x ;

			f << x << "\t" << y << "\t" << 0.1f << std::endl ;
		}
	}
	else if(k == 2)
	{
		f << "#DIM 2 1" << std::endl ;
		for(int i=0; i<100; ++i)
			for(int j=0; j<100; ++j)
			{
				const float x = i / (float)100.0f ;
				const float y = j / (float)100.0f ; 
				const float z = 100.0f * exp(-10.0 * x*x) * y*y - 0.01 *y*x*y ;
			
				f << x << "\t" << y << "\t" << z << "\t" << 0.1f << std::endl ;
			}
	}

	return 0 ;
}
