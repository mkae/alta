#include <iostream>
#include <fstream>
#include <cmath>

#include <core/args.h>

int main(int argc, char** argv)
{
	std::ofstream f("input.gnuplot") ;
	arguments args(argc, argv);	

	int nbx = 100;
	int nby = 100;
	int nbz = 100;
	if(args.is_defined("nbx"))
		nbx = args.get_int("nbx", 100) ;
	if(args.is_defined("nby"))
		nby = args.get_int("nby", 100) ;

	const int k = args.get_int("f", 1) ;
	if(k == 1)
	{
		f << "#DIM 1 1" << std::endl ;
		for(int i=0; i<nbx; ++i)
		{
			const float x = i / (float)nbx ;
            const float y = 100.0f * exp(-10.0 * x*x) * x*x - 0.01 *x*x*x + 0.1 ;

			f << x << "\t" << y << "\t" << 0.1f << std::endl ;
		}
	}
    else if(k == 2)
    {
        f << "#DIM 1 1" << std::endl ;
        for(int i=0; i<nbx; ++i)
        {
            const float x = i / (float)nbx ;
            const float y = (1.0 + 7.0*x - 10.5*x*x) / (1.0 + 7.0 * x) ;

            f << x << "\t" << y << "\t" << 0.1f << std::endl ;
        }
    }
    else if(k == 3)
	{
		f << "#DIM 2 1" << std::endl ;
		for(int i=0; i<nbx; ++i)
			for(int j=0; j<nby; ++j)
			{
				const float x = i / (float)nbx ;				
				const float y = j / (float)nby ; 
				const float z = 1 + 0.1f*x;
			
				f << x << "\t" << y << "\t" << z << "\t" << 0.1f << std::endl ;
			}
	}
    else if(k == 4)
	{
		f << "#DIM 2 1" << std::endl ;
		for(int i=0; i<nbx; ++i)
			for(int j=0; j<nby; ++j)
			{
				const float x = i / (float)nbx ;
				const float y = j / (float)nby ; 
				const float z = exp(-10.0 * x*x) +  x*y ;
			
				f << x << "\t" << y << "\t" << z << "\t" << 0.1f << std::endl ;
			}
	}
   else if(k == 5)
	{
		f << "#DIM 1 3" << std::endl ;
		f << "#PARAM_IN COS_TH" << std::endl;
		for(int i=0; i<nbx; ++i)
		{
				const float x = i / (float)nbx ;
				const float z1 = 0.1 + 0.5 * std::pow(x, 1.5) ;
				const float z2 = 0.0 - 0.1 * std::pow(x, 4.0) ;
				const float z3 = 0.5 + 0.7 * std::pow(x, 1.0) ;
			
				f << x << "\t" << z1 << "\t" << z2 << "\t" << z3 << std::endl ;
		}
	}
   else if(k == 6)
	{
		f << "#DIM 1 1" << std::endl ;
		f << "#PARAM_IN COS_TH" << std::endl;
		for(int i=0; i<nbx; ++i)
		{
				const float x = i / (float)nbx ;
				const float z = 0.1 + 0.5 * std::pow(x, 1.5) ;
			
				f << x << "\t" << z << std::endl ;
		}
	}

	return 0 ;
}
