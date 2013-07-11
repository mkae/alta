#include <iostream>
#include <fstream>
#include <cmath>

#include <core/args.h>
#include <core/params.h>

int main(int argc, char** argv)
{
	std::ofstream f("input.gnuplot") ;
	arguments args(argc, argv);	

    std::cout.precision(10);

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
        f << "#PARAM_IN UNKNOWN" << std::endl;
        //f << "#VS 2" << std::endl;
		for(int i=0; i<nbx; ++i)
		{
			const float x = i / (float)nbx ;
            //const float y = 100.0f * exp(-10.0 * x*x) * x*x - 0.01 *x*x*x + 0.1 ;
            const float y = (1.0) / (1.0E-10 + x*x*x);

            f << x << "\t" << y << "\t" << y*0.9f << "\t" << y*1.1f << std::endl ;
		}
	}
    else if(k == 2)
    {
        f << "#DIM 1 1" << std::endl ;
        f << "#PARAM_IN UNKNOWN" << std::endl;
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
                const float z = 10 * x + 1.0;
			
                f << x << "\t" << y << "\t" << z << "\t" << z-0.1f << "\t" << z << std::endl ;
			}
	}
    else if(k == 4)
	{
		f << "#DIM 2 1" << std::endl ;
        f << "#PARAM_IN UNKNOWN" << std::endl;
		for(int i=0; i<nbx; ++i)
			for(int j=0; j<nby; ++j)
			{
				const float x = i / (float)nbx ;
				const float y = j / (float)nby ; 
                const float z = x*y / (1.0E-3 + x*x*x) + 10.;

			
                f << x << "\t" << y << "\t" << z << std::endl ;
			}
	}
   else if(k == 5)
	{
		f << "#DIM 1 1" << std::endl ;
		f << "#PARAM_IN COS_TH" << std::endl;
		for(int i=0; i<nbx; ++i)
		{
                const double   x = i / (float)nbx ;
//				const float d[3] = {0.1, 0.0, 0.5};
				const float d[3] = {0.0, 0.0, 0.0};
				const float z1 = d[0] + 0.2 * std::pow(x, 1.5) ;
//				const float z2 = d[1] - 0.1 * std::pow(x, 4.0) ;
//				const float z3 = d[2] + 0.7 * std::pow(x, 1.0) ;
			
				f << x << "\t" << z1 ;
//				f << "\t" << z2 ;
//				f << "\t" << z3 ;
				f << std::endl ;
		}
	}
   else if(k == 6)
	{
		f << "#DIM 1 1" << std::endl ;
		f << "#PARAM_IN COS_TH" << std::endl;
		for(int i=0; i<nbx; ++i)
		{
                const double x = i / (float)nbx ;
                const double z = 0.1 + 0.5 * std::pow(x, 1.5) ;
			
				f << x << "\t" << z << std::endl ;
		}
	}
    // Lafortune fitting
    // Single lobe (0.86, 0.77, 18.6)
    else if(k == 7)
     {
         f << "#DIM 2 1" << std::endl ;
         f << "#PARAM_IN RUSIN_TH_TD" << std::endl;
         for(int i=0; i<nbx; ++i)
         {
             for(int j=0; j<nby; ++j)
             {
                 double in_r[2], in_c[6];
                 in_r[0] = M_PI * 0.5 * i / (float)nbx ;
                 in_r[1] = M_PI * 0.5 * j / (float)nby ;

                 params::convert(in_r, params::RUSIN_TH_TD, params::CARTESIAN, in_c);

                 const double Cx =0.86;
                 const double Cz =0.77;
                 const double n = 18.6;

                 const double cos = Cx * (in_c[0]*in_c[3] + in_c[1]*in_c[4]) + Cz*in_c[2]*in_c[5];

                 if(cos > 0.0)
                 {
                    const double z = std::pow(cos, n) ;

                    f << in_r[0] << "\t" << in_r[1] << "\t" << z << std::endl ;
                 }
             }
         }
     }
    // Lafortune fitting
    // Triple lobe (0.86, 0.77, 18.6), (-0.41, 0.018, 2.58), (-1.03, 0.7, 63.8)
    else if(k == 8)
     {
         f << "#DIM 2 1" << std::endl ;
         f << "#PARAM_IN RUSIN_TH_TD" << std::endl;
         for(int i=0; i<nbx; ++i)
         {
             for(int j=0; j<nby; ++j)
             {
                 double in_r[2], in_c[6];
                 in_r[0] = M_PI * 0.5 * i / (float)nbx ;
                 in_r[1] = M_PI * 0.5 * j / (float)nby ;

                 params::convert(in_r, params::RUSIN_TH_TD, params::CARTESIAN, in_c);

                 const double Cx[3] = {0.86, -0.410, -1.03};
                 const double Cz[3] = {0.77,  0.018,  0.70};
                 const double n[3]  = {18.6,  2.580,  63.8};

                 double z = 0.0;
                 for(int k=0; k<3; ++k)
                 {
                 const double cos = Cx[k] * (in_c[0]*in_c[3] + in_c[1]*in_c[4]) + Cz[k]*in_c[2]*in_c[5];

                 if(cos > 0.0)
                 {
                    z += std::pow(cos, n[k]) ;

                 }
                 }
                   f << in_r[0] << "\t" << in_r[1] << "\t" << z << std::endl ;
             }
         }
    }

	return 0 ;
}
