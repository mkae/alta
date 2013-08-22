#include <core/args.h>
#include <core/data.h>
#include <core/params.h>
#include <core/function.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>

#include <QApplication>

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <cmath>

int main(int argc, char** argv)
{
    QApplication app(argc, argv, false);
    arguments args(argc, argv) ;

    plugins_manager manager(args) ;

    if(args.is_defined("help")) {
        std::cout << "<<HELP>> data2moments --input data.file --output gnuplot.file --data loader.so" << std::endl ;
        std::cout << " - input, output and data are mandatory parameters" << std::endl ;
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
    if(! args.is_defined("data")) {
        std::cerr << "<<ERROR>> the data loader is not defined" << std::endl ;
        return 1 ;
    }

    // Import data
    data* d = NULL ;
    d = manager.get_data(args["data"]) ;
    d->load(args["input"]);

    // Create output file
    std::ofstream file(args["output"].c_str(), std::ios_base::trunc);

    if(d != NULL)
    {
        // Data parametrization
        params::input data_param = d->parametrization();
        int d_size = params::dimension(data_param);

        double in_angle[4] = {0.0, 0.0, 0.0, 0.0} ;

        // Sample every degree
        double dtheta = 0.5*M_PI / 90.0;

        // Moments
        vec rawm0(d->dimY());
        vec rawm1(d->dimY());

        for(int theta_in=0; theta_in<90; theta_in++)
        {
            in_angle[0] = theta_in * 0.5*M_PI / 90.0;

            // Integrate over the light hemisphere
            for(int theta_out=0; theta_out<90; theta_out++)
            {
                in_angle[2] = theta_out * 0.5*M_PI / 90.0;
                for(int phi_out=0; phi_out<180; phi_out++)
                {
                    in_angle[3] = phi_out * 0.5*M_PI / 180.0;

                    vec in(d_size);
                    params::convert(in_angle, params::SPHERICAL_TL_PL_TV_PV, data_param, &in[0]);

                    // Copy the input vector
                    vec x = d->value(in);

                    for(int i=0; i<d->dimY(); ++i)
                    {
                        double val = x[i] * cos(in_angle[2]);

                        rawm0[i] += val;
                        rawm1[i] += theta_out * val;
                    }
                }
            }

            for(int i=0; i<d->dimY(); ++i)
            {
                rawm0[i] /= 180.0*90.0;
                rawm1[i] /= 180.0*90.0 * rawm0[i];
            }

            // Output the value into the file
            file << in_angle[0] << "\t";

            for(int i=0; i<rawm0.size(); ++i)
                file << rawm0[i] << "\t";

            for(int i=0; i<rawm1.size(); ++i)
                file << rawm1[i] << "\t";
            file << std::endl;

        }

        file.close();
    }
    else
    {
        std::cerr << "<<ERROR>> load file \"" << args["input"] << "\"" << std::endl ;
    }

    return 0 ;
}
