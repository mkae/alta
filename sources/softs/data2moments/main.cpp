#include <core/args.h>
#include <core/data.h>
#include <core/params.h>
#include <core/function.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>

#include <QCoreApplication>

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <cmath>

int main(int argc, char** argv)
{
    QCoreApplication app(argc, argv, false);
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
		int y_size = d->dimY();

        double in_angle[4] = {0.0, 0.0, 0.0, 0.0} ;

        // Sample every degree
        double dtheta = 0.5*M_PI / 90.0;

        // Moments
        vec mean(d->dimY());
		vec raw_mnt_0(d->dimY());
		vec raw_mnt_1(d->dimY());

		// Number of elements for the integration and reconstruction
		// of the moment.
		int nb_theta_in  = 90;
		int nb_theta_out = 90;
		int nb_phi_out   = 360;

		double nb_elements = (double)(nb_phi_out*nb_theta_out);

		for(int theta_in=0; theta_in<nb_theta_in; theta_in++)
        {
            in_angle[0] = theta_in * 0.5*M_PI / (double)nb_theta_in;

            // Integrate over the light hemisphere
            for(int theta_out=0; theta_out<nb_theta_out; theta_out++)
            {
                in_angle[2] = theta_out * 0.5*M_PI / (double)nb_theta_out;
                for(int phi_out=0; phi_out<nb_phi_out; phi_out++)
                {
                    in_angle[3] = phi_out * 2.0*M_PI / (double)nb_phi_out - M_PI;

                    vec in(d_size);
                    params::convert(&in_angle[0], params::SPHERICAL_TL_PL_TV_PV, data_param, &in[0]);

                    // Copy the input vector
                    vec x = d->value(in);

                    for(int i=0; i<y_size; ++i)
					{
                        raw_mnt_0[i] += x[i] * sin(in_angle[2]);
						raw_mnt_1[i] += ((abs(in_angle[3]) < 0.5*M_PI) ? 1.0 : -1.0)*mean[i] * in_angle[2];
					}
                }
            }

			// Normalize and center the moments before export
            for(int i=0; i<y_size; ++i)
			{
				raw_mnt_0[i] /= nb_elements;
				raw_mnt_1[i] /= raw_mnt_0[i] * nb_elements;
			}

            // Output the value into the file
            file << in_angle[0] << "\t";
            
			for(int i=0; i<raw_mnt_0.size(); ++i)
                file << raw_mnt_0[i] << "\t";

			for(int i=0; i<raw_mnt_1.size(); ++i)
                file << raw_mnt_1[i] << "\t";

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
