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
        vec in_angle(4);

		// Number of elements for the integration and reconstruction
		// of the moment.
		int nb_theta_in  = 90;
		int nb_theta_out = 90;
		int nb_phi_out   = 360; // Bug here if there is a dense sampling of this domain

		// Moments
		vec raw_mnt_0(d->dimY());
		vec raw_mnt_1(d->dimY());
		vec raw_mnt_2(d->dimY());
		vec raw_mnt_3(d->dimY());

		double nb_elements = (double)(nb_phi_out*nb_theta_out);
		double weight = M_PI*M_PI/nb_elements;

		for(int theta_in=0; theta_in<nb_theta_in; ++theta_in)
        {
            in_angle[0] = theta_in * 0.5*M_PI / (double)nb_theta_in;

			// Init value for mnts
            for(int i=0; i<y_size; ++i)
			{
				raw_mnt_0[i] = 0.0;
				raw_mnt_1[i] = 0.0;
				raw_mnt_2[i] = 0.0;
				raw_mnt_3[i] = 0.0;
			}

            // Integrate over the light hemisphere
            for(int theta_out=0; theta_out<nb_theta_out; ++theta_out)
            {
                in_angle[2] = theta_out * 0.5*M_PI / (double)nb_theta_out;

                for(int phi_out=0; phi_out<nb_phi_out; ++phi_out)
                {
                    in_angle[3] = phi_out * 2.0*M_PI / (double)nb_phi_out;
					
					//! \todo Do not compute like Pascal! It does not make sense to do moments
					//! in the hemispherical parametrization. Use vectors instead (but they might
					//! degenerate to zero).
					double signed_theta_out = ((abs(in_angle[3]-M_PI) < 0.5*M_PI) ? -1.0 : 1.0) * in_angle[3];

					// Change the parametrization from SPHERICAL to BRDF's param
                    vec in(d_size);
					try
					{
						//std::cout << in_angle << std::endl;
						params::convert(&in_angle[0], params::SPHERICAL_TL_PL_TV_PV, data_param, &in[0]);
					}
					catch(...)
					{
						std::cout << "<<DEBUG>> error during conversion of " << in_angle << " to " << params::get_name(data_param) << std::endl;
						return 1;
					}

#ifdef DEBUG
					std::cout << std::endl;
					std::cout << in_angle << " -> " << in << std::endl;
#endif

                    // Evaluate the BRDF
                    vec x = d->value(in);

					// Density of samples
					double pdf = sin(in_angle[2]);

                    for(int i=0; i<y_size; ++i)
					{
						raw_mnt_0[i] += x[i] * pdf * weight;
						raw_mnt_1[i] += signed_theta_out * x[i] * pdf * weight;
						raw_mnt_2[i] += signed_theta_out*signed_theta_out * x[i] * pdf * weight;
						raw_mnt_3[i] += signed_theta_out*signed_theta_out*signed_theta_out * x[i] * pdf * weight;
					}
                }
            }
			/*
			double norm = 0.0;
			for(int i=0; i<y_size; ++i)
				norm += raw_mnt_0[i] / 3.0;
			*/
			// Normalize and center the moments before export
            for(int i=0; i<y_size; ++i)
			{
				raw_mnt_1[i] /= raw_mnt_0[i];
				raw_mnt_2[i] /= raw_mnt_0[i];
				raw_mnt_3[i] /= raw_mnt_0[i];
			}

            // Output the value into the file
            file << in_angle[0] << "\t";
            
			for(int i=0; i<raw_mnt_0.size(); ++i)
                file << raw_mnt_0[i] << "\t";

			for(int i=0; i<raw_mnt_1.size(); ++i)
                file << raw_mnt_1[i] << "\t";

			for(int i=0; i<raw_mnt_2.size(); ++i)
                file << raw_mnt_2[i] << "\t";

			for(int i=0; i<raw_mnt_3.size(); ++i)
                file << raw_mnt_3[i] << "\t";

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
