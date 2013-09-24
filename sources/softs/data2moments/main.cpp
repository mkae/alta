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
        const int nX = d->dimX();
        const int nY = d->dimY();

        // Options
        bool with_cosine = args.is_defined("cosine");
        bool use_angular = args.is_defined("angular-moments");

        // Normalisation factor for the integration
        const int nb_theta_int = 90;
        const int nb_phi_int   = (use_angular) ? 2 : 360;
        const double normalization = ((use_angular) ? M_PI : M_PI*M_PI) / (double)(nb_phi_int*nb_theta_int);

        double in_angle[4] = {0.0, 0.0, 0.0, 0.0} ;


        // Raw moments
        vec m_0(nY);
        vec m_x(nY);
        vec m_y(nY);
        vec m_xx(nY);
        vec m_xy(nY);
        vec m_yy(nY);

        // Cumulants
        vec k_x(nY);
        vec k_y(nY);
        vec k_xx(nY);
        vec k_xy(nY);
        vec k_yy(nY);

        // The X and Y directions to compute the moments. This is an argument of the command
        // line. The different directions can be: using stereographic coordinates, using the
        // theta of the classical parametrization (the second coordinate is then 0).
        vec xy(2);

        for(int theta_in=0; theta_in<90; theta_in++)
        {
            in_angle[0] = theta_in * 0.5*M_PI / 90.0;

            m_0  = vec::Zero(nY);
            m_x  = vec::Zero(nY);
            m_y  = vec::Zero(nY);
            m_xx = vec::Zero(nY);
            m_xy = vec::Zero(nY);
            m_yy = vec::Zero(nY);

            // Theta angle in [0 .. PI / 2]
            for(int theta_out=0; theta_out<nb_theta_int; theta_out++)
            {
                in_angle[2] = theta_out * 0.5*M_PI / (double)nb_theta_int;

                // Azimutal angle in [-PI .. PI]
                for(int phi_out=0; phi_out<nb_phi_int; phi_out++)
                {
                    in_angle[3] = phi_out * 2.0*M_PI / (double)nb_phi_int - M_PI;

                    vec in(nX), stereographics(4);
                    params::convert(in_angle, params::SPHERICAL_TL_PL_TV_PV, data_param, &in[0]);

                    if(use_angular)
                    {
                        xy[0] = (std::abs(in_angle[3]) < 0.5*M_PI) ? in_angle[2] : -in_angle[2];
                        xy[1] = 0.0;
                    }
                    else
                    {
                        params::convert(in_angle, params::SPHERICAL_TL_PL_TV_PV, params::STEREOGRAPHIC, &stereographics[0]);
                        xy[0] = stereographics[2];
                        xy[1] = stereographics[3];
                    }

                    // Copy the input vector
                    vec x = d->value(in);


                    for(int i=0; i<nY; ++i)
                    {
                        double val = x[i] * normalization;
                        if(!use_angular) {
                            val *= sin(in_angle[2]);
                        }
                        if(with_cosine)
                        {
                            val *= cos(in_angle[2]);
                        }

                        m_0[i]  += val ;
                        m_x[i]  += val * xy[0];
                        m_y[i]  += val * xy[1];
                        m_xx[i] += val * xy[0] * xy[0];
                        m_xy[i] += val * xy[0] * xy[1];
                        m_yy[i] += val * xy[1] * xy[1];
                    }
                }
            }


            for(int i=0; i<nY; ++i)
            {
                m_x[i]  /= m_0[i];
                m_y[i]  /= m_0[i];
                m_xx[i] /= m_0[i];
                m_xy[i] /= m_0[i];
                m_yy[i] /= m_0[i];
            }

            // compute cumulants
            k_x  = m_x;
            k_y  = m_y;
            k_xx = m_xx - product(m_x,m_x);
            k_xy = m_xy - product(m_x,m_y);
            k_yy = m_yy - product(m_y,m_y);

            // Output the value into the file
            file << in_angle[0] << "\t";

            for(int i=0; i<nY; ++i)
                file << m_0[i] << "\t";

            for(int i=0; i<nY; ++i)
                file << k_x[i] << "\t";

            if(!use_angular) {
                for(int i=0; i<nY; ++i) {
                    file << k_y[i] << "\t";
                }
            }

            for(int i=0; i<nY; ++i)
                file << k_xx[i] << "\t";

            if(!use_angular) {
                for(int i=0; i<nY; ++i) {
                    file << k_xy[i] << "\t";
                }

                for(int i=0; i<nY; ++i) {
                    file << k_yy[i] << "\t";
                }
            }
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
