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

#define  EPSILON 1.0E-5

vec coord(vec V, vec L, vec X, vec Y, vec N)
{
    vec pV = V-dot(V,N)*N;
    vec vCoord(2);
    vCoord[0] = dot(pV,X);
    vCoord[1] = dot(pV,Y);
    vCoord /= (1.0+dot(V,N));

    vec pL = L-dot(L,N)*N;
    vec lCoord(2);
    lCoord[0] = dot(pL,X);
    lCoord[1] = dot(pL,Y);
    lCoord /= (1.0+dot(L,N));

    if (norm(lCoord)>EPSILON)
    {
        vec lDir = normalize(lCoord);

        vec temp(2);
        temp[0] = lDir[0]*vCoord[0] + lDir[1]*vCoord[1];
        temp[1] = lDir[0]*vCoord[1] - lDir[1]*vCoord[0];
        vCoord = temp;
    }

    return vCoord;
}

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

        // Normalisation factor for the integration
        const int nb_theta_int = 90;
        const int nb_phi_int   = 360;
        const double normalization = M_PI*M_PI / (double)(nb_phi_int*nb_theta_int);

        double in_angle[4] = {0.0, 0.0, 0.0, 0.0} ;

        vec m_0(nY);
        vec m_x(nY);
        vec m_y(nY);
        vec m_xx(nY);
        vec m_xy(nY);
        vec m_yy(nY);

        // compute cumulants
        vec k_x(nY);
        vec k_y(nY);
        vec k_xx(nY);
        vec k_xy(nY);
        vec k_yy(nY);

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
                    params::convert(in_angle, params::SPHERICAL_TL_PL_TV_PV, params::STEREOGRAPHIC, &stereographics[0]);

                    // Copy the input vector
                    vec x = d->value(in);


                    for(int i=0; i<nY; ++i)
                    {
                        double val = x[i] * /* cos(in_angle[2]) */ normalization;

                        m_0[i]  += val ;
                        m_x[i]  += val * stereographics[2];
                        m_y[i]  += val * stereographics[3];
                        m_xx[i] += val * stereographics[2] * stereographics[2];
                        m_xy[i] += val * stereographics[2] * stereographics[3];
                        m_yy[i] += val * stereographics[3] * stereographics[3];
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

            for(int i=0; i<nY; ++i)
                file << k_y[i] << "\t";

            for(int i=0; i<nY; ++i)
                file << k_xx[i] << "\t";

            for(int i=0; i<nY; ++i)
                file << k_xy[i] << "\t";

            for(int i=0; i<nY; ++i)
                file << k_yy[i] << "\t";

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
