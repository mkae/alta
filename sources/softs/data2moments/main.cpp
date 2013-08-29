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
		  const int nb_phi_int   = 180;
		  const double normalization = M_PI*M_PI / (double)(nb_phi_int*nb_theta_int);

        double in_angle[4] = {0.0, 0.0, 0.0, 0.0} ;

		  // Static data
		  vec X(3); X[0] = 1; X[1] = 0; X[2] = 0;
		  vec Y(3); Y[0] = 0; Y[1] = 1; Y[2] = 0;
		  vec N(3); N[0] = 0; N[1] = 0; N[2] = 1;


        for(int theta_in=0; theta_in<90; theta_in++)
        {
            in_angle[0] = theta_in * 0.5*M_PI / 90.0;

            vec m_0(nY);  m_0 = vec::Zero(nY);
            vec m_x(nY);  m_0 = vec::Zero(nY);
            vec m_y(nY);  m_0 = vec::Zero(nY);
            vec m_xx(nY); m_0 = vec::Zero(nY);
            vec m_xy(nY); m_0 = vec::Zero(nY);
            vec m_yy(nY); m_0 = vec::Zero(nY);

            // Integrate over the light hemisphere
            for(int theta_out=0; theta_out<nb_theta_int; theta_out++)
            {
                in_angle[2] = theta_out * 0.5*M_PI / (double)nb_theta_int;
                for(int phi_out=0; phi_out<nb_phi_int; phi_out++)
                {
                    in_angle[3] = phi_out * 2.0*M_PI / nb_phi_int;

                    vec in(nX), cart(6), L(3), V(3);
                    params::convert(in_angle, params::SPHERICAL_TL_PL_TV_PV, data_param, &in[0]);
                    params::convert(in_angle, params::SPHERICAL_TL_PL_TV_PV, params::CARTESIAN, &cart[0]);
						  L[0] = cart[0];
						  L[1] = cart[1];
						  L[2] = cart[2];
						  V[0] = cart[3];
						  V[1] = cart[4];
						  V[2] = cart[5];

                    // Copy the input vector
                    vec x = d->value(in);

                    // Get the projected 2D coordinate
                    vec xy = coord(V, L, X, Y, N);

                    for(int i=0; i<nY; ++i)
                    {
                        double val = x[i] /* cos(in_angle[2])*/ * normalization;

                        m_0[i]  += val;
                        m_x[i]  += val * xy[0];
                        m_y[i]  += val * xy[1];
                        m_xx[i] += val * xy[0] * xy[0];
                        m_xy[i] += val * xy[0] * xy[1];
                        m_yy[i] += val * xy[1] * xy[1];
                    }
                }
            }

/*
            for(int i=0; i<nY; ++i)
            {
                m_x[i]  /= m_0[i];
                m_y[i]  /= m_0[i];
                m_xx[i] /= m_0[i];
                m_xy[i] /= m_0[i];
                m_yy[i] /= m_0[i];
            }
*/
            // compute cumulants
            vec k_x  = m_x;
            vec k_y  = m_y;
            vec k_xx = m_xx;// - product(m_x,m_x);
            vec k_xy = m_xy;// - product(m_x,m_y);
            vec k_yy = m_yy;// - product(m_y,m_y);

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
