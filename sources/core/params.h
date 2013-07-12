#pragma once

#include <string>
#include <map>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <iostream>

/*! \brief a static class allowing to change from one parametrization
 *  to another.
 *  \ingroup core
 *
 *  \details
 *  Any \typedef function object or \typedef data object should have
 *  an associated parametrization.
 */
class params
{
    public: // data

		 //! \brief list of all supported parametrization for the input space.
		 //! An unsupported parametrization will go under the name
		 //! <em>unknown</em>.
		 enum input
		 {
			 ROMEIRO_TH_TD,
			 RUSIN_TH_TD,
			 RUSIN_TH_PH_TD,
			 RUSIN_TH_TD_PD,
			 RUSIN_TH_PH_TD_PD,
			 COS_TH,
			 COS_TH_TD,
			 ISOTROPIC_TV_TL_DPHI,
			 ISOTROPIC_TD_PD, // Difference between two directions such as R and H
			 CARTESIAN,
			 SPHERICAL_TL_PL_TV_PV,
			 UNKNOWN_INPUT
		 };

		 //! \brief list of all supported parametrization for the output space.
		 //! An unsupported parametrization will go under the name
		 //! <em>unknown</em>.
		 enum output
		 {
			 INV_STERADIAN,
			 ENERGY,
			 RGB_COLOR,
			 XYZ_COLOR,
			 UNKNOWN_OUTPUT
		 };

    public: // methods

        //! \brief parse a string to provide a parametrization type.
        static params::input parse_input(const std::string& txt);

		  //! \brief look for the string associated with a parametrization
		  //! type.
		  static std::string get_name(const params::input param);

        //! \brief parse a string to provide a parametrization type.
        static params::output parse_output(const std::string& txt)
        {
            if(txt == std::string("ENERGY"))
            {
                return params::ENERGY;
            }
            else
            {
                return params::UNKNOWN_OUTPUT;
            }
        }

        //! \brief static function for input type convertion. This
        //! function allocate the resulting vector.
        static double* convert(const double* invec, params::input intype,
                               params::input outtype)
        {
            int dim = dimension(outtype); // Get the size of the output vector

            if(dim > 0)
            {
                double* outvec = new double[dim];
                double  temvec[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // Temp CARTESIAN vectors
                to_cartesian(invec, intype, temvec);
                from_cartesian(temvec, outtype, outvec);

                return outvec;
            }
            else
            {
                return NULL;
            }
        }

        //! \brief static function for input type convertion. The outvec
        //! resulting vector should be allocated with the correct
        //! output size.
        static void convert(const double* invec, params::input intype,
                            params::input outtype, double* outvec)
        {
            double  temvec[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // Temp CARTESIAN vectors
            to_cartesian(invec, intype, temvec);
#ifdef DEBUG
            std::cout << "<<DEBUG>> temp vec = ["
                      << temvec[0] << ", " << temvec[1] << ", " << temvec[2] << "] => ["
                      << temvec[3] << ", " << temvec[4] << ", " << temvec[5] << "]" << std::endl;
#endif
            from_cartesian(temvec, outtype, outvec);
        }

        //! \brief convert a input vector in a given parametrization to an
        //! output vector in a cartesian parametrization, that is two 3d
        //! vectors concatenated.
        static void to_cartesian(const double* invec, params::input intype,
                                 double* outvec);

        //! \brief convert a input CARTESIAN vector, that is two 3d vectors
        //! concatenated  to an output vector in a given parametrization.
        static void from_cartesian(const double* invec, params::input outtype,
                                   double* outvec);

        //! \brief provide a dimension associated with a parametrization
        static int  dimension(params::input t);

        //! \brief provide a dimension associated with a parametrization
        static int  dimension(params::output t)
        {
            switch(t)
            {
                // 1D Parametrizations
                case params::INV_STERADIAN:
                case params::ENERGY:
                    return 1;
                    break;

                // 3D Parametrization
                case params::RGB_COLOR:
                case params::XYZ_COLOR:
                    return 3;
                    break;

                default:
                    assert(false);
                    return -1;
                    break;
            }
        }

        //! \brief from the 4D definition of a half vector parametrization,
        //! export the cartesian coordinates.
        static void half_to_cartesian(double theta_h, double phi_h,
                                      double theta_d, double phi_d, double* out)
        {
            // Calculate the half vector
            double half[3];
            half[0] = sin(theta_h)*cos(phi_h);
            half[1] = sin(theta_h)*sin(phi_h);
            half[2] = cos(theta_h);

            // Compute the light vector using the rotation formula.
            out[0] = sin(theta_d)*cos(phi_d);
            out[1] = sin(theta_d)*sin(phi_d);
            out[2] = cos(theta_d);

				//! \todo investigate here, the rotation along N should be
				//1 of phi_h not theta_h !
            rotate_binormal(out, theta_h);
            rotate_normal(out, phi_h);

            // Compute the out vector from the in vector and the half
            // vector.
            const double dot = out[0]*half[0] + out[1]*half[1] + out[2]*half[2];
            out[3] = -out[0] + (dot+1.0) * half[0];
            out[4] = -out[1] + (dot+1.0) * half[1];
            out[5] = -out[2] + (dot+1.0) * half[2];

				assert(out[5] >= 0.0);
        }
			
        //! \brief from the 4D definition of a classical vector parametrization,
        //! export the cartesian coordinates.
		  static void classical_to_cartesian(double theta_l, double phi_l, 
		                                     double theta_v, double phi_v, double* out)
		  {
			  out[0] = cos(phi_l)*sin(theta_l);
			  out[1] = sin(phi_l)*sin(theta_l);
			  out[2] = cos(theta_l);
			  out[3] = cos(phi_v)*sin(theta_v);
			  out[4] = sin(phi_v)*sin(theta_v);
			  out[5] = cos(theta_v);
		  }

        //! \brief rotate a cartesian vector with respect to the normal of
        //! theta degrees.
        static void rotate_normal(double* vec, double theta)
        {
            const double cost = cos(theta);
            const double sint = sin(theta);

            vec[0] = cost * vec[0] - sint * vec[1];
            vec[1] = sint * vec[0] + cost * vec[1];
        }

        //! \brief rotate a cartesian vector with respect to the bi-normal of
        //! theta degrees.
        static void rotate_binormal(double* vec, double theta)
        {
            const double cost = cos(theta);
            const double sint = sin(theta);

            vec[1] = cost * vec[1] - sint * vec[2];
            vec[2] = sint * vec[1] + cost * vec[2];
        }

		  static void print_input_params();

};
