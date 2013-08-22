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
            half[0] = sin(theta_h)*sin(phi_h);
            half[1] = sin(theta_h)*cos(phi_h);
            half[2] = cos(theta_h);

            // Compute the light vector using the rotation formula.
            out[0] = sin(theta_d)*sin(phi_d);
            out[1] = sin(theta_d)*cos(phi_d);
            out[2] = cos(theta_d);

				// Rotate the diff vector to get the output vector
            rotate_binormal(out, theta_h);
            rotate_normal(out, phi_h);

            // Compute the out vector from the in vector and the half
            // vector.
            const double dot = out[0]*half[0] + out[1]*half[1] + out[2]*half[2];
            out[3] = -out[0] + 2.0*dot * half[0];
            out[4] = -out[1] + 2.0*dot * half[1];
            out[5] = -out[2] + 2.0*dot * half[2];

#ifdef DEBUG
				assert(out[2] >= 0.0 && out[5] >= 0.0);
#endif
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

			const double temp = cost * vec[0] + sint * vec[1];

            vec[1] = cost * vec[1] - sint * vec[0];
            vec[0] = temp;
        }

        //! \brief rotate a cartesian vector with respect to the bi-normal of
        //! theta degrees.
        static void rotate_binormal(double* vec, double theta)
        {
            const double cost = cos(theta);
            const double sint = sin(theta);

			const double temp = cost * vec[1] + sint * vec[2];

#ifdef DEBUG
			std::cout << acos(vec[2]) << std::endl;
#endif
            vec[2] = cost * vec[2] - sint * vec[1];
#ifdef DEBUG
			std::cout << acos(vec[2]) << std::endl;
#endif
            vec[1] = temp;
        }

		  static void print_input_params();

};

/*! \brief A parametrized object. Allow to define function object (either data
 *  or functions that are defined over an input space and output space. This
 *  Object allowas to change the parametrization of the input or output space.
 */
class parametrized
{
	public:
		parametrized() : _in_param(params::UNKNOWN_INPUT), 
		                 _out_param(params::UNKNOWN_OUTPUT) { }

		//! \brief provide the input parametrization of the object.
		virtual params::input parametrization() const
		{
			return _in_param;
		}
		
		//! \brief provide the input parametrization of the object.
		virtual params::input input_parametrization() const
		{
			return _in_param;
		}
		
		//! \brief provide the outout parametrization of the object.
		virtual params::output output_parametrization() const
		{
			return _out_param;
		}

		//! \brief can set the input parametrization of a non-parametrized
		//! object. Print an error if it is already defined.
		virtual void setParametrization(params::input new_param)
		{
			//! \todo Here is something strange happening. The equality between
			//! those enums is not correct for UNKNOWN_INPUT
			if(_in_param == new_param)
			{
				return;
			}
			else if(_in_param == params::UNKNOWN_INPUT)
			{
				_in_param = new_param;
			}
			else
			{
				std::cout << "<<ERROR>> an input parametrization is already defined: " << params::get_name(_in_param) << std::endl;
				std::cout << "<<ERROR>> trying to change to: " << params::get_name(new_param) << std::endl;
			}
		}
		
		//! \brief can set the output parametrization of a non-parametrized
		//! function. Throw an exception if it tries to erase a previously
		//! defined one.
		virtual void setParametrization(params::output new_param)
		{
			if(_out_param == new_param)
			{
				return;
			}
            else if(_out_param == params::UNKNOWN_OUTPUT)
			{
				_out_param = new_param;
			}
			else
			{
				std::cout << "<<ERROR>> an output parametrization is already defined: " << std::endl;
			}
		}

	protected:
		// Input and output parametrization
		params::input  _in_param ;
		params::output _out_param ;
};
