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
        static params::input parse_input(const std::string& txt)
        {
				for(std::map<params::input, const params::param_info>::const_iterator it=input_map.begin(); it != input_map.end(); ++it)
				{
					if(txt.compare(it->second.name) == 0)
					{
						return it->first;
					}
				}

				return params::UNKNOWN_INPUT;

			  /*
            if(txt == std::string("COS_TH"))
            {
                return params::COS_TH;
            }
            else if(txt == std::string("RUSIN_TH_TD"))
            {
                return params::RUSIN_TH_TD;
            }
            else if(txt == std::string("RUSIN_TH_PH_TD_PD"))
            {
                return params::RUSIN_TH_PH_TD_PD;
            }
            else
            {
                return params::UNKNOWN_INPUT;
            }
				*/
        }

		  static std::string get_name(const params::input param)
		  {
			  std::map<params::input, const params::param_info>::const_iterator it = input_map.find(param);
			  if(it != input_map.end())
			  {
				  return it->second.name;
			  }

			  return std::string();
		  }

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
                                 double* outvec)
        {
            switch(intype)
            {
                // 1D Parametrizations
                case params::COS_TH:
#ifndef USE_HALF
						 half_to_cartesian(acos(invec[0]), 0.0, 0.0, 0.0, outvec);
#else
						 outvec[0] = sqrt(1.0 - invec[0]*invec[0]);
						 outvec[1] = 0;
						 outvec[2] = invec[0];							
						 outvec[3] = sqrt(1.0 - invec[0]*invec[0]);
						 outvec[4] = 0;
						 outvec[5] = invec[0];							
#endif
                    break;

                // 2D Parametrizations
                case params::COS_TH_TD:
                    half_to_cartesian(acos(invec[0]), 0.0, acos(invec[1]), 0.0, outvec);
                    break;

                case params::RUSIN_TH_TD:
                    half_to_cartesian(invec[0], 0.0, invec[1], 0.0, outvec);
                    break;

                // 3D Parametrization
                case params::RUSIN_TH_PH_TD:
                    half_to_cartesian(invec[0], invec[1], invec[2], 0.0, outvec);
                    break;
                case params::RUSIN_TH_TD_PD:
                    half_to_cartesian(invec[0], 0.0, invec[1], invec[2], outvec);
                    break;

                // 4D Parametrization
                case params::RUSIN_TH_PH_TD_PD:
                    half_to_cartesian(invec[0], invec[1], invec[2], invec[3], outvec);
                    break;

                case params::SPHERICAL_TL_PL_TV_PV:
                    outvec[0] = cos(invec[1])*sin(invec[0]);
                    outvec[1] = sin(invec[1])*sin(invec[0]);
						  outvec[2] = cos(invec[0]);
                    outvec[3] = cos(invec[3])*sin(invec[2]);
                    outvec[4] = sin(invec[3])*sin(invec[2]);
						  outvec[5] = cos(invec[2]);
                    break;

                // 6D Parametrization
                case params::CARTESIAN:
                    memcpy(outvec, invec, 6*sizeof(double));
                    break;

                default:
                    throw("Transformation not implemented, params::to_cartesian");
                    break;
            }

        }

        //! \brief convert a input CARTESIAN vector, that is two 3d vectors
        //! concatenated  to an output vector in a given parametrization.
        static void from_cartesian(const double* invec, params::input outtype,
                                   double* outvec)
        {
            // Compute the half vector
            double half[3] ;
            half[0] = invec[0] + invec[3];
            half[1] = invec[1] + invec[4];
            half[2] = invec[2] + invec[5];
            double half_norm = sqrt(half[0]*half[0] + half[1]*half[1] + half[2]*half[2]);
            half[0] /= half_norm;
            half[1] /= half_norm;
            half[2] /= half_norm;
            
				// Difference vector 
				double diff[3];

            switch(outtype)
            {
                // 1D Parametrizations
                case params::COS_TH:
                    outvec[0] = half[2];
                    break;

                // 2D Parametrizations
                case params::COS_TH_TD:
                    outvec[0] = half[2];
                    outvec[1] = half[0]*outvec[0] + half[1]*outvec[1] + half[2]*outvec[2];
                    break;
                case params::RUSIN_TH_TD:
                    outvec[0] = acos(half[2]);
                    outvec[2] = acos(half[0]*outvec[0] + half[1]*outvec[1] + half[2]*outvec[2]);
                    break;

                // 3D Parametrization
                case params::RUSIN_TH_PH_TD:
                    outvec[0] = acos(half[2]);
                    outvec[1] = atan2(half[0], half[1]);
                    outvec[2] = acos(half[0]*outvec[0] + half[1]*outvec[1] + half[2]*outvec[2]);
                    break;
                
                // 4D Parametrization
                case params::RUSIN_TH_PH_TD_PD:
						  outvec[0] = acos(half[2]);
                    outvec[1] = atan2(half[0], half[1]);
						  
						  // Compute the diff vector
						  diff[0] = invec[0];
						  diff[1] = invec[1];
						  diff[2] = invec[2];
						  rotate_normal(diff, -outvec[1]);
						  rotate_binormal(diff, -outvec[0]);

						  outvec[2] = acos(diff[2]);
						  outvec[3] = atan2(diff[0], diff[1]);
                    break;

					 case params::SPHERICAL_TL_PL_TV_PV:
						  outvec[0] = acos(invec[2]);
						  outvec[1] = atan2(invec[0], invec[1]);
						  outvec[2] = acos(invec[5]);
						  outvec[3] = atan2(invec[3], invec[4]);
						  break;

                // 6D Parametrization
                case params::CARTESIAN:
                    memcpy(outvec, invec, 6*sizeof(double));
                    break;

                default:
                    assert(false);
                    break;
            }
        }

        //! \brief provide a dimension associated with a parametrization
        static int  dimension(params::input t)
        {
			  std::map<params::input, const params::param_info>::const_iterator it = params::input_map.find(t);
			  if(it != params::input_map.end())
			  {
				  return it->second.dimension;
			  }
			  else
			  {
				  return -1;
			  }
			  /*
            switch(t)
            {
                // 1D Parametrizations
                case params::COS_TH:
                    return 1;
                    break;

                // 2D Parametrizations
                case params::ISOTROPIC_TD_PD:
                case params::RUSIN_TH_TD:
                case params::ROMEIRO_TH_TD:
                case params::COS_TH_TD:
                    return 2;
                    break;

                // 3D Parametrization
                case params::RUSIN_TH_PH_TD:
                case params::RUSIN_TH_TD_PD:
                case params::ISOTROPIC_TV_TL_DPHI:
                    return 3;
                    break;

                // 4D Parametrization
                case params::RUSIN_TH_PH_TD_PD:
                case params::SPHERICAL_TL_PL_TV_PV:
                    return 4;
                    break;

                // 6D Parametrization
                case params::CARTESIAN:
                    return 6;
                    break;

                default:
                    assert(false);
                    return -1;
                    break;
            }
				*/
        }

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
            rotate_normal(out, phi_h);
            rotate_binormal(out, theta_h);

            // Compute the out vector from the in vector and the half
            // vector.
            const double dot = out[0]*half[0] + out[1]*half[1] + out[2]*half[2];
            out[3] = -out[0] + (dot+1.0) * half[0];
            out[4] = -out[1] + (dot+1.0) * half[1];
            out[5] = -out[2] + (dot+1.0) * half[2];

				assert(out[5] >= 0.0);
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

            vec[0] = cost * vec[0] - sint * vec[2];
            vec[2] = sint * vec[0] + cost * vec[2];
        }

		  static void print_input_params()
		  {
				for(std::map<params::input, const params::param_info>::const_iterator it=input_map.begin(); it != input_map.end(); ++it)
				{
					std::cout << it->second.name << std::endl;
				}
		  }


	 protected:

		  struct param_info
		  {
				param_info(std::string n, int d, std::string i) : 
					name(n), dimension(d), info(i) { };

			  std::string name;
			  int dimension;
			  std::string info;
		  };

		  static const std::map<params::input, const params::param_info> input_map;
};
