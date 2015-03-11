/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2014 CNRS
   Copyright (C) 2013, 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

#include <string>
#include <map>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <iostream>

#include "common.h"

/*! \class params
 *  \ingroup core
 *  \brief a static class allowing to change from one parametrization
 *  to another.
 *
 *  Any function object or data object should have an associated
 *  parametrization.
 *
 *  We use the following convention to defined the tangent, normal and
 *  bi-normal of the surface:
 *   * The normal is the upper vector (0, 0, 1)
 *   * The tangent direction is along x direction (1, 0, 0)
 *   * The bi-normal is along the y direction (0, 1, 0)
 */
class params
{
    public: // data

		 //! \brief list of all supported parametrization for the input space.
		 //! An unsupported parametrization will go under the name
		 //! <em>unknown</em>. We use the following notations:
		 //!   * The view vector is \f$\vec{v}\f$
		 //!   * The light vector is \f$\vec{l}\f$
		 //!   * The normal vector is \f$\vec{n}\f$
		 //!   * The reflected vector is \f$\vec{r} = 2\mbox{dot}(\vec{v}, \vec{n})\vec{n} - \vec{v}\f$
		 enum input
		 {
       RUSIN_TH_PH_TD_PD,     /*!< Half-angle parametrization as described in [Rusinkiewicz'98] */
       RUSIN_TH_PH_TD,
       RUSIN_TH_TD_PD,
       RUSIN_TH_TD,           /*!< Half-angle parametrization with no azimutal information */
       RUSIN_VH_VD,           /*!< Half-angle parametrization in vector format. Coordinates are:
                               [\f$\vec{h}_x, \vec{h}_y, \vec{h}_z, \vec{d}_x, \vec{d}_y, 
  									  \vec{d}_z \f$].*/
       RUSIN_VH,              /*!< Half-angle parametrization with no difference direction in 
  								     vector format. Coordinates are: [\f$\vec{h}_x, \vec{h}_y, 
  									  \vec{h}_z\f$]. */
       COS_TH_TD,      /*!< Cosine of the RUSIN_TH_TD parametrization: Coordinates are in 
                         $[\cos_\theta_h,\cos_\theta_d]$. */
       COS_TH,

       SCHLICK_TK_PK,         /*!< Schlick's back vector parametrization */
       SCHLICK_VK,            /*!< Schlick's back vector */
       SCHLICK_TL_TK_PROJ_DPHI,/*!< 3D Parametrization where the phi component is projected and
                              the parametrization is centered around the back direction.
  									 \f$[\theta_L, x, y] = [\theta_L, \theta_K \cos(\phi_K), \theta_K \sin(\phi_K)]\f$*/
       COS_TK,                /*!< Schlick's back vector dot product with the normal */


       RETRO_TL_TVL_PROJ_DPHI,/*!< 2D Parametrization where the phi component is projected and
                               the parametrization is centered around the retro direction
  									  \f$[x, y] = [\theta_{VL} \cos(\Delta\phi), \theta_{VL} 
  									  \sin(\Delta\phi)]\f$.*/

       STEREOGRAPHIC,         /*!< Stereographic projection of the Light and View vectors */


       SPHERICAL_TL_PL_TV_PV, /*!< Light and View vectors represented in spherical coordinates */
       COS_TLV,               /*!< Dot product between the Light and View vector */
       COS_TLR,               /*!< Dot product between the Light and Reflected vector */
       ISOTROPIC_TV_TL,       /*!< Light and View vectors represented in spherical coordinates, */
       ISOTROPIC_TV_TL_DPHI,  /*!< Light and View vectors represented in spherical coordinates,
                                   with the difference of azimutal coordinates in the last component  */
       ISOTROPIC_TV_PROJ_DPHI,/*!< 2D Parametrization where the phi component is projected.
                               Coordinates are: [\f$\theta_v \cos(\Delta\phi), \theta_v 
  									  \sin(\Delta\phi).\f$]*/
       ISOTROPIC_TL_TV_PROJ_DPHI,/*!< 3D Parametrization where the phi component is projected.
                                  Coordinates are: [\f$\theta_l, \theta_v \cos(\Delta\phi), 
  										  \theta_v \sin(\Delta\phi).\f$]*/
       ISOTROPIC_TD_PD,       /*!< Difference between two directions such as R and H */

       BARYCENTRIC_ALPHA_SIGMA, /*!< Barycentric parametrization defined in Stark et al. [2004].
                                 Coordinates are: \f$[\alpha, \sigma] = [{1\over 2}(1 - \vec{l}\vec{v}), 
  										 (1-(\vec{h}.\vec{n})^2)(1 - \alpha)]\f$ */

       CARTESIAN,             /*!< View and Light vectors represented in cartesian coordinates.
                               We always pack the view vector first: \f$\vec{c} = [v.x, v.y, 
  									  v.z, l.x, l.y, l.z] \f$*/

       UNKNOWN_INPUT          /*!< Default behaviour. Only use this is you do not fit BRDF data */
		 };

		 //! \brief list of all supported parametrization for the output space.
		 //! An unsupported parametrization will go under the name
		 //! <em>unknown</em>.
		 enum output
		 {
			 INV_STERADIAN,                /*!< Output values in inverse steradian (sr-1). 
														   This is the standard definition for a BRDF. */
			 INV_STERADIAN_COSINE_FACTOR,  /*!< Output values in inverse steradian (sr-1)
			                                    weighted by the cosine factor of the output
															direction. */
			 ENERGY,
			 RGB_COLOR,
			 XYZ_COLOR,
			 UNKNOWN_OUTPUT
		 };

    public: // methods

        //! \brief parse a string to provide a parametrization type.
        static params::input parse_input(const std::string& txt);

        //! \brief parse a string to provide a parametrization type.
        static params::output parse_output(const std::string& txt);

		  //! \brief look for the string associated with a parametrization
		  //! type.
		  static std::string get_name(const params::input param);
		  
		  //! \brief look for the string associated with a parametrization
		  //! type.
		  //! \todo Finish this implementation. It requires another static
		  //! object.
		  static std::string get_name(const params::output)
		  {
			  return std::string("UNKNOWN_OUTPUT");
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
        //! TODO:  RP this function is weird the outtype and intype ARE NEVER USED
        static void convert(const double* invec, params::output intype, int indim,
                            params::output outtype, int outdim, double* outvec)
        {
        	// The convertion is done using the cartesian parametrization as
    			// an intermediate one. If the two parametrizations are equals
    			// there is no need to perform the conversion.
    			if(outdim == indim)
    			{
    				for(int i=0; i<outdim; ++i) { outvec[i] = invec[i]; }
    			}
    			// If the dimension of the output data is bigger than the 
    			// dimensions of the input domain, and the input domain is of
    			// dimension one, spread the data over all dimensions.
    			else if(indim == 1)
    			{
    				for(int i=0; i<outdim; ++i) { outvec[i] = invec[0]; }
    			}
    			// If the output dimension is one, compute the average of the
    			// input vector values.
    			else if(outdim == 1)
    			{
    				for(int i=0; i<indim; ++i) 
            { 
                outvec[0] += invec[i]; 
            }	
            outvec[0] /= static_cast<double>(indim);
    			}
    			else
    			{
    				NOT_IMPLEMENTED();
    			}
        }

        //! \brief static function for input type convertion. The outvec
        //! resulting vector should be allocated with the correct
        //! output size.
        static void convert(const double* invec, params::input intype,
                            params::input outtype, double* outvec)
        {
			  // The convertion is done using the cartesian parametrization as
			  // an intermediate one. If the two parametrizations are equals
			  // there is no need to perform the conversion.
			  if(intype == outtype)
			  {
				  int dim = dimension(outtype);
				  for(int i=0; i<dim; ++i) { outvec[i] = invec[i]; }
			  }
			  // If the input parametrization is the CARTESIAN param, then 
			  // there is no need to transform the input data.
			  else if(intype == params::CARTESIAN)
			  {
				  from_cartesian(invec, outtype, outvec);
			  }
			  // If the output parametrization is the CARTESIAN param, then
			  // there is no need to convert back to another param.
			  else if(outtype == params::CARTESIAN)
			  {
				  to_cartesian(invec, intype, outvec);
			  }
			  else
			  {
				  // temporary CARTESIAN vector
				  double  temvec[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

				  to_cartesian(invec, intype, temvec);
				  from_cartesian(temvec, outtype, outvec);
			  }
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

		  //! \brief Is the value stored weighted by a cosine factor
		  static bool is_cosine_weighted(params::output t)
		  {
			  switch(t)
			  {
				  case params::INV_STERADIAN_COSINE_FACTOR:
					  return true;
					  break;

				  case params::INV_STERADIAN:
				  case params::ENERGY:
				  case params::RGB_COLOR:
				  case params::XYZ_COLOR:
				  default:
					  return false;
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

			  const double temp = cost * vec[0] + sint * vec[2];

			  vec[2] = cost * vec[2] - sint * vec[0];
			  vec[0] = temp;
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
		parametrized(params::input in_param, params::output out_param) {
			_in_param  = in_param;
			_out_param = out_param;
			_nX = params::dimension(_in_param);
			_nY = params::dimension(_out_param);
		}
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
				std::cout << "<<ERROR>> changing to: " << params::get_name(new_param) << std::endl;
				_in_param = new_param;
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


		/* DIMENSION OF THE INPUT AND OUTPUT DOMAIN */

		//! Provide the dimension of the input space of the function
		virtual int dimX() const { return _nX ; }
		//! Provide the dimension of the output space of the function
		virtual int dimY() const { return _nY ; }

		//! Set the dimension of the input space of the function
		virtual void setDimX(int nX) { _nX = nX ; }
		//! Set the dimension of the output space of the function
		virtual void setDimY(int nY) { _nY = nY ; }


		/* DEFINITION DOMAIN OF THE FUNCTION */

		//! \brief Set the minimum value the input can take
		virtual void setMin(const vec& min) ;

		//! \brief Set the maximum value the input can take
		virtual void setMax(const vec& max) ;

		//! \brief Get the minimum value the input can take
		virtual vec min() const ;

		//! \brief Get the maximum value the input can take
		virtual vec max() const ;


	protected:
		// Input and output parametrization
		params::input  _in_param ;
		params::output _out_param ;

		// Dimension of the function & domain of definition.
		int _nX, _nY ;
		vec _min, _max ;
};
