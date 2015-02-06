/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

// Include STL
#include <vector>
#include <string>

// Interface
#include <core/function.h>
#include <core/rational_function.h>
#include <core/data.h>
#include <core/fitter.h>
#include <core/args.h>
#include <core/common.h>

/*! \brief A isotropic_lafortune lobe class. It is provided for testing with the nonlinear
 *  fitting algorithms.
 *
 *  \details
 *  A isotropic_lafortune lobe is defined as \f$k_d + (L^T M V)^n\f$. We fit the restricted
 *  version where the M matrix is diagonal of coefficients \f$(Cx, Cy, Cz)\f$
 *
 *  Options: the <em>bootstrap</em> function will read the arguments to set the first value 
 *  of the function. By default, each lobe is defined as a forward lobe [-1,-1,1] using the 
 *  number of the lobe as the exponent.
 *
 *  Another initialization method <em>--booststrap</em> will put either a forward lobe, 
 *  a retro-reflective lobe [1,1,1] or the dot product [0,0,1]. The exponent will also 
 *  be the number of the lobe.
 */
class isotropic_lafortune_function : public nonlinear_function
{

	public: // methods

		isotropic_lafortune_function() : _n(1) 
		{ 
			setParametrization(params::CARTESIAN);
			setDimX(6);
		}

		// Overload the function operator
		virtual vec operator()(const vec& x) const ;
		virtual vec value(const vec& x) const ;
		virtual vec value(const vec& x, const vec& p) const;

		//! \brief Load function specific files
        virtual bool load(std::istream& in) ;

		//! \brief Export function
		virtual void save_call(std::ostream& out, const arguments& args) const;
      virtual void save_body(std::ostream& out, const arguments& args) const;

		//! \brief Boostrap the function by defining the diffuse term
		//!
		//! \details
		virtual void bootstrap(const ptr<data> d, const arguments& args);

		//! \brief Number of parameters to this non-linear function
		virtual int nbParameters() const ;

		//! \brief Get the vector of parameters for the function
		virtual vec parameters() const ;

		//! \brief get the min values for the parameters
		virtual vec getParametersMin() const;

		//! \brief Update the vector of parameters for the function
		virtual void setParameters(const vec& p) ;

		//! \brief Obtain the derivatives of the function with respect to the
		//! parameters. 
		virtual vec parametersJacobian(const vec& x) const ;

		//! \brief Provide the dimension of the input space of the function
		virtual int dimX() const
		{
			return 6;
		}
		
		//! \brief Set the number of output dimensions
		void setDimY(int nY);

		//! \brief Set the number of lobes to be used in the fit
		void setNbLobes(int N);

	private: // methods

		//! \brief Provide the coefficient of the monochromatic lobe number
		//! n for the color channel number c.
		inline void getCurrentLobe(int n, int c, double& Cx, double& Cz, double& N) const 
		{
			Cx = _C[(n*_nY + c)*2 + 0];
			Cz = _C[(n*_nY + c)*2 + 1];
			N  = _N[n*_nY + c];
		}


	private: // data

		//! \brief The isotropic_lafortune lobe data
		int _n; // Number of lobes
		vec _N, _C; // Lobes data
#ifdef WITH_DIFFUSE
		vec _kd; // Diffuse term
#endif
} ;

