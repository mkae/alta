/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2016 Inria

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
#include <core/data.h>
#include <core/fitter.h>
#include <core/args.h>
#include <core/common.h>

using namespace alta;

/*! 
 *	 \class nonlinear_function_ward
 *  \ingroup functions
 *  \ingroup plugins
 *  \brief Ward's anisotropic BRDF model.
 *
 *  \details
 *  A ward lobe is defined as \f${k_s \over 4 \pi \alpha_x \alpha_y 
 *  \sqrt((i.n) (o.n))} exp(- ({(h.x) / \alpha_x)^2 + ((h.y) / 
 *  \alpha_y)^2 \over (h.n)^2}\f$. Where \f$(\alpha_x, alpha_y)\f$ is 
 *  the anisotropic roughness control.
 *
 *  <h3>Plugin parameters</h3>
 *
 *  <ul>
 *     <li><b>bootstrap</b> function will set the lobe values as 1/li>
 *  </ul>
 */
class ward_function : public nonlinear_function
{

	public: // methods

    ward_function(const alta::parameters& params);

		// Overload the function operator
		virtual vec operator()(const vec& x) const ;
		virtual vec value(const vec& x) const ;

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

		//! \brief Set the number of output dimensions
		void setDimY(int nY);

	private: // data

		vec _ks, _ax, _ay; // Lobes data

        //! Allows to set the lobe to be isotropic or not
        bool isotropic;

    ward_function() {};
} ;

