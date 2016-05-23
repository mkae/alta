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
 *  \class nonlinear_function_retrobeckmann
 *  \ingroup functions
 *  \ingroup plugins
 *  \brief An empirical retro-reflecting BRDF model based on a Beckmann distribution 
 *  (see [Belcour et al. [2014]](https://hal.inria.fr/hal-01083366)).
 *
 *  \details
 *  This BRDF model uses the Beckmann distribution with the back vector 
 *  \f$ K = \frac{L - V}{||L - V||} \f$ in place of the standard Half vector.
 *
 *  \author Laurent Belcour \<laurent.belcour@umontreal.ca\>
 */
class beckmann_function : public nonlinear_function
{

	public: // methods

		beckmann_function():
        nonlinear_function(alta::parameters(6, 0,
                                            params::CARTESIAN,
                                            params::UNKNOWN_OUTPUT)),
       _use_back_param(true)
    {}

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

		vec _ks, _a; // Lobes data
		bool _use_back_param;
} ;

