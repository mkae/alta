/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2015 CNRS
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
#include <core/args.h>
#include <core/common.h>

using namespace alta;


/*! 
 *  \class nonlinear_function_blinn
 *  \ingroup functions
 *	 \brief A Blinn-Phong lobe.
 *
 *  \details
 *  A Blinn-Phong lobe has two parameters and is defined as:
 *  <center>
 *  \f$\rho(L, V)k_s |N.H|^a\f$, where \f$ H = \frac{L+V}{|| L+V ||}\f$.
 *  </center>
 *
 *  \author Laurent Belcour \<laurent.belcour@umontreal.ca\>
 */
class blinn_function : public nonlinear_function
{

	public: // methods

		blinn_function()
		{
        _parameters = alta::parameters(1, 1, params::COS_TH,
                                       params::UNKNOWN_OUTPUT);
		}

		// Overload the function operator
		virtual vec operator()(const vec& x) const ;
		virtual vec value(const vec& x) const ;


		//! \brief Boostrap the function by defining the diffuse term
		virtual void bootstrap(const ptr<data> d, const arguments& args);

		//! \brief Load function specific files
        virtual bool load(std::istream& filename) ;

		//! \brief Number of parameters to this non-linear function
		virtual int nbParameters() const ;

		//! \brief Get the vector of parameters for the function
		virtual vec parameters() const ;
		
		//! \brief The minimum parameter vector for this BRDF is the zero
		//! vector. The specular intensity cannot be negative and the
		//! exponent should not be either.
		virtual vec getParametersMin() const
		{
			return vec::Zero(_parameters.dimY()*2);
		}

		//! \brief Update the vector of parameters for the function
		virtual void setParameters(const vec& p) ;

		//! \brief Obtain the derivatives of the function with respect to the
		//! parameters. 
		virtual vec parametersJacobian(const vec& x) const ;

		void setDimY(int nY)
		{
			//CODE CONSISTENCY WITH THE OTHER PLUGIN
			//_nY = nY ;
			function::setDimY(nY);

			// Update the length of the vectors
			_ks.resize(_parameters.dimY()) ;
			_N.resize(_parameters.dimY()) ;
		}

		void save_call(std::ostream& out, const arguments& args) const;
		void save_body(std::ostream& out, const arguments& args) const;

	private: // lobe parameters

		//! \brief The blinn lobe parameters
		vec _ks, _N;
} ;

