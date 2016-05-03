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
#include <core/rational_function.h>
#include <core/data.h>
#include <core/fitter.h>
#include <core/args.h>
#include <core/common.h>

using namespace alta;


/*! \brief A retroblinn lobe class. It is provided for testing with the 
 *  nonlinear fitting algorithms.
 *
 *  \details
 *  A retroblinn lobe is defined as \f$k_d + k_s |L.V|^a\f$
 */
class retroblinn_function : public nonlinear_function
{

    public: // methods

		 // Overload the function operator
		 virtual vec operator()(const vec& x) const ;
		 virtual vec value(const vec& x) const ;

		 //! \brief Boostrap the function. This initialize the lobe width
		 //! and the specular coefficient to one.
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
       _parameters = alta::parameters(_parameters.dimX(), nY,
                                      _parameters.input_parametrization(),
                                      _parameters.output_parametrization());

			 // Update the length of the vectors
			 _ks.resize(_parameters.dimY()) ;
			 _N.resize(_parameters.dimY()) ;
		 }

 		 //! \brief Export function
		virtual void save_call(std::ostream& out, 
		                       const arguments& args) const;
      virtual void save_body(std::ostream& out, 
		                       const arguments& args) const;

 private: // data

		 //! \brief The retroblinn lobe data
		 vec _ks, _N;
} ;

