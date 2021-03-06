/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2014, 2016 Inria

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

class shifted_gamma_function : public nonlinear_function
{
	public: // methods

    shifted_gamma_function(const alta::parameters& params);

		// Overload the function operator
		virtual vec operator()(const vec& x) const ;
		virtual vec value(const vec& x) const ;

		//! \brief Export function
		virtual void save_call(std::ostream& out, const arguments& args) const 
    {
			NOT_IMPLEMENTED();
		}
		
    virtual void save_body(std::ostream& out, const arguments& args) const 
    {
			NOT_IMPLEMENTED();
		}
		
		//! Number of parameters to this non-linear function
		virtual int nbParameters() const ;

		//! Get the vector of parameters for the function
		virtual vec parameters() const ;

		//! Update the vector of parameters for the function
		virtual void setParameters(const vec& p) ;

		//! Obtain the derivatives of the function with respect to the 
		//! parameters. 
		virtual vec parametersJacobian(const vec& x) const ;

    //! Load BRDF parameters from .brdf file
    virtual bool load(std::istream& in);

    //! Save the Function 
    virtual void save(const std::string& filename, const arguments& args) const;

    friend std::ostream& operator<<(std::ostream& out, const shifted_gamma_function & sgd);

	private:
		//! Fresnel term of the microfacet distribution
		vec Fresnel(const vec& F0, const vec& F1, double V_H) const;

		//! Shifted gamma distribution function
		vec D(const vec& _alpha, const vec& _p, double cos_h, const vec& _K) const;

		//! Visibility function
		vec G1(double theta) const;

		//! Parameters of the visibility function
		vec sh_c, sh_theta0, sh_k, sh_lambda, p;

		//! Parameters of the Fresnel function
		vec F_0, F_1;

		//! Parameters of the microfacets function
		vec K_ap;

		//! Color parameters
		vec rho_d, rho_s, alpha;

    shifted_gamma_function() {};
} ;


std::ostream& operator<<(std::ostream& out, const shifted_gamma_function & sgd)
{
  out << " Shifted Gamma Function: Parameters = { rho _d = " << sgd.rho_d
      << " rho_s = " << sgd.rho_s << " alpha = " << sgd.alpha << std::endl
      << " K_alpha_p = " << sgd.K_ap << " F_0 = " << sgd.F_0 << " F_1 = " << sgd.F_1 << std::endl
      << " sh_c = " << sgd.sh_c << " sh_thetao = " << sgd.sh_theta0 << std::endl
      << " sh_k = " << sgd.sh_k << " sh_lambda = " << sgd.sh_lambda << " power = " << sgd.p 
      << " } " << std::endl; 
  
  return out;
}


