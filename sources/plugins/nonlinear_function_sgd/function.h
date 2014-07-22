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

class shifted_gamma_function : public nonlinear_function
{
	public: // methods

		// Overload the function operator
		virtual vec operator()(const vec& x) const ;
		virtual vec value(const vec& x) const ;

		//! Load function specific files
		virtual void load(const std::string& filename) ;

		//! Save the current function to a specific file type
		virtual void save(const std::string& filename, const arguments& args) const ;
		
		//! \brief Export function
		virtual void save_call(std::ostream& out, const arguments& args) const {
			NOT_IMPLEMENTED();
		}
		virtual void save_body(std::ostream& out, const arguments& args) const {
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
} ;

