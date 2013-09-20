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

/*! \brief A spherical Gaussian lobe class.
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

		ward_function()
		{ 
			setParametrization(params::CARTESIAN);
			setDimX(6);
		}

		// Overload the function operator
		virtual vec operator()(const vec& x) const ;
		virtual vec value(const vec& x) const ;

		//! \brief Load function specific files
		virtual void load(std::istream& in) ;

		//! \brief Export function
		virtual void save_call(std::ostream& out, const arguments& args) const;
      virtual void save_body(std::ostream& out, const arguments& args) const;

		//! \brief Boostrap the function by defining the diffuse term
		//!
		//! \details
		virtual void bootstrap(const data* d, const arguments& args);

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

	private: // data

		vec _ks, _ax, _ay; // Lobes data
} ;

