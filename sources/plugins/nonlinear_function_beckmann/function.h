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

/*! \brief A Microfacet model using a Gaussian distribution for the
 *  micro-gemeotry distribution parametrized by the half-angle.
 *
 *  \details
 *
 *  <h3>Plugin parameters</h3>
 *
 */
class beckmann_function : public nonlinear_function
{

	public: // methods

		beckmann_function()
		{
			setParametrization(params::CARTESIAN);
			setDimX(6);
		}

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

		vec _ks, _a; // Lobes data
} ;

