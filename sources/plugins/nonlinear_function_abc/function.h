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

/*!
 * \ingroup functions
 *	\brief  The ABC distribution for micro-facets models.
 *  \details
 *  Follows the implementation of Löw <i>et al.</i> [2012].
 *
 *  <h3>Plugin parameters</h3>
 *
 */
class abc_function : public nonlinear_function
{

	public: // methods

		abc_function()
		{
			setParametrization(params::COS_TH);
			setDimX(1);
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

		vec _a, _b, _c; // Lobes data
} ;

