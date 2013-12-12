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


class schlick : public nonlinear_function//fresnel
{

	public: // methods

		schlick()
		{
			setParametrization(params::CARTESIAN);
			setDimX(6);
		}

		//! \brief Load function specific files
		virtual bool load(std::istream& in) ;

		virtual void save_call(std::ostream& out, const arguments& args) const;
		virtual void save_body(std::ostream& out, const arguments& args) const;

	protected: // methods

		virtual vec operator()(const vec& x) const { return value(x); }
		virtual vec value(const vec& x) const;

		//! \brief Number of parameters to this non-linear function
		virtual int nbParameters() const ;

		//! \brief Get the vector of parameters for the function
		virtual vec parameters() const ;

		//! \brief Update the vector of parameters for the function
		virtual void setParameters(const vec& p) ;

		//! Get the vector of min parameters for the function
		virtual vec getParametersMin() const;

		//! \brief Obtain the derivatives of the function with respect to the
		//! parameters. 
		virtual vec parametersJacobian(const vec& x) const ;

		//! \brief Boostrap the function by defining the diffuse term
		virtual void bootstrap(const data* d, const arguments& args);

		//! \brief resize the parameter vector
		virtual void setDimY(int nY)
		{
			function::setDimY(nY);
			w.resize(nY);
		}

	private: // data

		//! Fresnel reflectance at theta = 0
		vec w;
} ;

