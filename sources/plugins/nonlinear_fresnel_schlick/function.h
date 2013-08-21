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


class schlick : public nonlinear_function
{

	public: // methods

		// Overload the function operator
		virtual vec operator()(const vec& x) const ;
		virtual vec value(const vec& x) const ;


		//! \brief Boostrap the function by defining the diffuse term
		virtual void bootstrap(const data* d, const arguments& args);

		//! \brief Load function specific files
		virtual void load(const std::string& filename) ;

		//! \brief Number of parameters to this non-linear function
		virtual int nbParameters() const ;

		//! \brief Get the vector of parameters for the function
		virtual vec parameters() const ;

		//! \brief Update the vector of parameters for the function
		virtual void setParameters(const vec& p) ;

		//! \brief Obtain the derivatives of the function with respect to the
		//! parameters. 
		virtual vec parametersJacobian(const vec& x) const ;

		//! \brief Provide the dimension of the input space of the function
		inline virtual int dimX() const
		{
			return 1;
		}
		inline virtual int dimY() const
		{
			return 1;
		}

		//! \brief Provide the parametrization of the input space of the function.
		//! For this one, we fix that the parametrization is in THETAD_PHID
		virtual params::input parametrization() const;
		virtual void setParametrization(params::input new_param);

	protected: // methods

	private: // data

		//! Unidimensional Fresnel reflectance at theta = 0
		double R;
} ;

