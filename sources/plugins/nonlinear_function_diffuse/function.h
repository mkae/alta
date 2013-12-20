#pragma once

// Include STL
#include <vector>
#include <string>

// Interface
#include <core/function.h>
#include <core/data.h>
#include <core/args.h>
#include <core/common.h>


/*!
 * \ingroup functions
 * \brief A Lambertian component class.
 */
class diffuse_function : public nonlinear_function
{

	public: // methods

        // Set the input parametrization to CARTESIAN to reduce the number
        // of transformations in a compound object.
        diffuse_function()
        {
            setParametrization(params::CARTESIAN);
				setDimX(6);
        }

		// Overload the function operator
		virtual vec operator()(const vec& x) const ;
		virtual vec value(const vec& x) const ;


		//! \brief Boostrap the function by defining the diffuse term
		virtual void bootstrap(const data* d, const arguments& args);

		//! \brief Load function specific files
        virtual bool load(std::istream& in) ;
		virtual void save_call(std::ostream& out, const arguments& args) const;

		//! \brief Number of parameters to this non-linear function
		virtual int nbParameters() const ;

		//! \brief Get the vector of parameters for the function
		virtual vec parameters() const ;

		//! \brief Update the vector of parameters for the function
		virtual void setParameters(const vec& p) ;

		//! \brief Obtain the derivatives of the function with respect to the
		//! parameters. 
		virtual vec parametersJacobian(const vec& x) const ;

		virtual void setDimY(int nY)
		{
			_nY = nY;
			_kd.resize(nY);
		}

	private: // data

		vec _kd;
} ;

