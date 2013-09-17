#pragma once

// Include STL
#include <vector>
#include <string>

// Interface
#include <core/function.h>
#include <core/data.h>
#include <core/args.h>
#include <core/common.h>


/*! \brief A blinn lobe class. It is provided for testing with the nonlinear
 *  fitting algorithms.
 *
 *  \details
 *  A blinn lobe is defined as \f$k_s |N.H|^a\f$
 *  \todo Finish implementation
 */
class blinn_function : public nonlinear_function
{

	public: // methods

		// Overload the function operator
		virtual vec operator()(const vec& x) const ;
		virtual vec value(const vec& x) const ;


		//! \brief Boostrap the function by defining the diffuse term
		virtual void bootstrap(const data* d, const arguments& args);

		//! \brief Load function specific files
		virtual void load(std::istream& filename) ;

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
		virtual int dimX() const
		{
			return 1 ;
		}

		//! \brief Provide the parametrization of the input space of the 
		//! function.
		virtual params::input parametrization() const
		{
			return params::COS_TH ;
		}
		virtual void setParametrization(params::input new_param)
		{
			std::cerr << "Cannot change the ouput parametrization " 
				       << __FILE__ << ":" << __LINE__ << std::endl;
			throw;
		}

		void setDimY(int nY)
		{
			_nY = nY ;

			// Update the length of the vectors
			_ks.resize(_nY) ;
			_N.resize(_nY) ;
		}

		void save_call(std::ostream& out, const arguments& args) const;
		void save_body(std::ostream& out, const arguments& args) const;

	private: // data

		//! \brief The blinn lobe data
		vec _ks, _N;
} ;

