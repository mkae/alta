#pragma once

// Include STL
#include <vector>
#include <string>

// Interface
#include <QObject>
#include <core/function.h>
#include <core/rational_function.h>
#include <core/data.h>
#include <core/fitter.h>
#include <core/args.h>
#include <core/common.h>

//#define ADAPT_TO_PARAM
//#define FIT_DIFFUSE

/*! \brief A lafortune lobe class. It is provided for testing with the nonlinear
 *  fitting algorithms.
 *
 *  \details
 *  A lafortune lobe is defined as \f$k_d + (L^T M V)^n\f$. We fit the restricted
 *  version where the M matrix is diagonal of coefficients \f$(Cx, Cy, Cz)\f$
 *  \todo Fitting the diffuse part is not stable
 */
class lafortune_function : public nonlinear_function, public QObject
{

	//    Q_OBJECT
	Q_INTERFACES(function)

	public: // methods

        lafortune_function() : _n(1) { }

	// Overload the function operator
		virtual vec operator()(const vec& x) const ;
		virtual vec value(const vec& x) const ;
		virtual vec value(const vec& x, const vec& p) const;

		//! \brief Load function specific files
		virtual void load(const std::string& filename) ;

		//! \brief Boostrap the function by defining the diffuse term
		virtual void bootstrap(const data* d, const arguments& args);

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
#ifdef ADAPT_TO_PARAM
			return _nX;
#else
			return 6;
#endif
		}

		//! \brief Provide the parametrization of the input space of the function.
		//! For this one, we fix that the parametrization is in THETAD_PHID
		virtual params::input parametrization() const
		{
#ifdef ADAPT_TO_PARAM
			return _in_param;
#else
			return params::CARTESIAN ;
#endif
		}
		virtual void setParametrization(params::input new_param)
		{
#ifdef ADAPT_TO_PARAM
			_in_param = new_param;
#else
			std::cerr << "<<ERROR>> Cannot change the ouput parametrization " << __FILE__ << ":" << __LINE__ << std::endl;
			throw;
#endif
		}

        //! \brief Set the number of output dimensions
        void setDimY(int nY);

        //! \brief Set the number of lobes to be used in the fit
        void setNbLobes(int N);

	protected: // methods

		virtual void save(const std::string& filename) const;


		//! \brief Output the function using a BRDF Explorer formating.
		virtual void save_brdfexplorer(const std::string& filename,
				const arguments& args) const;

	private: // methods

		//! \brief Provide the coefficient of the monochromatic lobe number
		//! n for the color channel number c.
		void getCurrentLobe(int n, int c, double& Cx, double& Cy, double& Cz, double& N) const 
		{
			Cx = _C[(n*_nY + c)*3 + 0];
			Cy = _C[(n*_nY + c)*3 + 1];
			Cz = _C[(n*_nY + c)*3 + 2];
			N  = _N[n*_nY + c];
		}


	private: // data

		//! \brief The lafortune lobe data
		int _n; // Number of lobes
		vec _N, _C; // Lobes data
		vec _kd; // Diffuse term
} ;

