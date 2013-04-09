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


/*! \brief A phong lobe class. It is provided for testing with the nonlinear 
 *  fitting algorithms.
 *
 *  \details
 *  A phong lobe is defined as \f$k_d + k_s |R.H|^N\f$
 *  \todo Finish implementation
 */
class phong_function : public nonlinear_function, public QObject
{

//    Q_OBJECT
    Q_INTERFACES(function)

	public: // methods

		// Overload the function operator
		virtual vec operator()(const vec& x) const ;
		virtual vec value(const vec& x) const ;

        //! \brief Load function specific files
		virtual void load(const std::string& filename) ;

        //! \brief Save the current function to a specific file type
		virtual void save(const std::string& filename, const arguments& args) const ;
		
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
            return 2 ;
        }

        //! \brief Provide the parametrization of the input space of the function.
        //! For this one, we fix that the parametrization is in THETAD_PHID
        virtual params::type parametrization() const
        {
            return params::ISOTROPIC_TD_PD ;
        }


	private: // data

        //! \brief The phong lobe data
		vec _kd, _ks, _N;
} ;

