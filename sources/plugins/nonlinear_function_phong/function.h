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
 *  A phong lobe is defined as \f$k_d + k_s |N.H|^a\f$
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
        virtual vec value(const vec& x, const vec& p) const
        {
            // Test input parameters for correct size
            assert(p.size() == nbParameters());

            vec res(dimY());
            for(int i=0; i<dimY(); ++i)
            {
                const double kd = p[i*3 + 0];
                const double ks = p[i*3 + 1];
                const double N  = p[i*3 + 2];
                res[i] = kd + ks * std::pow(x[0], N);
            }

            return res;
        }


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
        virtual int dimX() const
        {
            return 1 ;
        }

        //! \brief Provide the parametrization of the input space of the function.
        //! For this one, we fix that the parametrization is in THETAD_PHID
        virtual params::input parametrization() const
        {
            return params::COS_TH ;
        }
        virtual void setParametrization(params::input new_param)
        {
            std::cerr << "Cannot change the ouput parametrization " << __FILE__ << ":" << __LINE__ << std::endl;
            throw;
        }

        void setDimY(int nY)
        {
            _nY = nY ;

            // Update the length of the vectors
            _kd.resize(_nY) ;
            _ks.resize(_nY) ;
            _N.resize(_nY) ;
        }

    protected: // methods

        virtual void save(const std::string& filename) const
        {
            std::ofstream file(filename.c_str(), std::ios_base::trunc);
            file << "#DIM " << _nX << " " << _nY << std::endl ;
            file << "#kd" << std::endl;
            for(int i=0; i<_nY; ++i)
            {
                file << _kd[i] << std::endl;
            }
            file << std::endl;

            file << "#ks" << std::endl;
            for(int i=0; i<_nY; ++i)
            {
                file << _ks[i] << std::endl;
            }
            file << std::endl;

            file << "#N" << std::endl;
            for(int i=0; i<_nY; ++i)
            {
                file << _N[i] << std::endl;
            }
            file << std::endl;

        }

        //! \brief Output the function using a BRDF Explorer formating.
        virtual void save_brdfexplorer(const std::string& filename,
                                       const arguments& args) const;

	private: // data

        //! \brief The phong lobe data
		vec _kd, _ks, _N;
} ;

