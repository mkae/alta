#pragma once

#include <functional>
#include <string>
#include <fstream>

#include <QtPlugin>

#include "common.h"
#include "args.h"
#include "params.h"
#include "data.h"

/*! \brief A representation of an analytical function.
 *  \ingroup core
 *
 *  \details
 *  function are functors with a domain of definition specified by a vector 
 *  interval \f$[\vec{min} .. \vec{max}]\f$ where \f$\vec{min}\f$ and 
 *  \f$\vec{max}\f$ have the size of the input domain.
 *
 *  Any function used by the fitting algorithm should overload publicly this
 *  interface.
 */
class function 
{
	public: // methods

        // Constructor
        function() : _in_param(params::UNKNOWN_INPUT),
                     _out_param(params::UNKNOWN_OUTPUT) { }

		// Overload the function operator
		virtual vec operator()(const vec& x) const = 0 ;
		virtual vec value(const vec& x) const = 0 ;

		//! Load function specific files
		virtual void load(const std::string& filename) = 0 ;

		//! \brief Save the current function to a specific file type, args can 
		//! be used to differenciate the type of export.
		//!
		//! \see rational_function.cpp for an example
		virtual void save(const std::string& filename, const arguments& args) const
		{
			if(args.is_defined("export"))
			{
				if(args["export"].compare("c++") == 0)
				{
					std::cout << "<<INFO>> will export in C++ format" << std::endl;
					save_cpp(filename, args);
				}
				else if(args["export"].compare("matlab") == 0)
				{
					std::cout << "<<INFO>> will export in matlab format" << std::endl;
					save_matlab(filename, args);
				}
				else
				{
					std::cerr << "<<ERROR>> the export format is unknown" << std::endl ;
				}
			}
			else
			{
				save(filename) ;
			}
		}

		//! Provide the dimension of the input space of the function
		virtual int dimX() const { return _nX ; }
		//! Provide the dimension of the output space of the function
		virtual int dimY() const { return _nY ; }

		virtual void setDimX(int nX) { _nX = nX ; }
		virtual void setDimY(int nY) { _nY = nY ; }

		// Acces to the domain of definition of the function
		virtual void setMin(const vec& min) 
		{
#ifdef DEBUG
			assert(min.size() == _nX) ;
#endif
		  	_min = min ; 
		}
		virtual void setMax(const vec& max) 
		{
#ifdef DEBUG
			assert(max.size() == _nX) ;
#endif
			_max = max ; 
		}
		virtual vec getMin() const { return _min ; }
		virtual vec getMax() const { return _max ; }

        //! \brief provide the parametrization of the function.
        //! \note some function type can modify the parametrization to adapt
        //! to the data.
        virtual params::input parametrization() const
        {
            return _in_param;
        }

        //! \brief can set the input parametrization of a non-parametrized
        //! function. Throw an exception if it tries to erase a previously
        //! defined one.
        virtual void setParametrization(params::input new_param)
        {
            if(_in_param != params::UNKNOWN_INPUT)
                throw("A parametrization is already defined");

            _in_param = new_param;
        }
		
	protected: // function

		//! \brief Standard saving function.
		virtual void save(const std::string& filename) const 
		{
			NOT_IMPLEMENTED();
		}

		//! \brief Output the rational function as a gnuplot file. It requires
		//! the data object to output the function at the input location only.
		virtual void save_gnuplot(const std::string& filename, const data* d, 
		                          const arguments& args) const
        {
#ifndef OLD
            std::ofstream file(filename.c_str(), std::ios_base::trunc);
            for(int i=0; i<d->size(); ++i)
            {
                vec v = d->get(i) ;
        //		vec y1 ; y1.assign(d->dimY(), 0.0) ;
        //		for(int k=0; k<d->dimY(); ++k) { y1[k] = v[d->dimX() + k] ; }

                vec y2 = value(v) ;
                for(int u=0; u<d->dimX(); ++u)
                    file << v[u] << "\t" ;

                for(int u=0; u<d->dimY(); ++u)
                    file << y2[u] << "\t" ;

                file << std::endl ;
            }
            file.close();
#else
            NOT_IMPLEMENTED();
#endif
        }
		
		//! \brief Output the rational function using a C++ function formating.
		virtual void save_cpp(const std::string& filename, const arguments& args) const 

		{
			NOT_IMPLEMENTED();
		}
		//! \brief Output the rational function using a C++ function formating.
		virtual void save_matlab(const std::string& filename, const arguments& args) const 
		{
			NOT_IMPLEMENTED();
		}


	protected: // data

		// Dimension of the function & domain of definition.
		int _nX, _nY ;
		vec _min, _max ;

        // Input and output parametrization
        params::input  _in_param ;
        params::output _out_param ;
};

/*! \brief Non-linear function interface
 *  \ingroup core
 *
 * \details
 * Provide a way to obtain the dérivative of the function with respect to its
 * parameters. If the function \f$f(\vec{x})\f$ is defined for a vector of
 * parameters \f$\vec{a}\f$, the resulting vector is \f$df_i = {df \over 
 * da_i}\f$. 
 *
 * \note It is not necessary to have an analytical formulation
 * of the derivative and a numerical evaluation of it can be provided.
 */
class nonlinear_function: public function
{
	public: // methods

		//! Number of parameters to this non-linear function
		virtual int nbParameters() const = 0;

		//! Get the vector of parameters for the function
		virtual vec parameters() const = 0;

		//! Update the vector of parameters for the function
		virtual void setParameters(const vec& p) = 0;

		//! \brief Obtain the derivatives of the function with respect to the 
		//! parameters. 
		//
		// The x input of this function is the position in the input space and 
		// has size dimX(), the resulting vector has the size of the parameters
		// times the size of the output domain.
		//
		// The result vector should be orderer as res[i + dimY()*j], output
		// dimension first, then parameters.
		virtual vec parametersJacobian(const vec& x) const = 0;
};

Q_DECLARE_INTERFACE(function, "Fitter.Function") 
