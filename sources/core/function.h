#pragma once

#include <functional>
#include <string>

#include <QtPlugin>

#include "common.h"
#include "args.h"

class data ;

/*! \brief A representation of an analytical function.
 *
 */
class function 
{
	public: // methods

		// Overload the function operator
		virtual vec operator()(const vec& x) const = 0 ;
		virtual vec value(const vec& x) const = 0 ;

		// IO function to text files
		virtual void load(const std::string& filename) = 0 ;
		virtual void save(const std::string& filename, const arguments& args) const = 0 ;

		/*! \brief Provide the dimension of the input space of the function
		 */
		virtual int dimX() const { return _nX ; }
		/*! \brief Provide the dimension of the output space of the function
		 */
		virtual int dimY() const { return _nY ; }

		virtual void setDimX(int nX) { _nX = nX ; }
		virtual void setDimY(int nY) { _nY = nY ; }

		// Acces to the domain of definition of the function
		virtual void setMin(const vec& min) 
		{
			assert(min.size() == _nX) ;
		  	_min = min ; 
		}
		virtual void setMax(const vec& max) 
		{
			assert(max.size() == _nX) ;
			_max = max ; 
		}
		virtual vec getMin() const { return _min ; }
		virtual vec getMax() const { return _max ; }


	protected: //data

		// Dimension of the function & domain of
		// definition.
		int _nX, _nY ;
		vec _min, _max ;
};

class nonlinear_function: public function
{
	public: // methods

		// Set the vector of parameters for the function
		virtual vec parameters() const = 0;
		virtual void setParameters(const vec& p) = 0;

		// Obtain the derivatives of the function with respect
		// to the parameters. The x input of this function is 
		// the position in the input space and has size dimX(),
		// the resulting vector has the size of the parameters:
		// [df/dp1, ..., df/dpn]
		virtual vec parameters_derivatives(const vec& x) const = 0;
};

Q_DECLARE_INTERFACE(function, "Fitter.Function") 
