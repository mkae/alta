#pragma once

#include <functional>
#include <string>

#include <QtPlugin>

#include "common.h"
#include "args.h"

class data ;

class function 
{
	public: // methods

		// Overload the function operator
		virtual vec operator()(const vec& x) const = 0 ;
		virtual vec value(const vec& x) const = 0 ;

		// IO function to text files
		virtual void load(const std::string& filename) = 0 ;
		virtual void save(const std::string& filename, const arguments& args) const = 0 ;

		virtual int dimX() const { return _nX ; }
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
} ;

//Q_DECLARE_INTERFACE(function, "Fitter.Function") 
