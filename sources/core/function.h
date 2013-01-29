#pragma once

#include <functional>
#include <string>

#include <QtPlugin>

#include "common.h"
#include "args.h"

class function 
{
	public: // methods

		// Overload the function operator
		virtual vec operator()(const vec& x) const = 0 ;
		virtual vec value(const vec& x) const = 0 ;
		
		// IO function to text files
		virtual void load(const std::string& filename) = 0 ;
		virtual void save(const std::string& filename, const arguments& args) const = 0 ;

		virtual int dimX() const { return _nX ; } ;
		virtual int dimY() const { return _nY ; } ;

		virtual void setDimX(int nX) { _nX = nX ; } ;
		virtual void setDimY(int nY) { _nY = nY ; } ;

	protected:
		// Dimension of the function
		int _nX, _nY ;
} ;

Q_DECLARE_INTERFACE(function, "Fitter.Function") 
