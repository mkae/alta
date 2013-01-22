#pragma once

#include <string>
#include <pair>

template<class X, class Y> class data
{
		public: // methods

		// Load data from a file
		virtual void load(const std::string& filename) = 0 ;
		virtual void load(const std::string& filename, X min, X max) = 0 ;

		// Acces to data
		virtual bool get(int i, X& x, Y& y) const = 0 ;
		virtual const std::pair<X, Y>& operator[](int i) const = 0 ;

		// Get data size
		virtual int size() const = 0 ;

		// Get min and max input parameters
		virtual X min() const = 0 ;
		virtual X max() const = 0 ;
} ;
