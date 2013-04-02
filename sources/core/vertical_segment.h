#pragma once

// Include STL
#include <vector>
#include <string>

// Interface
#include "function.h"
#include "data.h"
#include "fitter.h"
#include "args.h"

class vertical_segment : public data
{
	public: // methods

		// Load data from a file
		virtual void load(const std::string& filename) ;
		virtual void load(const std::string& filename, const arguments& args) ;

		// Acces to data
        virtual void get(int i, vec &x, vec &yl, vec &yu) const ;
		virtual vec get(int i) const ;		
		virtual void get(int i, vec& yl, vec& yu) const ;		
		virtual vec operator[](int i) const ;

		// Get data size
		virtual int size() const ;

		// Get min and max input parameters
		virtual vec min() const ;
		virtual vec max() const ; 

	private: // data

		// Store for each point of data, the upper
		// and lower value
		std::vector<vec> _data ;

		// Store the min and max value on the input
		// domain
		vec _min, _max ;
} ;

