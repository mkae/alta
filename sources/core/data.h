#pragma once

#include <string>
#include <utility>

class data
{
		public: // methods

		// Load data from a file
		virtual void load(const std::string& filename) = 0 ;
		virtual void load(const std::string& filename, double min, double max) = 0 ;

		// Acces to data
		virtual bool get(int i, double& x, double& y, double& t) const = 0 ;
		virtual const std::vector<double>& operator[](int i) const = 0 ;

		// Get data size
		virtual int size() const = 0 ;

		// Get min and max input parameters
		virtual double min() const = 0 ;
		virtual double max() const = 0 ;
} ;
