#pragma once

#include <string>
#include <utility>

#include <QtPlugin>

class data
{
		public: // methods

		// Load data from a file
		virtual void load(const std::string& filename) = 0 ;
		virtual void load(const std::string& filename, double min, double max) = 0 ;

		// Acces to data
		virtual bool get(int i, double& x, double& y, double& t) const = 0 ;
		virtual const std::vector<double>& get(int i) const = 0 ;
		virtual const std::vector<double>& operator[](int i) const = 0 ;

		// Get data size, e.g. the number of samples to fit
		virtual int size() const = 0 ;

		// Get the dimension of the input and output space
		// WARNING: this dimension is defined after loading
		// the data!
		virtual int input_dimension() const = 0;
		virtual int output_dimension() const = 0;

		// Get min and max input space values
		virtual double min() const = 0 ;
		virtual double max() const = 0 ;
} ;
	 
Q_DECLARE_INTERFACE(data, "Fitter.Data") 
