#pragma once

// Include STL
#include <functional>
#include <vector>
#include <string>
#include <tuple>

// Interface
#include <QObject>
#include <core/function.h>
#include <core/data.h>
#include <core/fitter.h>
#include <core/args.h>

class rational_data : public QObject, public data
{
	Q_OBJECT
	Q_INTERFACES(data)

	public: // methods

		// Load data from a file
		virtual void load(const std::string& filename) ;
		virtual void load(const std::string& filename, double min, double max) ;

		// Acces to data
		virtual bool get(int i, double& x, double& yl, double& yu) const ;
		virtual const std::vector<double>& get(int i) const ;		
		virtual const std::vector<double>& operator[](int i) const ;

		// Get data size
		virtual int size() const ;

		// Get min and max input parameters
		virtual double min() const ;
		virtual double max() const ; 
	
		// Get the dimension of the input and output space
		// WARNING: this dimension is defined after loading
		// the data!
		virtual int input_dimension() const ;
		virtual int output_dimension() const ;

	private: // data

		// Store for each point of data, the upper
		// and lower value
		std::vector<std::vector<double> > _data ;

		// Store the min and max value on the input
		// domain
		double _min, _max ;
} ;

