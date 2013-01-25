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

class rational_1d : public function
{
	public: // methods

		rational_1d() ;
		rational_1d(const std::vector<double>& a, const std::vector<double>& b) ;
		virtual ~rational_1d() ;

		// Overload the function operator
		virtual double operator()(double x) const ;

		// Get the p_i and q_j function
		virtual double p(double x, int i) const ;
		virtual double q(double x, int j) const ;

		// IO function to text files
		void load(const std::string& filename) ;
		void save() const ;

		// STL stream ouput
		friend std::ostream& operator<< (std::ostream& out, const rational_1d& r) ;

	private: // data

		// Store the coefficients for the moment, I assume
		// the functions to be polynomials.
		std::vector<double> a ;
		std::vector<double> b ;
} ;

class rational_1d_data : public data
{
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

class rational_1d_fitter : public fitter
{
	public: // methods

		// Fitting a data object
		virtual bool fit_data(const data* d, function*& f) ;

		// Fitting a data object using np elements in the numerator and nq 
		// elements in the denominator
		virtual bool fit_data(const data* d, int np, int nq, function*& f) = 0 ;

		virtual void set_parameters(const arguments& args) ;

	protected: // data

		// min and Max usable np and nq values for the fitting
		int _max_np, _max_nq ;
		int _min_np, _min_nq ;
} ;

