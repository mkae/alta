#pragma once

// Include STL
#include <functional>
#include <vector>
#include <string>
#include <tuple>

// Personal include
//#include <fitting_algorithm.h>

class rational_1d : public std::function<float(float)> 
{
	public: // methods

		rational_1d() ;
		rational_1d(const std::vector<float>& a, const std::vector<float>& b) ;
		virtual ~rational_1d() ;

		// Overload the function operator
		virtual float operator()(float x) const ;

		// Get the p_i and q_j function
		virtual float p(float x, int i) const ;
		virtual float q(float x, int j) const ;

		// IO function to text files
		void load(const std::string& filename) ;
		void save() const ;

		// STL stream ouput
		friend std::ostream& operator<< (std::ostream& out, const rational_1d& r) ;

	private: // data

		// Store the coefficients for the moment, I assume
		// the functions to be polynomials.
		std::vector<float> a ;
		std::vector<float> b ;
} ;

class rational_1d_data // : public fitting_data
{
	public: // methods

		// Load data from a file
		void load(const std::string& filename) ;

		// Acces to data
		bool get(int i, float& x, float& yl, float& yu) const ;
		const std::vector<float>& operator[](int i) const ;

		// Get data size
		int size() const ;

		// Get min and max input parameters
		float min() const ;
		float max() const ; 

	private: // data

		// Store for each point of data, the upper
		// and lower value
		std::vector<std::vector<float> > _data ;

		// Store the min and max value on the input
		// domain
		float _min, _max ;
} ;

class rational_1d_fitter // : public fitting_algorithm
{
	public: // methods

		// Fitting a data object
		rational_1d fit_data(const rational_1d_data& data) ;

		// Fitting a data object using np elements
		// in the numerator and nq elements in the
		// denominator
		rational_1d fit_data(const rational_1d_data& data, int np, int nq) ;
} ;

