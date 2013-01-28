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
#include <core/common.h>

class rational_function : public QObject, public function
{
	Q_OBJECT
	Q_INTERFACES(function)

	public: // methods

		rational_function() ;
		rational_function(const std::vector<double>& a, const std::vector<double>& b) ;
		virtual ~rational_function() ;

		// Overload the function operator
		virtual vec operator()(const vec& x) const ;

		// Get the p_i and q_j function
		virtual double p(const vec& x, int i) const ;
		virtual double q(const vec& x, int j) const ;

		// IO function to text files
		void load(const std::string& filename) ;
		void save() const ;

		// STL stream ouput
		friend std::ostream& operator<< (std::ostream& out, const rational_function& r) ;

	private: // data

		// Store the coefficients for the moment, I assume
		// the functions to be polynomials.
		std::vector<double> a ;
		std::vector<double> b ;
} ;

