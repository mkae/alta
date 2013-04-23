#pragma once

// Include STL
#include <vector>
#include <string>

// Interface
#include <QObject>
#include "function.h"
#include "data.h"
#include "fitter.h"
#include "args.h"
#include "common.h"

class rational_function : public QObject, public function
{
	Q_OBJECT
	Q_INTERFACES(function)

	public: // methods

		rational_function() ;
        rational_function(int np, int nq) ;
		rational_function(const std::vector<double>& a, const std::vector<double>& b) ;
		virtual ~rational_function() ;

		// Overload the function operator
		virtual vec value(const vec& x) const ;
		virtual vec operator()(const vec& x) const { return value(x) ; }
		
		// Get the numerator (p) and denominator (q) functions
		virtual vec p(const vec& x) const ;
		virtual vec q(const vec& x) const ;

		// Get the p_i and q_j function
		virtual double p(const vec& x, int i) const ;
		virtual double q(const vec& x, int j) const ;

		// IO function to text files
		virtual void load(const std::string& filename) ;
	
		// Update the function
		virtual void update(const std::vector<double>& in_a, 
		                    const std::vector<double>& in_b) ;

		// Get the coefficients
		virtual double getP(int i) const { return a[i] ; }
		virtual double getQ(int i) const { return b[i] ; }

		// STL stream ouput
		friend std::ostream& operator<< (std::ostream& out, const rational_function& r) ;

        static int estimate_dk(int k, int d);
        static void populate(std::vector<int>& vec, int N, int M, int j);

	protected: // functions
		
        //! Convert a 1D index into a vector of degree for a
        //! multinomial coeffcient. The resulting vector v should
        //! be used as prod_k x[k]^v[k] for the monomial basis
		std::vector<int> index2degree(int i) const ;
		
		//! \brief Save the rational function to the rational format (see \ref formating).
		virtual void save(const std::string& filename) const ;

		//! \brief Output the rational function as a gnuplot file. It requires
		//! the data object to output the function at the input location only.
		virtual void save_gnuplot(const std::string& filename, const data* d, const arguments& args) const ;

		//! \brief Output the rational function using a C++ function formating.
		virtual void save_cpp(const std::string& filename, const arguments& args) const ;

		//! \brief Output the rational function using a C++ function formating.
		virtual void save_matlab(const std::string& filename, const arguments& args) const ;

	protected: // data

		// Store the coefficients for the moment, I assume
		// the functions to be polynomials.
		std::vector<double> a ;
		std::vector<double> b ;
} ;
