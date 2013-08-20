#pragma once

// Include STL
#include <vector>
#include <string>
#include <sstream>

// Interface
#include "function.h"
#include "data.h"
#include "fitter.h"
#include "args.h"
#include "common.h"


class rational_function_1d : public function
{
	public: // methods

		rational_function_1d() ;
		rational_function_1d(int np, int nq) ;
		rational_function_1d(const std::vector<double>& a,
		                     const std::vector<double>& b) ;
		virtual ~rational_function_1d() {}

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

		// Resize the polynomial
		virtual void resize(int np, int nq);

		// Get the coefficients
		virtual double getP(int i) const { return a[i]; }
		virtual double getQ(int i) const { return b[i]; }

		virtual std::vector<double> getP() const { return a; }
		virtual std::vector<double> getQ() const { return b; }

		// STL stream ouput
		friend std::ostream& operator<< (std::ostream& out,
		                                 const rational_function_1d& r) ;


		//! Convert a 1D index into a vector of degree for a
		//! multinomial coeffcient. The resulting vector v should
		//! be used as prod_k x[k]^v[k] for the monomial basis
		std::vector<int> index2degree(int i) const ;

	protected: // functions

		static int estimate_dk(int k, int d);
		static void populate(std::vector<int>& vec, int N, int M, int j);


	protected: // data

		// Store the coefficients for the moment, I assume
		// the functions to be polynomials.
		std::vector<double> a ;
		std::vector<double> b ;
} ;

class rational_function : public function
{
	public: // methods

		rational_function() ;
		rational_function(int np, int nq) ;
		virtual ~rational_function() ;

		// Overload the function operator
		virtual vec value(const vec& x) const ;
		virtual vec operator()(const vec& x) const { return value(x) ; }

		// IO function to text files
		virtual void load(const std::string& filename) ;

		// Update the function
		virtual void update(int i, rational_function_1d* r) ;

		//! Get the 1D function associated with color channel i. If no one exist, this
		//! function allocates a new element. If i > nY, it returns NULL.
		virtual rational_function_1d* get(int i) ;
		virtual rational_function_1d* get(int i) const ;

		// STL stream ouput
		friend std::ostream& operator<<(std::ostream& out, rational_function& r) ;

		//! \brief Output the rational function as a gnuplot file. It requires
		//! the data object to output the function at the input location only.
		virtual void save_gnuplot(const std::string& filename, 
		                          const data* d, const arguments& args) const ;

		//! Set the dimension of the output space of the function. This function will update
		//! the size of the rs vector size.
		virtual void setDimY(int nY) 
		{ 
			_nY = nY ;
			rs.resize(nY);
		}

	protected: // functions

		//! \brief Save the rational function to the rational format (see \ref formating).
		virtual void save(const std::string& filename) const ;

		//! \brief Output the rational function using a C++ function formating.
		virtual void save_cpp(const std::string& filename, const arguments& args) const ;

		//! \brief Output the rational function using a C++ function formating.
		virtual void save_matlab(const std::string& filename, const arguments& args) const ;

	protected: // data

		// Store the y \in R rational functions. Each channel is a distinct polynomial
		// and should be fitted separately.
		std::vector<rational_function_1d*> rs ;

		// Size of the polynomials
		//! \todo Change it by a more adaptive scheme, with different np, nq per color 
		//!channel?
		int np, nq;
} ;
