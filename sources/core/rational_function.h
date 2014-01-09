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


/*! \brief A one dimensional rational function class. A rational function has
 *  the form \f$r(x) = \dfrac{\sum_{i} a_i p_i(x)}{b_i q_i(x)}$\f.
 */
class rational_function_1d : public function
{
	public: // methods

		rational_function_1d() ;
		rational_function_1d(int np, int nq, bool separable = false) ;
        rational_function_1d(const vec& a, const vec& b) ;
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
		virtual bool load(std::istream& in);

		// Update the function
		virtual void update(const vec& in_a, 
		                    const vec& in_b) ;

		// Resize the polynomial
		virtual void resize(int np, int nq);

		// Get the coefficients
		virtual double getP(int i) const { return a[i]; }
		virtual double getQ(int i) const { return b[i]; }

		virtual vec getP() const { return a; }
		virtual vec getQ() const { return b; }

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
        vec a, b ;

		  //! Is the function separable with respect to its input dimensions?
		  //! \todo Make possible to have only part of the dimensions
		  //! separable.
		  bool _separable;
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
		virtual bool load(std::istream& in) ;

		// Update the function
		virtual void update(int i, rational_function_1d* r) ;

		//! Get the 1D function associated with color channel i. If no one exist, 
		//! this function allocates a new element. If i > nY, it returns NULL.
		virtual rational_function_1d* get(int i) ;
		virtual rational_function_1d* get(int i) const ;

		// STL stream ouput
		friend std::ostream& operator<<(std::ostream& out, rational_function& r) ;

		//! Set the dimension of the output space of the function. This function 
		//! will update the size of the rs vector size.
		virtual void setDimY(int nY) 
		{ 
			_nY = nY ;
			rs.resize(nY);
		}

		//! \brief Set the size of the rational function. Any newly created 1D 
		//! function will have np and nq fixed as defined.
		virtual void setSize(int np, int nq)
		{
			this->np = np;
			this->nq = nq;
			clear();
		}

		//! \brief Clear the vector of 1D rational functions.
		virtual void clear()
		{
			rs.clear();
			rs.resize(_nY);
		}

		virtual void setMin(const vec& min)
		{
			function::setMin(min);
		}

		virtual void setMax(const vec& max)
		{
			function::setMax(max);
		}

		//! \brief Save the rational function to the rational format (see 
		//! \ref formating).
		virtual void save_call(std::ostream& out, const arguments& args) const ;

	protected: // data

		//! Store the y \in R rational functions. Each channel is a distinct 
		//! polynomial and should be fitted separately.
		std::vector<rational_function_1d*> rs ;

		//! Size of the polynomials
		//! \todo Change it by a more adaptive scheme, with different np, nq per 
		//! color channel?
		int np, nq;
} ;
