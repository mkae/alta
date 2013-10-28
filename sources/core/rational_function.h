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
		rational_function_1d(int np, int nq, bool separable = true) ;
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

#ifdef TODO
template<class T> class rational_function_t : public function
{
	public: // methods

		rational_function_t() : np(0), nq(0) {}
		rational_function_t(int np, int nq) : np(np), nq(nq) {}
		virtual ~rational_function_t() {}

		// Overload the function operator
		virtual vec value(const vec& x) const 
		{
			vec res(_nY) ;

			for(int k=0; k<_nY; ++k)
			{
				res[k] = rs[k]->value(x)[0] ;
			}
			return res ;
		}
		virtual vec operator()(const vec& x) const { return value(x) ; }

		// IO function to text files
		virtual bool load(std::istream& in) 
		{

			// Parse line until the next comment
			while(in.peek() != '#')
			{
				char line[256];
				in.getline(line, 256);
			}

			// Checking for the comment line #FUNC nonlinear_function_lafortune
			std::string token;
			in >> token;
			if(token.compare("#FUNC") != 0) 
			{ 
				std::cerr << "<<ERROR>> parsing the stream. The #FUNC is not the next line defined." << std::endl; 
				return false;
			}

			in >> token;
			if(token.compare("rational_function") != 0) 
			{
				std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl; 
				return false;
			}

			int _np, _nq;
			// Shoudl have the #NP [int]
			in >> token >> _np;

			// Shoudl have the #NQ [int]
			in >> token >> _nq;
			setSize(_np, _nq);

			// Check for the MIN and MAX vector
			vec min(dimX()), max(dimX());
			in >> token;
			if(token.compare("#MIN") != 0) 
			{
				std::cerr << "<<ERROR>> the min value for the input space is not defined." << std::endl; 
				return false;
			}
			for(int k=0; k<dimX(); ++k) {in >> min[k];}
			setMin(min);

			in >> token;
			if(token.compare("#MAX") != 0) 
			{
				std::cerr << "<<ERROR>> the max value for the input space is not defined." << std::endl; 
				return false;
			}
			for(int k=0; k<dimX(); ++k) {in >> max[k]; }
			setMax(max);


			// Check for the polynomial basis type
			in >> token;
			if(token.compare("#BASIS") != 0) 
			{
				std::cerr << "<<ERROR>> the file is not specifying the polynomial basis." << std::endl; 
				return false;
			}
			in >> token;
			if(token.compare("LEGENDRE") != 0) 
			{
				std::cerr << "<<ERROR>> the basis is different than LEGENDRE." << std::endl; 
			}

			vec a(_np), b(_nq);
			for(int i=0; i<_nY; ++i)
			{
				// Parse the p_i coefficients
				for(int j=0; j<_np; ++j)
				{
					in >> token >> a[j];
				}

				// Parse the q_i coefficients
				for(int j=0; j<_nq; ++j)
				{
					in >> token >> b[j];
				}

				std::cout << a << std::endl;
				std::cout << b << std::endl;

				// Update the i_th color channel
				get(i)->update(a, b);
			}
			return true;
		}

#ifdef OLD
		// Update the function
		virtual void update(int i, T* r) ;
#endif

		//! Get the 1D function associated with color channel i. If no one exist, 
		//! this function allocates a new element. If i > nY, it returns NULL.
		T* get(int i) 
		{
			// Check for consistency in the index of color channel
			if(i < _nY)
			{
				if(rs[i] == NULL)
				{
					rs[i] = new T(np, nq);
					rs[i]->setDimX(dimX());
					rs[i]->setDimY(dimY());
					rs[i]->setMin(getMin()) ;
					rs[i]->setMax(getMax()) ;
				}
				return rs[i];
			}
			else
			{
				std::cout << "<<ERROR>> tried to access out of bound 1D RF" 
					<< std::endl;
				return NULL;
			}
		}
		T* get(int i) const 
		{
			// Check for consistency in the index of color channel
			if(i < _nY)
			{
				return rs[i];
			}
			else
			{
				std::cout << "<<ERROR>> tried to access out of bound 1D RF" 
					<< std::endl;
				return NULL;
			}
		}

#ifdef OLD
		// STL stream ouput
		friend std::ostream& operator<<(std::ostream& out, rational_function_t& r) 
		{
			for(int i=0; i<r.dimY(); ++i)
			{
				T* rf = r.get(i);
				out << "dimension " << i << ": ";
				if(rf != NULL)
				{
					out << *rf << std::endl;
				}
				else
				{
					out << "[NULL]" << std::endl;
				}
			}

			return out ;
		}
#endif

		//! Set the dimension of the output space of the function. This function 
		//! will update the size of the rs vector size.
		virtual void setDimY(int nY) 
		{ 
			_nY = nY ;
			rs.resize(nY);
		}

		//! \brief Set the size of the rational function. Any newly created 1D 
		//! function  will have np and nq fixed as defined.
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

		//! \brief Save the rational function to the rational format (see \ref 
		//! formating).
		virtual void save_call(std::ostream& out, const arguments& args) const 
		{
			out << "#FUNC rational_function" << std::endl;
			out << "#NP " << np << std::endl ;
			out << "#NQ " << nq << std::endl ;
			out << "#MIN "; for(int k=0; k<_nX; ++k) { out << _min[k] << " "; } 
			out << std::endl;
			out << "#MAX "; for(int k=0; k<_nX; ++k) { out << _max[k] << " "; } 
			out << std::endl; 
			out << "#BASIS LEGENDRE" << std::endl ;

			for(int k=0; k<_nY; ++k)
			{
				T* rf = get(k);
				vec a = rf->getP();
				vec b = rf->getQ();

				for(int i=0; i<np; ++i)
				{
					std::vector<int> index = rf->index2degree(i) ;
					for(unsigned int j=0; j<index.size(); ++j)
					{
						out << index[j] << "\t" ;
					}
					out << a[i] << std::endl ;
				}

				for(int i=0; i<nq; ++i)
				{
					std::vector<int> index = rf->index2degree(i) ;
					for(unsigned int j=0; j<index.size(); ++j)
					{
						out << index[j] << "\t" ;
					}
					out << b[i] << std::endl ;
				}
			}
		}

	protected: // data

		//! Store the y \in R rational functions. Each channel is a distinct 
		//! polynomial and should be fitted separately.
		std::vector<T*> rs ;

		//! Size of the polynomials
		//! \todo Change it by a more adaptive scheme, with different np, nq per 
		//! color channel?
		int np, nq;
}

typedef rational_function rational_function_t<rational_function_1d>;

#else
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
#endif
