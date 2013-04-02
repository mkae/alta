#pragma once

#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>

/*! \brief A core implementation of a vector of double. 
 *
 * \details
 * This type is used for any transmission of vector data with unfixed 
 * dimension. It allows to have a generic fitter working for 
 * n-Dimensional data.
 */
class vec : public std::vector<double>
{
	public:
		// Constructor & Destructors
		//
		vec() : std::vector<double>()
		{
		}
		vec(int dim) : std::vector<double>(dim)
		{
		assign(dim, 0.0) ;
		} ;
		virtual ~vec() 
		{
		} ;
		
		// Mathematical operators
		//
		friend vec operator-(const vec& a)
		{
			vec b(a.size()) ;
			for(unsigned int i=0; i<a.size(); ++i)
			{
				b[i] = -a[i] ;
			}
			return b ;
		}
		friend vec operator-(const vec& a, const vec& b)
		{
#ifdef DEBUG
			assert(a.size() == b.size()) ;
#endif
			vec c(a.size()) ;
			for(unsigned int i=0; i<a.size(); ++i)
			{
				c[i] = a[i] - b[i];
			}
			return c ;
		}
		friend vec operator+(const vec& a, const vec& b)
		{
#ifdef DEBUG
			assert(a.size() == b.size()) ;
#endif
			vec c(a.size()) ;
			for(unsigned int i=0; i<a.size(); ++i)
			{
				c[i] = a[i] + b[i];
			}
			return c ;
		}
		friend vec operator*(const vec& a, const vec& b)
		{
#ifdef DEBUG
			assert(a.size() == b.size()) ;
#endif
			vec c(a.size()) ;
			for(unsigned int i=0; i<a.size(); ++i)
			{
				c[i] = a[i] * b[i];
			}
			return c ;
		}
		friend vec operator*(const vec& a, double b)
		{
			vec c(a.size()) ;
			for(unsigned int i=0; i<a.size(); ++i)
			{
				c[i] = a[i] * b;
			}
			return c ;
		}
		friend vec operator*(double a, const vec& b)
		{
			vec c(b.size()) ;
			for(unsigned int i=0; i<b.size(); ++i)
			{
				c[i] = a * b[i];
			}
			return c ;
		}
		friend vec operator/(const vec& a, const vec& b)
		{
#ifdef DEBUG
			assert(a.size() == b.size()) ;
#endif
			vec c(a.size()) ;
			for(unsigned int i=0; i<a.size(); ++i)
			{
				c[i] = a[i] / b[i];
			}
			return c ;
		}
		friend vec operator/(const vec& a, double b)
		{
			vec c(a.size()) ;
			for(unsigned int i=0; i<a.size(); ++i)
			{
				c[i] = a[i] / b;
			}
			return c ;
		}
		friend vec operator/(double a, const vec& b)
		{
			vec c(b.size()) ;
			for(unsigned int i=0; i<b.size(); ++i)
			{
				c[i] = a / b[i];
			}
			return c ;
		}
		friend vec normalize(const vec& a)
		{
			vec b(a.size());
			double norm = 0.0 ;
			for(unsigned int i=0; i<a.size(); ++i)
			{
				norm += a[i]*a[i];
			}
			norm = sqrt(norm);
			
			for(unsigned int i=0; i<a.size(); ++i)
			{
				b[i] = a[i]/norm;
			}
			return b;
		}
		friend double dot(const vec& a, const vec& b)
		{
#ifdef DEBUG
			assert(a.size() == b.size());
#endif
			double res = 0.0;
			for(unsigned int i=0; i<a.size(); ++i)
			{
				res += a[i]*b[i];
			}

			return res;
		}
        friend bool operator<(const vec& a, const vec& b)
        {
            bool lessthan = true ;
            for(int i=0; i<a.size(); ++i)
                lessthan &= a[i] < b[i];
            return lessthan;
        }
        friend bool operator>(const vec& a, const vec& b)
        {
            bool greatthan = true ;
            for(int i=0; i<a.size(); ++i)
                greatthan &= a[i] > b[i];
            return greatthan;
        }

/*
		friend double norm(const vec& a)
		{

		}
*/
		//! \brief IO Functions
		//
		friend std::ostream& operator<< (std::ostream& out, const vec& v)
		{
			out << "[" ;
			for(unsigned int i=0; i<v.size(); ++i)
			{
				if(i > 0) out << ", " ;
				out << v[i] ;
			}
			out << "]" ;

			return out ;
		} ;


} ;
