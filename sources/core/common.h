/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2014 CNRS
   Copyright (C) 2013, 2014, 2015 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

#include <vector>
#include <iostream>
#include <cassert>

// Math Constants exist in Windows but they need to be included
// by defining the constant _USE_MATH_DEFINES
// Cf https://msdn.microsoft.com/en-us/library/4hwaceh6%28v=vs.110%29.aspx
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#else
#include <cmath>
#endif

#include <cstring>
#include <algorithm>

#ifdef OLD
/*! \brief A core implementation of a vector of double.
 *  \ingroup core
 *  \internal
 *
 *  \details
 *  This type is used for any transmission of vector data with unfixed
 *  dimension. It allows to have a generic fitter working for
 *  n-Dimensional data.
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
		}

		//! \brief Return a pointer to the underlying memory region.
		const double *data() const
		{
			// std::vector<T>::at returns a reference, and since the std::vector
			// storage is guaranteed to be contiguous, we can do this.
			return &at(0);
		}

		//! \brief get a subpart of the vector
        vec subvector(int start, int n) const
        {
            vec res(n);
            memcpy(&res[0], &this->at(start), n*sizeof(double));
            return res;
        }

		//! \brief copy operator. It resize the left operand to the size of the 
		//! right operand.
    //! RP: WHy isn't the return type vec& ? 
		vec operator=(const vec& a)
		{
			this->resize(a.size());
			for(unsigned int i=0; i<a.size(); ++i)
			{
				this->at(i) = a[i];
			}

			return *this ;
		}

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
        friend double norm(const vec& a)
        {
            double norm = 0.0 ;
            for(unsigned int i=0; i<a.size(); ++i)
            {
                norm += a[i]*a[i];
            }
            return sqrt(norm);
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
            for(unsigned int i=0; i<a.size(); ++i)
                lessthan &= a[i] < b[i];
            return lessthan;
        }
        friend bool operator>(const vec& a, const vec& b)
        {
            bool greatthan = true ;
            for(unsigned int i=0; i<a.size(); ++i)
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
        }


} ;
#endif 

#include <Eigen/Core>
typedef Eigen::VectorXd vec;

double norm(const vec& a);

vec normalize(const vec& a);

double dot(const vec& a, const vec& b);

vec product(const vec& a, const vec& b);

std::ostream& operator<<(std::ostream& out, const vec& v);

//! \brief locate the first index of value v in vector vec. Complexity in
//! O(n) is the worst case.
template<typename T> int is_in(std::vector<T> ve, T v)
{
    int res = -1;
    for(unsigned int i=0; i<ve.size(); ++i)
    {
        if(ve[i] == v)
            return i;
    }

    return res;
}

template<typename T> T clamp(T x, T a, T b)
{
	return std::max<T>(std::min<T>(x, b), a);
}

#ifdef _MSC_VER
#define NOT_IMPLEMENTED() \
std::cerr << "<<ERROR>> not implemented " << __FUNCDNAME__ << " in file " << __FILE__ \
          << ":" << __LINE__ << std::endl; \
throw
#else
#define NOT_IMPLEMENTED() \
std::cerr << "<<ERROR>> not implemented " << __PRETTY_FUNCTION__ << " in file " << __FILE__ \
          << ":" << __LINE__ << std::endl; \
throw
#endif

// Mathematical definition not provided on the Window plateform
#ifdef _MSC_VER

#if (_MSC_VER < 1800)
template<typename T> bool isnan(T x)
{
	return x==std::numeric_limits<T>::signaling_NaN();
}
#endif

template<typename T> T erf(T x)
{
    // constants
    const double a1 =  0.254829592;
    const double a2 = -0.284496736;
    const double a3 =  1.421413741;
    const double a4 = -1.453152027;
    const double a5 =  1.061405429;
    const double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0) {
        sign = -1;
	 }
    x = fabs(x);

    // A&S formula 7.1.26
    const double t = 1.0/(1.0 + p*x);
    const double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return sign*y;
}

template<typename T> T erfc( T x )
{
  return 1 - erf(x);
}

#endif

#ifdef _WIN32
#define ALTA_DLL_EXPORT extern "C" __declspec(dllexport)
#else
#define ALTA_DLL_EXPORT extern "C"
#endif

/*! \brief The timer class: a cross plateform timing solution
 */
class timer
{
    public:

        //! \brief Constructor
        timer();

        //! \brief Start the timer
        void start();

        //! \brief Stop the timer
        void stop();

        //! \brief Reset the timer
        void reset();

        //! \brief Compute the time taken since the start, accounting for
        //! all pauses.
        int elapsed() const;

        //! \brief Print the time to the standard output
        void print(std::ostream& out) const;

        //! \brief ostream compliant operator
        friend std::ostream& operator<<(std::ostream& out, const timer& t);

    private:

        //! \brief get the current time in seconds. This is a system dependant
        //! function.
        //! \internal
        unsigned int current_time() const;

        unsigned int _start, _stop;
        unsigned int _elapsed;
};


// I/O error handling in user interfaces.

#include <iostream>
#include <cstdlib>

#ifndef _WIN32

# define CATCH_FILE_IO_ERROR(file)								\
		catch (std::ios_base::failure& e) {						\
				std::cerr << "<<ERROR>> failed to load '"	\
									<< (file) << "'"								\
									<< ": " << strerror(errno)			\
									<< std::endl;										\
				exit(EXIT_FAILURE);												\
		}

#else

// We cannot rely on 'errno' on Windows.

# define CATCH_FILE_IO_ERROR(file)								\
		catch (std::ios_base::failure& e) {						\
				std::cerr << "<<ERROR>> failed to load '"	\
									<< (file) << "'"								\
									<< std::endl;										\
				exit(EXIT_FAILURE);												\
		}

#endif	/* _WIN32 */
