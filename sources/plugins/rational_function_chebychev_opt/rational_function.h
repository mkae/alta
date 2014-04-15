#pragma once

// Include STL
#include <vector>
#include <string>

// Interface
#include <core/function.h>
#include <core/rational_function.h>
#include <core/data.h>
#include <core/fitter.h>
#include <core/args.h>
#include <core/common.h>

class rational_function_chebychev_1d : public rational_function_1d
{
public: // methods

    rational_function_chebychev_1d() ;
    rational_function_chebychev_1d(int nX, int np, int nq) ;
    virtual ~rational_function_chebychev_1d() {}

    // Get the p_i and q_j function
    virtual double p(const vec& x, int i) const ;
    virtual double q(const vec& x, int j) const ;

	 // Overload the function operator
	 virtual vec value(const vec& x) const
	 {
		 vec res(1) ;

		 unsigned int const np = _p_coeffs.size();
		 unsigned int const nq = _q_coeffs.size();

		 double p = 0.0f ;
		 double q = 0.0f ;

		 for(unsigned int i=0; i<np; ++i)
		 {
			 p += _p_coeffs[i].a * this->p(x, i) ;
		 }

		 for(unsigned int i=0; i<nq; ++i)
		 {
			 q += _q_coeffs[i].a * this->q(x, i) ;
		 }

		 res[0] = p/q ;
		 return res ;
	 }
	 		
	 // Update the function
	virtual void update(const vec& in_a, 
	                    const vec& in_b)
	{
		// Get the size of the input vector
		const int np = in_a.size();
		const int nq = in_b.size();

		for(int k=0; k<np; ++k)
		{
			_p_coeffs[k].a = in_a[k];
		}
		
		for(int k=0; k<nq; ++k)
		{
			_q_coeffs[k].a = in_b[k];
		}
	}

	 // Resize the polynomial
	 virtual void resize(int np, int nq)
	 {
		if(dimX() == 0) return;

		 if(_p_coeffs.size() != np)
		 {
			 _p_coeffs.resize(np);
			 for(int k=0; k<np; ++k)
			 {
				 std::vector<int> deg = index2degree(k);
				 _p_coeffs[k].deg = deg;
			 }
		 }

		 if(_q_coeffs.size() != nq)
		 {
			 _q_coeffs.resize(nq);
			 for(int k=0; k<nq; ++k)
			 {
				 std::vector<int> deg = index2degree(k);
				 _q_coeffs[k].deg = deg;
			 }
		 }
	 }

	 // Get the coefficients
	 virtual double getP(int i) const { return _p_coeffs[i].a; }
	 virtual double getQ(int i) const { return _q_coeffs[i].a; }

	 virtual vec getP() const 
	 {
		 const int np = _p_coeffs.size();
		 vec t(np);
		 for(int i=0; i<np; ++i) {t[i] = _p_coeffs[i].a; }
		 return t; 
	 }
	 virtual vec getQ() const 
	 { 
		 const int nq = _q_coeffs.size();
		 vec t(nq);
		 for(int i=0; i<nq; ++i) {t[i] = _q_coeffs[i].a; }
		 return t; 
	 }
	 

protected:  // data

	 struct coeff
	 {
		 coeff() {}
		 coeff(double a, std::vector<int> deg) :
			 a(a), deg(deg) { }

		 double a;
		 std::vector<int> deg;
	 };

	 // Table of coefficients and indices, sorted with respect
	 // to the indices.
	 std::vector<coeff> _p_coeffs;
	 std::vector<coeff> _q_coeffs;

} ;

class rational_function_chebychev : public rational_function
{
	public: // methods

		rational_function_chebychev() ;
		virtual ~rational_function_chebychev() ;
		
		//! Get the 1D function associated with color channel i. If no one exist, 
		//! this function allocates a new element. If i > nY, it returns NULL.
		virtual rational_function_1d* get(int i) ;

		//! Update the y-1D function for the ith dimension.
		//! \note It will test if the 1D function provided is of the dynamic type
		//! \name rational_function_chebychev_1d
		virtual void update(int i, rational_function_1d* r)
		{
			if(dynamic_cast<rational_function_chebychev_1d*>(r) != NULL)
			{
				rational_function::update(i, r);
			}
			else
			{
#ifdef DEBUG
				std::cerr << "<<ERROR>> the function provided is not of type \"rational_function_chebychev\"" << std::endl;
#endif
			}
		}
} ;

