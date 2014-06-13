#include "rational_function.h"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

rational_function_1d::rational_function_1d()
{
}

rational_function_1d::rational_function_1d(int nX, unsigned int np, unsigned int nq, 
                                           bool separable) 
{
	setDimX(nX);
	setDimY(1);

	resize(np, nq);
	_separable = separable;
}

bool rational_function_1d::load(std::istream&)
{
	return true;
}

void rational_function_1d::save_body(std::ostream& out, const arguments& args) const
{
    bool is_matlab = args["export"] == "matlab";
    bool is_alta   = !args.is_defined("export") || args["export"] == "alta";
		 
	 const unsigned int np = _p_coeffs.size();
	 const unsigned int nq = _q_coeffs.size();

	 if(is_alta)
	 {
		 for(unsigned int i=0; i<np; ++i)
		 {
			 for(int j=0; j<dimX(); ++j)
			 {
				 out << _p_coeffs[i].deg[j] << "\t" ;
			 }
			 out << _p_coeffs[i].a << std::endl ;
		 }

		 for(unsigned int i=0; i<nq; ++i)
		 {
			 for(int j=0; j<dimX(); ++j)
			 {
				 out << _q_coeffs[i].deg[j] << "\t" ;
			 }
			 out << _q_coeffs[i].a << std::endl ;
		 }
	 }
	 else if(is_matlab)
	 {

		 out << "(";
		 for(unsigned int i=0; i<np; ++i)
		 {
			 out << _p_coeffs[i].a;
			 for(int k=0; k<dimX(); ++k)
			 {
				 if(k != dimX()-1) { out << ".*"; }
				 out << "x(k).^" << _q_coeffs[i].deg[k];
			 }
			 if(i != np-1) { out << " + "; }
		 }
		 out << ") / (";

		 for(unsigned int i=0; i<nq; ++i)
		 {
			 out << _p_coeffs[i].a << "x.^" << i;
			 for(int k=0; k<dimX(); ++k)
			 {
				 if(k != dimX()-1) { out << ".*"; }
				 out << "x(k).^" << _q_coeffs[i].deg[k];
			 }
			 if(i != nq-1) { out << " + "; }
		 }
		 out << ")";
	 }
	 else
	 {
		 NOT_IMPLEMENTED();
	 }
}

void rational_function_1d::update(const vec& in_a,
                                  const vec& in_b)
{
	// Get the size of the input vector
	const int np = in_a.size();
	const int nq = in_b.size();

#define NORMALIZE

#ifdef NORMALIZE
    const double b0 = (std::abs(in_b[0]) > 1.0E-16) ? in_b[0] : 1.0;
#else
    const double b0 = 1.0;
#endif

	for(int k=0; k<np; ++k)
	{
		_p_coeffs[k].a = in_a[k] / b0;
	}

	for(int k=0; k<nq; ++k)
	{
		_q_coeffs[k].a = in_b[k] / b0;
	}
}
		
void rational_function_1d::update(const rational_function_1d* r)
{
	// Get the size of the input vector
	const unsigned int np = r->_p_coeffs.size();
	const unsigned int nq = r->_q_coeffs.size();

	// Resize the rational function
	resize(np, nq);

	for(unsigned int k=0; k<np; ++k) {
		_p_coeffs[k].a = r->_p_coeffs[k].a;
	}

	for(unsigned int k=0; k<nq; ++k) {
		_q_coeffs[k].a = r->_q_coeffs[k].a;
	}

}

void rational_function_1d::resize(unsigned int np, unsigned int nq)
{
	// It is not possible to resize the multi-dimensional vector
	// if the input dimensions size is not defined. This can
	// happen at the creation of the rational function object.
	if(dimX() == 0) { return; }

	// Resize the numerator
	if(_p_coeffs.size() != np)
	{
		_p_coeffs.resize(np);
		for(unsigned int k=0; k<np; ++k)
		{
			std::vector<int> deg = index2degree(k);
			_p_coeffs[k].deg = deg;
		}
	}

	// Resize the denominator
	if(_q_coeffs.size() != nq)
	{
		_q_coeffs.resize(nq);
		for(unsigned int k=0; k<nq; ++k)
		{
			std::vector<int> deg = index2degree(k);
			_q_coeffs[k].deg = deg;
		}
	}
}



// Get the p_i and q_j function
vec rational_function_1d::p(const vec& x) const
{
	vec res(1) ;

	const unsigned int np = _p_coeffs.size() ;
	for(unsigned int i=0; i<np; ++i)
	{
		res[0] += _p_coeffs[i].a * this->p(x, i) ;
	}

	return res ;
}
vec rational_function_1d::q(const vec& x) const 
{
	vec res(1) ;

	const unsigned int np = _q_coeffs.size() ;
	for(unsigned int i=0; i<np; ++i)
	{
		res[0] += _q_coeffs[i].a * this->q(x, i) ;
	}

	return res ;
}

// Estimate the number of configuration for an indice
// vector of dimension d with maximum element value
// being k.
int rational_function_1d::estimate_dk(int k, int d)
{
	if(d == 1)
	{
		return 1;
	}
	else if(d ==2)
	{
		return k+1;
	}
	else
	{
		int res = 0;
		for(int i=0; i<=k; ++i)
		{
			res += estimate_dk(k-i, d-1);
		}
		return res;
	}
}

// Populate a vector of degrees of dimension N using a
// maximum degree of M. The index at the current level
// is j
void rational_function_1d::populate(std::vector<int>& vec, int N, int M, int j)
{
	// For each dimension, estimate the current level
	// based on the number of configurations in the
	// other dimensions
	int current_M = M ;
	int nb_conf = 0;
	for(int d=0; d<N-1; ++d)
	{
		int k;
		for(k=0; k<=current_M; ++k)
		{
			int oracle = estimate_dk(current_M-k, N-(d+1));
			if(nb_conf <= j && j < nb_conf+oracle)
			{
				break;
			}
			nb_conf += oracle;
		}

		vec[N-1 - d] = k ;	
		current_M -= k ;
	}
	vec[0] = current_M;
}

std::vector<int> rational_function_1d::index2degree(int i) const
{
	std::vector<int> deg ; deg.assign(dimX(), 0) ;

	if(i == 0)
		return deg ;
	
	// The case of one dimensional signals is trivial
	if(dimX() == 1)
	{
		deg[0] = i;
		return deg;
	}
#ifdef NOT_WORKING
	// Case of two or more dimension signal, we differ the treatment of
	// separable signals and non separable signals.
	if(_separable)
	{
		const int index = i-1;
		const int base  = index / dimX() + 1;
		for(int k=0; k<dimX(); ++k)
		{
			deg[k] = (index % dimX() == k) ? base : 0;
		}
	}
	else
#endif
	{
		if(dimX() == 2)
		{
			int Nk = 1 ;
			int k  = 1 ;
			while(!(i >= Nk && i < Nk+k+1))
			{
				Nk += k+1 ;
				++k ;
			}

			int r = i-Nk ;
			deg[0] = k-r;
			deg[1] = r;
		}
		else
		{
			int Nk = 1 ;
			int k  = 1 ;
			int dk = estimate_dk(k, dimX()) ;
			while(!(i >= Nk && i < Nk+dk))
			{
				Nk += dk ;
				++k ;
				dk = estimate_dk(k, dimX()) ;
			}

			// Populate the vector from front to back
			int j = i-Nk ;
			populate(deg, dimX(), k, j) ;
		}
	}

	return deg ;

}

// Get the p_i and q_j function
double rational_function_1d::p(const vec& x, int i) const
{
	double res = 1.0;
	for(int k=0; k<dimX(); ++k)
	{
		const double xp = 2.0*((x[k] - _min[k]) / (_max[k]-_min[k]) - 0.5);
		res *= pow(xp, _p_coeffs[i].deg[k]) ;
	}

	return res ;
}
double rational_function_1d::q(const vec& x, int i) const 
{
	double res = 1.0;
	for(int k=0; k<dimX(); ++k)
	{
		const double xp = 2.0*((x[k] - _min[k]) / (_max[k]-_min[k]) - 0.5);
		res *= pow(xp, _q_coeffs[i].deg[k]) ;
	}

	return res ;
}

// Overload the function operator
vec rational_function_1d::value(const vec& x) const 
{
	vec res(1) ;

	unsigned int const np = _p_coeffs.size();
	unsigned int const nq = _q_coeffs.size();

	double p = 0.0f ;
	double q = 0.0f ;

	for(unsigned int i=0; i<np; ++i)
	{
		p += _p_coeffs[i].a*this->p(x, i) ;
	}

	for(unsigned int i=0; i<nq; ++i)
	{
		q += _q_coeffs[i].a*this->q(x, i) ;
	}

	res[0] = p/q ;
	return res ;
}


std::ostream& operator<< (std::ostream& out, const rational_function_1d& r) 
{
	const unsigned int np = r._p_coeffs.size();
	const unsigned int nq = r._q_coeffs.size();

	std::cout << "p = [" ;
	for(unsigned int i=0; i<np; ++i)
	{
		if(i != 0)
		{
			std::cout << ", " ;
		}
		std::cout << r._p_coeffs[i].a ;
	}
	std::cout << "]" << std::endl ;

	std::cout << "q = [" ;
	for(unsigned int i=0; i<nq; ++i)
	{
		if(i != 0)
		{
			std::cout << ", " ;
		}
		std::cout << r._q_coeffs[i].a ;
	}
	std::cout << "]" << std::endl ;

	return out ;
}

#ifndef TODO
rational_function::rational_function() : np(0), nq(0)
{
}


rational_function::rational_function(int np, int nq) : np(np), nq(nq)
{
}

//! \todo clean memory here
rational_function::~rational_function()
{
}


void rational_function::update(int i, rational_function_1d* r)
{
	rs[i]->update(r);
}


rational_function_1d* rational_function::get(int i)
{
	// Check for consistency in the index of color channel
	if(i < _nY)
	{
		if(rs[i] == NULL)
		{
			rs[i] = new rational_function_1d(dimX(), np, nq);

			// Test if the input domain is not empty. If one dimension of
			// the input domain is a point, I manually inflate this dimension
			// to avoid numerical issues.
			vec _min = min();
			vec _max = max();
			for(int k=0; k<dimX(); ++k)
			{
				if(_min[k] == _max[k]) 
				{
					_min[k] -= 0.1;
					_max[k] += 0.1;
				}
			}

			rs[i]->setMin(_min) ;
			rs[i]->setMax(_max) ;
		}
		return rs[i];
	}
	else
	{
		std::cout << "<<ERROR>> tried to access out of bound 1D RF" << std::endl;
		return NULL;
	}
}

rational_function_1d* rational_function::get(int i) const
{
	// Check for consistency in the index of color channel
	if(i < _nY)
	{
		return rs[i];
	}
	else
	{
		std::cout << "<<ERROR>> tried to access out of bound 1D RF" << std::endl;
		return NULL;
	}
}

// Overload the function operator

vec rational_function::value(const vec& x) const
{
	vec res(_nY) ;

	for(int k=0; k<_nY; ++k)
	{
		res[k] = rs[k]->value(x)[0] ;
	}
	return res ;
}

// IO function to text files
bool rational_function::load(std::istream& in)
{

	// Parse line until the next comment
	while(in.peek() != '#')
	{
		char line[256];
		in.getline(line, 256);

		// If we cross the end of the file, or the badbit is
		// set, the file cannot be loaded
		if(!in.good())
			return false;
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

	vec a(_np), b(_nq);
	for(int i=0; i<_nY; ++i)
	{
		// Parse the p_i coefficients
		for(int j=0; j<_np; ++j)
		{
			for(int k=0; k<_nX; ++k)
			{
				in >> token;
			}
			in >> a[j];
		}
		
		// Parse the q_i coefficients
		for(int j=0; j<_nq; ++j)
		{
			for(int k=0; k<_nX; ++k)
			{
				in >> token;
			}
			in >> b[j];
		}
#ifdef NONE
        std::cout << "p_" << i << " = " << a << std::endl;
        std::cout << "q_" << i << " = " << b << std::endl;
#endif
		// Update the i_th color channel
		get(i)->update(a, b);
	}

	return true;
}


void rational_function::save_call(std::ostream& out, const arguments& args) const
{
	out.precision(64);
	out << std::scientific;
	out << "#FUNC rational_function" << std::endl;
	out << "#NP " << np << std::endl ;
	out << "#NQ " << nq << std::endl ;
	out << "#MIN "; for(int k=0; k<_nX; ++k) { out << _min[k] << " "; } out << std::endl;
	out << "#MAX "; for(int k=0; k<_nX; ++k) { out << _max[k] << " "; } out << std::endl; 

	for(int k=0; k<_nY; ++k)
	{
        const rational_function_1d* rf = get(k);
        rf->save_body(out, args);
	}
}

std::ostream& operator<< (std::ostream& out, rational_function& r)
{
	for(int i=0; i<r.dimY(); ++i)
	{
		rational_function_1d* rf = r.get(i);
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
