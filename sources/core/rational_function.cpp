#include "rational_function.h"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>


rational_function::rational_function() : a(), b()
{
}

rational_function::rational_function(int np, int nq) : a(np), b(nq)
{
}

rational_function::rational_function(const std::vector<double>& a,
                         const std::vector<double>& b) :
	a(a), b(b)
{
}
rational_function::~rational_function()
{
}
		
void rational_function::update(const std::vector<double>& in_a,
                               const std::vector<double>& in_b)
{
	a.reserve(in_a.size()) ;
	b.reserve(in_b.size()) ;
	a = in_a ;
	b = in_b ;
}

// Overload the function operator
vec rational_function::value(const vec& x) const 
{
	vec res(_nY) ;

    unsigned int const np = a.size() / _nY ;
    unsigned int const nq = b.size() / _nY ;

	for(int k=0; k<_nY; ++k)
	{
		double p = 0.0f ;
		double q = 0.0f ;
		
		for(unsigned int i=0; i<np; ++i)
		{
			p += a[k*_nY + i]*this->p(x, i) ;
		}

		for(unsigned int i=0; i<nq; ++i)
		{
			q += b[k*_nY + i]*this->q(x, i) ;
		}

		res[k] = p/q ;
	}
	return res ;
}

// Get the p_i and q_j function
vec rational_function::p(const vec& x) const
{
	vec res(_nY) ;

    unsigned int const np = a.size() / _nY ;

	for(int k=0; k<_nY; ++k)
	{
		double p = 0.0f ;
		
		for(unsigned int i=0; i<np; ++i)
		{
			p += a[k*_nY + i]*this->p(x, i) ;
		}

		res[k] = p ;
	}
	return res ;
}
vec rational_function::q(const vec& x) const 
{
	vec res(_nY) ;

    unsigned int const nq = b.size() / _nY ;

	for(int k=0; k<_nY; ++k)
	{
		double q = 0.0f ;

		for(unsigned int i=0; i<nq; ++i)
		{
			q += b[k*_nY + i]*this->q(x, i) ;
		}

		res[k] = q ;
	}
	return res ;
}

// Estimate the number of configuration for an indice
// vector of dimension d with maximum element value
// being k.
int rational_function::estimate_dk(int k, int d)
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
void rational_function::populate(std::vector<int>& vec, int N, int M, int j)
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

std::vector<int> rational_function::index2degree(int i) const
{
	std::vector<int> deg ; deg.assign(dimX(), 0) ;

	if(i == 0)
		return deg ;
		
	if(dimX() == 1)
	{
	    deg[0] = i;
	}
    else if(dimX() == 2)
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

	return deg ;

}

double legendre(double x, int i)
{
	if(i == 0)
	{
		return 1;
	}
	else if(i == 1)
	{
		return x;
	}
	else
	{
		return ((2*i-1)*x*legendre(x, i-1) - (i-1)*legendre(x, i-2)) / (double)i ;
	}
}

//#define POLYNOMIALS

// Get the p_i and q_j function
double rational_function::p(const vec& x, int i) const
{
	std::vector<int> deg = index2degree(i);
	double res = 1.0;
	for(int k=0; k<dimX(); ++k)
	{
#ifdef POLYNOMIALS
		res *= pow(x[k], deg[k]) ;
#else // LEGENDRE
		res *= legendre(2.0*((x[k] - _min[k]) / (_max[k]-_min[k]) - 0.5), deg[k]);
#endif
	}

	return res ;
}
double rational_function::q(const vec& x, int i) const 
{
	std::vector<int> deg = index2degree(i);
	double res = 1.0; 
	for(int k=0; k<dimX(); ++k)
	{
#ifdef POLYNOMIALS
		res *= pow(x[k], deg[k]) ;
#else // LEGENDRE
		res *= legendre(2.0*((x[k] - _min[k]) / (_max[k]-_min[k]) - 0.5), deg[k]);
#endif
	}

	return res ;
}

// IO function to text files
void rational_function::load(const std::string& filename)
{
	std::ifstream file(filename.c_str()) ;
	if(!file.is_open())
	{
		std::cerr << "<<ERROR>> unable to open file \"" << filename << "\"" << std::endl ;
		throw ;
	}

	int np, nq ;
	int nX, nY ;
	vec xmin, xmax ;
	int i = 0, j = 0;
	while(file.good())
	{
		std::string line ;
		std::getline(file, line) ;
		std::stringstream linestream(line) ;
		
		// Discard incorrect lines
		if(linestream.peek() == '#')
		{
			linestream.ignore(1) ;

			std::string comment ;
			linestream >> comment ;

			if(comment == std::string("DIM"))
			{
				linestream >> nX >> nY ;
				setDimX(nX) ;
				setDimY(nY) ;

				xmin.resize(nX) ;
				xmax.resize(nX) ;
				for(int k=0; k<nX; ++k)
					xmax[k] = 1.0;

				setMin(xmin) ;
				setMax(xmax) ;
			}
			else if(comment == std::string("NP"))
			{
				linestream >> np ;
				a.resize(np*nY);
			}
			else if(comment == std::string("NQ"))
			{
				linestream >> nq ;
				b.resize(nq*nY);
			}
			else if(comment == std::string("INPUT_PARAM"))
			{
				std::string param;
				linestream >> param ;
				setParametrization(params::parse_input(param));
			}
			continue ;
		} 
		else if(line.empty())
		{
			continue ;
		}
		else if(j < nY)
		{
			int index ; double val ;
			
			// Accessing the index
			for(int k=0; k<nX; ++k) {
				linestream >> index ;
			}

			// Accessing the value
			linestream >> val ;
			
			if(i < np)
			{
				a[i + np*j] = val ;
			}
			else
			{
				b[i-np + nq*j] = val ;
			}

			if(i < np+nq) {
				++i ;
			} else {
				i = 0 ;
				++j ;
			}
		}
	}	
/*
	for(int i=0; i<a.size(); ++i) {
		std::cout << a[i] << "\t" ;
	}
	for(int i=0; i<b.size(); ++i) {
		std::cout << b[i] << "\t" ;
	}
*/
}

//! \todo it should handle parametrization
void rational_function::save_matlab(const std::string& filename, const arguments& args) const
{
    unsigned int np = a.size() / _nY ;
    unsigned int nq = b.size() / _nY ;

    std::ofstream file(filename.c_str(), std::ios_base::trunc);


    file << "function y = brdf(x)" << std::endl;
    file << std::endl;
    file << "\ts = [";
    for(int i=0; i<dimX(); ++i)
    {
        file << 1.0 / (_max[i]-_min[i]);
        if(i < dimX()-1)
        {
            file << ", ";
        }
    }
    file << "];" << std::endl;
    file << "\tc = [";
    for(int i=0; i<dimX(); ++i)
    {
        file << _min[i];
        if(i < dimX()-1)
        {
            file << ", ";
        }
    }
    file << "];" << std::endl;
    file << std::endl ;

    // Export each color channel independantly
    for(int j=0; j<dimY(); ++j)
    {
        // Export the numerator of the jth color channel
        file << "\tp(" << j+1 << ",:) = ";
        for(unsigned int i=0; i<np; ++i)
        {
            if(i > 0 && a[np*j + i] >= 0.0)
                file << " + ";
            else if(a[np*j + i] < 0.0)
                file << " " ;

            file << a[np*j + i];

            std::vector<int> degree = index2degree(i);
            for(unsigned int k=0; k<degree.size(); ++k)
            {
               file << ".*legendrepoly(" << degree[k] << ", 2.0*((x(" << k+1 << ",:)"
					     << "-c(" << k+1 << "))*s(" << k+1 << ") - 0.5))" ;
            }
        }
        file << ";" << std::endl;

        // Export the denominator of the jth color channel
        file << "\tq(" << j+1 << ",:) = ";
        for(unsigned int i=0; i<nq; ++i)
        {
            if(i > 0 && b[np*j + i] >= 0.0)
                file << " + ";
            else if(b[np*j + i] < 0.0)
                file << " " ;

            file << b[np*j + i] ;

            std::vector<int> degree = index2degree(i);
            for(unsigned int k=0; k<degree.size(); ++k)
            {
               file << ".*legendrepoly(" << degree[k] << ", 2.0*((x(" << k+1 << ",:)" 
					     << "-c(" << k+1 << "))*s(" << k+1 << ") - 0.5))" ;
            }
        }
        file << ";" << std::endl;

        file << "\ty(" << j+1 << ",:) = p./q;" << std::endl;
        if(j < dimY()-1)
        {
            file << std::endl;
        }
    }


    file << "endfunction" << std::endl;

    file.close() ;
}

//! \todo it should handle parametrization
void rational_function::save_cpp(const std::string& filename, const arguments& args) const
{
    unsigned int np = a.size() / _nY ;
    unsigned int nq = b.size() / _nY ;

    std::ofstream file(filename.c_str(), std::ios_base::trunc);

    file << "double s[" << dimX() << "] = {";
    for(int i=0; i<dimX(); ++i)
    {
        file << 1.0 / (_max[i]-_min[i]);
        if(i < dimX()-1)
        {
            file << ", ";
        }
    }
    file << "};" << std::endl;
    file << "double c[" << dimX() << "] = {";
    for(int i=0; i<dimX(); ++i)
    {
        file << _min[i];
        if(i < dimX()-1)
        {
            file << ", ";
        }
    }
    file << "};" << std::endl;
    file << std::endl ;

    file << "// The Legendre polynomial of order i evaluated in x" << std::endl;
    file << "double l(double x, int i)" << std::endl;
    file << "{" << std::endl;
    file << "    if(i == 0)" << std::endl;
    file << "    {" << std::endl;
    file << "        return 1;" << std::endl;
    file << "    }" << std::endl;
    file << "    else if(i == 1)" << std::endl;
    file << "    {" << std::endl;
    file << "        return x;" << std::endl;
    file << "    }" << std::endl;
    file << "    else" << std::endl;
    file << "    {" << std::endl;
    file << "        return ((2*i-1)*x*l(x, i-1) - (i-1)*l(x, i-2)) / (double)i ;" << std::endl;
    file << "    }" << std::endl;
    file << "}" << std::endl;
    file << std::endl;

    file << "void brdf(double* x, double* y)" << std::endl;
    file << "{" << std::endl;
    file << "\tdouble p, q;" << std::endl;

    // Export each color channel independantly
    for(int j=0; j<dimY(); ++j)
    {
        // Export the numerator of the jth color channel
        file << "\tp = ";
        for(unsigned int i=0; i<np; ++i)
        {
            if(i > 0 && a[np*j + i] >= 0.0)
            {
                file << " + ";
            }
            file << a[np*j + i];

            std::vector<int> degree = index2degree(i);
            for(unsigned int k=0; k<degree.size(); ++k)
            {
                file << "*l(2.0*((x\[" << k << "\]-c[" << k << "])*s[" << k << "] - 0.5), " << degree[k] << ")" ;
            }
        }
        file << ";" << std::endl;

        // Export the denominator of the jth color channel
        file << "\tq = ";
        for(unsigned int i=0; i<nq; ++i)
        {
            if(i > 0 && b[np*j + i] >= 0.0)
                file << " + ";
            else if(b[np*j + i] < 0.0)
                file << " " ;

            file << b[np*j + i] ;

            std::vector<int> degree = index2degree(i);
            for(unsigned int k=0; k<degree.size(); ++k)
            {
                file << "*l(2.0*((x\[" << k << "\]-c[" << k << "])*s[" << k << "] - 0.5), " << degree[k] << ")" ;
            }
        }
        file << ";" << std::endl;

        file << "\ty[" << j << "] = p/q;" << std::endl;
        if(j < dimY()-1)
        {
            file << std::endl;
        }
    }


    file << "}" << std::endl;

    file.close() ;
}

void rational_function::save_gnuplot(const std::string& filename, const data* d, const arguments& args) const 
{
	std::ofstream file(filename.c_str(), std::ios_base::trunc);
	for(int i=0; i<d->size(); ++i)
	{
		vec v = d->get(i) ;
//		vec y1 ; y1.assign(d->dimY(), 0.0) ;
//		for(int k=0; k<d->dimY(); ++k) { y1[k] = v[d->dimX() + k] ; }

		vec y2 = value(v) ;
		for(int u=0; u<d->dimX(); ++u)
			file << v[u] << "\t" ;

		for(int u=0; u<d->dimY(); ++u)
			file << y2[u] << "\t" ;

		file << std::endl ;
	}
	file.close();
}

void rational_function::save(const std::string& filename) const 
{
	std::ofstream file(filename.c_str(), std::ios_base::trunc);
	file << "#DIM " << _nX << " " << _nY << std::endl ;
	file << "#NP " << a.size() / _nY << std::endl ;
	file << "#NQ " << b.size() / _nY << std::endl ;
	file << "#BASIS LEGENDRE" << std::endl ;
	file << "#INPUT_PARAM " << params::get_name(this->parametrization()) << std::endl;

	unsigned int np = a.size() / _nY ;
    unsigned int nq = b.size() / _nY ;
	for(int k=0; k<_nY; ++k)
	{
		for(unsigned int i=0; i<np; ++i)
		{
			std::vector<int> index = index2degree(i) ;
			for(unsigned int j=0; j<index.size(); ++j)
			{
				file << index[j] << "\t" ;
			}
			file << a[i+np*k] << std::endl ;
		}

		for(unsigned int i=0; i<nq; ++i)
		{
			std::vector<int> index = index2degree(i) ;
			for(unsigned int j=0; j<index.size(); ++j)
			{
				file << index[j] << "\t" ;
			}
			file << b[i+nq*k] << std::endl ;
		}
	}

}

std::ostream& operator<< (std::ostream& out, const rational_function& r) 
{
	std::cout << "p = [" ;
	for(unsigned int i=0; i<r.a.size(); ++i)
	{
		if(i != 0)
		{
			std::cout << ", " ;
		}
		std::cout << r.a[i] ;
	}
	std::cout << "]" << std::endl ;

	std::cout << "q = [" ;
	for(unsigned int i=0; i<r.b.size(); ++i)
	{
		if(i != 0)
		{
			std::cout << ", " ;
		}
		std::cout << r.b[i] ;
	}
	std::cout << "]" << std::endl ;

	return out ;
}


//Q_EXPORT_PLUGIN2(rational_function, rational_function)
