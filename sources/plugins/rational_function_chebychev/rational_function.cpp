#include "rational_function.h"

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <core/common.h>

ALTA_DLL_EXPORT function* provide_function()
{
    return new rational_function_chebychev();
}

rational_function_chebychev::rational_function_chebychev() : rational_function() 
{
}

rational_function_chebychev::rational_function_chebychev(const std::vector<double>& a,
                         const std::vector<double>& b) : rational_function(a, b)
{
}
rational_function_chebychev::~rational_function_chebychev()
{
}


rational_function_chebychev_1d::rational_function_chebychev_1d()
{
}

rational_function_chebychev_1d::rational_function_chebychev_1d(int np, int nq) :
    rational_function_1d(np, nq)
{
}

rational_function_chebychev_1d::rational_function_chebychev_1d(const std::vector<double>& a,
                                                               const std::vector<double>& b) :
    rational_function_1d(a, b)
{
}
		

double chebychev(double x, int i)
{
#ifdef RECURSIVE_FORM
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
		return 2.0*x*chebychev(x, i-1) - chebychev(x, i-2) ;
	}
#else
    return cos(i * acos(x));
#endif
}

// Get the p_i and q_j function
double rational_function_chebychev_1d::p(const vec& x, int i) const
{
	std::vector<int> deg = index2degree(i);
	double res = 1.0;
	for(int k=0; k<dimX(); ++k)
	{
		double xk = 2.0*((x[k] - _min[k]) / (_max[k]-_min[k]) - 0.5);
		res *= chebychev(xk, deg[k]);
	}

	return res ;
}
double rational_function_chebychev_1d::q(const vec& x, int i) const
{
	std::vector<int> deg = index2degree(i);
	double res = 1.0; 
	for(int k=0; k<dimX(); ++k)
	{
		double xk = 2.0*((x[k] - _min[k]) / (_max[k]-_min[k]) - 0.5);
		res *= chebychev(xk, deg[k]);
	}

	return res ;
}

//! \todo it should handle parametrization
void rational_function_chebychev::save_matlab(const std::string& filename, const arguments& args) const
{
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
        rational_function_1d* rf = get(j);
        std::vector<double> a = rf->getP();
        std::vector<double> b = rf->getQ();

        const unsigned int np = a.size();
        const unsigned int nq = b.size();

        // Export the numerator of the jth color channel
        file << "\tp(" << j+1 << ",:) = ";
        for(unsigned int i=0; i<np; ++i)
        {
            if(i > 0 && a[np*j + i] >= 0.0)
                file << " + ";
            else if(a[np*j + i] < 0.0)
                file << " " ;

            file << a[np*j + i];

            std::vector<int> degree = rf->index2degree(i);
            for(unsigned int k=0; k<degree.size(); ++k)
            {
               file << ".*chebychevpoly(1," << degree[k] << ", 2.0*((x(" << k+1 << ",:)"
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

            std::vector<int> degree = rf->index2degree(i);
            for(unsigned int k=0; k<degree.size(); ++k)
            {
               file << ".*chebychevpoly(1," << degree[k] << ", 2.0*((x(" << k+1 << ",:)"
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
void rational_function_chebychev::save_cpp(const std::string& filename, const arguments& args) const
{
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

    file << "// The Chebychev polynomial of order i evaluated in x" << std::endl;
    file << "double l(double x, int i)" << std::endl;
    file << "{" << std::endl;
    file << "\treturn cos(i * acos(x));" << std::endl;
    file << "}" << std::endl;
    file << std::endl;

    file << "void brdf(double* x, double* y)" << std::endl;
    file << "{" << std::endl;
    file << "\tdouble p, q;" << std::endl;

    // Export each color channel independantly
    for(int j=0; j<dimY(); ++j)
    {
        rational_function_1d* rf = get(j);
        std::vector<double> a = rf->getP();
        std::vector<double> b = rf->getQ();

        const unsigned int np = a.size();
        const unsigned int nq = b.size();

        // Export the numerator of the jth color channel
        file << "\tp = ";
        for(unsigned int i=0; i<np; ++i)
        {
            if(i > 0 && a[np*j + i] >= 0.0)
            {
                file << " + ";
            }
            file << a[np*j + i];

            std::vector<int> degree = rf->index2degree(i);
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

            std::vector<int> degree = rf->index2degree(i);
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

void rational_function_chebychev::save(const std::string& filename) const
{
    std::ofstream file(filename.c_str(), std::ios_base::trunc);
    file << "#DIM " << _nX << " " << _nY << std::endl ;
    file << "#NP " << np << std::endl ;
    file << "#NQ " << nq << std::endl ;
    file << "#BASIS CHEBYCHEV" << std::endl ;

    for(int k=0; k<_nY; ++k)
    {
        rational_function_1d* rf = get(k);
        std::vector<double> a = rf->getP();
        std::vector<double> b = rf->getQ();

        const unsigned int np = a.size();
        const unsigned int nq = b.size();

        for(unsigned int i=0; i<np; ++i)
        {
            std::vector<int> index = rf->index2degree(i) ;
            for(unsigned int j=0; j<index.size(); ++j)
            {
                file << index[j] << "\t" ;
            }
            file << a[i+np*k] << std::endl ;
        }

        for(unsigned int i=0; i<nq; ++i)
        {
            std::vector<int> index = rf->index2degree(i) ;
            for(unsigned int j=0; j<index.size(); ++j)
            {
                file << index[j] << "\t" ;
            }
            file << b[i+nq*k] << std::endl ;
        }
    }

}



//Q_EXPORT_PLUGIN2(rational_function_chebychev, rational_function_chebychev)
