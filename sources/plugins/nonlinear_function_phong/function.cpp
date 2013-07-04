#include "function.h"

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <core/common.h>

ALTA_DLL_EXPORT function* provide_function()
{
    return new phong_function();
}

// Overload the function operator
vec phong_function::operator()(const vec& x) const 
{
	return value(x);
}
vec phong_function::value(const vec& x) const 
{
    vec res(dimY());
    for(int i=0; i<dimY(); ++i)
    {
        res[i] = _kd[i] + _ks[i] * std::pow(x[0], _N[i]);
    }

    return res;
}

//! Load function specific files
void phong_function::load(const std::string& filename) 
{
    std::cerr << "Not implemented " << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
}

//! Number of parameters to this non-linear function
int phong_function::nbParameters() const 
{
#ifdef FIT_DIFFUSE
    return 3*dimY();
#else
    return 2*dimY();
#endif
}

//! Get the vector of parameters for the function
vec phong_function::parameters() const 
{
#ifdef FIT_DIFFUSE
    vec res(3*dimY());
    for(int i=0; i<dimY(); ++i)
    {
        res[i*3 + 0] = _kd[i];
        res[i*3 + 1] = _ks[i];
        res[i*3 + 2] = _N[i];
    }
#else
    vec res(2*dimY());
    for(int i=0; i<dimY(); ++i)
    {
        res[i*2 + 0] = _ks[i];
        res[i*2 + 1] = _N[i];
    }
#endif

    return res;
}

//! Update the vector of parameters for the function
void phong_function::setParameters(const vec& p) 
{
    for(int i=0; i<dimY(); ++i)
    {
#ifdef FIT_DIFFUSE
        _kd[i] = p[i*3 + 0];
        _ks[i] = p[i*3 + 1];
        _N[i]  = p[i*3 + 2];
#else
        _ks[i] = p[i*2 + 0];
        _N[i]  = p[i*2 + 1];
#endif
    }
}

//! Obtain the derivatives of the function with respect to the 
//! parameters. 
vec phong_function::parametersJacobian(const vec& x) const 
{
    vec jac(dimY()*nbParameters());
	 for(int i=0; i<dimY(); ++i)
		 for(int j=0; j<dimY(); ++j)
		 {
			 if(i == j)
			 {
#ifdef FIT_DIFFUSE
				 // df / dk_d
				 jac[i*nbParameters() + j*3+0] = 1.0;

				 // df / dk_s
				 jac[i*nbParameters() + j*3+1] = std::pow(x[0], _N[j]);

				 // df / dN
				 if(x[0] == 0.0)
					 jac[i*nbParameters() + j*3+2] = 0.0;
				 else
					 jac[i*nbParameters() + j*3+2] = _ks[j] * std::log(x[0]) * std::pow(x[0], _N[j]);
#else

                 // df / dk_s
                 jac[i*nbParameters() + j*2+0] = std::pow(x[0], _N[j]);

                 // df / dN
                 if(x[0] == 0.0)
                     jac[i*nbParameters() + j*2+1] = 0.0;
                 else
                     jac[i*nbParameters() + j*2+1] = _ks[j] * std::log(x[0]) * std::pow(x[0], _N[j]);
#endif
             }
			 else
			 {
#ifdef FIT_DIFFUSE
				 jac[i*nbParameters() + j*3+0] = 0.0;
				 jac[i*nbParameters() + j*3+1] = 0.0;
				 jac[i*nbParameters() + j*3+2] = 0.0;
#else
                 jac[i*nbParameters() + j*2+0] = 0.0;
                 jac[i*nbParameters() + j*2+1] = 0.0;
#endif
			 }
		 }

    return jac;
}


void phong_function::boostrap(const data* d, const arguments& args)
{

    vec x0 = d->get(0);
    for(int i=0; i<d->dimY(); ++i)
        _kd[i] = x0[d->dimX() + i];

    for(int i=1; i<d->size(); ++i)
    {
        vec xi = d->get(i);
        for(int j=0; j<d->dimY(); ++j)
            _kd[j] = std::min(xi[d->dimX() + j], _kd[j]);
    }
}

std::ofstream& type_affectation(std::ofstream& out, 
                                const std::string& name, 
                                const vec& x, int  nY)
{
	if(nY == 1)
		out << "float " ;
	else
		out << "vec" << nY << " ";
	
	out << name << " = ";

	if(nY != 1)
		out << "vec" << nY << "(";

	for(int i=0; i<nY; ++i)
	{
		if(i != 0) out << ", ";
		out << x[i];
	}

	if(nY != 1)
		out << ")";

	out << ";" << std::endl;

	return out;
}

//! \brief Output the function using a BRDF Explorer formating.
void phong_function::save_brdfexplorer(const std::string& filename,
                                       const arguments& args) const
{
    std::ofstream file(filename.c_str(), std::ios_base::trunc);
    file << "analytic" << std::endl;

    file << std::endl;
//    file << "::begin parameters" << std::endl;
//    file << "::end parameters" << std::endl;
//    file << std::endl;
    file << std::endl;
    file << "::begin shader" << std::endl;
	 
	 type_affectation(file, std::string("n"),  _N,  _nY);
	 type_affectation(file, std::string("kd"), _kd, _nY);
	 type_affectation(file, std::string("ks"), _ks, _nY);

    file << std::endl;
    file << "const float PI = 3.14159265358979323846;" << std::endl;
    file << std::endl;
    file << "vec3 BRDF( vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y )" << std::endl;
    file << "{" << std::endl;
    file << "    vec3 H = normalize(L+V);" << std::endl;
    if(_nY == 1)
    {
        file << "    float D = kd + ks * pow(max(0, dot(N,H)),n);" << std::endl;
    }
    else
    {
        file << "    vec" << _nY << " D = kd + ks * pow(vec" << _nY << "(max(0.0, dot(N,H))),n);" << std::endl;
    }
	 file << "    return vec3(D);" << std::endl;
//    file << "    if (normalized)" << std::endl;
//    file << "        D *= (2+n) / (2*PI);" << std::endl;
    file << "}" << std::endl;
    file << std::endl;
    file << "::end shader" << std::endl;
    file.close();
}


Q_EXPORT_PLUGIN2(phong_function, phong_function)
