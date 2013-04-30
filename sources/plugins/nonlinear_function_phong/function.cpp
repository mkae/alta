#include "function.h"

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

function* provide_function()
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
    return 3*dimY();
}

//! Get the vector of parameters for the function
vec phong_function::parameters() const 
{
    vec res(3*dimY());
    for(int i=0; i<dimY(); ++i)
    {
        res[i]          = _kd[i];
        res[i+  dimY()] = _ks[i];
        res[i+2*dimY()] = _N[i];
    }

    return res;
}

//! Update the vector of parameters for the function
void phong_function::setParameters(const vec& p) 
{
    for(int i=0; i<dimY(); ++i)
    {
        _kd[i] = p[i];
        _ks[i] = p[i+  dimY()];
        _N[i]  = p[i+2*dimY()];
    }
}

//! Obtain the derivatives of the function with respect to the 
//! parameters. 
vec phong_function::parametersJacobian(const vec& x) const 
{
    vec jac(dimY()*nbParameters());
    for(int i=0; i<dimY(); ++i)
    {
        for(int j=0; j<dimY(); ++j)
        {
            if(i == j)
            {
                // df / dk_d
                jac[i+0*dimY() + j*nbParameters()] = 1.0;

                // df / dk_s
                jac[i+1*dimY() + j*nbParameters()] = std::pow(x[0], _N[j]);

                // df / dN
                jac[i+2*dimY() + j*nbParameters()] = _N[j] * _ks[i] * std::pow(x[0], _N[j]);
            }
            else
            {
                jac[i+0*dimY() + j*nbParameters()] = 0.0;
                jac[i+1*dimY() + j*nbParameters()] = 0.0;
                jac[i+2*dimY() + j*nbParameters()] = 0.0;
            }
        }
    }

    return jac;
}

Q_EXPORT_PLUGIN2(phong_function, phong_function)
