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
    return new diffuse_function();
}

// Overload the function operator
vec diffuse_function::operator()(const vec& x) const 
{
	return value(x);
}
vec diffuse_function::value(const vec& x) const 
{
    vec res(dimY());
    for(int i=0; i<dimY(); ++i)
    {
        res[i] = _kd[i];
    }

    return res;
}

//! Load function specific files
void diffuse_function::load(const std::string& filename) 
{
	NOT_IMPLEMENTED();
}
		
void diffuse_function::save_alta(const std::string& filename, const arguments& args) const
{
	NOT_IMPLEMENTED();
}

//! Number of parameters to this non-linear function
int diffuse_function::nbParameters() const 
{
#ifdef FIT_DIFFUSE
    return dimY();
#else
    return 0;
#endif
}

//! Get the vector of parameters for the function
vec diffuse_function::parameters() const 
{
#ifdef FIT_DIFFUSE
    vec res(dimY());
    for(int i=0; i<dimY(); ++i)
    {
        res[i*3 + 0] = _kd[i];
    }
#else
    vec res(0);
#endif

    return res;
}

//! Update the vector of parameters for the function
void diffuse_function::setParameters(const vec& p) 
{
#ifdef FIT_DIFFUSE
    for(int i=0; i<dimY(); ++i)
    {
        _kd[i] = p[i];
    }
#endif
}

//! Obtain the derivatives of the function with respect to the 
//! parameters. 
vec diffuse_function::parametersJacobian(const vec& x) const 
{
#ifdef FIT_DIFFUSE
	vec jac(dimY());
	for(int i=0; i<dimY(); ++i)
		for(int j=0; j<dimY(); ++j)
		{
			if(i == j)
			{
				// df / dk_d
				jac[i*dimY() + j] = 1.0;

			}
			else
			{
				jac[i*dimY() + j] = 0.0;
			}
		}
#else
	vec jac(0);
#endif

    return jac;
}


void diffuse_function::bootstrap(const data* d, const arguments& args)
{
	// Set the diffuse component
	vec x0 = d->get(0);
	for(int i=0; i<d->dimY(); ++i)
		_kd[i] = x0[d->dimX() + i];

	for(int i=1; i<d->size(); ++i)
	{
		vec xi = d->get(i);
		for(int j=0; j<d->dimY(); ++j)
			_kd[j] = std::min(xi[d->dimX() + j], _kd[j]);
	}
	std::cout << "<<INFO>> found diffuse: " << _kd << std::endl;

}

