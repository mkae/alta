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
    return new schlick();
}

//! Load function specific files
void schlick::load(const std::string& filename) 
{
    std::cerr << "Cannot load a Schlick file." << std::endl;
    throw;
}

vec schlick::fresnelValue(const vec& x) const
{
	vec res(_nY);
	for(int i=0; i<_nY; ++i)
	{
		res[i] = R + (1.0 - R) * pow(1.0 - clamp(x[0], 0.0, 1.0), 5.0);
	}

	return res;
}

//! \brief Number of parameters to this non-linear function
int schlick::nbFresnelParameters() const 
{
	return dimY();
}

//! \brief Get the vector of parameters for the function
vec schlick::getFresnelParameters() const 
{
	vec p(1);
	p[0] = R;
	return p;
}

//! \brief Update the vector of parameters for the function
void schlick::setFresnelParameters(const vec& p) 
{
	R = p[0];
}

//! \brief Obtain the derivatives of the function with respect to the
//! parameters. 
vec schlick::getFresnelParametersJacobian(const vec& x) const 
{
	const int nY = dimY();

	vec jac(nY);
	for(int i=0; i<nY; ++i)
	{
		jac[i] = 1.0 - pow(1.0 - clamp(x[0], 0.0, 1.0), 5.0);
	}

	return jac;
}


void schlick::bootstrap(const data* d, const arguments& args)
{
	R = 1.0;
}
