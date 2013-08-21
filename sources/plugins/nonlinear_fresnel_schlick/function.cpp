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

// Overload the function operator
vec schlick::operator()(const vec& x) const 
{
	return value(x);
}
vec schlick::value(const vec& x) const 
{
	vec res(_nY);
	for(int i=0; i<_nY; ++i)
	{
		res[i] = R + (1.0 - R) * pow(1.0 - clamp(x[0], 0.0, 1.0), 5.0);
	}

	return res;
}

//! Load function specific files
void schlick::load(const std::string& filename) 
{
    std::cerr << "Cannot load a Schlick file." << std::endl;
    throw;
}

//! Number of parameters to this non-linear function
int schlick::nbParameters() const 
{
	return 0;
}

//! Get the vector of parameters for the function
vec schlick::parameters() const 
{
	vec r(0);
	return r;
}

//! Update the vector of parameters for the function
void schlick::setParameters(const vec& p) 
{
}

//! Obtain the derivatives of the function with respect to the 
//! parameters. 
vec schlick::parametersJacobian(const vec& x) const 
{
	vec r(0);
	return r;
}


void schlick::bootstrap(const data* d, const arguments& args)
{
}
