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
    res = _kd + _ks * pow(x[0], _N);

    return res;
}

//! Load function specific files
void phong_function::load(const std::string& filename) 
{
}

//! Save the current function to a specific file type
void phong_function::save(const std::string& filename, const arguments& args) const 
{
}

//! Number of parameters to this non-linear function
int phong_function::nbParameters() const 
{
}

//! Get the vector of parameters for the function
vec phong_function::parameters() const 
{
}

//! Update the vector of parameters for the function
void phong_function::setParameters(const vec& p) 
{
}

//! Obtain the derivatives of the function with respect to the 
//! parameters. 
vec phong_function::parametersJacobian(const vec& x) const 
{
}

Q_EXPORT_PLUGIN2(phong_function, phong_function)
