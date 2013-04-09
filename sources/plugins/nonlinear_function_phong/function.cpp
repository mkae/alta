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
    /*
    // shading
    vec lv(3); lv[0] = x[0]; lv[1] = x[1]; lv[2] = x[2];
    vec n(3); n[0] = 0.0; n[1] = 0.0; n[2] = 1.0;
    vec ev(3); ev[0] = x[3]; ev[1] = x[4]; ev = x[5];
    vec halfVector = normalize(lv + ev);

    double v_h = dot(ev, halfVector);
    double n_h = dot(n, halfVector);
    double n_l = dot(n, lv);
    double inLight = 1.0;
    if (n_l < 0.0) inLight = 0.0;
    double n_v = dot(n, ev);

    return one_pi * inLight * (n_l * rho_d + rho_s *
          D(alpha, p, n_h, K_ap) * G1(n_l) * G1 (n_v) *
             Fresnel(F_0, F_1, v_h));
    */
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
