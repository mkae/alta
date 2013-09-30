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
void schlick::load(std::istream& in)
{
    fresnel::load(in);

    // Parse line until the next comment
    while(in.peek() != '#')
    {
        char line[256];
        in.getline(line, 256);
    }

    // Checking for the comment line #FUNC nonlinear_fresnel_schlick
    std::string token;
    in >> token;
    if(token != "#FUNC") { std::cerr << "<<ERROR>> parsing the stream. The #FUNC is not the next line defined." << std::endl; }

    in >> token;
    if(token != "nonlinear_fresnel_schlick") { std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl; }

    // R [double]
    for(int i=0; i<dimY(); ++i)
    {
        in >> token >> R[i];
    }
}

void schlick::save_call(std::ostream& out, const arguments& args) const
{
    bool is_alta   = !args.is_defined("export") || args["export"] == "alta";

    if(is_alta)
    {
       f->save_call(out, args);
    }
    else
    {
        out << "(";	f->save_call(out, args); out << ")";
    }

    if(is_alta)
    {
        out << "#FUNC nonlinear_fresnel_schlick" << std::endl ;
        for(int i=0; i<dimY(); ++i)
        {
            out << "R " << R[i] << std::endl;
        }
        out << std::endl;
    }
    else
    {
        out << " * schlick_fresnel(L, V, N, X, Y, vec3(";
        for(int i=0; i<dimY(); ++i)
        {
            out << R[i];
            if(i < _nY-1) { out << ", "; }
        }
        out << "))";
    }
}

void schlick::save_body(std::ostream& out, const arguments& args) const
{
    f->save_body(out, args);
    bool is_shader = args["export"] == "shader" || args["export"] == "explorer";

    if(is_shader)
    {
        out << std::endl;
        out << "vec3 schlick_fresnel(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, vec3 R)" << std::endl;
        out << "{" << std::endl;
        out << "\tvec3 H = normalize(L + V);" << std::endl;
        out << "\treturn R + (vec3(1.0) - R) * pow(1.0f - clamp(dot(H,V), 0.0f, 1.0f), 5);" << std::endl;
        out << "}" << std::endl;
        out << std::endl;
    }

}

vec schlick::fresnelValue(const vec& x) const
{
    double xp[3], cart[6];
    params::convert(&x[0], input_parametrization(), params::RUSIN_VH, xp);
    params::convert(&x[0], input_parametrization(), params::CARTESIAN, cart);

    const double dotVH = xp[0]*cart[0] + xp[1]*cart[1] + xp[2]*cart[2];

    vec res(_nY);
    for(int i=0; i<_nY; ++i)
    {
        res[i] = R[i] + (1.0 - R[i]) * pow(1.0 - clamp(dotVH, 0.0, 1.0), 5.0);
    }

    return res;
}

//! \brief Number of parameters to this non-linear function
int schlick::nbFresnelParameters() const
{
    return dimY();
}

vec schlick::getFresnelParametersMin() const
{
    vec m(dimY());
    for(int i=0; i<dimY(); ++i) { m[i] = 0.0; }
    return m;
}

//! Get the vector of min parameters for the function
vec schlick::getFresnelParametersMax() const
{
    vec M(dimY());
    for(int i=0; i<dimY(); ++i) { M[i] = 1.0; }
    return M;
}

//! \brief Get the vector of parameters for the function
vec schlick::getFresnelParameters() const
{
    vec p(dimY());
    for(int i=0; i<dimY(); ++i) { p[i] = R[i]; }
    return p;
}

//! \brief Update the vector of parameters for the function
void schlick::setFresnelParameters(const vec& p)
{
    for(int i=0; i<dimY(); ++i) { R[i] = p[i]; }
}

//! \brief Obtain the derivatives of the function with respect to the
//! parameters.
vec schlick::getFresnelParametersJacobian(const vec& x) const
{
    const int nY = dimY();
    double xp[3], cart[6];
    params::convert(&x[0], input_parametrization(), params::RUSIN_VH, xp);
    params::convert(&x[0], input_parametrization(), params::CARTESIAN, cart);

    const double dotVH = xp[0]*cart[0] + xp[1]*cart[1] + xp[2]*cart[2];

    vec jac(nY*nY);
    for(int i=0; i<nY; ++i)
        for(int j=0; j<nY; ++j)
        {
            if(i == j) {
                jac[j*dimY() + i] = 1.0 - pow(1.0 - clamp(dotVH, 0.0, 1.0), 5.0);
            } else {
                jac[j*dimY() + i] = 0.0;
            }
        }

    return jac;
}


void schlick::fresnelBootstrap(const data* d, const arguments& args)
{
    for(int i=0; i<dimY(); ++i) { R[i] = 1.0; }
}
