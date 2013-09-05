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
    in >> token >> R;
}
		
void schlick::save_call(std::ostream& out, const arguments& args) const
{
    out << "(";	f->save_call(out, args); out << ")";
    bool is_alta   = !args.is_defined("export") || args["export"] == "alta";

    if(is_alta)
    {
        out << "#FUNC nonlinear_fresnel_schlick" << std::endl ;
        out << "R " << R << std::endl;
        out << std::endl;
    }
    else
    {
        out << " * schlick_fresnel(L, V, N, X, Y, " << R << ")";
    }
}

void schlick::save_body(std::ostream& out, const arguments& args) const
{
    f->save_body(out, args);
    bool is_shader = args["export"] == "shader" || args["export"] == "explorer";

    if(is_shader)
    {
        out << std::endl;
        out << "vec3 schlick_fresnel(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, float R)" << std::endl;
        out << "{" << std::endl;
        out << "\tvec3 H = normalize(L + V);" << std::endl;
        out << "\treturn vec3(R + (1.0f - R) * pow(1.0f - clamp(dot(H,L), 0.0f, 1.0f), 5));" << std::endl;
        out << "}" << std::endl;
    }

}


vec schlick::fresnelValue(const vec& x) const
{
	double xp[2];
	params::convert(&x[0], input_parametrization(), params::COS_TH_TD, xp);

	vec res(_nY);
	for(int i=0; i<_nY; ++i)
	{
		res[i] = R + (1.0 - R) * pow(1.0 - clamp(xp[1], 0.0, 1.0), 5.0);
	}

	return res;
}

//! \brief Number of parameters to this non-linear function
int schlick::nbFresnelParameters() const 
{
	return 1;
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
	double xp[2];
	params::convert(&x[0], input_parametrization(), params::COS_TH_TD, xp);

	vec jac(nY);
	for(int i=0; i<nY; ++i)
	{
		jac[i] = 1.0 - pow(1.0 - clamp(xp[1], 0.0, 1.0), 5.0);
	}

	return jac;
}


void schlick::fresnelBootstrap(const data* d, const arguments& args)
{
    R = 0.5;
}
