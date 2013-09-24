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
    return new blinn_function();
}

// Overload the function operator
vec blinn_function::operator()(const vec& x) const 
{
	return value(x);
}
vec blinn_function::value(const vec& x) const 
{
    vec res(dimY());
    for(int i=0; i<dimY(); ++i)
    {
        res[i] = _ks[i] * std::pow(x[0], _N[i]);
    }

    return res;
}

//! Load function specific files
void blinn_function::load(std::istream& in) 
{
	    // Parse line until the next comment
    while(in.peek() != '#')
    {
        char line[256];
        in.getline(line, 256);
    }

    // Checking for the comment line #FUNC nonlinear_function_blinn
    std::string token;
    in >> token;
    if(token != "#FUNC") { std::cerr << "<<ERROR>> parsing the stream. The #FUNC is not the next line defined." << std::endl; }

    in >> token;
    if(token != "nonlinear_function_blinn") { std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl; }

    // ks [double]
	 // N  [double]
    for(int i=0; i<dimY(); ++i)
    {
        in >> token >> _ks[i];
        in >> token >>  _N[i];
    }

    std::cout << "<<DEBUG>> load parameters " << parameters() << std::endl;
}

//! Number of parameters to this non-linear function
int blinn_function::nbParameters() const 
{
    return 2*dimY();
}

//! Get the vector of parameters for the function
vec blinn_function::parameters() const 
{
    vec res(2*dimY());
    for(int i=0; i<dimY(); ++i)
    {
        res[i*2 + 0] = _ks[i];
        res[i*2 + 1] = _N[i];
    }

    return res;
}

//! Update the vector of parameters for the function
void blinn_function::setParameters(const vec& p) 
{
    for(int i=0; i<dimY(); ++i)
    {
        _ks[i] = p[i*2 + 0];
        _N[i]  = p[i*2 + 1];
    }
}

//! Obtain the derivatives of the function with respect to the 
//! parameters. 
vec blinn_function::parametersJacobian(const vec& x) const 
{
    vec jac(dimY()*nbParameters());
	 for(int i=0; i<dimY(); ++i)
		 for(int j=0; j<dimY(); ++j)
		 {
			 if(i == j)
			 {

				 // df / dk_s
				 jac[i*nbParameters() + j*2+0] = std::pow(x[0], _N[j]);

				 // df / dN
				 if(x[0] == 0.0)
					 jac[i*nbParameters() + j*2+1] = 0.0;
				 else
					 jac[i*nbParameters() + j*2+1] = _ks[j] * std::log(x[0]) * std::pow(x[0], _N[j]);
			 }
			 else
			 {
				 jac[i*nbParameters() + j*2+0] = 0.0;
				 jac[i*nbParameters() + j*2+1] = 0.0;
			 }
		 }

    return jac;
}


void blinn_function::bootstrap(const data* d, const arguments& args)
{
    for(int i=0; i<dimY(); ++i)
    {
        _ks[i] = 1.0;
        _N[i]  = 1.0;
    }
}

void blinn_function::save_call(std::ostream& out, 
                               const arguments& args) const
{
    bool is_alta   = !args.is_defined("export") || args["export"] == "alta";

    if(is_alta)
    {
		out << "#FUNC nonlinear_function_blinn" << std::endl ;

		 for(int i=0; i<_nY; ++i)
		 {
			 out << "Ks " << _ks[i] << std::endl;
			 out << "N  " <<  _N[i] << std::endl;
		 }

		 out << std::endl;
	 }
	 else
	 {
		 out << "blinn(L, V, N, X, Y, vec3(";
		 for(int i=0; i<_nY; ++i)
		 {
			 out << _ks[i];
			 if(i < _nY-1) { out << ", "; }
		 }

		 out << "), vec3(";
		 for(int i=0; i<_nY; ++i)
		 {
			 out << _N[i];
			 if(i < _nY-1) { out << ", "; }
		 }

		 out << "))";
	 }

}

void blinn_function::save_body(std::ostream& out, 
                               const arguments& args) const
{
    bool is_shader = args["export"] == "shader" || args["export"] == "explorer";

    if(is_shader)
    {
        out << "vec3 blinn(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, vec3 ks, vec3 Nl)" << std::endl;
        out << "{" << std::endl;
		  out << "\tvec3 H = normalize(L + V);" << std::endl;
        out << "\tvec3 ext_dot = vec3(dot(H,N));" << std::endl;
        out << "\treturn ks * pow(max(ext_dot, vec3(0,0,0)), Nl);" << std::endl;
        out << "}" << std::endl;
    }

}


