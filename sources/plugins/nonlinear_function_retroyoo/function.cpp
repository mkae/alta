#include "function.h"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <core/common.h>

ALTA_DLL_EXPORT function* provide_function()
{
    return new yoo_function();
}

// Overload the function operator
vec yoo_function::operator()(const vec& x) const 
{
	return value(x);
}
vec yoo_function::value(const vec& x) const 
{
	vec res(dimY());

	const double factor = 3.0 / (16.0*M_PI);
	for(int i=0; i<dimY(); ++i)
	{
		const double q = 2.0 * M_PI * x[0] / 633.0;
		res[i] *= _kr[i] * factor * (2.4 + 1.0 / (1.0 + q*_lt[i]*_lt[i])) * (1.0 + (1.0 - exp(- 2.0 * q * 0.7 * _lt[i]) / (q*_lt[i])));
	}
	return res;
}

// Reset the output dimension
void yoo_function::setDimY(int nY)
{
    _nY = nY ;

    // Update the length of the vectors
    _lt.resize(_nY) ;
    _kr.resize(_nY) ;
}

//! Number of parameters to this non-linear function
int yoo_function::nbParameters() const 
{
	return 2*dimY();
}

//! Get the vector of parameters for the function
vec yoo_function::parameters() const 
{
	vec res(3*dimY());
	for(int i=0; i<dimY(); ++i)
	{
		res[i*2 + 0] = _kr[i];
		res[i*2 + 1] = _lt[i];
	}
	return res;
}

//! \brief get the min values for the parameters
vec yoo_function::getParametersMin() const
{
	vec res(2*dimY());
	for(int i=0; i<dimY(); ++i)
	{
		res[i*2 + 0] = 0.0;
		res[i*2 + 1] = 0.0;
	}

	return res;
}

//! Update the vector of parameters for the function
void yoo_function::setParameters(const vec& p) 
{
	for(int i=0; i<dimY(); ++i)
	{
		_kr[i] = p[i*2 + 0];
		_lt[i] = p[i*2 + 1];
	}
}

//! Obtain the derivatives of the function with respect to the 
//! parameters
//! \todo finish. 
vec yoo_function::parametersJacobian(const vec& x) const 
{
	const double factor = 3.0 / (16.0*M_PI);

	vec jac(dimY()*nbParameters());
	for(int i=0; i<dimY(); ++i)
	{
		for(int j=0; j<dimY(); ++j)
		{
			if(i == j)
			{
				const double q = 2.0 * M_PI * x[0] / 633.0;
				const double base = factor * (2.4 + 1.0 / (1.0 + q*_lt[i]*_lt[i])) * (1.0 + (1.0 - exp(- 2.0 * q * 0.7 * _lt[i]) / (q*_lt[i])));
				// df / dk_r
				jac[i*nbParameters() + j*2+0] = base;

				// df / dl_t
				// \todo
				jac[i*nbParameters() + j*2+1] = _kr[i] * base;
			}
			else
			{
				jac[i*nbParameters() + j*2+0] = 0.0;
				jac[i*nbParameters() + j*2+1] = 0.0;
			}
		}
	}

	return jac;
}
		
void yoo_function::bootstrap(const ptr<data>, const arguments&)
{
	for(int i=0; i<dimY(); ++i)
	{
		_kr[i] = 1.0;
		_lt[i] = 1.0;
	}
}

//! Load function specific files
bool yoo_function::load(std::istream& in)
{
	// Parse line until the next comment
	while(in.peek() != '#')
	{
		char line[256];
		in.getline(line, 256);

		// If we cross the end of the file, or the badbit is
		// set, the file cannot be loaded
		if(!in.good())
			return false;
	}

    // Checking for the comment line #FUNC nonlinear_function_lafortune
	std::string token;
	in >> token;
	if(token.compare("#FUNC") != 0) 
	{ 
		std::cerr << "<<ERROR>> parsing the stream. The #FUNC is not the next line defined." << std::endl; 
        return false;
	}

	in >> token;
   if(token.compare("nonlinear_function_retroyoo") != 0) 
	{
		std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl; 
        return false;
	}

	// Parse the lobe
	for(int i=0; i<_nY; ++i)
	{

		in >> token >> _kr[i];
		in >> token >> _lt[i];
	}
    return true;
}


void yoo_function::save_call(std::ostream& out, const arguments& args) const
{
	bool is_alta   = !args.is_defined("export") || args["export"] == "alta";

	if(is_alta)
	{
		out << "#FUNC nonlinear_function_retroyoo" << std::endl ;
		out << "#TYPE ";

		for(int i=0; i<_nY; ++i)
		{
			out << "Kr " << _kr[i] << std::endl;
			out << "Lt " << _lt[i]  << std::endl;
		}

		out << std::endl;
	}
	else
	{
		out << "retroyoo(L, V, N, X, Y, vec3(";
		for(int i=0; i<_nY; ++i)
		{
			out << _kr[i];
			if(i < _nY-1) { out << ", "; }
		}

		out << "), vec3(";
		for(int i=0; i<_nY; ++i)
		{
			out << _lt[i];
			if(i < _nY-1) { out << ", "; }
		}
		out << "))";
	}
}

void yoo_function::save_body(std::ostream& out, const arguments& args) const
{
    bool is_shader = args["export"] == "shader" || args["export"] == "explorer";

    if(is_shader)
    {
        out << "vec3 retroyoo(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, vec3 kr, vec3 lt)" << std::endl;
        out << "{" << std::endl;
        out << "\treturn ks / (4 * " << M_PI << " * a*a * ln*vn) * exp((bn*bn - 1.0) / (a*a*bn*bn)) * g_yoo(L,N,a) * g_yoo(V,N,a);" << std::endl;
        out << "}" << std::endl;
    }
}
