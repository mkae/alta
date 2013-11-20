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
    return new beckmann_function();
}
		
vec beckmann_function::G(const vec& x) const
{
	vec res(dimY());

	for(int i=0; i<dimY(); ++i)
	{
		const double cl = x[2] / (_a[i] * sqrt(1 - x[2]*x[2]));
		const double cv = x[5] / (_a[i] * sqrt(1 - x[5]*x[5]));

		res[i] = 1.0;
		if(cl < 1.6)
		{
			res[i] *= (3.535*cl + 2.181*cl*cl) / (1 + 2.276*cl + 2.577*cl*cl);
		}
		if(cv < 1.6)
		{
			res[i] *= (3.535*cv + 2.181*cv*cv) / (1 + 2.276*cv + 2.577*cv*cv);
		}
	}
	return res;
}

// Overload the function operator
vec beckmann_function::operator()(const vec& x) const 
{
	return value(x);
}
vec beckmann_function::value(const vec& x) const 
{
	double h[3];
	params::convert(&x[0], params::CARTESIAN, params::RUSIN_VH, &h[0]);

	// Compute the Shadow term to init res
	vec res = vec::Zero(dimY()); //G(x);

	for(int i=0; i<dimY(); ++i)
	{
		const double a    = _a[i];
		const double a2   = a*a;
		const double dh2  = h[2]*h[2];
		const double expo = exp((dh2 - 1.0) / (a2 * dh2));

		if(h[2] > 0.0 && x[2]*x[5]>0.0)
		{
			res[i] = _ks[i] / (4.0 /* x[2]*x[5] */* M_PI * a2 * dh2*dh2) * expo;
		}
		else
		{
			res[i] = 0.0; 
		}
	}
	return res;
}

// Reset the output dimension
void beckmann_function::setDimY(int nY)
{
    _nY = nY ;

    // Update the length of the vectors
    _a.resize(_nY) ;
    _ks.resize(_nY) ;
}

//! Number of parameters to this non-linear function
int beckmann_function::nbParameters() const 
{
	return 2*dimY();
}

//! Get the vector of parameters for the function
vec beckmann_function::parameters() const 
{
	vec res(2*dimY());
	for(int i=0; i<dimY(); ++i)
	{
		res[i*2 + 0] = _ks[i];
		res[i*2 + 1] = _a[i];
	}
	return res;
}

//! \brief get the min values for the parameters
vec beckmann_function::getParametersMin() const
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
void beckmann_function::setParameters(const vec& p) 
{
	for(int i=0; i<dimY(); ++i)
	{
		_ks[i] = p[i*2 + 0];
		_a[i]  = p[i*2 + 1];
	}
}

//! Obtain the derivatives of the function with respect to the 
//! parameters
//! \todo finish. 
vec beckmann_function::parametersJacobian(const vec& x) const 
{
	double h[3];
	params::convert(&x[0], params::CARTESIAN, params::RUSIN_VH, h);

	// Get the geometry term
	//vec g = G(x);

    vec jac(dimY()*nbParameters());
	 for(int i=0; i<dimY(); ++i)
	 {
		 for(int j=0; j<dimY(); ++j)
		 {
			 if(i == j && h[2]>0.0 && x[2]*x[5]>0.0)
			 {
				 const double a    = _a[i];
				 const double a2   = a*a;
				 const double dh2  = h[2]*h[2];
				 const double expo = exp((dh2 - 1.0) / (a2 * dh2));
				 const double fac  = (4.0 /* x[2]*x[5] */* M_PI * a2 * dh2*dh2);

				 // df / dk_s
				 jac[i*nbParameters() + j*2+0] = /*g[i] */ expo / fac;

				 // df / da_x
				 jac[i*nbParameters() + j*2+1] = -/* g[i] */ _ks[i] * (expo/(4.0/*x[2]*x[5]*/)) * ((2* a * h[2])/(M_PI*a2*a2*dh2)) * (1 + (dh2 - 1.0)*h[2]/(a2*dh2*h[2]));
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
		
void beckmann_function::bootstrap(const data* d, const arguments& args)
{
	for(int i=0; i<dimY(); ++i)
	{
		_ks[i] = 1.0;
		_a[i]  = 1.0;
	}
}

//! Load function specific files
bool beckmann_function::load(std::istream& in)
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
		std::cerr << "<<ERROR>> got \"" << token << "\" instead." << std::endl; 
      return false;
	}

	in >> token;
   if(token.compare("nonlinear_function_beckmann") != 0) 
	{
		std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl; 
		std::cerr << "<<ERROR>> got \"" << token << "\" instead." << std::endl; 
      return false;
	}

	// Parse the lobe
	for(int i=0; i<_nY; ++i)
	{

		in >> token >> _ks[i];
		in >> token >> _a[i];
	}
    return true;
}


void beckmann_function::save_call(std::ostream& out, const arguments& args) const
{
    bool is_alta   = !args.is_defined("export") || args["export"] == "alta";

    if(is_alta)
    {
		out << "#FUNC nonlinear_function_beckmann" << std::endl ;

		 for(int i=0; i<_nY; ++i)
		 {
			 out << "Ks " << _ks[i] << std::endl;
			 out << "a  " << _a[i]  << std::endl;
		 }
	
		 out << std::endl;
	 }
	 else
	 {
		 out << "beckmann(L, V, N, X, Y, vec3(";
		 for(int i=0; i<_nY; ++i)
		 {
			 out << _ks[i];
			 if(i < _nY-1) { out << ", "; }
		 }

		 out << "), vec3(";
		 for(int i=0; i<_nY; ++i)
		 {
			 out << _a[i];
			 if(i < _nY-1) { out << ", "; }
		 }
		 out << "))";
	 }

}

void beckmann_function::save_body(std::ostream& out, const arguments& args) const
{
    bool is_shader = args["export"] == "shader" || args["export"] == "explorer";

    if(is_shader)
    {
        out << "vec3 g_beckmann(vec3 M, vec3 N, vec3 a)" << std::endl;
		  out << "{" << std::endl;
        out << "\tfloat d = dot(M,N);" << std::endl;
		  out << "\tvec3 c = d / (a * sqrt(1.0f-d));" << std::endl;
		  out << "\tvec3 r;" << std::endl;
		  out << "\tif(c.x < 1.6f) {" << std::endl;
		  out << "\t\tr.x = (3.535*c.x + 2.181*c.x*c.x) / (1 + 2.276*c.x + 2.577*c.x*c.x);" << std::endl;
		  out << "\t} else {" << std::endl;
		  out << "\t\tr.x = 1.0f;" << std::endl;
		  out << "\t}" << std::endl;
		  out << "\tif(c.y < 1.6f) {" << std::endl;
		  out << "\t\tr.y = (3.535*c.y + 2.181*c.y*c.y) / (1 + 2.276*c.y + 2.577*c.y*c.y);" << std::endl;
		  out << "\t} else {" << std::endl;
		  out << "\t\tr.y = 1.0f;" << std::endl;
		  out << "\t}" << std::endl;
		  out << "\tif(c.z < 1.6f) {" << std::endl;
		  out << "\t\tr.z = (3.535*c.z + 2.181*c.z*c.z) / (1 + 2.276*c.z + 2.577*c.z*c.z);" << std::endl;
		  out << "\t} else {" << std::endl;
		  out << "\t\tr.z = 1.0f;" << std::endl;
		  out << "\t}" << std::endl;
		  out << "\treturn r;" << std::endl;
        out << "}" << std::endl;
		  out << std::endl;
        out << "vec3 beckmann(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, vec3 ks, vec3 a)" << std::endl;
        out << "{" << std::endl;
        out << "\tvec3  H   = normalize(L + V);" << std::endl;
        out << "\tfloat hn  = dot(H,N);" << std::endl;
		  out << "\tfloat ln  = dot(L,N);" << std::endl;
		  out << "\tfloat vn  = dot(V,N);" << std::endl;
		  out << "\t" << std::endl;
        out << "\treturn ks / (4 * " << M_PI << " * a*a * ln*vn) * exp((hn*hn - 1.0) / (a*a*hn*hn)) * g_beckmann(L,N,a) * g_beckmann(V,N,a);" << std::endl;
        out << "}" << std::endl;
    }
}
