/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2014, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "function.h"

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#define one_pi 0.31830988618

using namespace alta;

ALTA_DLL_EXPORT function* provide_function(const alta::parameters& params)
{
    return new shifted_gamma_function(params);
}

shifted_gamma_function::shifted_gamma_function(const alta::parameters& params)
    : nonlinear_function(params.set_input(6, params::CARTESIAN))
{
    auto nY = params.dimY();

    // Update the length of the vectors
    sh_c      = vec::Zero(nY);
    sh_theta0 = vec::Zero(nY);
    sh_k      = vec::Zero(nY);
    sh_lambda = vec::Zero(nY);
    p         = vec::Zero(nY);
    F_0       = vec::Zero(nY);
    F_1       = vec::Zero(nY);
    K_ap      = vec::Zero(nY);
    rho_d     = vec::Zero(nY);
    rho_s     = vec::Zero(nY);
    alpha     = vec::Zero(nY); 

    alpha.fill(1.0);
}


// Overload the function operator
vec shifted_gamma_function::operator()(const vec& x) const 
{
	return value(x);
}
vec shifted_gamma_function::value(const vec& x) const 
{
	// shading
	vec lv(3); 
	lv[0] = x[0]; lv[1] = x[1]; lv[2] = x[2];

	vec n(3); 
	n[0] = 0.0; n[1] = 0.0; n[2] = 1.0;

	vec ev(3); 
	ev[0] = x[3]; ev[1] = x[4]; ev[2] = x[5];
	
	vec const halfVector = normalize(lv + ev);
	
	double const v_h = dot(ev, halfVector);
	double const n_h = dot(n, halfVector);
	double const n_l = dot(n, lv); 
	
	double inLight = 1.0;
	if (n_l < 0.0)
	{
		inLight = 0.0;
	}

	double const n_v = dot(n, ev);

        // Tried it. It's worse. Back to previous one.
	// return inLight * (n_l * rho_d +  one_pi * rho_s.cwiseProduct(D(alpha, p, n_h, K_ap)).cwiseProduct(G1(n_l)).cwiseProduct(G1 (n_v)).cwiseProduct(Fresnel(F_0, F_1, v_h)));
	return one_pi * inLight * (rho_d + rho_s.cwiseProduct(D(alpha, p, n_h, K_ap)).cwiseProduct(G1(n_l)).cwiseProduct(G1 (n_v)).cwiseProduct(Fresnel(F_0, F_1, v_h)));
}

//! Number of parameters to this non-linear function
int shifted_gamma_function::nbParameters() const 
{
	return 11*_parameters.dimY();
}

//! Get the vector of parameters for the function
vec shifted_gamma_function::parameters() const 
{
	const int n = _parameters.dimY();

	vec res(nbParameters());
	for(int i=0; i<n; ++i) {
		res[i +  0*n] = sh_c[i];
		res[i +  1*n] = sh_theta0[i];
		res[i +  2*n] = sh_k[i];
		res[i +  3*n] = sh_lambda[i];
		res[i +  4*n] = p[i];
		res[i +  5*n] = F_0[i];
		res[i +  6*n] = F_1[i];
		res[i +  7*n] = K_ap[i];
		res[i +  8*n] = rho_d[i];
		res[i +  9*n] = rho_s[i];
		res[i + 10*n] = alpha[i];
	}

	return res;
}

//! Update the vector of parameters for the function
void shifted_gamma_function::setParameters(const vec& pi) 
{
	const int n = _parameters.dimY();
	for(int i=0; i<n; ++i) {
		sh_c[i]      = pi[i +  0*n];
		sh_theta0[i] = pi[i +  1*n];
		sh_k[i]      = pi[i +  2*n];
		sh_lambda[i] = pi[i +  3*n];
		p[i]         = pi[i +  4*n];
		F_0[i]       = pi[i +  5*n];
		F_1[i]       = pi[i +  6*n];
		K_ap[i]      = pi[i +  7*n];
		rho_d[i]     = pi[i +  8*n];
		rho_s[i]     = pi[i +  9*n];
		alpha[i]     = pi[i + 10*n];
	}

}

//! Obtain the derivatives of the function with respect to the 
//! parameters. \TODO
vec shifted_gamma_function::parametersJacobian(const vec& x) const {
	
	const int n = _parameters.dimY();
    vec jac(n*nbParameters());

	 for(int i=0; i<n; ++i) {
		 for(int j=0; j<nbParameters(); ++j) {
		 	jac[i + j*n] = 0.0;
		 }
	 }

    return jac;
}

		
vec shifted_gamma_function::Fresnel(const vec& F0, const vec& F1, double V_H) const
{
	vec F(_parameters.dimY());
	for(int i=0; i<_parameters.dimY(); ++i) {
		F[i] = F0[i] - V_H*F1[i] + (1.0 - F0[i])*pow(1.0 - V_H, 5.0);
	}
	return F;
}

vec shifted_gamma_function::D(const vec& _alpha, const vec& _p, 
                              double cos_h, const vec& _K) const
{
	double cos2 = cos_h*cos_h;
	double tan2 = (1.-cos2)/cos2;

	vec D(_parameters.dimY());
	for(int i=0; i<_parameters.dimY(); ++i) {
		const double ax      = _alpha[i] + tan2/_alpha[i];
		const double exp_pow = exp(-ax) / pow(ax, _p[i]);

		const double P22 = exp_pow;

		D[i] = P22 * (one_pi / (cos2 * cos2)) * _K[i];
	}
	return D;
}

vec shifted_gamma_function::G1(double theta) const
{
	vec G1(_parameters.dimY());
	for(int i=0; i<_parameters.dimY(); ++i)
	{
		const double exp_shc = exp(sh_c[i] * pow(std::max<double>(acos(theta) - sh_theta0[i],0.), sh_k[i]));
		G1[i] =  1.0 + sh_lambda[i] * (1.0 - exp_shc);
	}
	return G1;
}


bool shifted_gamma_function::load(std::istream& in)
{
	// Parse line until the next comment
	while(in.peek() != '#')
	{
		char line[256];
		in.getline(line, 256);

		if(!in.good())
		{
			return false;
		}
	}

	// Checking for the comment line #FUNC nonlinear_function_sgd
	std::string token;
	in >> token;
	if(token.compare("#FUNC") != 0) 
	{ 
		std::cerr << "<<ERROR>> parsing the stream. The #FUNC is not the next line defined." << std::endl; 
		std::cerr << "<<ERROR>> got \"" << token << "\" instead." << std::endl; 
      return false;
	}

	in >> token;
   if(token.compare("nonlinear_function_sgd") != 0) 
	{
		std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl; 
		std::cerr << "<<ERROR>> got \"" << token << "\" instead." << std::endl; 
    return false;
	}

	//Diffuse part first
	for( unsigned int i=0; i < _parameters.dimY(); i++)
	{
		in >> token >> rho_d[i];
	}

	//Specular part 
	for( unsigned int i=0; i < _parameters.dimY(); i++)
	{
		in >> token >> rho_s[i];
	}

	//Alpha
	for( unsigned int i=0; i < _parameters.dimY(); i++)
	{
		in >> token >> alpha[i];
	}

	//Power
	for( unsigned int i=0; i < _parameters.dimY(); i++)
	{
		in >> token >> p[i];
	}

	//F0
	for( unsigned int i=0; i < _parameters.dimY(); i++)
	{
		in >> token >> F_0[i];
	}

	//F1
	for( unsigned int i=0; i < _parameters.dimY(); i++)
	{
		in >> token >> F_1[i];
	}

	//K_alpha_p
	for( unsigned int i=0; i < _parameters.dimY(); i++)
	{
		in >> token >> K_ap[i];
	}

	//Lambda
	for( unsigned int i=0; i < _parameters.dimY(); i++)
	{
		in >> token >> sh_lambda[i];
	}

	// c
	for( unsigned int i=0; i < _parameters.dimY(); i++)
	{
		in >> token >> sh_c[i];
	}

	//k
	for( unsigned int i=0; i < _parameters.dimY(); i++)
	{
		in >> token >> sh_k[i];
	}

	//theta_o
	for( unsigned int i=0; i < _parameters.dimY(); i++)
	{
		in >> token >> sh_theta0[i];
	}

	return true;
}

void 
shifted_gamma_function::save(const std::string& filename, const arguments& args) const
{
	NOT_IMPLEMENTED();
}





