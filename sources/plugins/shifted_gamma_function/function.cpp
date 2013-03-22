#include "function.h"

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#define one_pi 0.31830988618

// Overload the function operator
vec shifted_gamma_function::operator()(const vec& x) const 
{
	return value(x);
}
vec shifted_gamma_function::value(const vec& x) const 
{
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
}

//! Load function specific files
void shifted_gamma_function::load(const std::string& filename) 
{
}

//! Save the current function to a specific file type
void shifted_gamma_function::save(const std::string& filename, const arguments& args) const 
{
}

//! Number of parameters to this non-linear function
int shifted_gamma_function::nbParameters() const 
{
}

//! Get the vector of parameters for the function
vec shifted_gamma_function::parameters() const 
{
}

//! Update the vector of parameters for the function
void shifted_gamma_function::setParameters(const vec& p) 
{
}

//! Obtain the derivatives of the function with respect to the 
//! parameters. 
vec shifted_gamma_function::parametersJacobian(const vec& x) const 
{
}

		
vec shifted_gamma_function::Fresnel(const vec& F0, const vec& F1, double V_H) const
{
	return F0 - V_H * F1  + (1. - F0)*pow(1. - V_H, 5.);
}

vec shifted_gamma_function::D(const vec& _alpha, const vec& _p, 
                              double cos_h, const vec& _K) const
{
	double cos2 = cos_h*cos_h;
	double tan2 = (1.-cos2)/cos2;
	vec ax = _alpha + tan2/_alpha;
	vec exp_pow(3) ;
	exp_pow[0] = exp(-ax[0]) / pow(ax[0], _p[0]);
	exp_pow[1] = exp(-ax[1]) / pow(ax[1], _p[1]);
	exp_pow[2] = exp(-ax[2]) / pow(ax[2], _p[2]);

	return (one_pi / (cos2 * cos2)) * _K;
}

vec shifted_gamma_function::G1(double theta) const
{
	vec exp_shc(3);
   exp_shc[0] = exp(sh_c[0] * pow(std::max<double>(acos(theta) - sh_theta0[0],0.), sh_k[0]));
   exp_shc[1] = exp(sh_c[1] * pow(std::max<double>(acos(theta) - sh_theta0[1],0.), sh_k[1]));
   exp_shc[2] = exp(sh_c[2] * pow(std::max<double>(acos(theta) - sh_theta0[2],0.), sh_k[2]));
	return 1.0 + sh_lambda * (1.0 - exp_shc);
}

Q_EXPORT_PLUGIN2(shifted_gamma_function, shifted_gamma_function)
