#include "function.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>
#include <cassert>

#include <core/common.h>



ALTA_DLL_EXPORT function* provide_function()
{
  return new WalterSmith();
}

WalterSmith::WalterSmith() : SQRT_PI( 1.772453850905516 )
{
  setParametrization(params::CARTESIAN);
  setDimX(6);
}


bool 
WalterSmith::load(std::istream& in) 
{
  assert(0);
  return false;
}

void 
WalterSmith::save_call(std::ostream& out, const arguments& args) const
{
  assert(0);
}

void 
WalterSmith::save_body(std::ostream& out, const arguments& args) const
{
  assert(0);
}


vec 
WalterSmith::value(const vec& x) const
{
  assert(0);
  
  double const LOCAL_EPSILON = 1e-6;


  //Recover light and view directions from x in Cartesian Parametrisation
  vec view(3);
  view[0] = x[0];   view[1] = x[1];   view[2] = x[2];
  
  vec light(3);
  light[0] = x[3];  light[1] = x[4]; light[2] = x[5];

  vec n(3); n[0] = 0.0; n[1] = 0.0; n[2] = 1.0;
  

  // Safe way to cpompute the half vector
  vec half_vector(3);
  half_vector = view + light;
  double const length  = half_vector.norm();

  if( length < LOCAL_EPSILON)
  {
    half_vector = n;    
  }
  else
  {
    half_vector /= length;   
  }
  
  vec res= vec::Zero( dimY() );


  //Now Compute Shadowing Quantities
  double const v_dot_h = dot(view, half_vector); 
  double const v_dot_n = dot(view, n );
  double const v_ratio = v_dot_h / v_dot_n;

  double const light_dot_h = dot(light, half_vector); 
  double const light_dot_n = dot(light, n );
  double const light_ratio = light_dot_h / light_dot_n;

  if (v_ratio <= 0 || light_ratio <=0 )  
  {
    return res;
  }
  else
  {
    for( unsigned int i=0; i < dimY(); i++ )
    {
      double const tan_theta_v = std::tan( std::acos( v_dot_n) );
      double const a_view = 1.0 / (_alpha[i] * tan_theta_v);
      double const geom_term_view = shadowTerm( a_view );

      double const tan_theta_l= std::tan( std::acos( light_dot_n) );
      double const a_light= 1.0 / (_alpha[i] * tan_theta_l);
      double const geom_term_light = shadowTerm( a_light );

      res[i] = geom_term_light * geom_term_view;      
      }
  }
  
  
  
  

}

int  
WalterSmith::nbParameters() const 
{
  return _alpha.size();
} 

vec  
WalterSmith::parameters() const 
{
  vec res;
  assert(0);
  return res;  
}

void 
WalterSmith::setParameters(const vec& p) 
{
  assert(0);
}

vec  
WalterSmith::getParametersMin() const
{
  vec res;
  assert(0);
  return res;

}

vec  
WalterSmith::parametersJacobian(const vec& x) const 
{
  vec res;
  assert(0);
  return res;

}

void 
WalterSmith::bootstrap(const ptr<data> d, const arguments& args)
{
  assert(0);

}

void 
WalterSmith:: setDimY(int nY)
{
  function::setDimY(nY);
  _alpha.resize(nY);
}


double  WalterSmith::shadowTerm( double a ) const 
{
  #define WALTER_SMITH_EXACT
  #ifdef WALTER_SMITH_EXACT
  return 2.0 / ( 1 + erf(a) + std::exp(-a*a) / ( a * SQRT_PI));
  #else //THIS IS THE RATIONAL APPROXIMATION. FASTER BUT LESS ACCURATE
    if( a >= 1.6 )
    {
      return 1.0;
    }
    else
    {
      double const num   = (3.535 + 2.181*a)*a;
      double const denom = (2.267 + 2.577*a)*a + 1.0;

      return num / denom;
    }
  #endif
}


