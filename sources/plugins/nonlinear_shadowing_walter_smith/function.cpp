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

WalterSmith::WalterSmith() 
: LOCAL_EPSILON (1e-6),
  FUNCTION_NAME_TOKEN("nonlinear_shadowing_walter_smith")
#ifdef WALTER_SMITH_EXACT
  ,
SQRT_PI( 1.772453850905516 )
#endif
{
  setParametrization(params::CARTESIAN);
  setDimX(6);
}


bool 
WalterSmith::load(std::istream& in) 
{
  using namespace std;

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

  // Checking for the comment line #FUNC nonlinear_function_blinn
  std::string token;
  in >> token;
  if(token != "#FUNC")
  {
      cerr << "<<ERROR>> parsing the stream. The #FUNC is not the next line defined." << endl;
      
      #ifndef DEBUG
      std::cerr << "<<ERROR>> got \"" << token << "\"" << std::endl;
      #endif      
      return false;
  }


  in >> token;
  if(token != FUNCTION_NAME_TOKEN)
  {
      std::cerr << "<<ERROR>> parsing the stream. function name is not the next token." << std::endl;
  #ifndef DEBUG
      std::cerr << "<<ERROR>> got \"" << token << "\"" << std::endl;
  #endif
      return false;
  }

  // alpha [double]
  for(int i=0; i<dimY(); ++i)
  {
      in >> token >> _alpha[i];
  }

  #ifdef DEBUG
  std::cout << "<<DEBUG>> load parameters " << parameters() << std::endl;
  #endif
  return true;


}

void 
WalterSmith::save_call(std::ostream& out, const arguments& args) const
{
  
  bool const is_alta   = !args.is_defined("export") || args["export"] == "alta";

  if(is_alta)
  {
    out << "#FUNC " << FUNCTION_NAME_TOKEN << std::endl ;

    for(int i=0; i< dimY(); ++i)
    {
      out << "alpha " << _alpha[i] << std::endl;
    }
    
    out << std::endl;
  }
  else
  {
   out << "walter_smith(L, V, N, X, Y, vec3(";
   for(int i=0; i<dimY(); ++i)
   {
     out << _alpha[i];
     if(i < dimY()-1) 
      { 
        out << ", "; 
      }
   }

   out << ") ) ";
   
  }

}

void 
WalterSmith::save_body(std::ostream& out, const arguments& args) const
{
  bool is_shader = args["export"] == "shader" || args["export"] == "explorer";

  if(is_shader)
  {
    out << "float shadow_approx_rational(float a )" << std::endl
        << "{" << std::endl
        << "\t if( a >= 1.6 ) " << std::endl
        << "\t {"  << std::endl
        << "\t\t return 1.0; " << std::endl
        << "\t } else { "  << std::endl
        << "\t\t " << std::endl
        << "\t\t float num   = (3.535 + 2.181*a)*a;"  << std::endl
        << "\t\t float denom = (2.267 + 2.577*a)*a + 1.0;" << std::endl
        << "\t\t return num/denom;" << std::endl
        << "\t }"
        << "}" << std::endl << std::endl;

    out << "vec3 walter_smith(vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, vec3 alpha )" << std::endl;
    out << "{" << std::endl;
    out << "\t vec3 H = normalize(L + V);" << std::endl;
    out << "\t float NdotL  = max( dot(N,L), 0.0 );" << std::endl;
    out << "\t float NdotV  = max( dot(N,V), 0.0 );" << std::endl;
    out << "\t float dot_NH = max( dot(N,H), 0.0 );" << std::endl;
    out << "\t float chi_v =  step( 0.0, dot_NH / NdotV ) ; " << std::endl;
    out << "\t if( chi_v <= 0.0 ) { return vec3(0.0, 0.0, 0.0); }" << std::endl;
    out << "\t float chi_l =  step( 0.0, dot_NH / NdotL ) ; " << std::endl;
    out << "\t if( chi_l <= 0.0 ) { return vec3(0.0, 0.0, 0.0); }" << std::endl
        << std::endl;

    out << "\t float tan_theta_v = tan( acos( NdotV ) ); " << std::endl
        << std::endl
        << "\t vec3 geom_term_view(0.0); " << std::endl
        << "\t geom_term_view[0] = shadow_approx_rational( 1.0 / (alpha[0]*tan_theta_v)) ;"<< std::endl
        << "\t geom_term_view[1] = shadow_approx_rational( 1.0 / (alpha[1]*tan_theta_v)) ;"<< std::endl
        << "\t geom_term_view[2] = shadow_approx_rational( 1.0 / (alpha[2]*tan_theta_v)) ;"<< std::endl
        << std::endl
        << "\t float tan_theta_l = tan( acos( NdotL ) ); " << std::endl
        << "\t vec3 geom_term_light(0.0); " << std::endl    
        << "\t geom_term_light[0] = shadow_approx_rational( 1.0 / (alpha[0]*tan_theta_l)) ;" << std::endl
        << "\t geom_term_light[1] = shadow_approx_rational( 1.0 / (alpha[1]*tan_theta_l)) ;" << std::endl
        << "\t geom_term_light[2] = shadow_approx_rational( 1.0 / (alpha[2]*tan_theta_l)) ;" << std::endl
        << "\t return geom_term_view * geom_term_light;" << std::endl
        << "}" << std::endl;

  }


}


vec 
WalterSmith::value(const vec& x) const
{   
  
   //Recover light and view directions from x in Cartesian Parametrisation
  vec view(3);   vec light(3);  
  //Normal Vector
  vec n(3); vec half_vector(3);

  getVectorsFromCartesianParam( x, view, light, n, half_vector);

  //Initialize result to zero by default
  vec res= vec::Zero( dimY() );

  if( partialShadowingTerm( x, view, light, n, half_vector) <= 0.0f )
  {
    return res;
  }
  else
  {
    //Now Compute Shadowing Quantities
    for(unsigned int i=0; i < dimY(); i++)
    {
      //double const tan_theta_v = std::tan( std::acos( v_dot_n) );
      double const tan_theta_v = std::sqrt(1.0 - view[2]*view[2]) / view[2];
      double const a_view = 1.0 / (_alpha[i] * tan_theta_v);
      double const geom_term_view = shadowTerm( a_view );

      //double const tan_theta_l= std::tan( std::acos( light_dot_n) );
      double const tan_theta_l = std::sqrt(1 - light[2]*light[2]) / light[2];
      double const a_light= 1.0 / (_alpha[i] * tan_theta_l);
      double const geom_term_light = shadowTerm( a_light );

      res[i] = geom_term_light * geom_term_view;      
    }
  }

  return res;    
}

int  
WalterSmith::nbParameters() const 
{
  return _alpha.size();
} 

vec  
WalterSmith::parameters() const 
{
  vec res( _alpha.size() );
  for(unsigned int i=0; i< _alpha.size(); i++)
  {
      res[i] = _alpha[i];
  }

  return res;
}

void 
WalterSmith::setParameters(const vec& p) 
{
  for(unsigned int i=0; i<_alpha.size(); i++)
  {
    _alpha[i] = p[i];
  }
}

vec  
WalterSmith::getParametersMin() const
{
  vec res( _alpha.size() );
 
  for( unsigned int i=0; i < _alpha.size(); i++)
  {
    res[i] = LOCAL_EPSILON;
  }

  return res;
}

vec  
WalterSmith::parametersJacobian(const vec& x) const 
{
  //A 3x3 matrix for the Jacobian
  vec jac = vec::Zero(dimY()*nbParameters());

  vec view(3);   vec light(3);  
  //Normal Vector
  vec n(3); 
  vec half_vector(3);

  getVectorsFromCartesianParam( x, view, light, n, half_vector);

  //First step compute Step functions to check that they are not equal to zero
  if( partialShadowingTerm( x, view, light, n, half_vector) <= 0.0f )
  {
    return jac;
  }


  double const tan_theta_v = std::sqrt(1.0 - view[2]*view[2]) / view[2];
  
  double const tan_theta_l = std::sqrt(1 - light[2]*light[2]) / light[2];
       
  

  // if they are not compute derivative of the alpha parameters
  for(int i=0; i<dimY(); ++i)
  {
   for(int j=0; j<dimY(); ++j)
   {
      if(i == j)
      {
          double const a_theta_v = 1.0 / (_alpha[i] * tan_theta_v);
          double const a_theta_l = 1.0 / (_alpha[i] * tan_theta_l);

          jac[i*nbParameters() + j] =  partialDerivative(_alpha[i], tan_theta_v ) * shadowTerm( a_theta_l ) + shadowTerm( a_theta_v ) * partialDerivative(_alpha[i], tan_theta_l ) ;

      }

    }//end for j loop
  }//end for i loop

  return jac;
}


void 
WalterSmith::bootstrap(const ptr<data> d, const arguments& args)
{
  
  for(unsigned int i=0; i<_alpha.size(); i++) 
  { 
    _alpha[i] = 1.0; 
  }
}

void 
WalterSmith:: setDimY(int nY)
{
  function::setDimY(nY);
  _alpha.resize(nY);
}


double  WalterSmith::shadowTerm( double a ) const 
{
  //#define WALTER_SMITH_EXACT
  #ifdef WALTER_SMITH_EXACT
  return 2.0 / ( 1 + erf(a) + std::exp(-a*a) / ( a * SQRT_PI) );
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


