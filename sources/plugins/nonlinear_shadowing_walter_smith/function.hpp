#pragma once

// Include STL
#include <vector>
#include <string>

// Interface
#include <core/function.h>
#include <core/data.h>
#include <core/fitter.h>
#include <core/args.h>
#include <core/common.h>


/*! 
 *  \ingroup functions
 *   \brief An implementation of the Walter Shadowing Term
 *   \brief This approximation is the correct one and users shoud NOT 
 *   use the one introduced by Schlick, which is available in alta in the nonlinear_shadowing_schlick plugin
 *  
 *
 *  \details
 *  Rational Approximation of the Smith Term that has been introduced
 *  in the paper  "Microfacet Models for Refraction" by Walter et al. EGSR 2007
 *  
 *  The Shadowing Term G1 is approximated by the following rational functions
 *  \f$ \f$
 *  
 *  \author Romain Pacanowski romain.pacanowski@institutoptique.fr
 *  
 *  \todo Finish implementation
 */
class WalterSmith : public nonlinear_function
{

  public: // methods

    WalterSmith();

    //! \brief Load function specific files
    virtual bool load(std::istream& in) ;

    virtual void save_call(std::ostream& out, const arguments& args) const;
    virtual void save_body(std::ostream& out, const arguments& args) const;

  protected: // methods

    
    virtual vec operator()(const vec& x) const { return value(x); }
    virtual vec value(const vec& x) const;

    //! \brief Number of parameters to this non-linear function
    virtual int nbParameters() const ;

    //! \brief Get the vector of parameters for the function
    virtual vec parameters() const ;

    //! \brief Update the vector of parameters for the function
    virtual void setParameters(const vec& p) ;

    //! Get the vector of min parameters for the function
    virtual vec getParametersMin() const;

    //! \brief Obtain the derivatives of the function with respect to the
    //! parameters. 
    virtual vec parametersJacobian(const vec& x) const ;

    //! \brief Boostrap the function by defining the diffuse term
    virtual void bootstrap(const ptr<data> d, const arguments& args);

    //! \brief resize the parameter vector
    virtual void setDimY(int nY);

  private:
    double shadowTerm( double a ) const;
    inline void getVectorsFromCartesianParam( vec const & x,
                                              vec & view,
                                              vec & light, 
                                              vec & normal,
                                              vec & half_vector ) const;

    
    /**
     * Returns 0 if x <= 0 and 1 in other case
     * @param  x the positionat which the step function is to be evaluated
     * @return   value of the step function 
     */
    inline float stepFunc( float x ) const ;

    /**
     * Compute the product of the two step functions  $\chi_1$ and $\chi_2$
     * from the input vector $x$
     * $\chi_1$ is defined as $\chi_1( \frac{v \cdot h}{v \cdot n})$
     * $\chi_2$ is defined as $\chi_2( \frac{l \cdot h}{l \cdot n})$
     *
     * where $v$ is the view direction, $l$ the light direction and $n$
     * is the normal of the surface
     * 
     * @param  x the input vector in the cartesian parametrization
     * @param  view  the view vector
     * @param  light the light vector
     * @param  normal the normal of the surface
     * @param  half_vector the half direction between the light and the view directions
     * @return   the product of the two step function Chi_1 and Chi_2
     */
    inline float partialShadowingTerm( vec const & x,  
                                       vec const & view,
                                       vec const & light, 
                                       vec const & normal,
                                       vec const & half_vector ) const;


    /**
     * Compute the partial derivative of the Shadowing term, expressed
     * as a rational function, according to the parameter alpha and 
     * the given tangent value. 
     *
     * 
     * 
     * @param  alpha the value of the parameter
     * @param  tan_value The tangent value which either tan( theta_view) or tan(theta_light)
     * @return       the value of the partial derivate 
     */
    float partialDerivative( float alpha, float tangent_value) const;

  private: // parameters

  //TODO: CHECK that this is the correct parameter
  vec _alpha;

  double const LOCAL_EPSILON;
  std::string const  FUNCTION_NAME_TOKEN; 
  
  #ifdef WALTER_SMITH_EXACT
  double const SQRT_PI;
  #endif

} ;




inline 
void WalterSmith::getVectorsFromCartesianParam( vec const & x,
                                                vec & view,
                                                vec & light, 
                                                vec & normal,
                                                vec & half_vector ) const
{
  //Demultiplexing the view and light directions from x
  view[0] = x[0];   view[1] = x[1];   view[2] = x[2];
  
  light[0] = x[3];  light[1] = x[4]; light[2] = x[5];

  //Normal Vector
   normal[0] = 0.0; normal[1] = 0.0; normal[2] = 1.0;
  
  // Safe way to compute the half vector
  half_vector = view + light;
  double const length  = half_vector.norm();

  if( length < LOCAL_EPSILON )
  {
    half_vector = normal;    
  }
  else
  {
    half_vector /= length;   
  }
}


inline 
float  WalterSmith::stepFunc( float x ) const 
{
  if ( x <= 0.0 )
  {
    return 0.0f;
  }
  
  return 1.0f;
}


inline float 
WalterSmith::partialShadowingTerm(vec const & x,
                                  vec const & view,
                                  vec const & light, 
                                  vec const & normal,
                                  vec const & half_vector ) const
{

  //Now Compute Step functions 
  double const v_dot_h = dot(view, half_vector); 
  double const v_dot_n = dot(view, normal );

  if( stepFunc( v_dot_h / v_dot_n ) <= 0.0 )
  {
    return 0.0f;
  }

  double const light_dot_h = dot(light, half_vector); 
  double const light_dot_n = dot(light, normal );

  if(stepFunc( light_dot_h / light_dot_n) <= 0.0 )  
  {
    return 0.0f;
  }

  return 1.0f;
}
 

inline 
float
WalterSmith::partialDerivative( float alpha, float tangent_value) const
{
  //Exported from Maple 
  float const t1 = tangent_value;
  float const t2 = alpha * t1;
  float const t4 = alpha * alpha;
  float const t5 = t1 * t1;
  float const t6 = t4 * t5;
  float t13 = 0.1000e4 * t6 + 0.2276e4 * t2 + 0.2577e4;
  t13 *= t13;
  return  -0.1e1 * t1 * (0.4362000e7 * t2 - 0.4145739e7 + 0.3535000e7 * t6) / t13;
}


