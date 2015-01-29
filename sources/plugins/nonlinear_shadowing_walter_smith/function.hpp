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

    WalterSmith()
    {
      setParametrization(params::CARTESIAN);
      setDimX(6);
    }

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
  

  private: // parameters

  //TODO: CHECK that this is the correct parameter
  vec _alpha;

} ;

