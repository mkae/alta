#pragma once

// Include STL
#include <vector>
#include <string>

// Interface
#include <core/function.h>
#include <core/rational_function.h>
#include <core/data.h>
#include <core/fitter.h>
#include <core/args.h>
#include <core/common.h>

class rational_function_chebychev_1d : public rational_function_1d
{
public: // methods

    rational_function_chebychev_1d() ;
    rational_function_chebychev_1d(int np, int nq) ;
    rational_function_chebychev_1d(const std::vector<double>& a, const std::vector<double>& b) ;
    virtual ~rational_function_chebychev_1d() {}

    // Get the p_i and q_j function
    virtual double p(const vec& x, int i) const ;
    virtual double q(const vec& x, int j) const ;

protected:  // methods

} ;

class rational_function_chebychev : public rational_function
{
	public: // methods

		rational_function_chebychev() ;
		virtual ~rational_function_chebychev() ;

        //! Update the y-1D function for the ith dimension.
        //! \note It will test if the 1D function provided is of the dynamic type
        //! \name rational_function_chebychev_1d
        virtual void update(int i, rational_function_1d* r)
        {
            if(dynamic_cast<rational_function_chebychev_1d*>(r) != NULL)
            {
                rational_function::update(i, r);
            }
            else
            {
#ifdef DEBUG
                std::cerr << "<<ERROR>> the function provided is not of type \"rational_function_chebychev\"" << std::endl;
#endif
            }
        }

    protected:  // methods

        //! \brief Save the rational function to the rational format (see \ref formating).
        virtual void save(const std::string& filename) const ;

        //! \brief Output the rational function using a C++ function formating.
        virtual void save_cpp(const std::string& filename, const arguments& args) const ;

        //! \brief Output the rational function using a C++ function formating.
        virtual void save_matlab(const std::string& filename, const arguments& args) const ;
} ;

