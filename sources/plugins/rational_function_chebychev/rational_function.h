#pragma once

// Include STL
#include <vector>
#include <string>

// Interface
#include <QObject>
#include <core/function.h>
#include <core/rational_function.h>
#include <core/data.h>
#include <core/fitter.h>
#include <core/args.h>
#include <core/common.h>

class rational_function_chebychev : public rational_function
{
	public: // methods

		rational_function_chebychev() ;
		rational_function_chebychev(const std::vector<double>& a, const std::vector<double>& b) ;
		virtual ~rational_function_chebychev() ;

		// Get the p_i and q_j function
		virtual double p(const vec& x, int i) const ;
		virtual double q(const vec& x, int j) const ;

    protected:  // methods

        //! \brief Save the rational function to the rational format (see \ref formating).
        virtual void save(const std::string& filename) const ;

        //! \brief Output the rational function using a C++ function formating.
        virtual void save_cpp(const std::string& filename, const arguments& args) const ;

        //! \brief Output the rational function using a C++ function formating.
        virtual void save_matlab(const std::string& filename, const arguments& args) const ;
} ;

