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

class rational_function_chebychev : public QObject, public rational_function
{
	Q_OBJECT
	Q_INTERFACES(function)

	public: // methods

		rational_function_chebychev() ;
		rational_function_chebychev(const std::vector<double>& a, const std::vector<double>& b) ;
		virtual ~rational_function_chebychev() ;

		// Get the p_i and q_j function
		virtual double p(const vec& x, int i) const ;
		virtual double q(const vec& x, int j) const ;
} ;

