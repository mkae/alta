#pragma once

#include <functional>
#include <string>

#include <QtPlugin>

class function : public std::function<double(double)>
{
	public: // methods

		// Overload the function operator
		virtual double operator()(double x) const = 0 ;
		
		// IO function to text files
		virtual void load(const std::string& filename) = 0 ;
		virtual void save() const = 0 ;

} ;

Q_DECLARE_INTERFACE(function, "Fitter.Function") 
