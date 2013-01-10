#pragma once

/*
 * Fitting interface for generic fitting algorithms
 *
 */
class fitting_algorithm
{
	public:
		// The class function has to be defined for
		// each fitting algorithm. Of course fitting
		// algorithms can share the same function.
		// Also we require that the function class 
		// is a functor.
		//
		class function ;

		// Data class that contains the data to fit.
		// For a rational function it will be intervals
		// but can be anything of any dimension.
		//
		class data ;

		// Static function to fit a data set d with the
		// underling function class. Return the best
		// fit (along with fitting information ?)
		//
		virtual static function fit_data(const data& d) = 0 ;
} ;
