#pragma once

template<class X, class Y> class function : public std::function<X(Y)>
{
	public: // methods

		// Overload the function operator
		virtual Y operator()(float X) const = 0 ;
		
		// IO function to text files
		virtual void load(const std::string& filename) = 0 ;
		virtual void save() const = 0 ;

} ;
