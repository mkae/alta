#pragma once

#include <vector>
#include <iostream>

class vec : public std::vector<double>
{
	public:

		friend std::ostream& operator<< (std::ostream& out, const vec& v)
		{
			for(unsigned int i=0; i<v.size(); ++i)
			{
				if(i == 0) out << "[" ; else out << ", " ;
				out << v[i] ;
			}
			out << "]" ;

			return out ;
		} ;


} ;
