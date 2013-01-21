#pragma once

#include <string>
#include <map>
#include <cstdlib>
#include <iostream>

class arguments
{
	public: // functions

		// Constructor and destructor
		arguments(int argc, char** const argv)
		{
			std::string key ;
			std::string data ;
			for(int i=0; i<argc; ++i)
			{
				std::string temp(argv[i]) ;
				std::string data ;

				if(temp.compare(0, 2, "--") == 0)
				{
					key = temp.substr(2, temp.size()-2) ;
					if(i+1 < argc) 
					{
						std::string next(argv[i+1]) ;
						if(next.compare(0, 2, "--") != 0)
						{
							data = next ;
						}
					}
				}
				_map.insert(std::pair<std::string, std::string>(key, data)) ;
			}
		} ;
		~arguments() { } ;

		// Access elements
		bool is_defined(const std::string& key) const
		{
			if(_map.count(key) > 0)
			{
				return true ;
			}
			else
			{
				return false ;
			}
		} ;
		std::string operator[](const std::string& key) const
		{
			if(_map.count(key) > 0)
			{
				return _map.find(key)->second ;
			}
			else
			{
				std::cerr << "Underfined request to key : \"" << key << "\"" << std::endl ;
				return std::string() ;
			}
		} ;
		float get_float(const std::string& key, float default_value = 0.0f)
		{
			std::string value = _map[key] ;
			if(value.empty())
				return default_value ;
			else
				return atof(value.c_str()) ;
		} ;
		int get_int(const std::string& key, int default_value = 0)
		{
			std::string value = _map[key] ;
			if(value.empty())
				return default_value ;
			else
				return atoi(value.c_str()) ;
		} ;

	private: // data

		std::map<std::string, std::string> _map ;

} ;
