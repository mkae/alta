#pragma once

#include <string>
#include <map>
#include <cstdlib>
#include <iostream>

#include "common.h"

class arguments
{
	public: // functions

		// Constructor and destructor
		arguments()
		{
		} ;
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
		float get_float(const std::string& key, float default_value = 0.0f) const
		{
			if(_map.count(key) > 0)
				return atof(_map.at(key).c_str()) ;
			else
				return default_value ;
		}
		int get_int(const std::string& key, int default_value = 0) const
		{
			if(_map.count(key) > 0)
				return atoi(_map.at(key).c_str()) ;
			else
				return default_value ;
		} 
		vec get_vec(const std::string& key, float default_value) const
		{
			vec res;
			if(_map.count(key) > 0)
			{
				std::string s = _map.at(key);
				if(s[0] == '\[') // Is an array of type [a, b, c]
				{
					int i = 0;
					size_t pos = 0;
					while(pos != std::string::npos)
					{
						size_t ppos = s.find(",", pos);

						if(ppos != std::string::npos)
						{
							std::cout << s.substr(pos, ppos) << std::endl ;
							res[i] = atof(s.substr(pos, ppos).c_str());
							pos = ppos+1;
							++i;
						}
					}
					return res;
				}
			}

			float val = get_float(key, default_value);
			res.push_back(default_value);
			return res;
		}

	private: // data

		std::map<std::string, std::string> _map ;

} ;
