#pragma once

#include <string>
#include <sstream>
#include <map>
#include <cstdlib>
#include <iostream>

#include "common.h"

/*! \brief A useful class for storing the high-level arguments of a program
 *  or a function.
 *  \ingroup core
 *  \internal
 *
 *  The set of parameters are parsed from the command line using the 
 *  constructor. They are stored as std::string in a std::map.
 *  \todo add functionalities to provide new parameters values.
 */
class arguments
{
	public: // functions

		// Constructor and destructor
		arguments()
		{
		}
		arguments(int argc, char** const argv)
        {
            for(int i=0; i<argc; ++i)
            {
                std::string temp(argv[i]) ;
                std::string key, data ;

                if(temp.compare(0, 2, "--") == 0)
                {
//#define DEBUG_ARGS
                    key = temp.substr(2, temp.size()-2) ;
#ifdef DEBUG_ARGS
                    std::cout << "<<DEBUG>> (" << i << ")" << key << " -> [ ";
#endif
                    if(++i < argc)
                    {

                        std::string next(argv[i]) ;
                        if(next[0] == '-')
                        {
                            --i;
                            continue;
                        }

                        if(next[0] == '[' && next[next.size()-1] != ']')
                        {
                            data.append(next);

                            ++i;
                            while(i < argc && next[next.size()-1] != ']')
                            {
                                next = argv[i] ;
                                data.append(" ");
#ifdef DEBUG_ARGS
                                std::cout << " ";
#endif
                                data.append(next);
#ifdef DEBUG_ARGS
                                std::cout << "(" << i << ")" << next;
#endif

                                ++i;
                            }
                            --i;
                        }
                        else
                        {
#ifdef DEBUG_ARGS
                            std::cout << next;
#endif
                            data.append(next);
                        }
                    }
#ifdef DEBUG_ARGS
                    std::cout << "]" << std::endl;
#endif
                }
                _map.insert(std::pair<std::string, std::string>(key, data)) ;
            }
        }
		~arguments()
		{
		}

		//! \brief is the elements in the command line ?
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
        }

        //! \brief is the data at the given key in a vector format?
        //! No matter the type of the data, this function will test is the
        //! mapped string is of type "[ .... ]".
        //! It returns false if there is no associated entry.
        bool is_vec(const std::string& key) const
        {
            if(_map.count(key) > 0)
            {
                return _map.find(key)->second[0] == '[' ;
            }
            else
            {
                return false ;
            }
        }

		//! \brief access the element stored value
		std::string operator[](const std::string& key) const
		{
			if(_map.count(key) > 0)
			{
				return _map.find(key)->second ;
			}
			else
			{
                //std::cerr << "Underfined request to key : \"" << key << "\"" << std::endl ;
				return std::string() ;
			}
        }
		//! \brief acces to the string value associated with the parameter
		//! \a key.
		//!
		//! The \a default_value argument will be returned if the \a key
		//! has no associated value.
		std::string get_string(const std::string& key, std::string default_value = "") const
		{
			if(_map.count(key) > 0)
				return _map.find(key)->second.c_str() ;
			else
				return default_value ;
		}
		//! \brief acces to the float value associated with the parameter
		//! \a key.
		//!
		//! The \a default_value argument will be returned if the \a key
		//! has no associated value.
		float get_float(const std::string& key, float default_value = 0.0f) const
		{
			if(_map.count(key) > 0)
				return (float)atof(_map.find(key)->second.c_str()) ;
			else
				return default_value ;
		}
		//! \brief acces to the integer value associated with the parameter
		//! \a key.
		//!
		//! The \a default_value argument will be returned if the \a key
		//! has no associated value.
		int get_int(const std::string& key, int default_value = 0) const
		{
			if(_map.count(key) > 0)
				return atoi(_map.find(key)->second.c_str()) ;
			else
				return default_value ;
		} 
		//! \brief acces to a vector of float of size \a size associated with
		//! the parameter \a key.
		//!
		//! The \a default_value argument will be returned if the \a key
		//! has no associated value.
		vec get_vec(const std::string& key, int size, float default_value = 0.0f) const
		{
			vec res(size);
            for(int i=0; i<size; ++i)
                res[i] = default_value;

			if(_map.count(key) > 0)
			{
				std::string s = _map.find(key)->second;
				if(s[0] == '[') // Is an array of type [a, b, c]
				{
					int i = 0;
					size_t pos = 1;
					while(pos != std::string::npos && i<size)
					{
						size_t ppos = s.find(',', pos);

						if(ppos != std::string::npos)
						{
                            res[i] = atof(s.substr(pos, ppos-pos).c_str());
							pos = ppos+1;
							++i;
						}
						else
						{
                            std::string temp = s.substr(pos, std::string::npos);
                            res[i] = atof(temp.substr(0, temp.size()-1).c_str());
							pos = ppos;
							++i;
						}
					}
					return res;
				}
			}

			float val = get_float(key, default_value);
			for(int i=0; i<size; ++i)
			{
				res[i] = val;
			}
			return res;
		}

        std::vector<std::string> get_vec(const std::string& key) const
        {
            std::vector<std::string> res;
            if(_map.count(key) > 0)
            {
                std::string s = _map.find(key)->second;

                if(s[0] == '[') // Is an array of type [a, b, c]
                {
                    size_t pos = 1;
                    while(pos != std::string::npos)
                    {
                        size_t ppos = s.find(',', pos);

                        if(ppos != std::string::npos)
                        {
                            std::string temp = s.substr(pos, ppos-pos);
                            res.push_back(temp);
                            pos = ppos+1;
                        }
                        else
                        {
                            std::string temp = s.substr(pos, std::string::npos);
                            temp = temp.substr(0, temp.size()-1);
                            res.push_back(temp);
                            pos = ppos;
                        }
                    }
                }
            }

            return res;
        }

		  //! \brief get the reconstructed command line arguments (without
		  //! the executable name
		  std::string get_cmd() const
		  {
			  std::string cmd;

			  std::map<std::string, std::string>::const_iterator it = _map.begin();
			  for(;it != _map.end(); it++)
			  {
				  if(!it->first.empty())
				  {
					  cmd.append(" --");
					  cmd.append(it->first);
					  cmd.append(" ");
					  cmd.append(it->second);
				  }
			  }

			  return cmd;
		  }

		  //! \brief Create an argument object from a string describing a command
		  //! line.
		  static arguments create_arguments(const std::string& n)
		  {
			  std::vector<std::string> cmd_vec;
			  std::stringstream stream(n);
#ifdef DEBUG_ARGS
			  std::cout << "<<DEBUG>> create argument vector: [";
#endif
			  while(stream.good())
			  {
				  std::string temp;
				  stream >> temp;
#ifdef DEBUG_ARGS
				  std::cout << temp << ", ";
#endif

				  cmd_vec.push_back(temp);
			  }
#ifdef DEBUG_ARGS
			  std::cout << "]" << std::endl;
#endif

			  int argc = cmd_vec.size();
			  char** argv;
			  argv = new char*[argc];
			  for(int i=0; i<argc; ++i)
			  {
				  argv[i] = &cmd_vec[i][0];
			  }

			  arguments current_args(argc, argv);
			  delete argv;
			  return current_args;
		  }


	private: // data

		std::map<std::string, std::string> _map ;

} ;
