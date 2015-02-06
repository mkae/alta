/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2014, 2015 CNRS
   Copyright (C) 2013, 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include "plugins_manager.h"
#include "rational_function.h"
#include "vertical_segment.h"

#ifdef _WIN32
    #include <windows.h>
#else
    #include <dlfcn.h>
#endif
#include <stdio.h>
#include <list>

//! Add dynamic library extension (.so or .dll) to a dynamic object as well as
//! prepending 'lib' depending on the plateform.
std::string library_name(const std::string name) 
{
	std::string filename;

	const int n = name.length();
#if defined( _WIN32)
	if(name.substr(n-4, n) != std::string("dll")) {
		filename.append(name);
		filename.append(".dll");
	}
#elif defined(__APPLE__)
	if(name.substr(n-5, n) != std::string("dylib")) {
		filename.append("lib");
		filename.append(name);
		filename.append(".dylib");
	}
#else
	if(name.substr(n-2, n) != std::string("so")) {
		filename.append("lib");
		filename.append(name);
		filename.append(".so");
	}
#endif
	else {
		filename.append(name);
	}

	return filename;
}

void get_library_dirs(std::list<std::string>& dirs)
{
	dirs.push_back("");

#ifndef _WIN32
	std::string obj_str ;
	const char *env_str = getenv("ALTA_LIB") ;

	if(env_str == NULL) {
		obj_str = "" ;
	} else {
		obj_str = std::string(env_str)+":" ;
	}

	while(true)	{
		int lb = obj_str.find_first_not_of(":",0) ;
		int le = obj_str.find_first_of(":",lb) ;

		if(lb < le)	{
			dirs.push_back(obj_str.substr(lb,le-lb)) ;
			obj_str = obj_str.substr(le+1,obj_str.size()) ;
		} else {
			break ;
		}
	}
#endif
}

//! \brief Open a dynamic library file (.so or .dll) and extract the associated
//! provide function. The template argument is used to cast the library to a
//! specific type.
template<typename T> T open_library(const std::string& filename, const char* function)
{
	std::list<std::string> basename;
	get_library_dirs(basename);

	std::list<std::string>::iterator iter;
	for(iter = basename.begin(); iter != basename.end(); ++iter)
	{
	std::string libname = *iter;
	libname.append(library_name(filename));


#ifdef _WIN32
    HINSTANCE handle = LoadLibraryA(libname.c_str());
    if(handle != NULL)
    {
        T res = (T)GetProcAddress(handle, function);

        if(res == NULL)
        {
            std::cerr << "<<ERROR>> unable to load the symbol \"" << function << "\" from " << filename << std::endl;
            continue;
        }
#ifdef DEBUG_CORE
        std::cout << "<<DEBUG>> will provide a " << function << " for library \"" << filename << "\"" << std::endl;
#endif
        return res;
    }
    else
    {
        std::cerr << "<<ERROR>> unable to load the dynamic library file \"" << libname << "\"" << std::endl;
        std::cerr << "          cause: \"" << GetLastError() << "\"" << std::endl;
        continue;
    }
#else
    void* handle = dlopen(libname.c_str(), RTLD_GLOBAL | RTLD_LAZY);

	 if(handle != NULL)
	 {
		 void (*res)();
		 *(void **)(&res) = dlsym(handle, function);

        if(dlerror() != NULL)
        {
            std::cerr << "<<ERROR>> unable to load the symbol \"" << function << "\" from " << libname << std::endl;
            continue;
        }
#ifdef DEBUG_CORE
        std::cout << "<<DEBUG>> will provide a " << function << " for library \"" << libname << "\"" << std::endl;
#endif
        return (T)res;
    }
    else
    {
        std::cerr << "<<ERROR>> unable to load the dynamic library file \"" << libname << "\"" << std::endl;
		  std::cerr << "          cause: \"" << dlerror() << "\"" << std::endl;
        continue;
    }
#endif
	}
return NULL;
}

//! \brief load a function from the ALTA input file.
function* plugins_manager::get_function(const std::string& filename)
{
	std::ifstream file;
	file.open(filename.c_str()) ;
	if(!file.is_open())
	{
		std::cerr << "<<ERROR>> unable to open file \"" << filename << "\"" << std::endl ;
		return NULL;
	}

	// Set the precision of the input
	file.precision(10);

	// Parameters of the function object
	int nX, nY;
	params::input param_in   = params::UNKNOWN_INPUT;
	params::output param_out = params::UNKNOWN_OUTPUT;
	arguments args;

	// Test for the first line of the file. Should be a ALTA FUNC HEADER
	std::string line ;
	std::getline(file, line) ;
	if(line != "#ALTA FUNC HEADER")
	{
		std::cerr << "<<ERROR>> this is not a function file" << std::endl;
	}

	// Parse the header for the function command line and the dimension
	// of the function
	while(line != "#ALTA HEADER END")
	{
		std::getline(file, line) ;
		std::stringstream linestream(line) ;

		linestream.ignore(1) ;

		std::string comment ;
		linestream >> comment ;

		if(comment == std::string("DIM"))
		{
			linestream >> nX >> nY ;
		}
		else if(comment == std::string("PARAM_IN"))
		{
			std::string name;
			linestream >> name;
			std::cout << "<<DEBUG>> parsed input parametrization: " << name << std::endl;
			param_in = params::parse_input(name);
		}
		else if(comment == std::string("PARAM_OUT"))
		{
			std::string name;
			linestream >> name;
			param_out = params::parse_output(name);
		}
		else if(comment == std::string("CMD"))
      {
			args = arguments::create_arguments(line.substr(4, std::string::npos));
		}
	}

	// Create the function from the command line
	function* f = get_function(args);
	f->setDimX(nX);
	f->setDimY(nY);
    if(f->input_parametrization() == params::UNKNOWN_INPUT)
    {
        f->setParametrization(param_in);
    }
	f->setParametrization(param_out);

	// Load the function part from the file object
    f->load(file);

	return f;
}

//! Get an instance of the function selected based on the name <em>n</em>.
//! Return NULL if no one exist.
function* plugins_manager::get_function(const arguments& args)
{
    if(!args.is_defined("func"))
    {
#ifdef DEBUG
        std::cout << "<<DEBUG>> no function plugin specified, returning a rational function" << std::endl;
#endif
        return new rational_function();
    }

    // The function to be returned.
    function* func = NULL;

    if(args.is_vec("func"))
    {
        std::vector<std::string> args_vec = args.get_vec("func");

        // Treating the case []
        if(args_vec.size() == 0)
        {
            return NULL;
        }

        //! create a <em>compound</em> class to store multiple
        //! functions in it.
        compound_function* compound = new compound_function();

        //! For each args_vec element, create a function object and add
        //! it to the compound one.
        for(unsigned int i=0; i<args_vec.size(); ++i)
        {
          std::string n("--func ");
          n.append(args_vec[i]);
  #ifdef DEBUG
          std::cout << "<<DEBUG>> load function with args: " << n << std::endl;
  #endif
          arguments temp_args = arguments::create_arguments(n);
          function* f = get_function(temp_args);
          if(dynamic_cast<nonlinear_function*>(f) == NULL)
          {
              std::cerr << "<<ERROR>> only non-linear function care compatible with a compound" << std::endl;
          }
          else
          {
              compound->push_back(dynamic_cast<nonlinear_function*>(f), temp_args);
          }
        }

		  //! return the compound class
		  func = compound;
    }
    else
    {
        const std::string filename = args["func"];
        FunctionPrototype myFunction = open_library<FunctionPrototype>(filename, "provide_function");
        if(myFunction != NULL)
        {
#ifdef DEBUG
            std::cout << "<<DEBUG>> using function provider in file \"" << filename << "\"" << std::endl;
#endif
            func =  myFunction();
        }
        else
        {
            std::cerr << "<<ERROR>> no function provider found in file \"" << filename << "\"" << std::endl;
            return new rational_function() ;
        }
    }

    // Treat the case of the Fresnel
    if(args.is_defined("fresnel"))
    {
      // Cast into a non linear function, only those are permitted to use
      // a Fresnel term.
      nonlinear_function* nl_func = dynamic_cast<nonlinear_function*>(func);
      if(nl_func == NULL)
      {
       std::cerr << "<<ERROR>> only non-linear function are permitted to use Fresnel" << std::endl;
       return func;
      }

      std::cout << "<<DEBUG>> multiplying by a Fresnel term" << std::endl;

      std::string n("--func ");
      if(args.is_vec("fresnel"))
      {
       std::string temp = args["fresnel"];
       n.append(temp.substr(1, temp.size()-2));
      }
      else
      {
       std::string fname = args["fresnel"];
       if(fname.empty()) // Nothing to do except print error, no plugin defined
       {
      	 std::cerr << "<<ERROR>> Fresnel plugin not defined" << std::endl;
      	 std::cerr << "<<ERROR>> using --fresnel alone is not permitted" << std::endl;
      	 return func;
       }
       else // Case where the fresnel parameters is only the plugin filename
       {
      	 n.append(fname);
       }
      }

      nonlinear_function* func_fres = dynamic_cast<nonlinear_function*>(get_function(arguments::create_arguments(n)));
      if(func_fres != NULL)
      {
      return new product_function(nl_func, func_fres);
      }
      else
      {
       std::cerr << "<<ERROR>> the right part of the product is not a nonlinear function. Will use only the left part." << std::endl;
       return func;
      }
	 }

   
/*
	 // Correction of the data by 1/cosine(theta_L)
	 if(args.is_defined("data-correct-cosine"))
	 {
		 nonlinear_function* cosine = new cosine_function();
		 func = new product_function(cosine, dynamic_cast<nonlinear_function*>(func));
	 }
	 // End of correction
*/


    return func;
}
ptr<data> plugins_manager::get_data(const std::string& n)
{
    if(n.empty())
    {
#ifdef DEBUG
        std::cout << "<<DEBUG>> no data plugin specified, returning a vertial_segment loader" << std::endl;
#endif
        return new vertical_segment();
    }

	 DataPrototype myData = open_library<DataPrototype>(n, "provide_data");
    if(myData != NULL)
    {
#ifdef DEBUG
        std::cout << "<<DEBUG>> using function provider in file \"" << n << "\"" << std::endl;
#endif
        return myData();
    }
    else
    {
        std::cerr << "<<ERROR>> no data provider found in file \"" << n << "\"" << std::endl;
        return new vertical_segment() ;
    }
}
ptr<fitter> plugins_manager::get_fitter(const std::string& n)
{
    if(n.empty())
    {
#ifdef DEBUG
        std::cout << "<<DEBUG>> no fitter plugin specified, returning null" << std::endl;
#endif
        return NULL;
    }

    FitterPrototype myFitter = open_library<FitterPrototype>(n, "provide_fitter");
    if(myFitter != NULL)
    {
#ifdef DEBUG
        std::cout << "<<DEBUG>> using function provider in file \"" << n << "\"" << std::endl;
#endif
        return myFitter();
    }
    else
    {
        std::cerr << "<<ERROR>> no fitter provider found in file \"" << n << "\"" << std::endl;
        return NULL ;
    }
}

void plugins_manager::check_compatibility(ptr<data>& d, const ptr<function>& f,
                                          const arguments& args)
{
	if(d->input_parametrization() == params::UNKNOWN_INPUT &&
		f->input_parametrization() == params::UNKNOWN_INPUT)
	{
		std::cout << "<<WARNING>> both function and data objects have no parametrization" << std::endl;
	}
	else
	{
		if(d->input_parametrization() == params::UNKNOWN_INPUT)
		{
			std::cout << "<<WARNING>> unknown parametrization for data" << std::endl;
		}

		if(f->input_parametrization() == params::UNKNOWN_INPUT)
		{
			std::cout << "<<DEBUG>> function will take the parametrization of the data" << std::endl;
			f->setParametrization(d->input_parametrization());
		}
		else if(d->input_parametrization() != f->input_parametrization() && args.is_defined("change-param"))
		{
			std::cout << "<<INFO>> has to change the parametrization of the input data " << params::get_name(d->input_parametrization()) << std::endl;
            std::cout << "<<INFO>> to " << params::get_name(f->input_parametrization()) << std::endl;
			ptr<data_params> dd = new data_params(d, f->input_parametrization());
			d = dynamic_pointer_cast<data>(dd) ;
		}
		else
		{
			std::cout << "<<DEBUG>> no change was made to the parametrization" << std::endl;
		}
	}

	if(f->dimY() != d->dimY())
	{
        f->setDimY(d->dimY());
	}

	/*
	// Check is the data has to be clusterized
	if(args.is_defined("cluster-dim"))
	{
	clustering* cluster = new clustering(d, args);
	d = cluster;
	}
	*/
}

// \todo implement the Darwin (MACOS) version.
#ifdef WIN32
#include <windows.h>
size_t plugins_manager::get_system_memory()
{
	MEMORYSTATUSEX status;
	status.dwLength = sizeof(status);
	GlobalMemoryStatusEx(&status);
	return status.ullTotalPhys;
}
#elif __APPLE__
size_t plugins_manager::get_system_memory()
{
	return 0;
}
#else
#include <unistd.h>
size_t plugins_manager::get_system_memory()
{
	long pages = sysconf(_SC_PHYS_PAGES);
	long page_size = sysconf(_SC_PAGE_SIZE);
	return pages * page_size;
}
#endif
